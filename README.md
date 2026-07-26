# biq

A recording format for SDR IQ that is 2–4× smaller than the 16-bit WAV files
everyone uses, and still seeks to any sample in constant time.

`.biq` is a thin container around [block floating
point](https://github.com/ibelinp/bfp-iq-codec): a 64-byte header, an optional JSON
metadata blob, then a stream of fixed-size BFP blocks. That's the whole format.
[**SPEC.md**](SPEC.md) is the normative description and fits in one sitting.

## Why

Every IQ recording format in common use stores raw samples at a **fixed scale**.
`cu8` off an RTL dongle is 2 bytes a sample. `ci16` — WAV with a SpectraVue
`auxi` chunk, RF64, SigMF `ci16_le`, Rohde & Schwarz `.wv`, headerless `.cs16` —
is 4. GNU Radio's `cf32` is 8. They differ in depth, and in how they spell the
sample rate and the centre frequency, and in nothing else that matters: the
scale is the same whatever the signal is doing. A 24-hour capture at 5 Msps is
864 GB, 1.7 TB, or 3.5 TB depending which one you picked.

The obvious fix is to compress, and the obvious compressors are all wrong for
this. FLAC and gzip are variable-rate and stateful, so you can't seek without
an index and you can't predict the size. Dropping to 8 bits works but wastes
range: a scale sized for the loud case quantises the quiet case near the
bottom. Server-side AGC destroys the levels you wanted to record in the first
place.

Block floating point sidesteps all three. A block of N samples shares one
exponent and stores a `b`-bit mantissa per component, so the scale follows the
signal — and the exponent is *bookkeeping the decoder undoes*, not a gain
decision. At 6 bits it is 2.7× smaller than `ci16`; that same 24-hour capture is
650 GB.

## What you get for keeping the block size fixed

Holding `N` and `b` constant for a whole file is the design decision the rest of
the format hangs on. It costs nothing (change parameters at a segment boundary
instead) and it buys four things no compressed format usually has:

**Constant-time seek.** No index, no scan — the offset of any sample is
arithmetic:

```
offset(s) = data_start + (s / N) * block_stride
```

The file is memory-mappable and parallel-decodable, and any block decodes with
no reference to any other.

**A free amplitude envelope.** Every block carries its exponent at a fixed
stride, so reading *only the exponent bytes* — one strided pass, `filesize /
block_stride` bytes — gives you a log-domain envelope of the entire recording at
N-sample resolution without decoding a single mantissa. Roughly 20 kB per second
of recording. That's your scrub bar and your "where's the burst in this 40 GB
file", for free, with no sidecar index.

**Truncation that just works.** `data_bytes = 0` legally means "runs to EOF", so
a recorder killed by a crash or a full disk leaves a fully readable file minus
at most the final partial block. A WAV in the same situation has a bogus
`dataSize` and needs repair tooling.

**Time that survives damage.** Because the stride is fixed, a sample's index is
positional — so corruption in the middle of a file cannot desync the index, and
samples *after* the damage are still at the right timestamp. The time reference
itself sits at byte 0 and is written before the first sample, so truncation
never reaches it. And every segment carries its own clock reading rather than
counting forward from the first one ([SPEC §10.1](SPEC.md)), so a lost or
corrupt segment leaves the rest of the series correctly timed and error cannot
accumulate across a recording.

The failure this *doesn't* absorb is samples the source dropped without the
recorder noticing — the file then looks continuous while everything after the
gap runs late, which is worse than an obvious failure. The spec's answer is to
roll a new segment at the gap so the drop becomes an ordinary boundary with a
fresh, correct-by-construction timestamp, rather than a note in the metadata
that a dying recorder might never write.

## Size

Bytes per complex sample is `b/4 + 1/N`. At N = 256:

| profile | b | bytes/sample | vs. int16 | 25 s @ 5 Msps | 24 h @ 5 Msps | intended use |
|---|---|---|---|---|---|---|
| `survey` | 4 | 1.004 | 4.0× | 125 MB | 434 GB | wideband monitoring, long unattended captures |
| `general` | 6 | 1.504 | 2.7× | 188 MB | 650 GB | **default.** Weak-signal, CW, digital modes |
| `archive` | 8 | 2.004 | 2.0× | 250 MB | 866 GB | keep-forever captures |
| — | int16 | 4.0 | 1.0× | 500 MB | 1.7 TB | what you're writing today |

One mantissa width governs both I and Q — the two components of a complex
baseband sample have identical statistics, so splitting the budget between them
buys nothing. Peak SQNR is the textbook `6.02·b + 1.76 dB`, measured against the
block's own peak rather than a fixed system full-scale, which is exactly why a
weak signal keeps its precision.

## Metadata

The 64-byte header carries what you need to decode and to know what you're
looking at: sample rate and centre frequency as `f64` (the SpectraVue `auxi`
chunk's `u32` Hz dies above 4.29 GHz), UTC start time in nanoseconds, and
`ref_dbm_full_scale` so playback shows the same S-meter reading the live
receiver did. That last one is the field none of the existing formats have.

Everything else goes in the JSON blob, which **reuses SigMF's key names
verbatim** (`core:sample_rate`, `core:frequency`, `core:datetime`,
`core:geolocation`, …). Converting to and from a `.sigmf-meta` sidecar is
mechanical, and nobody has to learn a second metadata vocabulary. Retunes and
dropouts are sample-indexed entries in `biq:events`, never inline in the block
stream — keeping the stream pure is what preserves the fixed stride.

## Try it

```sh
cd reference && make        # builds biqconv, no dependencies
./biqconv capture.wav capture.biq -b 6
./biqconv info capture.biq
./biqconv capture.biq back.wav
```

Measured on a 25-second 5 Msps FM-band capture from an ADALM-Pluto:

```
capture.wav -> capture.biq
  125009920 samples  100 MHz  5 Msps  b=6 N=256
  500.04 MB -> 188.00 MB   2.66x smaller
```

Round-tripped back to 16-bit WAV, that file measures 29.2 dB SNR against the
original with a worst-case single-sample error of 4 LSB out of 32768 — which is
what 6-bit mantissas against this signal's crest factor predicts.

`biqconv info` reads the whole recording's amplitude envelope from the exponent
bytes alone. On the 188 MB file that is 488,320 bytes touched and no mantissas
decoded.

## Sample recordings

No receiver, or a WAV `biqconv` won't read? [`samples/`](samples/) holds the
same 5-second capture at all three profiles, so you can start from real `.biq`
files and convert nothing:

```sh
cd reference && make
./biqconv info ../samples/fm-band-95.5MHz-2Msps-b6-general.biq
./biqconv ../samples/fm-band-95.5MHz-2Msps-b6-general.biq back.wav
```

| file | profile | b | size | vs. int16 |
|---|---|---|---|---|
| [`…-b4-survey.biq`](samples/fm-band-95.5MHz-2Msps-b4-survey.biq) | `survey` | 4 | 10.00 MB | 3.98× |
| [`…-b6-general.biq`](samples/fm-band-95.5MHz-2Msps-b6-general.biq) | `general` | 6 | 14.98 MB | 2.66× |
| [`…-b8-archive.biq`](samples/fm-band-95.5MHz-2Msps-b8-archive.biq) | `archive` | 8 | 19.96 MB | 2.00× |

FM broadcast band off an ADALM-Pluto: 95.495516 MHz centre, 2 Msps complex,
9,958,400 samples in 38,900 blocks of 256, recorded 2026-07-26T19:14:42Z. The
span holds a dozen stations, and RDS decodes from all three — including
`survey`, where four mantissa bits still carry a 57 kHz subcarrier well enough
for a real decoder to lock. That is the claim the profile table makes, in a
file you can check yourself.

The [`reference/golden/`](reference/golden/) vectors pin down the bit packing;
these pin down everything around it — the seek formula, the exponent envelope,
a header written by a real encoder — against a file you can also listen to.
Being `.biq` already, none of the usual IQ-WAV variation (8-bit `cu8`, `cf32`,
mono, RF64, whichever chunk order your recorder felt like) sits between you and
a first decode.

## Reference implementation

[`reference/`](reference/) is dependency-free and self-testing.

| file | what it is |
|---|---|
| [`biq.hpp`](reference/biq.hpp) | header-only C++11 encoder, decoder, `Writer`, `Reader`. Drop it in. |
| [`biq.py`](reference/biq.py) | the spec in Python — literal, slow, good for notebooks |
| [`biqconv.cpp`](reference/biqconv.cpp) | the CLI converter |
| [`biq_test.cpp`](reference/biq_test.cpp) | 1369 checks: packing, exponents, headers, seek, truncation, envelope |
| [`golden/`](reference/golden/) | test vectors from an independent Rust implementation |

```sh
make test     # C++ and Python, both against the golden vectors
```

The golden vectors are the part worth caring about if you write your own port.
They were produced by a separate Rust encoder, and both reference
implementations reproduce them **byte for byte** at b = 2, 4, 5, 6, 8, and 16 —
covering byte-aligned widths, widths that straddle two bytes, and widths that
straddle three. Interop on this format means agreeing on LSB-first bit packing
and on which exponent the encoder picks; the vectors pin down both. If your
implementation matches `golden/`, it will read anyone's files.

(On the exponent: `ceil(log2(peak / m_max))` is the rule, but a libm that rounds
`log2` the other way at an exact power of two picks `e±1`. Both are legal —
decoders never depend on the encoder's choice — so all three implementations
settle it with exact float compares instead, and agree everywhere.)

## Status

Draft. Implementations work and agree; the spec is not frozen until it has been
used in anger.

- [x] `SPEC.md` — format definition
- [x] `reference/biq.hpp` — header-only C++ encoder/decoder
- [x] `reference/biq.py` — Python reference, matching byte for byte
- [x] `biqconv` — CLI converter: `wav ↔ biq`, `sigmf ↔ biq`, `cs16 ↔ biq`
- [ ] multi-channel (the `channels` byte is reserved, MUST be 1 today)
- [ ] a recorder writing `.biq` live rather than converting after the fact

## Relationship to bfp-iq-codec

[bfp-iq-codec](https://github.com/ibelinp/bfp-iq-codec) explains block floating
point itself — the technique, the exponent rule, how to choose `b`, and how it
compares to block scaling and to entropy coders. It is the reading for *why*.

This repo is the *file*: a container around blocks that are byte-identical to the
ones bfp-iq-codec defines for streaming. One codec, two transports. If you have
already implemented bfp-iq-codec, adding `.biq` is a header parser and a seek
formula.

## License

MIT. See [LICENSE](LICENSE). Use it, ship it, port it.

There is nothing here to own. Block floating point is decades-old textbook DSP,
and this is sixty-four bytes of header wrapped around it. The spec is
deliberately small enough to implement in an afternoon, the metadata keys are
SigMF's rather than mine because the world does not need a second vocabulary for
"sample rate", and the golden vectors are checked in so nobody has to take my
word for what the bytes should be. If the container doesn't suit you, take the
ideas out of it — the fixed stride, the exponent envelope, the header that
survives a truncated file — and put them somewhere better.

I wrote it down because IQ is the raw material of everything interesting in
radio, and we throw it away constantly. Captures deleted to free a disk. The
good five seconds kept and all the context around them lost. Band surveys never
run because a week of spectrum is a wall of terabytes. Being four times smaller
doesn't change what's possible in principle — it changes what people actually
keep. If this means one more recording survives, or somebody hands over a whole
capture instead of a screenshot of it, it did its job.
