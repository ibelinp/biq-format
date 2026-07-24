# .biq — Block Floating Point IQ recording format

**Version 1 (magic `BIQ1`).** Status: draft, not yet frozen.

A container for recording SDR IQ as block floating point. It is a 64-byte
binary header, an optional JSON metadata blob, and then a stream of BFP blocks
— the same blocks defined in [bfp-iq-codec](https://github.com/ibelinp/bfp-iq-codec), unchanged.

The design goal is a file that is 2–4× smaller than the int16 recordings
everyone uses today while keeping the two properties those recordings have and
compressed formats usually lose: **constant-time seek to any sample**, and
**decode with no state and no lookahead**.

All multi-byte fields are little-endian. Floats are IEEE 754 (`f64` =
binary64). Signed integers are two's complement.

---

## 1. File layout

```
offset                        contents
0                             header          (64 bytes, §2)
header_bytes                  metadata JSON   (meta_bytes, §3, may be 0)
header_bytes + meta_bytes     block stream    (§4)
```

`header_bytes` and `meta_bytes` come from the header. Readers MUST locate the
metadata and the block stream using those fields, never by assuming 64.

## 2. Header (64 bytes)

| off | size | type | field | notes |
|---|---|---|---|---|
| 0 | 4 | char[4] | `magic` | `BIQ1` (0x42 0x49 0x51 0x31) |
| 4 | 2 | u16 | `header_bytes` | 64 in this version |
| 6 | 1 | u8 | `flags` | bit 0 = real-valued input (§4.1); bits 1–7 reserved, MUST be 0 |
| 7 | 1 | u8 | `channels` | **reserved — MUST be 1 (§2.1)** |
| 8 | 1 | u8 | `mantissa_bits` | `b`, 2…16, applies to both I and Q (§5) |
| 9 | 1 | u8 | `fill_rule` | 0 = `ceil`, 1 = tuned fill. Provenance only; decoders ignore it |
| 10 | 2 | u16 | `block_size` | `N`, samples per block. MUST be a multiple of 8 |
| 12 | 4 | u32 | `block_stride` | bytes per block on disk, ≥ `block_bytes` (§4.2) |
| 16 | 8 | f64 | `sample_rate_hz` | complex samples/s. Fractional rates permitted |
| 24 | 8 | f64 | `center_freq_hz` | tuned centre of the recording |
| 32 | 8 | i64 | `start_time_unix_ns` | UTC nanoseconds of this file's first sample. 0 = unknown. Per segment, from the clock — see §10.1 |
| 40 | 8 | f64 | `ref_dbm_full_scale` | dBm of a full-scale sine. NaN = uncalibrated (§6) |
| 48 | 8 | u64 | `data_bytes` | length of the block stream. **0 = read to EOF (§7)** |
| 56 | 4 | u32 | `meta_bytes` | length of the JSON blob. MUST be a multiple of 8 |
| 60 | 4 | — | reserved | MUST be 0 |

A reader MUST reject the file if `magic` is wrong, `mantissa_bits` is outside
2…16, `block_size` is 0 or not a multiple of 8, `block_stride` is less than the
computed `block_bytes`, or any reserved field is non-zero.

### 2.1 The `channels` byte

The field exists so that coherent multi-channel captures can be added without a
format break. **In version 1 it MUST be written as 1, and readers MUST reject
any other value.** Rejecting is deliberate: a v1 reader that silently treated a
4-channel file as single-channel would produce plausible-looking garbage.

Non-normative, so the field's meaning is unambiguous when it is unreserved: the
intended layout is block-major, channel-minor — for block index `k`, the `C`
channels appear consecutively, each an independent BFP block with its own
exponent. Independent exponents are safe for phase-coherent work because decode
restores absolute scale exactly, up to quantisation.

## 3. Metadata JSON

`meta_bytes` of UTF-8 immediately after the header. It is a single JSON object,
padded to a multiple of 8 bytes with spaces (0x20) — trailing whitespace is
legal JSON, so a parser can be handed the slice as-is. `meta_bytes = 0` means
no metadata, which is legal; everything needed to decode is in the header.

**Keys reuse SigMF's vocabulary verbatim** wherever SigMF has one, so
translation to and from a `.sigmf-meta` sidecar is mechanical:

```json
{
  "core:sample_rate": 5000000.0,
  "core:frequency": 100000000.0,
  "core:datetime": "2026-07-17T21:10:08.000000Z",
  "core:author": "Pieter Ibelings",
  "core:hw": "ADALM-Pluto",
  "core:recorder": "SpectraWeb 0.2.0",
  "core:geolocation": { "type": "Point", "coordinates": [-80.1, 26.1, 3.0] },

  "biq:events": [
    { "sample": 12500000, "center_freq_hz": 101100000.0 },
    { "sample": 41000000, "gap_samples": 4096 }
  ]
}
```

Values duplicated between the header and the JSON (sample rate, frequency,
start time) MUST agree; where they disagree, **the binary header wins**. The
header is the machine-readable truth and the JSON is for humans and for fields
that have no header slot.

Retunes, dropouts and tags go in `biq:events` as sample-indexed entries, not
inline in the block stream. Keeping the block stream pure is what preserves the
fixed stride, and the fixed stride is the whole point.

Unknown keys MUST be ignored, not rejected.

## 4. Block stream

The block stream is a concatenation of BFP blocks, byte-identical to the
`IQ_BLOCKS` body of the bfp-iq-codec streaming format. One block is:

```
offset            size                field
0                 1 byte (i8)         exponent e
1                 ceil(N·b/8) bytes   I mantissas, N of them, b bits each
1 + plane_bytes   ceil(N·b/8) bytes   Q mantissas, N of them, b bits each
```

```
plane_bytes = ceil(N · b / 8)          # exact, since N is a multiple of 8
block_bytes = 1 + 2 · plane_bytes
```

- Mantissas are signed two's-complement, `b` bits wide, **packed LSB-first**:
  sample `i` occupies bits `[i·b, (i+1)·b)` of its plane, low bits first.
- Planes are **planar, not interleaved** — all I, then all Q.
- The exponent is a signed byte; `-128` is the sentinel for an all-zero block.
- Decode is `sample = mantissa · 2^exponent`. One multiply per component, or a
  shift in fixed point. No cross-block state.

Encoding, the exponent rule, and the choice of `b` are specified in the
[bfp-iq-codec](https://github.com/ibelinp/bfp-iq-codec) README and are not restated here.

### 4.1 Real-valued input

When `flags` bit 0 is set the source is real, not complex: a block is the
exponent byte plus **one** plane, and `block_bytes = 1 + plane_bytes`. `N` is
then real samples per block. This exists for direct-sampling HF front ends.

### 4.2 `block_stride`

`block_stride` is the on-disk spacing of blocks and MAY exceed `block_bytes`.
Padding bytes MUST be written as zero and MUST be ignored on read.

Encoders SHOULD set `block_stride = block_bytes` — the format exists to save
bytes. An encoder targeting mmap or SIMD readers MAY round up to a multiple of
8 (N=256, b=5 → 321 becomes 328, costing 2.2%). Decoders MUST honour whatever
is in the header.

## 5. Mantissa width, block size, and profiles

One `mantissa_bits` governs **both** I and Q. Asymmetric depths are not
permitted: the two components of a complex baseband sample are statistically
identical, so spending different precision on them buys nothing and costs the
packing its uniform stride.

`block_size` MUST be a multiple of 8, which makes `N·b` a whole number of bytes
for every legal `b` — no padding bits, no partial trailing byte, ever.

Three profiles are named so tools can agree on defaults. They are conventions,
not restrictions: any legal `b` and `N` is a conforming file.

| profile | b | N | bytes/sample | vs. int16 | intended use |
|---|---|---|---|---|---|
| `survey` | 4 | 256 | 1.004 | 4.0× | wideband monitoring, long unattended captures |
| `general` | 6 | 256 | 1.504 | 2.7× | **default.** Weak-signal, CW, digital modes |
| `archive` | 8 | 256 | 2.004 | 2.0× | keep-forever captures; ~50 dB display floor |

`b ∈ {4, 8, 16}` unpack with byte and nibble operations alone. Depths 5, 6, 7
need cross-byte shifting for roughly 15% more saving. A decoder MAY special-case
the aligned depths for speed but MUST support all of 2…16.

N = 256 is the default because past ~256 the exponent byte is already amortised
to nothing, while larger blocks track level worse: one transient inside a long
block lifts the exponent for every sample around it and crushes the quiet parts.

## 6. Levels

`ref_dbm_full_scale` is the absolute power of a full-scale sine at the
recording point, so playback can show the same S-meter reading the live receiver
showed. Decoded samples are nominally within ±1.0; a value of NaN means the
recording is uncalibrated and levels are relative only.

Nothing else in this format applies gain, removes DC, or otherwise alters the
samples. A `.biq` file is the receiver's output quantised, and nothing more.

## 7. Seeking, truncation, and the exponent envelope

**Seek** to complex sample `s` is arithmetic, with no index and no scan:

```
data_start = header_bytes + meta_bytes
offset(s)  = data_start + (s / N) * block_stride        # integer division
```

The file is therefore memory-mappable and trivially parallel-decodable, and any
block decodes without reference to any other.

**Truncation** is recoverable by construction. `data_bytes = 0` means "the
stream runs to EOF", so a recorder killed by a crash, a full disk, or a pulled
cable leaves a file that is fully readable — at most the final partial block is
lost. Writers MAY leave `data_bytes` at 0 permanently; writers that know the
length SHOULD backfill it at finalise, and readers MUST accept both.

**The exponent envelope** falls out for free. Reading only the byte at each
`offset(k·N)` — one strided pass over `filesize / block_stride` bytes — yields a
log-domain amplitude envelope of the entire recording at N-sample resolution,
without decoding a single mantissa. At N=256 and 5 Msps that is ~19.5k points
per second of recording, about 20 kB per second of audio-rate scrub data. It is
enough to draw a full-length overview bar or find the burst in a 40 GB capture
in one pass. No sidecar index is defined because none is needed.

## 8. Conformance

A conforming **decoder** MUST:

1. Validate `magic`, and reject `channels != 1` and any non-zero reserved field.
2. Locate metadata and data via `header_bytes` and `meta_bytes`.
3. Support `b` = 2…16 and any `N` that is a multiple of 8.
4. Honour `block_stride`, treating padding as opaque.
5. Treat `data_bytes = 0` as read-to-EOF, and tolerate a truncated final block.
6. Ignore unknown JSON keys, and prefer the header where the two disagree.

A conforming **encoder** MUST write a header that satisfies §2, MUST NOT vary
`b`, `N`, or `block_stride` within a file, and SHOULD start a new file when a
parameter needs to change. An encoder that writes a segment series MUST stamp
each segment's `start_time_unix_ns` from the clock rather than deriving it
(§10.1), and one that continues through a known sample gap MUST record it
(§10.2).

The exponent an encoder picks is **not** part of conformance: any exponent that
does not clip the block is legal, and decoders never depend on the choice. The
`ceil` rule in §4 is the recommended one. Note that computing it as
`ceil(log2(peak / m_max))` leaves the answer at the mercy of libm — a platform
that rounds `log2` the other way at an exact power of two picks `e±1` — so
implementations that want to reproduce the test vectors byte for byte should
settle it with exact float compares instead (see `reference/biq.hpp`).

### 8.1 Test vectors

`reference/golden/` contains `golden_input.f32` (3072 complex f32 samples,
12 blocks of 256, including an all-zero block, an exactly-full-scale block, a
block containing a single sample far over nominal full scale, and a block near
the bottom of the exponent range) together with the encoded output at
`b` = 2, 4, 5, 6, 8, 16.

Those widths are chosen to cover every packing case: 4, 8 and 16 are
byte-aligned, 6 straddles two bytes at some sample indices, and 2 and 5 straddle
three. An implementation that reproduces all six byte for byte, and decodes them
to within half an exponent step, has the bit packing and the exponent rule
right. That is the whole interop surface.

The vectors were generated by an implementation independent of the two in
`reference/`; all three agree.

## 9. Versioning

Additive changes that a v1 reader can safely ignore — new JSON keys, new
`fill_rule` values — do not change the magic. Anything that alters the block
stream or the meaning of an existing header field takes a new magic (`BIQ2`).
Readers MUST reject unknown magic rather than guess.

## 10. Segmentation and time

Long recordings are normally split into a series of segment files. Each segment
is a complete, independently valid `.biq` file with its own header.

### 10.1 Absolute time

Within a file, the UTC time of complex sample `s` is:

```
t(s) = start_time_unix_ns + round(s * 1e9 / sample_rate_hz)
```

`start_time_unix_ns` of each segment **MUST be sampled from the clock at the
moment that segment's own first sample was captured.** It MUST NOT be derived
from an earlier segment's start time plus an elapsed-sample count.

This is the requirement that makes time survive damage, and it is easy to get
wrong because the arithmetic version looks equivalent and is simpler to write.
It is not equivalent. When each segment is stamped independently:

- A segment that is lost, corrupted, or truncated leaves every other segment's
  timing untouched. Nothing is chained, so nothing propagates.
- Timing error cannot accumulate across a recording. Whatever a segment gets
  wrong is corrected at the next segment boundary.
- A truncated segment keeps correct time for the samples it does hold, because
  the header is written at open and the time reference sits at byte 0.
- Corruption inside the block stream does not desync time: the fixed stride
  makes sample index positional, so samples after the damage are still at the
  right offset and therefore the right timestamp.

Derived timestamps lose all four. One segment written with a wrong sample count
shifts every timestamp after it, and nothing in the file indicates it happened.

### 10.2 Dropped samples

The one failure this format cannot absorb on its own is samples the source
dropped without the recorder noticing — a USB overrun, a ring-buffer overflow.
The file then looks continuous while every timestamp after the gap runs late,
which is worse than an obvious failure because it is plausible.

A recorder that detects a discontinuity in its input **SHOULD finalise the
current segment and start a new one at the gap.** The drop then becomes an
ordinary segment boundary carrying a fresh clock reading, which is correct by
construction and requires no cooperation from anything downstream.

A recorder that continues writing through a known gap MUST record it as a
`biq:events` entry (§3) with `gap_samples`. Note that the JSON is written at
finalise, so a recorder that dies before finalising loses the gap record even
though the samples are intact — which is why rolling a segment is preferred: it
commits the timing fact to the filesystem immediately instead of promising to
write it down later.

Readers MUST NOT assume that consecutive segments are gapless. Compare
`start_time_unix_ns` against the previous segment's start plus its duration; a
difference beyond the clock's own resolution is a real gap, not a rounding
artefact.

### 10.3 File naming

Not normative, but recommended, matching the convention SpectraFlux already
uses for its segmented recordings:

```
IQ_20260717_211008_100000000Hz_5000000sps_001.biq
```

Segments of one recording share the base name and differ only in the trailing
counter. Sort by name, not by date, to group a series.

The timestamp in the name is a convenience, not a substitute for
`start_time_unix_ns` — but it does mean a segment whose header is damaged can
still be placed in the series to within a second.
