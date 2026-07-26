# Sample recordings

Real `.biq` files, so you can read the format without owning a receiver or
getting your WAV past a converter first.

All three hold the same capture at the three profiles — one signal, three
mantissa widths, so the size and quality tradeoff is something you measure
rather than take on faith.

| file | profile | b | size | vs. int16 |
|---|---|---|---|---|
| `fm-band-95.5MHz-2Msps-b4-survey.biq` | `survey` | 4 | 10.00 MB | 3.98× |
| `fm-band-95.5MHz-2Msps-b6-general.biq` | `general` | 6 | 14.98 MB | 2.66× |
| `fm-band-95.5MHz-2Msps-b8-archive.biq` | `archive` | 8 | 19.96 MB | 2.00× |

## The capture

FM broadcast band, off an ADALM-Pluto:

```
centre       95.495516 MHz
sample rate  2 Msps complex
samples      9,958,400 in 38,900 blocks of 256  (4.979 s)
recorded     2026-07-26T19:14:42Z
```

The 2 MHz span holds a dozen stations. RDS decodes from all three files —
including `survey`, where four mantissa bits still carry the 57 kHz subcarrier
well enough for a decoder to lock. Nothing here is calibrated, so
`ref_dbm_full_scale` is left unset.

## Try one

```sh
cd ../reference && make
./biqconv info ../samples/fm-band-95.5MHz-2Msps-b6-general.biq
./biqconv ../samples/fm-band-95.5MHz-2Msps-b6-general.biq back.wav
```

`info` prints the whole recording's amplitude envelope from the exponent bytes
alone — 38,900 bytes touched, no mantissas decoded. `back.wav` is an ordinary
16-bit IQ WAV that any SDR tool will open.

## Checksums

```
c4ee4bece105f795d9e5ca5a3fe779cf85c21778dd4b10404d905ee75b2f67d5  fm-band-95.5MHz-2Msps-b4-survey.biq
ad8f63a38b5f8a9a0a45236091dec82cb20e1636cc1c0d5bdfae0de50a4b29c7  fm-band-95.5MHz-2Msps-b6-general.biq
cbbe56664098dc365a4199fee19b600600f902f815f573405d59560153a3cd18  fm-band-95.5MHz-2Msps-b8-archive.biq
```

## License

MIT, same as the rest of the repo. Use them in test suites, ports, and
comparisons without asking.
