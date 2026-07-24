#!/usr/bin/env python3
"""Reference reader/writer for the .biq block-floating-point IQ format.

Dependency-free and deliberately literal — this is the spec in Python, not a
fast implementation. See SPEC.md. The C++ in biq.hpp is the one to port from if
you want speed; this one is for reading a file in a notebook and for checking
that another implementation agrees.

    python3 biq.py --test           # round-trip + golden-vector check
    python3 biq.py info cap.biq     # header, geometry, envelope
"""

import math
import struct
import sys

MAGIC = b"BIQ1"
HEADER_BYTES = 64
FLAG_REAL = 0x01
ZERO_BLOCK = -128

# --------------------------------------------------------------------------
# Block codec
# --------------------------------------------------------------------------


def plane_bytes(block_size, bits):
    return (block_size * bits + 7) // 8


def block_bytes(block_size, bits, real=False):
    return 1 + (1 if real else 2) * plane_bytes(block_size, bits)


def mantissa_max(bits):
    return (1 << (bits - 1)) - 1


def choose_exponent(peak, bits):
    """Smallest e with peak <= m_max * 2**e — the `ceil` fill rule.

    Done with exact compares rather than ceil(log2(...)) so every
    implementation lands on the same exponent at power-of-two boundaries.
    """
    if not peak > 0.0:
        return ZERO_BLOCK
    m_max = float(mantissa_max(bits))
    e = int(math.ceil(math.log2(peak / m_max)))
    e = max(-127, min(127, e))
    while e < 127 and peak > m_max * 2.0**e:
        e += 1
    while e > -127 and peak <= m_max * 2.0 ** (e - 1):
        e -= 1
    return e


def _pack(plane, idx, bits, value):
    raw = value & ((1 << bits) - 1)
    bitpos = idx * bits
    byte, shift = bitpos // 8, bitpos % 8
    word = raw << shift
    while word:
        plane[byte] |= word & 0xFF
        word >>= 8
        byte += 1


def _unpack(plane, idx, bits):
    bitpos = idx * bits
    byte, shift = bitpos // 8, bitpos % 8
    word = 0
    for k in range(3):
        if byte + k < len(plane):
            word |= plane[byte + k] << (8 * k)
    raw = (word >> shift) & ((1 << bits) - 1)
    sign = 1 << (bits - 1)
    return (raw ^ sign) - sign


def _quantize(v, m_max):
    # Round half away from zero, matching Rust f32::round and C++ std::round.
    # Python's round() is banker's rounding, so it cannot be used here.
    q = math.floor(v + 0.5) if v >= 0 else math.ceil(v - 0.5)
    return max(-m_max - 1, min(m_max, int(q)))


def _f32(x):
    """Round a Python float to f32, so we quantise the same values Rust/C++ do."""
    return struct.unpack("<f", struct.pack("<f", x))[0]


def encode_block(samples, bits, real=False):
    """samples: list of (i, q) tuples, or plain floats when real=True."""
    n = len(samples)
    if real:
        peak = max((abs(_f32(s)) for s in samples), default=0.0)
    else:
        peak = max((max(abs(_f32(i)), abs(_f32(q))) for i, q in samples), default=0.0)
    e = choose_exponent(peak, bits)
    pb = plane_bytes(n, bits)
    out = bytearray(1 + (1 if real else 2) * pb)
    out[0] = e & 0xFF
    if e == ZERO_BLOCK:
        return bytes(out)

    scale = 2.0**-e
    m_max = mantissa_max(bits)
    if real:
        plane = bytearray(pb)
        for k, s in enumerate(samples):
            _pack(plane, k, bits, _quantize(_f32(_f32(s) * scale), m_max))
        out[1:] = plane
    else:
        ip = bytearray(pb)
        qp = bytearray(pb)
        for k, (i, q) in enumerate(samples):
            _pack(ip, k, bits, _quantize(_f32(_f32(i) * scale), m_max))
            _pack(qp, k, bits, _quantize(_f32(_f32(q) * scale), m_max))
        out[1 : 1 + pb] = ip
        out[1 + pb :] = qp
    return bytes(out)


def decode_block(data, block_size, bits, real=False):
    e = struct.unpack("<b", data[0:1])[0]
    scale = 0.0 if e == ZERO_BLOCK else 2.0**e
    pb = plane_bytes(block_size, bits)
    if real:
        plane = data[1 : 1 + pb]
        return [_f32(_unpack(plane, k, bits) * scale) for k in range(block_size)]
    ip = data[1 : 1 + pb]
    qp = data[1 + pb : 1 + 2 * pb]
    return [
        (_f32(_unpack(ip, k, bits) * scale), _f32(_unpack(qp, k, bits) * scale))
        for k in range(block_size)
    ]


# --------------------------------------------------------------------------
# Container
# --------------------------------------------------------------------------


class Header:
    __slots__ = (
        "header_bytes flags channels mantissa_bits fill_rule block_size block_stride "
        "sample_rate_hz center_freq_hz start_time_unix_ns ref_dbm_full_scale "
        "data_bytes meta_bytes"
    ).split()

    def __init__(self, sample_rate_hz=0.0, center_freq_hz=0.0, mantissa_bits=6,
                 block_size=256, real=False):
        self.header_bytes = HEADER_BYTES
        self.flags = FLAG_REAL if real else 0
        self.channels = 1
        self.mantissa_bits = mantissa_bits
        self.fill_rule = 0
        self.block_size = block_size
        self.sample_rate_hz = sample_rate_hz
        self.center_freq_hz = center_freq_hz
        self.start_time_unix_ns = 0
        self.ref_dbm_full_scale = float("nan")
        self.data_bytes = 0
        self.meta_bytes = 0
        self.block_stride = block_bytes(block_size, mantissa_bits, real)

    @property
    def real(self):
        return bool(self.flags & FLAG_REAL)

    @property
    def block_bytes(self):
        return block_bytes(self.block_size, self.mantissa_bits, self.real)

    @property
    def data_start(self):
        return self.header_bytes + self.meta_bytes

    @property
    def bytes_per_sample(self):
        return self.block_stride / self.block_size

    def pack(self):
        return (
            MAGIC
            + struct.pack(
                "<HBBBBHIddqdQI4x",
                self.header_bytes, self.flags, self.channels, self.mantissa_bits,
                self.fill_rule, self.block_size, self.block_stride,
                self.sample_rate_hz, self.center_freq_hz, self.start_time_unix_ns,
                self.ref_dbm_full_scale, self.data_bytes, self.meta_bytes,
            )
        )

    @classmethod
    def unpack(cls, buf):
        if len(buf) < HEADER_BYTES:
            raise ValueError("file shorter than a header")
        if buf[:4] != MAGIC:
            raise ValueError("bad magic (not a BIQ1 file)")
        h = cls()
        (h.header_bytes, h.flags, h.channels, h.mantissa_bits, h.fill_rule,
         h.block_size, h.block_stride, h.sample_rate_hz, h.center_freq_hz,
         h.start_time_unix_ns, h.ref_dbm_full_scale, h.data_bytes,
         h.meta_bytes) = struct.unpack("<HBBBBHIddqdQI4x", buf[4:HEADER_BYTES])

        if h.header_bytes < HEADER_BYTES:
            raise ValueError("header_bytes below the v1 minimum")
        if h.channels != 1:
            raise ValueError("channels != 1 (reserved in v1; refusing to guess a layout)")
        if h.flags & ~FLAG_REAL:
            raise ValueError("unknown flag bits set")
        if not 2 <= h.mantissa_bits <= 16:
            raise ValueError("mantissa_bits outside 2..16")
        if h.block_size == 0 or h.block_size % 8:
            raise ValueError("block_size zero or not a multiple of 8")
        if h.block_stride < h.block_bytes:
            raise ValueError("block_stride smaller than the encoded block")
        if h.meta_bytes % 8:
            raise ValueError("meta_bytes not a multiple of 8")
        if struct.unpack("<I", buf[60:64])[0] != 0:
            raise ValueError("reserved field is non-zero")
        return h


def write(path, samples, header, meta=""):
    """Write all `samples` (list of (i,q), or floats if header.real) to `path`."""
    meta = meta + " " * (-len(meta) % 8)
    header.meta_bytes = len(meta)
    n = header.block_size
    blocks = (len(samples) + n - 1) // n
    pad = (0.0 if header.real else (0.0, 0.0))
    with open(path, "wb") as f:
        f.write(header.pack())
        f.write(meta.encode())
        written = 0
        for k in range(blocks):
            chunk = list(samples[k * n : (k + 1) * n])
            chunk += [pad] * (n - len(chunk))
            blk = encode_block(chunk, header.mantissa_bits, header.real)
            blk += b"\0" * (header.block_stride - len(blk))
            f.write(blk)
            written += len(blk)
        header.data_bytes = written
        f.seek(0)
        f.write(header.pack())
    return blocks


class Reader:
    """Random-access reader. Seeking is arithmetic — see SPEC §7."""

    def __init__(self, path):
        self.f = open(path, "rb")
        self.header = Header.unpack(self.f.read(HEADER_BYTES))
        self.f.seek(self.header.header_bytes)
        self.metadata = self.f.read(self.header.meta_bytes).decode("utf-8", "replace")
        self.f.seek(0, 2)
        self.file_bytes = self.f.tell()
        avail = max(0, self.file_bytes - self.header.data_start)
        declared = self.header.data_bytes or avail  # 0 = read to EOF
        self.data_bytes = min(declared, avail)
        self.blocks = self.data_bytes // self.header.block_stride
        self.samples = self.blocks * self.header.block_size

    @property
    def truncated(self):
        return bool(self.header.data_bytes) and \
            self.header.data_bytes > self.file_bytes - self.header.data_start

    @property
    def duration_s(self):
        return self.samples / self.header.sample_rate_hz if self.header.sample_rate_hz else 0.0

    def offset_of(self, sample):
        return self.header.data_start + (sample // self.header.block_size) * self.header.block_stride

    def read_block(self, k):
        """Decode block k. O(1) — no scan, no index."""
        self.f.seek(self.header.data_start + k * self.header.block_stride)
        raw = self.f.read(self.header.block_bytes)
        return decode_block(raw, self.header.block_size, self.header.mantissa_bits,
                            self.header.real)

    def read(self, start=0, count=None):
        if count is None:
            count = self.samples - start
        out = []
        s = start
        while len(out) < count and s < self.samples:
            blk = self.read_block(s // self.header.block_size)
            off = s % self.header.block_size
            take = min(len(blk) - off, count - len(out))
            out.extend(blk[off : off + take])
            s += take
        return out

    def envelope(self):
        """One exponent per block, no mantissas decoded — SPEC §7."""
        env = []
        for k in range(self.blocks):
            self.f.seek(self.header.data_start + k * self.header.block_stride)
            b = self.f.read(1)
            if not b:
                break
            env.append(struct.unpack("<b", b)[0])
        return env

    def envelope_dbfs(self):
        m = float(mantissa_max(self.header.mantissa_bits))
        return [
            float("-inf") if e == ZERO_BLOCK else 20.0 * math.log10(m * 2.0**e)
            for e in self.envelope()
        ]

    def close(self):
        self.f.close()


# --------------------------------------------------------------------------
# CLI / self-test
# --------------------------------------------------------------------------


def _info(path):
    r = Reader(path)
    h = r.header
    print(path)
    print(f"  geometry      b={h.mantissa_bits}  N={h.block_size}  stride={h.block_stride}"
          f"{' (exact)' if h.block_stride == h.block_bytes else ' (padded)'}"
          f"  {'real' if h.real else 'complex'}")
    print(f"  sample rate   {h.sample_rate_hz:.10g} Hz")
    print(f"  centre        {h.center_freq_hz:.10g} Hz")
    print(f"  samples       {r.samples} in {r.blocks} blocks  ({r.duration_s:.3f} s)")
    print(f"  size          {r.blocks * h.block_stride / 1e6:.2f} MB "
          f"({h.bytes_per_sample:.4f} bytes/sample, {4.0 / h.bytes_per_sample:.2f}x vs int16)"
          f"{'  TRUNCATED' if r.truncated else ''}")
    if h.meta_bytes:
        print(f"  metadata      {r.metadata.strip()}")
    db = [d for d in r.envelope_dbfs() if d != float("-inf")]
    if db:
        print(f"  envelope      {max(db):.0f} .. {min(db):.0f} dBFS "
              f"({len(db)} exponents, zero mantissas decoded)")
    r.close()


def _test():
    import os
    import tempfile

    fails = [0]

    def check(cond, msg):
        if not cond:
            fails[0] += 1
            print(f"  FAIL {msg}")

    print("bit packing")
    plane = bytearray(2)
    _pack(plane, 0, 5, 3)
    _pack(plane, 1, 5, 7)
    check(plane[0] == 0xE3, f"byte0 = 0x{plane[0]:02x}, want 0xE3")
    check(_unpack(plane, 0, 5) == 3 and _unpack(plane, 1, 5) == 7, "unpack mismatch")

    print("golden vectors match the Rust encoder")
    here = os.path.dirname(os.path.abspath(__file__))
    gdir = os.path.join(here, "golden")
    if not os.path.isdir(gdir):
        print("  (golden/ not present — skipped)")
    else:
        raw = open(os.path.join(gdir, "golden_input.f32"), "rb").read()
        flat = struct.unpack(f"<{len(raw) // 4}f", raw)
        samples = list(zip(flat[0::2], flat[1::2]))
        N = 256
        for bits in (2, 4, 5, 6, 8, 16):
            want = open(os.path.join(gdir, f"golden_b{bits}.bin"), "rb").read()
            got = b"".join(
                encode_block(samples[k * N : (k + 1) * N], bits)
                for k in range(len(samples) // N)
            )
            bad = next((i for i in range(min(len(got), len(want))) if got[i] != want[i]), None)
            check(len(got) == len(want) and bad is None,
                  f"b={bits}: byte {bad} differs" if bad is not None
                  else f"b={bits}: {len(got)} bytes, want {len(want)}")

    print("file round-trip, seek, envelope")
    n = 256 * 20 + 33
    sig = [
        (
            _f32(0.9 * 10 ** (-3.0 * i / n) * math.cos(2 * math.pi * 0.0100073 * i)),
            _f32(0.9 * 10 ** (-3.0 * i / n) * math.sin(2 * math.pi * 0.0100073 * i)),
        )
        for i in range(n)
    ]
    path = os.path.join(tempfile.gettempdir(), "biq_py_test.biq")
    h = Header(5e6, 100e6, mantissa_bits=8)
    write(path, sig, h, '{"core:author":"test"}')

    r = Reader(path)
    check(r.blocks == 21, f"blocks = {r.blocks}, want 21")
    check(r.header.sample_rate_hz == 5e6, "sample rate lost")
    check(not r.truncated, "clean file reported truncated")
    check('"core:author"' in r.metadata, "metadata lost")
    dec = r.read(0, n)
    e = sum((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 for a, b in zip(sig, dec))
    s = sum(a[0] ** 2 + a[1] ** 2 for a in sig)
    snr = 10 * math.log10(s / e)
    check(snr > 37.0, f"SQNR {snr:.1f} dB, expected > 37")
    for p in (0, 1, 255, 256, 1000, 4000):
        check(r.read(p, 3) == dec[p : p + 3], f"seek to {p} disagreed with a linear read")
    env = r.envelope()
    check(len(env) == r.blocks, "envelope length")
    check(all(env[i] <= env[i - 1] for i in range(1, len(env))), "envelope did not decay")
    r.close()
    os.remove(path)

    print("real-valued mode")
    xr = [_f32(0.8 * math.sin(2 * math.pi * 0.013 * i)) for i in range(256 * 4)]
    hr = Header(2e6, 0.0, mantissa_bits=8, real=True)
    check(hr.block_bytes == 257, f"real block = {hr.block_bytes} bytes, want 257")
    write(path, xr, hr)
    r = Reader(path)
    check(r.header.real, "real flag lost")
    dr = r.read(0, len(xr))
    e = sum((a - b) ** 2 for a, b in zip(xr, dr))
    s = sum(a * a for a in xr)
    check(10 * math.log10(s / e) > 37.0, "real-mode SQNR too low")
    r.close()
    os.remove(path)

    print(f"\n{'ok' if not fails[0] else str(fails[0]) + ' failed'}")
    return 1 if fails[0] else 0


def main(argv):
    if len(argv) > 1 and argv[1] == "--test":
        return _test()
    if len(argv) > 2 and argv[1] == "info":
        _info(argv[2])
        return 0
    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
