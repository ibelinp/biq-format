// biq.hpp — reference encoder/decoder for the .biq block-floating-point IQ
// recording format. See SPEC.md; this file is the normative behaviour in code.
//
// Header-only, C++11, no dependencies beyond the standard library. Drop it in.
//
//   biq::Writer w;
//   biq::Header h = biq::Header::make(5e6, 100e6, /*b=*/6, /*N=*/256);
//   w.open("cap.biq", h, R"({"core:author":"me"})");
//   w.write(iq, n);          // interleaved float I,Q
//   w.close();
//
//   biq::Reader r;
//   r.open("cap.biq");
//   r.seek(1'000'000);
//   r.read(buf, 4096);
//   std::vector<int8_t> env = r.envelope();   // no mantissas decoded
//
// Samples are interleaved float (I,Q,I,Q,...). std::complex<float> is
// layout-compatible with float[2], so `(float*)cplx_vec.data()` is legal.

#ifndef BIQ_HPP
#define BIQ_HPP

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace biq {

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

static const char    kMagic[4]    = {'B', 'I', 'Q', '1'};
static const uint16_t kHeaderBytes = 64;
static const uint8_t kFlagReal    = 0x01;  // real-valued input: one plane, not two
static const int8_t  kZeroBlock   = -128;  // exponent sentinel for an all-zero block

// ---------------------------------------------------------------------------
// Block geometry
// ---------------------------------------------------------------------------

/// Bytes in one packed mantissa plane. N is a multiple of 8, so this is exact.
inline size_t plane_bytes(uint16_t block_size, uint8_t mantissa_bits) {
    return (size_t(block_size) * mantissa_bits + 7) / 8;
}

/// Bytes in one encoded block: exponent + one plane (real) or two (complex).
inline size_t block_bytes(uint16_t block_size, uint8_t mantissa_bits, bool real) {
    return 1 + (real ? 1u : 2u) * plane_bytes(block_size, mantissa_bits);
}

inline int32_t mantissa_max(uint8_t mantissa_bits) {
    return (int32_t(1) << (mantissa_bits - 1)) - 1;
}

// ---------------------------------------------------------------------------
// Bit packing — LSB-first, two's complement, b bits per value (SPEC §4)
// ---------------------------------------------------------------------------

/// Write value `v` (b bits, two's complement) at index `idx`. The plane must be
/// zeroed first; this ORs bits in, exactly like the Rust reference.
inline void pack(uint8_t* plane, size_t idx, uint8_t bits, int32_t v) {
    uint32_t raw    = uint32_t(v) & ((uint32_t(1) << bits) - 1);
    size_t   bitpos = idx * bits;
    size_t   byte   = bitpos / 8;
    uint32_t word   = raw << (bitpos % 8);  // spans at most 3 bytes for b <= 16
    while (word != 0) {
        plane[byte++] |= uint8_t(word & 0xff);
        word >>= 8;
    }
}

/// Read the b-bit two's-complement value at index `idx`. `plane_len` bounds the
/// read so the last value in a plane doesn't run past the end.
inline int32_t unpack(const uint8_t* plane, size_t plane_len, size_t idx, uint8_t bits) {
    size_t   bitpos = idx * bits;
    size_t   byte   = bitpos / 8;
    uint32_t word   = 0;
    for (size_t k = 0; k < 3 && byte + k < plane_len; ++k)
        word |= uint32_t(plane[byte + k]) << (8 * k);
    uint32_t raw  = (word >> (bitpos % 8)) & ((uint32_t(1) << bits) - 1);
    uint32_t sign = uint32_t(1) << (bits - 1);
    return int32_t(raw ^ sign) - int32_t(sign);
}

// ---------------------------------------------------------------------------
// Block codec
// ---------------------------------------------------------------------------

/// Smallest exponent e with peak <= m_max * 2^e (the `ceil` fill rule).
///
/// Computed with exact float compares rather than ceil(log2(...)) so every
/// implementation agrees byte for byte — a libm that rounds log2 the other way
/// at a power-of-two boundary would pick e±1. Decoders never depend on the
/// encoder's choice, but golden-vector tests do.
inline int8_t choose_exponent(float peak, uint8_t mantissa_bits) {
    if (!(peak > 0.0f)) return kZeroBlock;  // also catches NaN
    float m_max = float(mantissa_max(mantissa_bits));
    int   e     = int(std::ceil(std::log2(peak / m_max)));  // starting guess
    if (e < -127) e = -127;
    if (e > 127) e = 127;
    while (e < 127 && peak > m_max * std::ldexp(1.0f, e)) ++e;
    while (e > -127 && peak <= m_max * std::ldexp(1.0f, e - 1)) --e;
    return int8_t(e);
}

inline int32_t quantize(float v, int32_t m_max) {
    // Round-half-away-from-zero, matching Rust's f32::round. The clamp only
    // fires on the rounding edge, which is why `ceil` is safe.
    float  r = std::round(v);
    int32_t q = int32_t(r);
    if (q > m_max) q = m_max;
    if (q < -m_max - 1) q = -m_max - 1;
    return q;
}

/// Encode `block_size` complex samples (interleaved I,Q) into `out`, which must
/// be block_bytes() long. Returns the chosen exponent.
inline int8_t encode_block(const float* iq, uint16_t block_size, uint8_t mantissa_bits,
                           uint8_t* out) {
    float peak = 0.0f;
    for (size_t i = 0; i < size_t(block_size) * 2; ++i) {
        float a = std::fabs(iq[i]);
        if (a > peak) peak = a;
    }
    int8_t e = choose_exponent(peak, mantissa_bits);
    out[0]   = uint8_t(e);

    size_t pb = plane_bytes(block_size, mantissa_bits);
    std::memset(out + 1, 0, 2 * pb);
    if (e == kZeroBlock) return e;

    float   scale = std::ldexp(1.0f, -e);
    int32_t mmax  = mantissa_max(mantissa_bits);
    uint8_t* ip   = out + 1;
    uint8_t* qp   = out + 1 + pb;
    for (uint16_t i = 0; i < block_size; ++i) {
        pack(ip, i, mantissa_bits, quantize(iq[2 * i + 0] * scale, mmax));
        pack(qp, i, mantissa_bits, quantize(iq[2 * i + 1] * scale, mmax));
    }
    return e;
}

/// Real-valued variant: one plane, `block_size` real samples.
inline int8_t encode_block_real(const float* x, uint16_t block_size, uint8_t mantissa_bits,
                                uint8_t* out) {
    float peak = 0.0f;
    for (uint16_t i = 0; i < block_size; ++i) {
        float a = std::fabs(x[i]);
        if (a > peak) peak = a;
    }
    int8_t e = choose_exponent(peak, mantissa_bits);
    out[0]   = uint8_t(e);

    size_t pb = plane_bytes(block_size, mantissa_bits);
    std::memset(out + 1, 0, pb);
    if (e == kZeroBlock) return e;

    float   scale = std::ldexp(1.0f, -e);
    int32_t mmax  = mantissa_max(mantissa_bits);
    for (uint16_t i = 0; i < block_size; ++i)
        pack(out + 1, i, mantissa_bits, quantize(x[i] * scale, mmax));
    return e;
}

/// Decode one block into `block_size` interleaved I,Q floats.
inline void decode_block(const uint8_t* in, uint16_t block_size, uint8_t mantissa_bits,
                         float* iq) {
    int8_t e     = int8_t(in[0]);
    float  scale = (e == kZeroBlock) ? 0.0f : std::ldexp(1.0f, e);
    size_t pb    = plane_bytes(block_size, mantissa_bits);
    const uint8_t* ip = in + 1;
    const uint8_t* qp = in + 1 + pb;
    for (uint16_t i = 0; i < block_size; ++i) {
        iq[2 * i + 0] = float(unpack(ip, pb, i, mantissa_bits)) * scale;
        iq[2 * i + 1] = float(unpack(qp, pb, i, mantissa_bits)) * scale;
    }
}

inline void decode_block_real(const uint8_t* in, uint16_t block_size, uint8_t mantissa_bits,
                              float* x) {
    int8_t e     = int8_t(in[0]);
    float  scale = (e == kZeroBlock) ? 0.0f : std::ldexp(1.0f, e);
    size_t pb    = plane_bytes(block_size, mantissa_bits);
    for (uint16_t i = 0; i < block_size; ++i)
        x[i] = float(unpack(in + 1, pb, i, mantissa_bits)) * scale;
}

// ---------------------------------------------------------------------------
// Little-endian scalar I/O
// ---------------------------------------------------------------------------

namespace detail {

inline void put_u16(uint8_t* p, uint16_t v) {
    p[0] = uint8_t(v);
    p[1] = uint8_t(v >> 8);
}
inline void put_u32(uint8_t* p, uint32_t v) {
    for (int i = 0; i < 4; ++i) p[i] = uint8_t(v >> (8 * i));
}
inline void put_u64(uint8_t* p, uint64_t v) {
    for (int i = 0; i < 8; ++i) p[i] = uint8_t(v >> (8 * i));
}
inline void put_f64(uint8_t* p, double v) {
    uint64_t bits;
    std::memcpy(&bits, &v, 8);
    put_u64(p, bits);
}
inline uint16_t get_u16(const uint8_t* p) { return uint16_t(p[0]) | uint16_t(p[1]) << 8; }
inline uint32_t get_u32(const uint8_t* p) {
    uint32_t v = 0;
    for (int i = 0; i < 4; ++i) v |= uint32_t(p[i]) << (8 * i);
    return v;
}
inline uint64_t get_u64(const uint8_t* p) {
    uint64_t v = 0;
    for (int i = 0; i < 8; ++i) v |= uint64_t(p[i]) << (8 * i);
    return v;
}
inline double get_f64(const uint8_t* p) {
    uint64_t bits = get_u64(p);
    double   v;
    std::memcpy(&v, &bits, 8);
    return v;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Header (SPEC §2)
// ---------------------------------------------------------------------------

struct Header {
    uint16_t header_bytes       = kHeaderBytes;
    uint8_t  flags              = 0;
    uint8_t  channels           = 1;  // reserved: MUST be 1 in v1
    uint8_t  mantissa_bits      = 6;
    uint8_t  fill_rule          = 0;  // 0 = ceil
    uint16_t block_size         = 256;
    uint32_t block_stride       = 0;  // 0 here means "exact" — filled in by make()
    double   sample_rate_hz     = 0.0;
    double   center_freq_hz     = 0.0;
    int64_t  start_time_unix_ns = 0;
    double   ref_dbm_full_scale = std::numeric_limits<double>::quiet_NaN();
    uint64_t data_bytes         = 0;  // 0 = read to EOF
    uint32_t meta_bytes         = 0;

    bool   real() const { return (flags & kFlagReal) != 0; }
    size_t block_bytes() const { return biq::block_bytes(block_size, mantissa_bits, real()); }
    size_t data_start() const { return size_t(header_bytes) + meta_bytes; }

    /// Bytes on disk per complex (or real) sample.
    double bytes_per_sample() const { return double(block_stride) / block_size; }

    static Header make(double sample_rate_hz, double center_freq_hz, uint8_t mantissa_bits,
                       uint16_t block_size = 256, bool real = false) {
        Header h;
        h.sample_rate_hz = sample_rate_hz;
        h.center_freq_hz = center_freq_hz;
        h.mantissa_bits  = mantissa_bits;
        h.block_size     = block_size;
        h.flags          = real ? kFlagReal : 0;
        h.block_stride   = uint32_t(h.block_bytes());
        return h;
    }

    /// Round the stride up to a multiple of `align` for mmap/SIMD readers.
    /// Costs a little size; the default exact stride is what you usually want.
    void pad_stride(size_t align) {
        size_t bb    = block_bytes();
        block_stride = uint32_t((bb + align - 1) / align * align);
    }

    void serialize(uint8_t out[kHeaderBytes]) const {
        std::memset(out, 0, kHeaderBytes);
        std::memcpy(out, kMagic, 4);
        detail::put_u16(out + 4, header_bytes);
        out[6] = flags;
        out[7] = channels;
        out[8] = mantissa_bits;
        out[9] = fill_rule;
        detail::put_u16(out + 10, block_size);
        detail::put_u32(out + 12, block_stride);
        detail::put_f64(out + 16, sample_rate_hz);
        detail::put_f64(out + 24, center_freq_hz);
        detail::put_u64(out + 32, uint64_t(start_time_unix_ns));
        detail::put_f64(out + 40, ref_dbm_full_scale);
        detail::put_u64(out + 48, data_bytes);
        detail::put_u32(out + 56, meta_bytes);
    }

    /// Parse and validate. Returns false and sets `err` on any violation of §2.
    bool parse(const uint8_t* in, size_t len, std::string* err) {
        const char* e = nullptr;
        if (len < kHeaderBytes)
            e = "file shorter than a header";
        else if (std::memcmp(in, kMagic, 4) != 0)
            e = "bad magic (not a BIQ1 file)";
        if (e) {
            if (err) *err = e;
            return false;
        }
        header_bytes       = detail::get_u16(in + 4);
        flags              = in[6];
        channels           = in[7];
        mantissa_bits      = in[8];
        fill_rule          = in[9];
        block_size         = detail::get_u16(in + 10);
        block_stride       = detail::get_u32(in + 12);
        sample_rate_hz     = detail::get_f64(in + 16);
        center_freq_hz     = detail::get_f64(in + 24);
        start_time_unix_ns = int64_t(detail::get_u64(in + 32));
        ref_dbm_full_scale = detail::get_f64(in + 40);
        data_bytes         = detail::get_u64(in + 48);
        meta_bytes         = detail::get_u32(in + 56);

        if (header_bytes < kHeaderBytes)
            e = "header_bytes below the v1 minimum";
        else if (channels != 1)
            e = "channels != 1 (reserved in v1; refusing to guess a layout)";
        else if (flags & ~kFlagReal)
            e = "unknown flag bits set";
        else if (mantissa_bits < 2 || mantissa_bits > 16)
            e = "mantissa_bits outside 2..16";
        else if (block_size == 0 || block_size % 8 != 0)
            e = "block_size zero or not a multiple of 8";
        else if (block_stride < block_bytes())
            e = "block_stride smaller than the encoded block";
        else if (meta_bytes % 8 != 0)
            e = "meta_bytes not a multiple of 8";
        else if (detail::get_u32(in + 60) != 0)
            e = "reserved field is non-zero";
        if (e) {
            if (err) *err = e;
            return false;
        }
        return true;
    }
};

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

class Writer {
public:
    ~Writer() { close(); }

    /// `meta` is a JSON object (or empty). It is padded to a multiple of 8 with
    /// spaces, which is legal trailing whitespace for any JSON parser.
    bool open(const std::string& path, const Header& h, const std::string& meta = "") {
        close();
        hdr_  = h;
        meta_ = meta;
        while (meta_.size() % 8 != 0) meta_.push_back(' ');
        hdr_.meta_bytes = uint32_t(meta_.size());
        if (hdr_.block_stride == 0) hdr_.block_stride = uint32_t(hdr_.block_bytes());
        hdr_.data_bytes = 0;  // backfilled at close; 0 stays valid if we crash

        f_ = std::fopen(path.c_str(), "wb");
        if (!f_) {
            err_ = "cannot open " + path + " for writing";
            return false;
        }
        uint8_t hb[kHeaderBytes];
        hdr_.serialize(hb);
        if (std::fwrite(hb, 1, kHeaderBytes, f_) != kHeaderBytes ||
            (!meta_.empty() && std::fwrite(meta_.data(), 1, meta_.size(), f_) != meta_.size())) {
            err_ = "short write on header";
            close_file();
            return false;
        }
        blockbuf_.assign(hdr_.block_stride, 0);
        pending_.clear();
        pending_.reserve(size_t(hdr_.block_size) * 2);
        data_bytes_ = 0;
        return true;
    }

    /// Append `n` complex samples (interleaved I,Q; or `n` reals if the header
    /// says real). Samples are buffered until a whole block is available.
    bool write(const float* samples, size_t n) {
        size_t comp = hdr_.real() ? 1 : 2;
        size_t full = size_t(hdr_.block_size) * comp;
        for (size_t i = 0; i < n * comp; ++i) {
            pending_.push_back(samples[i]);
            if (pending_.size() == full && !flush_block()) return false;
        }
        return true;
    }

    /// Finish the file: zero-pad any partial block, backfill `data_bytes`.
    bool close() {
        if (!f_) return true;
        bool ok = true;
        if (!pending_.empty()) {
            pending_.resize(size_t(hdr_.block_size) * (hdr_.real() ? 1 : 2), 0.0f);
            ok = flush_block();
        }
        if (ok) {
            hdr_.data_bytes = data_bytes_;
            uint8_t hb[kHeaderBytes];
            hdr_.serialize(hb);
            ok = std::fseek(f_, 0, SEEK_SET) == 0 &&
                 std::fwrite(hb, 1, kHeaderBytes, f_) == kHeaderBytes;
            if (!ok) err_ = "could not backfill data_bytes";
        }
        close_file();
        return ok;
    }

    const Header&      header() const { return hdr_; }
    const std::string& error() const { return err_; }
    uint64_t           blocks_written() const { return data_bytes_ / hdr_.block_stride; }

private:
    bool flush_block() {
        std::memset(blockbuf_.data(), 0, blockbuf_.size());
        if (hdr_.real())
            encode_block_real(pending_.data(), hdr_.block_size, hdr_.mantissa_bits,
                              blockbuf_.data());
        else
            encode_block(pending_.data(), hdr_.block_size, hdr_.mantissa_bits, blockbuf_.data());
        pending_.clear();
        if (std::fwrite(blockbuf_.data(), 1, blockbuf_.size(), f_) != blockbuf_.size()) {
            err_ = "short write on block";
            return false;
        }
        data_bytes_ += blockbuf_.size();
        return true;
    }
    void close_file() {
        if (f_) std::fclose(f_);
        f_ = nullptr;
    }

    std::FILE*           f_ = nullptr;
    Header               hdr_;
    std::string          meta_, err_;
    std::vector<float>   pending_;
    std::vector<uint8_t> blockbuf_;
    uint64_t             data_bytes_ = 0;
};

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

class Reader {
public:
    ~Reader() {
        if (f_) std::fclose(f_);
    }

    bool open(const std::string& path) {
        if (f_) std::fclose(f_);
        f_ = std::fopen(path.c_str(), "rb");
        if (!f_) {
            err_ = "cannot open " + path;
            return false;
        }
        uint8_t hb[kHeaderBytes];
        if (std::fread(hb, 1, kHeaderBytes, f_) != kHeaderBytes) {
            err_ = "file shorter than a header";
            return false;
        }
        if (!hdr_.parse(hb, kHeaderBytes, &err_)) return false;

        // Skip any header extension a future minor version added.
        if (hdr_.header_bytes > kHeaderBytes &&
            std::fseek(f_, long(hdr_.header_bytes), SEEK_SET) != 0) {
            err_ = "truncated header extension";
            return false;
        }
        meta_.assign(hdr_.meta_bytes, '\0');
        if (hdr_.meta_bytes &&
            std::fread(&meta_[0], 1, hdr_.meta_bytes, f_) != hdr_.meta_bytes) {
            err_ = "truncated metadata";
            return false;
        }

        std::fseek(f_, 0, SEEK_END);
        file_bytes_ = uint64_t(std::ftell(f_));
        uint64_t avail =
            file_bytes_ > hdr_.data_start() ? file_bytes_ - hdr_.data_start() : 0;
        // data_bytes == 0 means read-to-EOF (a recorder that died mid-write).
        // Trust the smaller of the two, so a truncated file just reads short.
        uint64_t declared = hdr_.data_bytes ? hdr_.data_bytes : avail;
        data_bytes_       = declared < avail ? declared : avail;
        blocks_           = data_bytes_ / hdr_.block_stride;  // partial tail block dropped

        blockbuf_.assign(hdr_.block_stride, 0);
        scratch_.assign(size_t(hdr_.block_size) * (hdr_.real() ? 1 : 2), 0.0f);
        cached_ = kNoBlock;
        pos_    = 0;
        return true;
    }

    const Header&      header() const { return hdr_; }
    const std::string& metadata() const { return meta_; }
    const std::string& error() const { return err_; }
    uint64_t           blocks() const { return blocks_; }
    uint64_t           samples() const { return blocks_ * hdr_.block_size; }
    double duration_s() const {
        return hdr_.sample_rate_hz > 0 ? double(samples()) / hdr_.sample_rate_hz : 0.0;
    }
    /// True if the file ends mid-block (recorder crashed, disk filled).
    bool truncated() const {
        return hdr_.data_bytes != 0 &&
               hdr_.data_bytes > file_bytes_ - hdr_.data_start();
    }

    /// Byte offset of the block containing sample `s` — SPEC §7. Pure arithmetic.
    uint64_t offset_of(uint64_t s) const {
        return hdr_.data_start() + (s / hdr_.block_size) * hdr_.block_stride;
    }

    /// Seek to an exact sample — O(1), no index, no scan.
    bool seek(uint64_t sample) {
        if (sample > samples()) {
            err_ = "seek past end";
            return false;
        }
        pos_ = sample;
        return true;
    }

    uint64_t tell() const { return pos_; }

    /// Read up to `n` samples into `out` (interleaved I,Q, so 2n floats — or n
    /// floats for a real file). Returns the number of samples actually read,
    /// which is short only at end of file.
    size_t read(float* out, size_t n) {
        size_t comp = hdr_.real() ? 1 : 2;
        size_t done = 0;
        while (done < n && pos_ < samples()) {
            uint64_t blk = pos_ / hdr_.block_size;
            if (!load_block(blk)) break;
            size_t off  = size_t(pos_ % hdr_.block_size);
            size_t take = hdr_.block_size - off;
            if (take > n - done) take = n - done;
            std::memcpy(out + done * comp, scratch_.data() + off * comp,
                        take * comp * sizeof(float));
            done += take;
            pos_ += take;
        }
        return done;
    }

    /// The free amplitude envelope (SPEC §7): one exponent per block, read with
    /// a strided pass that touches `blocks()` bytes and decodes no mantissas.
    /// Value kZeroBlock (-128) marks an all-zero block.
    std::vector<int8_t> envelope() {
        std::vector<int8_t> env;
        env.reserve(size_t(blocks_));
        for (uint64_t k = 0; k < blocks_; ++k) {
            if (std::fseek(f_, long(hdr_.data_start() + k * hdr_.block_stride), SEEK_SET) != 0)
                break;
            int c = std::fgetc(f_);
            if (c == EOF) break;
            env.push_back(int8_t(uint8_t(c)));
        }
        cached_ = kNoBlock;  // we moved the file pointer behind load_block's back
        return env;
    }

    /// Envelope in dBFS relative to full scale (±1.0), for drawing. Blocks whose
    /// peak is unknown are placed at their exponent's ceiling, which is what the
    /// exponent means: the block's peak is within (0.5, 1.0] of m_max * 2^e.
    std::vector<float> envelope_dbfs() {
        std::vector<int8_t> e = envelope();
        std::vector<float>  db(e.size());
        float               mmax = float(mantissa_max(hdr_.mantissa_bits));
        for (size_t i = 0; i < e.size(); ++i)
            db[i] = (e[i] == kZeroBlock)
                        ? -std::numeric_limits<float>::infinity()
                        : 20.0f * std::log10(mmax * std::ldexp(1.0f, e[i]));
        return db;
    }

private:
    static const uint64_t kNoBlock = ~uint64_t(0);

    /// Decode block `blk` into the scratch buffer, unless it's already there.
    bool load_block(uint64_t blk) {
        if (cached_ == blk) return true;
        if (std::fseek(f_, long(hdr_.data_start() + blk * hdr_.block_stride), SEEK_SET) != 0)
            return false;
        if (std::fread(blockbuf_.data(), 1, blockbuf_.size(), f_) != blockbuf_.size())
            return false;
        if (hdr_.real())
            decode_block_real(blockbuf_.data(), hdr_.block_size, hdr_.mantissa_bits,
                              scratch_.data());
        else
            decode_block(blockbuf_.data(), hdr_.block_size, hdr_.mantissa_bits, scratch_.data());
        cached_ = blk;
        return true;
    }

    std::FILE*           f_ = nullptr;
    Header               hdr_;
    std::string          meta_, err_;
    std::vector<uint8_t> blockbuf_;
    std::vector<float>   scratch_;
    uint64_t             cached_     = kNoBlock;
    uint64_t             file_bytes_ = 0, data_bytes_ = 0, blocks_ = 0, pos_ = 0;
};

}  // namespace biq

#endif  // BIQ_HPP
