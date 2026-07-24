// biq_test.cpp — self-test for biq.hpp.
//
//   make test
//
// The interesting part is `golden_vectors_match_rust`: the vectors in golden/
// were produced by an independent Rust implementation, so agreeing with them
// byte for byte is the actual interop guarantee. Everything else is container
// behaviour.

#include "biq.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------

static int g_checks = 0, g_fails = 0;
static const char* g_case = "";

#define CHECK(cond, ...)                                              \
    do {                                                              \
        ++g_checks;                                                   \
        if (!(cond)) {                                                \
            ++g_fails;                                                \
            std::printf("  FAIL %s:%d  %s\n       ", __FILE__, __LINE__, #cond); \
            std::printf(__VA_ARGS__);                                 \
            std::printf("\n");                                        \
        }                                                             \
    } while (0)

static void begin(const char* name) {
    g_case = name;
    std::printf("%s\n", name);
}

static std::string g_dir = "golden";

static std::vector<uint8_t> slurp(const std::string& path) {
    std::vector<uint8_t> out;
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) {
        std::printf("  FATAL cannot open %s\n", path.c_str());
        std::exit(2);
    }
    std::fseek(f, 0, SEEK_END);
    out.resize(size_t(std::ftell(f)));
    std::fseek(f, 0, SEEK_SET);
    if (!out.empty() && std::fread(out.data(), 1, out.size(), f) != out.size()) out.clear();
    std::fclose(f);
    return out;
}

static std::string tmp_path(const char* name) {
    const char* d = std::getenv("TMPDIR");
    std::string dir = d ? d : "/tmp/";
    if (!dir.empty() && dir[dir.size() - 1] != '/') dir += '/';
    return dir + name;
}

static double sqnr_db(const std::vector<float>& a, const std::vector<float>& b, size_t from,
                      size_t to) {
    double sig = 0, err = 0;
    for (size_t i = from; i < to; ++i) {
        sig += double(a[i]) * a[i];
        double d = double(a[i]) - b[i];
        err += d * d;
    }
    return err > 0 ? 10.0 * std::log10(sig / err) : 1e9;
}

// ---------------------------------------------------------------------------

static const uint16_t kBlock = 256;
static const uint8_t  kBits[] = {2, 4, 5, 6, 8, 16};

/// The one that matters: our encoder must produce the exact bytes the Rust
/// encoder produced, and decode Rust's bytes to the same floats.
static void golden_vectors_match_rust() {
    begin("golden vectors match the Rust encoder");

    std::vector<uint8_t> raw = slurp(g_dir + "/golden_input.f32");
    size_t n = raw.size() / 8;  // complex f32
    std::vector<float> in(n * 2);
    std::memcpy(in.data(), raw.data(), raw.size());
    CHECK(n == size_t(kBlock) * 12, "expected 3072 samples, got %zu", n);

    for (size_t bi = 0; bi < sizeof(kBits); ++bi) {
        uint8_t b = kBits[bi];
        std::vector<uint8_t> want = slurp(g_dir + "/golden_b" + std::to_string(b) + ".bin");
        size_t bb = biq::block_bytes(kBlock, b, false);
        size_t nb = n / kBlock;
        CHECK(want.size() == nb * bb, "b=%u golden is %zu bytes, expected %zu", b, want.size(),
              nb * bb);

        // Encode → compare bytes.
        std::vector<uint8_t> got(nb * bb, 0);
        for (size_t k = 0; k < nb; ++k)
            biq::encode_block(in.data() + k * kBlock * 2, kBlock, b, got.data() + k * bb);

        size_t first_bad = size_t(-1);
        for (size_t i = 0; i < got.size() && i < want.size(); ++i)
            if (got[i] != want[i]) {
                first_bad = i;
                break;
            }
        CHECK(first_bad == size_t(-1), "b=%u byte %zu: got 0x%02x want 0x%02x (block %zu)", b,
              first_bad, first_bad == size_t(-1) ? 0 : got[first_bad],
              first_bad == size_t(-1) ? 0 : want[first_bad], first_bad / bb);

        // Decode Rust's bytes → must round-trip within half an exponent step.
        std::vector<float> dec(n * 2);
        for (size_t k = 0; k < nb; ++k)
            biq::decode_block(want.data() + k * bb, kBlock, b, dec.data() + k * kBlock * 2);

        size_t worst = 0;
        double worst_err = 0;
        for (size_t k = 0; k < nb; ++k) {
            int8_t e = int8_t(want[k * bb]);
            if (e == biq::kZeroBlock) continue;
            double step = std::ldexp(1.0, e);
            for (size_t i = 0; i < size_t(kBlock) * 2; ++i) {
                size_t idx = k * kBlock * 2 + i;
                double d   = std::fabs(double(in[idx]) - dec[idx]);
                if (d / step > worst_err) {
                    worst_err = d / step;
                    worst     = idx;
                }
            }
        }
        CHECK(worst_err <= 0.5001, "b=%u: error %.4f exponent steps at sample %zu", b, worst_err,
              worst / 2);
    }
}

static void packing_is_lsb_first() {
    begin("bit packing is LSB-first two's complement");
    uint8_t plane[2] = {0, 0};
    biq::pack(plane, 0, 5, 3);
    biq::pack(plane, 1, 5, 7);
    // value0 in the low 5 bits, value1's low 3 bits above it.
    CHECK(plane[0] == 0xE3, "byte0 = 0x%02x, want 0xE3", plane[0]);
    CHECK(plane[1] == 0x00, "byte1 = 0x%02x, want 0x00", plane[1]);
    CHECK(biq::unpack(plane, 2, 0, 5) == 3, "unpack[0] = %d", biq::unpack(plane, 2, 0, 5));
    CHECK(biq::unpack(plane, 2, 1, 5) == 7, "unpack[1] = %d", biq::unpack(plane, 2, 1, 5));

    uint8_t p2[2] = {0, 0};
    biq::pack(p2, 0, 5, -13);
    CHECK(biq::unpack(p2, 2, 0, 5) == -13, "negative round-trip gave %d",
          biq::unpack(p2, 2, 0, 5));

    // Every width, every alignment, full value range — the last value in a
    // plane must not read past the end.
    for (uint8_t b = 2; b <= 16; ++b) {
        size_t pb = biq::plane_bytes(64, b);
        std::vector<uint8_t> plane2(pb, 0);
        int32_t lo = -(int32_t(1) << (b - 1)), hi = (int32_t(1) << (b - 1)) - 1;
        std::vector<int32_t> vals(64);
        for (size_t i = 0; i < 64; ++i) {
            vals[i] = lo + int32_t(i * 7919 % uint32_t(hi - lo + 1));
            biq::pack(plane2.data(), i, b, vals[i]);
        }
        bool ok = true;
        for (size_t i = 0; i < 64; ++i)
            if (biq::unpack(plane2.data(), pb, i, b) != vals[i]) ok = false;
        CHECK(ok, "b=%u round-trip failed across the plane", b);
    }
}

static void exponent_rule_is_exact() {
    begin("exponent is the smallest that does not clip");
    for (uint8_t b = 2; b <= 16; ++b) {
        double mmax = biq::mantissa_max(b);
        // Exact powers of two are where a log2-based rule misrounds.
        for (int e = -20; e <= 20; ++e) {
            float peak = float(mmax * std::ldexp(1.0, e));
            int8_t got = biq::choose_exponent(peak, b);
            CHECK(got == e, "b=%u peak=m_max*2^%d chose e=%d", b, e, got);
            // A hair over must step up exactly one.
            float over = std::nextafterf(peak, 1e30f);
            CHECK(biq::choose_exponent(over, b) == e + 1, "b=%u peak just over 2^%d chose %d", b,
                  e, biq::choose_exponent(over, b));
        }
    }
    CHECK(biq::choose_exponent(0.0f, 6) == biq::kZeroBlock, "zero peak must be the sentinel");
}

static void header_round_trips_and_validates() {
    begin("header serialises, parses, and rejects malformed files");
    biq::Header h = biq::Header::make(5e6, 100e6, 6, 256);
    h.start_time_unix_ns = 1784106608000000000LL;
    h.ref_dbm_full_scale = -10.5;
    h.meta_bytes         = 24;

    uint8_t buf[biq::kHeaderBytes];
    h.serialize(buf);
    biq::Header r;
    std::string err;
    CHECK(r.parse(buf, sizeof(buf), &err), "valid header rejected: %s", err.c_str());
    CHECK(r.sample_rate_hz == 5e6 && r.center_freq_hz == 100e6, "scalars did not survive");
    CHECK(r.start_time_unix_ns == h.start_time_unix_ns, "timestamp did not survive");
    CHECK(r.ref_dbm_full_scale == -10.5, "ref level did not survive");
    CHECK(r.block_stride == 385, "stride = %u, want 385", r.block_stride);
    CHECK(r.data_start() == 64 + 24, "data_start = %zu", r.data_start());

    struct Bad {
        const char* what;
        size_t      off;
        uint8_t     val;
    } bad[] = {
        {"bad magic", 0, 'X'},
        {"channels != 1", 7, 2},
        {"mantissa_bits = 17", 8, 17},
        {"mantissa_bits = 1", 8, 1},
        {"block_size not a multiple of 8", 10, 100},
        {"unknown flag bit", 6, 0x02},
        {"reserved non-zero", 60, 1},
    };
    for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); ++i) {
        uint8_t b2[biq::kHeaderBytes];
        std::memcpy(b2, buf, sizeof(b2));
        b2[bad[i].off] = bad[i].val;
        biq::Header  x;
        std::string  e2;
        CHECK(!x.parse(b2, sizeof(b2), &e2), "%s was accepted", bad[i].what);
    }
    // Stride below the encoded block size must be caught too.
    uint8_t b3[biq::kHeaderBytes];
    std::memcpy(b3, buf, sizeof(b3));
    b3[12] = 10;
    b3[13] = 0;
    biq::Header x3;
    std::string e3;
    CHECK(!x3.parse(b3, sizeof(b3), &e3), "undersized block_stride was accepted");
    // meta_bytes must be 8-aligned.
    uint8_t b4[biq::kHeaderBytes];
    std::memcpy(b4, buf, sizeof(b4));
    b4[56] = 25;
    biq::Header x4;
    std::string e4;
    CHECK(!x4.parse(b4, sizeof(b4), &e4), "unaligned meta_bytes was accepted");
}

/// Build a test signal: a tone that steps down in level so the exponent moves.
static std::vector<float> test_signal(size_t n) {
    std::vector<float> x(n * 2);
    for (size_t i = 0; i < n; ++i) {
        double amp = 0.9 * std::pow(10.0, -3.0 * double(i) / double(n));  // 0 → -60 dB
        double ph  = 2.0 * M_PI * 0.0100073 * double(i);
        x[2 * i]     = float(amp * std::cos(ph));
        x[2 * i + 1] = float(amp * std::sin(ph));
    }
    return x;
}

static void file_round_trip_and_seek() {
    begin("file round-trips, seeks exactly, and reports honest geometry");
    const size_t N = 256 * 40 + 77;  // deliberately not a whole number of blocks
    std::vector<float> in = test_signal(N);
    std::string path = tmp_path("biq_test_rt.biq");
    std::string meta = "{\"core:author\":\"test\",\"core:hw\":\"synthetic\"}";

    biq::Header h = biq::Header::make(5e6, 100e6, 8, 256);
    h.ref_dbm_full_scale = -13.0;
    biq::Writer w;
    CHECK(w.open(path, h, meta), "open failed: %s", w.error().c_str());
    CHECK(w.write(in.data(), N), "write failed: %s", w.error().c_str());
    CHECK(w.close(), "close failed: %s", w.error().c_str());

    biq::Reader r;
    CHECK(r.open(path), "reopen failed: %s", r.error().c_str());
    CHECK(r.header().mantissa_bits == 8 && r.header().block_size == 256, "geometry lost");
    CHECK(r.header().sample_rate_hz == 5e6, "sample rate lost");
    CHECK(r.header().ref_dbm_full_scale == -13.0, "ref level lost");
    CHECK(r.metadata().compare(0, meta.size(), meta) == 0, "metadata lost: [%s]",
          r.metadata().c_str());
    CHECK(r.header().meta_bytes % 8 == 0, "metadata was not padded to 8");
    CHECK(!r.truncated(), "clean file reported as truncated");
    // The partial trailing block is zero-padded and kept, so we get 41 blocks.
    CHECK(r.blocks() == 41, "blocks = %llu, want 41", (unsigned long long)r.blocks());

    std::vector<float> out(r.samples() * 2);
    CHECK(r.read(out.data(), r.samples()) == r.samples(), "short read");
    double s = sqnr_db(in, out, 0, N * 2);
    CHECK(s > 37.0, "8-bit SQNR %.1f dB, expected > 37", s);

    // Seek to arbitrary (non-block-aligned) offsets and confirm the samples are
    // the same ones a linear read produced.
    const uint64_t probes[] = {0, 1, 255, 256, 257, 1000, 4095, 10239};
    for (size_t i = 0; i < sizeof(probes) / sizeof(probes[0]); ++i) {
        uint64_t p = probes[i];
        if (p >= r.samples()) continue;
        CHECK(r.seek(p), "seek(%llu) failed", (unsigned long long)p);
        CHECK(r.tell() == p, "tell = %llu after seek(%llu)", (unsigned long long)r.tell(),
              (unsigned long long)p);
        float got[8];
        size_t k = r.read(got, 4);
        CHECK(k == 4 || p + 4 > r.samples(), "seek+read short at %llu",
              (unsigned long long)p);
        bool same = true;
        for (size_t j = 0; j < k * 2; ++j)
            if (got[j] != out[p * 2 + j]) same = false;
        CHECK(same, "seek(%llu) gave different samples than a linear read",
              (unsigned long long)p);
    }
    CHECK(!r.seek(r.samples() + 1), "seek past end was allowed");

    // Reads that straddle blocks must be seamless.
    CHECK(r.seek(250), "seek failed");
    std::vector<float> straddle(20);
    CHECK(r.read(straddle.data(), 10) == 10, "straddling read short");
    bool same = true;
    for (size_t j = 0; j < 20; ++j)
        if (straddle[j] != out[250 * 2 + j]) same = false;
    CHECK(same, "read across a block boundary lost samples");

    std::remove(path.c_str());
}

static void envelope_is_free_and_correct() {
    begin("exponent envelope matches the decoded peaks");
    const size_t N = 256 * 32;
    std::vector<float> in = test_signal(N);
    std::string path = tmp_path("biq_test_env.biq");

    biq::Writer w;
    CHECK(w.open(path, biq::Header::make(5e6, 100e6, 6, 256)), "open failed");
    w.write(in.data(), N);
    w.close();

    biq::Reader r;
    CHECK(r.open(path), "reopen failed");
    std::vector<int8_t> env = r.envelope();
    CHECK(env.size() == r.blocks(), "envelope has %zu entries for %llu blocks", env.size(),
          (unsigned long long)r.blocks());

    // Each block's true peak must sit in (m_max*2^(e-1), m_max*2^e] — that is
    // exactly what the exponent promises, and what makes the envelope usable.
    double mmax = biq::mantissa_max(6);
    bool   ok = true;
    for (size_t k = 0; k < env.size(); ++k) {
        double peak = 0;
        for (size_t i = 0; i < 256 * 2; ++i) {
            size_t idx = k * 512 + i;
            if (idx < in.size()) peak = std::fmax(peak, std::fabs(in[idx]));
        }
        double hi = mmax * std::ldexp(1.0, env[k]);
        double lo = mmax * std::ldexp(1.0, env[k] - 1);
        if (!(peak <= hi * 1.0000001 && peak > lo * 0.9999999)) ok = false;
    }
    CHECK(ok, "an exponent did not bracket its block's peak");

    // Envelope must be monotonically decaying for a decaying signal.
    bool decays = true;
    for (size_t k = 1; k < env.size(); ++k)
        if (env[k] > env[k - 1]) decays = false;
    CHECK(decays, "envelope of a decaying signal went up");

    std::vector<float> db = r.envelope_dbfs();
    CHECK(db.size() == env.size(), "dBFS envelope size mismatch");
    CHECK(db.front() > db.back() + 40.0f, "expected >40 dB of decay, got %.1f",
          db.front() - db.back());

    std::remove(path.c_str());
}

static void truncation_is_recoverable() {
    begin("a truncated file still reads (SPEC §7)");
    const size_t N = 256 * 20;
    std::vector<float> in = test_signal(N);
    std::string path = tmp_path("biq_test_trunc.biq");

    biq::Writer w;
    CHECK(w.open(path, biq::Header::make(5e6, 100e6, 6, 256)), "open failed");
    w.write(in.data(), N);
    w.close();

    std::vector<uint8_t> whole = slurp(path);

    // Case 1: the recorder died — data_bytes is still 0 and the tail is a
    // partial block. Everything before it must still decode.
    {
        std::vector<uint8_t> t(whole.begin(), whole.end() - 500);
        for (int i = 0; i < 8; ++i) t[48 + i] = 0;  // data_bytes = 0
        std::FILE* f = std::fopen(path.c_str(), "wb");
        std::fwrite(t.data(), 1, t.size(), f);
        std::fclose(f);

        biq::Reader r;
        CHECK(r.open(path), "truncated file rejected: %s", r.error().c_str());
        CHECK(r.blocks() == 18, "recovered %llu blocks, want 18",
              (unsigned long long)r.blocks());
        std::vector<float> out(r.samples() * 2);
        CHECK(r.read(out.data(), r.samples()) == r.samples(), "short read on recovery");
        CHECK(sqnr_db(in, out, 0, r.samples() * 2) > 25.0, "recovered data is wrong");
    }

    // Case 2: data_bytes says more than the file holds (finalised, then the
    // copy was cut short). The reader must clamp, not read garbage.
    {
        std::vector<uint8_t> t(whole.begin(), whole.end() - 1000);
        std::FILE* f = std::fopen(path.c_str(), "wb");
        std::fwrite(t.data(), 1, t.size(), f);
        std::fclose(f);

        biq::Reader r;
        CHECK(r.open(path), "clamped file rejected");
        CHECK(r.truncated(), "truncation not reported");
        CHECK(r.blocks() == 17, "clamped to %llu blocks, want 17",
              (unsigned long long)r.blocks());
    }

    std::remove(path.c_str());
}

static void padded_stride_round_trips() {
    begin("padded block_stride round-trips");
    const size_t N = 256 * 8;
    std::vector<float> in = test_signal(N);
    std::string path = tmp_path("biq_test_pad.biq");

    biq::Header h = biq::Header::make(5e6, 100e6, 5, 256);
    CHECK(h.block_stride == 321, "exact stride = %u, want 321", h.block_stride);
    h.pad_stride(8);
    CHECK(h.block_stride == 328, "padded stride = %u, want 328", h.block_stride);

    biq::Writer w;
    CHECK(w.open(path, h), "open failed");
    w.write(in.data(), N);
    w.close();

    biq::Reader r;
    CHECK(r.open(path), "reopen failed: %s", r.error().c_str());
    CHECK(r.header().block_stride == 328, "stride lost");
    CHECK(r.samples() == N, "sample count lost");
    std::vector<float> out(N * 2);
    r.read(out.data(), N);
    CHECK(sqnr_db(in, out, 0, N * 2) > 20.0, "padded file decoded wrong");
    CHECK(r.seek(1000) && r.tell() == 1000, "seek broken with padded stride");

    std::remove(path.c_str());
}

static void real_mode_round_trips() {
    begin("real-valued mode round-trips");
    const size_t N = 256 * 8;
    std::vector<float> x(N);
    for (size_t i = 0; i < N; ++i) x[i] = float(0.8 * std::sin(2.0 * M_PI * 0.013 * double(i)));
    std::string path = tmp_path("biq_test_real.biq");

    biq::Header h = biq::Header::make(2e6, 0.0, 8, 256, /*real=*/true);
    CHECK(h.block_bytes() == 257, "real block = %zu bytes, want 257", h.block_bytes());
    biq::Writer w;
    CHECK(w.open(path, h), "open failed");
    w.write(x.data(), N);
    w.close();

    biq::Reader r;
    CHECK(r.open(path), "reopen failed");
    CHECK(r.header().real(), "real flag lost");
    std::vector<float> out(N);
    CHECK(r.read(out.data(), N) == N, "short read");
    double sig = 0, err = 0;
    for (size_t i = 0; i < N; ++i) {
        sig += double(x[i]) * x[i];
        err += (double(x[i]) - out[i]) * (double(x[i]) - out[i]);
    }
    CHECK(10.0 * std::log10(sig / err) > 37.0, "real-mode SQNR too low");

    std::remove(path.c_str());
}

static void profiles_hit_their_size_and_sqnr() {
    begin("profiles hit their advertised size and SQNR");
    const size_t N = 256 * 64;
    std::vector<float> in = test_signal(N);
    struct P {
        const char* name;
        uint8_t     bits;
        double      bytes_per_sample;
        double      sqnr_floor;
    } profiles[] = {
        {"survey ", 4, 1.00390625, 13.0},
        {"general", 6, 1.50390625, 25.0},
        {"archive", 8, 2.00390625, 37.0},
    };
    for (size_t i = 0; i < 3; ++i) {
        std::string path = tmp_path("biq_test_profile.biq");
        biq::Header h = biq::Header::make(5e6, 100e6, profiles[i].bits, 256);
        CHECK(std::fabs(h.bytes_per_sample() - profiles[i].bytes_per_sample) < 1e-9,
              "%s: %.6f bytes/sample, want %.6f", profiles[i].name, h.bytes_per_sample(),
              profiles[i].bytes_per_sample);

        biq::Writer w;
        w.open(path, h);
        w.write(in.data(), N);
        w.close();

        biq::Reader r;
        r.open(path);
        std::vector<float> out(N * 2);
        r.read(out.data(), N);
        double s = sqnr_db(in, out, 0, N * 2);
        CHECK(s > profiles[i].sqnr_floor, "%s: SQNR %.1f dB below floor %.1f", profiles[i].name,
              s, profiles[i].sqnr_floor);
        std::printf("  %s  b=%2u  %.4f bytes/sample  %.1f dB SQNR  %.2fx vs int16\n",
                    profiles[i].name, profiles[i].bits, h.bytes_per_sample(), s,
                    4.0 / h.bytes_per_sample());
        std::remove(path.c_str());
    }
}

// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
    if (argc > 1) g_dir = argv[1];

    golden_vectors_match_rust();
    packing_is_lsb_first();
    exponent_rule_is_exact();
    header_round_trips_and_validates();
    file_round_trip_and_seek();
    envelope_is_free_and_correct();
    truncation_is_recoverable();
    padded_stride_round_trips();
    real_mode_round_trips();
    profiles_hit_their_size_and_sqnr();

    std::printf("\n%d checks, %d failed\n", g_checks, g_fails);
    return g_fails == 0 ? 0 : 1;
}
