// biqconv — convert between .biq and the int16 IQ formats everyone already has.
//
//   biqconv in.wav out.biq -b 6          # WAV (or RF64) -> biq
//   biqconv in.biq out.wav               # back again, losslessly re-expanded
//   biqconv in.sigmf-data out.biq        # SigMF ci16_le + sidecar
//   biqconv in.cs16 out.biq -r 5e6 -f 100e6
//   biqconv info cap.biq                 # header + free envelope, decodes nothing
//
// The point of this tool is that an existing WAV-based workflow can try the
// format in one command and lose nothing but the bits it chose to drop.
//
// Formats are picked from the extension: .wav .cs16 .sigmf-data .biq

#include "biq.hpp"

#include <cctype>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <limits>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

static void die(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    std::fprintf(stderr, "biqconv: ");
    std::vfprintf(stderr, fmt, ap);
    std::fprintf(stderr, "\n");
    va_end(ap);
    std::exit(1);
}

static std::string ext_of(const std::string& p) {
    size_t slash = p.find_last_of('/');
    size_t dot   = p.find_last_of('.');
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) return "";
    std::string e = p.substr(dot + 1);
    for (size_t i = 0; i < e.size(); ++i) e[i] = char(std::tolower(e[i]));
    return e;
}

static std::string with_ext(const std::string& p, const char* e) {
    size_t dot = p.find_last_of('.');
    return (dot == std::string::npos ? p : p.substr(0, dot)) + e;
}

static std::string human(double bytes) {
    const char* u[] = {"B", "kB", "MB", "GB", "TB"};
    int         i   = 0;
    while (bytes >= 1000.0 && i < 4) {
        bytes /= 1000.0;
        ++i;
    }
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.2f %s", bytes, u[i]);
    return buf;
}

// ---------------------------------------------------------------------------
// UTC time
//
// By hand rather than with timegm(): that is not in C++11 and not on Windows.
// These are the standard civil<->days algorithms, exact over the whole range
// the header's int64 nanoseconds can hold.
// ---------------------------------------------------------------------------

static const int64_t kNsPerSec = 1000000000LL;

static int64_t days_from_civil(int y, int m, int d) {
    y -= m <= 2;
    int64_t era = (y >= 0 ? y : y - 399) / 400;
    int64_t yoe = y - era * 400;                                   // [0, 399]
    int64_t doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;  // [0, 365]
    int64_t doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;           // [0, 146096]
    return era * 146097 + doe - 719468;
}

static void civil_from_days(int64_t z, int* y, int* m, int* d) {
    z += 719468;
    int64_t era = (z >= 0 ? z : z - 146096) / 146097;
    int64_t doe = z - era * 146097;                                       // [0, 146096]
    int64_t yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;  // [0, 399]
    int64_t doy = doe - (365 * yoe + yoe / 4 - yoe / 100);                // [0, 365]
    int64_t mp  = (5 * doy + 2) / 153;                                    // [0, 11]
    *d          = int(doy - (153 * mp + 2) / 5 + 1);
    *m          = int(mp + (mp < 10 ? 3 : -9));
    *y          = int(yoe + era * 400 + (*m <= 2));
}

/// Civil UTC -> Unix nanoseconds. False on a date that isn't one.
static bool utc_to_ns(int y, int mo, int d, int h, int mi, int s, int64_t frac_ns, int64_t* out) {
    if (mo < 1 || mo > 12 || d < 1 || d > 31 || h > 23 || mi > 59 || s > 60) return false;
    if (y < 1970 || y > 2500) return false;
    *out = ((days_from_civil(y, mo, d) * 24 + h) * 60 + mi) * 60 * kNsPerSec +
           int64_t(s) * kNsPerSec + frac_ns;
    return true;
}

static void ns_to_civil(int64_t ns, int* y, int* mo, int* d, int* h, int* mi, int* s,
                        int64_t* frac) {
    int64_t secs = ns / kNsPerSec;
    *frac        = ns % kNsPerSec;
    if (*frac < 0) {  // floor, so pre-1970 stamps come out right too
        *frac += kNsPerSec;
        --secs;
    }
    int64_t days = secs / 86400, rem = secs % 86400;
    if (rem < 0) {
        rem += 86400;
        --days;
    }
    civil_from_days(days, y, mo, d);
    *h  = int(rem / 3600);
    *mi = int(rem / 60 % 60);
    *s  = int(rem % 60);
}

static std::string iso8601_from_ns(int64_t ns) {
    int     y, mo, d, h, mi, s;
    int64_t frac;
    ns_to_civil(ns, &y, &mo, &d, &h, &mi, &s, &frac);
    char buf[64];
    if (frac)
        std::snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%02d.%09lldZ", y, mo, d, h, mi,
                      s, (long long)frac);
    else
        std::snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%02dZ", y, mo, d, h, mi, s);
    return buf;
}

/// A zone suffix: ±HH:MM, ±HHMM or ±HH. Offset east of UTC, so UTC = clock - offset.
static bool parse_utc_offset(const char* p, int64_t* out) {
    if (*p != '+' && *p != '-') return false;
    int sign = (*p == '-') ? -1 : 1;
    ++p;
    int n = 0;
    while (std::isdigit(static_cast<unsigned char>(p[n]))) ++n;
    int oh = 0, om = 0;
    if (n == 2) {
        oh = (p[0] - '0') * 10 + (p[1] - '0');
        if (p[2] == ':' && std::isdigit(static_cast<unsigned char>(p[3])) &&
            std::isdigit(static_cast<unsigned char>(p[4])))
            om = (p[3] - '0') * 10 + (p[4] - '0');
        else if (p[2])
            return false;
    } else if (n == 4) {
        oh = (p[0] - '0') * 10 + (p[1] - '0');
        om = (p[2] - '0') * 10 + (p[3] - '0');
        if (p[4]) return false;
    } else {
        return false;
    }
    if (oh > 23 || om > 59) return false;
    *out = sign * (int64_t(oh) * 3600 + int64_t(om) * 60) * kNsPerSec;
    return true;
}

/// RFC 3339, which is also what SigMF's `core:datetime` is:
/// YYYY-MM-DDTHH:MM:SS[.frac][Z|±HH:MM]. A bare stamp is read as UTC.
static bool parse_iso8601(const std::string& t, int64_t* out) {
    int  y, mo, d, h, mi, s;
    char rest[32] = {0};
    if (std::sscanf(t.c_str(), "%4d-%2d-%2dT%2d:%2d:%2d%31s", &y, &mo, &d, &h, &mi, &s, rest) < 6 &&
        std::sscanf(t.c_str(), "%4d-%2d-%2d %2d:%2d:%2d%31s", &y, &mo, &d, &h, &mi, &s, rest) < 6)
        return false;

    const char* p       = rest;
    int64_t     frac_ns = 0;
    if (*p == '.') {
        ++p;
        int64_t scale = kNsPerSec;
        while (std::isdigit(static_cast<unsigned char>(*p))) {
            scale /= 10;
            frac_ns += (*p++ - '0') * scale;  // digits past nanoseconds fall off
        }
    }
    int64_t offset_ns = 0;
    if (*p == '+' || *p == '-') {
        if (!parse_utc_offset(p, &offset_ns)) return false;
    } else if (*p && *p != 'Z' && *p != 'z') {
        return false;
    }
    if (!utc_to_ns(y, mo, d, h, mi, s, frac_ns, out)) return false;
    *out -= offset_ns;
    return true;
}

/// How to read a recorder's zone-less wall clock.
///
/// `auxi` StartTime is a bare Windows SYSTEMTIME, and the SpectraVue convention
/// it follows is the operator's local time — SpectraFlux writes local, as does
/// most Windows SDR software. The .biq header is defined as UTC, so something
/// has to supply the zone. Defaulting to the converting machine's is right
/// whenever a file is converted where it was recorded, which is the ordinary
/// case; `--wav-tz` covers the rest, and the resolved UTC time is always
/// printed so a bad assumption shows up immediately instead of years later.
struct WallClock {
    enum Mode { kLocal, kUtc, kFixed };
    Mode    mode      = kLocal;
    int64_t offset_ns = 0;  // kFixed: UTC = clock - offset

    bool to_ns(int y, int mo, int d, int h, int mi, int s, int64_t frac_ns, int64_t* out) const {
        if (!utc_to_ns(y, mo, d, h, mi, s, frac_ns, out)) return false;
        if (mode == kUtc) return true;
        if (mode == kFixed) {
            *out -= offset_ns;
            return true;
        }
        // Local: mktime reads the fields in the host zone. tm_isdst = -1 lets it
        // work out DST; the hour that repeats at a fall-back boundary is
        // genuinely ambiguous and mktime's choice is as good as any.
        std::tm tm;
        std::memset(&tm, 0, sizeof(tm));
        tm.tm_year  = y - 1900;
        tm.tm_mon   = mo - 1;
        tm.tm_mday  = d;
        tm.tm_hour  = h;
        tm.tm_min   = mi;
        tm.tm_sec   = s;
        tm.tm_isdst = -1;
        std::time_t t = std::mktime(&tm);
        if (t == std::time_t(-1)) return false;
        *out = int64_t(t) * kNsPerSec + frac_ns;
        return true;
    }
};

/// `--wav-tz`: "local", "utc", or a fixed offset such as "+02:00" / "-0700".
static bool parse_wav_tz(const std::string& a, WallClock* wc) {
    if (a == "local") {
        wc->mode = WallClock::kLocal;
        return true;
    }
    if (a == "utc" || a == "UTC" || a == "z" || a == "Z") {
        wc->mode = WallClock::kUtc;
        return true;
    }
    if (a.size() >= 2 && (a[0] == '+' || a[0] == '-') && parse_utc_offset(a.c_str(), &wc->offset_ns)) {
        wc->mode = WallClock::kFixed;
        return true;
    }
    return false;
}

/// Center frequency out of a SpectraFlux/SDR-style filename, e.g.
/// "IQ_20260717_211008_100000000Hz_5000000sps_001.wav" -> 100e6. 0 if absent.
static double freq_from_name(const std::string& path, const char* suffix) {
    size_t n = std::strlen(suffix);
    for (size_t i = 0; i + n <= path.size(); ++i) {
        if (path.compare(i, n, suffix) != 0) continue;
        size_t j = i;
        while (j > 0 && std::isdigit(static_cast<unsigned char>(path[j - 1]))) --j;
        if (j < i) return std::atof(path.substr(j, i - j).c_str());
    }
    return 0.0;
}

// Minimal JSON scalar lookup. Enough for a SigMF sidecar, which is a flat object
// of known keys — not a general parser, and it does not pretend to be one.
static bool json_number(const std::string& js, const char* key, double* out) {
    std::string k = std::string("\"") + key + "\"";
    size_t      p = js.find(k);
    if (p == std::string::npos) return false;
    p = js.find(':', p + k.size());
    if (p == std::string::npos) return false;
    *out = std::atof(js.c_str() + p + 1);
    return true;
}
static bool json_string(const std::string& js, const char* key, std::string* out) {
    std::string k = std::string("\"") + key + "\"";
    size_t      p = js.find(k);
    if (p == std::string::npos) return false;
    p = js.find(':', p + k.size());
    if (p == std::string::npos) return false;
    size_t a = js.find('"', p);
    if (a == std::string::npos) return false;
    size_t b = js.find('"', a + 1);
    if (b == std::string::npos) return false;
    *out = js.substr(a + 1, b - a - 1);
    return true;
}

// ---------------------------------------------------------------------------
// int16 IQ sources: WAV / RF64, raw .cs16, SigMF
// ---------------------------------------------------------------------------

struct IntSource {
    std::FILE* f            = nullptr;
    uint64_t   data_offset  = 0;
    uint64_t   data_bytes   = 0;
    double     sample_rate  = 0;
    double     center_freq  = 0;
    int64_t    start_ns     = 0;  // 0 = unknown, per SPEC §2
    std::string datetime;
    std::string meta;  // verbatim sidecar JSON, if there was one

    uint64_t samples() const { return data_bytes / 4; }  // 2 ch * int16
    ~IntSource() {
        if (f) std::fclose(f);
    }
};

static uint16_t rd_u16(const uint8_t* p) {
    return uint16_t(uint16_t(p[0]) | uint16_t(p[1]) << 8);
}
static uint32_t rd_u32(const uint8_t* p) {
    return uint32_t(p[0]) | uint32_t(p[1]) << 8 | uint32_t(p[2]) << 16 | uint32_t(p[3]) << 24;
}
static uint64_t rd_u64(const uint8_t* p) {
    uint64_t v = 0;
    for (int i = 0; i < 8; ++i) v |= uint64_t(p[i]) << (8 * i);
    return v;
}

/// RIFF/WAVE and RF64/BW64. Wants 16-bit PCM stereo (= complex int16).
static void open_wav(IntSource& s, const std::string& path, const WallClock& wc) {
    s.f = std::fopen(path.c_str(), "rb");
    if (!s.f) die("cannot open %s", path.c_str());

    uint8_t riff[12];
    if (std::fread(riff, 1, 12, s.f) != 12) die("%s: too short", path.c_str());
    bool rf64 = std::memcmp(riff, "RF64", 4) == 0 || std::memcmp(riff, "BW64", 4) == 0;
    if (!rf64 && std::memcmp(riff, "RIFF", 4) != 0) die("%s: not a RIFF/RF64 file", path.c_str());
    if (std::memcmp(riff + 8, "WAVE", 4) != 0) die("%s: not WAVE", path.c_str());

    uint64_t ds64_data = 0;
    uint16_t channels = 0, bits = 0, fmt = 0;

    for (;;) {
        uint8_t ch[8];
        if (std::fread(ch, 1, 8, s.f) != 8) break;
        std::string id(reinterpret_cast<char*>(ch), 4);
        uint64_t    sz  = rd_u32(ch + 4);
        long        pos = std::ftell(s.f);

        if (id == "ds64") {
            uint8_t b[24];
            if (std::fread(b, 1, 24, s.f) == 24) ds64_data = rd_u64(b + 8);
        } else if (id == "fmt ") {
            std::vector<uint8_t> b(sz < 16 ? 16 : size_t(sz), 0);
            std::fread(b.data(), 1, size_t(sz), s.f);
            fmt         = uint16_t(b[0] | b[1] << 8);
            channels    = uint16_t(b[2] | b[3] << 8);
            s.sample_rate = rd_u32(&b[4]);
            bits        = uint16_t(b[14] | b[15] << 8);
        } else if (id == "auxi") {
            // SpectraVue `auxi`: two SYSTEMTIMEs (32 bytes) then DWORD CenterFreq.
            std::vector<uint8_t> b(size_t(sz), 0);
            std::fread(b.data(), 1, size_t(sz), s.f);
            if (sz >= 36) {
                uint32_t cf = rd_u32(&b[32]);
                if (cf > 0) s.center_freq = cf;
            }
            // StartTime is a SYSTEMTIME: wYear, wMonth, wDayOfWeek, wDay,
            // wHour, wMinute, wSecond, wMilliseconds. wDayOfWeek is redundant
            // and some writers leave it zero, so it is ignored.
            if (sz >= 16) {
                int     y = rd_u16(&b[0]), mo = rd_u16(&b[2]), d = rd_u16(&b[6]);
                int     h = rd_u16(&b[8]), mi = rd_u16(&b[10]), sec = rd_u16(&b[12]);
                int64_t ms = rd_u16(&b[14]), ns;
                if (y && wc.to_ns(y, mo, d, h, mi, sec, ms * 1000000, &ns)) s.start_ns = ns;
            }
        } else if (id == "data") {
            s.data_offset = uint64_t(pos);
            s.data_bytes  = (rf64 && sz == 0xFFFFFFFFu) ? ds64_data : sz;
            break;
        }
        if (std::fseek(s.f, pos + long(sz) + long(sz & 1), SEEK_SET) != 0) break;
    }

    if (!s.data_bytes) die("%s: no data chunk", path.c_str());
    if (fmt != 1) die("%s: not PCM (format tag %u)", path.c_str(), fmt);
    if (channels != 2) die("%s: %u channels, expected 2 (I and Q)", path.c_str(), channels);
    if (bits != 16) die("%s: %u bits/sample, expected 16", path.c_str(), bits);

    // Clamp a declared size that runs past EOF — a WAV cut short by a crash.
    std::fseek(s.f, 0, SEEK_END);
    uint64_t avail = uint64_t(std::ftell(s.f)) - s.data_offset;
    if (s.data_bytes > avail) s.data_bytes = avail;
    std::fseek(s.f, long(s.data_offset), SEEK_SET);

    if (s.center_freq == 0) s.center_freq = freq_from_name(path, "Hz");
}

static void open_cs16(IntSource& s, const std::string& path, double rate, double freq) {
    s.f = std::fopen(path.c_str(), "rb");
    if (!s.f) die("cannot open %s", path.c_str());
    std::fseek(s.f, 0, SEEK_END);
    s.data_bytes  = uint64_t(std::ftell(s.f));
    s.data_offset = 0;
    std::fseek(s.f, 0, SEEK_SET);
    s.sample_rate = rate > 0 ? rate : freq_from_name(path, "sps");
    s.center_freq = freq > 0 ? freq : freq_from_name(path, "Hz");
    if (s.sample_rate <= 0)
        die("%s is headerless — pass -r <sample rate>", path.c_str());
}

static void open_sigmf(IntSource& s, const std::string& path) {
    std::string meta_path = with_ext(path, ".sigmf-meta");
    std::FILE*  m         = std::fopen(meta_path.c_str(), "rb");
    if (!m) die("cannot open sidecar %s", meta_path.c_str());
    std::fseek(m, 0, SEEK_END);
    s.meta.resize(size_t(std::ftell(m)));
    std::fseek(m, 0, SEEK_SET);
    if (!s.meta.empty()) std::fread(&s.meta[0], 1, s.meta.size(), m);
    std::fclose(m);

    std::string dt;
    if (json_string(s.meta, "core:datatype", &dt) && dt != "ci16_le" && dt != "ci16")
        die("%s: datatype %s — only ci16_le is supported", meta_path.c_str(), dt.c_str());
    json_number(s.meta, "core:sample_rate", &s.sample_rate);
    json_number(s.meta, "core:frequency", &s.center_freq);
    if (json_string(s.meta, "core:datetime", &s.datetime) && !s.datetime.empty() &&
        !parse_iso8601(s.datetime, &s.start_ns))
        std::fprintf(stderr,
                     "biqconv: warning: %s: core:datetime \"%s\" is not RFC 3339; the .biq "
                     "start time will be left unknown\n",
                     meta_path.c_str(), s.datetime.c_str());

    s.f = std::fopen(path.c_str(), "rb");
    if (!s.f) die("cannot open %s", path.c_str());
    std::fseek(s.f, 0, SEEK_END);
    s.data_bytes  = uint64_t(std::ftell(s.f));
    s.data_offset = 0;
    std::fseek(s.f, 0, SEEK_SET);
}

// ---------------------------------------------------------------------------
// int16 IQ sinks
// ---------------------------------------------------------------------------

static void put_u32(std::FILE* f, uint32_t v) {
    uint8_t b[4];
    for (int i = 0; i < 4; ++i) b[i] = uint8_t(v >> (8 * i));
    std::fwrite(b, 1, 4, f);
}
static void put_u16(std::FILE* f, uint16_t v) {
    uint8_t b[2] = {uint8_t(v), uint8_t(v >> 8)};
    std::fwrite(b, 1, 2, f);
}

/// Fill a Windows SYSTEMTIME at `out` from Unix nanoseconds, for `auxi`.
/// Written as UTC: the reader's default is local, so a `biq -> wav -> biq`
/// round trip only returns the same instant when both ends agree. `--wav-tz
/// utc` on the way back is what closes that loop.
static void put_systemtime(uint8_t* out, int64_t ns) {
    int     y, mo, d, h, mi, s;
    int64_t frac;
    ns_to_civil(ns, &y, &mo, &d, &h, &mi, &s, &frac);
    int64_t days = days_from_civil(y, mo, d);
    int     dow  = int(((days % 7) + 11) % 7);  // 1970-01-01 was a Thursday; Sunday = 0
    const uint16_t f[8] = {uint16_t(y),  uint16_t(mo), uint16_t(dow),
                           uint16_t(d),  uint16_t(h),  uint16_t(mi),
                           uint16_t(s),  uint16_t(frac / 1000000)};
    for (int i = 0; i < 8; ++i) {
        out[i * 2]     = uint8_t(f[i]);
        out[i * 2 + 1] = uint8_t(f[i] >> 8);
    }
}

/// 16-bit stereo WAV with a SpectraVue `auxi` chunk so the centre frequency and
/// start time survive the trip. Sizes are backfilled at close.
static std::FILE* open_wav_out(const std::string& path, double rate, double freq, int64_t start_ns,
                               uint64_t samples) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) die("cannot open %s for writing", path.c_str());
    std::fwrite("RIFF", 1, 4, f);
    put_u32(f, 0);  // backfilled
    std::fwrite("WAVE", 1, 4, f);

    std::fwrite("fmt ", 1, 4, f);
    put_u32(f, 16);
    put_u16(f, 1);                              // PCM
    put_u16(f, 2);                              // I and Q
    put_u32(f, uint32_t(rate));
    put_u32(f, uint32_t(rate) * 4);             // byte rate
    put_u16(f, 4);                              // block align
    put_u16(f, 16);                             // bits

    std::fwrite("auxi", 1, 4, f);
    put_u32(f, 64);
    uint8_t auxi[64];
    std::memset(auxi, 0, sizeof(auxi));
    if (start_ns) {
        put_systemtime(auxi, start_ns);
        double dur = rate > 0 ? double(samples) / rate : 0.0;
        put_systemtime(auxi + 16, start_ns + int64_t(dur * 1e9));
    }
    uint32_t cf = uint32_t(freq < 0 ? 0 : (freq > 4294967295.0 ? 0 : freq));
    for (int i = 0; i < 4; ++i) auxi[32 + i] = uint8_t(cf >> (8 * i));
    std::fwrite(auxi, 1, sizeof(auxi), f);

    std::fwrite("data", 1, 4, f);
    put_u32(f, 0);  // backfilled
    return f;
}

static void close_wav_out(std::FILE* f, uint64_t data_bytes) {
    if (data_bytes & 1) std::fputc(0, f);
    long total = std::ftell(f);
    std::fseek(f, 4, SEEK_SET);
    put_u32(f, uint32_t(total - 8));
    std::fseek(f, 12 + 8 + 16 + 8 + 64 + 4, SEEK_SET);  // the data chunk's size field
    put_u32(f, uint32_t(data_bytes));
    std::fclose(f);
    if (data_bytes > 0xFFFFFFFFull)
        std::fprintf(stderr,
                     "biqconv: warning: %llu bytes exceeds WAV's 4 GB limit; the size field "
                     "wrapped. Use RF64-aware tools or keep the .biq.\n",
                     (unsigned long long)data_bytes);
}

static void write_sigmf_meta(const std::string& data_path, const biq::Header& h,
                             const std::string& biq_meta) {
    std::string p = with_ext(data_path, ".sigmf-meta");
    std::FILE*  f = std::fopen(p.c_str(), "wb");
    if (!f) die("cannot write %s", p.c_str());
    std::string dt;
    json_string(biq_meta, "core:datetime", &dt);
    std::fprintf(f,
                 "{\n"
                 "  \"global\": {\n"
                 "    \"core:datatype\": \"ci16_le\",\n"
                 "    \"core:sample_rate\": %.10g,\n"
                 "    \"core:version\": \"1.0.0\",\n"
                 "    \"core:recorder\": \"biqconv\"\n"
                 "  },\n"
                 "  \"captures\": [\n"
                 "    {\n"
                 "      \"core:sample_start\": 0,\n"
                 "      \"core:frequency\": %.10g%s%s%s\n"
                 "    }\n"
                 "  ],\n"
                 "  \"annotations\": []\n"
                 "}\n",
                 h.sample_rate_hz, h.center_freq_hz, dt.empty() ? "" : ",\n      \"core:datetime\": \"",
                 dt.empty() ? "" : dt.c_str(), dt.empty() ? "" : "\"");
    std::fclose(f);
}

// ---------------------------------------------------------------------------
// Conversions
// ---------------------------------------------------------------------------

static const size_t kChunk = 1 << 16;  // samples per I/O chunk

static void to_biq(const std::string& in, const std::string& out, uint8_t bits,
                   uint16_t block, size_t align, double rate, double freq, int64_t start_ns,
                   const WallClock& wc) {
    IntSource   s;
    std::string e = ext_of(in);
    if (e == "wav" || e == "rf64" || e == "bw64")
        open_wav(s, in, wc);
    else if (e == "sigmf-data")
        open_sigmf(s, in);
    else
        open_cs16(s, in, rate, freq);

    if (rate > 0) s.sample_rate = rate;
    if (freq > 0) s.center_freq = freq;
    if (start_ns) s.start_ns = start_ns;

    biq::Header h = biq::Header::make(s.sample_rate, s.center_freq, bits, block);
    if (align > 1) h.pad_stride(align);
    h.start_time_unix_ns = s.start_ns;

    // core:datetime is regenerated from the header field rather than passed
    // through, so the JSON and byte 32 can never disagree.
    std::string meta = "{\"core:sample_rate\":" + std::to_string(s.sample_rate) +
                       ",\"core:frequency\":" + std::to_string(s.center_freq);
    if (s.start_ns) meta += ",\"core:datetime\":\"" + iso8601_from_ns(s.start_ns) + "\"";
    meta += ",\"core:recorder\":\"biqconv\"}";

    biq::Writer w;
    if (!w.open(out, h, meta)) die("%s", w.error().c_str());

    std::vector<int16_t> ibuf(kChunk * 2);
    std::vector<float>   fbuf(kChunk * 2);
    uint64_t             left = s.samples(), done = 0;
    while (left) {
        size_t n = size_t(left < kChunk ? left : kChunk);
        if (std::fread(ibuf.data(), 4, n, s.f) != n) break;
        for (size_t i = 0; i < n * 2; ++i) fbuf[i] = float(ibuf[i]) * (1.0f / 32768.0f);
        if (!w.write(fbuf.data(), n)) die("%s", w.error().c_str());
        left -= n;
        done += n;
    }
    if (!w.close()) die("%s", w.error().c_str());

    double src = double(done) * 4.0, dst = double(done) * w.header().bytes_per_sample();
    std::printf("%s -> %s\n", in.c_str(), out.c_str());
    std::printf("  %llu samples  %.6g MHz  %.6g Msps  b=%u N=%u\n", (unsigned long long)done,
                s.center_freq / 1e6, s.sample_rate / 1e6, bits, block);
    std::printf("  %s -> %s   %.2fx smaller\n", human(src).c_str(), human(dst).c_str(),
                src / dst);
    std::printf("  start time %s\n",
                s.start_ns ? iso8601_from_ns(s.start_ns).c_str() : "unknown");
}

static void from_biq(const std::string& in, const std::string& out) {
    biq::Reader r;
    if (!r.open(in)) die("%s: %s", in.c_str(), r.error().c_str());
    if (r.header().real()) die("%s holds real samples; int16 IQ output expects complex", in.c_str());
    if (r.truncated())
        std::fprintf(stderr, "biqconv: warning: %s is truncated; converting what is there\n",
                     in.c_str());

    std::string e = ext_of(out);
    bool        wav = (e == "wav");
    std::FILE*  f;
    if (wav)
        f = open_wav_out(out, r.header().sample_rate_hz, r.header().center_freq_hz,
                         r.header().start_time_unix_ns, r.samples());
    else if ((f = std::fopen(out.c_str(), "wb")) == nullptr)
        die("cannot open %s for writing", out.c_str());
    if (e == "sigmf-data") write_sigmf_meta(out, r.header(), r.metadata());

    std::vector<float>   fbuf(kChunk * 2);
    std::vector<int16_t> ibuf(kChunk * 2);
    uint64_t             done = 0, bytes = 0;
    for (;;) {
        size_t n = r.read(fbuf.data(), kChunk);
        if (!n) break;
        for (size_t i = 0; i < n * 2; ++i) {
            float v = fbuf[i] * 32768.0f;
            if (v > 32767.0f) v = 32767.0f;
            if (v < -32768.0f) v = -32768.0f;
            ibuf[i] = int16_t(std::lrintf(v));
        }
        std::fwrite(ibuf.data(), 4, n, f);
        done += n;
        bytes += n * 4;
    }
    if (wav)
        close_wav_out(f, bytes);
    else
        std::fclose(f);

    std::printf("%s -> %s\n", in.c_str(), out.c_str());
    std::printf("  %llu samples  %.6g MHz  %.6g Msps  b=%u\n", (unsigned long long)done,
                r.header().center_freq_hz / 1e6, r.header().sample_rate_hz / 1e6,
                r.header().mantissa_bits);
    std::printf("  %s -> %s\n",
                human(double(done) * r.header().bytes_per_sample()).c_str(),
                human(double(bytes)).c_str());
}

static void info(const std::string& path) {
    biq::Reader r;
    if (!r.open(path)) die("%s: %s", path.c_str(), r.error().c_str());
    const biq::Header& h = r.header();

    std::printf("%s\n", path.c_str());
    std::printf("  geometry      b=%u  N=%u  stride=%u%s  %s\n", h.mantissa_bits, h.block_size,
                h.block_stride, h.block_stride == h.block_bytes() ? " (exact)" : " (padded)",
                h.real() ? "real" : "complex");
    std::printf("  sample rate   %.10g Hz\n", h.sample_rate_hz);
    std::printf("  centre        %.10g Hz\n", h.center_freq_hz);
    if (h.start_time_unix_ns)
        std::printf("  start time    %s  (%lld ns)\n", iso8601_from_ns(h.start_time_unix_ns).c_str(),
                    (long long)h.start_time_unix_ns);
    else
        std::printf("  start time    unknown\n");
    if (h.ref_dbm_full_scale == h.ref_dbm_full_scale)
        std::printf("  full scale    %.2f dBm\n", h.ref_dbm_full_scale);
    else
        std::printf("  full scale    uncalibrated\n");
    std::printf("  samples       %llu in %llu blocks  (%.3f s)\n",
                (unsigned long long)r.samples(), (unsigned long long)r.blocks(),
                r.duration_s());
    std::printf("  size          %s  (%.4f bytes/sample, %.2fx vs int16)%s\n",
                human(double(r.blocks()) * h.block_stride).c_str(), h.bytes_per_sample(),
                4.0 / h.bytes_per_sample(), r.truncated() ? "  TRUNCATED" : "");
    if (h.meta_bytes) std::printf("  metadata      %s\n", r.metadata().c_str());

    // The free envelope: one strided pass, no mantissas decoded.
    std::vector<float> db = r.envelope_dbfs();
    if (db.empty()) return;
    const int W = 64;
    float lo = 1e30f, hi = -1e30f;
    for (size_t i = 0; i < db.size(); ++i) {
        if (db[i] == -std::numeric_limits<float>::infinity()) continue;
        if (db[i] < lo) lo = db[i];
        if (db[i] > hi) hi = db[i];
    }
    if (lo > hi) return;
    static const char* bars[] = {" ", "▁", "▂", "▃", "▄", "▅", "▆", "▇", "█"};
    std::printf("  envelope      %.0f dBFS ", hi);
    for (int c = 0; c < W; ++c) {
        size_t a = db.size() * size_t(c) / W, b = db.size() * size_t(c + 1) / W;
        if (b <= a) b = a + 1;
        float peak = -1e30f;
        for (size_t i = a; i < b && i < db.size(); ++i)
            if (db[i] > peak) peak = db[i];
        int level = (hi > lo) ? int((peak - lo) / (hi - lo) * 8.0f + 0.5f) : 8;
        if (level < 0) level = 0;
        if (level > 8) level = 8;
        std::printf("%s", bars[level]);
    }
    std::printf(" %.0f dBFS\n", lo);
    std::printf("                %llu exponent bytes read, zero mantissas decoded\n",
                (unsigned long long)r.blocks());
}

// ---------------------------------------------------------------------------

static void usage() {
    std::printf(
        "biqconv — convert between .biq and int16 IQ formats\n"
        "\n"
        "  biqconv <in> <out> [options]\n"
        "  biqconv info <file.biq>\n"
        "\n"
        "Formats are chosen by extension:\n"
        "  .wav .rf64 .bw64   16-bit stereo PCM, SpectraVue `auxi` centre frequency\n"
        "  .sigmf-data        SigMF ci16_le with a .sigmf-meta sidecar\n"
        "  .cs16              headerless interleaved int16 (GNU Radio convention)\n"
        "  .biq               block floating point\n"
        "\n"
        "Options (converting to .biq):\n"
        "  -b <bits>   mantissa width, 2..16          (default 6)\n"
        "  -N <n>      block size, multiple of 8      (default 256)\n"
        "  --pad <n>   round the block stride up to n bytes for mmap readers\n"
        "  -r <hz>     sample rate   (required for headerless input)\n"
        "  -f <hz>     centre frequency\n"
        "  -t <time>   UTC start time, RFC 3339: \"2026-07-26T19:14:42Z\"\n"
        "  --wav-tz <z>  zone of the WAV `auxi` clock: local (default), utc,\n"
        "              or a fixed offset like -04:00\n"
        "\n"
        "Start time comes from the `auxi` chunk of a WAV or a SigMF\n"
        "`core:datetime`; -t overrides either. Without one the header records\n"
        "\"unknown\", which is legal. `auxi` carries no zone and SpectraVue-style\n"
        "recorders write local time, so it is read in the host's zone unless\n"
        "--wav-tz says otherwise; the resolved UTC time is always printed.\n"
        "\n"
        "Profiles: -b 4 survey (4.0x), -b 6 general (2.7x), -b 8 archive (2.0x)\n");
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage();
        return 1;
    }
    if (std::strcmp(argv[1], "info") == 0) {
        if (argc < 3) die("info needs a file");
        info(argv[2]);
        return 0;
    }
    if (std::strcmp(argv[1], "-h") == 0 || std::strcmp(argv[1], "--help") == 0) {
        usage();
        return 0;
    }
    if (argc < 3) {
        usage();
        return 1;
    }

    std::string in = argv[1], out = argv[2];
    uint8_t     bits  = 6;
    uint16_t    block = 256;
    size_t      align = 1;
    double      rate = 0, freq = 0;
    int64_t     start_ns = 0;
    WallClock   wc;
    for (int i = 3; i < argc; ++i) {
        std::string a = argv[i];
        const char* v = (i + 1 < argc) ? argv[i + 1] : nullptr;
        if (a == "-b" && v)
            bits = uint8_t(std::atoi(argv[++i]));
        else if (a == "-N" && v)
            block = uint16_t(std::atoi(argv[++i]));
        else if (a == "--pad" && v)
            align = size_t(std::atoi(argv[++i]));
        else if (a == "-r" && v)
            rate = std::atof(argv[++i]);
        else if (a == "-f" && v)
            freq = std::atof(argv[++i]);
        else if (a == "-t" && v) {
            if (!parse_iso8601(argv[++i], &start_ns))
                die("-t wants RFC 3339, e.g. 2026-07-26T19:14:42Z (got %s)", argv[i]);
            if (start_ns == 0) die("-t 1970-01-01T00:00:00Z is how a .biq header spells \"unknown\"");
        } else if (a == "--wav-tz" && v) {
            if (!parse_wav_tz(argv[++i], &wc))
                die("--wav-tz wants local, utc, or an offset like -04:00 (got %s)", argv[i]);
        } else
            die("unknown option %s", a.c_str());
    }
    if (bits < 2 || bits > 16) die("-b must be 2..16");
    if (block == 0 || block % 8) die("-N must be a non-zero multiple of 8");

    if (ext_of(in) == "biq")
        from_biq(in, out);
    else if (ext_of(out) == "biq")
        to_biq(in, out, bits, block, align, rate, freq, start_ns, wc);
    else
        die("one side must be a .biq file");
    return 0;
}
