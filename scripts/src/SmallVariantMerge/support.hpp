#pragma once
#include <algorithm>
#include <array>
#include <charconv>
#include <cctype>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <zlib.h>

namespace sm {
using U32 = uint32_t;
using U64 = uint64_t;
using Str = std::string;
using View = std::string_view;
namespace fs = std::filesystem;
inline void require(bool ok, const char *what) {
    if (!ok)
        throw std::runtime_error(what);
}
inline void require(bool ok, const Str &what) {
    if (!ok)
        throw std::runtime_error(what);
}
inline U32 index32(U64 n) {
    require(n <= UINT32_MAX, "coordinate/index outside uint32");
    return U32(n);
}
inline int64_t integer(View s) {
    int64_t n = 0;
    if (!s.empty() && s.front() == '+') {
        s.remove_prefix(1);
        require(!s.empty() && s.front() >= '0' && s.front() <= '9', "invalid integer sign");
    }
    require(!s.empty(), "invalid empty integer");
    auto r = std::from_chars(s.data(), s.data() + s.size(), n);
    if (r.ec != std::errc() || r.ptr != s.data() + s.size())
        throw std::runtime_error("invalid integer: " + Str(s));
    return n;
}
inline U32 coordinate(View s) {
    auto n = integer(s);
    require(n >= 0, "negative coordinate/index");
    return index32(n);
}
inline std::vector<View> split(View s, char sep) {
    std::vector<View> out;
    size_t start = 0;
    for (;;) {
        auto end = s.find(sep, start);
        if (end == View::npos) {
            out.push_back(s.substr(start));
            return out;
        }
        out.push_back(s.substr(start, end - start));
        start = end + 1;
    }
}
inline Str join(const std::vector<Str> &a, View sep) {
    Str s;
    for (size_t i = 0; i < a.size(); ++i) {
        if (i)
            s += sep;
        s += a[i];
    }
    return s;
}
inline Str replace(Str s, View from, View to) {
    size_t p = 0;
    while ((p = s.find(from, p)) != Str::npos) {
        s.replace(p, from.size(), to);
        p += to.size();
    }
    return s;
}
inline Str unescape(View s) {
    if (s == ".")
        return "";
    Str r(s);
    for (auto pair : {std::pair<View, View>{"@0A", "\n"},
                      {"@09", "\t"},
                      {"@2C", ","},
                      {"@3B", ";"},
                      {"@3A", ":"},
                      {"@40", "@"}})
        r = replace(std::move(r), pair.first, pair.second);
    return r;
}
inline void append_escaped(Str &out, View s) {
    if (s.empty()) {
        out += '.';
        return;
    }
    for (char c : s) {
        switch (c) {
        case '@':
            out += "@40";
            break;
        case ':':
            out += "@3A";
            break;
        case ';':
            out += "@3B";
            break;
        case ',':
            out += "@2C";
            break;
        case '\t':
            out += "@09";
            break;
        case '\n':
            out += "@0A";
            break;
        default:
            out += c;
        }
    }
}
inline Str escape(View s) {
    Str out;
    append_escaped(out, s);
    return out;
}
inline Str upper(Str s) {
    for (auto &c : s)
        if (c >= 'a' && c <= 'z')
            c -= 32;
    return s;
}
inline Str quote(View s) {
    // Python json.dumps(..., ensure_ascii=True), including non-BMP names.
    Str out = "\"";
    auto hex = [&](U32 n) {
        const char *h = "0123456789abcdef";
        out += "\\u";
        for (int shift = 12; shift >= 0; shift -= 4)
            out += h[(n >> shift) & 15];
    };
    for (size_t i = 0; i < s.size(); ++i) {
        unsigned char c = s[i];
        if (c == '"' || c == '\\') {
            out += '\\';
            out += char(c);
        } else if (c == '\n')
            out += "\\n";
        else if (c == '\r')
            out += "\\r";
        else if (c == '\t')
            out += "\\t";
        else if (c == '\b')
            out += "\\b";
        else if (c == '\f')
            out += "\\f";
        else if (c < 32)
            hex(c);
        else if (c < 127)
            out += char(c);
        else if (c == 127)
            hex(c);
        else {
            int extra = c < 224 ? 1 : c < 240 ? 2 : 3;
            U32 cp = c & ((1u << (6 - extra)) - 1);
            for (int j = 0; j < extra; ++j) {
                require(++i < s.size() && (uint8_t(s[i]) & 192) == 128, "invalid UTF-8");
                cp = (cp << 6) | (uint8_t(s[i]) & 63);
            }
            if (cp <= 65535)
                hex(cp);
            else {
                cp -= 65536;
                hex(0xd800 + (cp >> 10));
                hex(0xdc00 + (cp & 1023));
            }
        }
    }
    return out + '"';
}
inline Str json_strings(const std::vector<Str> &values) {
    Str out = "[";
    for (size_t i = 0; i < values.size(); ++i) {
        if (i)
            out += ", ";
        out += quote(values[i]);
    }
    return out + ']';
}
inline Str sha256(View text) {
    static constexpr U32 k[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
    std::array<U32, 8> h{0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                         0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    Str data(text);
    data += char(128);
    while (data.size() % 64 != 56)
        data += '\0';
    U64 bits = U64(text.size()) * 8;
    for (int i = 56; i >= 0; i -= 8)
        data += char(bits >> i);
    auto rr = [](U32 n, int b) { return (n >> b) | (n << (32 - b)); };
    for (size_t offset = 0; offset < data.size(); offset += 64) {
        U32 w[64];
        for (int i = 0; i < 16; ++i) {
            w[i] = 0;
            for (int j = 0; j < 4; ++j)
                w[i] = (w[i] << 8) | uint8_t(data[offset + i * 4 + j]);
        }
        for (int i = 16; i < 64; ++i)
            w[i] = w[i - 16] + (rr(w[i - 15], 7) ^ rr(w[i - 15], 18) ^ (w[i - 15] >> 3)) + w[i - 7] +
                   (rr(w[i - 2], 17) ^ rr(w[i - 2], 19) ^ (w[i - 2] >> 10));
        auto a = h;
        for (int i = 0; i < 64; ++i) {
            U32 t1 = a[7] + (rr(a[4], 6) ^ rr(a[4], 11) ^ rr(a[4], 25)) + ((a[4] & a[5]) ^ (~a[4] & a[6])) +
                     k[i] + w[i];
            U32 t2 =
                (rr(a[0], 2) ^ rr(a[0], 13) ^ rr(a[0], 22)) + ((a[0] & a[1]) ^ (a[0] & a[2]) ^ (a[1] & a[2]));
            a = {t1 + t2, a[0], a[1], a[2], a[3] + t1, a[4], a[5], a[6]};
        }
        for (int i = 0; i < 8; ++i)
            h[i] += a[i];
    }
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (auto v : h)
        out << std::setw(8) << v;
    return out.str();
}
struct Reader {
    std::ifstream f;
    explicit Reader(const Str &path) : f(path, std::ios::binary) { require(bool(f), "cannot open " + path); }
    void bytes(char *p, size_t n) {
        f.read(p, n);
        require(bool(f), "truncated native merge request/shard");
    }
    U32 u32() {
        unsigned char b[4];
        bytes(reinterpret_cast<char *>(b), 4);
        return U32(b[0]) | (U32(b[1]) << 8) | (U32(b[2]) << 16) | (U32(b[3]) << 24);
    }
    U64 u64() {
        U64 lo = u32();
        return lo | (U64(u32()) << 32);
    }
    Str str() {
        U32 n = u32();
        Str s(n, '\0');
        bytes(s.data(), n);
        return s;
    }
    std::vector<Str> strings() {
        U32 n = u32();
        std::vector<Str> a;
        a.reserve(n);
        while (n--)
            a.push_back(str());
        return a;
    }
};
inline void put32(std::ostream &f, U32 n) {
    char b[4];
    for (int i = 0; i < 4; ++i)
        b[i] = char(n >> (i * 8));
    f.write(b, 4);
}
inline void putstr(std::ostream &f, View s) {
    put32(f, index32(s.size()));
    f.write(s.data(), s.size());
}
struct Lines {
    gzFile f = nullptr;
    std::array<char, 1 << 16> buffer{};
    size_t pos = 0, end = 0;
    explicit Lines(const Str &path) {
        f = gzopen(path.c_str(), "rb");
        require(f, "cannot read " + path);
        gzbuffer(f, 1 << 20);
    }
    ~Lines() {
        if (f)
            gzclose(f);
    }
    bool next(Str &line) {
        line.clear();
        for (;;) {
            if (pos == end) {
                int n = gzread(f, buffer.data(), buffer.size());
                require(n >= 0, "compressed input read failed");
                if (!n) {
                    int error = Z_OK;
                    gzerror(f, &error);
                    require(error == Z_OK || error == Z_STREAM_END, "truncated compressed input");
                    if (!line.empty() && line.back() == '\r')
                        line.pop_back();
                    return !line.empty();
                }
                pos = 0;
                end = n;
            }
            auto p = static_cast<const char *>(memchr(buffer.data() + pos, '\n', end - pos));
            size_t stop = p ? size_t(p - buffer.data()) : end;
            line.append(buffer.data() + pos, stop - pos);
            pos = stop + (p ? 1 : 0);
            if (p) {
                if (!line.empty() && line.back() == '\r')
                    line.pop_back();
                return true;
            }
        }
    }
};
struct Intern {
    std::vector<Str> values;
    std::unordered_map<Str, U32> ids;
    U32 add(View s) {
        auto it = ids.find(Str(s));
        if (it != ids.end())
            return it->second;
        U32 id = index32(values.size());
        values.emplace_back(s);
        ids.emplace(values.back(), id);
        return id;
    }
};
inline void checked_close(std::ofstream &f) {
    f.close();
    require(!f.fail(), "failed to write merge output (disk full?)");
}
inline std::unordered_map<Str, Str> info(View s) {
    std::unordered_map<Str, Str> out;
    for (auto p : split(s, ';'))
        if (!p.empty() && p != ".") {
            auto i = p.find('=');
            out[Str(p.substr(0, i))] = i == View::npos ? "" : Str(p.substr(i + 1));
        }
    return out;
}
inline Str get(const std::unordered_map<Str, Str> &m, const Str &key, const Str &fallback = "") {
    auto i = m.find(key);
    return i == m.end() ? fallback : i->second;
}
inline Str trim(View s) {
    while (!s.empty() && std::isspace(uint8_t(s.front())))
        s.remove_prefix(1);
    while (!s.empty() && std::isspace(uint8_t(s.back())))
        s.remove_suffix(1);
    return Str(s);
}
inline std::unordered_map<Str, Str> structured(View line, View tag) {
    while (!line.empty() && std::isspace(uint8_t(line.back())))
        line.remove_suffix(1);
    Str prefix = "##" + Str(tag) + "=<";
    std::unordered_map<Str, Str> out;
    if (line.substr(0, prefix.size()) != prefix || line.back() != '>')
        return out;
    line = line.substr(prefix.size(), line.size() - prefix.size() - 1);
    size_t p = 0;
    while (p < line.size()) {
        auto e = line.find('=', p);
        if (e == View::npos)
            break;
        Str key = trim(line.substr(p, e - p)), value;
        p = e + 1;
        if (p < line.size() && line[p] == '"') {
            ++p;
            while (p < line.size()) {
                char c = line[p++];
                if (c == '\\' && p < line.size())
                    value += line[p++];
                else if (c == '"')
                    break;
                else
                    value += c;
            }
        } else {
            e = line.find(',', p);
            if (e == View::npos)
                e = line.size();
            value = trim(line.substr(p, e - p));
            p = e;
        }
        out[key] = value;
        if (p < line.size() && line[p] == ',')
            ++p;
    }
    return out;
}
} // namespace sm
