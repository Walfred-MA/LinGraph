// largeSingletonFilter: drop large singleton insertions from a merged cohort grVCF.
//
// An insertion row (INFO/SVTYPE=INS) of the main VCF (cohort.sv.vcf) is
// removed when one haplotype carries it (INFO/NSUP=1; a row without NSUP
// counts the sample columns whose GT is not 0 or .) and its size
// (INFO/MAXSIZE, else |SVLEN|) is at least --min-size. Deletions and other
// types are kept.
// Rows nested inside a removed row (CHROM = the removed row's ID, at any
// depth) are removed too, so no row is left on a sequence that is gone; with
// -c the same is done in the companion files (cohort.indel.vcf: nested small
// indels); -l is for a companion whose rows never hold other rows
// (cohort.snp.vcf: insertion SNPs), which is then read once instead of twice. Header lines are kept, except the
// ##contig lines of removed rows (an insertion's sequence is declared as a
// contig for the rows on it), plus one ##largeSingletonFilter line. Inputs
// may be plain or .gz; an output ending in .gz is written through gzip.
//
// Usage:
//   largeSingletonFilter -i cohort.sv.vcf -o filtered.sv.vcf [-s 10000]
//       [-c cohort.indel.vcf filtered.indel.vcf] [-l cohort.snp.vcf filtered.snp.vcf]
//       [-r removed.tsv] [-x removed.vcf]
//
// The cohort pipeline does not use this tool: its merge publishes these rows
// to cohort.sv.10Ksingleton.vcf directly. With -x, a merge published before
// that gets the same layout (cohort.sv.vcf + the removed rows as a VCF).

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

struct Input {
    std::string in, out;
    bool leaf = false;     // its rows never hold other rows: no first pass
    long long rows = 0, removed = 0;
};

constexpr size_t STREAM_BUFFER = 16 << 20;   // large requests to shared file systems

double seconds_since(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

void report_pass(const char *pass, const std::string &path, unsigned long long bytes, double seconds) {
    std::fprintf(stderr, "[largeSingletonFilter] %s %s: %.1f GB in %.0f s (%.0f MB/s)\n", pass, path.c_str(),
                 bytes / 1e9, seconds, seconds > 0 ? bytes / 1e6 / seconds : 0.0);
}

bool ends_with(const std::string &text, const char *suffix) {
    size_t n = std::strlen(suffix);
    return text.size() >= n && text.compare(text.size() - n, n, suffix) == 0;
}

std::string quoted(const std::string &path) {
    std::string out = "'";
    for (char c : path) {
        if (c == '\'') out += "'\\''";
        else out += c;
    }
    return out + "'";
}

// Plain files through fopen, .gz through gzip (no zlib dependency).
struct Stream {
    FILE *file = nullptr;
    bool pipe = false;
    std::string path;

    void open_read(const std::string &name) {
        path = name;
        if (ends_with(name, ".gz") || ends_with(name, ".bgz")) {
            file = popen(("gzip -cd -- " + quoted(name)).c_str(), "r");
            pipe = true;
        } else {
            file = std::fopen(name.c_str(), "r");
        }
        if (!file) fail("cannot read");
        std::setvbuf(file, nullptr, _IOFBF, STREAM_BUFFER);
    }
    void open_write(const std::string &name) {
        path = name;
        if (ends_with(name, ".gz")) {
            file = popen(("gzip -c > " + quoted(name)).c_str(), "w");
            pipe = true;
        } else {
            file = std::fopen(name.c_str(), "w");
        }
        if (!file) fail("cannot write");
        std::setvbuf(file, nullptr, _IOFBF, STREAM_BUFFER);
    }
    void close() {
        if (!file) return;
        int status = pipe ? pclose(file) : std::fclose(file);
        file = nullptr;
        if (status != 0) fail(pipe ? "gzip failed" : "close failed");
    }
    [[noreturn]] void fail(const char *what) const {
        std::fprintf(stderr, "[largeSingletonFilter] ERROR: %s: %s\n", what, path.c_str());
        std::exit(1);
    }
};

// Line reader for very long lines (one column per haplotype).
struct Lines {
    char *buffer = nullptr;
    size_t capacity = 0;
    ~Lines() { std::free(buffer); }
    bool next(FILE *file, std::string_view &line) {
        ssize_t length = getline(&buffer, &capacity, file);
        if (length < 0) return false;
        while (length > 0 && (buffer[length - 1] == '\n' || buffer[length - 1] == '\r')) --length;
        line = std::string_view(buffer, static_cast<size_t>(length));
        return true;
    }
};

// The first N tab-separated columns; returns how many were found.
size_t columns(std::string_view line, std::string_view *out, size_t n) {
    size_t found = 0, start = 0;
    while (found < n) {
        size_t tab = line.find('\t', start);
        if (tab == std::string_view::npos) {
            out[found++] = line.substr(start);
            break;
        }
        out[found++] = line.substr(start, tab - start);
        start = tab + 1;
    }
    return found;
}

// Rows nested in another row sit on that row's ID: I_/D_ (INS_/DEL_ in older
// runs) and other SV types' prefixes.
bool nested_id(std::string_view id) {
    static const char *prefixes[] = {"I_", "D_", "INS_", "DEL_", "SUB_", "DUP_", "INV_"};
    for (const char *prefix : prefixes)
        if (id.substr(0, std::strlen(prefix)) == prefix) return true;
    return false;
}

std::string_view info_value(std::string_view info, std::string_view key) {
    size_t start = 0;
    while (start <= info.size()) {
        size_t end = info.find(';', start);
        if (end == std::string_view::npos) end = info.size();
        std::string_view item = info.substr(start, end - start);
        if (item.size() > key.size() && item.substr(0, key.size()) == key && item[key.size()] == '=')
            return item.substr(key.size() + 1);
        start = end + 1;
    }
    return {};
}

long long row_size(std::string_view info) {
    std::string_view value = info_value(info, "MAXSIZE");
    if (value.empty() || value == ".") value = info_value(info, "SVLEN");
    if (!value.empty() && value[0] == '-') value.remove_prefix(1);
    long long size = 0;
    for (char c : value) {
        if (c < '0' || c > '9') return 0;
        size = size * 10 + (c - '0');
    }
    return size;
}

// Carrier columns of a row (stops counting at 2); FIRST = first carrier's column.
int carriers(std::string_view line, size_t &first) {
    size_t column = 0, start = 0;
    int count = 0;
    while (start <= line.size()) {
        size_t end = line.find('\t', start);
        if (end == std::string_view::npos) end = line.size();
        if (column >= 9) {
            std::string_view field = line.substr(start, end - start);
            std::string_view gt = field.substr(0, field.find(':'));
            if (!gt.empty() && gt != "0" && gt != ".") {
                if (count == 0) first = column;
                if (++count > 1) return count;
            }
        }
        ++column;
        start = end + 1;
    }
    return count;
}

void usage() {
    std::fprintf(stderr,
        "usage: largeSingletonFilter -i cohort.sv.vcf -o filtered.sv.vcf [-s 10000]\n"
        "           [-c IN OUT]... [-l IN OUT]... [-r removed.tsv]\n"
        "  -i, --input FILE     merged cohort SV VCF (plain or .gz)\n"
        "  -o, --output FILE    filtered VCF (.gz: compressed with gzip)\n"
        "  -s, --min-size BP    remove singleton insertions of at least BP (default 10000;\n"
        "                       deletions and other types are kept)\n"
        "  -c, --companion IN OUT  also drop the removed rows' nested rows from IN\n"
        "                       (cohort.indel.vcf); repeatable\n"
        "  -l, --leaf IN OUT    the same for a file whose rows never hold other rows\n"
        "                       (cohort.snp.vcf): read once instead of twice; repeatable\n"
        "  -r, --removed FILE   write the removed rows (ID, CHROM, POS, SVTYPE, size, carrier, reason)\n"
        "  -x, --removed-vcf FILE  write the main file's removed rows as a VCF (same header)\n");
}

}  // namespace

int main(int argc, char **argv) {
    std::vector<Input> inputs(1);
    long long min_size = 10000;
    std::string removed_path, removed_vcf;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto value = [&](int count = 1) {
            if (i + count >= argc) { usage(); std::exit(2); }
        };
        if (arg == "-i" || arg == "--input") { value(); inputs[0].in = argv[++i]; }
        else if (arg == "-o" || arg == "--output") { value(); inputs[0].out = argv[++i]; }
        else if (arg == "-s" || arg == "--min-size") {
            value();
            char *end = nullptr;
            min_size = std::strtoll(argv[++i], &end, 10);
            if (!end || *end || min_size < 1) {
                std::fprintf(stderr, "[largeSingletonFilter] ERROR: --min-size must be a positive integer\n");
                return 2;
            }
        }
        else if (arg == "-c" || arg == "--companion" || arg == "-l" || arg == "--leaf") {
            value(2);
            Input companion;
            companion.leaf = arg == "-l" || arg == "--leaf";
            companion.in = argv[++i];
            companion.out = argv[++i];
            inputs.push_back(companion);
        }
        else if (arg == "-r" || arg == "--removed") { value(); removed_path = argv[++i]; }
        else if (arg == "-x" || arg == "--removed-vcf") { value(); removed_vcf = argv[++i]; }
        else if (arg == "-h" || arg == "--help") { usage(); return 0; }
        else { usage(); return 2; }
    }
    if (inputs[0].in.empty() || inputs[0].out.empty()) { usage(); return 2; }
    for (const Input &input : inputs)
        if (input.in == input.out) {
            std::fprintf(stderr, "[largeSingletonFilter] ERROR: output equals input: %s\n", input.in.c_str());
            return 2;
        }

    // Pass 1: large singletons of the main file; nested rows of every file
    // under their parent row's ID.
    std::unordered_set<std::string> removed;
    std::unordered_map<std::string, std::vector<std::string>> children;
    std::unordered_map<std::string, std::string> reason;
    std::vector<std::string> samples;
    Stream removed_out;
    if (!removed_path.empty()) {
        removed_out.open_write(removed_path);
        std::fprintf(removed_out.file, "ID\tCHROM\tPOS\tSVTYPE\tSIZE\tCARRIER\tREASON\n");
    }
    long long seeds = 0;
    for (size_t index = 0; index < inputs.size(); ++index) {
        if (inputs[index].leaf) continue;
        Stream in;
        in.open_read(inputs[index].in);
        Lines lines;
        std::string_view line;
        auto started = std::chrono::steady_clock::now();
        unsigned long long bytes = 0;
        while (lines.next(in.file, line)) {
            bytes += line.size() + 1;
            if (line.empty()) continue;
            if (line[0] == '#') {
                if (index == 0 && line.substr(0, 6) == "#CHROM") {
                    std::string_view rest = line;
                    size_t column = 0, start = 0;
                    while (start <= rest.size()) {
                        size_t end = rest.find('\t', start);
                        if (end == std::string_view::npos) end = rest.size();
                        if (column >= 9) samples.emplace_back(rest.substr(start, end - start));
                        ++column;
                        start = end + 1;
                    }
                }
                continue;
            }
            std::string_view field[8];
            if (columns(line, field, 8) < 8) continue;
            if (nested_id(field[0]) && nested_id(field[2]))
                children[std::string(field[0])].emplace_back(field[2]);
            if (index != 0 || info_value(field[7], "SVTYPE") != "INS") continue;
            long long size = row_size(field[7]);
            if (size < min_size) continue;
            // NSUP: merged rows count their supporting haplotypes; the
            // sample columns are read only to name the carrier.
            std::string_view nsup = info_value(field[7], "NSUP");
            size_t first = 0;
            if (!nsup.empty() && nsup != ".") {
                if (nsup != "1") continue;
                if (removed_out.file) carriers(line, first);
            } else if (carriers(line, first) != 1) {
                continue;
            }
            std::string id(field[2]);
            if (!removed.insert(id).second) continue;
            ++seeds;
            if (removed_out.file) {
                std::string svtype(info_value(field[7], "SVTYPE"));
                std::string carrier = first >= 9 && first - 9 < samples.size() ? samples[first - 9] : "";
                std::fprintf(removed_out.file, "%s\t%.*s\t%.*s\t%s\t%lld\t%s\tsingleton\n", id.c_str(),
                             (int)field[0].size(), field[0].data(), (int)field[1].size(), field[1].data(),
                             svtype.c_str(), size, carrier.c_str());
            }
        }
        in.close();
        report_pass("pass 1 read", inputs[index].in, bytes, seconds_since(started));
    }

    // Everything nested under a removed row, at any depth.
    std::vector<std::string> stack(removed.begin(), removed.end());
    long long nested = 0;
    while (!stack.empty()) {
        std::string parent = std::move(stack.back());
        stack.pop_back();
        auto found = children.find(parent);
        if (found == children.end()) continue;
        for (const std::string &child : found->second)
            if (removed.insert(child).second) {
                ++nested;
                reason[child] = parent;
                stack.push_back(child);
            }
    }
    if (removed_out.file) {
        for (const auto &item : reason)
            std::fprintf(removed_out.file, "%s\t\t\t\t\t\tnested_in:%s\n", item.first.c_str(), item.second.c_str());
        removed_out.close();
    }

    // Pass 2: write every file without the removed rows and the rows on them.
    std::string header = "##largeSingletonFilter=<MinSize=" + std::to_string(min_size)
        + ",Singletons=" + std::to_string(seeds) + ",NestedRemoved=" + std::to_string(nested)
        + ",Description=\"haplotype-singleton insertions of at least MinSize bp removed, with the rows nested in them\">\n";
    for (size_t index = 0; index < inputs.size(); ++index) {
        Input &input = inputs[index];
        Stream in, out, extra;
        in.open_read(input.in);
        out.open_write(input.out);
        // The main file's removed rows, as a VCF with the full header.
        if (index == 0 && !removed_vcf.empty()) extra.open_write(removed_vcf);
        Lines lines;
        std::string_view line;
        auto started = std::chrono::steady_clock::now();
        unsigned long long bytes = 0;
        while (lines.next(in.file, line)) {
            bytes += line.size() + 1;
            if (!line.empty() && line[0] == '#') {
                if (extra.file) {
                    if (line.substr(0, 6) == "#CHROM") std::fputs(header.c_str(), extra.file);
                    std::fwrite(line.data(), 1, line.size(), extra.file);
                    std::fputc('\n', extra.file);
                }
                if (line.substr(0, 6) == "#CHROM") std::fputs(header.c_str(), out.file);
                if (line.substr(0, 13) == "##contig=<ID=") {
                    std::string_view id = line.substr(13);
                    id = id.substr(0, id.find_first_of(",>"));
                    if (removed.count(std::string(id))) continue;
                }
            } else if (!line.empty()) {
                std::string_view field[3];
                if (columns(line, field, 3) >= 3) {
                    ++input.rows;
                    if (!removed.empty() && (removed.count(std::string(field[2]))
                                             || removed.count(std::string(field[0])))) {
                        ++input.removed;
                        if (extra.file) {
                            std::fwrite(line.data(), 1, line.size(), extra.file);
                            std::fputc('\n', extra.file);
                        }
                        continue;
                    }
                }
            }
            std::fwrite(line.data(), 1, line.size(), out.file);
            std::fputc('\n', out.file);
        }
        in.close();
        out.close();
        extra.close();
        report_pass("pass 2 read+write", input.in, bytes, seconds_since(started));
        std::fprintf(stderr, "[largeSingletonFilter] %s: %lld row(s), %lld removed -> %s\n",
                     input.in.c_str(), input.rows, input.removed, input.out.c_str());
    }
    std::fprintf(stderr, "[largeSingletonFilter] %lld singleton insertion(s) >= %lld bp; %lld nested row(s) removed with them\n",
                 seeds, min_size, nested);
    return 0;
}
