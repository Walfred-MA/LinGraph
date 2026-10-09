// GrvcfGraph: exact merged grVCF -> rGFA (segments, links, paths). Assembly-free:
// it reads the merged grVCFs and the backbone/alternative/template/catalog
// FASTAs, never an assembly.
//
//   GrvcfGraph -v VCF [VCF ...] -r BACKBONE.fa [-a ALT.fa ...]
//              [--local-reference-templates TEMPLATES.fa ...] [--local-path-fasta CATALOG.fa ...]
//              [--reference-haplotype CHM13_h1] -o OUT.gfa [-t THREADS]
//              [--max-node-length 1024] [--placements FILE] [--roots-table FILE]
//              [--gaf FOLDER [-s SAMPLE_VCF_OR_HEADERS ...] [--query-lengths TSV]]
//
// 1. Scan: VCF columns 1-8 only (sample columns are skipped); every row is a
//    fixed 64-byte record, allele bases come from ALT / INFO/SEQ.
// 2. Paths: the names kept rows are placed on that are not rows are resolved
//    to FASTA records exactly as the Python converter does (_index_roots,
//    resolve_local_sources, template_lift, add_template_roots): -r/-a records,
//    lifted templates (own path, or an alias of a sequence-identical included
//    interval), lookup-only catalog names (aliases) and --exact DUP_ templates.
// 3. Segments: every kept row (PASS rows and what they are placed on) is
//    placed on the sequence carrying its parent's bases; each sequence is cut
//    at those breakpoints, at emitted alias and lift points, and every
//    --max-node-length bases, one S line per piece. Sequences are numbered
//    paths first, then row alleles by nesting level and VCF order.
// 4. Links: consecutive segments of a sequence; each row's allele entered
//    from the base before it on its parent and left to the base after it
//    (through the parent's own flanks at its ends; deletions link across);
//    lifted templates attached at their lift points. These are the Python
//    converter's link walks with a 1-base flank; carrier junctions come later
//    from the GAF walks.
// 5. Paths: one P line per emitted path (backbone, alternatives, novels,
//    templates, DUP_ templates); an alias walks the segments it aliases.
//    Variant alleles get no P line (SN/SO tags name them).
// 6. GAF (--gaf): one SAMPLE.gaf per sample with ##pseudoLinearMapping lines
//    (from -s per-sample VCFs or cohort.samples.headers.gz, else the merged
//    headers), walked as gfa_sample_gaf.py does from the merged sample columns.
//    Links the walks need are added to the GFA (L lines without SR). Column 2
//    is --query-lengths, else the largest mapped query coordinate. Also
//    SAMPLE.gaf.breaks/unplaced/lengthdebug.tsv and gaf_stats.tsv.
//
// Scale: the VCFs are read once, uncompressed inputs in parallel byte ranges.
// The scan writes the samples' alleles (20 bytes each) to per-worker spill
// files in --spill-dir (default: the GAF folder), grouped by sample, so memory
// does not grow with the cohort; each sample is read back by its walker.
// --gaf-threads samples are walked at once (default --threads), each holding
// about 400 bytes per carried allele.
//
// Coordinates follow gfa_stable_coords.StableResolver._breakpoints: a SNP
// replaces [POS-1, POS); other rows start at POS on a backbone/alternative
// path and at POS-1 on a row or DUP_ template (nested POS is 1-based), end at
// END (start for a pure insertion), a DEL with SVLEN at start+|SVLEN|. Rows
// nested in a DEL replace that deletion's interval. A --exact full-locus
// duplication site (PACLASS=fulllocusdup, EXTENDGRAPHCIGAR >DUP_..:N=) walks
// its shared DUP_ template; rows nested in the site are placed on it.

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cerrno>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sys/mman.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <thread>
#include <tuple>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::string;
using std::string_view;
using std::vector;

namespace {

enum Kind : uint8_t { INS = 0, DEL = 1, SNP = 2, SUB = 3 };
enum Flag : uint8_t { PASS = 1, SITE = 2 };
const char *KIND_NAMES[] = {"insertion", "deletion", "snp", "substitution"};

// One VCF data row. Offsets point into the names / seq arenas.
struct Row {
    uint64_t name;       // row ID
    uint64_t chrom;      // the name it is placed on (a row ID or a FASTA path)
    uint64_t seq;        // allele bases in the seq arena; SITE: DUP_ template name in names
    uint64_t line;       // data line index, numbered from 0 across the inputs
    uint32_t name_len, chrom_len, seq_len;
    uint32_t pos, end;   // VCF POS, INFO/END (POS when absent)
    uint32_t svlen;      // |SVLEN| (0 when absent)
    uint32_t size;       // allele length: INS SVLEN, SUB/SNP bases, DEL 0
    uint8_t kind, flags;
    uint8_t pad[2];
};
static_assert(sizeof(Row) == 64, "Row layout");

[[noreturn]] void fail(const string &message) { throw std::runtime_error(message); }

double peak_gib() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
#ifdef __APPLE__
    return usage.ru_maxrss / double(1ull << 30);
#else
    return usage.ru_maxrss / double(1ull << 20);
#endif
}

auto started = std::chrono::steady_clock::now();
void log(const string &message) {
    double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    fprintf(stderr, "[GrvcfGraph %.1fs] %s; peak RSS=%.2f GiB\n", seconds, message.c_str(), peak_gib());
    fflush(stderr);
}

uint64_t hash_name(string_view text) {
    uint64_t value = 1469598103934665603ull;
    for (unsigned char c : text) value = (value ^ c) * 1099511628211ull;
    value ^= value >> 33;
    value *= 0xff51afd7ed558ccdull;
    value ^= value >> 33;
    return value;
}

bool parse_u64(string_view text, uint64_t &value) {
    if (text.empty()) return false;
    auto result = std::from_chars(text.data(), text.data() + text.size(), value);
    return result.ec == std::errc() && result.ptr == text.data() + text.size();
}

bool parse_i64(string_view text, int64_t &value) {
    if (text.empty()) return false;
    auto result = std::from_chars(text.data(), text.data() + text.size(), value);
    return result.ec == std::errc() && result.ptr == text.data() + text.size();
}

vector<string> split(string_view text, char separator);
bool starts_with(string_view text, string_view prefix);
bool all_digits(string_view text);

uint32_t to_u32(int64_t value, const char *what) {
    if (value < 0 || value > int64_t(UINT32_MAX)) fail(string(what) + " is outside 0.." + std::to_string(UINT32_MAX));
    return uint32_t(value);
}

bool is_alpha(string_view text) {
    if (text.empty()) return false;
    for (unsigned char c : text)
        if (!((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'))) return false;
    return true;
}

// gfa_query_anchors._unescape
string unescape(string_view text) {
    if (text.find('@') == string_view::npos) return string(text);
    string out;
    out.reserve(text.size());
    for (size_t i = 0; i < text.size(); ++i) {
        if (text[i] == '@' && i + 2 < text.size()) {
            string_view code = text.substr(i + 1, 2);
            char c = 0;
            if (code == "0A") c = '\n';
            else if (code == "09") c = '\t';
            else if (code == "2C") c = ',';
            else if (code == "3B") c = ';';
            else if (code == "3A") c = ':';
            else if (code == "40") c = '@';
            if (c) {
                out.push_back(c);
                i += 2;
                continue;
            }
        }
        out.push_back(text[i]);
    }
    return out;
}

// INFO/SEQ as plain bases. A graph-encoded SEQ is accepted only when all its
// query bases are literal payloads (I/X/S); copied (=/M) bases are an error.
string literal_sequence(string_view raw, uint32_t size) {
    if (raw.empty() || raw == ".") fail("INS/SUB requires INFO/SEQ");
    string text = unescape(raw);
    if (text[0] != '<' && text[0] != '>') {
        if (!is_alpha(text)) fail("INFO/SEQ is not plain bases");
        if (text.size() != size) fail("INFO/SEQ length " + std::to_string(text.size()) + " differs from SVLEN " + std::to_string(size));
        return text;
    }
    string out;
    size_t position = 0;
    while (position < text.size()) {
        if (text[position] != '<' && text[position] != '>') fail("INFO/SEQ: text outside an encoded graph chunk");
        size_t end = text.find_first_of("<>", position + 1);
        if (end == string::npos) end = text.size();
        size_t body = position + 1;
        size_t colon = text.find(':', body);
        bool target = colon != string::npos && colon < end;
        if (target) body = colon + 1;
        size_t cursor = body;
        while (cursor < end) {
            size_t digits = cursor;
            while (digits < end && text[digits] >= '0' && text[digits] <= '9') ++digits;
            if (digits == cursor || digits >= end) fail("INFO/SEQ: malformed graph mapping operation");
            uint64_t length = 0;
            if (!parse_u64(string_view(text).substr(cursor, digits - cursor), length))
                fail("INFO/SEQ: graph mapping operation length is out of range");
            char op = text[digits];
            if (!strchr("=MXIDHS", op)) fail("INFO/SEQ: malformed graph mapping operation");
            cursor = digits + 1;
            size_t payload = cursor;
            while (payload < end && isalpha((unsigned char)text[payload])) ++payload;
            if (op == 'X' && length && !target) fail("INFO/SEQ: unnamed X operation has no graph target");
            if ((op == '=' || op == 'M') && length)
                fail("INFO/SEQ copies bases from another path (" + std::to_string(length) + string(1, op) +
                     "); the exact grVCF must carry literal insertion bases");
            if ((op == 'I' || op == 'X' || op == 'S') && length) {
                if (payload - cursor != length) fail("INFO/SEQ: missing or short literal payload");
                out.append(text, cursor, length);
            }
            cursor = payload;
        }
        position = end;
    }
    if (out.size() != size) fail("decoded INFO/SEQ length " + std::to_string(out.size()) + " differs from SVLEN " + std::to_string(size));
    return out;
}

// DUP_<locus>_<start>_<end> (gfa_interval_metadata.TEMPLATE_PATH, greedy locus).
bool template_name(const string &name, string *locus, uint64_t *start, uint64_t *end) {
    if (name.compare(0, 4, "DUP_") != 0) return false;
    size_t last = name.rfind('_');
    if (last == string::npos || last < 6) return false;
    size_t middle = name.rfind('_', last - 1);
    if (middle == string::npos || middle < 5) return false;
    string_view first = string_view(name).substr(middle + 1, last - middle - 1);
    string_view second = string_view(name).substr(last + 1);
    uint64_t a = 0, b = 0;
    for (string_view digits : {first, second})
        if (digits.empty() || digits.find_first_not_of("0123456789") != string_view::npos) return false;
    if (!parse_u64(first, a) || !parse_u64(second, b)) return false;
    if (locus) *locus = name.substr(4, middle - 4);
    if (start) *start = a;
    if (end) *end = b;
    return true;
}

// ------------------------------------------------------------------ scan

// Lines of a file, or of the lines starting in byte range [begin, end) of an
// uncompressed file (begin is a line start).
class LineReader {
  public:
    explicit LineReader(const string &path, uint64_t begin = 0, uint64_t end = UINT64_MAX)
        : path_(path), buffer_(16u << 20), position_(begin), limit_(end) {
        if (path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0) {
            if (begin) fail("internal: a .gz input cannot be read from an offset");
            string quoted;
            for (char c : path) quoted += c == '\'' ? string("'\\''") : string(1, c);
            string command = "gzip -dc '" + quoted + "'";
            file_ = popen(command.c_str(), "r");
            piped_ = true;
        } else {
            file_ = fopen(path.c_str(), "rb");
            if (file_ && begin && fseeko(file_, off_t(begin), SEEK_SET)) fail("cannot seek in " + path);
        }
        if (!file_) fail("cannot open " + path + ": " + strerror(errno));
    }
    // File offset of the line next() returned last.
    uint64_t offset() const { return line_start_; }
    ~LineReader() {
        if (file_) piped_ ? pclose(file_) : fclose(file_);
    }
    bool next(string_view &line) {
        if (position_ >= limit_) return false;  // the next line belongs to the next range
        for (;;) {
            if (begin_ < end_) {
                const char *start = buffer_.data() + begin_;
                const void *newline = memchr(start, '\n', end_ - begin_);
                if (newline) {
                    size_t length = (const char *)newline - start;
                    begin_ += length + 1;
                    line_start_ = position_;
                    position_ += length + 1;
                    if (length && start[length - 1] == '\r') --length;
                    line = string_view(start, length);
                    return true;
                }
            }
            if (eof_) {
                if (begin_ < end_) {
                    size_t length = end_ - begin_;
                    line_start_ = position_;
                    position_ += length;
                    if (buffer_[begin_ + length - 1] == '\r') --length;
                    line = string_view(buffer_.data() + begin_, length);
                    begin_ = end_;
                    return true;
                }
                return false;
            }
            if (begin_) {
                memmove(buffer_.data(), buffer_.data() + begin_, end_ - begin_);
                end_ -= begin_;
                begin_ = 0;
            }
            if (end_ == buffer_.size()) buffer_.resize(buffer_.size() * 2);
            size_t got = fread(buffer_.data() + end_, 1, buffer_.size() - end_, file_);
            if (got == 0) {
                if (ferror(file_)) fail("read error in " + path_);
                eof_ = true;
                if (piped_) {
                    // A truncated or corrupt .gz must not pass as a shorter file.
                    int status = pclose(file_);
                    file_ = nullptr;
                    if (status != 0) fail("gzip -dc failed on " + path_ + " (truncated or corrupt file?)");
                }
            }
            end_ += got;
        }
    }

  private:
    string path_;
    FILE *file_ = nullptr;
    bool piped_ = false, eof_ = false;
    vector<char> buffer_;
    size_t begin_ = 0, end_ = 0;
    uint64_t position_ = 0, limit_ = UINT64_MAX, line_start_ = 0;
};

struct Scanned {
    vector<Row> rows;
    string names, seq;
    uint64_t lines = 0;           // data lines, including skipped rows
    vector<uint64_t> file_starts;  // first data line of each input
    string error;
    // GAF: query contigs of the collected alleles, numbered within one range
    vector<string> contigs;
    std::unordered_map<string, uint32_t> contig_ids;
    // joined scan: each range's contigs and the index of its first row
    vector<vector<string>> chunk_contigs;
    vector<uint64_t> chunk_row_base;
};

// ---- GAF alleles, collected by the scan (gfa_sample_gaf.build_shards)
// Empty or whitespace-only line (skipped wherever VCF lines are read).
bool blank_line(string_view line) {
    for (char c : line)
        if (!isspace((unsigned char)c)) return false;
    return true;
}

// A sample's carried allele as the scan saw it: contig and row numbered within
// the scanned range `chunk`; mapped to global numbers when its sample is walked.
struct DiskObs {
    uint32_t contig, row, qs, qe;
    uint16_t chunk;
    uint8_t strand, flags;
};
static_assert(sizeof(DiskObs) == 20, "DiskObs layout");
enum : uint8_t {
    OBS_NOT_SNP_TYPE = 1,    // TYPE names this observation something else than SNP
    OBS_COORD_OVERFLOW = 2,  // QUERYCOORD beyond 2^32: an error only if the row is walked
};

// Records of one sample in one spill file: [offset, offset + count) records.
struct ObsSegment {
    uint32_t sample, count;
    uint64_t offset;
};

// One scan worker's spill file: blocks of records grouped by sample.
struct ObsSpill {
    string path;
    FILE *file = nullptr;
    int fd = -1;
    uint64_t written = 0;  // records
    vector<DiskObs> buffer;
    vector<uint32_t> owner;  // sample of each buffered record
    vector<ObsSegment> segments;
};

// Every sample's alleles spilled to disk during the scan, so memory does not
// grow with the cohort: each worker appends blocks grouped by sample to its own
// file, and a sample's records are read back by its walker.
class ObsCollector {
  public:
    // `prefix` + 6 random characters: unique also across hosts sharing the folder.
    ObsCollector(const string &prefix, size_t samples) : folder_(prefix + "XXXXXX"), samples_(samples) {
        if (!mkdtemp(&folder_[0])) fail("cannot create " + folder_ + ": " + strerror(errno));
        if (const char *value = getenv("GRVCFGRAPH_SPILL_RECORDS"))  // tests: force small blocks
            if (!parse_u64(value, block_) || !block_) fail("GRVCFGRAPH_SPILL_RECORDS must be a positive integer");
    }
    ~ObsCollector() {
        for (auto &spill : spills)
            if (spill.file) fclose(spill.file);
        rmdir(folder_.c_str());
    }
    std::unordered_map<string, size_t> slot;  // sample name -> index
    vector<vector<int64_t>> columns;          // per input: #CHROM column -> sample index (-1)
    vector<ObsSpill> spills;                  // per scan worker

    void start(size_t workers) {
        spills.resize(workers);
        for (size_t w = 0; w < workers; ++w) {
            spills[w].path = folder_ + "/worker_" + std::to_string(w) + ".obs";
            spills[w].file = fopen(spills[w].path.c_str(), "w+b");
            if (!spills[w].file) fail("cannot write " + spills[w].path + ": " + strerror(errno));
            unlink(spills[w].path.c_str());  // the open file stays readable; the kernel frees it at exit
        }
        rmdir(folder_.c_str());
    }
    void add(ObsSpill &spill, size_t sample, const DiskObs &record) {
        spill.buffer.push_back(record);
        spill.owner.push_back(uint32_t(sample));
    }
    bool full(const ObsSpill &spill) const { return spill.buffer.size() >= block_; }
    // Append the buffer as one block grouped by sample (append order kept within a sample).
    void flush(ObsSpill &spill) {
        size_t count = spill.buffer.size();
        if (!count) return;
        vector<uint64_t> start(samples_ + 1, 0);
        for (uint32_t sample : spill.owner) ++start[sample + 1];
        for (size_t k = 0; k < samples_; ++k) start[k + 1] += start[k];
        vector<uint64_t> next(start.begin(), start.end() - 1);
        vector<DiskObs> grouped(count);
        for (size_t i = 0; i < count; ++i) grouped[next[spill.owner[i]]++] = spill.buffer[i];
        if (fwrite(grouped.data(), sizeof(DiskObs), count, spill.file) != count) fail("write failed: " + spill.path);
        for (size_t k = 0; k < samples_; ++k)
            if (start[k + 1] > start[k])
                spill.segments.push_back({uint32_t(k), uint32_t(start[k + 1] - start[k]), spill.written + start[k]});
        spill.written += count;
        spill.buffer.clear();
        spill.owner.clear();
    }
    // After the scan: flush, reopen for reading, index every sample's segments.
    void finish() {
        by_sample_.assign(samples_, {});
        for (size_t w = 0; w < spills.size(); ++w) {
            ObsSpill &spill = spills[w];
            flush(spill);
            vector<DiskObs>().swap(spill.buffer);
            vector<uint32_t>().swap(spill.owner);
            if (fflush(spill.file) || ferror(spill.file)) fail("write failed: " + spill.path);
            spill.fd = fileno(spill.file);  // read back with pread; closed with spill.file
            for (auto &segment : spill.segments) by_sample_[segment.sample].push_back({uint32_t(w), segment});
            total_ += spill.written;
        }
    }
    uint64_t total() const { return total_; }
    // One sample's records, in spill order (thread-safe: pread).
    vector<DiskObs> load(size_t sample) const {
        uint64_t count = 0;
        for (auto &entry : by_sample_[sample]) count += entry.second.count;
        vector<DiskObs> out(count);
        char *at = reinterpret_cast<char *>(out.data());
        for (auto &[worker, segment] : by_sample_[sample]) {
            size_t size = size_t(segment.count) * sizeof(DiskObs), done = 0;
            off_t offset = off_t(segment.offset * sizeof(DiskObs));
            while (done < size) {
                ssize_t got = pread(spills[worker].fd, at + done, size - done, offset + off_t(done));
                if (got <= 0) fail("read failed: " + spills[worker].path);
                done += size_t(got);
            }
            at += size;
        }
        return out;
    }

  private:
    string folder_;
    size_t samples_;
    uint64_t total_ = 0, block_ = 4u << 20;  // records buffered per worker before a block is written
    vector<vector<std::pair<uint32_t, ObsSegment>>> by_sample_;
};

// FORMAT key positions, kept while consecutive rows share their FORMAT.
struct FormatKeys {
    string text;
    int64_t gt = -1, contig = -1, coord = -1, type = -1;
    size_t keys = 0;
};

// The sample-column alleles of one scanned row for the GAF: every carrier
// column's ASSEMBLYCONTIG/QUERYCOORD pairs, for samples that get a GAF.
void collect_alleles(string_view line, uint32_t row, bool snp, uint16_t chunk, const vector<int64_t> &columns,
                     Scanned &part, ObsCollector &collector, ObsSpill &spill, FormatKeys &format,
                     vector<string_view> &values, vector<string_view> &names, vector<string_view> &coords,
                     vector<string_view> &kinds) {
    auto split_view = [](string_view text, char separator, vector<string_view> &out) {
        out.clear();
        size_t start = 0;
        for (;;) {
            size_t stop = text.find(separator, start);
            out.push_back(text.substr(start, stop == string_view::npos ? string_view::npos : stop - start));
            if (stop == string_view::npos) return;
            start = stop + 1;
        }
    };
    // FORMAT is column 9 (index 8); sample columns follow.
    size_t position = 0;
    for (int i = 0; i < 8; ++i) {
        position = line.find('\t', position);
        if (position == string_view::npos) return;
        ++position;
    }
    size_t stop = line.find('\t', position);
    if (stop == string_view::npos) return;  // fewer than 10 columns
    string_view format_text = line.substr(position, stop - position);
    if (format_text != format.text) {
        format = FormatKeys();
        format.text = string(format_text);
        split_view(format_text, ':', values);
        format.keys = values.size();
        for (size_t k = 0; k < values.size(); ++k) {
            if (values[k] == "GT" && format.gt < 0) format.gt = int64_t(k);
            else if (values[k] == "ASSEMBLYCONTIG" && format.contig < 0) format.contig = int64_t(k);
            else if (values[k] == "QUERYCOORD" && format.coord < 0) format.coord = int64_t(k);
            else if (values[k] == "TYPE" && format.type < 0) format.type = int64_t(k);
        }
    }
    if (format.gt < 0 || format.contig < 0 || format.coord < 0) return;
    position = stop + 1;
    for (size_t c = 9; position <= line.size(); ++c) {
        size_t end = line.find('\t', position);
        if (end == string_view::npos) end = line.size();
        string_view text = line.substr(position, end - position);
        position = end + 1;
        if (c >= columns.size() || columns[c] < 0 || text.empty() || text[0] != '1') continue;
        split_view(text, ':', values);
        if (values.size() != format.keys || values[format.gt] != "1") continue;
        split_view(values[format.contig], ',', names);
        split_view(values[format.coord], ',', coords);
        if (format.type >= 0) split_view(values[format.type], ',', kinds);
        else kinds.clear();
        for (size_t k = 0; k < coords.size(); ++k) {
            // START[-END][+|-]
            string_view coord = coords[k];
            char strand = 0;
            if (!coord.empty() && (coord.back() == '+' || coord.back() == '-')) {
                strand = coord.back();
                coord.remove_suffix(1);
            }
            size_t dash = coord.find('-');
            string_view first = coord.substr(0, dash), second = dash == string_view::npos ? string_view() : coord.substr(dash + 1);
            if (!all_digits(first) || (dash != string_view::npos && !all_digits(second))) continue;
            size_t which = names.size() > 1 ? k : 0;
            if (which >= names.size()) continue;
            uint64_t a = 0, b = 0;
            bool overflow = !parse_u64(first, a) || a > UINT32_MAX ||
                            (dash != string_view::npos && (!parse_u64(second, b) || b > UINT32_MAX));
            int64_t low = 0, high = 0;
            if (!overflow) {
                low = int64_t(a);
                if (snp && dash == string_view::npos && strand == '-') low -= 1;  // '-' points are traversal boundaries
                high = dash != string_view::npos ? int64_t(b) : low + (snp ? 1 : 0);
                if (low > high) std::swap(low, high);
                if (low < 0) continue;
                if (high > int64_t(UINT32_MAX)) overflow = true;
            }
            if (overflow) low = high = 0;  // rows that are not walked may carry any coordinate
            auto found = part.contig_ids.find(string(names[which]));
            uint32_t contig;
            if (found == part.contig_ids.end()) {
                contig = uint32_t(part.contigs.size());
                part.contig_ids.emplace(string(names[which]), contig);
                part.contigs.emplace_back(names[which]);
            } else {
                contig = found->second;
            }
            uint8_t flags = (format.type >= 0 && kinds.size() == coords.size() && kinds[k] != "SNP") ? OBS_NOT_SNP_TYPE : 0;
            if (overflow) flags |= OBS_COORD_OVERFLOW;
            collector.add(spill, size_t(columns[c]),
                          {contig, row, uint32_t(low), uint32_t(high), chunk, uint8_t(strand == '-'), flags});
        }
    }
}

void add_name(string &arena, string_view text, uint64_t &offset, uint32_t &length) {
    offset = arena.size();
    length = uint32_t(text.size());
    arena.append(text);
}

// One data line -> a Row (false: skipped, as gfa_interval_metadata skips SVLEN=0 INS).
bool parse_row(string_view line, Scanned &out, Row &row) {
    string_view field[8];
    size_t position = 0;
    for (int i = 0; i < 8; ++i) {
        size_t tab = line.find('\t', position);
        if (tab == string_view::npos) {
            if (i < 7) fail("malformed VCF record (fewer than 8 columns)");
            tab = line.size();
        }
        field[i] = line.substr(position, tab - position);
        position = tab + 1;
    }
    string_view svtype, svlen, end, seq, paclass, extend;
    bool has_end = false;
    string_view info = field[7];
    size_t cursor = 0;
    while (cursor <= info.size()) {
        size_t semicolon = info.find(';', cursor);
        if (semicolon == string_view::npos) semicolon = info.size();
        string_view item = info.substr(cursor, semicolon - cursor);
        size_t equals = item.find('=');
        if (equals != string_view::npos) {
            string_view key = item.substr(0, equals), value = item.substr(equals + 1);
            if (key == "SVTYPE") svtype = value;
            else if (key == "SVLEN") svlen = value;
            else if (key == "END") end = value, has_end = true;
            else if (key == "SEQ") seq = value;
            else if (key == "PACLASS") paclass = value;
            else if (key == "EXTENDGRAPHCIGAR") extend = value;
        }
        cursor = semicolon + 1;
    }
    memset(&row, 0, sizeof(row));
    int64_t pos = 0;
    if (!parse_i64(field[1], pos)) fail("invalid POS");
    row.pos = to_u32(pos, "POS");
    row.end = row.pos;
    if (field[6] == "PASS" || field[6] == ".") row.flags |= PASS;
    // END and SVLEN are read only for the types that use them (SNPs ignore both).
    auto read_end = [&] {
        if (!has_end) return;
        int64_t value = 0;
        if (!parse_i64(end, value)) fail("invalid INFO/END");
        row.end = to_u32(value, "END");
    };
    auto read_svlen = [&] {
        if (svlen.empty() || svlen == ".") return;
        int64_t value = 0;
        if (!parse_i64(svlen, value)) fail("invalid SVLEN");
        if (value == INT64_MIN) fail("SVLEN is outside the supported range");
        row.svlen = to_u32(value < 0 ? -value : value, "SVLEN");
    };
    string allele;
    if (svtype == "INS") {
        row.kind = INS;
        if (svlen.empty() || svlen == ".") fail("INS requires a numeric SVLEN");
        read_svlen();
        if (row.svlen == 0) return false;
        read_end();
        row.size = row.svlen;
        if (seq.empty() || seq == ".") fail("INS requires INFO/SEQ");
        string site;
        if (paclass == "fulllocusdup" && !extend.empty()) {
            // >DUP_<locus>_<start>_<end>:<size>=
            string text = unescape(extend);
            size_t colon = text.rfind(':');
            uint64_t length = 0;
            if (text.size() > 6 && text.compare(0, 5, ">DUP_") == 0 && text.back() == '=' && colon != string::npos &&
                parse_u64(string_view(text).substr(colon + 1, text.size() - colon - 2), length) && length == row.size &&
                template_name(text.substr(1, colon - 1), nullptr, nullptr, nullptr))
                site = text.substr(1, colon - 1);
        }
        if (!site.empty()) {
            row.flags |= SITE;
            uint64_t offset;
            add_name(out.names, site, offset, row.seq_len);
            row.seq = offset;
        } else {
            allele = literal_sequence(seq, row.size);
        }
    } else {
        string_view type = svtype;
        if (type.empty() && field[3].size() == 1 && field[4].size() == 1 && is_alpha(field[4])) type = "SNP";
        if (type == "SNP") {
            row.kind = SNP;
            if (!is_alpha(field[4])) fail("SNP requires a literal ALT allele");
            allele = string(field[4]);
            row.end = row.pos;
        } else if (type == "DEL") {
            row.kind = DEL;
            if (!has_end) fail("DEL requires INFO/END");
            read_svlen();
            read_end();
            if (!svlen.empty() && svlen != "." && row.svlen == 0) fail("deletion SVLEN must be nonzero");
            if (row.svlen == 0 && row.end <= row.pos) fail("deletion END must be greater than POS");
        } else if (type == "SUB") {
            row.kind = SUB;
            if (!has_end) fail("SUB requires INFO/END");
            read_end();
            allele = unescape(seq);
            if (!is_alpha(allele)) fail("SUB requires a literal INFO/SEQ allele");
        } else {
            fail("unsupported variant type " + string(type.empty() ? field[4] : type));
        }
    }
    if (!allele.empty()) {
        row.seq = out.seq.size();
        row.seq_len = uint32_t(allele.size());
        row.size = row.seq_len;
        out.seq.append(allele);
    }
    add_name(out.names, field[2], row.name, row.name_len);
    // Consecutive rows usually share their parent: store it once.
    if (!out.rows.empty()) {
        const Row &previous = out.rows.back();
        if (previous.chrom_len == field[0].size() &&
            out.names.compare(previous.chrom, previous.chrom_len, field[0].data(), field[0].size()) == 0) {
            row.chrom = previous.chrom;
            row.chrom_len = previous.chrom_len;
            return true;
        }
    }
    add_name(out.names, field[0], row.chrom, row.chrom_len);
    return true;
}

// A byte range of one input: whole lines, scanned by one thread.
struct Chunk {
    size_t file;
    uint64_t begin, end;  // end UINT64_MAX: to the end of the file
};

void scan_chunk(const string &path, const Chunk &chunk, uint16_t chunk_index, Scanned &out,
                ObsCollector *collector, ObsSpill *spill) {
    try {
        LineReader reader(path, chunk.begin, chunk.end);
        string_view line;
        uint64_t number = 0;
        FormatKeys format;
        vector<string_view> values, names, coords, kinds;
        while (reader.next(line)) {
            ++number;
            if (line.empty() || line[0] == '#' || blank_line(line)) continue;
            Row row;
            try {
                if (parse_row(line, out, row)) {
                    row.line = out.lines;
                    if (collector) {
                        collect_alleles(line, uint32_t(out.rows.size()), row.kind == SNP, chunk_index,
                                        collector->columns[chunk.file], out, *collector, *spill, format,
                                        values, names, coords, kinds);
                        if (collector->full(*spill)) collector->flush(*spill);
                    }
                    out.rows.push_back(row);
                }
            } catch (const std::exception &error) {
                fail(path + (chunk.begin ? " (line at byte " + std::to_string(reader.offset()) + ")"
                                         : ":" + std::to_string(number)) + ": " + error.what());
            }
            ++out.lines;
        }
    } catch (const std::exception &error) {
        out.error = error.what();
    }
}

string read_file(const string &path) {
    FILE *file = fopen(path.c_str(), "rb");
    if (!file) fail("cannot read " + path + ": " + strerror(errno));
    string data;
    struct stat info;
    if (fstat(fileno(file), &info) == 0) data.resize(info.st_size);
    if (!data.empty() && fread(&data[0], 1, data.size(), file) != data.size()) fail("read failed: " + path);
    fclose(file);
    return data;
}

// Sorted (hash, row) pairs: name -> row lookups without a string hash table.
struct NameIndex {
    vector<std::pair<uint64_t, uint32_t>> entries;
    const string *names = nullptr;
    const vector<Row> *rows = nullptr;

    void build(const vector<Row> &all, const string &arena) {
        rows = &all;
        names = &arena;
        entries.resize(all.size());
        for (size_t i = 0; i < all.size(); ++i)
            entries[i] = {hash_name(string_view(arena).substr(all[i].name, all[i].name_len)), uint32_t(i)};
        std::sort(entries.begin(), entries.end());
    }
    string_view name(uint32_t row) const {
        return string_view(*names).substr((*rows)[row].name, (*rows)[row].name_len);
    }
    int64_t find(string_view text) const {
        uint64_t key = hash_name(text);
        auto it = std::lower_bound(entries.begin(), entries.end(), std::make_pair(key, uint32_t(0)));
        for (; it != entries.end() && it->first == key; ++it)
            if (name(it->second) == text) return it->second;
        return -1;
    }
    // A row ID used twice.
    string duplicate() const {
        for (size_t i = 0; i < entries.size(); ++i)
            for (size_t j = i + 1; j < entries.size() && entries[j].first == entries[i].first; ++j)
                if (name(entries[i].second) == name(entries[j].second)) return string(name(entries[i].second));
        return "";
    }
};

// First line start at or after byte `at` of an open file (size when none).
uint64_t line_start(int fd, uint64_t at, uint64_t size) {
    if (at == 0) return 0;
    vector<char> block(1u << 20);
    for (uint64_t position = at - 1; position < size;) {
        ssize_t got = pread(fd, block.data(), block.size(), off_t(position));
        if (got <= 0) break;
        const void *newline = memchr(block.data(), '\n', size_t(got));
        if (newline) return position + ((const char *)newline - block.data()) + 1;
        position += uint64_t(got);
    }
    return size;
}

// Uncompressed inputs are cut at line starts into ranges of at most 256 MB
// (up to 4 per thread); a .gz input is one range.
vector<Chunk> plan_chunks(const vector<string> &inputs, unsigned threads) {
    uint64_t piece = 256ull << 20;
    if (const char *value = getenv("GRVCFGRAPH_CHUNK_BYTES"))  // tests: force small ranges
        if (!parse_u64(value, piece) || !piece) fail("GRVCFGRAPH_CHUNK_BYTES must be a positive integer");
    vector<Chunk> chunks;
    for (size_t f = 0; f < inputs.size(); ++f) {
        const string &path = inputs[f];
        struct stat info;
        bool gz = path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
        if (gz || stat(path.c_str(), &info) || !S_ISREG(info.st_mode) || info.st_size == 0) {
            chunks.push_back({f, 0, UINT64_MAX});
            continue;
        }
        uint64_t size = uint64_t(info.st_size);
        uint64_t count = std::min<uint64_t>((size + piece - 1) / piece, uint64_t(std::max(1u, threads)) * 4);
        int fd = open(path.c_str(), O_RDONLY);
        if (fd < 0) fail("cannot open " + path + ": " + strerror(errno));
        uint64_t previous = 0;
        for (uint64_t k = 1; k < count; ++k) {
            uint64_t boundary = line_start(fd, size * k / count, size);
            if (boundary > previous && boundary < size) {
                chunks.push_back({f, previous, boundary});
                previous = boundary;
            }
        }
        close(fd);
        chunks.push_back({f, previous, UINT64_MAX});
    }
    return chunks;
}

// Scan the inputs' ranges in parallel and join them in file order, numbering
// data lines across the inputs.
Scanned scan_inputs(const vector<string> &inputs, unsigned threads, ObsCollector *collector = nullptr) {
    if (inputs.empty()) fail("-v VCF required");
    vector<Chunk> chunks = plan_chunks(inputs, threads);
    if (chunks.size() > 65536) fail("too many scan ranges");
    unsigned workers = std::max(1u, std::min<unsigned>(threads, chunks.size()));
    if (collector) {
        // Which sample column of each input belongs to which GAF sample.
        for (auto &path : inputs) {
            vector<int64_t> columns;
            LineReader reader(path);
            string_view line;
            while (reader.next(line)) {
                if (blank_line(line)) continue;
                if (line[0] != '#') break;
                if (!starts_with(line, "#CHROM")) continue;
                auto fields = split(line, '\t');
                columns.assign(fields.size(), -1);
                for (size_t c = 9; c < fields.size(); ++c) {
                    auto it = collector->slot.find(fields[c]);
                    if (it != collector->slot.end()) columns[c] = int64_t(it->second);
                }
                break;
            }
            collector->columns.push_back(std::move(columns));
        }
        collector->start(workers);
    }
    vector<Scanned> parts(chunks.size());
    {
        std::atomic<size_t> next{0};
        vector<std::thread> pool;
        for (unsigned t = 0; t < workers; ++t)
            pool.emplace_back([&, t] {
                ObsSpill *spill = collector ? &collector->spills[t] : nullptr;
                for (size_t i; (i = next++) < chunks.size();)
                    scan_chunk(inputs[chunks[i].file], chunks[i], uint16_t(i), parts[i], collector, spill);
            });
        for (auto &thread : pool) thread.join();
    }
    for (auto &part : parts)
        if (!part.error.empty()) fail(part.error);
    if (collector) collector->finish();
    Scanned all;
    size_t total = 0;
    for (auto &part : parts) total += part.rows.size();
    if (total >= UINT32_MAX) fail("more than 2^32 rows");
    all.rows.reserve(total);
    for (size_t c = 0; c < parts.size(); ++c) {
        Scanned &part = parts[c];
        if (c == 0 || chunks[c].file != chunks[c - 1].file) all.file_starts.push_back(all.lines);
        all.chunk_row_base.push_back(all.rows.size());
        all.chunk_contigs.push_back(std::move(part.contigs));
        std::unordered_map<string, uint32_t>().swap(part.contig_ids);
        for (Row row : part.rows) {
            row.name += all.names.size();
            row.chrom += all.names.size();
            row.seq += (row.flags & SITE) ? all.names.size() : all.seq.size();
            row.line += all.lines;
            all.rows.push_back(row);
        }
        all.names += part.names;
        all.seq += part.seq;
        all.lines += part.lines;
        vector<Row>().swap(part.rows);
        string().swap(part.names);
        string().swap(part.seq);
    }
    log("scanned " + std::to_string(all.lines) + " data lines: " + std::to_string(all.rows.size()) + " rows, " +
        std::to_string(all.seq.size()) + " allele bases");
    return all;
}

// ------------------------------------------------------------------ FASTA paths

vector<string> split(string_view text, char separator) {
    vector<string> fields;
    size_t start = 0;
    for (;;) {
        size_t stop = text.find(separator, start);
        fields.emplace_back(text.substr(start, stop == string_view::npos ? string_view::npos : stop - start));
        if (stop == string_view::npos) return fields;
        start = stop + 1;
    }
}

vector<string> split_space(string_view text) {
    vector<string> fields;
    size_t i = 0;
    while (i < text.size()) {
        while (i < text.size() && isspace((unsigned char)text[i])) ++i;
        size_t start = i;
        while (i < text.size() && !isspace((unsigned char)text[i])) ++i;
        if (i > start) fields.emplace_back(text.substr(start, i - start));
    }
    return fields;
}

bool starts_with(string_view text, string_view prefix) { return text.substr(0, prefix.size()) == prefix; }

bool all_digits(string_view text) {
    if (text.empty()) return false;
    for (char c : text)
        if (c < '0' || c > '9') return false;
    return true;
}

string real_path(const string &path) {
    char *resolved = realpath(path.c_str(), nullptr);
    if (!resolved) return path;
    string out(resolved);
    free(resolved);
    return out;
}

// gfa_query_anchors.reverse_complement, IUPAC codes included.
string reverse_complement(string text) {
    static const string from = "ACGTRYMKBDHVNacgtrymkbdhvn", to = "TGCAYRKMVHDBNtgcayrkmvhdbn";
    std::reverse(text.begin(), text.end());
    for (char &c : text) {
        size_t at = from.find(c);
        if (at != string::npos) c = to[at];
    }
    return text;
}

struct FaiRecord {
    uint64_t length, offset, line_bases, line_width;
};

class Fasta {
  public:
    explicit Fasta(const string &path) : path_(path) {
        fd_ = open(path.c_str(), O_RDONLY);
        if (fd_ < 0) fail("cannot open " + path + ": " + strerror(errno));
        struct stat info;
        if (stat((path + ".fai").c_str(), &info)) fail(path + " requires an adjacent index " + path + ".fai (samtools faidx)");
        string data = read_file(path + ".fai");
        for (auto &line : split(data, '\n')) {
            auto fields = split(line, '\t');
            if (fields.size() < 5) continue;
            FaiRecord record;
            if (!parse_u64(fields[1], record.length) || !parse_u64(fields[2], record.offset) ||
                !parse_u64(fields[3], record.line_bases) || !parse_u64(fields[4], record.line_width) || !record.line_bases)
                fail(path + ".fai: malformed line");
            if (!index_.count(fields[0])) order_.push_back(fields[0]);
            index_[fields[0]] = record;
        }
    }
    ~Fasta() { close(fd_); }
    const string &path() const { return path_; }
    const vector<string> &names() const { return order_; }
    bool has(const string &name) const { return index_.count(name) > 0; }
    uint64_t length(const string &name) const { return record(name).length; }

    // Whitespace-separated tokens of the record's header line (read just before its bases).
    vector<string> header(const string &name) const {
        uint64_t offset = record(name).offset;
        for (uint64_t width = 4096;; width *= 2) {
            uint64_t start = offset > width ? offset - width : 0;
            string block = read(start, offset - start);
            while (!block.empty() && (block.back() == '\n' || block.back() == '\r')) block.pop_back();
            size_t marker = block.rfind("\n>");
            string header = marker != string::npos ? block.substr(marker + 2)
                            : (starts_with(block, ">") ? block.substr(1) : "");
            auto fields = split_space(header);
            if (!fields.empty() && fields[0] == name && header.find('\n') == string::npos) return fields;
            if (start == 0 || width >= (1u << 20)) fail(path_ + ": cannot read the header of " + name + "; check the .fai");
        }
    }

    // Bases [start, end) of one record.
    void fetch(const string &name, uint64_t start, uint64_t end, string &out) const {
        const FaiRecord &r = record(name);
        if (end > r.length || start > end) fail(path_ + ": " + name + ":" + std::to_string(start) + "-" + std::to_string(end) + " is out of bounds");
        out.clear();
        if (start == end) return;
        uint64_t first = r.offset + start / r.line_bases * r.line_width + start % r.line_bases;
        uint64_t last = r.offset + (end - 1) / r.line_bases * r.line_width + (end - 1) % r.line_bases + 1;
        string raw = read(first, last - first);
        out.reserve(end - start);
        for (char c : raw)
            if (c != '\n' && c != '\r') out.push_back(c);
        if (out.size() != end - start) fail(path_ + ": incomplete FASTA interval of " + name);
    }

    // Upper-case bases, reverse-complemented for strand '-'.
    string sequence(const string &name, uint64_t start, uint64_t end, char strand) const {
        string out;
        fetch(name, start, end, out);
        for (char &c : out) c = char(toupper((unsigned char)c));
        return strand == '-' ? reverse_complement(out) : out;
    }

  private:
    const FaiRecord &record(const string &name) const {
        auto it = index_.find(name);
        if (it == index_.end()) fail(path_ + ": record " + name + " is not in the .fai");
        return it->second;
    }
    string read(uint64_t offset, uint64_t size) const {
        string data(size, '\0');
        size_t done = 0;
        while (done < size) {
            ssize_t got = pread(fd_, &data[done], size - done, offset + done);
            if (got <= 0) fail("read failed: " + path_);
            done += got;
        }
        return data;
    }
    string path_;
    int fd_ = -1;
    std::unordered_map<string, FaiRecord> index_;
    vector<string> order_;
};

const uint32_t NONE = UINT32_MAX;

struct Lift {
    int64_t root = -1;  // backbone path the template end attaches to
    uint64_t point = 0;
    char strand = '+';
};

// A FASTA path, or a name whose bases are an interval of one (alias).
struct Root {
    string name, kind, path, record;  // record: the FASTA record holding the bases
    bool emit = false;
    uint64_t base = 0, length = 0;    // bases [base, base+length) of record
    int64_t order = -1;               // output order of emitted paths
    int64_t alias = -1;               // root index carrying this name's bases
    uint64_t alias_start = 0, alias_end = 0;
    char alias_strand = '+';
    string lift_status;
    Lift left, right;
    uint32_t sequence = NONE;
};

// gfa_source_catalog.SOURCE: SAMPLE:CONTIG:START-END[+-]
struct Source {
    string sample, contig;
    uint64_t start = 0, end = 0;
    char strand = '+';
};

bool parse_source(string value, Source &out) {
    out.strand = '+';
    if (!value.empty() && (value.back() == '+' || value.back() == '-')) {
        out.strand = value.back();
        value.pop_back();
    }
    size_t colon = value.rfind(':');
    if (colon == string::npos) return false;
    string range = value.substr(colon + 1), left = value.substr(0, colon);
    size_t dash = range.find('-');
    if (dash == string::npos || !all_digits(range.substr(0, dash)) || !all_digits(range.substr(dash + 1))) return false;
    size_t first = left.find(':');
    if (first == string::npos || first == 0 || first + 1 >= left.size()) return false;
    out.sample = left.substr(0, first);
    out.contig = left.substr(first + 1);
    return parse_u64(range.substr(0, dash), out.start) && parse_u64(range.substr(dash + 1), out.end);
}

std::unordered_map<string, string> header_attributes(const vector<string> &fields) {
    std::unordered_map<string, string> attrs;
    for (size_t i = 1; i < fields.size(); ++i) {
        size_t equals = fields[i].find('=');
        if (equals != string::npos) attrs[fields[i].substr(0, equals)] = fields[i].substr(equals + 1);
    }
    return attrs;
}

bool has_token(const vector<string> &fields, const string &token) {
    return std::find(fields.begin(), fields.end(), token) != fields.end();
}

// gfa_source_catalog.source_interval
bool source_interval(const vector<string> &fields, uint64_t length, Source &out) {
    string value;
    for (size_t i = 1; i < fields.size() && value.empty(); ++i)
        if (starts_with(fields[i], "source=")) value = fields[i].substr(7);
    auto attrs = header_attributes(fields);
    if (value.empty() && attrs.count("source_haplotype") && attrs.count("source_contig") &&
        attrs.count("source_start") && attrs.count("source_end") && attrs.count("source_strand"))
        value = attrs["source_haplotype"] + ":" + attrs["source_contig"] + ":" + attrs["source_start"] + "-" +
                attrs["source_end"] + attrs["source_strand"];
    if (value.empty() && fields.size() > 1) value = fields[1];
    if (value.empty() || !parse_source(value, out)) return false;
    if (out.end < out.start || out.end - out.start != length)
        fail(fields[0] + ": source interval length differs from FASTA length " + std::to_string(length));
    return true;
}

// assembly_contigs.accepted_contig
bool accepted_contig(const string &sample, string contig) {
    if (sample != "HG38_h1") return true;
    if (starts_with(contig, "HG38#1#")) contig = contig.substr(7);
    if (contig == "chrX" || contig == "chrY") return true;
    if (!starts_with(contig, "chr") || !all_digits(contig.substr(3))) return false;
    string number = contig.substr(3);
    if (number[0] == '0') return false;
    uint64_t value = 0;
    parse_u64(number, value);
    return value >= 1 && value <= 22;
}

// graphvcfmerge._safe_variant_token
string safe_token(const string &text) {
    if (text.empty()) return "locus";
    string out;
    bool replaced = false;
    for (char c : text) {
        bool ok = isalnum((unsigned char)c) || c == '_' || c == '.' || c == '-';
        if (ok) {
            out.push_back(c);
            replaced = false;
        } else if (!replaced) {
            out.push_back('_');
            replaced = true;
        }
    }
    return out;
}

class Roots {
  public:
    vector<Root> roots;
    std::unordered_map<string, uint32_t> index;

    Fasta &fasta(const string &path) {
        auto it = fastas_.find(path);
        if (it == fastas_.end()) it = fastas_.emplace(path, std::make_unique<Fasta>(path)).first;
        return *it->second;
    }
    bool has(const string &name) const { return index.count(name) > 0; }
    uint32_t add(Root root) {
        index[root.name] = uint32_t(roots.size());
        roots.push_back(std::move(root));
        return uint32_t(roots.size() - 1);
    }
    int64_t next_order() const {
        int64_t order = -1;
        for (auto &root : roots) order = std::max(order, root.order);
        return order + 1;
    }
    string bases(uint32_t root, uint64_t start, uint64_t end, char strand) {
        const Root &r = roots[root];
        return fasta(r.path).sequence(r.record, r.base + start, r.base + end, strand);
    }

    // gfa_interval_pipeline._index_roots (rGFA): -r and -a records.
    void index_inputs(const string &reference, const vector<string> &alternatives, string &backbone) {
        vector<std::pair<string, string>> inputs{{reference, "reference"}};
        for (auto &path : alternatives) inputs.push_back({path, "alternative"});
        std::unordered_set<string> seen;
        int64_t emitted = 0;
        for (auto &[value, input_kind] : inputs) {
            string path = real_path(value);
            if (!seen.insert(path + "\t" + input_kind).second) continue;
            Fasta &reader = fasta(path);
            for (auto &name : reader.names()) {
                uint64_t length = reader.length(name);
                vector<string> fields = length ? reader.header(name) : vector<string>();
                bool fixed = has_token(fields, "sequence_role=imported_alternative");
                Source source;
                bool sourced = !fields.empty() && source_interval(fields, length, source);
                bool is_reference = input_kind == "reference";
                if (is_reference && !fixed && sourced && backbone.empty()) backbone = source.sample;
                string sample = sourced ? source.sample : (is_reference ? backbone : "");
                string contig = sourced ? source.contig : name;
                if (!fixed && !accepted_contig(sample, contig)) continue;
                bool emit = length > 0;
                bool novel = starts_with(name, "novel");
                string kind = (input_kind == "alternative" && novel) ? "novel" : input_kind;
                if (fixed || (is_reference && sourced && source.sample != backbone)) kind = novel ? "novel" : "alternative";
                auto previous = index.find(name);
                if (previous != index.end()) {
                    Root &prior = roots[previous->second];
                    if (prior.path != path || prior.length != length) fail("FASTA root " + name + " is defined more than once");
                    if (emit && !prior.emit) {
                        prior.emit = true;
                        prior.kind = kind;
                        prior.order = emitted++;
                    }
                    continue;
                }
                Root root;
                root.name = root.record = name;
                root.kind = kind;
                root.path = path;
                root.length = length;
                root.emit = emit;
                root.order = emit ? emitted++ : -1;
                add(root);
            }
        }
    }

    // gfa_source_catalog.template_lift
    void template_lift(const vector<string> &fields, const string &backbone, Root *out) {
        auto attrs = header_attributes(fields);
        string status = attrs.count("lift_status") ? attrs["lift_status"] : "";
        if (attrs["backbone"].empty() ||
            (status != "mapped" && status != "one_sided" && status != "both" && status != "unmapped"))
            fail(fields[0] + ": fallback template must be lifted before export; run lift_local_templates.py");
        if (!backbone.empty() && attrs["backbone"] != backbone)
            fail(fields[0] + ": stale template lift to " + attrs["backbone"] + ", expected " + backbone);
        if (out) out->lift_status = status;
        if (status != "mapped" && status != "one_sided") return;
        // TARGET:START-END:STRAND
        string placement = attrs["reference"];
        size_t strand_colon = placement.rfind(':');
        bool ok = strand_colon != string::npos && strand_colon + 2 == placement.size() &&
                  (placement.back() == '+' || placement.back() == '-');
        string target, range;
        if (ok) {
            string head = placement.substr(0, strand_colon);
            size_t colon = head.rfind(':');
            ok = colon != string::npos && colon > 0;
            if (ok) {
                target = head.substr(0, colon);
                range = head.substr(colon + 1);
            }
        }
        size_t dash = range.find('-');
        uint64_t start = 0, end = 0;
        ok = ok && dash != string::npos && all_digits(range.substr(0, dash)) && all_digits(range.substr(dash + 1));
        if (!ok) fail(fields[0] + ": malformed lifted reference placement");
        if (!parse_u64(range.substr(0, dash), start) || !parse_u64(range.substr(dash + 1), end))
            fail(fields[0] + ": malformed lifted reference placement");
        char strand = placement.back();
        auto it = index.find(target);
        if (it == index.end() || roots[it->second].kind != "reference" || start > end || end > roots[it->second].length)
            fail(fields[0] + ": lift target " + target + ":" + std::to_string(start) + "-" + std::to_string(end) +
                 " is outside the selected backbone");
        Lift left{int64_t(it->second), strand == '+' ? start : end, strand};
        Lift right{int64_t(it->second), strand == '+' ? end : start, strand};
        if (status == "mapped") {
            if (out) out->left = left, out->right = right;
        } else if (start != end) {
            fail(fields[0] + ": one-sided lift must identify one breakpoint");
        } else if (attrs["lift_note"] == "one_side_up") {
            if (out) out->left = left;
        } else if (attrs["lift_note"] == "one_side_down") {
            if (out) out->right = right;
        } else {
            fail(fields[0] + ": one-sided lift lacks side metadata");
        }
    }

    // gfa_source_catalog.resolve_local_sources: templates get their own path or
    // alias a sequence-identical included interval; catalog names only alias.
    void resolve_local(const std::set<string> &targets, const vector<string> &templates,
                       const vector<string> &catalogs, const string &backbone) {
        std::set<string> needed;
        for (auto &name : targets)
            if (!has(name)) needed.insert(name);
        if ((needed.empty() || catalogs.empty()) && templates.empty()) return;
        std::unordered_map<string, vector<std::pair<uint32_t, Source>>> included;
        auto key = [](const string &sample, const string &contig) { return sample + '\x01' + contig; };
        for (uint32_t r = 0; r < roots.size(); ++r) {
            Source source;
            if (source_interval(fasta(roots[r].path).header(roots[r].name), roots[r].length, source)) {
                included[key(source.sample, source.contig)].push_back({r, source});
            } else if (!backbone.empty() && roots[r].kind == "reference" &&
                       included[key(backbone, roots[r].name)].empty()) {
                Source whole{backbone, roots[r].name, 0, roots[r].length, '+'};
                included[key(backbone, roots[r].name)].push_back({r, whole});
            }
        }
        vector<std::pair<string, bool>> roles;
        std::unordered_set<string> listed;
        for (auto &path : templates)
            if (listed.insert(path + "\t1").second) roles.push_back({path, true});
        for (auto &path : catalogs)
            if (listed.insert(path + "\t0").second) roles.push_back({path, false});
        int64_t order = next_order();
        uint64_t added = 0, reused = 0, aliased = 0;
        for (auto &[catalog, is_template] : roles) {
            Fasta &reader = fasta(catalog);
            vector<string> names;
            for (auto &name : reader.names())
                if (is_template || needed.count(name)) names.push_back(name);
            std::sort(names.begin(), names.end());
            names.erase(std::unique(names.begin(), names.end()), names.end());
            for (auto &name : names) {
                uint64_t length = reader.length(name);
                auto existing = index.find(name);
                if (existing != index.end()) {
                    if (is_template) template_lift(reader.header(name), backbone, nullptr);
                    const Root &prior = roots[existing->second];
                    uint32_t target = prior.alias >= 0 ? uint32_t(prior.alias) : existing->second;
                    uint64_t start = prior.alias >= 0 ? prior.alias_start : 0;
                    uint64_t end = prior.alias >= 0 ? prior.alias_end : prior.length;
                    char strand = prior.alias >= 0 ? prior.alias_strand : '+';
                    if (length != end - start ||
                        reader.sequence(name, 0, length, '+') != bases(target, start, end, strand))
                        fail(name + ": conflicting sequence definitions in FASTA catalogs");
                    continue;
                }
                vector<string> fields = reader.header(name);
                Source source;
                if (!source_interval(fields, length, source))
                    fail(name + ": local path lacks source SAMPLE:CONTIG:START-END metadata");
                if (!has_token(fields, "sequence_role=imported_alternative") && !accepted_contig(source.sample, source.contig)) {
                    if (is_template) continue;
                    fail(name + ": excluded HG38 non-primary source contig " + source.contig);
                }
                Root root;
                root.name = root.record = name;
                root.length = length;
                if (is_template) template_lift(fields, backbone, &root);
                vector<std::tuple<uint32_t, uint64_t, uint64_t, char>> matching;
                string local;
                for (auto &[target, origin] : included[key(source.sample, source.contig)]) {
                    if (!(origin.start <= source.start && source.end <= origin.end)) continue;
                    uint64_t a = origin.strand == '+' ? source.start - origin.start : origin.end - source.end;
                    uint64_t b = origin.strand == '+' ? source.end - origin.start : origin.end - source.start;
                    char strand = source.strand == origin.strand ? '+' : '-';
                    if (local.empty() && length) local = reader.sequence(name, 0, length, '+');
                    if (local == bases(target, a, b, strand)) {
                        auto candidate = std::make_tuple(target, a, b, strand);
                        if (std::find(matching.begin(), matching.end(), candidate) == matching.end())
                            matching.push_back(candidate);
                    }
                }
                if (matching.empty()) {
                    if (!is_template)
                        fail(name + ": no sequence-identical covering source in the -r/-a FASTAs or templates for " +
                             source.sample + ":" + source.contig + ":" + std::to_string(source.start) + "-" +
                             std::to_string(source.end) + source.strand +
                             "; the local catalog is lookup-only and cannot add this path");
                    root.kind = "alternative";
                    root.path = catalog;
                    root.emit = true;
                    root.order = order++;
                    uint32_t added_root = add(root);
                    included[key(source.sample, source.contig)].push_back({added_root, source});
                    ++added;
                    continue;
                }
                if (matching.size() != 1) fail(name + ": ambiguous source mapping to multiple included paths");
                auto [target, a, b, strand] = matching[0];
                root.alias = target;
                root.alias_start = a;
                root.alias_end = b;
                root.alias_strand = strand;
                root.path = roots[target].path;
                root.kind = is_template ? "alternative" : roots[target].kind;
                root.emit = is_template;
                root.order = is_template ? order++ : -1;
                add(root);
                ++(is_template ? reused : aliased);
            }
        }
        log("templates: " + std::to_string(added) + " own paths, " + std::to_string(reused) +
            " reusing included sequence; " + std::to_string(aliased) + " catalog names aliased");
    }

    // gfa_duplications.add_template_roots: DUP_<locus>_<start>_<end> copies [start, end) of locus.
    void add_templates(const std::set<string> &targets) {
        int64_t order = next_order();
        for (auto &name : targets) {
            string locus;
            uint64_t start = 0, end = 0;
            if (has(name) || !template_name(name, &locus, &start, &end)) continue;
            if (!has(locus)) {
                vector<string> known;
                for (auto &root : roots)
                    if (safe_token(root.name) == locus) known.push_back(root.name);
                std::sort(known.begin(), known.end());
                if (known.size() > 1) fail(name + ": locus " + locus + " matches several paths");
                if (known.empty())
                    fail(name + ": template locus " + locus + " is not in -r/-a, the local reference templates or the catalog");
                locus = known[0];
            }
            const Root &named = roots[index[locus]];
            uint32_t target = index[locus];
            uint64_t low = start, high = end;
            char strand = '+';
            if (named.alias >= 0) {
                target = uint32_t(named.alias);
                low = named.alias_strand == '+' ? named.alias_start + start : named.alias_end - end;
                high = named.alias_strand == '+' ? named.alias_start + end : named.alias_end - start;
                strand = named.alias_strand;
            }
            const Root base = roots[target];
            if (!(low < high && high <= base.length))
                fail(name + ": template interval " + locus + ":" + std::to_string(start) + "-" + std::to_string(end) +
                     " is not in -r/-a or the local reference templates");
            if (strand != '+') fail(name + ": template maps in reverse onto " + base.name + "; not supported");
            Root root;
            root.name = name;
            root.kind = "duplication";
            root.path = base.path;
            root.record = base.record;
            root.base = base.base + low;
            root.length = high - low;
            root.emit = true;
            root.order = order++;
            add(root);
        }
    }

    void write_table(const string &path) const {
        FILE *file = fopen(path.c_str(), "w");
        if (!file) fail("cannot write " + path);
        fprintf(file, "#name\tkind\temit\torder\tfasta\trecord\tbase\tlength\talias_target\talias_start\talias_end\t"
                      "alias_strand\tlift_status\tlift_left\tlift_right\n");
        auto lift = [&](const Lift &value) {
            return value.root < 0 ? string(".")
                                  : roots[value.root].name + ":" + std::to_string(value.point) + ":" + value.strand;
        };
        for (auto &root : roots) {
            fprintf(file, "%s\t%s\t%d\t%lld\t%s\t%s\t%llu\t%llu\t", root.name.c_str(), root.kind.c_str(), int(root.emit),
                    (long long)root.order, root.path.c_str(), root.record.c_str(), (unsigned long long)root.base,
                    (unsigned long long)root.length);
            if (root.alias >= 0)
                fprintf(file, "%s\t%llu\t%llu\t%c\t", roots[root.alias].name.c_str(), (unsigned long long)root.alias_start,
                        (unsigned long long)root.alias_end, root.alias_strand);
            else
                fprintf(file, ".\t.\t.\t.\t");
            fprintf(file, "%s\t%s\t%s\n", root.lift_status.empty() ? "." : root.lift_status.c_str(), lift(root.left).c_str(),
                    lift(root.right).c_str());
        }
        fclose(file);
    }

  private:
    std::unordered_map<string, std::unique_ptr<Fasta>> fastas_;
};

struct Options {
    vector<string> vcfs, alternatives, templates, catalogs, sample_vcfs;
    string reference, backbone, output, placements, roots_table, gaf_folder, query_lengths, spill_folder;
    unsigned threads = 1, gaf_threads = 0;  // gaf_threads 0: as --threads
    uint64_t max_node = 1024;
};

struct Sequence {
    bool is_row;
    uint32_t owner;   // root or row index
    uint64_t length;
    uint16_t rank;
};

struct Placement {
    uint32_t sequence = NONE;  // canonical parent sequence (NONE: inside a deletion's interval)
    uint32_t start = 0, end = 0;
    char strand = '+';
};

void parallel_for(size_t count, unsigned threads, const std::function<void(size_t)> &body) {
    std::atomic<size_t> next{0};
    std::mutex lock;
    string error;
    vector<std::thread> pool;
    for (unsigned t = 0; t < std::max(1u, threads); ++t)
        pool.emplace_back([&] {
            try {
                for (size_t i; (i = next++) < count;) body(i);
            } catch (const std::exception &e) {
                std::lock_guard<std::mutex> guard(lock);
                if (error.empty()) error = e.what();
                next = count;
            }
        });
    for (auto &thread : pool) thread.join();
    if (!error.empty()) fail(error);
}



// A buffered output file that is closed before its buffer is freed, also when
// a block fails part-way (fclose would otherwise flush freed memory at exit).
class Writer {
  public:
    explicit Writer(const string &path) : path_(path), buffer_(1u << 20) {
        file_ = fopen(path.c_str(), "wb");
        if (!file_) fail("cannot write " + path + ": " + strerror(errno));
        setvbuf(file_, buffer_.data(), _IOFBF, buffer_.size());
    }
    ~Writer() {
        if (file_) fclose(file_);
    }
    void write(const string &text) {
        if (fwrite(text.data(), 1, text.size(), file_) != text.size()) fail("write failed: " + path_);
    }
    void close() {
        FILE *file = file_;
        file_ = nullptr;
        if (fclose(file)) fail("write failed: " + path_);
    }

  private:
    string path_;
    vector<char> buffer_;
    FILE *file_ = nullptr;
};

void append_u64(string &out, uint64_t value) {
    char number[24];
    out.append(number, std::to_chars(number, number + sizeof(number), value).ptr);
}

// Sorted chunks in parallel, then pairwise merges.
template <class T, class Less>
void parallel_sort(vector<T> &items, unsigned threads, Less less) {
    size_t parts = std::min<size_t>(std::max(1u, threads), std::max<size_t>(1, items.size() / 100000));
    if (parts <= 1) {
        std::sort(items.begin(), items.end(), less);
        return;
    }
    vector<size_t> bounds(parts + 1);
    for (size_t i = 0; i <= parts; ++i) bounds[i] = items.size() * i / parts;
    parallel_for(parts, threads, [&](size_t i) { std::sort(items.begin() + bounds[i], items.begin() + bounds[i + 1], less); });
    for (size_t width = 1; width < parts; width *= 2) {
        vector<size_t> starts;
        for (size_t i = 0; i + width < parts; i += 2 * width) starts.push_back(i);
        parallel_for(starts.size(), threads, [&](size_t k) {
            size_t i = starts[k];
            std::inplace_merge(items.begin() + bounds[i], items.begin() + bounds[i + width],
                               items.begin() + bounds[std::min(parts, i + 2 * width)], less);
        });
    }
}

// ------------------------------------------------------------------ GAF
// One GAF per sample, as gfa_sample_gaf.write_sample_gaf builds it: the
// sample's ##pseudoLinearMapping location lines are chained into runs (query
// and reference contiguous on the same strands), the sample's alleles (its
// QUERYCOORD observations in the merged sample columns) are spliced into each
// run at their breakpoints (nested rows into their parent allele, a
// deletion's nested rows as its allele), and consecutive runs of a contig are
// joined into one record while the query is contiguous. Joins the graph does
// not link are kept and their links added to the GFA.

// A link between two handles; its rank sits in the low 16 bits under the
// target, so (from, key) orders by from, target, rank. 16 bytes.
struct Link {
    uint64_t from, key;
    static Link make(uint64_t from, uint64_t to, uint32_t rank) {
        if (to >> 48) fail("internal: handle beyond 2^48");
        if (rank != UINT32_MAX && rank >= 0xFFFF) fail("nesting rank beyond 65534");
        return {from, to << 16 | (rank == UINT32_MAX ? 0xFFFF : rank)};
    }
    uint64_t to() const { return key >> 16; }
    uint32_t rank() const { return (key & 0xFFFF) == 0xFFFF ? UINT32_MAX : uint32_t(key & 0xFFFF); }  // UINT32_MAX: added by the GAF walks
    bool operator<(const Link &other) const { return from != other.from ? from < other.from : key < other.key; }
    bool same(const Link &other) const { return from == other.from && to() == other.to(); }
};
static_assert(sizeof(Link) == 16, "Link layout");

const uint64_t NO_KEY = UINT64_MAX;

// One ##pseudoLinearMapping location line (Path=pri).
struct Mapping {
    string contig;
    uint64_t qs = 0, qe = 0;
    char qstrand = '+';
    string category;
    bool has_reference = false;
    string reference;
    uint64_t rlow = 0, rhigh = 0;
    char rstrand = '+';
};

// gfa_sample_gaf._META: key="value" pairs, backslash escapes decoded.
std::unordered_map<string, string> meta_values(string_view line) {
    std::unordered_map<string, string> values;
    size_t i = 0;
    while (i < line.size()) {
        if (!(isalnum((unsigned char)line[i]) || line[i] == '_')) {
            ++i;
            continue;
        }
        size_t key_start = i;
        while (i < line.size() && (isalnum((unsigned char)line[i]) || line[i] == '_')) ++i;
        if (i + 1 >= line.size() || line[i] != '=' || line[i + 1] != '"') continue;
        string key(line.substr(key_start, i - key_start)), value;
        i += 2;
        bool closed = false;
        while (i < line.size()) {
            char c = line[i++];
            if (c == '\\' && i < line.size()) {
                char e = line[i++];
                value.push_back(e == 'n' ? '\n' : e == 't' ? '\t' : e);
            } else if (c == '"') {
                closed = true;
                break;
            } else {
                value.push_back(c);
            }
        }
        if (closed) values[key] = value;
    }
    return values;
}

// gfa_sample_gaf._INTERVAL: NAME:START-END(+|-)
bool parse_interval(const string &text, string &name, uint64_t &low, uint64_t &high, char &strand) {
    if (text.size() < 5 || (text.back() != '+' && text.back() != '-')) return false;
    strand = text.back();
    string body = text.substr(0, text.size() - 1);
    size_t colon = body.rfind(':');
    if (colon == string::npos || colon == 0) return false;
    string range = body.substr(colon + 1);
    size_t dash = range.find('-');
    if (dash == string::npos || !all_digits(range.substr(0, dash)) || !all_digits(range.substr(dash + 1))) return false;
    name = body.substr(0, colon);
    return parse_u64(range.substr(0, dash), low) && parse_u64(range.substr(dash + 1), high);
}

// Mapping lines of every sample (gfa_sample_gaf.read_pseudolinear). A file with
// one sample column holds that sample's lines; otherwise (a merged header, or
// cohort.samples.headers.gz with SAMPLE<TAB>LINE rows) each line's Sample names it.
void read_mappings(const string &path, std::map<string, vector<Mapping>> &samples,
                   std::map<string, std::set<std::tuple<string, string, string>>> &seen) {
    LineReader reader(path);
    string_view line;
    vector<std::pair<string, string>> lines;  // (sample filter or "", line)
    string single;
    while (reader.next(line)) {
        if (starts_with(line, "#CHROM")) {
            auto columns = split(line, '\t');
            if (columns.size() == 10) single = columns[9];
            break;
        }
        size_t at = line.find("##pseudoLinearMapping=<");
        if (at == string_view::npos) {
            if (!line.empty() && line[0] != '#' && line.find('\t') == string_view::npos) break;
            continue;
        }
        lines.push_back({string(line.substr(0, at > 0 ? at - 1 : 0)), string(line.substr(at))});
    }
    if (!single.empty()) samples[single];  // its sample gets a GAF even without mapping lines
    for (auto &[prefix, text] : lines) {
        auto values = meta_values(text);
        string sample = !single.empty() ? single : !prefix.empty() ? prefix : values["Sample"];
        if (sample.empty() || sample == ".") continue;
        if (values["Path"] == "alt") continue;  // a sequence source places no bases
        Mapping m;
        if (!parse_interval(values["Query"], m.contig, m.qs, m.qe, m.qstrand)) continue;
        m.category = values.count("Category") ? values["Category"] : ".";
        // A deletion (query point) claims no bases: several at one point are
        // distinct reference spans, all kept; repeated intervals keep the first.
        string reference = values["Reference"];
        if (!seen[sample].insert({values["Query"], m.category, m.qs == m.qe ? reference : string()}).second) continue;
        m.has_reference = parse_interval(reference, m.reference, m.rlow, m.rhigh, m.rstrand);
        samples[sample].push_back(m);
    }
}

void sort_mappings(vector<Mapping> &intervals) {
    // Ties at one query point in walking order along the reference.
    auto tie = [](const Mapping &m) -> int64_t {
        if (!m.has_reference) return 0;
        return m.qstrand == m.rstrand ? int64_t(m.rlow) : -int64_t(m.rhigh);
    };
    std::stable_sort(intervals.begin(), intervals.end(), [&](const Mapping &a, const Mapping &b) {
        if (a.contig != b.contig) return a.contig < b.contig;
        if (a.qs != b.qs) return a.qs < b.qs;
        if (a.qe != b.qe) return a.qe < b.qe;
        return tie(a) < tie(b);
    });
}

// One carried allele of a sample: query interval, variant row (20 bytes;
// query contigs are shorter than 2^32).
struct Obs {
    uint32_t contig, row, qs, qe;
    uint8_t strand;
};

// Read-only graph and variant tables shared by the sample walkers.
struct GafGraph {
    const vector<Sequence> &sequences;
    const vector<uint64_t> &offsets, &unique, &node_first;
    const vector<uint32_t> &cuts, &gap_nodes;
    uint64_t max_node;
    const vector<Link> &links;
    const vector<Row> &rows;  // rows[i].line: VCF data line (variants.tsv vcf_index)
    // Per kept row (gfa_interval_pipeline._write_variant_index): parent and
    // allele keys (a sequence, or sequences.size() + row for a deletion;
    // UINT32_MAX: none), breakpoints on the parent, kind, and whether the
    // parent maps in reverse.
    vector<uint32_t> var_parent, var_path;
    vector<uint32_t> var_start, var_end;
    vector<uint8_t> var_kind, var_reverse;
    std::function<string(uint64_t)> key_name;

    struct Bounds {
        uint64_t node, start, end;
    };
    uint32_t sequence_of(uint64_t node) const {
        return uint32_t(std::upper_bound(node_first.begin(), node_first.end(), node) - node_first.begin() - 1);
    }
    // The segment holding position p (< length) of sequence s, and its extent.
    Bounds locate(uint32_t s, uint64_t p) const {
        const uint32_t *cut = cuts.data() + offsets[s];
        uint64_t k = std::upper_bound(cut, cut + unique[s], uint32_t(p)) - cut - 1;
        if (p >= sequences[s].length || k + 1 >= unique[s]) fail("internal: position outside its sequence");
        uint64_t j = (p - cut[k]) / max_node, start = cut[k] + j * max_node;
        return {node_first[s] + gap_nodes[offsets[s] + k] + j, start, std::min<uint64_t>(cut[k + 1], start + max_node)};
    }
    uint64_t node_length(uint64_t node) const {
        uint32_t s = sequence_of(node);
        uint64_t local = node - node_first[s];
        const uint32_t *first = gap_nodes.data() + offsets[s];
        uint64_t k = std::upper_bound(first, first + unique[s] - 1, uint32_t(local)) - first - 1;
        const uint32_t *cut = cuts.data() + offsets[s];
        uint64_t start = cut[k] + (local - first[k]) * max_node;
        return std::min<uint64_t>(cut[k + 1], start + max_node) - start;
    }
    uint64_t length(uint64_t key) const { return key < sequences.size() ? sequences[key].length : 0; }
    // Whether the GFA links handle a to handle b.
    bool linked(uint64_t a, uint64_t b) const {
        uint64_t from = a, to = b;
        if (std::make_pair(b ^ 1, a ^ 1) < std::make_pair(from, to)) from = b ^ 1, to = a ^ 1;
        if (!(from & 1) && to == from + 2 && sequence_of(from >> 1) == sequence_of(to >> 1)) return true;
        auto it = std::lower_bound(links.begin(), links.end(), std::make_pair(from, to),
                                   [](const Link &link, const std::pair<uint64_t, uint64_t> &key) {
                                       return std::make_pair(link.from, link.to()) < key;
                                   });
        return it != links.end() && it->from == from && it->to() == to;
    }
};

// Walk piece: consecutive segments of one sequence, first..last in walking
// order (descending and each reversed when reverse), with the bases trimmed
// off its first and last segment.
struct Piece {
    uint64_t first, last;
    bool reverse;
    uint32_t offset, trim, bases;  // within one sequence (< 2^32 bases)
};
static_assert(sizeof(Piece) == 32, "Piece layout");

Piece reverse_piece(const Piece &p) { return {p.last, p.first, !p.reverse, p.trim, p.offset, p.bases}; }

void reverse_pieces(vector<Piece> &pieces) {
    std::reverse(pieces.begin(), pieces.end());
    for (auto &piece : pieces) piece = reverse_piece(piece);
}

uint64_t walked_length(const vector<Piece> &pieces) {
    uint64_t total = 0;
    for (auto &p : pieces) total += p.bases - p.offset - p.trim;
    return total;
}

using Stats = std::map<string, int64_t>;

// allele_edges.edge_owned: which zero-length nested rows on an allele's query
// edge belong to it.
vector<size_t> edge_owned(vector<std::tuple<uint64_t, uint64_t, size_t>> starts,
                          vector<std::tuple<uint64_t, uint64_t, size_t>> ends,
                          vector<std::pair<uint64_t, uint64_t>> taken, int64_t current, int64_t span,
                          uint64_t length) {
    vector<size_t> owned;
    auto fits = [&](uint64_t first, uint64_t last) {
        if (current - int64_t(last - first) < span) return false;
        for (auto &[a, b] : taken)
            if (a < b && !(last <= a || b <= first)) return false;
        return true;
    };
    std::stable_sort(starts.begin(), starts.end(), [](auto &x, auto &y) {
        return std::make_pair(std::get<0>(x), std::get<1>(x)) < std::make_pair(std::get<0>(y), std::get<1>(y));
    });
    uint64_t covered = 0;
    for (auto &[first, last, key] : starts)
        if (first == covered && fits(first, last)) {
            owned.push_back(key);
            taken.push_back({first, last});
            current -= int64_t(last - first);
            covered = last;
        }
    std::stable_sort(ends.begin(), ends.end(), [](auto &x, auto &y) {
        return std::make_pair(std::get<1>(x), std::get<0>(x)) > std::make_pair(std::get<1>(y), std::get<0>(y));
    });
    covered = length;
    for (auto &[first, last, key] : ends)
        if (last == covered && fits(first, last)) {
            owned.push_back(key);
            taken.push_back({first, last});
            current -= int64_t(last - first);
            covered = first;
        }
    return owned;
}

enum { K_INS = 0, K_SNP = 1, K_SUB = 2, K_DEL = 3 };  // gfa_sample_gaf.KINDS
const char *GAF_KIND_NAMES[] = {"insertion", "snp", "substitution", "deletion"};

// One sample's alleles sorted by contig and query start (gfa_sample_gaf.Sample).
class SampleWalker {
  public:
    SampleWalker(const GafGraph &graph, vector<Obs> records, Stats &stats)
        : g(graph), rec(std::move(records)), used(rec.size(), 0), stats(stats) {
        std::stable_sort(rec.begin(), rec.end(), [](const Obs &a, const Obs &b) {
            return std::make_tuple(a.contig, a.qs, a.qe) < std::make_tuple(b.contig, b.qs, b.qe);
        });
    }
    const GafGraph &g;
    vector<Obs> rec;
    vector<uint8_t> used;
    std::unordered_set<size_t> skipped;  // positions skipped for overlapping another variant
    Stats &stats;

    std::pair<size_t, size_t> contig_range(uint32_t contig) const {
        auto low = std::lower_bound(rec.begin(), rec.end(), contig, [](const Obs &o, uint32_t c) { return o.contig < c; });
        auto high = std::upper_bound(rec.begin(), rec.end(), contig, [](uint32_t c, const Obs &o) { return c < o.contig; });
        return {size_t(low - rec.begin()), size_t(high - rec.begin())};
    }
    size_t first_at(size_t left, size_t right, uint64_t low) const {
        return std::lower_bound(rec.begin() + left, rec.begin() + right, low,
                                [](const Obs &o, uint64_t value) { return o.qs < value; }) - rec.begin();
    }
    // Unused records with query span inside [low, high] on contig.
    vector<size_t> between(uint32_t contig, uint64_t low, uint64_t high) const {
        auto [left, right] = contig_range(contig);
        vector<size_t> out;
        for (size_t p = first_at(left, right, low); p < right && rec[p].qs <= high; ++p)
            if (rec[p].qe <= high && !used[p]) out.push_back(p);
        return out;
    }
    // Used flags of records starting in [low, high] (restore()).
    std::pair<size_t, vector<uint8_t>> snapshot(uint32_t contig, uint64_t low, uint64_t high) const {
        auto [left, right] = contig_range(contig);
        size_t first = first_at(left, right, low);
        size_t last = std::upper_bound(rec.begin() + left, rec.begin() + right, high,
                                       [](uint64_t value, const Obs &o) { return value < o.qs; }) - rec.begin();
        last = std::max(first, last);  // a reversed interval (low > high) has no records
        return {first, vector<uint8_t>(used.begin() + first, used.begin() + last)};
    }
    void restore(const std::pair<size_t, vector<uint8_t>> &saved) {
        std::copy(saved.second.begin(), saved.second.end(), used.begin() + saved.first);
    }

    // Forward piece covering [low, high) of a sequence (none when empty).
    bool span(uint64_t key, uint64_t low, uint64_t high, vector<Piece> &out) const {
        if (high <= low) return false;
        if (key >= g.sequences.size() || high > g.sequences[key].length) fail("internal: span outside its sequence");
        auto a = g.locate(uint32_t(key), low), b = g.locate(uint32_t(key), high - 1);
        out.push_back({a.node, b.node, false, uint32_t(low - a.start), uint32_t(b.end - high), uint32_t(b.end - a.start)});
        return true;
    }

    // Nested candidates of one allele without the edge points it does not own.
    vector<size_t> edge_filter(const vector<size_t> &nested, uint64_t low, uint64_t high, bool descending,
                               bool has_length, uint64_t length, uint64_t path) const {
        if (!has_length || low == high) return nested;
        uint64_t first = descending ? high : low;
        vector<size_t> kept;
        vector<std::tuple<uint64_t, uint64_t, size_t>> starts, ends;
        vector<std::pair<uint64_t, uint64_t>> spans;
        int64_t current = int64_t(length);
        for (size_t candidate : nested) {
            uint64_t point = rec[candidate].qs;
            uint32_t variant = rec[candidate].row;
            bool child = g.var_parent[variant] == path;
            uint64_t a = g.var_start[variant], b = g.var_end[variant];
            if (child && point == rec[candidate].qe && (point == low || point == high)) {
                (point == first ? starts : ends).push_back({a, b, candidate});
                continue;
            }
            kept.push_back(candidate);
            if (child) {
                spans.push_back({a, b});
                current += int64_t(rec[candidate].qe - point) - int64_t(b - a);
            }
        }
        for (size_t key : edge_owned(starts, ends, spans, current, int64_t(high - low), length)) kept.push_back(key);
        return kept;
    }

    // Pieces for [low, high) of a sequence with the candidates on it spliced in.
    vector<Piece> splice(uint64_t path, uint64_t low, uint64_t high, const vector<size_t> &candidates, bool descending) {
        vector<std::tuple<uint64_t, uint64_t, int64_t, uint32_t>> chosen;
        for (size_t position : candidates) {
            uint32_t variant = rec[position].row;
            if (g.var_parent[variant] != path) continue;
            uint64_t a = g.var_start[variant], b = g.var_end[variant];
            if (low <= a && a <= b && b <= high)
                chosen.push_back({a, b, descending ? -int64_t(position) : int64_t(position), variant});
        }
        std::sort(chosen.begin(), chosen.end());
        // A SNP inside one of this sample's SVs here is redundant: the SV replaces that base.
        vector<uint64_t> sv_starts, sv_reach;
        for (auto &[a, b, position, variant] : chosen)
            if (b > a && g.var_kind[variant] != K_SNP) {
                sv_starts.push_back(a);
                sv_reach.push_back(std::max(b, sv_reach.empty() ? b : sv_reach.back()));
            }
        vector<Piece> pieces;
        uint64_t cursor = low;
        for (auto &[a, b, signed_position, variant] : chosen) {
            size_t position = size_t(signed_position < 0 ? -signed_position : signed_position);
            if (g.var_kind[variant] == K_SNP && !sv_starts.empty()) {
                int64_t covering = int64_t(std::lower_bound(sv_starts.begin(), sv_starts.end(), a + 1) - sv_starts.begin()) - 1;
                if (covering >= 0 && sv_reach[covering] >= b) {
                    used[position] = 1;
                    ++stats["snp_inside_sv"];
                    continue;
                }
            }
            if (a < cursor) {
                ++stats["variant_overlap_skipped"];
                skipped.insert(position);
                continue;
            }
            span(path, cursor, a, pieces);
            vector<Piece> inner = allele(position, variant);
            pieces.insert(pieces.end(), inner.begin(), inner.end());
            used[position] = 1;
            ++stats["variants_spliced"];
            cursor = b;
        }
        span(path, cursor, high, pieces);
        return pieces;
    }

    // Allele walk in its parent's orientation, with nested variants.
    vector<Piece> allele(size_t position, uint32_t variant) {
        uint64_t path = g.var_path[variant];
        uint64_t low = rec[position].qs, high = rec[position].qe;
        bool descending = rec[position].strand;
        bool deletion = g.var_kind[variant] == K_DEL;
        uint64_t length = deletion ? 0 : g.length(path);
        vector<size_t> nested;
        for (size_t candidate : between(rec[position].contig, low, high))
            if (candidate != position) nested.push_back(candidate);
        nested = edge_filter(nested, low, high, descending, !deletion, length, path);
        used[position] = 1;
        vector<Piece> pieces;
        if (deletion) {
            // A deletion has no bases; rows nested in it (--exact: a member's
            // kept bases) are its allele, by offset and then query order.
            vector<std::pair<uint64_t, int64_t>> chosen;
            for (size_t candidate : nested)
                if (g.var_parent[rec[candidate].row] == path)
                    chosen.push_back({g.var_start[rec[candidate].row], descending ? -int64_t(candidate) : int64_t(candidate)});
            std::sort(chosen.begin(), chosen.end());
            for (auto &[offset, signed_candidate] : chosen) {
                size_t candidate = size_t(signed_candidate < 0 ? -signed_candidate : signed_candidate);
                if (used[candidate]) continue;
                vector<Piece> inner = allele(candidate, rec[candidate].row);
                pieces.insert(pieces.end(), inner.begin(), inner.end());
                ++stats["variants_spliced"];
            }
        } else {
            pieces = splice(path, 0, length, nested, descending);
        }
        if (g.var_reverse[variant]) reverse_pieces(pieces);
        return pieces;
    }
};

// One GAF record under construction: segment ranges plus edge trims.
struct Walk {
    struct Part {
        uint64_t first, last;
        bool reverse;
        uint32_t bases;
    };
    uint64_t qstart = 0, qend = 0;
    vector<Part> parts;
    uint64_t last = NO_KEY;  // last handle; NO_KEY while empty
    uint64_t start_offset = 0, end_trim = 0, variants = 0;
    vector<std::pair<uint64_t, uint64_t>> pending;  // links added by the extend() in progress
};

using LinkSet = std::set<std::pair<uint64_t, uint64_t>>;

std::pair<uint64_t, uint64_t> link_key(uint64_t a, uint64_t b) {
    return std::min(std::make_pair(a, b), std::make_pair(b ^ 1, a ^ 1));
}

bool walk_append(const GafGraph &g, Walk &walk, const Piece &piece, Stats &stats, LinkSet &added) {
    uint64_t first = piece.first << 1 | piece.reverse, last = piece.last << 1 | piece.reverse;
    if (walk.last == NO_KEY) {
        walk.parts.push_back({piece.first, piece.last, piece.reverse, piece.bases});
        walk.last = last;
        walk.start_offset = piece.offset;
        walk.end_trim = piece.trim;
        return true;
    }
    uint64_t length = g.node_length(walk.last >> 1);
    if (walk.last == first && length - walk.end_trim == piece.offset) {
        // Continue inside the same segment: the piece without its first segment.
        if (piece.first != piece.last) {
            uint64_t next = piece.reverse ? piece.first - 1 : piece.first + 1;
            walk.parts.push_back({next, piece.last, piece.reverse, uint32_t(piece.bases - g.node_length(piece.first))});
            walk.last = last;
        }
        walk.end_trim = piece.trim;
        return true;
    }
    if (walk.end_trim || piece.offset) {
        ++stats["break_mid_node"];
        return false;
    }
    if (!g.linked(walk.last, first)) {
        auto key = link_key(walk.last, first);
        if (added.insert(key).second) {
            walk.pending.push_back(key);
            ++stats["link_added"];
        }
    }
    walk.parts.push_back({piece.first, piece.last, piece.reverse, piece.bases});
    walk.last = last;
    walk.end_trim = piece.trim;
    return true;
}

// Append all pieces or none (rolls back on a failed join).
bool walk_extend(const GafGraph &g, Walk &walk, const vector<Piece> &pieces, Stats &stats, LinkSet &added) {
    size_t count = walk.parts.size();
    uint64_t last = walk.last, start_offset = walk.start_offset, end_trim = walk.end_trim;
    walk.pending.clear();
    for (auto &piece : pieces)
        if (!walk_append(g, walk, piece, stats, added)) {
            walk.parts.resize(count);
            walk.last = last;
            walk.start_offset = start_offset;
            walk.end_trim = end_trim;
            for (auto &key : walk.pending) added.erase(key);
            stats["link_added"] -= int64_t(walk.pending.size());
            walk.pending.clear();
            return false;
        }
    return true;
}

struct GafInput {
    vector<Mapping> intervals;
    vector<Obs> records;
};

// Canonical reference of a mapping line: a sequence key, interval and strand.
struct Canon {
    bool ok = false;
    uint64_t key = 0, low = 0, high = 0;
    char strand = '+';
};

struct Run {
    string contig;
    uint64_t qs, qe;
    char qstrand;
    Canon run;  // ok: an anchored run
    string category;
    struct Member {
        uint64_t qs, qe;
        Canon canonical;
        string category;
    };
    vector<Member> members;
};

// gfa_sample_gaf._runs: contiguous pseudo-linear intervals on one sequence.
vector<Run> make_runs(const vector<Mapping> &intervals, const std::function<Canon(const Mapping &)> &canonical,
                      Stats &stats) {
    vector<Run> runs;
    std::optional<Run> current;
    for (auto &m : intervals) {
        ++stats["intervals"];
        Canon c;
        if (m.category == "UNMAPPED") {
            ++stats["break_unmapped"];
        } else if (m.has_reference) {
            c = canonical(m);
            if (!c.ok) ++stats["break_reference_not_in_gfa"];
        } else if (m.category != "INSERTION" && m.category != "DELETION") {
            ++stats["break_reference_not_in_gfa"];
        }
        Run::Member member{m.qs, m.qe, c, m.category};
        if (c.ok && current) {
            Canon &r = current->run;
            if (m.contig == current->contig && m.qstrand == current->qstrand && c.key == r.key &&
                c.strand == r.strand && m.qs == current->qe) {
                if (m.qstrand == c.strand && c.low == r.high) {
                    r.high = c.high;
                    current->qe = m.qe;
                    current->category = "RUN";
                    current->members.push_back(member);
                    continue;
                }
                if (m.qstrand != c.strand && c.high == r.low) {
                    r.low = c.low;
                    current->qe = m.qe;
                    current->category = "RUN";
                    current->members.push_back(member);
                    continue;
                }
            }
        }
        if (current) {
            ++stats["runs"];
            runs.push_back(*current);
            current.reset();
        }
        if (!c.ok) runs.push_back(Run{m.contig, m.qs, m.qe, m.qstrand, Canon(), m.category, {member}});
        else current = Run{m.contig, m.qs, m.qe, m.qstrand, c, m.category, {member}};
    }
    if (current) {
        ++stats["runs"];
        runs.push_back(*current);
    }
    return runs;
}

struct GafResult {
    Stats stats;
    LinkSet added;
};

// gfa_sample_gaf.write_sample_gaf for one sample.
GafResult write_sample_gaf(const GafGraph &g, const string &sample, GafInput input,
                           const vector<string> &contig_names,
                           const std::unordered_map<string, uint32_t> &contig_ids,
                           const std::unordered_map<string, uint64_t> &lengths,
                           const std::function<Canon(const Mapping &)> &canonical, const string &folder) {
    GafResult result;
    Stats &stats = result.stats;
    LinkSet &added = result.added;
    SampleWalker sw(g, std::move(input.records), stats);
    string base = folder + "/" + sample + ".gaf";
    Writer out(base + ".tmp"), breaks(base + ".breaks.tsv.tmp"), debug(base + ".lengthdebug.tsv.tmp");
    breaks.write("contig\tquery_start\tquery_end\tbases\tcategory\toutcome\treason\n");
    debug.write("contig\tquery_start\tquery_end\tcategory\trun_path\twalk_minus_query\tcandidates\t"
                "vcf_index:kind:start-end:query:state\n");
    int64_t written = 0;
    std::optional<Walk> walk;
    string current_contig;
    string line;

    auto count = [&](const char *key) -> int64_t {
        auto it = stats.find(key);
        return it == stats.end() ? 0 : it->second;
    };
    auto close = [&] {
        if (walk && !walk->parts.empty() && walk->qend > walk->qstart) {
            uint64_t path_length = 0;
            for (auto &part : walk->parts) path_length += part.bases;
            uint64_t path_start = walk->start_offset, path_end = path_length - walk->end_trim;
            uint64_t query_span = walk->qend - walk->qstart;
            uint64_t block = std::max(query_span, path_end - path_start), matches = std::min(query_span, path_end - path_start);
            auto known = lengths.find(current_contig);
            uint64_t length = known != lengths.end() ? known->second : walk->qend;
            line.assign(current_contig);
            line.push_back('\t');
            append_u64(line, length);
            line.push_back('\t');
            append_u64(line, walk->qstart);
            line.push_back('\t');
            append_u64(line, walk->qend);
            line.append("\t+\t");
            for (auto &part : walk->parts) {
                for (uint64_t n = part.first;; part.reverse ? --n : ++n) {
                    line.push_back(part.reverse ? '<' : '>');
                    append_u64(line, n);
                    if (n == part.last) break;
                }
                if (line.size() > (8u << 20)) {
                    out.write(line);
                    line.clear();
                }
            }
            line.push_back('\t');
            append_u64(line, path_length);
            line.push_back('\t');
            append_u64(line, path_start);
            line.push_back('\t');
            append_u64(line, path_end);
            line.push_back('\t');
            append_u64(line, matches);
            line.push_back('\t');
            append_u64(line, block);
            line.append("\t255\tnv:i:");
            append_u64(line, walk->variants);
            line.push_back('\n');
            out.write(line);
            ++written;
        }
        walk.reset();
    };

    auto report = [&](const string &contig, uint64_t qs, uint64_t qe, const string &category, const string &outcome,
                      const Stats &before) {
        string reason;
        for (auto &[key, value] : stats) {
            if (key.compare(0, 6, "break_") != 0) continue;
            auto it = before.find(key);
            if (value > (it == before.end() ? 0 : it->second)) reason += (reason.empty() ? "" : ",") + key;
        }
        breaks.write(contig + "\t" + std::to_string(qs) + "\t" + std::to_string(qe) + "\t" +
                     std::to_string(int64_t(qe) - int64_t(qs)) +
                     "\t" + category + "\t" + outcome + "\t" + (reason.empty() ? "." : reason) + "\n");
    };

    auto length_report = [&](uint64_t qs, uint64_t qe, const string &category, const Canon &run,
                             const vector<size_t> &candidates, int64_t delta) {
        string items;
        for (size_t position : candidates) {
            uint32_t variant = sw.rec[position].row;
            uint64_t parent = g.var_parent[variant], a = g.var_start[variant], b = g.var_end[variant];
            string state = sw.used[position] ? "placed"
                           : sw.skipped.count(position) ? "overlap_skipped"
                           : (!run.ok || parent != run.key) ? "other_path:" + g.key_name(parent)
                           : !(run.low <= a && a <= b && b <= run.high) ? "outside" : "unused";
            if (!items.empty()) items.push_back(';');
            items += std::to_string(g.rows[variant].line) + ":" + GAF_KIND_NAMES[g.var_kind[variant]] + ":" + std::to_string(a) + "-" +
                     std::to_string(b) + ":q" + std::to_string(sw.rec[position].qs) + "-" +
                     std::to_string(sw.rec[position].qe) + ":" + state;
        }
        debug.write(current_contig + "\t" + std::to_string(qs) + "\t" + std::to_string(qe) + "\t" + category + "\t" +
                    (run.ok ? g.key_name(run.key) : string(".")) + ":" + std::to_string(run.ok ? run.low : 0) + "-" +
                    std::to_string(run.ok ? run.high : 0) + "\t" + std::to_string(delta) + "\t" +
                    std::to_string(candidates.size()) + "\t" + items + "\n");
    };

    // The walk pieces of one run (or anchorless interval); false when it cannot be walked.
    auto place = [&](uint64_t qs, uint64_t qe, char qstrand, const Canon &run, const string &category,
                     int64_t contig_id, vector<Piece> &pieces) -> bool {
        vector<size_t> candidates = contig_id >= 0 ? sw.between(uint32_t(contig_id), qs, qe) : vector<size_t>();
        pieces.clear();
        bool ok = false;
        if (run.ok) {
            pieces = sw.splice(run.key, run.low, run.high, candidates, run.strand != qstrand);
            if (run.strand != qstrand) reverse_pieces(pieces);
            ok = true;
        } else if (category == "INSERTION" || category == "DELETION") {
            // No usable anchor: the variants carried inside this query interval.
            std::unordered_set<uint64_t> inner;
            for (size_t position : candidates) inner.insert(g.var_path[sw.rec[position].row]);
            for (size_t position : candidates) {
                uint32_t variant = sw.rec[position].row;
                if (sw.used[position] || inner.count(g.var_parent[variant])) continue;
                if (g.var_kind[variant] == K_DEL) continue;  // belongs to the anchored run next to it
                vector<Piece> allele = sw.allele(position, variant);
                if (sw.rec[position].strand) reverse_pieces(allele);
                pieces.insert(pieces.end(), allele.begin(), allele.end());
                ++stats["variants_spliced"];
            }
            ok = true;
            if (qe > qs && pieces.empty()) {
                ++stats["break_insertion_without_variant"];
                ok = false;
            }
        }
        if (ok && walked_length(pieces) != qe - qs) {
            // A variant was not placed here: never write a walk that does not spell the query.
            ++stats["break_length_mismatch"];
            length_report(qs, qe, category, run, candidates, int64_t(walked_length(pieces)) - int64_t(qe - qs));
            ok = false;
        }
        return ok;
    };

    auto walk_on = [&](uint64_t qs, uint64_t qe, const vector<Piece> &pieces, int64_t spliced) -> bool {
        // A record continues only where the query is contiguous.
        if (!walk || qs != walk->qend || !walk_extend(g, *walk, pieces, stats, added)) {
            Walk fresh;
            fresh.qstart = fresh.qend = qs;
            if (!walk_extend(g, fresh, pieces, stats, added)) return false;
            close();
            walk = std::move(fresh);
        }
        walk->qend = qe;
        walk->variants += uint64_t(spliced);
        ++stats["runs_walked"];
        return true;
    };

    vector<Piece> pieces;
    for (auto &run : make_runs(input.intervals, canonical, stats)) {
        if (run.contig != current_contig) {
            close();
            current_contig = run.contig;
        }
        auto found = contig_ids.find(run.contig);
        int64_t contig_id = found == contig_ids.end() ? -1 : int64_t(found->second);
        Stats saved = stats;
        std::pair<size_t, vector<uint8_t>> snapshot{0, {}};
        if (contig_id >= 0) snapshot = sw.snapshot(uint32_t(contig_id), run.qs, run.qe);
        int64_t spliced_start = count("variants_spliced");
        bool placed = place(run.qs, run.qe, run.qstrand, run.run, run.category, contig_id, pieces);
        if (placed && walk_on(run.qs, run.qe, pieces, count("variants_spliced") - spliced_start)) continue;
        report(run.contig, run.qs, run.qe, run.category, "run_failed", saved);
        // Variants consumed by the failed attempt stay available to the neighbouring runs.
        if (contig_id >= 0) sw.restore(snapshot);
        if (run.members.size() > 1) {
            // Retry interval by interval, so only intervals that cannot be walked are lost.
            stats = saved;
            if (contig_id >= 0) sw.restore(snapshot);
            ++stats["run_split"];
            for (auto &member : run.members) {
                spliced_start = count("variants_spliced");
                Stats before = stats;
                std::pair<size_t, vector<uint8_t>> member_snapshot{0, {}};
                if (contig_id >= 0) member_snapshot = sw.snapshot(uint32_t(contig_id), member.qs, member.qe);
                bool member_placed = place(member.qs, member.qe, run.qstrand, member.canonical, member.category,
                                           contig_id, pieces);
                if (!member_placed ||
                    !walk_on(member.qs, member.qe, pieces, count("variants_spliced") - spliced_start)) {
                    if (contig_id >= 0) sw.restore(member_snapshot);
                    stats["run_dropped"] += member_placed;
                    report(run.contig, member.qs, member.qe, member.category,
                           member_placed ? "interval_dropped" : "interval_unwalked", before);
                    close();
                }
            }
            continue;
        }
        stats["run_dropped"] += placed;
        close();
    }
    close();
    out.close();
    breaks.close();
    debug.close();
    for (const char *suffix : {"", ".breaks.tsv", ".lengthdebug.tsv"})
        if (rename((base + suffix + ".tmp").c_str(), (base + suffix).c_str())) fail("cannot rename " + base + suffix);
    int64_t unplaced = 0;
    Writer report_unplaced(base + ".unplaced.tsv.tmp");
    report_unplaced.write("contig\tquery_start\tquery_end\tstrand\tvcf_index\tparent\tparent_start\tparent_end\treason\n");
    for (size_t position = 0; position < sw.rec.size(); ++position) {
        if (sw.used[position]) continue;
        ++unplaced;
        const Obs &o = sw.rec[position];
        report_unplaced.write(contig_names[o.contig] + "\t" + std::to_string(o.qs) + "\t" + std::to_string(o.qe) + "\t" +
                              (o.strand ? "-" : "+") + "\t" + std::to_string(g.rows[o.row].line) + "\t" +
                              g.key_name(g.var_parent[o.row]) + "\t" + std::to_string(g.var_start[o.row]) + "\t" +
                              std::to_string(g.var_end[o.row]) + "\t" +
                              (sw.skipped.count(position) ? "overlap_skipped" : "unused") + "\n");
    }
    report_unplaced.close();
    if (rename((base + ".unplaced.tsv.tmp").c_str(), (base + ".unplaced.tsv").c_str())) fail("cannot rename unplaced report");
    stats["variants_unplaced"] = unplaced;
    stats["gaf_records"] = written;
    return result;
}

int build(Options &opt) {
    // GAF samples come first: the scan collects their alleles from the sample
    // columns as it reads the rows (one pass over the merged VCFs).
    vector<string> sample_names;
    vector<vector<Mapping>> sample_intervals;
    std::unique_ptr<ObsCollector> collector;
    if (!opt.gaf_folder.empty()) {
        if (mkdir(opt.gaf_folder.c_str(), 0777) && errno != EEXIST)
            fail("cannot create " + opt.gaf_folder + ": " + strerror(errno));
        // Each sample's mapping lines: per-sample VCFs or a headers file (-s),
        // else the merged VCF headers.
        std::map<string, vector<Mapping>> mapping;
        {
            std::map<string, std::set<std::tuple<string, string, string>>> seen;
            for (auto &path : opt.sample_vcfs.empty() ? opt.vcfs : opt.sample_vcfs) read_mappings(path, mapping, seen);
        }
        if (mapping.empty())
            fail("--gaf: no samples; give the per-sample VCFs or cohort.samples.headers.gz with -s");
        for (auto &[name, intervals] : mapping) {
            sort_mappings(intervals);
            sample_names.push_back(name);
            sample_intervals.push_back(std::move(intervals));
        }
        string spill = (opt.spill_folder.empty() ? opt.gaf_folder : opt.spill_folder) + "/.alleles.";
        collector = std::make_unique<ObsCollector>(spill, sample_names.size());
        for (size_t k = 0; k < sample_names.size(); ++k) collector->slot[sample_names[k]] = k;
    }
    Scanned scan = scan_inputs(opt.vcfs, opt.threads, collector.get());
    vector<Row> &rows = scan.rows;
    const string &names = scan.names;
    const string &allele_bases = scan.seq;
    string_view arena(names);
    auto row_name = [&](uint32_t i) { return arena.substr(rows[i].name, rows[i].name_len); };
    NameIndex index;
    index.build(rows, names);
    string duplicate = index.duplicate();
    if (!duplicate.empty()) fail("duplicate row ID " + duplicate);

    // Kept rows: PASS rows and the rows they are nested in
    // (gfa_interval_pipeline._reachable_events with include_parents).
    vector<int64_t> row_parent(rows.size());
    parallel_for(rows.size(), opt.threads, [&](size_t i) {
        row_parent[i] = index.find(arena.substr(rows[i].chrom, rows[i].chrom_len));
    });
    decltype(index.entries)().swap(index.entries);  // not needed any more
    vector<uint8_t> keep(rows.size(), 0);
    for (size_t i = 0; i < rows.size(); ++i)
        if (rows[i].flags & PASS)
            for (int64_t r = int64_t(i); r >= 0 && !keep[r]; r = row_parent[r]) keep[r] = 1;
    // FASTA path names the kept rows are placed on, and their DUP_ templates
    // (gfa_interval_pipeline._needed_targets): only these are resolved.
    std::set<string> targets;
    {
        uint64_t previous = UINT64_MAX;
        for (size_t i = 0; i < rows.size(); ++i) {
            if (!keep[i]) continue;
            if (row_parent[i] < 0 && rows[i].chrom != previous) {
                previous = rows[i].chrom;
                targets.emplace(arena.substr(rows[i].chrom, rows[i].chrom_len));
            }
            if (rows[i].flags & SITE) targets.emplace(arena.substr(rows[i].seq, rows[i].seq_len));
        }
    }
    Roots table;
    string backbone = opt.backbone;
    table.index_inputs(opt.reference, opt.alternatives, backbone);
    table.resolve_local(targets, opt.templates, opt.catalogs, backbone);
    table.add_templates(targets);
    if (!opt.roots_table.empty()) table.write_table(opt.roots_table);
    vector<Root> &roots = table.roots;
    const std::unordered_map<string, uint32_t> &root_index = table.index;
    // Emitted paths in output order.
    vector<uint32_t> root_order;
    for (uint32_t r = 0; r < roots.size(); ++r)
        if (roots[r].emit) root_order.push_back(r);
    std::stable_sort(root_order.begin(), root_order.end(), [&](uint32_t a, uint32_t b) { return roots[a].order < roots[b].order; });
    log("resolved " + std::to_string(targets.size()) + " target names; " + std::to_string(roots.size()) +
        " FASTA paths and aliases (backbone " + (backbone.empty() ? string("unknown") : backbone) + ")");
    const unsigned threads = opt.threads;
    const uint64_t max_node = opt.max_node;
    const string &output = opt.output;
    const string &placements_path = opt.placements;
    const string work = output + ".parts." + std::to_string(getpid());

    // Parents: a row (index) or a root (index | ROOT_BIT).
    const uint64_t ROOT_BIT = 1ull << 32, MISSING = UINT64_MAX;
    vector<uint64_t> parent(rows.size(), MISSING);
    vector<uint32_t> site(rows.size(), NONE);  // SITE rows: their DUP_ template root
    parallel_for(rows.size(), threads, [&](size_t i) {
        if (!keep[i]) return;
        if (row_parent[i] >= 0) {
            parent[i] = uint64_t(row_parent[i]);
        } else {
            auto it = root_index.find(string(arena.substr(rows[i].chrom, rows[i].chrom_len)));
            if (it != root_index.end()) parent[i] = ROOT_BIT | it->second;
        }
        if (rows[i].flags & SITE) {
            auto it = root_index.find(string(arena.substr(rows[i].seq, rows[i].seq_len)));
            if (it != root_index.end()) site[i] = it->second;
        }
    });
    vector<int64_t>().swap(row_parent);
    for (size_t r = 0; r < rows.size(); ++r) {
        if (!keep[r]) continue;
        if (parent[r] == MISSING)
            fail(string(row_name(r)) + ": graph target " + string(arena.substr(rows[r].chrom, rows[r].chrom_len)) +
                 " has no row or FASTA definition; supply it with -r/-a, the templates or the catalog");
        if ((rows[r].flags & SITE) && site[r] == NONE)
            fail(string(row_name(r)) + ": DUP_ template " + string(arena.substr(rows[r].seq, rows[r].seq_len)) +
                 " has no FASTA definition");
        if (root_index.count(string(row_name(r)))) fail(string(row_name(r)) + ": row ID also names a FASTA root");
        if (rows[r].flags & SITE) {
            const Root &tmpl = roots[site[r]];
            if (tmpl.alias >= 0) fail(tmpl.name + ": DUP_ template must have its own sequence");
            if (tmpl.length != rows[r].size)
                fail(string(row_name(r)) + ": SVLEN differs from the length of DUP_ template " + tmpl.name);
        }
    }
    for (auto &root : roots)
        if (root.length > UINT32_MAX) fail(root.name + ": paths of 2^32 bases or more are not supported");

    // Nesting level (1 on a root) and rank (SR: reference 0, other roots 1, rows parent + 1).
    auto canonical_root = [&](uint32_t r) { return roots[r].alias >= 0 ? uint32_t(roots[r].alias) : r; };
    auto root_rank = [&](uint32_t r) -> uint16_t { return roots[canonical_root(r)].kind == "reference" ? 0 : 1; };
    vector<uint16_t> level(rows.size(), 0), rank(rows.size(), 0);
    {
        vector<uint32_t> stack;
        for (size_t i = 0; i < rows.size(); ++i) {
            if (!keep[i] || level[i]) continue;
            // Walk up to a placed parent, then resolve from the outermost row inward.
            stack.clear();
            for (uint64_t r = i;;) {
                stack.push_back(uint32_t(r));
                if (stack.size() > rows.size()) fail("cyclic row/parent graph involving " + string(row_name(i)));
                uint64_t p = parent[r];
                if ((p & ROOT_BIT) || level[p]) break;
                r = p;
            }
            for (size_t k = stack.size(); k-- > 0;) {
                uint32_t row = stack[k];
                uint64_t p = parent[row];
                if (p & ROOT_BIT) {
                    level[row] = 1;
                    rank[row] = root_rank(uint32_t(p)) + 1;
                } else {
                    if (level[p] == UINT16_MAX) fail("nesting deeper than 65535 levels");
                    level[row] = level[p] + 1;
                    rank[row] = rank[p] + 1;
                }
            }
        }
    }

    // Sequences: emitted roots with their own bases, then row alleles by (level, VCF order).
    vector<Sequence> sequences;
    for (uint32_t r : root_order)
        if (roots[r].alias < 0 && roots[r].length) {
            roots[r].sequence = uint32_t(sequences.size());
            sequences.push_back({false, r, roots[r].length, root_rank(r)});
        }
    vector<uint32_t> own(rows.size(), NONE);
    {
        uint16_t deepest = 0;
        for (size_t i = 0; i < rows.size(); ++i)
            if (keep[i]) deepest = std::max(deepest, level[i]);
        vector<vector<uint32_t>> by_level(deepest + 1);
        for (size_t i = 0; i < rows.size(); ++i)
            if (keep[i] && rows[i].kind != DEL && !(rows[i].flags & SITE) && rows[i].size)
                by_level[level[i]].push_back(uint32_t(i));
        for (auto &group : by_level) {
            std::stable_sort(group.begin(), group.end(), [&](uint32_t a, uint32_t b) { return rows[a].line < rows[b].line; });
            for (uint32_t i : group) {
                own[i] = uint32_t(sequences.size());
                sequences.push_back({true, i, rows[i].size, rank[i]});
            }
        }
    }
    if (sequences.size() >= NONE) fail("too many sequences");
    log("levels resolved: " + std::to_string(sequences.size()) + " sequences with bases");

    // Placement of every kept row on its canonical parent sequence.
    auto coordinate_length = [&](uint64_t p) -> uint64_t {
        if (p & ROOT_BIT) return roots[uint32_t(p)].length;
        return rows[p].kind == DEL ? 0 : rows[p].size;
    };
    vector<Placement> place(rows.size());
    vector<uint8_t> on_deletion(rows.size(), 0);
    // Rows in level order so a deletion parent is placed before its children.
    vector<uint32_t> order;
    for (size_t i = 0; i < rows.size(); ++i)
        if (keep[i]) order.push_back(uint32_t(i));
    std::stable_sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) { return level[a] < level[b]; });
    for (uint32_t i : order) {
        const Row &row = rows[i];
        uint64_t p = parent[i];
        string name(row_name(i));
        bool parent_is_row = !(p & ROOT_BIT);
        int64_t start, end;
        if (row.kind == SNP) {
            start = int64_t(row.pos) - 1;
            end = row.pos;
            if (start < 0 || uint64_t(end) > coordinate_length(p))
                fail(name + ": SNP coordinate is outside its parent");
        } else if (parent_is_row && rows[p].kind == DEL) {
            // --exact: a member's kept bases replace the deletion's interval.
            if (row.kind != INS) fail(name + ": only insertions can be nested in a deletion");
            int64_t offset = int64_t(row.pos) - 1;
            int64_t deleted = rows[p].svlen ? rows[p].svlen : int64_t(rows[p].end) - rows[p].pos;
            if (offset < 0 || offset > std::max<int64_t>(deleted, int64_t(rows[p].end) - rows[p].pos))
                fail(name + ": offset is outside its deletion");
            place[i] = place[p];
            on_deletion[i] = 1;
            continue;
        } else {
            bool nested = parent_is_row || roots[uint32_t(p)].kind == "duplication";
            start = int64_t(row.pos) - (nested ? 1 : 0);
            end = row.end == row.pos ? start : int64_t(row.end);
            if (row.kind == DEL && row.svlen) end = start + row.svlen;
            if (start < 0 || end < start || uint64_t(end) > coordinate_length(p))
                fail(name + ": parent interval " + std::to_string(start) + "-" + std::to_string(end) +
                     " is outside 0-" + std::to_string(coordinate_length(p)));
        }
        Placement &out = place[i];
        if (parent_is_row) {
            if (rows[p].flags & SITE) {
                const Root &tmpl = roots[site[p]];
                if (tmpl.length != rows[p].size) fail(string(row_name(p)) + ": DUP_ template length differs from SVLEN");
                out.sequence = roots[canonical_root(site[p])].sequence;
                if (tmpl.alias >= 0) fail(tmpl.name + ": DUP_ template must have its own sequence");
            } else {
                out.sequence = own[p];
            }
            out.start = uint32_t(start);
            out.end = uint32_t(end);
        } else {
            const Root &root = roots[uint32_t(p)];
            if (root.alias >= 0) {
                if (root.alias_strand == '+') {
                    out.start = uint32_t(root.alias_start + start);
                    out.end = uint32_t(root.alias_start + end);
                } else {
                    out.start = uint32_t(root.alias_end - end);
                    out.end = uint32_t(root.alias_end - start);
                    out.strand = '-';
                }
                out.sequence = roots[root.alias].sequence;
            } else {
                out.sequence = root.sequence;
                out.start = uint32_t(start);
                out.end = uint32_t(end);
            }
        }
        if (out.sequence == NONE) fail(name + ": parent " + string(arena.substr(row.chrom, row.chrom_len)) + " has no emitted sequence");
    }
    log("placed " + std::to_string(order.size()) + " kept rows");

    // Cut points per sequence (CSR): both ends, row breakpoints, lift attachments.
    vector<uint64_t> offsets(sequences.size() + 1, 0);
    auto each_cut = [&](const std::function<void(uint32_t, uint64_t)> &emit) {
        for (uint32_t s = 0; s < sequences.size(); ++s) {
            emit(s, 0);
            emit(s, sequences[s].length);
        }
        for (uint32_t i : order) {
            if (on_deletion[i]) continue;
            emit(place[i].sequence, place[i].start);
            if (place[i].end != place[i].start) emit(place[i].sequence, place[i].end);
        }
        // An emitted alias (a template reusing included bases) walks that interval.
        for (auto &root : roots)
            if (root.emit && root.alias >= 0 && roots[root.alias].sequence != NONE) {
                emit(roots[root.alias].sequence, root.alias_start);
                emit(roots[root.alias].sequence, root.alias_end);
            }
        // Lifted template ends attach to the backbone there (lift targets are reference paths).
        for (auto &root : roots)
            for (const Lift *lift : {&root.left, &root.right}) {
                if (lift->root < 0) continue;
                const Root &named = roots[lift->root];
                uint64_t point = named.alias < 0 ? lift->point
                                 : named.alias_strand == '+' ? named.alias_start + lift->point : named.alias_end - lift->point;
                const Root &target = roots[canonical_root(uint32_t(lift->root))];
                if (target.sequence != NONE) emit(target.sequence, point);
            }
    };
    each_cut([&](uint32_t s, uint64_t) { ++offsets[s + 1]; });
    for (size_t s = 0; s < sequences.size(); ++s) offsets[s + 1] += offsets[s];
    vector<uint32_t> cuts(offsets.back());
    {
        vector<uint64_t> fill(offsets.begin(), offsets.end() - 1);
        each_cut([&](uint32_t s, uint64_t point) {
            if (point > sequences[s].length) fail("cut point outside a sequence");
            cuts[fill[s]++] = uint32_t(point);
        });
    }
    vector<uint64_t> unique(sequences.size()), node_first(sequences.size() + 1, 0);
    vector<uint32_t> gap_nodes(cuts.size());  // segments before each gap, within its sequence
    parallel_for(sequences.size(), threads, [&](size_t s) {
        uint32_t *begin = cuts.data() + offsets[s], *end = cuts.data() + offsets[s + 1];
        std::sort(begin, end);
        unique[s] = std::unique(begin, end) - begin;
        uint64_t count = 0;
        for (uint64_t k = 0; k + 1 < unique[s]; ++k) {
            gap_nodes[offsets[s] + k] = uint32_t(count);
            count += (uint64_t(begin[k + 1]) - begin[k] + max_node - 1) / max_node;
        }
        if (count > UINT32_MAX) fail("more than 2^32 segments in one sequence");
        node_first[s + 1] = count;
    });
    node_first[0] = 1;
    for (size_t s = 0; s < sequences.size(); ++s) node_first[s + 1] += node_first[s];
    uint64_t nodes = node_first.back() - 1;
    log("cut " + std::to_string(sequences.size()) + " sequences into " + std::to_string(nodes) + " segments");

    // The segment holding position p (< length) of sequence s.
    auto node_at = [&](uint32_t s, uint64_t p) -> uint64_t {
        const uint32_t *cut = cuts.data() + offsets[s];
        uint64_t k = std::upper_bound(cut, cut + unique[s], uint32_t(p)) - cut - 1;
        if (p >= sequences[s].length || k + 1 >= unique[s]) fail("internal: position outside its sequence");
        return node_first[s] + gap_nodes[offsets[s] + k] + (p - cut[k]) / max_node;
    };
    auto sequence_of = [&](uint64_t node) {
        return uint32_t(std::upper_bound(node_first.begin(), node_first.end(), node) - node_first.begin() - 1);
    };

    // ---------------------------------------------------------- links
    // A handle is an oriented segment, segment << 1 | reverse. Links follow
    // gfa_stable_coords: each allele is entered from the base before it on its
    // parent and left to the base after it; at an end of the parent the walk
    // continues through the parent's own parent (a deletion's nested rows use
    // the deletion's interval, a duplication site's rows the site's flanks),
    // and a lifted template's ends through its lift attachments.
    const uint64_t NO_HANDLE = UINT64_MAX;
    auto flip = [&](uint64_t handle) { return handle == NO_HANDLE ? handle : handle ^ 1; };
    std::function<uint64_t(uint32_t, uint64_t, bool)> outside;
    // Leaving root r on one side through its lift attachment (none: path end).
    auto through_lift = [&](uint32_t r, bool left) -> uint64_t {
        const Lift &lift = left ? roots[r].left : roots[r].right;
        if (lift.root < 0) return NO_HANDLE;
        const Root &named = roots[lift.root];
        uint64_t point = lift.point;
        char strand = lift.strand;
        if (named.alias >= 0) {
            point = named.alias_strand == '+' ? named.alias_start + point : named.alias_end - point;
            if (named.alias_strand == '-') strand = strand == '+' ? '-' : '+';
        }
        uint32_t s = roots[canonical_root(uint32_t(lift.root))].sequence;
        if (s == NONE) return NO_HANDLE;
        return strand == '+' ? outside(s, point, left) : flip(outside(s, point, !left));
    };
    // The segment just before (left) / at (right) position p of sequence s, walking it forward.
    outside = [&](uint32_t s, uint64_t p, bool left) -> uint64_t {
        if (left ? p > 0 : p < sequences[s].length) return node_at(s, left ? p - 1 : p) << 1;
        if (sequences[s].is_row) return NO_HANDLE;
        return through_lift(sequences[s].owner, left);
    };
    vector<Link> links;
    links.reserve(order.size() * 2);
    auto add_link = [&](uint64_t from, uint64_t to, uint32_t link_rank) {
        if (from == NO_HANDLE || to == NO_HANDLE) return;
        // A link and its reverse complement are one link: keep the smaller form.
        uint64_t reverse_from = to ^ 1, reverse_to = from ^ 1;
        if (std::make_pair(reverse_from, reverse_to) < std::make_pair(from, to)) from = reverse_from, to = reverse_to;
        // Consecutive segments of one sequence are written with the sequence.
        if (!(from & 1) && to == from + 2 && sequence_of(from >> 1) == sequence_of(to >> 1)) return;
        links.push_back(Link::make(from, to, link_rank));
    };
    // Rows nested in q sit on q's allele, or on its DUP_ template for a site row.
    auto allele_of = [&](uint32_t q) -> uint32_t {
        return (rows[q].flags & SITE) ? roots[site[q]].sequence : own[q];
    };
    auto path_rank = [&](uint32_t r) -> uint32_t { return roots[r].kind == "reference" ? 0 : 1; };
    {
        vector<uint64_t> anchor_left(rows.size(), NO_HANDLE), anchor_right(rows.size(), NO_HANDLE);
        for (uint32_t i : order) {  // parents before children
            const Placement &pl = place[i];
            uint64_t p = parent[i], left, right;
            if (p & ROOT_BIT) {
                uint32_t r = uint32_t(p);
                const Root &root = roots[r];
                // A template alias leaves through its own lift at its ends; a
                // lookup-only alias continues on the path it aliases.
                bool lifted = root.alias >= 0 && !root.lift_status.empty();
                uint64_t local_start = pl.start;
                if (root.alias >= 0)
                    local_start = root.alias_strand == '+' ? pl.start - root.alias_start : root.alias_end - pl.end;
                uint64_t local_end = local_start + (pl.end - pl.start);
                if (lifted && local_start == 0) left = through_lift(r, true);
                else left = pl.strand == '+' ? outside(pl.sequence, pl.start, true) : flip(outside(pl.sequence, pl.end, false));
                if (lifted && local_end == root.length) right = through_lift(r, false);
                else right = pl.strand == '+' ? outside(pl.sequence, pl.end, false) : flip(outside(pl.sequence, pl.start, true));
            } else if (rows[p].kind == DEL) {
                left = anchor_left[p];
                right = anchor_right[p];
            } else {
                uint32_t s = allele_of(uint32_t(p));
                left = pl.start > 0 ? node_at(s, pl.start - 1) << 1 : anchor_left[p];
                right = pl.end < sequences[s].length ? node_at(s, pl.end) << 1 : anchor_right[p];
            }
            anchor_left[i] = left;
            anchor_right[i] = right;
            uint32_t allele = rows[i].kind == DEL ? NONE : allele_of(i);
            if (allele == NONE) {
                add_link(left, right, rank[i]);
            } else {
                add_link(left, node_first[allele] << 1, rank[i]);
                add_link((node_first[allele + 1] - 1) << 1, right, rank[i]);
            }
        }
    }
    // A path's first and last handle (an alias walks its target interval).
    auto walk_ends = [&](uint32_t r, uint64_t &first, uint64_t &last) -> bool {
        const Root &root = roots[r];
        if (root.alias >= 0) {
            uint32_t t = roots[root.alias].sequence;
            if (t == NONE || root.alias_end <= root.alias_start) return false;
            uint64_t a = node_at(t, root.alias_start), b = node_at(t, root.alias_end - 1);
            first = root.alias_strand == '+' ? a << 1 : b << 1 | 1;
            last = root.alias_strand == '+' ? b << 1 : a << 1 | 1;
            return true;
        }
        if (root.sequence == NONE) return false;
        first = node_first[root.sequence] << 1;
        last = (node_first[root.sequence + 1] - 1) << 1;
        return true;
    };
    for (uint32_t r : root_order) {
        uint64_t first, last;
        if (roots[r].lift_status.empty() || !walk_ends(r, first, last)) continue;
        add_link(through_lift(r, true), first, path_rank(r));
        add_link(last, through_lift(r, false), path_rank(r));
    }
    parallel_sort(links, threads, [](const Link &a, const Link &b) { return a < b; });
    links.erase(std::unique(links.begin(), links.end(), [](const Link &a, const Link &b) { return a.same(b); }),
                links.end());
    uint64_t internal = nodes - sequences.size();
    log("links: " + std::to_string(internal) + " within sequences, " + std::to_string(links.size()) + " between them");

    // ---------------------------------------------------------- GAF
    if (!opt.gaf_folder.empty()) {
        GafGraph g{sequences, offsets, unique, node_first, cuts, gap_nodes, max_node, links, rows,
                   {}, {}, {}, {}, {}, {}, {}};
        // Variant tables (gfa_interval_pipeline._write_variant_index): a row
        // nested in another row sits on that row's allele (a site's DUP_
        // template; a deletion's own key, by offset); a row on a FASTA path on
        // its canonical sequence, reversed when that path aliases in reverse.
        const uint64_t deletion_base = sequences.size();
        if (deletion_base + rows.size() >= UINT32_MAX) fail("too many sequences and rows for the GAF tables");
        auto allele_key = [&](uint32_t q) -> uint32_t {
            if (rows[q].kind == DEL) return uint32_t(deletion_base + q);
            if (rows[q].flags & SITE) return roots[site[q]].sequence;
            return own[q];
        };
        size_t row_count = rows.size();
        g.var_parent.assign(row_count, UINT32_MAX);
        g.var_path.assign(row_count, UINT32_MAX);
        g.var_start.assign(row_count, 0);
        g.var_end.assign(row_count, 0);
        g.var_kind.assign(row_count, 0);
        g.var_reverse.assign(row_count, 0);
        for (uint32_t i : order) {
            const Row &row = rows[i];
            g.var_kind[i] = row.kind == INS ? K_INS : row.kind == SNP ? K_SNP : row.kind == SUB ? K_SUB : K_DEL;
            g.var_path[i] = allele_key(i);
            uint64_t p = parent[i];
            if (p & ROOT_BIT) {
                g.var_parent[i] = place[i].sequence;
                g.var_start[i] = place[i].start;
                g.var_end[i] = place[i].end;
                g.var_reverse[i] = place[i].strand == '-';
            } else if (rows[p].kind == DEL) {
                g.var_parent[i] = allele_key(uint32_t(p));
                g.var_start[i] = g.var_end[i] = row.pos - 1;  // rows nested in a deletion keep their offset
            } else {
                g.var_parent[i] = allele_key(uint32_t(p));
                g.var_start[i] = place[i].start;
                g.var_end[i] = place[i].end;
            }
        }
        g.key_name = [&](uint64_t key) -> string {
            if (key >= UINT32_MAX) return ".";
            if (key >= deletion_base) return string(row_name(uint32_t(key - deletion_base)));
            const Sequence &s = sequences[key];
            return s.is_row ? string(row_name(s.owner)) : roots[s.owner].name;
        };

        // The scan collected every sample's alleles (gfa_sample_gaf.build_shards);
        // only kept rows count, and a top-level SNP row's non-SNP observations
        // (INS_SNP restating an insertion) are skipped.
        vector<uint8_t> top_level(row_count, 0);
        for (uint32_t i : order) {
            string_view chrom = arena.substr(rows[i].chrom, rows[i].chrom_len);
            bool nested_name = false;
            for (const char *prefix : {"INS_", "DEL_", "SUB_", "DUP_", "I_", "D_"}) nested_name |= starts_with(chrom, prefix);
            uint64_t p = parent[i];
            top_level[i] = rows[i].kind == SNP && (p & ROOT_BIT) && roots[uint32_t(p)].alias < 0 && !nested_name;
        }
        // Query contigs numbered by first appearance over the scanned ranges in order.
        vector<string> contig_names;
        std::unordered_map<string, uint32_t> contig_ids;
        vector<vector<uint32_t>> contig_remap(scan.chunk_contigs.size());
        for (size_t c = 0; c < scan.chunk_contigs.size(); ++c)
            for (auto &name : scan.chunk_contigs[c]) {
                auto inserted = contig_ids.emplace(name, uint32_t(contig_names.size()));
                if (inserted.second) contig_names.push_back(name);
                contig_remap[c].push_back(inserted.first->second);
            }
        const std::unordered_map<string, size_t> &slot = collector->slot;

        // GAF column 2: --query-lengths (SAMPLE CONTIG LENGTH, or CONTIG LENGTH for
        // every sample), else the largest mapped query coordinate of the contig.
        vector<std::unordered_map<string, uint64_t>> lengths(sample_names.size());
        if (!opt.query_lengths.empty()) {
            LineReader reader(opt.query_lengths);
            string_view line;
            while (reader.next(line)) {
                if (line.empty() || line[0] == '#') continue;
                auto fields = split_space(line);
                uint64_t length = 0;
                if (fields.size() == 3 && slot.count(fields[0]) && parse_u64(fields[2], length)) {
                    lengths[slot.at(fields[0])][fields[1]] = length;
                } else if (fields.size() == 2 && parse_u64(fields[1], length)) {
                    for (auto &table : lengths) table[fields[0]] = length;
                } else if (fields.size() != 3) {
                    fail(opt.query_lengths + ": expected SAMPLE CONTIG LENGTH or CONTIG LENGTH");
                }
            }
        }

        // A mapping line's reference on its canonical sequence (an alias maps
        // onto the interval it reuses, reversed for a '-' alias).
        std::function<Canon(const Mapping &)> canonical = [&](const Mapping &m) -> Canon {
            Canon c;
            auto it = root_index.find(m.reference);
            if (it == root_index.end()) return c;
            const Root &root = roots[it->second];
            if (m.rlow > m.rhigh || m.rhigh > root.length) return c;
            if (root.alias >= 0) {
                const Root &target = roots[root.alias];
                if (target.sequence == NONE) return c;
                c.key = target.sequence;
                if (root.alias_strand == '+') {
                    c.low = root.alias_start + m.rlow;
                    c.high = root.alias_start + m.rhigh;
                    c.strand = m.rstrand;
                } else {
                    c.low = root.alias_end - m.rhigh;
                    c.high = root.alias_end - m.rlow;
                    c.strand = m.rstrand == '+' ? '-' : '+';
                }
            } else {
                if (root.sequence == NONE) return c;
                c.key = root.sequence;
                c.low = m.rlow;
                c.high = m.rhigh;
                c.strand = m.rstrand;
            }
            c.ok = true;
            return c;
        };

        vector<GafResult> results(sample_names.size());
        std::atomic<uint64_t> observations{0};
        // Samples walked at once: each holds about 400 bytes per carried allele.
        parallel_for(sample_names.size(), opt.gaf_threads ? opt.gaf_threads : threads, [&](size_t k) {
            // This sample's alleles: global rows and contigs, kept rows only, in
            // VCF order within one query interval.
            vector<Obs> records;
            {
                vector<DiskObs> raw = collector->load(k);
                records.reserve(raw.size());
                for (const DiskObs &d : raw) {
                    uint64_t row = scan.chunk_row_base[d.chunk] + d.row;
                    if (!keep[row] || ((d.flags & OBS_NOT_SNP_TYPE) && top_level[row])) continue;
                    if (d.flags & OBS_COORD_OVERFLOW)
                        fail("QUERYCOORD beyond 2^32 (sample " + sample_names[k] + ", row " + string(row_name(uint32_t(row))) + ")");
                    records.push_back({contig_remap[d.chunk][d.contig], uint32_t(row), d.qs, d.qe, d.strand});
                }
            }
            std::stable_sort(records.begin(), records.end(), [](const Obs &a, const Obs &b) {
                return std::make_tuple(a.contig, a.qs, a.qe, a.row) < std::make_tuple(b.contig, b.qs, b.qe, b.row);
            });
            observations += records.size();
            if (opt.query_lengths.empty()) {
                // GAF column 2: the largest mapped query coordinate of the contig.
                for (auto &m : sample_intervals[k]) lengths[k][m.contig] = std::max(lengths[k][m.contig], m.qe);
                for (auto &o : records)
                    lengths[k][contig_names[o.contig]] = std::max<uint64_t>(lengths[k][contig_names[o.contig]], o.qe);
            }
            GafInput input{std::move(sample_intervals[k]), std::move(records)};
            results[k] = write_sample_gaf(g, sample_names[k], std::move(input), contig_names, contig_ids, lengths[k],
                                          canonical, opt.gaf_folder);
        });
        collector.reset();  // removes the spill files
        log("GAF: " + std::to_string(sample_names.size()) + " samples, " + std::to_string(observations.load()) +
            " carried alleles on " + std::to_string(contig_names.size()) + " query contigs");
        // Links the walks need that the graph lacks become part of the graph.
        uint64_t before = links.size();
        for (auto &result : results)
            for (auto &[from, to] : result.added) links.push_back(Link::make(from, to, UINT32_MAX));
        parallel_sort(links, threads, [](const Link &a, const Link &b) { return a < b; });
        links.erase(std::unique(links.begin(), links.end(), [](const Link &a, const Link &b) { return a.same(b); }),
                    links.end());
        {
            std::set<string> keys;
            for (auto &result : results)
                for (auto &item : result.stats) keys.insert(item.first);
            Writer table(opt.gaf_folder + "/gaf_stats.tsv");
            string text = "sample";
            for (auto &key : keys) text += "\t" + key;
            text += "\n";
            for (size_t k = 0; k < sample_names.size(); ++k) {
                text += sample_names[k];
                for (auto &key : keys) {
                    auto it = results[k].stats.find(key);
                    text += "\t" + std::to_string(it == results[k].stats.end() ? 0 : it->second);
                }
                text += "\n";
            }
            table.write(text);
            table.close();
        }
        int64_t records = 0;
        for (auto &result : results) records += result.stats["gaf_records"];
        log("GAF: wrote " + std::to_string(records) + " records for " + std::to_string(sample_names.size()) +
            " samples to " + opt.gaf_folder + "; walks added " + std::to_string(links.size() - before) + " links");
    }

    // ---------------------------------------------------------- output
    // S, L and P lines are written in parallel blocks to temporary files,
    // then joined in order.
    for (auto &root : roots)
        if (root.sequence != NONE) table.fasta(root.path);
    vector<string> parts;
    const string temporary_output = output + ".tmp." + std::to_string(getpid());
    auto part_path = [&](char tag, size_t number) {
        char name[32];
        snprintf(name, sizeof(name), "/%c_part.%06zu", tag, number);
        return work + name;
    };
    struct Cleanup {
        std::function<void()> run;
        ~Cleanup() { run(); }
    } cleanup{[&] {
        for (auto &path : parts) unlink(path.c_str());
        rmdir(work.c_str());
        unlink(temporary_output.c_str());
    }};
    if (mkdir(work.c_str(), 0777) && errno != EEXIST) fail("cannot create " + work + ": " + strerror(errno));
    // Blocks of consecutive items with about equal weight.
    auto make_blocks = [&](size_t count, const std::function<uint64_t(size_t)> &weight) {
        uint64_t total = 0;
        for (size_t i = 0; i < count; ++i) total += weight(i);
        uint64_t target = std::max<uint64_t>(1, total / (uint64_t(std::max(1u, threads)) * 8));
        vector<std::pair<size_t, size_t>> blocks;
        size_t first = 0;
        uint64_t sum = 0;
        for (size_t i = 0; i < count; ++i) {
            sum += weight(i);
            if (sum >= target) {
                blocks.push_back({first, i + 1});
                first = i + 1;
                sum = 0;
            }
        }
        if (first < count) blocks.push_back({first, count});
        return blocks;
    };
    auto write_blocks = [&](char tag, const vector<std::pair<size_t, size_t>> &blocks,
                            const std::function<void(size_t, size_t, Writer &)> &body) {
        size_t base = parts.size();
        for (size_t b = 0; b < blocks.size(); ++b) parts.push_back(part_path(tag, b));
        parallel_for(blocks.size(), threads, [&](size_t b) {
            Writer writer(parts[base + b]);
            body(blocks[b].first, blocks[b].second, writer);
            writer.close();
        });
    };

    // S lines: whole sequences per block, balanced by bases.
    uint64_t bases = 0;
    for (auto &s : sequences) bases += s.length;
    const uint64_t WINDOW = 8u << 20;
    write_blocks('S', make_blocks(sequences.size(), [&](size_t s) { return sequences[s].length; }),
                 [&](size_t first_sequence, size_t last_sequence, Writer &writer) {
        string window, line;
        for (size_t s = first_sequence; s < last_sequence; ++s) {
            const Sequence &sequence = sequences[s];
            string name, kind;
            const Root *root = nullptr;
            const char *bases_of_row = nullptr;
            if (sequence.is_row) {
                const Row &row = rows[sequence.owner];
                name = string(row_name(sequence.owner));
                kind = KIND_NAMES[row.kind];
                if (row.seq + row.size > allele_bases.size()) fail(name + ": allele bases are outside the scan");
                bases_of_row = allele_bases.data() + row.seq;
            } else {
                root = &roots[sequence.owner];
                name = root->name;
                kind = root->kind;
            }
            const uint32_t *cut = cuts.data() + offsets[s];
            uint64_t id = node_first[s], window_start = 0;
            window.clear();
            for (uint64_t k = 0; k + 1 < unique[s]; ++k) {
                for (uint64_t low = cut[k]; low < cut[k + 1]; low += max_node) {
                    uint64_t high = std::min<uint64_t>(cut[k + 1], low + max_node);
                    const char *piece;
                    if (bases_of_row) {
                        piece = bases_of_row + low;
                    } else {
                        if (low < window_start || high > window_start + window.size()) {
                            window_start = low;
                            uint64_t stop = std::min<uint64_t>(sequence.length, std::max<uint64_t>(high, low + WINDOW));
                            table.fasta(root->path).fetch(root->record, root->base + low, root->base + stop, window);
                        }
                        piece = window.data() + (low - window_start);
                    }
                    line.assign("S\t");
                    append_u64(line, id);
                    line.push_back('\t');
                    size_t at = line.size();
                    line.append(piece, high - low);
                    for (size_t c = at; c < line.size(); ++c) {
                        char base = line[c] & ~0x20;  // upper case
                        if (base != 'A' && base != 'C' && base != 'G' && base != 'T' && base != 'N')
                            fail(name + ": sequence contains bases outside ACGTN; VG would change them on import");
                        line[c] = base;
                    }
                    line.append("\tLN:i:");
                    append_u64(line, high - low);
                    line.append("\tSN:Z:");
                    line.append(name);
                    line.append("\tSO:i:");
                    append_u64(line, low);
                    line.append("\tSR:i:");
                    append_u64(line, sequence.rank);
                    line.append("\tTP:Z:");
                    line.append(kind);
                    line.push_back('\n');
                    writer.write(line);
                    ++id;
                }
            }
            if (id != node_first[s + 1]) fail("internal: segment numbering mismatch for " + name);
        }
    });

    // L lines in segment order: consecutive segments of a sequence merged with
    // the links between sequences, both sorted by (from, to).
    const size_t node_blocks = std::max<size_t>(1, std::min<uint64_t>(nodes, uint64_t(std::max(1u, threads)) * 8));
    vector<std::pair<size_t, size_t>> link_blocks;
    for (size_t b = 0; b < node_blocks && nodes; ++b)
        link_blocks.push_back({1 + nodes * b / node_blocks, 1 + nodes * (b + 1) / node_blocks});
    write_blocks('L', link_blocks, [&](size_t low, size_t high, Writer &writer) {
        string line;
        auto emit = [&](uint64_t from, uint64_t to, uint32_t link_rank) {
            line.assign("L\t");
            append_u64(line, from >> 1);
            line.append(from & 1 ? "\t-\t" : "\t+\t");
            append_u64(line, to >> 1);
            line.append(to & 1 ? "\t-\t0M" : "\t+\t0M");
            if (link_rank != UINT32_MAX) {  // walk-added links have no rank
                line.append("\tSR:i:");
                append_u64(line, link_rank);
            }
            line.push_back('\n');
            writer.write(line);
        };
        auto next = std::lower_bound(links.begin(), links.end(), uint64_t(low) << 1,
                                     [](const Link &link, uint64_t from) { return link.from < from; });
        auto stop = std::lower_bound(links.begin(), links.end(), uint64_t(high) << 1,
                                     [](const Link &link, uint64_t from) { return link.from < from; });
        uint32_t s = sequence_of(low);
        for (uint64_t n = low; n < high; ++n) {
            while (n >= node_first[s + 1]) ++s;
            // Links from n+ that sort before n+ -> (n+1)+, then that link, then n- links.
            uint64_t from = n << 1;
            while (next != stop && (next->from < from || (next->from == from && next->to() < from + 2)))
                emit(next->from, next->to(), next->rank()), ++next;
            if (n + 1 < node_first[s + 1]) emit(from, from + 2, sequences[s].rank);
            while (next != stop && next->from <= (from | 1)) emit(next->from, next->to(), next->rank()), ++next;
        }
        while (next != stop) emit(next->from, next->to(), next->rank()), ++next;
    });

    // P lines for the emitted paths, in output order.
    auto path_segments = [&](uint32_t r) -> uint64_t {
        uint64_t first, last;
        if (!walk_ends(r, first, last)) return 0;
        return (std::max(first, last) >> 1) - (std::min(first, last) >> 1) + 1;
    };
    write_blocks('P', make_blocks(root_order.size(), [&](size_t k) { return path_segments(root_order[k]) + 1; }),
                 [&](size_t first_path, size_t last_path, Writer &writer) {
        string line;
        for (size_t k = first_path; k < last_path; ++k) {
            uint32_t r = root_order[k];
            const Root &root = roots[r];
            uint64_t first, last;
            if (!walk_ends(r, first, last)) continue;  // an empty record has no segments
            line.assign("P\t");
            line.append(root.name);
            line.push_back('\t');
            bool reverse = first & 1;
            uint64_t a = first >> 1, b = last >> 1;
            for (uint64_t n = a;; reverse ? --n : ++n) {
                append_u64(line, n);
                line.append(reverse ? "-" : "+");
                if (n == b) break;
                line.push_back(',');
                if (line.size() >= (8u << 20)) {
                    writer.write(line);
                    line.clear();
                }
            }
            line.append("\t*\tTP:Z:");
            line.append(root.kind);
            if (!root.lift_status.empty()) {
                line.append("\tLS:Z:");
                line.append(root.lift_status);
                if (root.lift_status != "mapped" && root.lift_status != "one_sided") line.append("\tUP:Z:unplaced");
            }
            line.push_back('\n');
            writer.write(line);
        }
    });
    {
        const string &temporary = temporary_output;
        Writer out(temporary);
        out.write("H\tVN:Z:1.0\tTS:Z:GrvcfGraph\n");
        vector<char> buffer(16u << 20);
        for (auto &path : parts) {
            FILE *in = fopen(path.c_str(), "rb");
            if (!in) fail("cannot read " + path);
            for (size_t got; (got = fread(buffer.data(), 1, buffer.size(), in)) > 0;)
                out.write(string(buffer.data(), got));
            bool failed = ferror(in);
            fclose(in);
            if (failed) fail("read failed: " + path);
            unlink(path.c_str());
        }
        out.close();
        if (rename(temporary.c_str(), output.c_str())) fail("cannot rename " + temporary);
    }
    log("wrote " + std::to_string(nodes) + " segments (" + std::to_string(bases) + " bases), " +
        std::to_string(internal + links.size()) + " links and " + std::to_string(root_order.size()) + " paths to " + output);

    if (!placements_path.empty()) {
        Writer file(placements_path);
        file.write("id\tkind\tlevel\trank\tparent\tsequence\tstart\tend\tstrand\tinside_deletion\tallele\n");
        string line;
        auto sequence_name = [&](uint32_t s) -> string {
            if (s == NONE) return ".";
            return sequences[s].is_row ? string(row_name(sequences[s].owner)) : roots[sequences[s].owner].name;
        };
        for (uint32_t i : order) {
            const Row &row = rows[i];
            string allele = own[i] != NONE ? sequence_name(own[i])
                            : (row.flags & SITE) ? roots[site[i]].name : ".";
            line.assign(row_name(i));
            line += string("\t") + KIND_NAMES[row.kind] + "\t" + std::to_string(level[i]) + "\t" + std::to_string(rank[i]) +
                    "\t" + string(arena.substr(row.chrom, row.chrom_len)) + "\t" + sequence_name(place[i].sequence) + "\t" +
                    std::to_string(place[i].start) + "\t" + std::to_string(place[i].end) + "\t" + place[i].strand + "\t" +
                    std::to_string(int(on_deletion[i])) + "\t" + allele + "\n";
            file.write(line);
        }
        file.close();
    }
    return 0;
}


void usage() {
    fprintf(stderr,
            "usage: GrvcfGraph -v VCF [VCF ...] -r BACKBONE.fa [-a ALT.fa ...]\n"
            "                  [--local-reference-templates FASTA ...] [--local-path-fasta FASTA ...]\n"
            "                  [--reference-haplotype NAME] -o OUT.gfa [-t THREADS]\n"
            "                  [--max-node-length 1024] [--placements FILE] [--roots-table FILE]\n"
            "                  [--gaf FOLDER [-s SAMPLE_VCF_OR_HEADERS ...] [--query-lengths TSV]\n"
            "                   [--spill-dir DIR] [--gaf-threads N]]\n");
}

}  // namespace

int main(int argc, char **argv) {
    Options opt;
    try {
        for (int i = 1; i < argc; ++i) {
            string option = argv[i];
            auto value = [&]() -> string {
                if (i + 1 >= argc) fail(option + " requires a value");
                return argv[++i];
            };
            auto values = [&](vector<string> &out) {
                size_t before = out.size();
                while (i + 1 < argc && argv[i + 1][0] != '-') out.push_back(argv[++i]);
                if (out.size() == before) fail(option + " requires a value");
            };
            if (option == "-h" || option == "--help") {
                usage();
                return 0;
            } else if (option == "-v" || option == "--vcf") {
                values(opt.vcfs);
            } else if (option == "-r" || option == "--reference-fasta") {
                opt.reference = value();
            } else if (option == "-a" || option == "--alternatives-fasta") {
                values(opt.alternatives);
            } else if (option == "--local-reference-templates") {
                values(opt.templates);
            } else if (option == "--local-path-fasta") {
                values(opt.catalogs);
            } else if (option == "--reference-haplotype") {
                opt.backbone = value();
            } else if (option == "-o" || option == "--output") {
                opt.output = value();
            } else if (option == "--placements") {
                opt.placements = value();
            } else if (option == "--roots-table") {
                opt.roots_table = value();
            } else if (option == "--gaf") {
                opt.gaf_folder = value();
            } else if (option == "-s" || option == "--sample-vcf") {
                values(opt.sample_vcfs);
            } else if (option == "--query-lengths") {
                opt.query_lengths = value();
            } else if (option == "--spill-dir") {
                opt.spill_folder = value();
            } else if (option == "--gaf-threads") {
                opt.gaf_threads = unsigned(std::max(1, atoi(value().c_str())));
            } else if (option == "-t" || option == "--threads") {
                opt.threads = unsigned(std::max(1, atoi(value().c_str())));
            } else if (option == "--max-node-length") {
                if (!parse_u64(value(), opt.max_node) || !opt.max_node || opt.max_node > UINT32_MAX)
                    fail("--max-node-length must be between 1 and 2^32-1");
            } else {
                fail("unknown option " + option);
            }
        }
        if (opt.vcfs.empty() || opt.reference.empty() || opt.output.empty()) {
            usage();
            return 2;
        }
        return build(opt);
    } catch (const std::exception &error) {
        fprintf(stderr, "GrvcfGraph: error: %s\n", error.what());
        return 1;
    }
}
