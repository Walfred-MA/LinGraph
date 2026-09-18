#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <condition_variable>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/types.h>
#include <unistd.h>

namespace {

constexpr std::size_t KMER_SIZE = 31;
constexpr std::size_t TOP_SHIFT = 2 * (KMER_SIZE - 1);
constexpr std::uint64_t FORWARD_MASK =
    (std::uint64_t{1} << (2 * (KMER_SIZE - 1))) - 1;
constexpr std::uint64_t DEFAULT_QUERY_MERGE_GAP = 100;
constexpr std::int64_t MIN_DELETE_SCORE = 100;
constexpr std::uint64_t MIN_RESIDUAL_SIZE = 100;
constexpr std::uint32_t MIN_MASKED_REFERENCE_RUN = 100;
constexpr std::uint64_t MASKED_REGION_HARD_GAP = 50;
constexpr std::size_t FASTA_WRAP = 80;

struct Options {
    std::string input_path;
    std::string input_list_path;
    std::string reference;
    std::string output_path;
    std::string output_list_path;
    std::string covered_bed_path;
    std::size_t threads = 8;
    std::uint64_t query_merge_gap = DEFAULT_QUERY_MERGE_GAP;
    std::uint64_t minimum_unmasked_bases = 0;
    bool reference_given = false;
    bool exclude_all = false;
    bool mask_covered = false;
};

struct FastaEntry {
    std::string name;
    std::string header;
    std::uint64_t length = 0;
    std::uint64_t offset = 0;
    std::uint64_t line_bases = 0;
    std::uint64_t line_width = 0;
};

struct Interval {
    std::uint64_t start = 0;
    std::uint64_t end = 0;
};

class KmerSearcherPresenceIndex {
public:
    void initialize(std::uint64_t maximum_items) {
        if (initialized_) {
            throw std::runtime_error(
                "internal error: KmerSearcher index initialized twice");
        }
        initialized_ = true;
        // KmerSearcher uses 30 bucket bits, a uint32 suffix, a flat suffix
        // array, and one byte of bucket capacity. Preserve that layout for
        // whole genomes. Smaller inputs use fewer buckets and uint64 suffixes
        // so --exclude-all does not have a fixed 5-GiB minimum.
        bucket_bits_ = 8;
        while (bucket_bits_ < 30 &&
               (maximum_items / 8 + (maximum_items % 8 != 0)) >
                   (std::uint64_t{1} << bucket_bits_)) {
            ++bucket_bits_;
        }
        if (bucket_bits_ < 30) {
            const __uint128_t adaptive_bytes =
                static_cast<__uint128_t>(maximum_items) * 8 +
                (static_cast<__uint128_t>(5) << bucket_bits_);
            const __uint128_t kmer_searcher_bytes =
                static_cast<__uint128_t>(maximum_items) * 4 +
                (static_cast<__uint128_t>(5) << 30);
            if (kmer_searcher_bytes < adaptive_bytes) bucket_bits_ = 30;
        }
        bucket_count_ = std::uint64_t{1} << bucket_bits_;
        bucket_mask_ = bucket_count_ - 1;
        bucket_capacities_.assign(
            static_cast<std::size_t>(bucket_count_), 0);
    }

    void preview(std::uint64_t kmer) {
        if (!initialized_ || allocated_) {
            throw std::runtime_error(
                "internal error: invalid KmerSearcher preview phase");
        }
        std::uint8_t& capacity = bucket_capacities_[bucket_for(kmer)];
        if (capacity != std::numeric_limits<std::uint8_t>::max()) {
            ++capacity;
        }
    }

    void allocate_suffixes() {
        if (!initialized_ || allocated_) {
            throw std::runtime_error(
                "internal error: invalid KmerSearcher allocation phase");
        }
        allocated_ = true;
        std::uint64_t slot_count = 0;
        for (const std::uint8_t capacity : bucket_capacities_) {
            if (capacity > std::numeric_limits<std::uint64_t>::max() -
                    slot_count) {
                throw std::runtime_error(
                    "reference KmerSearcher suffix count exceeds UINT64_MAX");
            }
            slot_count += capacity;
        }
        slot_count_ = slot_count;
        use_64_bit_offsets_ =
            slot_count > std::numeric_limits<std::uint32_t>::max();
        if (use_64_bit_offsets_) {
            bucket_offsets_64_.resize(
                static_cast<std::size_t>(bucket_count_));
        } else {
            bucket_offsets_32_.resize(
                static_cast<std::size_t>(bucket_count_));
        }

        std::uint64_t offset = 0;
        for (std::uint64_t bucket = 0; bucket < bucket_count_; ++bucket) {
            if (use_64_bit_offsets_) {
                bucket_offsets_64_[static_cast<std::size_t>(bucket)] = offset;
            } else {
                bucket_offsets_32_[static_cast<std::size_t>(bucket)] =
                    static_cast<std::uint32_t>(offset);
            }
            offset += bucket_capacities_[static_cast<std::size_t>(bucket)];
        }

        if (bucket_bits_ == 30) {
            if (slot_count > suffixes_32_.max_size()) {
                throw std::runtime_error(
                    "reference KmerSearcher suffix array is too large");
            }
            suffixes_32_.assign(static_cast<std::size_t>(slot_count), 0);
        } else {
            if (slot_count > suffixes_64_.max_size()) {
                throw std::runtime_error(
                    "reference KmerSearcher suffix array is too large");
            }
            suffixes_64_.assign(static_cast<std::size_t>(slot_count), 0);
        }
    }

    void insert(std::uint64_t kmer) {
        if (!allocated_ || finalized_) {
            throw std::runtime_error(
                "internal error: invalid KmerSearcher insertion phase");
        }
        const std::uint64_t bucket = bucket_for(kmer);
        const std::uint64_t suffix = suffix_for(kmer);
        if (suffix == 0) {
            throw std::runtime_error(
                "internal error: canonical 31-mer has a zero suffix");
        }
        const std::uint64_t start = offset_for(bucket);
        const std::uint64_t end = start +
            bucket_capacities_[static_cast<std::size_t>(bucket)];
        for (std::uint64_t index = start; index < end; ++index) {
            std::uint64_t stored = bucket_bits_ == 30
                ? suffixes_32_[static_cast<std::size_t>(index)]
                : suffixes_64_[static_cast<std::size_t>(index)];
            if (stored == suffix) return;
            if (stored == 0) {
                if (bucket_bits_ == 30) {
                    suffixes_32_[static_cast<std::size_t>(index)] =
                        static_cast<std::uint32_t>(suffix);
                } else {
                    suffixes_64_[static_cast<std::size_t>(index)] = suffix;
                }
                ++size_;
                return;
            }
        }
        // The original KmerSearcher silently drops a 256th distinct suffix.
        // Preserve exact cleaner semantics by retaining such keys separately.
        if (overflow_.insert(kmer).second) ++size_;
    }

    void finalize() {
        if (!allocated_ || finalized_) {
            throw std::runtime_error(
                "internal error: invalid KmerSearcher finalize phase");
        }
        for (std::uint64_t bucket = 0; bucket < bucket_count_; ++bucket) {
            const std::uint8_t capacity =
                bucket_capacities_[static_cast<std::size_t>(bucket)];
            if (capacity < 8) continue;
            const std::uint64_t start = offset_for(bucket);
            if (bucket_bits_ == 30) {
                auto begin = suffixes_32_.begin() +
                    static_cast<std::ptrdiff_t>(start);
                std::sort(begin, begin + capacity);
            } else {
                auto begin = suffixes_64_.begin() +
                    static_cast<std::ptrdiff_t>(start);
                std::sort(begin, begin + capacity);
            }
        }
        finalized_ = true;
    }

    bool contains(std::uint64_t kmer) const {
        if (!finalized_) return false;
        const std::uint64_t bucket = bucket_for(kmer);
        const std::uint64_t suffix = suffix_for(kmer);
        const std::uint8_t capacity =
            bucket_capacities_[static_cast<std::size_t>(bucket)];
        const std::uint64_t start = offset_for(bucket);
        bool found = false;
        if (capacity < 8) {
            for (std::uint64_t index = start;
                 index < start + capacity; ++index) {
                const std::uint64_t stored = bucket_bits_ == 30
                    ? suffixes_32_[static_cast<std::size_t>(index)]
                    : suffixes_64_[static_cast<std::size_t>(index)];
                if (stored == suffix) {
                    found = true;
                    break;
                }
            }
        } else if (bucket_bits_ == 30) {
            const auto begin = suffixes_32_.begin() +
                static_cast<std::ptrdiff_t>(start);
            found = std::binary_search(
                begin, begin + capacity, static_cast<std::uint32_t>(suffix));
        } else {
            const auto begin = suffixes_64_.begin() +
                static_cast<std::ptrdiff_t>(start);
            found = std::binary_search(begin, begin + capacity, suffix);
        }
        if (found || overflow_.empty()) return found;
        return overflow_.find(kmer) != overflow_.end();
    }

    std::uint64_t size() const { return size_; }
    std::uint64_t slot_count() const { return slot_count_; }
    std::uint64_t overflow_size() const { return overflow_.size(); }
    std::uint8_t bucket_bits() const { return bucket_bits_; }

    long double allocated_gib() const {
        const __uint128_t bytes =
            static_cast<__uint128_t>(bucket_capacities_.size()) +
            static_cast<__uint128_t>(bucket_offsets_32_.size()) *
                sizeof(std::uint32_t) +
            static_cast<__uint128_t>(bucket_offsets_64_.size()) *
                sizeof(std::uint64_t) +
            static_cast<__uint128_t>(suffixes_32_.size()) *
                sizeof(std::uint32_t) +
            static_cast<__uint128_t>(suffixes_64_.size()) *
                sizeof(std::uint64_t);
        return static_cast<long double>(bytes) /
            (1024.0L * 1024.0L * 1024.0L);
    }

private:
    bool initialized_ = false;
    bool allocated_ = false;
    bool finalized_ = false;
    bool use_64_bit_offsets_ = false;
    std::uint8_t bucket_bits_ = 0;
    std::uint64_t bucket_count_ = 0;
    std::uint64_t bucket_mask_ = 0;
    std::uint64_t slot_count_ = 0;
    std::uint64_t size_ = 0;
    std::vector<std::uint8_t> bucket_capacities_;
    std::vector<std::uint32_t> bucket_offsets_32_;
    std::vector<std::uint64_t> bucket_offsets_64_;
    std::vector<std::uint32_t> suffixes_32_;
    std::vector<std::uint64_t> suffixes_64_;
    std::unordered_set<std::uint64_t> overflow_;

    static std::uint64_t mix64(std::uint64_t value) {
        value += 0x9e3779b97f4a7c15ULL;
        value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
        value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
        return value ^ (value >> 31);
    }

    static std::uint32_t rotate_right_32(
        std::uint32_t value,
        unsigned shift
    ) {
        return (value >> shift) | (value << ((32 - shift) & 31));
    }

    static std::uint32_t mix_32_to_30(std::uint32_t value) {
        return (rotate_right_32(value, 7) ^
                rotate_right_32(value, 13) ^
                rotate_right_32(value, 19)) & ((1U << 30) - 1);
    }

    std::uint64_t suffix_for(std::uint64_t kmer) const {
        return kmer >> bucket_bits_;
    }

    std::uint64_t bucket_for(std::uint64_t kmer) const {
        const std::uint64_t suffix = suffix_for(kmer);
        const std::uint64_t low = kmer & bucket_mask_;
        if (bucket_bits_ == 30) {
            return low ^ mix_32_to_30(static_cast<std::uint32_t>(suffix));
        }
        return low ^ (mix64(suffix) & bucket_mask_);
    }

    std::uint64_t offset_for(std::uint64_t bucket) const {
        return use_64_bit_offsets_
            ? bucket_offsets_64_[static_cast<std::size_t>(bucket)]
            : bucket_offsets_32_[static_cast<std::size_t>(bucket)];
    }
};

struct ReferenceKmers {
    // One low-copy reference k-mer may remove every copy of the same
    // high-copy query k-mer, so counts are recorded but lookup uses presence.
    std::unordered_map<std::uint64_t, std::uint32_t> counts;
    // Sparse side table: the longest fully lowercase A/C/G/T run containing
    // this canonical k-mer in any reference record.
    std::unordered_map<std::uint64_t, std::uint32_t> masked_run_lengths;
    // --exclude-all needs presence only. Use KmerSearcher's compact
    // bucket/flat-suffix layout rather than a node-based count map.
    KmerSearcherPresenceIndex searcher_presence;
    bool use_searcher_presence = false;
    std::uint64_t total_occurrences = 0;

    bool contains(std::uint64_t kmer) const {
        return use_searcher_presence
            ? searcher_presence.contains(kmer)
            : counts.find(kmer) != counts.end();
    }

    std::uint64_t distinct_size() const {
        return use_searcher_presence
            ? searcher_presence.size()
            : counts.size();
    }
};

struct CleanResult {
    std::string sequence;
    std::vector<Interval> covered_intervals;
    std::vector<Interval> fragments;
    std::uint64_t deleted_bases = 0;
    std::uint64_t ambiguous_bases = 0;
    std::uint64_t short_bases = 0;
    std::uint64_t masked_bases = 0;
    std::uint64_t unmasked_bases = 0;
    std::uint64_t hit_intervals = 0;
    std::uint64_t deletion_intervals = 0;
};

struct Totals {
    std::uint64_t query_records = 0;
    std::uint64_t query_bases = 0;
    std::uint64_t output_records = 0;
    std::uint64_t output_bases = 0;
    std::uint64_t deleted_bases = 0;
    std::uint64_t ambiguous_bases = 0;
    std::uint64_t short_bases = 0;
    std::uint64_t masked_bases = 0;
    std::uint64_t unmasked_bases = 0;
    std::uint64_t filtered_records = 0;
    std::uint64_t filtered_bases = 0;
};

struct SampleJob {
    std::string sample;
    std::string input_path;
    std::string fai_path;
    std::string output_path;
};

struct ReferenceSummary {
    std::uint64_t records = 0;
    std::uint64_t bases = 0;
};

void usage(const char* program) {
    std::cerr
        << "Usage:\n  " << program
        << " -i input.fa -o output.fa [-r reference.fa|name1,name2] "
           "[-t threads] [--query-merge-gap INT] [--exclude-all] "
           "[--covered-bed hits.bed] [--mask-covered "
           "--min-unmasked-bases INT]\n  "
        << program
        << " -I Sample.list -O Output.list [-r reference.fa|name1,name2] "
           "[-t threads] [--query-merge-gap INT] [--exclude-all]\n\n"
        << "Remove query sequence represented by canonical 31-mers in one "
           "fixed reference.\n"
        << "Both input.fa and an external reference.fa require adjacent .fai "
           "indexes.\n"
        << "Without -r, the first input record is the reference. A non-path -r "
           "is a\n"
        << "comma-separated list of exact input FASTA record names (the first "
           "header\n"
        << "field). The first occurrence of each requested name is selected.\n"
        << "Sample.list rows are SAMPLE FASTA [FAI]; Output.list rows are "
           "SAMPLE FASTA.\n"
        << "--exclude-all disables masked-repeat eligibility and uses every "
           "reference k-mer.\n"
        << "--covered-bed writes the raw union of matching query 31-mer "
           "windows as\n0-based, half-open BED3 intervals (single -i/-o "
           "mode only).\n"
        << "--mask-covered lowercases those raw covered intervals instead of "
           "deleting\nthem; --min-unmasked-bases filters whole records by "
           "remaining uppercase\nACGT bases (default: 0).\n"
        << "--query-merge-gap is the hard query hit-gap break (default: 100).\n"
        << "Coordinates in generated fragment names are 0-based, half-open.\n";
}

std::string require_value(int argc, char** argv, int& index) {
    if (index + 1 >= argc) {
        throw std::runtime_error(
            std::string("missing value after ") + argv[index]);
    }
    return argv[++index];
}

std::uint64_t parse_uint64(
    const std::string& text,
    const std::string& context
) {
    if (text.empty() || text.front() == '-') {
        throw std::runtime_error(context + ": invalid nonnegative integer: " + text);
    }
    std::size_t consumed = 0;
    std::uint64_t value = 0;
    try {
        value = std::stoull(text, &consumed);
    } catch (...) {
        throw std::runtime_error(context + ": invalid nonnegative integer: " + text);
    }
    if (consumed != text.size()) {
        throw std::runtime_error(context + ": invalid nonnegative integer: " + text);
    }
    return value;
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "-h" || argument == "--help") {
            usage(argv[0]);
            std::exit(EXIT_SUCCESS);
        } else if (argument == "-i" || argument == "--input") {
            options.input_path = require_value(argc, argv, index);
        } else if (argument == "-I" || argument == "--input-list") {
            options.input_list_path = require_value(argc, argv, index);
        } else if (argument == "-r" || argument == "--reference") {
            if (options.reference_given) {
                throw std::runtime_error("-r/--reference may be supplied only once");
            }
            options.reference_given = true;
            options.reference = require_value(argc, argv, index);
            if (options.reference.empty()) {
                throw std::runtime_error("-r/--reference requires a file or name");
            }
        } else if (argument == "-o" || argument == "--output") {
            options.output_path = require_value(argc, argv, index);
        } else if (argument == "-O" || argument == "--output-list") {
            options.output_list_path = require_value(argc, argv, index);
        } else if (argument == "--exclude-all") {
            options.exclude_all = true;
        } else if (argument == "--covered-bed") {
            options.covered_bed_path = require_value(argc, argv, index);
            if (options.covered_bed_path.empty()) {
                throw std::runtime_error("--covered-bed requires a file path");
            }
        } else if (argument == "--mask-covered") {
            options.mask_covered = true;
        } else if (argument == "--min-unmasked-bases") {
            options.minimum_unmasked_bases = parse_uint64(
                require_value(argc, argv, index), argument);
        } else if (argument == "--query-merge-gap" ||
                   argument == "--merge-size") {
            options.query_merge_gap = parse_uint64(
                require_value(argc, argv, index), argument);
        } else if (argument == "-t" || argument == "--threads" ||
                   argument == "-c" || argument == "--cores") {
            const std::uint64_t parsed = parse_uint64(
                require_value(argc, argv, index), argument);
            if (parsed == 0 || parsed > std::numeric_limits<std::size_t>::max()) {
                throw std::runtime_error("thread count must be at least 1");
            }
            options.threads = static_cast<std::size_t>(parsed);
        } else {
            throw std::runtime_error("unknown argument: " + argument);
        }
    }
    const bool have_single =
        !options.input_path.empty() || !options.output_path.empty();
    const bool have_batch =
        !options.input_list_path.empty() || !options.output_list_path.empty();
    if (have_single && have_batch) {
        throw std::runtime_error(
            "single -i/-o mode and batch -I/-O mode are mutually exclusive");
    }
    if (!have_single && !have_batch) {
        throw std::runtime_error("provide either -i/-o or -I/-O");
    }
    if (have_single &&
        (options.input_path.empty() || options.output_path.empty())) {
        throw std::runtime_error("-i/--input and -o/--output must be used together");
    }
    if (have_batch &&
        (options.input_list_path.empty() || options.output_list_path.empty())) {
        throw std::runtime_error(
            "-I/--input-list and -O/--output-list must be used together");
    }
    if (have_batch && !options.covered_bed_path.empty()) {
        throw std::runtime_error(
            "--covered-bed is available only with single -i/-o mode");
    }
    if (!options.mask_covered && options.minimum_unmasked_bases != 0) {
        throw std::runtime_error(
            "--min-unmasked-bases requires --mask-covered");
    }
    std::error_code error;
    const std::string source_path = have_single
        ? options.input_path : options.input_list_path;
    if (!std::filesystem::is_regular_file(source_path, error)) {
        throw std::runtime_error(
            std::string(have_single ? "input FASTA" : "input list") +
            " is not a regular file: " + source_path);
    }
    if (have_batch &&
        !std::filesystem::is_regular_file(options.output_list_path, error)) {
        throw std::runtime_error(
            "output list is not a regular file: " + options.output_list_path);
    }
    if (have_single) {
        const auto input_absolute =
            std::filesystem::absolute(options.input_path, error).lexically_normal();
        error.clear();
        const auto output_absolute =
            std::filesystem::absolute(options.output_path, error).lexically_normal();
        if (input_absolute == output_absolute) {
            throw std::runtime_error("input and output FASTA paths must differ");
        }
    }
    return options;
}

class IndexedFasta {
public:
    explicit IndexedFasta(
        std::string path,
        std::string index_path = {}
    ) : path_(std::move(path)) {
        if (path_.size() >= 3 && path_.substr(path_.size() - 3) == ".gz") {
            throw std::runtime_error(
                "compressed FASTA is unsupported for .fai random access: " + path_);
        }
        const std::string fai_path =
            index_path.empty() ? path_ + ".fai" : std::move(index_path);
        std::ifstream index(fai_path);
        if (!index) {
            throw std::runtime_error(
                "FASTA index not found: " + fai_path +
                "\nPlease create it with: samtools faidx " + path_);
        }
        fd_ = ::open(path_.c_str(), O_RDONLY);
        if (fd_ < 0) {
            throw std::runtime_error(
                "cannot open FASTA " + path_ + ": " + std::strerror(errno));
        }
        try {
            read_index(index, fai_path);
            if (entries_.empty()) {
                throw std::runtime_error("FASTA index has no records: " + fai_path);
            }
            for (FastaEntry& entry : entries_) {
                entry.header = read_header_before(entry.offset);
                const std::string header_name = first_header_field(entry.header);
                if (header_name != entry.name) {
                    throw std::runtime_error(
                        "FASTA and index disagree near record " + entry.name +
                        "; recreate the stale index with: samtools faidx " + path_);
                }
            }
        } catch (...) {
            ::close(fd_);
            fd_ = -1;
            throw;
        }
    }

    IndexedFasta(const IndexedFasta&) = delete;
    IndexedFasta& operator=(const IndexedFasta&) = delete;

    ~IndexedFasta() {
        if (fd_ >= 0) ::close(fd_);
    }

    const std::string& path() const { return path_; }
    const std::vector<FastaEntry>& entries() const { return entries_; }

    std::string sequence(const FastaEntry& entry) const {
        if (entry.length > std::numeric_limits<std::size_t>::max()) {
            throw std::runtime_error(
                "FASTA record is too large for this process: " + entry.name);
        }
        std::string result(static_cast<std::size_t>(entry.length), '\0');
        if (entry.length == 0) return result;
        if (entry.line_bases == 0 || entry.line_width < entry.line_bases) {
            throw std::runtime_error(
                "invalid line widths in FASTA index for record: " + entry.name);
        }
        std::uint64_t position = 0;
        while (position < entry.length) {
            const std::uint64_t line_index = position / entry.line_bases;
            const std::uint64_t within_line = position % entry.line_bases;
            const std::uint64_t count = std::min(
                entry.length - position, entry.line_bases - within_line);
            const std::uint64_t byte_offset = entry.offset +
                line_index * entry.line_width + within_line;
            pread_exact(
                &result[static_cast<std::size_t>(position)],
                static_cast<std::size_t>(count), byte_offset,
                "sequence for " + entry.name);
            position += count;
        }
        return result;
    }

private:
    std::string path_;
    int fd_ = -1;
    std::vector<FastaEntry> entries_;

    static std::vector<std::string> split_tabs(const std::string& line) {
        std::vector<std::string> fields;
        std::size_t start = 0;
        while (true) {
            const std::size_t tab = line.find('\t', start);
            fields.push_back(line.substr(start, tab - start));
            if (tab == std::string::npos) break;
            start = tab + 1;
        }
        return fields;
    }

    static std::string first_header_field(const std::string& header) {
        std::size_t end = 0;
        while (end < header.size() &&
               !std::isspace(static_cast<unsigned char>(header[end]))) {
            ++end;
        }
        return header.substr(0, end);
    }

    void read_index(std::ifstream& input, const std::string& fai_path) {
        std::string line;
        std::size_t line_number = 0;
        while (std::getline(input, line)) {
            ++line_number;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            const std::vector<std::string> fields = split_tabs(line);
            if (fields.size() < 5 || fields[0].empty()) {
                throw std::runtime_error(
                    fai_path + ":" + std::to_string(line_number) +
                    ": expected five tab-delimited .fai fields");
            }
            const std::string context =
                fai_path + ":" + std::to_string(line_number);
            FastaEntry entry;
            entry.name = fields[0];
            entry.length = parse_uint64(fields[1], context);
            entry.offset = parse_uint64(fields[2], context);
            entry.line_bases = parse_uint64(fields[3], context);
            entry.line_width = parse_uint64(fields[4], context);
            entries_.push_back(std::move(entry));
        }
    }

    void pread_exact(
        char* destination,
        std::size_t count,
        std::uint64_t offset,
        const std::string& context
    ) const {
        std::size_t completed = 0;
        while (completed < count) {
            const auto current_offset = static_cast<off_t>(offset + completed);
            if (current_offset < 0 ||
                static_cast<std::uint64_t>(current_offset) != offset + completed) {
                throw std::runtime_error("FASTA offset is too large while reading " + context);
            }
            const ssize_t result = ::pread(
                fd_, destination + completed, count - completed, current_offset);
            if (result < 0) {
                if (errno == EINTR) continue;
                throw std::runtime_error(
                    "cannot read " + context + " from " + path_ + ": " +
                    std::strerror(errno));
            }
            if (result == 0) {
                throw std::runtime_error(
                    "unexpected end of FASTA while reading " + context +
                    "; recreate its .fai index");
            }
            completed += static_cast<std::size_t>(result);
        }
    }

    char pread_byte(std::uint64_t offset, const std::string& context) const {
        char value = '\0';
        pread_exact(&value, 1, offset, context);
        return value;
    }

    std::string read_header_before(std::uint64_t sequence_offset) const {
        if (sequence_offset == 0) {
            throw std::runtime_error(
                "invalid zero sequence offset in FASTA index: " + path_ + ".fai");
        }
        std::uint64_t position = sequence_offset;
        while (position > 0) {
            const char value = pread_byte(position - 1, "FASTA header");
            if (value != '\n' && value != '\r') break;
            --position;
        }
        const std::uint64_t end = position;
        while (position > 0) {
            const char value = pread_byte(position - 1, "FASTA header");
            if (value == '\n' || value == '\r') break;
            --position;
        }
        if (end <= position || end - position > std::numeric_limits<std::size_t>::max()) {
            throw std::runtime_error("cannot recover FASTA header before indexed sequence");
        }
        std::string raw(static_cast<std::size_t>(end - position), '\0');
        pread_exact(&raw[0], raw.size(), position, "FASTA header");
        if (raw.empty() || raw.front() != '>') {
            throw std::runtime_error(
                "indexed sequence is not preceded by a FASTA header in " + path_ +
                "; recreate its .fai index");
        }
        raw.erase(raw.begin());
        if (raw.empty()) {
            throw std::runtime_error("empty FASTA header in " + path_);
        }
        return raw;
    }
};

int base_to_int(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

bool is_lower_acgt(char base) {
    return base == 'a' || base == 'c' || base == 'g' || base == 't';
}

template <class Emit>
void emit_kmers(
    const std::string& sequence,
    std::uint64_t begin,
    std::uint64_t end,
    Emit emit
) {
    std::uint64_t forward = 0;
    std::uint64_t reverse = 0;
    std::size_t have = 0;
    std::array<std::uint8_t, KMER_SIZE> lowercase{};
    std::size_t lowercase_count = 0;
    std::size_t ring_position = 0;

    for (std::uint64_t position = begin; position < end; ++position) {
        const char base = sequence[static_cast<std::size_t>(position)];
        const int converted = base_to_int(base);
        if (converted < 0) {
            forward = reverse = 0;
            have = lowercase_count = ring_position = 0;
            lowercase.fill(0);
            continue;
        }
        const std::uint8_t is_lower =
            static_cast<std::uint8_t>(is_lower_acgt(base));
        if (have < KMER_SIZE) {
            forward = (forward << 2) | static_cast<std::uint64_t>(converted);
            reverse |= static_cast<std::uint64_t>(3 - converted) << (2 * have);
            lowercase[ring_position] = is_lower;
            lowercase_count += is_lower;
            ring_position = (ring_position + 1) % KMER_SIZE;
            ++have;
        } else {
            lowercase_count -= lowercase[ring_position];
            lowercase[ring_position] = is_lower;
            lowercase_count += is_lower;
            ring_position = (ring_position + 1) % KMER_SIZE;
            forward = ((forward & FORWARD_MASK) << 2) |
                static_cast<std::uint64_t>(converted);
            reverse = (reverse >> 2) |
                (static_cast<std::uint64_t>(3 - converted) << TOP_SHIFT);
        }
        if (have >= KMER_SIZE) {
            emit(
                std::max(forward, reverse),
                position + 1 - KMER_SIZE,
                position + 1,
                lowercase_count == KMER_SIZE);
        }
    }
}

std::int64_t interval_length_signed(const Interval& interval) {
    const std::uint64_t length = interval.end - interval.start;
    if (length > static_cast<std::uint64_t>(
            std::numeric_limits<std::int64_t>::max())) {
        throw std::runtime_error("sequence interval exceeds INT64_MAX bases");
    }
    return static_cast<std::int64_t>(length);
}

std::int64_t saturating_score_add(
    std::int64_t score,
    std::int64_t contribution
) {
    if (contribution > std::numeric_limits<std::int64_t>::max() - score) {
        return std::numeric_limits<std::int64_t>::max();
    }
    return score + contribution;
}

std::vector<std::uint8_t> bidirectional_score_breaks(
    const std::vector<Interval>& intervals,
    std::uint64_t hard_gap,
    bool hard_gap_is_inclusive
) {
    std::vector<std::uint8_t> breaks(
        intervals.size() > 1 ? intervals.size() - 1 : 0, 0);
    if (intervals.size() < 2) return breaks;

    const auto is_hard_break = [hard_gap, hard_gap_is_inclusive](
        std::uint64_t gap
    ) {
        return hard_gap_is_inclusive ? gap >= hard_gap : gap > hard_gap;
    };

    std::int64_t score = interval_length_signed(intervals.front());
    for (std::size_t index = 0; index + 1 < intervals.size(); ++index) {
        const std::uint64_t gap =
            intervals[index + 1].start - intervals[index].end;
        if (gap > static_cast<std::uint64_t>(
                std::numeric_limits<std::int64_t>::max())) {
            throw std::runtime_error("sequence gap exceeds INT64_MAX bases");
        }
        score -= static_cast<std::int64_t>(gap);
        if (score < 0 || is_hard_break(gap)) {
            breaks[index] = 1;
            score = 0;
        }
        score = saturating_score_add(
            score, interval_length_signed(intervals[index + 1]));
    }

    score = interval_length_signed(intervals.back());
    for (std::size_t right = intervals.size() - 1; right > 0; --right) {
        const std::size_t gap_index = right - 1;
        const std::uint64_t gap =
            intervals[right].start - intervals[gap_index].end;
        score -= static_cast<std::int64_t>(gap);
        if (score < 0 || is_hard_break(gap)) {
            breaks[gap_index] = 1;
            score = 0;
        }
        score = saturating_score_add(
            score, interval_length_signed(intervals[gap_index]));
    }
    return breaks;
}

std::vector<Interval> lowercase_intervals(const std::string& sequence) {
    std::vector<Interval> intervals;
    std::uint64_t position = 0;
    while (position < sequence.size()) {
        if (!is_lower_acgt(sequence[static_cast<std::size_t>(position)])) {
            ++position;
            continue;
        }
        const std::uint64_t start = position;
        do {
            ++position;
        } while (position < sequence.size() &&
                 is_lower_acgt(sequence[static_cast<std::size_t>(position)]));
        intervals.push_back({start, position});
    }
    return intervals;
}

void add_reference_sequence(
    const std::string& sequence,
    ReferenceKmers& reference,
    bool exclude_all
) {
    if (exclude_all) {
        emit_kmers(sequence, 0, sequence.size(), [&reference](
            std::uint64_t kmer,
            std::uint64_t,
            std::uint64_t,
            bool
        ) {
            reference.searcher_presence.insert(kmer);
        });
        return;
    }

    emit_kmers(sequence, 0, sequence.size(), [&reference](
        std::uint64_t kmer,
        std::uint64_t,
        std::uint64_t,
        bool
    ) {
        std::uint32_t& count = reference.counts[kmer];
        if (count != std::numeric_limits<std::uint32_t>::max()) ++count;
        if (reference.total_occurrences !=
            std::numeric_limits<std::uint64_t>::max()) {
            ++reference.total_occurrences;
        }
    });

    const std::vector<Interval> masked = lowercase_intervals(sequence);
    const std::vector<std::uint8_t> breaks = bidirectional_score_breaks(
        masked, MASKED_REGION_HARD_GAP, true);
    std::size_t group_start = 0;
    while (group_start < masked.size()) {
        std::size_t group_end = group_start;
        while (group_end + 1 < masked.size() && !breaks[group_end]) {
            ++group_end;
        }
        const std::uint64_t merged_length =
            masked[group_end].end - masked[group_start].start;
        const std::uint32_t stored_length = static_cast<std::uint32_t>(
            std::min<std::uint64_t>(
                merged_length, std::numeric_limits<std::uint32_t>::max()));
        for (std::size_t index = group_start; index <= group_end; ++index) {
            emit_kmers(
                sequence, masked[index].start, masked[index].end,
                [&reference, stored_length](
                    std::uint64_t kmer,
                    std::uint64_t,
                    std::uint64_t,
                    bool
                ) {
                    std::uint32_t& old = reference.masked_run_lengths[kmer];
                    old = std::max(old, stored_length);
                });
        }
        group_start = group_end + 1;
    }
}

void preview_reference_sequence(
    const std::string& sequence,
    ReferenceKmers& reference
) {
    emit_kmers(sequence, 0, sequence.size(), [&reference](
        std::uint64_t kmer,
        std::uint64_t,
        std::uint64_t,
        bool
    ) {
        reference.searcher_presence.preview(kmer);
        if (reference.total_occurrences !=
            std::numeric_limits<std::uint64_t>::max()) {
            ++reference.total_occurrences;
        }
    });
}

void append_merged_hit(
    std::vector<Interval>& intervals,
    std::uint64_t start,
    std::uint64_t end
) {
    if (!intervals.empty() && start <= intervals.back().end) {
        intervals.back().end = std::max(intervals.back().end, end);
    } else {
        intervals.push_back({start, end});
    }
}

std::vector<Interval> query_hit_intervals(
    const std::string& sequence,
    const ReferenceKmers& reference,
    bool exclude_all
) {
    std::vector<Interval> hits;
    emit_kmers(sequence, 0, sequence.size(), [&reference, &hits, exclude_all](
        std::uint64_t kmer,
        std::uint64_t start,
        std::uint64_t end,
        bool fully_masked
    ) {
        if (!reference.contains(kmer)) return;
        if (fully_masked && !exclude_all) {
            const auto run = reference.masked_run_lengths.find(kmer);
            if (run == reference.masked_run_lengths.end() ||
                run->second <= MIN_MASKED_REFERENCE_RUN) {
                return;
            }
        }
        append_merged_hit(hits, start, end);
    });
    return hits;
}

std::vector<Interval> polish_deletion_intervals(
    const std::vector<Interval>& hits,
    std::uint64_t query_merge_gap
) {
    if (hits.empty()) return {};
    const std::vector<std::uint8_t> breaks = bidirectional_score_breaks(
        hits, query_merge_gap, false);

    std::vector<Interval> deletions;
    std::size_t group_start = 0;
    while (group_start < hits.size()) {
        std::size_t group_end = group_start;
        std::uint64_t hit_bases =
            hits[group_start].end - hits[group_start].start;
        while (group_end + 1 < hits.size() && !breaks[group_end]) {
            ++group_end;
            hit_bases += hits[group_end].end - hits[group_end].start;
        }
        const std::uint64_t span =
            hits[group_end].end - hits[group_start].start;
        const std::uint64_t gap_bases = span - hit_bases;
        const __int128 final_score = static_cast<__int128>(hit_bases) -
            2 * static_cast<__int128>(gap_bases);
        if (final_score >= MIN_DELETE_SCORE) {
            deletions.push_back({hits[group_start].start, hits[group_end].end});
        }
        group_start = group_end + 1;
    }
    return deletions;
}

void collect_valid_fragments(
    const std::string& sequence,
    const Interval& residual,
    CleanResult& result
) {
    std::uint64_t position = residual.start;
    while (position < residual.end) {
        if (base_to_int(sequence[static_cast<std::size_t>(position)]) < 0) {
            const std::uint64_t start = position;
            do {
                ++position;
            } while (position < residual.end &&
                     base_to_int(sequence[static_cast<std::size_t>(position)]) < 0);
            result.ambiguous_bases += position - start;
            continue;
        }
        const std::uint64_t start = position;
        do {
            ++position;
        } while (position < residual.end &&
                 base_to_int(sequence[static_cast<std::size_t>(position)]) >= 0);
        if (position - start >= MIN_RESIDUAL_SIZE) {
            result.fragments.push_back({start, position});
        } else {
            result.short_bases += position - start;
        }
    }
}

CleanResult clean_query(
    std::string sequence,
    const ReferenceKmers& reference,
    bool exclude_all,
    std::uint64_t query_merge_gap
) {
    CleanResult result;
    result.sequence = std::move(sequence);
    result.covered_intervals =
        query_hit_intervals(result.sequence, reference, exclude_all);
    const std::vector<Interval> deletions = polish_deletion_intervals(
        result.covered_intervals, query_merge_gap);
    result.hit_intervals = result.covered_intervals.size();
    result.deletion_intervals = deletions.size();

    std::uint64_t cursor = 0;
    for (const Interval& deletion : deletions) {
        if (cursor < deletion.start) {
            collect_valid_fragments(
                result.sequence, {cursor, deletion.start}, result);
        }
        result.deleted_bases += deletion.end - deletion.start;
        cursor = std::max(cursor, deletion.end);
    }
    if (cursor < result.sequence.size()) {
        collect_valid_fragments(
            result.sequence,
            {cursor, static_cast<std::uint64_t>(result.sequence.size())},
            result);
    }
    return result;
}

bool is_upper_acgt(char base) {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

CleanResult mask_covered_query(
    std::string sequence,
    const ReferenceKmers& reference,
    bool exclude_all,
    std::uint64_t minimum_unmasked_bases
) {
    CleanResult result;
    result.sequence = std::move(sequence);
    result.covered_intervals =
        query_hit_intervals(result.sequence, reference, exclude_all);
    result.hit_intervals = result.covered_intervals.size();
    for (const Interval& interval : result.covered_intervals) {
        result.masked_bases += interval.end - interval.start;
        for (std::uint64_t position = interval.start;
             position < interval.end; ++position) {
            char& base = result.sequence[static_cast<std::size_t>(position)];
            if (is_upper_acgt(base)) {
                base = static_cast<char>(
                    std::tolower(static_cast<unsigned char>(base)));
            }
        }
    }
    for (const char base : result.sequence) {
        if (is_upper_acgt(base)) ++result.unmasked_bases;
    }
    if (result.unmasked_bases >= minimum_unmasked_bases) {
        result.fragments.push_back({
            0, static_cast<std::uint64_t>(result.sequence.size()),
        });
    }
    return result;
}

std::string header_suffix(const FastaEntry& entry) {
    const std::size_t position = entry.name.size();
    if (entry.header.compare(0, position, entry.name) != 0) return {};
    return entry.header.substr(position);
}

void write_wrapped_sequence(
    std::ofstream& output,
    const std::string& sequence,
    std::uint64_t start,
    std::uint64_t end
) {
    for (std::uint64_t position = start; position < end; position += FASTA_WRAP) {
        const std::uint64_t count = std::min<std::uint64_t>(FASTA_WRAP, end - position);
        output.write(
            sequence.data() + static_cast<std::size_t>(position),
            static_cast<std::streamsize>(count));
        output.put('\n');
    }
}

void write_complete_record(
    std::ofstream& output,
    const FastaEntry& entry,
    const std::string& sequence
) {
    output << '>' << entry.header << '\n';
    write_wrapped_sequence(output, sequence, 0, sequence.size());
}

void write_query_result(
    std::ofstream& output,
    const FastaEntry& entry,
    const CleanResult& result,
    Totals& totals
) {
    const bool complete =
        result.fragments.size() == 1 && result.fragments[0].start == 0 &&
        result.fragments[0].end == result.sequence.size();
    const std::string suffix = header_suffix(entry);
    for (const Interval& fragment : result.fragments) {
        if (complete) {
            output << '>' << entry.header << '\n';
        } else {
            output << '>' << entry.name << ':' << fragment.start << '-'
                   << fragment.end << suffix << '\n';
        }
        write_wrapped_sequence(
            output, result.sequence, fragment.start, fragment.end);
        ++totals.output_records;
        totals.output_bases += fragment.end - fragment.start;
    }
    ++totals.query_records;
    totals.query_bases += result.sequence.size();
    totals.deleted_bases += result.deleted_bases;
    totals.ambiguous_bases += result.ambiguous_bases;
    totals.short_bases += result.short_bases;
    totals.masked_bases += result.masked_bases;
    totals.unmasked_bases += result.unmasked_bases;
    if (result.fragments.empty() && result.sequence.size() != 0) {
        ++totals.filtered_records;
        totals.filtered_bases += result.sequence.size();
    }
}

void write_covered_bed(
    std::ofstream& output,
    const FastaEntry& entry,
    const CleanResult& result
) {
    for (const Interval& interval : result.covered_intervals) {
        output << entry.name << '\t' << interval.start << '\t'
               << interval.end << '\n';
    }
}

bool is_existing_path(const std::string& value) {
    if (value.empty()) return false;
    std::error_code error;
    return std::filesystem::exists(value, error);
}

std::string trim_ascii_space(const std::string& text) {
    std::size_t begin = 0;
    while (begin < text.size() &&
           std::isspace(static_cast<unsigned char>(text[begin]))) {
        ++begin;
    }
    std::size_t end = text.size();
    while (end > begin &&
           std::isspace(static_cast<unsigned char>(text[end - 1]))) {
        --end;
    }
    return text.substr(begin, end - begin);
}

std::vector<std::string> split_reference_names(const std::string& text) {
    std::vector<std::string> names;
    std::size_t begin = 0;
    while (true) {
        const std::size_t comma = text.find(',', begin);
        const std::string name = trim_ascii_space(
            text.substr(begin, comma == std::string::npos
                ? std::string::npos : comma - begin));
        if (name.empty()) {
            throw std::runtime_error(
                "-r/--reference contains an empty record name: " + text);
        }
        names.push_back(name);
        if (comma == std::string::npos) break;
        begin = comma + 1;
    }
    return names;
}

std::vector<std::string> split_whitespace(const std::string& line) {
    std::istringstream input(line);
    std::vector<std::string> fields;
    std::string field;
    while (input >> field) fields.push_back(field);
    return fields;
}

std::string resolve_list_path(
    const std::string& list_path,
    const std::string& value
) {
    std::error_code error;
    std::filesystem::path path(value);
    if (path.is_relative()) {
        const std::filesystem::path absolute_list =
            std::filesystem::absolute(list_path, error);
        if (error) {
            throw std::runtime_error(
                "cannot resolve list path " + list_path + ": " + error.message());
        }
        path = absolute_list.parent_path() / path;
    }
    return path.lexically_normal().string();
}

std::string absolute_normalized_path(const std::string& value) {
    std::error_code error;
    const std::filesystem::path path = std::filesystem::absolute(value, error);
    if (error) {
        throw std::runtime_error(
            "cannot resolve path " + value + ": " + error.message());
    }
    return path.lexically_normal().string();
}

void require_regular_file(
    const std::string& path,
    const std::string& context
) {
    std::error_code error;
    if (!std::filesystem::is_regular_file(path, error)) {
        throw std::runtime_error(context + ": file not found: " + path);
    }
}

std::vector<SampleJob> read_batch_jobs(
    const std::string& input_list_path,
    const std::string& output_list_path
) {
    std::ifstream input_list(input_list_path);
    if (!input_list) {
        throw std::runtime_error("cannot open input list: " + input_list_path);
    }
    std::vector<SampleJob> jobs;
    std::unordered_set<std::string> input_samples;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input_list, line)) {
        ++line_number;
        const std::string trimmed = trim_ascii_space(line);
        if (trimmed.empty() || trimmed.front() == '#') continue;
        const std::vector<std::string> fields = split_whitespace(trimmed);
        if (fields.size() < 2 || fields.size() > 3) {
            throw std::runtime_error(
                input_list_path + ":" + std::to_string(line_number) +
                ": expected SAMPLE FASTA [FAI]");
        }
        if (!input_samples.insert(fields[0]).second) {
            throw std::runtime_error(
                input_list_path + ":" + std::to_string(line_number) +
                ": duplicate sample: " + fields[0]);
        }
        SampleJob job;
        job.sample = fields[0];
        job.input_path = resolve_list_path(input_list_path, fields[1]);
        job.fai_path = fields.size() == 3
            ? resolve_list_path(input_list_path, fields[2])
            : job.input_path + ".fai";
        require_regular_file(
            job.input_path,
            input_list_path + ":" + std::to_string(line_number) +
                ": input FASTA");
        require_regular_file(
            job.fai_path,
            input_list_path + ":" + std::to_string(line_number) +
                ": FASTA index");
        jobs.push_back(std::move(job));
    }
    if (jobs.empty()) {
        throw std::runtime_error(input_list_path + ": no sample FASTA sources");
    }

    std::ifstream output_list(output_list_path);
    if (!output_list) {
        throw std::runtime_error("cannot open output list: " + output_list_path);
    }
    std::unordered_map<std::string, std::string> outputs;
    std::unordered_set<std::string> output_paths;
    line_number = 0;
    while (std::getline(output_list, line)) {
        ++line_number;
        const std::string trimmed = trim_ascii_space(line);
        if (trimmed.empty() || trimmed.front() == '#') continue;
        const std::vector<std::string> fields = split_whitespace(trimmed);
        if (fields.size() != 2) {
            throw std::runtime_error(
                output_list_path + ":" + std::to_string(line_number) +
                ": expected SAMPLE OUTPUT_FASTA");
        }
        if (outputs.find(fields[0]) != outputs.end()) {
            throw std::runtime_error(
                output_list_path + ":" + std::to_string(line_number) +
                ": duplicate sample: " + fields[0]);
        }
        const std::string output =
            resolve_list_path(output_list_path, fields[1]);
        if (!output_paths.insert(output).second) {
            throw std::runtime_error(
                output_list_path + ":" + std::to_string(line_number) +
                ": duplicate output path: " + output);
        }
        outputs.emplace(fields[0], output);
    }
    if (outputs.size() != jobs.size()) {
        throw std::runtime_error(
            "input and output lists contain different sample sets");
    }
    for (SampleJob& job : jobs) {
        const auto found = outputs.find(job.sample);
        if (found == outputs.end()) {
            throw std::runtime_error(
                "sample is missing from output list: " + job.sample);
        }
        job.output_path = found->second;
        if (std::filesystem::path(job.input_path).lexically_normal() ==
            std::filesystem::path(job.output_path).lexically_normal()) {
            throw std::runtime_error(
                "input and output FASTA paths must differ for sample " +
                job.sample);
        }
        outputs.erase(found);
    }
    if (!outputs.empty()) {
        throw std::runtime_error(
            "sample is missing from input list: " + outputs.begin()->first);
    }
    return jobs;
}

std::string temporary_output_path(const std::string& output_path) {
    return output_path + ".tmp." + std::to_string(static_cast<long long>(::getpid()));
}

ReferenceSummary build_reference_from_fasta(
    const IndexedFasta& fasta,
    ReferenceKmers& reference,
    bool exclude_all
) {
    if (exclude_all) {
        std::uint64_t maximum_kmers = 0;
        for (const FastaEntry& entry : fasta.entries()) {
            const std::uint64_t count = entry.length >= KMER_SIZE
                ? entry.length - KMER_SIZE + 1 : 0;
            if (count > std::numeric_limits<std::uint64_t>::max() -
                    maximum_kmers) {
                throw std::runtime_error(
                    "reference 31-mer upper bound exceeds UINT64_MAX");
            }
            maximum_kmers += count;
        }
        reference.use_searcher_presence = true;
        reference.searcher_presence.initialize(maximum_kmers);
    }
    ReferenceSummary summary;
    for (const FastaEntry& entry : fasta.entries()) {
        std::string sequence = fasta.sequence(entry);
        summary.bases += sequence.size();
        ++summary.records;
        if (exclude_all) {
            preview_reference_sequence(sequence, reference);
        } else {
            add_reference_sequence(sequence, reference, false);
        }
    }
    if (exclude_all) {
        reference.searcher_presence.allocate_suffixes();
        std::cerr
            << "[shared] Allocated exact KmerSearcher-style --exclude-all index "
            << "with "
            << static_cast<unsigned>(reference.searcher_presence.bucket_bits())
            << " bucket bits and "
            << reference.searcher_presence.slot_count() << " suffix slots ("
            << static_cast<double>(reference.searcher_presence.allocated_gib())
            << " GiB).\n";
        for (const FastaEntry& entry : fasta.entries()) {
            add_reference_sequence(fasta.sequence(entry), reference, true);
        }
        reference.searcher_presence.finalize();
        if (reference.searcher_presence.overflow_size() != 0) {
            std::cerr << "[shared] Preserved "
                      << reference.searcher_presence.overflow_size()
                      << " k-mer(s) from buckets exceeding 255 distinct "
                         "suffixes.\n";
        }
    }
    return summary;
}

void log_reference_summary(
    const std::string& prefix,
    const ReferenceSummary& summary,
    const ReferenceKmers& reference,
    bool exclude_all
) {
    std::cerr
        << prefix << "Indexed " << summary.records << " reference record(s), "
        << summary.bases << " bp, " << reference.total_occurrences
        << " 31-mer occurrences, " << reference.distinct_size()
        << " distinct canonical 31-mers";
    if (exclude_all) {
        std::cerr << "; masked-repeat eligibility disabled by --exclude-all.\n";
    } else {
        std::cerr << ", and " << reference.masked_run_lengths.size()
                  << " fully masked k-mers.\n";
    }
}

void process_sample(
    const SampleJob& job,
    const Options& options,
    const ReferenceKmers* shared_reference,
    bool external_reference
) {
    IndexedFasta input(job.input_path, job.fai_path);
    std::vector<std::size_t> input_reference_indices;
    std::vector<std::uint8_t> is_input_reference(input.entries().size(), 0);

    ReferenceKmers local_reference;
    const ReferenceKmers* reference = shared_reference;
    if (!external_reference) {
        if (options.reference_given) {
            for (const std::string& name :
                 split_reference_names(options.reference)) {
                const auto& entries = input.entries();
                const auto found = std::find_if(
                    entries.begin(), entries.end(),
                    [&name](const FastaEntry& entry) {
                        return entry.name == name;
                    });
                if (found == entries.end()) {
                    throw std::runtime_error(
                        "[" + job.sample + "] reference name not found in "
                        "input FASTA: " + name);
                }
                const std::size_t index =
                    static_cast<std::size_t>(found - entries.begin());
                is_input_reference[index] = 1;
            }
        } else {
            is_input_reference[0] = 1;
        }
        for (std::size_t index = 0; index < is_input_reference.size(); ++index) {
            if (is_input_reference[index]) input_reference_indices.push_back(index);
        }
        ReferenceSummary summary;
        if (options.exclude_all) {
            std::uint64_t maximum_kmers = 0;
            for (const std::size_t index : input_reference_indices) {
                const std::uint64_t length = input.entries()[index].length;
                const std::uint64_t count = length >= KMER_SIZE
                    ? length - KMER_SIZE + 1 : 0;
                if (count > std::numeric_limits<std::uint64_t>::max() -
                        maximum_kmers) {
                    throw std::runtime_error(
                        "reference 31-mer upper bound exceeds UINT64_MAX");
                }
                maximum_kmers += count;
            }
            local_reference.use_searcher_presence = true;
            local_reference.searcher_presence.initialize(maximum_kmers);
        }
        for (const std::size_t index : input_reference_indices) {
            std::string sequence = input.sequence(input.entries()[index]);
            summary.bases += sequence.size();
            ++summary.records;
            if (options.exclude_all) {
                preview_reference_sequence(sequence, local_reference);
            } else {
                add_reference_sequence(sequence, local_reference, false);
            }
        }
        if (options.exclude_all) {
            local_reference.searcher_presence.allocate_suffixes();
            for (const std::size_t index : input_reference_indices) {
                add_reference_sequence(
                    input.sequence(input.entries()[index]), local_reference,
                    true);
            }
            local_reference.searcher_presence.finalize();
        }
        reference = &local_reference;
        log_reference_summary(
            "[" + job.sample + "] ", summary, local_reference,
            options.exclude_all);
    }
    if (reference == nullptr) {
        throw std::runtime_error("internal error: missing reference k-mer index");
    }

    std::vector<const FastaEntry*> queries;
    queries.reserve(input.entries().size());
    for (std::size_t index = 0; index < input.entries().size(); ++index) {
        if (!external_reference && is_input_reference[index]) continue;
        queries.push_back(&input.entries()[index]);
    }

    const std::string temporary_path = temporary_output_path(job.output_path);
    const std::string temporary_bed_path = options.covered_bed_path.empty()
        ? std::string()
        : temporary_output_path(options.covered_bed_path);
    std::remove(temporary_path.c_str());
    if (!temporary_bed_path.empty()) {
        std::remove(temporary_bed_path.c_str());
    }
    try {
        std::ofstream output(temporary_path, std::ios::binary | std::ios::trunc);
        output.exceptions(std::ios::badbit | std::ios::failbit);
        std::unique_ptr<std::ofstream> covered_bed;
        if (!temporary_bed_path.empty()) {
            covered_bed = std::make_unique<std::ofstream>(
                temporary_bed_path, std::ios::binary | std::ios::trunc);
            covered_bed->exceptions(std::ios::badbit | std::ios::failbit);
        }

        std::uint64_t copied_reference_records = 0;
        std::uint64_t copied_reference_bases = 0;
        if (!external_reference) {
            for (const std::size_t index : input_reference_indices) {
                const FastaEntry& entry = input.entries()[index];
                const std::string sequence = input.sequence(entry);
                write_complete_record(output, entry, sequence);
                ++copied_reference_records;
                copied_reference_bases += sequence.size();
            }
        }

        Totals totals;
        std::atomic<std::size_t> next_task{0};
        std::atomic<bool> stop{false};
        std::mutex write_mutex;
        std::condition_variable write_condition;
        std::size_t next_write = 0;
        std::exception_ptr worker_error;

        auto record_error = [&](std::exception_ptr error) {
            {
                std::lock_guard<std::mutex> lock(write_mutex);
                if (!worker_error) worker_error = std::move(error);
                stop.store(true, std::memory_order_release);
            }
            write_condition.notify_all();
        };

        auto worker = [&]() {
            while (!stop.load(std::memory_order_acquire)) {
                const std::size_t query_index =
                    next_task.fetch_add(1, std::memory_order_relaxed);
                if (query_index >= queries.size()) return;
                CleanResult result;
                try {
                    if (options.mask_covered) {
                        result = mask_covered_query(
                            input.sequence(*queries[query_index]), *reference,
                            options.exclude_all,
                            options.minimum_unmasked_bases);
                    } else {
                        result = clean_query(
                            input.sequence(*queries[query_index]), *reference,
                            options.exclude_all, options.query_merge_gap);
                    }
                } catch (...) {
                    record_error(std::current_exception());
                    return;
                }

                std::unique_lock<std::mutex> lock(write_mutex);
                write_condition.wait(lock, [&]() {
                    return stop.load(std::memory_order_acquire) ||
                        query_index == next_write;
                });
                if (stop.load(std::memory_order_acquire)) return;
                try {
                    if (covered_bed) {
                        write_covered_bed(
                            *covered_bed, *queries[query_index], result);
                    }
                    write_query_result(
                        output, *queries[query_index], result, totals);
                    ++next_write;
                } catch (...) {
                    if (!worker_error) worker_error = std::current_exception();
                    stop.store(true, std::memory_order_release);
                    lock.unlock();
                    write_condition.notify_all();
                    return;
                }
                lock.unlock();
                write_condition.notify_all();
            }
        };

        const std::size_t thread_count = std::min(
            options.threads, std::max<std::size_t>(1, queries.size()));
        std::vector<std::thread> threads;
        threads.reserve(thread_count);
        for (std::size_t index = 0; index < thread_count; ++index) {
            threads.emplace_back(worker);
        }
        for (std::thread& thread : threads) thread.join();
        if (worker_error) std::rethrow_exception(worker_error);

        output.close();
        if (covered_bed) covered_bed->close();
        if (!temporary_bed_path.empty() &&
            std::rename(
                temporary_bed_path.c_str(),
                options.covered_bed_path.c_str()) != 0) {
            throw std::runtime_error(
                "cannot move completed BED to " + options.covered_bed_path +
                ": " + std::strerror(errno));
        }
        if (std::rename(temporary_path.c_str(), job.output_path.c_str()) != 0) {
            if (!options.covered_bed_path.empty()) {
                std::remove(options.covered_bed_path.c_str());
            }
            throw std::runtime_error(
                "cannot move completed output to " + job.output_path + ": " +
                std::strerror(errno));
        }
        if (options.mask_covered) {
            std::cerr
                << "[" << job.sample << "] Masked " << totals.masked_bases
                << " reference-covered bp across " << totals.query_records
                << " query record(s); wrote "
                << copied_reference_records + totals.output_records
                << " record(s), "
                << copied_reference_bases + totals.output_bases << " bp to "
                << job.output_path << ". Filtered " << totals.filtered_records
                << " record(s), " << totals.filtered_bases
                << " bp, below " << options.minimum_unmasked_bases
                << " unmasked uppercase A/C/G/T bases.\n";
        } else {
            std::cerr
                << "[" << job.sample << "] Cleaned " << totals.query_records
                << " query record(s), " << totals.query_bases << " bp; wrote "
                << copied_reference_records + totals.output_records
                << " record(s), "
                << copied_reference_bases + totals.output_bases << " bp to "
                << job.output_path << ". Removed " << totals.deleted_bases
                << " redundant bp, " << totals.ambiguous_bases
                << " ambiguous bp, and " << totals.short_bases
                << " bp in residual fragments shorter than 100.\n";
        }
        if (!options.covered_bed_path.empty()) {
            std::cerr << "[" << job.sample
                      << "] Wrote raw reference-covered 31-mer intervals to "
                      << options.covered_bed_path << ".\n";
        }
    } catch (...) {
        std::remove(temporary_path.c_str());
        if (!temporary_bed_path.empty()) {
            std::remove(temporary_bed_path.c_str());
        }
        throw;
    }
}

int run(const Options& options) {
    std::vector<SampleJob> jobs;
    if (!options.input_path.empty()) {
        jobs.push_back({
            "single", options.input_path, options.input_path + ".fai",
            options.output_path,
        });
    } else {
        jobs = read_batch_jobs(
            options.input_list_path, options.output_list_path);
    }

    const bool external_reference =
        options.reference_given && is_existing_path(options.reference);
    std::unordered_set<std::string> protected_paths;
    for (const SampleJob& job : jobs) {
        protected_paths.insert(absolute_normalized_path(job.input_path));
        protected_paths.insert(absolute_normalized_path(job.fai_path));
    }
    if (!options.input_list_path.empty()) {
        protected_paths.insert(
            absolute_normalized_path(options.input_list_path));
        protected_paths.insert(
            absolute_normalized_path(options.output_list_path));
    }
    if (external_reference) {
        protected_paths.insert(absolute_normalized_path(options.reference));
        protected_paths.insert(
            absolute_normalized_path(options.reference + ".fai"));
    }
    std::unordered_set<std::string> output_paths;
    for (const SampleJob& job : jobs) {
        const std::string output = absolute_normalized_path(job.output_path);
        if (protected_paths.find(output) != protected_paths.end()) {
            throw std::runtime_error(
                "output path would overwrite an input, index, reference, or "
                "list file for sample " + job.sample + ": " + output);
        }
        if (!output_paths.insert(output).second) {
            throw std::runtime_error("duplicate output path: " + output);
        }
    }
    if (!options.covered_bed_path.empty()) {
        const std::string bed =
            absolute_normalized_path(options.covered_bed_path);
        if (protected_paths.find(bed) != protected_paths.end()) {
            throw std::runtime_error(
                "--covered-bed would overwrite an input, index, reference, "
                "or list file: " + bed);
        }
        if (!output_paths.insert(bed).second) {
            throw std::runtime_error(
                "--covered-bed must differ from every FASTA output: " + bed);
        }
    }

    ReferenceKmers shared_reference;
    if (external_reference) {
        std::error_code error;
        if (!std::filesystem::is_regular_file(options.reference, error)) {
            throw std::runtime_error(
                "external reference is not a regular file: " + options.reference);
        }
        IndexedFasta external(options.reference);
        const ReferenceSummary summary = build_reference_from_fasta(
            external, shared_reference, options.exclude_all);
        log_reference_summary(
            "[shared] ", summary, shared_reference, options.exclude_all);
    }

    for (const SampleJob& job : jobs) {
        process_sample(
            job, options,
            external_reference ? &shared_reference : nullptr,
            external_reference);
    }
    if (jobs.size() > 1) {
        std::cerr << "Completed " << jobs.size() << " sample/output jobs.\n";
    }
    return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(parse_options(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
