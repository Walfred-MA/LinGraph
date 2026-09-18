#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr std::size_t KMER_SIZE = 31;
constexpr std::size_t TOP_SHIFT = 2 * (KMER_SIZE - 1);
constexpr std::uint64_t FORWARD_MASK =
    (std::uint64_t{1} << (2 * (KMER_SIZE - 1))) - 1;
// A Gram-matrix feature is multiplied on both sides.  Using sqrt(0.1) as
// the feature amplitude therefore gives a lowercase-overlapping occurrence
// exactly 0.1 weight in diagonal variance.
using Score = float;
static_assert(sizeof(Score) == 4, "KmerMatch requires 32-bit float scores");
static_assert(
    std::numeric_limits<Score>::is_iec559 &&
        std::numeric_limits<Score>::digits == 24,
    "KmerMatch requires IEEE-754 binary32 scores");

constexpr Score MASKED_AMPLITUDE = 0.31622776601683793320F;

struct FastaRecord {
    std::string name;
    std::string sequence;
};

struct Interval {
    std::uint64_t start = 0;
    std::uint64_t end = 0;
};

struct ReferenceView {
    const FastaRecord* record = nullptr;
    std::vector<Interval> intervals;
};

struct KmerOccurrence {
    std::uint8_t occurrence = 0;
    std::uint8_t masked_occurrence = 0;
};

using KmerCounts = std::unordered_map<std::uint64_t, KmerOccurrence>;

struct KmerIndexEntry {
    std::uint32_t index = 0;
    // Largest count observed in one indexed reference record.  Counts are
    // saturated because the dense per-record matrix also stores uint8_t.
    std::uint8_t occurrence = 0;
};

using KmerIndex = std::unordered_map<std::uint64_t, KmerIndexEntry>;

struct Options {
    std::string input_path;
    std::string reference_path;
    std::string output_path;
    std::string bed_path;
    bool masking = true;
    bool repeat_score = false;
    bool paired = false;
};

void usage(const char* program) {
    std::cerr
        << "Usage: " << program
        << " -i input.fa -r ref.fa -o report.txt [-m 0|1] [-b regions.bed] [--repeat-score] [--paired]\n\n"
        << "Compare every input FASTA record (A) with every selected reference\n"
        << "FASTA record (B) using occurrence-weighted canonical 31-mers. BED\n"
        << "coordinates are\n"
        << "0-based half-open. With -m 1 (default), 31-mers overlapping any\n"
        << "lowercase base receive 0.1 weight in variance.\n"
        << "--repeat-score uses reference occurrence weight 1/n instead of 1/n^2.\n"
        << "--paired compares input record i only with reference record i.\n";
}

std::string require_value(int argc, char** argv, int& index) {
    if (index + 1 >= argc) {
        throw std::runtime_error(std::string("missing value after ") + argv[index]);
    }
    return argv[++index];
}

bool parse_bool(const std::string& text) {
    if (text == "1" || text == "true" || text == "TRUE") return true;
    if (text == "0" || text == "false" || text == "FALSE") return false;
    throw std::runtime_error("masking must be 0/1 or false/true, found: " + text);
}

bool ends_with_case_insensitive(
    const std::string& text,
    const std::string& suffix
) {
    if (text.size() < suffix.size()) return false;
    return std::equal(
        suffix.rbegin(), suffix.rend(), text.rbegin(),
        [](char left, char right) {
            return std::tolower(static_cast<unsigned char>(left)) ==
                std::tolower(static_cast<unsigned char>(right));
        });
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
        } else if (argument == "-r" || argument == "--reference") {
            options.reference_path = require_value(argc, argv, index);
        } else if (argument == "-o" || argument == "--output") {
            options.output_path = require_value(argc, argv, index);
        } else if (argument == "-m" || argument == "--masking") {
            options.masking = parse_bool(require_value(argc, argv, index));
        } else if (argument == "--repeat-score") {
            options.repeat_score = true;
        } else if (argument == "--paired") {
            options.paired = true;
        } else if (argument == "-b" || argument == "--bed") {
            options.bed_path = require_value(argc, argv, index);
        } else {
            throw std::runtime_error("unknown argument: " + argument);
        }
    }
    if (options.input_path.empty() || options.reference_path.empty() ||
        options.output_path.empty()) {
        throw std::runtime_error("-i, -r, and -o are required");
    }
    if (ends_with_case_insensitive(options.input_path, ".bed")) {
        throw std::runtime_error("-i/--input must be FASTA, not BED");
    }
    return options;
}

std::string first_token(const std::string& text) {
    std::istringstream stream(text);
    std::string token;
    stream >> token;
    return token;
}

std::vector<FastaRecord> read_fasta(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open FASTA: " + path);

    std::vector<FastaRecord> records;
    std::unordered_set<std::string> names;
    std::string line;
    FastaRecord current;
    bool have_record = false;

    auto finish = [&]() {
        if (!have_record) return;
        if (!names.insert(current.name).second) {
            throw std::runtime_error(
                path + ": duplicate FASTA record name: " + current.name);
        }
        records.push_back(std::move(current));
        current = FastaRecord{};
        have_record = false;
    };

    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line.front() == '>') {
            finish();
            current.name = first_token(line.substr(1));
            if (current.name.empty()) {
                throw std::runtime_error(
                    path + ":" + std::to_string(line_number) +
                    ": empty FASTA record name");
            }
            have_record = true;
            continue;
        }
        if (!have_record) {
            throw std::runtime_error(
                path + ":" + std::to_string(line_number) +
                ": sequence before first FASTA header");
        }
        for (char base : line) {
            if (!std::isspace(static_cast<unsigned char>(base))) {
                current.sequence.push_back(base);
            }
        }
    }
    finish();
    if (records.empty()) throw std::runtime_error("FASTA has no records: " + path);
    return records;
}

std::uint64_t parse_coordinate(
    const std::string& text,
    const std::string& context
) {
    if (text.empty() || text.front() == '-') {
        throw std::runtime_error(context + ": invalid coordinate: " + text);
    }
    std::size_t consumed = 0;
    std::uint64_t value = 0;
    try {
        value = std::stoull(text, &consumed);
    } catch (...) {
        throw std::runtime_error(context + ": invalid coordinate: " + text);
    }
    if (consumed != text.size()) {
        throw std::runtime_error(context + ": invalid coordinate: " + text);
    }
    return value;
}

std::unordered_map<std::string, std::vector<Interval>> read_bed(
    const std::string& path
) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open BED: " + path);
    std::unordered_map<std::string, std::vector<Interval>> intervals;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line.front() == '#' ||
            line.rfind("track", 0) == 0 || line.rfind("browser", 0) == 0) {
            continue;
        }
        std::istringstream stream(line);
        std::string name;
        std::string start_text;
        std::string end_text;
        if (!(stream >> name >> start_text >> end_text)) {
            throw std::runtime_error(
                path + ":" + std::to_string(line_number) +
                ": expected at least three BED columns");
        }
        const std::string context = path + ":" + std::to_string(line_number);
        const std::uint64_t start = parse_coordinate(start_text, context);
        const std::uint64_t end = parse_coordinate(end_text, context);
        if (end < start) {
            throw std::runtime_error(context + ": BED end is smaller than start");
        }
        intervals[name].push_back({start, end});
    }
    if (intervals.empty()) throw std::runtime_error("BED has no intervals: " + path);
    return intervals;
}

std::vector<ReferenceView> select_references(
    const std::vector<FastaRecord>& references,
    const std::string& bed_path
) {
    std::vector<ReferenceView> selected;
    if (bed_path.empty()) {
        selected.reserve(references.size());
        for (const auto& reference : references) {
            selected.push_back({
                &reference,
                {{0, static_cast<std::uint64_t>(reference.sequence.size())}},
            });
        }
        return selected;
    }

    auto bed = read_bed(bed_path);
    for (const auto& reference : references) {
        auto found = bed.find(reference.name);
        if (found == bed.end()) continue;
        for (const Interval& interval : found->second) {
            if (interval.end > reference.sequence.size()) {
                throw std::runtime_error(
                    bed_path + ": interval for " + reference.name +
                    " extends past FASTA length " +
                    std::to_string(reference.sequence.size()));
            }
        }
        selected.push_back({&reference, std::move(found->second)});
        bed.erase(found);
    }
    if (!bed.empty()) {
        throw std::runtime_error(
            bed_path + ": BED contig is absent from reference FASTA: " +
            bed.begin()->first);
    }
    if (selected.empty()) {
        throw std::runtime_error("BED selected no reference FASTA records");
    }
    return selected;
}

int base_to_int(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

template <class Emit>
void emit_interval_kmers(
    const std::string& sequence,
    const Interval& interval,
    bool downweight_lowercase,
    Emit emit
) {
    std::uint64_t forward = 0;
    std::uint64_t reverse = 0;
    std::size_t have = 0;
    std::array<std::uint8_t, KMER_SIZE> lowercase{};
    std::size_t lowercase_count = 0;
    std::size_t ring_position = 0;

    for (std::uint64_t position = interval.start; position < interval.end; ++position) {
        const char base = sequence[static_cast<std::size_t>(position)];
        const int converted = base_to_int(base);
        if (converted < 0) {
            forward = reverse = 0;
            have = lowercase_count = ring_position = 0;
            lowercase.fill(0);
            continue;
        }
        const std::uint8_t is_lower = static_cast<std::uint8_t>(
            std::islower(static_cast<unsigned char>(base)) != 0);
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
                downweight_lowercase && lowercase_count != 0);
        }
    }
}

template <class Emit>
void emit_reference_kmers(
    const ReferenceView& reference,
    bool masking,
    Emit emit
) {
    for (const Interval& interval : reference.intervals) {
        emit_interval_kmers(reference.record->sequence, interval, masking, emit);
    }
}

void increment_count(KmerOccurrence& count, bool masked) {
    if (count.occurrence == std::numeric_limits<std::uint8_t>::max()) return;
    ++count.occurrence;
    if (masked) ++count.masked_occurrence;
}

KmerCounts reference_kmer_counts(
    const ReferenceView& reference,
    bool masking
) {
    KmerCounts kmers;
    emit_reference_kmers(reference, masking, [&](
        std::uint64_t kmer, bool masked
    ) {
        increment_count(kmers[kmer], masked);
    });
    return kmers;
}

KmerCounts query_kmer_counts(const FastaRecord& query, bool masking) {
    KmerCounts kmers;
    emit_interval_kmers(
        query.sequence,
        {0, static_cast<std::uint64_t>(query.sequence.size())},
        masking,
        [&](std::uint64_t kmer, bool masked) {
            increment_count(kmers[kmer], masked);
        });
    return kmers;
}

Score effective_count(const KmerOccurrence& count) {
    return static_cast<Score>(
        count.occurrence - count.masked_occurrence) +
        MASKED_AMPLITUDE * count.masked_occurrence;
}

Score standalone_variance(const KmerOccurrence& count, bool repeat_score) {
    if (count.occurrence == 0) return 0.0F;
    if (repeat_score) {
        const Score value = effective_count(count);
        return value * value / static_cast<Score>(count.occurrence);
    }
    const Score value = effective_count(count) /
        static_cast<Score>(count.occurrence);
    return value * value;
}

Score occurrence_weight(std::uint8_t occurrence, bool repeat_score) {
    const Score count = static_cast<Score>(occurrence);
    return repeat_score ? 1.0F / count : 1.0F / (count * count);
}

struct RawKmer {
    std::uint64_t kmer = 0;
    bool masked = false;
};

using CondensedKmers =
    std::vector<std::pair<std::uint64_t, KmerOccurrence>>;

template <class Generate>
CondensedKmers collect_condensed_kmers(
    std::size_t maximum_kmer_count,
    Generate generate
) {
    std::vector<RawKmer> raw;
    raw.reserve(maximum_kmer_count);
    generate([&](std::uint64_t kmer, bool masked) {
        raw.push_back({kmer, masked});
    });
    std::sort(raw.begin(), raw.end(), [](const RawKmer& left,
                                         const RawKmer& right) {
        return left.kmer < right.kmer;
    });
    CondensedKmers condensed;
    condensed.reserve(raw.size());
    std::size_t begin = 0;
    while (begin < raw.size()) {
        std::size_t end = begin;
        KmerOccurrence count;
        while (end < raw.size() && raw[end].kmer == raw[begin].kmer) {
            increment_count(count, raw[end].masked);
            ++end;
        }
        condensed.emplace_back(raw[begin].kmer, count);
        begin = end;
    }
    return condensed;
}

CondensedKmers condensed_query_kmers(
    const FastaRecord& query, bool masking
) {
    const std::size_t maximum_kmer_count =
        query.sequence.size() >= KMER_SIZE
            ? query.sequence.size() - KMER_SIZE + 1
            : 0;
    return collect_condensed_kmers(maximum_kmer_count, [&](auto emit) {
        emit_interval_kmers(
            query.sequence,
            {0, static_cast<std::uint64_t>(query.sequence.size())},
            masking, emit);
    });
}

CondensedKmers condensed_reference_kmers(
    const ReferenceView& reference, bool masking
) {
    std::size_t maximum_kmer_count = 0;
    for (const Interval& interval : reference.intervals) {
        const std::uint64_t length = interval.end - interval.start;
        if (length >= KMER_SIZE) {
            maximum_kmer_count += static_cast<std::size_t>(
                length - KMER_SIZE + 1);
        }
    }
    return collect_condensed_kmers(maximum_kmer_count, [&](auto emit) {
        emit_reference_kmers(reference, masking, emit);
    });
}

void compare_pairwise_sorted(
    const FastaRecord& query,
    const ReferenceView& reference,
    bool masking,
    bool repeat_score,
    std::ofstream& output
) {
    const CondensedKmers query_kmers =
        condensed_query_kmers(query, masking);
    const CondensedKmers reference_kmers =
        condensed_reference_kmers(reference, masking);
    Score common = 0.0F;
    Score only_query = 0.0F;
    Score only_reference = 0.0F;
    Score query_variance = 0.0F;
    Score reference_variance = 0.0F;
    std::size_t query_index = 0;
    std::size_t reference_index = 0;
    while (query_index < query_kmers.size() ||
           reference_index < reference_kmers.size()) {
        if (reference_index == reference_kmers.size() ||
            (query_index < query_kmers.size() &&
             query_kmers[query_index].first <
                 reference_kmers[reference_index].first)) {
            const Score variance = standalone_variance(
                query_kmers[query_index].second, repeat_score);
            only_query += variance;
            query_variance += variance;
            ++query_index;
            continue;
        }
        if (query_index == query_kmers.size() ||
            reference_kmers[reference_index].first <
                query_kmers[query_index].first) {
            const Score variance = standalone_variance(
                reference_kmers[reference_index].second, repeat_score);
            only_reference += variance;
            reference_variance += variance;
            ++reference_index;
            continue;
        }
        const KmerOccurrence& query_count =
            query_kmers[query_index].second;
        const KmerOccurrence& reference_count =
            reference_kmers[reference_index].second;
        const Score weight = occurrence_weight(
            reference_count.occurrence, repeat_score);
        const Score query_value = effective_count(query_count);
        const Score reference_value = effective_count(reference_count);
        common += query_value * reference_value * weight;
        query_variance += query_value * query_value * weight;
        reference_variance += reference_value * reference_value * weight;
        ++query_index;
        ++reference_index;
    }
    output
        << query.name << '\t' << reference.record->name << '\t'
        << common << '\t' << only_query << '\t'
        << only_reference << '\t' << query_variance << '\t'
        << reference_variance << '\n';
}

void write_header(std::ofstream& output) {
    output
        << "seqAname\tseqBname\tcommon_kmers\t"
        << "kmers_only_in_A\tkmers_only_in_B\t"
        << "variance_A\tvariance_B\n";
}

void compare_single(
    const FastaRecord& query,
    const ReferenceView& reference,
    bool masking,
    bool repeat_score,
    std::ofstream& output
) {
    const KmerCounts reference_counts =
        reference_kmer_counts(reference, masking);
    const KmerCounts query_counts = query_kmer_counts(query, masking);
    Score common = 0.0F;
    Score only_query = 0.0F;
    Score only_reference = 0.0F;
    Score query_variance = 0.0F;
    Score reference_variance = 0.0F;
    for (const auto& [kmer, query_count] : query_counts) {
        const auto found = reference_counts.find(kmer);
        if (found == reference_counts.end()) {
            // There is no reference-derived denominator for this feature.
            // Use its own count as denominator; an entirely uppercase-only
            // feature contributes 1 (or its count with --repeat-score);
            // an entirely masked feature receives a factor of 0.1.
            const Score variance = standalone_variance(query_count, repeat_score);
            only_query += variance;
            query_variance += variance;
            continue;
        }
        const Score weight = occurrence_weight(found->second.occurrence, repeat_score);
        const Score query_value = effective_count(query_count);
        const Score reference_value = effective_count(found->second);
        common += query_value * reference_value * weight;
        query_variance += query_value * query_value * weight;
    }
    for (const auto& [kmer, reference_count] : reference_counts) {
        const Score reference_value = effective_count(reference_count);
        reference_variance += reference_value * reference_value *
            occurrence_weight(reference_count.occurrence, repeat_score);
        if (query_counts.find(kmer) == query_counts.end()) {
            only_reference += standalone_variance(reference_count, repeat_score);
        }
    }
    output
        << query.name << '\t' << reference.record->name << '\t'
        << common << '\t' << only_query << '\t'
        << only_reference << '\t' << query_variance << '\t'
        << reference_variance << '\n';
}

void compare_multiple(
    const std::vector<FastaRecord>& queries,
    const std::vector<ReferenceView>& references,
    bool masking,
    bool repeat_score,
    std::ofstream& output
) {
    KmerIndex kmer_index;
    std::vector<std::vector<std::pair<std::uint32_t, KmerOccurrence>>>
        reference_counts(references.size());

    for (std::size_t reference_index = 0;
         reference_index < references.size(); ++reference_index) {
        const KmerCounts counts = reference_kmer_counts(
            references[reference_index], masking);
        auto& indexed_counts = reference_counts[reference_index];
        indexed_counts.reserve(counts.size());
        for (const auto& [kmer, occurrence] : counts) {
            auto found = kmer_index.find(kmer);
            if (found == kmer_index.end()) {
                if (kmer_index.size() >
                    std::numeric_limits<std::uint32_t>::max()) {
                    throw std::runtime_error(
                        "more than UINT32_MAX distinct reference 31-mers");
                }
                const auto index = static_cast<std::uint32_t>(
                    kmer_index.size());
                kmer_index.emplace(
                    kmer, KmerIndexEntry{index, occurrence.occurrence});
                indexed_counts.emplace_back(index, occurrence);
            } else {
                found->second.occurrence = std::max(
                    found->second.occurrence, occurrence.occurrence);
                indexed_counts.emplace_back(found->second.index, occurrence);
            }
        }
    }

    const std::size_t kmer_count = kmer_index.size();
    std::vector<Score> weights(kmer_count, 0.0F);
    for (const auto& [kmer, entry] : kmer_index) {
        (void)kmer;
        weights[entry.index] = occurrence_weight(entry.occurrence, repeat_score);
    }

    // Keep only nonzero reference occurrences.  The former dense matrix used
    // references.size() * kmer_count cells and could require tens of GiB for
    // one repeat-rich graph even though nearly every cell was zero.
    using Posting = std::pair<std::uint32_t, KmerOccurrence>;
    // CSR-style offsets avoid one std::vector allocation/object per distinct
    // k-mer, which is itself costly for graphs with millions of unique k-mers.
    std::vector<std::size_t> posting_offsets(kmer_count + 1, 0);
    std::vector<Score> reference_variance(references.size(), 0.0F);
    for (std::size_t index = 0; index < references.size(); ++index) {
        for (const auto& [index_in_matrix, occurrence] :
             reference_counts[index]) {
            ++posting_offsets[index_in_matrix + 1];
            const Score value = effective_count(occurrence);
            reference_variance[index] +=
                value * value * weights[index_in_matrix];
        }
    }
    for (std::size_t index = 1; index < posting_offsets.size(); ++index) {
        posting_offsets[index] += posting_offsets[index - 1];
    }
    std::vector<Posting> reference_postings(posting_offsets.back());
    std::vector<std::size_t> posting_positions = posting_offsets;
    for (std::size_t index = 0; index < references.size(); ++index) {
        for (const auto& [index_in_matrix, occurrence] :
             reference_counts[index]) {
            reference_postings[posting_positions[index_in_matrix]++] = {
                static_cast<std::uint32_t>(index), occurrence,
            };
        }
    }
    posting_positions.clear();
    posting_positions.shrink_to_fit();
    // Release the construction copy before query processing.  Each nonzero
    // reference occurrence now exists exactly once in reference_postings.
    reference_counts.clear();
    reference_counts.shrink_to_fit();

    // Self variances need only one float32 value per record. Query k-mer maps
    // are still streamed so the much larger per-k-mer state is never retained
    // for all queries at once.
    std::vector<Score> query_variance(queries.size(), 0.0F);

    // Stream one query at a time.  The three Gram accumulators scale with the
    // number of reference records, not queries * references or
    // queries * distinct-kmers.  Sorting indexed query k-mers preserves the
    // old k-mer-index accumulation order and therefore its numeric behavior.
    for (std::size_t query_index = 0;
         query_index < queries.size(); ++query_index) {
        const KmerCounts counts = query_kmer_counts(
            queries[query_index], masking);
        std::vector<std::pair<std::uint32_t, KmerOccurrence>> indexed_query;
        indexed_query.reserve(counts.size());
        Score& this_query_variance = query_variance[query_index];
        for (const auto& [kmer, occurrence] : counts) {
            const auto found = kmer_index.find(kmer);
            if (found == kmer_index.end()) {
                this_query_variance += standalone_variance(occurrence, repeat_score);
            } else {
                indexed_query.emplace_back(found->second.index, occurrence);
                const Score value = effective_count(occurrence);
                this_query_variance +=
                    value * value * weights[found->second.index];
            }
        }
        std::sort(indexed_query.begin(), indexed_query.end(), [](
            const auto& left, const auto& right
        ) { return left.first < right.first; });

        std::vector<Score> common(references.size(), 0.0F);
        std::vector<Score> shared_query_variance(references.size(), 0.0F);
        std::vector<Score> shared_reference_variance(references.size(), 0.0F);
        for (const auto& [index_in_matrix, query_count] : indexed_query) {
            const Score query_value = effective_count(query_count);
            const Score query_contribution =
                query_value * query_value * weights[index_in_matrix];
            for (std::size_t posting = posting_offsets[index_in_matrix];
                 posting < posting_offsets[index_in_matrix + 1]; ++posting) {
                const auto& [reference_index, reference_count] =
                    reference_postings[posting];
                const Score reference_value = effective_count(reference_count);
                common[reference_index] +=
                    query_value * reference_value *
                    weights[index_in_matrix];
                shared_query_variance[reference_index] += query_contribution;
                shared_reference_variance[reference_index] +=
                    reference_value * reference_value *
                    weights[index_in_matrix];
            }
        }

        for (std::size_t reference_index = 0;
             reference_index < references.size(); ++reference_index) {
            const Score only_query = std::max(
                0.0F,
                this_query_variance -
                    shared_query_variance[reference_index]);
            const Score only_reference = std::max(
                0.0F,
                reference_variance[reference_index] -
                    shared_reference_variance[reference_index]);
            output
                << queries[query_index].name << '\t'
                << references[reference_index].record->name << '\t'
                << common[reference_index] << '\t'
                << only_query << '\t'
                << only_reference << '\t' << this_query_variance << '\t'
                << reference_variance[reference_index] << '\n';
        }
    }
}

void run(const Options& options) {
    const auto queries = read_fasta(options.input_path);
    std::vector<FastaRecord> separate_reference_records;
    const std::vector<FastaRecord>* reference_records = &queries;
    if (options.reference_path != options.input_path) {
        separate_reference_records = read_fasta(options.reference_path);
        reference_records = &separate_reference_records;
    }
    const auto references = select_references(
        *reference_records, options.bed_path);

    std::ofstream output(options.output_path);
    if (!output) {
        throw std::runtime_error("cannot open output report: " + options.output_path);
    }
    output << std::setprecision(std::numeric_limits<Score>::max_digits10);
    write_header(output);
    if (options.paired) {
        if (queries.size() != references.size()) {
            throw std::runtime_error(
                "--paired requires equal input and reference record counts");
        }
        for (std::size_t index = 0; index < queries.size(); ++index) {
            compare_pairwise_sorted(
                queries[index], references[index], options.masking,
                options.repeat_score, output);
        }
    } else if (queries.size() == 1 && references.size() == 1) {
        compare_single(queries.front(), references.front(), options.masking, options.repeat_score, output);
    } else {
        compare_multiple(queries, references, options.masking, options.repeat_score, output);
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        run(parse_options(argc, argv));
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        usage(argv[0]);
        return EXIT_FAILURE;
    }
}
