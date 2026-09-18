#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <regex>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr std::uint32_t kKmerLength = 31;
constexpr std::uint32_t kSharedThreshold = 1000;
constexpr std::uint32_t kDefaultGap = 300;

struct Options {
    std::string fasta;
    std::string target_fasta;
    std::string bed;
    std::string output_bed;
    std::string groups_tsv;
    std::string pairs_tsv;
    std::uint32_t max_partitions = 20;
    std::uint32_t threads = std::max(1U, std::thread::hardware_concurrency());
    std::uint32_t gap = kDefaultGap;
    std::uint32_t shared_threshold = kSharedThreshold;
    bool alignment_only = false;
    bool only_kmer = false;
};

struct Interval {
    std::uint32_t start = 0;
    std::uint32_t end = 0;
};

struct Partition {
    std::string name;
    std::unordered_map<std::string, std::vector<Interval>> intervals;
    std::string source_contig;
    Interval source_interval;
};

struct Group {
    std::vector<std::uint32_t> original_members;
    // Lifted BED intervals remain separate because only they participate in
    // the aligned-reference overlap rule and the output BED.
    std::unordered_map<std::string, std::vector<Interval>> intervals;
    // Lifted half of the k-mer input. The other half is each original member's
    // complete -T target record, validated against its decoded source span.
    std::unordered_map<std::string, std::vector<Interval>> kmer_intervals;
};

class DisjointSet {
public:
    explicit DisjointSet(std::size_t size) : parent_(size) {
        std::iota(parent_.begin(), parent_.end(), 0U);
    }

    std::uint32_t find(std::uint32_t value) {
        std::uint32_t root = value;
        while (parent_[root] != root) root = parent_[root];
        while (parent_[value] != value) {
            const std::uint32_t next = parent_[value];
            parent_[value] = root;
            value = next;
        }
        return root;
    }

    void unite(std::uint32_t left, std::uint32_t right) {
        left = find(left);
        right = find(right);
        if (left == right) return;
        const std::uint32_t root = std::min(left, right);
        parent_[left] = root;
        parent_[right] = root;
    }

private:
    std::vector<std::uint32_t> parent_;
};

#pragma pack(push, 1)
struct KmerAssociation {
    std::uint64_t kmer;
    std::uint16_t partition;
};
#pragma pack(pop)
static_assert(sizeof(KmerAssociation) == 10, "association must remain packed");

[[noreturn]] void usage(const char* program, const std::string& error = "") {
    if (!error.empty()) std::cerr << program << ": error: " << error << '\n';
    std::cerr
        << "Usage: " << program
        << " -f reference.fa -T block_targets.fa -b mappings.filtered.bed"
        << " -o merged.bed [options]\n"
        << "  -f PATH       FASTA containing the lifted BED column-1 contigs\n"
        << "  -T PATH       Block-target FASTA; enumerates targets with no BED lift and\n"
        << "                supplies each column-4 region's exact sequence\n"
        << "  -b PATH       16-column filtered mapping BED; column 4 is partition\n"
        << "  -o PATH       Output merged BED4\n"
        << "  -g PATH       Output original-to-merged TSV (default: OUTPUT.groups.tsv)\n"
        << "  -p PATH       Optional shared unmasked-k-mer pair report\n"
        << "  -m INT        Ignore k-mers in more than INT partitions (default: 20)\n"
        << "  -t INT        K-mer extraction threads (default: hardware threads)\n"
        << "  --gap INT     Merge region gaps strictly below INT (default: 300)\n"
        << "  --shared INT  Strict shared-base/k-mer merge threshold (default: 1000)\n"
        << "  --alignment-only\n"
        << "                Merge only by aligned interval overlap; skip all k-mer work\n"
        << "  --only-kmer\n"
        << "                Merge only by shared unmasked k-mers; skip aligned-base merging\n";
    std::exit(error.empty() ? 0 : 1);
}

std::uint32_t parse_u32(const std::string& text, const std::string& option,
                        bool allow_zero = false) {
    std::size_t consumed = 0;
    const unsigned long long value = std::stoull(text, &consumed);
    if (consumed != text.size() || value > std::numeric_limits<std::uint32_t>::max() ||
        (!allow_zero && value == 0)) {
        throw std::runtime_error("invalid value for " + option + ": " + text);
    }
    return static_cast<std::uint32_t>(value);
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        auto value = [&]() -> std::string {
            if (++index >= argc) throw std::runtime_error("missing value after " + argument);
            return argv[index];
        };
        if (argument == "-f") options.fasta = value();
        else if (argument == "-T" || argument == "--target-fasta") {
            options.target_fasta = value();
        }
        else if (argument == "-b") options.bed = value();
        else if (argument == "-o") options.output_bed = value();
        else if (argument == "-g") options.groups_tsv = value();
        else if (argument == "-p") options.pairs_tsv = value();
        else if (argument == "-m") options.max_partitions = parse_u32(value(), argument);
        else if (argument == "-t") options.threads = parse_u32(value(), argument);
        else if (argument == "--gap") options.gap = parse_u32(value(), argument, true);
        else if (argument == "--shared") {
            options.shared_threshold = parse_u32(value(), argument);
        } else if (argument == "--alignment-only") {
            options.alignment_only = true;
        } else if (argument == "--only-kmer") {
            options.only_kmer = true;
        } else if (argument == "-h" || argument == "--help") {
            usage(argv[0]);
        } else {
            throw std::runtime_error("unknown option: " + argument);
        }
    }
    if (options.fasta.empty() || options.target_fasta.empty() ||
        options.bed.empty() || options.output_bed.empty()) {
        throw std::runtime_error("-f, -T, -b, and -o are required");
    }
    if (options.groups_tsv.empty()) options.groups_tsv = options.output_bed + ".groups.tsv";
    if (options.max_partitions > std::numeric_limits<std::uint16_t>::max()) {
        throw std::runtime_error("-m cannot exceed 65535");
    }
    if (options.shared_threshold >= std::numeric_limits<std::uint16_t>::max()) {
        throw std::runtime_error("--shared must be below 65535 for uint16 counters");
    }
    if (options.alignment_only && options.only_kmer) {
        throw std::runtime_error(
            "--alignment-only and --only-kmer are mutually exclusive");
    }
    if (options.alignment_only && !options.pairs_tsv.empty()) {
        throw std::runtime_error("-p cannot be used with --alignment-only");
    }
    return options;
}

std::vector<std::string> split_tab(const std::string& line) {
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

std::uint32_t checked_coordinate(const std::string& text,
                                 const std::string& context) {
    std::size_t consumed = 0;
    const unsigned long long value = std::stoull(text, &consumed);
    if (consumed != text.size() || value > std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error(context + ": coordinate exceeds uint32: " + text);
    }
    return static_cast<std::uint32_t>(value);
}

std::vector<std::string> split_underscore(const std::string& text) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (true) {
        const std::size_t separator = text.find('_', start);
        fields.push_back(text.substr(start, separator - start));
        if (separator == std::string::npos) break;
        start = separator + 1;
    }
    return fields;
}

std::pair<std::string, Interval> source_interval_from_partition_name(
    const std::string& name, const std::string& context) {
    // Required format:
    //   gN_SAMPLE_HAP_CONTIG[_CONTIG...]_START_END
    // Remove gN, SAMPLE, and HAP; the last two fields are coordinates and all
    // intervening fields, joined by underscores, are the true contig name.
    const auto fields = split_underscore(name);
    if (fields.size() < 6 || fields[0].size() < 2 || fields[0][0] != 'g' ||
        !std::all_of(fields[0].begin() + 1, fields[0].end(), [](char value) {
            return std::isdigit(static_cast<unsigned char>(value));
        }) || fields[1].empty() || fields[2].empty()) {
        throw std::runtime_error(
            context + ": column-4 partition must be "
            "gN_SAMPLE_HAP_CONTIG_START_END, got " + name);
    }
    std::string contig;
    for (std::size_t index = 3; index + 2 < fields.size(); ++index) {
        if (fields[index].empty()) {
            throw std::runtime_error(context + ": empty contig field in " + name);
        }
        if (!contig.empty()) contig.push_back('_');
        contig += fields[index];
    }
    if (contig.empty()) {
        throw std::runtime_error(context + ": empty source contig in " + name);
    }
    const std::uint32_t start = checked_coordinate(
        fields[fields.size() - 2], context + " column-4 source start");
    const std::uint32_t end = checked_coordinate(
        fields.back(), context + " column-4 source end");
    if (end <= start) {
        throw std::runtime_error(context + ": invalid column-4 source interval " + name);
    }
    return {contig, Interval{start, end}};
}

std::unordered_map<std::string, std::string> read_fasta(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open FASTA " + path);
    std::unordered_map<std::string, std::string> sequences;
    std::string current;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line.front() == '>') {
            const std::size_t end = line.find_first_of(" \t", 1);
            current = line.substr(1, end == std::string::npos ? end : end - 1);
            if (current.empty()) throw std::runtime_error(path + ": empty FASTA name");
            if (!sequences.emplace(current, std::string()).second) {
                throw std::runtime_error(path + ": duplicate FASTA record " + current);
            }
            continue;
        }
        if (current.empty()) {
            if (!line.empty()) throw std::runtime_error(path + ": sequence before header");
            continue;
        }
        auto& sequence = sequences[current];
        for (const char base : line) {
            if (!std::isspace(static_cast<unsigned char>(base))) sequence.push_back(base);
        }
    }
    if (sequences.empty()) throw std::runtime_error(path + ": no FASTA records");
    return sequences;
}

std::vector<std::string> read_fasta_names(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open FASTA " + path);
    std::vector<std::string> names;
    std::unordered_set<std::string> seen;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line.front() != '>') continue;
        const std::size_t end = line.find_first_of(" \t", 1);
        const std::string name = line.substr(
            1, end == std::string::npos ? end : end - 1);
        if (name.empty()) throw std::runtime_error(path + ": empty FASTA name");
        if (!seen.insert(name).second) {
            throw std::runtime_error(path + ": duplicate FASTA record " + name);
        }
        names.push_back(name);
    }
    if (names.empty()) throw std::runtime_error(path + ": no FASTA records");
    return names;
}

std::vector<Partition> read_bed(
    const std::string& path, const std::vector<std::string>& target_names) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open BED " + path);
    std::vector<Partition> partitions;
    std::unordered_map<std::string, std::uint32_t> indexes;
    partitions.reserve(target_names.size());
    for (const std::string& name : target_names) {
        const auto source = source_interval_from_partition_name(
            name, "target FASTA record " + name);
        indexes.emplace(name, static_cast<std::uint32_t>(partitions.size()));
        partitions.push_back(Partition{
            name, {}, source.first, source.second,
        });
    }
    std::string line;
    std::uint64_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line.front() == '#') continue;
        const auto fields = split_tab(line);
        if (fields.size() < 16) {
            throw std::runtime_error(path + ":" + std::to_string(line_number) +
                                     ": expected 16-column filtered BED");
        }
        const std::uint32_t start = checked_coordinate(
            fields[1], path + ":" + std::to_string(line_number));
        const std::uint32_t end = checked_coordinate(
            fields[2], path + ":" + std::to_string(line_number));
        if (end <= start || fields[0].empty() || fields[3].empty()) {
            throw std::runtime_error(path + ":" + std::to_string(line_number) +
                                     ": invalid BED interval or target");
        }
        auto found = indexes.find(fields[3]);
        std::uint32_t partition_index;
        if (found == indexes.end()) {
            throw std::runtime_error(
                path + ":" + std::to_string(line_number) +
                ": column-4 target absent from -T FASTA: " + fields[3]);
        } else {
            partition_index = found->second;
        }
        partitions[partition_index].intervals[fields[0]].push_back({start, end});
    }
    if (partitions.empty()) throw std::runtime_error(path + ": no mapping rows");
    if (partitions.size() > std::numeric_limits<std::uint16_t>::max() - 1ULL) {
        throw std::runtime_error("partition count exceeds 65534");
    }
    return partitions;
}

std::vector<Interval> merge_intervals(std::vector<Interval> intervals,
                                      std::uint32_t gap) {
    if (intervals.empty()) return {};
    std::sort(intervals.begin(), intervals.end(), [](const Interval& left,
                                                      const Interval& right) {
        return left.start < right.start ||
               (left.start == right.start && left.end < right.end);
    });
    std::vector<Interval> merged;
    merged.reserve(intervals.size());
    merged.push_back(intervals.front());
    for (std::size_t index = 1; index < intervals.size(); ++index) {
        Interval& current = merged.back();
        const Interval next = intervals[index];
        const std::uint64_t allowed_end = static_cast<std::uint64_t>(current.end) + gap;
        if (next.start <= current.end || (gap > 0 && next.start < allowed_end)) {
            current.end = std::max(current.end, next.end);
        } else {
            merged.push_back(next);
        }
    }
    return merged;
}

void merge_partition_intervals(std::vector<Partition>& partitions,
                               std::uint32_t gap) {
    for (auto& partition : partitions) {
        for (auto& entry : partition.intervals) {
            entry.second = merge_intervals(std::move(entry.second), gap);
        }
    }
}

std::uint64_t pair_key(std::uint32_t left, std::uint32_t right) {
    if (left > right) std::swap(left, right);
    return (static_cast<std::uint64_t>(left) << 32U) | right;
}

std::pair<std::uint32_t, std::uint32_t> decode_pair(std::uint64_t key) {
    return {static_cast<std::uint32_t>(key >> 32U),
            static_cast<std::uint32_t>(key)};
}

void merge_by_shared_unmasked_regions(
    const std::vector<Partition>& partitions,
    const std::unordered_map<std::string, std::string>& sequences,
    std::uint32_t threshold,
    DisjointSet& groups) {
    struct EventInterval {
        Interval interval;
        std::uint32_t partition;
    };
    std::unordered_map<std::string, std::vector<EventInterval>> by_contig;
    for (std::uint32_t partition = 0; partition < partitions.size(); ++partition) {
        for (const auto& entry : partitions[partition].intervals) {
            for (const Interval interval : entry.second) {
                by_contig[entry.first].push_back({interval, partition});
            }
        }
    }

    std::unordered_map<std::uint64_t, std::uint32_t> shared;
    for (auto& contig_entry : by_contig) {
        const auto sequence_found = sequences.find(contig_entry.first);
        if (sequence_found == sequences.end()) {
            throw std::runtime_error("BED contig absent from FASTA: " + contig_entry.first);
        }
        const std::string& sequence = sequence_found->second;
        std::vector<std::uint32_t> prefix(sequence.size() + 1, 0);
        for (std::size_t position = 0; position < sequence.size(); ++position) {
            const char base = sequence[position];
            const bool unmasked = base == 'A' || base == 'C' || base == 'G' || base == 'T';
            prefix[position + 1] = prefix[position] + static_cast<std::uint32_t>(unmasked);
        }
        auto& intervals = contig_entry.second;
        std::sort(intervals.begin(), intervals.end(), [](const EventInterval& left,
                                                         const EventInterval& right) {
            if (left.interval.start != right.interval.start) {
                return left.interval.start < right.interval.start;
            }
            return left.interval.end < right.interval.end;
        });
        std::vector<EventInterval> active;
        for (const EventInterval current : intervals) {
            if (current.interval.end > sequence.size()) {
                throw std::runtime_error("BED interval exceeds FASTA contig " +
                                         contig_entry.first);
            }
            active.erase(std::remove_if(active.begin(), active.end(),
                                        [&](const EventInterval& prior) {
                return prior.interval.end <= current.interval.start;
            }), active.end());
            for (const EventInterval prior : active) {
                if (prior.partition == current.partition) continue;
                const std::uint32_t start = std::max(prior.interval.start,
                                                     current.interval.start);
                const std::uint32_t end = std::min(prior.interval.end,
                                                   current.interval.end);
                if (end <= start) continue;
                const std::uint32_t amount = prefix[end] - prefix[start];
                if (amount == 0) continue;
                std::uint32_t& count = shared[pair_key(prior.partition,
                                                       current.partition)];
                if (count <= threshold) {
                    count = std::min<std::uint32_t>(threshold + 1, count + amount);
                }
            }
            active.push_back(current);
        }
    }
    for (const auto& entry : shared) {
        if (entry.second > threshold) {
            const auto pair = decode_pair(entry.first);
            groups.unite(pair.first, pair.second);
        }
    }
}

std::vector<Group> build_groups(const std::vector<Partition>& partitions,
                                DisjointSet& membership,
                                std::uint32_t gap) {
    std::unordered_map<std::uint32_t, std::uint32_t> root_to_group;
    std::vector<Group> groups;
    for (std::uint32_t partition = 0; partition < partitions.size(); ++partition) {
        const std::uint32_t root = membership.find(partition);
        auto found = root_to_group.find(root);
        std::uint32_t group_index;
        if (found == root_to_group.end()) {
            group_index = static_cast<std::uint32_t>(groups.size());
            root_to_group.emplace(root, group_index);
            groups.push_back(Group{});
        } else {
            group_index = found->second;
        }
        Group& group = groups[group_index];
        group.original_members.push_back(partition);
        for (const auto& entry : partitions[partition].intervals) {
            auto& destination = group.intervals[entry.first];
            destination.insert(destination.end(), entry.second.begin(), entry.second.end());
            auto& kmer_destination = group.kmer_intervals[entry.first];
            kmer_destination.insert(
                kmer_destination.end(), entry.second.begin(), entry.second.end());
        }
    }
    for (auto& group : groups) {
        for (auto& entry : group.intervals) {
            entry.second = merge_intervals(std::move(entry.second), gap);
        }
        for (auto& entry : group.kmer_intervals) {
            entry.second = merge_intervals(std::move(entry.second), gap);
        }
    }
    return groups;
}

bool encode_base(char base, std::uint64_t& value, bool& masked) {
    masked = base >= 'a' && base <= 'z';
    const char upper = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    switch (upper) {
        case 'A': value = 0; return true;
        case 'C': value = 1; return true;
        case 'G': value = 2; return true;
        case 'T': value = 3; return true;
        default: return false;
    }
}

void collect_interval_kmers(const std::string& sequence, Interval interval,
                            std::unordered_set<std::uint64_t>& output) {
    if (interval.end > sequence.size()) {
        throw std::runtime_error("merged BED interval exceeds FASTA sequence");
    }
    constexpr std::uint64_t sequence_mask = (1ULL << (2U * kKmerLength)) - 1ULL;
    constexpr std::uint64_t lowercase_mask = (1ULL << kKmerLength) - 1ULL;
    std::uint64_t forward = 0;
    std::uint64_t reverse = 0;
    std::uint64_t masked_window = 0;
    std::uint32_t valid = 0;
    for (std::uint32_t position = interval.start; position < interval.end; ++position) {
        std::uint64_t encoded = 0;
        bool masked = false;
        if (!encode_base(sequence[position], encoded, masked)) {
            forward = reverse = masked_window = 0;
            valid = 0;
            continue;
        }
        forward = ((forward << 2U) | encoded) & sequence_mask;
        reverse = (reverse >> 2U) |
                  ((3ULL - encoded) << (2U * (kKmerLength - 1U)));
        masked_window = ((masked_window << 1U) | static_cast<std::uint64_t>(masked)) &
                        lowercase_mask;
        if (valid < kKmerLength) ++valid;
        // A soft-masked base invalidates every 31-mer window containing it.
        // Retain only windows composed entirely of uppercase A/C/G/T.
        if (valid >= kKmerLength && masked_window == 0) {
            output.insert(std::max(forward, reverse));
        }
    }
}

std::vector<KmerAssociation> collect_unmasked_kmers(
    const std::vector<Group>& groups,
    const std::vector<Partition>& partitions,
    const std::unordered_map<std::string, std::string>& sequences,
    const std::unordered_map<std::string, std::string>& target_sequences,
    std::uint32_t threads) {
    if (groups.size() > std::numeric_limits<std::uint16_t>::max() - 1ULL) {
        throw std::runtime_error("merged partition count exceeds 65534");
    }
    threads = std::max<std::uint32_t>(1, std::min<std::uint32_t>(
        threads, static_cast<std::uint32_t>(groups.size())));
    std::atomic<std::uint32_t> next{0};
    std::vector<std::vector<KmerAssociation>> local(threads);
    std::vector<std::thread> workers;
    std::mutex error_mutex;
    std::exception_ptr error;
    for (std::uint32_t worker = 0; worker < threads; ++worker) {
        workers.emplace_back([&, worker]() {
            try {
                while (true) {
                    const std::uint32_t group_index = next.fetch_add(1);
                    if (group_index >= groups.size()) break;
                    std::unordered_set<std::uint64_t> unique;
                    for (const auto& entry : groups[group_index].kmer_intervals) {
                        const auto found = sequences.find(entry.first);
                        if (found == sequences.end()) {
                            throw std::runtime_error(
                                "lifted k-mer source contig absent from -f FASTA: " +
                                entry.first);
                        }
                        for (const Interval interval : entry.second) {
                            collect_interval_kmers(found->second, interval, unique);
                        }
                    }
                    for (const std::uint32_t member :
                         groups[group_index].original_members) {
                        const std::string& target = partitions[member].name;
                        const auto found = target_sequences.find(target);
                        if (found == target_sequences.end()) {
                            throw std::runtime_error(
                                "partition absent from -T FASTA: " + target);
                        }
                        if (found->second.size() >
                            std::numeric_limits<std::uint32_t>::max()) {
                            throw std::runtime_error(
                                "target sequence exceeds uint32 length: " + target);
                        }
                        collect_interval_kmers(
                            found->second,
                            Interval{0, static_cast<std::uint32_t>(found->second.size())},
                            unique);
                    }
                    auto& output = local[worker];
                    output.reserve(output.size() + unique.size());
                    for (const std::uint64_t kmer : unique) {
                        output.push_back({kmer, static_cast<std::uint16_t>(group_index)});
                    }
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!error) error = std::current_exception();
                next.store(static_cast<std::uint32_t>(groups.size()));
            }
        });
    }
    for (auto& worker : workers) worker.join();
    if (error) std::rethrow_exception(error);
    std::size_t total = 0;
    for (const auto& output : local) total += output.size();
    std::vector<KmerAssociation> associations;
    associations.reserve(total);
    for (auto& output : local) {
        associations.insert(associations.end(), output.begin(), output.end());
        std::vector<KmerAssociation>().swap(output);
    }
    return associations;
}

std::uint64_t triangular_index(std::uint32_t left, std::uint32_t right) {
    if (left > right) std::swap(left, right);
    return static_cast<std::uint64_t>(right) * (right + 1ULL) / 2ULL + left;
}

std::vector<std::uint16_t> count_shared_kmers(
    std::vector<KmerAssociation>& associations,
    std::size_t group_count,
    std::uint32_t max_partitions) {
    std::sort(associations.begin(), associations.end(),
              [](const KmerAssociation& left, const KmerAssociation& right) {
        return left.kmer < right.kmer ||
               (left.kmer == right.kmer && left.partition < right.partition);
    });
    const std::uint64_t pair_count = group_count * (group_count + 1ULL) / 2ULL;
    if (pair_count > std::numeric_limits<std::size_t>::max()) {
        throw std::runtime_error("partition pair matrix exceeds address space");
    }
    std::vector<std::uint16_t> counts(static_cast<std::size_t>(pair_count), 0);
    std::size_t begin = 0;
    while (begin < associations.size()) {
        std::size_t end = begin + 1;
        while (end < associations.size() &&
               associations[end].kmer == associations[begin].kmer) ++end;
        const std::size_t association_count = end - begin;
        if (association_count <= max_partitions) {
            for (std::size_t left = begin; left < end; ++left) {
                for (std::size_t right = left + 1; right < end; ++right) {
                    std::uint16_t& count = counts[triangular_index(
                        associations[left].partition,
                        associations[right].partition)];
                    if (count != std::numeric_limits<std::uint16_t>::max()) ++count;
                }
            }
        }
        begin = end;
    }
    return counts;
}

std::string short_label(const std::string& name) {
    static const std::regex leading_group("^(g[0-9]+)(?:_|$)");
    std::smatch match;
    if (std::regex_search(name, match, leading_group)) return match[1].str();
    if (!name.empty() && std::all_of(name.begin(), name.end(), [](char value) {
            return std::isdigit(static_cast<unsigned char>(value));
        })) return "g" + name;
    return name;
}

std::string merged_name(const std::vector<std::uint32_t>& members,
                        const std::vector<Partition>& partitions) {
    if (members.size() == 1) return partitions[members.front()].name;
    std::string name = "merged_";
    for (const std::uint32_t member : members) name += short_label(partitions[member].name);
    return name;
}

void write_outputs(const Options& options,
                   const std::vector<Partition>& partitions,
                   const std::vector<Group>& first_groups,
                   DisjointSet& final_membership,
                   const std::vector<std::uint16_t>& shared_counts) {
    std::unordered_map<std::uint32_t, std::uint32_t> root_to_final;
    std::vector<Group> final_groups;
    for (std::uint32_t first = 0; first < first_groups.size(); ++first) {
        const std::uint32_t root = final_membership.find(first);
        auto found = root_to_final.find(root);
        std::uint32_t final_index;
        if (found == root_to_final.end()) {
            final_index = static_cast<std::uint32_t>(final_groups.size());
            root_to_final.emplace(root, final_index);
            final_groups.push_back(Group{});
        } else {
            final_index = found->second;
        }
        Group& destination = final_groups[final_index];
        destination.original_members.insert(destination.original_members.end(),
            first_groups[first].original_members.begin(),
            first_groups[first].original_members.end());
        for (const auto& entry : first_groups[first].intervals) {
            auto& intervals = destination.intervals[entry.first];
            intervals.insert(intervals.end(), entry.second.begin(), entry.second.end());
        }
    }
    for (auto& group : final_groups) {
        std::sort(group.original_members.begin(), group.original_members.end());
        for (auto& entry : group.intervals) {
            entry.second = merge_intervals(std::move(entry.second), options.gap);
        }
    }
    std::sort(final_groups.begin(), final_groups.end(), [](const Group& left,
                                                           const Group& right) {
        return left.original_members.front() < right.original_members.front();
    });

    std::ofstream bed(options.output_bed);
    if (!bed) throw std::runtime_error("cannot write " + options.output_bed);
    std::ofstream mapping(options.groups_tsv);
    if (!mapping) throw std::runtime_error("cannot write " + options.groups_tsv);
    mapping << "partition\tmerged_partition\n";
    for (const Group& group : final_groups) {
        const std::string name = merged_name(group.original_members, partitions);
        for (const std::uint32_t member : group.original_members) {
            mapping << partitions[member].name << '\t' << name << '\n';
        }
        std::vector<std::string> contigs;
        contigs.reserve(group.intervals.size());
        for (const auto& entry : group.intervals) contigs.push_back(entry.first);
        std::sort(contigs.begin(), contigs.end());
        for (const std::string& contig : contigs) {
            for (const Interval interval : group.intervals.at(contig)) {
                bed << contig << '\t' << interval.start << '\t' << interval.end
                    << '\t' << name << '\n';
            }
        }
    }
    if (!bed || !mapping) throw std::runtime_error("failed while writing outputs");

    if (!options.pairs_tsv.empty()) {
        std::ofstream pairs(options.pairs_tsv);
        if (!pairs) throw std::runtime_error("cannot write " + options.pairs_tsv);
        pairs << "partition_a\tpartition_b\tshared_unmasked_kmers\n";
        for (std::uint32_t right = 1; right < first_groups.size(); ++right) {
            for (std::uint32_t left = 0; left < right; ++left) {
                const std::uint16_t count = shared_counts[triangular_index(left, right)];
                if (count == 0) continue;
                pairs << merged_name(first_groups[left].original_members, partitions)
                      << '\t'
                      << merged_name(first_groups[right].original_members, partitions)
                      << '\t' << count << '\n';
            }
        }
    }
}

int run(const Options& options) {
    std::cerr << "Loading FASTA...\n";
    const auto sequences = read_fasta(options.fasta);
    std::cerr << "Loading block-target FASTA...\n";
    const auto target_names = read_fasta_names(options.target_fasta);
    const auto target_sequences = read_fasta(options.target_fasta);
    std::cerr << "Loading filtered BED...\n";
    auto partitions = read_bed(options.bed, target_names);
    for (const Partition& partition : partitions) {
        const auto found = target_sequences.find(partition.name);
        if (found == target_sequences.end()) {
            throw std::runtime_error(
                "partition absent from -T FASTA: " + partition.name);
        }
        const std::uint64_t source_length =
            static_cast<std::uint64_t>(partition.source_interval.end) -
            partition.source_interval.start;
        if (found->second.size() != source_length) {
            throw std::runtime_error(
                "column-4 source span differs from -T sequence length for " +
                partition.name + ": span=" + std::to_string(source_length) +
                " sequence=" + std::to_string(found->second.size()));
        }
    }
    merge_partition_intervals(partitions, options.gap);

    DisjointSet region_membership(partitions.size());
    if (options.only_kmer) {
        std::cerr << "K-mer-only mode: skipping aligned-reference-base merge rule\n";
    } else {
        std::cerr << "Merging partitions sharing more than " << options.shared_threshold
                  << " unmasked reference bases...\n";
        merge_by_shared_unmasked_regions(partitions, sequences,
                                         options.shared_threshold, region_membership);
    }
    auto first_groups = build_groups(partitions, region_membership, options.gap);
    if (options.only_kmer) {
        std::cerr << "K-mer input: " << first_groups.size()
                  << " original partitions\n";
    } else {
        std::cerr << "Region rule: " << partitions.size() << " partitions -> "
                  << first_groups.size() << " groups\n";
    }

    DisjointSet final_membership(first_groups.size());
    std::vector<std::uint16_t> shared_counts;
    if (options.alignment_only) {
        std::cerr << "Alignment-only mode: skipping unmasked k-mer collection and merge rule\n";
    } else {
        std::cerr << "Collecting canonical unmasked 31-mers with " << options.threads
                  << " threads...\n";
        auto associations = collect_unmasked_kmers(
            first_groups, partitions, sequences, target_sequences, options.threads);
        std::cerr << "Sorting " << associations.size()
                  << " k-mer/partition associations...\n";
        shared_counts = count_shared_kmers(associations, first_groups.size(),
                                           options.max_partitions);
        std::vector<KmerAssociation>().swap(associations);

        for (std::uint32_t right = 1; right < first_groups.size(); ++right) {
            for (std::uint32_t left = 0; left < right; ++left) {
                if (shared_counts[triangular_index(left, right)] >
                    options.shared_threshold) {
                    final_membership.unite(left, right);
                }
            }
        }
    }
    write_outputs(options, partitions, first_groups, final_membership, shared_counts);
    std::cerr << "Wrote " << options.output_bed << " and " << options.groups_tsv << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(parse_options(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << "kmerpartitionmerger: error: " << error.what() << '\n';
        return 1;
    }
}
