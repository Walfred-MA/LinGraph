#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
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

#include "htslib/faidx.h"

using std::size_t;
using std::string;
using std::uint32_t;
using std::uint64_t;
using std::vector;

namespace {

constexpr int KMER_SIZE = 31;
constexpr uint64_t KMER_MASK = (uint64_t{1} << (2 * KMER_SIZE)) - 1;
constexpr size_t DEFAULT_CHUNK_SIZE = 4ULL * 1024ULL * 1024ULL;
constexpr size_t RDATA_FLUSH_SIZE = 1'000'000;

struct Options {
    string reference_fasta;
    string reference_bed;
    string query_fasta;
    string query_bed;
    string output_bed;
    int threads = std::max(1u, std::thread::hardware_concurrency());
    size_t hotspot_kmers = 500;
    uint64_t hotspot_window = 1000;
    uint64_t merge_distance = 500;
    size_t chunk_size = DEFAULT_CHUNK_SIZE;
};

struct Region {
    string chrom;
    hts_pos_t start = 0;   // BED: 0-based, half-open
    hts_pos_t end = 0;
    string name;
    uint32_t index = 0;

    uint64_t length() const {
        return static_cast<uint64_t>(end - start);
    }
};

struct Hotspot {
    hts_pos_t start = 0;
    hts_pos_t end = 0;     // BED: half-open
    vector<uint64_t> scores;
};

struct BedOutput {
    string chrom;
    hts_pos_t start = 0;
    hts_pos_t end = 0;
    string name;
    uint32_t best_region = 0;
};

struct FaidxDeleter {
    void operator()(faidx_t *fai) const noexcept {
        if (fai != nullptr) fai_destroy(fai);
    }
};
using FaidxPtr = std::unique_ptr<faidx_t, FaidxDeleter>;

void log(const string &message) {
    std::cerr << "[KmerRegionHotspot] " << message << '\n';
}

[[noreturn]] void usage(const char *program, int exit_code = 1) {
    std::ostream &out = exit_code == 0 ? std::cout : std::cerr;
    out
        << "Usage: " << program << " -r REF.fa -rb REF.bed -q QUERY.fa -o OUT.bed [options]\n\n"
        << "Required:\n"
        << "  -r FILE                 faidx-indexed reference FASTA\n"
        << "  -rb FILE                reference BED3+; BED4 name is used when present\n"
        << "  -q FILE                 de novo assembly FASTA (index is created if absent)\n"
        << "  -o FILE                 output BED6\n\n"
        << "Optional:\n"
        << "  -qb FILE                restrict query scanning to BED3+ intervals\n"
        << "  -t, --threads INT       worker threads [hardware concurrency]\n"
        << "  --hotspot-kmers INT     matching 31-mers required for a seed [500]\n"
        << "  --hotspot-window INT    seed window in bp [1000]\n"
        << "  --merge-distance INT    merge hotspot seeds within this many bp [500]\n"
        << "  --chunk-size INT        faidx fetch chunk size [4194304]\n"
        << "  -h, --help              show this help\n\n"
        << "K-mers are canonicalized exactly as in KmerSeacher.hpp by selecting\n"
        << "max(forward_2bit, reverse_complement_2bit). Lowercase A/C/G/T are kept;\n"
        << "non-ACGT characters break the rolling k-mer.\n";
    std::exit(exit_code);
}

uint64_t parse_u64(const string &text, const string &label) {
    size_t used = 0;
    unsigned long long value = 0;
    try {
        value = std::stoull(text, &used);
    } catch (...) {
        throw std::runtime_error("invalid " + label + ": " + text);
    }
    if (used != text.size()) throw std::runtime_error("invalid " + label + ": " + text);
    return static_cast<uint64_t>(value);
}

Options parse_args(int argc, char **argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];
        auto require_value = [&](const string &flag) -> string {
            if (i + 1 >= argc) throw std::runtime_error("missing value after " + flag);
            return argv[++i];
        };

        if (arg == "-r") opt.reference_fasta = require_value(arg);
        else if (arg == "-rb") opt.reference_bed = require_value(arg);
        else if (arg == "-q") opt.query_fasta = require_value(arg);
        else if (arg == "-qb") opt.query_bed = require_value(arg);
        else if (arg == "-o") opt.output_bed = require_value(arg);
        else if (arg == "-t" || arg == "--threads") {
            const uint64_t value = parse_u64(require_value(arg), "thread count");
            if (value == 0 || value > static_cast<uint64_t>(std::numeric_limits<int>::max()))
                throw std::runtime_error("--threads must be >= 1");
            opt.threads = static_cast<int>(value);
        } else if (arg == "--hotspot-kmers") {
            const uint64_t value = parse_u64(require_value(arg), "hotspot k-mer count");
            if (value == 0) throw std::runtime_error("--hotspot-kmers must be >= 1");
            opt.hotspot_kmers = static_cast<size_t>(value);
        } else if (arg == "--hotspot-window") {
            opt.hotspot_window = parse_u64(require_value(arg), "hotspot window");
            if (opt.hotspot_window == 0) throw std::runtime_error("--hotspot-window must be >= 1");
        } else if (arg == "--merge-distance") {
            opt.merge_distance = parse_u64(require_value(arg), "merge distance");
        } else if (arg == "--chunk-size") {
            const uint64_t value = parse_u64(require_value(arg), "chunk size");
            if (value == 0) throw std::runtime_error("--chunk-size must be >= 1");
            opt.chunk_size = static_cast<size_t>(value);
        } else if (arg == "-h" || arg == "--help") {
            usage(argv[0], 0);
        } else {
            throw std::runtime_error("unknown argument: " + arg);
        }
    }

    if (opt.reference_fasta.empty() || opt.reference_bed.empty() ||
        opt.query_fasta.empty() || opt.output_bed.empty()) {
        usage(argv[0]);
    }
    return opt;
}

FaidxPtr load_faidx(const string &fasta, bool allow_create) {
    const string fai_path = fasta + ".fai";
    FaidxPtr fai(allow_create
                     ? fai_load(fasta.c_str())
                     : fai_load3(fasta.c_str(), fai_path.c_str(), nullptr, 0));
    if (!fai) {
        if (allow_create) {
            throw std::runtime_error(
                "cannot load/create FASTA index " + fai_path +
                "; run: samtools faidx " + fasta);
        }
        throw std::runtime_error("cannot load required FASTA index: " + fai_path);
    }
    return fai;
}

vector<string> split_ws(const string &line) {
    std::istringstream in(line);
    vector<string> fields;
    string value;
    while (in >> value) fields.push_back(value);
    return fields;
}

vector<Region> read_bed_regions(
    const string &path,
    faidx_t *fai,
    bool require_reference_names,
    const string &generated_prefix)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open BED: " + path);

    vector<Region> regions;
    string line;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("track", 0) == 0 || line.rfind("browser", 0) == 0) continue;
        const vector<string> fields = split_ws(line);
        if (fields.size() < 3) {
            throw std::runtime_error(path + ":" + std::to_string(line_number) + ": expected BED3+");
        }

        const hts_pos_t contig_length = faidx_seq_len64(fai, fields[0].c_str());
        if (contig_length < 0) {
            throw std::runtime_error(path + ":" + std::to_string(line_number) +
                                     ": contig absent from FASTA: " + fields[0]);
        }

        const uint64_t start_u = parse_u64(fields[1], "BED start");
        const uint64_t end_u = parse_u64(fields[2], "BED end");
        if (end_u <= start_u) {
            throw std::runtime_error(path + ":" + std::to_string(line_number) +
                                     ": BED end must be greater than start");
        }
        if (end_u > static_cast<uint64_t>(contig_length)) {
            throw std::runtime_error(path + ":" + std::to_string(line_number) +
                                     ": BED interval exceeds contig length");
        }

        Region region;
        region.chrom = fields[0];
        region.start = static_cast<hts_pos_t>(start_u);
        region.end = static_cast<hts_pos_t>(end_u);
        region.index = static_cast<uint32_t>(regions.size());
        if (require_reference_names && fields.size() >= 4 && !fields[3].empty() && fields[3] != ".") {
            region.name = fields[3];
        } else if (require_reference_names) {
            region.name = generated_prefix + std::to_string(region.index);
        } else if (fields.size() >= 4 && !fields[3].empty()) {
            region.name = fields[3];
        } else {
            region.name = generated_prefix + std::to_string(region.index);
        }
        regions.push_back(std::move(region));
    }

    if (regions.empty()) throw std::runtime_error("no usable intervals in BED: " + path);
    return regions;
}

vector<Region> whole_contig_regions(faidx_t *fai) {
    const int nseq = faidx_nseq(fai);
    if (nseq < 0) throw std::runtime_error("cannot enumerate query FASTA contigs");
    vector<Region> regions;
    regions.reserve(static_cast<size_t>(nseq));
    for (int i = 0; i < nseq; ++i) {
        const char *name = faidx_iseq(fai, i);
        if (!name) throw std::runtime_error("invalid FASTA index entry");
        const hts_pos_t len = faidx_seq_len64(fai, name);
        if (len <= 0) continue;
        Region region;
        region.chrom = name;
        region.start = 0;
        region.end = len;
        region.name = name;
        region.index = static_cast<uint32_t>(regions.size());
        regions.push_back(std::move(region));
    }
    return regions;
}

inline bool base_to_int(char base, uint64_t &value) {
    switch (base) {
        case 'A': case 'a': value = 0; return true;
        case 'C': case 'c': value = 1; return true;
        case 'G': case 'g': value = 2; return true;
        case 'T': case 't': value = 3; return true;
        default: return false;
    }
}

struct RollingKmer31 {
    int current_size = 0;
    uint64_t forward = 0;
    uint64_t reverse = 0;

    void reset() {
        current_size = 0;
        forward = 0;
        reverse = 0;
    }

    bool push(char base, uint64_t &canonical) {
        uint64_t converted = 0;
        if (!base_to_int(base, converted)) {
            reset();
            return false;
        }

        ++current_size;
        forward = ((forward << 2) & KMER_MASK) | converted;
        reverse >>= 2;
        reverse |= (uint64_t{3} - converted) << (2 * KMER_SIZE - 2);

        if (current_size < KMER_SIZE) return false;
        // KmerSeacher.hpp uses the numerically larger canonical orientation.
        canonical = forward >= reverse ? forward : reverse;
        return true;
    }
};

template <class Consume>
void scan_region_kmers(
    faidx_t *fai,
    const Region &region,
    size_t chunk_size,
    Consume consume)
{
    RollingKmer31 rolling;
    for (hts_pos_t chunk_start = region.start;
         chunk_start < region.end;
         chunk_start += static_cast<hts_pos_t>(chunk_size)) {
        const hts_pos_t chunk_end_exclusive = std::min<hts_pos_t>(
            region.end, chunk_start + static_cast<hts_pos_t>(chunk_size));
        const hts_pos_t fetch_end = chunk_end_exclusive - 1;
        hts_pos_t fetched_length = 0;
        std::unique_ptr<char, decltype(&std::free)> sequence(
            faidx_fetch_seq64(fai, region.chrom.c_str(), chunk_start,
                              fetch_end, &fetched_length),
            &std::free);
        const hts_pos_t expected = chunk_end_exclusive - chunk_start;
        if (!sequence || fetched_length != expected) {
            throw std::runtime_error("cannot fetch " + region.chrom + ":" +
                                     std::to_string(static_cast<long long>(chunk_start)) + "-" +
                                     std::to_string(static_cast<long long>(chunk_end_exclusive)));
        }

        for (hts_pos_t offset = 0; offset < fetched_length; ++offset) {
            uint64_t kmer = 0;
            if (!rolling.push(sequence.get()[offset], kmer)) continue;
            const hts_pos_t end_position = chunk_start + offset; // inclusive base position
            const hts_pos_t kmer_start = end_position - (KMER_SIZE - 1);
            if (kmer_start < region.start) continue;
            consume(kmer_start, kmer);
        }
    }
}

template <class Work>
void parallel_for_regions(size_t n, int threads, Work work) {
    if (n == 0) return;
    const size_t worker_count = std::min<size_t>(n, static_cast<size_t>(std::max(1, threads)));
    std::atomic_size_t next{0};
    std::atomic_bool stop{false};
    std::mutex error_mutex;
    std::exception_ptr error;

    auto worker = [&](size_t worker_id) {
        try {
            while (!stop.load(std::memory_order_relaxed)) {
                const size_t index = next.fetch_add(1, std::memory_order_relaxed);
                if (index >= n) break;
                work(worker_id, index);
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!error) error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
        }
    };

    vector<std::thread> pool;
    pool.reserve(worker_count);
    for (size_t i = 0; i < worker_count; ++i) pool.emplace_back(worker, i);
    for (auto &thread : pool) thread.join();
    if (error) std::rethrow_exception(error);
}

std::unordered_map<uint64_t, uint32_t> build_kmer_index(
    const Options &opt,
    const vector<Region> &regions)
{
    const size_t worker_count = std::min<size_t>(
        regions.size(), static_cast<size_t>(std::max(1, opt.threads)));
    vector<std::unordered_set<uint64_t>> local_sets(worker_count);

    log("reference pass 1/2: collecting unique canonical 31-mers from " +
        std::to_string(regions.size()) + " BED regions using " +
        std::to_string(worker_count) + " threads");

    // Each worker owns both its faidx handle and unordered_set. No locking is
    // needed while scanning, and the final global unordered_map is immutable
    // during all subsequent parallel work.
    vector<FaidxPtr> worker_fai;
    worker_fai.reserve(worker_count);
    for (size_t i = 0; i < worker_count; ++i)
        worker_fai.push_back(load_faidx(opt.reference_fasta, false));

    parallel_for_regions(regions.size(), static_cast<int>(worker_count),
        [&](size_t worker_id, size_t region_index) {
            auto &local = local_sets[worker_id];
            scan_region_kmers(worker_fai[worker_id].get(), regions[region_index],
                              opt.chunk_size,
                              [&](hts_pos_t, uint64_t kmer) { local.insert(kmer); });
        });

    uint64_t rough_unique = 0;
    for (const auto &set : local_sets) rough_unique += set.size();
    std::unordered_map<uint64_t, uint32_t> kmer_index;
    kmer_index.reserve(static_cast<size_t>(rough_unique));

    uint64_t next_index = 0;
    for (auto &set : local_sets) {
        for (const uint64_t kmer : set) {
            auto inserted = kmer_index.emplace(kmer, 0);
            if (!inserted.second) continue;
            if (next_index > std::numeric_limits<uint32_t>::max())
                throw std::runtime_error("more than 2^32 unique reference k-mers; uint32 index exhausted");
            inserted.first->second = static_cast<uint32_t>(next_index++);
        }
        std::unordered_set<uint64_t>().swap(set);
    }

    log("unique reference 31-mers: " + std::to_string(kmer_index.size()));
    return kmer_index;
}

struct RDataBuilder {
    vector<uint32_t> data;
    std::mutex mutex;

    void append(vector<uint32_t> &buffer) {
        if (buffer.empty()) return;
        std::lock_guard<std::mutex> lock(mutex);
        const size_t needed = data.size() + buffer.size();
        if (needed > data.capacity()) {
            const size_t new_capacity = ((needed + RDATA_FLUSH_SIZE - 1) / RDATA_FLUSH_SIZE) * RDATA_FLUSH_SIZE;
            data.reserve(new_capacity);
        }
        data.insert(data.end(), buffer.begin(), buffer.end());
        buffer.clear();
    }
};

vector<uint32_t> build_rdata(
    const Options &opt,
    const vector<Region> &regions,
    const std::unordered_map<uint64_t, uint32_t> &kmer_index)
{
    const uint64_t num_regions = regions.size();
    const uint64_t num_kmers = kmer_index.size();
    const uint64_t uint32_space = uint64_t{1} << 32;
    if (num_regions == 0 || num_kmers > uint32_space / num_regions) {
        throw std::runtime_error(
            "packed rdata key kmer_index*num_regions+region_index does not fit uint32: "
            "unique_kmers=" + std::to_string(num_kmers) +
            ", regions=" + std::to_string(num_regions) +
            ". Use a uint64 rdata representation for this dataset.");
    }

    const size_t worker_count = std::min<size_t>(
        regions.size(), static_cast<size_t>(std::max(1, opt.threads)));
    vector<FaidxPtr> worker_fai;
    worker_fai.reserve(worker_count);
    for (size_t i = 0; i < worker_count; ++i)
        worker_fai.push_back(load_faidx(opt.reference_fasta, false));

    RDataBuilder builder;
    log("reference pass 2/2: recording k-mer occurrences in 1,000,000-entry rdata blocks");

    parallel_for_regions(regions.size(), static_cast<int>(worker_count),
        [&](size_t worker_id, size_t region_index) {
            vector<uint32_t> buffer;
            buffer.reserve(RDATA_FLUSH_SIZE);
            scan_region_kmers(worker_fai[worker_id].get(), regions[region_index],
                              opt.chunk_size,
                              [&](hts_pos_t, uint64_t kmer) {
                const auto found = kmer_index.find(kmer);
                if (found == kmer_index.end())
                    throw std::runtime_error("internal error: pass-2 k-mer absent from index");
                const uint64_t packed = static_cast<uint64_t>(found->second) * num_regions + region_index;
                buffer.push_back(static_cast<uint32_t>(packed));
                if (buffer.size() == RDATA_FLUSH_SIZE) builder.append(buffer);
            });
            builder.append(buffer);
        });

    log("rdata occurrences: " + std::to_string(builder.data.size()) + "; sorting by k-mer index");
    std::sort(builder.data.begin(), builder.data.end());
    return std::move(builder.data);
}

vector<uint64_t> build_rdata_offsets(
    const vector<uint32_t> &rdata,
    size_t num_kmers,
    size_t num_regions)
{
    vector<uint64_t> offsets(num_kmers + 1, 0);
    size_t pos = 0;
    for (size_t kmer_index = 0; kmer_index < num_kmers; ++kmer_index) {
        offsets[kmer_index] = pos;
        const uint64_t upper = (static_cast<uint64_t>(kmer_index) + 1) * num_regions;
        while (pos < rdata.size() && static_cast<uint64_t>(rdata[pos]) < upper) ++pos;
    }
    offsets[num_kmers] = rdata.size();
    if (pos != rdata.size()) throw std::runtime_error("internal error while indexing rdata");
    return offsets;
}

vector<Hotspot> find_hotspots(
    faidx_t *fai,
    const Region &query_region,
    const Options &opt,
    const std::unordered_map<uint64_t, uint32_t> &kmer_index)
{
    vector<Hotspot> hotspots;
    std::deque<hts_pos_t> last_hits;

    scan_region_kmers(fai, query_region, opt.chunk_size,
                      [&](hts_pos_t kmer_start, uint64_t kmer) {
        if (kmer_index.find(kmer) == kmer_index.end()) return;

        last_hits.push_back(kmer_start);
        if (last_hits.size() > opt.hotspot_kmers) last_hits.pop_front();
        if (last_hits.size() < opt.hotspot_kmers) return;

        const hts_pos_t oldest = last_hits.front();
        const uint64_t span = static_cast<uint64_t>(kmer_start - oldest);
        if (span >= opt.hotspot_window) return;

        const hts_pos_t seed_start = oldest;
        const hts_pos_t seed_end = std::min<hts_pos_t>(
            query_region.end, kmer_start + KMER_SIZE);

        if (!hotspots.empty() &&
            static_cast<uint64_t>(seed_start - hotspots.back().end > 0 ?
                                  seed_start - hotspots.back().end : 0) <= opt.merge_distance) {
            hotspots.back().end = std::max(hotspots.back().end, seed_end);
        } else {
            Hotspot hotspot;
            hotspot.start = seed_start;
            hotspot.end = seed_end;
            hotspots.push_back(std::move(hotspot));
        }
    });

    return hotspots;
}

void score_hotspots(
    faidx_t *fai,
    const Region &query_region,
    const Options &opt,
    const std::unordered_map<uint64_t, uint32_t> &kmer_index,
    const vector<uint32_t> &rdata,
    const vector<uint64_t> &rdata_offsets,
    size_t num_regions,
    vector<Hotspot> &hotspots)
{
    if (hotspots.empty()) return;
    for (auto &hotspot : hotspots) hotspot.scores.assign(num_regions, 0);

    size_t hotspot_index = 0;
    scan_region_kmers(fai, query_region, opt.chunk_size,
                      [&](hts_pos_t kmer_start, uint64_t kmer) {
        while (hotspot_index < hotspots.size() && kmer_start >= hotspots[hotspot_index].end)
            ++hotspot_index;
        if (hotspot_index >= hotspots.size()) return;
        Hotspot &hotspot = hotspots[hotspot_index];
        if (kmer_start < hotspot.start) return;

        const auto found = kmer_index.find(kmer);
        if (found == kmer_index.end()) return;
        const uint32_t kidx = found->second;
        const uint64_t begin = rdata_offsets[kidx];
        const uint64_t end = rdata_offsets[static_cast<size_t>(kidx) + 1];
        for (uint64_t i = begin; i < end; ++i) {
            const uint32_t region_index = rdata[static_cast<size_t>(i)] % static_cast<uint32_t>(num_regions);
            ++hotspot.scores[region_index];
        }
    });
}

string bed_name_for_hotspot(
    const Hotspot &hotspot,
    const vector<Region> &reference_regions,
    uint32_t &best_region)
{
    struct Candidate {
        uint32_t region = 0;
        uint64_t raw_score = 0;
        double similarity = 0.0;
    };

    vector<Candidate> candidates;
    candidates.reserve(reference_regions.size());
    double best_similarity = -1.0;
    best_region = 0;
    uint64_t best_raw = 0;

    for (size_t i = 0; i < reference_regions.size(); ++i) {
        const uint64_t raw = hotspot.scores[i];
        const uint64_t denom = reference_regions[i].length();
        const double similarity = denom == 0 ? 0.0 : static_cast<double>(raw) / static_cast<double>(denom);
        if (similarity > best_similarity ||
            (similarity == best_similarity && raw > best_raw)) {
            best_similarity = similarity;
            best_raw = raw;
            best_region = static_cast<uint32_t>(i);
        }
        candidates.push_back(Candidate{static_cast<uint32_t>(i), raw, similarity});
    }

    if (best_similarity <= 0.0) return {};
    const double threshold = 0.5 * best_similarity;
    candidates.erase(
        std::remove_if(candidates.begin(), candidates.end(),
                       [&](const Candidate &c) { return !(c.similarity > threshold); }),
        candidates.end());

    std::sort(candidates.begin(), candidates.end(), [](const Candidate &a, const Candidate &b) {
        if (a.similarity != b.similarity) return a.similarity > b.similarity;
        if (a.raw_score != b.raw_score) return a.raw_score > b.raw_score;
        return a.region < b.region;
    });

    std::ostringstream name;
    name << std::fixed << std::setprecision(6);
    for (size_t j = 0; j < candidates.size(); ++j) {
        if (j) name << ',';
        const Candidate &c = candidates[j];
        name << reference_regions[c.region].name
             << "|similarity=" << c.similarity
             << "|score=" << c.raw_score;
    }
    return name.str();
}

vector<BedOutput> process_query_region(
    faidx_t *fai,
    const Region &query_region,
    const Options &opt,
    const vector<Region> &reference_regions,
    const std::unordered_map<uint64_t, uint32_t> &kmer_index,
    const vector<uint32_t> &rdata,
    const vector<uint64_t> &rdata_offsets)
{
    vector<Hotspot> hotspots = find_hotspots(fai, query_region, opt, kmer_index);
    if (hotspots.empty()) return {};

    score_hotspots(fai, query_region, opt, kmer_index, rdata, rdata_offsets,
                   reference_regions.size(), hotspots);

    vector<BedOutput> output;
    output.reserve(hotspots.size());
    for (const Hotspot &hotspot : hotspots) {
        uint32_t best_region = 0;
        string name = bed_name_for_hotspot(hotspot, reference_regions, best_region);
        if (name.empty()) continue;
        output.push_back(BedOutput{query_region.chrom, hotspot.start, hotspot.end,
                                   std::move(name), best_region});
    }
    return output;
}

void write_bed6(const string &path, const vector<vector<BedOutput>> &results) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open output: " + path);
    for (const auto &task_results : results) {
        for (const auto &record : task_results) {
            out << record.chrom << '\t'
                << record.start << '\t'
                << record.end << '\t'
                << record.name << '\t'
                << record.best_region << '\t'
                << ".\n";
        }
    }
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options opt = parse_args(argc, argv);
        if (!std::filesystem::is_regular_file(opt.reference_fasta))
            throw std::runtime_error("reference FASTA not found: " + opt.reference_fasta);
        if (!std::filesystem::is_regular_file(opt.reference_bed))
            throw std::runtime_error("reference BED not found: " + opt.reference_bed);
        if (!std::filesystem::is_regular_file(opt.query_fasta))
            throw std::runtime_error("query FASTA not found: " + opt.query_fasta);
        if (!opt.query_bed.empty() && !std::filesystem::is_regular_file(opt.query_bed))
            throw std::runtime_error("query BED not found: " + opt.query_bed);

        FaidxPtr reference_fai = load_faidx(opt.reference_fasta, false);
        vector<Region> reference_regions = read_bed_regions(
            opt.reference_bed, reference_fai.get(), true, "region_");
        log("reference BED regions: " + std::to_string(reference_regions.size()));

        auto kmer_index = build_kmer_index(opt, reference_regions);
        vector<uint32_t> rdata = build_rdata(opt, reference_regions, kmer_index);
        vector<uint64_t> rdata_offsets = build_rdata_offsets(
            rdata, kmer_index.size(), reference_regions.size());

        FaidxPtr query_metadata = load_faidx(opt.query_fasta, true);
        vector<Region> query_regions;
        if (opt.query_bed.empty()) {
            query_regions = whole_contig_regions(query_metadata.get());
        } else {
            query_regions = read_bed_regions(
                opt.query_bed, query_metadata.get(), false, "query_region_");
        }
        log("query scan targets: " + std::to_string(query_regions.size()));

        const size_t worker_count = std::min<size_t>(
            query_regions.size(), static_cast<size_t>(std::max(1, opt.threads)));
        vector<FaidxPtr> query_fai;
        query_fai.reserve(worker_count);
        for (size_t i = 0; i < worker_count; ++i)
            query_fai.push_back(load_faidx(opt.query_fasta, false));

        vector<vector<BedOutput>> results(query_regions.size());
        log("scanning query with " + std::to_string(worker_count) +
            " threads; hotspot seed=" + std::to_string(opt.hotspot_kmers) +
            " matching k-mers/" + std::to_string(opt.hotspot_window) +
            " bp; merge distance=" + std::to_string(opt.merge_distance) + " bp");

        parallel_for_regions(query_regions.size(), static_cast<int>(worker_count),
            [&](size_t worker_id, size_t query_index) {
                results[query_index] = process_query_region(
                    query_fai[worker_id].get(), query_regions[query_index], opt,
                    reference_regions, kmer_index, rdata, rdata_offsets);
            });

        write_bed6(opt.output_bed, results);

        uint64_t hotspot_count = 0;
        for (const auto &entries : results) hotspot_count += entries.size();
        log("wrote " + std::to_string(hotspot_count) + " hotspots to " + opt.output_bed);
        log("BED score column is the zero-based index of the best reference BED region");
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
