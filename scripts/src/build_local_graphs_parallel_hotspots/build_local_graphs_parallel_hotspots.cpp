// Build per-partition hotspot FASTAs with shared-memory C++ threads.
// Compile: g++ -O3 -std=c++17 -pthread -o build_local_graphs_parallel_hotspots build_local_graphs_parallel_hotspots.cpp

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace fs = std::filesystem;
using Coordinate = uint32_t;
using Interval = std::pair<Coordinate, Coordinate>;
using ContigIntervals = std::unordered_map<std::string, std::vector<Interval>>;

struct SourceContig {
    std::string sample;
    std::string contig;

    bool operator==(const SourceContig& other) const {
        return sample == other.sample && contig == other.contig;
    }
};

struct SourceContigHash {
    size_t operator()(const SourceContig& value) const {
        size_t seed = std::hash<std::string>{}(value.sample);
        seed ^= std::hash<std::string>{}(value.contig) +
                0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        return seed;
    }
};

using PartitionIntervals =
    std::unordered_map<SourceContig, std::vector<Interval>, SourceContigHash>;
using HotspotList = std::vector<PartitionIntervals>;

struct Options {
    std::string hotspot_dir;
    std::string query_paths;
    std::string target_list;
    std::string output_dir;
    std::string materialize_header;
    std::string output_fasta;
    unsigned jobs = 16;
    Coordinate kmer_length = 31;
    Coordinate merge_distance = 30000;
    Coordinate anchor = 15000;
    bool resume = false;
};

struct FaiRecord {
    Coordinate length;
    uint64_t offset;
    Coordinate line_bases;
    Coordinate line_width;
};

struct AssemblyIndex {
    size_t order;
    std::string name;
    std::string fasta;
    std::string fai;
    std::unordered_map<std::string, FaiRecord> records;
};

std::mutex log_mutex;

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void log(const std::string& message) {
    std::lock_guard<std::mutex> lock(log_mutex);
    std::cerr << message << '\n';
}

std::vector<std::string> split_ws(const std::string& line) {
    std::istringstream input(line);
    std::vector<std::string> fields;
    std::string field;
    while (input >> field) fields.push_back(field);
    return fields;
}

uint64_t parse_u64(const std::string& text, const std::string& context) {
    size_t used = 0;
    unsigned long long value = 0;
    try {
        value = std::stoull(text, &used);
    } catch (...) {
        fail(context + ": invalid unsigned integer " + text);
    }
    if (used != text.size()) fail(context + ": invalid unsigned integer " + text);
    return static_cast<uint64_t>(value);
}

Coordinate coordinate(const std::string& text, const std::string& context) {
    uint64_t value = parse_u64(text, context);
    if (value > std::numeric_limits<Coordinate>::max()) {
        fail(context + ": coordinate exceeds uint32: " + text);
    }
    return static_cast<Coordinate>(value);
}

bool ends_with(const std::string& value, const std::string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string partition_name_from_path(std::string path) {
    while (!path.empty() && path.back() == '/') path.pop_back();
    std::string filename = fs::path(path).filename().string();
    std::string lowered = filename;
    std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    static const std::vector<std::string> suffixes = {
        ".fasta_kmer.list", ".fna_kmer.list", ".fa_kmer.list",
        "_kmer.list", ".fasta", ".fna", ".fa", ".list"
    };
    std::string name = filename;
    for (const auto& suffix : suffixes) {
        if (ends_with(lowered, suffix)) {
            name.resize(filename.size() - suffix.size());
            break;
        }
    }
    std::string parent = fs::path(path).parent_path().filename().string();
    return !parent.empty() && parent == name ? parent : name;
}

std::vector<std::string> read_targets(const std::string& path) {
    std::ifstream input(path);
    if (!input) fail("cannot open target list: " + path);
    std::vector<std::string> targets(1, "");  // 1-based names.
    std::unordered_map<std::string, bool> seen;
    std::string line;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        std::istringstream row(line);
        std::string target_path, extra;
        if (!(row >> target_path) || (row >> extra)) {
            fail(path + ":" + std::to_string(line_number) +
                 ": expected exactly one target path");
        }
        std::string name = partition_name_from_path(target_path);
        if (name.empty() || name == "." || name == ".." || seen.count(name)) {
            fail(path + ":" + std::to_string(line_number) +
                 ": empty/duplicate partition name " + name);
        }
        seen[name] = true;
        targets.push_back(std::move(name));
    }
    if (targets.size() == 1) fail(path + ": no target paths");
    return targets;
}

std::vector<std::string> discover_hotspots(const std::string& directory) {
    std::vector<std::string> paths;
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.is_regular_file() && ends_with(entry.path().filename().string(), "_hotspot.txt")) {
            paths.push_back(entry.path().string());
        }
    }
    std::sort(paths.begin(), paths.end());
    if (paths.empty()) fail(directory + ": no *_hotspot.txt files");
    return paths;
}

template <typename Function>
void parallel_indices(size_t count, unsigned jobs, Function function) {
    if (count == 0) return;
    std::atomic<size_t> next{0};
    std::mutex error_mutex;
    std::exception_ptr error;
    const unsigned workers = std::min<unsigned>(jobs, static_cast<unsigned>(count));
    std::vector<std::thread> threads;
    threads.reserve(workers);
    for (unsigned worker = 0; worker < workers; ++worker) {
        threads.emplace_back([&]() {
            try {
                while (true) {
                    size_t index = next.fetch_add(1);
                    if (index >= count) break;
                    function(index);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!error) error = std::current_exception();
                next.store(count);
            }
        });
    }
    for (auto& thread : threads) thread.join();
    if (error) std::rethrow_exception(error);
}

struct HotspotStats {
    std::atomic<uint64_t> rows{0};
    std::atomic<uint64_t> mitochondrial{0};
};

std::string infer_sample(const std::string& contig);

void load_hotspots(
    const std::vector<std::string>& files,
    HotspotList& hotspots,
    std::vector<std::mutex>& partition_mutexes,
    const std::unordered_map<std::string, size_t>& sample_ids,
    Coordinate kmer_length,
    unsigned jobs,
    HotspotStats& stats
) {
    parallel_indices(files.size(), jobs, [&](size_t file_index) {
        const std::string& path = files[file_index];
        const std::string filename = fs::path(path).filename().string();
        const std::string suffix = "_hotspot.txt";
        std::string source_sample = filename.substr(0, filename.size() - suffix.size());
        if (!sample_ids.count(source_sample)) source_sample.clear();
        std::ifstream input(path);
        if (!input) fail("cannot open hotspot file: " + path);
        std::string line;
        size_t line_number = 0;
        while (std::getline(input, line)) {
            ++line_number;
            if (line.empty() || line[0] == '#') continue;
            auto fields = split_ws(line);
            std::string context = path + ":" + std::to_string(line_number);
            if (fields.size() < 4) fail(context + ": expected CONTIG TARGET START END");
            const std::string& contig = fields[0];
            uint64_t partition64 = parse_u64(fields[1], context);
            Coordinate raw_start = coordinate(fields[2], context);
            Coordinate raw_end = coordinate(fields[3], context);
            if (partition64 == 0 || partition64 >= hotspots.size() ||
                raw_start == 0 || raw_end < raw_start) {
                fail(context + ": invalid hotspot row");
            }
            if (contig == "chrM") {
                ++stats.mitochondrial;
                continue;
            }
            Coordinate start = raw_start > kmer_length ? raw_start - kmer_length : 0;
            Coordinate end = raw_end;
            if (start < end) {
                size_t partition = static_cast<size_t>(partition64);
                std::string sample = source_sample.empty() ?
                    infer_sample(contig) : source_sample;
                std::lock_guard<std::mutex> lock(partition_mutexes[partition]);
                hotspots[partition][SourceContig{sample, contig}].emplace_back(start, end);
                ++stats.rows;
            }
        }
    });
}

AssemblyIndex load_fai(
    size_t order, const std::string& name, const std::string& fasta,
    const std::string& fai
) {
    if (ends_with(fasta, ".gz") || ends_with(fasta, ".bgz") || ends_with(fasta, ".bgzf")) {
        fail("FAI extraction requires uncompressed FASTA: " + fasta);
    }
    std::ifstream input(fai);
    if (!input) fail("cannot open FAI: " + fai);
    AssemblyIndex index{order, name, fs::absolute(fasta).string(), fs::absolute(fai).string(), {}};
    std::string line;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty()) continue;
        auto fields = split_ws(line);
        std::string context = fai + ":" + std::to_string(line_number);
        if (fields.size() < 5) fail(context + ": expected FAI5+");
        FaiRecord record{
            coordinate(fields[1], context), parse_u64(fields[2], context),
            coordinate(fields[3], context), coordinate(fields[4], context)
        };
        if (!record.line_bases || record.line_width < record.line_bases ||
            !index.records.emplace(fields[0], record).second) {
            fail(context + ": invalid/duplicate FAI record");
        }
    }
    if (index.records.empty()) fail(fai + ": empty FAI");
    return index;
}

std::vector<AssemblyIndex> read_query_indexes(const std::string& path, unsigned jobs) {
    std::ifstream input(path);
    if (!input) fail("cannot open query paths: " + path);
    struct Row { std::string name, fasta, fai; };
    std::vector<Row> rows;
    std::unordered_map<std::string, bool> seen;
    std::string line;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') continue;
        auto fields = split_ws(line);
        if (fields.size() < 2) fail(path + ":" + std::to_string(line_number) + ": expected NAME FASTA [FAI]");
        if (seen.count(fields[0])) fail(path + ":" + std::to_string(line_number) + ": duplicate sample " + fields[0]);
        seen[fields[0]] = true;
        rows.push_back({fields[0], fields[1], fields.size() > 2 ? fields[2] : fields[1] + ".fai"});
    }
    if (rows.empty()) fail(path + ": no assemblies");
    std::vector<AssemblyIndex> indexes(rows.size());
    parallel_indices(rows.size(), jobs, [&](size_t i) {
        indexes[i] = load_fai(i, rows[i].name, rows[i].fasta, rows[i].fai);
    });
    return indexes;
}

std::string infer_sample(const std::string& contig) {
    size_t first = contig.find('#');
    if (first != std::string::npos) {
        size_t second = contig.find('#', first + 1);
        std::string sample = contig.substr(0, first);
        std::string hap = contig.substr(first + 1, second == std::string::npos ?
                                        std::string::npos : second - first - 1);
        if (sample.empty() || hap.empty()) fail("cannot infer sample from contig " + contig);
        if (hap[0] == 'h') hap.erase(0, 1);
        return sample + "_h" + hap;
    }
    if (contig.rfind("NC_0609", 0) == 0) return "CHM13_h1";
    return "HG38_h1";
}

std::vector<Interval> merge_anchor(
    std::vector<Interval> intervals, Coordinate length,
    Coordinate merge_distance, Coordinate anchor
) {
    std::sort(intervals.begin(), intervals.end());
    std::vector<Interval> merged;
    for (auto interval : intervals) {
        interval.second = std::min(interval.second, length);
        if (interval.first >= interval.second) continue;
        bool separated = !merged.empty() && interval.first > merged.back().second &&
                         interval.first - merged.back().second >= merge_distance;
        if (merged.empty() || separated) merged.push_back(interval);
        else merged.back().second = std::max(merged.back().second, interval.second);
    }
    for (auto& interval : merged) {
        interval.first = interval.first > anchor ? interval.first - anchor : 0;
        uint64_t extended = static_cast<uint64_t>(interval.second) + anchor;
        interval.second = static_cast<Coordinate>(std::min<uint64_t>(length, extended));
    }
    return merged;
}

std::string fetch_sequence(int descriptor, const FaiRecord& record,
                           Coordinate start, Coordinate end,
                           const std::string& context) {
    if (start >= end || end > record.length) fail(context + ": invalid interval");
    std::string sequence(end - start, '\0');
    Coordinate position = start;
    size_t destination = 0;
    while (position < end) {
        Coordinate line = position / record.line_bases;
        Coordinate column = position % record.line_bases;
        Coordinate size = std::min<Coordinate>(end - position, record.line_bases - column);
        uint64_t file_offset = record.offset + static_cast<uint64_t>(line) * record.line_width + column;
        ssize_t got = ::pread(descriptor, sequence.data() + destination, size, static_cast<off_t>(file_offset));
        if (got != static_cast<ssize_t>(size)) {
            fail(context + ": short FASTA read: " + std::strerror(errno));
        }
        position += size;
        destination += size;
    }
    return sequence;
}

constexpr size_t PARTITION_BUFFER_LIMIT = size_t{1} << 30;  // 1 GiB.
constexpr size_t GLOBAL_BUFFER_LIMIT = size_t{2} << 30;     // 2 GiB.
constexpr size_t PARTITION_BUFFER_BLOCK = size_t{16} << 20; // 16 MiB.
constexpr Coordinate SEQUENCE_READ_CHUNK = Coordinate{16} << 20;
std::atomic<size_t> global_buffered_bytes{0};

class AtomicBufferedWriter {
public:
    explicit AtomicBufferedWriter(const fs::path& final_path)
        : final_path_(final_path),
          temporary_path_(
              final_path.string() + ".tmp." + std::to_string(::getpid()) + "." +
              std::to_string(std::hash<std::thread::id>{}(std::this_thread::get_id()))),
          output_(temporary_path_, std::ios::binary) {
        if (!output_) fail("cannot create " + temporary_path_);
    }

    AtomicBufferedWriter(const AtomicBufferedWriter&) = delete;
    AtomicBufferedWriter& operator=(const AtomicBufferedWriter&) = delete;

    ~AtomicBufferedWriter() {
        if (!committed_) {
            output_.close();
            discard_buffer();
            std::error_code ignored;
            fs::remove(temporary_path_, ignored);
        }
    }

    void append(const char* data, size_t size) {
        while (size) {
            if (blocks_.empty() || blocks_.back().size() == PARTITION_BUFFER_BLOCK) {
                blocks_.emplace_back();
                blocks_.back().reserve(PARTITION_BUFFER_BLOCK);
            }
            std::string& block = blocks_.back();
            size_t amount = std::min(size, PARTITION_BUFFER_BLOCK - block.size());
            block.append(data, amount);
            data += amount;
            size -= amount;
            buffered_ += amount;
            size_t global_size =
                global_buffered_bytes.fetch_add(amount, std::memory_order_relaxed) + amount;
            if (buffered_ == PARTITION_BUFFER_LIMIT ||
                global_size >= GLOBAL_BUFFER_LIMIT) {
                flush();
            }
        }
    }

    void append(const std::string& data) {
        append(data.data(), data.size());
    }

    void append(char value) {
        append(&value, 1);
    }

    void commit() {
        flush();
        output_.close();
        if (!output_) fail("failed closing " + temporary_path_);
        std::error_code error;
        fs::rename(temporary_path_, final_path_, error);
        if (error) fail("cannot replace " + final_path_.string() + ": " + error.message());
        committed_ = true;
    }

private:
    void flush() {
        if (!buffered_) return;
        for (const std::string& block : blocks_) {
            output_.write(block.data(), static_cast<std::streamsize>(block.size()));
        }
        if (!output_) fail("failed writing " + temporary_path_);
        discard_buffer();
    }

    void discard_buffer() noexcept {
        const size_t released = buffered_;
        std::vector<std::string>().swap(blocks_);
        buffered_ = 0;
        if (released) {
            global_buffered_bytes.fetch_sub(released, std::memory_order_relaxed);
        }
    }

    fs::path final_path_;
    std::string temporary_path_;
    std::ofstream output_;
    std::vector<std::string> blocks_;
    size_t buffered_ = 0;
    bool committed_ = false;
};

void append_wrapped_chunk(
    AtomicBufferedWriter& output,
    const std::string& sequence,
    size_t& line_column
) {
    constexpr size_t width = 80;
    size_t start = 0;
    while (start < sequence.size()) {
        size_t amount = std::min(width - line_column, sequence.size() - start);
        output.append(sequence.data() + start, amount);
        start += amount;
        line_column += amount;
        if (line_column == width) {
            output.append('\n');
            line_column = 0;
        }
    }
}

std::string record_prefix(const std::string& partition) {
    if (partition.size() > 2 && partition[0] == 'g') {
        size_t end = 1;
        while (end < partition.size() && std::isdigit(
                   static_cast<unsigned char>(partition[end]))) ++end;
        if (end > 1 && end < partition.size() && partition[end] == '_') {
            return partition.substr(0, end);
        }
    }
    return partition;
}

bool output_has_partition_prefix(
    const fs::path& output,
    const std::string& partition
) {
    if (!fs::is_regular_file(output) || fs::file_size(output) == 0) return false;
    std::ifstream input(output, std::ios::binary);
    std::string first_header;
    if (!std::getline(input, first_header)) return false;
    const std::string expected = ">" + record_prefix(partition) + "_";
    return first_header.compare(0, expected.size(), expected) == 0;
}

struct BuildStats {
    std::atomic<uint64_t> partitions{0}, records{0}, bases{0}, resumed{0};
};

struct HeaderRequest {
    std::string header;
    std::string record_id;
    std::string sample;
    std::string contig;
    Coordinate start;
    Coordinate end;
};

HeaderRequest parse_header_request(
    const std::string& line,
    size_t line_number,
    const std::string& path,
    const std::vector<AssemblyIndex>& assemblies
) {
    const std::string context = path + ":" + std::to_string(line_number);
    if (line.empty() || line[0] != '>') {
        fail(context + ": expected a FASTA header beginning with >");
    }
    const auto fields = split_ws(line);
    if (fields.size() < 2 || fields[0].size() < 2) {
        fail(context + ": header lacks a record ID or source coordinates");
    }
    const std::string record_id = fields[0].substr(1);
    std::string sample;
    std::string contig;
    Coordinate start = 0, end = 0;
    bool found_coordinate = false;
    for (size_t index = 1; index < fields.size(); ++index) {
        const std::string& field = fields[index];
        if (field.rfind("sample=", 0) == 0) {
            sample = field.substr(7);
            continue;
        }
        if (found_coordinate) continue;
        const size_t colon = field.rfind(':');
        if (colon == std::string::npos || colon == 0) continue;
        std::string interval = field.substr(colon + 1);
        if (!interval.empty() && (interval.back() == '+' || interval.back() == '-')) {
            interval.pop_back();
        }
        const size_t dash = interval.find('-');
        if (dash == std::string::npos || dash == 0 || dash + 1 >= interval.size()) continue;
        try {
            start = coordinate(interval.substr(0, dash), context);
            end = coordinate(interval.substr(dash + 1), context);
        } catch (const std::exception&) {
            continue;
        }
        if (end <= start) fail(context + ": invalid source interval " + field);
        contig = field.substr(0, colon);
        found_coordinate = true;
    }
    if (!found_coordinate) fail(context + ": no valid CONTIG:START-END source coordinate");
    if (sample.empty()) {
        for (const AssemblyIndex& assembly : assemblies) {
            const std::string token = "_" + assembly.name + "_";
            if (record_id.find(token) != std::string::npos) {
                if (!sample.empty()) {
                    fail(context + ": record ID matches more than one query sample; add sample=NAME");
                }
                sample = assembly.name;
            }
        }
    }
    if (sample.empty()) {
        fail(context + ": cannot determine source sample; add sample=NAME to the header");
    }
    return HeaderRequest{line, record_id, sample, contig, start, end};
}

void materialize_one_header(
    const Options& options,
    const std::vector<AssemblyIndex>& assemblies,
    const std::unordered_map<std::string, size_t>& sample_ids
) {
    std::ifstream input(options.materialize_header);
    if (!input) fail("cannot open header manifest: " + options.materialize_header);
    fs::path output_path(options.output_fasta);
    if (!output_path.parent_path().empty()) fs::create_directories(output_path.parent_path());
    AtomicBufferedWriter writer(output_path);
    uint64_t records = 0, bases = 0;
    std::vector<std::vector<HeaderRequest>> requests_by_assembly(assemblies.size());
    std::string line;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) fail(options.materialize_header + ":" + std::to_string(line_number) + ": empty header line");
        HeaderRequest request = parse_header_request(
            line, line_number, options.materialize_header, assemblies
        );
        auto sample_it = sample_ids.find(request.sample);
        if (sample_it == sample_ids.end()) {
            fail(options.materialize_header + ":" + std::to_string(line_number) +
                 ": sample is absent from query paths: " + request.sample);
        }
        const size_t assembly_index = sample_it->second;
        const AssemblyIndex& assembly = assemblies[assembly_index];
        auto record_it = assembly.records.find(request.contig);
        if (record_it == assembly.records.end()) {
            fail(options.materialize_header + ":" + std::to_string(line_number) +
                 ": contig " + request.contig + " is absent from " + assembly.fai);
        }
        if (request.end > record_it->second.length) {
            fail(options.materialize_header + ":" + std::to_string(line_number) +
                 ": source interval exceeds contig length");
        }
        requests_by_assembly[assembly_index].push_back(std::move(request));
        ++records;
    }
    if (!records) fail(options.materialize_header + ": no header records");

    for (size_t assembly_index = 0; assembly_index < assemblies.size(); ++assembly_index) {
        const auto& requests = requests_by_assembly[assembly_index];
        if (requests.empty()) continue;
        const AssemblyIndex& assembly = assemblies[assembly_index];
        int descriptor = ::open(assembly.fasta.c_str(), O_RDONLY);
        if (descriptor < 0) {
            fail("cannot open FASTA " + assembly.fasta + ": " + std::strerror(errno));
        }
        try {
            for (const HeaderRequest& request : requests) {
                const auto record_it = assembly.records.find(request.contig);
                if (record_it == assembly.records.end()) {
                    fail("internal error: validated contig disappeared: " + request.contig);
                }
                writer.append(request.header);
                writer.append('\n');
                size_t line_column = 0;
                Coordinate position = request.start;
                while (position < request.end) {
                    const Coordinate amount = std::min(
                        static_cast<Coordinate>(request.end - position),
                        SEQUENCE_READ_CHUNK
                    );
                    const Coordinate chunk_end = position + amount;
                    std::string sequence = fetch_sequence(
                        descriptor, record_it->second, position, chunk_end,
                        assembly.fasta + ":" + request.contig
                    );
                    append_wrapped_chunk(writer, sequence, line_column);
                    position = chunk_end;
                }
                if (line_column) writer.append('\n');
                bases += static_cast<uint64_t>(request.end - request.start);
            }
        } catch (...) {
            ::close(descriptor);
            throw;
        }
        ::close(descriptor);
    }
    writer.commit();
    log("Materialized " + std::to_string(records) + " records (" +
        std::to_string(bases) + " bp) from " + options.materialize_header +
        " to " + options.output_fasta);
}

void build_partition(
    size_t partition_index,
    const std::vector<std::string>& targets,
    PartitionIntervals& partition_hotspots,
    const std::vector<AssemblyIndex>& assemblies,
    const std::unordered_map<std::string, size_t>& sample_ids,
    const Options& options,
    BuildStats& stats
) {
    const std::string& partition = targets[partition_index];
    fs::path directory = fs::path(options.output_dir) / partition;
    fs::path output = directory / (partition + "_samples.fasta");
    if (options.resume && output_has_partition_prefix(output, partition)) {
        ++stats.resumed;
        return;
    }

    std::vector<ContigIntervals> by_sample(assemblies.size());
    for (auto& item : partition_hotspots) {
        const std::string& sample = item.first.sample;
        const std::string& contig = item.first.contig;
        auto sample_it = sample_ids.find(sample);
        if (sample_it == sample_ids.end()) fail(partition + ": inferred sample absent from query paths: " + sample);
        const AssemblyIndex& assembly = assemblies[sample_it->second];
        auto record = assembly.records.find(contig);
        if (record == assembly.records.end()) fail(partition + ": contig absent from " + assembly.fai + ": " + contig);
        auto intervals = merge_anchor(
            std::move(item.second), record->second.length,
            options.merge_distance, options.anchor
        );
        if (!intervals.empty()) by_sample[assembly.order].emplace(contig, std::move(intervals));
    }

    uint64_t made_records = 0, made_bases = 0;
    std::string prefix = record_prefix(partition);
    fs::create_directories(directory);
    AtomicBufferedWriter writer(output);
    for (const AssemblyIndex& assembly : assemblies) {
        auto& contigs = by_sample[assembly.order];
        if (contigs.empty()) continue;
        std::vector<std::string> names;
        names.reserve(contigs.size());
        for (const auto& item : contigs) names.push_back(item.first);
        std::sort(names.begin(), names.end());
        int descriptor = ::open(assembly.fasta.c_str(), O_RDONLY);
        if (descriptor < 0) fail("cannot open FASTA " + assembly.fasta + ": " + std::strerror(errno));
        try {
            uint64_t sequence_index = 0;
            for (const std::string& contig : names) {
                const FaiRecord& record = assembly.records.at(contig);
                for (const Interval& interval : contigs.at(contig)) {
                    ++sequence_index;
                    writer.append(">" + prefix + "_" + assembly.name + "_" +
                                  std::to_string(sequence_index) + "\t" + contig + ":" +
                                  std::to_string(interval.first) + "-" +
                                  std::to_string(interval.second) + "\tpartition=" + partition + "\n");
                    size_t line_column = 0;
                    Coordinate position = interval.first;
                    while (position < interval.second) {
                        Coordinate remaining = interval.second - position;
                        Coordinate amount = std::min(remaining, SEQUENCE_READ_CHUNK);
                        Coordinate chunk_end = position + amount;
                        std::string sequence = fetch_sequence(
                            descriptor, record, position, chunk_end,
                            assembly.fasta + ":" + contig
                        );
                        append_wrapped_chunk(writer, sequence, line_column);
                        position = chunk_end;
                    }
                    if (line_column) writer.append('\n');
                    ++made_records;
                    made_bases += static_cast<uint64_t>(interval.second - interval.first);
                }
            }
        } catch (...) {
            ::close(descriptor);
            throw;
        }
        ::close(descriptor);
    }
    if (!made_records) return;
    writer.commit();
    ++stats.partitions;
    stats.records += made_records;
    stats.bases += made_bases;
}

Options parse_args(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        std::string argument = argv[i];
        auto value = [&]() -> std::string {
            if (++i >= argc) fail("missing value after " + argument);
            return argv[i];
        };
        if (argument == "--hotspot-dir") options.hotspot_dir = value();
        else if (argument == "-q" || argument == "--query-paths") options.query_paths = value();
        else if (argument == "-T" || argument == "--target-list") options.target_list = value();
        else if (argument == "-o" || argument == "--output-dir") options.output_dir = value();
        else if (argument == "--materialize-header") options.materialize_header = value();
        else if (argument == "--output-fasta") options.output_fasta = value();
        else if (argument == "-j" || argument == "--jobs") options.jobs = static_cast<unsigned>(parse_u64(value(), argument));
        else if (argument == "--hotspot-kmer-length") options.kmer_length = coordinate(value(), argument);
        else if (argument == "--hotspot-merge-distance") options.merge_distance = coordinate(value(), argument);
        else if (argument == "--hotspot-anchor") options.anchor = coordinate(value(), argument);
        else if (argument == "--resume") options.resume = true;
        else if (argument == "-h" || argument == "--help") {
            std::cout
                << "Usage:\n  " << argv[0]
                << " --hotspot-dir DIR -q QUERY_PATHS -T TARGET_LIST -o DIR [-j N] [--resume]\n  "
                << argv[0]
                << " --materialize-header PARTITION.header -q QUERY_PATHS --output-fasta FILE\n";
            std::exit(0);
        } else fail("unknown argument: " + argument);
    }
    const bool materialize = !options.materialize_header.empty() || !options.output_fasta.empty();
    if (options.query_paths.empty()) fail("-q/--query-paths is required");
    if (materialize) {
        if (options.materialize_header.empty() || options.output_fasta.empty()) {
            fail("--materialize-header and --output-fasta must be used together");
        }
        if (!options.hotspot_dir.empty() || !options.target_list.empty() || !options.output_dir.empty()) {
            fail("header materialization cannot be combined with hotspot build arguments");
        }
    } else if (options.hotspot_dir.empty() || options.target_list.empty() || options.output_dir.empty()) {
        fail("missing required hotspot-build argument; use --help");
    }
    if (!options.jobs || !options.kmer_length) fail("jobs and kmer length must be positive");
    return options;
}

int main(int argc, char** argv) {
    try {
        Options options = parse_args(argc, argv);
        options.query_paths = fs::absolute(options.query_paths).string();
        if (!options.materialize_header.empty()) {
            options.materialize_header = fs::absolute(options.materialize_header).string();
            options.output_fasta = fs::absolute(options.output_fasta).string();
            auto assemblies = read_query_indexes(options.query_paths, options.jobs);
            std::unordered_map<std::string, size_t> sample_ids;
            for (size_t i = 0; i < assemblies.size(); ++i) {
                sample_ids.emplace(assemblies[i].name, i);
            }
            materialize_one_header(options, assemblies, sample_ids);
            return 0;
        }
        options.hotspot_dir = fs::absolute(options.hotspot_dir).string();
        options.target_list = fs::absolute(options.target_list).string();
        options.output_dir = fs::absolute(options.output_dir).string();
        fs::create_directories(options.output_dir);

        std::vector<std::string> targets = read_targets(options.target_list);
        const size_t target_count = targets.size() - 1;

        log("Loading all query FAI files into shared RAM");
        auto assemblies = read_query_indexes(options.query_paths, options.jobs);
        std::unordered_map<std::string, size_t> sample_ids;
        uint64_t contigs = 0;
        for (size_t i = 0; i < assemblies.size(); ++i) {
            sample_ids.emplace(assemblies[i].name, i);
            contigs += assemblies[i].records.size();
        }
        log("Loaded " + std::to_string(assemblies.size()) + " FAI catalogs with " +
            std::to_string(contigs) + " contigs");

        HotspotList hotspots(2 * target_count);
        std::vector<std::mutex> partition_mutexes(hotspots.size());
        auto hotspot_files = discover_hotspots(options.hotspot_dir);
        log("Loading " + std::to_string(hotspot_files.size()) + " hotspot files with " +
            std::to_string(options.jobs) + " shared-memory threads into " +
            std::to_string(hotspots.size()) + " partition slots");
        HotspotStats hotspot_stats;
        load_hotspots(hotspot_files, hotspots, partition_mutexes,
                      sample_ids, options.kmer_length, options.jobs, hotspot_stats);
        log("Loaded " + std::to_string(hotspot_stats.rows.load()) +
            " hotspots; excluded " + std::to_string(hotspot_stats.mitochondrial.load()) + " chrM rows");

        size_t unnamed = 0;
        for (size_t i = target_count + 1; i < hotspots.size(); ++i) {
            unnamed += !hotspots[i].empty();
            PartitionIntervals().swap(hotspots[i]);
        }
        if (unnamed) log("WARNING: " + std::to_string(unnamed) +
                         " valid partition IDs exceed the target-list row count and cannot be named/materialized");

        BuildStats build_stats;
        std::vector<size_t> active;
        size_t total_active = 0;
        for (size_t i = 1; i <= target_count; ++i) {
            if (hotspots[i].empty()) continue;
            ++total_active;
            const fs::path partition_output = fs::path(options.output_dir) /
                targets[i] / (targets[i] + "_samples.fasta");
            if (options.resume && output_has_partition_prefix(
                    partition_output, targets[i])) {
                ++build_stats.resumed;
                PartitionIntervals().swap(hotspots[i]);
            } else {
                if (options.resume && fs::is_regular_file(partition_output) &&
                    fs::file_size(partition_output) > 0) {
                    log("Rebuilding " + targets[i] +
                        ": existing sample FASTA has a legacy/mismatched prefix");
                }
                active.push_back(i);
            }
        }
        log("Resume scan released hotspot data for " +
            std::to_string(build_stats.resumed.load()) + " completed partitions; " +
            std::to_string(active.size()) + " partitions remain");
        std::atomic<size_t> completed_partitions{build_stats.resumed.load()};
        const size_t progress_step = std::max<size_t>(1, total_active / 100);
        log("Starting partition write phase: " + std::to_string(total_active) +
            " active partitions (" + std::to_string(active.size()) + " unfinished) with " +
            std::to_string(options.jobs) + " threads");
        parallel_indices(active.size(), options.jobs, [&](size_t i) {
            const size_t partition_index = active[i];
            build_partition(
                partition_index, targets, hotspots[partition_index], assemblies,
                sample_ids, options, build_stats
            );
            PartitionIntervals().swap(hotspots[partition_index]);
            size_t completed = completed_partitions.fetch_add(1) + 1;
            if (completed == total_active || completed % progress_step == 0) {
                double percent = total_active == 0 ? 100.0 :
                    100.0 * static_cast<double>(completed) /
                    static_cast<double>(total_active);
                std::ostringstream message;
                message.setf(std::ios::fixed);
                message.precision(1);
                message << "Partition progress: " << completed << "/"
                        << total_active << " (" << percent << "%); written="
                        << build_stats.partitions.load() << ", resumed="
                        << build_stats.resumed.load() << ", sequences="
                        << build_stats.records.load() << ", bases="
                        << build_stats.bases.load();
                log(message.str());
            }
        });
        log("Complete: " + std::to_string(build_stats.partitions.load()) + " partitions, " +
            std::to_string(build_stats.records.load()) + " sequences, " +
            std::to_string(build_stats.bases.load()) + " bp, " +
            std::to_string(build_stats.resumed.load()) + " resumed under " + options.output_dir);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
