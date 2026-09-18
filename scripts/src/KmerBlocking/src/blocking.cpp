#include "blocking.hpp"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstring>
#include <deque>
#include <exception>
#include <fstream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <unordered_set>
#include <utility>

namespace kmer_blocking {
namespace {

struct Contig {
    std::string name;
    std::uint64_t length = 0;
    std::uint64_t global_offset = 0;
    std::vector<std::uint64_t> kmers;
};

struct Block {
    std::uint64_t start = 0;
    std::uint64_t end = 0;
};

struct Interval {
    std::uint64_t start = 0;
    std::uint64_t end = 0;
};

struct Extensions {
    std::uint64_t start = 0;
    std::uint64_t end = 0;
    std::uint64_t hotspots = 0;
};

struct ContigResult {
    std::vector<Block> blocks;
    std::uint64_t initial_blocks = 0;
    std::uint64_t hotspots = 0;
    std::uint64_t moved_breakpoints = 0;
    std::uint64_t absorbed_blocks = 0;
};

std::uint64_t parse_u64(const std::string& text, const std::string& label,
                        std::size_t line_number) {
    if (text.empty() || text[0] == '-') {
        throw std::runtime_error("invalid " + label + " on .fai line " +
                                 std::to_string(line_number) + ": " + text);
    }
    std::size_t used = 0;
    unsigned long long value = 0;
    try {
        value = std::stoull(text, &used, 10);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid " + label + " on .fai line " +
                                 std::to_string(line_number) + ": " + text);
    }
    if (used != text.size()) {
        throw std::runtime_error("invalid " + label + " on .fai line " +
                                 std::to_string(line_number) + ": " + text);
    }
    return static_cast<std::uint64_t>(value);
}

std::vector<Contig> read_fai(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open FASTA index " + path);
    }

    std::vector<Contig> contigs;
    std::string line;
    std::uint64_t global_offset = 0;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        const std::size_t first_tab = line.find('\t');
        const std::size_t second_tab =
            first_tab == std::string::npos
                ? std::string::npos
                : line.find('\t', first_tab + 1);
        if (first_tab == std::string::npos || first_tab == 0) {
            throw std::runtime_error(
                "expected tab-separated contig name and length on .fai line " +
                std::to_string(line_number));
        }
        const std::string name = line.substr(0, first_tab);
        const std::string length_text = line.substr(
            first_tab + 1,
            second_tab == std::string::npos ? std::string::npos
                                             : second_tab - first_tab - 1);
        const std::uint64_t length =
            parse_u64(length_text, "contig length", line_number);
        if (length > std::numeric_limits<std::uint32_t>::max() - global_offset) {
            throw std::runtime_error(
                "total assembly length must be less than 2^32 bases");
        }
        Contig contig;
        contig.name = name;
        contig.length = length;
        contig.global_offset = global_offset;
        if (length > std::numeric_limits<std::size_t>::max()) {
            throw std::runtime_error("contig is too long for this system: " + name);
        }
        contig.kmers.resize(static_cast<std::size_t>(length), 0);
        contigs.push_back(std::move(contig));
        global_offset += length;
    }
    if (contigs.empty()) {
        throw std::runtime_error("FASTA index is empty: " + path);
    }
    return contigs;
}

class BinaryReader {
public:
    explicit BinaryReader(const std::string& path)
        : path_(path), stream_(path, std::ios::binary) {
        if (!stream_) {
            throw std::runtime_error("cannot open k-mer index " + path + ": " +
                                     std::strerror(errno));
        }
    }

    bool next_u64(std::uint64_t& value) {
        char bytes[8];
        stream_.read(bytes, 8);
        const std::streamsize count = stream_.gcount();
        if (count == 0 && stream_.eof()) {
            return false;
        }
        if (count != 8) {
            throw std::runtime_error("truncated uint64 in k-mer index " + path_);
        }
        value = 0;
        for (int i = 0; i < 8; ++i) {
            value |= static_cast<std::uint64_t>(
                         static_cast<unsigned char>(bytes[i]))
                     << (8 * i);
        }
        return true;
    }

    std::uint32_t u32() {
        char bytes[4];
        stream_.read(bytes, 4);
        if (stream_.gcount() != 4) {
            throw std::runtime_error("truncated uint32 in k-mer index " + path_);
        }
        std::uint32_t value = 0;
        for (int i = 0; i < 4; ++i) {
            value |= static_cast<std::uint32_t>(
                         static_cast<unsigned char>(bytes[i]))
                     << (8 * i);
        }
        return value;
    }

private:
    std::string path_;
    std::ifstream stream_;
};

std::size_t contig_for_position(const std::vector<Contig>& contigs,
                                std::uint32_t position) {
    const auto it = std::upper_bound(
        contigs.begin(), contigs.end(), static_cast<std::uint64_t>(position),
        [](std::uint64_t value, const Contig& contig) {
            return value < contig.global_offset;
        });
    if (it == contigs.begin()) {
        throw std::runtime_error("index location precedes the first contig");
    }
    const std::size_t index = static_cast<std::size_t>(it - contigs.begin() - 1);
    if (position >= contigs[index].global_offset + contigs[index].length) {
        throw std::runtime_error("index location falls outside all contigs: " +
                                 std::to_string(position));
    }
    return index;
}

void store_location(std::vector<Contig>& contigs, std::uint64_t key,
                    std::uint32_t global_position) {
    const std::size_t contig_index =
        contig_for_position(contigs, global_position);
    Contig& contig = contigs[contig_index];
    const std::uint64_t local_position =
        global_position - contig.global_offset;
    if (local_position + kKmerLength > contig.length) {
        throw std::runtime_error(
            "index contains a 31-mer start beyond the usable part of contig " +
            contig.name + ": " + std::to_string(local_position));
    }
    std::uint64_t& slot = contig.kmers[static_cast<std::size_t>(local_position)];
    if (slot != 0) {
        throw std::runtime_error("index contains duplicate location " +
                                 std::to_string(global_position));
    }
    slot = key;
}

void load_index(const std::string& path, std::vector<Contig>& contigs,
                BlockingStats& stats) {
    std::uint64_t largest_start = 0;
    bool has_usable_bases = false;
    for (const Contig& contig : contigs) {
        if (contig.length >= kKmerLength) {
            largest_start = std::max(
                largest_start,
                contig.global_offset + contig.length - kKmerLength);
            has_usable_bases = true;
        }
    }

    BinaryReader reader(path);
    std::uint64_t key = 0;
    std::uint64_t previous_key = 0;
    bool have_previous_key = false;
    while (reader.next_u64(key)) {
        if (key == 0) {
            throw std::runtime_error("k-mer index contains reserved key zero");
        }
        if (have_previous_key && key <= previous_key) {
            throw std::runtime_error(
                "k-mer index keys are not in strictly increasing order");
        }
        previous_key = key;
        have_previous_key = true;
        ++stats.indexed_kmers;

        const std::uint32_t value = reader.u32();
        if (has_usable_bases && value <= largest_start) {
            store_location(contigs, key, value);
            ++stats.indexed_locations;
            continue;
        }

        const std::uint32_t count =
            std::numeric_limits<std::uint32_t>::max() - value;
        if (count < 2 || count > std::numeric_limits<std::uint16_t>::max()) {
            throw std::runtime_error(
                "invalid location/count marker " + std::to_string(value) +
                " for k-mer " + std::to_string(key));
        }
        for (std::uint32_t i = 0; i < count; ++i) {
            const std::uint32_t position = reader.u32();
            if (!has_usable_bases || position > largest_start) {
                throw std::runtime_error("invalid indexed location " +
                                         std::to_string(position));
            }
            store_location(contigs, key, position);
            ++stats.indexed_locations;
        }
    }
}

std::vector<Block> make_initial_blocks(std::uint64_t length,
                                       std::uint64_t block_size) {
    std::vector<Block> blocks;
    for (std::uint64_t start = 0; start < length;) {
        const std::uint64_t end =
            std::min(length, start + std::min(block_size, length - start));
        blocks.push_back({start, end});
        start = end;
    }
    if (blocks.size() > 1 &&
        blocks.back().end - blocks.back().start < kMinimumRemainder) {
        blocks[blocks.size() - 2].end = blocks.back().end;
        blocks.pop_back();
    }
    return blocks;
}

std::vector<Interval> find_hotspots(
    const std::vector<std::uint64_t>& container, std::uint64_t begin,
    std::uint64_t end, const std::unordered_set<std::uint64_t>& query_keys) {
    std::vector<Interval> hotspots;
    bool in_component = false;
    bool qualifies = false;
    std::uint64_t component_start = 0;
    std::uint64_t previous_hit = 0;
    std::deque<std::uint64_t> consecutive_hits;

    const auto finish_component = [&] {
        if (in_component && qualifies) {
            hotspots.push_back(
                {component_start,
                 std::min<std::uint64_t>(container.size(),
                                         previous_hit + kKmerLength)});
        }
    };

    for (std::uint64_t position = begin; position < end; ++position) {
        const std::uint64_t key = container[static_cast<std::size_t>(position)];
        if (key == 0 || query_keys.find(key) == query_keys.end()) {
            continue;
        }

        if (!in_component || position - previous_hit >= kConnectDistance) {
            finish_component();
            in_component = true;
            qualifies = false;
            component_start = position;
            consecutive_hits.clear();
        }
        previous_hit = position;
        consecutive_hits.push_back(position);
        if (consecutive_hits.size() > kHotspotHitCount) {
            consecutive_hits.pop_front();
        }
        if (consecutive_hits.size() == kHotspotHitCount &&
            consecutive_hits.back() - consecutive_hits.front() <=
                kHotspotSpan) {
            qualifies = true;
        }
    }
    finish_component();
    return hotspots;
}

Extensions find_extensions(const std::vector<std::uint64_t>& container,
                           const Block& block,
                           std::uint64_t search_distance) {
    Extensions extensions{block.start, block.end, 0};
    std::unordered_set<std::uint64_t> query_keys;
    query_keys.reserve(static_cast<std::size_t>(block.end - block.start));
    for (std::uint64_t position = block.start; position < block.end; ++position) {
        const std::uint64_t key = container[static_cast<std::size_t>(position)];
        if (key != 0) {
            query_keys.insert(key);
        }
    }
    if (query_keys.empty()) {
        return extensions;
    }

    const std::uint64_t left_begin =
        block.start > search_distance ? block.start - search_distance : 0;
    const std::vector<Interval> left =
        find_hotspots(container, left_begin, block.start, query_keys);
    for (const Interval& hotspot : left) {
        extensions.start = std::min(extensions.start, hotspot.start);
    }

    const std::uint64_t right_end =
        block.end > container.size() -
                        std::min<std::uint64_t>(search_distance, container.size())
            ? container.size()
            : block.end + search_distance;
    const std::vector<Interval> right =
        find_hotspots(container, block.end, right_end, query_keys);
    for (const Interval& hotspot : right) {
        extensions.end = std::max(extensions.end, hotspot.end);
    }
    extensions.hotspots = left.size() + right.size();
    return extensions;
}

ContigResult block_contig(const Contig& contig, const Options& options) {
    ContigResult result;
    result.blocks = make_initial_blocks(contig.length, options.block_size);
    result.initial_blocks = result.blocks.size();

    for (std::size_t index = 0; index < result.blocks.size(); ++index) {
        const Extensions extensions = find_extensions(
            contig.kmers, result.blocks[index], options.search_distance);
        result.hotspots += extensions.hotspots;

        if (extensions.start < result.blocks[index].start) {
            const std::uint64_t target = extensions.start;
            bool changed = false;
            while (index > 0 && target < result.blocks[index].start) {
                const Block previous = result.blocks[index - 1];
                if (target <= previous.start) {
                    result.blocks[index].start = previous.start;
                    result.blocks.erase(result.blocks.begin() +
                                        static_cast<std::ptrdiff_t>(index - 1));
                    --index;
                    ++result.absorbed_blocks;
                    changed = true;
                    continue;
                }
                if (target - previous.start < kMinimumRemainder) {
                    result.blocks[index].start = previous.start;
                    result.blocks.erase(result.blocks.begin() +
                                        static_cast<std::ptrdiff_t>(index - 1));
                    --index;
                    ++result.absorbed_blocks;
                } else {
                    result.blocks[index - 1].end = target;
                    result.blocks[index].start = target;
                }
                changed = true;
                break;
            }
            if (changed) {
                ++result.moved_breakpoints;
            }
        }

        if (extensions.end > result.blocks[index].end) {
            const std::uint64_t target = extensions.end;
            bool changed = false;
            while (index + 1 < result.blocks.size() &&
                   target > result.blocks[index].end) {
                const Block next = result.blocks[index + 1];
                if (target >= next.end) {
                    result.blocks[index].end = next.end;
                    result.blocks.erase(result.blocks.begin() +
                                        static_cast<std::ptrdiff_t>(index + 1));
                    ++result.absorbed_blocks;
                    changed = true;
                    continue;
                }
                if (next.end - target < kMinimumRemainder) {
                    result.blocks[index].end = next.end;
                    result.blocks.erase(result.blocks.begin() +
                                        static_cast<std::ptrdiff_t>(index + 1));
                    ++result.absorbed_blocks;
                } else {
                    result.blocks[index].end = target;
                    result.blocks[index + 1].start = target;
                }
                changed = true;
                break;
            }
            if (changed) {
                ++result.moved_breakpoints;
            }
        }
    }
    return result;
}

void validate_partition(const Contig& contig,
                        const std::vector<Block>& blocks) {
    if (contig.length == 0) {
        if (!blocks.empty()) {
            throw std::runtime_error("internal error: empty contig has blocks");
        }
        return;
    }
    if (blocks.empty() || blocks.front().start != 0 ||
        blocks.back().end != contig.length) {
        throw std::runtime_error("internal error: incomplete block partition for " +
                                 contig.name);
    }
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        if (blocks[i].start >= blocks[i].end ||
            (i != 0 && blocks[i - 1].end != blocks[i].start)) {
            throw std::runtime_error("internal error: invalid block partition for " +
                                     contig.name);
        }
    }
}

}  // namespace

BlockingStats create_blocks(const Options& options) {
    if (options.block_size < kMinimumRemainder) {
        throw std::runtime_error("--block-size must be at least 10000");
    }
    if (options.threads == 0) {
        throw std::runtime_error("--threads must be at least 1");
    }
    std::vector<Contig> contigs = read_fai(options.fai_path);
    BlockingStats stats;
    stats.contigs = contigs.size();
    for (const Contig& contig : contigs) {
        stats.bases += contig.length;
    }

    load_index(options.index_path, contigs, stats);

    std::vector<ContigResult> results(contigs.size());
    std::atomic<std::size_t> next_contig{0};
    std::atomic<bool> stop{false};
    std::exception_ptr worker_error;
    std::mutex error_mutex;
    const std::size_t thread_count = std::max<std::size_t>(
        1, std::min<std::size_t>(options.threads, contigs.size()));
    std::vector<std::thread> workers;
    workers.reserve(thread_count);
    for (std::size_t worker = 0; worker < thread_count; ++worker) {
        workers.emplace_back([&] {
            try {
                while (!stop.load(std::memory_order_relaxed)) {
                    const std::size_t index =
                        next_contig.fetch_add(1, std::memory_order_relaxed);
                    if (index >= contigs.size()) {
                        break;
                    }
                    results[index] = block_contig(contigs[index], options);
                    validate_partition(contigs[index], results[index].blocks);
                }
            } catch (...) {
                stop.store(true, std::memory_order_relaxed);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!worker_error) {
                    worker_error = std::current_exception();
                }
            }
        });
    }
    for (std::thread& worker : workers) {
        worker.join();
    }
    if (worker_error) {
        std::rethrow_exception(worker_error);
    }

    std::ofstream bed(options.output_path);
    if (!bed) {
        throw std::runtime_error("cannot open BED output " + options.output_path +
                                 ": " + std::strerror(errno));
    }
    for (std::size_t i = 0; i < contigs.size(); ++i) {
        stats.initial_blocks += results[i].initial_blocks;
        stats.final_blocks += results[i].blocks.size();
        stats.hotspots += results[i].hotspots;
        stats.moved_breakpoints += results[i].moved_breakpoints;
        stats.absorbed_blocks += results[i].absorbed_blocks;
        for (const Block& block : results[i].blocks) {
            bed << contigs[i].name << '\t' << block.start << '\t' << block.end
                << '\n';
        }
    }
    bed.flush();
    if (!bed) {
        throw std::runtime_error("failed while writing BED output " +
                                 options.output_path);
    }
    return stats;
}

}  // namespace kmer_blocking
