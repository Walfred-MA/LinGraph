#include "block_hash.hpp"

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <fcntl.h>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

namespace kmer_blocking {
namespace {

constexpr std::uint32_t kMarkerMinimum = 0xffff0000U;
constexpr std::uint32_t kMaximumUint32 = 0xffffffffU;
constexpr std::uint64_t kBoundaryMarker =
    std::numeric_limits<std::uint64_t>::max();
constexpr std::uint32_t kSearcherHash1Size = 0x40000000U;
constexpr std::uint32_t kSearcherHashMask = 0x3fffffffU;
constexpr std::uint32_t kSearcherValueBuckets = 1U << 16;
constexpr std::uint32_t kMaximumSearcherBlockIndex = 0xfffeU;
constexpr std::uint64_t kSearcherHeaderSize =
    sizeof(std::uint64_t) * 4 + sizeof(std::uint32_t);

#pragma pack(push, 1)
struct KmerPosition {
    std::uint64_t kmer;
    std::uint32_t position;
};

struct HashRecord {
    std::uint32_t hash1;
    std::uint32_t hash2;
    std::uint32_t block_index;
};
#pragma pack(pop)

static_assert(sizeof(KmerPosition) == 12,
              "KmerPosition must occupy exactly 12 bytes");
static_assert(sizeof(HashRecord) == 12,
              "HashRecord must occupy exactly 12 bytes");

struct Contig {
    std::string name;
    std::uint64_t length = 0;
    std::uint64_t global_offset = 0;
};

struct Block {
    std::string contig;
    std::uint64_t start = 0;
    std::uint64_t end = 0;
    std::uint64_t global_start = 0;
    std::uint64_t global_end = 0;
    std::uint32_t index = 0;
    std::uint64_t unique_kmers = 0;
};

std::uint32_t read_u32(const unsigned char* bytes) {
    return static_cast<std::uint32_t>(bytes[0]) |
           (static_cast<std::uint32_t>(bytes[1]) << 8) |
           (static_cast<std::uint32_t>(bytes[2]) << 16) |
           (static_cast<std::uint32_t>(bytes[3]) << 24);
}

std::uint64_t read_u64(const unsigned char* bytes) {
    std::uint64_t value = 0;
    for (int i = 0; i < 8; ++i) {
        value |= static_cast<std::uint64_t>(bytes[i]) << (8 * i);
    }
    return value;
}

std::uint32_t rotate_right_32(std::uint32_t value, unsigned amount) {
    return (value >> amount) | (value << ((32 - amount) & 31));
}

std::uint32_t mix_32_to_30(std::uint32_t value) {
    const std::uint32_t mixed =
        rotate_right_32(value, 7) ^ rotate_right_32(value, 13) ^
        rotate_right_32(value, 19);
    return mixed & kSearcherHashMask;
}

void searcher_hash(std::uint64_t kmer, std::uint32_t& hash1,
                   std::uint32_t& hash2) {
    hash1 = static_cast<std::uint32_t>(kmer) & kSearcherHashMask;
    hash2 = static_cast<std::uint32_t>(kmer >> 30);
    hash1 ^= mix_32_to_30(hash2);
}

void store_u16(unsigned char* destination, std::uint16_t value) {
    destination[0] = static_cast<unsigned char>(value & 0xff);
    destination[1] = static_cast<unsigned char>((value >> 8) & 0xff);
}

void store_u32(unsigned char* destination, std::uint32_t value) {
    for (int i = 0; i < 4; ++i) {
        destination[i] =
            static_cast<unsigned char>((value >> (8 * i)) & 0xff);
    }
}

void store_u64(unsigned char* destination, std::uint64_t value) {
    for (int i = 0; i < 8; ++i) {
        destination[i] =
            static_cast<unsigned char>((value >> (8 * i)) & 0xff);
    }
}

std::uint64_t checked_add(std::uint64_t left, std::uint64_t right,
                          const char* description) {
    if (right > std::numeric_limits<std::uint64_t>::max() - left) {
        throw std::runtime_error(std::string("size overflow for ") +
                                 description);
    }
    return left + right;
}

std::uint64_t checked_multiply(std::uint64_t left, std::uint64_t right,
                               const char* description) {
    if (left != 0 &&
        right > std::numeric_limits<std::uint64_t>::max() / left) {
        throw std::runtime_error(std::string("size overflow for ") +
                                 description);
    }
    return left * right;
}

std::vector<Contig> read_fai(
    const std::string& path,
    std::unordered_map<std::string, std::size_t>& by_name,
    std::uint64_t& total_bases) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open FASTA index " + path);
    }
    std::vector<Contig> contigs;
    total_bases = 0;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        std::istringstream fields(line);
        Contig contig;
        fields >> contig.name >> contig.length;
        if (!fields || contig.name.empty() || by_name.count(contig.name)) {
            throw std::runtime_error("invalid or duplicate .fai line: " + line);
        }
        if (contig.length > kMaximumUint32 - total_bases) {
            throw std::runtime_error(
                "total assembly length must be less than 2^32 bases");
        }
        contig.global_offset = total_bases;
        by_name.emplace(contig.name, contigs.size());
        contigs.push_back(contig);
        total_bases += contig.length;
    }
    if (contigs.empty()) {
        throw std::runtime_error("FASTA index is empty: " + path);
    }
    return contigs;
}

std::vector<Block> read_bed(
    const std::string& path, const std::vector<Contig>& contigs,
    const std::unordered_map<std::string, std::size_t>& by_name) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open block BED " + path);
    }
    std::vector<Block> blocks;
    std::string line;
    std::uint64_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream fields(line);
        Block block;
        fields >> block.contig >> block.start >> block.end;
        if (!fields || block.start >= block.end) {
            throw std::runtime_error("invalid BED coordinates on line " +
                                     std::to_string(line_number));
        }
        const auto found = by_name.find(block.contig);
        if (found == by_name.end()) {
            throw std::runtime_error("BED contig is absent from .fai: " +
                                     block.contig);
        }
        const Contig& contig = contigs[found->second];
        if (block.end > contig.length) {
            throw std::runtime_error("BED interval exceeds contig on line " +
                                     std::to_string(line_number));
        }
        block.global_start = contig.global_offset + block.start;
        block.global_end = contig.global_offset + block.end;
        blocks.push_back(std::move(block));
    }
    if (blocks.empty()) {
        throw std::runtime_error("block BED contains no intervals");
    }

    std::sort(blocks.begin(), blocks.end(),
              [](const Block& left, const Block& right) {
                  if (left.global_start != right.global_start) {
                      return left.global_start < right.global_start;
                  }
                  return left.global_end < right.global_end;
              });
    if (blocks.size() >= kMaximumUint32) {
        throw std::runtime_error("too many blocks for uint32 block indexes");
    }
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        if (i != 0 && blocks[i].global_start < blocks[i - 1].global_end) {
            throw std::runtime_error("block BED contains overlapping intervals");
        }
        blocks[i].index = static_cast<std::uint32_t>(i + 1);
    }
    return blocks;
}

std::vector<KmerPosition> load_unique_kmers(const std::string& path) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) {
        throw std::runtime_error("cannot open k-mer index " + path);
    }
    const std::streamoff file_size = input.tellg();
    if (file_size < 0) {
        throw std::runtime_error("cannot determine k-mer index size");
    }
    const std::size_t size = static_cast<std::size_t>(file_size);
    input.seekg(0, std::ios::beg);
    std::vector<unsigned char> data(size);
    input.read(reinterpret_cast<char*>(data.data()),
               static_cast<std::streamsize>(size));
    if (input.gcount() != static_cast<std::streamsize>(size)) {
        throw std::runtime_error("failed while reading k-mer index " + path);
    }

    std::vector<KmerPosition> unique;
    unique.reserve(size / 12);
    std::size_t offset = 0;
    while (offset < size) {
        if (size - offset < 12) {
            throw std::runtime_error("truncated k-mer record header");
        }
        const std::uint64_t kmer = read_u64(data.data() + offset);
        const std::uint32_t value = read_u32(data.data() + offset + 8);
        offset += 12;
        if (value < kMarkerMinimum) {
            unique.push_back({kmer, value});
            continue;
        }
        const std::uint32_t count = kMaximumUint32 - value;
        const std::size_t location_bytes = static_cast<std::size_t>(count) * 4;
        if (location_bytes > size - offset) {
            throw std::runtime_error("truncated multiple-location k-mer record");
        }
        offset += location_bytes;
    }
    return unique;
}

void assign_blocks(std::vector<KmerPosition>& events,
                   std::vector<Block>& blocks,
                   BlockHashStats& stats) {
    const std::size_t kmer_count = events.size();
    events.reserve(events.size() + blocks.size());
    for (const Block& block : blocks) {
        events.push_back({kBoundaryMarker,
                          static_cast<std::uint32_t>(block.global_start)});
    }
    std::sort(events.begin(), events.end(),
              [](const KmerPosition& left, const KmerPosition& right) {
                  if (left.position != right.position) {
                      return left.position < right.position;
                  }
                  return left.kmer > right.kmer;
              });

    std::size_t next_block = 0;
    const Block* active = nullptr;
    std::size_t write_index = 0;
    for (KmerPosition event : events) {
        if (event.kmer == kBoundaryMarker) {
            if (next_block >= blocks.size() ||
                blocks[next_block].global_start != event.position) {
                throw std::runtime_error("internal block boundary marker mismatch");
            }
            active = &blocks[next_block++];
            continue;
        }
        if (active == nullptr || event.position < active->global_start ||
            event.position >= active->global_end) {
            ++stats.unique_kmers_outside_blocks;
            continue;
        }
        ++blocks[active->index - 1].unique_kmers;
        event.position = active->index;
        events[write_index++] = event;
    }
    events.resize(write_index);
    stats.unique_kmers_loaded = kmer_count;
    stats.unique_kmers_assigned = write_index;
}

void write_annotated_bed(const std::string& path,
                         const std::vector<Block>& blocks) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("cannot open annotated BED output " + path);
    }
    for (const Block& block : blocks) {
        output << block.contig << '\t' << block.start << '\t' << block.end << '\t'
               << block.index << '\t' << block.unique_kmers << '\n';
    }
    output.flush();
    if (!output) {
        throw std::runtime_error("failed while writing annotated BED " + path);
    }
}

void write_hash(const std::string& path,
                const std::vector<KmerPosition>& assigned,
                std::uint64_t block_count, std::uint32_t threads) {
    if (block_count > kMaximumSearcherBlockIndex) {
        throw std::runtime_error(
            "KmerSearcher-compatible hashes support at most 65534 blocks");
    }
    if (assigned.size() >
        static_cast<std::size_t>(kMaximumUint32 - 1U)) {
        throw std::runtime_error(
            "too many k-mers for KmerSearcher's uint32 hash redirects");
    }

    std::vector<HashRecord> records(assigned.size());
    const std::size_t thread_count = std::max<std::size_t>(
        1, std::min<std::size_t>(threads, assigned.empty() ? 1 : assigned.size()));
    std::vector<std::thread> workers;
    workers.reserve(thread_count);
    for (std::size_t thread_index = 0; thread_index < thread_count;
         ++thread_index) {
        workers.emplace_back([&, thread_index] {
            const std::size_t begin = assigned.size() * thread_index / thread_count;
            const std::size_t end =
                assigned.size() * (thread_index + 1) / thread_count;
            for (std::size_t i = begin; i < end; ++i) {
                const KmerPosition& kmer = assigned[i];
                std::uint32_t hash1 = 0;
                std::uint32_t hash2 = 0;
                searcher_hash(kmer.kmer, hash1, hash2);
                records[i] = {hash1, hash2, kmer.position};
            }
        });
    }
    for (std::thread& worker : workers) {
        worker.join();
    }
    if (std::any_of(records.begin(), records.end(),
                    [](const HashRecord& record) {
                        return record.hash1 >= kSearcherHash1Size;
                    })) {
        throw std::runtime_error(
            "a k-mer maps to the out-of-range bucket in KmerSearcher's "
            "current 30-bit hash table");
    }
    std::sort(records.begin(), records.end(),
              [](const HashRecord& left, const HashRecord& right) {
                  if (left.hash1 != right.hash1) {
                      return left.hash1 < right.hash1;
                  }
                  return left.hash2 < right.hash2;
              });

    for (std::size_t begin = 0; begin < records.size();) {
        std::size_t end = begin + 1;
        while (end < records.size() &&
               records[end].hash1 == records[begin].hash1) {
            ++end;
        }
        if (end - begin > std::numeric_limits<std::uint8_t>::max()) {
            throw std::runtime_error(
                "a KmerSearcher hash bucket contains more than 255 k-mers");
        }
        begin = end;
    }

    // Match Kmer_hash::savehash exactly. Empty regions are left as sparse
    // zero-filled ranges, avoiding KmerSearcher's approximately 5 GiB fixed
    // table cost in resident memory while producing the same bytes when read.
    const std::uint64_t hash2_entries =
        checked_add(records.size(), 11, "Hash2 entries");
    const std::uint64_t total_kmers =
        checked_add(records.size(), 1, "total k-mers");
    const std::uint64_t total_targets =
        checked_add(block_count, 1, "total targets");

    const std::uint64_t hash1_offset = kSearcherHeaderSize;
    const std::uint64_t hash2_offset = checked_add(
        hash1_offset,
        checked_multiply(kSearcherHash1Size, sizeof(std::uint32_t), "Hash1"),
        "Hash2 offset");
    const std::uint64_t hash2_size_offset = checked_add(
        hash2_offset,
        checked_multiply(hash2_entries, sizeof(std::uint32_t), "Hash2"),
        "Hash2_size offset");
    const std::uint64_t values_offset = checked_add(
        hash2_size_offset, kSearcherHash1Size, "values offset");
    const std::uint64_t values_mul_sizes_offset = checked_add(
        values_offset,
        checked_multiply(hash2_entries, sizeof(std::uint16_t), "values"),
        "multi-value sizes offset");
    const std::uint64_t values_mul_offset = checked_add(
        values_mul_sizes_offset,
        checked_multiply(kSearcherValueBuckets, sizeof(std::uint32_t),
                         "multi-value sizes"),
        "multi-value buckets offset");
    const std::uint64_t output_size = checked_add(
        values_mul_offset,
        checked_multiply(kSearcherValueBuckets, sizeof(std::uint64_t),
                         "empty multi-value buckets"),
        "output size");
    if (output_size >
        static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()) ||
        output_size >
            static_cast<std::uint64_t>(std::numeric_limits<off_t>::max())) {
        throw std::runtime_error(
            "KmerSearcher-compatible hash is too large for this platform");
    }

    const std::string temporary_path = path + ".tmp";
    const int descriptor =
        ::open(temporary_path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
    if (descriptor < 0) {
        throw std::runtime_error("cannot open block hash output " +
                                 temporary_path + ": " + std::strerror(errno));
    }
    if (::ftruncate(descriptor, static_cast<off_t>(output_size)) != 0) {
        const std::string message = std::strerror(errno);
        ::close(descriptor);
        ::unlink(temporary_path.c_str());
        throw std::runtime_error("cannot size block hash output " +
                                 temporary_path + ": " + message);
    }

    void* mapped =
        ::mmap(nullptr, static_cast<std::size_t>(output_size),
               PROT_READ | PROT_WRITE, MAP_SHARED, descriptor, 0);
    if (mapped == MAP_FAILED) {
        const std::string message = std::strerror(errno);
        ::close(descriptor);
        ::unlink(temporary_path.c_str());
        throw std::runtime_error("cannot map block hash output " +
                                 temporary_path + ": " + message);
    }
    unsigned char* bytes = static_cast<unsigned char*>(mapped);

    store_u64(bytes, hash2_entries);
    store_u64(bytes + 8, kSearcherHash1Size);
    store_u64(bytes + 16, total_kmers);
    store_u64(bytes + 24, total_targets);
    store_u32(bytes + 32, 0);  // KmerSearcher's bucketcollision.

    std::uint32_t hash2_index = 1;
    for (std::size_t begin = 0; begin < records.size();) {
        std::size_t end = begin + 1;
        while (end < records.size() &&
               records[end].hash1 == records[begin].hash1) {
            ++end;
        }
        const std::uint32_t bucket_size =
            static_cast<std::uint32_t>(end - begin);
        store_u32(bytes + hash1_offset +
                      static_cast<std::uint64_t>(records[begin].hash1) * 4,
                  hash2_index);
        bytes[hash2_size_offset + records[begin].hash1] =
            static_cast<unsigned char>(bucket_size);
        for (std::size_t i = begin; i < end; ++i, ++hash2_index) {
            store_u32(bytes + hash2_offset +
                          static_cast<std::uint64_t>(hash2_index) * 4,
                      records[i].hash2);
            store_u16(bytes + values_offset +
                          static_cast<std::uint64_t>(hash2_index) * 2,
                      static_cast<std::uint16_t>(records[i].block_index));
        }
        begin = end;
    }

    if (::msync(mapped, static_cast<std::size_t>(output_size), MS_SYNC) != 0) {
        const std::string message = std::strerror(errno);
        ::munmap(mapped, static_cast<std::size_t>(output_size));
        ::close(descriptor);
        ::unlink(temporary_path.c_str());
        throw std::runtime_error("cannot flush block hash output " +
                                 temporary_path + ": " + message);
    }
    if (::munmap(mapped, static_cast<std::size_t>(output_size)) != 0) {
        const std::string message = std::strerror(errno);
        ::close(descriptor);
        ::unlink(temporary_path.c_str());
        throw std::runtime_error("cannot unmap block hash output " +
                                 temporary_path + ": " + message);
    }
    if (::close(descriptor) != 0) {
        const std::string message = std::strerror(errno);
        ::unlink(temporary_path.c_str());
        throw std::runtime_error("cannot close block hash output " +
                                 temporary_path + ": " + message);
    }
    if (::rename(temporary_path.c_str(), path.c_str()) != 0) {
        const std::string message = std::strerror(errno);
        ::unlink(temporary_path.c_str());
        throw std::runtime_error("cannot finalize block hash output " + path +
                                 ": " + message);
    }
}

}  // namespace

BlockHashStats create_block_hash(const BlockHashOptions& options) {
    if (options.bed_path.empty() || options.index_path.empty() ||
        options.fai_path.empty()) {
        throw std::runtime_error("-b, -k, and -f are required in block-hash mode");
    }
    if (options.annotated_bed_path.empty() && options.hash_path.empty()) {
        throw std::runtime_error("block-hash mode requires -ob and/or -ok");
    }

    std::unordered_map<std::string, std::size_t> contig_by_name;
    std::uint64_t total_bases = 0;
    const std::vector<Contig> contigs =
        read_fai(options.fai_path, contig_by_name, total_bases);
    if (total_bases >= kMarkerMinimum) {
        throw std::runtime_error(
            "assembly coordinates overlap the assumed index marker range");
    }
    std::vector<Block> blocks =
        read_bed(options.bed_path, contigs, contig_by_name);
    std::fprintf(stderr, "Loading count-1 k-mers from %s...\n",
                 options.index_path.c_str());
    std::vector<KmerPosition> kmers = load_unique_kmers(options.index_path);

    BlockHashStats stats;
    stats.blocks = blocks.size();
    std::fprintf(stderr,
                 "Sorting %llu count-1 k-mers plus %llu block markers...\n",
                 static_cast<unsigned long long>(kmers.size()),
                 static_cast<unsigned long long>(blocks.size()));
    assign_blocks(kmers, blocks, stats);
    if (!options.annotated_bed_path.empty()) {
        write_annotated_bed(options.annotated_bed_path, blocks);
    }
    if (!options.hash_path.empty()) {
        std::fprintf(stderr, "Building and writing block hash to %s...\n",
                     options.hash_path.c_str());
        write_hash(options.hash_path, kmers, blocks.size(), options.threads);
    }
    return stats;
}

}  // namespace kmer_blocking
