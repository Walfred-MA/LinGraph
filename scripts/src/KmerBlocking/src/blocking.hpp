#ifndef KMER_BLOCKING_HPP
#define KMER_BLOCKING_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace kmer_blocking {

constexpr std::uint32_t kKmerLength = 31;
constexpr std::uint32_t kConnectDistance = 150;
constexpr std::uint32_t kHotspotHitCount = 100;
constexpr std::uint32_t kHotspotSpan = 500;
constexpr std::uint32_t kMinimumRemainder = 10'000;

struct Options {
    std::string index_path;
    std::string fai_path;
    std::string output_path;
    std::uint64_t block_size = 50'000;
    std::uint64_t search_distance = 10'000;
    std::uint32_t threads = 1;
};

struct BlockingStats {
    std::uint64_t contigs = 0;
    std::uint64_t bases = 0;
    std::uint64_t indexed_kmers = 0;
    std::uint64_t indexed_locations = 0;
    std::uint64_t initial_blocks = 0;
    std::uint64_t final_blocks = 0;
    std::uint64_t hotspots = 0;
    std::uint64_t moved_breakpoints = 0;
    std::uint64_t absorbed_blocks = 0;
};

BlockingStats create_blocks(const Options& options);

}  // namespace kmer_blocking

#endif
