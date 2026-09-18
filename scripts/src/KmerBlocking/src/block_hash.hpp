#ifndef KMER_BLOCK_HASH_HPP
#define KMER_BLOCK_HASH_HPP

#include <cstdint>
#include <string>

namespace kmer_blocking {

struct BlockHashOptions {
    std::string bed_path;
    std::string index_path;
    std::string fai_path;
    std::string annotated_bed_path;
    std::string hash_path;
    std::uint32_t threads = 1;
};

struct BlockHashStats {
    std::uint64_t blocks = 0;
    std::uint64_t unique_kmers_loaded = 0;
    std::uint64_t unique_kmers_assigned = 0;
    std::uint64_t unique_kmers_outside_blocks = 0;
};

BlockHashStats create_block_hash(const BlockHashOptions& options);

}  // namespace kmer_blocking

#endif
