#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr std::uint32_t kNoBlock = std::numeric_limits<std::uint32_t>::max();
constexpr std::uint32_t kMarkerMinimum = 0xffff0000U;
constexpr std::uint32_t kMaximumUint32 = 0xffffffffU;
constexpr std::uint64_t kMinimumBlockLength = 5'000;
constexpr std::uint64_t kLowExclusiveLimit = 1'000;

struct Contig {
    std::string name;
    std::uint64_t length = 0;
    std::uint64_t global_offset = 0;
};

struct Block {
    std::string name;
    std::string contig;
    std::uint64_t start = 0;
    std::uint64_t end = 0;
    std::uint64_t global_start = 0;
    std::uint64_t exclusive_kmers = 0;
};

std::vector<Contig> read_fai(
    const std::string& path,
    std::unordered_map<std::string, std::size_t>& contig_by_name,
    std::uint64_t& total_bases) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open FAI: " + path);
    }

    std::vector<Contig> contigs;
    std::string line;
    total_bases = 0;
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        std::istringstream fields(line);
        Contig contig;
        fields >> contig.name >> contig.length;
        if (!fields || contig_by_name.count(contig.name)) {
            throw std::runtime_error("invalid or duplicate FAI contig: " + line);
        }
        contig.global_offset = total_bases;
        contig_by_name.emplace(contig.name, contigs.size());
        contigs.push_back(contig);
        total_bases += contig.length;
    }
    if (contigs.empty()) {
        throw std::runtime_error("FAI contains no contigs");
    }
    if (total_bases >= kMarkerMinimum) {
        throw std::runtime_error(
            "strong marker assumption failed: assembly reaches 0xffff0000");
    }
    return contigs;
}

std::vector<Block> read_blocks(
    const std::string& path, const std::vector<Contig>& contigs,
    const std::unordered_map<std::string, std::size_t>& contig_by_name,
    std::uint64_t& discarded_small) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open BED: " + path);
    }

    std::vector<Block> blocks;
    std::vector<std::uint64_t> previous_end(contigs.size(), 0);
    std::vector<bool> seen(contigs.size(), false);
    std::string line;
    std::uint64_t line_number = 0;
    discarded_small = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream fields(line);
        Block block;
        fields >> block.contig >> block.start >> block.end;
        if (!fields || block.start >= block.end) {
            throw std::runtime_error("invalid BED line " +
                                     std::to_string(line_number));
        }
        const auto found = contig_by_name.find(block.contig);
        if (found == contig_by_name.end()) {
            throw std::runtime_error("BED contig is absent from FAI: " +
                                     block.contig);
        }
        const std::size_t contig_index = found->second;
        const Contig& contig = contigs[contig_index];
        if (block.end > contig.length ||
            (seen[contig_index] && block.start < previous_end[contig_index])) {
            throw std::runtime_error("overlapping or out-of-range BED block at line " +
                                     std::to_string(line_number));
        }
        seen[contig_index] = true;
        previous_end[contig_index] = block.end;
        if (block.end - block.start < kMinimumBlockLength) {
            ++discarded_small;
            continue;
        }
        if (!(fields >> block.name)) {
            block.name = block.contig + ":" + std::to_string(block.start) + "-" +
                         std::to_string(block.end);
        }
        block.global_start = contig.global_offset + block.start;
        blocks.push_back(std::move(block));
    }
    return blocks;
}

class MappedFile {
public:
    explicit MappedFile(const std::string& path) {
        descriptor_ = open(path.c_str(), O_RDONLY);
        if (descriptor_ < 0) {
            throw std::runtime_error("cannot open index: " + path + ": " +
                                     std::strerror(errno));
        }
        struct stat status {};
        if (fstat(descriptor_, &status) != 0 || status.st_size <= 0) {
            close(descriptor_);
            descriptor_ = -1;
            throw std::runtime_error("cannot stat or empty index: " + path);
        }
        size_ = static_cast<std::size_t>(status.st_size);
        void* mapped = mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, descriptor_, 0);
        if (mapped == MAP_FAILED) {
            close(descriptor_);
            descriptor_ = -1;
            throw std::runtime_error("cannot mmap index: " +
                                     std::string(std::strerror(errno)));
        }
        data_ = static_cast<const unsigned char*>(mapped);
#ifdef MADV_SEQUENTIAL
        madvise(const_cast<unsigned char*>(data_), size_, MADV_SEQUENTIAL);
#endif
    }

    ~MappedFile() {
        if (data_ != nullptr) {
            munmap(const_cast<unsigned char*>(data_), size_);
        }
        if (descriptor_ >= 0) {
            close(descriptor_);
        }
    }

    const unsigned char* data() const { return data_; }
    std::size_t size() const { return size_; }

private:
    int descriptor_ = -1;
    const unsigned char* data_ = nullptr;
    std::size_t size_ = 0;
};

std::uint32_t load_u32(const unsigned char* source) {
    std::uint32_t value;
    std::memcpy(&value, source, sizeof(value));
    return value;
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc != 5 && argc != 6) {
            std::fprintf(
                stderr,
                "Usage: %s INDEX.bin ASSEMBLY.fai BLOCKS.bed BLOCK_BINS.tsv "
                "[THREADS]\n",
                argv[0]);
            return EXIT_FAILURE;
        }
        const int threads = argc == 6 ? std::max(1, std::atoi(argv[5])) : 64;
#ifdef _OPENMP
        omp_set_num_threads(threads);
#else
        (void)threads;
#endif

        std::unordered_map<std::string, std::size_t> contig_by_name;
        std::uint64_t total_bases = 0;
        const std::vector<Contig> contigs =
            read_fai(argv[2], contig_by_name, total_bases);
        std::uint64_t discarded_small = 0;
        std::vector<Block> blocks = read_blocks(
            argv[3], contigs, contig_by_name, discarded_small);
        if (blocks.size() >= kNoBlock) {
            throw std::runtime_error("too many retained BED blocks");
        }

        std::fprintf(stderr,
                     "Allocating %.2f GiB position-to-block table for %llu "
                     "bases...\n",
                     static_cast<double>(total_bases * sizeof(std::uint32_t)) /
                         1073741824.0,
                     static_cast<unsigned long long>(total_bases));
        std::vector<std::uint32_t> position_to_block(
            static_cast<std::size_t>(total_bases), kNoBlock);

#pragma omp parallel for schedule(dynamic, 1)
        for (std::int64_t i = 0;
             i < static_cast<std::int64_t>(blocks.size()); ++i) {
            const Block& block = blocks[static_cast<std::size_t>(i)];
            std::fill(position_to_block.begin() +
                          static_cast<std::ptrdiff_t>(block.global_start),
                      position_to_block.begin() +
                          static_cast<std::ptrdiff_t>(block.global_start +
                                                      block.end - block.start),
                      static_cast<std::uint32_t>(i));
        }

        std::vector<std::uint64_t> block_counts(blocks.size(), 0);
        const MappedFile index(argv[1]);
        const unsigned char* data = index.data();
        const std::size_t size = index.size();
        std::size_t offset = 0;
        std::uint64_t index_kmers = 0;
        std::uint64_t exclusive_kmers = 0;
        std::uint32_t maximum_occurrence = 1;

        std::fprintf(stderr,
                     "Pass 1/2: selecting blocks from %.2f GiB index...\n",
                     static_cast<double>(size) / 1073741824.0);
        while (offset < size) {
            const std::uint32_t value = load_u32(data + offset + 8);
            offset += 12;
            ++index_kmers;
            if (value >= kMarkerMinimum) {
                const std::uint32_t occurrence = kMaximumUint32 - value;
                maximum_occurrence = std::max(maximum_occurrence, occurrence);
                offset += static_cast<std::size_t>(occurrence) * 4;
                continue;
            }
            ++exclusive_kmers;
            const std::uint32_t block_index = position_to_block[value];
            if (block_index != kNoBlock) {
                ++block_counts[block_index];
            }
        }

        std::vector<std::uint32_t> low_index_for_block(blocks.size(), kNoBlock);
        std::vector<std::size_t> low_blocks;
        for (std::size_t i = 0; i < blocks.size(); ++i) {
            blocks[i].exclusive_kmers = block_counts[i];
            if (block_counts[i] < kLowExclusiveLimit) {
                low_index_for_block[i] =
                    static_cast<std::uint32_t>(low_blocks.size());
                low_blocks.push_back(i);
            }
        }

        const std::uint64_t bins_per_block =
            static_cast<std::uint64_t>(maximum_occurrence) + 1;
        if (!low_blocks.empty() &&
            bins_per_block > std::numeric_limits<std::size_t>::max() /
                                 low_blocks.size()) {
            throw std::runtime_error("per-block bin table is too large");
        }
        std::vector<std::uint64_t> per_block_bins(
            static_cast<std::size_t>(bins_per_block * low_blocks.size()), 0);

        std::fprintf(stderr,
                     "Pass 2/2: counting bins 1..%u for %zu selected blocks...\n",
                     maximum_occurrence, low_blocks.size());
        offset = 0;
        while (offset < size) {
            const std::uint32_t value = load_u32(data + offset + 8);
            offset += 12;
            if (value < kMarkerMinimum) {
                const std::uint32_t block_index = position_to_block[value];
                if (block_index != kNoBlock) {
                    const std::uint32_t low_index =
                        low_index_for_block[block_index];
                    if (low_index != kNoBlock) {
                        ++per_block_bins[static_cast<std::size_t>(low_index) *
                                             bins_per_block +
                                         1];
                    }
                }
                continue;
            }

            const std::uint32_t occurrence = kMaximumUint32 - value;
            for (std::uint32_t i = 0; i < occurrence; ++i) {
                const std::uint32_t position = load_u32(data + offset + i * 4);
                const std::uint32_t block_index = position_to_block[position];
                if (block_index == kNoBlock) {
                    continue;
                }
                const std::uint32_t low_index =
                    low_index_for_block[block_index];
                if (low_index != kNoBlock) {
                    ++per_block_bins[static_cast<std::size_t>(low_index) *
                                         bins_per_block +
                                     occurrence];
                }
            }
            offset += static_cast<std::size_t>(occurrence) * 4;
        }

        std::ofstream report(argv[4]);
        if (!report) {
            throw std::runtime_error("cannot open block-bin report: " +
                                     std::string(argv[4]));
        }
        report << "block_name\tbin\tkmer_positions\n";
        for (std::size_t low_index = 0; low_index < low_blocks.size();
             ++low_index) {
            const Block& block = blocks[low_blocks[low_index]];
            const std::size_t base =
                static_cast<std::size_t>(low_index * bins_per_block);
            for (std::uint32_t occurrence = 1;
                 occurrence <= maximum_occurrence; ++occurrence) {
                const std::uint64_t count = per_block_bins[base + occurrence];
                if (count != 0) {
                    report << block.name << "\tbin_" << occurrence << '\t'
                           << count << '\n';
                }
            }
        }
        std::fprintf(
            stderr,
            "Index k-mers: %llu; exclusive index k-mers: %llu; retained BED "
            "blocks: %zu; discarded blocks under 5000 bp: %llu; blocks under "
            "1000 exclusive k-mers: %llu\n",
            static_cast<unsigned long long>(index_kmers),
            static_cast<unsigned long long>(exclusive_kmers), blocks.size(),
            static_cast<unsigned long long>(discarded_small),
            static_cast<unsigned long long>(low_blocks.size()));
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "temp_exclusive_bins: error: %s\n", error.what());
        return EXIT_FAILURE;
    }
}
