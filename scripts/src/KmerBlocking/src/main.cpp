#include "blocking.hpp"
#include "block_hash.hpp"

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>

namespace {

void usage(std::ostream& output, const char* program) {
    output
        << "Usage: " << program
        << " --index assembly.k31.idx --fai assembly.fa.fai --output blocks.bed [options]\n"
        << "\n"
        << "Required:\n"
        << "  -i, --index PATH           KmerIndex binary file\n"
        << "  -f, --fai PATH             FASTA .fai used to build the index\n"
        << "  -o, --output PATH          BED3 block output\n"
        << "\n"
        << "Options:\n"
        << "  -S, --block-size N         Initial block size (default: 50000)\n"
        << "  -D, --search-distance N    External flank search distance\n"
        << "                               (default: 10000)\n"
        << "  -t, --threads N            Contig worker threads (default: CPU count)\n"
        << "  -h, --help                 Show this help\n";
    output
        << "\nIndependent unique-k-mer/block-hash mode:\n"
        << "  -b PATH                    Input block BED\n"
        << "  -k PATH                    KmerIndex binary file\n"
        << "  -f PATH                    Matching FASTA .fai\n"
        << "  -ob PATH                   Optional BED5 with unique-k-mer score\n"
        << "  -ok PATH                   Optional KmerSearcher-compatible cache\n"
        << "                             (block indexes are target IDs)\n";
}

std::uint64_t parse_number(const std::string& text,
                           const std::string& option) {
    if (text.empty() || text[0] == '-') {
        throw std::runtime_error("invalid value for " + option + ": " + text);
    }
    std::size_t used = 0;
    unsigned long long value = 0;
    try {
        value = std::stoull(text, &used, 10);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid value for " + option + ": " + text);
    }
    if (used != text.size()) {
        throw std::runtime_error("invalid value for " + option + ": " + text);
    }
    return static_cast<std::uint64_t>(value);
}

std::string next_value(int argc, char* argv[], int& index,
                       const std::string& option) {
    if (++index >= argc) {
        throw std::runtime_error("missing value after " + option);
    }
    return argv[index];
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        kmer_blocking::Options options;
        kmer_blocking::BlockHashOptions hash_options;
        const unsigned hardware_threads = std::thread::hardware_concurrency();
        options.threads = hardware_threads == 0 ? 1 : hardware_threads;
        hash_options.threads = options.threads;

        for (int i = 1; i < argc; ++i) {
            const std::string argument = argv[i];
            if (argument == "-h" || argument == "--help") {
                usage(std::cout, argv[0]);
                return EXIT_SUCCESS;
            }
            if (argument == "-i" || argument == "--index") {
                options.index_path = next_value(argc, argv, i, argument);
            } else if (argument == "-f" || argument == "--fai") {
                options.fai_path = next_value(argc, argv, i, argument);
                hash_options.fai_path = options.fai_path;
            } else if (argument == "-o" || argument == "--output") {
                options.output_path = next_value(argc, argv, i, argument);
            } else if (argument == "-b") {
                hash_options.bed_path = next_value(argc, argv, i, argument);
            } else if (argument == "-k") {
                hash_options.index_path = next_value(argc, argv, i, argument);
            } else if (argument == "-ob") {
                hash_options.annotated_bed_path =
                    next_value(argc, argv, i, argument);
            } else if (argument == "-ok") {
                hash_options.hash_path = next_value(argc, argv, i, argument);
            } else if (argument == "-S" || argument == "--block-size") {
                options.block_size = parse_number(
                    next_value(argc, argv, i, argument), argument);
            } else if (argument == "-D" ||
                       argument == "--search-distance") {
                options.search_distance = parse_number(
                    next_value(argc, argv, i, argument), argument);
            } else if (argument == "-t" || argument == "--threads") {
                const std::uint64_t value = parse_number(
                    next_value(argc, argv, i, argument), argument);
                if (value == 0 || value > std::numeric_limits<std::uint32_t>::max()) {
                    throw std::runtime_error("--threads must be between 1 and " +
                                             std::to_string(
                                                 std::numeric_limits<std::uint32_t>::max()));
                }
                options.threads = static_cast<std::uint32_t>(value);
                hash_options.threads = options.threads;
            } else {
                throw std::runtime_error("unknown option: " + argument);
            }
        }

        const bool block_hash_mode = !hash_options.bed_path.empty() ||
                                     !hash_options.index_path.empty() ||
                                     !hash_options.annotated_bed_path.empty() ||
                                     !hash_options.hash_path.empty();
        if (block_hash_mode) {
            if (!options.index_path.empty() || !options.output_path.empty()) {
                throw std::runtime_error(
                    "do not mix -i/-o blocking options with -b/-k hash mode");
            }
            const kmer_blocking::BlockHashStats stats =
                kmer_blocking::create_block_hash(hash_options);
            std::fprintf(
                stderr,
                "Blocks: %llu\n"
                "Count-1 k-mers loaded: %llu\n"
                "Count-1 k-mers assigned: %llu\n"
                "Count-1 k-mers outside blocks: %llu\n",
                static_cast<unsigned long long>(stats.blocks),
                static_cast<unsigned long long>(stats.unique_kmers_loaded),
                static_cast<unsigned long long>(stats.unique_kmers_assigned),
                static_cast<unsigned long long>(
                    stats.unique_kmers_outside_blocks));
            return EXIT_SUCCESS;
        }

        if (options.index_path.empty() || options.fai_path.empty() ||
            options.output_path.empty()) {
            std::fprintf(stderr,
                         "Usage: %s -i INDEX -f FASTA.fai -o BLOCKS.bed "
                         "[-t THREADS]\n",
                         argv[0]);
            throw std::runtime_error("--index, --fai, and --output are required");
        }

        const kmer_blocking::BlockingStats stats =
            kmer_blocking::create_blocks(options);
        std::fprintf(
            stderr,
            "Loaded %llu indexed 31-mers (%llu locations)\n"
            "Processed %llu contigs (%llu bases)\n"
            "Initial blocks: %llu\n"
            "Detected flank hotspots: %llu\n"
            "Moved breakpoints: %llu\n"
            "Absorbed blocks: %llu\n"
            "Final BED blocks: %llu\n",
            static_cast<unsigned long long>(stats.indexed_kmers),
            static_cast<unsigned long long>(stats.indexed_locations),
            static_cast<unsigned long long>(stats.contigs),
            static_cast<unsigned long long>(stats.bases),
            static_cast<unsigned long long>(stats.initial_blocks),
            static_cast<unsigned long long>(stats.hotspots),
            static_cast<unsigned long long>(stats.moved_breakpoints),
            static_cast<unsigned long long>(stats.absorbed_blocks),
            static_cast<unsigned long long>(stats.final_blocks));
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "kmerblocking: error: %s\n", error.what());
        return EXIT_FAILURE;
    }
}
