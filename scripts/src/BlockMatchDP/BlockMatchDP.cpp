// BlockMatchDP: KmerNorm uint16 matches between cells of the same graph block,
// then a glocal Needleman-Wunsch alignment of two ordered cell sequences.
//
// Cells are FASTA records named ...|>NODE or ...|<NODE (one graph-path step, as
// written by tools/alignment_cells_fasta.py); their group ("block kind") is
// NODE. Graph blocks given with -g (records graph|NODE or NODE) join their own
// group. Query (-q) and target (-t) cells are kept in file order.
//
// Step 1, match (KmerNorm calculatenorm, the uint16 mode, ifweight=1), over all
// selected cells of -q, -t and -g together, using canonical 31-mers:
//   * a k-mer present once in exactly one cell adds 0.05 to that cell's
//     singleton weight; (uint)(weight + 0.5) is added to its self-match;
//   * every other k-mer adds n_i * n_j to each pair of cells containing it
//     (KmerNorm's +1 and offsite paths are the n = 1 case).
// Values are summed exactly and saturated at 65535 (KmerNorm's uint16 matrix
// wraps on its += path). Only pairs of cells in the same group are stored.
// --max-occurrence N drops k-mers seen >= N times in the query or in the
// target (--occurrence-source target: in the target only, so the k-mers used do
// not depend on how often the query repeats them); --top-cells N keeps, per sequence, the N cells with the most kept
// k-mer positions (order is kept).
//
// Step 2, glocal NW: the query is aligned end to end; the target has free
// leading and trailing gaps. Pairing two cells scores their match (0 when
// their groups differ; --pair-score log uses -ln(1 - similarity) times the
// smaller self-match, with similarity capped at 0.999). Unpaired cells use affine costs: the first cell of a
// run costs --gap-open-factor (default: equal to --gap-factor) times its
// self-match, each further cell --gap-factor (default 1) times its self-match,
// rounded to the nearest integer. --gap-open-mean F adds F times the mean
// self-match of the selected query and target cells to every run's first cell,
// so the unpaired cells of an expansion are kept together in few runs.
//
// Outputs (PREFIX):
//   .cells.tsv   every selected cell: source, name, group, self-match
//   .groups.tsv  group sizes and offsets into .pairs.u16
//   .pairs.u16   per group, the upper triangle (row-major) of uint16 matches
//   .dp.tsv      the alignment, one column per row
//   .dp.summary  score and counts
//   .best.tsv    per query cell, its best target cell anywhere in the same group
//                by normalized match m/sqrt(self_q*self_t), and the runner-up
//
// Build: c++ -O3 -std=c++17 -pthread -o BlockMatchDP BlockMatchDP.cpp -lz

#include <zlib.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

constexpr int K = 31;
constexpr std::uint64_t MASK = (std::uint64_t{1} << (2 * K)) - 1;
constexpr int REVERSE_SHIFT = 2 * (K - 1);
constexpr double SINGLETON_WEIGHT = 0.05;
constexpr std::uint64_t MAX_U16 = 65535;

struct KmerHash {
    std::size_t operator()(std::uint64_t x) const {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return static_cast<std::size_t>(x ^ (x >> 31));
    }
};

[[noreturn]] void fail(const std::string& message) {
    std::cerr << "BlockMatchDP: " << message << "\n";
    std::exit(1);
}

struct Options {
    std::string query, target, graph, output;
    std::string pair_score = "match";  // match | log
    std::uint64_t max_occurrence = 0;
    std::string occurrence_source = "both";  // both | target
    std::size_t top_cells = 0;
    int threads = 1;
    double gap_factor = 1.0;
    double gap_open_factor = -1.0;  // default: equal to gap_factor (linear gaps)
    double gap_open_mean = 0.0;     // extra cost per gap run, in mean self-matches
};

void usage() {
    std::cerr <<
        "Usage: BlockMatchDP -q query.cells.fa -t target.cells.fa -o PREFIX [-g graph.FA]\n"
        "                    [--max-occurrence N] [--occurrence-source both|target]\n"
        "                    [--top-cells N] [--gap-factor F]\n"
        "                    [--gap-open-factor F] [--gap-open-mean F]\n"
        "                    [--pair-score match|log] [-T threads]\n";
}

Options parse(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc) fail("missing value for " + a);
            return argv[++i];
        };
        if (a == "-q") o.query = value();
        else if (a == "-t") o.target = value();
        else if (a == "-g") o.graph = value();
        else if (a == "-o") o.output = value();
        else if (a == "-T" || a == "--threads") o.threads = std::max(1, std::atoi(value().c_str()));
        else if (a == "--max-occurrence") o.max_occurrence = std::strtoull(value().c_str(), nullptr, 10);
        else if (a == "--occurrence-source") o.occurrence_source = value();
        else if (a == "--top-cells") o.top_cells = std::strtoull(value().c_str(), nullptr, 10);
        else if (a == "--gap-factor") o.gap_factor = std::strtod(value().c_str(), nullptr);
        else if (a == "--gap-open-factor") o.gap_open_factor = std::strtod(value().c_str(), nullptr);
        else if (a == "--gap-open-mean") o.gap_open_mean = std::strtod(value().c_str(), nullptr);
        else if (a == "--pair-score") o.pair_score = value();
        else if (a == "-h" || a == "--help") { usage(); std::exit(0); }
        else fail("unknown option " + a);
    }
    if (o.query.empty() || o.target.empty() || o.output.empty()) { usage(); std::exit(1); }
    if (!(o.gap_factor >= 0)) fail("--gap-factor must be >= 0");
    if (o.gap_open_factor < 0) o.gap_open_factor = o.gap_factor;
    if (!(o.gap_open_mean >= 0)) fail("--gap-open-mean must be >= 0");
    if (o.pair_score != "match" && o.pair_score != "log") fail("--pair-score must be match or log");
    if (o.occurrence_source != "both" && o.occurrence_source != "target")
        fail("--occurrence-source must be both or target");
    return o;
}

inline int base_code(char base) {
    if (base >= 'a') base -= 32;
    switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

struct Cell {
    std::string name, group;
    int source = 0;  // 0 query, 1 target, 2 graph
    std::uint64_t positions = 0, kept = 0;
};

std::string group_of(const std::string& name, int source) {
    if (source == 2) {
        return name.rfind("graph|", 0) == 0 ? name.substr(6) : name;
    }
    auto bar = name.rfind('|');
    if (bar == std::string::npos || bar + 2 > name.size())
        fail("cell name has no |>NODE or |<NODE suffix: " + name);
    std::string node = name.substr(bar + 1);
    if (node[0] == '>' || node[0] == '<') node = node.substr(1);
    return node;
}

}  // namespace

int main(int argc, char** argv) {
    const Options opt = parse(argc, argv);
    std::unordered_map<std::uint64_t, std::uint32_t, KmerHash> dictionary;
    dictionary.reserve(1 << 22);
    std::vector<Cell> cells;
    struct Entry { std::uint32_t kmer, cell, count; };
    std::vector<Entry> entries;

    // Read every cell; k-mer ids follow first-seen order.
    std::vector<std::uint32_t> current;
    auto flush = [&]() {
        if (cells.empty()) return;
        std::sort(current.begin(), current.end());
        const auto cell = static_cast<std::uint32_t>(cells.size() - 1);
        for (std::size_t a = 0; a < current.size();) {
            std::size_t b = a;
            while (b < current.size() && current[b] == current[a]) ++b;
            entries.push_back({current[a], cell, static_cast<std::uint32_t>(b - a)});
            a = b;
        }
        current.clear();
    };
    const std::string inputs[3] = {opt.query, opt.target, opt.graph};
    for (int source = 0; source < 3; ++source) {
        if (inputs[source].empty()) continue;
        gzFile file = gzopen(inputs[source].c_str(), "rb");
        if (!file) fail("cannot open " + inputs[source]);
        gzbuffer(file, 1 << 20);
        std::string line;
        std::vector<char> buffer(1 << 16);
        std::uint64_t forward = 0, reverse = 0;
        int size = 0;
        auto next = [&]() {
            line.clear();
            while (gzgets(file, buffer.data(), static_cast<int>(buffer.size()))) {
                line += buffer.data();
                if (!line.empty() && line.back() == '\n') return true;
            }
            return !line.empty();
        };
        bool in_record = false;
        while (next()) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                flush();
                std::size_t end = line.find_first_of(" \t\r\n", 1);
                Cell cell;
                cell.name = line.substr(1, end == std::string::npos ? std::string::npos : end - 1);
                cell.source = source;
                cell.group = group_of(cell.name, source);
                cells.push_back(cell);
                in_record = true;
                size = 0;
                continue;
            }
            if (!in_record) continue;
            for (char base : line) {
                if (base == '\n' || base == ' ' || base == '\r') continue;
                int code = base_code(base);
                if (code < 0) { size = 0; forward = reverse = 0; continue; }
                forward = ((forward << 2) | static_cast<std::uint64_t>(code)) & MASK;
                reverse = (reverse >> 2) | (static_cast<std::uint64_t>(3 - code) << REVERSE_SHIFT);
                if (++size < K) continue;
                const std::uint64_t kmer = forward >= reverse ? forward : reverse;
                auto found = dictionary.emplace(kmer, static_cast<std::uint32_t>(dictionary.size()));
                current.push_back(found.first->second);
                ++cells.back().positions;
            }
        }
        flush();
        in_record = false;
        gzclose(file);
    }
    const std::size_t kmer_count = dictionary.size();
    decltype(dictionary)().swap(dictionary);
    std::cerr << "[BlockMatchDP] cells read: " << cells.size() << "; 31-mers: " << kmer_count << "\n";

    // K-mer filter: occurrence within the query and within the target.
    std::vector<std::uint64_t> occurrence(2 * kmer_count, 0);
    for (const auto& e : entries) {
        const int source = cells[e.cell].source;
        if (source < 2) occurrence[2 * e.kmer + source] += e.count;
    }
    std::vector<char> keep(kmer_count, 1);
    if (opt.max_occurrence) {
        const bool target_only = opt.occurrence_source == "target";
        for (std::size_t k = 0; k < kmer_count; ++k)
            keep[k] = (target_only || occurrence[2 * k] < opt.max_occurrence) &&
                      occurrence[2 * k + 1] < opt.max_occurrence;
    }
    std::vector<std::uint64_t>().swap(occurrence);
    for (const auto& e : entries) if (keep[e.kmer]) cells[e.cell].kept += e.count;

    // Cell selection per sequence (graph blocks are always kept).
    std::vector<char> selected(cells.size(), 1);
    for (int source = 0; source < 2; ++source) {
        std::vector<std::uint32_t> members;
        for (std::uint32_t c = 0; c < cells.size(); ++c) if (cells[c].source == source) members.push_back(c);
        if (opt.top_cells && opt.top_cells < members.size()) {
            std::vector<std::uint32_t> ranked = members;
            std::stable_sort(ranked.begin(), ranked.end(), [&](std::uint32_t a, std::uint32_t b) {
                return cells[a].kept > cells[b].kept;
            });
            for (std::size_t r = opt.top_cells; r < ranked.size(); ++r) selected[ranked[r]] = 0;
        }
    }

    // Groups and member positions.
    std::unordered_map<std::string, std::uint32_t> group_index;
    std::vector<std::string> group_names;
    std::vector<std::vector<std::uint32_t>> group_members;
    std::vector<std::int64_t> group_of_cell(cells.size(), -1), position(cells.size(), -1);
    for (std::uint32_t c = 0; c < cells.size(); ++c) {
        if (!selected[c]) continue;
        auto found = group_index.emplace(cells[c].group, static_cast<std::uint32_t>(group_names.size()));
        if (found.second) { group_names.push_back(cells[c].group); group_members.emplace_back(); }
        const auto g = found.first->second;
        group_of_cell[c] = g;
        position[c] = static_cast<std::int64_t>(group_members[g].size());
        group_members[g].push_back(c);
    }
    std::size_t total_selected = 0;
    for (auto s : selected) total_selected += s;

    // Per k-mer lists of (group, position, count), grouped by group.
    struct Posting { std::uint32_t group, pos; std::uint16_t count; };
    std::vector<std::uint64_t> offset(kmer_count + 1, 0);
    for (const auto& e : entries) if (keep[e.kmer] && selected[e.cell]) ++offset[e.kmer + 1];
    for (std::size_t k = 0; k < kmer_count; ++k) offset[k + 1] += offset[k];
    std::vector<Posting> postings(offset[kmer_count]);
    {
        std::vector<std::uint64_t> fill(offset.begin(), offset.end() - 1);
        for (const auto& e : entries) {
            if (!keep[e.kmer] || !selected[e.cell]) continue;
            postings[fill[e.kmer]++] = {static_cast<std::uint32_t>(group_of_cell[e.cell]),
                                        static_cast<std::uint32_t>(position[e.cell]),
                                        static_cast<std::uint16_t>(std::min<std::uint32_t>(e.count, 65535))};
        }
    }
    std::vector<Entry>().swap(entries);
    for (std::size_t k = 0; k < kmer_count; ++k)
        std::sort(postings.begin() + offset[k], postings.begin() + offset[k + 1],
                  [](const Posting& a, const Posting& b) {
                      return a.group != b.group ? a.group < b.group : a.pos < b.pos;
                  });

    // Singleton weights (KmerNorm order: k-mer id order, float += 0.05).
    std::vector<float> singleton(cells.size(), 0.0F);
    std::vector<char> is_singleton(kmer_count, 0);
    for (std::size_t k = 0; k < kmer_count; ++k) {
        if (offset[k + 1] - offset[k] == 1 && postings[offset[k]].count == 1) {
            is_singleton[k] = 1;
            const auto& p = postings[offset[k]];
            singleton[group_members[p.group][p.pos]] += SINGLETON_WEIGHT;
        }
    }

    // Exact match sums per group (uint64), split by row ownership across threads.
    std::vector<std::vector<std::uint64_t>> sums(group_members.size());
    for (std::size_t g = 0; g < group_members.size(); ++g) {
        const std::uint64_t m = group_members[g].size();
        sums[g].assign(m * (m + 1) / 2, 0);
    }
    auto row_base = [](std::uint64_t m, std::uint64_t i) { return i * m - i * (i + 1) / 2; };
    const int threads = opt.threads;
    auto work = [&](int t) {
        for (std::size_t k = 0; k < kmer_count; ++k) {
            if (is_singleton[k]) continue;
            const auto begin = offset[k], end = offset[k + 1];
            for (auto a = begin; a < end;) {
                auto b = a;
                while (b < end && postings[b].group == postings[a].group) ++b;
                const auto g = postings[a].group;
                const std::uint64_t m = group_members[g].size();
                auto& sum = sums[g];
                for (auto x = a; x < b; ++x) {
                    const std::uint32_t i = postings[x].pos;
                    if (static_cast<int>(i % threads) != t) continue;
                    const std::uint64_t base = row_base(m, i);
                    const std::uint64_t ni = postings[x].count;
                    for (auto y = x; y < b; ++y) sum[base + postings[y].pos] += ni * postings[y].count;
                }
                a = b;
            }
        }
    };
    {
        std::vector<std::thread> pool;
        for (int t = 0; t < threads; ++t) pool.emplace_back(work, t);
        for (auto& thread : pool) thread.join();
    }

    // Saturate to uint16; add singleton weights to the self-matches.
    std::vector<std::vector<std::uint16_t>> match(group_members.size());
    std::vector<std::uint64_t> self(cells.size(), 0);
    std::uint64_t saturated = 0;
    for (std::size_t g = 0; g < group_members.size(); ++g) {
        const std::uint64_t m = group_members[g].size();
        match[g].resize(sums[g].size());
        for (std::uint64_t i = 0; i < m; ++i) {
            const auto c = group_members[g][i];
            sums[g][row_base(m, i) + i] += static_cast<unsigned>(singleton[c] + 0.5);
            for (std::uint64_t j = i; j < m; ++j) {
                const auto v = sums[g][row_base(m, i) + j];
                if (v > MAX_U16) ++saturated;
                match[g][row_base(m, i) + j] = static_cast<std::uint16_t>(std::min(v, MAX_U16));
            }
            self[c] = match[g][row_base(m, i) + i];
        }
        std::vector<std::uint64_t>().swap(sums[g]);
    }
    auto score = [&](std::uint32_t a, std::uint32_t b) -> std::uint64_t {
        const auto g = group_of_cell[a];
        if (g != group_of_cell[b]) return 0;
        const std::uint64_t m = group_members[g].size();
        std::uint64_t i = position[a], j = position[b];
        if (i > j) std::swap(i, j);
        return match[g][row_base(m, i) + j];
    };
    std::cerr << "[BlockMatchDP] selected cells: " << total_selected << "; groups: " << group_members.size()
              << "; values saturated at 65535: " << saturated << "\n";

    // Step 1 outputs.
    {
        std::FILE* list = std::fopen((opt.output + ".cells.tsv").c_str(), "w");
        if (!list) fail("cannot write " + opt.output + ".cells.tsv");
        std::fprintf(list, "cell\tsource\tname\tgroup\tgroup_position\tselected\tkmer_positions\tkept_positions\tself_match\n");
        const char* names[3] = {"query", "target", "graph"};
        for (std::size_t c = 0; c < cells.size(); ++c)
            std::fprintf(list, "%zu\t%s\t%s\t%s\t%lld\t%d\t%llu\t%llu\t%llu\n", c, names[cells[c].source],
                         cells[c].name.c_str(), cells[c].group.c_str(), static_cast<long long>(position[c]),
                         static_cast<int>(selected[c]), static_cast<unsigned long long>(cells[c].positions),
                         static_cast<unsigned long long>(cells[c].kept), static_cast<unsigned long long>(self[c]));
        std::fclose(list);
        std::FILE* groups = std::fopen((opt.output + ".groups.tsv").c_str(), "w");
        std::FILE* pairs = std::fopen((opt.output + ".pairs.u16").c_str(), "wb");
        if (!groups || !pairs) fail("cannot write group outputs");
        std::fprintf(groups, "group\tmembers\tquery\ttarget\tgraph\toffset_values\tvalues\n");
        std::uint64_t at = 0;
        for (std::size_t g = 0; g < group_members.size(); ++g) {
            int counts[3] = {0, 0, 0};
            for (auto c : group_members[g]) ++counts[cells[c].source];
            std::fprintf(groups, "%s\t%zu\t%d\t%d\t%d\t%llu\t%zu\n", group_names[g].c_str(),
                         group_members[g].size(), counts[0], counts[1], counts[2],
                         static_cast<unsigned long long>(at), match[g].size());
            if (std::fwrite(match[g].data(), sizeof(std::uint16_t), match[g].size(), pairs) != match[g].size())
                fail("write failed: " + opt.output + ".pairs.u16");
            at += match[g].size();
        }
        std::fclose(groups);
        std::fclose(pairs);
    }

    std::vector<std::uint32_t> query, target;
    for (std::uint32_t c = 0; c < cells.size(); ++c) {
        if (!selected[c]) continue;
        if (cells[c].source == 0) query.push_back(c);
        else if (cells[c].source == 1) target.push_back(c);
    }

    // Best target cell per query cell (any order), for non-collinear placement.
    {
        std::vector<std::int64_t> target_rank(cells.size(), -1);
        for (std::size_t r = 0; r < target.size(); ++r) target_rank[target[r]] = static_cast<std::int64_t>(r);
        std::FILE* best = std::fopen((opt.output + ".best.tsv").c_str(), "w");
        if (!best) fail("cannot write " + opt.output + ".best.tsv");
        std::fprintf(best, "query_cell\tquery_name\tbest_target_cell\tbest_target_name\ttarget_rank\t"
                           "match\tsimilarity\tsecond_similarity\tcandidates\n");
        for (const auto a : query) {
            const auto g = group_of_cell[a];
            double first = 0.0, second = 0.0;
            std::int64_t best_cell = -1;
            std::uint64_t best_match = 0, candidates = 0;
            for (const auto b : group_members[g]) {
                if (cells[b].source != 1) continue;
                ++candidates;
                const std::uint64_t m = score(a, b);
                const double sim = static_cast<double>(m) /
                    std::sqrt(std::max(1.0, static_cast<double>(self[a])) * std::max(1.0, static_cast<double>(self[b])));
                if (sim > first) { second = first; first = sim; best_cell = b; best_match = m; }
                else if (sim > second) second = sim;
            }
            if (best_cell < 0) {
                std::fprintf(best, "%u\t%s\t.\t.\t.\t0\t0\t0\t%llu\n", a, cells[a].name.c_str(),
                             static_cast<unsigned long long>(candidates));
            } else {
                std::fprintf(best, "%u\t%s\t%lld\t%s\t%lld\t%llu\t%.4f\t%.4f\t%llu\n", a, cells[a].name.c_str(),
                             static_cast<long long>(best_cell), cells[best_cell].name.c_str(),
                             static_cast<long long>(target_rank[best_cell]),
                             static_cast<unsigned long long>(best_match), first, second,
                             static_cast<unsigned long long>(candidates));
            }
        }
        std::fclose(best);
    }

    // Step 2: glocal NW (query global, target free ends).
    const std::size_t nq = query.size(), nt = target.size();
    std::cerr << "[BlockMatchDP] DP: " << nq << " x " << nt << "\n";
    // Pair score. match: the raw uint16 match. log: with similarity
    // s = min(match / sqrt(self_a * self_b), 0.999) and the smaller
    // self-match m, -ln(1 - s) * m (rounded); 0 when match is 0.
    auto pair_score = [&](std::uint32_t a, std::uint32_t b) -> std::int64_t {
        const std::uint64_t m = score(a, b);
        if (opt.pair_score == "match" || m == 0) return static_cast<std::int64_t>(m);
        const double sa = static_cast<double>(self[a]), sb = static_cast<double>(self[b]);
        const double similarity = static_cast<double>(m) / std::sqrt(std::max(1.0, sa) * std::max(1.0, sb));
        const double smaller = std::min(sa, sb);
        return static_cast<std::int64_t>(std::llround(-std::log(1.0 - std::min(similarity, 0.999)) * smaller));
    };
    // Affine (Gotoh) gaps: the first unpaired cell of a run costs
    // gap_open_factor x self-match, every further cell gap_factor x self-match.
    // With equal factors this is the linear scheme.
    auto scaled = [&](double factor, std::uint32_t c) -> std::int64_t {
        return static_cast<std::int64_t>(std::llround(factor * static_cast<double>(self[c])));
    };
    // --gap-open-mean: a fixed extra cost per gap run, F x the mean self-match.
    double self_sum = 0;
    for (auto c : query) self_sum += static_cast<double>(self[c]);
    for (auto c : target) self_sum += static_cast<double>(self[c]);
    const std::int64_t open_extra = (opt.gap_open_mean > 0 && (!query.empty() || !target.empty()))
        ? static_cast<std::int64_t>(std::llround(opt.gap_open_mean * self_sum /
                                                 static_cast<double>(query.size() + target.size())))
        : 0;
    if (open_extra) std::cerr << "[BlockMatchDP] gap open: +" << open_extra << " per gap run\n";
    auto open_cost = [&](std::uint32_t c) { return scaled(opt.gap_open_factor, c) + open_extra; };
    auto extend_cost = [&](std::uint32_t c) { return scaled(opt.gap_factor, c); };
    constexpr std::int64_t NEG = std::numeric_limits<std::int64_t>::min() / 4;
    enum State : std::uint8_t { PAIR = 0, QGAP = 1, TGAP = 2 };
    // Per cell: predecessor state of PAIR (bits 0-1), QGAP (bits 2-3), TGAP (bits 4-5).
    std::vector<std::uint8_t> trace((nq + 1) * (nt + 1), 0);
    std::vector<std::int64_t> pm(nt + 1, 0), px(nt + 1, NEG), py(nt + 1, NEG);  // row 0: free target start
    std::vector<std::int64_t> cm(nt + 1), cx(nt + 1), cy(nt + 1);
    auto argmax3 = [](std::int64_t m, std::int64_t x, std::int64_t y, std::uint8_t& who) {
        who = PAIR;
        std::int64_t best = m;
        if (x > best) { best = x; who = QGAP; }
        if (y > best) { best = y; who = TGAP; }
        return best;
    };
    for (std::size_t i = 1; i <= nq; ++i) {
        const auto a = query[i - 1];
        const std::int64_t open_a = open_cost(a), extend_a = extend_cost(a);
        cm[0] = NEG; cy[0] = NEG;
        cx[0] = (i == 1) ? pm[0] - open_a : px[0] - extend_a;
        trace[i * (nt + 1)] = static_cast<std::uint8_t>(((i == 1) ? PAIR : QGAP) << 2);
        for (std::size_t j = 1; j <= nt; ++j) {
            const auto b = target[j - 1];
            std::uint8_t from_m, from_x, from_y;
            const std::int64_t diag = argmax3(pm[j - 1], px[j - 1], py[j - 1], from_m);
            cm[j] = diag + pair_score(a, b);
            cx[j] = argmax3(pm[j] - open_a, px[j] - extend_a, py[j] - open_a, from_x);
            const std::int64_t open_b = open_cost(b);
            cy[j] = argmax3(cm[j - 1] - open_b, cx[j - 1] - open_b, cy[j - 1] - extend_cost(b), from_y);
            // argmax3 order is (pair, query gap, target gap); for TGAP the middle
            // argument is the QGAP predecessor and the last the TGAP one.
            trace[i * (nt + 1) + j] = static_cast<std::uint8_t>(from_m | (from_x << 2) | (from_y << 4));
        }
        std::swap(pm, cm); std::swap(px, cx); std::swap(py, cy);
    }
    std::size_t end_j = 0;
    std::uint8_t end_state = PAIR;
    std::int64_t total = NEG;
    for (std::size_t j = 0; j <= nt; ++j) {  // free trailing target
        std::uint8_t who;
        const std::int64_t v = argmax3(pm[j], px[j], py[j], who);
        if (v > total) { total = v; end_j = j; end_state = who; }
    }

    struct Column { int op; std::int64_t a, b; };
    std::vector<Column> columns;
    std::size_t i = nq, j = end_j;
    std::uint8_t state = end_state;
    while (i > 0) {
        const std::uint8_t bits = trace[i * (nt + 1) + j];
        if (state == PAIR) {
            columns.push_back({0, static_cast<std::int64_t>(query[i - 1]), static_cast<std::int64_t>(target[j - 1])});
            state = bits & 3; --i; --j;
        } else if (state == QGAP) {
            columns.push_back({1, static_cast<std::int64_t>(query[i - 1]), -1});
            state = (bits >> 2) & 3; --i;
        } else {
            columns.push_back({2, -1, static_cast<std::int64_t>(target[j - 1])});
            state = (bits >> 4) & 3; --j;
        }
    }
    std::reverse(columns.begin(), columns.end());
    const std::size_t start_j = j;

    std::FILE* dp = std::fopen((opt.output + ".dp.tsv").c_str(), "w");
    if (!dp) fail("cannot write " + opt.output + ".dp.tsv");
    std::fprintf(dp, "op\tquery_cell\tquery_name\ttarget_cell\ttarget_name\tsame_group\tmatch\tgap_cost\tpair_score\n");
    std::uint64_t paired = 0, same_group = 0, query_gaps = 0, target_gaps = 0, identical = 0;
    int previous_op = 0;
    for (const auto& col : columns) {
        const bool opens = col.op != previous_op;
        previous_op = col.op;
        if (col.op == 0) {
            const auto a = static_cast<std::uint32_t>(col.a), b = static_cast<std::uint32_t>(col.b);
            const bool same = group_of_cell[a] == group_of_cell[b];
            ++paired; same_group += same; identical += cells[a].name == cells[b].name;
            std::fprintf(dp, "pair\t%u\t%s\t%u\t%s\t%d\t%llu\t0\t%lld\n", a, cells[a].name.c_str(), b,
                         cells[b].name.c_str(), same ? 1 : 0, static_cast<unsigned long long>(score(a, b)),
                         static_cast<long long>(pair_score(a, b)));
        } else if (col.op == 1) {
            const auto a = static_cast<std::uint32_t>(col.a);
            ++query_gaps;
            std::fprintf(dp, "query_gap\t%u\t%s\t.\t.\t.\t0\t%lld\t0\n", a, cells[a].name.c_str(),
                         static_cast<long long>(opens ? open_cost(a) : extend_cost(a)));
        } else {
            const auto b = static_cast<std::uint32_t>(col.b);
            ++target_gaps;
            std::fprintf(dp, "target_gap\t.\t.\t%u\t%s\t.\t0\t%lld\t0\n", b, cells[b].name.c_str(),
                         static_cast<long long>(opens ? open_cost(b) : extend_cost(b)));
        }
    }
    std::fclose(dp);
    std::FILE* summary = std::fopen((opt.output + ".dp.summary").c_str(), "w");
    if (!summary) fail("cannot write " + opt.output + ".dp.summary");
    std::fprintf(summary,
                 "score\t%lld\npair_score\t%s\ngap_factor\t%g\ngap_open_factor\t%g\nquery_cells\t%zu\ntarget_cells\t%zu\ntarget_span\t%zu-%zu\n"
                 "paired\t%llu\npaired_same_group\t%llu\npaired_same_name\t%llu\n"
                 "query_unpaired\t%llu\ntarget_unpaired_inside_span\t%llu\nsaturated_values\t%llu\n",
                 static_cast<long long>(total), opt.pair_score.c_str(), opt.gap_factor, opt.gap_open_factor, nq, nt, start_j, end_j,
                 static_cast<unsigned long long>(paired), static_cast<unsigned long long>(same_group),
                 static_cast<unsigned long long>(identical), static_cast<unsigned long long>(query_gaps),
                 static_cast<unsigned long long>(target_gaps), static_cast<unsigned long long>(saturated));
    if (opt.occurrence_source != "both") std::fprintf(summary, "occurrence_source\t%s\n", opt.occurrence_source.c_str());
    if (open_extra) std::fprintf(summary, "gap_open_mean\t%g\ngap_open_extra\t%lld\n", opt.gap_open_mean,
                                 static_cast<long long>(open_extra));
    std::fclose(summary);
    std::cerr << "[BlockMatchDP] score " << total << "; paired " << paired << " (same group " << same_group
              << ", same name " << identical << "); query unpaired " << query_gaps << "\n";
    return 0;
}
