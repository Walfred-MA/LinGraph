// Native small-variant worker. Python owns CLI, headers, manifests and atomic
// publication; this process owns record parsing, sorting, grouping and output.
#include "support.hpp"
#include <regex>
#include <atomic>

using namespace sm;
static const Str SNP_FORMAT = "GT:TYPE:SIZE:BASE:ASSEMBLYCONTIG:QUERYCOORD:ALLELENAME:LABEL_H";
static const Str SV_FORMAT =
    "GT:TYPE:SIZE:EXTENDGRAPHCIGAR:ASSEMBLYCONTIG:QUERYCOORD:TEMPLATEOFFSET:LIFTCATEGORY:ALLELENAME:LABEL_H";
static const Str SV_LEGACY = "GT:TYPE:SIZE:EXTENDGRAPHCIGAR:ASSEMBLYCONTIG:QUERYCOORD:ALLELENAME:LABEL_H";
static const std::array<Str, 10> EVENT{
    "GT",         "TYPE",       "SIZE",    "EXTENDGRAPHCIGAR", "ASSEMBLYCONTIG",
    "QUERYCOORD", "ALLELENAME", "LABEL_H", "TEMPLATEOFFSET",   "LIFTCATEGORY"};
static const Str BASES = "ACGTRYSWKMBDHVN";
static bool snp_format(View s) {
    return s == "HSNP" || s == SNP_FORMAT || s == View(SNP_FORMAT).substr(0, SNP_FORMAT.rfind(":ALLELENAME"));
}
static bool sv_format(View s) {
    return s == "HSV" || s == SV_FORMAT || s == SV_LEGACY ||
           s == View(SV_FORMAT).substr(0, SV_FORMAT.rfind(":ALLELENAME")) ||
           s == View(SV_LEGACY).substr(0, SV_LEGACY.rfind(":ALLELENAME"));
}
static bool no_call(View s) { return s.empty() || s == "." || s == "./." || s == "0" || s == "0/0"; }
static bool signed_digits(View &s) {
    if (!s.empty() && (s.front() == '+' || s.front() == '-'))
        s.remove_prefix(1);
    size_t n = 0;
    while (n < s.size() && s[n] >= '0' && s[n] <= '9')
        ++n;
    s.remove_prefix(n);
    return n != 0;
}
static bool label_ok(View s) {
    auto coordinate_h = [](View &value) {
        if (!signed_digits(value) || value.empty() || value.front() != 'H')
            return false;
        value.remove_prefix(1);
        return true;
    };
    for (;;) {
        if (s.empty())
            return false;
        if (s.front() == '.')
            s.remove_prefix(1);
        else {
            if (!coordinate_h(s))
                return false;
            if (!s.empty() && s.front() != '&') {
                if (s.front() == '_')
                    s.remove_prefix(1);
                if (!coordinate_h(s))
                    return false;
            }
        }
        if (s.empty())
            return true;
        if (s.front() != '&')
            return false;
        s.remove_prefix(1);
    }
}
static bool offset_ok(View s) { return s == "." || (signed_digits(s) && s.empty()); }
using Allele = std::array<Str, 10>;
static Allele parse_legacy(Str s) {
    if (s.find(':') == Str::npos)
        s = unescape(replace(std::move(s), "@7C", "|"));
    auto fields = split(s, ':');
    require(fields.size() >= 7, "malformed observation");
    size_t suffix = 3, n = fields.size();
    if (n >= 10 && label_ok(fields[n - 3]) && offset_ok(fields[n - 2]) &&
        (fields.back() == "Pri" || fields.back() == "Dup" || fields.back() == "DivergeDup" ||
         fields.back() == "."))
        suffix = 6;
    else if (n >= 9 && label_ok(fields[n - 2]) && offset_ok(fields.back()))
        suffix = 5;
    else if (n >= 8 && label_ok(fields.back()))
        suffix = 4;
    Allele a;
    a.fill(".");
    for (size_t i = 0; i < 3; ++i)
        a[i] = fields[i];
    size_t end = n - suffix;
    a[3] = Str(fields[3]);
    for (size_t i = 4; i < end; ++i)
        a[3] += ':' + Str(fields[i]);
    for (size_t i = 0; i < suffix; ++i)
        a[4 + i] = fields[end + i];
    return a;
}
static Str allele_text(const Allele &a) { return join(std::vector<Str>(a.begin(), a.end()), ":"); }
static std::vector<Allele> observations(View field, View fmt, bool compact) {
    if (no_call(field.substr(0, field.find(':'))))
        return {};
    std::vector<Allele> result;
    if (fmt == "HSNP" || fmt == "HSV") {
        for (auto part : split(field, '|')) {
            if (no_call(part))
                continue;
            Str value(part);
            if (part.find(':') == View::npos)
                value = unescape(replace(std::move(value), "@7C", "|"));
            auto a = parse_legacy(value);
            if (compact)
                for (auto &x : a)
                    if (x != ".")
                        x = unescape(x);
            result.push_back(std::move(a));
        }
        return result;
    }
    auto names = split(fmt, ':'), cols = split(field, ':');
    require(cols.size() == names.size(), "malformed allele: FORMAT column count");
    std::vector<std::vector<View>> vectors;
    for (size_t i = 1; i < cols.size(); ++i)
        vectors.push_back(split(cols[i], ','));
    require(!vectors.empty(), "malformed allele");
    size_t count = vectors[0].size();
    for (auto &v : vectors)
        require(v.size() == count, "malformed allele: FORMAT observation count");
    for (size_t e = 0; e < count; ++e) {
        Allele a;
        a.fill(".");
        a[0] = cols[0];
        for (size_t i = 1; i < names.size(); ++i) {
            for (size_t j = 1; j < EVENT.size(); ++j)
                if (names[i] == EVENT[j] || (j == 3 && names[i] == "BASE")) {
                    a[j] = vectors[i - 1][e] == "." ? "." : unescape(vectors[i - 1][e]);
                    break;
                }
        }
        result.push_back(compact ? std::move(a) : parse_legacy(allele_text(a)));
    }
    return result;
}
static uint8_t bases(View ref, View alt) {
    Str r = upper(Str(ref)), a = upper(Str(alt));
    auto ri = BASES.find(r), ai = BASES.find(a);
    require(r.size() == 1 && a.size() == 1 && ri != Str::npos && ai != Str::npos,
            "SNP requires single IUPAC REF/ALT bases");
    return uint8_t((ri << 4) | ai);
}

// Output columns have equal-width reference/missing fields. Coverage updates
// patch this reusable byte string; each site replaces only its touched samples.
// Copying the dense VCF bytes is unavoidable, but no per-site sample scan is used.
struct Columns {
    Intern samples;
    std::vector<std::vector<U32>> slots;
    Str defaults, missing, reference;
    size_t width = 0;
    struct Interval {
        U64 start, end;
        U32 sample;
    };
    std::vector<std::vector<std::pair<U64, U64>>> coverage;
    std::vector<Interval> events;
    using End = std::pair<U64, U32>;
    std::multimap<U64, U32> active;
    std::vector<std::multimap<U64, U32>::iterator> where;
    std::vector<bool> in_active;
    size_t cursor = 0;
    std::vector<std::pair<U32, const Str *>> replacements, field_refs;
    Str row_buffer;
    Columns(const std::vector<Str> &names, const Str &fmt) {
        size_t n = split(fmt, ':').size();
        missing = ".";
        for (size_t i = 1; i < n; ++i)
            missing += ":.";
        reference = missing;
        reference[0] = '0';
        width = missing.size() + 1;
        for (U32 column = 0; column < names.size(); ++column) {
            U32 s = sample(names[column]);
            slots[s].push_back(column);
            defaults += '\t' + missing;
        }
    }
    U32 sample(const Str &name) {
        U32 s = samples.add(name);
        if (slots.size() <= s) {
            slots.resize(size_t(s) + 1);
            coverage.resize(size_t(s) + 1);
        }
        return s;
    }
    void covered(const Str &name, U64 start, U64 end) { coverage[sample(name)].emplace_back(start, end); }
    void prepare() {
        for (U32 s = 0; s < coverage.size(); ++s) {
            auto &values = coverage[s];
            std::sort(values.begin(), values.end());
            std::vector<std::pair<U64, U64>> merged;
            for (auto x : values)
                if (x.second > x.first) {
                    if (!merged.empty() && x.first <= merged.back().second)
                        merged.back().second = std::max(merged.back().second, x.second);
                    else
                        merged.push_back(x);
                }
            values = std::move(merged);
            for (auto x : values)
                events.push_back({x.first, x.second, s});
        }
        std::sort(events.begin(), events.end(), [](auto a, auto b) { return a.start < b.start; });
        where.resize(slots.size());
        in_active.assign(slots.size(), false);
    }
    void set(U32 s, char state) {
        for (auto c : slots[s])
            defaults[size_t(c) * width + 1] = state;
    }
    void advance(U64 start) {
        while (cursor < events.size() && events[cursor].start <= start) {
            auto e = events[cursor++];
            if (in_active[e.sample])
                active.erase(where[e.sample]);
            where[e.sample] = active.emplace(e.end, e.sample);
            in_active[e.sample] = true;
            set(e.sample, '0');
        }
        while (!active.empty() && active.begin()->first <= start) {
            U32 s = active.begin()->second;
            in_active[s] = false;
            set(s, '.');
            active.erase(active.begin());
        }
    }
    void write(std::ostream &out, U64 end, const std::vector<std::pair<U32, Str>> &overrides) {
        field_refs.clear();
        for (const auto &x : overrides)
            field_refs.emplace_back(x.first, &x.second);
        write_refs(out, end, field_refs);
    }
    void write_refs(std::ostream &out, U64 end, const std::vector<std::pair<U32, const Str *>> &overrides) {
        replacements.clear();
        // Rare boundary-crossing indels require a missing default even when
        // their first base is covered. Explicit calls still take precedence.
        for (auto it = active.begin(); it != active.end() && it->first < end; ++it)
            for (auto c : slots[it->second])
                replacements.emplace_back(c, &missing);
        for (const auto &x : overrides)
            for (auto c : slots[x.first])
                replacements.emplace_back(c, x.second);
        std::stable_sort(replacements.begin(), replacements.end(),
                         [](auto a, auto b) { return a.first < b.first; });
        size_t from = 0;
        row_buffer.clear();
        for (size_t i = 0; i < replacements.size();) {
            size_t j = i + 1;
            while (j < replacements.size() && replacements[j].first == replacements[i].first)
                ++j;
            auto [column, text] = replacements[j - 1];
            size_t offset = size_t(column) * width;
            row_buffer.append(defaults.data() + from, offset - from);
            row_buffer += '\t';
            row_buffer += *text;
            from = offset + width;
            i = j;
        }
        row_buffer.append(defaults.data() + from, defaults.size() - from);
        row_buffer += '\n';
        out.write(row_buffer.data(), row_buffer.size());
    }
};

struct Record {
    U32 pos, qpos, query, allele, label;
    uint8_t base;
};
static_assert(sizeof(Record) == 24, "unexpected record alignment");
static U64 key(const Record &r) { return (U64(r.pos) << 8) | r.base; }
static bool record_less(const Record &a, const Record &b) { return key(a) < key(b); }
static U32 le32(const unsigned char *p) {
    return U32(p[0]) | (U32(p[1]) << 8) | (U32(p[2]) << 16) | (U32(p[3]) << 24);
}
static Record unpack(const unsigned char *p) {
    return {le32(p), le32(p + 5), le32(p + 9), le32(p + 13), le32(p + 17), p[4]};
}
static void pack(std::ostream &out, const Record &r) {
    put32(out, r.pos);
    out.put(char(r.base));
    put32(out, r.qpos);
    put32(out, r.query);
    put32(out, r.allele);
    put32(out, r.label);
}
static void radix_sort(std::vector<Record> &a, U32 threads) {
    if (a.size() < 32768) {
        std::sort(a.begin(), a.end(), record_less);
        return;
    }
    threads = std::max<U32>(1, std::min<U32>(threads, U32(a.size() / 32768)));
    std::vector<Record> scratch(a.size());
    using Counts = std::array<size_t, 256>;
    std::vector<Counts> counts(threads);
    auto parallel = [&](auto fn) {
        std::vector<std::thread> workers;
        for (U32 t = 1; t < threads; ++t)
            workers.emplace_back(fn, t);
        fn(0);
        for (auto &w : workers)
            w.join();
    };
    for (unsigned shift = 0; shift < 40; shift += 8) {
        parallel([&](U32 t) {
            counts[t].fill(0);
            for (size_t i = a.size() * t / threads; i < a.size() * (t + 1) / threads; ++i)
                ++counts[t][(key(a[i]) >> shift) & 255];
        });
        size_t offset = 0;
        for (unsigned b = 0; b < 256; ++b)
            for (U32 t = 0; t < threads; ++t) {
                size_t n = counts[t][b];
                counts[t][b] = offset;
                offset += n;
            }
        parallel([&](U32 t) {
            auto positions = counts[t];
            for (size_t i = a.size() * t / threads; i < a.size() * (t + 1) / threads; ++i)
                scratch[positions[(key(a[i]) >> shift) & 255]++] = a[i];
        });
        a.swap(scratch);
    }
}

// Sorted spill runs have a bounded fan-in, so large cohorts never need one
// open descriptor per sample. All scratch files belong to this invocation.
template <class T, class Less, class Save, class Load> struct Runs {
    Less less;
    Save save;
    Load load;
    Str prefix;
    std::vector<Str> paths, owned;
    U64 serial = 0;
    Runs(Less l, Save s, Load r, Str p) : less(l), save(s), load(r), prefix(std::move(p)) {}
    ~Runs() {
        for (auto &p : owned) {
            std::error_code ec;
            fs::remove(p, ec);
        }
    }
    Str path() {
        Str p = prefix + ".run." + std::to_string(serial++);
        owned.push_back(p);
        return p;
    }
    void store(const std::vector<T> &a) {
        Str p = path();
        std::ofstream f(p, std::ios::binary);
        require(bool(f), "cannot create spill run");
        for (const auto &x : a)
            save(f, x);
        checked_close(f);
        paths.push_back(p);
    }
    template <class Emit> void merge(const std::vector<Str> &inputs, Emit emit) {
        std::vector<std::unique_ptr<Reader>> files;
        std::vector<T> rows(inputs.size());
        auto cmp = [&](U32 a, U32 b) { return less(rows[b], rows[a]); };
        std::priority_queue<U32, std::vector<U32>, decltype(cmp)> heap(cmp);
        for (U32 i = 0; i < inputs.size(); ++i) {
            files.emplace_back(new Reader(inputs[i]));
            if (load(*files.back(), rows[i]))
                heap.push(i);
        }
        while (!heap.empty()) {
            U32 i = heap.top();
            heap.pop();
            emit(rows[i]);
            if (load(*files[i], rows[i]))
                heap.push(i);
        }
    }
    template <class Emit> void finish(Emit emit) {
        while (paths.size() > 64) {
            std::vector<Str> next;
            for (size_t i = 0; i < paths.size(); i += 64) {
                std::vector<Str> batch(paths.begin() + i, paths.begin() + std::min(paths.size(), i + 64));
                Str p = path();
                std::ofstream out(p, std::ios::binary);
                merge(batch, [&](const T &row) { save(out, row); });
                checked_close(out);
                next.push_back(p);
                for (auto &old : batch)
                    fs::remove(old);
            }
            paths = std::move(next);
        }
        merge(paths, emit);
    }
};
static void save_record(std::ostream &out, const Record &r) { pack(out, r); }
static bool load_record(Reader &f, Record &r) {
    if (f.f.peek() == EOF)
        return false;
    unsigned char b[21];
    f.bytes(reinterpret_cast<char *>(b), 21);
    r = unpack(b);
    return true;
}

struct Query {
    U32 sample, contig, kind, filter;
    uint8_t mode;
    char state, strand;
};
struct Dictionaries {
    Columns &columns;
    Intern contigs, kinds, filters, alleles, groups;
    std::vector<Query> queries;
    explicit Dictionaries(Columns &c) : columns(c) {}
    U32 query(const std::vector<Str> &q) {
        require(q.size() == 7, "invalid SNP query dictionary");
        auto mode = coordinate(q[5]);
        require(mode <= 2 && q[2].size() == 1 && q[4].size() == 1, "invalid SNP query metadata");
        Query value{columns.sample(q[0]), contigs.add(q[1]), kinds.add(q[3]), filters.add(q[6]),
                    uint8_t(mode),        q[4][0],           q[2][0]};
        U32 id = index32(queries.size());
        queries.push_back(value);
        return id;
    }
    Str label(const Record &r) const {
        auto mode = queries[r.query].mode;
        if (mode == 2) {
            require(r.label < groups.values.size(), "invalid SNP LABEL_H dictionary index");
            return groups.values[r.label];
        }
        if (!mode)
            return ".";
        int64_t n = int64_t(r.label >> 1) ^ -int64_t(r.label & 1);
        return std::to_string(n) + 'H';
    }
};
struct SnpEmitter {
    Columns &columns;
    Dictionaries &dict;
    std::ostream &out;
    Str chrom, token;
    struct Observation {
        U32 query, allele;
        Str coordinate, label;
    };
    struct Sample {
        bool touched = false;
        char state = 0;
        std::vector<Observation> observations;
        Str field;
    };
    std::vector<Sample> samples;
    std::vector<U32> touched, filter_ids;
    std::vector<bool> seen_filter;
    std::vector<std::vector<Str>> filter_parts;
    std::vector<std::pair<U32, const Str *>> fields;
    Record first{};
    bool have_site = false;
    U64 count = 0;
    SnpEmitter(Columns &c, Dictionaries &d, std::ostream &o, Str name, Str t)
        : columns(c), dict(d), out(o), chrom(std::move(name)), token(std::move(t)),
          samples(c.samples.values.size()), seen_filter(d.filters.values.size(), false) {
        for (const auto &filter : d.filters.values) {
            std::vector<Str> parts;
            for (auto part : split(filter, ';'))
                parts.emplace_back(part);
            filter_parts.push_back(std::move(parts));
        }
    }
    auto order(const Observation &value) const {
        const auto &q = dict.queries[value.query];
        // GT, SIZE and ALT are constant within a site. These five remaining
        // fields reproduce Python's raw eight-field tuple order exactly.
        return std::tie(dict.kinds.values[q.kind], dict.contigs.values[q.contig], value.coordinate,
                        dict.alleles.values[value.allele], value.label);
    }
    void add(const Record &r) {
        if (have_site && key(first) != key(r))
            flush();
        if (!have_site) {
            first = r;
            have_site = true;
        }
        const auto &q = dict.queries[r.query];
        if (!seen_filter[q.filter]) {
            seen_filter[q.filter] = true;
            filter_ids.push_back(q.filter);
        }
        auto &sample = samples[q.sample];
        if (!sample.touched) {
            sample.touched = true;
            touched.push_back(q.sample);
        }
        if (q.state != '1') {
            if (q.state == '.' || !sample.state)
                sample.state = q.state;
            return;
        }
        Observation value{r.query, r.allele, std::to_string(r.qpos) + q.strand, dict.label(r)};
        // Adjacent duplicates are common in overlapping source shards. The
        // final sort/unique also handles duplicates separated in the input.
        if (sample.observations.empty() || order(sample.observations.back()) != order(value))
            sample.observations.push_back(std::move(value));
    }
    void flush() {
        if (!have_site)
            return;
        fields.clear();
        U32 supporters = 0;
        for (U32 id : touched) {
            auto &sample = samples[id];
            auto &values = sample.observations;
            if (values.empty()) {
                if (sample.state == '.' || sample.state == '0')
                    fields.emplace_back(id, sample.state == '0' ? &columns.reference : &columns.missing);
                continue;
            }
            ++supporters;
            std::sort(values.begin(), values.end(),
                      [&](const auto &a, const auto &b) { return order(a) < order(b); });
            values.erase(std::unique(values.begin(), values.end(),
                                     [&](const auto &a, const auto &b) { return order(a) == order(b); }),
                         values.end());
            Str &field = sample.field;
            field = "1";
            for (size_t column = 1; column < 8; ++column) {
                field += ':';
                for (size_t i = 0; i < values.size(); ++i) {
                    if (i)
                        field += ',';
                    const auto &v = values[i];
                    const auto &q = dict.queries[v.query];
                    switch (column) {
                    case 1:
                        append_escaped(field, dict.kinds.values[q.kind]);
                        break;
                    case 2:
                        field += '1';
                        break;
                    case 3:
                        field += BASES[first.base & 15];
                        break;
                    case 4:
                        append_escaped(field, dict.contigs.values[q.contig]);
                        break;
                    case 5:
                        append_escaped(field, v.coordinate);
                        break;
                    case 6:
                        append_escaped(field, dict.alleles.values[v.allele]);
                        break;
                    case 7:
                        append_escaped(field, v.label);
                        break;
                    }
                }
            }
            fields.emplace_back(id, &field);
        }
        if (supporters) {
            columns.advance(U64(first.pos) - 1);
            std::set<Str> filters;
            bool pass = false;
            for (U32 id : filter_ids) {
                const auto &parts = filter_parts[id];
                if (std::find(parts.begin(), parts.end(), "PASS") != parts.end()) {
                    pass = true;
                    break;
                }
                filters.insert(parts.begin(), parts.end());
            }
            Str filter = "PASS";
            if (!pass) {
                filters.erase(".");
                filters.erase("");
                filter = join(std::vector<Str>(filters.begin(), filters.end()), ";");
                if (filter.empty())
                    filter = ".";
            }
            out << chrom << '\t' << first.pos << "\tS_" << token << '_' << first.pos << '_'
                << BASES[first.base & 15] << '\t' << BASES[first.base >> 4] << '\t' << BASES[first.base & 15]
                << "\t.\t" << filter << "\tNSUP=" << supporters << '\t' << SNP_FORMAT;
            columns.write_refs(out, first.pos, fields);
            ++count;
        }
        // Reset only the touched samples; retained capacities are reused by
        // later sites. No O(cohort size) initialization occurs inside this loop.
        for (U32 id : touched) {
            auto &sample = samples[id];
            sample.touched = false;
            sample.state = 0;
            sample.observations.clear();
            // A rare high-copy site must not permanently enlarge every
            // sample's cache over the rest of the chromosome.
            if (sample.observations.capacity() > 64)
                std::vector<Observation>().swap(sample.observations);
            if (sample.field.capacity() > 4096)
                Str().swap(sample.field);
        }
        for (U32 id : filter_ids)
            seen_filter[id] = false;
        touched.clear();
        filter_ids.clear();
        have_site = false;
    }
};

using DropKey = std::tuple<U32, U32, U32, U32, uint8_t, char>;
struct DropHash {
    size_t operator()(const DropKey &k) const {
        auto [sample, contig, pos, qpos, alt, strand] = k;
        U64 h = sample;
        for (U64 v : {U64(contig), U64(pos), U64(qpos), U64(alt), U64(uint8_t(strand))})
            h = (h ^ (v + 0x9e3779b9 + (h << 6) + (h >> 2)));
        return size_t(h);
    }
};
static void compact(Reader &request, const Str &output) {
    Str chrom = request.str(), token = request.str();
    U64 budget = request.u64();
    U32 threads = request.u32();
    Columns columns(request.strings(), SNP_FORMAT);
    Dictionaries dict(columns);
    std::unordered_map<DropKey, bool, DropHash> drops;
    U32 drop_count = request.u32();
    while (drop_count--) {
        Str sample = request.str();
        U32 pos = request.u32();
        Str alt = upper(request.str()), contig = request.str();
        U32 qpos = request.u32();
        Str strand = request.str();
        require(alt.size() == 1 && strand.size() == 1, "invalid realignment drop");
        auto b = BASES.find(alt);
        require(b != Str::npos, "invalid realignment drop base");
        DropKey k{columns.sample(sample), dict.contigs.add(contig), pos, qpos, uint8_t(b), strand[0]};
        require(drops.emplace(k, false).second, "duplicate realignment drop not found after first removal");
    }
    struct Add {
        Str sample, ref, alt, contig, strand;
        U32 pos, qpos;
    };
    std::vector<Add> adds;
    U32 add_count = request.u32();
    while (add_count--) {
        Add a;
        a.sample = request.str();
        a.pos = request.u32();
        a.ref = request.str();
        a.alt = request.str();
        a.contig = request.str();
        a.qpos = request.u32();
        a.strand = request.str();
        adds.push_back(std::move(a));
    }
    U64 total = request.u64();
    // Reserve half the advertised budget for dictionaries, source metadata,
    // output and the parent Python process. Radix sorting uses two arrays.
    U64 capacity = std::max<U64>(1, std::min<U64>(total + adds.size(), budget / (4 * sizeof(Record))));
    std::vector<Record> records;
    records.reserve(size_t(capacity));
    Runs<Record, decltype(&record_less), decltype(&save_record), decltype(&load_record)> runs(
        record_less, save_record, load_record, output);
    auto flush = [&] {
        radix_sort(records, threads);
        runs.store(records);
        records.clear();
    };
    auto append = [&](Record r, bool realign) {
        require(r.pos > 0 && (r.base >> 4) < BASES.size() && (r.base & 15) < BASES.size(),
                "invalid packed SNP position/base");
        if (realign && !drops.empty()) {
            auto q = dict.queries[r.query];
            if (q.state == '1' && dict.kinds.values[q.kind] == "SNP") {
                auto it = drops.find({q.sample, q.contig, r.pos, r.qpos, uint8_t(r.base & 15), q.strand});
                if (it != drops.end()) {
                    it->second = true;
                    return;
                }
            }
        }
        if (records.size() == capacity)
            flush();
        records.push_back(r);
    };
    U32 source_count = request.u32();
    U64 loaded = 0;
    while (source_count--) {
        U32 nq = request.u32();
        std::vector<U32> qm;
        qm.reserve(nq);
        while (nq--)
            qm.push_back(dict.query(request.strings()));
        auto source_alleles = request.strings(), source_groups = request.strings();
        std::vector<U32> am, gm;
        for (auto &a : source_alleles)
            am.push_back(dict.alleles.add(a));
        for (auto &g : source_groups)
            gm.push_back(dict.groups.add(g));
        Str binary = request.str();
        U32 spans = request.u32();
        std::unique_ptr<Reader> data;
        if (spans)
            data = std::make_unique<Reader>(binary);
        std::vector<unsigned char> buffer(21 * 65536);
        while (spans--) {
            U64 offset = request.u64(), count = request.u64();
            data->f.seekg(offset);
            require(bool(data->f), "invalid SNP section offset");
            while (count) {
                size_t n = size_t(std::min<U64>(count, 65536));
                data->bytes(reinterpret_cast<char *>(buffer.data()), n * 21);
                for (size_t i = 0; i < n; ++i) {
                    auto r = unpack(buffer.data() + 21 * i);
                    require(r.query < qm.size() && r.allele < am.size(), "invalid SNP dictionary index");
                    r.query = qm[r.query];
                    r.allele = am[r.allele];
                    if (dict.queries[r.query].mode == 2) {
                        require(r.label < gm.size(), "invalid SNP LABEL_H dictionary index");
                        r.label = gm[r.label];
                    }
                    append(r, true);
                }
                count -= n;
                loaded += n;
            }
        }
        U32 nc = request.u32();
        while (nc--) {
            Str s = request.str();
            U64 start = request.u64(), end = request.u64();
            columns.covered(s, start, end);
        }
    }
    require(loaded == total, "inconsistent SNP observation count");
    for (const auto &d : drops)
        require(d.second, "realignment drop not found");
    for (const auto &a : adds) {
        U32 q = dict.query({a.sample, a.contig, a.strand, "SNP", "1", "0", "PASS"});
        append({a.pos, a.qpos, q, dict.alleles.add("."), 0, bases(a.ref, a.alt)}, false);
    }
    columns.prepare();
    std::ofstream out(output, std::ios::binary);
    require(bool(out), "cannot create SNP part");
    SnpEmitter emit(columns, dict, out, chrom, token);
    if (runs.paths.empty()) {
        radix_sort(records, threads);
        for (const auto &r : records)
            emit.add(r);
    } else {
        if (!records.empty())
            flush();
        std::vector<Record>().swap(records);
        runs.finish([&](const Record &r) { emit.add(r); });
    }
    emit.flush();
    checked_close(out);
}

struct SmallKey {
    bool snp;
    U32 pos;
    Str ref, alt;
    int64_t end;
    Str kind, sequence;
    auto tuple() const { return std::tie(snp, pos, ref, alt, end, kind, sequence); }
    bool operator<(const SmallKey &b) const { return tuple() < b.tuple(); }
    bool operator==(const SmallKey &b) const { return tuple() == b.tuple(); }
    Str json() const {
        Str s = "[" + std::to_string(pos) + ", " + quote(ref) + ", " + quote(alt);
        if (!snp)
            s += ", " + std::to_string(end) + ", " + quote(kind) + ", " + quote(sequence);
        return s + "]";
    }
};
struct SmallHash {
    size_t operator()(const SmallKey &k) const {
        size_t h = k.pos;
        auto mix = [&](size_t v) { h ^= v + 0x9e3779b9 + (h << 6) + (h >> 2); };
        mix(k.snp);
        mix(std::hash<Str>{}(k.ref));
        mix(std::hash<Str>{}(k.alt));
        mix(size_t(k.end));
        mix(std::hash<Str>{}(k.kind));
        mix(std::hash<Str>{}(k.sequence));
        return h;
    }
};
static bool small_key(const std::vector<View> &row, U32 cutoff, SmallKey &k) {
    require(row.size() >= 9, "VCF record has fewer than 9 columns");
    k.snp = snp_format(row[8]);
    k.pos = coordinate(row[1]);
    k.ref = row[3];
    k.alt = row[4];
    k.end = 0;
    k.kind.clear();
    k.sequence.clear();
    if (k.snp)
        return true;
    require(sv_format(row[8]), "unsupported graphreftovcf FORMAT");
    auto values = info(row[7]);
    int64_t size = integer(get(values, "MAXSIZE", get(values, "SVLEN", "0")));
    require(size != INT64_MIN, "invalid indel size");
    if (size < 0)
        size = -size;
    if (U64(size) >= cutoff)
        return false;
    require(size > 0, "indel needs positive MAXSIZE or nonzero SVLEN");
    k.end = integer(get(values, "END", std::to_string(k.pos)));
    Str fallback = k.alt;
    while (!fallback.empty() && (fallback.front() == '<' || fallback.front() == '>'))
        fallback.erase(0, 1);
    while (!fallback.empty() && (fallback.back() == '<' || fallback.back() == '>'))
        fallback.pop_back();
    k.kind = get(values, "SVTYPE", fallback);
    k.sequence = unescape(get(values, "SEQ", get(values, "SVINSSEQ")));
    if (k.kind == "INS" || k.kind == "SUB")
        require(!k.sequence.empty() &&
                    std::all_of(k.sequence.begin(), k.sequence.end(),
                                [](char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'); }),
                "small INS/SUB requires plain INFO/SEQ for exact allele matching");
    k.sequence = upper(std::move(k.sequence));
    return true;
}
static Str nsup_info(const Str &text, size_t nsup) {
    auto values = info(text);
    std::vector<Str> order;
    std::unordered_set<Str> seen;
    for (auto item : split(text, ';')) {
        Str name(item.substr(0, item.find('=')));
        if (name != "PAMAP" && name != "PAPATH" && !name.empty() && seen.insert(name).second)
            order.push_back(name);
    }
    values.erase("PAMAP");
    values.erase("PAPATH");
    values["NSUP"] = std::to_string(nsup);
    if (seen.insert("NSUP").second)
        order.push_back("NSUP");
    std::vector<Str> items;
    for (auto &name : order)
        if (values.count(name))
            items.push_back(name + (values[name].empty() ? "" : "=" + values[name]));
    return join(items, ";");
}
static void shift_offset(Allele &a) {
    if (a[8] != ".")
        return;
    View cigar = a[3];
    int64_t shift = 0;
    while (cigar.size() > 2 && (cigar.front() == '>' || cigar.front() == '<')) {
        size_t i = 1;
        while (i < cigar.size() && cigar[i] >= '0' && cigar[i] <= '9')
            ++i;
        if (i == 1 || i + 1 >= cigar.size() || cigar[i] != 'H' ||
            (cigar[i + 1] != '>' && cigar[i + 1] != '<'))
            break;
        int64_t n = integer(cigar.substr(1, i - 1));
        shift += cigar.front() == '>' ? n : -n;
        cigar.remove_prefix(i + 1);
    }
    a[3] = Str(cigar);
    a[8] = std::to_string(shift);
}
struct SmallRecord {
    U32 site, sample;
    char state;
    Str allele;
};
static void save_small(std::ostream &out, const SmallRecord &r) {
    put32(out, r.site);
    put32(out, r.sample);
    out.put(r.state);
    putstr(out, r.allele);
}
static bool load_small(Reader &f, SmallRecord &r) {
    if (f.f.peek() == EOF)
        return false;
    r.site = f.u32();
    r.sample = f.u32();
    f.bytes(&r.state, 1);
    r.allele = f.str();
    return true;
}
struct SmallSite {
    SmallKey key;
    Str qual, filter, info;
};
static void small_merge(Reader &request, const Str &snp_output, const Str &indel_output) {
    Str chrom = request.str(), token = request.str();
    U32 cutoff = request.u32();
    U64 budget = request.u64();
    auto names = request.strings();
    Columns snp_columns(names, SNP_FORMAT), indel_columns(names, SV_FORMAT);
    std::vector<SmallSite> sites;
    std::unordered_map<SmallKey, U32, SmallHash> index;
    std::vector<SmallRecord> records;
    U64 bytes = 0, site_bytes = 0;
    auto less = [&](const SmallRecord &a, const SmallRecord &b) {
        if (a.site != b.site)
            return sites[a.site].key < sites[b.site].key;
        return std::tie(a.sample, a.state, a.allele) < std::tie(b.sample, b.state, b.allele);
    };
    Runs<SmallRecord, decltype(less), decltype(&save_small), decltype(&load_small)> runs(
        less, save_small, load_small, snp_output);
    auto flush = [&] {
        std::sort(records.begin(), records.end(), less);
        runs.store(records);
        std::vector<SmallRecord>().swap(records);
        bytes = 0;
    };
    auto append = [&](SmallRecord r) {
        bytes += sizeof(SmallRecord) * 2 + r.allele.size() + 32;
        records.push_back(std::move(r));
        if (bytes > budget / 2)
            flush();
    };
    U32 files = request.u32();
    while (files--) {
        Str path = request.str();
        auto samples = request.strings();
        std::vector<U32> sample_ids;
        for (auto &name : samples) {
            U32 s = snp_columns.sample(name);
            require(indel_columns.sample(name) == s, "inconsistent sample map");
            sample_ids.push_back(s);
        }
        Lines input(path);
        Str line;
        U64 number = 0;
        while (input.next(line)) {
            ++number;
            try {
                if (line.rfind("##referenceCoverage=<", 0) == 0) {
                    auto v = structured(line, "referenceCoverage");
                    Str s = get(v, "Sample");
                    U64 start = coordinate(get(v, "Start")), end = coordinate(get(v, "End"));
                    snp_columns.covered(s, start, end);
                    indel_columns.covered(s, start, end);
                    continue;
                }
                auto row = split(line, '\t');
                SmallKey k;
                require(small_key(row, cutoff, k), "unexpected large variant in small shard");
                require(row.size() == 9 + samples.size(), "sample column count differs from header");
                auto it = index.find(k);
                U32 id;
                if (it == index.end()) {
                    id = index32(sites.size());
                    site_bytes += sizeof(SmallSite) + sizeof(SmallKey) + 256 +
                                  2 * (k.ref.size() + k.alt.size() + k.kind.size() + k.sequence.size()) +
                                  row[5].size() + row[6].size() + row[7].size();
                    require(site_bytes < budget / 4,
                            "small-variant site metadata exceeds RAM budget; increase GRAPHVCFMERGE_RAM_GB");
                    index.emplace(k, id);
                    sites.push_back({std::move(k), Str(row[5]), Str(row[6]), Str(row[7])});
                } else
                    id = it->second;
                if (row[6] == "PASS")
                    sites[id].filter = "PASS";
                for (size_t i = 0; i < samples.size(); ++i) {
                    auto field = row[9 + i];
                    auto gt = field.substr(0, field.find(':'));
                    auto values = observations(field, row[8], false);
                    if (values.empty()) {
                        require(no_call(gt), "malformed allele");
                        append({id, sample_ids[i], char(gt == "0" || gt == "0/0" ? '0' : '.'), {}});
                    }
                    for (auto &a : values) {
                        if (!sites[id].key.snp)
                            shift_offset(a);
                        append({id, sample_ids[i], '1', allele_text(a)});
                    }
                }
            } catch (const std::exception &e) {
                throw std::runtime_error(path + ":" + std::to_string(number) + ": " + e.what());
            }
        }
    }
    // Source sample coverage can include names absent from the output header.
    // The lookup table is needed only while reading. Free its duplicate keys
    // before sorting and formatting the observation records.
    index.clear();
    index.rehash(0);
    snp_columns.prepare();
    indel_columns.prepare();
    std::ofstream snp(snp_output, std::ios::binary), indel(indel_output, std::ios::binary);
    require(bool(snp) && bool(indel), "cannot create small-variant parts");
    U32 current = 0, sample = 0, supporters = 0;
    bool have = false, have_sample = false;
    char state = 0;
    Str last_allele;
    std::vector<Allele> parsed;
    std::vector<std::pair<U32, Str>> fields;
    // Runs already order each site's records by sample, state and allele.
    // Consume each sample once; adjacent equal alleles need no map/set nodes
    // and only unique observations need parsing for FORMAT output.
    auto finish_sample = [&] {
        if (!have_sample)
            return;
        const auto &k = sites[current].key;
        auto &columns = k.snp ? snp_columns : indel_columns;
        if (parsed.empty())
            fields.emplace_back(sample, state == '.' ? columns.missing : columns.reference);
        else {
            ++supporters;
            Str field = parsed[0][0];
            static constexpr size_t snp_order[]{1, 2, 3, 4, 5, 6, 7};
            static constexpr size_t indel_order[]{1, 2, 3, 4, 5, 8, 9, 6, 7};
            const auto *order = k.snp ? snp_order : indel_order;
            for (size_t col = 0; col < (k.snp ? 7 : 9); ++col) {
                field += ':';
                for (size_t i = 0; i < parsed.size(); ++i) {
                    if (i)
                        field += ',';
                    append_escaped(field, parsed[i][order[col]]);
                }
            }
            fields.emplace_back(sample, std::move(field));
        }
        parsed.clear();
        if (parsed.capacity() > 64)
            std::vector<Allele>().swap(parsed);
        last_allele.clear();
        state = 0;
        have_sample = false;
    };
    auto emit = [&] {
        finish_sample();
        if (!have || !supporters)
            return;
        const auto &site = sites[current];
        const auto &k = site.key;
        auto &columns = k.snp ? snp_columns : indel_columns;
        auto &out = k.snp ? snp : indel;
        U64 start = k.snp ? (k.pos ? U64(k.pos) - 1 : 0) : k.pos;
        U64 end = k.snp ? start + 1 : std::max<int64_t>(int64_t(start) + 1, k.end);
        columns.advance(start);
        Str id = (k.snp ? "SNP_" : "INDEL_") + token + '_' + std::to_string(k.pos) + '_' +
                 sha256("[" + quote(chrom) + ", " + k.json() + "]").substr(0, 20);
        out << chrom << '\t' << k.pos << '\t' << id << '\t' << k.ref << '\t' << k.alt << '\t' << site.qual
            << '\t' << site.filter << '\t' << nsup_info(site.info, supporters) << '\t'
            << (k.snp ? SNP_FORMAT : SV_FORMAT);
        columns.write(out, end, fields);
    };
    auto accept = [&](const SmallRecord &r) {
        if (!have || r.site != current) {
            emit();
            current = r.site;
            have = true;
            fields.clear();
            supporters = 0;
        }
        if (have_sample && r.sample != sample)
            finish_sample();
        sample = r.sample;
        have_sample = true;
        if (r.state == '1') {
            if (parsed.empty() || r.allele != last_allele) {
                parsed.push_back(parse_legacy(r.allele));
                last_allele = r.allele;
            }
        } else {
            if (!state || r.state == '.')
                state = r.state;
        }
    };
    if (runs.paths.empty()) {
        std::sort(records.begin(), records.end(), less);
        for (auto &r : records)
            accept(r);
    } else {
        if (!records.empty())
            flush();
        runs.finish(accept);
    }
    emit();
    checked_close(snp);
    checked_close(indel);
}

static std::pair<U32, uint8_t> label_value(const Str &text, Intern &groups) {
    auto scalar = [](View s) -> std::pair<U32, uint8_t> {
        if (s == ".")
            return {0, 0};
        require(s.size() > 1 && s.back() == 'H', "invalid SNP LABEL_H");
        int64_t n = integer(s.substr(0, s.size() - 1));
        require(n >= INT32_MIN && n <= INT32_MAX, "SNP LABEL_H outside signed 32-bit encoding");
        return {U32((U64(n) << 1) ^ U64(n >> 31)), 1};
    };
    if (text.find('&') == Str::npos)
        return scalar(text);
    for (auto p : split(text, '&'))
        scalar(p);
    return {groups.add(text), 2};
}
static bool alternative_record(View text, const std::unordered_set<U32> &entries) {
    if (text.find("PAPATH=alt") != View::npos)
        return true;
    if (text.find("PAMAP=") == View::npos)
        return false;
    static const std::regex pamap(R"((?:^|;)PAMAP=([0-9,]+)(?:;|$))");
    std::match_results<View::const_iterator> m;
    if (!std::regex_search(text.begin(), text.end(), m, pamap))
        return false;
    Str ids = m[1].str();
    for (auto x : split(ids, ','))
        if (!entries.count(coordinate(x)))
            return false;
    return true;
}
struct Scan {
    Str directory;
    bool compact;
    U32 cutoff;
    std::ofstream binary;
    std::map<Str, std::vector<std::pair<U64, U64>>> sections;
    std::map<Str, U64> counts;
    std::vector<Str> chrom_order, metadata, all_samples;
    std::unordered_set<Str> seen_samples;
    std::unordered_map<Str, Str> keys, names;
    std::map<Str, std::vector<std::tuple<Str, U32, U32>>> coverage;
    Intern alleles, groups, query_keys;
    std::vector<std::vector<Str>> queries;
    U64 offset = 0;
    Str last;
    // At most 32 shard handles. This limit applies independently of cohort size.
    std::map<Str, std::pair<U64, std::unique_ptr<std::ofstream>>> handles;
    U64 tick = 0;
    Scan(Str dir, bool mode, U32 limit) : directory(std::move(dir)), compact(mode), cutoff(limit) {
        fs::create_directories(directory);
        if (compact) {
            binary.open(directory + "/records.bin", std::ios::binary);
            require(bool(binary), "cannot create SNP records");
        }
    }
    Str chrom(const Str &name) {
        auto it = keys.find(name);
        if (it != keys.end())
            return it->second;
        Str key = sha256(name);
        keys[name] = key;
        names[key] = name;
        chrom_order.push_back(key);
        return key;
    }
    void text(const Str &name, const Str &line) {
        Str key = chrom(name);
        auto it = handles.find(key);
        if (it == handles.end()) {
            if (handles.size() >= 32) {
                auto old = std::min_element(handles.begin(), handles.end(), [](const auto &a, const auto &b) {
                    return a.second.first < b.second.first;
                });
                checked_close(*old->second.second);
                handles.erase(old);
            }
            auto f = std::make_unique<std::ofstream>(directory + "/" + key + ".vcf",
                                                     std::ios::app | std::ios::binary);
            require(bool(*f), "cannot create small shard");
            it = handles.emplace(key, std::make_pair(tick++, std::move(f))).first;
        }
        it->second.first = tick++;
        *it->second.second << line << '\n';
    }
    void add(const std::vector<View> &row, const Str &sample, const Allele *a, char state) {
        U32 pos = coordinate(row[1]);
        require(pos > 0, "SNP VCF position must be positive");
        Str contig = ".", kind = "SNP", allele = ".", strand = "+";
        U32 qpos = 0, label = 0;
        uint8_t mode = 0;
        if (a) {
            require((*a)[0] == "1" && (*a)[1].find("SNP") != Str::npos && (*a)[2] == "1",
                    "malformed SNP observation");
            require(upper((*a)[3]) == upper(Str(row[4])), "SNP observation/base disagreement");
            const auto &q = (*a)[5];
            require(q.size() > 1 && (q.back() == '+' || q.back() == '-') &&
                        std::all_of(q.begin(), q.end() - 1, [](char c) { return c >= '0' && c <= '9'; }),
                    "invalid SNP QUERYCOORD");
            qpos = coordinate(View(q).substr(0, q.size() - 1));
            strand = q.back();
            contig = (*a)[4];
            kind = (*a)[1];
            allele = (*a)[6];
            std::tie(label, mode) = label_value((*a)[7], groups);
        }
        std::vector<Str> query{sample,     contig, strand, kind, Str(1, state), std::to_string(mode),
                               Str(row[6])};
        // Length prefixes avoid collisions when custom values contain separators.
        Str query_key;
        for (auto &part : query)
            query_key += std::to_string(part.size()) + ":" + part;
        U32 qi = query_keys.add(query_key);
        if (qi == queries.size())
            queries.push_back(std::move(query));
        Record r{pos, qpos, qi, alleles.add(allele), label, bases(row[3], row[4])};
        pack(binary, r);
        Str key = chrom(Str(row[0]));
        if (key != last) {
            sections[key].emplace_back(offset, 0);
            last = key;
        }
        ++sections[key].back().second;
        ++counts[key];
        offset += 21;
    }
    void file(const Str &path) {
        Lines input(path);
        Str line;
        U64 number = 0;
        std::vector<Str> samples;
        std::unordered_set<U32> alternative;
        while (input.next(line)) {
            ++number;
            try {
                if (line.rfind("##referenceCoverage=<", 0) == 0) {
                    auto v = structured(line, "referenceCoverage");
                    Str sample = get(v, "Sample");
                    if (!sample.empty()) {
                        Str name = get(v, "Chrom");
                        chrom(name);
                        if (compact)
                            coverage[name].emplace_back(sample, coordinate(get(v, "Start")),
                                                        coordinate(get(v, "End")));
                        else
                            text(name, line);
                    }
                } else if (line.rfind("##", 0) == 0) {
                    if (line.rfind("##pseudoLinearMapping=<", 0) == 0) {
                        auto v = structured(line, "pseudoLinearMapping");
                        if (get(v, "Path") == "alt" && !get(v, "ID").empty())
                            alternative.insert(coordinate(v["ID"]));
                    }
                    if (compact || line.rfind("##referenceCoverage", 0) != 0)
                        metadata.push_back(trim(line));
                } else if (line.rfind("#CHROM\t", 0) == 0) {
                    auto fields = split(line, '\t');
                    require(fields.size() > 9, "expected unique sample names in the VCF header");
                    samples.clear();
                    std::unordered_set<Str> unique;
                    for (size_t i = 9; i < fields.size(); ++i) {
                        Str s(fields[i]);
                        require(unique.insert(s).second, "expected unique sample names in the VCF header");
                        samples.push_back(s);
                    }
                } else if (!trim(line).empty() && line[0] != '#') {
                    require(!samples.empty(), "missing #CHROM sample header");
                    auto row = split(line, '\t');
                    require(row.size() >= 9, "VCF record has fewer than 9 columns");
                    if (alternative_record(row[7], alternative))
                        continue;
                    require(row.size() == 9 + samples.size(), "sample column count differs from header");
                    if (compact) {
                        if (!snp_format(row[8]))
                            continue;
                        for (size_t i = 0; i < samples.size(); ++i) {
                            auto field = row[9 + i], gt = field.substr(0, field.find(':'));
                            auto values = observations(field, row[8], true);
                            if (values.empty()) {
                                require(no_call(gt), "malformed allele");
                                add(row, samples[i], nullptr, gt == "0" || gt == "0/0" ? '0' : '.');
                            }
                            for (auto &a : values)
                                add(row, samples[i], &a, '1');
                        }
                    } else {
                        SmallKey k;
                        if (small_key(row, cutoff, k))
                            text(Str(row[0]), line);
                    }
                }
            } catch (const std::exception &e) {
                throw std::runtime_error(path + ":" + std::to_string(number) + ": " + e.what());
            }
        }
        require(!samples.empty(), path + ": missing #CHROM sample header");
        for (auto &sample : samples)
            if (seen_samples.insert(sample).second)
                all_samples.push_back(sample);
    }
    void finish(const Str &output) {
        if (compact)
            checked_close(binary);
        for (auto &h : handles)
            checked_close(*h.second.second);
        handles.clear();
        Str chroms = "{";
        for (size_t i = 0; i < chrom_order.size(); ++i) {
            if (i)
                chroms += ",";
            auto &k = chrom_order[i];
            chroms += quote(k) + ":" + quote(names[k]);
        }
        chroms += '}';
        if (compact) {
            std::ofstream source(directory + "/source.json");
            source << "{\"protocol\":4,\"directory\":" << quote(fs::absolute(directory).string())
                   << ",\"queries\":[";
            for (size_t i = 0; i < queries.size(); ++i) {
                if (i)
                    source << ',';
                const auto &q = queries[i];
                source << '[';
                for (size_t j = 0; j < q.size(); ++j) {
                    if (j)
                        source << ',';
                    if (j == 5)
                        source << q[j];
                    else
                        source << quote(q[j]);
                }
                source << ']';
            }
            source << "],\"alleles\":" << json_strings(alleles.values)
                   << ",\"label_groups\":" << json_strings(groups.values) << ",\"chroms\":" << chroms
                   << ",\"coverage\":{";
            bool comma = false;
            for (auto &[name, values] : coverage) {
                if (comma)
                    source << ',';
                comma = true;
                source << quote(name) << ":[";
                for (size_t i = 0; i < values.size(); ++i) {
                    if (i)
                        source << ',';
                    auto &[s, a, b] = values[i];
                    source << '[' << quote(s) << ',' << a << ',' << b << ']';
                }
                source << ']';
            }
            source << "},\"lengths\":{},\"counts\":{";
            comma = false;
            for (auto &[k, n] : counts) {
                if (comma)
                    source << ',';
                comma = true;
                source << quote(k) << ':' << n;
            }
            source << "},\"sections\":{";
            comma = false;
            for (auto &[k, spans] : sections) {
                if (comma)
                    source << ',';
                comma = true;
                source << quote(k) << ":[";
                for (size_t i = 0; i < spans.size(); ++i) {
                    if (i)
                        source << ',';
                    source << '[' << spans[i].first << ',' << spans[i].second << ']';
                }
                source << ']';
            }
            source << "},\"metadata\":[]}";
            checked_close(source);
        }
        std::ofstream result(output);
        result << "{\"samples\":" << json_strings(all_samples) << ",\"chroms\":" << chroms
               << ",\"metadata\":" << json_strings(metadata)
               << ",\"source\":" << quote(directory + "/source.json") << '}';
        checked_close(result);
    }
};
static void scan(Reader &request, const Str &output, bool compact) {
    Str dir = request.str();
    U32 cutoff = request.u32();
    Scan scanner(dir, compact, cutoff);
    for (auto &path : request.strings())
        scanner.file(path);
    scanner.finish(output);
}

int main(int argc, char **argv) {
    try {
        require(argc >= 4,
                "usage: SmallVariantMerge compact|small|scan-snp|scan-small REQUEST OUTPUT [INDEL_OUTPUT]");
        Reader request(argv[2]);
        require(request.str() == "small-variant-merge-v1", "unsupported native request protocol");
        Str command = argv[1];
        if (command == "compact")
            compact(request, argv[3]);
        else if (command == "small") {
            require(argc == 5, "missing indel output");
            small_merge(request, argv[3], argv[4]);
        } else if (command == "scan-snp" || command == "scan-small")
            scan(request, argv[3], command == "scan-snp");
        else
            throw std::runtime_error("unknown native merge command");
        require(request.f.peek() == EOF, "trailing data in native merge request");
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "SmallVariantMerge: " << e.what() << '\n';
        return 1;
    }
}
