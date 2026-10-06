"""Insertion SNPs finished by SV-merge chromosome jobs concat byte-identically."""
import json
from pathlib import Path
import shutil
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import graphvcfmerge as vcf  # noqa: E402
import graphvcfmerge_insertion_store as store  # noqa: E402
import graphvcfmerge_snp as snp  # noqa: E402
import graphvcfmerge_snp_compact as compact  # noqa: E402

SAMPLES = ["b", "a", "c"]  # SNP manifest order: input VCF order


def write_vcf(path, sample, rows):
    text = ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=200>\n##contig=<ID=chr2,length=200>\n"
            f"##referenceCoverage=<Chrom=chr1,Sample={sample},Start=0,End=100>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n" + "".join(rows))
    path.write_text(text)
    return path


def snp_row(chrom, pos, sample, alt="C"):
    field = f"1:SNP:1:{alt}:{sample}:2+:allele:0H"
    return f"{chrom}\t{pos}\told_snp\tA\t{alt}\t.\tPASS\tNSUP=1\t{vcf.SNP_FORMAT}\t{field}\n"


def original_parts(tmp):
    paths = [write_vcf(tmp / f"{s}.vcf", s, [snp_row("chr2", 9, s), snp_row("chr1", 5, s)]) for s in SAMPLES]
    root = tmp / "snps"
    manifest = compact.scan(paths, root)
    for key in manifest["chroms"]:
        compact.merge_chrom(root, key)
    return root


# (insertion name, chromosome job, template, carriers: [(sample index, member sequence)])
TEMPLATE = "ACGTACGTACGTACGTACGT"
INSERTIONS = [
    ("INS_chr1_10", "chr1", [(1, "ACGTACCTACGTACGTACGT"), (0, "ACGTACGTACGAACGTACGT")]),
    ("INS_chr1_40", "chr1", [(2, "TCGTACGTACGTACGTACGA"), (1, "TCGTACGTACGTACGTACGT")]),
    ("INS_chr2_15", "chr2", [(0, "ACGTACGTACGTACGTACGG")]),
    ("INS_chr2_90", "chr2", [(1, "ACGTACGTTCGTACGTACGT"), (2, "ACGTACGTTCGTACGTACGT"),
                              (0, "ACGTACGTACGTACGTACGT")]),
]


def body(member):
    # Equal lengths: one '=' or 'X' run per position change.
    out, run, op = [], 0, None
    for ref, alt in zip(TEMPLATE, member):
        current = "=" if ref == alt else "X"
        if current != op and run:
            out.append(f"{run}{op}")
            run = 0
        op, run = current, run + 1
    out.append(f"{run}{op}")
    return "".join(out)


def save_store(insertions_root):
    for name, chrom, carriers in INSERTIONS:
        sequences = [TEMPLATE] + [seq for _i, seq in carriers]
        bodies = [""] + [body(seq) for _i, seq in carriers]
        sources = [(0, "rep:ctg", "allele", ".", 100, 120, "+")] + [
            (index, f"{SAMPLES[index]}:ctg", "allele", ".", 200, 220, "+") for index, _seq in carriers]
        compact.save_insertion(insertions_root, name, sequences, bodies, sources, 0, SAMPLES, owner=chrom)
    # A nested --exact realignment adds one observation to an insertion.
    compact.save_insertion_realignment(insertions_root, "INS_chr2_90", {
        "chrom": "INS_chr2_90", "drops": [], "adds": [["c", 3, "G", "T", "c:ctg", 203, "+"]]})
    store.flush()
    compact.finalize_insertions(insertions_root, insertion_ids=[name for name, _c, _s in INSERTIONS])


def sv_part(directory, chrom):
    path = directory / f"{chrom}.part"
    rows = [f"{chrom}\t{10 + i}\t{name}\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=20;END={10 + i}\tGT\t1\n"
            for i, (name, job, _c) in enumerate(INSERTIONS) if job == chrom]
    path.write_text("".join(rows))
    return str(path)


def concat(root, insertions, out):
    snp.concat(root, out, insertion_snps=str(insertions))
    return out.read_bytes()


@pytest.fixture
def cohort(tmp_path):
    root = original_parts(tmp_path)
    insertions = tmp_path / "insertions"
    save_store(insertions)
    parts = {chrom: sv_part(tmp_path, chrom) for chrom in ("chr1", "chr2")}
    expected = concat(root, insertions, tmp_path / "expected.vcf")
    assert b"INS_chr2_90" in expected and b"INS_chr1_10" in expected
    return tmp_path, root, insertions, parts, expected


def finish(insertions, parts, chrom, samples=SAMPLES):
    compact.finish_chrom_insertions(insertions, [parts[chrom], None], samples, f"locus_{chrom}")


def finished_count(capsys):
    text = capsys.readouterr().err
    marker = [line for line in text.splitlines() if "finished by chromosome jobs" in line]
    return int(marker[-1].split("(")[-1].split()[0]) if marker else 0


LABELS = ["locus_chr1", "locus_chr2"]


def manifest_of(insertions):
    return json.loads((insertions / "manifest.json").read_text())


def test_all_finished_streams_parts_byte_identically(cohort, capsys):
    """Finished parts + a manifest naming only the chromosomes reproduce the
    per-insertion concat exactly; nothing per insertion is kept."""
    tmp, root, insertions, parts, expected = cohort
    finish(insertions, parts, "chr1")
    finish(insertions, parts, "chr2")
    assert sorted(p.name for p in (insertions / "finished").iterdir()) == [
        f"{label}{suffix}" for label in LABELS for suffix in (".contigs", ".json", ".part")]
    info = json.loads((insertions / "finished" / "locus_chr1.json").read_text())
    assert set(info) == {"protocol", "samples", "insertions", "sites"} and info["insertions"] == 2
    compact.write_finished_manifest(insertions, LABELS)
    assert manifest_of(insertions) == {"protocol": compact.PROTOCOL, "layout": compact.MANIFEST_LAYOUT,
                                       "chromosomes": LABELS}
    capsys.readouterr()
    assert concat(root, insertions, tmp / "new.vcf") == expected
    assert finished_count(capsys) == 4


def test_manifest_needs_every_chromosome_finished(cohort, capsys):
    tmp, root, insertions, parts, expected = cohort
    finish(insertions, parts, "chr2")
    with pytest.raises(ValueError, match="not finished"):
        compact.write_finished_manifest(insertions, LABELS)
    # The per-insertion manifest still merges everything itself.
    capsys.readouterr()
    assert concat(root, insertions, tmp / "old.vcf") == expected
    assert finished_count(capsys) == 0


def test_rerun_chromosome_replaces_its_part(cohort, capsys):
    """A chromosome saved again drops its old part first; the new part holds
    the new data."""
    tmp, root, insertions, parts, expected = cohort
    finish(insertions, parts, "chr1")
    finish(insertions, parts, "chr2")
    compact.discard_finished(insertions, "locus_chr1")
    assert compact.finished_info(insertions, "locus_chr1") is None
    assert not any(p.name.startswith("locus_chr1") for p in (insertions / "finished").iterdir())
    with pytest.raises(ValueError, match="not finished"):
        compact.write_finished_manifest(insertions, LABELS)
    compact.save_insertion(insertions, "INS_chr1_40", [TEMPLATE, "ACGTACGTACGTACGTACCT"],
                           ["", body("ACGTACGTACGTACGTACCT")],
                           [(0, "rep:ctg", "allele", ".", 1, 21, "+"), (0, "b:ctg", "allele", ".", 5, 25, "+")],
                           0, SAMPLES, owner="chr1")
    store.flush()
    fresh_dir = tmp / "fresh"
    shutil.copytree(insertions, fresh_dir)
    shutil.rmtree(fresh_dir / "finished")
    compact.finalize_insertions(fresh_dir, insertion_ids=[name for name, _c, _s in INSERTIONS])
    fresh = concat(root, fresh_dir, tmp / "fresh.vcf")
    assert fresh != expected  # the re-saved insertion really changed
    finish(insertions, parts, "chr1")
    compact.write_finished_manifest(insertions, LABELS)
    capsys.readouterr()
    assert concat(root, insertions, tmp / "rerun.vcf") == fresh
    assert finished_count(capsys) == 4


def test_part_in_another_sample_order_is_reordered(cohort, capsys):
    tmp, root, insertions, parts, expected = cohort
    finish(insertions, parts, "chr1", samples=["a", "b", "c"])
    finish(insertions, parts, "chr2")
    compact.write_finished_manifest(insertions, LABELS)
    capsys.readouterr()
    assert concat(root, insertions, tmp / "order.vcf") == expected
    assert finished_count(capsys) == 4


def test_name_used_by_two_chromosomes_fails_concat(cohort):
    """One store entry per name: a name in two chromosomes lost one of them."""
    tmp, root, insertions, parts, _expected = cohort
    with open(parts["chr2"], "a") as handle:
        handle.write("chr2\t99\tINS_chr1_10\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=20;END=99\tGT\t1\n")
    finish(insertions, parts, "chr1")
    finish(insertions, parts, "chr2")
    compact.write_finished_manifest(insertions, LABELS)
    with pytest.raises(ValueError, match="used by two chromosomes"):
        concat(root, insertions, tmp / "twice.vcf")


def test_unsaved_insertion_fails_its_chromosome(cohort):
    _tmp, _root, insertions, parts, _expected = cohort
    with open(parts["chr1"], "a") as handle:
        handle.write("chr1\t99\tINS_never_saved\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=20;END=99\tGT\t1\n")
    with pytest.raises(ValueError, match="no saved SNP entry"):
        finish(insertions, parts, "chr1")


def bundle_entries(path):
    """Every committed (key, offset), all generations, by walking the bundle."""
    with open(path, "rb") as handle:
        return store._walk(handle.fileno(), 0, path.stat().st_size)[0]


def index_records(path):
    data = store._index_path(path).read_bytes()
    return [store._IDX.unpack_from(data, offset) for offset in range(0, len(data) - len(data) % 16, 16)]


def test_store_index_is_16_bytes_per_entry(cohort):
    _tmp, _root, insertions, _parts, _expected = cohort
    bundles = sorted(insertions.glob("insertions-*.bundle"))
    assert bundles
    for bundle in bundles:
        entries = bundle_entries(bundle)
        assert store._index_path(bundle).stat().st_size == 16 * len(entries)
        assert index_records(bundle) == [(int.from_bytes(key[:8], "little"), offset)
                                         for key, offset in entries]
    # Lookups through the index find the same latest generations as a walk.
    keys = [compact.chrom_key(name) for name, _c, _s in INSERTIONS]
    indexed = store.find_many(insertions, keys)
    for path in insertions.glob("insertions-*.idx"):
        path.unlink()
    store._TABLES.clear()
    walked = store.find_many(insertions, keys)
    assert {k: v["offset"] for k, v in indexed.items()} == {k: v["offset"] for k, v in walked.items()}
    assert len(indexed) == len(INSERTIONS)


def test_store_index_recovers_a_killed_writer(cohort):
    """Entries committed without index records (a writer killed in between)
    are still found, and the next append to that bundle indexes them."""
    _tmp, _root, insertions, _parts, _expected = cohort
    keys = [compact.chrom_key(name) for name, _c, _s in INSERTIONS]
    before = {k: v["offset"] for k, v in store.find_many(insertions, keys).items()}
    bundle = store._path(insertions, keys[0])
    index = store._index_path(bundle)
    data = index.read_bytes()
    index.write_bytes(data[:-16] + b"partial")      # last record lost, a torn one left
    store._TABLES.clear()
    assert {k: v["offset"] for k, v in store.find_many(insertions, keys).items()} == before
    # Append to the same bundle: the index is repaired, then extended.
    new_key = bytes([int(keys[0][:2], 16)]) .hex() + "ab" * 31
    assert store._path(insertions, new_key) == bundle
    store.put(insertions, new_key, {"x": 1})
    entries = bundle_entries(bundle)
    assert index_records(bundle) == [(int.from_bytes(key[:8], "little"), offset) for key, offset in entries]
    store._TABLES.clear()
    assert store.find_many(insertions, [new_key])[new_key]["offset"] == entries[-1][1]
    assert {k: v["offset"] for k, v in store.find_many(insertions, keys).items()} == before


def ins_row(chrom, pos, sequence, sample):
    span = len(sequence)
    field = f"1:INS:{span}:>{span}I:{sample}:2-{2 + span}+:0:.:allele:0H0H"
    return (f"{chrom}\t{pos}\tINS_{chrom}_{pos}_1\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN={span};MAXSIZE={span};"
            f"END={pos};SEQ={sequence};NSUP=1\t{vcf.SV_FORMAT}\t{field}\n")


def test_backfill_stage_finishes_old_run_byte_identically(tmp_path, capsys):
    """A run whose chromosome jobs predate finishing: backfill, then concat."""
    import random
    import cohort_vcf_merge as cohort
    rng = random.Random(11)
    rows = {sample: [] for sample in SAMPLES}
    for chrom, pos in (("chr1", 40), ("chr1", 120), ("chr2", 60)):
        template = "".join(rng.choice("ACGT") for _ in range(150))
        for number, sample in enumerate(SAMPLES):
            member = list(template)
            if number:  # each later carrier differs at its own position
                index = 20 * number + pos % 7
                member[index] = next(b for b in "ACGT" if b != template[index])
            rows[sample].append(ins_row(chrom, pos, "".join(member), sample))
    paths = [write_vcf(tmp_path / f"{s}.vcf", s, sorted(rows[s], key=lambda r: (r.split()[0], int(r.split()[1]))))
             for s in SAMPLES]
    cohort.run(tmp_path / "cohort", "all", paths=paths, processes=2, var_in_insert=0, keep_merge_tmpdir=True)
    root, = (tmp_path / "cohort" / "tmp" / "merge_only").iterdir()
    insertion_root = root / "insertion_snps"
    finished = insertion_root / "finished"
    expected = (root / "snp.vcf").read_bytes()
    # A current run: chromosome jobs finished, the manifest names only them.
    assert manifest_of(insertion_root)["layout"] == compact.MANIFEST_LAYOUT
    assert "sources" not in manifest_of(insertion_root)
    assert sorted(p.suffix for p in finished.iterdir()) == [".contigs"] * 2 + [".json"] * 2 + [".part"] * 2
    assert expected.count(b"\nINS") or b"SNP_" in expected

    # An old run: no finished parts, one manifest listing every insertion.
    shutil.rmtree(finished)
    compact.finalize_insertions(insertion_root, [p for p in (root / "sv.vcf", root / "sv.vcf.small.vcf")
                                                  if p.is_file()])
    assert "sources" in manifest_of(insertion_root)
    assert snp.concat(root / "snp_shards", tmp_path / "old.vcf") is None
    assert (tmp_path / "old.vcf").read_bytes() == expected
    # One chromosome alone (a SLURM job each), then the rest.
    first = (root / "sv_shards" / "chroms.txt").read_text().split()[0]
    vcf.main(["--sv", "--stage", "finish-insertions", "--shards-dir", str(root / "sv_shards"),
              "--insertion-snps", str(root / "insertion_snps"), "--chrom", first])
    assert [p.name for p in finished.glob("*.json")] == [vcf._safe_locus_key(first) + ".json"]
    assert "sources" in manifest_of(insertion_root)     # not every chromosome finished yet
    capsys.readouterr()
    assert vcf.main(["--sv", "--stage", "finish-insertions", "--shards-dir", str(root / "sv_shards"),
                     "--insertion-snps", str(root / "insertion_snps"), "--chrom", "no_such_chrom"]) not in (0, None)
    assert "no_such_chrom" in capsys.readouterr().err
    assert vcf.main(["--sv", "--stage", "finish-insertions", "--shards-dir", str(root / "sv_shards"),
                     "--insertion-snps", str(root / "insertion_snps"), "--processes", "2"]) in (0, None)
    assert sorted(p.suffix for p in finished.iterdir()) == [".contigs"] * 2 + [".json"] * 2 + [".part"] * 2
    log = capsys.readouterr().err
    assert "reading entries at indexed offsets from manifest.json" in log
    assert "manifest.json now lists the 2 finished chromosome(s)" in log
    assert manifest_of(insertion_root)["layout"] == compact.MANIFEST_LAYOUT
    snp.concat(root / "snp_shards", tmp_path / "backfilled.vcf")
    assert (tmp_path / "backfilled.vcf").read_bytes() == expected
    line = [l for l in capsys.readouterr().err.splitlines() if "finished by chromosome jobs" in l][-1]
    loci = int(line.split(" insertion loci")[0].split()[-1])
    assert loci >= 1 and f"({loci} finished by chromosome jobs)" in line  # every locus from a part

    # Rerunning skips chromosomes already finished.
    vcf.main(["--sv", "--stage", "finish-insertions", "--shards-dir", str(root / "sv_shards"),
              "--insertion-snps", str(root / "insertion_snps")])
    assert "0 chromosome(s) to finish" in capsys.readouterr().err



@pytest.mark.parametrize("batch", [1, 128])
def test_indexed_direct_read_matches_locked_lookup(cohort, monkeypatch, batch):
    """Backfill reads (manifest offsets, open unlocked bundles, file-order
    batches) finish every chromosome byte-identically, realignment included."""
    _tmp, _root, insertions, parts, _expected = cohort
    monkeypatch.setattr(compact, "DIRECT_READ_BATCH", batch)
    index = compact.backfill_index(insertions)
    assert index is not None and "INS_chr2_90" not in index[1]  # realign keys are digests
    assert len(index[1]) == 1
    finished = insertions / "finished"
    for chrom in ("chr1", "chr2"):
        finish(insertions, parts, chrom)
        locked = [(finished / f"locus_{chrom}{suffix}").read_bytes() for suffix in (".part", ".contigs", ".json")]
        with monkeypatch.context() as patch:  # no bundle lock on the indexed path
            patch.setattr(store.fcntl, "flock", lambda *_a: pytest.fail("bundle locked"))
            compact.finish_chrom_insertions(insertions, [parts[chrom], None], SAMPLES, f"locus_{chrom}",
                                            index=index)
        direct = [(finished / f"locus_{chrom}{suffix}").read_bytes() for suffix in (".part", ".contigs", ".json")]
        assert direct == locked and locked[0]
