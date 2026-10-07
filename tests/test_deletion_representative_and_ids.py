"""Deletion groups take the smallest member as representative; row IDs use I_/D_."""
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import graphvcfmerge as vcf  # noqa: E402

TOOLS = Path(__file__).resolve().parents[1] / "tools"


def member(index, pos, size):
    return (pos, pos + size, size, None, index)


def refine(svtype, members):
    return vcf._fast_refine_one_cluster(svtype, members, None, {"size_similarity": 0.7})


def test_deletion_representative_is_the_smallest_member():
    members = [member(0, 1002, 75), member(1, 1000, 60), member(2, 995, 70)]
    (representative, chosen), = refine("DEL", members)
    assert representative[2] == 60
    assert sorted(item[4] for item in chosen) == [0, 1, 2]


def test_deletion_ties_go_to_the_leftmost():
    (representative, _chosen), = refine("DEL", [member(0, 1005, 60), member(1, 1000, 60)])
    assert representative[0] == 1000


def test_smallest_seed_absorbs_up_to_its_size_window():
    # 60 / 0.7 = 85.7: the 90 bp deletion forms its own group.
    groups = refine("DEL", [member(0, 1000, 60), member(1, 1000, 80), member(2, 1000, 90)])
    assert [(rep[2], sorted(item[2] for item in chosen)) for rep, chosen in groups] == [
        (60, [60, 80]), (90, [90])]


def test_other_types_keep_the_largest_representative():
    (representative, _chosen), = refine("DUP", [member(0, 1000, 60), member(1, 1000, 75)])
    assert representative[2] == 75


def test_row_ids_use_short_prefixes_and_shorten_older_ones():
    assert vcf._prefixed_variant_id("INS", "INS_chr1_100_3") == "I_chr1_100_3"
    assert vcf._prefixed_variant_id("DEL", "DEL_chr1_100_3") == "D_chr1_100_3"
    assert vcf._prefixed_variant_id("INS", "I_chr1_100_3") == "I_chr1_100_3"
    assert vcf._prefixed_variant_id("DEL", "chr1_100") == "D_chr1_100"
    assert vcf._prefixed_variant_id("DUP", "chr1_100") == "DUP_chr1_100"
    assert vcf._prefixed_variant_id("DUP", "DUP_chr1_100") == "DUP_chr1_100"


def test_nested_rows_of_both_id_styles_are_recognized():
    for path in (TOOLS / "grvcf_to_vcf.py", TOOLS / "vcf_chromfix.py",
                 Path(vcf.__file__).with_name("gfa_sample_gaf.py")):
        pattern = re.search(r"NESTED_CHROM = re\.compile\(r'([^']+)'\)", path.read_text()).group(1)
        for chrom in ("I_chr1_100_1", "D_chr1_100_1", "INS_chr1_100_1", "DEL_chr1_100_1"):
            assert re.match(pattern, chrom), (path.name, chrom)
        assert not re.match(pattern, "chr1") and not re.match(pattern, "NC_060925.1")


def test_count_sv_skips_nested_rows_of_both_id_styles(tmp_path):
    rows = ["chr1\t100\tD_chr1_100_1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-60;END=160",
            "D_chr1_100_1\t5\tI_D_chr1_100_1_5_1\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=70;END=5",
            "DEL_chr1_200_1\t5\tINS_DEL_chr1_200_1_5_1\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=70;END=5"]
    path = tmp_path / "rows.vcf"
    path.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "\n".join(rows) + "\n")
    out = subprocess.run(["bash", str(TOOLS / "count_sv.sh"), str(path)],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    assert out[1].split("\t")[1:] == ["0", "1", "0", "1"]


def test_repeat_numbering_skips_ids_already_taken():
    used = {}
    assert [vcf._unique_variant_id(used, base) for base in ("X", "X_2", "X", "X")] == [
        "X", "X_2", "X_3", "X_4"]


def test_top_level_ids_carry_their_chromosome():
    assert vcf._chrom_variant_id("DEL", "DEL_chr1_100_3", "chr1") == "D_chr1_100_3"
    assert vcf._chrom_variant_id("INS", "old_indel", "chr1") == "I_chr1_old_indel"
    assert vcf._chrom_variant_id("INS", "I_old_indel", "chr2") == "I_chr2_old_indel"
    token = vcf._safe_variant_token("I_HG19#1#chr1_100_3")
    assert vcf._chrom_variant_id("INS", token, "HG19#1#chr1") == "I_HG19_1_chr1_100_3"
    # chr10's token is not chr1's
    assert vcf._chrom_variant_id("INS", "I_chr10_5_1", "chr1") == "I_chr1_chr10_5_1"


def test_duplicate_row_ids_in_a_chromosome_fail(tmp_path):
    main, small = tmp_path / "c.part", tmp_path / "c.small.part"
    main.write_text("c\t10\tD_c_10_1\tN\t<DEL>\nD_c_10_1\t1\tI_D_c_10_1_1_1\tN\t<INS>\n")
    small.write_text("c\t50\tD_c_50_1\tN\t<DEL>\n")
    vcf._check_unique_row_ids([str(main), str(small), None])
    small.write_text("c\t50\tD_c_10_1\tN\t<DEL>\n")
    import pytest
    with pytest.raises(ValueError, match="duplicate merged row ID 'D_c_10_1'"):
        vcf._check_unique_row_ids([str(main), str(small)])


def test_inputs_sharing_one_id_give_unique_merged_ids(tmp_path):
    """Every input row has ID 'same': merged IDs are still unique within and
    across chromosomes."""
    import cohort_vcf_merge as cohort
    fmt = vcf.SV_FORMAT
    paths = []
    for sample in ("a", "b"):
        rows = []
        for chrom, pos, size in (("chr1", 100, 60), ("chr1", 900, 80), ("chr2", 100, 60)):
            rows.append(f"{chrom}\t{pos}\tsame\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-{size};"
                        f"END={pos + size};MAXSIZE={size};NSUP=1\t{fmt}\t"
                        f"1:DEL:{size}:>{size}D:{sample}ctg:{pos}+:0:.:allele:0H0H\n")
        path = tmp_path / f"{sample}.vcf"
        path.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=2000>\n"
                        "##contig=<ID=chr2,length=2000>\n"
                        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
                        + "".join(rows))
        paths.append(path)
    targets = cohort.run(tmp_path / "out", "svonly", paths=paths, processes=1)
    ids = [line.split("\t")[2] for line in Path(targets[0]).read_text().splitlines()
           if not line.startswith("#")]
    assert len(ids) == 3 and len(set(ids)) == 3
    assert sorted(ids) == ["D_chr1_same", "D_chr1_same_2", "D_chr2_same"]


def test_id_letters_count_in_bijective_base_26():
    assert [vcf._id_letters(n) for n in (1, 2, 26, 27, 52, 53, 702, 703)] == [
        "a", "b", "z", "aa", "az", "ba", "zz", "aaa"]


def test_nested_ids_are_position_and_size():
    used = {}
    parent = "I_chr1_1000_3"
    assert vcf._nested_variant_id(used, "DEL", parent, 57, 60) == "D_I_chr1_1000_3_57_60"
    # Insertions of different bases at one position and size: a, b, ...
    assert [vcf._nested_variant_id(used, "INS", parent, 12, 100) for _ in range(4)] == [
        "I_I_chr1_1000_3_12_100", "I_I_chr1_1000_3_12_100a",
        "I_I_chr1_1000_3_12_100b", "I_I_chr1_1000_3_12_100c"]
    assert vcf._nested_variant_id(used, "INS", parent, 12, 101) == "I_I_chr1_1000_3_12_101"


def test_nested_rows_are_named_by_position_and_size(tmp_path):
    """A 340 bp insertion; three 320 bp copies lack 40 bp of it and carry a
    20 bp insertion: nested rows are named <type>_<parent>_<pos>_<size>."""
    import random
    import cohort_vcf_merge as cohort
    rng = random.Random(5)
    seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    T, E, X = seq(300), seq(40), seq(20)
    alleles = {"a": T[:200] + E + T[200:], "b": T[:150] + X + T[150:], "c": T[:150] + X + T[150:]}
    paths = []
    for sample, allele in alleles.items():
        span = len(allele)
        field = f"1:INS:{span}:>{span}I:{sample}ctg:1000-{1000 + span}+:0:.:allele:0H0H"
        path = tmp_path / f"{sample}.vcf"
        path.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=5000>\n"
                        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
                        f"chr1\t1000\tI_chr1_1000_{sample}\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN={span};"
                        f"MAXSIZE={span};END=1000;SEQ={allele};NSUP=1\t{vcf.SV_FORMAT}\t{field}\n")
        paths.append(path)
    targets = cohort.run(tmp_path / "out", "all", paths=paths, processes=1)
    rows = [line.split("\t") for target in targets for line in Path(target).read_text().splitlines()
            if not line.startswith("#") and "SVTYPE=" in line]
    named = {}
    for chrom, pos, row_id, *rest in rows:
        info = dict(item.split("=", 1) for item in rest[4].split(";") if "=" in item)
        if chrom != "chr1":     # nested: <I|D>_<parent>_<POS>_<size>
            assert row_id == f"{info['SVTYPE'][0]}_{chrom}_{pos}_{abs(int(info['SVLEN']))}"
        named[row_id] = (chrom, info["SVTYPE"], abs(int(info["SVLEN"])))
    assert named["I_chr1_1000_a"] == ("chr1", "INS", 340)
    assert ("I_chr1_1000_a", "DEL", 40) in named.values()
    assert ("I_chr1_1000_a", "INS", 20) in named.values()
    ids = [row[2] for row in rows]
    assert len(ids) == len(set(ids))
