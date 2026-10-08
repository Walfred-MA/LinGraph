#!/usr/bin/env python3
"""Build the versioned grVCF format guide beside this source file.

Requires reportlab. Run from any directory: python docs/build_grvcf_guide.py
"""

from __future__ import annotations

from html import escape
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (
    HRFlowable, KeepTogether, PageBreak, Paragraph, Preformatted,
    SimpleDocTemplate, Spacer, Table, TableStyle,
)


DESTINATION = Path(__file__).with_name("LinGraph_Graph_Recursive_VCF.pdf")
NAVY = colors.HexColor("#123047")
TEAL = colors.HexColor("#007E87")
PALE = colors.HexColor("#EAF4F4")
LIGHT = colors.HexColor("#F2F5F7")
INK = colors.HexColor("#243240")
MUTED = colors.HexColor("#52616B")
LINE = colors.HexColor("#CCD9DC")


def styles():
    s = getSampleStyleSheet()
    s.add(ParagraphStyle(
        name="TitleGr", fontName="Helvetica-Bold", fontSize=25, leading=28,
        textColor=NAVY, spaceAfter=11,
    ))
    s.add(ParagraphStyle(
        name="SubtitleGr", fontName="Helvetica", fontSize=12.8, leading=17,
        textColor=TEAL, spaceAfter=16,
    ))
    s.add(ParagraphStyle(
        name="H1Gr", fontName="Helvetica-Bold", fontSize=15, leading=19,
        textColor=NAVY, spaceBefore=0, spaceAfter=9,
    ))
    s.add(ParagraphStyle(
        name="H2Gr", fontName="Helvetica-Bold", fontSize=10.5, leading=14,
        textColor=TEAL, spaceBefore=11, spaceAfter=5,
    ))
    s.add(ParagraphStyle(
        name="BodyGr", fontName="Helvetica", fontSize=9.1, leading=13.5,
        textColor=INK, spaceAfter=7,
    ))
    s.add(ParagraphStyle(
        name="SmallGr", fontName="Helvetica", fontSize=8.05, leading=11.3,
        textColor=INK, spaceAfter=5,
    ))
    s.add(ParagraphStyle(
        name="TinyGr", fontName="Helvetica", fontSize=7.35, leading=10.1,
        textColor=INK,
    ))
    s.add(ParagraphStyle(
        name="THGr", fontName="Helvetica-Bold", fontSize=8.1, leading=11,
        textColor=colors.white,
    ))
    s.add(ParagraphStyle(
        name="CalloutGr", fontName="Helvetica-Bold", fontSize=9.4, leading=14,
        textColor=NAVY,
    ))
    s.add(ParagraphStyle(
        name="CodeGr", fontName="Courier", fontSize=7.6, leading=10.6,
        textColor=INK,
    ))
    s.add(ParagraphStyle(
        name="CenterGr", fontName="Helvetica-Bold", fontSize=9.1,
        leading=12.5, alignment=TA_CENTER, textColor=NAVY,
    ))
    return s


S = styles()
STORY = []


def p(text, style="BodyGr"):
    STORY.append(Paragraph(text, S[style]))


def h1(text):
    p(text, "H1Gr")


def h2(text):
    p(text, "H2Gr")


def bullet(text):
    p("<font color='#007E87'><b>\u2022</b></font>  " + text, "BodyGr")


def callout(text):
    box = Table([[Paragraph(text, S["CalloutGr"])]], colWidths=[174 * mm])
    box.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), PALE),
        ("BOX", (0, 0), (-1, -1), 0.5, LINE),
        ("LEFTPADDING", (0, 0), (-1, -1), 12),
        ("RIGHTPADDING", (0, 0), (-1, -1), 12),
        ("TOPPADDING", (0, 0), (-1, -1), 10),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 10),
    ]))
    STORY.extend([Spacer(1, 4 * mm), box, Spacer(1, 4 * mm)])


def table(headers, rows, widths, *, compact=False):
    body_style = S["TinyGr" if compact else "SmallGr"]
    data = [[Paragraph(escape(item), S["THGr"]) for item in headers]]
    data.extend([[Paragraph(item, body_style) for item in row] for row in rows])
    tab = Table(data, colWidths=[x * mm for x in widths], repeatRows=1, hAlign="LEFT")
    tab.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), NAVY),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, LIGHT]),
        ("GRID", (0, 0), (-1, -1), 0.35, LINE),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 7 if not compact else 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 7 if not compact else 5),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ]))
    STORY.extend([tab, Spacer(1, 3 * mm)])


def code(lines):
    block = Table([[Preformatted(lines.strip("\n"), S["CodeGr"])]], colWidths=[174 * mm])
    block.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), LIGHT),
        ("LINEBEFORE", (0, 0), (0, 0), 3, TEAL),
        ("LEFTPADDING", (0, 0), (-1, -1), 10),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
    ]))
    STORY.extend([block, Spacer(1, 3 * mm)])


def next_page():
    STORY.append(PageBreak())


def header_footer(canvas, doc):
    canvas.saveState()
    w, h = A4
    canvas.setStrokeColor(LINE)
    canvas.line(18 * mm, h - 18 * mm, w - 18 * mm, h - 18 * mm)
    canvas.line(18 * mm, 17 * mm, w - 18 * mm, 17 * mm)
    canvas.setFont("Helvetica", 7.4)
    canvas.setFillColor(MUTED)
    canvas.drawString(18 * mm, 12.5 * mm, "grVCF  |  Graph Recursive VCF  |  7 October 2026")
    canvas.drawRightString(w - 18 * mm, 12.5 * mm, str(doc.page))
    canvas.restoreState()


# 1 / Purpose and mode distinction.
p("Graph Recursive VCF (grVCF)", "TitleGr")
p("A graph-aware VCF for lossless assembly information and unambiguous cohort variants", "SubtitleGr")
p("LinGraph is a variants-guided <b>pangenome graph builder</b>. It records variant sequences and assembly coordinates in grVCF and can export cohort calls as a pangenome graph. grVCF uses VCF-like rows plus graph-aware fields and parent-child relationships.")
callout("Two representations of a merged cohort share one goal: retain enough sequence, mapping, and variant information to reconstruct the represented assembly regions. The recursive (exact) representation additionally makes every carrier of one row use the same sequence and the same reference location.")
h2("Choose the representation")
table(["Mode", "Where differences live", "Same-row guarantee"], [
    ["<b>CIGAR</b><br/>Default merged output", "A row groups a representative allele with each carrier's own breakpoint, size, and alignment in FORMAT; TEMPLATEOFFSET can differ.", "No. A shared row can cover distinct member sequences and locations."],
    ["<b>Recursive / exact</b><br/><font size='8'>--exact</font>", "Carriers use the row's allele and placement. Any residual differences are separate child, flank, or SNP rows; children can nest again.", "Yes, for the row operation. A complete haplotype allele can still differ after its carried children are applied."],
], [33, 85, 56])
h2("How the hierarchy encodes graph topology")
code("reference left flank  -->  Insertion A  -->  reference right flank\n                              |\n                              +--> Child B: nested variant on Insertion A\n                                      +--> another nested change, if needed")
p("A child row's CHROM names its parent allele or shared template, so its coordinates live on that sequence. This allows variants on insertions, deletions, and duplication templates to become explicit graph branches. A graph is a downstream export of the complete grVCF, not its defining purpose.", "SmallGr")
p("<b>Scope of losslessness.</b> The representation preserves mapped assembly sequence when all variant classes, nested dependencies, sequence sources, and placement metadata are retained. Unmapped or filtered regions require separate accounting; an SV-only subset or standard VCF projection cannot reproduce all represented bases.", "SmallGr")


# 2 / Envelope, fields, and coordinates.
next_page()
h1("01  Read a grVCF record")
p("grVCF currently uses a VCF 4.2-style envelope and a .vcf extension. INFO describes a representative allele; FORMAT describes each haplotype observation and its assembly location. Read the actual FORMAT column, because caller options and merge stage can change which optional fields appear.")
table(["Column or field", "Meaning"], [
    ["CHROM / POS / ID", "A reference or alternative locus, or a parent variant ID. On a nested row POS is parent-relative. ID may be the sequence target for descendants."],
    ["REF / ALT", "SNPs use bases; structural rows can use symbolic &lt;INS&gt;, &lt;DEL&gt;, or &lt;SUB&gt;. INFO/SEQ carries representative inserted or replacement bases."],
    ["INFO: SVTYPE, END, SVLEN", "Event type, endpoint, and representative signed length. NSUP counts supporting haplotype columns, not sequencing reads."],
    ["INFO: SEQ, EXTENDGRAPHCIGAR", "Representative sequence and graph alignment. SEQ may be plain bases or a named graph encoding; keep referenced targets."],
    ["FORMAT: GT", "Haploid state: 1 supports the row, 0 is callable reference, . is missing. Do not read . as 0."],
    ["FORMAT: TYPE, SIZE, EXTENDGRAPHCIGAR", "Per-observation event, size, and alignment to the representative. CIGAR mode uses these to retain member differences."],
    ["FORMAT: ASSEMBLYCONTIG, QUERYCOORD", "Source assembly contig and zero-based query coordinate or range, with strand."],
    ["FORMAT: TEMPLATEOFFSET", "Per-observation breakpoint displacement from row POS. In exact shared placement it is 0; in CIGAR mode it may differ."],
    ["FORMAT: BASE, ALLELENAME, LABEL_H", "SNP alternate base and optional source/provenance fields. Multiple observations appear as corresponding comma-separated values."],
], [54, 120], compact=True)
h2("Coordinates and names")
p("SNP POS is one-based. For symbolic root events, the graph converter interprets replacement intervals from POS and END; a pure insertion has END=POS. A nested insertion starts at parent offset POS-1. A nested deletion removes from that offset for abs(SVLEN) bases. Interpret coordinates in the parent named by CHROM, not automatically on the reference chromosome.", "SmallGr")
p("New merged IDs use I_, D_, and S_ for insertions, deletions, and SNPs; older files may use INS_, DEL_, and SNP_. A nested insertion can be named I_&lt;parent&gt;_&lt;pos&gt;_&lt;size&gt;, with letter suffixes for distinct sequences at the same position and size. Preserve IDs when copying dependent rows.", "SmallGr")


# 3 / CIGAR mode.
next_page()
h1("02  CIGAR representation: lossless member information")
p("The default cohort merge writes <b>##graphvcfmergeVersion=cigar</b>. It chooses one representative allele for a row. Each carrier retains its own size, breakpoint, assembly coordinate, and alignment to that representative. This is a compact, lossless way to describe members without requiring one exact allele per row.")
table(["Carrier", "Row shown at chr1:100", "Own observation saved in FORMAT"], [
    ["h1", "Representative insertion A = ACGTTGCA", "8 bases at chr1:100; TEMPLATEOFFSET=0; assembly query range identifies the bases."],
    ["h2", "Same representative row A", "10 bases at chr1:105; TEMPLATEOFFSET=5; extended graph CIGAR describes GG after the fourth base."],
], [20, 62, 92])
p("Illustrative observations above are not literal VCF records. h1 and h2 support one displayed row but differ in both placement and sequence. The FORMAT arrays retain which observation belongs to which haplotype; do not replace them with the representative sequence alone.")
h2("What must be kept")
bullet("The representative SEQ or its named sequence source, all per-observation EXTENDGRAPHCIGAR fields, TYPE, SIZE, TEMPLATEOFFSET, ASSEMBLYCONTIG, QUERYCOORD, and matching comma-list entries.")
bullet("SNP and small-indel files, nested variants and insertion SNPs, plus the reference and local sequence catalogs needed by encoded alleles.")
bullet("The caller's ##pseudoLinearMapping and ##referenceCoverage lines. Merged output stores the per-sample lines in cohort.samples.headers.gz beside the cohort VCFs.")
h2("Practical effect")
p("CIGAR mode reconstructs an individual observation from its row plus the observation alignment and placement metadata. The row itself is <b>not</b> a universal allele at one breakpoint. Use the recursive / exact merge if downstream graph interpretation needs that stronger row-level guarantee.")
callout("Lossless sequence reconstruction and same-row identity are separate properties. CIGAR mode targets the first; recursive / exact mode targets both.")


# 4 / Recursive exact mode.
next_page()
h1("03  Recursive representation: one row, one operation")
p("An exact merge writes <b>##graphvcfmergeVersion=exact</b>. It realigns members to a representative placement when possible. A shared row then denotes the same sequence change at the same parent-relative position for every carrier. Differences that remain become separate nested, flank (_F), separated (_S), or SNP rows.")
p("The next abbreviated example omits FORMAT details. Row A is an insertion at chr1:100 with representative sequence ACGTTGCA. Child B is a <b>nested variant on Insertion A</b>; its CHROM is A, and it inserts GG at parent offset 4 (POS=5). S and D are other changes on A.")
table(["CHROM", "POS", "ID", "ALT / INFO", "Meaning"], [
    ["chr1", "100", "A", "&lt;INS&gt;; SEQ=ACGTTGCA", "Parent insertion, 8 bases."],
    ["A", "5", "B", "&lt;INS&gt;; SEQ=GG", "Child insertion after base 4 of A."],
    ["A", "3", "S", "G to A", "Child SNP at base 3 of A."],
    ["A", "3", "D", "&lt;DEL&gt;; SVLEN=-2", "Child deletion of bases 3-4 of A."],
], [22, 13, 17, 63, 59], compact=True)
table(["Haplotype", "Carried rows", "Complete inserted sequence"], [
    ["h1", "A", "ACGTTGCA"],
    ["h2", "A + B", "ACGTGGTGCA"],
    ["h3", "A + S", "ACATTGCA"],
    ["h4", "A + D", "ACTGCA"],
], [34, 37, 103], compact=True)
p("All carriers of A agree on A's eight-base sequence and chr1 placement. Their <i>complete</i> inserted sequences differ because they carry different child rows. If several haplotypes carry B, B itself has one GG insertion at the same offset on A. A row with CHROM=B can represent a further variant on that inserted GG sequence.")
h2("Reconstruction rule")
p("Start from the reference placement, apply the carried top-level row, then apply carried children on the unchanged parent coordinate axis. Recurse into any child allele before placing it into its parent. Sibling coordinates should not be shifted by edits already applied to the string. Strand and QUERYCOORD map the result back to the assembly.")


# 5 / Exact merge and validation.
next_page()
h1("04  What the exact merge changes")
p("The merge reads indexed assemblies and reference FASTAs with --exact. It compares shifted observations against the chosen row placement; compatible alleles are realigned there and their remaining bases are represented by child operations. If a move cannot be represented exactly, the member is kept on a separate row.")
table(["Situation", "Representation"], [
    ["Different bases within a shared insertion", "Nested INS/DEL/SUB or insertion SNP rows on the insertion ID; deeper changes recurse on child IDs."],
    ["Member deletes fewer bases than the shared deletion", "Bases retained by that member are expressed as child insertions on the deletion's coordinate space."],
    ["Shifted member, or difference at a row flank", "Realign when exact; otherwise use a separate row or _F/_S pieces so one row does not silently mix distinct placements."],
    ["Full-locus duplicated sequence", "A shared DUP_ template can carry nested variants while preserving copy-specific provenance."],
], [68, 106])
h2("Validation has two independent questions")
table(["Question", "Current checker output"], [
    ["Can the represented mapped assembly regions be reconstructed?", "tools/check_merge_lossless.py reports lossless=yes when its region reconstruction checks pass."],
    ["Does each carrier of each row agree with that row's allele and breakpoint?", "The same checker reports unified=yes; mismatches are written to line_check.tsv, including nonzero TEMPLATEOFFSET."],
], [82, 92])
p("A passing check applies to the files, assemblies, reference catalogs, and mapped regions supplied. It does not assert that unmapped bases or filtered-away variants are represented. The per-sample checker supports plain SEQ; retain encoded sequence catalogs when using other sequence modes.", "SmallGr")
h2("Sequence encoding versus variant hierarchy")
p("INFO/SEQ may reuse a named target through graph-CIGAR pieces. That is a way to store the <i>bases</i> compactly; it does not move the row or replace its child records. FORMAT/EXTENDGRAPHCIGAR describes an observation relative to a representative. The two encodings have different roles even when both contain CIGAR-like operations.")


# 6 / Commands and round trip.
next_page()
h1("05  Produce, inspect, and convert")
p("Run these examples from the repository root. vcfs.list has one individual grVCF path per line. query_paths.txt has NAME FASTA [FAI] per line, with names matching haplotype columns.")
h2("Merge individual grVCFs")
code("# CIGAR mode (default)\npython tools/merge_grvcfs.py -I vcfs.list -O merged_cigar -t 16\n\n# Recursive / exact mode: indexed assemblies and reference are required\npython tools/merge_grvcfs.py -I vcfs.list -O merged_exact -t 16 \\\n  --exact query_paths.txt --reference-fasta reference.fa")
p("The cohort merge writes cohort.sv.vcf, cohort.indel.vcf, cohort.snp.vcf, and cohort.samples.headers.gz. Read ##graphvcfmergeVersion=cigar|exact from each merged VCF. A direct LinGraph graph run also uses --exact for its default merged cohort calls.", "SmallGr")
h2("Check losslessness and shared-row identity")
code("python tools/check_merge_lossless.py \\\n  -v merged_exact/cohort.sv.vcf merged_exact/cohort.indel.vcf \\\n     merged_exact/cohort.snp.vcf \\\n  -s samples/*/*.vcf -r reference.fa -q query_paths.txt -o checks")
p("Inspect checks/summary.tsv for lossless and unified. The sample VCFs supply their mapping lines; the query FASTAs supply the expected assembly bases. Add -r for each local reference template FASTA used by the merge.", "SmallGr")
h2("Split or change representations")
code("python tools/convert_merged_grvcf.py split \\\n  -v merged_exact/cohort.sv.vcf merged_exact/cohort.indel.vcf \\\n     merged_exact/cohort.snp.vcf -r reference.fa -o individual\n\npython tools/convert_merged_grvcf.py to-exact \\\n  -v merged_cigar/cohort.sv.vcf merged_cigar/cohort.indel.vcf \\\n     merged_cigar/cohort.snp.vcf -r reference.fa -q query_paths.txt \\\n  --reference-fasta reference.fa -O merged_exact")
p("The converter also provides to-cigar. Keep cohort.samples.headers.gz next to merged files (or pass --sample-headers) so split VCFs regain ##pseudoLinearMapping and ##referenceCoverage. Equivalent placements in repetitive sequence may yield different row positions on a round trip even when reconstructed alleles agree.", "SmallGr")


# 7 / Downstream uses and references.
next_page()
h1("06  Use the complete grVCF")
h2("Build a pangenome graph")
p("Feed the complete merged SV, indel, and SNP grVCFs to the graph exporter, with the reference and required sequence sources. Parent and child alleles become graph paths and branches. The rGFA exporter includes all variant sizes by default; a size-filtered or SV-only graph omits parts of the assembly representation.")
code("python scripts/merged_vcf_to_gfa.py \\\n  -v merged_exact/cohort.sv.vcf merged_exact/cohort.indel.vcf \\\n     merged_exact/cohort.snp.vcf \\\n  -q query_paths.txt --graph-folder graph -o cohort.gfa -t 8")
p("Use the actual graph's local reference templates and alternative sequence catalog when required by your cohort. For LinGraph's full cohort workflow, --make-graph exports cohort.gfa and per-sample GAFs; --gfa-only exports from existing merged calls.", "SmallGr")
h2("Export a conventional VCF")
code("python tools/grvcf_to_vcf.py \\\n  -i merged_exact/cohort.sv.vcf -o cohort.sv.standard.vcf")
p("This is a <b>representative-allele projection</b>. The current streaming exporter removes nested rows and POS=0 records, keeps GT and a small INFO subset, and converts an &lt;INS&gt; allele with plain INFO/SEQ bases to explicit REF/ALT. It does not reconstruct each carrier's complete allele, normalize against a reference, or emit an unplaced sidecar. Graph-encoded insertion SEQ stays symbolic. Keep the original grVCFs for lossless work and graph building.", "SmallGr")
h2("Files and implementation to consult")
table(["Purpose", "Repository source"], [
    ["Caller schema and per-sample sequence encoding", "scripts/graphreftovcf.py; tools/check_vcf_lossless.py"],
    ["Cohort CIGAR / exact layout, nested rows", "scripts/graphvcfmerge.py; scripts/graphvcfmerge_snp_compact.py"],
    ["Merge and convert cohorts", "tools/merge_grvcfs.py; tools/convert_merged_grvcf.py"],
    ["Reconstruction, graph, standard export", "tools/check_merge_lossless.py; scripts/merged_vcf_to_gfa.py; tools/grvcf_to_vcf.py"],
], [63, 111], compact=True)
p("These are LinGraph grVCF conventions, not an independent VCF standard. Read the file's own header as its field schema and retain the producer version and reference identity when sharing data.", "SmallGr")


def main():
    doc = SimpleDocTemplate(
        str(DESTINATION), pagesize=A4, leftMargin=18 * mm,
        rightMargin=18 * mm, topMargin=23 * mm, bottomMargin=22 * mm,
        title="Graph Recursive VCF (grVCF): Lossless Assembly Information and Exact Cohort Variants",
        author="LinGraph", subject="Graph-aware VCF format guide",
    )
    doc.build(STORY, onFirstPage=header_footer, onLaterPages=header_footer)
    print(DESTINATION)


if __name__ == "__main__":
    main()
