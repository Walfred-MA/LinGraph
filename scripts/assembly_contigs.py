"""Contig eligibility shared by cohort calling and graph export."""
HG38_MAIN_CONTIGS = frozenset([*(f'chr{i}' for i in range(1, 23)), 'chrX', 'chrY'])


def accepted_contig(sample, contig):
    if sample != 'HG38_h1':
        return True
    contig = contig or ''
    if contig.startswith('HG38#1#'):
        contig = contig[len('HG38#1#'):]
    return contig in HG38_MAIN_CONTIGS
