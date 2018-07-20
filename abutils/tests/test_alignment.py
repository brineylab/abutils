#!usr/env/python
# filename: test_alignment.py

import os
from tempfile import NamedTemporaryFile

from Bio import SeqIO

from ..core.sequence import Sequence
from ..utils import alignment


QUERY = Sequence('ATGCATGC', id='query')
ALN_QUERY = 'ATGCATGC'
TARGET = Sequence('ATGATGC', id='target')
ALN_TARGET = 'ATG-ATGC'


def test_muscle():
    aln = alignment.muscle(sequences=[QUERY, TARGET])
    aln_query = [a for a in aln if a.id == 'query'][0]
    aln_target = [a for a in aln if a.id == 'target'][0]
    assert str(aln_query.seq) == ALN_QUERY and str(aln_target.seq) == ALN_TARGET


def test_muscle_as_file():
    aln_file = NamedTemporaryFile(delete=False)
    aln_file.close()
    aln = alignment.muscle(sequences=[QUERY, TARGET], as_file=True, alignment_file=aln_file.name)
    with open(aln) as f:
        seqs = [s for s in SeqIO.parse(f, 'fasta')]
    os.unlink(aln)
    aln_query = [a for a in seqs if a.id == 'query'][0]
    aln_target = [a for a in seqs if a.id == 'target'][0]
    assert str(aln_query.seq) == ALN_QUERY and str(aln_target.seq) == ALN_TARGET


# def test_mafft():
#     aln = alignment.mafft(sequences=[QUERY, TARGET], print_stdout=True, print_stderr=True)
#     aln_query = [a for a in aln if a.id == 'query'][0]
#     aln_target = [a for a in aln if a.id == 'target'][0]
#     assert str(aln_query.seq).upper() == '---ATGCATGC' and str(aln_target.seq).upper() == 'ATGATGC----'


def test_local_alignment():
    aln = alignment.local_alignment(QUERY, TARGET, gap_open=-1)
    assert aln.aligned_query == ALN_QUERY and aln.aligned_target == ALN_TARGET


def test_global_alignment():
    aln = alignment.global_alignment(QUERY, TARGET, gap_open=-1)
    assert aln.aligned_query == ALN_QUERY and aln.aligned_target == ALN_TARGET
