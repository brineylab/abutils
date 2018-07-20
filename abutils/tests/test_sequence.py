#!usr/env/python
# filename: test_alignment.py

import os
from ..core.sequence import Sequence, read_fasta, read_json


def test_read_fasta():
    seqs = read_fasta('data/PG9.fasta')
    assert len(seqs) == 1 and seqs[0].id == 'PG9'


def test_read_json():
    seqs = read_json(os.path.abspath('tests/data/PG9.json'))
    assert len(seqs) == 1 and seqs[0].id == 'PG9'


