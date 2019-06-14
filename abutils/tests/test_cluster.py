#!usr/env/python
# filename: test_cluster.py

import os

from ..core.sequence import Sequence
from ..utils import cluster


CLUSTER_SEQS = [Sequence('ATGCATGCATGCATGCATGCATGC', id='seq1'),
                Sequence('ATGCATGCATGCATGCATGGATGC', id='seq2'),
                Sequence('ATGCATCCATGCATGCATGCATGC', id='seq3')]
CONSENSUS = 'ATGCATGCATGCATGCATGCATGC'
CENTROID = 'ATGCATGCATGCATGCATGCATGC'
TEMP = os.path.abspath('./')


# def test_cdhit():
#     ofile, cfile = cluster.cdhit(CLUSTER_SEQS, make_db=False, debug=True)
#     errors = []
#     if not os.path.getsize(ofile):
#         errors.append('CDHIT output file ({}) is empty'.format(ofile))
#     if not os.path.getsize(cfile):
#         errors.append('CDHIT cluster file ({}) is empty'.format(cfile))
#     os.unlink(ofile)
#     os.unlink(cfile)
#     assert len(errors) == 0, '\n'.join(errors)


# def test_parse_clusters():
#     seq_dict = {s.id: s for s in CLUSTER_SEQS}
#     ofile, cfile = cluster.cdhit(CLUSTER_SEQS, make_db=False, threshold=0.9, debug=True)
#     cdhit_result = cluster.parse_clusters(ofile, cfile, seq_dict=seq_dict)
#     assert len(cdhit_result) == 1


def test_cluster_centroid():
    clu = cluster.cluster(CLUSTER_SEQS, threshold=0.9, debug=True)
    c = clu[0]
    assert c.centroid.sequence == CENTROID


def test_cluster_consensus():
    clu = cluster.cluster(CLUSTER_SEQS, threshold=0.9, debug=True)
    c = clu[0]
    assert c.consensus.sequence == CONSENSUS
