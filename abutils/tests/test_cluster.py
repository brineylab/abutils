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


def test_cdhit():
    ofile, cfile = cluster.cdhit(CLUSTER_SEQS, make_db=False)
    files_exist = False
    if all([os.path.getsize(ofile), os.path.getsize(cfile)]):
        files_exist = True
    os.unlink(ofile)
    os.unlink(cfile)
    assert files_exist


def test_parse_clusters():
    seq_dict = {s.id: s for s in CLUSTER_SEQS}
    ofile, cfile = cluster.cdhit(CLUSTER_SEQS, make_db=False, threshold=0.9)
    cdhit_result = cluster.parse_clusters(ofile, cfile, seq_dict=seq_dict)
    assert len(cdhit_result) == 1


def test_cluster_centroid():
    cdhit_result = cluster.cluster(CLUSTER_SEQS, make_db=False, threshold=0.9)
    c = cdhit_result.clusters[0]
    assert c.centroid.sequence == CENTROID


def test_cluster_consensus():
    cdhit_result = cluster.cluster(CLUSTER_SEQS, make_db=False, threshold=0.9)
    c = cdhit_result.clusters[0]
    assert c.consensus.sequence == CONSENSUS
