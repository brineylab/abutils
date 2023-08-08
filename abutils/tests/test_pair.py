#!/usr/bin/env python
# filename: test_pair.py

import csv
import os
import pytest
import tempfile

from ..core.pair import assign_pairs, Pair
from ..core.sequence import Sequence


@pytest.fixture
def bcr_hk_seqs():
    return [
        Sequence({"locus": "IGH", "sequence": "AAAAA", "sequence_id": "seq1"}),
        Sequence({"locus": "IGK", "sequence": "CCCCC", "sequence_id": "seq2"}),
    ]


@pytest.fixture
def bcr_hl_seqs():
    return [
        Sequence({"locus": "IGH", "sequence": "AAAAA", "sequence_id": "seq1"}),
        Sequence({"locus": "IGL", "sequence": "CCCCC", "sequence_id": "seq2"}),
    ]


@pytest.fixture
def tcr_ab_seqs():
    return [
        Sequence({"locus": "TRA", "sequence": "AAAAA", "sequence_id": "seq1"}),
        Sequence({"locus": "TRB", "sequence": "CCCCC", "sequence_id": "seq2"}),
    ]


@pytest.fixture
def tcr_dg_seqs():
    return [
        Sequence({"locus": "TRD", "sequence": "AAAAA", "sequence_id": "seq1"}),
        Sequence({"locus": "TRG", "sequence": "CCCCC", "sequence_id": "seq2"}),
    ]


def test_init(bcr_hk_seqs, bcr_hl_seqs):
    pair = Pair(bcr_hk_seqs + bcr_hl_seqs)
    assert pair._seqs == bcr_hk_seqs + bcr_hl_seqs
    assert pair._receptor is None
    assert pair._heavy is None
    assert pair._light is None
    assert pair._alpha is None
    assert pair._beta is None
    assert pair._delta is None
    assert pair._gamma is None
    assert pair._heavies is None
    assert pair._lights is None
    assert pair._alphas is None
    assert pair._betas is None
    assert pair._deltas is None
    assert pair._gammas is None
    assert pair._name is None
    assert pair._fasta is None
    assert pair._sample is None
    assert pair._subject is None
    assert pair._group is None
    assert pair._experiment is None
    assert pair._timepoint is None
    assert pair._is_pair is None
    assert pair._lineage is None
    assert pair._select_chain == pair._chain_selector


def test_receptor():
    # with "chain" annotation
    p = Pair(
        [
            Sequence({"chain": "heavy", "sequence": "AAAAA", "sequence_id": "seq1"}),
            Sequence({"chain": "kappa", "sequence": "CCCCC", "sequence_id": "seq2"}),
        ]
    )
    assert p.receptor == "bcr"
    p = Pair(
        [
            Sequence({"chain": "heavy", "sequence": "AAAAA", "sequence_id": "seq1"}),
            Sequence({"chain": "lambda", "sequence": "CCCCC", "sequence_id": "seq2"}),
        ]
    )
    assert p.receptor == "bcr"
    p = Pair(
        [
            Sequence({"chain": "alpha", "sequence": "AAAAA", "sequence_id": "seq1"}),
            Sequence({"chain": "beta", "sequence": "CCCCC", "sequence_id": "seq2"}),
        ]
    )
    assert p.receptor == "tcr"
    p = Pair(
        [
            Sequence({"chain": "delta", "sequence": "AAAAA", "sequence_id": "seq1"}),
            Sequence({"chain": "gamma", "sequence": "CCCCC", "sequence_id": "seq2"}),
        ]
    )
    assert p.receptor == "tcr"
    p = Pair(
        [
            Sequence({"chain": "heavy", "sequence": "AAAAA", "sequence_id": "seq1"}),
            Sequence({"chain": "alpha", "sequence": "CCCCC", "sequence_id": "seq2"}),
        ]
    )
    assert p.receptor == "unknown"
    # with "locus" annotation
    p = Pair(
        [
            Sequence({"locus": "IGH", "sequence": "AAAAA", "sequence_id": "seq1"}),
            Sequence({"locus": "IGK", "sequence": "CCCCC", "sequence_id": "seq2"}),
        ]
    )
    assert p.receptor == "bcr"
    p = Pair(
        [
            Sequence({"locus": "IGH", "sequence": "AAAAA", "sequence_id": "seq1"}),
            Sequence({"locus": "IGL", "sequence": "CCCCC", "sequence_id": "seq2"}),
        ]
    )
    assert p.receptor == "bcr"
    p = Pair(
        [
            Sequence({"locus": "TRA", "sequence": "AAAAA", "sequence_id": "seq1"}),
            Sequence({"locus": "TRB", "sequence": "CCCCC", "sequence_id": "seq2"}),
        ]
    )
    assert p.receptor == "tcr"
    p = Pair(
        [
            Sequence({"locus": "TRD", "sequence": "AAAAA", "sequence_id": "seq1"}),
            Sequence({"locus": "TRA", "sequence": "CCCCC", "sequence_id": "seq2"}),
        ]
    )
    assert p.receptor == "tcr"
    p = Pair(
        [
            Sequence({"locus": "IGH", "sequence": "AAAAA", "sequence_id": "seq1"}),
            Sequence({"locus": "TRA", "sequence": "CCCCC", "sequence_id": "seq2"}),
        ]
    )
    assert p.receptor == "unknown"


def test_heavy(bcr_hk_seqs, bcr_hl_seqs):
    # heavy/kappa pair
    pair = Pair(bcr_hk_seqs)
    assert pair.heavy is not None
    assert pair.heavy["sequence_id"] == "seq1"
    assert pair.heavy["sequence"] == "AAAAA"
    assert pair.heavy["locus"] == "IGH"
    # heavy/lambda pair
    pair = Pair(bcr_hl_seqs)
    assert pair.heavy is not None
    assert pair.heavy["sequence_id"] == "seq1"
    assert pair.heavy["sequence"] == "AAAAA"
    assert pair.heavy["locus"] == "IGH"


def test_light(bcr_hk_seqs, bcr_hl_seqs):
    # heavy/kappa pair
    pair = Pair(bcr_hk_seqs)
    assert pair.light is not None
    assert pair.light["sequence_id"] == "seq2"
    assert pair.light["sequence"] == "CCCCC"
    assert pair.light["locus"] == "IGK"
    # heavy/lambda pair
    pair = Pair(bcr_hl_seqs)
    assert pair.light is not None
    assert pair.light["sequence_id"] == "seq2"
    assert pair.light["sequence"] == "CCCCC"
    assert pair.light["locus"] == "IGL"


def test_alpha(tcr_ab_seqs):
    pair = Pair(tcr_ab_seqs)
    assert pair.alpha is not None
    assert pair.alpha["sequence_id"] == "seq1"
    assert pair.alpha["sequence"] == "AAAAA"
    assert pair.alpha["locus"] == "TRA"


def test_beta(tcr_ab_seqs):
    pair = Pair(tcr_ab_seqs)
    assert pair.beta is not None
    assert pair.beta["sequence_id"] == "seq2"
    assert pair.beta["sequence"] == "CCCCC"
    assert pair.beta["locus"] == "TRB"


def test_delta(tcr_dg_seqs):
    pair = Pair(tcr_dg_seqs)
    assert pair.delta is not None
    assert pair.delta["sequence_id"] == "seq1"
    assert pair.delta["sequence"] == "AAAAA"
    assert pair.delta["locus"] == "TRD"


def test_gamma(tcr_dg_seqs):
    pair = Pair(tcr_dg_seqs)
    assert pair.gamma is not None
    assert pair.gamma["sequence_id"] == "seq2"
    assert pair.gamma["sequence"] == "CCCCC"
    assert pair.gamma["locus"] == "TRG"


def test_is_pair(bcr_hk_seqs, bcr_hl_seqs, tcr_ab_seqs, tcr_dg_seqs):
    # actual pairs
    pair = Pair(bcr_hk_seqs)
    assert pair.is_pair is True
    pair = Pair(bcr_hl_seqs)
    assert pair.is_pair is True
    pair = Pair(tcr_ab_seqs)
    assert pair.is_pair is True
    pair = Pair(tcr_dg_seqs)
    assert pair.is_pair is True
    # mixed pairs
    pair = Pair(
        [
            Sequence({"locus": "IGH", "sequence": "AAAAA", "sequence_id": "seq1"}),
            Sequence({"locus": "TRA", "sequence": "CCCCC", "sequence_id": "seq2"}),
        ]
    )
    assert pair.is_pair is False
    pair = Pair(
        [
            Sequence({"locus": "TRD", "sequence": "AAAAA", "sequence_id": "seq1"}),
            Sequence({"locus": "TRA", "sequence": "CCCCC", "sequence_id": "seq2"}),
        ]
    )
    assert pair.is_pair is False
    # singletons
    pair = Pair(
        [
            Sequence({"locus": "IGH", "sequence": "AAAAA", "sequence_id": "seq1"}),
        ]
    )
    assert pair.is_pair is False
    pair = Pair(
        [
            Sequence({"locus": "TRB", "sequence": "AAAAA", "sequence_id": "seq1"}),
        ]
    )
    assert pair.is_pair is False


def test_name():
    pair = Pair(
        [
            Sequence({"locus": "IGH", "sequence": "AAAAA", "sequence_id": "bcr-pair"}),
            Sequence({"locus": "IGK", "sequence": "CCCCC", "sequence_id": "bcr-pair"}),
        ]
    )
    assert pair.name == "bcr-pair"
    pair = Pair(
        [
            Sequence({"locus": "TRA", "sequence": "AAAAA", "sequence_id": "tcr-pair"}),
            Sequence({"locus": "TRB", "sequence": "CCCCC", "sequence_id": "tcr-pair"}),
        ]
    )
    assert pair.name == "tcr-pair"


def test_select_chain():
    def umi_selector(seq_list):
        if len(seq_list) == 1:
            return seq_list[0]
        sorted_seqs = sorted(seq_list, key=lambda x: x["umi"], reverse=True)
        return sorted_seqs[0]

    pair = Pair(
        [
            Sequence(
                {
                    "locus": "IGH",
                    "sequence": "AAAAA",
                    "sequence_id": "heavy1",
                    "umi": 100,
                }
            ),
            Sequence(
                {
                    "locus": "IGH",
                    "sequence": "AAAAA",
                    "sequence_id": "heavy2",
                    "umi": 10,
                }
            ),
            Sequence(
                {
                    "locus": "IGK",
                    "sequence": "CCCCC",
                    "sequence_id": "kappa1",
                    "umi": 20,
                }
            ),
            Sequence(
                {
                    "locus": "IGK",
                    "sequence": "CCCCC",
                    "sequence_id": "kappa2",
                    "umi": 200,
                }
            ),
        ],
        chain_selection_func=umi_selector,
    )
    assert len(pair.heavies) == 2
    assert pair.heavy is not None
    assert pair.heavy.id == "heavy1"
    assert len(pair.lights) == 2
    assert pair.light is not None
    assert pair.light.id == "kappa2"


# -------------------
#    assign_pairs
# -------------------


def test_assign_pairs_with_no_delimiter():
    seqs = [
        Sequence({"sequence_id": "seq1", "sequence": "ATCG", "locus": "IGH"}),
        Sequence({"sequence_id": "seq1", "sequence": "ATCG", "locus": "IGK"}),
        Sequence({"sequence_id": "seq2", "sequence": "ATCG", "locus": "IGH"}),
        Sequence({"sequence_id": "seq2", "sequence": "ATCG", "locus": "IGL"}),
        Sequence({"sequence_id": "seq3", "sequence": "ATCG", "locus": "TRA"}),
        Sequence({"sequence_id": "seq3", "sequence": "ATCG", "locus": "TRB"}),
        Sequence({"sequence_id": "seq4", "sequence": "ATCG", "locus": "TRD"}),
        Sequence({"sequence_id": "seq4", "sequence": "ATCG", "locus": "TRG"}),
    ]
    pairs = assign_pairs(seqs)
    assert isinstance(pairs, list)
    assert len(pairs) == 4
    for pair in pairs:
        assert pair.name in ["seq1", "seq2", "seq3", "seq4"]


def test_assign_pairs_with_delimiter():
    seqs = [
        Sequence({"sequence_id": "seq1_heavy", "sequence": "ATCG", "locus": "IGH"}),
        Sequence({"sequence_id": "seq1_kappa", "sequence": "ATCG", "locus": "IGK"}),
        Sequence({"sequence_id": "seq2_heavy", "sequence": "ATCG", "locus": "IGH"}),
        Sequence({"sequence_id": "seq2_lambda", "sequence": "ATCG", "locus": "IGL"}),
        Sequence({"sequence_id": "seq3_alpha", "sequence": "ATCG", "locus": "TRA"}),
        Sequence({"sequence_id": "seq3_beta", "sequence": "ATCG", "locus": "TRB"}),
        Sequence({"sequence_id": "seq4_delta", "sequence": "ATCG", "locus": "TRD"}),
        Sequence({"sequence_id": "seq4_gamma", "sequence": "ATCG", "locus": "TRG"}),
    ]
    pairs = assign_pairs(seqs, delim="_", delim_occurance=1)
    assert isinstance(pairs, list)
    assert len(pairs) == 4
    for pair in pairs:
        assert pair.name in ["seq1", "seq2", "seq3", "seq4"]


def test_assign_pairs_with_pairs_only():
    seqs = [
        Sequence({"sequence_id": "seq1", "sequence": "ATCG", "locus": "IGH"}),
        Sequence({"sequence_id": "seq1", "sequence": "ATCG", "locus": "IGK"}),
        Sequence({"sequence_id": "seq2", "sequence": "ATCG", "locus": "IGH"}),
        Sequence({"sequence_id": "seq2", "sequence": "ATCG", "locus": "IGL"}),
        Sequence({"sequence_id": "seq3", "sequence": "ATCG", "locus": "TRA"}),
        Sequence({"sequence_id": "seq3", "sequence": "ATCG", "locus": "TRB"}),
        Sequence({"sequence_id": "seq4", "sequence": "ATCG", "locus": "TRD"}),
    ]
    pairs = assign_pairs(seqs, pairs_only=True)
    assert isinstance(pairs, list)
    assert len(pairs) == 3
    for pair in pairs:
        assert pair.name in ["seq1", "seq2", "seq3"]


def test_assign_pairs_with_tenx_annot_file():
    tenx_annot_file = tempfile.NamedTemporaryFile(delete=False)
    tenx_annot_file.close()
    with open(tenx_annot_file.name, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["contig_id", "umis"])
        writer.writerow(["seq1_heavy", "10"])
        writer.writerow(["seq1_kappa", "100"])
        writer.writerow(["seq2_heavy", "20"])
        writer.writerow(["seq2_kappa", "200"])

    seqs = [
        Sequence({"sequence_id": "seq1_heavy", "sequence": "ATCG", "locus": "IGH"}),
        Sequence({"sequence_id": "seq1_kappa", "sequence": "ATCG", "locus": "IGK"}),
        Sequence({"sequence_id": "seq2_heavy", "sequence": "ATCG", "locus": "IGH"}),
        Sequence({"sequence_id": "seq2_kappa", "sequence": "ATCG", "locus": "IGK"}),
    ]
    pairs = assign_pairs(
        seqs, tenx_annot_file=tenx_annot_file.name, delim="_", delim_occurance=1
    )
    assert isinstance(pairs, list)
    assert len(pairs) == 2
    for pair in pairs:
        assert isinstance(pair, Pair)
        if pair.name == "seq1":
            assert pair.heavy["umis"] == 10
            assert pair.light["umis"] == 100
        if pair.name == "seq2":
            assert pair.heavy["umis"] == 20
            assert pair.light["umis"] == 200
    os.unlink(tenx_annot_file.name)
