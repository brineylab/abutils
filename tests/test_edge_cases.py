"""Edge case tests for core data models: Sequence, Pair, Lineage."""

import pytest

from abutils.core.lineage import Lineage
from abutils.core.pair import Pair
from abutils.core.sequence import Sequence


# ============================
#   Sequence edge cases
# ============================


class TestSequenceEdgeCases:
    def test_empty_string_sequence(self):
        seq = Sequence("")
        assert seq.sequence == ""
        assert len(seq) == 0

    def test_none_id_generates_uuid(self):
        seq = Sequence("ATCG")
        assert seq.id is not None
        assert len(seq.id) > 0

    def test_dict_with_missing_id_key(self):
        seq = Sequence({"sequence": "ATCG"})
        assert seq.sequence == "ATCG"
        # id should be None when no id key is found
        assert seq.id is None

    def test_dict_with_missing_sequence_key(self):
        seq = Sequence({"sequence_id": "seq1"})
        assert seq.id == "seq1"
        assert seq.sequence is None

    def test_getitem_with_annotation_key(self):
        seq = Sequence({"sequence_id": "seq1", "sequence": "ATCG", "v_call": "IGHV1-2*02"})
        assert seq["v_call"] == "IGHV1-2*02"

    def test_getitem_missing_key_returns_none(self):
        seq = Sequence({"sequence_id": "seq1", "sequence": "ATCG"})
        assert seq["nonexistent_key"] is None

    def test_contains_with_dict_checks_keys(self):
        seq = Sequence({"sequence_id": "seq1", "sequence": "ATCG"})
        assert "sequence_id" in seq
        assert "nonexistent" not in seq

    def test_contains_with_string_checks_motif(self):
        seq = Sequence("ATCGATCG")
        assert "ATCG" in seq
        assert "GGGG" not in seq

    def test_sequence_length(self):
        seq = Sequence("ATCGATCG")
        assert len(seq) == 8

    def test_sequence_equality(self):
        seq1 = Sequence("ATCG", id="seq1")
        seq2 = Sequence("ATCG", id="seq1")
        seq3 = Sequence("ATCG", id="seq2")
        assert seq1 == seq2
        assert seq1 != seq3

    def test_reverse_complement_palindrome(self):
        seq = Sequence("AATT")
        assert seq.reverse_complement == "AATT"

    def test_reverse_complement_single_base(self):
        seq = Sequence("A")
        assert seq.reverse_complement == "T"

    def test_translate_stop_codon(self):
        seq = Sequence("TAA")
        assert seq.translate() == "*"

    def test_translate_incomplete_codon(self):
        """Translating a sequence not divisible by 3."""
        seq = Sequence("ATCGA")
        # should handle gracefully (translate the full codons)
        result = seq.translate()
        assert isinstance(result, str)

    def test_fasta_property(self):
        seq = Sequence("ATCG", id="test_seq")
        assert seq.fasta == ">test_seq\nATCG"

    def test_fasta_property_cached(self):
        seq = Sequence("ATCG", id="test_seq")
        fasta1 = seq.fasta
        fasta2 = seq.fasta
        assert fasta1 is fasta2  # same object (cached)

    def test_annotations_dict_access(self):
        data = {
            "sequence_id": "PG9",
            "sequence": "ATCG",
            "v_call": "IGHV1-2*02",
            "j_call": "IGHJ4*02",
            "cdr3_aa": "CARDRGSGYYDFWSG",
        }
        seq = Sequence(data)
        assert seq["v_call"] == "IGHV1-2*02"
        assert seq["j_call"] == "IGHJ4*02"
        assert seq["cdr3_aa"] == "CARDRGSGYYDFWSG"

    def test_sequence_from_seqrecord(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        rec = SeqRecord(Seq("ATCGATCG"), id="bio_seq", description="test record")
        seq = Sequence(rec)
        assert seq.id == "bio_seq"
        assert seq.sequence == "ATCGATCG"

    def test_sequence_from_another_sequence(self):
        orig = Sequence("ATCG", id="original")
        copy = Sequence(orig)
        assert copy.sequence == "ATCG"
        assert copy.id == "original"

    def test_sequence_from_list(self):
        seq = Sequence(["my_id", "GCTAGCTA"])
        assert seq.id == "my_id"
        assert seq.sequence == "GCTAGCTA"


# ============================
#   Sequence with real data
# ============================


class TestSequenceWithRealData:
    """Tests using the PG9 antibody data from conftest."""

    def test_pg9_has_vdj_sequence(self, pg9_sequence):
        assert pg9_sequence.sequence is not None
        assert len(pg9_sequence.sequence) > 200

    def test_pg9_annotation_access(self, pg9_sequence):
        assert pg9_sequence["chain"] == "lambda"
        assert pg9_sequence["prod"] == "yes"

    def test_pg9_v_gene_info(self, pg9_sequence):
        v_gene = pg9_sequence["v_gene"]
        assert v_gene["full"] == "IGLV2-14*01"
        assert v_gene["fam"] == "IGLV2"
        assert v_gene["gene"] == "IGLV2-14"

    def test_pg9_j_gene_info(self, pg9_sequence):
        j_gene = pg9_sequence["j_gene"]
        assert j_gene["full"] == "IGLJ3*02"
        assert j_gene["gene"] == "IGLJ3"

    def test_pg9_cdr3(self, pg9_sequence):
        assert pg9_sequence["cdr3_aa"] == "KSLTSTRRRV"
        assert pg9_sequence["cdr3_len"] == 10

    def test_pg9_mutation_counts(self, pg9_sequence):
        assert pg9_sequence["mut_count_nt"] == 31
        assert pg9_sequence["var_muts_nt"]["num"] == 25
        assert pg9_sequence["join_muts_nt"]["num"] == 6

    def test_pg9_identity(self, pg9_sequence):
        nt_ident = pg9_sequence["nt_identity"]
        assert 0 < nt_ident["v"] < 100
        assert 0 < nt_ident["j"] < 100

    def test_pg9_regions(self, pg9_sequence):
        """Verify all FR/CDR regions are present."""
        for region in ["fr1_aa", "cdr1_aa", "fr2_aa", "cdr2_aa", "fr3_aa", "cdr3_aa", "fr4_aa"]:
            assert pg9_sequence[region] is not None
            assert len(pg9_sequence[region]) > 0


# ============================
#   Pair edge cases
# ============================


class TestPairEdgeCases:
    def test_heavy_only_pair(self):
        heavy = Sequence({"locus": "IGH", "sequence": "ATCG", "sequence_id": "h1"})
        pair = Pair([heavy])
        assert pair.heavy is not None
        assert pair.light is None
        assert pair.is_pair is False

    def test_light_only_pair(self):
        light = Sequence({"locus": "IGK", "sequence": "GCTA", "sequence_id": "l1"})
        pair = Pair([light])
        assert pair.heavy is None
        assert pair.light is not None
        assert pair.is_pair is False

    def test_bcr_pair_detection(self):
        h = Sequence({"locus": "IGH", "sequence": "ATCG", "sequence_id": "h"})
        l = Sequence({"locus": "IGK", "sequence": "GCTA", "sequence_id": "l"})
        pair = Pair([h, l])
        assert pair.receptor == "bcr"
        assert pair.is_pair is True

    def test_tcr_pair_detection(self):
        a = Sequence({"locus": "TRA", "sequence": "ATCG", "sequence_id": "a"})
        b = Sequence({"locus": "TRB", "sequence": "GCTA", "sequence_id": "b"})
        pair = Pair([a, b])
        assert pair.receptor == "tcr"
        assert pair.is_pair is True

    def test_tcr_gamma_delta_detection(self):
        g = Sequence({"locus": "TRG", "sequence": "ATCG", "sequence_id": "g"})
        d = Sequence({"locus": "TRD", "sequence": "GCTA", "sequence_id": "d"})
        pair = Pair([g, d])
        assert pair.receptor == "tcr"
        assert pair.is_pair is True
        assert pair.gamma is not None
        assert pair.delta is not None

    def test_locus_from_v_call(self):
        """When locus field is missing, infer from v_call."""
        h = Sequence({"v_call": "IGHV1-2*02", "sequence": "ATCG", "sequence_id": "h"})
        l = Sequence({"v_call": "IGKV3-20*01", "sequence": "GCTA", "sequence_id": "l"})
        pair = Pair([h, l])
        assert pair.receptor == "bcr"

    def test_pair_name_from_heavy(self):
        h = Sequence({"locus": "IGH", "sequence": "ATCG", "sequence_id": "my_heavy"})
        l = Sequence({"locus": "IGK", "sequence": "GCTA", "sequence_id": "my_light"})
        pair = Pair([h, l])
        assert pair.name == "my_heavy"

    def test_pair_name_from_light_when_no_heavy(self):
        l = Sequence({"locus": "IGK", "sequence": "GCTA", "sequence_id": "my_light"})
        pair = Pair([l])
        assert pair.name == "my_light"

    def test_pair_explicit_name(self):
        h = Sequence({"locus": "IGH", "sequence": "ATCG", "sequence_id": "h"})
        pair = Pair([h], name="custom_name")
        assert pair.name == "custom_name"

    def test_heavies_property(self):
        h1 = Sequence({"locus": "IGH", "sequence": "ATCG", "sequence_id": "h1"})
        h2 = Sequence({"locus": "IGH", "sequence": "GCTA", "sequence_id": "h2"})
        l = Sequence({"locus": "IGK", "sequence": "TTTT", "sequence_id": "l1"})
        pair = Pair([h1, h2, l])
        assert len(pair.heavies) == 2
        assert len(pair.lights) == 1

    def test_lights_property_kappa_and_lambda(self):
        k = Sequence({"locus": "IGK", "sequence": "AAAA", "sequence_id": "k"})
        l = Sequence({"locus": "IGL", "sequence": "CCCC", "sequence_id": "l"})
        pair = Pair([k, l])
        assert len(pair.lights) == 2

    def test_fasta_output(self):
        h = Sequence({"locus": "IGH", "sequence": "ATCGATCG", "sequence_id": "h1"})
        l = Sequence({"locus": "IGK", "sequence": "GCTAGCTA", "sequence_id": "l1"})
        pair = Pair([h, l])
        fasta = pair.fasta()
        assert ">h1_heavy" in fasta
        assert "ATCGATCG" in fasta
        assert ">l1_light" in fasta
        assert "GCTAGCTA" in fasta

    def test_fasta_no_chain_append(self):
        h = Sequence({"locus": "IGH", "sequence": "ATCG", "sequence_id": "h1"})
        pair = Pair([h])
        fasta = pair.fasta(append_chain=False)
        assert ">h1\n" in fasta
        assert "_heavy" not in fasta

    def test_chain_selection_with_umis(self):
        """When multiple heavy chains exist, select by highest UMI count."""
        h1 = Sequence(
            {"locus": "IGH", "sequence": "AAAA", "sequence_id": "h1", "umis": 5}
        )
        h2 = Sequence(
            {"locus": "IGH", "sequence": "CCCC", "sequence_id": "h2", "umis": 10}
        )
        pair = Pair([h1, h2])
        assert pair.heavy.id == "h2"  # higher UMI count


# ============================
#   Pair with real data
# ============================


class TestPairWithRealData:
    def test_bcr_pairs_are_valid(self, bcr_pairs):
        assert len(bcr_pairs) == 5
        for pair in bcr_pairs:
            assert pair.is_pair is True
            assert pair.receptor == "bcr"
            assert pair.heavy is not None
            assert pair.light is not None

    def test_bcr_pair_v_genes(self, bcr_pairs):
        """Verify V gene annotations are accessible through pairs."""
        for pair in bcr_pairs:
            assert pair.heavy["v_call"] is not None
            assert pair.heavy["v_call"].startswith("IGH")
            assert pair.light["v_call"] is not None
            assert pair.light["v_call"][:3] in ("IGK", "IGL")

    def test_bcr_pair_cdr3(self, bcr_pairs):
        for pair in bcr_pairs:
            assert pair.heavy["cdr3_aa"] is not None
            assert len(pair.heavy["cdr3_aa"]) > 0

    def test_tcr_pairs_are_valid(self, tcr_pairs):
        assert len(tcr_pairs) == 3
        for pair in tcr_pairs:
            assert pair.is_pair is True
            assert pair.receptor == "tcr"
            assert pair.alpha is not None
            assert pair.beta is not None

    def test_tcr_pair_v_genes(self, tcr_pairs):
        for pair in tcr_pairs:
            assert pair.alpha["v_call"].startswith("TRA")
            assert pair.beta["v_call"].startswith("TRB")


# ============================
#   Lineage edge cases
# ============================


class TestLineageEdgeCases:
    def test_single_pair_lineage(self):
        h = Sequence({"locus": "IGH", "sequence": "ATCG", "sequence_id": "h1"})
        l = Sequence({"locus": "IGK", "sequence": "GCTA", "sequence_id": "l1"})
        lineage = Lineage([Pair([h, l])])
        assert lineage.size() == 1

    def test_heavy_only_lineage(self):
        pairs = [
            Pair([Sequence({"locus": "IGH", "sequence": "ATCG", "sequence_id": f"h{i}"})])
            for i in range(3)
        ]
        lineage = Lineage(pairs)
        assert lineage.size() == 3
        assert lineage.size(pairs_only=True) == 0
        assert len(lineage.heavies) == 3
        assert len(lineage.lights) == 0

    def test_lineage_iteration(self):
        pairs = [
            Pair(
                [
                    Sequence({"locus": "IGH", "sequence": "ATCG", "sequence_id": f"h{i}"}),
                    Sequence({"locus": "IGK", "sequence": "GCTA", "sequence_id": f"l{i}"}),
                ]
            )
            for i in range(4)
        ]
        lineage = Lineage(pairs)
        collected = list(lineage)
        assert len(collected) == 4

    def test_lineage_contains_and_getitem(self):
        h = Sequence({"locus": "IGH", "sequence": "ATCG", "sequence_id": "h1"})
        l = Sequence({"locus": "IGK", "sequence": "GCTA", "sequence_id": "l1"})
        pair = Pair([h, l])
        lineage = Lineage([pair])
        assert pair.name in lineage
        retrieved = lineage[pair.name]
        assert retrieved.name == pair.name


# ============================
#   Lineage with real data
# ============================


class TestLineageWithRealData:
    def test_lineage_from_bcr_pairs(self, bcr_lineage):
        assert bcr_lineage.size() == 3
        assert bcr_lineage.size(pairs_only=True) == 3

    def test_lineage_heavies_have_correct_vgene(self, bcr_lineage):
        for pair in bcr_lineage.heavies:
            assert pair.heavy["v_call"] == "IGHV1-2*02"

    def test_lineage_just_pairs_all_paired(self, bcr_lineage):
        for pair in bcr_lineage.just_pairs:
            assert pair.is_pair is True

    def test_lineage_indel_flags(self, bcr_lineage):
        # These sequences don't have indel annotations
        assert bcr_lineage.has_insertion is False
        assert bcr_lineage.has_deletion is False
        assert bcr_lineage.has_indel is False
