"""Tests for abutils.utils.codons."""

from abutils.utils.codons import codon_lookup


class TestCodonLookup:
    def test_has_64_codons(self):
        assert len(codon_lookup) == 64

    def test_all_codons_present(self):
        bases = "ACGT"
        for b1 in bases:
            for b2 in bases:
                for b3 in bases:
                    assert f"{b1}{b2}{b3}" in codon_lookup

    def test_stop_codons(self):
        assert codon_lookup["TAA"] == "*"
        assert codon_lookup["TAG"] == "*"
        assert codon_lookup["TGA"] == "*"

    def test_methionine_start(self):
        assert codon_lookup["ATG"] == "M"

    def test_tryptophan_unique(self):
        assert codon_lookup["TGG"] == "W"

    def test_all_20_amino_acids_represented(self):
        amino_acids = set(codon_lookup.values()) - {"*"}
        assert len(amino_acids) == 20
        expected = set("ACDEFGHIKLMNPQRSTVWY")
        assert amino_acids == expected

    def test_degenerate_codons(self):
        """Multiple codons should encode the same amino acid."""
        # Leucine has 6 codons
        leu_codons = [k for k, v in codon_lookup.items() if v == "L"]
        assert len(leu_codons) == 6
        # Serine has 6 codons
        ser_codons = [k for k, v in codon_lookup.items() if v == "S"]
        assert len(ser_codons) == 6
        # Isoleucine has 3 codons
        ile_codons = [k for k, v in codon_lookup.items() if v == "I"]
        assert len(ile_codons) == 3

    def test_specific_codons(self):
        """Spot-check several codons against the standard genetic code."""
        assert codon_lookup["TTT"] == "F"
        assert codon_lookup["GGG"] == "G"
        assert codon_lookup["AAA"] == "K"
        assert codon_lookup["CCC"] == "P"
        assert codon_lookup["GAT"] == "D"
        assert codon_lookup["GAG"] == "E"
        assert codon_lookup["CGT"] == "R"
        assert codon_lookup["ACT"] == "T"
        assert codon_lookup["GTG"] == "V"
        assert codon_lookup["GCT"] == "A"
        assert codon_lookup["CAT"] == "H"
        assert codon_lookup["CAG"] == "Q"
        assert codon_lookup["AAT"] == "N"
        assert codon_lookup["TGT"] == "C"
        assert codon_lookup["TAT"] == "Y"
