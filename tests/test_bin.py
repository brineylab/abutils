"""Tests for abutils.bin (binary path resolution)."""

import os
from unittest import mock

import pytest

from abutils.bin import BIN_DIR, MACHINE, SYSTEM, get_binary_directory, get_path


class TestGetBinaryDirectory:
    def test_returns_string(self):
        result = get_binary_directory()
        assert isinstance(result, str)

    def test_path_contains_binaries(self):
        result = get_binary_directory()
        assert result.endswith("binaries")

    def test_directory_exists(self):
        result = get_binary_directory()
        assert os.path.isdir(result)


class TestGetPath:
    def test_invalid_binary_raises(self):
        with pytest.raises(ValueError, match="not available"):
            get_path("nonexistent_tool")

    def test_name_normalization_mmseqs2(self):
        """mmseqs2 should be normalized to mmseqs."""
        with mock.patch("os.path.exists", return_value=True):
            path = get_path("mmseqs2")
        assert "mmseqs_" in os.path.basename(path)
        assert "mmseqs2_" not in os.path.basename(path)

    def test_name_normalization_cd_hit(self):
        """cd-hit should be normalized to cdhit."""
        with mock.patch("os.path.exists", return_value=True):
            path = get_path("cd-hit")
        assert "cdhit_" in os.path.basename(path)

    def test_name_normalization_muscle_v3(self):
        """muscle_v3 should be normalized to muscle3."""
        with mock.patch("os.path.exists", return_value=True):
            path = get_path("muscle_v3")
        assert "muscle3_" in os.path.basename(path)

    def test_name_normalization_minimap(self):
        """minimap should be normalized to minimap2."""
        with mock.patch("os.path.exists", return_value=True):
            path = get_path("minimap")
        assert "minimap2_" in os.path.basename(path)

    def test_case_insensitive(self):
        with mock.patch("os.path.exists", return_value=True):
            path = get_path("FASTTREE")
        assert "fasttree_" in os.path.basename(path)

    def test_system_only_binaries(self):
        """fastp, muscle3, blastn, makeblastdb use system-only paths (no machine arch)."""
        for binary in ["fastp", "muscle_v3", "blastn", "makeblastdb"]:
            with mock.patch("os.path.exists", return_value=True):
                path = get_path(binary)
            basename = os.path.basename(path)
            assert f"_{SYSTEM}" in basename
            # Should NOT have machine arch suffix for these
            assert not basename.endswith(f"_{MACHINE}")

    def test_system_machine_binaries(self):
        """Most binaries include both system and machine in the path."""
        for binary in ["cdhit", "fasttree", "mmseqs", "muscle", "vsearch", "minimap2"]:
            with mock.patch("os.path.exists", return_value=True):
                path = get_path(binary)
            basename = os.path.basename(path)
            assert f"_{SYSTEM}_{MACHINE}" in basename

    def test_mafft_has_subdirectory(self):
        """mafft binary lives in a subdirectory with mafft.bat."""
        with mock.patch("os.path.exists", return_value=True):
            path = get_path("mafft")
        assert path.endswith("mafft.bat")
        assert f"mafft_{SYSTEM}" in path

    @pytest.mark.parametrize(
        "name",
        [
            "blastn",
            "makeblastdb",
            "cdhit",
            "mafft",
            "minimap2",
            "mmseqs",
            "muscle",
            "muscle_v3",
            "vsearch",
            "fastp",
            "sickle",
            "fasttree",
        ],
    )
    def test_all_available_binaries_resolve(self, name):
        """Every documented binary name resolves without ValueError."""
        with mock.patch("os.path.exists", return_value=True):
            path = get_path(name)
        assert path.startswith(BIN_DIR)
