# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

abutils is a Python library for working with adaptive immune receptor repertoire (AIRR) data, including antibody and T-cell receptor sequences. It provides models, functions, and visualization tools for bulk AIRR datasets and serves as a core component of the ab[x] toolkit.

## Build and Test Commands

```bash
# Install dependencies
pip install -r requirements.txt

# Run all tests
pytest

# Run a specific test file
pytest abutils/tests/test_sequence.py

# Run a specific test
pytest abutils/tests/test_sequence.py::test_function_name -v
```

## Code Architecture

### Module Namespace Pattern (scanpy-style)
The library uses a scanpy-inspired namespace pattern with short module aliases:
- `abutils.tl` - Tools: clustering, alignment, phylogenetics, clonify (lineage assignment)
- `abutils.pl` - Plots: bar, scatter, heatmap, kde, ridge, donut visualizations
- `abutils.io` - I/O: reading/writing FASTA, FASTQ, CSV, AIRR, Parquet formats
- `abutils.pp` - Preprocessing: data preparation utilities
- `abutils.cl` / `abutils.color` - Color utilities and palettes

### Core Data Models (`abutils/core/`)
Three primary classes represent AIRR data:
- **Sequence**: Container for individual antibody/TCR sequences with annotations. Accepts raw strings, dicts, BioPython SeqRecords, or iterables. Annotations stored in `.annotations` dict are accessible via bracket notation (`seq['v_gene']`).
- **Pair**: Represents paired heavy/light (BCR) or alpha/beta (TCR) chains. Contains multiple Sequence objects with automatic chain type detection via locus or V-gene.
- **Lineage**: Groups related Pair objects into clonal lineages for phylogenetic analysis.

### Tools Module (`abutils/tools/`)
- `alignment.py`: Multiple sequence alignment (FAMSA, MAFFT, MUSCLE) and pairwise alignment (local, global, semiglobal via parasail)
- `cluster.py`: Sequence clustering with CD-HIT, MMseqs2, or VSEARCH backends
- `clonify.py`: Lineage/clonotype assignment using CDR3 similarity with shared mutation bonus
- `phylo.py`: Phylogenetic tree construction with FastTree, visualization with baltic

### Bundled Binaries (`abutils/binaries/`)
External tools are bundled for convenience: CD-HIT, FastTree, MAFFT, MUSCLE, MMseqs2, VSEARCH, minimap2, fastp. Binaries are platform-specific (darwin/linux, amd64/arm64) and downloaded on first use if not included in the package. Access via `abutils.bin.get_path("tool_name")`.

### I/O Patterns
The `io` module provides consistent read/write functions:
- `read_*` / `to_*` functions for each format (fasta, fastq, csv, airr, parquet)
- `from_polars` / `to_polars` convert between Sequence/Pair objects and Polars DataFrames
- Pair data uses column suffixes `:0` (heavy) and `:1` (light) in tabular formats

### DataFrame Support
Both Pandas and Polars are supported, with Polars preferred internally for performance. Lazy evaluation (`pl.LazyFrame`) is used where possible.

## Key Conventions

- AIRR format uses tab-separated values with standard field names (`sequence_id`, `sequence`, `v_call`, `j_call`, `cdr3`, etc.)
- Sequence IDs default to `sequence_id`, sequences to `sequence` field
- Test files are in `abutils/tests/` and follow `test_*.py` naming
- Python 3.10+ required
