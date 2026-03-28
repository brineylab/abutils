# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

abutils is a Python library for working with adaptive immune receptor repertoire (AIRR) data, including antibody and T-cell receptor sequences. It provides models, functions, and visualization tools for bulk AIRR datasets and serves as a core component of the ab[x] toolkit. Used by [`scab`](https://github.com/briney/scab) for single cell AIRR analysis.

## Build and Test Commands

```bash
# Install for development (editable install)
pip install -e .

# Install dependencies only
pip install -r requirements.txt

# Run all tests
pytest

# Run a specific test file
pytest abutils/tests/test_sequence.py

# Run a specific test
pytest abutils/tests/test_sequence.py::test_function_name -v
```

Build uses `setup.py` with setuptools (no pyproject.toml). CI runs pytest on Python 3.10-3.13 across Ubuntu versions.

## Key Version Constraints

- Python 3.10+
- `numpy<2` (hard requirement)
- `polars>=1.6`
- `pyfamsa<0.7.0`
- `seaborn>=0.11`

## Code Architecture

### Namespace Pattern (scanpy-style)

The public API uses short module aliases defined as top-level re-export modules:
- `abutils.tl` (`tl.py`) - Tools: re-exports from `tools/*` and select `core` functions (e.g., `assign_pairs`, `group_lineages`, `codon_optimize`)
- `abutils.pl` (`pl.py`) - Plots: re-exports from `plots/*`
- `abutils.io` (`io.py`) - Dual role: contains substantial I/O implementation code AND re-exports from `core` modules
- `abutils.pp` (`pp.py`) - Preprocessing: re-exports from `tools/preprocessing.py`
- `abutils.cl` / `abutils.color` (`cl.py`) - Color utilities: re-exports from `utils/color.py`

Core classes `Sequence`, `Pair`, `Lineage` are directly accessible from `abutils`.

### Core Data Models (`abutils/core/`)

- **Sequence** (`sequence.py`, ~2000 lines): Container for individual sequences. Accepts raw strings, dicts, BioPython SeqRecords, or iterables. Annotations in `._annotations` dict accessible via bracket notation (`seq['v_gene']`). Auto-generates UUIDs when no ID provided. Includes static I/O methods (`read_fasta`, `read_fastq`, `read_airr`, `read_csv`, `read_parquet`, `from_mongodb`).
- **Pair** (`pair.py`, ~1500 lines): Paired heavy/light (BCR) or alpha/beta (TCR) chains. Auto-detects receptor type from locus/V-gene annotations. Chain access via properties: `.heavy`, `.light`, `.alpha`, `.beta`.
- **Lineage** (`lineage.py`, ~1500 lines): Groups related Pairs into clonal lineages. Key methods: `phylogeny()` (tree construction), `donut()` (visualization), `cluster()`. Also defines `LineageSummary` class.

### Tools Module (`abutils/tools/`)

Implementation of core algorithms:
- `alignment.py`: MSA (FAMSA, MAFFT, MUSCLE) and pairwise alignment (local, global, semiglobal via parasail)
- `cluster.py`: Sequence clustering with CD-HIT, MMseqs2, or VSEARCH backends
- `clonify.py`: Lineage/clonotype assignment using CDR3 similarity with shared mutation bonus
- `phylo.py`: Phylogenetic tree construction (FastTree) and visualization (baltic)
- `preprocessing.py`: FASTQFile classes, deduplication, quality filtering, adapter trimming
- `similarity.py`: Repertoire comparison metrics (Morisita-Horn, Jensen-Shannon, Jaccard, etc.)

### Plots Module (`abutils/plots/`)

Uses a composable pattern with base classes + mixins:
- `base.py` / `models/data.py` - Base plotting infrastructure and `PlotData` class for flexible input processing
- `mixins/` - Modular behaviors (e.g., `hue.py` for color handling)
- Individual plot types: `bar.py`, `scatter.py`, `heatmap.py`, `kde.py`, `ridge.py`, `donut.py`, `lineage.py`

### Binary Management (`bin.py` + `abutils/binaries/`)

External tools (CD-HIT, FastTree, MAFFT, MUSCLE, MMseqs2, VSEARCH, minimap2, fastp) are bundled as platform-specific binaries. `bin.get_path("tool_name")` resolves the correct binary for the current platform (darwin/linux, amd64/arm64) and downloads from S3 on first use if not included in the package.

### Utils Module (`abutils/utils/`)

Support functions organized by domain: `alignment.py`, `cluster.py`, `phylogeny.py`, `color.py`, `germlines.py`, `codons.py`, `seqio.py`, `database.py`, `path.py`, etc. Also contains `decorators.py` which provides `@lazy_property` used throughout core classes for caching expensive computed values.

## Key Conventions

- AIRR format: tab-separated values with standard fields (`sequence_id`, `sequence`, `v_call`, `j_call`, `cdr3`, etc.)
- DataFrame column suffixes for paired data: `:0` (heavy/alpha chain) and `:1` (light/beta chain)
- Polars preferred internally over Pandas; lazy evaluation (`pl.LazyFrame`) used where possible
- Test data lives in `abutils/tests/data/` (e.g., `PG9.json`, `PG9.fasta`)
- Tests are in `abutils/tests/` following `test_*.py` naming
