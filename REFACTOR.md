# abutils Refactoring Plan

Comprehensive analysis of the abutils codebase for reorganization, modernization, and cleanup.
**Generated:** 2026-03-27 | **Current version:** 0.5.4 | **Branch:** refactor

---

## Table of Contents

1. [Project Structure & Packaging](#1-project-structure--packaging)
2. [Dead Code & Commented-Out Blocks](#2-dead-code--commented-out-blocks)
3. [Code Duplication](#3-code-duplication)
4. [Test Suite Assessment](#4-test-suite-assessment)
5. [Code Organization Issues](#5-code-organization-issues)
6. [Outdated Python Syntax](#6-outdated-python-syntax)
7. [Implementation Plan](#7-implementation-plan)

---

## 1. Project Structure & Packaging

### 1.1 Current State

The project uses legacy `setup.py` + `setup.cfg` packaging:

```
abutils/
├── setup.py              # setuptools-based, some commented-out code
├── setup.cfg             # minimal (just description-file)
├── MANIFEST.in           # explicit binary/data includes
├── requirements.txt      # 31 dependencies
├── abutils/              # package code (no src/ layout)
│   ├── __init__.py
│   ├── version.py        # __version__ = "0.5.4"
│   ├── bin.py, cl.py, io.py, log.py, pl.py, pp.py, tl.py
│   ├── core/             # Sequence, Pair, Lineage
│   ├── tools/            # analysis algorithms
│   ├── plots/            # visualization
│   ├── utils/            # support utilities
│   ├── preprocess/       # entirely commented out
│   ├── binaries/         # platform-specific executables
│   └── tests/            # pytest tests
└── .github/workflows/
    ├── pytest.yml         # CI: Python 3.10-3.13, Ubuntu
    └── pythonpublish.yml  # PyPI: setuptools sdist+wheel
```

### 1.2 Issues

| Issue | Severity | Description |
|-------|----------|-------------|
| No `pyproject.toml` | HIGH | Still using `setup.py` -- PEP 517/518/621 non-compliant |
| No `src/` layout | HIGH | Package lives at repo root, not under `src/` |
| `version.py` pattern | MEDIUM | Should use dynamic versioning (hatchling, setuptools-scm) |
| `MANIFEST.in` + `include_package_data` | MEDIUM | Fragile binary inclusion; hatchling handles this better |
| `pythonpublish.yml` uses `setup.py sdist bdist_wheel` | MEDIUM | Should use `python -m build` with hatchling |
| `requirements.txt` as dependency source | LOW | Dependencies should live in `pyproject.toml` |
| `setup.cfg` is vestigial | LOW | Only contains `description-file = README.md` |

### 1.3 Plan: Migrate to Modern Packaging

**Target structure:**
```
abutils/
├── pyproject.toml           # single source of truth
├── src/
│   └── abutils/
│       ├── __init__.py
│       ├── _version.py      # or use hatch-vcs
│       ├── bin.py, cl.py, io.py, ...
│       ├── core/
│       ├── tools/
│       ├── plots/
│       ├── utils/
│       ├── binaries/
│       └── py.typed         # PEP 561 marker
├── tests/                   # outside src/
│   ├── conftest.py
│   ├── data/                # test data (moved from abutils/tests/data/)
│   ├── test_alignment.py
│   └── ...
├── .github/workflows/
│   ├── pytest.yml
│   └── pythonpublish.yml    # use hatchling
└── docs/
```

**Steps:**

1. Create `pyproject.toml` with hatchling backend:
   - Move all metadata from `setup.py`
   - Move dependencies from `requirements.txt`
   - Configure binary inclusion via `[tool.hatch.build]`
   - Add `[tool.ruff]`, `[tool.pytest.ini_options]` sections
   - Pin `python_requires = ">=3.11"`

2. Create `src/abutils/` layout:
   - `mkdir -p src && mv abutils src/`
   - Update all internal imports (they should remain relative, so mostly unchanged)
   - Update `BINARY_DIR` in `__init__.py`

3. Move tests outside the package:
   - `mv src/abutils/tests tests/`
   - Create `tests/conftest.py` with shared fixtures
   - Move `tests/data/` alongside tests
   - Update pytest config in `pyproject.toml`

4. Update CI workflows:
   - `pythonpublish.yml`: replace `setup.py sdist bdist_wheel` with `python -m build`
   - Install hatchling as build dependency
   - Keep trusted publisher OIDC (already correct)

5. Remove legacy files:
   - Delete `setup.py`, `setup.cfg`, `MANIFEST.in`
   - Delete `requirements.txt` (deps now in pyproject.toml; optionally keep for dev convenience)
   - Delete `abutils/version.py` (use `importlib.metadata` or hatch-vcs)

---

## 2. Dead Code & Commented-Out Blocks

### 2.1 Entirely Commented-Out Files

These files contain nothing but license headers and commented-out code. All functionality has been migrated to `tools/`.

| File | Lines | Status |
|------|-------|--------|
| `utils/alignment.py` | 2,293 | 100% commented out -- migrated to `tools/alignment.py` |
| `utils/cluster.py` | 1,079 | 100% commented out -- migrated to `tools/cluster.py` |
| `preprocess/preprocess.py` | 590 | 100% commented out -- `deduplicate()`, `reduction()` and helpers |

**Action:** Delete all three files. Remove `preprocess/` package entirely. Content preserved in git history.

**Lines recovered: ~3,962**

### 2.2 Large Commented-Out Blocks in Active Files

| File | Lines | Description |
|------|-------|-------------|
| `utils/path.py` | ~139 lines (161-299) | Legacy `initialize()`, alternative `list_files()`/`make_dir()`, `print_splash()` |
| `setup.py` | ~55 lines (1-56) | Old config dict and imports |
| `tl.py` | ~10 lines (37-46) | Old `utils.phylogeny` and `utils.alignment` imports |
| `__init__.py` | several blocks | Commented-out import groups from older API |
| `io.py` | several blocks | Commented-out imports from `core.sequence` |
| `pp.py` | several lines | Old `preprocess.preprocess` imports |

**Action:** Remove all commented-out code blocks. They add noise and are preserved in git history.

### 2.3 TODO/FIXME Comments

| File | Line | Comment |
|------|------|---------|
| `tests/test_alignment.py` | 273 | `# TODO: fix famsa alignment file and alignment string tests` |
| `plots/scatter.py` | 451 | `# TODO: allow different markers that correspond to marker "categories"` |
| `plots/kde.py` | 430 | `# TODO: allow different markers that correspond to marker "categories"` |
| `plots/kde.py` | 1332 | `# TODO` (bare, no description) |
| `tools/phylo.py` | 593 | `# TODO: add support for branch color orders` |

**Action:** Resolve or remove. The bare `# TODO` in `kde.py` should definitely be removed. The scatter/kde marker TODOs are duplicates suggesting a shared concern that should be tracked in an issue.

### 2.4 Empty/Vestigial Files

Several `utils/` files appear to be empty stubs (just license headers, no active code):

| File | Lines | Notes |
|------|-------|-------|
| `utils/circos.py` | 33 | Minimal -- verify if used |
| `utils/codons.py` | 91 | Codon table lookup -- verify if used |
| `utils/seqio.py` | 332 | Verify overlap with `io.py` |
| `preprocess/__init__.py` | ~0 | Package init for dead module |

**Action:** Verify whether each is imported anywhere. Delete if unused; consolidate if functionality overlaps with another module.

---

## 3. Code Duplication

### 3.1 Critical: Duplicate Function Implementations

The following functions exist as near-identical copies in two locations:

#### File Operations (io.py vs utils/path.py)

| Function | `io.py` | `utils/path.py` |
|----------|---------|-----------------|
| `make_dir()` | lines 80-91 | lines 38-49 |
| `list_files()` | lines 94-165 | lines 52-123 |
| `rename_file()` | lines 168-181 | lines 126-139 |
| `delete_files()` | lines 184-200 | lines 142-158 |

**Action:** Keep implementations in `io.py` (the public API surface). Delete `utils/path.py` entirely and update any internal imports.

#### Logging Infrastructure (tools/log.py vs utils/log.py)

| Function | `tools/log.py` | `utils/log.py` |
|----------|----------------|----------------|
| `setup_logging()` | lines 254-309 | lines 75-130 |
| `null_logger()` | lines 312-316 | lines 133-137 |
| `get_logger()` | lines 319-351 | lines 140-172 |
| `_add_stream_handler()` | lines 354-380 | lines 175-201 |

Additionally, `tools/log.py` has extra classes (`LoggingMixin`, `SimpleLogger`, `NotebookLogger`) not in `utils/log.py`.

**Action:** Delete `utils/log.py`. The top-level `log.py` already re-exports from `tools/log.py`.

#### Sequence I/O (io.py vs core/sequence.py)

| Function | `io.py` | `core/sequence.py` |
|----------|---------|---------------------|
| `read_csv()` | lines 977-1046 | lines 898-967 |
| `read_airr()` | lines 1147-1207 | lines 969-1000 |
| `read_parquet()` | lines 1289-1352 | lines 1002-1065 |
| `to_airr()` | lines 1209-1287 | lines 1948-2026 |

**Action:** Keep canonical implementations in `core/sequence.py` as static methods. Have `io.py` delegate to them rather than re-implementing.

#### Lineage Plotting Helpers

| Function | `core/lineage.py` | `plots/lineage.py` |
|----------|-------------------|---------------------|
| `_get_donut_colors()` | lines 1370-1373 | lines 197-200 |
| `_get_monochrome_colors()` | lines 1376-1388 | lines 203-211 |

**Action:** Keep in `plots/lineage.py` only. Extract plotting methods from `core/lineage.py` entirely (see Section 5.2).

### 3.2 Duplicate Data Processing Classes

Two parallel plot data processing implementations:

- `plots/data.py` (310 lines): `InputData` class + `process_input_data()` function
- `plots/models/data.py` (815 lines): `PlotData` class -- comprehensive, used by `pl.py`

**Action:** Determine which plot modules use which. If `plots/data.py` is unused or legacy, delete it. If both are needed, clarify their roles and rename to avoid confusion.

---

## 4. Test Suite Assessment

### 4.1 Coverage Overview

**14 test files, ~275 tests total.** Test data consists of only 2 real files (`PG9.fasta`, `PG9.json`); everything else is synthetic.

#### Modules With NO Tests

| Module | Lines | Public Functions/Classes |
|--------|-------|------------------------|
| `bin.py` | ~100 | `get_path()`, `copy_to_binary_directory()` |
| `utils/database.py` | 395 | `SQLiteDatabase`, `KeyValueStore` |
| `utils/decorators.py` | 121 | `@lazy_property` |
| `utils/germlines.py` | 544 | `germline_names()`, `get_germline()`, `germlines()` |
| `utils/path.py` | 298 | `make_dir()`, `list_files()`, etc. (duplicated in io.py) |
| `utils/jobs.py` | 173 | Job monitoring functions |
| `utils/convert.py` | 74 | `abi_to_fasta()` |
| `utils/s3.py` | 310 | S3 upload/download |
| `utils/ssh_tunnel.py` | 232 | SSH tunneling |
| `utils/mongodb.py` | 492 | MongoDB query functions |
| `utils/seqio.py` | 332 | Sequence I/O utilities |
| `tools/anarci.py` | 93 | `number()`, `assign_numbering()` |
| `tools/search.py` | 229 | `search()` MMseqs2 wrapper |
| `plots/ring.py` | 145 | Ring plot |
| `plots/summary.py` | ~100 | Summary plot |
| `plots/lineage.py` | 210 | Lineage plot |
| `plots/models/data.py` | 815 | `PlotData` class |

#### Partially Tested Modules (significant gaps)

| Module | Tested | Untested |
|--------|--------|----------|
| `core/sequence.py` | `__init__`, `__contains__`, `__getitem__`, `reverse_complement`, `translate`, `codon_optimize`, `region`, `as_fasta`, `read_fasta`, `read_fastq`, `read_airr`, `read_csv`, `read_parquet`, `from_mongodb` | `kmers()`, `mutation_count()`, `junction()`, `write_fasta()`, `write_fastq()`, `write_json()`, `fasta` property, `sequence_aa`, `qualscores`, `header`, `name` |
| `core/pair.py` | `__init__`, chain properties, `is_pair`, `name`, `select_chain` | Plural chain properties (`heavies`, `lights`, `alphas`, etc.), chain assignment edge cases |
| `core/lineage.py` | Basic construction, pair access | `phylogeny()`, `cluster()`, `donut()`, UCA assignment |
| `tools/alignment.py` | MSA functions (FAMSA, MAFFT, MUSCLE), pairwise alignment | `GlobalAlignment`, `LocalAlignment`, `SemiGlobalAlignment` classes directly, `consensus()` |
| `tools/cluster.py` | `cluster()`, backend-specific functions | `Cluster.consensus`, `Clusters` iteration, centroid access |
| `tools/preprocessing.py` | `FASTQFile`, `IlluminaFile`, `ElementFile` init | Actual merge/trim execution (all mocked), `deduplicate()`, `reduction()` |
| `io.py` | File read/write basics | `concatenate_*()`, `split_*()`, `from_pandas()`, `from_polars()` |

### 4.2 Test Quality Issues

#### A. Over-Reliance on Synthetic Data

Almost every test constructs trivial synthetic sequences (e.g., `"AAAAA"`, `"CCCCC"`, repeated p53 sequences). This misses:
- Real antibody sequences with complex V(D)J annotations
- Sequences with ambiguous bases, gaps, or quality issues
- Realistic paired heavy/light data with proper annotation fields
- Diverse CDR3 lengths and germline families

**Recommendation:** Create a `tests/data/` directory (at project root after migration) with:
- A small set of real paired BCR sequences (10-20 pairs) with full AIRR annotations
- A small AIRR TSV file with realistic annotation fields
- A representative FASTA/FASTQ with real antibody sequences
- At least one TCR dataset for alpha/beta pair testing

#### B. Shallow Assertions

Many tests only verify "doesn't crash" or check output types:

```python
# Typical weak assertion pattern found throughout:
ax = bar(x="category", data=categorical_data)
assert isinstance(ax, matplotlib.axes.Axes)  # doesn't verify content

result = cluster(sequences, threshold=0.9)
assert isinstance(result, Clusters)  # doesn't verify cluster contents
```

Tests should verify actual computed values where possible (like `test_similarity.py` does well with `pytest.approx`).

#### C. Mocked External Tools in Preprocessing Tests

`test_preprocessing.py` mocks all subprocess calls to fastp, FLASH2, sickle, etc. This means the tests never verify that the actual tool invocations produce correct results. These should be integration tests with real (small) FASTQ files.

#### D. Missing Edge Case Coverage

- Empty input handling (empty lists, empty strings, None values)
- Sequences with non-standard characters (ambiguous bases, gaps)
- Unicode in sequence IDs
- Very long sequences (stress testing)
- Malformed FASTA/FASTQ files
- Pairs with missing chains or conflicting annotations

### 4.3 Unnecessary/Redundant Tests

| Test | File | Issue |
|------|------|-------|
| `test_nested_dict_to_dataframe` | `test_utilities.py` | Only test in file; tests a trivial utility function |
| Smoke tests in `test_plots.py` | `test_plots.py` | 11 tests that only check `isinstance(ax, Axes)` -- low value |

### 4.4 Plan: Test Improvements

**Phase 1 -- Real data and fixtures:**
- Curate a small real-world test dataset with diverse antibody sequences
- Create shared `conftest.py` fixtures for common test data
- Replace synthetic sequences in existing tests where feasible

**Phase 2 -- Coverage gaps:**
- Add tests for all completely untested modules (prioritize `utils/decorators.py`, `bin.py`, `utils/germlines.py`)
- Add tests for untested methods on `Sequence`, `Pair`, `Lineage`
- Add `PlotData` unit tests

**Phase 3 -- Quality improvements:**
- Strengthen assertions in existing tests (verify values, not just types)
- Add edge case tests for core data models
- Replace mocked subprocess calls in preprocessing with real integration tests (mark `@pytest.mark.slow`)

---

## 5. Code Organization Issues

### 5.1 Plotting Logic in Core Data Models

`core/lineage.py` contains ~200 lines of plotting code that belongs in `plots/`:

- `donut()` method (~70 lines) -- creates donut chart visualization
- `_get_donut_colors()` helper
- `_get_monochrome_colors()` helper
- Various matplotlib imports used only by plotting methods

**Action:** Move all plotting methods to `plots/lineage.py`. The `Lineage` class should be a pure data model. Plotting functions can accept a `Lineage` object as input.

### 5.2 Overlapping Module Responsibilities

The `utils/` vs `tools/` split is unclear. The intended pattern appears to be:
- `tools/` = user-facing analysis functions (the "active" layer)
- `utils/` = internal support functions

But several `utils/` modules overlap with `tools/`:

| `utils/` module | `tools/` counterpart | Status |
|-----------------|---------------------|--------|
| `utils/alignment.py` | `tools/alignment.py` | utils is commented out -- **DELETE** |
| `utils/cluster.py` | `tools/cluster.py` | utils is commented out -- **DELETE** |
| `utils/phylogeny.py` | `tools/phylo.py` | Both active -- unclear boundary |
| `utils/log.py` | `tools/log.py` | Duplicate -- **DELETE utils version** |
| `utils/seqio.py` | `io.py` | Overlap -- needs consolidation |

**Action for `utils/phylogeny.py` vs `tools/phylo.py`:** These serve different purposes -- `tools/phylo.py` wraps FastTree and provides the `Phylogeny` class, while `utils/phylogeny.py` contains lower-level tree manipulation utilities. If `utils/phylogeny.py` is only used by `tools/phylo.py`, make it a private submodule of tools or merge.

### 5.3 `io.py` Is Too Large and Has Mixed Responsibilities

At 1,519 lines, `io.py` serves three distinct roles:
1. **File system operations** (`make_dir`, `list_files`, `delete_files`, `rename_file`)
2. **Sequence format I/O** (`read_fasta`, `read_csv`, `read_airr`, `to_parquet`, etc.)
3. **Format conversion and splitting** (`split_fasta`, `concatenate_airr`, etc.)

**Action:** If `io.py` is converted to a subpackage (`io/`), split into:
- `io/_files.py` -- file system operations
- `io/_sequences.py` -- sequence read/write
- `io/_formats.py` -- format conversion, concatenation, splitting
- `io/__init__.py` -- re-exports (preserves `abutils.io.read_fasta()` API)

### 5.4 Very Large Files

| File | Lines | Recommendation |
|------|-------|----------------|
| `tools/alignment.py` | 2,086 | Split: MSA classes/functions vs pairwise alignment classes |
| `core/sequence.py` | 2,005 | Acceptable (single class), but consider extracting static I/O methods |
| `plots/kde.py` | 1,806 | Split: core KDE logic vs marginal/bivariate variants |
| `core/lineage.py` | 1,566 | Remove plotting methods (~200 lines); rest is acceptable |
| `io.py` | 1,519 | Split into subpackage (see 5.3) |
| `core/pair.py` | 1,473 | Acceptable (single class) |
| `tools/preprocessing.py` | 1,411 | Acceptable (single domain) |
| `tools/cluster.py` | 1,220 | Acceptable (single domain) |
| `tools/clonify.py` | 1,185 | Acceptable (single domain) |

### 5.5 Naming Issues

#### Confusing Duplicate Names

- Three different `donut()` functions: `core/lineage.py`, `plots/donut.py`, and accessible via `pl.donut()`
- Two `data.py` files in plots: `plots/data.py` and `plots/models/data.py`
- Two `lineage.py` files: `core/lineage.py` and `plots/lineage.py` (this is fine -- different subpackages)
- Two `alignment.py` files: `tools/alignment.py` and `utils/alignment.py` (utils one should be deleted)
- Two `cluster.py` files: `tools/cluster.py` and `utils/cluster.py` (utils one should be deleted)
- Two `log.py` files: `tools/log.py` and `utils/log.py` (utils one should be deleted)

#### Non-Descriptive Names

- `utils/utilities.py` -- generic name, only contains `nested_dict_to_dataframe()`
- `utils/jobs.py` -- unclear what "jobs" means without context
- `utils/click.py` -- could be confused with the Click CLI library itself

---

## 6. Outdated Python Syntax

With Python 3.11+ as the minimum, the following patterns should be modernized:

### 6.1 Type Annotation Modernization

| Pattern | Count | Files | Replacement |
|---------|-------|-------|-------------|
| `Optional[X]` | **520** | 31 files | `X \| None` |
| `Union[X, Y]` | **392** | 24 files | `X \| Y` |
| `from typing import Optional, Union, ...` | **28** active | 28 files | Remove; keep only `Iterable`, `Callable` if needed |
| `List[X]` | 1 | `tools/phylo.py` | `list[X]` |
| `Tuple[X, Y]` | 1 | `plots/heatmap.py` | `tuple[X, Y]` |

**Top files by occurrence count:**

| File | `Optional[]` | `Union[]` | Total |
|------|-------------|-----------|-------|
| `io.py` | 54 | 30 | 84 |
| `core/sequence.py` | 52 | 36 | 88 |
| `tools/alignment.py` | 52 | 26 | 78 |
| `plots/models/data.py` | 47 | 52 | 99 |
| `core/pair.py` | 44 | 14 | 58 |
| `plots/heatmap.py` | 28 | 15 | 43 |
| `plots/scatter.py` | 24 | 23 | 47 |
| `plots/kde.py` | 24 | 20 | 44 |
| `plots/ring.py` | 24 | 4 | 28 |
| `tools/preprocessing.py` | 17 | 3 | 20 |
| `tools/phylo.py` | 13 | 35 | 48 |

**Note:** After converting all `Optional[X]` to `X | None` and `Union[X, Y]` to `X | Y`, most `from typing import` statements can be reduced to just `Iterable` and `Callable` (which are still needed from `collections.abc` or `typing`). With 3.11+, `Iterable` and `Callable` can be imported from `collections.abc` instead.

### 6.2 Class Definitions

| Pattern | Count | Location | Fix |
|---------|-------|----------|-----|
| `class Sequence(object):` | 1 | `core/sequence.py:52` | `class Sequence:` |
| `class Pair(object):` | 1 | `core/pair.py:37` | `class Pair:` |

### 6.3 Super Calls

| Pattern | Location | Fix |
|---------|----------|-----|
| `super(SQLiteDatabase, self).__init__()` | `utils/database.py:49` | `super().__init__()` |
| `super(KeyValueStore, self).__init__(...)` | `utils/database.py:223` | `super().__init__(...)` |

### 6.4 `from __future__` Imports (all no-ops in 3.11+)

| File | Import |
|------|--------|
| `utils/convert.py:26` | `from __future__ import absolute_import, division, print_function, unicode_literals` |
| `utils/decorators.py:26` | same |
| `utils/germlines.py:16` | same |
| `utils/progbar.py:26` | same |
| `plots/base.py:26` | `from __future__ import print_function` |
| `plots/summary.py:26` | `from __future__ import print_function` |
| `plots/__init__.py:1` | `from __future__ import absolute_import` |

**Action:** Remove all. These are Python 2 compatibility artifacts.

### 6.5 String Formatting

**171 occurrences** of `.format()` across 24 files, plus a few `%`-style format strings (e.g., in `utils/color.py`).

**Top files:**
- `utils/phylogeny.py` (20), `utils/database.py` (16), `tools/alignment.py` (14), `utils/s3.py` (11), `utils/cluster.py` (11), `utils/mongodb.py` (10)

**Action:** Convert to f-strings during any file touched for other reasons. Not worth a dedicated pass for the non-commented files since `.format()` is still valid Python.

### 6.6 `os.path` Usage

**292 occurrences** across 28 files.

**Heaviest users:** `tools/preprocessing.py` (33), `utils/cluster.py` (23, but commented out), `utils/path.py` (14), `io.py` (21), `bin.py` (19)

**Action:** Migrate to `pathlib.Path` incrementally. Prioritize `bin.py`, `io.py`, and `tools/preprocessing.py`. Accept `os.path` in commented-out files that will be deleted.

### 6.7 Unnecessary `.keys()` Calls

Patterns like `if item in self.pair_dict.keys():` found in:
- `core/lineage.py:152`
- `core/sequence.py:196,221`
- Various test files

**Action:** Replace `x in dict.keys()` with `x in dict`.

---

## 7. Implementation Plan

### Phase 1: Dead Code Removal (Low Risk, High Impact)

**Estimated scope:** Delete ~4,500 lines, zero functional changes.

| Step | Action | Lines Removed |
|------|--------|---------------|
| 1.1 | Delete `utils/alignment.py` (entirely commented out) | 2,293 |
| 1.2 | Delete `utils/cluster.py` (entirely commented out) | 1,079 |
| 1.3 | Delete `preprocess/preprocess.py` and `preprocess/__init__.py` | ~600 |
| 1.4 | Delete `utils/log.py` (duplicate of `tools/log.py`) | 201 |
| 1.5 | Remove commented-out blocks in `utils/path.py`, `tl.py`, `pp.py`, `__init__.py`, `io.py`, `setup.py` | ~250 |
| 1.6 | Remove all `from __future__ import` statements | 7 lines across 7 files |
| 1.7 | Remove bare/orphan TODO comments | 5 |
| 1.8 | Verify no imports reference deleted modules; fix any that do | -- |

**Validation:** Run `pytest` after each deletion to catch broken imports.

### Phase 2: Deduplication (Medium Risk)

**Estimated scope:** Consolidate duplicated functions, net reduction ~800 lines.

| Step | Action |
|------|--------|
| 2.1 | Delete `utils/path.py`; ensure `io.py` has the canonical file-ops implementations; update internal imports |
| 2.2 | Make `io.py` I/O functions delegate to `core/sequence.py` static methods instead of reimplementing (`read_csv`, `read_airr`, `read_parquet`, `to_airr`) |
| 2.3 | Remove `donut()`, `_get_donut_colors()`, `_get_monochrome_colors()` from `core/lineage.py`; ensure `plots/lineage.py` and `plots/donut.py` are the sole owners |
| 2.4 | Audit `plots/data.py` vs `plots/models/data.py` -- delete the unused one |

**Validation:** Run full test suite. Check that `abutils.io.read_csv()`, `abutils.tl.*`, `abutils.pl.*` all still work.

### Phase 3: Project Structure Migration (High Risk, High Impact)

| Step | Action |
|------|--------|
| 3.1 | Create `pyproject.toml` with hatchling backend, migrate all metadata and dependencies |
| 3.2 | Create `src/` layout: move `abutils/` under `src/abutils/` |
| 3.3 | Move `tests/` to project root (out of package) |
| 3.4 | Update `BINARY_DIR` and any hardcoded paths |
| 3.5 | Update CI workflows (`pytest.yml`, `pythonpublish.yml`) |
| 3.6 | Delete `setup.py`, `setup.cfg`, `MANIFEST.in` |
| 3.7 | Test editable install: `pip install -e .` |
| 3.8 | Test wheel build: `python -m build` |
| 3.9 | Verify binary inclusion in built wheel |

**Validation:** Full CI matrix (3.10-3.13, Ubuntu). Verify `pip install` from built wheel works.

### Phase 4: Syntax Modernization (Low Risk)

| Step | Action | Scope |
|------|--------|-------|
| 4.1 | Replace all `Optional[X]` with `X \| None` | 520 occurrences, 31 files |
| 4.2 | Replace all `Union[X, Y]` with `X \| Y` | 392 occurrences, 24 files |
| 4.3 | Replace `List[X]`, `Tuple[X, Y]` with `list[X]`, `tuple[X, Y]` | 2 occurrences |
| 4.4 | Update `from typing import` to import only `Iterable`, `Callable` from `collections.abc` | 28 files |
| 4.5 | Fix `class X(object):` -> `class X:` | 2 occurrences |
| 4.6 | Fix `super(ClassName, self).__init__()` -> `super().__init__()` | 2 occurrences |
| 4.7 | Remove `.keys()` from `in dict.keys()` | ~10 active occurrences |
| 4.8 | Convert `.format()` to f-strings in files touched for other reasons | opportunistic |
| 4.9 | Migrate `os.path` -> `pathlib.Path` in core files | opportunistic |

**Validation:** `ruff check`, then `pytest`.

### Phase 5: Code Organization Improvements (Medium Risk)

| Step | Action |
|------|--------|
| 5.1 | Extract plotting methods from `core/lineage.py` to `plots/lineage.py` |
| 5.2 | Convert `io.py` to `io/` subpackage if warranted by size after dedup |
| 5.3 | Clarify `utils/phylogeny.py` vs `tools/phylo.py` boundary |
| 5.4 | Audit remaining `utils/` modules for usage; delete truly unused ones |
| 5.5 | Consider splitting `tools/alignment.py` (2,086 lines) into MSA + pairwise modules |

### Phase 6: Test Suite Overhaul (ongoing)

| Step | Action |
|------|--------|
| 6.1 | Create shared `tests/conftest.py` with real-data fixtures |
| 6.2 | Curate real-world test datasets (BCR pairs, TCR pairs, AIRR TSV) |
| 6.3 | Add tests for untested modules (prioritize: `decorators.py`, `bin.py`, `germlines.py`, `PlotData`) |
| 6.4 | Strengthen assertions in existing tests (verify values, not just types) |
| 6.5 | Add edge case tests: empty inputs, malformed data, ambiguous bases |
| 6.6 | Replace mocked subprocess calls in preprocessing tests with integration tests |
| 6.7 | Mark integration tests with `@pytest.mark.slow` |

---

## Appendix: File Inventory

### By Size (active code only, excluding commented-out files)

| File | Lines | Category |
|------|-------|----------|
| `tools/alignment.py` | 2,086 | tools |
| `core/sequence.py` | 2,005 | core |
| `plots/kde.py` | 1,806 | plots |
| `core/lineage.py` | 1,566 | core |
| `io.py` | 1,519 | top-level |
| `core/pair.py` | 1,473 | core |
| `tools/preprocessing.py` | 1,411 | tools |
| `tools/cluster.py` | 1,220 | tools |
| `tools/clonify.py` | 1,185 | tools |
| `tools/phylo.py` | 1,066 | tools |
| `utils/phylogeny.py` | 994 | utils |
| `tools/similarity.py` | 834 | tools |
| `plots/models/data.py` | 815 | plots |
| `plots/scatter.py` | 738 | plots |
| `utils/germlines.py` | 544 | utils |
| `plots/heatmap.py` | 509 | plots |
| `plots/donut.py` | 501 | plots |
| `utils/mongodb.py` | 492 | utils |
| `utils/color.py` | 410 | utils |
| `utils/database.py` | 395 | utils |
| `tools/log.py` | 380 | tools |
| `plots/bar.py` | 377 | plots |
| `utils/seqio.py` | 332 | utils |
| `plots/data.py` | 310 | plots |
| `utils/s3.py` | 310 | utils |
| `utils/path.py` | 298 | utils (duplicate) |
| `plots/ridge.py` | 290 | plots |
| `utils/ssh_tunnel.py` | 232 | utils |
| `tools/search.py` | 229 | tools |
| `plots/lineage.py` | 210 | plots |
| `tools/cloning.py` | 204 | tools |
| `utils/jobs.py` | 173 | utils |
| `plots/ring.py` | 145 | plots |
| `utils/decorators.py` | 121 | utils |
| `utils/utilities.py` | 110 | utils |
| `utils/progbar.py` | 102 | utils |
| `plots/base.py` | 100 | plots |
| `tools/anarci.py` | 93 | tools |
| `utils/codons.py` | 91 | utils |
| `utils/convert.py` | 74 | utils |
| `utils/click.py` | 56 | utils |
| `plots/utils.py` | 54 | plots |
| `tl.py` | 47 | top-level |
| `utils/circos.py` | 33 | utils |
| `plots/mixins/hue.py` | 32 | plots |
| `pl.py` | 31 | top-level |
| `pp.py` | 28 | top-level |
| `cl.py` | 28 | top-level |
| **Active total** | **~27,000** | |
| **Dead code** | **~4,500** | utils/alignment, utils/cluster, preprocess, commented blocks |
| **Grand total** | **~35,250** | |

### Files to Delete (Phase 1-2)

| File | Lines | Reason |
|------|-------|--------|
| `utils/alignment.py` | 2,293 | 100% commented out |
| `utils/cluster.py` | 1,079 | 100% commented out |
| `preprocess/preprocess.py` | 590 | 100% commented out |
| `preprocess/__init__.py` | ~5 | Dead package |
| `utils/log.py` | 201 | Exact duplicate of tools/log.py |
| `utils/path.py` | 298 | Duplicate of io.py functions |
| `plots/data.py` | 310 | Likely superseded by plots/models/data.py (verify first) |
| **Total** | **~4,776** | |
