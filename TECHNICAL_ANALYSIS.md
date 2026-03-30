# Technical Analysis of `abutils`

Note: the request said "oplm", but the repository analyzed here is `/home/briney/git/abutils`, and the requested concerns
(AIRR utilities, packaged binaries, plotting/tools API) clearly match `abutils`.

## Scope and method

- Read the package structure, docs, CI workflow, public namespace shims, and representative core/tool/plot modules.
- Ran the full test suite locally on March 30, 2026:
  - `349 passed, 1 skipped, 1 xfailed, 1 xpassed in 2.87s`
- Performed a few runtime API smoke checks for documented import paths and edge-case behaviors.
- Compared the current API shape against the scanpy public-API model:
  - Scanpy API docs: https://scanpy.readthedocs.io/en/1.9.x/api.html
  - Scanpy repo/README public API statement: https://github.com/scverse/scanpy

## Executive summary

The codebase is in materially better shape than many research libraries: the suite is fast, green, CPU-friendly, and it does execute several real packaged binaries. The core models also now have meaningful real-data fixtures, which is a real improvement over purely synthetic smoke tests.

The main problem is not basic test count. The main problem is that the package's documented public API and its actual import surface have drifted apart badly, and the current tests do not assert that contract. That is why obvious breakages still exist:

- `abutils.pl` is documented as the plotting namespace but is empty (`src/abutils/pl.py`).
- `abutils.tl` is documented as the main tools namespace, but important functions the docs advertise there are missing, including `clonify` and preprocessing helpers.
- `abutils.tools` and `abutils.plots` are used in docs/examples as callable namespaces, but both package `__init__` files are empty, so those examples are also misleading.
- There are several verified behavioral bugs in important code paths that currently have no test coverage.

If I had to prioritize, I would do this first:

1. Fix the public namespace contract and add import-surface tests.
2. Add a small set of genuine end-to-end workflow tests built on the real AIRR fixtures already in the repo.
3. Replace string-built `shell=True` subprocess wrappers with structured argv-based wrappers.
4. Split monolithic modules and pull optional dependencies into extras.

## 1. Test harness: current state

### What is strong today

- The suite is genuinely CI-friendly.
  - The entire suite completes in under 3 seconds on a CPU-only runner.
  - GitHub Actions runs across Python 3.10-3.13 on two Linux images (`.github/workflows/pytest.yml:1-31`).
- There is real binary execution coverage for some important paths:
  - MSA wrappers call packaged binaries in tests for MAFFT, MUSCLE 5, MUSCLE 3, and FAMSA (`tests/test_alignment.py:66-176`, `tests/test_alignment.py:185-266`, `tests/test_alignment.py:274-343`, `tests/test_alignment.py:352-442`).
  - Clustering wrappers call MMseqs2, CD-HIT, and VSEARCH (`tests/test_cluster.py:120-161`).
  - FastTree and `phylogeny()` are exercised with real inputs (`tests/test_phylo.py:87-160`, `tests/test_phylo.py:169-272`).
  - `merge_fastqs_fastp()` and `merge_fastqs()` run real code paths (`tests/test_preprocessing.py:273-304`).
- Real-world AIRR fixtures exist and are used:
  - PG9 JSON/FASTA and BCR/TCR AIRR TSV fixtures are centralized in `tests/conftest.py:20-174`.
  - Core model tests use those fixtures meaningfully (`tests/test_edge_cases.py:134-174`, `tests/test_edge_cases.py:291-324`, `tests/test_edge_cases.py:379-396`).

### Where the current harness is still weak

#### 1. API-contract testing is almost nonexistent

The docs promise one API; the tests validate a different one.

- Docs say plotting functions are accessible through `abutils.pl` (`docs/source/index.rst:85-102`).
- `src/abutils/pl.py` is empty (`src/abutils/pl.py:1-24`).
- Plot tests bypass the public namespace entirely and import leaf modules directly, e.g. `from abutils.plots.bar import bar` (`tests/test_plots.py:92-236`).

This same pattern exists for tools:

- Docs say preprocessing is exposed as `abutils.tl.merge_fastqs()` (`docs/source/tools/preprocessing.rst:10-21`, `docs/source/tools/preprocessing.rst:35-97`).
- Docs say clonify is exposed as `abutils.tl.clonify()` (`docs/source/tools/clonify.rst:10-20`, `docs/source/tools/clonify.rst:35-79`).
- `src/abutils/tl.py` only star-imports alignment, anarci, cloning, cluster, phylo, and search (`src/abutils/tl.py:25-30`), so `clonify` and preprocessing functions are not actually there.
- The clonify tests also bypass the documented namespace and import directly from `abutils.tools.clonify` (`tests/test_clonify.py:5-6`).

This is the single biggest testing deficiency because it allows silent breakage of the user-facing API.

#### 2. The suite has real binary coverage, but not full packaged-binary coverage

The user requirement here is stricter than "some binaries are tested". The repository packages a larger binary set than the suite truly exercises.

`src/abutils/bin.py` advertises these binaries (`src/abutils/bin.py:91-172`):

- `blastn`
- `makeblastdb`
- `cdhit`
- `mafft`
- `minimap2`
- `mmseqs`
- `muscle`
- `muscle3`
- `vsearch`
- `fastp`
- `sickle`
- `fasttree`

What is actually exercised end-to-end:

- `mafft`, `muscle`, `muscle3`, `mmseqs`, `cdhit`, `vsearch`, `fasttree`, `fastp`

What is not meaningfully exercised:

- `blastn`
- `makeblastdb`
- `minimap2`
- `sickle`
- ANARCI/HMMER-related paths
- `mmseqs_search()` itself

Also, `tests/test_bin.py` is only a path-construction test. It patches `os.path.exists` to always return `True` (`tests/test_bin.py:30-105`), so it does not validate that a packaged binary actually exists, is executable, or works on the current platform.

#### 3. Important preprocessing entry points are mostly untested

The high-level preprocessing functions are exactly the kind of workflow code that benefits from end-to-end tests, but most of that surface is not covered.

- `merge_fastqs_vsearch()` is mocked rather than run against the packaged binary (`tests/test_preprocessing.py:307-335`).
- `merge_fastqs_fastp()` only checks the returned path, not output content or generated reports (`tests/test_preprocessing.py:290-304`).
- `deduplicate()` test is commented out (`tests/test_preprocessing.py:356-377`).
- `reduction()` test is commented out (`tests/test_preprocessing.py:521`).
- `_process_chains()` is tested with mocked clustering and consensus generation, not with real AIRR data and real binaries (`tests/test_preprocessing.py:381-...`).

The result is that the lowest-level helper logic is partially covered, but the actual user workflows are not.

#### 4. Known-broken plotting paths are tolerated instead of driven to green

- `heatmap()` is still xfailed because it "requires specific data format via PlotData" (`tests/test_plots.py:170-177`).
- `show_palettes()` is skipped because the test author was unsure how to test it (`tests/test_color.py:169-172`).
- FAMSA file/string output tests are commented out due a known failure (`tests/test_alignment.py:283-299`).

This is not catastrophic, but it does mean the green suite currently overstates functional completeness.

### What a stronger, still-fast harness should look like

I would organize the tests into three layers:

#### A. Public API contract tests

Small import-surface tests that assert the documented namespace exists:

- `abutils.tl.mafft`
- `abutils.tl.cluster`
- `abutils.tl.clonify`
- `abutils.pp.merge_fastqs` if preprocessing stays in `pp`
- `abutils.pl.bar`, `abutils.pl.scatter`, `abutils.pl.kde`, `abutils.pl.donut`
- top-level model availability: `abutils.Sequence`, `abutils.Pair`, `abutils.Lineage`

These tests would be extremely fast and would have caught the current API drift immediately.

#### B. Binary smoke tests

For every packaged binary on Linux:

- resolve path with `get_path()`
- assert the file exists and is executable
- run a very cheap command such as `--help`, `-h`, or version output

This should be parametrized across the binary table and kept separate from higher-level workflow tests.

#### C. Small real-world workflow tests

These should use the existing real AIRR fixtures and tiny FASTQ fixtures:

- `read_airr` or `from_polars` -> `Sequence` -> `Pair` -> `clonify`
- `clonify` result -> `group_lineages` -> `phylogeny(cluster=True)`
- `merge_fastqs_fastp` and `merge_fastqs_vsearch` on tiny paired FASTQs with output validation
- `mmseqs_search` on a tiny query/target set with known hit expectations
- a plot smoke suite through the public namespace (`abutils.pl.*`) using the real-data objects where sensible

All of these can remain CPU-only and CI-friendly if the datasets stay tiny.

### Test-harness verdict

The current suite is good enough to support refactoring inside already-tested code paths. It is not yet good enough to guarantee the documented user experience or all packaged-binary workflows.

## 2. API structure and usability

### What scanpy gets right

The useful part of scanpy is not just the names `tl` and `pl`. It is the discipline around them.

From the scanpy API docs:

- `tl` is the analysis/tools namespace.
- `pl` is the plotting namespace.
- `pl` intentionally parallels `tl` for many workflows.

From the scanpy repo/README:

- Scanpy explicitly distinguishes public API from internal API and documents only the supported surface.

That combination matters. Users learn one mental model and the docs, examples, and imports all agree.

### What `abutils` does today

The repository appears to want a similar model:

- `abutils.tl`
- `abutils.pl`
- `abutils.pp`
- `abutils.cl`
- `abutils.io`

But the implementation is inconsistent.

#### 1. `abutils.pl` is not usable

- Docs explicitly say all plotting functions are accessible via `abutils.pl` (`docs/source/index.rst:85-102`).
- `src/abutils/pl.py` is empty (`src/abutils/pl.py:1-24`).
- A runtime check confirms `abutils.pl.bar` raises `AttributeError`.

That is a hard API break, not just a style concern.

#### 2. `abutils.tl` is incomplete relative to the docs

- Docs place preprocessing under `abutils.tl` (`docs/source/tools/preprocessing.rst:10-21`).
- Docs place clonify under `abutils.tl` (`docs/source/tools/clonify.rst:10-20`).
- `src/abutils/tl.py` does not import preprocessing or clonify (`src/abutils/tl.py:25-30`).
- Runtime checks confirm `abutils.tl.clonify`, `abutils.tl.pairwise_distance`, `abutils.tl.merge_fastqs`, and `abutils.tl.reduction` are absent.

This is a more serious usability problem than naming taste, because users cannot trust the docs.

#### 3. `abutils.pp` exists but is only half-committed

- `src/abutils/pp.py` does re-export preprocessing (`src/abutils/pp.py:25`).
- But `abutils.__init__` does not import `pp` (`src/abutils/__init__.py:4-9`).
- A plain `import abutils` therefore does not expose `abutils.pp` until `abutils.pp` is explicitly imported later.

This is fixable, but it reinforces the sense that the public API was never fully curated.

#### 4. `abutils.tools` and `abutils.plots` are not public namespaces despite docs/examples using them

- `src/abutils/tools/__init__.py` is empty.
- `src/abutils/plots/__init__.py` is empty.
- Runtime checks confirm attributes like `abutils.tools.mmseqs_search`, `abutils.tools.clonify`, `abutils.tools.merge_fastqs`, `abutils.plots.bar`, and `abutils.plots.scatter` are absent.

That means docs examples using `abutils.tools.*` or `abutils.plots.*` are also misleading.

### The biggest design issue: namespace curation is implicit, not explicit

`abutils.tl`, `abutils.pp`, and `abutils.cl` are built with wildcard re-exports:

- `src/abutils/tl.py:25-30`
- `src/abutils/pp.py:25`
- `src/abutils/cl.py:26`

That has several bad effects:

- internal helpers leak into the public surface
- the exported API depends on whether underlying modules remembered to define `__all__`
- static analysis and docs become harder to trust
- it becomes easy for docs and code to drift

This is already visible in practice: `abutils.tl` exposes helpers like `sp` and `get_binary_path`, which are implementation details, while failing to expose some truly public algorithms.

### Recommended API shape

If the goal is scanpy-like ergonomics for interactive and pipeline use, I would make the public surface explicit and boring:

- `abutils.io.*`
  - file readers/writers and dataframe conversions
- `abutils.pp.*`
  - preprocessing and data-reduction steps
  - e.g. `merge_fastqs`, `deduplicate`, `reduction`
- `abutils.tl.*`
  - analysis/transformation algorithms
  - e.g. alignment, clustering, search, clonify, phylogeny
- `abutils.pl.*`
  - plotting only
- `abutils.cl.*`
  - color/palette helpers only
- top-level models only
  - `Sequence`, `Pair`, `Lineage`

Then add explicit `__all__` lists in those namespace modules and make the docs/examples test against that exact surface.

### API verdict

The high-level intended API is good. The actual current API is not coherent enough to serve as a stable notebook/pipeline interface. The structure is close to something good, but it needs explicit curation and contract tests.

## 3. Verified red-flag bugs and design problems

These are higher priority than style cleanup because they either break behavior today or make future change riskier than necessary.

### A. Verified behavioral bugs

#### 1. `Cluster.consensus` is broken

In `src/abutils/tools/cluster.py:103-108`, the property returns `self._consensus()` instead of `self._consensus`.

Impact:

- Accessing `cluster.consensus` raises `TypeError: 'Sequence' object is not callable`.
- The old consensus assertion in `tests/test_cluster.py` is commented out, so the bug is undetected.

#### 2. `Phylogeny.fasta_string` breaks when clustering is disabled

In `src/abutils/tools/phylo.py:281-288`, `sequences` is only initialized inside `if self.do_clustering:`. If clustering is disabled and root is `None` or `False`, `to_fasta(sequences)` uses an unbound local variable.

I verified this with a runtime smoke test: `phylogeny(..., cluster=False, root=False).fasta_string` raises `UnboundLocalError`.

Impact:

- Disabling clustering is not safe, even though the API advertises it.
- This is exactly the kind of flexibility regression that hurts pipeline users later.

#### 3. The ANARCI wrapper is effectively broken

`src/abutils/tools/anarci.py` has multiple issues:

- The header is mislabeled `# filename: phylo.py` (`src/abutils/tools/anarci.py:1-2`).
- `numbering()` passes `debug` positionally into `anarci_wrap()`, which means the boolean is bound to `anarci_path`, not `debug` (`src/abutils/tools/anarci.py:50-56`).
- `anarci_wrap()` asks `get_binary_path("hmmer")`, but `hmmer` is not an available binary in `src/abutils/bin.py:133-150`.
- It also sets the ANARCI executable path to `get_binary_path("hmmer")` (`src/abutils/tools/anarci.py:69-70`), which is clearly wrong.

Impact:

- This surface is almost certainly unusable in its current form.
- There are no tests covering it.

#### 4. `mmseqs_search()` default argument handling is likely wrong

`src/abutils/tools/search.py:193` always appends `--search-type {search_type}` to the command, but the default value is `None`.

Impact:

- The default call likely emits `--search-type None`, which MMseqs2 should reject.
- There are no tests for `mmseqs_search()` at all.

### B. Architecture and maintainability problems

#### 1. Monolithic modules still dominate the codebase

The repo still has very large mixed-responsibility files:

- `src/abutils/core/sequence.py`: 1729 lines
- `src/abutils/io.py`: 1497 lines
- `src/abutils/tools/preprocessing.py`: 1411 lines
- `src/abutils/core/pair.py`: 1279 lines
- `src/abutils/tools/msa.py`: 1230 lines
- `src/abutils/tools/cluster.py`: 1207 lines
- `src/abutils/tools/phylo.py`: 1055 lines

This matters because the problem is not just size. The abstractions are mixed:

- `Sequence` is both a data model and a home for translation/codon-optimization logic and file parsing/writing (`src/abutils/core/sequence.py:49-140`, `src/abutils/core/sequence.py:1091-1450`).
- `io.py` then re-imports and re-surfaces those `Sequence`-level readers and converters (`src/abutils/io.py:35-51`).
- `preprocessing.py` mixes filename parsing, FASTQ merging, AIRR-table transforms, clustering-driven reduction, and console UI (`src/abutils/tools/preprocessing.py` throughout).

This increases coupling and makes future refactors harder than they need to be.

#### 2. CLI wrappers rely heavily on string-built `shell=True` subprocesses

Representative examples:

- MAFFT: `src/abutils/tools/msa.py:516-525`
- MUSCLE: `src/abutils/tools/msa.py:655-658`
- MUSCLE v3: `src/abutils/tools/msa.py:939-955`
- VSEARCH clustering: `src/abutils/tools/cluster.py:530-540`
- MMseqs clustering: `src/abutils/tools/cluster.py:695-723`
- CD-HIT: `src/abutils/tools/cluster.py:840-843`
- FastTree: `src/abutils/tools/phylo.py:122-127`
- `mmseqs_search`: `src/abutils/tools/search.py:193-214`
- FASTQ merge wrappers: `src/abutils/tools/preprocessing.py:488-499`, `src/abutils/tools/preprocessing.py:612-655`

Impact:

- path quoting is inconsistent
- spaces/shell metacharacters are harder to reason about
- debugging and error reporting are weaker than they need to be
- commands are harder to unit test cleanly

For a library that intentionally wraps bundled binaries, this is one of the most important technical debts to pay down.

#### 3. Logging is still largely `print()`-driven

There are many library-level `print()` calls across clustering, phylogeny, preprocessing, search, pairwise alignment, and plotting helper code.

Impact:

- noisy notebook output
- poor integration into larger pipelines
- hard to redirect or structure logs
- inconsistent quiet/debug behavior

This is especially awkward for a library intended for both notebooks and automated pipelines.

#### 4. The public API boundary is not documented in code

Scanpy's README explicitly tells users to rely on documented public APIs, not internal modules. `abutils` currently does the opposite in practice: users often have to reach into internal module paths because the documented public namespaces are incomplete.

That is a maintainability trap. It guarantees more downstream breakage during future refactors.

## 4. Other relevant observations

### Dependency footprint is heavier than it needs to be

`pyproject.toml` includes plotting and test dependencies in the main runtime dependency set:

- `matplotlib`, `seaborn`, `python-circos`
- `pytest`

See `pyproject.toml:27-55`.

Impact:

- headless/pipeline users pay for plotting and test dependencies
- `import abutils` eagerly imports `cl`, which imports matplotlib/seaborn (`src/abutils/__init__.py:4-9`, `src/abutils/utils/color.py:28-31`)
- packaging is less modular than it could be

I would strongly consider extras such as:

- `abutils[plot]`
- `abutils[test]`
- possibly `abutils[phylo]` / `abutils[search]` if future optionalization is desired

### CI does not exercise macOS, despite shipping Darwin binaries

The CI matrix is Linux-only (`.github/workflows/pytest.yml:12-31`), but the repository ships Darwin-targeted binaries.

Impact:

- packaging or permission regressions for macOS binaries will not be caught in CI
- a package that bundles platform-specific executables should ideally test at least one job per supported platform

### Repository hygiene still has some noise

The tree includes committed built docs and legacy binary directories:

- `docs/source/_build/...`
- `src/abutils/binaries/_old/...`

These may not ship in the wheel because the build config is constrained, but they still add review noise and maintenance overhead.

## 5. Recommended next steps

### Highest priority

1. Fix the broken public API surface.
   - populate `abutils.pl`
   - decide whether preprocessing lives in `pp` or `tl`, then make docs and code agree
   - export `clonify` from the intended namespace
   - stop documenting `abutils.tools.*` and `abutils.plots.*` as public callable surfaces unless those packages are intentionally curated

2. Add import-contract tests.
   - These will catch most of the current user-facing breakages immediately and cheaply.

3. Fix the verified bugs.
   - `Cluster.consensus`
   - `Phylogeny.fasta_string` when `cluster=False`
   - ANARCI wrapper
   - `mmseqs_search()` default CLI construction

### Next priority

4. Add a real integration test tier.
   - tiny real AIRR workflow
   - tiny FASTQ merge workflow
   - tiny MMseqs search workflow
   - public plotting smoke tests through `abutils.pl`

5. Convert binary wrappers to `subprocess.run([...], check=True, ...)`.
   - preserve debug logging
   - return clearer exceptions
   - reduce shell-quoting brittleness

6. Split monolithic modules by responsibility.
   - `Sequence` model vs sequence IO vs sequence transforms
   - preprocessing FASTQ merge vs AIRR reduction
   - plotting data normalization vs rendering

### Nice to have

7. Move plotting and test dependencies into extras.
8. Add at least one macOS CI job if Darwin binaries are considered supported.
9. Remove committed built-doc artifacts and clearly archive or delete legacy `_old` binaries.

## Bottom line

`abutils` already has a respectable fast test suite and some genuinely useful binary-backed integration coverage. The immediate risk is not "there are no tests"; the immediate risk is that the tests are not enforcing the public API that the docs promise, and a few untested code paths are clearly broken today.

If this library is meant to be a stable notebook-and-pipeline utility layer for AIRR analysis, the right move is to treat the public namespace as a product, not an accident of star imports. Once that contract is explicit and tested, the rest of the refactor work becomes much safer.
