# provlock

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23037252.svg)](https://doi.org/10.5281/zenodo.23037252)

**Checksum-locked, seed-locked, changelog-tracked provenance for spatial transcriptomics and omics pipelines.**

`provlock` grew out of real reproducibility failures encountered during spatial
transcriptomics research: metadata silently drifting between reruns, group
labels changing without notice, and — later — multi-gigabyte AnnData
checkpoints that made naive provenance-locking storage-prohibitive at scale.
It is designed to sit alongside Scanpy/Squidpy, not replace them: they handle
the analysis, `provlock` handles making sure you can trust and reproduce your
own intermediate results.

## What it does

- **`save_data_locked` / `load_data_verified`** — checksum-locked CSV/DataFrame
  checkpoints. Re-saving identical content is a no-op; genuine changes are
  archived (never deleted) and logged with a timestamp.
- **`save_adata_locked` / `load_adata_verified`** — the same guarantee for
  AnnData `.h5ad` checkpoints, with size-aware defaults: gzip compression and
  a lossless int32 downcast of a sparse `counts` layer (raw counts are
  integers — this alone cuts checkpoint size substantially with zero
  precision loss).
- **`validate_groups`** — confirms group/category labels match what you
  expect before any statistical test runs on them.
- **`verify_lock`** — a standalone check: "has this file changed since I last
  locked it?" Use it as a first step when debugging an unexpected result —
  it tells you whether to look at your data or your code.
- **`set_seed` / `freeze_requirements`** — reproducible randomness and a
  frozen library-version snapshot per project.

Every function returns a `LockResult` (or a plain dict for `validate_groups`)
so results can be checked programmatically, not just read off printed text.

## Install

Not yet on PyPI. For now, install directly from source:

```bash
pip install -e .
```

or, in Colab/Jupyter without a local clone:

```python
import sys
sys.path.insert(0, "/path/to/provlock/src")
```

(A running notebook kernel does not pick up a fresh `pip install -e` without
a restart — pointing `sys.path` directly at the source avoids that.)

For AnnData/h5ad support:

```bash
pip install -e ".[spatial]"
```

## Quick start

The example below is excerpted from
[`examples/demo_lymph_node.ipynb`](examples/demo_lymph_node.ipynb), run
against a real public dataset (`V1_Human_Lymph_Node`, 10x Genomics Visium,
4,035 spots) — every number here is actual output, not illustrative.

```python
from provlock import ProvenanceLock, validate_groups

pl = ProvenanceLock(base_dir="my_project/provenance")
pl.set_seed(42)
```
```
Random seed locked to 42. Pass random_state=42 to every stochastic call
(e.g. sc.tl.umap, sc.tl.leiden, sklearn estimators).
```

Lock a raw AnnData checkpoint:

```python
result = pl.save_adata_locked(adata, "lymph_node_raw.h5ad")
```
```
Saved: lymph_node_demo/provenance/lymph_node_raw.h5ad (103.06 MB, compression=gzip)
```

Validate group labels before running statistics on them:

```python
validate_groups(adata.obs, "zone", ["GC_high", "GC_low"])
```
```
Group labels OK: {'GC_high': 2013, 'GC_low': 2012}
```

Catch a silent external change to a locked file:

```python
check = pl.verify_lock("lymph_node_scored.h5ad")
```
```
'lymph_node_scored.h5ad': CHANGED. Investigate data before touching code.
```

## Why the `counts` layer gets special treatment

Raw AnnData checkpoints from spatial transcriptomics pipelines can reach
1.5–2 GB at full resolution — naively archiving a new copy on every change
makes provenance-locking storage-prohibitive fast. `provlock`'s default fix
targets the actual source of the bloat: if a `counts` layer is present and
sparse, its data is cast to `int32` (raw counts are integers — this is exact,
not an approximation), combined with gzip compression at level 4. This was
validated in production on a real multi-section Visium dataset before being
generalized into this package; `X` itself and other layers are left
untouched by default, since blanket downcasting isn't always appropriate for
every workflow.

```python
pl.save_adata_locked(adata, "checkpoint.h5ad")               # defaults: int32 counts + gzip
pl.save_adata_locked(adata, "checkpoint.h5ad", downcast_float32=True)  # opt-in, broader
```

## Design notes

- No hardcoded paths, no Colab-only imports — works in Colab, local Jupyter,
  Kaggle, HPC, or plain scripts. Pass `base_dir` explicitly per project.
- Archive-on-change, never delete-on-change.
- `freeze_requirements()` uses `sys.executable -m pip freeze`, not a bare
  `pip` call, so it always resolves to the active environment.

## Known issues

- Installing alongside some Colab environments may print a pip dependency
  resolver warning (e.g. a `pandas` version conflict with `google-colab`'s
  pinned version). This has not caused functional failures in testing, but
  if you see downstream pandas errors, pin your pandas version explicitly.

## AI-assisted development disclosure

Portions of this package's code and documentation were developed with the
assistance of generative AI tools. All code was reviewed, tested, and
verified by the author before inclusion.

## Testing

```bash
pip install -e ".[dev]"
pytest tests/
```

## How to cite

If you use `provlock` in your own work, please cite it as:

> Venkatesh, S. (2026). provlock: Checksum-locked, seed-locked,
> changelog-tracked provenance for spatial transcriptomics and omics
> pipelines (v0.1.0) [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.23037252

A machine-readable citation is also available in [`CITATION.cff`](CITATION.cff).

## Author

S. Venkatesh, Integrative Cancer and Aging Research Laboratory (ICARL),
Department of Physiology, Saveetha Medical College and Hospital, Saveetha
Institute of Medical and Technical Sciences (SIMATS), Saveetha University,
Thandalam - 602105, Chennai, Tamil Nadu, India.

## License

MIT
