---
title: 'provlock: Checksum-locked, seed-locked, changelog-tracked provenance for spatial transcriptomics and omics pipelines'
tags:
  - Python
  - reproducibility
  - provenance
  - spatial transcriptomics
  - single-cell genomics
  - bioinformatics
authors:
  - name: S. Venkatesh
    affiliation: 1
affiliations:
  - name: "Integrative Cancer and Aging Research Laboratory (ICARL), Department of Physiology, Saveetha Medical College and Hospital, Saveetha Institute of Medical and Technical Sciences (SIMATS), Saveetha University, Thandalam - 602105, Chennai, Tamil Nadu, India"
    index: 1
date: 08 September 2026
bibliography: paper.bib
---

# Summary

`provlock` adds a reproducibility layer to spatial transcriptomics and
single-cell omics pipelines built on tools such as Scanpy [@wolf2018scanpy]
and Squidpy [@palla2022squidpy]. The package checksum-locks tabular data
(CSV / `pandas.DataFrame`) and `AnnData` (`.h5ad`) checkpoints, locks random
seeds for stochastic pipeline steps, validates group labels before
statistical testing, and records a changelog of when and why an
intermediate result changed. `provlock` does not replace Scanpy, Squidpy,
or similar analysis tools; it wraps the save and load steps that already
exist in a typical pipeline.

# Statement of need

Computational omics pipelines fail to reproduce for reasons that are often
subtle rather than catastrophic. A metadata file can change silently
between analysis sessions. A group label can shift without record. A
checkpoint can be overwritten with no note of what changed or why. These
problems typically surface at peer review, when a reviewer requests a
rerun and the original intermediate state can no longer be reconstructed
with confidence.

`provlock` originated from two such failures encountered during composite
spatial-index research on public Visium spatial transcriptomics data.
First, a pseudoreplication issue traced back to uncertainty over whether a
metadata file had changed between sessions. Second, `AnnData` checkpoints
reaching 1.5-2 GB per file, which made naive checksum-and-archive locking
storage-prohibitive at pipeline scale. The package's `AnnData`-locking
function addresses the second problem directly: when a sparse `counts`
layer is present, its values are downcast to `int32` before saving. Raw
sequencing counts are integers, so this downcast is lossless, not an
approximation, and combined with gzip compression it substantially reduces
checkpoint size. This exact approach was used in production on a real
multi-section Visium dataset before being generalized into this package.

Workflow managers such as Snakemake and Nextflow orchestrate pipeline
execution and dependency graphs. Data version control tools such as DVC
version large files against a remote store. Both solve a different problem
at a different scale. `provlock` targets the single-analyst or small-team
notebook workflow — Colab, Jupyter, or plain scripts — where a full
workflow manager or DVC setup is more infrastructure than the analysis
needs, but a silently changed metadata file or checkpoint still needs to
be caught before it reaches a result.

# Functionality

- `ProvenanceLock.save_data_locked` / `load_data_verified` lock CSV/
  DataFrame checkpoints by checksum. Identical content is not re-saved. A
  genuine change archives the previous version, never deletes it, and
  writes a timestamped changelog entry.
- `ProvenanceLock.save_adata_locked` / `load_adata_verified` provide the
  same guarantee for `AnnData` `.h5ad` checkpoints, with gzip compression
  and a lossless `int32` downcast of a sparse `counts` layer applied by
  default.
- `validate_groups` checks that group or category labels in a DataFrame
  match an expected set, and optionally expected per-group counts, before
  a statistical test runs on them.
- `ProvenanceLock.verify_lock` reports whether a previously locked file has
  changed outside the tracked workflow. It serves as a first diagnostic
  step when a result changes unexpectedly.
- `ProvenanceLock.set_seed` and `freeze_requirements` lock random seeds and
  snapshot installed library versions for a project.

Each function returns a result object, or a dictionary in the case of
`validate_groups`, rather than only a printed message. This allows the
checks to run inside a test suite or automated pipeline, not only
interactively in a notebook.

# Usage example

The example notebook (`examples/demo_lymph_node.ipynb`) runs the full
workflow on a public 10x Genomics Visium dataset, `V1_Human_Lymph_Node`
[@tenxgenomics2020lymphnode]
(4,035 spots), which is unrelated to the dataset used in the motivating
manuscript work. The raw checkpoint locks at 103.06 MB after gzip
compression. A set of established germinal-centre B-cell marker genes
[@desilva2015dynamics] — CD19 as a pan-B-cell marker, BCL6 and AICDA as
core regulators of the germinal-centre reaction, and MME, RGS13, and CD83
as commonly used germinal-centre B-cell markers — is scored across spots
as a non-novel demonstration signature, not as a biological finding. The
resulting zone labels are checked with `validate_groups`
(GC-high: 2,013 spots, GC-low: 2,012 spots). A simulated external
modification to the locked, scored checkpoint is correctly flagged by
`verify_lock`.

# Acknowledgements

Portions of this package's code and documentation were developed with the
assistance of generative AI tools. The author reviewed, tested, and
verified all code before release.

# References
