"""
provlock.core
=============
Reproducibility / data-provenance toolkit for spatial transcriptomics and
omics pipelines. Originally developed as a set of Colab-notebook helper
functions during composite spatial-index research (checksum-locked data,
seed-locked stochastic steps, changelog-tracked reruns); generalized here
into a platform-agnostic, installable package.

Design goals:
- No hardcoded paths, no Colab-only imports (works in Colab, local Jupyter,
  Kaggle, HPC, or plain scripts)
- Every function returns a status object (bool / dict), not just prints,
  so results can be checked programmatically or in a test suite
- Archive-on-change, never delete-on-change
- AnnData (.h5ad) support with default size-reduction and a size warning
  before archiving, since raw checkpoints can be multi-GB
"""

import hashlib
import os
import random
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional

import numpy as np
import pandas as pd

try:
    import anndata as ad
    _HAS_ANNDATA = True
except ImportError:
    _HAS_ANNDATA = False

try:
    from scipy import sparse as sp_sparse
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False


@dataclass
class LockResult:
    """Return type for save/load/verify operations — check .ok programmatically."""
    ok: bool
    message: str
    path: Optional[str] = None
    size_mb: Optional[float] = None
    extra: dict = field(default_factory=dict)

    def __bool__(self):
        return self.ok


class ProvenanceLock:
    """
    One instance per project. Handles checksum-locked saving/loading of
    tabular data (CSV/DataFrame) and AnnData (.h5ad) checkpoints, with
    archive-on-change and a plain-text changelog.

    Example
    -------
    >>> pl = ProvenanceLock(base_dir="my_project/provenance")
    >>> pl.set_seed(42)
    >>> pl.save_data_locked(metadata_df, "metadata.csv")
    >>> df = pl.load_data_verified("metadata.csv")
    """

    def __init__(self, base_dir: str = "provenance", verbose: bool = True):
        self.base_dir = base_dir
        self.archive_dir = os.path.join(base_dir, "archive")
        self.checksum_file = os.path.join(base_dir, "checksums.txt")
        self.changelog_file = os.path.join(base_dir, "CHANGELOG.txt")
        self.verbose = verbose

        os.makedirs(self.base_dir, exist_ok=True)
        os.makedirs(self.archive_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # internal helpers
    # ------------------------------------------------------------------
    def _print(self, msg):
        if self.verbose:
            print(msg)

    def _checksum_of_file(self, filepath):
        h = hashlib.sha256()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()

    def _checksum_of_df(self, df):
        return hashlib.sha256(df.to_csv(index=False).encode()).hexdigest()

    def _read_checksum(self, name):
        if not os.path.exists(self.checksum_file):
            return None
        with open(self.checksum_file) as f:
            for line in f:
                if line.startswith(name + "="):
                    return line.strip().split("=", 1)[1]
        return None

    def _write_checksum(self, name, checksum):
        lines = []
        found = False
        if os.path.exists(self.checksum_file):
            with open(self.checksum_file) as f:
                lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith(name + "="):
                lines[i] = f"{name}={checksum}\n"
                found = True
                break
        if not found:
            lines.append(f"{name}={checksum}\n")
        with open(self.checksum_file, "w") as f:
            f.writelines(lines)

    def _log(self, message):
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(self.changelog_file, "a") as f:
            f.write(f"[{ts}] {message}\n")

    def _archive_existing(self, canonical_path, filename, size_warning_mb=500):
        """Copy old file to archive before overwrite. Warns if large."""
        size_mb = os.path.getsize(canonical_path) / (1024 * 1024)
        if size_mb > size_warning_mb:
            self._print(
                f"NOTE: archiving '{filename}' will copy {size_mb:.0f} MB "
                f"to {self.archive_dir}. Large checkpoints accumulate fast — "
                f"consider periodically pruning old archives if disk is limited."
            )
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        archive_name = f"{filename}.{ts}.bak"
        shutil.copy(canonical_path, os.path.join(self.archive_dir, archive_name))
        return archive_name, size_mb

    # ------------------------------------------------------------------
    # STEP A1: reproducible randomness
    # ------------------------------------------------------------------
    def set_seed(self, seed: int = 42) -> LockResult:
        np.random.seed(seed)
        random.seed(seed)
        msg = (f"Random seed locked to {seed}. Pass random_state={seed} to every "
               f"stochastic call (e.g. sc.tl.umap, sc.tl.leiden, sklearn estimators).")
        self._print(msg)
        return LockResult(ok=True, message=msg, extra={"seed": seed})

    # ------------------------------------------------------------------
    # Tabular data (CSV / DataFrame)
    # ------------------------------------------------------------------
    def save_data_locked(self, df: pd.DataFrame, filename: str) -> LockResult:
        """Save a DataFrame under a fixed name. No-ops if content is identical.
        Archives (never deletes) the old version if content changed, and logs it."""
        path = os.path.join(self.base_dir, filename)
        new_sum = self._checksum_of_df(df)
        old_sum = self._read_checksum(filename)

        if old_sum == new_sum:
            self._print(f"'{filename}': unchanged, not re-saved.")
            return LockResult(ok=True, message="unchanged", path=path)

        if os.path.exists(path):
            archive_name, size_mb = self._archive_existing(path, filename)
            self._log(f"'{filename}' changed -> archived as {archive_name}. NOTE WHY.")
            self._print(f"WARNING: '{filename}' changed since last save. "
                        f"Old version archived as {archive_name}.")
        else:
            self._log(f"'{filename}' created for the first time.")

        df.to_csv(path, index=False)
        self._write_checksum(filename, new_sum)
        size_mb = os.path.getsize(path) / (1024 * 1024)
        self._print(f"Saved: {path} ({size_mb:.2f} MB)")
        return LockResult(ok=True, message="saved", path=path, size_mb=size_mb)

    def load_data_verified(self, filename: str) -> pd.DataFrame:
        """Load a CSV previously saved with save_data_locked(); warns loudly on mismatch."""
        path = os.path.join(self.base_dir, filename)
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"'{filename}' not found in {self.base_dir}. "
                f"Save it with save_data_locked() first."
            )
        df = pd.read_csv(path)
        recomputed = self._checksum_of_df(df)
        saved = self._read_checksum(filename)
        if recomputed == saved:
            self._print(f"'{filename}' verified — matches last locked version.")
        else:
            self._print(f"WARNING: '{filename}' does NOT match the last locked "
                        f"version. Something changed outside this system. STOP "
                        f"and investigate before continuing analysis.")
        return df

    # ------------------------------------------------------------------
    # AnnData (.h5ad) — the new, size-aware extension
    # ------------------------------------------------------------------
    def save_adata_locked(
        self,
        adata,
        filename: str,
        compress: bool = True,
        compression_opts: int = 4,
        downcast_counts_int32: bool = True,
        downcast_float32: bool = False,
        sparsify: bool = False,
        drop_raw: bool = False,
    ) -> LockResult:
        """
        Save an AnnData checkpoint with checksum-lock + archive-on-change,
        same semantics as save_data_locked(), but sized for large omics data.

        Default size-reduction mirrors a verified production fix (brought a
        ~2GB TNBC Visium checkpoint down to MB scale with no precision loss
        on the counts data):
          - compress + compression_opts: gzip at the given level (default 4)
          - downcast_counts_int32: if a 'counts' layer exists and is sparse,
            cast its .data to int32 — raw counts are integers, so this is a
            lossless 4x reduction vs the float64 default, not an approximation
          - downcast_float32 / sparsify: OFF by default — these are blanket,
            lossy-in-spirit transforms not part of the verified fix; opt in
            explicitly if your data doesn't have a 'counts' layer to target
          - drop_raw: drop .raw before saving (off by default — opt in)
        """
        if not _HAS_ANNDATA:
            raise ImportError(
                "anndata is required for save_adata_locked(). "
                "Install with: pip install provlock[spatial]"
            )

        path = os.path.join(self.base_dir, filename)
        a = adata.copy()

        if drop_raw and a.raw is not None:
            a.raw = None

        if downcast_counts_int32 and "counts" in a.layers:
            c = a.layers["counts"]
            if _HAS_SCIPY and sp_sparse.issparse(c):
                c.data = c.data.astype(np.int32)
                a.layers["counts"] = c
            elif hasattr(c, "astype"):
                a.layers["counts"] = c.astype(np.int32)

        if downcast_float32:
            if hasattr(a.X, "dtype") and a.X.dtype != np.float32:
                a.X = a.X.astype(np.float32)
            for layer in list(a.layers.keys()):
                if layer == "counts":
                    continue  # already handled above as int32, don't re-cast to float
                if hasattr(a.layers[layer], "dtype") and a.layers[layer].dtype != np.float32:
                    a.layers[layer] = a.layers[layer].astype(np.float32)

        if sparsify and _HAS_SCIPY and not sp_sparse.issparse(a.X):
            density = np.count_nonzero(a.X) / a.X.size if a.X.size else 1
            if density < 0.5:
                a.X = sp_sparse.csr_matrix(a.X)

        compression = "gzip" if compress else None

        # checksum on the file bytes after write (content-based, matches CSV logic
        # in spirit: we hash what's actually on disk)
        tmp_path = path + ".tmp"
        if compression:
            a.write_h5ad(tmp_path, compression=compression, compression_opts=compression_opts)
        else:
            a.write_h5ad(tmp_path)
        new_sum = self._checksum_of_file(tmp_path)
        old_sum = self._read_checksum(filename)
        new_size_mb = os.path.getsize(tmp_path) / (1024 * 1024)

        if old_sum == new_sum:
            os.remove(tmp_path)
            self._print(f"'{filename}': unchanged, not re-saved.")
            return LockResult(ok=True, message="unchanged", path=path, size_mb=new_size_mb)

        if os.path.exists(path):
            archive_name, old_size_mb = self._archive_existing(path, filename)
            self._log(f"'{filename}' changed -> archived as {archive_name}. NOTE WHY.")
            self._print(f"WARNING: '{filename}' changed since last save. "
                        f"Old version ({old_size_mb:.0f} MB) archived as {archive_name}.")
        else:
            self._log(f"'{filename}' created for the first time.")

        shutil.move(tmp_path, path)
        self._write_checksum(filename, new_sum)
        self._print(f"Saved: {path} ({new_size_mb:.2f} MB, "
                    f"compression={'gzip' if compress else 'none'})")
        return LockResult(ok=True, message="saved", path=path, size_mb=new_size_mb)

    def load_adata_verified(self, filename: str):
        """Load an .h5ad checkpoint previously saved with save_adata_locked();
        warns loudly on checksum mismatch."""
        if not _HAS_ANNDATA:
            raise ImportError("anndata is required for load_adata_verified().")

        path = os.path.join(self.base_dir, filename)
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"'{filename}' not found in {self.base_dir}. "
                f"Save it with save_adata_locked() first."
            )
        current = self._checksum_of_file(path)
        saved = self._read_checksum(filename)
        if current == saved:
            self._print(f"'{filename}' verified — matches last locked version.")
        else:
            self._print(f"WARNING: '{filename}' does NOT match the last locked "
                        f"version. STOP and investigate before continuing analysis.")
        return ad.read_h5ad(path)

    # ------------------------------------------------------------------
    # STEP F: standalone debug check
    # ------------------------------------------------------------------
    def verify_lock(self, filename: str) -> LockResult:
        """Quick check: has this file silently changed since last save?
        Works for any file type (checksum is over raw file bytes)."""
        path = os.path.join(self.base_dir, filename)
        if not os.path.exists(path):
            msg = f"'{filename}' does not exist yet."
            self._print(msg)
            return LockResult(ok=False, message=msg, path=path)
        current = self._checksum_of_file(path)
        saved = self._read_checksum(filename)
        if current == saved:
            msg = f"'{filename}': OK, unchanged. Bug is likely in code, not data."
            self._print(msg)
            return LockResult(ok=True, message=msg, path=path)
        else:
            msg = f"'{filename}': CHANGED. Investigate data before touching code."
            self._print(msg)
            return LockResult(ok=False, message=msg, path=path)

    # ------------------------------------------------------------------
    # STEP G: environment freezing
    # ------------------------------------------------------------------
    def freeze_requirements(self) -> LockResult:
        """Snapshot exact installed library versions. Uses sys.executable
        (not a bare 'pip' call) so it always resolves to the active env."""
        path = os.path.join(self.base_dir, "requirements.txt")
        result = subprocess.run(
            [sys.executable, "-m", "pip", "freeze"],
            capture_output=True, text=True
        )
        with open(path, "w") as f:
            f.write(result.stdout)
        msg = f"Library versions frozen to {path}."
        self._print(msg)
        return LockResult(ok=True, message=msg, path=path)


# ------------------------------------------------------------------------
# STEP D: group-label validation (standalone — doesn't need file locking,
# kept as a plain function rather than a class method)
# ------------------------------------------------------------------------
def validate_groups(df, group_column, expected_groups, expected_counts=None, verbose=True):
    """
    Validate that a DataFrame's group column matches expected categories
    (and optionally expected per-group counts) before analysis proceeds.

    Returns a dict of actual counts. Raises ValueError on mismatch of the
    group SET (not counts — count mismatches are warned, not fatal, since
    a slightly different n is sometimes legitimate).
    """
    actual = set(df[group_column].dropna().unique())
    expected = set(expected_groups)

    if actual != expected:
        if verbose:
            print(f"MISMATCH — expected {expected}, found {actual}")
        raise ValueError("Group labels do not match expected groups. STOP.")

    counts = df[group_column].value_counts().to_dict()
    if verbose:
        print(f"Group labels OK: {counts}")

    if expected_counts:
        for g, n in expected_counts.items():
            actual_n = counts.get(g, 0)
            if verbose:
                status = "OK" if actual_n == n else "<-- MISMATCH"
                print(f"  '{g}': expected {n}, found {actual_n} {status}")

    n_missing = df[group_column].isna().sum()
    if n_missing > 0 and verbose:
        print(f"WARNING: {n_missing} rows missing a group label.")

    return counts
