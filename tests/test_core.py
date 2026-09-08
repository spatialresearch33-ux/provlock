import os
import numpy as np
import pandas as pd
import pytest

from provlock.core import ProvenanceLock, validate_groups


@pytest.fixture
def tmp_lock(tmp_path):
    return ProvenanceLock(base_dir=str(tmp_path / "provenance"), verbose=False)


def test_set_seed_reproducible(tmp_lock):
    tmp_lock.set_seed(42)
    a = np.random.rand(5)
    tmp_lock.set_seed(42)
    b = np.random.rand(5)
    assert np.allclose(a, b)


def test_save_data_locked_creates_file(tmp_lock):
    df = pd.DataFrame({"a": [1, 2, 3]})
    result = tmp_lock.save_data_locked(df, "test.csv")
    assert result.ok
    assert os.path.exists(result.path)


def test_save_data_locked_noop_on_identical_content(tmp_lock):
    df = pd.DataFrame({"a": [1, 2, 3]})
    r1 = tmp_lock.save_data_locked(df, "test.csv")
    mtime1 = os.path.getmtime(r1.path)
    r2 = tmp_lock.save_data_locked(df.copy(), "test.csv")
    mtime2 = os.path.getmtime(r2.path)
    assert r2.message == "unchanged"
    assert mtime1 == mtime2  # file was not rewritten


def test_save_data_locked_archives_on_change(tmp_lock):
    df1 = pd.DataFrame({"a": [1, 2, 3]})
    df2 = pd.DataFrame({"a": [1, 2, 999]})  # content genuinely changed
    tmp_lock.save_data_locked(df1, "test.csv")
    tmp_lock.save_data_locked(df2, "test.csv")
    archived = os.listdir(tmp_lock.archive_dir)
    assert len(archived) == 1
    assert archived[0].startswith("test.csv.")


def test_load_data_verified_matches(tmp_lock):
    df = pd.DataFrame({"a": [1, 2, 3]})
    tmp_lock.save_data_locked(df, "test.csv")
    loaded = tmp_lock.load_data_verified("test.csv")
    pd.testing.assert_frame_equal(df, loaded)


def test_load_data_verified_missing_file_raises(tmp_lock):
    with pytest.raises(FileNotFoundError):
        tmp_lock.load_data_verified("does_not_exist.csv")


def test_verify_lock_detects_external_tampering(tmp_lock):
    df = pd.DataFrame({"a": [1, 2, 3]})
    result = tmp_lock.save_data_locked(df, "test.csv")
    # simulate silent external modification of the saved file
    with open(result.path, "a") as f:
        f.write("\ntampered,row")
    check = tmp_lock.verify_lock("test.csv")
    assert check.ok is False


def test_verify_lock_ok_when_untouched(tmp_lock):
    df = pd.DataFrame({"a": [1, 2, 3]})
    tmp_lock.save_data_locked(df, "test.csv")
    check = tmp_lock.verify_lock("test.csv")
    assert check.ok is True


def test_validate_groups_passes_on_match():
    df = pd.DataFrame({"group": ["TNBC"] * 5 + ["Normal"] * 5})
    counts = validate_groups(df, "group", ["TNBC", "Normal"], verbose=False)
    assert counts == {"TNBC": 5, "Normal": 5}


def test_validate_groups_raises_on_mismatch():
    df = pd.DataFrame({"group": ["TNBC"] * 5 + ["Unexpected"] * 5})
    with pytest.raises(ValueError):
        validate_groups(df, "group", ["TNBC", "Normal"], verbose=False)


def test_freeze_requirements_creates_file(tmp_lock):
    result = tmp_lock.freeze_requirements()
    assert result.ok
    assert os.path.exists(result.path)
    with open(result.path) as f:
        content = f.read()
    assert len(content) > 0  # some packages are installed in any real env


# ---------------------------------------------------------------
# AnnData tests — skipped automatically if anndata isn't installed
# ---------------------------------------------------------------
anndata = pytest.importorskip("anndata")
import scipy.sparse as sp  # noqa: E402


def _make_test_adata(n_obs=50, n_vars=20, sparse=False, with_counts_layer=False):
    X = np.random.rand(n_obs, n_vars).astype(np.float64)
    if sparse:
        X = sp.csr_matrix(X)
    a = anndata.AnnData(X=X)
    if with_counts_layer:
        counts = sp.csr_matrix(
            np.random.poisson(5, size=(n_obs, n_vars)).astype(np.float64)
        )
        a.layers["counts"] = counts
    return a


def test_save_adata_locked_creates_file(tmp_lock):
    a = _make_test_adata()
    result = tmp_lock.save_adata_locked(a, "test.h5ad")
    assert result.ok
    assert os.path.exists(result.path)


def test_save_adata_locked_downcasts_counts_layer_to_int32(tmp_lock):
    """Verifies the actual production fix: sparse 'counts' layer -> int32,
    lossless size reduction, not a blanket float32 downcast."""
    a = _make_test_adata(with_counts_layer=True)
    tmp_lock.save_adata_locked(a, "test.h5ad", downcast_counts_int32=True)
    loaded = tmp_lock.load_adata_verified("test.h5ad")
    assert loaded.layers["counts"].dtype == np.int32
    # X itself should be untouched (default downcast_float32=False)
    assert loaded.X.dtype == np.float64


def test_save_adata_locked_optional_float32_downcast(tmp_lock):
    a = _make_test_adata()
    tmp_lock.save_adata_locked(a, "test.h5ad", downcast_float32=True)
    loaded = tmp_lock.load_adata_verified("test.h5ad")
    assert loaded.X.dtype == np.float32


def test_save_adata_locked_noop_on_identical_content(tmp_lock):
    a = _make_test_adata()
    r1 = tmp_lock.save_adata_locked(a, "test.h5ad")
    mtime1 = os.path.getmtime(r1.path)
    r2 = tmp_lock.save_adata_locked(a.copy(), "test.h5ad")
    assert r2.message == "unchanged"
    assert os.path.getmtime(r2.path) == mtime1


def test_load_adata_verified_roundtrip(tmp_lock):
    a = _make_test_adata()
    tmp_lock.save_adata_locked(a, "test.h5ad")
    loaded = tmp_lock.load_adata_verified("test.h5ad")
    assert loaded.n_obs == a.n_obs
    assert loaded.n_vars == a.n_vars
