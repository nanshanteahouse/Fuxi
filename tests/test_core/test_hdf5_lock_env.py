"""Test HDF5 file-locking env var defense-in-depth.

Ensures ``os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")``
sets the default but does NOT overwrite an existing value.
All tests use ``monkeypatch`` to avoid leaking env changes.
"""

import os


def test_env_var_set_default(monkeypatch):
    """setdefault should set the value when the env var is absent."""
    monkeypatch.delenv("HDF5_USE_FILE_LOCKING", raising=False)
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    assert os.environ["HDF5_USE_FILE_LOCKING"] == "FALSE"


def test_env_var_not_overwritten(monkeypatch):
    """setdefault must NOT overwrite an already-set value."""
    monkeypatch.setenv("HDF5_USE_FILE_LOCKING", "TRUE")
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    assert os.environ["HDF5_USE_FILE_LOCKING"] == "TRUE"
