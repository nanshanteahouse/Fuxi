"""Validate that requirements/*.txt files are thin wrappers (not stale pin lists).

A "thin wrapper" is a .txt file whose only non-comment line is ``-e .[xxx]``,
delegating to pyproject.toml as the single source of truth.  If someone
accidentally restores an old explicit-dependency requirements file, this test
catches it.
"""

from pathlib import Path

REQUIREMENTS_DIR = Path("requirements")

EXPECTED_THIN_WRAPPERS = [
    "requirements.txt",
    "requirements/base.txt",
    "requirements/all.txt",
    "requirements/rna.txt",
    "requirements/atac.txt",
    "requirements/spatial.txt",
    "requirements/bulk.txt",
    "requirements/paper.txt",
    "requirements/dev.txt",
]


def _is_thin_wrapper(path: Path) -> bool:
    """Return True if *path* exists and its only non-comment line starts with ``-e ``."""
    if not path.exists():
        return False
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    non_comment = [line for line in lines if not line.startswith("#") and line.strip()]
    return len(non_comment) == 1 and non_comment[0].strip().startswith("-e ")


class TestThinWrappers:
    """Quick structural checks — pure file I/O, <0.1 s."""

    def test_all_expected_thin_wrappers_exist(self):
        """Every expected thin-wrapper file exists and is non-empty."""
        for fname in EXPECTED_THIN_WRAPPERS:
            path = Path(fname)
            assert path.exists(), f"Missing: {fname}"
            assert path.stat().st_size > 0, f"Empty: {fname}"

    def test_all_thin_wrappers_have_correct_format(self):
        """Each thin-wrapper contains exactly one ``-e ...`` line."""
        for fname in EXPECTED_THIN_WRAPPERS:
            path = Path(fname)
            assert _is_thin_wrapper(path), (
                f"{fname} is not a thin wrapper — expected a single -e line"
            )
