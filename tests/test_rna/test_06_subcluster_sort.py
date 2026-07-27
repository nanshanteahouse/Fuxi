"""Test crash-safe cluster ID sort key for 06_subcluster.py."""


# The sort key used in 06_subcluster.py (extracted for testing)
def _safe_sort_key(x: str) -> tuple:
    """Cluster ID sort key: numeric by value if ASCII-digit, else push to end."""
    return (len(x), x) if (x.isascii() and x.isdigit()) else (999, x)


def test_numeric_cluster_ids_sort_correctly():
    """Numeric cluster IDs sort by numeric value, not lexicographic."""
    ids = ["0", "1", "10", "2"]
    result = sorted(ids, key=_safe_sort_key)
    assert result == ["0", "1", "2", "10"], f"Expected numeric sort, got {result}"


def test_non_numeric_cluster_ids_no_crash():
    """Non-numeric cluster IDs do not crash and sort to the end."""
    ids = ["0", "1", "abc"]
    result = sorted(ids, key=_safe_sort_key)
    # Numeric ones come first in order, non-numeric at end
    assert result[:2] == ["0", "1"], f"Numeric prefix wrong: {result}"
    assert result[2] == "abc", f"Non-numeric should be last: {result}"
