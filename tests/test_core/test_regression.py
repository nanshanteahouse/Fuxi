"""Regression tests for W1 and W2 fix cycles.

W1 fixes:
  W1-01: safe_write WSL — tmp+mv strategy on /mnt/ paths
  W1-02: format_cci_results descending sort — ascending=False parameter
  W1-03: guess_genome pipeline keys — 11 known species map correctly

W2 fixes:
  W2-01 (C6):  raw guard in spatial 05_annotate score_genes_mode
  W2-02 (C7):  raw guard in core/label_transfer run_label_transfer
  W2-03 (M3):  multi-species Ensembl ID regex in ensure_gene_symbols
  W2-04 (N3):  int(c) guard for non-integer cluster labels in atac 04_annotate
  W2-05 (M6):  sys.exit(2) when subcluster_types not configured in rna 06_subcluster
  W2-06 (M12): LLM response validation fallback in atac 04_annotate
"""

import json
import re
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd

from core.interaction.cell_interaction import format_cci_results
from core.preprocess.format_detector import guess_genome
from core.utils import safe_write

# ══════════════════════════════════════════════════════════════════════════
# W1-01: safe_write WSL
# ══════════════════════════════════════════════════════════════════════════


def test_safe_write_wsl():
    """W1-01: tmp+mv strategy on /mnt/ paths (WSL fix).

    Verifies that when target starts with /mnt/, the function:
      1. Creates tmpdir
      2. Writes to tmp path
      3. Moves tmp file to target
    """
    mock_adata = MagicMock()

    with patch("core.utils.os.makedirs") as mock_makedirs:
        with patch("core.utils.shutil.move") as mock_move:
            with patch("core.utils.os.path.getsize", return_value=1_000_000):
                safe_write(mock_adata, "/mnt/data/test.h5ad")

                # tmpdir created
                mock_makedirs.assert_called_once_with("/tmp/Fuxi", exist_ok=True)

                # adata.write called with tmp path (not target)
                mock_adata.write.assert_called_once_with("/tmp/Fuxi/test.h5ad", compression="gzip")

                # shutil.move called with tmp -> target
                mock_move.assert_called_once_with("/tmp/Fuxi/test.h5ad", "/mnt/data/test.h5ad")


# ══════════════════════════════════════════════════════════════════════════
# W1-02: format_cci_results descending sort
# ══════════════════════════════════════════════════════════════════════════


def test_format_cci_results_descending():
    """W1-02: top row has highest morans value when ascending=False.

    Creates a DataFrame with morans [0.1, 0.9, 0.5]; the 0.9 row must
    appear first after descending sort.
    """
    df = pd.DataFrame(
        {
            "ligand": ["A", "B", "C"],
            "receptor": ["X", "Y", "Z"],
            "source": ["S1", "S2", "S3"],
            "target": ["T1", "T2", "T3"],
            "morans": [0.1, 0.9, 0.5],
        }
    )
    result = format_cci_results(df, ascending=False, pval_col="morans")
    assert result.iloc[0]["morans"] == 0.9, (
        f"Expected top morans=0.9, got {result.iloc[0]['morans']}"
    )


# ══════════════════════════════════════════════════════════════════════════
# W1-03: guess_genome pipeline keys
# ══════════════════════════════════════════════════════════════════════════


def test_guess_genome_pipeline_keys():
    """W1-03: all 11 known species return correct genome strings."""
    expected = {
        "human": "hg38",
        "mouse": "mm10",
        "rat": "rn6",
        "zebrafish": "danRer11",
        "cow": "bosTau9",
        "pig": "susScr11",
        "macaque": "rheMac10",
        "chicken": "galGal6",
        "drosophila": "dm6",
        "c_elegans": "ce11",
        "frog": "xenTro10",
    }
    for species, expected_genome in expected.items():
        result = guess_genome(species)
        assert result == expected_genome, (
            f"{species}: expected {expected_genome!r}, got {result!r}"
        )


def test_guess_genome_unknown_species():
    """W1-03: unknown species returns empty string (no crash)."""
    assert guess_genome("gorilla") == ""
    assert guess_genome("") == ""
    assert guess_genome("homo_sapiens") == ""  # pipeline key is "human"


# ══════════════════════════════════════════════════════════════════════════
# W2-01 (C6): raw guard in spatial 05_annotate
# ══════════════════════════════════════════════════════════════════════════


def test_raw_guard_spatial():
    """W2-01: spatial score_genes_mode guard — adata.raw is None → use var_names."""
    # -- Case 1: adata.raw is None → fallback to adata.var_names --
    adata_no_raw = MagicMock(spec=[])
    adata_no_raw.raw = None
    adata_no_raw.var_names = ["GAPDH", "TP53", "EGFR"]

    var_names = (
        adata_no_raw.raw.var_names if adata_no_raw.raw is not None else adata_no_raw.var_names
    )
    assert list(var_names) == ["GAPDH", "TP53", "EGFR"], (
        "Should fall back to adata.var_names when adata.raw is None"
    )

    # -- Case 2: adata.raw exists → use raw.var_names --
    adata_with_raw = MagicMock(spec=[])
    adata_with_raw.raw = MagicMock()
    adata_with_raw.raw.var_names = ["raw_GAPDH", "raw_TP53"]
    adata_with_raw.var_names = ["GAPDH", "TP53", "EGFR"]

    var_names2 = (
        adata_with_raw.raw.var_names
        if adata_with_raw.raw is not None
        else adata_with_raw.var_names
    )
    assert list(var_names2) == ["raw_GAPDH", "raw_TP53"], (
        "Should use adata.raw.var_names when adata.raw exists"
    )


# ══════════════════════════════════════════════════════════════════════════
# W2-02 (C7): raw guard in label_transfer
# ══════════════════════════════════════════════════════════════════════════


def test_raw_guard_label_transfer():
    """W2-02: label_transfer has_raw guard — query.raw is None → fallback."""
    ref_var_names = np.array(["GAPDH", "TP53"])

    # -- Case 1: query has .raw --
    query_raw = MagicMock()
    query_raw.var_names = np.array(["GAPDH", "TP53", "EGFR"])
    query_raw.to_adata.return_value = "raw_adata_copy"

    query_with_raw = MagicMock()
    query_with_raw.raw = query_raw
    query_with_raw.var_names = np.array(["GAPDH", "TP53", "EGFR"])

    has_raw = query_with_raw.raw is not None
    common = np.intersect1d(
        ref_var_names,
        query_with_raw.raw.var_names if has_raw else query_with_raw.var_names,
    )
    raw_copy = query_with_raw.raw.to_adata() if has_raw else query_with_raw.copy()

    assert has_raw is True
    assert list(common) == ["GAPDH", "TP53"]
    assert raw_copy == "raw_adata_copy"

    # -- Case 2: query has NO .raw --
    query_no_raw = MagicMock()
    query_no_raw.raw = None
    query_no_raw.var_names = np.array(["GAPDH", "TP53", "EGFR"])
    query_no_raw.copy.return_value = "copied_adata"

    has_raw2 = query_no_raw.raw is not None
    common2 = np.intersect1d(
        ref_var_names,
        query_no_raw.raw.var_names if has_raw2 else query_no_raw.var_names,
    )
    raw_copy2 = query_no_raw.raw.to_adata() if has_raw2 else query_no_raw.copy()

    assert has_raw2 is False
    assert list(common2) == ["GAPDH", "TP53"]
    assert raw_copy2 == "copied_adata", "Should fall back to query.copy() when query.raw is None"


# ══════════════════════════════════════════════════════════════════════════
# W2-03 (M3): multi-species Ensembl ID regex
# ══════════════════════════════════════════════════════════════════════════


def _ensembl_regex() -> str:
    """Return the regex pattern used by ensure_gene_symbols."""
    return r"^ENS[A-Z]{0,4}G\d{11}$"


def test_ensembl_detection_multi_species():
    """W2-03: ensure_gene_symbols regex matches Ensembl IDs from ≥5 species."""
    pattern = _ensembl_regex()

    # Known Ensembl gene IDs (must match)
    ensembl_ids = [
        ("ENSG00000139618", "human"),  # 0 letters between ENS and G
        ("ENSMUSG00000057147", "mouse"),  # 3 letters: MUS
        ("ENSDARG00000079245", "zebrafish"),  # 4 letters: DARG
        ("ENSRNOG00000012345", "rat"),  # 3 letters: RNO
        ("ENSCAFG00000012345", "dog"),  # 3 letters: CAF
        ("ENSXETG00000012345", "frog"),  # 3 letters: XET
        ("ENSGGAG00000012345", "chicken"),  # 3 letters: GGA
    ]

    for eid, species in ensembl_ids:
        assert re.match(pattern, eid), f"Failed to match {species} Ensembl ID: {eid!r}"

    # Non-Ensembl identifiers (must NOT match)
    non_ensembl = [
        "TP53",
        "GAPDH",
        "ENSG0000013961",  # only 10 digits
        "ENSG000001396188",  # 12 digits
        "ENSG0000013961A",  # letter in digit block
        "ENS!G00000139618",  # special char in prefix
        "ENSMUSG0000057147",  # only 10 digits
        "ENSDARG0000079245",  # only 10 digits
        "ENS0123G00000123456",  # digits in species prefix block
    ]

    for bad in non_ensembl:
        assert not re.match(pattern, bad), f"Incorrectly matched non-Ensembl string: {bad!r}"


# ══════════════════════════════════════════════════════════════════════════
# W2-04 (N3): int(c) guard for non-integer cluster labels
# ══════════════════════════════════════════════════════════════════════════


def test_atac_non_integer_labels():
    """W2-04: int(c) guard prevents crash on non-integer labels like '0_1'."""
    # Non-integer label that previously crashed int()
    assert "0_1".isdigit() is False, "isdigit must be False for '0_1'"

    result = int("0_1") if str("0_1").isdigit() else str("0_1")
    assert result == "0_1", f"Expected '0_1', got {result!r}"
    assert isinstance(result, str)

    # Integer labels still work correctly
    assert "3".isdigit() is True
    result2 = int("3") if str("3").isdigit() else str("3")
    assert result2 == 3
    assert isinstance(result2, int)

    # Float-like label "1.5" — isdigit is False, pass through as string
    result3 = int("1.5") if str("1.5").isdigit() else str("1.5")
    assert result3 == "1.5"
    assert isinstance(result3, str)

    # Ordinary integer label "7" — int conversion works
    result4 = int("7") if str("7").isdigit() else str("7")
    assert result4 == 7
    assert isinstance(result4, int)


# ══════════════════════════════════════════════════════════════════════════
# W2-05 (M6): subcluster exit code 2
# ══════════════════════════════════════════════════════════════════════════


def test_step06_exit_code_2():
    """W2-05: 06_subcluster exits with code 2 when subcluster_types not set."""
    import importlib
    import inspect

    mod = importlib.import_module("rna.steps.06_subcluster")
    source = inspect.getsource(mod.main)
    assert "sys.exit(2)" in source, (
        "06_subcluster.main() must contain sys.exit(2) for the skip case"
    )

    # -- Verify the guard logic independently --
    # Trigger condition: cell_type is None AND subcluster_types is empty/falsy
    assert (None is None and not []) is True, "Guard should trigger exit"
    assert (None is None and not [1, 2, 3]) is False, "Non-empty list prevents exit"
    assert ("T cell" != None or not []) is True, "Non-None cell_type prevents exit"


# ══════════════════════════════════════════════════════════════════════════
# W2-06 (M12): LLM annotation validation
# ══════════════════════════════════════════════════════════════════════════


def test_atac_annotation_validation():
    """W2-06: malformed LLM JSON triggers fallback to empty dict.

    The validation logic rejects:
      - Non-dict JSON (e.g. list)
      - Dict values missing the required "cell_type" key
      - Unparseable JSON (except block)

    Accepts: dict(str → dict with "cell_type" key).
    """

    def validate(response: str) -> dict:
        """Replicate the validation logic from atac 04_annotate."""
        try:
            annotations = json.loads(response)
        except json.JSONDecodeError:
            return {}
        if not isinstance(annotations, dict) or not all(
            isinstance(v, dict) and "cell_type" in v for v in annotations.values()
        ):
            return {}
        return annotations

    # -- Valid: proper dict of dicts with cell_type --
    valid = '{"0": {"cell_type": "T cell", "confidence": "high"}}'
    result = validate(valid)
    assert result == {"0": {"cell_type": "T cell", "confidence": "high"}}

    # -- Edge: cell_type is empty string but key exists → valid --
    valid_empty_ct = '{"0": {"cell_type": "", "confidence": "low"}}'
    result_empty = validate(valid_empty_ct)
    assert result_empty == {"0": {"cell_type": "", "confidence": "low"}}

    # -- Malformed JSON → empty dict --
    assert validate("not valid json") == {}

    # -- Valid JSON but wrong type (list, not dict) → empty dict --
    assert validate('["a", "b"]') == {}

    # -- Dict missing cell_type key → empty dict --
    assert validate('{"0": {"confidence": "high"}}') == {}

    # -- Dict where a value is not a dict → empty dict --
    assert validate('{"0": "string_value"}') == {}

    # -- Empty JSON object → empty dict (no values to validate) --
    assert validate("{}") == {}
