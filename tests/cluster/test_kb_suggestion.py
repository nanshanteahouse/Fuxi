"""Tests for KB-informed cluster count suggestion: ``suggest_target_n_clusters``.

Strategy C: species-priority with all-source median fallback.
Uses real KB data — no monkeypatching or mock data.
"""

from core.kb import suggest_target_n_clusters


class TestSuggestTargetNClusters:
    """KB-informed cluster count suggestion."""

    # ── Species-filtered median ──────────────────────────────────────

    def test_retina_mouse_species_filtered(self) -> None:
        """retina + mouse returns median of mouse-only sources (~39)."""
        result = suggest_target_n_clusters("retina", "mouse")
        # median(macosko2015=39, shekhar2016=15, tran2019=46)
        assert isinstance(result, int)
        assert 35 <= result <= 45

    def test_retina_zebrafish_species_filtered(self) -> None:
        """retina + zebrafish returns median of zebrafish-only sources (~45)."""
        result = suggest_target_n_clusters("retina", "zebrafish")
        # median(kolsch2021=30, hahn2023=60)
        assert isinstance(result, int)
        assert result == 45

    # ── Species fallback (all-source median) ─────────────────────────

    def test_retina_unknown_species_fallback(self) -> None:
        """retina + unknown species falls back to all-source tissue median (~60)."""
        result = suggest_target_n_clusters("retina", "unknown_species")
        # All-source median dominated by hahn2023's 13 multi-species entries at 60
        assert isinstance(result, int)
        assert result == 60

    # ── Unknown tissue ───────────────────────────────────────────────

    def test_unknown_tissue_returns_none(self) -> None:
        """Unknown tissue returns None regardless of species."""
        result = suggest_target_n_clusters("unknown_tissue", "mouse")
        assert result is None

    # ── Edge: no species at all ──────────────────────────────────────

    def test_retina_no_species_fallback(self) -> None:
        """retina with species=None returns all-source median (~60)."""
        result = suggest_target_n_clusters("retina")
        assert isinstance(result, int)
        assert result == 60

    def test_unknown_tissue_no_species_returns_none(self) -> None:
        """Unknown tissue with species=None returns None."""
        result = suggest_target_n_clusters("unknown_tissue")
        assert result is None
