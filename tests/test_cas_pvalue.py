"""Unit tests for calcCombinedAffinityPvalue — the CAS p-value lookup.

Consumes the per-TF slice of cas_pvalues_dict built by
data_loader.species_specific_data from a species' CAS threshold TSV.
Per-TF entry shape: {"scores": ndarray (ascending), "pvalues": ndarray
(descending, parallel to scores)}.
"""
from __future__ import annotations

import numpy as np

from tfbs_footprinter3.scoring import calcCombinedAffinityPvalue


def _build(scores_descending_by_p: list[tuple[float, float]]):
    """Build the per-TF slice calcCombinedAffinityPvalue expects.

    Input convention: p-value descending / score ascending rows (the TSV
    native order). Stored as parallel ndarrays sorted ascending by score.
    """
    arr = sorted(scores_descending_by_p, key=lambda pp: pp[1])  # by score asc
    return {
        "scores": np.array([s for _, s in arr], dtype=np.float64),
        "pvalues": np.array([p for p, _ in arr], dtype=np.float64),
    }


class TestCalcCombinedAffinityPvalue:
    def test_exact_hit_returns_its_pvalue(self):
        tf = _build([(1.0, 1.0), (0.5, 2.0), (0.1, 3.0), (0.01, 4.0)])
        assert calcCombinedAffinityPvalue(3.0, tf) == "0.1"

    def test_above_max_returns_tightest_pvalue(self):
        tf = _build([(1.0, 1.0), (0.5, 2.0), (0.1, 3.0), (0.01, 4.0)])
        # 5.0 is above the largest threshold (4.0); still get the p of the
        # tightest observable bin.
        assert calcCombinedAffinityPvalue(5.0, tf) == "0.01"

    def test_below_min_returns_prefixed_max_pvalue(self):
        tf = _build([(1.0, 1.0), (0.5, 2.0), (0.1, 3.0), (0.01, 4.0)])
        # 0.5 is below the smallest threshold (1.0): no row applies; emit
        # ">" + pvalue-at-first-row (here ">1.0") so callers can tell this
        # apart from a real empirical p-value.
        assert calcCombinedAffinityPvalue(0.5, tf) == ">1.0"

    def test_between_bins_returns_conservative_lower_bin(self):
        """This is the big semantic shift vs the legacy dense-JSON lookup.

        With a sparse TSV cutoff table, the supportable claim for a
        between-bin query is the LOOSER (larger) p-value of the lower
        cutoff. Old JSON-based lookup returned the next-higher bin's
        tighter p, which overstated significance.
        """
        tf = _build([(1.0, 1.0), (0.5, 2.0), (0.1, 3.0), (0.01, 4.0)])
        # 2.5 is above 2.0 (p=0.5) but below 3.0 (p=0.1). The strictest
        # claim supported by observed data at score=2.5 is p<=0.5.
        assert calcCombinedAffinityPvalue(2.5, tf) == "0.5"

    def test_exact_at_min_threshold(self):
        tf = _build([(1.0, 1.0), (0.5, 2.0), (0.1, 3.0), (0.01, 4.0)])
        # Hitting the smallest threshold exactly: that row applies.
        assert calcCombinedAffinityPvalue(1.0, tf) == "1.0"
