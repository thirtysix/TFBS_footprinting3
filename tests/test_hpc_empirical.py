"""Unit tests for hpc/build_cas_pvalues.py empirical p-value calculation.

This is the producer of the JSON that test_cas_pvalue.py's lookup consumes.
Keeps Stage D numerically pinned so we can evolve the pipeline without
silently shifting the p-value math.
"""
from __future__ import annotations

import pytest

from hpc.build_cas_pvalues import empirical_pvalues


class TestEmpiricalPvalues:
    def test_empty_input_returns_empty(self):
        assert empirical_pvalues([]) == []

    def test_monotonic_decreasing(self):
        """P(X >= s) is non-increasing as s increases."""
        scores = [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 5.0]
        result = empirical_pvalues(scores)
        pvs = [p for _, p in result]
        for a, b in zip(pvs[:-1], pvs[1:], strict=True):
            assert a >= b

    def test_sorted_ascending_by_score(self):
        """Output must be sorted by score ascending for the bisect_left lookup."""
        scores = [5.0, 1.0, 3.0, 2.0, 4.0]
        result = empirical_pvalues(scores)
        scores_out = [s for s, _ in result]
        assert scores_out == sorted(scores_out)

    def test_hand_worked_distribution(self):
        """Exact values for the 8-score fixture used in the CLI sanity check."""
        scores = [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 5.0]
        result = empirical_pvalues(scores)
        expected = {
            1.0: 8 / 8,  # all 8 are >= 1
            2.0: 7 / 8,  # 7 are >= 2
            3.0: 5 / 8,  # 5 are >= 3
            4.0: 4 / 8,  # 4 are >= 4
            5.0: 1 / 8,  # 1 is  >= 5
        }
        got = dict(result)
        for score, p in expected.items():
            assert got[score] == pytest.approx(p)

    def test_rounding_to_two_decimals(self):
        """Scores are rounded to 2 decimals so Stage D's output aligns with
        tfbs_footprinter3.py:1073 where CAS itself is rounded to 2 decimals."""
        # 1.234 and 1.236 both round to 1.23 and 1.24 respectively -> distinct
        scores = [1.234, 1.236, 1.241, 1.244]
        result = empirical_pvalues(scores)
        scores_out = [s for s, _ in result]
        # 1.234 -> 1.23; 1.236 -> 1.24; 1.241 -> 1.24; 1.244 -> 1.24
        assert 1.23 in scores_out
        assert 1.24 in scores_out
        # 1.23 should occur only once, 1.24 three times collapsed into one entry
        assert len([s for s, _ in result if s == 1.24]) == 1

    def test_single_value(self):
        """Degenerate case: one score -> p=1.0 at that score."""
        assert empirical_pvalues([2.5]) == [(2.5, 1.0)]

    def test_all_identical(self):
        """All-same distribution collapses to one entry with p=1.0."""
        result = empirical_pvalues([3.0, 3.0, 3.0, 3.0])
        assert result == [(3.0, 1.0)]
