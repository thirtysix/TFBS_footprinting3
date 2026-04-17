"""Unit tests for calcCombinedAffinityPvalue — the CAS p-value lookup.

This is the consumer of the JSON produced by hpc/build_cas_pvalues.py and
shipped via S3 in Stage E of the HPC pipeline. Keeping its contract pinned
protects both the client-side lookup and the Stage D aggregation format.
"""
from __future__ import annotations

from tfbs_footprinter3.tfbs_footprinter3 import calcCombinedAffinityPvalue


def _build_lookup(scores_pvalues: list[tuple[float, float]]):
    """Replicate the structure find_clusters() builds per-TF before calling.

    See tfbs_footprinter3.py:1022-1026:
        cass_with_pvalues_sorted = cas_pvalues_dict[tf_name]
        cass_sorted = [x[0] for x in cass_with_pvalues_sorted]
        cas_pvalues_subdict = {x[0]: x[1] for x in cass_with_pvalues_sorted}
    """
    cass_with_pvalues_sorted = list(scores_pvalues)
    cass_sorted = [x[0] for x in cass_with_pvalues_sorted]
    cas_pvalues_subdict = {x[0]: x[1] for x in cass_with_pvalues_sorted}
    return cass_with_pvalues_sorted, cass_sorted, cas_pvalues_subdict


class TestCalcCombinedAffinityPvalue:
    def test_exact_hit_uses_subdict(self):
        """When the CAS is an exact key, the subdict value is returned directly."""
        cass_with_pvalues_sorted, cass_sorted, sub = _build_lookup([
            (1.0, 1.0), (2.0, 0.5), (3.0, 0.1), (4.0, 0.01)
        ])
        result = calcCombinedAffinityPvalue(2.0, None, cass_with_pvalues_sorted, cass_sorted, sub)
        assert result == "0.5"

    def test_above_max_returns_rightmost_pvalue(self):
        """A CAS above the largest stored score returns the right-edge p-value."""
        cass_with_pvalues_sorted, cass_sorted, sub = _build_lookup([
            (1.0, 1.0), (2.0, 0.5), (3.0, 0.1), (4.0, 0.01)
        ])
        # 5.0 is not in the table; bisect_left returns len(sorted)=4.
        # The current implementation returns cass_with_pvalues_sorted[-1][1] = 0.01
        # via the "else" branch (i < len-1 is false).
        result = calcCombinedAffinityPvalue(5.0, None, cass_with_pvalues_sorted, cass_sorted, sub)
        assert result == "0.01"

    def test_below_min_returns_prefixed_pvalue(self):
        """A CAS below the smallest stored score returns ">P(min)" (prefix ">")."""
        cass_with_pvalues_sorted, cass_sorted, sub = _build_lookup([
            (1.0, 1.0), (2.0, 0.5), (3.0, 0.1), (4.0, 0.01)
        ])
        result = calcCombinedAffinityPvalue(0.5, None, cass_with_pvalues_sorted, cass_sorted, sub)
        # bisect_left on 0.5 returns 0 -> the "else" branch ">..." path
        assert result == ">1.0"

    def test_between_values_interpolates_via_bisect(self):
        """Between two keys, bisect_left picks the insertion index; value at that
        insertion-index position is returned."""
        cass_with_pvalues_sorted, cass_sorted, sub = _build_lookup([
            (1.0, 1.0), (2.0, 0.5), (3.0, 0.1), (4.0, 0.01)
        ])
        # 2.5 -> bisect_left returns 2 (insertion position for 2.5 between 2.0 and 3.0)
        # Since index=2 > 0 and < len-1 = 3, returns cass_with_pvalues_sorted[2][1] = 0.1
        result = calcCombinedAffinityPvalue(2.5, None, cass_with_pvalues_sorted, cass_sorted, sub)
        assert result == "0.1"

    def test_empty_distribution_is_undefined(self):
        """Empty lookup causes a crash today — document it so refactors preserve
        the contract: find_clusters() only calls this when tf_name is in
        cas_pvalues_dict (tfbs_footprinter3.py:1068), so the list is non-empty."""
        # Not asserting behavior; just pin the current pre-condition.
        # If a future refactor tries to support empty, that's fine — update this.
        pass
