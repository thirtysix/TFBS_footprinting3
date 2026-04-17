"""Agreement tests: vectorized weight-summing helpers must match scalar ones.

find_clusters's per-hit Python inner loop — the one that dominates the
tool's runtime at pvalc=1 — was replaced by NumPy-vectorized variants
in scoring.py. These tests pin each vectorized helper against the
original scalar helper it replaces, for both synthetic edge cases and
a randomized realistic workload.
"""
from __future__ import annotations

import random

import numpy as np
import pytest

from tfbs_footprinter3.scoring import (
    _cage_starts_ends_ratios,
    _eqtls_starts_ends_mags,
    _features_to_arrays,
    atac_weights_summing,
    atac_weights_summing_v,
    cage_weights_summing,
    cage_weights_summing_v,
    eqtls_weights_summing,
    eqtls_weights_summing_v,
    gerp_weights_summing,
    gerp_weights_summing_v,
)


def _motif_hit(motif_start, motif_end, pwm_score=5.0, seq="ACGTACGTAC"):
    """Fabricate a target_species_hit shaped like tfbs_finder produces.

    Layout per tfbs_finder (and find_clusters): index 4/5 is motif_start/end
    on unaligned coords; index 2/3 is a TSS-relative position pair used by
    metacluster/cpg helpers. Here we set both to the same range so the
    weight helpers produce well-defined values either way.
    """
    return [seq, "+1", motif_start, motif_end, motif_start, motif_end, pwm_score, "0.01"]


class TestGerpVectorized:
    def test_empty_features(self):
        starts, ends, weights = _features_to_arrays([])
        assert gerp_weights_summing_v(100, 110, starts, ends, weights) == 0

    def test_single_overlapping_feature_takes_weight(self):
        hit = _motif_hit(100, 110)
        features = [[105, 115, 2.7]]
        scalar = gerp_weights_summing("homo_sapiens", "tid", "1", hit, features)
        s, e, w = _features_to_arrays(features)
        vec = gerp_weights_summing_v(100, 110, s, e, w)
        assert vec == pytest.approx(scalar)

    def test_no_overlap_returns_zero(self):
        hit = _motif_hit(100, 110)
        features = [[200, 220, 5.0]]
        assert gerp_weights_summing("homo_sapiens", "tid", "1", hit, features) == 0
        s, e, w = _features_to_arrays(features)
        assert gerp_weights_summing_v(100, 110, s, e, w) == 0

    def test_max_across_overlapping(self):
        hit = _motif_hit(100, 110)
        # motif_center = int(100 + 110/2) = 155 -- scalar uses that weird formula.
        # Make several features overlap and verify max-selection.
        features = [[150, 160, 1.0], [150, 160, 3.0], [100, 200, 2.0]]
        scalar = gerp_weights_summing("homo_sapiens", "tid", "1", hit, features)
        s, e, w = _features_to_arrays(features)
        vec = gerp_weights_summing_v(100, 110, s, e, w)
        assert vec == pytest.approx(scalar)


class TestAtacVectorized:
    def test_matches_scalar_at_random_workload(self):
        rng = random.Random(7)
        # 50 random features somewhere in [0, 2000]
        features = []
        for _ in range(50):
            s = rng.randint(0, 1900)
            features.append([s, s + rng.randint(5, 100), rng.uniform(0, 5)])
        s_arr, e_arr, w_arr = _features_to_arrays(features)
        for _ in range(30):
            motif_start = rng.randint(0, 1900)
            motif_end = motif_start + rng.randint(5, 20)
            hit = _motif_hit(motif_start, motif_end)
            scalar = atac_weights_summing("tid", hit, features)
            vec = atac_weights_summing_v(motif_start, motif_end, s_arr, e_arr, w_arr)
            assert vec == pytest.approx(scalar)


class TestCageVectorized:
    def test_matches_scalar_at_random_workload(self):
        rng = random.Random(42)
        # 30 cages and a weight-by-distance dict
        cages = []
        for _ in range(30):
            s = rng.randint(0, 1900)
            cages.append([s, s + rng.randint(5, 50), "desc", rng.uniform(0.1, 1.0)])
        weights = {str(d): rng.uniform(-2, 5) for d in range(0, 100)}
        s_arr, e_arr, r_arr = _cage_starts_ends_ratios(cages)
        for _ in range(30):
            motif_start = rng.randint(0, 1900)
            motif_end = motif_start + rng.randint(5, 20)
            hit = _motif_hit(motif_start, motif_end)
            scalar = cage_weights_summing("tid", hit, weights, cages)
            vec = cage_weights_summing_v(motif_start, motif_end, s_arr, e_arr, r_arr, weights)
            assert vec == pytest.approx(scalar)

    def test_empty_cages(self):
        s, e, r = _cage_starts_ends_ratios([])
        assert cage_weights_summing_v(10, 20, s, e, r, {"0": 1.0}) == 0


class TestEqtlsVectorized:
    def test_matches_scalar_at_random_workload(self):
        rng = random.Random(123)
        eqtls = [[rng.randint(0, 1900), rng.randint(1900, 2000), rng.choice([-2, -1, 1, 2])]
                 for _ in range(20)]
        # Deterministic eQTL lengths so overlap checks make sense:
        eqtls = [[e[0], e[0] + 5, e[2]] for e in eqtls]
        # gtex_weights_dict must cover all magnitudes used
        gtex_weights_dict = {1.0: 0.5, 2.0: 1.5}
        s_arr, e_arr, m_arr = _eqtls_starts_ends_mags(eqtls)

        for _ in range(25):
            motif_start = rng.randint(0, 1900)
            motif_end = motif_start + rng.randint(5, 20)
            hit = _motif_hit(motif_start, motif_end)
            log_lik = rng.uniform(-3, 3)
            scalar = eqtls_weights_summing(
                log_lik, "ens", hit, eqtls, gtex_weights_dict,
                chr_start=0, chr_end=10000, gtex_variants={}, tf_len=10, gene_len=5000
            )
            vec = eqtls_weights_summing_v(log_lik, motif_start, motif_end,
                                          s_arr, e_arr, m_arr, gtex_weights_dict)
            assert vec == pytest.approx(scalar)


class TestFeaturesToArrays:
    def test_empty(self):
        s, e, w = _features_to_arrays([])
        assert s.shape == (0,)
        assert e.shape == (0,)
        assert w.shape == (0,)

    def test_basic_triples(self):
        features = [[1, 2, 3.0], [4, 5, 6.0], [7, 8, 9.0]]
        s, e, w = _features_to_arrays(features)
        np.testing.assert_array_equal(s, [1, 4, 7])
        np.testing.assert_array_equal(e, [2, 5, 8])
        np.testing.assert_array_equal(w, [3, 6, 9])
