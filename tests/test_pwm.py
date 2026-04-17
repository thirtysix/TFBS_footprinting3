"""Unit tests for PWM_scorer and pwm_maker.

These functions sit in tfbs_footprinter3's hottest loop (tfbs_finder calls
PWM_scorer once per position x per TF). They are the primary target of the
upcoming NumPy vectorization pass — these tests are the safety net that
lets that refactor land without behavior drift.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from tfbs_footprinter3.pwm import pwm_scan_sliding, seq_to_int_array
from tfbs_footprinter3.tfbs_footprinter3 import PWM_scorer, pwm_maker

# ---------- PWM_scorer ----------


class TestPWMScorer:
    def test_all_zero_pwm(self, pwm_dict):
        """A zero PWM scores any sequence as 0.0."""
        motif_length = 4
        pwm = [[0.0] * motif_length for _ in range(4)]
        assert PWM_scorer("ACGT", pwm, pwm_dict, "fake") == 0.0

    def test_identity_contribution(self, pwm_dict):
        """Each position contributes the matrix value at (row=nuc, col=pos)."""
        # Design a PWM where only position 0 has non-zero values
        pwm = [
            [1.0, 0.0, 0.0],  # A row
            [2.0, 0.0, 0.0],  # C row
            [3.0, 0.0, 0.0],  # G row
            [4.0, 0.0, 0.0],  # T row
        ]
        assert PWM_scorer("AAA", pwm, pwm_dict, "fake") == 1.0
        assert PWM_scorer("CAA", pwm, pwm_dict, "fake") == 2.0
        assert PWM_scorer("GAA", pwm, pwm_dict, "fake") == 3.0
        assert PWM_scorer("TAA", pwm, pwm_dict, "fake") == 4.0

    def test_additive_over_positions(self, pwm_dict):
        """Total score is the sum over position contributions."""
        pwm = [
            [1.0, 0.5, 0.25],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
        assert PWM_scorer("AAA", pwm, pwm_dict, "fake") == pytest.approx(1.75)

    def test_mixed_sequence(self, pwm_dict):
        pwm = [
            [1.0, 2.0, 3.0],  # A
            [0.1, 0.2, 0.3],  # C
            [0.01, 0.02, 0.03],  # G
            [-1.0, -2.0, -3.0],  # T
        ]
        # ACG -> A@0=1.0 + C@1=0.2 + G@2=0.03 = 1.23
        assert PWM_scorer("ACG", pwm, pwm_dict, "fake") == pytest.approx(1.23)

    def test_single_position(self, pwm_dict):
        pwm = [[-5.0], [10.0], [0.0], [0.0]]
        assert PWM_scorer("C", pwm, pwm_dict, "fake") == 10.0
        assert PWM_scorer("A", pwm, pwm_dict, "fake") == -5.0


# ---------- pwm_maker ----------


class TestPwmMaker:
    """pwm_maker converts a PFM (counts matrix) to a PWM (log2 likelihood)."""

    def test_balanced_motif_with_uniform_bg(self, uniform_bg):
        """A perfectly balanced column (equal counts) at uniform bg -> 0 weight."""
        motif_length = 2
        # Each column has 25 of each nt -> 100 total -> equal probability 0.25
        tf_motif = [
            [25, 25],  # A
            [25, 25],  # C
            [25, 25],  # G
            [25, 25],  # T
        ]
        pwm = pwm_maker("+1", motif_length, tf_motif, uniform_bg)
        # (25 + 0.2) / (100 + 0.8) = 25.2/100.8 = 0.25 exactly -> log2(0.25/0.25) = 0
        for row in pwm:
            for val in row:
                assert val == pytest.approx(0.0, abs=1e-12)

    def test_shape_matches_motif(self, uniform_bg):
        motif_length = 5
        tf_motif = [[10] * motif_length for _ in range(4)]
        pwm = pwm_maker("+1", motif_length, tf_motif, uniform_bg)
        assert len(pwm) == 4
        for row in pwm:
            assert len(row) == motif_length

    def test_skewed_column_log_likelihood(self, uniform_bg):
        """Column with 100% A (plus pseudocount) produces expected log2 weights."""
        motif_length = 1
        tf_motif = [
            [100],  # all A
            [0],
            [0],
            [0],
        ]
        pwm = pwm_maker("+1", motif_length, tf_motif, uniform_bg)
        # With pseudo_count = 0.8 spread uniformly (0.2 per base):
        #   A_prob = (100 + 0.2) / (100 + 0.8) = 100.2 / 100.8
        #   C/G/T_prob = (0 + 0.2) / (100 + 0.8) = 0.2 / 100.8
        pseudo_count = 0.8
        N = 100
        expected_A = math.log((100 + pseudo_count / 4) / (N + pseudo_count) / 0.25, 2)
        expected_other = math.log((0 + pseudo_count / 4) / (N + pseudo_count) / 0.25, 2)
        assert pwm[0][0] == pytest.approx(expected_A)
        assert pwm[1][0] == pytest.approx(expected_other)
        assert pwm[2][0] == pytest.approx(expected_other)
        assert pwm[3][0] == pytest.approx(expected_other)

    def test_pwm_scorer_roundtrip(self, pwm_dict, uniform_bg):
        """PFM -> pwm_maker -> PWM_scorer is consistent with the log2 definition."""
        motif_length = 3
        tf_motif = [
            [80, 10, 30],  # A
            [10, 70, 30],  # C
            [5, 10, 30],  # G
            [5, 10, 10],  # T
        ]
        pwm = pwm_maker("+1", motif_length, tf_motif, uniform_bg)
        # Score the "best" motif (pick the argmax nucleotide at each col)
        # Col 0: A=80 (max), Col 1: C=70, Col 2: A=C=G=30 (tied; pick A)
        score_ACA = PWM_scorer("ACA", pwm, pwm_dict, "fake")
        # Sanity: best-case score should be >= score of a random sequence
        score_TTT = PWM_scorer("TTT", pwm, pwm_dict, "fake")
        assert score_ACA > score_TTT


# ---------- pwm_scan_sliding (NumPy vectorized) ----------


class TestPwmScanSliding:
    """Vectorized scanner must agree exactly with the scalar PWM_scorer
    at every window position — it's the hot path that replaces the
    per-position inner Python loop in tfbs_finder.
    """

    def test_empty_when_seq_shorter_than_motif(self, pwm_dict):
        seq_int = seq_to_int_array("AC", pwm_dict)
        pwm_array = np.array([[1.0, 2.0, 3.0]] * 4)
        scores = pwm_scan_sliding(seq_int, pwm_array)
        assert scores.shape == (0,)

    def test_single_window_matches_scalar(self, pwm_dict):
        pwm = [
            [1.0, 2.0, 3.0],
            [0.1, 0.2, 0.3],
            [0.01, 0.02, 0.03],
            [-1.0, -2.0, -3.0],
        ]
        pwm_array = np.array(pwm, dtype=np.float64)
        seq = "ACG"
        expected = PWM_scorer(seq, pwm, pwm_dict, "fake")
        scores = pwm_scan_sliding(seq_to_int_array(seq, pwm_dict), pwm_array)
        assert scores.shape == (1,)
        assert scores[0] == pytest.approx(expected)

    def test_matches_scalar_at_every_position(self, pwm_dict):
        """Run a realistic scan: random 200-nt sequence, 10-wide PWM.

        Every vectorized score must equal PWM_scorer on the same window
        to within float64 epsilon.
        """
        import random
        rng = random.Random(42)
        seq = "".join(rng.choice("ACGT") for _ in range(200))
        motif_length = 10
        pwm = [[rng.uniform(-3, 3) for _ in range(motif_length)] for _ in range(4)]
        pwm_array = np.array(pwm, dtype=np.float64)

        scores = pwm_scan_sliding(seq_to_int_array(seq, pwm_dict), pwm_array)

        # Compare against scalar PWM_scorer at each window
        expected = [
            PWM_scorer(seq[i:i + motif_length], pwm, pwm_dict, "fake")
            for i in range(len(seq) - motif_length + 1)
        ]
        assert len(scores) == len(expected)
        np.testing.assert_allclose(scores, expected, rtol=1e-12, atol=1e-12)

    def test_seq_to_int_array_maps_correctly(self, pwm_dict):
        arr = seq_to_int_array("ACGTACGT", pwm_dict)
        np.testing.assert_array_equal(arr, [0, 1, 2, 3, 0, 1, 2, 3])

    def test_seq_to_int_array_marks_unknown_chars(self, pwm_dict):
        arr = seq_to_int_array("AXN-", pwm_dict)
        assert arr[0] == 0  # 'A'
        assert arr[1] == -1  # 'X' not in dict
        assert arr[2] == -1  # 'N' not in dict
        assert arr[3] == -1  # '-' not in dict
