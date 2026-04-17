"""Position weight matrix (PWM) construction and scoring.

Extracted from tfbs_footprinter3.py. These are the hottest pure
functions in the tool — tfbs_finder() previously called PWM_scorer
once per (position × TF). Pinned by the unit tests in tests/test_pwm.py;
any refactor here must keep those tests green.

`pwm_scan_sliding` is the NumPy-based fast path introduced in the
vectorization PR. It computes PWM scores at every sliding-window
position of a sequence in one shot, replacing the O(seq_length ×
motif_length) inner Python loop in tfbs_finder. `PWM_scorer` is kept
as a drop-in scalar reference (and is what the unit tests exercise).
"""
from __future__ import annotations

import math

import numpy as np


def pwm_maker(strand, motif_length, tf_motif, bg_nuc_freq_dict):
    """
    Make a PWM from a nucleotide frequency table.
    """
    pwm = [[], [], [], []]
    nuc_list = 'ACGT'
    # PWM according to http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647310/
    for i in range(0, motif_length):
        col = [tf_motif[0][i], tf_motif[1][i], tf_motif[2][i], tf_motif[3][i]]

        # number of sequences
        N = sum(col)

        # for each position (col) in the PFM.
        for j in range(0, len(tf_motif)):
            # nuc count at this position.
            nuc_count = tf_motif[j][i]

            # pseudo-count = sqrt(total number of samples).
            pseudo_count = 0.8

            # background frequency for this nucleotide in the promoter.
            nuc_bg = bg_nuc_freq_dict[nuc_list[j]]

            # probability of nuc
            nuc_probability = (nuc_count + pseudo_count / 4) / (N + pseudo_count)
            nuc_weight = math.log((nuc_probability / nuc_bg), 2)
            pwm[j].append(nuc_weight)

    pwm = pwm[:]

    return pwm


def PWM_scorer(seq, pwm, pwm_dict, pwm_type):
    """
    Generate score for current seq given a pwm.
    """

    seq_score = 0.0
    for i in range(0, len(seq)):
        seq_score += pwm[pwm_dict[seq[i:i + 1]]][i]

    return seq_score


def seq_to_int_array(seq: str, pwm_dict: dict) -> np.ndarray:
    """Map an ACGT sequence string to an int8 array of row indices.

    Positions with a character that isn't in `pwm_dict` are set to -1;
    callers should avoid scoring windows that contain -1 since such
    windows produce garbage PWM rows.
    """
    arr = np.full(len(seq), -1, dtype=np.int8)
    for char, idx in pwm_dict.items():
        arr[np.frombuffer(seq.encode("ascii"), dtype=np.uint8) == ord(char)] = idx
    return arr


def pwm_scan_sliding(seq_int: np.ndarray, pwm_array: np.ndarray) -> np.ndarray:
    """Return PWM scores at every sliding-window start position.

    scores[i] = sum_j pwm_array[seq_int[i+j], j]   for j = 0 .. motif_length - 1

    Equivalent to calling PWM_scorer on seq[i:i+motif_length] for every
    i in range(0, len(seq) - motif_length + 1), but does the sum in
    NumPy in a single fancy-indexing + sum-along-axis operation.

    `seq_int` must be integer row indices (use seq_to_int_array).
    `pwm_array` must be shape (4, motif_length), dtype float64.
    """
    motif_length = pwm_array.shape[1]
    n_positions = len(seq_int) - motif_length + 1
    if n_positions <= 0:
        return np.empty(0, dtype=np.float64)

    # Sliding windows over seq_int: shape (n_positions, motif_length)
    windows = np.lib.stride_tricks.sliding_window_view(seq_int, motif_length)
    # Fancy-index pwm_array at (row = windows[i, j], col = j) for all i, j.
    col_idx = np.arange(motif_length)
    # pwm_array[windows, col_idx] has shape (n_positions, motif_length);
    # NumPy broadcasts col_idx across the windows axis.
    return pwm_array[windows, col_idx].sum(axis=1)
