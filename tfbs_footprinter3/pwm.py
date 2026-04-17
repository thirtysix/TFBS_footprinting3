"""Position weight matrix (PWM) construction and scoring.

Extracted from tfbs_footprinter3.py (previously lines 772-822). These are
the hottest pure functions in the tool — tfbs_finder() calls PWM_scorer
once per (position × TF). Pinned by the unit tests in tests/test_pwm.py;
any refactor here must keep those tests green.
"""
from __future__ import annotations

import math


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
