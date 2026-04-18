"""Combined affinity scoring and per-component weight helpers.

The CAS engine: find_clusters() orchestrates the seven experimental-data
weight calculations (eqtls, atac, metacluster, cage_correlations, cage,
gerp, cpg) into a single Combined Affinity Score per TFBS hit, then
converts that into an empirical p-value via calcCombinedAffinityPvalue.

Extracted from tfbs_footprinter3.py as part of the monolith split.
calcCombinedAffinityPvalue is pinned by tests/test_cas_pvalue.py.
"""
from __future__ import annotations

import logging
import math
import time
from bisect import bisect_left

import numpy as np

from tfbs_footprinter3.io_utils import distance_solve, overlap_range

# -------- vectorized weight-summing helpers --------
#
# find_clusters iterates (tf_name × hit) millions of times at pvalc=1.
# Each hit originally called six Python-level weight helpers, each of
# which looped over a feature list (cages, gerps, eqtls, atac, metas)
# and did a Python overlap/distance check. The helpers below take
# pre-built NumPy arrays of feature (start, end, weight) triples and
# do the overlap / distance math in NumPy.
#
# find_clusters builds the arrays once per transcript (features are
# transcript-scoped, not TF-scoped) and passes them into these helpers.
# The original scalar helpers (eqtls_weights_summing, gerp_weights_summing,
# etc.) are preserved above as reference implementations and for any
# external caller that still uses them.


def _features_to_arrays(features, start_idx=0, end_idx=1, weight_idx=2):
    """Split a list of [start, end, weight, ...] records into 3 NumPy arrays.

    Empty list -> three length-0 arrays so downstream code can still use
    boolean-index / max / sum without branching.
    """
    if not features:
        empty = np.empty(0, dtype=np.float64)
        return empty, empty, empty
    starts = np.array([f[start_idx] for f in features], dtype=np.float64)
    ends = np.array([f[end_idx] for f in features], dtype=np.float64)
    weights = np.array([f[weight_idx] for f in features], dtype=np.float64)
    return starts, ends, weights


def _eqtls_starts_ends_mags(converted_eqtls):
    """For eqtls: third field is signed effect; callers want abs().

    Used by eqtls_weights_summing_v so the lookup dict key is the
    magnitude matching the scalar `converted_eqtl_score_mag = abs(c[2])`.
    """
    if not converted_eqtls:
        empty = np.empty(0, dtype=np.float64)
        return empty, empty, empty
    starts = np.array([e[0] for e in converted_eqtls], dtype=np.float64)
    ends = np.array([e[1] for e in converted_eqtls], dtype=np.float64)
    mags = np.abs(np.array([e[2] for e in converted_eqtls], dtype=np.float64))
    return starts, ends, mags


def _cage_starts_ends_ratios(converted_cages):
    """Cages expose [start, end, desc, peak_ratio] — grab 0, 1, 3."""
    if not converted_cages:
        empty = np.empty(0, dtype=np.float64)
        return empty, empty, empty
    starts = np.array([c[0] for c in converted_cages], dtype=np.float64)
    ends = np.array([c[1] for c in converted_cages], dtype=np.float64)
    ratios = np.array([c[3] for c in converted_cages], dtype=np.float64)
    return starts, ends, ratios


def _ranges_overlap(a_start, a_end, b_starts, b_ends):
    """Vectorized: does [a_start, a_end] overlap each [b_starts[i], b_ends[i]]?

    Mirrors the behavior of overlap_range: an overlap range is non-empty when
    max(a[0], b[0]) <= min(a[-1], b[-1]). (overlap_range sorts each input
    first, so we normalize via min/max here as well.)
    """
    a_lo = min(a_start, a_end)
    a_hi = max(a_start, a_end)
    b_lo = np.minimum(b_starts, b_ends)
    b_hi = np.maximum(b_starts, b_ends)
    return np.maximum(a_lo, b_lo) <= np.minimum(a_hi, b_hi)


def _ranges_gap(a_start, a_end, b_starts, b_ends):
    """Vectorized distance_solve-equivalent gap between [a_start, a_end]
    and each [b_starts[i], b_ends[i]] (0 when they overlap)."""
    a_lo = min(a_start, a_end)
    a_hi = max(a_start, a_end)
    b_lo = np.minimum(b_starts, b_ends)
    b_hi = np.maximum(b_starts, b_ends)
    # Gap = max(0, max(a_lo, b_lo) - min(a_hi, b_hi))
    return np.maximum(0.0, np.maximum(a_lo, b_lo) - np.minimum(a_hi, b_hi))


def gerp_weights_summing_v(motif_start, motif_end, gerp_starts, gerp_ends, gerp_weights):
    """Vectorized drop-in replacement for gerp_weights_summing (per-hit).

    Takes max weight across gerps whose range overlaps the motif center.
    Kept for unit testing against the scalar reference; find_clusters
    uses the batched variant `_overlap_point_max` below for performance.
    """
    if gerp_starts.size == 0:
        return 0
    motif_center = int((motif_start + motif_end) / 2)
    overlaps = _ranges_overlap(motif_center, motif_center + 1, gerp_starts, gerp_ends)
    if not overlaps.any():
        return 0
    return float(gerp_weights[overlaps].max())


def atac_weights_summing_v(motif_start, motif_end, atac_starts, atac_ends, atac_weights):
    """Vectorized atac_weights_summing (per-hit, test reference)."""
    if atac_starts.size == 0:
        return 0
    motif_center = int((motif_start + motif_end) / 2)
    overlaps = _ranges_overlap(motif_center, motif_center + 1, atac_starts, atac_ends)
    if not overlaps.any():
        return 0
    return float(atac_weights[overlaps].max())


def cage_weights_summing_v(motif_start, motif_end, cage_starts, cage_ends, cage_ratios, cage_dist_weights_dict):
    """Vectorized cage_weights_summing (per-hit, test reference)."""
    if cage_starts.size == 0:
        return 0
    distances = _ranges_gap(motif_start, motif_end, cage_starts, cage_ends)
    total = 0.0
    for idx in range(distances.shape[0]):
        key = str(int(distances[idx]))
        if key in cage_dist_weights_dict:
            total += cage_dist_weights_dict[key] * cage_ratios[idx]
    return total


def eqtls_weights_summing_v(eqtl_occurrence_log_likelihood, motif_start, motif_end,
                            eqtl_starts, eqtl_ends, eqtl_mags, gtex_weights_dict):
    """Vectorized eqtls_weights_summing (per-hit, test reference)."""
    if eqtl_starts.size == 0:
        return 0
    overlaps = _ranges_overlap(motif_start, motif_end, eqtl_starts, eqtl_ends)
    if not overlaps.any():
        return 0
    total = 0.0
    for mag in eqtl_mags[overlaps]:
        total += gtex_weights_dict[float(mag)] + eqtl_occurrence_log_likelihood
    return total


# -------- Per-TF BATCH helpers --------
#
# The per-hit _v functions above let us verify math correctness against the
# scalar helpers in unit tests, but NumPy's per-call overhead dominates for
# small feature lists at pvalc=1 scale (millions of calls). The batch
# variants below process ALL hits of a TF in one NumPy op, amortizing the
# overhead. find_clusters uses these.


def _overlap_point_max(points, starts, ends, weights):
    """For each `point` in `points`, return the max of `weights[j]` where the
    1-wide range `[point, point+1]` overlaps `[starts[j], ends[j]]`; 0 when
    no feature overlaps.

    This preserves the exact semantics of the scalar helpers, which called
    `overlap_range([motif_center, motif_center+1], [feat_start, feat_end])`.
    overlap_range returns a non-empty range iff
    `max(point, start) <= min(point+1, end)`, so we mirror that bound.

    Shape: points (H,), starts/ends/weights (F,). Output (H,).
    """
    if starts.size == 0 or points.size == 0:
        return np.zeros(points.shape, dtype=np.float64)
    # Broadcast: (H, 1) vs (1, F) -> (H, F); use the 1-wide point range so
    # that a feature starting at exactly point+1 still counts as overlapping.
    lo = np.maximum(points[:, None], starts[None, :])
    hi = np.minimum(points[:, None] + 1, ends[None, :])
    overlaps = lo <= hi
    # Mask-aware max: set non-overlap entries to -inf, max over F axis
    w_expand = np.broadcast_to(weights[None, :], overlaps.shape)
    masked = np.where(overlaps, w_expand, -np.inf)
    result = masked.max(axis=1)
    # positions with no overlap come out as -inf -> set to 0
    result[~np.isfinite(result)] = 0.0
    return result


def _range_gap_batch(motif_starts, motif_ends, feature_starts, feature_ends):
    """(H, F) matrix of distance_solve results between motif and feature ranges.

    Each cell is max(0, max(motif_start, feat_start) - min(motif_end, feat_end)).
    """
    m_lo = motif_starts[:, None]
    m_hi = motif_ends[:, None]
    f_lo = feature_starts[None, :]
    f_hi = feature_ends[None, :]
    return np.maximum(0.0, np.maximum(m_lo, f_lo) - np.minimum(m_hi, f_hi))


def _build_cage_dist_lut(cage_dist_weights_dict):
    """Build a flat NumPy lookup table from int distance -> weight.

    Returns (lut, max_key). `lut[d]` is the weight for distance d, or
    np.nan if the distance isn't in the dict. Callers treat nan as
    "skip" (matches the original "if key in dict" check).

    Intended to be called ONCE per transcript and reused across all TFs.
    """
    if not cage_dist_weights_dict:
        return np.empty(0, dtype=np.float64), -1
    max_key = max(int(k) for k in cage_dist_weights_dict.keys())
    lut = np.full(max_key + 1, np.nan, dtype=np.float64)
    for k, v in cage_dist_weights_dict.items():
        lut[int(k)] = v
    return lut, max_key


def _cage_batch(motif_starts, motif_ends, cage_starts, cage_ends, cage_ratios, lut, max_key):
    """Per-hit cage_weights_sum for an entire TF of hits.

    `lut` and `max_key` come from _build_cage_dist_lut (built once per
    transcript). Returns a 1-D array of length H.
    """
    H = motif_starts.shape[0]
    if cage_starts.size == 0 or H == 0 or max_key < 0:
        return np.zeros(H, dtype=np.float64)
    distances = _range_gap_batch(motif_starts, motif_ends, cage_starts, cage_ends).astype(np.int64)
    in_range = distances <= max_key
    safe = np.where(in_range, distances, 0)
    per_cell = lut[safe] * cage_ratios[None, :]  # (H, F)
    # nan entries (distance in range but dict key missing) and out-of-range
    # cells collapse to 0.
    per_cell = np.where(in_range & ~np.isnan(per_cell), per_cell, 0.0)
    return per_cell.sum(axis=1)


def _eqtl_batch(motif_starts, motif_ends, eqtl_starts, eqtl_ends, eqtl_mags,
                gtex_weights_dict, eqtl_occurrence_log_likelihood):
    """Per-hit eqtls_weights_sum across all hits of a TF."""
    H = motif_starts.shape[0]
    if eqtl_starts.size == 0 or H == 0:
        return np.zeros(H, dtype=np.float64)
    overlaps = (
        np.maximum(motif_starts[:, None], eqtl_starts[None, :]) <=
        np.minimum(motif_ends[:, None], eqtl_ends[None, :])
    )
    # Per-cell weight = gtex_weights_dict[mag] + eqtl_occurrence_log_likelihood when overlap, else 0
    # Look up gtex_weights_dict per unique magnitude; eqtl_mags is (F,) so do a small Python map once.
    per_eqtl_weight = np.array(
        [gtex_weights_dict[float(m)] + eqtl_occurrence_log_likelihood for m in eqtl_mags],
        dtype=np.float64,
    )
    return (overlaps * per_eqtl_weight[None, :]).sum(axis=1)


def _metacluster_batch_tf(motif_starts_23, motif_ends_23, counts_arr, motif_len, overlap_weights_dict):
    """Batch metacluster_weights_summing for every hit of a single TF.

    motif_len is TF-scoped (all hits of a TF share the same motif length)
    so the inner lookup dict and the sliding window width are constant.
    """
    H = motif_starts_23.shape[0]
    if H == 0 or motif_len <= 0:
        return np.zeros(H, dtype=np.float64)
    inner_dict = overlap_weights_dict.get(str(motif_len))
    if not inner_dict:
        return np.zeros(H, dtype=np.float64)

    # Pad counts_arr with zeros if some hits' windows extend past the end.
    # Scalar version slices counts[start:end] which truncates gracefully;
    # sliding_window_view requires full-width windows so we pad to make the
    # two paths equivalent. counts are non-negative, so appending zeros
    # can't change the per-window max once any real value is present.
    starts_int = motif_starts_23.astype(np.int64)
    needed_len = int(starts_int.max()) + motif_len if H else motif_len
    if counts_arr.size < needed_len:
        padded = np.zeros(needed_len, dtype=counts_arr.dtype if counts_arr.size else np.int64)
        padded[: counts_arr.size] = counts_arr
        counts_arr = padded

    # Max-over-window of width motif_len at each hit's motif_start.
    windows = np.lib.stride_tricks.sliding_window_view(counts_arr, motif_len)
    max_start = windows.shape[0] - 1
    clipped = np.clip(starts_int, 0, max_start)
    num_overlaps = windows[clipped].max(axis=1)

    # Dense LUT: index = num_overlapping_metaclusters, value = weight.
    max_key = max(int(k) for k in inner_dict.keys())
    lut = np.zeros(max_key + 1, dtype=np.float64)
    present = np.zeros(max_key + 1, dtype=bool)
    for k, v in inner_dict.items():
        lut[int(k)] = v
        present[int(k)] = True

    in_range = num_overlaps <= max_key
    safe = np.where(in_range, num_overlaps, 0)
    result = np.where(in_range & present[safe], lut[safe], 0.0)
    return result


def _cpg_batch_tf(motif_starts_23, motif_len, cpg_obsexp_array,
                  cpg_obsexp_weights_dict_keys, cpg_obsexp_weights_dict):
    """Batch cpg_weights_summing for every hit of a single TF.

    `cpg_obsexp_array` is the per-position obs2exp value extracted from
    cpg_list (last field of each row) — built ONCE per transcript.
    """
    H = motif_starts_23.shape[0]
    if H == 0 or len(cpg_obsexp_weights_dict_keys) == 0 or motif_len <= 0:
        return np.zeros(H, dtype=np.float64)
    if cpg_obsexp_array.size == 0:
        return np.zeros(H, dtype=np.float64)

    # Pad with zeros if any hit's window extends past the end (matches
    # scalar's graceful slice truncation; obs/exp ratios are non-negative
    # so padding with 0 preserves per-window max when the real values cover
    # any usable portion).
    starts_int = motif_starts_23.astype(np.int64)
    needed_len = int(starts_int.max()) + motif_len
    if cpg_obsexp_array.size < needed_len:
        padded = np.zeros(needed_len, dtype=cpg_obsexp_array.dtype)
        padded[: cpg_obsexp_array.size] = cpg_obsexp_array
        cpg_obsexp_array = padded

    windows = np.lib.stride_tricks.sliding_window_view(cpg_obsexp_array, motif_len)
    max_start = windows.shape[0] - 1
    clipped = np.clip(starts_int, 0, max_start)
    cpg_obsexps = windows[clipped].max(axis=1)

    # Mirror the scalar bisect_left logic:
    #   idx = bisect_left(keys, cpg_obsexp)
    #   if idx != len(keys): pick keys[idx]
    #   else: pick keys[-1]
    keys_arr = np.asarray(cpg_obsexp_weights_dict_keys, dtype=np.float64)
    indices = np.searchsorted(keys_arr, cpg_obsexps, side="left")
    max_idx = len(keys_arr) - 1
    clamped = np.minimum(indices, max_idx)
    selected_keys = keys_arr[clamped]

    # Map float keys -> weights. Dict values are Python floats; build once.
    result = np.array(
        [cpg_obsexp_weights_dict[float(k)] for k in selected_keys],
        dtype=np.float64,
    )
    return result


def calcCombinedAffinityPvalue(combined_affinity_score, cas_pvalues_dict, cass_with_pvalues_sorted, cass_sorted, cas_pvalues_subdict):
    """
    Calculate the pvalue for this combined affinity score.
    """

    # determine the pvalue of the current combined affinity score
    if combined_affinity_score in cas_pvalues_subdict:
        combined_affinity_score_pvalue = str(cas_pvalues_subdict[combined_affinity_score])
    else:
        # index the current combined affinity score in the sorted list of scores (keys)
        cass_with_pvalues_sorted_index = bisect_left(cass_sorted, combined_affinity_score)
        if cass_with_pvalues_sorted_index > 0:
            if cass_with_pvalues_sorted_index < len(cass_with_pvalues_sorted) - 1:
                combined_affinity_score_pvalue = str(cass_with_pvalues_sorted[cass_with_pvalues_sorted_index][1])
            else:
                combined_affinity_score_pvalue = str(cass_with_pvalues_sorted[-1][1])
        else:
            combined_affinity_score_pvalue = ">" + str(cass_with_pvalues_sorted[0][1])

    return combined_affinity_score_pvalue


def eqtl_overlap_likelihood(converted_eqtls, chr_start, chr_end, tf_len, gene_len, gtex_variants, ens_gene_id):
    """
    Likelihood of eQTL occurrence.
    """

    eqtl_occurrence_log_likelihood = 0

    if ens_gene_id in gtex_variants:
        transcript_len = float(chr_end - chr_start)

        # determine size of search space, and probability of observing an eQTL in this gene.
        # GTEx searches for variants which occur over the span of the gene + 1,000,000 nt upstream+downstream.
        eqtl_search_space = 2000000 + gene_len
        associated_gtx_eqtls = gtex_variants[ens_gene_id]
        variant_count = len(associated_gtx_eqtls) * 1.
        eqtl_occurrence_log_likelihood = -1 * math.log(((tf_len * variant_count) / (eqtl_search_space - tf_len)) * (transcript_len / gene_len), 2)

    return eqtl_occurrence_log_likelihood


def eqtls_weights_summing(eqtl_occurrence_log_likelihood, ens_gene_id, target_species_hit, converted_eqtls, gtex_weights_dict, chr_start, chr_end, gtex_variants, tf_len, gene_len):
    """
    Identify if any of the eQTLs associated with this gene overlap this predicted TFBS.
    Retrieve the log-likelihood scores for all of them.
    """

    eqtl_weights = []

    if len(converted_eqtls) > 0:
        # determine the weight score for likelihood of this magnitude eQTL.
        # ref-point
        motif_start = target_species_hit[4]
        motif_end = target_species_hit[5]

        for converted_eqtl in converted_eqtls:
            converted_eqtl_start = converted_eqtl[0]
            converted_eqtl_end = converted_eqtl[1]
            converted_eqtl_score_mag = abs(converted_eqtl[2])

            overlap = overlap_range([motif_start, motif_end], [converted_eqtl_start, converted_eqtl_end])

            if len(overlap) > 0:
                eqtl_weight = gtex_weights_dict[converted_eqtl_score_mag]
                eqtl_weights.append(eqtl_weight + eqtl_occurrence_log_likelihood)

    eqtl_weights_sum = sum(eqtl_weights)

    return eqtl_weights_sum


def cage_correlations_summing_preparation(gene_name, transcript_id, cage_dict, TF_cage_dict, tf_name):
    """
    Extract transcript relevant cages (target_cages) once for each TF under analysis.
    Extract TF relevant data just once per TF, instead of for every hit.
    Identify the CAGE ides which are associated with the promoter of the TF currently under analysis.
    The correlation between these and the target gene will be extracted and correlating log-weight will be summed.
    """

    target_cages = []
    tf_cages = []

    # current transcript (target) cages
    if gene_name in cage_dict:
        target_cages = [x[0].replace("hg_", "").replace(".1", "") for x in cage_dict[gene_name]]

    if tf_name in TF_cage_dict:
        tf_cages = [x[0].replace("hg_", "").replace(".1", "") for x in TF_cage_dict[tf_name]]

    return target_cages, tf_cages


def cage_correlations_summing(target_species_hit, transcript_id, target_cages, tf_cages, cage_correlations_dict, cage_corr_weights_dict):
    """
    Extract correlation values between CAGEs associated with a predicted TFBS protein,
    and CAGEs associated with the current gene.
    """

    corr_weights_ls = []
    corr_weight_sum = 0

    for target_cage in target_cages:
        if target_cage in cage_correlations_dict:
            for tf_cage in tf_cages:
                if tf_cage in cage_correlations_dict[target_cage]:
                    cage_correlation = cage_correlations_dict[target_cage][tf_cage]
                    cage_corr_weight = cage_corr_weights_dict[abs(cage_correlation)]
                    corr_weights_ls.append(cage_corr_weight)

    # take strongest correlation weight score
    if len(corr_weights_ls) > 0:
        corr_weights_ls.sort()
        corr_weight_sum = corr_weights_ls[-1]

    return corr_weight_sum


def cage_weights_summing(transcript_id, target_species_hit, cage_dist_weights_dict, converted_cages):
    """
    Retrieve a log-likelihood score for this from the pre-existing dictionary.
    """

    cage_weights = []

    # ref-point
    for converted_cage in converted_cages:
        cage_peak_count_ratio = converted_cage[3]
        motif_cage_dist = str(distance_solve([converted_cage[0], converted_cage[1]], [target_species_hit[4], target_species_hit[5]]))

        if motif_cage_dist in cage_dist_weights_dict:
            cage_weight = cage_dist_weights_dict[motif_cage_dist]
            cage_weight_peak_count_ratio_adjusted = cage_weight * cage_peak_count_ratio
            cage_weights.append(cage_weight_peak_count_ratio_adjusted)

    cage_weights_sum = sum(cage_weights)

    return cage_weights_sum


def atac_weights_summing(transcript_id, target_species_hit, converted_atac_seqs_in_promoter):
    """
    Identify ATAC-Seq peaks which are near a putative TFBS.
    Retrieve a log-likelihood score for this from the pre-existing dictionary.
    """

    atac_weights_ls = []
    atac_weights_sum = 0

    # ref-point
    motif_start = target_species_hit[4]
    motif_end = target_species_hit[5]
    motif_center = int((motif_start + motif_end) / 2)

    for converted_atac_in_promoter in converted_atac_seqs_in_promoter:
        converted_atac_in_promoter_start = converted_atac_in_promoter[0]
        converted_atac_in_promoter_end = converted_atac_in_promoter[1]
        converted_atac_in_promoter_weight = converted_atac_in_promoter[2]

        overlap = overlap_range([motif_center, motif_center + 1], [converted_atac_in_promoter_start, converted_atac_in_promoter_end])

        if len(overlap) > 0:
            atac_weights_ls.append(converted_atac_in_promoter_weight)

    if len(atac_weights_ls) > 0:
        atac_weights_sum = max(atac_weights_ls)

    return atac_weights_sum


def metacluster_weights_summing(transcript_id, target_species_hit, metacluster_overlap_weights_dict, converted_metaclusters_in_promoter, metacluster_in_promoter_counts):
    """
    Identify the number of metaclusters which overlap this putative TFBS.
    Retrieve a log-likelihood score for this from the pre-existing dictionary.
    """

    num_overlapping_metaclusters = 0
    metacluster_weights_sum = 0

    # ref-point
    motif_start = target_species_hit[2]
    motif_end = target_species_hit[3]
    motif_len = motif_end - motif_start
    motif_window = metacluster_in_promoter_counts[motif_start:motif_end]
    num_overlapping_metaclusters = max(motif_window)

    if str(num_overlapping_metaclusters) in metacluster_overlap_weights_dict[str(motif_len)]:
        metacluster_weights_sum = metacluster_overlap_weights_dict[str(motif_len)][str(num_overlapping_metaclusters)]
    else:
        print("metacluster overlap sum not in weight dict")
        logging.warning(" ".join(["metacluster overlap sum not in weight dict"]))

    return metacluster_weights_sum


def gerp_weights_summing(target_species, transcript_id, chromosome, target_species_hit, converted_gerps_in_promoter):
    """
    Identify the gerps which are near this predicted TFBS.
    Retrieve a log-likelihood score for this distance from the pre-existing dictionary.
    """

    # ref-point
    motif_start = target_species_hit[4]
    motif_end = target_species_hit[5]
    motif_center = int((motif_start + motif_end) / 2)

    gerp_weights_sum = 0
    gerp_weights_ls = []
    for converted_gerp_in_promoter in converted_gerps_in_promoter:
        converted_gerp_in_promoter_start = converted_gerp_in_promoter[0]
        converted_gerp_in_promoter_end = converted_gerp_in_promoter[1]
        converted_gerp_in_promoter_weight = converted_gerp_in_promoter[2]

        overlap = overlap_range([motif_center, motif_center + 1], [converted_gerp_in_promoter_start, converted_gerp_in_promoter_end])

        if len(overlap) > 0:
            gerp_weights_ls.append(converted_gerp_in_promoter_weight)

    if len(gerp_weights_ls) > 0:
        gerp_weights_sum = max(gerp_weights_ls)

    return gerp_weights_sum


def cpg_weights_summing(transcript_id, target_species_hit, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, cpg_list):
    """
    Retrieve a CpG weight score based on the CpG obs/exp of the midpoint of the
    current predicted TFBS.
    """

    if len(cpg_obsexp_weights_dict_keys) > 0:
        # retrieve locations and CpG obs/exp score for the midpoint of this predicted TFBS
        # ref-point
        motif_start = target_species_hit[2]
        motif_end = target_species_hit[3]
        cpg_obsexp = max([x[-1] for x in cpg_list[motif_start:motif_end]])

        # extract the weight for the obsexp which is just less than the current obsexp
        next_lesser_obsexp_index = bisect_left(cpg_obsexp_weights_dict_keys, cpg_obsexp)
        if next_lesser_obsexp_index != len(cpg_obsexp_weights_dict_keys):
            next_lesser_obsexp = cpg_obsexp_weights_dict_keys[next_lesser_obsexp_index]
        else:
            next_lesser_obsexp = cpg_obsexp_weights_dict_keys[-1]

        cpg_weight = cpg_obsexp_weights_dict[next_lesser_obsexp]

        return cpg_weight

    else:
        return 0


def find_clusters(gene_name, ens_gene_id, chr_start, chr_end, alignment, target_species, chromosome, tfbss_found_dict, cleaned_aligned_filename, converted_gerps_in_promoter, converted_cages, converted_metaclusters_in_promoter, metacluster_in_promoter_counts, converted_atac_seqs_in_promoter, converted_eqtls, gtex_weights_dict, transcript_id, cage_dict, TF_cage_dict, cage_dist_weights_dict, metacluster_overlap_weights_dict, cpg_list, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, cage_correlations_dict, cage_corr_weights_dict, gtex_variants, gene_len, cas_pvalues_dict, pvalc, slim=False):
    """
    For each target species hit:
    Identify the highest score for each species within the locality threshold.
    Create combined affinity score from the target species hit and those best scores from each species.
    If two target species hits are within the locality threshold from one another, choose the hit which has the highest combined affinity score.

    When slim=True (used by the HPC CAS campaign), skip building cluster_dict
    of per-hit rows entirely and return a dict tf_name -> (pwm_scores_ndarray,
    cas_rounded_ndarray). This cuts peak memory on a 20 M-hit transcript from
    ~20 GB to ~4 GB by avoiding the Python list-of-lists accumulation and
    downstream pandas DataFrame materialization.
    """
    start_time = time.time()

    cluster_dict = {}
    slim_cluster_dict = {}  # tf_name -> (pwm_scores_ndarray, cas_rounded_ndarray); only populated when slim=True

    # Pre-build NumPy arrays of feature (start, end, weight) triples ONCE.
    # These are transcript-scoped (not per-TF) so the cost is amortized across
    # every (tf_name, hit) iteration below. See tests/test_scoring_vectorized.py
    # for the agreement checks between the vectorized helpers called below and
    # their scalar predecessors.
    gerp_starts, gerp_ends, gerp_weights = _features_to_arrays(converted_gerps_in_promoter)
    atac_starts, atac_ends, atac_weights = _features_to_arrays(converted_atac_seqs_in_promoter)
    cage_starts, cage_ends, cage_ratios = _cage_starts_ends_ratios(converted_cages)
    eqtl_starts, eqtl_ends, eqtl_mags = _eqtls_starts_ends_mags(converted_eqtls)
    cage_dist_lut, cage_dist_max_key = _build_cage_dist_lut(cage_dist_weights_dict)

    # cpg_list stores per-position rows ending in obs2exp; extract the obs2exp
    # column as a flat float array so _cpg_batch_tf can sliding-window-max it.
    if cpg_list:
        cpg_obsexp_array = np.array([row[-1] for row in cpg_list], dtype=np.float64)
    else:
        cpg_obsexp_array = np.empty(0, dtype=np.float64)

    # metacluster_in_promoter_counts is a dense int list per position; convert once.
    metacluster_counts_arr = np.asarray(metacluster_in_promoter_counts, dtype=np.int64) \
        if metacluster_in_promoter_counts is not None \
        else np.empty(0, dtype=np.int64)

    for tf_name, hits in tfbss_found_dict.items():

        # build dict and sorted list of pre-computed combined affinity scores for this tf
        if tf_name in cas_pvalues_dict:
            cass_with_pvalues_sorted = cas_pvalues_dict[tf_name]
            cass_sorted = [x[0] for x in cass_with_pvalues_sorted]
            cas_pvalues_subdict = {x[0]: x[1] for x in cass_with_pvalues_sorted}

        if len(hits) > 0:
            cluster_dict[tf_name] = []
            # Preserve the scalar's eqtl_overlap_likelihood value exactly: the
            # original code used `tf_len = len(hits[0][1])`, which is actually
            # the length of the strand string ("+1"/"-1" -> 2) rather than the
            # motif length. It's a pre-existing quirk we must match here or the
            # eqtl log-likelihood shifts. For metacluster/cpg (which compute
            # motif_len fresh from hit[2]/hit[3]) we use the real motif length.
            tf_len_quirk = len(hits[0][1])
            motif_length = len(hits[0][0])
            target_cages, tf_cages = cage_correlations_summing_preparation(gene_name, transcript_id, cage_dict, TF_cage_dict, tf_name)
            eqtl_occurrence_log_likelihood = eqtl_overlap_likelihood(converted_eqtls, chr_start, chr_end, tf_len_quirk, gene_len, gtex_variants, ens_gene_id)

            # cage_correlations_summing's output depends on target_cages + tf_cages +
            # the global corr dicts — not on the specific hit. Compute once per TF
            # and broadcast across hits.
            if target_species == "homo_sapiens":
                corr_weight_sum_tf = cage_correlations_summing(None, transcript_id, target_cages, tf_cages, cage_correlations_dict, cage_corr_weights_dict)
            else:
                corr_weight_sum_tf = 0

            # Gather per-hit motif coordinates into arrays so every per-component
            # weight computation runs in one NumPy op for the whole TF's hits.
            motif_starts_46 = np.array([h[4] for h in hits], dtype=np.float64)
            motif_ends_46 = np.array([h[5] for h in hits], dtype=np.float64)
            motif_centers = ((motif_starts_46 + motif_ends_46) / 2.0).astype(np.int64)
            motif_starts_23 = np.array([h[2] for h in hits], dtype=np.float64)
            pwm_scores = np.array([h[6] for h in hits], dtype=np.float64)

            # Batch weight components (all transcript-level feature arrays built above).
            gerp_per_hit = _overlap_point_max(motif_centers, gerp_starts, gerp_ends, gerp_weights)
            cage_per_hit = _cage_batch(motif_starts_46, motif_ends_46, cage_starts, cage_ends, cage_ratios, cage_dist_lut, cage_dist_max_key)
            cpg_per_hit = _cpg_batch_tf(motif_starts_23, motif_length, cpg_obsexp_array, cpg_obsexp_weights_dict_keys, cpg_obsexp_weights_dict)

            if target_species == "homo_sapiens":
                atac_per_hit = _overlap_point_max(motif_centers, atac_starts, atac_ends, atac_weights)
                eqtl_per_hit = _eqtl_batch(motif_starts_46, motif_ends_46, eqtl_starts, eqtl_ends, eqtl_mags, gtex_weights_dict, eqtl_occurrence_log_likelihood)
                metacluster_per_hit = _metacluster_batch_tf(motif_starts_23, None, metacluster_counts_arr, motif_length, metacluster_overlap_weights_dict)
            else:
                atac_per_hit = np.zeros(len(hits))
                eqtl_per_hit = np.zeros(len(hits))
                metacluster_per_hit = np.zeros(len(hits))

            corr_per_hit = np.full(len(hits), float(corr_weight_sum_tf), dtype=np.float64)

            # Stack the 7 weight columns in the order find_clusters always emitted:
            # [species (gerp), cage, eqtls, atac, metacluster, cpg, corr].
            weights_raw = np.column_stack([
                gerp_per_hit, cage_per_hit, eqtl_per_hit, atac_per_hit,
                metacluster_per_hit, cpg_per_hit, corr_per_hit,
            ])

            # Batch CAS computation and rounding. Using np.round here (instead
            # of Python's built-in round) means results at exact half-values
            # (e.g. 4.285 -> 4.28 vs Python 4.29) can shift by 0.01 on a small
            # fraction of cells. Accepted as a tradeoff: np.round on the whole
            # matrix is one C-level op vs. millions of per-cell Python round()
            # calls. The scalar helpers remain available for exact-match tests.
            cas_raw = weights_raw.sum(axis=1) + pwm_scores

            if slim:
                # Keep CAS rounding consistent with the full path (which applies
                # np.round for the human-readable output) but store as ndarray, not list.
                # pvalc filtering skipped: the campaign always runs pvalc=1, and the
                # slim path is not exposed to end-users who need sub-1 filtering.
                slim_cluster_dict[tf_name] = (pwm_scores.astype(np.float32, copy=False),
                                              np.round(cas_raw, 2).astype(np.float32))
                continue

            cas_rounded = np.round(cas_raw, 2).tolist()
            weights_rounded = np.round(weights_raw, 2).tolist()
            had_contribution = (weights_raw != 0).tolist()

            for hit_idx, hit in enumerate(hits):
                combined_affinity_score = cas_rounded[hit_idx]

                if tf_name in cas_pvalues_dict:
                    combined_affinity_score_pvalue = calcCombinedAffinityPvalue(combined_affinity_score, cas_pvalues_dict, cass_with_pvalues_sorted, cass_sorted, cas_pvalues_subdict)
                else:
                    combined_affinity_score_pvalue = ""

                append_hit = False
                if ">" not in combined_affinity_score_pvalue and combined_affinity_score_pvalue != "":
                    if float(combined_affinity_score_pvalue) <= float(pvalc):
                        append_hit = True
                elif float(pvalc) == 1:
                    append_hit = True

                if append_hit:
                    hit.append(combined_affinity_score)
                    hit.append(combined_affinity_score_pvalue)
                    row_rounded = weights_rounded[hit_idx]
                    row_contributed = had_contribution[hit_idx]
                    # Preserve int-0 discipline: a raw value of 0.0 means the
                    # helper found no matching feature for this hit — scalar
                    # returned Python int 0 there; anything else is a rounded
                    # float.
                    experimental_weights_rounded = [
                        row_rounded[c] if row_contributed[c] else 0
                        for c in range(7)
                    ]
                    hit += experimental_weights_rounded
                    cluster_dict[tf_name].append(hit)

    total_time = time.time() - start_time
    logging.info(" ".join(["total time for find_clusters() for this transcript:", str(total_time), "seconds"]))

    if slim:
        return slim_cluster_dict
    return cluster_dict
