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

from tfbs_footprinter3.io_utils import distance_solve, overlap_range


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
    motif_center = int(motif_start + motif_end / 2)

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
    motif_center = int(motif_start + motif_end / 2)

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


def find_clusters(gene_name, ens_gene_id, chr_start, chr_end, alignment, target_species, chromosome, tfbss_found_dict, cleaned_aligned_filename, converted_gerps_in_promoter, converted_cages, converted_metaclusters_in_promoter, metacluster_in_promoter_counts, converted_atac_seqs_in_promoter, converted_eqtls, gtex_weights_dict, transcript_id, cage_dict, TF_cage_dict, cage_dist_weights_dict, metacluster_overlap_weights_dict, cpg_list, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, cage_correlations_dict, cage_corr_weights_dict, gtex_variants, gene_len, cas_pvalues_dict, pvalc):
    """
    For each target species hit:
    Identify the highest score for each species within the locality threshold.
    Create combined affinity score from the target species hit and those best scores from each species.
    If two target species hits are within the locality threshold from one another, choose the hit which has the highest combined affinity score.
    """
    start_time = time.time()

    cluster_dict = {}

    for tf_name, hits in tfbss_found_dict.items():

        # build dict and sorted list of pre-computed combined affinity scores for this tf
        if tf_name in cas_pvalues_dict:
            cass_with_pvalues_sorted = cas_pvalues_dict[tf_name]
            cass_sorted = [x[0] for x in cass_with_pvalues_sorted]
            cas_pvalues_subdict = {x[0]: x[1] for x in cass_with_pvalues_sorted}

        if len(hits) > 0:
            cluster_dict[tf_name] = []
            tf_len = len(hits[0][1])
            target_cages, tf_cages = cage_correlations_summing_preparation(gene_name, transcript_id, cage_dict, TF_cage_dict, tf_name)
            eqtl_occurrence_log_likelihood = eqtl_overlap_likelihood(converted_eqtls, chr_start, chr_end, tf_len, gene_len, gtex_variants, ens_gene_id)

            for hit in hits:
                # ref-point
                combined_affinity_score = 0
                target_species_hit = hit
                target_species_pwm_score = target_species_hit[6]
                species_weights_sum = 0
                cage_weights_sum = 0
                eqtls_weights_sum = 0
                atac_weights_sum = 0
                metacluster_weights_sum = 0
                corr_weight_sum = 0
                tf_len = len(hit[0])

                # datasets only available for homo sapiens
                # todo: build within function checking
                if target_species == "homo_sapiens":
                    eqtls_weights_sum = eqtls_weights_summing(eqtl_occurrence_log_likelihood, ens_gene_id, target_species_hit, converted_eqtls, gtex_weights_dict, chr_start, chr_end, gtex_variants, tf_len, gene_len)
                    atac_weights_sum = atac_weights_summing(transcript_id, target_species_hit, converted_atac_seqs_in_promoter)
                    metacluster_weights_sum = metacluster_weights_summing(transcript_id, target_species_hit, metacluster_overlap_weights_dict, converted_metaclusters_in_promoter, metacluster_in_promoter_counts)
                    corr_weight_sum = cage_correlations_summing(target_species_hit, transcript_id, target_cages, tf_cages, cage_correlations_dict, cage_corr_weights_dict)

                cage_weights_sum = cage_weights_summing(transcript_id, target_species_hit, cage_dist_weights_dict, converted_cages)
                species_weights_sum = gerp_weights_summing(target_species, transcript_id, chromosome, target_species_hit, converted_gerps_in_promoter)
                cpg_weight = cpg_weights_summing(transcript_id, target_species_hit, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, cpg_list)

                # calculate the complete score (combined affinity)
                experimental_weights = [species_weights_sum, cage_weights_sum, eqtls_weights_sum, atac_weights_sum, metacluster_weights_sum, cpg_weight, corr_weight_sum]
                combined_affinity_score += sum(experimental_weights) + target_species_pwm_score
                combined_affinity_score = round(combined_affinity_score, 2)

                if tf_name in cas_pvalues_dict:
                    combined_affinity_score_pvalue = calcCombinedAffinityPvalue(combined_affinity_score, cas_pvalues_dict, cass_with_pvalues_sorted, cass_sorted, cas_pvalues_subdict)
                else:
                    combined_affinity_score_pvalue = ""

                if ">" not in combined_affinity_score_pvalue and combined_affinity_score_pvalue != "":
                    if float(combined_affinity_score_pvalue) <= float(pvalc):
                        # append the combined affinity score and its pvalue
                        hit.append(combined_affinity_score)
                        hit.append(combined_affinity_score_pvalue)

                        # round all of the experimental weights to two places and append to hit
                        experimental_weights_rounded = [round(x, 2) for x in experimental_weights]
                        hit += experimental_weights_rounded

                        cluster_dict[tf_name].append(hit)
                elif float(pvalc) == 1:
                    # append the combined affinity score and its pvalue
                    hit.append(combined_affinity_score)
                    hit.append(combined_affinity_score_pvalue)

                    # round all of the experimental weights to two places and append to hit
                    experimental_weights_rounded = [round(x, 2) for x in experimental_weights]
                    hit += experimental_weights_rounded

                    cluster_dict[tf_name].append(hit)

    total_time = time.time() - start_time
    logging.info(" ".join(["total time for find_clusters() for this transcript:", str(total_time), "seconds"]))

    return cluster_dict
