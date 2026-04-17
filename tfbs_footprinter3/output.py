"""Post-scoring output and ranking helpers.

Collects the handful of small functions that turn find_clusters' raw
output into ranked hits for plotting and into a CSV for the user.

Extracted from tfbs_footprinter3.py as part of the monolith split.
"""
from __future__ import annotations

import csv
from decimal import Decimal
from operator import itemgetter


def target_species_hits_table_writer(sorted_clusters_target_species_hits_list, output_table_name):
    """
    Write to table results sorted by combined affinity score.
    """

    with open(output_table_name, 'w') as output_table:
        writerUS = csv.writer(output_table)
        writerUS.writerow(['binding prot.', 'motif', 'strand', 'start', 'end', 'TSS-relative start', 'TSS-relative end', 'PWM score', 'p-value', 'combined\naffinity\nscore', 'combined\naffinity\nscore\np-value', 'species\nweights\nsum', 'cage\nweights\nsum', 'eqtls\nweights\nsum', 'atac\nweights\nsum', 'metacluster\nweights\nsum', 'cpg\nweight', 'corr.\nweight\nsum'])

        # for all results which have passed thresholds, write full result to .csv
        # ref-point

        if len(sorted_clusters_target_species_hits_list) > 0:
            for hit in sorted_clusters_target_species_hits_list:
                frame_score_pval_str = hit[8]
                combined_affinity_score_pval_str = hit[10]

                if ">" not in frame_score_pval_str and frame_score_pval_str != "":
                    if float(frame_score_pval_str) <= 0.0001:
                        hit[8] = f"{Decimal(frame_score_pval_str):.3e}"
                if ">" not in combined_affinity_score_pval_str and combined_affinity_score_pval_str != "":
                    if float(combined_affinity_score_pval_str) <= 0.0001:
                        hit[10] = f"{Decimal(combined_affinity_score_pval_str):.3e}"

                writerUS.writerow([str(x) for x in hit])


def sort_target_species_hits(cluster_dict):
    """
    Sort target_species hits which are part of a cluster by combined affinity score.
    """
    sorted_clusters_target_species_hits_list = []

    for tf_name, hits in cluster_dict.items():
        for hit in hits:
            sorted_clusters_target_species_hits_list.append([tf_name] + hit)

    # ref-point
    if len(sorted_clusters_target_species_hits_list) > 0:
        sorted_clusters_target_species_hits_list = sorted(sorted_clusters_target_species_hits_list, key=itemgetter(10), reverse=False)

    return sorted_clusters_target_species_hits_list


def top_x_greatest_hits(sorted_clusters_target_species_hits_list, top_x_tfs_count):
    """
    Identify the best scoring hits up to some threshold of number of tfs.
    Allows plotting more than one instance of a top tf, without increasing the total tf used count.
    e.g. 3 instances of KLF4 will count as only one tf used towards the top_x_tfs_count threshold.
    """

    # to keep track of how many tfs have been added
    top_x_tfs = []
    # to store the x greatest tfs and their locations
    top_x_greatest_hits_dict = {}

    # add all hits to single pool so top hits can be identified
    for sorted_clusters_target_species_hit in sorted_clusters_target_species_hits_list:
        # ref-point
        tf_name = sorted_clusters_target_species_hit[0]
        if len(top_x_tfs) < top_x_tfs_count:

            # keep track of what & how many tfs have been added
            if tf_name not in top_x_tfs:
                top_x_tfs.append(tf_name)

            # add the hit to the top hits if the count threshold has not been met
            if tf_name in top_x_greatest_hits_dict:
                top_x_greatest_hits_dict[tf_name].append(sorted_clusters_target_species_hit)

            else:
                top_x_greatest_hits_dict[tf_name] = [sorted_clusters_target_species_hit]

    return top_x_greatest_hits_dict


def top_greatest_hits(sorted_clusters_target_species_hits_list, top_x_tfs_count):
    """
    Identify the best scoring hits up to some threshold of number of tfs.
    Allows plotting more than one instance of a top tf, without increasing the total tf used count.
    e.g. 3 instances of KLF4 will count as only one tf used towards the top_x_tfs_count threshold.
    """

    # to keep track of how many tfs have been added
    top_x_tfs = []
    # to store the x greatest tfs and their locations
    top_x_greatest_hits_dict = {}

    added_count = 0

    # add all hits to single pool so top hits can be identified
    for sorted_clusters_target_species_hit in sorted_clusters_target_species_hits_list:
        # ref-point
        tf_name = sorted_clusters_target_species_hit[0]

        if added_count < 1000:
            # keep track of what & how many tfs have been added
            if tf_name not in top_x_tfs:
                top_x_tfs.append(tf_name)

            # add the hit to the top hits if the count threshold has not been met
            if tf_name in top_x_greatest_hits_dict:
                top_x_greatest_hits_dict[tf_name].append(sorted_clusters_target_species_hit)

            else:
                top_x_greatest_hits_dict[tf_name] = [sorted_clusters_target_species_hit]

            added_count += 1

    return top_x_greatest_hits_dict
