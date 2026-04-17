"""Post-scoring output and ranking helpers.

Collects the handful of small functions that turn find_clusters' raw
output into ranked hits for plotting and into a CSV for the user.

Extracted from tfbs_footprinter3.py as part of the monolith split.
"""
from __future__ import annotations

import csv
from decimal import Decimal
from operator import itemgetter

_SMALL_PVAL_THRESHOLD = 0.0001
_CSV_HEADER = [
    'binding prot.', 'motif', 'strand', 'start', 'end',
    'TSS-relative start', 'TSS-relative end', 'PWM score', 'p-value',
    'combined\naffinity\nscore', 'combined\naffinity\nscore\np-value',
    'species\nweights\nsum', 'cage\nweights\nsum', 'eqtls\nweights\nsum',
    'atac\nweights\nsum', 'metacluster\nweights\nsum', 'cpg\nweight',
    'corr.\nweight\nsum',
]

# Parquet column names mirror the CSV header but drop the embedded newlines
# that only exist for display width in the rendered CSV. Keeping them here
# would force consumers to spell "combined\naffinity\nscore" to address the
# column, which is annoying in DataFrame queries.
_PARQUET_HEADER = [name.replace('\n', ' ') for name in _CSV_HEADER]
# Column indices of the 7 experimental-weight columns + the CAS score + PWM
# score. Cast to float64 before writing so Parquet gets a uniform numeric
# dtype (the rows mix int-0 and float-nonzero for weights).
_FLOAT_COLUMN_INDICES = (7, 9, 11, 12, 13, 14, 15, 16, 17)
_INT_COLUMN_INDICES = (3, 4, 5, 6)  # start, end, TSS-rel start/end


def _scientific_pvalue_if_small(s):
    """Format a p-value string as scientific notation iff it's <= 1e-4.

    Mirrors the in-place transformation the original per-row code did on
    hit[8] and hit[10] — factored out so the hot loop is simpler.
    """
    if s == "" or ">" in s:
        return s
    if float(s) <= _SMALL_PVAL_THRESHOLD:
        return f"{Decimal(s):.3e}"
    return s


def _apply_sci_pval_formatting(rows):
    """In-place scientific-notation reformat for the two p-value columns.

    Cached per unique input string — the tool reuses each PWM p-value
    threshold across millions of hits per TF, and the CAS p-value bins
    are a small discrete set, so the cache hit rate is near 100%.
    """
    sci_cache = {}

    def to_sci(s):
        cached = sci_cache.get(s)
        if cached is not None:
            return cached
        sci_cache[s] = out = _scientific_pvalue_if_small(s)
        return out

    for hit in rows:
        hit[8] = to_sci(hit[8])
        hit[10] = to_sci(hit[10])


def target_species_hits_table_writer_parquet(sorted_clusters_target_species_hits_list, output_table_name):
    """Write the sorted-clusters table as a Parquet file instead of CSV.

    Parquet is ~10x smaller on disk and ~10x faster to write at the
    ~2M-row scale we hit at pvalc=1. The HPC pipeline (hpc/cas_only.py)
    prefers Parquet when available: faster load + smaller on-disk
    footprint × 5000 transcripts × 123 species.

    Requires pyarrow; imported lazily so the CSV path has no extra dep
    surface at startup.
    """
    import pandas as pd  # noqa: PLC0415 — lazy to keep --help cold-start fast

    if not sorted_clusters_target_species_hits_list:
        pd.DataFrame(columns=_PARQUET_HEADER).to_parquet(
            output_table_name, engine="pyarrow", compression="snappy", index=False
        )
        return

    _apply_sci_pval_formatting(sorted_clusters_target_species_hits_list)

    df = pd.DataFrame(sorted_clusters_target_species_hits_list, columns=_PARQUET_HEADER)
    # Force uniform numeric dtypes so pyarrow writes fixed-width binary
    # columns rather than an object union; int-0 vs float-0.0 is only a
    # CSV-serialization concern and irrelevant to consumers that read the
    # columns via df['species weights sum'].
    for col_idx in _INT_COLUMN_INDICES:
        col = _PARQUET_HEADER[col_idx]
        df[col] = df[col].astype("int64")
    for col_idx in _FLOAT_COLUMN_INDICES:
        col = _PARQUET_HEADER[col_idx]
        df[col] = df[col].astype("float64")

    df.to_parquet(output_table_name, engine="pyarrow", compression="snappy", index=False)


def target_species_hits_table_writer(sorted_clusters_target_species_hits_list, output_table_name):
    """
    Write to table results sorted by combined affinity score.
    """
    with open(output_table_name, 'w') as output_table:
        writerUS = csv.writer(output_table)
        writerUS.writerow(_CSV_HEADER)

        if not sorted_clusters_target_species_hits_list:
            return

        # Hot loop: ~2M rows at pvalc=1. See _apply_sci_pval_formatting for
        # the cached Decimal conversion; csv.writer.writerows internally
        # calls str() per cell so we hand it the raw rows.
        _apply_sci_pval_formatting(sorted_clusters_target_species_hits_list)
        writerUS.writerows(sorted_clusters_target_species_hits_list)


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
