"""Post-scoring output and ranking helpers.

Collects the handful of small functions that turn find_clusters' raw
output into ranked hits for plotting and into a CSV for the user.

Extracted from tfbs_footprinter3.py as part of the monolith split.
"""
from __future__ import annotations

import csv
from decimal import Decimal
from operator import itemgetter

import numpy as np

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

    Pass through the prefixed sentinels:
      ">..."  -- query below min observed (less significant than anything
                 in the empirical table)
      "<..."  -- query above max observed (more significant than the
                 tightest observed p)
      ""      -- no p-value available
    """
    if s == "" or s.startswith(">") or s.startswith("<"):
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


_PYARROW_MISSING_HINT = (
    "The 'parquet' output format requires pyarrow. Install it with:\n"
    "    pip install 'tfbs_footprinting3[parquet]'\n"
    "or\n"
    "    pip install pyarrow\n"
    "Alternatively pass -of csv (the default) to keep the CSV output."
)


def _require_pyarrow():
    """Raise a helpful ImportError if pyarrow is missing.

    pyarrow is an optional dependency — only callers that want Parquet
    output need to pay the ~60 MB install cost. We check here (rather
    than at module import) so CSV users never see the error.
    """
    try:
        import pyarrow  # noqa: F401 — just checking availability
    except ImportError as exc:
        raise ImportError(_PYARROW_MISSING_HINT) from exc


def target_species_hits_table_writer_parquet(sorted_clusters_target_species_hits_list, output_table_name):
    """Write the sorted-clusters table as a Parquet file instead of CSV.

    Parquet is ~10x smaller on disk and ~10x faster to write at the
    ~2M-row scale we hit at pvalc=1. The HPC pipeline (hpc/cas_only.py)
    prefers Parquet when available: faster load + smaller on-disk
    footprint × 5000 transcripts × 123 species.

    pyarrow is an OPTIONAL dependency — install with
    `pip install tfbs_footprinting3[parquet]` to enable this output
    path. A missing pyarrow produces a pointed error from this function
    rather than a cryptic pandas traceback.
    """
    _require_pyarrow()
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


def target_species_hits_table_writer_slim_parquet(slim_cluster_dict, output_table_name):
    """Write a 3-column (tf_name, PWM score, combined affinity score) parquet.

    Purpose: the non-human CAS campaign (see hpc/puhti/) only needs the CAS
    score distribution per TF to build empirical p-value tables. The full
    18-column `sortedclusters` table is ~320 MB per transcript and
    requires building a pandas DataFrame + python list-of-lists of 20 M
    rows, which peaks at ~20 GB of process memory. The slim path skips
    that entire accumulation: it receives per-TF numpy arrays directly
    from find_clusters and streams them TF-by-TF into a pyarrow parquet
    via ParquetWriter, so peak memory is bounded by one TF's arrays.

    slim_cluster_dict: dict {tf_name: (pwm_scores_ndarray, cas_rounded_ndarray)}
    output_table_name: path to the output parquet.

    Schema matches the subset of _PARQUET_HEADER that hpc/puhti/build_cas_distributions.py
    actually reads (`binding prot.`, `combined affinity score`), plus PWM
    score for QC. No p-value columns (would all be empty in campaign mode).
    """
    _require_pyarrow()
    import pyarrow as pa  # noqa: PLC0415
    import pyarrow.parquet as pq  # noqa: PLC0415

    schema = pa.schema([
        ("binding prot.", pa.string()),
        ("PWM score", pa.float32()),
        ("combined affinity score", pa.float32()),
    ])

    writer = pq.ParquetWriter(output_table_name, schema, compression="snappy")
    try:
        if not slim_cluster_dict:
            # Write a single empty batch so readers see a valid empty parquet.
            writer.write_table(pa.Table.from_arrays(
                [pa.array([], type=pa.string()),
                 pa.array([], type=pa.float32()),
                 pa.array([], type=pa.float32())],
                schema=schema,
            ))
            return

        for tf_name, (pwm_arr, cas_arr) in slim_cluster_dict.items():
            n = int(pwm_arr.size)
            if n == 0:
                continue
            # pa.array copies the numpy buffer once; cheaper than pandas
            # roundtrip and doesn't double-materialize into a list.
            table_chunk = pa.Table.from_arrays(
                [
                    pa.array([tf_name] * n, type=pa.string()),
                    pa.array(pwm_arr.astype(np.float32), type=pa.float32()),
                    pa.array(cas_arr.astype(np.float32), type=pa.float32()),
                ],
                schema=schema,
            )
            writer.write_table(table_chunk)
    finally:
        writer.close()


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

    # Sort most-significant first. CAS p-value (index 10) is a monotone
    # transform of CAS score (index 9), so sorting by score descending is
    # biologically identical -- and avoids the lexical-vs-numeric string
    # issues with scientific notation and "<p" / ">p" prefix sentinels
    # that make a direct sort by the p-value column unreliable.
    if len(sorted_clusters_target_species_hits_list) > 0:
        sorted_clusters_target_species_hits_list = sorted(sorted_clusters_target_species_hits_list, key=itemgetter(9), reverse=True)

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
