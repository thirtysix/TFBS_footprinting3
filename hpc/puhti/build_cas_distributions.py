#!/usr/bin/env python3
"""Aggregate tfbs_footprinter3's native per-transcript Parquet output into
the per-species threshold artifacts.

Reads: `<results_dir>/*/TFBSs_found.sortedclusters.parquet`, the layout
tfbs_footprinter3 produces with `-of parquet` or `-of slim-parquet`.

Writes per species:

  1. <out_dir>/<species>/<species>.CAS_thresholds.jaspar_2026.tsv.gz

     The CAS p-value lookup table. Consumed by the tool's
     calcCombinedAffinityPvalue (scoring.py). Grid: 0.01..1.0 step 0.01
     + decades 1e-3..1e-12 + a 0.1x subgrid.

  2. <out_dir>/<species>/<species>.tfs_thresholds.jaspar_2026.tsv.gz

     The per-PWM p-value lookup table. Consumed by the tool's PWM hit
     filter / p-value assignment at pipeline.py. Grid: 0.01..1.0 step
     0.01 + decades 1e-3..1e-9 (narrower than CAS; raw PWM scores
     don't resolve beyond that at typical sampling depth).

Both files use the same (tf_name, p_value, score) TSV schema and drop
duplicate-score rows below the statistical floor 1/N.

Usage:
    python hpc/puhti/build_cas_distributions.py \\
        --species acanthochromis_polyacanthus \\
        --results-dir runs/acanthochromis_polyacanthus/tfbs_results \\
        --out-dir hpc/artifacts
"""
from __future__ import annotations

import argparse
import gzip
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Column names in the tool's native parquet output (see output.py).
_TF_COL = "binding prot."
_CAS_COL = "combined affinity score"
_PWM_COL = "PWM score"
_ROUND_DECIMALS = 2


def _cas_target_p_values() -> list[float]:
    """CAS p-value target grid: coarse 0.01..1.0 + decades 1e-3..1e-12
    + a 0.1x subgrid between decades. ~200 values; the fine subgrid is
    worth it for CAS because the sum-of-contributions distribution
    supports deeper tail resolution than raw PWM scores do."""
    coarse = [x / 100 for x in range(1, 101)]  # 0.01..1.00
    decades = [10 ** (-x) for x in range(3, 13)]  # 1e-3..1e-12
    fine = [y * (z / 10) for y in decades for z in range(1, 11)]  # 0.1..1.0 x each decade
    return sorted(set(coarse + decades + fine), reverse=True)


def _pwm_target_p_values() -> list[float]:
    """PWM p-value target grid: coarse 0.01..1.0 + decades 1e-3..1e-9.
    107 values, matching the reference
    `{species}.tfs_thresholds.jaspar_2018.tsv.gz` files. PWM score
    distributions are narrower than CAS (no signal sum), so the
    statistical floor shows up ~1e-6 for 2M-hit samples; going
    deeper than 1e-9 wastes rows on uninformative saturated tails."""
    coarse = [x / 100 for x in range(1, 101)]  # 0.01..1.00
    decades = [10 ** (-x) for x in range(3, 10)]  # 1e-3..1e-9
    return sorted(set(coarse + decades), reverse=True)


def _load_values_per_tf(results_dir: Path, value_col: str) -> dict[str, np.ndarray]:
    """Concatenate `value_col` across every per-transcript parquet in
    `results_dir`, grouped by TF name.

    Returns a dict tf_name -> 1-D float64 ndarray of all that TF's
    values across every hit across every transcript.

    Stream-concat per-file reads of just the two columns we need so
    peak memory stays bounded by one parquet's rows (~200 MB) plus the
    per-TF accumulator for the value column (~4-8 GB for 1019 TFs x
    ~2M hits).
    """
    parquet_files = sorted(results_dir.rglob("TFBSs_found.sortedclusters.parquet"))
    if not parquet_files:
        raise FileNotFoundError(f"no parquet files under {results_dir}")
    logging.info("aggregating %d parquet files for column %r from %s",
                 len(parquet_files), value_col, results_dir)

    per_tf_chunks: dict[str, list[np.ndarray]] = {}
    for i, path in enumerate(parquet_files, 1):
        df = pd.read_parquet(path, columns=[_TF_COL, value_col])
        for tf_name, sub in df.groupby(_TF_COL, sort=False):
            per_tf_chunks.setdefault(tf_name, []).append(sub[value_col].to_numpy())
        if i % 10 == 0 or i == len(parquet_files):
            logging.info("  %d/%d parquets read", i, len(parquet_files))

    return {tf: np.concatenate(chunks) for tf, chunks in per_tf_chunks.items()}


def _empirical_survival_pvalues(scores: np.ndarray) -> list[tuple[float, float]]:
    """Return sorted-ascending (score, P(CAS>=score)) tuples.

    Uses the rounded-to-2-decimal convention matching
    tfbs_footprinter3.py:1073. Survival: P(X >= s) = (# observations
    >= s) / N.
    """
    if scores.size == 0:
        return []
    rounded = np.round(scores, _ROUND_DECIMALS)
    # Sorted unique scores + their counts
    sorted_scores = np.sort(rounded)
    unique_scores, first_idx = np.unique(sorted_scores, return_index=True)
    # For each unique score s, # observations >= s = N - first_idx[s]
    n = sorted_scores.size
    survival_counts = n - first_idx.astype(np.float64)
    p_values = survival_counts / n
    return [(float(s), float(p)) for s, p in zip(unique_scores, p_values, strict=True)]


def _scores_at_target_pvalues(
    score_pvalue_pairs: list[tuple[float, float]], targets: list[float]
) -> list[tuple[float, float]]:
    """Interpolate the score at each target p-value.

    score_pvalue_pairs is sorted ascending by score, so p-values are
    descending. For each target p, pick the smallest score whose
    empirical p <= target. If target is smaller than any observed p
    (tail too thin), emit the largest observed score.
    """
    if not score_pvalue_pairs:
        return []
    scores = np.array([s for s, _ in score_pvalue_pairs])
    pvalues = np.array([p for _, p in score_pvalue_pairs])
    # pvalues are descending in score order. For each target t, find the
    # first (smallest-score) index where p <= t.
    out = []
    for t in targets:
        mask = pvalues <= t
        if mask.any():
            i = int(np.argmax(mask))  # first True
            out.append((t, float(scores[i])))
        else:
            # Target tail is tighter than the distribution resolves; report
            # the extreme tail score so downstream tools see the right-most
            # available reference.
            out.append((t, float(scores[-1])))
    return out


def _build_threshold_tsv(
    per_tf_values: dict[str, np.ndarray],
    target_p_values: list[float],
    min_hits: int,
    tsv_path: Path,
) -> tuple[int, int, int]:
    """Build and gzip-write a threshold TSV from a pre-loaded per-TF
    value dict.

    Returns (tfs_kept, tfs_dropped, rows_written).
    """
    tsv_rows: list[tuple[str, float, float]] = []
    tfs_kept = tfs_dropped = 0
    for tf_name in sorted(per_tf_values.keys()):
        scores = per_tf_values[tf_name]
        if scores.size < min_hits:
            tfs_dropped += 1
            continue
        pairs = _empirical_survival_pvalues(scores)
        # Drop consecutive rows that re-emit the same score under a
        # tighter target p-value (below the statistical floor of 1/N).
        # Targets iterate descending in p, so keep the FIRST (loosest) row
        # for each distinct score.
        seen_scores: set[float] = set()
        for t, s in _scores_at_target_pvalues(pairs, target_p_values):
            if s in seen_scores:
                continue
            seen_scores.add(s)
            tsv_rows.append((tf_name, t, s))
        tfs_kept += 1

    with gzip.open(tsv_path, "wt") as f:
        f.write("tf_name\tp_value\tscore\n")
        for tf_name, p, s in tsv_rows:
            f.write(f"{tf_name}\t{p}\t{s}\n")
    logging.info("wrote %s (%d rows, %d TFs kept, %d dropped)",
                 tsv_path, len(tsv_rows), tfs_kept, tfs_dropped)
    return tfs_kept, tfs_dropped, len(tsv_rows)


def build_for_species(
    species: str,
    results_dir: Path,
    out_dir: Path,
    min_hits: int,
) -> dict:
    out_species_dir = out_dir / species
    out_species_dir.mkdir(parents=True, exist_ok=True)

    # CAS threshold table -- tool loads via data_loader.py:306 for the
    # combined-affinity-score p-value lookup.
    cas_per_tf = _load_values_per_tf(results_dir, _CAS_COL)
    cas_tsv = out_species_dir / f"{species}.CAS_thresholds.jaspar_2026.tsv.gz"
    cas_kept, cas_dropped, cas_rows = _build_threshold_tsv(
        cas_per_tf, _cas_target_p_values(), min_hits, cas_tsv,
    )
    del cas_per_tf  # free ~8 GB before the PWM pass

    # Per-PWM threshold table -- tool loads via data_loader.py (glob
    # `.tfs_thresholds.*.tsv.gz`) for the PWM-score p-value lookup at
    # each individual hit.
    pwm_per_tf = _load_values_per_tf(results_dir, _PWM_COL)
    pwm_tsv = out_species_dir / f"{species}.tfs_thresholds.jaspar_2026.tsv.gz"
    pwm_kept, pwm_dropped, pwm_rows = _build_threshold_tsv(
        pwm_per_tf, _pwm_target_p_values(), min_hits, pwm_tsv,
    )

    return {
        "species": species,
        "status": "ok",
        "min_hits": min_hits,
        "cas": {"tfs_kept": cas_kept, "tfs_dropped": cas_dropped, "rows": cas_rows, "tsv": str(cas_tsv)},
        "pwm": {"tfs_kept": pwm_kept, "tfs_dropped": pwm_dropped, "rows": pwm_rows, "tsv": str(pwm_tsv)},
    }


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species", required=True, help="Species slug (e.g. acanthochromis_polyacanthus).")
    p.add_argument("--results-dir", type=Path, required=True,
                   help="The tool's tfbs_results/ directory (contains per-transcript subdirs with parquet shards).")
    p.add_argument("--out-dir", type=Path, default=Path("hpc/artifacts"),
                   help="Where per-species artifacts are written (default: hpc/artifacts).")
    p.add_argument("--min-hits", type=int, default=100,
                   help="Drop TFs with fewer than this many total hits (default: 100).")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    summary = build_for_species(args.species, args.results_dir, args.out_dir, args.min_hits)
    logging.info("summary: %s", summary)
    return 0 if summary["status"] == "ok" else 1


if __name__ == "__main__":
    sys.exit(main())
