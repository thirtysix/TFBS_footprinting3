#!/usr/bin/env python3
"""Aggregate tfbs_footprinter3's native per-transcript Parquet output into
the per-species CAS distribution artifacts.

Reads: `<results_dir>/*/TFBSs_found.sortedclusters.parquet`, the layout
tfbs_footprinter3 produces when driven with `-of parquet`.

Writes two artifacts per species:

  1. <out_dir>/<species>/CAS_pvalues.0.1.tf_ls.json

     The runtime lookup file `species_specific_data()` loads from the S3
     tarball at `tfbs_footprinter3/tfbs_footprinter3.py:755`. For each
     TF, a sorted list of (score, empirical-p-value) pairs.

  2. <out_dir>/<species>/<species>.CAS_thresholds.jaspar_2018.tsv.gz

     The human-readable threshold table mirroring the published
     `homo_sapiens.CAS_thresholds.jaspar_2018.tsv.gz`. For each TF, one
     row per target p-value in the same grid the reference pipeline used:
     [0.01..1.0 in 0.01 steps] + [1e-3..1e-12 decades] + a 0.1x grid
     between them.

Usage:
    python hpc/puhti/build_cas_distributions.py \\
        --species acanthochromis_polyacanthus \\
        --results-dir runs/acanthochromis_polyacanthus/tfbs_results \\
        --out-dir hpc/artifacts
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Column names in the tool's native parquet output (see output.py).
_TF_COL = "binding prot."
_CAS_COL = "combined affinity score"
_ROUND_DECIMALS = 2

# NOTE on the "0.1" in the output filename: this is a historical artifact.
# Early tool releases truncated the JSON to p <= 0.1 to save space, then
# bisected the remaining sorted list for lookup. The tool still loads the
# same filename (`tfbs_footprinter3/data_loader.py:306`) and still uses
# bisect_left on the sorted score list, but there's no reason to keep the
# truncation: emitting the full distribution lets the tool report p-values
# for non-significant hits too (p between 0.1 and 1.0), which
# `pwm_results_summaries/*.CAS_pvalue_score.tsv` was designed to provide.
# We keep the literal filename for backward compatibility with the loader.


def _target_p_values() -> list[float]:
    """Reproduce the grid used in the reference
    results_extraction.005.py for the human CAS run."""
    coarse = [x / 100 for x in range(1, 101)]  # 0.01..1.00
    decades = [10 ** (-x) for x in range(3, 13)]  # 1e-3..1e-12
    fine = [y * (z / 10) for y in decades for z in range(1, 11)]  # 0.1..1.0 x each decade
    return sorted(set(coarse + decades + fine), reverse=True)


def _load_cas_scores_per_tf(results_dir: Path) -> dict[str, np.ndarray]:
    """Concatenate CAS scores across every per-transcript parquet in
    `results_dir`, grouped by TF.

    Returns a dict tf_name -> 1-D float64 ndarray of all that TF's CAS
    scores across every hit across every transcript.
    """
    parquet_files = sorted(results_dir.rglob("TFBSs_found.sortedclusters.parquet"))
    if not parquet_files:
        raise FileNotFoundError(f"no parquet files under {results_dir}")
    logging.info("aggregating %d parquet files from %s", len(parquet_files), results_dir)

    # Stream-concat per-file reads of just the two columns we need so peak
    # memory stays proportional to one parquet file (~200 MB) instead of
    # the whole 17 GB species total.
    per_tf_chunks: dict[str, list[np.ndarray]] = {}
    for i, path in enumerate(parquet_files, 1):
        df = pd.read_parquet(path, columns=[_TF_COL, _CAS_COL])
        for tf_name, sub in df.groupby(_TF_COL, sort=False):
            per_tf_chunks.setdefault(tf_name, []).append(sub[_CAS_COL].to_numpy())
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


def build_for_species(
    species: str,
    results_dir: Path,
    out_dir: Path,
    min_hits: int,
) -> dict:
    out_species_dir = out_dir / species
    out_species_dir.mkdir(parents=True, exist_ok=True)

    per_tf_scores = _load_cas_scores_per_tf(results_dir)

    tf_pvalues_json: dict[str, list[list[float]]] = {}
    tsv_rows: list[tuple[str, float, float]] = []
    tfs_kept = tfs_dropped = 0
    targets = _target_p_values()

    for tf_name in sorted(per_tf_scores.keys()):
        scores = per_tf_scores[tf_name]
        if scores.size < min_hits:
            tfs_dropped += 1
            continue
        pairs = _empirical_survival_pvalues(scores)
        tf_pvalues_json[tf_name] = [[s, p] for s, p in pairs]
        for t, s in _scores_at_target_pvalues(pairs, targets):
            tsv_rows.append((tf_name, t, s))
        tfs_kept += 1

    json_path = out_species_dir / "CAS_pvalues.0.1.tf_ls.json"
    json_path.write_text(json.dumps(tf_pvalues_json))
    logging.info("wrote %s (%d TFs)", json_path, tfs_kept)

    tsv_path = out_species_dir / f"{species}.CAS_thresholds.jaspar_2018.tsv.gz"
    with gzip.open(tsv_path, "wt") as f:
        f.write("tf_name\tp_value\tscore\n")
        for tf_name, p, s in tsv_rows:
            f.write(f"{tf_name}\t{p}\t{s}\n")
    logging.info("wrote %s (%d rows)", tsv_path, len(tsv_rows))

    return {
        "species": species,
        "status": "ok",
        "tfs_kept": tfs_kept,
        "tfs_dropped": tfs_dropped,
        "min_hits": min_hits,
        "json": str(json_path),
        "tsv": str(tsv_path),
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
