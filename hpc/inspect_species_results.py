#!/usr/bin/env python3
"""Systematic sanity check on a species' Puhti CAS-campaign output.

Given a directory containing per-transcript parquet shards (e.g.
`hpc/cas_scores/acanthochromis_polyacanthus/`, populated by the
`rsync … tasks/` step of the Puhti pipeline), reports:

  - Completion rate vs the expected transcript count
  - Per-transcript row/TF/file-size summary (sortable)
  - Nonzero-rate for each contextual signal (GERP, CAGE, CpG, ...)
  - Cross-transcript CAS distribution stats
  - Per-TF hit counts (top + bottom) to flag TFs with anomalously few/many hits
  - Outliers that warrant a closer look

Schema assumed is the one emitted by `output.py:_PARQUET_HEADER`.

Usage:
    python hpc/inspect_species_results.py \\
        --results-dir hpc/cas_scores/acanthochromis_polyacanthus \\
        [--expected-transcripts 100]
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_CONTEXT_COLS = [
    "species weights sum",      # GERP
    "cage weights sum",
    "cpg weight",
    "atac weights sum",
    "eqtls weights sum",
    "metacluster weights sum",
    "corr. weight sum",
]
_SCORE_COLS = ["PWM score", "combined affinity score"]
_TF_COL = "binding prot."


def _find_parquets(results_dir: Path) -> list[Path]:
    return sorted(results_dir.rglob("TFBSs_found.sortedclusters.parquet"))


def _summarize_transcript(path: Path) -> dict:
    """Compact per-transcript summary. Opens the parquet once."""
    df = pd.read_parquet(path)
    # Extract transcript ID from the "<ENSAPOT…>(5000_5000)_1.0" parent dir
    transcript_id = path.parent.name.split("(")[0]

    rec: dict = {
        "transcript_id": transcript_id,
        "path": path,
        "rows": len(df),
        "tfs": int(df[_TF_COL].nunique()),
        "file_mb": round(path.stat().st_size / (1024 * 1024), 1),
        "cas_mean": float(df["combined affinity score"].mean()),
        "cas_max": float(df["combined affinity score"].max()),
        "pwm_max": float(df["PWM score"].max()),
    }
    for col in _CONTEXT_COLS:
        nz = int((df[col] != 0).sum())
        rec[f"nz_{col.replace(' ', '_')}_pct"] = round(100 * nz / len(df), 1) if len(df) else 0.0
    return rec


def _render_transcript_table(rows: list[dict]) -> str:
    if not rows:
        return "  (no transcripts)"
    df = pd.DataFrame(rows).drop(columns=["path"])
    # Trim to the most important columns for display
    cols = ["transcript_id", "rows", "tfs", "file_mb",
            "cas_mean", "cas_max", "pwm_max",
            "nz_species_weights_sum_pct",  # GERP
            "nz_cage_weights_sum_pct",
            "nz_cpg_weight_pct"]
    df = df[cols]
    return df.to_string(index=False)


def _aggregate_cas(paths: list[Path], tf_col: str = _TF_COL, cas_col: str = "combined affinity score") -> pd.DataFrame:
    """Stream-concat every transcript's CAS scores, grouped by TF.

    Peak memory is bounded by one transcript at a time; the output is a
    small DataFrame of per-TF counts + quantile summaries.
    """
    per_tf_stats: dict[str, dict] = {}
    for i, p in enumerate(paths, 1):
        df = pd.read_parquet(p, columns=[tf_col, cas_col])
        for tf, sub in df.groupby(tf_col, sort=False):
            s = sub[cas_col].to_numpy()
            rec = per_tf_stats.setdefault(tf, {"count": 0, "sum": 0.0, "sq_sum": 0.0, "min": np.inf, "max": -np.inf})
            rec["count"] += s.size
            rec["sum"] += float(s.sum())
            rec["sq_sum"] += float((s * s).sum())
            rec["min"] = min(rec["min"], float(s.min()))
            rec["max"] = max(rec["max"], float(s.max()))
        if i % 10 == 0 or i == len(paths):
            logging.info("  aggregated %d/%d parquets", i, len(paths))

    rows = []
    for tf, rec in per_tf_stats.items():
        n = rec["count"]
        mean = rec["sum"] / n if n else 0
        var = max(0.0, rec["sq_sum"] / n - mean * mean) if n else 0
        rows.append({
            "tf_name": tf,
            "hits": n,
            "mean": round(mean, 3),
            "std": round(var ** 0.5, 3),
            "min": round(rec["min"], 2),
            "max": round(rec["max"], 2),
        })
    return pd.DataFrame(rows).sort_values("hits", ascending=False).reset_index(drop=True)


def inspect(results_dir: Path, expected_transcripts: int | None) -> int:
    parquets = _find_parquets(results_dir)

    print(f"== Results dir: {results_dir} ==")
    print(f"Parquets found: {len(parquets)}"
          + (f" / expected {expected_transcripts}" if expected_transcripts else ""))
    if not parquets:
        print("Nothing to inspect.")
        return 1

    completion_pct = (100 * len(parquets) / expected_transcripts) if expected_transcripts else None
    if completion_pct is not None:
        print(f"Completion: {completion_pct:.1f}%")
        if len(parquets) < expected_transcripts:
            print(f"  ⚠ {expected_transcripts - len(parquets)} transcripts missing")
    print()

    # Per-transcript pass
    print("== Per-transcript summary ==")
    logging.info("scanning per-transcript...")
    records = [_summarize_transcript(p) for p in parquets]
    print(_render_transcript_table(records))
    print()

    # Quick health flags at transcript level
    print("== Transcript-level health flags ==")
    small = [r for r in records if r["tfs"] < 1000]
    zero_gerp = [r for r in records if r["nz_species_weights_sum_pct"] == 0.0]
    zero_cpg = [r for r in records if r["nz_cpg_weight_pct"] == 0.0]
    if small:
        print(f"  ⚠ {len(small)} transcripts with <1000 unique TFs (expected ~1019):"
              + ", ".join(r["transcript_id"] for r in small[:5]) + (" …" if len(small) > 5 else ""))
    if zero_gerp:
        print(f"  ⚠ {len(zero_gerp)} transcripts with 0% nonzero GERP "
              f"(likely regions of genome with no conservation data): "
              + ", ".join(r["transcript_id"] for r in zero_gerp[:5]) + (" …" if len(zero_gerp) > 5 else ""))
    if zero_cpg:
        print(f"  ⚠ {len(zero_cpg)} transcripts with 0% CpG weight")
    if not (small or zero_gerp or zero_cpg):
        print("  ✓ no transcript-level anomalies")
    print()

    # Cross-transcript CAS distribution
    print("== Cross-transcript CAS aggregation (per-TF) ==")
    per_tf_df = _aggregate_cas(parquets)
    print(f"  Total hits:    {per_tf_df['hits'].sum():,}")
    print(f"  Unique TFs:    {len(per_tf_df)}")
    print(f"  Median hits/TF: {int(per_tf_df['hits'].median()):,}")
    print(f"  Min hits/TF:   {int(per_tf_df['hits'].min()):,} ({per_tf_df.iloc[-1]['tf_name']})")
    print(f"  Max hits/TF:   {int(per_tf_df['hits'].max()):,} ({per_tf_df.iloc[0]['tf_name']})")
    print()
    print("  Top 5 TFs by hit count:")
    print(per_tf_df.head(5).to_string(index=False))
    print()
    print("  Bottom 5 TFs by hit count:")
    print(per_tf_df.tail(5).to_string(index=False))
    print()

    # CAS-score health
    global_max_cas = per_tf_df["max"].max()
    global_min_cas = per_tf_df["min"].min()
    print("== CAS score range ==")
    print(f"  Global max:    {global_max_cas:.2f}  (TF: {per_tf_df.loc[per_tf_df['max'].idxmax(), 'tf_name']})")
    print(f"  Global min:    {global_min_cas:.2f}")

    under_min_hits = per_tf_df[per_tf_df["hits"] < 100]
    if not under_min_hits.empty:
        print()
        print(f"  ⚠ {len(under_min_hits)} TFs with <100 hits — below the typical"
              f" empirical-p-value floor. build_cas_distributions defaults to"
              f" dropping these.")

    print()
    print("== Done ==")
    return 0


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-dir", type=Path, required=True,
                   help="Directory containing per-transcript parquet shards "
                        "(e.g. hpc/cas_scores/<species>/).")
    p.add_argument("--expected-transcripts", type=int, default=None,
                   help="Expected transcript count (e.g. 100 for the pilot). "
                        "Drives the completion-percent line.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    return inspect(args.results_dir, args.expected_transcripts)


if __name__ == "__main__":
    sys.exit(main())
