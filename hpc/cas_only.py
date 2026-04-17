#!/usr/bin/env python3
"""Stage C: CAS-only computation harness for a chunk of transcripts.

Designed to be launched as a SLURM array task (see hpc/slurm/submit.sh).
Each invocation processes a slice of a species' transcript list, drives
tfbs_footprinter3 through its normal CLI path with:

    - ensemblrest() monkey-patched to read from a shared sqlite cache
      (populated by Stage B: hpc/prefetch.py)
    - pvalc=1.0 so every hit is emitted (we need the full empirical
      distribution, not just significant hits)
    - plotting disabled (--nofig) to strip ~seconds per transcript

After the run, parses the per-transcript `TFBSs_found.sortedclusters.csv`
output files and emits a single Parquet with (transcript_id, tf_name,
cas_score) per hit. These shards get concatenated by Stage D
(hpc/build_cas_pvalues.py) into per-species empirical p-value tables.

Usage:
    python -m hpc.cas_only \\
        --list hpc/transcript_lists/mus_musculus.txt \\
        --range 0:10 \\
        --cache hpc/cache/ensembl.sqlite \\
        --output hpc/cas_scores/mus_musculus/task_000.parquet

    # Benchmark mode: run one transcript, print timing, don't write
    python -m hpc.cas_only --list hpc/transcript_lists/mus_musculus.txt \\
        --range 0:1 --cache hpc/cache/ensembl.sqlite --benchmark
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
import tempfile
import time
from pathlib import Path

# CAS_OUTPUT_COLUMNS — output Parquet schema
# Keep as simple dataclass-less dict-list so pandas/pyarrow picks it up cleanly
# and Stage D can read without any HPC-local imports.
CAS_OUTPUT_COLUMNS = ["transcript_id", "tf_name", "cas_score"]

# Columns present in TFBSs_found.sortedclusters.csv (see README.md section 4.2)
# We only need `binding prot` and `combined affinity score`.
CSV_TF_COL = "binding prot"
CSV_CAS_COL = "combined affinity score"


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--list", type=Path, required=True, help="Transcript list file (one Ensembl ID per line).")
    p.add_argument("--range", type=str, required=True, help="Slice 'START:END' (0-indexed, end exclusive).")
    p.add_argument("--cache", type=Path, required=True, help="Sqlite cache written by prefetch.py (read-only).")
    p.add_argument("--output", type=Path, help="Destination Parquet file (required unless --benchmark).")
    p.add_argument("--promoter-before", type=int, default=900)
    p.add_argument("--promoter-after", type=int, default=100)
    p.add_argument("--pval", type=float, default=0.01, help="PWM p-value threshold (passed to tool).")
    p.add_argument("--top-x", type=int, default=10, help="--tx flag forwarded to tool (affects figure only).")
    p.add_argument("--benchmark", action="store_true",
                   help="Run the slice, print per-transcript timing, don't write Parquet or allow live REST.")
    p.add_argument("--allow-live-rest", action="store_true",
                   help="Let ensemblrest() fall back to live REST if cache miss. Default off.")
    return p.parse_args(argv)


def slice_transcripts(path: Path, rng: str) -> list[str]:
    start_s, end_s = rng.split(":", 1)
    start = int(start_s)
    end = int(end_s)
    ids = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    return ids[start:end]


def write_tool_csv(transcripts: list[str], out_csv: Path, args: argparse.Namespace) -> None:
    """Emit the CSV format tfbs_footprinter3 accepts via -t <file>.

    Column order per parse_transcript_ids() / file_to_datalist() in the
    monolith: transcript_id, tf_ids_file, promoter_before, promoter_after,
    top_x, pval, pvalc. Set pvalc=1.0 so every hit passes the CAS filter.
    """
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        for tid in transcripts:
            w.writerow([tid, "", args.promoter_before, args.promoter_after,
                        args.top_x, args.pval, 1.0])


def collect_cas_rows(results_dir: Path, transcripts: list[str], args: argparse.Namespace) -> list[tuple]:
    """Parse each transcript's TFBSs_found.sortedclusters.csv, emit hit rows."""
    # Tool names its per-transcript subdir as {tid}_(before_after)_{pval}
    suffix = f"({args.promoter_before}_{args.promoter_after})_{args.pval}"
    rows: list[tuple] = []
    missing = []
    for tid in transcripts:
        csv_path = results_dir / f"{tid}_{suffix}" / "TFBSs_found.sortedclusters.csv"
        if not csv_path.exists():
            missing.append(tid)
            continue
        with csv_path.open() as f:
            reader = csv.DictReader(f)
            for r in reader:
                try:
                    score = float(r[CSV_CAS_COL])
                except (KeyError, ValueError, TypeError):
                    continue
                rows.append((tid, r[CSV_TF_COL], score))
    if missing:
        logging.warning("no output for %d transcript(s): %s", len(missing), missing[:5])
    return rows


def write_parquet(rows: list[tuple], out_path: Path) -> None:
    import pandas as pd  # noqa: PLC0415  local import to keep --help fast
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows, columns=CAS_OUTPUT_COLUMNS)
    df.to_parquet(out_path, index=False)
    logging.info("wrote %d rows to %s", len(df), out_path)


def run_main_inprocess(tool_csv: Path, workdir: Path) -> None:
    """Drive tfbs_footprinter3.main() after chdir + curdir patch."""
    original_cwd = os.getcwd()
    os.chdir(workdir)
    try:
        from tfbs_footprinter3 import tfbs_footprinter3 as tff
        # curdir was captured at import-time (tfbs_footprinter3.py:64); re-bind.
        tff.curdir = str(workdir)
        saved_argv = sys.argv[:]
        sys.argv = ["tfbs_footprinter3", "-t", str(tool_csv), "-no"]
        try:
            tff.main()
        finally:
            sys.argv = saved_argv
    finally:
        os.chdir(original_cwd)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    if not args.benchmark and not args.output:
        sys.exit("--output is required unless --benchmark is set")

    transcripts = slice_transcripts(args.list, args.range)
    if not transcripts:
        logging.warning("empty slice: %s %s", args.list, args.range)
        return 0
    logging.info("processing %d transcripts (%s)", len(transcripts), args.range)

    # Apply the cache monkey-patch BEFORE importing tfbs_footprinter3 main
    from hpc.ensembl_cache import CachedEnsemblClient, patch_tfbs_footprinter3
    client = CachedEnsemblClient(args.cache, read_only=not args.allow_live_rest)
    patch_tfbs_footprinter3(client)

    with tempfile.TemporaryDirectory(prefix="tfbs_cas_only_") as tmp:
        workdir = Path(tmp)
        tool_csv = workdir / "transcripts.csv"
        write_tool_csv(transcripts, tool_csv, args)

        t0 = time.time()
        run_main_inprocess(tool_csv, workdir)
        elapsed = time.time() - t0
        logging.info("tool run complete in %.1fs (%.1fs/transcript avg)",
                     elapsed, elapsed / len(transcripts))

        results_dir = workdir / "tfbs_results"
        rows = collect_cas_rows(results_dir, transcripts, args)
        logging.info("collected %d CAS hit rows", len(rows))

        if args.benchmark:
            print(f"BENCHMARK  transcripts={len(transcripts)}  elapsed={elapsed:.1f}s  "
                  f"avg={elapsed/len(transcripts):.1f}s/tid  cas_rows={len(rows)}")
            return 0

        write_parquet(rows, args.output)
    client.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
