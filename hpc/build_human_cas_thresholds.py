#!/usr/bin/env python3
"""Concatenate the per-TF human CAS threshold TSVs from the research archive
into the single gzipped artifact the tool will load at runtime.

Reads:  <summaries_dir>/*.CAS_pvalue_score.tsv
Writes: <out_dir>/homo_sapiens/homo_sapiens.CAS_thresholds.jaspar_2024.tsv.gz

Each input TSV already carries `tf_name\tp_value\tscore` rows across the
same p-value grid that hpc/puhti/build_cas_distributions.py emits for
non-human species, so this script is a straight concatenation: strip the
header from files 2..N, gzip the result.

Usage:
    python hpc/build_human_cas_thresholds.py \\
        --summaries-dir /path/to/pwm_results_summaries \\
        --out-dir hpc/artifacts
"""
from __future__ import annotations

import argparse
import gzip
import logging
import sys
from pathlib import Path


def build(summaries_dir: Path, out_dir: Path) -> dict:
    tsv_files = sorted(summaries_dir.glob("*.CAS_pvalue_score.tsv"))
    if not tsv_files:
        raise FileNotFoundError(f"no *.CAS_pvalue_score.tsv under {summaries_dir}")
    logging.info("found %d per-TF TSVs", len(tsv_files))

    out_species_dir = out_dir / "homo_sapiens"
    out_species_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_species_dir / "homo_sapiens.CAS_thresholds.jaspar_2024.tsv.gz"

    rows_written = 0
    with gzip.open(out_path, "wt") as out:
        out.write("tf_name\tp_value\tscore\n")
        for i, path in enumerate(tsv_files, 1):
            with open(path) as f:
                header = f.readline()
                if header.strip() != "tf_name\tp_value\tscore":
                    raise ValueError(f"unexpected header in {path}: {header!r}")
                for line in f:
                    out.write(line)
                    rows_written += 1
            if i % 100 == 0 or i == len(tsv_files):
                logging.info("  %d/%d TFs merged", i, len(tsv_files))
    logging.info("wrote %s (%d rows across %d TFs)", out_path, rows_written, len(tsv_files))
    return {"out": str(out_path), "tfs": len(tsv_files), "rows": rows_written}


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summaries-dir", type=Path, required=True,
                   help="Directory containing per-TF *.CAS_pvalue_score.tsv files.")
    p.add_argument("--out-dir", type=Path, default=Path("hpc/artifacts"),
                   help="Where to write the aggregated artifact (default: hpc/artifacts).")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    build(args.summaries_dir, args.out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
