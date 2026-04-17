#!/usr/bin/env python3
"""Stage D: aggregate Parquet shards from Stage C into CAS_pvalues JSON per species.

For each species:
1. Concatenate all hpc/cas_scores/{species}/*.parquet shards into a single
   dataframe of (transcript_id, tf_name, cas_score) rows.
2. For each TF with at least MIN_HITS observations, build the empirical
   survival p-value function:
       p(s) = P(CAS >= s)
   computed on the unique rounded scores (2 decimals, matching
   tfbs_footprinter3.py:1073).
3. Emit {tf_name: [[score, p_value], ...]} sorted ascending by score —
   that's the exact shape consumed by calcCombinedAffinityPvalue's
   bisect_left lookup at tfbs_footprinter3.py:1001.
4. Write hpc/artifacts/{species}/CAS_pvalues.0.1.tf_ls.json.

Stage E (not yet scripted) injects each JSON into the per-species S3
tarball at s3://tfbssexperimentaldata/{species}.tar.gz. After re-upload,
tfbs_footprinter3 -update on user machines will download the new tarball
and species_specific_data() will pick the JSON up automatically with no
code changes (tfbs_footprinter3.py:755).

Usage:
    python hpc/build_cas_pvalues.py --species mus_musculus
    python hpc/build_cas_pvalues.py --species all
    python hpc/build_cas_pvalues.py --species mus_musculus --min-hits 50
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

DEFAULT_MIN_HITS = 100
ROUND_DECIMALS = 2


def empirical_pvalues(scores: list[float]) -> list[tuple[float, float]]:
    """Return sorted-ascending [(score, p_value)] tuples.

    p_value at each unique rounded score = P(CAS >= score) across the
    input distribution. The p-value at the minimum score is 1.0 by
    definition; p-value at the maximum is 1/N.
    """
    if not scores:
        return []
    import pandas as pd
    s = pd.Series(scores).round(ROUND_DECIMALS).sort_values().reset_index(drop=True)
    n = len(s)
    # Count-by-unique, ascending
    counts = s.value_counts().sort_index(ascending=True)
    # Survival: P(X >= score). Since scores are sorted, the running
    # "ranks from the top" is: N - cumulative_count_up_to_before_this_score.
    # Use cumcount via cumsum of counts shifted up by one.
    cumulative = counts.cumsum()
    # Number of observations strictly less than each score
    strictly_less = cumulative.shift(1).fillna(0)
    p_values = (n - strictly_less) / n
    # Build sorted (score, p) tuples
    return [(float(score), float(pv)) for score, pv in zip(p_values.index, p_values.values, strict=True)]


def build_one_species(species: str, shard_dir: Path, out_dir: Path, min_hits: int) -> dict:
    """Aggregate shards, compute p-values per TF, write JSON, return summary."""
    import pandas as pd
    shards = sorted(shard_dir.glob("*.parquet"))
    if not shards:
        logging.warning("%s: no shards under %s", species, shard_dir)
        return {"species": species, "status": "no_shards"}

    logging.info("%s: concatenating %d shards", species, len(shards))
    df = pd.concat([pd.read_parquet(p) for p in shards], ignore_index=True)
    if df.empty:
        return {"species": species, "status": "empty"}

    logging.info("%s: %d total hits across %d TFs", species, len(df), df["tf_name"].nunique())

    tf_pvalues: dict[str, list[list[float]]] = {}
    tf_counts: dict[str, int] = {}
    dropped: list[tuple[str, int]] = []
    for tf_name, grp in df.groupby("tf_name"):
        count = len(grp)
        tf_counts[tf_name] = count
        if count < min_hits:
            dropped.append((tf_name, count))
            continue
        # Store as list of 2-element lists so the JSON matches existing format
        tf_pvalues[tf_name] = [[score, pv] for score, pv in empirical_pvalues(grp["cas_score"].tolist())]

    out_path = out_dir / species / "CAS_pvalues.0.1.tf_ls.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(tf_pvalues))
    logging.info("%s: wrote %d TFs (%d dropped below min_hits=%d) to %s",
                 species, len(tf_pvalues), len(dropped), min_hits, out_path)

    return {
        "species": species,
        "status": "ok",
        "total_hits": int(len(df)),
        "tfs_kept": len(tf_pvalues),
        "tfs_dropped": len(dropped),
        "min_hits": min_hits,
        "output": str(out_path),
    }


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species", type=str, required=True,
                   help="Species slug (e.g. 'mus_musculus') or 'all' to iterate every species under --shard-root.")
    p.add_argument("--shard-root", type=Path, default=Path("hpc/cas_scores"))
    p.add_argument("--out-root", type=Path, default=Path("hpc/artifacts"))
    p.add_argument("--min-hits", type=int, default=DEFAULT_MIN_HITS)
    p.add_argument("--manifest", type=Path, default=Path("hpc/artifacts/build_manifest.json"))
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    if args.species == "all":
        species_list = sorted(p.name for p in args.shard_root.iterdir() if p.is_dir())
    else:
        species_list = [args.species]

    manifest = []
    failed = 0
    for sp in species_list:
        summary = build_one_species(sp, args.shard_root / sp, args.out_root, args.min_hits)
        manifest.append(summary)
        if summary.get("status") != "ok":
            failed += 1

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2))
    logging.info("manifest: %s (%d ok, %d failed)",
                 args.manifest, len(species_list) - failed, failed)
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
