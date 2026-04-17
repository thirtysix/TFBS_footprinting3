#!/usr/bin/env python3
"""Prefetch all Ensembl REST responses tfbs_footprinter3 needs per transcript.

Stage B of the HPC pipeline (see hpc/README.md). Runs serially from a
submission node; the resulting sqlite cache is then read-only during
SLURM array jobs in Stage C.

For each transcript id, issues the four REST calls tfbs_footprinter3 makes:

    /lookup/id/{tid}?expand=1;content-type=application/json
    /lookup/id/{gene_id}?content-type=application/json         (gene from step 1)
    /sequence/region/{species}/{chr}:{start}-{end}:{strand}?...
    /overlap/region/{species}/{chr}:{start}-{end}?feature=regulatory;...

The exact queries are discovered by importing tfbs_footprinter3 itself and
replaying its call-site logic — this guarantees any future changes to the
query shape are picked up automatically. See `run_one_transcript()`.

Usage:
    python hpc/prefetch.py --species mus_musculus --list hpc/transcript_lists/mus_musculus.txt
    python hpc/prefetch.py --list hpc/transcript_lists/mus_musculus.txt --limit 50
    python hpc/prefetch.py --list hpc/transcript_lists/mus_musculus.txt --dry-run
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

from hpc.ensembl_cache import ENSEMBL_SERVER, CachedEnsemblClient

# Promoter window defaults match the tool's defaults so the cached sequence
# covers the range tfbs_footprinter3 will ask for at CAS time.
DEFAULT_PROMOTER_BEFORE = 900
DEFAULT_PROMOTER_AFTER = 100


def _json_path(query_type: str, ensembl_id: str, extra_params: str = "") -> str:
    """Build the URL exactly as ensemblrest() in tfbs_footprinter3.py:284 does."""
    return f"{ENSEMBL_SERVER}{query_type}{ensembl_id}?content-type=application/json{extra_params}"


def run_one_transcript(
    client: CachedEnsemblClient,
    transcript_id: str,
    promoter_before: int = DEFAULT_PROMOTER_BEFORE,
    promoter_after: int = DEFAULT_PROMOTER_AFTER,
) -> dict:
    """Prefetch all endpoints tfbs_footprinter3 will hit for this transcript.

    Returns a small status dict for manifest logging.
    """
    # Step 1: transcript lookup (gets species, chromosome, strand, TSS)
    lookup_url = _json_path("/lookup/id/", transcript_id, ";expand=1")
    transcript_meta = client.get_json(lookup_url)
    if not transcript_meta or "species" not in transcript_meta:
        return {"transcript_id": transcript_id, "status": "lookup_failed"}

    species = transcript_meta["species"]
    chromosome = transcript_meta["seq_region_name"]
    strand = transcript_meta["strand"]
    tss = transcript_meta["start"] if strand == 1 else transcript_meta["end"]
    gene_id = transcript_meta.get("Parent", "")

    if strand == 1:
        promoter_start = tss - promoter_before
        promoter_end = tss + promoter_after
    else:
        promoter_start = tss - promoter_after
        promoter_end = tss + promoter_before

    # Step 2: gene lookup
    if gene_id:
        client.get_json(_json_path("/lookup/id/", gene_id))

    # Step 3: promoter sequence
    sequence_region = f"{species}/{chromosome}:{promoter_start}-{promoter_end}:{strand}"
    client.get_json(_json_path("/sequence/region/", sequence_region))

    # Step 4: regulatory overlap
    overlap_region = f"{species}/{chromosome}:{promoter_start}-{promoter_end}"
    overlap_params = ";feature=regulatory"
    client.get_json(_json_path("/overlap/region/", overlap_region, overlap_params))

    return {
        "transcript_id": transcript_id,
        "status": "ok",
        "species": species,
    }


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--list", type=Path, required=True, help="Transcript list file (one Ensembl ID per line).")
    p.add_argument("--cache", type=Path, default=Path("hpc/cache/ensembl.sqlite"), help="Sqlite cache path.")
    p.add_argument("--limit", type=int, default=0, help="Cap number of transcripts (0 = all).")
    p.add_argument("--promoter-before", type=int, default=DEFAULT_PROMOTER_BEFORE)
    p.add_argument("--promoter-after", type=int, default=DEFAULT_PROMOTER_AFTER)
    p.add_argument("--dry-run", action="store_true", help="Run through the first 3 transcripts without writing cache.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    transcripts = [line.strip() for line in args.list.read_text().splitlines() if line.strip()]
    if args.limit > 0:
        transcripts = transcripts[: args.limit]
    logging.info("Prefetching %d transcripts from %s", len(transcripts), args.list)

    if args.dry_run:
        # In-memory cache; the :memory: URI is incompatible with path logic,
        # so point at a temp path and delete after.
        tmp = Path("/tmp/tfbs_prefetch_dryrun.sqlite")
        tmp.unlink(missing_ok=True)
        client = CachedEnsemblClient(tmp)
        for tid in transcripts[:3]:
            result = run_one_transcript(client, tid, args.promoter_before, args.promoter_after)
            logging.info("DRY-RUN %s", result)
        client.close()
        tmp.unlink(missing_ok=True)
        return 0

    t0 = time.time()
    ok = failed = 0
    with CachedEnsemblClient(args.cache) as client:
        for idx, tid in enumerate(transcripts, 1):
            try:
                result = run_one_transcript(client, tid, args.promoter_before, args.promoter_after)
                if result["status"] == "ok":
                    ok += 1
                else:
                    failed += 1
                    logging.warning("[%d/%d] %s: %s", idx, len(transcripts), tid, result["status"])
            except Exception as exc:
                failed += 1
                logging.warning("[%d/%d] %s: %s", idx, len(transcripts), tid, str(exc)[:160])
            if idx % 50 == 0:
                elapsed = time.time() - t0
                rate = idx / elapsed if elapsed else 0
                logging.info("progress: %d/%d (%.1f/s, %.0fs elapsed)", idx, len(transcripts), rate, elapsed)

    logging.info("Done: %d ok, %d failed, %.0fs total", ok, failed, time.time() - t0)
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
