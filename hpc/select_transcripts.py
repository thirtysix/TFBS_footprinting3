#!/usr/bin/env python3
"""Select ~N canonical protein-coding transcript IDs per Ensembl species.

Outputs one file per species at hpc/transcript_lists/{species}.txt, one
Ensembl transcript ID per line, consumable by tfbs_footprinter3 -t.

Uses Ensembl BioMart bulk TSV (one HTTP POST per species) rather than REST,
so 124 species = 124 requests total. Run this once on a workstation; the
resulting lists are checked into git-ignored hpc/transcript_lists/.

Usage:
    python hpc/select_transcripts.py                       # all species
    python hpc/select_transcripts.py --species mus_musculus danio_rerio
    python hpc/select_transcripts.py --max 5000 --out-dir hpc/transcript_lists
    python hpc/select_transcripts.py --dry-run             # fetch one species, don't write
"""
from __future__ import annotations

import argparse
import io
import json
import logging
import random
import sys
import time
import urllib.request
from pathlib import Path

ENSEMBL_REST = "https://rest.ensembl.org"
BIOMART_URL = "https://www.ensembl.org/biomart/martservice"

BIOMART_QUERY = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "{dataset}" interface = "default" >
        <Filter name = "biotype" value = "protein_coding"/>
        <Filter name = "transcript_is_canonical" excluded = "0"/>
        <Attribute name = "ensembl_transcript_id" />
    </Dataset>
</Query>
"""


def get_species_list() -> list[str]:
    """Fetch the canonical Ensembl species slugs (underscore form, e.g. 'mus_musculus')."""
    req = urllib.request.Request(
        f"{ENSEMBL_REST}/info/species?division=EnsemblVertebrates",
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        data = json.load(resp)
    return sorted(s["name"] for s in data["species"])


def biomart_dataset_name(species: str) -> str:
    """Convert 'mus_musculus' -> 'mmusculus_gene_ensembl' (BioMart convention).

    The convention is <first-letter-of-genus><full-species>_gene_ensembl.
    Examples: homo_sapiens -> hsapiens_gene_ensembl, danio_rerio -> drerio_gene_ensembl.
    """
    parts = species.split("_")
    if len(parts) < 2:
        raise ValueError(f"Unrecognized species slug: {species!r}")
    genus, *rest = parts
    return f"{genus[0]}{''.join(rest)}_gene_ensembl"


def fetch_transcripts(species: str, timeout: int = 300, max_attempts: int = 4) -> list[str]:
    """POST the BioMart query and return a list of transcript IDs (all canonical protein-coding).

    BioMart occasionally returns 'Could not connect to mysql database' errors;
    retry with exponential backoff.
    """
    dataset = biomart_dataset_name(species)
    body = f"query={BIOMART_QUERY.format(dataset=dataset)}".encode()
    last_err: Exception | None = None
    for attempt in range(max_attempts):
        try:
            req = urllib.request.Request(BIOMART_URL, data=body, method="POST")
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                payload = resp.read().decode("utf-8", errors="replace")
            if payload.startswith("Query ERROR") or "<h3>Error</h3>" in payload:
                raise RuntimeError(f"BioMart error: {payload[:200]}")
            ids = [line.strip() for line in io.StringIO(payload) if line.strip()]
            return ids
        except Exception as exc:
            last_err = exc
            if attempt < max_attempts - 1:
                delay = 5 * (2 ** attempt)  # 5s, 10s, 20s
                logging.info("Retry %s in %ds: %s", species, delay, str(exc)[:120])
                time.sleep(delay)
    raise RuntimeError(f"BioMart error for {species} ({dataset}) after {max_attempts} attempts: {last_err}")


def sample_ids(ids: list[str], max_count: int, seed: int) -> list[str]:
    """Return up to max_count IDs; if smaller, return all. Deterministic via seed."""
    if len(ids) <= max_count:
        return sorted(ids)
    rng = random.Random(seed)
    return sorted(rng.sample(ids, max_count))


def write_list(ids: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(ids) + "\n")


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species", nargs="*", help="Subset of species slugs to fetch (default: all).")
    p.add_argument("--max", type=int, default=5000, help="Cap transcripts per species (default: 5000).")
    p.add_argument("--out-dir", type=Path, default=Path("hpc/transcript_lists"))
    p.add_argument("--manifest", type=Path, default=Path("hpc/transcript_lists/manifest.json"))
    p.add_argument("--seed", type=int, default=20260417, help="RNG seed for subsampling (default: today's date).")
    p.add_argument("--sleep", type=float, default=1.0, help="Seconds to sleep between species to be polite to BioMart.")
    p.add_argument("--dry-run", action="store_true", help="Fetch one species (mus_musculus) and print, don't write files.")
    p.add_argument("--resume", action="store_true", help="Skip species whose output file already exists.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    if args.dry_run:
        species = "mus_musculus"
        logging.info("DRY RUN: fetching %s", species)
        ids = fetch_transcripts(species)
        sample = sample_ids(ids, args.max, args.seed)
        logging.info("%s: %d total, %d sampled, first 5: %s", species, len(ids), len(sample), sample[:5])
        return 0

    species_list = args.species if args.species else get_species_list()
    logging.info("Targeting %d species", len(species_list))

    manifest = {}
    if args.resume and args.manifest.exists():
        manifest = json.loads(args.manifest.read_text())

    for idx, species in enumerate(species_list, 1):
        out_path = args.out_dir / f"{species}.txt"
        if args.resume and out_path.exists():
            logging.info("[%d/%d] %s: skipping (already present)", idx, len(species_list), species)
            continue
        try:
            ids = fetch_transcripts(species)
        except Exception as exc:
            logging.warning("[%d/%d] %s: FAILED (%s)", idx, len(species_list), species, exc)
            manifest[species] = {"status": "failed", "error": str(exc)}
            continue
        sample = sample_ids(ids, args.max, args.seed)
        write_list(sample, out_path)
        manifest[species] = {
            "status": "ok",
            "total_canonical": len(ids),
            "sampled": len(sample),
            "reduced": len(sample) < args.max,
        }
        logging.info(
            "[%d/%d] %s: %d total -> %d sampled%s",
            idx, len(species_list), species, len(ids), len(sample),
            " (REDUCED)" if len(sample) < args.max else "",
        )
        time.sleep(args.sleep)

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    logging.info("Wrote manifest: %s", args.manifest)
    failed = [s for s, v in manifest.items() if v.get("status") != "ok"]
    if failed:
        logging.warning("%d species failed: %s", len(failed), failed)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
