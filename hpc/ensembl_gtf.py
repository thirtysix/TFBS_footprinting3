#!/usr/bin/env python3
"""Download and parse Ensembl release-113 GTF files, pick candidate transcripts.

Used by Stage A of the non-human CAS-distribution campaign:

  1. For each target species (slug in Ensembl's lowercase convention),
     download the release-113 GTF once and cache it under
     hpc/ensembl_gtfs/<species>.gtf.gz.
  2. Parse the GTF for protein-coding transcripts that carry
     'tag "Ensembl_canonical"'. Record (transcript_id, chromosome).
  3. Sample N transcripts (default 100) spread across up to K distinct
     chromosomes (default 10). The emitted order is shuffled -- the
     downstream pre-sort in cli.py clusters them again by chromosome
     before the main loop runs, so transcript-file order only affects
     reproducibility, not wall time.

Why GTF instead of Ensembl BioMart (used by the earlier
select_transcripts.py): BioMart throws transient 'Could not connect to
mysql database' errors that force retry loops and slow the full 316-
species fetch. GTF is a single FTP download per species, no retry logic
needed, and pins us to a specific release (important because the tool's
experimental-data bucket is built against Ensembl 113).

Usage:
    # Download + select for one species
    python hpc/ensembl_gtf.py --species acanthochromis_polyacanthus \\
        --n 100 --max-chromosomes 10

    # Batch mode: all 316 species
    python hpc/ensembl_gtf.py --species-file /tmp/species_316.txt --n 100
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import random
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

ENSEMBL_RELEASE = 113
ENSEMBL_FTP_BASE = f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/gtf"

# GTF attribute extractors. The attribute column is a semicolon-separated
# string of key "value"; pairs. Using specific regex (faster than a full
# parse) for just the fields we need.
_ATTR_BIOTYPE = re.compile(r'biotype "([^"]+)"')
_ATTR_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')
_ATTR_TAG_CANONICAL = re.compile(r'tag "Ensembl_canonical"')


def _ensembl_gtf_url(species: str) -> str:
    """Find the GTF URL for a species by listing the Ensembl FTP directory.

    Ensembl's filename embeds the assembly name (e.g.
    Acanthochromis_polyacanthus.ASM210954v1.113.gtf.gz), which varies by
    species. We hit the directory index and regex out the expected file.
    """
    dir_url = f"{ENSEMBL_FTP_BASE}/{species}/"
    req = urllib.request.Request(dir_url, headers={"User-Agent": "tfbs_footprinter3-hpc/0.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        html = resp.read().decode("utf-8", errors="replace")
    species_cap = species[0].upper() + species[1:]
    # Exclude .abinitio, .chr, .chr_patch_hapl_scaff variants -- we want the plain species.assembly.113.gtf.gz
    pattern = re.compile(rf'"({re.escape(species_cap)}\.[A-Za-z0-9_.\-]+\.{ENSEMBL_RELEASE}\.gtf\.gz)"')
    matches = [m.group(1) for m in pattern.finditer(html)]
    # Prefer the shortest filename (which is the plain .113.gtf.gz, not .chr_patch_hapl_scaff.113.gtf.gz etc.)
    if not matches:
        raise RuntimeError(f"no Ensembl {ENSEMBL_RELEASE} GTF found for {species!r} at {dir_url}")
    matches.sort(key=len)
    return dir_url + matches[0]


def download_gtf(species: str, cache_dir: Path) -> Path:
    """Download the Ensembl 113 GTF for `species` (cached)."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_path = cache_dir / f"{species}.gtf.gz"
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path
    url = _ensembl_gtf_url(species)
    logging.info("downloading %s -> %s", url, out_path)
    with urllib.request.urlopen(url, timeout=600) as resp, open(out_path, "wb") as f:
        while True:
            chunk = resp.read(1 << 16)
            if not chunk:
                break
            f.write(chunk)
    return out_path


def iter_canonical_protein_coding_transcripts(gtf_path: Path):
    """Yield (transcript_id, chromosome) for canonical protein-coding transcripts.

    'Canonical' here means the GTF line carries `tag "Ensembl_canonical"`.
    Ensembl designates exactly one canonical transcript per gene (since
    release ~107); this matches what `hpc/select_transcripts.py`'s
    BioMart query used to request.
    """
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature_type = fields[2]
            if feature_type != "transcript":
                continue
            attrs = fields[8]
            # Cheap early-out on biotype (most rows aren't protein_coding).
            m_biotype = _ATTR_BIOTYPE.search(attrs)
            if not m_biotype or m_biotype.group(1) != "protein_coding":
                continue
            if not _ATTR_TAG_CANONICAL.search(attrs):
                continue
            m_tid = _ATTR_TRANSCRIPT_ID.search(attrs)
            if not m_tid:
                continue
            yield m_tid.group(1), fields[0]  # (transcript_id, chromosome)


def select_transcripts(
    candidates: list[tuple[str, str]],
    n: int,
    max_chromosomes: int | None,
    seed: int,
) -> list[str]:
    """Sample `n` transcripts spread across up to `max_chromosomes` chromosomes.

    Strategy:
      * Group candidates by chromosome.
      * Pick up to `max_chromosomes` chromosomes randomly (from those with
        the most candidates, to avoid running out on tiny contigs).
      * Round-robin draw from those buckets until we have n, or fewer if
        the species genuinely has fewer canonical protein-coding transcripts.

    Output is the resulting transcript_ids as a plain list (order shuffled
    by the RNG; downstream pre-sort in cli.py re-clusters by chromosome).
    """
    rng = random.Random(seed)
    by_chr: dict[str, list[str]] = defaultdict(list)
    for tid, chrom in candidates:
        by_chr[chrom].append(tid)

    # Pick chromosomes. Prefer ones with more candidates so a 10-chromosome
    # spread on a small genome doesn't fail from a single-transcript contig.
    ranked_chroms = sorted(by_chr.keys(), key=lambda c: -len(by_chr[c]))
    if max_chromosomes is None or max_chromosomes >= len(ranked_chroms):
        picked_chroms = ranked_chroms
    else:
        # Take the top `max_chromosomes` by size, then shuffle to decorrelate
        # chromosome choice from genome assembly order.
        picked_chroms = ranked_chroms[:max_chromosomes]

    rng.shuffle(picked_chroms)
    for chrom in picked_chroms:
        rng.shuffle(by_chr[chrom])

    # Round-robin pull to achieve even spread
    out: list[str] = []
    idx = [0] * len(picked_chroms)
    while len(out) < n:
        progressed = False
        for i, chrom in enumerate(picked_chroms):
            if idx[i] < len(by_chr[chrom]):
                out.append(by_chr[chrom][idx[i]])
                idx[i] += 1
                progressed = True
                if len(out) >= n:
                    break
        if not progressed:
            break
    return out


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    target = p.add_mutually_exclusive_group(required=True)
    target.add_argument("--species", help="Single species slug (e.g. acanthochromis_polyacanthus).")
    target.add_argument("--species-file", type=Path, help="Newline-separated species slugs.")
    p.add_argument("--n", type=int, default=100, help="Transcripts per species (default: 100).")
    p.add_argument("--max-chromosomes", type=int, default=10,
                   help="Cap distinct chromosomes sampled per species (default: 10).")
    p.add_argument("--cache-dir", type=Path, default=Path("hpc/ensembl_gtfs"),
                   help="Where to cache downloaded GTFs.")
    p.add_argument("--out-dir", type=Path, default=Path("hpc/transcript_lists"),
                   help="Where to write per-species transcript-id files.")
    p.add_argument("--seed", type=int, default=20260417,
                   help="RNG seed for reproducible selection.")
    p.add_argument("--resume", action="store_true",
                   help="Skip species whose transcript list already exists.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    if args.species:
        species_list = [args.species]
    else:
        species_list = [s.strip() for s in args.species_file.read_text().splitlines() if s.strip()]
    logging.info("processing %d species", len(species_list))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest: list[dict] = []

    for i, species in enumerate(species_list, 1):
        out_path = args.out_dir / f"{species}.txt"
        if args.resume and out_path.exists():
            logging.info("[%d/%d] %s: skipping (already present)", i, len(species_list), species)
            continue
        try:
            gtf_path = download_gtf(species, args.cache_dir)
            candidates = list(iter_canonical_protein_coding_transcripts(gtf_path))
            picked = select_transcripts(candidates, args.n, args.max_chromosomes, args.seed)
            out_path.write_text("\n".join(picked) + "\n")
            entry = {
                "species": species,
                "status": "ok",
                "total_canonical_protein_coding": len(candidates),
                "picked": len(picked),
                "reduced": len(picked) < args.n,
            }
            logging.info("[%d/%d] %s: %d total -> %d picked%s",
                         i, len(species_list), species, len(candidates), len(picked),
                         " (REDUCED)" if entry["reduced"] else "")
        except Exception as exc:
            entry = {"species": species, "status": "error", "error": str(exc)[:240]}
            logging.warning("[%d/%d] %s: %s", i, len(species_list), species, exc)
        manifest.append(entry)

    (args.out_dir / "gtf_selection_manifest.json").write_text(json.dumps(manifest, indent=2))
    failed = sum(1 for e in manifest if e.get("status") != "ok")
    logging.info("done: %d ok, %d failed", len(manifest) - failed, failed)
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
