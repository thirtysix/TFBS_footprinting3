#!/usr/bin/env python3
"""Stage E: inject CAS_pvalues JSON into per-species S3 tarballs and re-upload.

For each species, download s3://tfbssexperimentaldata/{species}.tar.gz,
add Stage D's CAS_pvalues.0.1.tf_ls.json into its {species}/ subdir
(without disturbing the existing experimental data files), and re-upload.

After re-upload, any user running `tfbs_footprinter3 -update` will pull
the new tarball, and the species_specific_data() loader at
tfbs_footprinter3.py:755 will find the CAS_pvalues JSON and start
emitting combined-affinity-score p-values for non-human analyses with
no client-side code change.

Uses the `aws` CLI via subprocess so there's no boto3 dependency to add
to the package. Assumes `aws` is on PATH and credentials are configured
(env vars, ~/.aws/credentials, or an instance role — any of the usual).

Usage:
    # Dry run: download, inspect, repack locally, don't upload
    python hpc/publish_to_s3.py --species mus_musculus --dry-run

    # Upload for one species
    python hpc/publish_to_s3.py --species mus_musculus

    # Batch upload for every species with an artifact on disk
    python hpc/publish_to_s3.py --species all
"""
from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

BUCKET = "tfbssexperimentaldata"
CAS_FILENAME = "CAS_pvalues.0.1.tf_ls.json"


def run(cmd: list[str]) -> None:
    """Run a subprocess, stream output, raise on non-zero exit."""
    logging.info("$ %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def s3_key(species: str) -> str:
    return f"s3://{BUCKET}/{species}.tar.gz"


def process_one(species: str, artifact_root: Path, dry_run: bool) -> dict:
    """Download species tarball, inject CAS JSON, re-upload (unless dry run)."""
    cas_json_path = artifact_root / species / CAS_FILENAME
    if not cas_json_path.exists():
        return {"species": species, "status": "no_artifact", "path": str(cas_json_path)}

    with tempfile.TemporaryDirectory(prefix=f"tfbs_publish_{species}_") as tmp:
        tmpdir = Path(tmp)
        local_tar = tmpdir / f"{species}.tar.gz"
        extract_dir = tmpdir / "extract"
        extract_dir.mkdir()

        # 1. Download the existing tarball
        run(["aws", "s3", "cp", s3_key(species), str(local_tar)])

        # 2. Extract
        with tarfile.open(local_tar, "r:gz") as tar:
            # Safety: refuse absolute paths / .. traversal
            for member in tar.getmembers():
                if member.name.startswith(("/", "..")) or ".." in Path(member.name).parts:
                    raise RuntimeError(f"unsafe path in tarball: {member.name}")
            tar.extractall(extract_dir)

        # 3. Locate the species directory inside the tar
        species_dir = extract_dir / species
        if not species_dir.is_dir():
            # Some tarballs may have just the files at root; tolerate both
            candidates = [p for p in extract_dir.iterdir() if p.is_dir() and p.name == species]
            if not candidates:
                raise RuntimeError(f"species subdir {species!r} not found in {local_tar}")
            species_dir = candidates[0]

        # 4. Inject the CAS JSON
        dst = species_dir / CAS_FILENAME
        dst.write_bytes(cas_json_path.read_bytes())
        logging.info("injected %s into %s", CAS_FILENAME, species_dir)

        # 5. Re-tar
        new_tar = tmpdir / f"{species}.new.tar.gz"
        with tarfile.open(new_tar, "w:gz") as tar:
            tar.add(species_dir, arcname=species)
        logging.info("rebuilt tarball: %d bytes", new_tar.stat().st_size)

        if dry_run:
            # Preserve for inspection
            preserved = artifact_root / f"{species}.dry-run.tar.gz"
            preserved.write_bytes(new_tar.read_bytes())
            logging.info("DRY-RUN: preserved tarball at %s", preserved)
            return {"species": species, "status": "dry_run", "local": str(preserved)}

        # 6. Upload (overwrite)
        run(["aws", "s3", "cp", str(new_tar), s3_key(species)])
        return {"species": species, "status": "uploaded", "bytes": new_tar.stat().st_size}


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species", type=str, required=True,
                   help="Species slug (or 'all' for every directory under --artifacts).")
    p.add_argument("--artifacts", type=Path, default=Path("hpc/artifacts"),
                   help="Directory with {species}/CAS_pvalues.0.1.tf_ls.json files from Stage D.")
    p.add_argument("--dry-run", action="store_true",
                   help="Rebuild tarballs locally, skip the s3 upload.")
    p.add_argument("--manifest", type=Path, default=Path("hpc/artifacts/publish_manifest.json"))
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    if args.species == "all":
        species_list = sorted(
            p.name for p in args.artifacts.iterdir()
            if p.is_dir() and (p / CAS_FILENAME).exists()
        )
    else:
        species_list = [args.species]

    results = []
    failed = 0
    for sp in species_list:
        try:
            results.append(process_one(sp, args.artifacts, args.dry_run))
        except Exception as exc:
            failed += 1
            logging.error("%s: %s", sp, exc)
            results.append({"species": sp, "status": "error", "error": str(exc)})

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(results, indent=2))
    logging.info("manifest: %s (%d ok, %d failed)",
                 args.manifest, len(species_list) - failed, failed)
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
