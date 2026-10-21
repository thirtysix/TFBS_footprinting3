#!/usr/bin/env python3
"""Pin the SHA-256 of every data tarball in the public bucket into the package.

Writes tfbs_footprinter3/data/experimental_data_sha256.json, which data_loader checks each
download against before extracting it. The checksums ship inside the released package, not in
the bucket, so a write to the bucket alone cannot change both a tarball and its checksum.

Run it after every publish_to_s3.py upload, commit the JSON, and release the package. Users on an
older release then get a clear "no pinned checksum" error for new species (or can opt in with
TFBS_FOOTPRINTER3_ALLOW_UNPINNED=1) instead of silently extracting unverified data.

Reads the bucket anonymously over HTTPS (it is public-read), streaming each tarball through
SHA-256 without keeping it on disk: about 2 GB of downloads for the current ~317 tarballs.

Usage:
    python hpc/build_data_manifest.py                  # all top-level tarballs
    python hpc/build_data_manifest.py --jobs 4
    python hpc/build_data_manifest.py --only homo_sapiens.tar.gz mus_musculus.tar.gz   # update some
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import sys
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BUCKET_URL = "https://s3.us-east-2.amazonaws.com/tfbssexperimentaldata"
OUT = Path(__file__).resolve().parents[1] / "tfbs_footprinter3" / "data" / "experimental_data_sha256.json"


def list_tarballs() -> list[str]:
    """Top-level *.tar.gz keys (backups under backup/ are left out)."""
    keys, token = [], None
    while True:
        url = f"{BUCKET_URL}/?list-type=2" + (f"&continuation-token={urllib.parse.quote(token)}" if token else "")
        xml = urllib.request.urlopen(url, timeout=60).read().decode()
        keys += [k for k in re.findall(r"<Key>([^<]+)</Key>", xml) if k.endswith(".tar.gz") and "/" not in k]
        m = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", xml)
        if not m:
            return sorted(keys)
        token = m[1]


def hash_one(key: str) -> tuple[str, str, int]:
    digest, size = hashlib.sha256(), 0
    with urllib.request.urlopen(f"{BUCKET_URL}/{urllib.parse.quote(key)}", timeout=300) as r:
        for chunk in iter(lambda: r.read(1 << 20), b""):
            digest.update(chunk)
            size += len(chunk)
    return key, digest.hexdigest(), size


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--jobs", type=int, default=4, help="parallel downloads (default 4)")
    p.add_argument("--only", nargs="+", help="re-hash just these keys, keeping the rest of the manifest")
    a = p.parse_args(argv)

    current = json.loads(OUT.read_text()) if OUT.exists() else {"files": {}}
    keys = a.only or list_tarballs()
    files = dict(current.get("files", {})) if a.only else {}
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for i, (key, sha, size) in enumerate(pool.map(hash_one, keys), 1):
            files[key] = {"sha256": sha, "size": size}
            print(f"[{i}/{len(keys)}] {key} {sha[:12]} {size:,}", file=sys.stderr)
    OUT.write_text(json.dumps({"bucket": BUCKET_URL, "generated": dt.date.today().isoformat(),
                               "files": dict(sorted(files.items()))}, indent=1) + "\n")
    print(f"wrote {OUT} ({len(files)} tarballs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
