"""SQLite-backed cache for Ensembl REST JSON responses.

Used by Stage B (prefetch.py) to populate the cache, and by Stage C
(cas_only.py) to replace tfbs_footprinter3.ensemblrest with a cache-aware
version at runtime, so SLURM array jobs don't hammer the Ensembl REST API.

Cache format: a single sqlite file with one table `responses(url PRIMARY
KEY, fetched_at, body)`. The body is a JSON string — we parse at read time.
This keeps the cache introspectable (sqlite CLI, DB Browser) and rsync-able
between a submission node and compute nodes.
"""
from __future__ import annotations

import contextlib
import json
import logging
import sqlite3
import time
from pathlib import Path

import httplib2

ENSEMBL_SERVER = "http://rest.ensembl.org"
DEFAULT_RETRY_STATUSES = {429, 500, 502, 503, 504}
DEFAULT_MAX_ATTEMPTS = 4
DEFAULT_BACKOFF_BASE = 5.0


class CachedEnsemblClient:
    """Sqlite-backed cache for Ensembl REST JSON calls.

    Use get_json(url) to fetch; the client serves from cache when present
    and falls back to httplib2 otherwise, storing every success.
    """

    def __init__(
        self,
        db_path: Path | str,
        *,
        read_only: bool = False,
        max_attempts: int = DEFAULT_MAX_ATTEMPTS,
        backoff_base: float = DEFAULT_BACKOFF_BASE,
        rate_limit_sleep: float = 0.07,
    ):
        self.db_path = Path(db_path)
        self.read_only = read_only
        self.max_attempts = max_attempts
        self.backoff_base = backoff_base
        self.rate_limit_sleep = rate_limit_sleep
        self._http = httplib2.Http()
        self._conn = self._open_conn()

    def _open_conn(self) -> sqlite3.Connection:
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(
            f"file:{self.db_path}?mode={'ro' if self.read_only else 'rwc'}",
            uri=True,
            isolation_level=None,
            check_same_thread=False,
        )
        if not self.read_only:
            conn.execute("PRAGMA journal_mode=WAL;")
            conn.execute(
                "CREATE TABLE IF NOT EXISTS responses ("
                "url TEXT PRIMARY KEY, "
                "fetched_at REAL NOT NULL, "
                "body TEXT NOT NULL"
                ");"
            )
        return conn

    def get_cached(self, url: str) -> dict | list | None:
        row = self._conn.execute(
            "SELECT body FROM responses WHERE url = ?", (url,)
        ).fetchone()
        if row is None:
            return None
        return json.loads(row[0])

    def store(self, url: str, body: dict | list) -> None:
        if self.read_only:
            raise RuntimeError("store() called on a read-only client")
        self._conn.execute(
            "INSERT OR REPLACE INTO responses (url, fetched_at, body) VALUES (?, ?, ?)",
            (url, time.time(), json.dumps(body)),
        )

    def get_json(self, url: str) -> dict | list:
        """Return cached JSON for `url`, fetching via REST if absent."""
        cached = self.get_cached(url)
        if cached is not None:
            return cached
        if self.read_only:
            raise KeyError(f"URL not in cache: {url}")
        body = self._fetch(url)
        self.store(url, body)
        time.sleep(self.rate_limit_sleep)  # polite throttle
        return body

    def _fetch(self, url: str) -> dict | list:
        last_exc: Exception | None = None
        for attempt in range(self.max_attempts):
            try:
                resp, raw = self._http.request(url, method="GET")
                status = int(resp.status)
                if status in DEFAULT_RETRY_STATUSES:
                    raise RuntimeError(f"HTTP {status} (retryable)")
                if status >= 400:
                    raise RuntimeError(f"HTTP {status}: {raw[:200]!r}")
                return json.loads(raw)
            except Exception as exc:
                last_exc = exc
                if attempt < self.max_attempts - 1:
                    delay = self.backoff_base * (2 ** attempt)
                    logging.info("Retry %s in %.0fs: %s", url, delay, str(exc)[:120])
                    time.sleep(delay)
        raise RuntimeError(f"Failed after {self.max_attempts} attempts: {url} ({last_exc})")

    def close(self) -> None:
        with contextlib.suppress(Exception):
            self._conn.close()

    def __enter__(self) -> CachedEnsemblClient:
        return self

    def __exit__(self, *exc) -> None:
        self.close()


# ---- monkey-patch helper ----


def patch_tfbs_footprinter3(client: CachedEnsemblClient) -> None:
    """Replace the Ensembl REST entry point with a cache-aware wrapper.

    Rebinds `tfbs_footprinter3.ensembl.ensemblrest` (the canonical
    location) AND `tfbs_footprinter3.tfbs_footprinter3.ensemblrest`
    (the legacy re-export). Pipeline-side callers in
    `tfbs_footprinter3.pipeline` use attribute access
    (`ensembl.ensemblrest(...)`) so they pick up the patch transparently.

    All known call sites use output_type='json'; the 'fasta' branch in
    ensemblrest is dead as of v0.0.7.
    """
    from tfbs_footprinter3 import ensembl
    from tfbs_footprinter3 import tfbs_footprinter3 as tff

    def _cached_ensemblrest(query_type, options, output_type, ensembl_id=None, log=False):
        full_query = ENSEMBL_SERVER + query_type + (ensembl_id or "") + options
        if log:
            logging.info("Ensembl REST (cached): %s", full_query)
        return client.get_json(full_query)

    ensembl.ensemblrest = _cached_ensemblrest
    tff.ensemblrest = _cached_ensemblrest
