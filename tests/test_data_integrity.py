"""Downloaded data tarballs: pinned SHA-256 and safe extraction (2026-10-21).

The bucket's write key could replace a tarball; before this, the tool extracted whatever it
downloaded with a bare extractall(), so a crafted tarball could write anywhere the user can.
"""
from __future__ import annotations

import hashlib
import io
import json
import tarfile

import pytest

from tfbs_footprinter3 import data_loader as dl


def make_tar(path, members):
    """members: list of (name, bytes) for files, or (name, ('symlink', target))."""
    with tarfile.open(path, "w:gz") as tar:
        for name, body in members:
            info = tarfile.TarInfo(name)
            if isinstance(body, tuple):
                info.type, info.linkname = tarfile.SYMTYPE, body[1]
                tar.addfile(info)
            else:
                info.size = len(body)
                tar.addfile(info, io.BytesIO(body))
    return path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_a_pinned_tarball_that_matches_is_accepted(tmp_path):
    t = make_tar(tmp_path / "mus_musculus.tar.gz", [("mus_musculus/a.json", b"{}")])
    dl.verify_download(str(t), "mus_musculus.tar.gz", {"mus_musculus.tar.gz": {"sha256": sha(t)}})


def test_a_changed_tarball_is_refused(tmp_path):
    t = make_tar(tmp_path / "mus_musculus.tar.gz", [("mus_musculus/a.json", b"{}")])
    with pytest.raises(dl.DataIntegrityError, match="does not match its pinned checksum"):
        dl.verify_download(str(t), "mus_musculus.tar.gz", {"mus_musculus.tar.gz": {"sha256": "0" * 64}})


def test_an_unpinned_tarball_needs_an_explicit_opt_in(tmp_path, monkeypatch):
    t = make_tar(tmp_path / "new_species.tar.gz", [("new_species/a.json", b"{}")])
    monkeypatch.delenv(dl.ALLOW_UNPINNED_ENV, raising=False)
    with pytest.raises(dl.DataIntegrityError, match="no pinned checksum"):
        dl.verify_download(str(t), "new_species.tar.gz", {})
    monkeypatch.setenv(dl.ALLOW_UNPINNED_ENV, "1")
    dl.verify_download(str(t), "new_species.tar.gz", {})
    with pytest.raises(dl.DataIntegrityError):                   # the opt-in never excuses a mismatch
        dl.verify_download(str(t), "new_species.tar.gz", {"new_species.tar.gz": {"sha256": "0" * 64}})


@pytest.mark.parametrize("kind", ["dotdot", "absolute", "symlink"])
@pytest.mark.parametrize("builtin_filter", [True, False], ids=["data-filter", "fallback"])
def test_nothing_lands_outside_the_data_dir(tmp_path, monkeypatch, kind, builtin_filter):
    """Refused, or (the built-in filter, for an absolute path) rewritten to land inside: either
    way nothing is written outside the destination and no link out is created."""
    if not builtin_filter:                       # interpreters without tarfile's "data" filter
        monkeypatch.delattr(dl.tarfile, "data_filter", raising=False)
    outside = tmp_path / "outside"
    members = {"dotdot": [("../outside/escaped.txt", b"x")],
               "absolute": [(str(outside / "escaped.txt"), b"x")],
               "symlink": [("mus_musculus/link", ("symlink", "/etc/passwd"))]}[kind]
    t = make_tar(tmp_path / "evil.tar.gz", members)
    dest = tmp_path / "data"
    dest.mkdir()
    try:
        dl.safe_extract(str(t), str(dest))
        refused = False
    except dl.DataIntegrityError:
        refused = True
    assert not (outside / "escaped.txt").exists()
    assert not (dest / "mus_musculus" / "link").is_symlink()
    assert refused or kind == "absolute"


@pytest.mark.parametrize("builtin_filter", [True, False], ids=["data-filter", "fallback"])
def test_a_normal_tarball_extracts(tmp_path, monkeypatch, builtin_filter):
    if not builtin_filter:
        monkeypatch.delattr(dl.tarfile, "data_filter", raising=False)
    t = make_tar(tmp_path / "ok.tar.gz", [("mus_musculus/a.json", b"{}")])
    dest = tmp_path / "data"
    dest.mkdir()
    dl.safe_extract(str(t), str(dest))
    assert (dest / "mus_musculus" / "a.json").read_bytes() == b"{}"


def test_a_failed_check_deletes_the_download_and_extracts_nothing(tmp_path, monkeypatch):
    src = make_tar(tmp_path / "src.tar.gz", [("mus_musculus/a.json", b"{}")])
    def fake_retrieve(url, path):
        with open(path, "wb") as fh:
            fh.write(src.read_bytes())
    monkeypatch.setattr(dl.urllib.request, "urlretrieve", fake_retrieve)
    monkeypatch.setattr(dl, "_pinned_hashes", lambda: {"mus_musculus.tar.gz": {"sha256": "0" * 64}})
    dest = tmp_path / "data"
    dest.mkdir()
    tar_path = dest / "mus_musculus.tar.gz"
    with pytest.raises(dl.DataIntegrityError):
        dl.fetch_and_extract("https://example.invalid/x", str(tar_path), str(dest))
    assert not tar_path.exists() and not (dest / "mus_musculus").exists()


def test_the_shipped_manifest_pins_the_main_species():
    files = json.load(open(dl.DATA_MANIFEST))["files"]
    for name in ("homo_sapiens.tar.gz", "mus_musculus.tar.gz", "danio_rerio.tar.gz"):
        assert len(files[name]["sha256"]) == 64
