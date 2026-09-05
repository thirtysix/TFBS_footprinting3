"""Regression tests for per-chromosome data-file discovery.

Motivation: `data_loader.species_specific_data` built the GTEx eQTL filename
from a hard-coded `"gtex_v7"` release string, while every shipped tarball has
carried `gtex_v8` since the data was rebuilt. The lookup was a bare
`os.path.exists` with no `else`, so the mismatch silently zeroed the eQTL
weight for every hit rather than raising -- one of the seven advertised
evidence layers contributed nothing, undetectably.

These tests pin the replacement: discovery by species + chromosome + a stable
token, version-agnostic.
"""
from __future__ import annotations

import pytest

from tfbs_footprinter3.data_loader import find_chromosome_data_file


def _make_gtex_dir(tmp_path, release="gtex_v8", chromosomes=("1", "2", "10", "19", "X")):
    d = tmp_path / "gtex_data"
    d.parent.mkdir(parents=True, exist_ok=True)
    d.mkdir()
    for c in chromosomes:
        (d / f"homo_sapiens.{release}.Chr{c}.min_unique.eqtls.grch38.msg").touch()
    # the weights file lives in the same directory and must never be picked up
    (d / f"homo_sapiens.{release}.gtex_weights.json").touch()
    return d


def test_finds_file_for_any_release(tmp_path):
    """The whole point: a release bump must not disable the layer."""
    for release in ("gtex_v7", "gtex_v8", "gtex_v10"):
        d = _make_gtex_dir(tmp_path / release, release=release)
        found = find_chromosome_data_file(d, "homo_sapiens", "1", "eqtls")
        assert found is not None, f"release {release} not discovered"
        assert found.endswith(f"homo_sapiens.{release}.Chr1.min_unique.eqtls.grch38.msg")


def test_chr1_does_not_match_chr10_or_chr19(tmp_path):
    """`.Chr1.` must not prefix-match Chr10..Chr19."""
    d = _make_gtex_dir(tmp_path)
    found = find_chromosome_data_file(d, "homo_sapiens", "1", "eqtls")
    assert found.endswith(".Chr1.min_unique.eqtls.grch38.msg")


@pytest.mark.parametrize("query,expected", [("10", "Chr10"), ("x", "ChrX"), ("X", "ChrX")])
def test_chromosome_is_upcased(tmp_path, query, expected):
    d = _make_gtex_dir(tmp_path)
    found = find_chromosome_data_file(d, "homo_sapiens", query, "eqtls")
    assert f".{expected}." in found


def test_weights_file_is_not_matched(tmp_path):
    """gtex_weights.json shares the directory but has no .Chr<N>. segment."""
    d = _make_gtex_dir(tmp_path)
    for c in ("1", "2", "10", "19", "X"):
        assert "gtex_weights" not in find_chromosome_data_file(d, "homo_sapiens", c, "eqtls")


def test_newest_release_wins(tmp_path):
    d = _make_gtex_dir(tmp_path, release="gtex_v7", chromosomes=("1",))
    (d / "homo_sapiens.gtex_v8.Chr1.min_unique.eqtls.grch38.msg").touch()
    found = find_chromosome_data_file(d, "homo_sapiens", "1", "eqtls")
    assert "gtex_v8" in found


def test_missing_chromosome_and_missing_dir_return_none(tmp_path):
    d = _make_gtex_dir(tmp_path)
    assert find_chromosome_data_file(d, "homo_sapiens", "Y", "eqtls") is None
    assert find_chromosome_data_file(d, "mus_musculus", "1", "eqtls") is None
    assert find_chromosome_data_file(tmp_path / "absent", "homo_sapiens", "1", "eqtls") is None


def test_shipped_human_data_is_discoverable():
    """End-to-end against the real installed data, if present.

    This is the assertion that would have caught the original bug.
    """
    import os

    from tfbs_footprinter3.data_loader import script_dir

    d = os.path.join(script_dir, "data", "homo_sapiens", "gtex_data")
    if not os.path.exists(d):
        pytest.skip("human experimental data not downloaded")

    for chromosome in ("1", "10", "X"):
        found = find_chromosome_data_file(d, "homo_sapiens", chromosome, "eqtls")
        assert found is not None, f"no eQTL file discovered for chr{chromosome}"
        assert os.path.exists(found)
