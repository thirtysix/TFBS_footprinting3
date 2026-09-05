"""Regression tests for TF-name resolution into the CAGE-correlation data.

Since 0.1.0 the TF identifier carried through scoring is the composite
`{tf_name}__{matrix_id}`, upper-cased by `cli.py`. But
`{species}.CAGE.jasparTFs.dict.json` still holds the 575 *bare* JASPAR 2018
names, so the membership test `tf_name in TF_cage_dict` matched none of the
1019 JASPAR 2026 motifs and the expression-correlation weight was identically
zero for every hit -- one of the seven advertised evidence layers, silently
contributing nothing.

Measured against the shipped human data before the fix: 0 / 1019 motifs
matched. After stripping the matrix id and matching case-insensitively:
597 / 1019.
"""
from __future__ import annotations

import pytest

from tfbs_footprinter3.scoring import (
    cage_correlations_summing_preparation,
    resolve_tf_cage_key,
)


@pytest.fixture
def tf_cage():
    """Bare JASPAR 2018 style keys, as the shipped file actually has them."""
    return {
        "SRF": [["hg_CAGE_1", 1.0]],
        "MNT": [["hg_CAGE_2", 1.0]],
        "NR1H3::RXRA": [["hg_CAGE_3", 1.0]],
        "Arnt": [["hg_CAGE_4", 1.0]],
    }


def test_composite_key_resolves_to_bare_name(tf_cage):
    assert resolve_tf_cage_key("SRF__MA0083.3", tf_cage) == "SRF"


def test_uppercased_composite_resolves(tf_cage):
    """cli.py upper-cases the matrix dict keys, so ARNT must find 'Arnt'."""
    assert resolve_tf_cage_key("ARNT__MA0004.1", tf_cage) == "Arnt"


def test_dimer_names_survive_the_split(tf_cage):
    """'::' must not be confused with the '__' matrix-id separator."""
    assert resolve_tf_cage_key("NR1H3::RXRA__MA0074.1", tf_cage) == "NR1H3::RXRA"


def test_plain_name_still_works(tf_cage):
    assert resolve_tf_cage_key("MNT", tf_cage) == "MNT"


def test_unknown_tf_returns_none(tf_cage):
    assert resolve_tf_cage_key("NOTATF__MA9999.1", tf_cage) is None


def test_preparation_returns_cages_for_a_composite_name(tf_cage):
    cage_dict = {"MYGENE": [["hg_CAGE_9", 1.0]]}
    target, tf_cages = cage_correlations_summing_preparation(
        "MYGENE", "ENST0", cage_dict, tf_cage, "SRF__MA0083.3")
    assert target == ["CAGE_9"]
    assert tf_cages == ["CAGE_1"], "composite TF name must still find its CAGE peaks"


def test_preparation_was_broken_for_composite_names_before(tf_cage):
    """Pin the failure mode: a bare membership test finds nothing."""
    assert "SRF__MA0083.3" not in tf_cage
