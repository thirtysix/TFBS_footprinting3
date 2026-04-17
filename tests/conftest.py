"""Shared pytest fixtures.

Keep this file minimal — unit tests should pull their own inputs; fixtures
here are for things that multiple suites really need.
"""
from __future__ import annotations

import pytest


@pytest.fixture
def pwm_dict() -> dict[str, int]:
    """The ACGT -> row-index lookup used throughout tfbs_footprinter3.

    Matches the mapping the tool establishes when calling PWM_scorer
    (nucleotide letter -> pwm row index).
    """
    return {"A": 0, "C": 1, "G": 2, "T": 3}


@pytest.fixture
def uniform_bg() -> dict[str, float]:
    """Equal 0.25 background frequency for each base."""
    return {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
