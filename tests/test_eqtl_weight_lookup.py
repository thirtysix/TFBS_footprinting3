"""Regression tests for eQTL effect-magnitude weight lookup.

`gtex_weights_dict` is keyed on a 0.001-resolution grid of effect magnitudes
(2,617 keys spanning 0.0-6.023 in the shipped human table), but the magnitude
stored on each variant carries full precision -- 0.392323, 0.116223, and so on.
The old code did an exact dict lookup, `gtex_weights_dict[float(m)]`, which
raises KeyError for any magnitude not landing exactly on a grid point.

That crash was unreachable in 0.1.4 and earlier: `data_loader` looked for a
`gtex_v7` filename that no shipped tarball contains, so `converted_eqtls` was
always empty and the lookup never ran. Fixing the filename revived the layer
and immediately hit this on real data.

Resolved the way the CpG layer resolves the identical problem: bisect onto the
grid instead of requiring an exact hit.
"""
from __future__ import annotations

import numpy as np
import pytest

from tfbs_footprinter3.scoring import _eqtl_batch, _gtex_weight_lookup


@pytest.fixture
def grid():
    """A 0.001-step magnitude grid, like the shipped GTEx weights table."""
    return {round(0.001 * i, 3): float(i) for i in range(0, 1001)}


def test_exact_grid_values_are_unchanged(grid):
    got = _gtex_weight_lookup([0.0, 0.034, 0.5, 1.0], grid)
    assert list(got) == [0.0, 34.0, 500.0, 1000.0]


@pytest.mark.parametrize("mag", [0.392323, 0.388412, 0.116223, 0.0005, 0.9999])
def test_off_grid_magnitudes_do_not_raise(grid, mag):
    """These exact values crashed the previous implementation on real data."""
    got = _gtex_weight_lookup([mag], grid)
    assert np.isfinite(got).all()


def test_snaps_to_the_grid_key_at_or_above_the_query(grid):
    # 0.392323 sits between 0.392 and 0.393; bisect_left semantics take 0.393,
    # matching cpg_weights_summing.
    assert _gtex_weight_lookup([0.392323], grid)[0] == pytest.approx(393.0)
    assert _gtex_weight_lookup([0.392], grid)[0] == pytest.approx(392.0)


def test_magnitude_above_the_table_clamps_to_the_top(grid):
    assert _gtex_weight_lookup([12.5], grid)[0] == pytest.approx(1000.0)


def test_negative_magnitudes_use_absolute_value(grid):
    assert _gtex_weight_lookup([-0.25], grid)[0] == _gtex_weight_lookup([0.25], grid)[0]


def test_empty_table_gives_zero_weights(grid):
    assert list(_gtex_weight_lookup([0.1, 0.2], {})) == [0.0, 0.0]


def test_eqtl_batch_survives_off_grid_magnitudes(grid):
    """End-to-end through the vectorized path that actually crashed."""
    motif_starts = np.array([10, 40], dtype=np.int64)
    motif_ends = np.array([20, 50], dtype=np.int64)
    eqtl_starts = np.array([12, 45], dtype=np.int64)
    eqtl_ends = np.array([13, 46], dtype=np.int64)
    eqtl_mags = np.array([0.392323, 0.116223], dtype=np.float64)

    out = _eqtl_batch(motif_starts, motif_ends, eqtl_starts, eqtl_ends,
                      eqtl_mags, grid, 1.5)
    assert out.shape == (2,)
    assert np.isfinite(out).all()
    # first motif overlaps only the first eQTL, second only the second
    assert out[0] == pytest.approx(393.0 + 1.5)
    assert out[1] == pytest.approx(117.0 + 1.5)


def test_scalar_reference_agrees_with_batch_off_grid(grid):
    """The scalar test-reference must not diverge from the production path.

    `eqtls_weights_summing_v` is what tests/test_scoring_vectorized.py checks
    `_eqtl_batch` against, so if only one of them snaps to the grid the
    agreement test silently stops meaning anything.
    """
    from tfbs_footprinter3.scoring import eqtls_weights_summing_v

    motif_start, motif_end = 10, 20
    eqtl_starts = np.array([12, 15], dtype=np.int64)
    eqtl_ends = np.array([13, 16], dtype=np.int64)
    eqtl_mags = np.array([0.392323, 0.116223], dtype=np.float64)

    scalar = eqtls_weights_summing_v(1.5, motif_start, motif_end,
                                     eqtl_starts, eqtl_ends, eqtl_mags, grid)
    batch = _eqtl_batch(np.array([motif_start]), np.array([motif_end]),
                        eqtl_starts, eqtl_ends, eqtl_mags, grid, 1.5)[0]
    assert scalar == pytest.approx(batch)
