"""Round-trip test for the Parquet output writer.

Ensures target_species_hits_table_writer_parquet produces a Parquet file
that round-trips back to the TF name and CAS score values hpc/cas_only.py
expects (those two columns are the HPC pipeline's only real dependency
on the tool's output).
"""
from __future__ import annotations

from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")
pytest.importorskip("pyarrow")

from tfbs_footprinter3.output import (  # noqa: E402
    _PARQUET_HEADER,
    target_species_hits_table_writer_parquet,
)


def _sample_hit(tf_name, pwm, cas):
    """Build a hit row in the same shape find_clusters + sort produce.

    Layout: [tf_name, seq, strand, start, end, TSSstart, TSSend,
             PWM_score, PWMp, CAS, CASp, species, cage, eqtls,
             atac, metaclust, cpg, corr]
    """
    return [tf_name, "ACGTACGTAC", "+1", 100, 110, -800, -790,
            pwm, "0.01", cas, "0.001",
            0, 4.33, 0, 0, 3.69, 4.29, 0]


def test_empty_list_writes_header_only(tmp_path):
    path = tmp_path / "empty.parquet"
    target_species_hits_table_writer_parquet([], str(path))
    df = pd.read_parquet(path)
    assert list(df.columns) == _PARQUET_HEADER
    assert len(df) == 0


def test_roundtrip_preserves_tf_and_cas(tmp_path):
    path = tmp_path / "rt.parquet"
    rows = [
        _sample_hit("CTCF", 15.2, 42.7),
        _sample_hit("SP1", 14.1, 39.9),
        _sample_hit("SP1", 13.5, 28.8),
    ]
    target_species_hits_table_writer_parquet(rows, str(path))
    df = pd.read_parquet(path)

    assert list(df["binding prot."]) == ["CTCF", "SP1", "SP1"]
    assert list(df["combined affinity score"]) == pytest.approx([42.7, 39.9, 28.8])


def test_int_and_float_columns_are_typed(tmp_path):
    path = tmp_path / "typed.parquet"
    rows = [_sample_hit("CTCF", 15.2, 42.7)]
    target_species_hits_table_writer_parquet(rows, str(path))
    df = pd.read_parquet(path)

    # Coordinates typed as integer columns (not object / string)
    for col in ("start", "end", "TSS-relative start", "TSS-relative end"):
        assert df[col].dtype.kind in "iu", f"{col} should be integer, got {df[col].dtype}"
    # Numeric score/weight columns are float
    for col in ("PWM score", "combined affinity score",
                "species weights sum", "cage weights sum"):
        assert df[col].dtype.kind == "f", f"{col} should be float, got {df[col].dtype}"


def test_pvalue_scientific_conversion(tmp_path):
    path = tmp_path / "sci.parquet"
    # Small p-value triggers scientific notation in the writer; carries
    # through to the parquet string cell.
    row = _sample_hit("CTCF", 15.2, 42.7)
    row[8] = "0.00001"   # PWM pval
    row[10] = "0.00002"  # CAS pval
    target_species_hits_table_writer_parquet([row], str(path))
    df = pd.read_parquet(path)

    assert df["p-value"].iloc[0] == "1.000e-5"
    assert df["combined affinity score p-value"].iloc[0] == "2.000e-5"


def test_parquet_file_is_smaller_than_equivalent_csv(tmp_path):
    """Sanity: a meaningful-sized write should show Parquet's size advantage.

    Not a strict contract (compression ratio varies with data), but a
    100-row batch of uniform hits should fit well under an equivalent
    CSV written via csv.writer.writerows.
    """
    from tfbs_footprinter3.output import target_species_hits_table_writer
    rows_csv = [_sample_hit(f"TF_{i}", 10.0 + i * 0.01, 30.0 + i * 0.01) for i in range(1000)]
    rows_parquet = [r.copy() for r in rows_csv]

    csv_path = tmp_path / "out.csv"
    parquet_path = tmp_path / "out.parquet"

    target_species_hits_table_writer(rows_csv, str(csv_path))
    target_species_hits_table_writer_parquet(rows_parquet, str(parquet_path))

    assert Path(parquet_path).stat().st_size < Path(csv_path).stat().st_size


# --- Slim parquet writer -----------------------------------------------------

from tfbs_footprinter3.output import (  # noqa: E402
    target_species_hits_table_writer_slim_parquet,
)


def _slim_dict(items):
    """Build the (pwm_arr, cas_arr) dict shape find_clusters produces in slim mode."""
    import numpy as np  # noqa: PLC0415
    return {tf: (np.asarray(pwm, dtype=np.float32),
                 np.asarray(cas, dtype=np.float32))
            for tf, pwm, cas in items}


def test_slim_parquet_empty_dict(tmp_path):
    path = tmp_path / "slim_empty.parquet"
    target_species_hits_table_writer_slim_parquet({}, str(path))
    df = pd.read_parquet(path)
    assert list(df.columns) == ["binding prot.", "PWM score", "combined affinity score"]
    assert len(df) == 0


def test_slim_parquet_roundtrip(tmp_path):
    path = tmp_path / "slim.parquet"
    slim = _slim_dict([
        ("CTCF__MA1930.2", [15.2, 14.8, 14.6], [42.7, 42.1, 41.9]),
        ("SP1__MA0079.6", [14.1, 13.5], [39.9, 28.8]),
    ])
    target_species_hits_table_writer_slim_parquet(slim, str(path))
    df = pd.read_parquet(path)

    # Preserves per-TF grouping (streamed TF-by-TF, no sort)
    assert list(df["binding prot."]) == ["CTCF__MA1930.2", "CTCF__MA1930.2", "CTCF__MA1930.2",
                                          "SP1__MA0079.6", "SP1__MA0079.6"]
    assert list(df["PWM score"]) == pytest.approx([15.2, 14.8, 14.6, 14.1, 13.5], rel=1e-5)
    assert list(df["combined affinity score"]) == pytest.approx([42.7, 42.1, 41.9, 39.9, 28.8], rel=1e-5)


def test_slim_parquet_typed(tmp_path):
    path = tmp_path / "slim_typed.parquet"
    slim = _slim_dict([("CTCF", [15.2], [42.7])])
    target_species_hits_table_writer_slim_parquet(slim, str(path))
    df = pd.read_parquet(path)

    # pandas 2.x defaults string columns to StringDtype; older pandas uses
    # object dtype. Both are valid — we only care that the column is
    # string-typed, not numeric or something weird.
    assert pd.api.types.is_string_dtype(df["binding prot."])
    assert df["PWM score"].dtype.kind == "f"
    assert df["combined affinity score"].dtype.kind == "f"
