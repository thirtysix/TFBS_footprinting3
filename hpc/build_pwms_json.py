#!/usr/bin/env python3
"""Convert a JASPAR CORE PFM-plain-text file to tfbs_footprinter3's pwms.json format.

Input (JASPAR "jaspar" format, one record per 5 lines):

    >MA0004.1       Arnt
    A  [  4 19  0  0  0  0 ]
    C  [ 16  0 20  0  0  0 ]
    G  [  0  1  0 20  0 20 ]
    T  [  0  0  0  0 20  0 ]

Output JSON schema (matches the existing tfbs_footprinter3/data/pwms.json):

    { "<tf_name>__<matrix_id>": [[A counts...], [C counts...], [G counts...], [T counts...]] }

The composite key (`tf_name__matrix_id`) disambiguates TFs whose name is shared
across multiple JASPAR matrices (e.g. multiple MA-IDs for the same TF). Double
underscore is the separator because single underscore appears naturally inside
some TF names; double underscore or asterisk are safe delimiters per the tool
author's convention.

Usage:
    python hpc/build_pwms_json.py \\
        --input hpc/jaspar_2026/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt \\
        --output tfbs_footprinter3/data/JASPAR_2026_pwms.json
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from pathlib import Path

_HEADER_RE = re.compile(r"^>(\S+)\s+(.+?)\s*$")
_COUNTS_RE = re.compile(r"\[\s*(.*?)\s*\]")
_ROW_ORDER = ("A", "C", "G", "T")


def parse_jaspar_pfm(text: str) -> dict[str, list[list[int]]]:
    """Parse JASPAR plain-text PFM into {composite_key: [[A],[C],[G],[T]]}."""
    out: dict[str, list[list[int]]] = {}
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if not line.strip():
            i += 1
            continue
        header = _HEADER_RE.match(line)
        if not header:
            raise ValueError(f"line {i+1}: expected '>MAxxxx.v<TAB>name', got {line!r}")
        matrix_id = header.group(1)
        tf_name = header.group(2)
        counts_by_base: dict[str, list[int]] = {}
        for offset, base in enumerate(_ROW_ORDER, start=1):
            row_line = lines[i + offset]
            if not row_line.lstrip().startswith(base):
                raise ValueError(
                    f"line {i+1+offset}: expected '{base} [ ... ]' row, got {row_line!r}"
                )
            m = _COUNTS_RE.search(row_line)
            if not m:
                raise ValueError(f"line {i+1+offset}: no '[ counts ]' block in {row_line!r}")
            counts_by_base[base] = [int(x) for x in m.group(1).split()]
        # Sanity: all four rows must have equal width.
        widths = {len(counts_by_base[b]) for b in _ROW_ORDER}
        if len(widths) != 1:
            raise ValueError(
                f"{matrix_id} {tf_name}: A/C/G/T rows have mismatched widths {widths}"
            )
        composite_key = f"{tf_name}__{matrix_id}"
        if composite_key in out:
            raise ValueError(f"duplicate composite key {composite_key!r}")
        out[composite_key] = [counts_by_base[b] for b in _ROW_ORDER]
        i += 5
    return out


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, required=True,
                   help="Path to JASPAR CORE vertebrates pfm-jaspar plain-text file.")
    p.add_argument("--output", type=Path, required=True,
                   help="Where to write the JSON file (e.g. tfbs_footprinter3/data/JASPAR_2026_pwms.json).")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    text = args.input.read_text()
    pwms = parse_jaspar_pfm(text)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(pwms))
    logging.info("wrote %s (%d PWMs)", args.output, len(pwms))
    return 0


if __name__ == "__main__":
    sys.exit(main())
