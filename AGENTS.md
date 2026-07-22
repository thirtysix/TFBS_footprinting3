# AGENTS.md

Orientation for AI coding agents (and humans) arriving at this repo: what it is, how to run it, and the conventions worth knowing before making changes.

## What this is
TFBSFootprinter3: a Python tool that predicts conserved transcription-factor binding sites (TFBSs) in vertebrate promoters. Candidate sites are first scored against the JASPAR 2026 CORE vertebrates non-redundant motif set (1019 PWMs), then reweighted by seven contextual data types: sequence conservation across homologous species, proximity to CAGE-supported TSSs (FANTOM), TF-to-target expression correlation, ChIP-seq TFBS metaclusters (GTRD), eQTLs (GTEx), CpG content, and ATAC-seq open chromatin (ENCODE). Analyses are keyed on Ensembl transcript IDs (up to 316 species in Ensembl release 113). Output is a predictions table plus an SVG promoter figure.

It ships as an installable package (conda and PyPI) exposing the `tfbs_footprinter3` console command. Per-species reference and experimental data are downloaded on demand from AWS S3.

## Install
Conda (recommended), Python 3.9-3.12, 3.12 recommended:
```bash
conda create -n tfbs python=3.12
conda activate tfbs
conda install -c thirtysix -c conda-forge tfbs-footprinting3
```
PyPI:
```bash
pip install tfbs_footprinting3
# optional parquet output (-of parquet): brings in pyarrow
pip install 'tfbs_footprinting3[parquet]'
```

## Run
```bash
# analysis driven by a CSV of transcript IDs plus per-transcript arguments
tfbs_footprinter3 -t sample_csv.csv
# or a plain text file of Ensembl transcript IDs with shared CLI arguments
tfbs_footprinter3 -t sample_ids.txt -tfs jaspar_tf_ids_selected.txt -pb 900 -pa 100 -p 0.01
# download / refresh the experimental data files
tfbs_footprinter3 -update
```
Sample inputs live in `sample_analysis/`. Runtime deps: biopython, numpy, matplotlib, pandas, httplib2, msgpack, wget (pyarrow optional).

## Test and lint
Unit tests are offline and numeric (PWM scoring, CAS p-value lookup, scalar-vs-vectorized scoring agreement, parquet round-trip). Some tests import the top-level `hpc` package, so run pytest from the repo root:
```bash
pytest tests/
ruff check .
```

## Layout
- `tfbs_footprinter3/` - the package. The former monolith is split into `cli`, `pipeline`, `scoring`, `plotting`, `pwm`, `data_loader`, `ensembl`, `translators`, and others. `tfbs_footprinter3.py` is a backward-compatibility shim that re-exports the public API and is the console-script entry point (`main`).
- `tfbs_footprinter3/data/` - bundled JSON (JASPAR PWMs, background frequencies).
- `hpc/` - helper scripts for large runs: S3 publishing, CAS p-value / threshold building, and SLURM / CSC Puhti submission wrappers.
- `conda-recipe/meta.yaml`, `pyproject.toml`, `setup.py` (thin shim) - packaging. Build steps are in `Build_TFBSFootprinter.txt`.
- `CHANGELOG.md` tracks releases; version lives in `tfbs_footprinter3/__init__.py`.

## Gotchas
- The distribution name (`TFBS_footprinting3` / conda `tfbs-footprinting3`) differs from the import package name `tfbs_footprinter3` (note the `-ter` spelling). The console script is `tfbs_footprinter3`.
- First run (or `-update`) downloads experimental data from S3 into a local `data/` directory; the tool needs network access and more than 6 GB of RAM.
- Analyses call the Ensembl REST API and are subject to its rate limits. Most events (queries, rate-limit warnings, unknown transcript IDs, download errors) are written to `TFBS_footprinter3.log` in the working directory rather than to stdout; check it first when troubleshooting.
- Since 0.1.0 the TF name in output is a `{tf_name}__{matrix_id}` composite (e.g. `ARNT__MA0004.1`) so TFs with multiple JASPAR matrices no longer collide.
- Parquet output (`-of parquet` / `-of slim-parquet`) requires the `[parquet]` extra; pyarrow is intentionally not a core dependency.
- The Docker image is deprecated and runs Python 2.7. Prefer the conda or PyPI install.
