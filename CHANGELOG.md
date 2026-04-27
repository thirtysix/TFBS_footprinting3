# Changelog

## 0.1.3 — 2026-04-27

Hotfix that pins the Ensembl REST endpoint, recovering 17 species
that silently broke when Ensembl rolled forward to assembly versions
newer than release 113.

- `tfbs_footprinter3/ensembl.py` now defaults its REST server to
  `https://oct2024.rest.ensembl.org` (the Ensembl 113 archive) instead
  of `http://rest.ensembl.org` (always-latest). The shipped species
  data tarballs and the JASPAR 2026 threshold artifacts are all keyed
  on Ensembl 113 transcript IDs, but the live REST endpoint had moved
  on to newer assemblies (e.g. ovis_aries went from `ARS-UI_Ramb_v2.0`
  to `v3.0`; the inbred-mouse `MGP_*` IDs and `ENSFCAT*` cat IDs are
  absent on the live endpoint as a result). Pinning fixes:
  - sheep (`ovis_aries`)
  - cat (`felis_catus`) — also resolves the species-name remap that
    previously sent the tool looking for `felis_catus_abyssinian`
  - 16 mouse strain genomes (the `MGP_*` lab strains plus `mus_spretus`)
- Override the endpoint with `TFBS_FOOTPRINTER3_ENSEMBL_REST=<url>`
  if you need a different release.
- 18 new species' tarballs uploaded to S3 with JASPAR 2026 thresholds
  (sheep + 17 recovered above), bringing the runtime-supported species
  count to 314 / 316 vertebrates.

## 0.1.2 — 2026-04-21

Hotfix for a plotting crash reported by users after 0.1.1.

- `plot_promoter` crashed with `ValueError: high <= 0` when the figure
  contained more than ~20 unique TFs (common with composite JASPAR
  2026 TF names, since each motif variant is a distinct `tf_name`).
  Two bugs compounded: the color palette was drained via
  `color_series.remove(...)` on every new TF, and the pick call
  `numpy.random.randint(0, len(color_series) - 1)` errored once the
  palette reached size 1 AND excluded the last color entry even when
  it didn't error. Replaced with deterministic modulo cycling
  (`color_series[len(color_dict) % len(color_series)]`) -- stable
  colors across re-runs, no palette exhaustion, no crash on any TF
  count.

## 0.1.1 — 2026-04-20

Hotfix on top of 0.1.0.

- The 0.1.0 conda package listed `wget` in run requirements, but
  conda-forge's `wget` is the GNU CLI binary -- not the Python `wget`
  package required by `data_loader.py`. Fresh conda installs crashed
  on first import with `ModuleNotFoundError: No module named 'wget'`.
- Fixed by replacing the `wget.download` calls in `data_loader.py`
  with `urllib.request.urlretrieve` (stdlib, no new dep). The Python
  `wget` package dependency is now dropped from both pyproject.toml
  and the conda recipe. PyPI installs were unaffected (they got the
  Python wget from pip) but are also now slimmer.

## 0.1.0 — 2026-04-20

Major release: JASPAR 2026 motif catalog, new threshold table format,
several long-standing scoring bugs fixed.

### JASPAR 2026 motif upgrade

- Replaces the 575-motif JASPAR 2018 catalog with the 1019-motif
  JASPAR 2026 CORE vertebrates non-redundant set (`data/JASPAR_2026_pwms.json`).
- Composite TF naming `{tf_name}__{matrix_id}` (e.g. `ARNT__MA0004.1`)
  disambiguates TFs that have multiple JASPAR matrices; same naming
  threads through all downstream output columns and data files.
- Per-species S3 tarballs now carry matching `*.CAS_thresholds.jaspar_2026.tsv.gz`
  and `*.tfs_thresholds.jaspar_2026.tsv.gz`.

### Bug fixes (pre-existing, now resolved)

- **GERP conservation was silently zero on every hit**. A motif-center
  calculation had an operator-precedence bug (`start + end/2` instead
  of `(start + end)/2`), which made the overlap check with GERP/ATAC
  ranges always miss. Fix restores ~16 % nonzero conservation scores
  on forward-strand transcripts.
- **Reverse-strand GERP/ATAC was doubly broken**: the translator
  emitted `(start, end)` pairs with `start > end` after TSS-relative
  conversion, which made every overlap check on reverse-strand genes
  fail silently. Fixed so all reverse-strand transcripts now contribute
  conservation properly.
- **PWM p-value output sentinels**: the boundary logic was mixing
  "query below empirical minimum" and "query above empirical maximum"
  into the same `">1.0"` output. Now returns `">{max_p}"` for below-min
  (legitimately non-significant) and `"<{min_p}"` for above-max (more
  significant than anything observed in the sample).

### Output format changes

- CAS p-value table switched from a ~300 MB per-species JSON
  (`CAS_pvalues.0.1.tf_ls.json`) to a ~800 KB gzipped TSV
  (`{species}.CAS_thresholds.jaspar_2026.tsv.gz`). The old JSON files
  are preserved inside the new tarballs for backward compatibility
  with 0.0.7 and earlier.
- Output tables now sort by CAS score (descending) instead of by
  CAS p-value string. Sorts are numerically correct for scientific
  notation and for the new `<p`/`>p` sentinel strings.
- Matplotlib is lazily imported, so `tfbs_footprinter3 --help`
  no longer pays the ~1.5 s matplotlib startup cost.

### HPC pipeline

New `hpc/` subtree for reproducing the per-species threshold tables:

- `hpc/ensembl_gtf.py` — pick N canonical protein-coding transcripts
  per species from a cached Ensembl GTF.
- `hpc/puhti/array_submit.sh` — CSC Puhti SLURM array, 10 transcripts
  per task via subprocess-per-transcript so memory doesn't accumulate.
- `hpc/puhti/aggregate_submit.sh` — chained aggregation job, fits in
  8 GB / 4 h for species up to 1000 transcripts.
- `hpc/puhti/build_cas_distributions.py` — streaming histogram
  aggregator (int64-encoded `(tf, score)` keys); memory stays flat at
  ~4 GB regardless of hit count.
- `hpc/build_release_tarballs.sh` — merges new threshold TSVs into
  each species' existing S3 tarball for release.

### Performance

- `-of slim-parquet` output format for the HPC campaign: writes only
  `(binding prot., PWM score, combined affinity score)`. Peak memory
  ~4 GB on a 20 M-hit transcript (vs ~20 GB for full parquet). Users
  who want the full contextual columns continue to use `-of parquet`.
- Vectorized scoring path (batched weight components per TF) replaces
  a scalar inner loop.

### Breaking changes

- TF names in all output tables are now composite keys
  (`ARNT__MA0004.1` rather than `ARNT`). Downstream code that joins on
  the `binding prot.` column needs to account for the `__MA...` suffix.
- `calcCombinedAffinityPvalue` signature changed to accept a single
  per-TF slice `{"scores": ndarray, "pvalues": ndarray}` rather than
  the old multi-argument decomposition. External callers are unlikely
  but would need updating.

## 0.0.7

Prior baseline. See git history for fine-grained changes.
