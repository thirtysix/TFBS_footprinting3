# Changelog

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
