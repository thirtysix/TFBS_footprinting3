# Changelog

## 0.1.5 — 2026-09-05

Fixes a silent regression that disabled the eQTL evidence layer entirely.

- `data_loader.species_specific_data` built the per-chromosome GTEx
  filename from a hard-coded `"gtex_v7"` release string, but every
  shipped tarball has carried `gtex_v8` since the data was rebuilt
  (verified: `homo_sapiens.tar.gz` contains 24 `gtex_v8` files and no
  `gtex_v7`). The lookup was a bare `os.path.exists` with no `else`, so
  the mismatch did not raise -- `gtex_variants` stayed empty and the
  `eqtls weights sum` column was identically zero for every hit. One of
  the seven advertised evidence layers had been contributing nothing,
  undetectably. Confirmed against the committed sample output: 0 of
  ~500,000 rows had a nonzero eQTL weight.
- Replaced with `data_loader.find_chromosome_data_file()`, which matches
  on species + `.Chr<N>.` + a stable token and takes the
  lexicographically last hit, so a future release bump cannot silently
  disable the layer the same way. The chromosome is matched with
  surrounding dots so `Chr1` does not also match `Chr10`..`Chr19`.
- New `tests/test_eqtl_data_discovery.py` pins the behaviour, including
  an end-to-end assertion against the installed human data.

Reviving that layer immediately exposed a second, latent crash in the
same code path:

- `_eqtl_batch` (and its scalar counterpart `eqtls_weights_summing`)
  looked the eQTL effect magnitude up in `gtex_weights_dict` by exact
  float key. That table is a 0.001-resolution grid -- 2,617 keys from
  0.0 to 6.023 in the shipped human data -- while the magnitude stored
  on each variant carries full precision, so any value not landing on a
  grid point raised `KeyError` (observed: 0.392323, 0.388412, 0.116223).
  Unreachable before, because `converted_eqtls` was always empty.
- Now bisects onto the grid, the same resolution `cpg_weights_summing`
  already uses for the identical problem, clamping above the top of the
  table and taking the absolute magnitude. Pinned by
  `tests/test_eqtl_weight_lookup.py`.
- `eqtls_weights_summing_v`, the scalar reference that
  `tests/test_scoring_vectorized.py` checks `_eqtl_batch` against, was
  updated the same way. Left as it was, the agreement test would have
  compared a snapping implementation against a non-snapping one and
  quietly stopped testing anything.

### Expression-correlation layer was dead since 0.1.0

The composite TF naming introduced with JASPAR 2026 broke the lookup into
the CAGE expression-correlation data. Scoring carries
`{tf_name}__{matrix_id}` (upper-cased by `cli.py`), but
`{species}.CAGE.jasparTFs.dict.json` was never rebuilt and still holds the
575 *bare* JASPAR 2018 names, so `if tf_name in TF_cage_dict` matched
**0 of 1019** motifs and `corr. weight sum` was identically zero for every
hit -- a second one of the seven advertised evidence layers contributing
nothing, undetectably.

- `resolve_tf_cage_key()` now strips the matrix id and falls back to a
  case-insensitive match, recovering **597 of 1019** motifs (the remainder
  are TFs the 2018-era correlation data does not cover at all).
- Pinned by `tests/test_tf_cage_key_resolution.py`, including the dimer
  case (`NR1H3::RXRA__MA0074.1`), where `::` must not be confused with the
  `__` matrix-id separator.
- Rebuilding that dictionary for the full JASPAR 2026 catalog remains the
  real fix; this restores the coverage the data can support.

### Metacluster weights: a silent gap for new motif widths

The GTRD metacluster weight table is keyed by motif length and was built
against JASPAR 2018, whose 575 motifs it covered completely. JASPAR 2026
introduced widths the table does not carry (4, 22, 23 and 30 nt, covering
5 of 1019 motifs). The batch scoring path returned 0 for those without
comment, and the scalar path raised `KeyError` on the same input.

- Both paths now return 0 and log a warning once per unseen motif length,
  so the gap is discoverable rather than silent, and the two paths agree.
- This does not fix the underlying data gap: those 5 motifs still score 0
  for one of the strongest evidence layers. Rebuilding the table over the
  JASPAR 2026 width range is the real fix.

Note for anyone comparing against published results: the eQTL term was
zero in the version used for the benchmarks in
doi:10.1080/21541264.2025.2521764, so eQTL feature-contribution figures
from that paper describe a layer that was not loading.


## 0.1.4 — 2026-05-22

Hotfix for a false-negative connectivity check that locked out users
behind Google-blocking firewalls.

- `io_utils.is_online()` previously probed `www.google.com:80`, which
  has nothing to do with what the tool actually downloads. Users on
  networks where Google is unreachable but the S3 experimental-data
  bucket and the Ensembl REST endpoint are both reachable (e.g. a
  CentOS 7 user in mainland China who reported `curl` against the
  bucket returning `200 OK` while `tfbs_footprinter3 --version`
  exited with "System does not appear to be connected to the
  internet") were unable to run the tool at all.
- The new probe targets the hosts the tool actually depends on
  (`s3.us-east-2.amazonaws.com` + the pinned Ensembl REST endpoint),
  uses HTTPS/443 (port 80 is increasingly firewalled), and returns
  True if *any* target is reachable. `TFBS_FOOTPRINTER3_ENSEMBL_REST`
  is honored if set, so users on networks that need an Ensembl mirror
  can point the probe at the same URL.

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
