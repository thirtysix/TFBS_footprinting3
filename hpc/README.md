# HPC pipeline for non-human CAS distributions

Scripts that build per-species empirical Combined Affinity Score p-value
distributions (`CAS_pvalues.0.1.tf_ls.json`) for the 123 non-human species
supported by tfbs_footprinter3.

## Scope

`PWM + GERP + CAGE + CpG` — the four components already computed by
`find_clusters()` for non-human species at
`tfbs_footprinter3/tfbs_footprinter3.py:1015`.

## Stages

| Stage | Script | Runs where | Output |
|-------|--------|------------|--------|
| A | `select_transcripts.py` | workstation, one-time | `transcript_lists/{species}.txt` |
| B | `prefetch.py` + `ensembl_cache.py` | workstation, one-time | `cache/ensembl.sqlite` |
| C | `cas_only.py` via `slurm/submit.sh` | SLURM job array | `cas_scores/{species}/*.parquet` |
| D | `build_cas_pvalues.py` | workstation | `artifacts/{species}/CAS_pvalues.0.1.tf_ls.json` |
| E | S3 upload (todo) | workstation | `s3://tfbssexperimentaldata/{species}.tar.gz` |

## Stages B–D usage

```bash
# Stage B: populate the Ensembl REST cache (single-threaded on submission node)
python hpc/prefetch.py --list hpc/transcript_lists/mus_musculus.txt \
    --cache hpc/cache/ensembl.sqlite

# Stage C benchmark: measure per-transcript time before scheduling full run
python -m hpc.cas_only \
    --list hpc/transcript_lists/mus_musculus.txt --range 0:5 \
    --cache hpc/cache/ensembl.sqlite --benchmark

# Stage C full run: SLURM array of ~500 tasks per species
export SPECIES=mus_musculus
export TOTAL_TRANSCRIPTS=5000
export CHUNK_SIZE=10
sbatch --array=0-$((TOTAL_TRANSCRIPTS / CHUNK_SIZE - 1)) hpc/slurm/submit.sh

# Stage D: per-species aggregation into CAS_pvalues JSON
python hpc/build_cas_pvalues.py --species mus_musculus
# or to do every species present under hpc/cas_scores/
python hpc/build_cas_pvalues.py --species all
```

## Key design decisions

- **CAS scope**: PWM + GERP + CAGE + CpG — the four components
  `find_clusters()` already computes for non-human species
  (`tfbs_footprinter3.py:1055`). No scoring code changes.
- **Cache-first REST**: `ensembl_cache.py` monkey-patches
  `tfbs_footprinter3.ensemblrest` to read from a sqlite file populated in
  Stage B. SLURM workers open the cache read-only by default, so N
  parallel tasks never re-fetch from Ensembl. `--allow-live-rest` opts
  a worker into live-fetch fallback on cache miss.
- **Rounding**: p-values are computed on scores rounded to 2 decimals
  to match the rounding applied at CAS computation time
  (`tfbs_footprinter3.py:1073`).
- **Min-hits threshold**: TFs with fewer than 100 hits across the 5000
  transcripts are dropped from the p-value table (configurable via
  `--min-hits`). Too-sparse distributions are misleading.

## Stage A — transcript selection

```bash
# All species, 5000 canonical protein-coding transcripts each (default)
python hpc/select_transcripts.py

# Subset
python hpc/select_transcripts.py --species mus_musculus danio_rerio

# Resume after partial failures (BioMart occasionally returns transient db errors)
python hpc/select_transcripts.py --resume

# Dry run
python hpc/select_transcripts.py --dry-run
```

Output: `transcript_lists/{species}.txt` (gitignored) + `manifest.json`
recording per-species success/fail and actual counts (flags species with
<5000 canonical protein-coding transcripts as `reduced`).

## Directory layout

```
hpc/
├── README.md                    (this file)
├── select_transcripts.py        (Stage A)
├── ensembl_cache.py             (shared cache client + monkey-patch helper)
├── prefetch.py                  (Stage B)
├── cas_only.py                  (Stage C worker)
├── build_cas_pvalues.py         (Stage D)
├── slurm/
│   └── submit.sh                (Stage C SLURM array wrapper)
├── transcript_lists/            (gitignored; Stage A output)
├── cache/                       (gitignored; Stage B output)
├── cas_scores/                  (gitignored; Stage C output)
└── artifacts/                   (gitignored; Stage D output)
```
