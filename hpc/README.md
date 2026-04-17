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
| B | `prefetch.py` (todo) | workstation, one-time | `cache/ensembl.sqlite` |
| C | `cas_only.py` (todo) | SLURM job array | `cas_scores/{species}/*.parquet` |
| D | `build_cas_pvalues.py` (todo) | workstation | `artifacts/{species}/CAS_pvalues.0.1.tf_ls.json` |
| E | S3 upload (todo) | workstation | `s3://tfbssexperimentaldata/{species}.tar.gz` |

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
├── prefetch.py                  (Stage B, todo)
├── cas_only.py                  (Stage C, todo)
├── build_cas_pvalues.py         (Stage D, todo)
├── slurm/                       (SLURM submit scripts)
│   └── submit.sh                (todo)
├── transcript_lists/            (gitignored; Stage A output)
├── cache/                       (gitignored; Stage B output)
├── cas_scores/                  (gitignored; Stage C output)
└── artifacts/                   (gitignored; Stage D output)
```
