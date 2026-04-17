# CSC Puhti pipeline for the non-human CAS campaign

End-to-end recipe for running the CAS empirical-p-value build on
`barker@puhti.csc.fi` (CSC Puhti). Pilot scope: one species,
100 canonical protein-coding transcripts, promoter window
`-pb 5000 -pa 5000`, `pvalc=1` (unfiltered), parquet output.

**JASPAR version**: this campaign targets JASPAR 2026 (1019 vertebrates
non-redundant PWMs). The tool on the `feat/jaspar-2026` branch ships
`JASPAR_2026_pwms.json` and uses composite TF keys like `ARNT__MA0004.1`
throughout. The existing S3 species tarballs are still keyed on 2018 TF
names — the purpose of this Puhti campaign is to regenerate that data
against JASPAR 2026 motifs, after which the tarballs get re-uploaded.

Project dir: `/scratch/project_2001307/ensembl_genomes_CAS_scoring_2026`
(referred to as `$PROJECT` below).

Workflow is **workstation → Puhti → workstation**; the GTF parsing and
final aggregation happen on the workstation (small data, quick loops),
the per-transcript compute runs on Puhti.

## One-time setup on Puhti

```bash
ssh barker@puhti.csc.fi
cd /scratch/project_2001307/ensembl_genomes_CAS_scoring_2026

# Clone the tool on the JASPAR 2026 branch (1019 vertebrates non-redundant PWMs)
git clone -b feat/jaspar-2026 https://github.com/thirtysix/TFBS_footprinting3.git
cd TFBS_footprinting3

# Load a suitable Python (adjust to Puhti's current selection)
module load python-data/3.10

# Virtualenv with parquet extras
python -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -e '.[parquet]'

# Sanity check
tfbs_footprinter3 --help | head -2
```

## Per-species workflow

### 1. Select transcripts (on workstation)

```bash
# In the repo on your workstation:
python hpc/ensembl_gtf.py \
    --species acanthochromis_polyacanthus \
    --n 100 --max-chromosomes 10

# Produces: hpc/transcript_lists/acanthochromis_polyacanthus.txt
```

### 2. Copy the transcript list to Puhti

```bash
SPECIES=acanthochromis_polyacanthus
PROJECT=/scratch/project_2001307/ensembl_genomes_CAS_scoring_2026

rsync hpc/transcript_lists/${SPECIES}.txt \
    barker@puhti.csc.fi:${PROJECT}/runs/${SPECIES}/transcripts.txt
```

### 3. Submit the SLURM job on Puhti

```bash
ssh barker@puhti.csc.fi
cd /scratch/project_2001307/ensembl_genomes_CAS_scoring_2026

# One-time: create the per-species run dir
mkdir -p runs/${SPECIES}

# Launch
sbatch --export=SPECIES=${SPECIES} \
    TFBS_footprinting3/hpc/puhti/pilot_submit.sh
```

**Edit `pilot_submit.sh`'s `#SBATCH` header first** for the billing
project and partition on your account. The defaults in this repo are
placeholders.

The job:

* creates `$PROJECT/runs/$SPECIES/` as the working directory
* runs `tfbs_footprinter3` on the 100 transcripts with
  `-pb 5000 -pa 5000 -p 1 -pc 1 -no -of parquet`
* experimental data for the species is auto-downloaded on first run
  from `s3://tfbssexperimentaldata/${SPECIES}.tar.gz` — no manual copy
  required

### 4. Aggregate results (on workstation)

```bash
# Pull the parquet shards back
rsync -a --include='*.parquet' --include='*/' --exclude='*' \
    barker@puhti.csc.fi:${PROJECT}/runs/${SPECIES}/tfbs_results/ \
    hpc/cas_scores/${SPECIES}/

# Roll up into the CAS p-value JSON + thresholds TSV
python hpc/build_cas_pvalues.py --species ${SPECIES}
```

Produces `hpc/artifacts/${SPECIES}/CAS_pvalues.0.1.tf_ls.json` in the
exact format consumed by `tfbs_footprinter3.py:755`
(`species_specific_data` loads it on startup for that species).

## Scaling to all 316 species

Once the pilot for one species validates end-to-end, the same recipe
drives all 316 — just run step 1 in batch mode (`--species-file
/path/to/species_list.txt`) and step 3 once per species (or wrap in
a SLURM job array that iterates over species).

Projected cost (based on local benchmarking, with pvalc=1, 10kb window,
parquet output, Ensembl prefetch cache):

| | per species | all 316 |
|---|---|---|
| CPU-time | ~1 h (100 transcripts × 35 s) | ~316 h |
| Wall @ 10-way SLURM array per species | ~6 min | hours |
| Disk (parquet) | ~5 GB | ~1.6 TB |
| CAS_pvalues JSON uploaded to S3 | ~MB | ~GB total |

## Notes

- The SLURM submit script expects `SPECIES=<slug>` in the env. Without
  it the job aborts with a clear error.
- The tool auto-downloads per-species tarballs from S3 on first run.
  On compute nodes this requires outbound HTTP access, which Puhti's
  standard partitions allow. If a partition blocks outbound, stage the
  tarball manually:
  `aws s3 cp s3://tfbssexperimentaldata/${SPECIES}.tar.gz ${PROJECT}/TFBS_footprinting3/tfbs_footprinter3/data/`.
- `-p 1 -pc 1` keeps every hit (no filtering). At 100 transcripts × 10kb
  window this produces ~1M hits per TF per species, the statistical
  floor for stable empirical p-values to ~p=1e-4.
