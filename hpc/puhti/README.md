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

# rsync 3.1.x will not auto-create missing parent dirs, so create the
# per-species run dir on Puhti first. (The -p in `mkdir -p` means this
# is also safe to re-run across species and creates $PROJECT itself
# on first use.)
ssh barker@puhti.csc.fi "mkdir -p ${PROJECT}/runs/${SPECIES}"

rsync hpc/transcript_lists/${SPECIES}.txt \
    barker@puhti.csc.fi:${PROJECT}/runs/${SPECIES}/transcripts.txt
```

### 3. Submit the SLURM array + chained aggregator on Puhti

```bash
ssh barker@puhti.csc.fi
cd /scratch/project_2001307/ensembl_genomes_CAS_scoring_2026

# One task per transcript (20-way concurrent, 8 GB each) plus an
# aggregation job that fires automatically once the array succeeds.
# Only the ~800 KB threshold TSV needs to leave Puhti afterwards.
ARRAY_ID=$(sbatch --parsable --export=SPECIES=${SPECIES} \
    TFBS_footprinting3/hpc/puhti/array_submit.sh)
sbatch --dependency=afterok:${ARRAY_ID} --export=SPECIES=${SPECIES} \
    TFBS_footprinting3/hpc/puhti/aggregate_submit.sh
```

**Why a batched array**: memory grows with each transcript processed
in a single Python process (observed: 32.8 GB after 23 transcripts
in-process). Two mitigations stacked:

  1. Each array task processes N transcripts (default
     `TRANSCRIPTS_PER_TASK=10`), invoking `tfbs_footprinter3` as a
     fresh subprocess per transcript so memory fully releases between
     them. Per-transcript peak stays at ~4-5 GB under slim-parquet.
  2. The array size drops from 100 to 10 tasks per species. At CSC's
     ~200-slot AssocMaxSubmit quota, the campaign can run ~20
     species concurrently (10 tasks × 20 species = 200) instead of
     ~2 — a ~5x speedup for the full 316-species run.

The array:

* creates `$PROJECT/runs/$SPECIES/tasks/NNN/` per task
* picks 10 consecutive transcripts from `transcripts.txt` (task 1:
  lines 1-10; task 2: lines 11-20; …)
* runs `tfbs_footprinter3 -t single_transcript.txt -pb 5000 -pa 5000 -p 1 -pc 1 -no -of slim-parquet`
  once per transcript in a bash loop
* species experimental data is pre-staged by submit_batch.sh on the
  Puhti login node (one wget + tar-extract) before the array fans
  out, eliminating the concurrent-download race that would otherwise
  corrupt in-flight msgpack files

Defaults are `--array=1-10%10`, `--mem=8G`, `--time=01:30:00`.
Override for a different layout — e.g. 5 transcripts per task:
`sbatch --array=1-20%10 --export=SPECIES=X,TRANSCRIPTS_PER_TASK=5 …`.

For one-off single-transcript tests, `pilot_submit.sh` still runs the
whole list serially; kept as a reference but not used for the campaign.

### 4. Pull the threshold TSV back to the workstation

Aggregation happens on Puhti (step 3's second sbatch); only the small
artifact needs to come back:

```bash
mkdir -p hpc/artifacts/${SPECIES}
rsync barker@puhti.csc.fi:${PROJECT}/artifacts/${SPECIES}/${SPECIES}.CAS_thresholds.jaspar_2026.tsv.gz \
    hpc/artifacts/${SPECIES}/
```

That TSV is the single runtime artifact the tool consumes. It can be
dropped into a species' S3 tarball (`tfbs_footprinter3/data/${SPECIES}/`
at install time) and `data_loader.py:306` will load it automatically.

For QA / manual inspection of the raw parquet shards, the older
"pull everything back" flow still works:

```bash
mkdir -p hpc/cas_scores/${SPECIES}
rsync -a --include='*.parquet' --include='*/' --exclude='*' \
    barker@puhti.csc.fi:${PROJECT}/runs/${SPECIES}/tasks/ \
    hpc/cas_scores/${SPECIES}/
python hpc/inspect_species_results.py \
    --results-dir hpc/cas_scores/${SPECIES} --expected-transcripts 100
```

## Batch mode (multiple species at once)

A small wrapper sbatches array + aggregate pairs for every species in a
list. Intended for scaling up: first a 10-species test, then the full
316 once that's validated.

1. **Workstation**: generate all transcript lists in one call
   ```bash
   python hpc/ensembl_gtf.py --species-file hpc/species_lists/test10.txt \
       --n 100 --max-chromosomes 10 --resume
   ```

2. **Workstation → Puhti**: push each species' transcripts.txt into its
   run dir (both sides get created on first use). `ssh -n` is required
   inside the loop — without it, ssh inherits the loop's stdin and
   swallows the rest of the species list after the first iteration.
   ```bash
   while IFS= read -r SP; do
       [[ -z "${SP}" || "${SP:0:1}" == "#" ]] && continue
       ssh -n barker@puhti.csc.fi "mkdir -p ${PROJECT}/runs/${SP}"
       rsync hpc/transcript_lists/${SP}.txt \
           barker@puhti.csc.fi:${PROJECT}/runs/${SP}/transcripts.txt
   done < hpc/species_lists/test10.txt
   ```

3. **Puhti**: submit all species' paired jobs via the batch wrapper
   ```bash
   ssh barker@puhti.csc.fi
   cd /scratch/project_2001307/ensembl_genomes_CAS_scoring_2026
   cd TFBS_footprinting3 && git pull && cd ..
   bash TFBS_footprinting3/hpc/puhti/submit_batch.sh \
       TFBS_footprinting3/hpc/species_lists/test10.txt
   ```
   Each species gets its own array + chained aggregate pair. With
   `%20` in-array concurrency × N species, peak concurrent tasks is
   `20 × N` — fine on Puhti's `small` partition for the 10-species
   test; drop the `%20` in array_submit.sh if you want to be gentler
   for the 316 run.

4. **Monitor**: `squeue -u $USER`. A completed species produces
   `${PROJECT}/artifacts/${SPECIES}/${SPECIES}.CAS_thresholds.jaspar_2026.tsv.gz`.

5. **Pull artifacts back (single rsync for everything)**
   ```bash
   mkdir -p hpc/artifacts && rsync -a \
       barker@puhti.csc.fi:${PROJECT}/artifacts/ hpc/artifacts/
   ```

Projected cost for the full 316-species run, assuming this 10-species
test confirms ~6-min/task wall and ~4-GB/task peak:

| | per species | all 316 |
|---|---|---|
| CPU-time | ~10 h (100 tasks × ~6 min) | ~3,200 h |
| Wall @ 20-way in-array concurrency | ~30 min | ~1-2 days (many species in parallel) |
| Puhti scratch (parquet + logs) | ~9 GB | ~2.9 TB |
| **Artifact back to workstation** | **~800 KB TSV** | **~250 MB total** |

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
