#!/bin/bash
# CSC Puhti SLURM pilot submission for the non-human CAS build.
#
# Runs tfbs_footprinter3 on 100 canonical protein-coding transcripts of
# a single species, promoter window -5kb to +5kb, pvalc=1, parquet
# output. Intended for one pilot species end-to-end before fanning out
# to the full 316-species campaign.
#
# Usage:
#   sbatch --export=SPECIES=<species_slug> hpc/puhti/pilot_submit.sh
#
# EDIT the #SBATCH lines below for your actual account/partition. The
# defaults below are placeholders that should work on Puhti's `small`
# partition for CSC project 2001307.

#SBATCH --job-name=tfbs_cas_pilot
#SBATCH --account=project_2001307
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:30:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

: "${SPECIES:?SPECIES env var is required (e.g. acanthochromis_polyacanthus)}"
: "${PROJECT:=/scratch/project_2001307/ensembl_genomes_CAS_scoring_2026}"
: "${REPO_ROOT:=${PROJECT}/TFBS_footprinting3}"
: "${RUN_DIR:=${PROJECT}/runs/${SPECIES}}"

echo "[$(date -Iseconds)] species=${SPECIES}"
echo "[$(date -Iseconds)] run_dir=${RUN_DIR}"
mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"

# Require the transcripts file to have been rsync'd ahead of time.
if [[ ! -s transcripts.txt ]]; then
    echo "FATAL: ${RUN_DIR}/transcripts.txt is missing or empty." >&2
    echo "Run hpc/ensembl_gtf.py on the workstation and rsync the result first." >&2
    exit 1
fi

# Python env set up by hpc/puhti/README.md one-time setup
module load python-data/3.10
source "${REPO_ROOT}/venv/bin/activate"

echo "[$(date -Iseconds)] tool version:"
tfbs_footprinter3 --help | head -1

# The heavy lift. -pb/-pa 5000: 10 kb window. -p -pc 1: keep every hit.
# -no: no figure (we only need the sortedclusters table). -of parquet:
# compact binary output, ~10x smaller than CSV.
echo "[$(date -Iseconds)] launching tfbs_footprinter3"
tfbs_footprinter3 \
    -t transcripts.txt \
    -pb 5000 -pa 5000 \
    -p 1 -pc 1 \
    -no \
    -of parquet

echo "[$(date -Iseconds)] done"
ls tfbs_results/ | head
echo "[$(date -Iseconds)] parquet file count: $(find tfbs_results -name '*.parquet' | wc -l)"
