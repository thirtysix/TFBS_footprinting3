#!/bin/bash -l
# CSC Puhti SLURM pilot submission for the non-human CAS build.
#
# The `-l` (login shell) is required so /etc/profile.d/* sources Lmod; the
# default non-login SLURM shell on Puhti has no `module` command otherwise.
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

# Ensure Lmod / `module` is available. sbatch on Puhti can spawn a shell
# where `module` is not a defined function, even with -l. Try the common
# init locations; the first one present wins. The `set +u` wrapper is
# because CSC's zz-csc-env.sh references LC_CTYPE (and possibly other
# unset vars) that would otherwise abort us under `set -u`.
if ! command -v module >/dev/null 2>&1; then
    set +u
    for _lmod_init in \
        /appl/profile/zz-csc-env.sh \
        /etc/profile.d/lmod.sh \
        /etc/profile.d/z00_lmod.sh \
        /usr/share/lmod/lmod/init/bash; do
        if [[ -f "${_lmod_init}" ]]; then
            # shellcheck disable=SC1090
            source "${_lmod_init}"
            break
        fi
    done
    set -u
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
