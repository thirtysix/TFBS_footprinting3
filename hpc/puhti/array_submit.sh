#!/bin/bash -l
# Per-transcript CSC Puhti SLURM submission for the CAS campaign.
#
# One SLURM task per transcript. Each task runs tfbs_footprinter3 on a
# single transcript with the campaign's standard settings:
#
#   -pb 5000 -pa 5000 -p 1 -pc 1 -no -of parquet
#
# Why per-transcript instead of one big serial job:
#
#   - Memory doesn't accumulate. Serial runs climb to 30+ GB over 20-odd
#     transcripts (observed: job 34055541 OOM'd at 102.5% of 32 GB after
#     23 transcripts). Per-transcript peak is ~7-9 GB, well inside the
#     12 GB allocation below.
#
#   - True wall-clock parallelism. With 20 concurrent tasks a 100-
#     transcript species completes in ~30 minutes instead of ~11 hours.
#
#   - Failure isolation. A bad transcript only kills that one task;
#     re-submit the failed array indices without losing the others.
#
# Usage:
#
#   # Default: 100 transcripts (1..100 in transcripts.txt)
#   sbatch --export=SPECIES=acanthochromis_polyacanthus \
#       TFBS_footprinting3/hpc/puhti/array_submit.sh
#
#   # Different species count: override --array on the command line
#   sbatch --array=1-316%20 --export=SPECIES=homo_sapiens \
#       TFBS_footprinting3/hpc/puhti/array_submit.sh
#
# The %20 suffix caps concurrency at 20 to be a polite scheduler citizen;
# raise to %50 if you want more throughput and the queue is empty.

#SBATCH --job-name=tfbs_cas
#SBATCH --account=project_2001307
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:45:00
#SBATCH --array=1-100%20
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -eo pipefail

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
fi
set -u

: "${SPECIES:?SPECIES env var is required (e.g. acanthochromis_polyacanthus)}"
: "${PROJECT:=/scratch/project_2001307/ensembl_genomes_CAS_scoring_2026}"
: "${REPO_ROOT:=${PROJECT}/TFBS_footprinting3}"
: "${RUN_DIR:=${PROJECT}/runs/${SPECIES}}"

# Require the transcripts file (one-time rsync from workstation)
if [[ ! -s "${RUN_DIR}/transcripts.txt" ]]; then
    echo "FATAL: ${RUN_DIR}/transcripts.txt is missing or empty." >&2
    exit 1
fi

# Pick this task's transcript ID (sed is 1-indexed, SLURM array starts at 1 by our default)
TRANSCRIPT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${RUN_DIR}/transcripts.txt")
if [[ -z "${TRANSCRIPT}" ]]; then
    echo "FATAL: no transcript on line ${SLURM_ARRAY_TASK_ID} of ${RUN_DIR}/transcripts.txt" >&2
    exit 1
fi

# Per-task working dir keeps concurrent TFBS_footprinter3.log writes
# from interleaving, and gives each task a clean cwd for tfbs_results/.
TASK_DIR="${RUN_DIR}/tasks/$(printf '%03d' "${SLURM_ARRAY_TASK_ID}")"
mkdir -p "${TASK_DIR}"
cd "${TASK_DIR}"
echo "${TRANSCRIPT}" > transcripts.txt

echo "[$(date -Iseconds)] array_id=${SLURM_ARRAY_TASK_ID} species=${SPECIES} transcript=${TRANSCRIPT}"

module load python-data/3.10
source "${REPO_ROOT}/venv/bin/activate"

tfbs_footprinter3 \
    -t transcripts.txt \
    -pb 5000 -pa 5000 \
    -p 1 -pc 1 \
    -no \
    -of parquet

echo "[$(date -Iseconds)] done"
ls tfbs_results/ | head
