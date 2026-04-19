#!/bin/bash -l
# Batched CSC Puhti SLURM submission for the CAS campaign.
#
# Each array task processes N consecutive transcripts
# (TRANSCRIPTS_PER_TASK, default 10) from transcripts.txt. Within the
# task, transcripts are run one at a time via separate tfbs_footprinter3
# subprocess invocations -- that guarantees memory is fully reclaimed
# by the OS between transcripts (Python's allocator can hold onto
# arenas across iterations in-process), keeping each task's peak RSS
# at the single-transcript level (~4-5 GB under -of slim-parquet).
#
# Why batch instead of one task per transcript:
#
#   - 10x fewer queue slots per species (10 tasks vs 100). At CSC's
#     ~200-slot user quota, the campaign can fan out to ~20 species
#     concurrently instead of ~2 -- a ~5x speedup for the full
#     316-species run.
#   - Failure isolation still works -- a bad transcript only kills one
#     subprocess; the bash loop continues. Re-running the task picks
#     up only the failed transcripts via the tool's results_files_exist
#     resume check.
#
# Usage:
#
#   # Default: 100 transcripts in 10 tasks of 10 transcripts each
#   sbatch --export=SPECIES=acanthochromis_polyacanthus \
#       TFBS_footprinting3/hpc/puhti/array_submit.sh
#
#   # Change batch size or total count
#   sbatch --export=SPECIES=X,TRANSCRIPTS_PER_TASK=5 --array=1-20%10 \
#       TFBS_footprinting3/hpc/puhti/array_submit.sh

#SBATCH --job-name=tfbs_cas
#SBATCH --account=project_2001307
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --array=1-10%10
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
: "${TRANSCRIPTS_PER_TASK:=10}"

# Require the transcripts file (one-time rsync from workstation)
if [[ ! -s "${RUN_DIR}/transcripts.txt" ]]; then
    echo "FATAL: ${RUN_DIR}/transcripts.txt is missing or empty." >&2
    exit 1
fi

# Pick this task's slice of transcripts (sed is 1-indexed)
START=$(( (SLURM_ARRAY_TASK_ID - 1) * TRANSCRIPTS_PER_TASK + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * TRANSCRIPTS_PER_TASK ))
TRANSCRIPTS=$(sed -n "${START},${END}p" "${RUN_DIR}/transcripts.txt")
if [[ -z "${TRANSCRIPTS}" ]]; then
    echo "FATAL: no transcripts on lines ${START}-${END} of ${RUN_DIR}/transcripts.txt" >&2
    exit 1
fi

# Per-task working dir keeps concurrent TFBS_footprinter3.log writes
# from interleaving, and gives each task a clean cwd for tfbs_results/.
# All N transcripts' output subdirs land under this one task dir.
TASK_DIR="${RUN_DIR}/tasks/$(printf '%03d' "${SLURM_ARRAY_TASK_ID}")"
mkdir -p "${TASK_DIR}"
cd "${TASK_DIR}"

echo "[$(date -Iseconds)] array_id=${SLURM_ARRAY_TASK_ID} species=${SPECIES} transcripts=$(echo "${TRANSCRIPTS}" | wc -l)"

module load python-data/3.10
source "${REPO_ROOT}/venv/bin/activate"

# Process transcripts ONE AT A TIME via separate subprocesses so each
# run starts with a clean memory footprint. If any single transcript
# fails, log it and continue with the rest -- the array task still
# succeeds as long as the bash loop exits cleanly. The resubmit
# strategy is to re-run the array; tfbs_footprinter3's
# results_files_exist check (cli.py:366) skips transcripts whose
# sortedclusters.parquet already exists, so only failures recompute.
FAILED_TRANSCRIPTS=()
SUCCESS_COUNT=0

while IFS= read -r TRANSCRIPT; do
    [[ -z "${TRANSCRIPT}" ]] && continue
    echo "[$(date -Iseconds)]   transcript ${TRANSCRIPT}"
    echo "${TRANSCRIPT}" > single_transcript.txt
    if tfbs_footprinter3 \
        -t single_transcript.txt \
        -pb 5000 -pa 5000 \
        -p 1 -pc 1 \
        -no \
        -of slim-parquet; then
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "  [warn] ${TRANSCRIPT} failed" >&2
        FAILED_TRANSCRIPTS+=("${TRANSCRIPT}")
    fi
done <<< "${TRANSCRIPTS}"

echo "[$(date -Iseconds)] done: ${SUCCESS_COUNT} succeeded, ${#FAILED_TRANSCRIPTS[@]} failed"
if [[ ${#FAILED_TRANSCRIPTS[@]} -gt 0 ]]; then
    printf '  failed: %s\n' "${FAILED_TRANSCRIPTS[@]}"
fi
ls tfbs_results/ 2>/dev/null | head
