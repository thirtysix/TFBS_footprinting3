#!/bin/bash -l
# Batch-submit paired (array + chained aggregate) jobs for a list of species.
# Runs on Puhti's login node; no compute is used here — it just sbatches.
#
# Usage (from Puhti login):
#   cd /scratch/project_2001307/ensembl_genomes_CAS_scoring_2026
#   bash TFBS_footprinting3/hpc/puhti/submit_batch.sh \
#       TFBS_footprinting3/hpc/species_lists/test10.txt
#
# For each species in the list:
#   * confirms $RUN_DIR/transcripts.txt exists (otherwise skips with a warning)
#   * sbatches array_submit.sh -> captures $ARRAY_ID
#   * sbatches aggregate_submit.sh with --dependency=afterok:$ARRAY_ID
#   * prints both IDs and the species name
#
# The 20-way concurrency cap inside array_submit.sh means 10 species x
# 20 concurrent = 200 tasks running simultaneously at peak. Drop the
# `%20` in array_submit.sh's #SBATCH --array line to a smaller number
# (or pass --array=1-100%N on the command line here) if you want to
# be gentler on the queue.

set -eo pipefail

SPECIES_LIST="${1:?usage: submit_batch.sh <species-list-file>}"

if [[ ! -s "${SPECIES_LIST}" ]]; then
    echo "FATAL: species list ${SPECIES_LIST} is missing or empty" >&2
    exit 1
fi

: "${PROJECT:=/scratch/project_2001307/ensembl_genomes_CAS_scoring_2026}"
: "${REPO_ROOT:=${PROJECT}/TFBS_footprinting3}"

declare -a SUBMITTED_ARRAYS=()
declare -a SUBMITTED_AGGS=()
declare -a SKIPPED=()

while IFS= read -r SPECIES; do
    # Tolerate blank lines and comments
    SPECIES="${SPECIES%%#*}"
    SPECIES="${SPECIES// /}"
    [[ -z "${SPECIES}" ]] && continue

    RUN_DIR="${PROJECT}/runs/${SPECIES}"
    if [[ ! -s "${RUN_DIR}/transcripts.txt" ]]; then
        echo "  [skip] ${SPECIES}: ${RUN_DIR}/transcripts.txt missing (rsync first)"
        SKIPPED+=("${SPECIES}")
        continue
    fi

    ARRAY_ID=$(sbatch --parsable \
        --export=SPECIES="${SPECIES}" \
        "${REPO_ROOT}/hpc/puhti/array_submit.sh")
    AGG_ID=$(sbatch --parsable \
        --dependency=afterok:"${ARRAY_ID}" \
        --export=SPECIES="${SPECIES}" \
        "${REPO_ROOT}/hpc/puhti/aggregate_submit.sh")

    echo "  [ok]   ${SPECIES}   array=${ARRAY_ID}   agg=${AGG_ID}"
    SUBMITTED_ARRAYS+=("${ARRAY_ID}")
    SUBMITTED_AGGS+=("${AGG_ID}")
done < "${SPECIES_LIST}"

echo
echo "submitted: ${#SUBMITTED_ARRAYS[@]} arrays, ${#SUBMITTED_AGGS[@]} aggregators"
echo "skipped:   ${#SKIPPED[@]}"
if [[ ${#SKIPPED[@]} -gt 0 ]]; then
    printf '  %s\n' "${SKIPPED[@]}"
fi
echo
echo "monitor:   squeue -u \$USER"
echo "once done, pull artifacts back:"
echo "  mkdir -p hpc/artifacts && rsync -a barker@puhti.csc.fi:${PROJECT}/artifacts/ hpc/artifacts/"
