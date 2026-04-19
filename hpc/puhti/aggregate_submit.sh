#!/bin/bash -l
# Run hpc/puhti/build_cas_distributions.py on Puhti after a species'
# array completes, so only the ~800 KB threshold TSV needs to be rsync'd
# back -- not the ~9 GB of parquet shards.
#
# Typical chained submission (so this kicks off automatically once the
# array succeeds):
#
#   ARRAY_ID=$(sbatch --parsable --export=SPECIES=${SPECIES} \
#       TFBS_footprinting3/hpc/puhti/array_submit.sh)
#   sbatch --dependency=afterok:${ARRAY_ID} --export=SPECIES=${SPECIES} \
#       TFBS_footprinting3/hpc/puhti/aggregate_submit.sh
#
# `afterok:` skips aggregation if any array task failed, matching the
# "only aggregate complete data" contract. If you want to aggregate a
# partial run (e.g. 99/100 tasks), submit this by itself without the
# dependency flag.

#SBATCH --job-name=tfbs_agg
#SBATCH --account=project_2001307
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -eo pipefail

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
: "${ARTIFACTS_DIR:=${PROJECT}/artifacts}"
: "${MIN_HITS:=100}"

mkdir -p "${RUN_DIR}/logs"
mkdir -p "${ARTIFACTS_DIR}"

module load python-data/3.10
source "${REPO_ROOT}/venv/bin/activate"

# Count completed parquets so we fail fast if the array didn't actually
# produce output (e.g. dependency bypassed via --dependency=afternotok
# by mistake, or a stale $RUN_DIR/tasks/ from a previous species).
PARQUET_COUNT=$(find "${RUN_DIR}/tasks" -name 'TFBSs_found.sortedclusters.parquet' 2>/dev/null | wc -l)
echo "[$(date -Iseconds)] species=${SPECIES} parquets_found=${PARQUET_COUNT}"
if [[ "${PARQUET_COUNT}" -lt 1 ]]; then
    echo "FATAL: no parquets under ${RUN_DIR}/tasks/ -- did the array run?" >&2
    exit 1
fi

python "${REPO_ROOT}/hpc/puhti/build_cas_distributions.py" \
    --species "${SPECIES}" \
    --results-dir "${RUN_DIR}/tasks" \
    --out-dir "${ARTIFACTS_DIR}" \
    --min-hits "${MIN_HITS}"

echo "[$(date -Iseconds)] done"
ls -lh "${ARTIFACTS_DIR}/${SPECIES}/"
