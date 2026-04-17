#!/bin/bash
# SLURM array submission for Stage C of the CAS distribution pipeline.
#
# Each array task processes CHUNK_SIZE transcripts from a single species'
# list, writing one Parquet shard. See hpc/cas_only.py for the per-task
# logic and hpc/README.md for the overall pipeline.
#
# Usage (from repo root on the submission node):
#
#     export SPECIES=mus_musculus
#     export TOTAL_TRANSCRIPTS=5000          # or whatever Stage A produced
#     export CHUNK_SIZE=10                   # -> 500 array tasks
#     sbatch --array=0-$((TOTAL_TRANSCRIPTS / CHUNK_SIZE - 1)) hpc/slurm/submit.sh
#
# You MUST edit the #SBATCH lines below for your cluster: partition name,
# account, time limit, memory. The defaults here are conservative starting
# points — calibrate after a benchmark run (see the --benchmark flag on
# hpc/cas_only.py).

#SBATCH --job-name=tfbs_cas
#SBATCH --output=hpc/logs/%x_%A_%a.out
#SBATCH --error=hpc/logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
# Uncomment and edit for your cluster:
# #SBATCH --partition=batch
# #SBATCH --account=YOUR_ACCOUNT

set -euo pipefail

: "${SPECIES:?SPECIES env var is required (e.g. mus_musculus)}"
: "${CHUNK_SIZE:=10}"
: "${REPO_ROOT:=$(pwd)}"
: "${CACHE:=${REPO_ROOT}/hpc/cache/ensembl.sqlite}"

cd "${REPO_ROOT}"

LIST="hpc/transcript_lists/${SPECIES}.txt"
OUT_DIR="hpc/cas_scores/${SPECIES}"
START=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))
END=$(( START + CHUNK_SIZE ))
TASK_ID_PADDED=$(printf "%04d" "${SLURM_ARRAY_TASK_ID}")
OUT="${OUT_DIR}/task_${TASK_ID_PADDED}.parquet"

mkdir -p "${OUT_DIR}" hpc/logs

echo "[$(date -Iseconds)] species=${SPECIES} range=${START}:${END} -> ${OUT}"

# Cache is read-only in worker tasks (prefetch.py populates it once on
# the submission node). If you want workers to fall back to live REST
# on cache miss, pass --allow-live-rest.
python -m hpc.cas_only \
    --list "${LIST}" \
    --range "${START}:${END}" \
    --cache "${CACHE}" \
    --output "${OUT}" \
    ${TFBS_CAS_EXTRA:-}

echo "[$(date -Iseconds)] done"
