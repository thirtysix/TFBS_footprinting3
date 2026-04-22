#!/bin/bash
#SBATCH --job-name=runs_inventory
#SBATCH --account=project_2001307
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/project_2001307/ensembl_genomes_CAS_scoring_2026/logs/runs_inventory_%j.out
#SBATCH --error=/scratch/project_2001307/ensembl_genomes_CAS_scoring_2026/logs/runs_inventory_%j.err

set -u

CAMP=/scratch/project_2001307/ensembl_genomes_CAS_scoring_2026
OUT=$CAMP/logs/runs_inventory_${SLURM_JOB_ID}.txt

{
  echo "=== $(date -u +%FT%TZ) runs/ inventory ==="
  echo

  echo "=== per-species size in runs/ (sorted, all) ==="
  du -sh "$CAMP"/runs/*/ 2>/dev/null | sort -h
  echo

  echo "=== species in artifacts/ (already aggregated locally) ==="
  ls -1 "$CAMP"/artifacts/ 2>/dev/null | sort
  echo

  echo "=== species in runs/ ==="
  ls -1d "$CAMP"/runs/*/ 2>/dev/null | xargs -n1 basename | sort
  echo

  echo "=== species in runs/ WITH artifacts (candidates for deletion after S3 verify) ==="
  comm -12 \
    <(ls -1d "$CAMP"/runs/*/ 2>/dev/null | xargs -n1 basename | sort) \
    <(ls -1 "$CAMP"/artifacts/ 2>/dev/null | sort)
  echo

  echo "=== species in runs/ WITHOUT artifacts (in-flight / failed) ==="
  comm -23 \
    <(ls -1d "$CAMP"/runs/*/ 2>/dev/null | xargs -n1 basename | sort) \
    <(ls -1 "$CAMP"/artifacts/ 2>/dev/null | sort)
  echo

  echo "=== artifact file sanity — one per species, CAS + tfs thresholds present? ==="
  for sp in $(ls -1 "$CAMP"/artifacts/ 2>/dev/null | sort); do
    cas=$(ls "$CAMP/artifacts/$sp"/*.CAS_thresholds.jaspar_2026.tsv.gz 2>/dev/null | head -1)
    tfs=$(ls "$CAMP/artifacts/$sp"/*.tfs_thresholds.jaspar_2026.tsv.gz 2>/dev/null | head -1)
    echo "$sp  CAS=$(test -n "$cas" && echo YES || echo NO)  tfs=$(test -n "$tfs" && echo YES || echo NO)"
  done
  echo

  echo "=== done $(date -u +%FT%TZ) ==="
} > "$OUT" 2>&1

echo "Wrote runs inventory to $OUT"
