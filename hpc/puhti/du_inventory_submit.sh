#!/bin/bash
#SBATCH --job-name=du_inventory
#SBATCH --account=project_2001307
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/project_2001307/ensembl_genomes_CAS_scoring_2026/logs/du_inventory_%j.out
#SBATCH --error=/scratch/project_2001307/ensembl_genomes_CAS_scoring_2026/logs/du_inventory_%j.err

set -u

SCRATCH=/scratch/project_2001307
OUT=/scratch/project_2001307/ensembl_genomes_CAS_scoring_2026/logs/du_inventory_${SLURM_JOB_ID}.txt
mkdir -p "$(dirname "$OUT")"

{
  echo "=== $(date -u +%FT%TZ) ==="
  echo
  echo "=== csc-workspaces ==="
  csc-workspaces 2>&1 || true
  echo
  echo "=== top-level ($SCRATCH/*) ==="
  du -sh "$SCRATCH"/* 2>/dev/null | sort -h
  echo
  CAMP="$SCRATCH/ensembl_genomes_CAS_scoring_2026"
  if [ -d "$CAMP" ]; then
    echo "=== campaign dir ($CAMP/*) ==="
    du -sh "$CAMP"/* 2>/dev/null | sort -h
    echo
    for sub in cas_scores artifacts ensembl_cache tarballs logs tmp outputs; do
      if [ -d "$CAMP/$sub" ]; then
        echo "=== $sub — top 20 children ==="
        du -sh "$CAMP/$sub"/* 2>/dev/null | sort -h | tail -20
        echo
      fi
    done
  fi
  echo "=== done $(date -u +%FT%TZ) ==="
} > "$OUT" 2>&1

echo "Wrote inventory to $OUT"
