#!/bin/bash
#SBATCH --job-name=cleanup_runs
#SBATCH --account=project_2001307
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/users/barker/cleanup_runs_%j.out
#SBATCH --error=/users/barker/cleanup_runs_%j.err

# Deletes runs/{species}/ for every species that has BOTH
#   artifacts/{species}/*.CAS_thresholds.jaspar_2026.tsv.gz
#   artifacts/{species}/*.tfs_thresholds.jaspar_2026.tsv.gz
#
# Default is DRY RUN. Pass --apply as the first argument to actually rm.

set -u

MODE="${1:-dry}"
CAMP=/scratch/project_2001307/ensembl_genomes_CAS_scoring_2026
OUT=/users/barker/cleanup_runs_${SLURM_JOB_ID}.txt

{
  echo "=== $(date -u +%FT%TZ) cleanup_runs mode=$MODE ==="
  echo

  total_bytes=0
  count=0
  skipped=0
  missing_artifact=0

  for run_dir in "$CAMP"/runs/*/; do
    sp=$(basename "$run_dir")
    art_dir="$CAMP/artifacts/$sp"

    cas_ok=0
    tfs_ok=0
    [ -n "$(ls "$art_dir"/*.CAS_thresholds.jaspar_2026.tsv.gz 2>/dev/null)" ] && cas_ok=1
    [ -n "$(ls "$art_dir"/*.tfs_thresholds.jaspar_2026.tsv.gz 2>/dev/null)" ] && tfs_ok=1

    if [ "$cas_ok" -eq 1 ] && [ "$tfs_ok" -eq 1 ]; then
      bytes=$(du -sb "$run_dir" 2>/dev/null | awk '{print $1}')
      human=$(du -sh "$run_dir" 2>/dev/null | awk '{print $1}')
      total_bytes=$(( total_bytes + bytes ))
      count=$(( count + 1 ))
      if [ "$MODE" = "--apply" ]; then
        rm -rf "$run_dir" && echo "DELETED  $human  $sp"
      else
        echo "WOULD DELETE  $human  $sp"
      fi
    else
      skipped=$(( skipped + 1 ))
      missing_artifact=$(( missing_artifact + 1 ))
    fi
  done

  echo
  echo "=== summary ==="
  echo "species targeted:  $count"
  echo "species skipped:   $skipped (no matching artifact pair)"
  printf "total reclaim:     %s (%d bytes)\n" \
    "$(numfmt --to=iec-i --suffix=B "$total_bytes" 2>/dev/null || echo "$total_bytes")" \
    "$total_bytes"
  echo "mode:              $MODE (use --apply to delete)"
  echo
  echo "=== post-state ==="
  du -sh "$CAMP"/runs 2>/dev/null
  du -sh "$CAMP"/artifacts 2>/dev/null
  echo
  echo "=== done $(date -u +%FT%TZ) ==="
} > "$OUT" 2>&1

echo "Wrote cleanup report to $OUT"
