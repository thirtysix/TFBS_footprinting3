#!/bin/bash
# Scan recent tfbs_agg_*.err logs for Parquet-corrupt failures and
# emit "<species>" one per line to /users/barker/failed_agg_species.txt.
set -u

LOGDIR=/scratch/project_2001307/ensembl_genomes_CAS_scoring_2026/logs
OUT=/users/barker/failed_agg_species.txt
: > "$OUT"

for f in "$LOGDIR"/tfbs_agg_*.err; do
  [ -s "$f" ] || continue
  grep -q "ArrowInvalid\|Parquet file size is 0\|magic bytes not found" "$f" || continue
  sp=$(grep -oE "runs/[a-z_0-9]+/" "$f" | head -1 | awk -F/ '{print $2}')
  [ -n "$sp" ] && echo "$sp"
done | sort -u > "$OUT"

echo "wrote $(wc -l < "$OUT") species to $OUT"
