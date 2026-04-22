#!/bin/bash
#SBATCH --job-name=scrub_parquets
#SBATCH --account=project_2001307
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/users/barker/scrub_parquets_%j.out
#SBATCH --error=/users/barker/scrub_parquets_%j.err

# Two modes:
#   (default)     Histogram parquet file sizes across runs/*/tasks/*.parquet
#                 so you can pick a sensible threshold.
#   --delete-below BYTES [--apply]
#                 Dry-run (without --apply) or delete every parquet file
#                 strictly smaller than BYTES. The regenerating SLURM
#                 array fills the gaps on resubmit.
#
# Examples:
#   sbatch scrub_bad_parquets_submit.sh
#   sbatch scrub_bad_parquets_submit.sh --delete-below 1024
#   sbatch scrub_bad_parquets_submit.sh --delete-below 1024 --apply

set -u

CAMP=/scratch/project_2001307/ensembl_genomes_CAS_scoring_2026
OUT=/users/barker/scrub_parquets_${SLURM_JOB_ID}.txt

MODE=histogram
THRESH=0
APPLY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --delete-below) MODE=delete; THRESH="$2"; shift 2 ;;
    --apply)        APPLY=1; shift ;;
    *)              echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

{
  echo "=== $(date -u +%FT%TZ) scrub_bad_parquets mode=$MODE thresh=$THRESH apply=$APPLY ==="
  echo
  echo "scanning $CAMP/runs/**/TFBSs_found.sortedclusters.parquet ..."

  # Collect sizes once (listing 200K+ small files on Lustre is the slow part)
  SIZES=$(mktemp)
  find "$CAMP"/runs -name 'TFBSs_found.sortedclusters.parquet' -printf '%s %p\n' > "$SIZES"
  total=$(wc -l < "$SIZES")
  echo "found $total parquet files"
  echo

  if [ "$MODE" = "histogram" ]; then
    echo "=== size histogram ==="
    awk '{
      s=$1
      if (s == 0)            b="0B"
      else if (s < 1024)     b="<1KB"
      else if (s < 10*1024)  b="<10KB"
      else if (s < 100*1024) b="<100KB"
      else if (s < 1024*1024)b="<1MB"
      else if (s < 10*1024*1024)b="<10MB"
      else if (s < 100*1024*1024)b="<100MB"
      else                   b=">=100MB"
      c[b]++
    } END {
      for (k in c) printf "  %-10s %d\n", k, c[k]
    }' "$SIZES" | sort -k1
    echo
    echo "=== samples per bucket (first 3) ==="
    awk '{
      s=$1
      if (s == 0)            b="0B"
      else if (s < 1024)     b="<1KB"
      else if (s < 10*1024)  b="<10KB"
      else if (s < 100*1024) b="<100KB"
      else if (s < 1024*1024)b="<1MB"
      else if (s < 10*1024*1024)b="<10MB"
      else if (s < 100*1024*1024)b="<100MB"
      else                   b=">=100MB"
      if (n[b]<3) { print b, $0; n[b]++ }
    }' "$SIZES" | sort

  elif [ "$MODE" = "delete" ]; then
    echo "=== files strictly below $THRESH bytes ==="
    BAD=$(mktemp)
    awk -v t="$THRESH" '$1 < t { print }' "$SIZES" > "$BAD"
    nbad=$(wc -l < "$BAD")
    sumbad=$(awk '{s+=$1} END{print s+0}' "$BAD")
    echo "candidates: $nbad files, total $sumbad bytes"
    echo
    echo "=== candidate species summary ==="
    awk '{n=split($2,a,"/"); print a[n-2]}' "$BAD" | sort | uniq -c | sort -rn | head -30
    echo
    if [ "$APPLY" = "1" ]; then
      echo "=== deleting ==="
      deleted=0
      while IFS= read -r line; do
        path="${line#* }"
        rm -f "$path" && deleted=$((deleted+1))
      done < "$BAD"
      echo "deleted $deleted files"
    else
      echo "(dry-run; re-submit with --apply to actually delete)"
      echo "=== first 10 candidates ==="
      head -10 "$BAD"
    fi
    rm -f "$BAD"
  fi

  rm -f "$SIZES"
  echo
  echo "=== done $(date -u +%FT%TZ) ==="
} > "$OUT" 2>&1

echo "Wrote scrub report to $OUT"
