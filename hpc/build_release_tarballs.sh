#!/bin/bash
# Build updated per-species tarballs for the JASPAR 2026 release.
#
# For each species in the input list:
#   1. Download the current tarball from S3 (the ground-truth archive)
#   2. Extract into a scratch workdir
#   3. Inject the new 2026 CAS + PWM threshold TSVs from hpc/artifacts/
#   4. Re-tar
#   5. Optionally upload back to S3 (with pre-upload backup)
#
# The old 2018-era files (tfs_thresholds.jaspar_2018.tsv.gz,
# CAS_pvalues.0.1.tf_ls.json if present) are preserved so older tool
# installs that glob for them continue to work. The new files sit
# alongside; the 2026 branch of tfbs_footprinter3 picks them up via
# data_loader's sort-and-take-latest glob logic.
#
# Usage:
#   bash hpc/build_release_tarballs.sh <species-list-file> [--upload]
#
# Without --upload, tarballs land in /tmp/release_tarballs/ for review.
# With --upload, they're also pushed to S3 after backing up each species'
# existing tarball to s3://tfbssexperimentaldata/backup/<date>/.

set -euo pipefail

SPECIES_LIST="${1:?usage: build_release_tarballs.sh <species-list-file> [--upload]}"
UPLOAD=false
[[ "${2:-}" == "--upload" ]] && UPLOAD=true

: "${ARTIFACTS_DIR:=hpc/artifacts}"
: "${OUT_DIR:=/tmp/release_tarballs}"
: "${S3_BUCKET_HTTP:=https://s3.us-east-2.amazonaws.com/tfbssexperimentaldata}"
: "${S3_BUCKET_AWS:=s3://tfbssexperimentaldata}"
: "${BACKUP_PREFIX:=backup/$(date +%Y%m%d)}"

mkdir -p "${OUT_DIR}"

total=0
uploaded=0
skipped=()

while IFS= read -r SPECIES; do
    SPECIES="${SPECIES%%#*}"
    SPECIES="${SPECIES// /}"
    [[ -z "${SPECIES}" ]] && continue

    total=$((total + 1))
    echo "=== ${SPECIES} ==="

    CAS_TSV="${ARTIFACTS_DIR}/${SPECIES}/${SPECIES}.CAS_thresholds.jaspar_2026.tsv.gz"
    PWM_TSV="${ARTIFACTS_DIR}/${SPECIES}/${SPECIES}.tfs_thresholds.jaspar_2026.tsv.gz"
    if [[ ! -s "${CAS_TSV}" ]] || [[ ! -s "${PWM_TSV}" ]]; then
        echo "  [skip] missing artifact TSV(s) in ${ARTIFACTS_DIR}/${SPECIES}/"
        skipped+=("${SPECIES}")
        continue
    fi

    WORK="${OUT_DIR}/_work_${SPECIES}"
    rm -rf "${WORK}"
    mkdir -p "${WORK}"

    # Pull current tarball from the public bucket endpoint
    if ! curl -fsSL -o "${WORK}/${SPECIES}.tar.gz" "${S3_BUCKET_HTTP}/${SPECIES}.tar.gz"; then
        echo "  [skip] could not download ${S3_BUCKET_HTTP}/${SPECIES}.tar.gz"
        skipped+=("${SPECIES}")
        rm -rf "${WORK}"
        continue
    fi

    tar -xzf "${WORK}/${SPECIES}.tar.gz" -C "${WORK}"
    rm "${WORK}/${SPECIES}.tar.gz"

    # Inject new TSVs alongside the existing files
    cp "${CAS_TSV}" "${WORK}/${SPECIES}/"
    cp "${PWM_TSV}" "${WORK}/${SPECIES}/"

    # Re-tar
    NEW_TARBALL="${OUT_DIR}/${SPECIES}.tar.gz"
    tar -czf "${NEW_TARBALL}" -C "${WORK}" "${SPECIES}"
    rm -rf "${WORK}"

    NEW_SIZE=$(du -h "${NEW_TARBALL}" | cut -f1)
    echo "  [ok] ${NEW_SIZE} ${NEW_TARBALL}"

    if ${UPLOAD}; then
        # Backup current S3 object to versioned path first
        BACKUP_KEY="${BACKUP_PREFIX}/${SPECIES}.tar.gz"
        echo "  [backup] ${S3_BUCKET_AWS}/${SPECIES}.tar.gz -> ${S3_BUCKET_AWS}/${BACKUP_KEY}"
        aws s3 cp "${S3_BUCKET_AWS}/${SPECIES}.tar.gz" "${S3_BUCKET_AWS}/${BACKUP_KEY}" \
            --only-show-errors
        echo "  [upload] ${NEW_TARBALL} -> ${S3_BUCKET_AWS}/${SPECIES}.tar.gz"
        # --acl public-read: the bucket has a global AllUsers:READ grant at
        # the bucket level, but new objects otherwise inherit private ACLs.
        # The tool downloads via unauthenticated HTTP from
        # s3.us-east-2.amazonaws.com, so every release tarball must be
        # explicitly set public-readable.
        aws s3 cp "${NEW_TARBALL}" "${S3_BUCKET_AWS}/${SPECIES}.tar.gz" --acl public-read --only-show-errors
        uploaded=$((uploaded + 1))
    fi
done < "${SPECIES_LIST}"

echo
echo "built:    ${total}"
echo "uploaded: ${uploaded}"
echo "skipped:  ${#skipped[@]}"
if [[ ${#skipped[@]} -gt 0 ]]; then
    printf '  %s\n' "${skipped[@]}"
fi
echo
echo "tarballs under ${OUT_DIR}/"
ls -lh "${OUT_DIR}"/*.tar.gz 2>/dev/null || true
