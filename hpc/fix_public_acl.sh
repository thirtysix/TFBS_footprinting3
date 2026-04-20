#!/bin/bash
# One-shot: set public-read ACL on already-uploaded tarballs that
# inherited private ACLs by mistake. Run once after adding
# s3:PutObjectAcl to the IAM policy.
#
# Usage:
#   bash hpc/fix_public_acl.sh <species-list-file>

set -euo pipefail

SPECIES_LIST="${1:?usage: fix_public_acl.sh <species-list-file>}"
: "${S3_BUCKET_AWS:=s3://tfbssexperimentaldata}"

while IFS= read -r SPECIES; do
    SPECIES="${SPECIES%%#*}"
    SPECIES="${SPECIES// /}"
    [[ -z "${SPECIES}" ]] && continue
    echo "  ${SPECIES}"
    aws s3api put-object-acl \
        --bucket tfbssexperimentaldata \
        --key "${SPECIES}.tar.gz" \
        --acl public-read
done < "${SPECIES_LIST}"

echo "done"
