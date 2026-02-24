#!/usr/bin/env bash
#
# Download TCGA GDC MAF files from cBioPortal/datahub for validation testing.
# Files are saved to testdata/tcga/ (~1.6GB total).
#
set -euo pipefail

STUDIES=(
  blca_tcga_gdc
  brca_tcga_gdc
  chol_tcga_gdc
  coad_tcga_gdc
  gbm_tcga_gdc
  luad_tcga_gdc
  skcm_tcga_gdc
)

BASE_URL="https://media.githubusercontent.com/media/cBioPortal/datahub/master/crdc/gdc"
DEST_DIR="$(cd "$(dirname "$0")/.." && pwd)/testdata/tcga"

mkdir -p "$DEST_DIR"

for study in "${STUDIES[@]}"; do
  dest="$DEST_DIR/${study}_data_mutations.txt"
  if [[ -f "$dest" ]]; then
    echo "Skipping $study (already exists)"
    continue
  fi
  url="$BASE_URL/$study/data_mutations.txt"
  echo "Downloading $study..."
  curl -fL --progress-bar -o "$dest" "$url"
done

echo "Done. Files saved to $DEST_DIR"
