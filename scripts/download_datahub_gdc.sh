#!/usr/bin/env bash
#
# Download all 33 GDC TCGA MAF files from cBioPortal/datahub for validation testing.
# Files are saved to testdata/datahub_gdc/ (~5GB total).
#
set -euo pipefail

STUDIES=(
  acc_tcga_gdc
  aml_tcga_gdc
  blca_tcga_gdc
  brca_tcga_gdc
  ccrcc_tcga_gdc
  cesc_tcga_gdc
  chol_tcga_gdc
  chrcc_tcga_gdc
  coad_tcga_gdc
  difg_tcga_gdc
  dlbclnos_tcga_gdc
  esca_tcga_gdc
  gbm_tcga_gdc
  hcc_tcga_gdc
  hgsoc_tcga_gdc
  hnsc_tcga_gdc
  luad_tcga_gdc
  lusc_tcga_gdc
  mnet_tcga_gdc
  nsgct_tcga_gdc
  paad_tcga_gdc
  plmeso_tcga_gdc
  prad_tcga_gdc
  prcc_tcga_gdc
  read_tcga_gdc
  skcm_tcga_gdc
  soft_tissue_tcga_gdc
  stad_tcga_gdc
  thpa_tcga_gdc
  thym_tcga_gdc
  ucec_tcga_gdc
  ucs_tcga_gdc
  um_tcga_gdc
)

BASE_URL="https://media.githubusercontent.com/media/cBioPortal/datahub/master/crdc/gdc"
DEST_DIR="$(cd "$(dirname "$0")/.." && pwd)/testdata/datahub_gdc"

mkdir -p "$DEST_DIR"

downloaded=0
skipped=0
failed=0

for study in "${STUDIES[@]}"; do
  dest="$DEST_DIR/${study}_data_mutations.txt"
  if [[ -f "$dest" ]]; then
    echo "Skipping $study (already exists)"
    skipped=$((skipped + 1))
    continue
  fi
  url="$BASE_URL/$study/data_mutations.txt"
  echo "Downloading $study..."
  if curl -fL --progress-bar -o "$dest" "$url"; then
    downloaded=$((downloaded + 1))
  else
    echo "WARNING: Failed to download $study (may not exist)"
    rm -f "$dest"
    failed=$((failed + 1))
  fi
done

echo ""
echo "Done. Downloaded: $downloaded, Skipped: $skipped, Failed: $failed"
echo "Files saved to $DEST_DIR"
