#!/usr/bin/env bash
#
# Download ALL datahub MAF files: GDC (33 studies) + public/ studies.
# Files are saved to testdata/datahub_all/ organized by assembly.
#
# Public studies are mixed assemblies (GRCh37/GRCh38). The script detects
# NCBI_Build from the first data line and organizes into subdirectories:
#   testdata/datahub_all/GRCh38/  (GDC + GRCh38 public)
#   testdata/datahub_all/GRCh37/  (GRCh37 public)
#   testdata/datahub_all/unknown/ (no NCBI_Build column)
#
set -euo pipefail

DEST_DIR="$(cd "$(dirname "$0")/.." && pwd)/testdata/datahub_all"
mkdir -p "$DEST_DIR/GRCh38" "$DEST_DIR/GRCh37" "$DEST_DIR/unknown"

downloaded=0
skipped=0
failed=0

# Helper: download a study and sort by assembly
download_study() {
  local study="$1"
  local url="$2"
  local filename="${study}_data_mutations.txt"

  # Check if already downloaded in any assembly dir
  for d in "$DEST_DIR/GRCh38" "$DEST_DIR/GRCh37" "$DEST_DIR/unknown"; do
    if [[ -f "$d/$filename" ]]; then
      return 1  # skipped
    fi
  done

  local tmpfile
  tmpfile=$(mktemp)
  if ! curl -fL --progress-bar -o "$tmpfile" "$url" 2>/dev/null; then
    rm -f "$tmpfile"
    return 2  # failed
  fi

  # Detect assembly from NCBI_Build column using a single awk pass.
  # Skips comment lines (#), finds NCBI_Build in header, reads first data line.
  local assembly
  assembly=$(awk -F'\t' '
    /^#/ { next }
    !hdr {
      hdr = 1
      for (i = 1; i <= NF; i++) {
        if ($i == "NCBI_Build") { col = i; break }
      }
      next
    }
    col > 0 {
      print $col
      exit
    }
  ' "$tmpfile")

  case "$assembly" in
    GRCh38|hg38) assembly="GRCh38" ;;
    GRCh37|hg19) assembly="GRCh37" ;;
    *) assembly="unknown" ;;
  esac

  mv "$tmpfile" "$DEST_DIR/$assembly/$filename"
  echo "  -> $assembly"
  return 0
}

# --- Part 1: GDC studies (all GRCh38) ---
echo "=== Downloading GDC TCGA studies ==="

GDC_STUDIES=(
  acc_tcga_gdc aml_tcga_gdc blca_tcga_gdc brca_tcga_gdc ccrcc_tcga_gdc
  cesc_tcga_gdc chol_tcga_gdc chrcc_tcga_gdc coad_tcga_gdc difg_tcga_gdc
  dlbclnos_tcga_gdc esca_tcga_gdc gbm_tcga_gdc hcc_tcga_gdc hgsoc_tcga_gdc
  hnsc_tcga_gdc luad_tcga_gdc lusc_tcga_gdc mnet_tcga_gdc nsgct_tcga_gdc
  paad_tcga_gdc plmeso_tcga_gdc prad_tcga_gdc prcc_tcga_gdc read_tcga_gdc
  skcm_tcga_gdc soft_tissue_tcga_gdc stad_tcga_gdc thpa_tcga_gdc thym_tcga_gdc
  ucec_tcga_gdc ucs_tcga_gdc um_tcga_gdc
)

GDC_BASE="https://media.githubusercontent.com/media/cBioPortal/datahub/master/crdc/gdc"

for study in "${GDC_STUDIES[@]}"; do
  echo "Downloading $study..."
  url="$GDC_BASE/$study/data_mutations.txt"
  download_study "$study" "$url"
  rc=$?
  case $rc in
    0) downloaded=$((downloaded + 1)) ;;
    1) echo "  Skipping (already exists)"; skipped=$((skipped + 1)) ;;
    2) echo "  WARNING: Failed to download"; failed=$((failed + 1)) ;;
  esac
done

# --- Part 2: Public studies ---
echo ""
echo "=== Enumerating public/ studies via GitHub API ==="

# Paginate through GitHub API to get all directories
page=1
all_studies=()
while true; do
  response=$(curl -fsSL "https://api.github.com/repos/cBioPortal/datahub/contents/public?per_page=100&page=$page" 2>/dev/null) || break

  # Extract directory names
  names=$(echo "$response" | python3 -c "
import json, sys
items = json.load(sys.stdin)
if not isinstance(items, list):
    sys.exit(0)
for item in items:
    if item.get('type') == 'dir':
        print(item['name'])
" 2>/dev/null) || break

  if [[ -z "$names" ]]; then
    break
  fi

  while IFS= read -r name; do
    all_studies+=("$name")
  done <<< "$names"

  # Check if we got a full page (need to continue)
  count=$(echo "$names" | wc -l)
  if [[ $count -lt 100 ]]; then
    break
  fi
  page=$((page + 1))
done

echo "Found ${#all_studies[@]} public studies"

PUBLIC_BASE="https://media.githubusercontent.com/media/cBioPortal/datahub/master/public"

for study in "${all_studies[@]}"; do
  echo "Downloading $study..."
  url="$PUBLIC_BASE/$study/data_mutations.txt"
  download_study "$study" "$url"
  rc=$?
  case $rc in
    0) downloaded=$((downloaded + 1)) ;;
    1) echo "  Skipping (already exists)"; skipped=$((skipped + 1)) ;;
    2) failed=$((failed + 1)) ;;  # Silently skip — many public studies lack mutations
  esac
done

echo ""
echo "=== Summary ==="
echo "Downloaded: $downloaded, Skipped: $skipped, Failed: $failed"
echo "GRCh38 files: $(ls "$DEST_DIR/GRCh38/" 2>/dev/null | wc -l)"
echo "GRCh37 files: $(ls "$DEST_DIR/GRCh37/" 2>/dev/null | wc -l)"
echo "Unknown assembly: $(ls "$DEST_DIR/unknown/" 2>/dev/null | wc -l)"
echo "Files saved to $DEST_DIR"
