#!/usr/bin/env bash
# Annotate a cohort through genome-nexus-annotation-pipeline against a running
# vibe-vep, and score the result against a reference MAF.
#
# This is the canonical-transcript half of the comparison: the pipeline writes
# one transcript per variant, so this checks what a MAF carries and nothing
# else. The other half is `go test -tags apidiff -run TestAnnotationFieldParity`,
# which compares every field of every transcript.
#
#   scripts/maf_compare.sh <tag> [input.maf] [reference.maf]
#
# Environment:
#   PIPELINE_JAR   annotationPipeline jar (required unless discoverable)
#   VIBE_BASE      vibe-vep base URL                (default http://127.0.0.1:38081)
#   CHUNKS         parallel pipeline chunks         (default 4)
#   OUT_DIR        where outputs are written        (default /home/xiang/data)
set -euo pipefail

TAG=${1:?usage: maf_compare.sh <tag> [input.maf] [reference.maf]}
INPUT=${2:-/home/xiang/data/mskimpact-50k.txt}
REFERENCE=${3:-/home/xiang/data/output_vibe_vep_v60.txt}
OUT_DIR=${OUT_DIR:-/home/xiang/data}
VIBE_BASE=${VIBE_BASE:-http://127.0.0.1:38081}
CHUNKS=${CHUNKS:-4}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

JAR=${PIPELINE_JAR:-}
if [[ -z "$JAR" ]]; then
  JAR=$(find "$HOME/Github/genome-nexus-annotation-pipeline" -name 'annotationPipeline-*.jar' \
        ! -name '*sources*' 2>/dev/null | head -1)
fi
[[ -n "$JAR" && -f "$JAR" ]] || { echo "set PIPELINE_JAR to the annotationPipeline jar" >&2; exit 1; }

curl -sf "$VIBE_BASE/health" >/dev/null || { echo "no vibe-vep at $VIBE_BASE" >&2; exit 1; }

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

echo "input      $INPUT"
echo "reference  $REFERENCE"
echo "pipeline   $(basename "$JAR")"

head -1 "$INPUT" > "$WORK/hdr"
tail -n +2 "$INPUT" > "$WORK/body"
LINES=$(( ($(wc -l < "$WORK/body") + CHUNKS - 1) / CHUNKS ))
split -l "$LINES" "$WORK/body" "$WORK/part."

for part in "$WORK"/part.*; do
  cat "$WORK/hdr" "$part" > "$part.maf"
  GENOMENEXUS_BASE="$VIBE_BASE" java -jar "$JAR" \
    --filename "$part.maf" --output-filename "$part.out" --isoform-override mskcc \
    > "$part.log" 2>&1 &
done
wait

OUT="$OUT_DIR/output_vibe_vep_$TAG.txt"
first=1
for part in "$WORK"/part.*.out; do
  if (( first )); then cat "$part" > "$OUT"; first=0
  else grep -v '^Hugo_Symbol' "$part" >> "$OUT" || true
  fi
done
echo "wrote $OUT ($(grep -vc '^Hugo_Symbol' "$OUT") rows)"

python3 -I "$HERE/maf_score.py" "$OUT" "$REFERENCE"
