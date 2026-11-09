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
#   CHUNKS         parallel pipeline chunks         (default 2)
#   JVM_HEAP       heap cap per chunk               (default 2g)
#   OUT_DIR        where outputs are written        (default /home/xiang/data)
set -euo pipefail

TAG=${1:?usage: maf_compare.sh <tag> [input.maf] [reference.maf]}
INPUT=${2:-/home/xiang/data/mskimpact-50k.txt}
REFERENCE=${3:-/home/xiang/data/output_vibe_vep_v60.txt}
OUT_DIR=${OUT_DIR:-/home/xiang/data}
VIBE_BASE=${VIBE_BASE:-http://127.0.0.1:38081}
# Two chunks of a capped heap, because this box has 7 GB and is shared: four
# uncapped JVMs plus vibe-vep is enough to get one of them OOM-killed, which
# leaves a short output file rather than an error.
CHUNKS=${CHUNKS:-2}
JVM_HEAP=${JVM_HEAP:-2g}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

JAR=${PIPELINE_JAR:-}
if [[ -z "$JAR" ]]; then
  JAR=$(find "$HOME/Github/genome-nexus-annotation-pipeline" -name 'annotationPipeline-*.jar' \
        ! -name '*sources*' 2>/dev/null | head -1)
fi
[[ -n "$JAR" && -f "$JAR" ]] || { echo "set PIPELINE_JAR to the annotationPipeline jar" >&2; exit 1; }

curl -sf "$VIBE_BASE/health" >/dev/null || { echo "no vibe-vep at $VIBE_BASE" >&2; exit 1; }

WORK=$(mktemp -d)
# Kept on failure: the chunk logs are the only record of why a chunk died.
trap '[[ -n "${KEEP_WORK:-}" ]] && echo "working files in $WORK" || rm -rf "$WORK"' EXIT

echo "input      $INPUT"
echo "reference  $REFERENCE"
echo "pipeline   $(basename "$JAR")"

head -1 "$INPUT" > "$WORK/hdr"
tail -n +2 "$INPUT" > "$WORK/body"
LINES=$(( ($(wc -l < "$WORK/body") + CHUNKS - 1) / CHUNKS ))
split -l "$LINES" "$WORK/body" "$WORK/part."

for part in "$WORK"/part.*; do
  cat "$WORK/hdr" "$part" > "$part.maf"
  # --mode extended is what adds the Additional_Transcripts and HGVS_Offset
  # columns. Without it the output is a different shape and scores against a
  # reference as if two columns had emptied.
  GENOMENEXUS_BASE="$VIBE_BASE" java -Xmx"$JVM_HEAP" -jar "$JAR" \
    --filename "$part.maf" --output-filename "$part.out" --isoform-override mskcc \
    --mode extended \
    > "$part.log" 2>&1 &
done
wait

# A chunk that was killed leaves no output, and concatenating what is left
# would score a short file as if it were whole.
for part in "$WORK"/part.*.maf; do
  [[ -s "${part%.maf}.out" ]] || {
    KEEP_WORK=1
    echo "chunk ${part##*/} produced no output; see ${part%.maf}.log" >&2
    exit 1
  }
done

OUT="$OUT_DIR/output_vibe_vep_$TAG.txt"
# Each chunk repeats the pipeline's own "#" preamble and the column header,
# so only the first chunk contributes them.
first=1
for part in "$WORK"/part.*.out; do
  if (( first )); then
    cat "$part" > "$OUT"
    first=0
  else
    grep -v -e '^#' -e '^Hugo_Symbol' "$part" >> "$OUT" || true
  fi
done
WROTE=$(grep -vc -e '^#' -e '^Hugo_Symbol' "$OUT")
EXPECTED=$(( $(wc -l < "$WORK/body") ))
echo "wrote $OUT ($WROTE rows)"
[[ "$WROTE" == "$EXPECTED" ]] || { KEEP_WORK=1; echo "expected $EXPECTED rows, got $WROTE" >&2; exit 1; }

# A column that is absent scores as a column that emptied, so the shape is
# checked before the contents are.
if ! diff -q <(grep -m1 '^Hugo_Symbol' "$OUT" | tr '\t' '\n') \
             <(grep -m1 '^Hugo_Symbol' "$REFERENCE" | tr '\t' '\n') >/dev/null; then
  echo "columns differ from the reference:" >&2
  diff <(grep -m1 '^Hugo_Symbol' "$OUT" | tr '\t' '\n') \
       <(grep -m1 '^Hugo_Symbol' "$REFERENCE" | tr '\t' '\n') >&2 || true
  exit 1
fi

python3 -I "$HERE/maf_score.py" "$OUT" "$REFERENCE"
