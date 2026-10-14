#!/usr/bin/env bash
# Derive the transcript -> RefSeq_mRNA mapping VEP uses, from Ensembl's GRCh37
# core database dumps.
#
# GENCODE ships its own metadata.RefSeq, but where the two disagree VEP follows
# Ensembl: on a VEP111 MSK-IMPACT MAF, GENCODE's ordering reproduces 87.33% of
# the RefSeq column and this mapping reproduces 100%. TP53 ENST00000269305 is
# the clearest case -- GENCODE lists NM_000546.5 first, Ensembl NM_001126118.1,
# and VEP reports the latter.
#
# The rule is: the first RefSeq_mRNA xref (external_db_id 1801) in object_xref_id
# order. Predicted accessions (XM_, external_db_id 1806) are excluded -- they
# often sort first and are not what VEP reports.
#
# Usage: build_ensembl_refseq_map.sh <outdir>   (writes ensembl_refseq_mrna.GRCh37.tsv.gz)
set -euo pipefail
out=${1:?usage: build_ensembl_refseq_map.sh <outdir>}
base=https://ftp.ensembl.org/pub/grch37/current/mysql/homo_sapiens_core_116_37
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
for f in transcript.txt.gz xref.txt.gz object_xref.txt.gz; do
  echo "fetching $f" >&2
  curl -sL --fail -o "$tmp/$f" "$base/$f"
done
python3 - "$tmp" "$out" <<'PY'
import gzip, sys, collections
tmp, out = sys.argv[1], sys.argv[2]
tx = {}
with gzip.open(f'{tmp}/transcript.txt.gz', 'rt', errors='replace') as fh:
    for line in fh:
        p = line.rstrip('\n').split('\t')
        for col in p:
            if col.startswith('ENST'):
                tx[p[0]] = col
                break
xr = {}
with gzip.open(f'{tmp}/xref.txt.gz', 'rt', errors='replace') as fh:
    for line in fh:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 4 and p[1] == '1801':      # RefSeq_mRNA only
            xr[p[0]] = p[3]                      # versioned accession
m = collections.defaultdict(list)
with gzip.open(f'{tmp}/object_xref.txt.gz', 'rt', errors='replace') as fh:
    for line in fh:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 4 and p[2] == 'Transcript' and p[3] in xr:
            t = tx.get(p[1])
            if t:
                m[t].append((int(p[0]), xr[p[3]]))
path = f'{out}/ensembl_refseq_mrna.GRCh37.tsv.gz'
with gzip.open(path, 'wt') as fh:
    for t in sorted(m):
        for _, acc in sorted(m[t]):              # object_xref_id order
            fh.write(f'{t}\t{acc}\n')
print(f'wrote {path}: {len(m)} transcripts', file=sys.stderr)
PY
