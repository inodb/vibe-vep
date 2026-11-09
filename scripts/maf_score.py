#!/usr/bin/env python3
"""Score two pipeline MAFs row for row, without holding either in memory.

Usage: maf_score.py <candidate> <reference>

Only the columns an annotator decides are compared; everything else in a MAF
is carried through from the input and would agree trivially.
"""
import sys

COLUMNS = [
    "Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Type",
    "Reference_Allele", "Tumor_Seq_Allele2", "Hugo_Symbol", "Entrez_Gene_Id",
    "dbSNP_RS", "Consequence", "Variant_Classification", "HGVSc", "HGVSp",
    "HGVSp_Short", "Transcript_ID", "RefSeq", "Protein_position", "Codons",
    "Exon_Number", "HGVS_Offset", "Annotation_Status",
]


def rows(path):
    header = None
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = {name: i for i, name in enumerate(fields)}
                continue
            yield header, fields


def main(candidate, reference):
    total = identical = 0
    differing = {name: 0 for name in COLUMNS}
    for (lh, left), (rh, right) in zip(rows(candidate), rows(reference)):
        total += 1
        same = True
        for name in COLUMNS:
            if name not in lh or name not in rh:
                continue
            if left[lh[name]] != right[rh[name]]:
                differing[name] += 1
                same = False
        identical += same

    print(f"rows={total}  identical={identical} ({100 * identical / max(total, 1):.4f}%)")
    print(f"  candidate: {candidate}")
    print(f"  reference: {reference}")
    for name in COLUMNS:
        n = differing[name]
        if n:
            print(f"  {name}: {100 * (total - n) / max(total, 1):.2f}% match ({n} differ)")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
