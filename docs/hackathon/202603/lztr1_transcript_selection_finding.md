# LZTR1 Transcript Selection Finding

**Date:** 2026-03-11  
**Analyst:** vibe-vep validation against MSK-IMPACT Solid Heme dataset  
**Affected variants:** 582 (all LZTR1 SNPs in dataset)

---

## Observation

In the v3 discrepancy report (`data_mutations_allsnps_vibevep_v3_discrepancies.tsv`), all 582 LZTR1
variants show a transcript ID mismatch:

| Field | MSK-IMPACT (input) | vibe-vep (output) |
|---|---|---|
| Transcript ID | `ENST00000215739` | `ENST00000646124` |
| HGVSc | ✅ identical | |
| HGVSp | ✅ identical (538/582) | |
| Variant Classification | ✅ identical (580/582) | |

The 44 records where HGVSp also differs are splice-site variants where vibe-vep omits the
`p.X<pos>_splice` notation that MSK-IMPACT records — an unrelated formatting difference.
The 2 VC differences (`Splice_Region → Silent`) reflect a known splice-boundary classification
difference, not a consequence of the transcript switch.

**In all 582 cases, the biological interpretation is equivalent.** The discrepancy is purely a
transcript identifier change between annotation versions.

---

## Root Cause

### 1. `ENST00000215739` was retired in GENCODE v46

MSK-IMPACT was annotated using an older Ensembl/GENCODE release in which `ENST00000215739` was
the canonical LZTR1 isoform. This transcript is **absent from GENCODE v46lift37**
(`gencode.v46lift37.annotation.gtf.gz`), the GTF currently used by vibe-vep for GRCh37 annotation.

```
$ zgrep "ENST00000215739" ~/.vibe-vep/grch37/gencode.v46lift37.annotation.gtf.gz
(no output — transcript not present)
```

### 2. Biomart canonical override is silently skipped

The Genome Nexus biomart file maps all three canonical columns to `ENST00000215739`:

```
LZTR1   ENST00000215739   ENST00000215739   ENST00000215739
        (Ensembl col 3)   (UniProt col 9)   (MSKCC col 12)
```

`applyCanonicalOverrides` in `internal/cache/gtf_loader.go` checks whether the override transcript
exists in the loaded GTF before applying it. Since `ENST00000215739` is not in the GTF, the
`!found → continue` guard silently skips the override for LZTR1.

### 3. vibe-vep selects the GTF's Ensembl canonical

With no override applied, transcript selection falls back to the GTF's `Ensembl_canonical` tag.
The `parseAttributes` fix (shipped in this branch) correctly reads all `tag` attributes, so
`ENST00000646124` — which carries `Ensembl_canonical`, `MANE_Select`, and `GENCODE_Primary` tags
in v46lift37 — is correctly identified as canonical.

```
ENST00000646124  tags: Ensembl_canonical, MANE_Select, GENCODE_Primary, appris_principal_1
```

This transcript is a GRCh38-native isoform **remapped to GRCh37** coordinates via the GENCODE
liftover pipeline (`remap_status "full_contig"`), and represents the current best-supported
LZTR1 isoform.

---

## This is Not a GRCh38/GRCh37 Assembly Mix-Up

`ENST00000646124` is present in **both** the GRCh37 liftover GTF and the GRCh38 GTF, but
with different genomic coordinates:

| Assembly | Coordinates |
|---|---|
| GRCh37 (v46lift37) | chr22:21,336,586–21,353,321 |
| GRCh38 (v46) | chr22:20,982,297–20,999,032 |

vibe-vep is using the GRCh37 coordinates from the GRCh37 GTF. The transcript ID is shared
across assemblies, which caused the initial impression of an assembly mix-up. The HGVSc
and HGVSp values computed against GRCh37 coordinates are consistent with those from MSK-IMPACT,
confirming the correct assembly is in use.

---

## Impact Assessment

| Category | Count |
|---|---|
| Total affected LZTR1 variants | 582 |
| Transcript ID only differs | 536 |
| Transcript ID + HGVSp formatting | 44 |
| Variant Classification difference | 2 |
| Variants with incorrect biology | **0** |

---

## Proposed Fix

Add `LZTR1 → ENST00000646124` to the canonical overrides, either:

1. **Update the Genome Nexus biomart file** to reference `ENST00000646124` for LZTR1 — this
   will be correct once the biomart source is updated to GENCODE v46.

2. **Add a manual override** in a local overrides layer (similar to the CHEK2 pattern) mapping
   `LZTR1 → ENST00000646124` explicitly. This avoids waiting on upstream biomart updates but
   aligns vibe-vep's output with the GENCODE v46 canonical.

Note: With option 2, the output transcript would differ from MSK-IMPACT's `ENST00000215739` by
design — vibe-vep would be reflecting the current GENCODE canonical rather than the legacy one.
A note in the validation report would be appropriate.

---

## Related Issue

This pattern is structurally identical to the CHEK2 case documented in
[`chek2_transcript_selection_root_cause.md`](./chek2_transcript_selection_root_cause.md):
a transcript present in the biomart overrides file has been retired in GENCODE v46, causing
the override to be silently skipped and the GTF canonical to be used instead. The CHEK2 case
differed in that the GTF canonical was a *different* isoform with different biology; for LZTR1
the biology is equivalent.
