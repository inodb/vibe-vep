# Exploration: vibe-vep Selecting ENST00000403642 Instead of ENST00000328354 for CHEK2

**Variant:** CHEK2 chr22:29091207 G>A  
**Expected transcript:** `ENST00000328354` (MSKCC/UniProt preferred isoform)  
**Actual transcript:** `ENST00000403642`

---

## Root Cause

**`ENST00000328354` does not exist in the GENCODE v46lift37 GTF** — it was retired/replaced in GENCODE v46. The canonical override file (`ensembl_biomart_canonical_transcripts_per_hgnc.txt`) correctly lists it as the MSKCC/UniProt preferred isoform, but `applyCanonicalOverrides` has a guard:

```go
if !found {
    continue  // override is silently skipped if transcript isn't in the GTF
}
```

So the override for CHEK2 is **never applied**, and the GTF's own `Ensembl_canonical` tag takes effect instead — which points to **`ENST00000404276`** (MANE Select, `appris_principal_2`).

**Why is `ENST00000403642` selected then?** It's not canonical either — for the specific variant at chr22:29091207, `ENST00000403642` is being selected over `ENST00000404276` likely due to **impact ranking**: `ENST00000404276` may produce a lower-impact consequence at that position (or the variant falls outside its CDS), while `ENST00000403642` produces a missense. Transcript selection in `SelectBestAnnotation` (`internal/output/compare.go`) prioritises canonical first, then impact — so a non-canonical higher-impact annotation can win over a canonical lower-impact one.

---

## CHEK2 Transcript Summary (GENCODE v46lift37)

| Transcript | In GTF | Ensembl_canonical | In Override File | Notes |
|---|---|---|---|---|
| `ENST00000328354` | ❌ **absent** | — | ✅ MSKCC/UniProt preferred | Retired in GENCODE v46 |
| `ENST00000404276` | ✅ | ✅ **MANE Select** | — | GTF canonical (`appris_principal_2`), but may not cover this variant position |
| `ENST00000403642` | ✅ | ❌ | — | `GENCODE_Primary`, selected by impact ranking |

---

## Code Path

1. **`internal/cache/canonical.go`** — `LoadCanonicalOverrides()` reads col 8 (`uniprot_canonical_transcript`) from the Genome Nexus biomart file, mapping `CHEK2 → ENST00000328354`
2. **`internal/cache/gtf_loader.go`** — `applyCanonicalOverrides()` (line 438) silently skips the override when the target transcript is not found in the GTF
3. **`internal/cache/gtf_loader.go`** (line 132) — `ENST00000404276` receives `IsCanonical = true` from the GTF `Ensembl_canonical` tag
4. **`internal/output/compare.go`** — `SelectBestAnnotation()` / `AnnotationBetter()` selects `ENST00000403642` over `ENST00000404276` via impact ranking at this position

---

## Suggested Fixes

**Option A — Update the override file mapping:**  
Map `CHEK2 → ENST00000404276` (the current GENCODE v46 MANE Select transcript) in the canonical overrides, replacing the retired `ENST00000328354`.

**Option B — Fall back gracefully when the override transcript is absent:**  
In `applyCanonicalOverrides`, instead of silently skipping when the target transcript is not found, log a warning and leave the GTF `Ensembl_canonical` assignment intact (current behaviour already does this, but it should be made explicit and observable).

**Option C — Prefer `Ensembl_canonical` over impact in `AnnotationBetter`:**  
Ensure that a transcript marked `IsCanonical` (`ENST00000404276`) always wins over a non-canonical transcript (`ENST00000403642`) regardless of impact. Review whether the current ranking in `AnnotationBetter` achieves this consistently.
