# vibe-vep Regression Report — Latest Binary vs Post-Fix Baseline

**Date:** 2026-03-11  
**Input:** `testdata/msksolidheme/data_mutations_allsnps.txt` (871,305 SNP records)  
**Baseline:** `data_mutations_allsnps_vibevep_grch37_fixed.txt` (commit prior to new updates)  
**Regressed:** `data_mutations_allsnps_vibevep_latest.txt` (binary with new code updates)

---

## Field-Level Match Rates

| Field | Latest Binary | Post-Fix Baseline | Δ |
|---|---|---|---|
| hugo_symbol | 97.3% | 97.3% | ±0.0% |
| consequence | 98.0% | 99.2% | **−1.2%** |
| variant_classification | 99.1% | 99.2% | −0.1% |
| transcript_id | 91.4% | 91.4% | ±0.0% |
| hgvsc | 92.8% | 96.8% | **−3.9%** |
| hgvsp | 89.7% | 93.4% | **−3.8%** |
| hgvsp_short | 92.6% | 96.5% | **−3.9%** |

---

## Variant Classification Differences

**Total:** 3,412 (vs 3,120 in post-fix baseline)

| input.Variant_Classification | vibe.Variant_Classification | Count |
|---|---|---|
| 5'Flank | IGR | 650 |
| Missense_Mutation | 3'UTR | 640 |
| Missense_Mutation | Nonsense_Mutation | 321 |
| Splice_Region | Silent | 197 |
| Missense_Mutation | Intron | 190 |
| Missense_Mutation | Nonstop_Mutation | 178 |
| 5'Flank | Missense_Mutation | 139 |
| Missense_Mutation | 5'UTR | 138 |
| 3'Flank | Missense_Mutation | 133 |
| 5'Flank | RNA | 131 |
| 5'UTR | Missense_Mutation | 94 |
| Intron | Missense_Mutation | 78 |
| Missense_Mutation | Silent | 73 |
| 5'UTR | Intron | 44 |
| Missense_Mutation | Splice_Site | 40 |
| 5'Flank | 5'UTR | 37 |
| Silent | Missense_Mutation | 33 |
| Nonsense_Mutation | Nonstop_Mutation | 30 |
| Missense_Mutation | Translation_Start_Site | 29 |
| Nonsense_Mutation | Missense_Mutation | 21 |
| *(all others)* | | ~268 |

---

## Key Regressions Identified

1. **`Missense_Mutation → Nonsense_Mutation` (321 cases):** New regression not present in post-fix baseline. Suggests a change in stop-codon detection or amino acid consequence logic is incorrectly classifying missense variants as nonsense.

2. **`Missense_Mutation → Nonstop_Mutation` (178 cases):** Increased substantially from ~30 cases in the post-fix baseline. Indicates a similar regression in stop-codon/readthrough logic.

3. **HGVSc/HGVSp −3.9% drop:** 15,600+ additional variants with differing HGVS notation. Strongly correlated with the above consequence classification changes.

4. **`Nonsense_Mutation → Nonstop_Mutation` (30 cases):** Also new/increased.

---

## Assessment

The new code updates introduced a regression in codon consequence prediction, primarily affecting stop-codon handling. The `transcript_id` match rate is unchanged, ruling out a transcript selection cause. The root cause is likely in `internal/annotate/` codon or consequence logic changes introduced in the latest commits.

The fix from the canonical tag parsing (`parseAttributes`) is confirmed to be included (transcript_id match rate stable at 91.4% vs the broken pre-fix baseline of ~88%).
