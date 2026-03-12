# vibe-vep Annotation Comparison Report

**Input:** `testdata/msksolidheme/data_mutations_snp500.txt` (500 SNP records from MSK-IMPACT Solid Heme)  
**Output:** `testdata/msksolidheme/data_mutations_snp500_vibevep.txt`  
**Assembly:** GRCh37 (GENCODE v46)  
**Matched variants (by chrom/pos/ref/alt):** 410 / 500

---

## Field-by-Field Agreement

| Field | Match | Mismatch | Match% | Notes |
|---|---|---|---|---|
| `hugo_symbol` | 407 | 3 | **99.3%** | 3 TERT variants not annotated (upstream/IGR boundary) |
| `consequence` | 397 | 13 | **96.8%** | See breakdown below |
| `variant_classification` | 403 | 7 | **98.3%** | Driven by consequence differences |
| `transcript_id` | 0 | 410 | **0%\*** | *Format only* — vibe includes version suffix (e.g. `ENST00000328354` → `ENST00000403642.5_6`); after stripping version: **90%** |
| `hgvsc` | 3 | 407 | **0.7%\*** | *Format only* — input prefixes transcript ID (`ENST...:c.1A>C`), vibe omits it; after stripping: **94.4%** |
| `hgvsp` | 369 | 41 | **90.0%** | Mostly driven by different transcript selection |
| `hgvsp_short` | 381 | 29 | **92.9%** | Same cause as `hgvsp` |

> \* These are formatting differences, not biological disagreements.

---

## Consequence Mismatch Breakdown (13 total)

| Input | vibe-vep | Count |
|---|---|---|
| `upstream_gene_variant` | `intergenic_variant` | 3 |
| `missense_variant` | `missense_variant,NMD_transcript_variant` | 3 |
| `splice_donor_region_variant,intron_variant` | `splice_region_variant,intron_variant` | 1 |
| `missense_variant` | `stop_gained` | 1 |
| `splice_donor_5th_base_variant,intron_variant` | `splice_region_variant,intron_variant` | 1 |
| `splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant` | `splice_region_variant,intron_variant` | 1 |
| `missense_variant` | `3_prime_UTR_variant,NMD_transcript_variant` | 1 |
| `splice_region_variant,synonymous_variant` | `synonymous_variant,splice_region_variant` | 1 *(term order only)* |
| `missense_variant` | `stop_lost` | 1 |

---

## Key Takeaways

- **Format differences** account for most `transcript_id` and `hgvsc` mismatches:
  - vibe-vep appends GENCODE version suffixes to transcript IDs
  - Input `HGVSc` values include the transcript prefix (`ENST...:c.35G>T`); vibe-vep emits only `c.35G>T`
- **Transcript selection** differences explain most `hgvsp`/`hgvsp_short` divergences — vibe-vep picks a different canonical transcript in some cases (e.g. `CHEK2` at chr22:29091207 uses `ENST00000403642` vs `ENST00000328354`)
- **Real biological differences (2 cases):** vibe-vep predicts `stop_gained` or `stop_lost` where the input records `missense_variant` — likely due to a different transcript context
- **SO term granularity:** vibe-vep uses broader Sequence Ontology terms (`splice_region_variant`) compared to finer Ensembl-specific terms in the input (`splice_donor_region_variant`, `splice_donor_5th_base_variant`, `splice_polypyrimidine_tract_variant`)
- **NMD transcripts:** vibe-vep annotates NMD (nonsense-mediated decay) transcripts explicitly; the input does not
