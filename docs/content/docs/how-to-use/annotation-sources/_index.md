---
title: Annotation Sources
description: Data sources available for variant annotation.
weight: 3
aliases:
  - /docs/annotation-sources/
---

## Core

- **GENCODE**: Gene annotations from [GENCODE](https://www.gencodegenes.org/) — transcripts, exons, CDS coordinates

## Annotation Sources (opt-in)

| Source | Match Level | Assembly | Data Size | Description |
|--------|------------|----------|-----------|-------------|
| **OncoKB** | Gene symbol | Any | ~50 KB | Cancer gene classification (ONCOGENE/TSG) from [OncoKB](https://www.oncokb.org/) |
| **AlphaMissense** | Genomic (chr:pos:ref:alt) | GRCh38 | ~643 MB | Missense pathogenicity scores from [AlphaMissense](https://github.com/google-deepmind/alphamissense) (Cheng et al., Science 2023). CC BY 4.0 |
| **ClinVar** | Genomic (chr:pos:ref:alt) | GRCh38 | ~182 MB | Clinical significance from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) (4.1M variants) |
| **Cancer Hotspots** | Protein position (transcript + AA pos) | Any | ~200 KB | Recurrent mutation hotspots from [cancerhotspots.org](https://www.cancerhotspots.org/) |
| **SIGNAL** | Genomic (chr:pos:ref:alt) | GRCh37 only | ~32 MB | Germline mutation frequencies from [SIGNAL](https://signal.mutationalsignatures.com/) |
| **SIFT** | Protein (peptide MD5) | Any | ~4.1 GB (shared) | SIFT missense prediction scores via [Ensembl](https://ftp.ensembl.org/pub/) variation database |
| **PolyPhen-2** | Protein (peptide MD5) | Any | ~4.1 GB (shared) | PolyPhen-2 HDIV missense prediction scores via [Ensembl](https://ftp.ensembl.org/pub/) variation database |
| **dbSNP** | Genomic (chr:pos:ref:alt) | GRCh38 | ~17 GB | RS identifiers from [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) |

## Match Levels

**Match level** determines how each source links to variants:

- **Genomic**: matches on exact chr:pos:ref:alt — assembly-specific (GRCh37 vs GRCh38)
- **Protein position**: matches on transcript + amino acid position — transcript-version sensitive. Hotspot positions are only annotated when the annotation's transcript matches the hotspot's transcript.
- **Protein (peptide MD5)**: matches on MD5 hash of the protein sequence + amino acid position + alternate residue. Assembly-independent since it operates on protein sequences, not genomic coordinates.
- **Gene symbol**: matches on gene name only — assembly-independent

## Storage

### Genomic index

Genomic annotation sources (AlphaMissense, ClinVar, SIGNAL, gnomAD, dbSNP) are merged into a single SQLite database (`genomic_annotations.sqlite`) with a `WITHOUT ROWID` clustered primary key on `(chrom, pos, ref, alt)`. This gives ~1-5us point lookups with near-zero Go heap via mmap — one DB lookup per variant instead of many.

All variants are stored with normalized coordinates:
- **Chromosome**: without "chr" prefix (e.g. "12", not "chr12")
- **Position and alleles**: canonical MAF-style (no anchor base for indels). VCF-style indels are normalized during build and lookup, so both formats match correctly.

### SIFT/PolyPhen-2 predictions

SIFT and PolyPhen-2 predictions are stored in a separate SQLite database (`ensembl_sift_polyphen.sqlite`), built from Ensembl's [variation database MySQL dumps](https://ftp.ensembl.org/pub/). The data contains pre-computed prediction matrices for every possible amino acid substitution in every Ensembl protein.

**Data source:** Ensembl runs SIFT 6.2.1 and PolyPhen-2 2.2.3 on their proteome and publishes the results as `protein_function_predictions.txt.gz` in the MySQL dump for each release. This is the same data that Ensembl VEP uses — no registration or external tools required.

**How it works:**
1. `vibe-vep download` fetches two files from the Ensembl FTP (~4.1 GB total):
   - `translation_md5.txt.gz` — maps translation IDs to peptide sequence MD5 hashes
   - `protein_function_predictions.txt.gz` — gzip-compressed binary prediction matrices
2. `vibe-vep prepare` (or first annotation run) builds a local SQLite keyed by `(md5, analysis)`
3. During annotation, for each missense variant, the CDS sequence is translated to protein, MD5-hashed, and used to look up the prediction matrix. The matrix is indexed by amino acid position and alternate residue to retrieve the SIFT/PolyPhen score and qualitative prediction.

**Binary matrix format:** Each prediction is a 16-bit little-endian value: top 2 bits encode the qualitative prediction, bottom 10 bits encode the score (0-1000, divided by 1000). A value of 0xFFFF indicates the reference amino acid (no prediction). Matrices have a 3-byte header followed by 20 entries per amino acid position (one per possible substitution in alphabetical order: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y).

**Version matching:** vibe-vep downloads predictions from the Ensembl release that matches the GENCODE gene models in use (GENCODE v45 = Ensembl 111, v19 = Ensembl 75). This mapping is maintained in `gencodeEnsemblMap` in `download.go`.

**Coverage:** 92.8% of GENCODE v45 protein-coding transcripts have matching predictions in Ensembl 111 (99.0% of canonical transcripts). The remaining ~7% are proteins where Ensembl did not compute SIFT/PolyPhen predictions (e.g. too short, non-standard). Variants in unmatched proteins simply receive no SIFT/PolyPhen annotation — the same behavior as Ensembl VEP itself.

**Output columns:**

| Column | Range | Description |
|--------|-------|-------------|
| `sift.score` | 0-1 | SIFT score (lower = more damaging) |
| `sift.prediction` | tolerated, deleterious, tolerated_low_confidence, deleterious_low_confidence | Qualitative SIFT prediction |
| `polyphen.score` | 0-1 | PolyPhen-2 HDIV score (higher = more damaging) |
| `polyphen.prediction` | probably_damaging, possibly_damaging, benign, unknown | Qualitative PolyPhen-2 prediction |

### Variant cache

The DuckDB variant cache (`variant_cache.duckdb`) is separate from both — used for `--save-results` / `--from-cache` / `export parquet` post-analysis.

## Enabling Sources

```bash
# AlphaMissense (GRCh38): download + prepare + enable
vibe-vep config set annotations.alphamissense true
vibe-vep download  # fetches ~643 MB
vibe-vep prepare   # builds SQLite index

# ClinVar (GRCh38): download + enable
vibe-vep config set annotations.clinvar true
vibe-vep download  # fetches ~182 MB clinvar.vcf.gz

# Hotspots: point to TSV file
vibe-vep config set annotations.hotspots /path/to/hotspots_v2_and_3d.txt

# SIGNAL (GRCh37 only): enable
vibe-vep config set annotations.signal true

# SIFT + PolyPhen-2 (via Ensembl): download + enable
vibe-vep config set annotations.sift true
vibe-vep config set annotations.polyphen true
vibe-vep download  # fetches ~4.1 GB from Ensembl FTP

# dbSNP RS IDs: download + enable
vibe-vep config set annotations.dbsnp true
vibe-vep download  # fetches ~17 GB dbsnp.vcf.gz
```

Use `vibe-vep version` to see which sources are loaded and `vibe-vep version --maf-columns` for the full column mapping.
