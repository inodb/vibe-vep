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
| **SIFT** | Genomic (chr:pos:ref:alt) | GRCh38 | via dbNSFP | SIFT missense prediction scores (0-1, lower = more damaging) via [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP) |
| **PolyPhen-2** | Genomic (chr:pos:ref:alt) | GRCh38 | via dbNSFP | PolyPhen-2 HDIV missense prediction scores (0-1, higher = more damaging) via [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP) |
| **dbSNP** | Genomic (chr:pos:ref:alt) | GRCh38 | ~17 GB | RS identifiers from [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) |

## Match Levels

**Match level** determines how each source links to variants:

- **Genomic**: matches on exact chr:pos:ref:alt — assembly-specific (GRCh37 vs GRCh38)
- **Protein position**: matches on transcript + amino acid position — transcript-version sensitive. Hotspot positions are only annotated when the annotation's transcript matches the hotspot's transcript.
- **Gene symbol**: matches on gene name only — assembly-independent

## Storage

Genomic annotation sources (AlphaMissense, ClinVar, SIGNAL, gnomAD, SIFT, PolyPhen-2, dbSNP) are merged into a single SQLite database (`genomic_annotations.sqlite`) with a `WITHOUT ROWID` clustered primary key on `(chrom, pos, ref, alt)`. This gives ~1-5μs point lookups with near-zero Go heap via mmap — one DB lookup per variant instead of many.

All variants are stored with normalized coordinates:
- **Chromosome**: without "chr" prefix (e.g. "12", not "chr12")
- **Position and alleles**: canonical MAF-style (no anchor base for indels). VCF-style indels are normalized during build and lookup, so both formats match correctly.

This is separate from the DuckDB variant cache (`variant_cache.duckdb`) used for `--save-results` / `--from-cache` / `export parquet` post-analysis.

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

# SIFT + PolyPhen-2 (via dbNSFP): enable
# Requires dbNSFP per-chromosome files in ~/.vibe-vep/{assembly}/dbnsfp/
vibe-vep config set annotations.sift true
vibe-vep config set annotations.polyphen true

# dbSNP RS IDs: download + enable
vibe-vep config set annotations.dbsnp true
vibe-vep download  # fetches ~17 GB dbsnp.vcf.gz
```

Use `vibe-vep version` to see which sources are loaded and `vibe-vep version --maf-columns` for the full column mapping.
