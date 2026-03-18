---
title: Development
description: Setting up a development environment and contributing to vibe-vep.
weight: 4
---

## Prerequisites

- Go 1.24+ with `CGO_ENABLED=1` (required for DuckDB and SQLite)
- Git

## Building

```bash
git clone https://github.com/inodb/vibe-vep.git
cd vibe-vep
CGO_ENABLED=1 go build -o vibe-vep ./cmd/vibe-vep
```

## Running Tests

```bash
# Run all tests (fast mode, ~1s)
go test ./... -short -count=1

# Run corner case tests (66 tests, 7ms)
go test ./internal/annotate/ -run 'TestDatahub_|TestEdge_|TestCorner_' -v

# Run benchmarks
go test ./internal/annotate/ -bench . -benchmem

# Run fast validation benchmark (GDC, ~2 min)
go test ./internal/output/ -run TestValidationBenchmarkGDC -v -count=1 -timeout 10m

# Run full validation benchmark (all assemblies, ~21 min)
go test ./internal/output/ -run 'TestValidationBenchmark(GDC|All|AllGRCh37)' -v -count=1 -timeout 120m

# Run annotation sources benchmark
go test ./internal/output/ -run TestAnnotationSourcesBenchmark -v -count=1 -timeout 30m
```

## Test Data

Three tiers of test data are available:

| Tier | Size | Studies | Command |
|------|------|---------|---------|
| `datahub_gdc` | ~5 GB | 32 GDC studies (GRCh38) | `make download-datahub-gdc` |
| `datahub_all` | ~50 GB | 431+ studies (GRCh38 + GRCh37) | `make download-datahub-all` |

For backward compatibility, `make download-testdata` still downloads the original 7-study TCGA subset.

## Project Structure

```
cmd/vibe-vep/       CLI entry point (annotate, download commands)
internal/
  annotate/         Consequence prediction (PredictConsequence, Annotator)
    corner_case_test.go       13 edge case tests
    vep_edge_cases_test.go    17 Ensembl VEP test patterns
    datahub_mismatch_test.go  36 datahub mismatch tests
  cache/            Transcript cache (GENCODE GTF/FASTA loader)
  duckdb/           DuckDB cache for transcripts and variant results
  genomicindex/     Unified SQLite index for annotation source lookups (AM, ClinVar, SIGNAL)
  maf/              MAF file parser
  output/           Output formatting and validation comparison
  vcf/              VCF file parser
testdata/
  cache/            Test transcript data (JSON)
  datahub_gdc/      32 GDC studies for GRCh38 validation (~5 GB)
  datahub_all/      All datahub studies (GRCh38 + GRCh37, ~50 GB)
```

## Roadmap

- [ ] **Feature parity for MAF annotation** — Match the annotation capabilities of the [genome-nexus-annotation-pipeline](https://github.com/genome-nexus/genome-nexus-annotation-pipeline)
  - [x] Consequence prediction (99.8% concordance)
  - [x] HGVSp/HGVSc notation
  - [x] Full MAF output format
  - [x] Cancer gene annotations (OncoKB)
  - [x] VCF to MAF conversion
  - [x] `--pick` / `--most-severe` annotation filtering
- [ ] **Additional annotation sources**
  - [x] AlphaMissense pathogenicity scores
  - [x] ClinVar clinical significance
  - [x] Cancer Hotspots
  - [x] SIGNAL germline frequencies
  - [x] SIFT predictions
  - [x] PolyPhen-2 predictions
  - [ ] gnomAD allele frequencies
- [x] **Re-annotate datahub GDC studies**
- [ ] **Replace genome-nexus-annotation-pipeline for datahub**

## License

MIT License
