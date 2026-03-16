# vibe-vep Datahub Benchmark Summary

Comprehensive validation of vibe-vep consequence prediction against cBioPortal datahub MAF annotations.

## Benchmark Tiers

| Tier | Assembly | Studies | Variants | Conseq Match | Time | Command |
|------|----------|---------|----------|-------------|------|---------|
| `datahub_gdc` | GRCh38 | 32 | 2,489,269 | 99.8% | ~2 min | `go test ./internal/output/ -run TestValidationBenchmarkGDC -v -timeout 10m` |
| `datahub_all` GRCh38 | GRCh38 | 54 | 4,432,166 | 99.6% | ~6 min | `go test ./internal/output/ -run TestValidationBenchmarkAll$ -v -timeout 30m` |
| `datahub_all` GRCh37 | GRCh37 | 433 | 16,271,072 | 99.7% | ~19 min | `go test ./internal/output/ -run TestValidationBenchmarkAllGRCh37 -v -timeout 120m` |
| **Total** | | **519** | **~23M** | **99.7%** | **~27 min** | |

## Quick Start

```bash
# Download test data
scripts/download_datahub_gdc.sh          # 32 GDC studies (~5 GB)
scripts/download_datahub_all.sh          # 487 studies (~50 GB, GDC + public)

# Download GENCODE reference
vibe-vep download --assembly GRCh38
vibe-vep download --assembly GRCh37      # only needed for GRCh37 benchmark

# Run fast benchmark (GDC only, ~2 min)
go test ./internal/output/ -run TestValidationBenchmarkGDC -v -count=1 -timeout 10m

# Run full benchmark (all assemblies, ~27 min)
go test ./internal/output/ -run 'TestValidationBenchmark(GDC|All|AllGRCh37)' -v -count=1 -timeout 120m
```

## Corner Case Tests (Fast Feedback Loop)

66 unit tests covering all mismatch patterns, runs in **7ms**:

```bash
go test ./internal/annotate/ -run 'TestDatahub_|TestEdge_|TestCorner_' -v
```

| File | Tests | Source |
|------|-------|--------|
| `corner_case_test.go` | 13 | Manual edge cases + GRCh37 bug fixes |
| `vep_edge_cases_test.go` | 17 | Ensembl VEP test patterns |
| `datahub_mismatch_test.go` | 36 | Datahub mismatch mining |

## Mismatch Collection

Collect all mismatches for analysis:

```bash
# GRCh38 mismatches
go test ./internal/output/ -run TestMismatchCollection -v -count=1 -timeout 60m
# Output: testdata/datahub_all/GRCh38/mismatches.json

# GRCh37 mismatches
go test ./internal/output/ -run TestMismatchCollectionGRCh37 -v -count=1 -timeout 120m
# Output: testdata/datahub_all/GRCh37/mismatches.json
```

## Reports

Each benchmark generates reports in its study directory:

- `validation_report.md` — Human-readable markdown with match rates, category breakdowns, cancer gene analysis, and performance
- `validation_report.json` — Machine-readable JSON for programmatic analysis

All reports include system hardware info (CPU, cores, memory, OS, kernel, Go version).

| Report | Location |
|--------|----------|
| GDC (fast) | `testdata/datahub_gdc/validation_report.md` |
| All GRCh38 | `testdata/datahub_all/GRCh38/validation_report.md` |
| All GRCh37 | `testdata/datahub_all/GRCh37/validation_report.md` |
| GRCh38 mismatches | `testdata/datahub_all/GRCh38/mismatches.json` |

## Iterative Fix Cycle

1. Run corner case tests: `go test ./internal/annotate/ -run 'TestDatahub_|TestEdge_|TestCorner_' -v` (7ms)
2. Fix bugs in `consequence.go` / `hgvsc.go` / `hgvsp.go`
3. Verify: `go test ./...` (all tests pass)
4. Run fast benchmark: `go test ./internal/output/ -run TestValidationBenchmarkGDC -v` (~2 min)
5. Collect new mismatches, mine for corner cases, repeat
