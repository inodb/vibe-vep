# Tasks: VCF Parsing and Variant Annotation

**Input**: Design documents from `/specs/001-vcf-variant-annotation/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, contracts/cli.md

**Last Updated**: 2026-02-05

## Status Summary

| Phase | Status | Notes |
|-------|--------|-------|
| Phase 1: Setup | ✅ Complete | Go module, directory structure |
| Phase 2: Foundational | ✅ Complete | Core types, test fixtures |
| Phase 3: User Story 1 (VCF Annotation) | ✅ Complete | Basic annotation working |
| Phase 4: User Story 2 (Output Formats) | ⏳ Partial | Tab output done, VCF output pending |
| Phase 5: User Story 3 (Error Handling) | ✅ Complete | Error messages implemented |
| Phase 6: Polish | ✅ Complete | README, help, validation |
| **Extended Features** | ✅ Complete | GENCODE, MAF, REST, validation mode |

---

## Phase 1: Setup (Shared Infrastructure) ✅

- [x] T001 Initialize Go module
- [x] T002 Create directory structure
- [x] T003 Add dependencies to go.mod
- [x] T004 Create .golangci.yml
- [x] T005 Create Makefile

---

## Phase 2: Foundational ✅

### Test Fixtures
- [x] T006 Create testdata/kras_g12c.vcf
- [x] T007 Create testdata/multi_variant.vcf
- [x] T008 Create testdata/cache/ with test data

### Core Types
- [x] T009 Implement Variant type
- [x] T010 Implement Gene type
- [x] T011 Implement Transcript type
- [x] T012 Implement Exon type
- [x] T013 Implement Annotation type
- [x] T014 Implement codon table and TranslateCodon()
- [x] T015 Implement ReverseComplement()
- [x] T016 Write codon translation tests
- [x] T017 Write Variant tests

---

## Phase 3: User Story 1 - VCF Annotation ✅

### Tests
- [x] T018 Test VCF parser
- [x] T019 Test cache loader
- [x] T020 Test consequence classifier
- [x] T021 Test amino acid change calculation
- [x] T022 Integration test: KRAS G12C

### Implementation
- [x] T023 Implement VCF parser
- [x] T024 Implement cache loader (JSON + Sereal)
- [x] T025 Implement transcript lookup (linear search, interval tree TODO)
- [x] T026 Implement GenomicToCDS()
- [x] T027 Implement CDSToCodonPosition()
- [x] T028 Implement consequence classifier
- [x] T029 Implement Annotator.Annotate()
- [x] T030 Implement tab output writer
- [x] T031 Implement CLI with `annotate` subcommand
- [x] T032 Add cache flags
- [x] T033 Handle multi-allelic sites
- [x] T034 Mark canonical transcript

---

## Phase 4: User Story 2 - Output Formats ⏳

- [x] T035 Tab-delimited output (default)
- [ ] T036 VCF output with CSQ INFO field
- [ ] T037 Test --output-format flag
- [ ] T038 Implement VCF output writer

---

## Phase 5: User Story 3 - Error Handling ✅

- [x] T042 File not found error with path
- [x] T043 Invalid format error with line number
- [x] T044 Unknown chromosome warning
- [x] T046 Error wrapping with context
- [x] T048 Cache not found error with hint
- [x] T049 Handle unknown chromosome gracefully
- [x] T050 Define exit codes

---

## Phase 6: Polish ✅

- [x] T051 Gzip support for .vcf.gz
- [x] T052 Stdin support (vibe-vep annotate -)
- [x] T053 Implement --help
- [x] T055 Add README.md
- [x] T056 Verify KRAS G12C output

---

## Extended Features (Beyond Original Spec) ✅

These features were added to improve usability and reduce setup complexity.

### GENCODE Cache (replaces VEP Sereal requirement)
- [x] E001 Implement GTF parser (`internal/cache/gtf_loader.go`)
- [x] E002 Implement FASTA loader (`internal/cache/fasta_loader.go`)
- [x] E003 Implement download command (`cmd/vibe-vep/download.go`)
- [x] E004 Auto-detect GENCODE cache in ~/.vibe-vep/
- [x] E005 Chromosome normalization (chr12 -> 12)
- [x] E006 GTF/FASTA unit tests

### MAF Support
- [x] E007 Implement MAF parser (`internal/maf/parser.go`)
- [x] E008 Implement VariantParser interface (`internal/vcf/interface.go`)
- [x] E009 Auto-detect input format (VCF vs MAF)
- [x] E010 MAF parser tests

### REST API Fallback
- [x] E011 Implement REST loader (`internal/cache/rest_loader.go`)
- [x] E012 Implement CacheWithLoader for on-demand loading
- [x] E013 Add --rest flag to CLI

### Validation Mode
- [x] E014 Implement validation writer (`internal/output/validation.go`)
- [x] E015 Add --validate flag
- [x] E016 Add --validate-all flag
- [x] E017 Generate validation summary statistics
- [x] E018 Validation report for TCGA PAAD

### DuckDB Cache (alternative backend)
- [x] E019 Implement DuckDB loader (`internal/cache/duckdb.go`)
- [x] E020 S3 URL support for remote caches

---

## Remaining Work

| Task | Priority | Description |
|------|----------|-------------|
| VCF output format | P2 | Add --output-format vcf with CSQ field |
| Interval tree | P3 | Optimize transcript lookup for large files |
| GRCh37 validation | P3 | Test with GRCh37 data |

---

## Validation Results

**TCGA PAAD (GRCh38)**: 24,849 variants
- Matches: 15,472 (62.3%)
- Mismatches: 9,377 (37.7%)

See `docs/validation_report_paad.md` for detailed analysis.
