# Tasks: VCF Parsing and Variant Annotation

**Input**: Design documents from `/specs/001-vcf-variant-annotation/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, contracts/cli.md

**Tests**: Included per Constitution Principle I (Test-First). TDD approach: write tests first, ensure they fail, then implement.

**Organization**: Tasks grouped by user story for independent implementation and testing.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (US1, US2, US3)
- All paths are relative to repository root

## Path Conventions

```text
cmd/vibe-vep/           # CLI entry point
internal/vcf/           # VCF parsing
internal/cache/         # VEP cache loading
internal/annotate/      # Annotation logic
internal/output/        # Output formatters
testdata/               # Test fixtures
```

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Project initialization and Go module setup

- [x] T001 Initialize Go module with `go mod init github.com/your-org/vibe-vep`
- [x] T002 [P] Create directory structure: cmd/vibe-vep/, internal/vcf/, internal/cache/, internal/annotate/, internal/output/, testdata/
- [x] T003 [P] Add dependencies to go.mod: github.com/brentp/vcfgo, github.com/biogo/hts, github.com/Sereal/Sereal/Go/sereal
- [x] T004 [P] Create .golangci.yml with project linting configuration
- [x] T005 [P] Create Makefile with build, test, lint targets

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Core types and infrastructure that ALL user stories depend on

**CRITICAL**: No user story work can begin until this phase is complete

### Test Fixtures

- [x] T006 [P] Create testdata/kras_g12c.vcf with single KRAS G12C variant (chr12:25245351 C>A, genomic coords for c.34G>T)
- [x] T007 [P] Create testdata/multi_variant.vcf with 5 variants including synonymous, missense, intergenic
- [x] T008 [P] Create testdata/cache/homo_sapiens/homo_sapiens_GRCh38/12/ with minimal KRAS transcript data (JSON format for testing)

### Core Types (shared across stories)

- [x] T009 [P] Implement Variant type with IsSNV(), IsIndel() methods in internal/vcf/variant.go
- [x] T010 [P] Implement Gene type in internal/cache/gene.go
- [x] T011 [P] Implement Transcript type with Exons slice in internal/cache/transcript.go
- [x] T012 [P] Implement Exon type in internal/cache/transcript.go
- [x] T013 [P] Implement Annotation type with consequence fields in internal/annotate/annotation.go
- [x] T014 Implement codon table and TranslateCodon() in internal/annotate/codon.go
- [x] T015 Implement ReverseComplement() for strand handling in internal/annotate/codon.go

### Core Tests for Foundational Types

- [x] T016 [P] Write table-driven tests for codon translation in internal/annotate/codon_test.go
- [x] T017 [P] Write tests for Variant.IsSNV(), IsIndel() in internal/vcf/variant_test.go

**Checkpoint**: Foundation ready - core types defined, test fixtures created

---

## Phase 3: User Story 1 - Parse Single VCF File (Priority: P1) MVP

**Goal**: Annotate variants in a VCF file with consequence predictions using KRAS G12C as validation

**Independent Test**: Run `vibe-vep annotate testdata/kras_g12c.vcf` and verify output shows KRAS, missense_variant, G12C

### Tests for User Story 1

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [x] T018 [P] [US1] Test VCF parser reads single variant correctly in internal/vcf/parser_test.go
- [x] T019 [P] [US1] Test VEP cache loader finds KRAS transcript in internal/cache/loader_test.go
- [x] T020 [P] [US1] Test consequence classifier returns missense_variant for G12C in internal/annotate/consequence_test.go
- [x] T021 [P] [US1] Test amino acid change calculation returns "G12C" in internal/annotate/consequence_test.go
- [x] T022 [US1] Integration test: annotate KRAS G12C end-to-end in internal/annotate/annotator_test.go

### Implementation for User Story 1

- [x] T023 [US1] Implement VCF parser wrapper using vcfgo in internal/vcf/parser.go
- [x] T024 [US1] Implement VEP cache loader with JSON support in internal/cache/loader.go (Sereal support TODO)
- [ ] T025 [US1] Implement transcript interval tree for position lookup in internal/cache/loader.go (using linear search for now)
- [x] T026 [US1] Implement GenomicToCDS() coordinate mapping in internal/annotate/consequence.go
- [x] T027 [US1] Implement CDSToCodonPosition() in internal/annotate/consequence.go
- [x] T028 [US1] Implement consequence classifier (missense, synonymous, stop_gained, frameshift) in internal/annotate/consequence.go
- [x] T029 [US1] Implement Annotator.Annotate() orchestration in internal/annotate/annotator.go
- [x] T030 [US1] Implement tab-delimited output writer (default format) in internal/output/tab.go
- [x] T031 [US1] Implement CLI main with `annotate` subcommand in cmd/vibe-vep/main.go
- [x] T032 [US1] Add --cache-dir, --species, --assembly flags in cmd/vibe-vep/main.go
- [x] T033 [US1] Handle multi-allelic sites by splitting ALT alleles in internal/vcf/parser.go
- [x] T034 [US1] Mark canonical transcript in annotation output in internal/annotate/annotator.go

**Checkpoint**: User Story 1 complete - `vibe-vep annotate` works with tab output, KRAS G12C validates correctly

---

## Phase 4: User Story 2 - Output Format Selection (Priority: P2)

**Goal**: Support VCF output format with CSQ INFO field in addition to tab-delimited

**Independent Test**: Run `vibe-vep annotate -f vcf testdata/kras_g12c.vcf` and verify valid VCF output with CSQ field

### Tests for User Story 2

- [ ] T035 [P] [US2] Test VCF output writer produces valid VCF with CSQ header in internal/output/vcf_test.go
- [ ] T036 [P] [US2] Test CSQ field format matches VEP specification in internal/output/vcf_test.go
- [ ] T037 [US2] Test --output-format flag switches between tab and vcf in cmd/vibe-vep/main_test.go

### Implementation for User Story 2

- [ ] T038 [US2] Implement VCF output writer with CSQ INFO field in internal/output/vcf.go
- [ ] T039 [US2] Implement Writer interface in internal/output/writer.go
- [ ] T040 [US2] Add --output-format flag to CLI in cmd/vibe-vep/main.go
- [ ] T041 [US2] Add --output flag for file output (default stdout) in cmd/vibe-vep/main.go

**Checkpoint**: User Story 2 complete - both output formats work independently

---

## Phase 5: User Story 3 - Error Handling for Invalid Input (Priority: P3)

**Goal**: Provide clear, actionable error messages for invalid VCF files and missing files

**Independent Test**: Run `vibe-vep annotate nonexistent.vcf` and verify error message includes file path and remediation hint

### Tests for User Story 3

- [ ] T042 [P] [US3] Test file not found error includes path in cmd/vibe-vep/main_test.go
- [ ] T043 [P] [US3] Test invalid VCF format error includes line number in internal/vcf/parser_test.go
- [ ] T044 [P] [US3] Test unknown chromosome warning continues processing in internal/annotate/annotator_test.go
- [ ] T045 [US3] Create testdata/invalid.vcf with malformed line for error testing

### Implementation for User Story 3

- [ ] T046 [US3] Implement error wrapping with context in internal/vcf/parser.go
- [ ] T047 [US3] Add line number tracking to VCF parse errors in internal/vcf/parser.go
- [ ] T048 [US3] Implement cache not found error with download hint in internal/cache/loader.go
- [ ] T049 [US3] Handle unknown chromosome gracefully (warn and continue) in internal/annotate/annotator.go
- [ ] T050 [US3] Define exit codes (0=success, 1=error, 2=usage) in cmd/vibe-vep/main.go

**Checkpoint**: User Story 3 complete - all error cases handled with actionable messages

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Final integration, documentation, and validation

- [ ] T051 [P] Add gzip/BGZF support for .vcf.gz files in internal/vcf/parser.go
- [ ] T052 [P] Add stdin support (vibe-vep annotate -) in cmd/vibe-vep/main.go
- [ ] T053 [P] Implement --help with usage examples in cmd/vibe-vep/main.go
- [ ] T054 Run quickstart.md validation: execute all examples and verify output
- [ ] T055 Add README.md with installation and usage instructions
- [ ] T056 Verify KRAS G12C output matches Ensembl VEP (SC-002 validation)

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies - start immediately
- **Foundational (Phase 2)**: Depends on Setup - BLOCKS all user stories
- **User Story 1 (Phase 3)**: Depends on Foundational - MVP deliverable
- **User Story 2 (Phase 4)**: Depends on Foundational + Writer interface from US1
- **User Story 3 (Phase 5)**: Depends on Foundational - can parallel with US1/US2
- **Polish (Phase 6)**: Depends on all user stories complete

### User Story Dependencies

- **User Story 1 (P1)**: No dependencies on other stories - core MVP
- **User Story 2 (P2)**: Uses Writer interface from US1 implementation
- **User Story 3 (P3)**: Independent of US1/US2 - can run in parallel

### Within Each User Story

1. Tests written FIRST and verified to FAIL
2. Models/types before services
3. Services before CLI integration
4. Core implementation before edge cases

### Parallel Opportunities

**Phase 1 (Setup):**
```
T002, T003, T004, T005 can run in parallel
```

**Phase 2 (Foundational):**
```
T006, T007, T008 (fixtures) in parallel
T009, T010, T011, T012, T013 (types) in parallel
T016, T017 (tests) in parallel
```

**Phase 3 (US1 Tests):**
```
T018, T019, T020, T021 can run in parallel
```

**Phase 4 (US2 Tests):**
```
T035, T036 can run in parallel
```

**Phase 5 (US3 Tests):**
```
T042, T043, T044 can run in parallel
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup
2. Complete Phase 2: Foundational
3. Complete Phase 3: User Story 1
4. **STOP and VALIDATE**: Run `vibe-vep annotate testdata/kras_g12c.vcf`
5. Verify output: KRAS, missense_variant, G12C, CANONICAL=YES

### Incremental Delivery

1. Setup + Foundational → Foundation ready
2. Add User Story 1 → MVP: basic annotation with tab output
3. Add User Story 2 → VCF output format option
4. Add User Story 3 → Robust error handling
5. Polish → Production ready

### Validation Checkpoints

| Checkpoint | Validation |
|------------|------------|
| After US1 | KRAS G12C annotates correctly |
| After US2 | Both output formats produce valid output |
| After US3 | Error cases return actionable messages |
| After Polish | quickstart.md examples all work |

---

## Notes

- [P] tasks = different files, no dependencies
- [US#] label maps task to specific user story
- Constitution requires TDD: tests fail before implementation
- KRAS G12C (chr12:25245351 C>A genomic = c.34G>T coding) is the primary validation test
- Commit after each task or logical group
- Stop at any checkpoint to validate independently
