# Implementation Plan: VCF Parsing and Variant Annotation

**Branch**: `001-vcf-variant-annotation` | **Date**: 2026-01-24 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/001-vcf-variant-annotation/spec.md`

## Summary

Implement basic VCF parsing and variant consequence annotation for vibe-vep. The system will
read VCF files, determine which genes/transcripts are affected by each variant, and classify
the consequence using Sequence Ontology terms. Uses KRAS G12C (chr12:25245350 G>C) as the
primary validation test case.

## Technical Context

**Language/Version**: Go 1.22+
**Primary Dependencies**:
  - github.com/brentp/vcfgo (VCF parsing)
  - github.com/biogo/hts (BGZF compression, tabix)
  - github.com/Sereal/Sereal/Go/sereal (VEP cache reading)
**Storage**: VEP cache files (read-only), no database
**Testing**: go test with table-driven tests
**Target Platform**: Linux, macOS (single static binary)
**Project Type**: Single CLI application
**Performance Goals**: Annotate 1000 variants/second (after cache loading)
**Constraints**: Memory-efficient streaming; no loading entire VCF into memory
**Scale/Scope**: Human genome, GRCh38 only (for this feature)

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Test-First | ✅ PASS | TDD with KRAS G12C as validation; table-driven tests |
| II. Simplicity | ✅ PASS | Single binary, stdin/stdout support, sensible defaults |
| III. Performance | ✅ PASS | Streaming VCF parsing, lazy loading, interval trees |
| IV. Compatibility | ✅ PASS | Reads VEP cache files, outputs VEP-compatible formats |
| V. Documentation | ✅ PASS | CLI help, quickstart.md, error messages with hints |

## Project Structure

### Documentation (this feature)

```text
specs/001-vcf-variant-annotation/
├── plan.md              # This file
├── spec.md              # Feature specification
├── research.md          # Phase 0 research findings
├── data-model.md        # Entity definitions
├── quickstart.md        # Usage documentation
├── contracts/
│   └── cli.md           # CLI interface contract
└── checklists/
    └── requirements.md  # Spec quality checklist
```

### Source Code (repository root)

```text
cmd/
└── vibe-vep/
    └── main.go              # CLI entry point

internal/
├── vcf/
│   ├── parser.go            # VCF parsing wrapper
│   ├── parser_test.go
│   └── variant.go           # Variant type
├── cache/
│   ├── loader.go            # VEP cache loader
│   ├── loader_test.go
│   ├── transcript.go        # Transcript type
│   └── gene.go              # Gene type
├── annotate/
│   ├── annotator.go         # Main annotation logic
│   ├── annotator_test.go
│   ├── consequence.go       # Consequence classification
│   ├── consequence_test.go
│   ├── codon.go             # Codon table and translation
│   └── codon_test.go
└── output/
    ├── writer.go            # Output interface
    ├── tab.go               # Tab-delimited formatter
    ├── vcf.go               # VCF/CSQ formatter
    └── writer_test.go

testdata/
├── kras_g12c.vcf            # Test VCF with KRAS G12C variant
├── multi_variant.vcf        # Multiple variants test
├── multi_allelic.vcf        # Multi-allelic site test
├── invalid.vcf              # Malformed VCF for error testing
└── cache/                   # Minimal VEP cache subset for testing
    └── homo_sapiens/
        └── test_GRCh38/
            └── 12/          # Chromosome 12 data for KRAS
```

**Structure Decision**: Single project layout following Go conventions. All packages under
`internal/` are private. CLI entry point in `cmd/vibe-vep/`. Test fixtures in `testdata/`.

## Complexity Tracking

> **No violations - design follows all Constitution principles.**

| Principle | Design Choice | Justification |
|-----------|---------------|---------------|
| Simplicity | 3 external deps only | vcfgo, hts, sereal - all essential, no bloat |
| Performance | Streaming + interval tree | Required for million-variant scale |
| Compatibility | Sereal format | VEP cache format, requires one-time conversion |

## Phase 0 Artifacts

- [research.md](research.md) - VCF library comparison, VEP cache format, consequence algorithm

## Phase 1 Artifacts

- [data-model.md](data-model.md) - Variant, Gene, Transcript, Exon, Annotation entities
- [contracts/cli.md](contracts/cli.md) - CLI interface specification
- [quickstart.md](quickstart.md) - User documentation and validation test

## Next Steps

Run `/speckit.tasks` to generate the implementation task list based on this plan.
