# Implementation Plan: Variant Annotation with GENCODE

**Branch**: `main` | **Updated**: 2026-02-05 | **Spec**: [spec.md](spec.md)
**Status**: Complete

## Summary

vibe-vep is a standalone variant effect predictor that uses GENCODE annotations instead of VEP's Sereal cache. This reduces setup from a 17GB download to 95MB while providing equivalent functionality.

## Technical Stack

| Component | Choice | Rationale |
|-----------|--------|-----------|
| Language | Go 1.22+ | Fast, single binary, no runtime deps |
| Annotations | GENCODE GTF + FASTA | Open format, primary data source |
| Fallback | Ensembl REST API | Zero-setup option |
| Alternative | DuckDB cache | Fast queries, S3 support |

## Project Structure

```
cmd/vibe-vep/
├── main.go              # CLI entry, annotate command
├── download.go          # Download GENCODE files
└── convert.go           # Convert VEP cache to DuckDB

internal/
├── vcf/
│   ├── parser.go        # VCF parsing
│   ├── variant.go       # Variant type
│   └── interface.go     # VariantParser interface
├── maf/
│   └── parser.go        # MAF parsing
├── cache/
│   ├── gtf_loader.go    # GENCODE GTF parser
│   ├── fasta_loader.go  # GENCODE FASTA parser
│   ├── rest_loader.go   # Ensembl REST API
│   ├── duckdb.go        # DuckDB backend
│   ├── loader.go        # VEP Sereal loader (legacy)
│   ├── cache.go         # Cache interface
│   └── transcript.go    # Transcript type
├── annotate/
│   ├── annotator.go     # Main annotation logic
│   ├── consequence.go   # Consequence classification
│   ├── codon.go         # Codon translation
│   └── annotation.go    # Annotation type
└── output/
    ├── tab.go           # Tab-delimited output
    └── validation.go    # Validation mode output

testdata/
├── sample.vcf           # Test VCF files
├── sample.maf           # Test MAF file
├── sample.gtf           # Test GTF (KRAS)
└── sample_cds.fa        # Test FASTA (KRAS CDS)

docs/
└── validation_report_paad.md  # TCGA validation results

specs/001-vcf-variant-annotation/
├── spec.md              # Feature specification
├── plan.md              # This file
├── tasks.md             # Implementation tasks
└── ...                  # Other design docs
```

## Architecture

```
┌─────────────┐     ┌──────────────┐     ┌─────────────┐
│ VCF/MAF     │────▶│  Annotator   │────▶│ Tab Output  │
│ Parser      │     │              │     │             │
└─────────────┘     └──────┬───────┘     └─────────────┘
                          │
                          ▼
              ┌───────────────────────┐
              │    Cache Layer        │
              ├───────────────────────┤
              │ GENCODE (GTF+FASTA)   │ ◀── Primary
              │ REST API              │ ◀── Fallback
              │ DuckDB                │ ◀── Alternative
              └───────────────────────┘
```

## Key Design Decisions

### 1. GENCODE over VEP Sereal

**Problem**: VEP cache is 17GB, uses Perl's Sereal format
**Solution**: Use GENCODE GTF + FASTA directly (95MB)
**Trade-off**: Slightly different transcript versions, 62% match rate

### 2. Chromosome Normalization

**Problem**: GENCODE uses "chr12", VCF/MAF often use "12"
**Solution**: Strip "chr" prefix when loading GTF
**Location**: `internal/cache/gtf_loader.go:normalizeChrom()`

### 3. Multiple Cache Backends

**Problem**: Different users have different needs
**Solution**: Support multiple backends with auto-detection
- GENCODE: Default, smallest download
- REST API: Zero setup, slower
- DuckDB: Fast queries, S3 support

### 4. Validation Mode

**Problem**: Users want to verify annotations against existing data
**Solution**: `--validate` flag compares MAF annotations to predictions
**Output**: Summary statistics + per-variant comparison

## Implementation Phases

### Phase 1: Core VCF Annotation ✅
- VCF parser with vcfgo
- Transcript lookup from test cache
- Consequence classification
- Tab output

### Phase 2: GENCODE Cache ✅
- GTF parser for transcript structures
- FASTA loader for CDS sequences
- Download command
- Auto-detection in ~/.vibe-vep/

### Phase 3: Extended Features ✅
- MAF parser
- REST API fallback
- Validation mode
- DuckDB backend

### Phase 4: Pending
- VCF output with CSQ INFO field
- Interval tree for faster lookups
- GRCh37 validation

## Performance

| Operation | Time | Notes |
|-----------|------|-------|
| Download | ~2 min | 95MB from GENCODE FTP |
| Load cache | ~9 sec | 254K transcripts |
| Annotate 100 variants | <1 sec | After cache loaded |
| Validate 25K variants | ~15 sec | TCGA PAAD dataset |

## Testing

```bash
# Unit tests
go test ./...

# Integration test
vibe-vep annotate testdata/sample.vcf

# Validation test
vibe-vep annotate --validate /path/to/data_mutations.txt
```

## Design Artifacts

- [spec.md](spec.md) - Feature specification (updated)
- [tasks.md](tasks.md) - Implementation tasks (updated)
- [data-model.md](data-model.md) - Entity definitions
- [contracts/cli.md](contracts/cli.md) - CLI interface
- [research.md](research.md) - Technical research
- [quickstart.md](quickstart.md) - User guide

## Next Steps

See [tasks.md](tasks.md) for remaining work:
- VCF output format (P2)
- Interval tree optimization (P3)
- GRCh37 validation (P3)
