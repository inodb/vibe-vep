# Feature Specification: Variant Annotation with GENCODE

**Feature Branch**: `main` (merged from `001-vcf-variant-annotation`)
**Created**: 2026-01-24
**Updated**: 2026-02-05
**Status**: Complete

## Summary

vibe-vep annotates genomic variants with predicted consequences using GENCODE annotations. It supports VCF and MAF input formats, provides a simple download command for annotation data, and includes validation mode for comparing against existing annotations.

## User Scenarios

### User Story 1 - Annotate Variants (Priority: P1) ✅

A researcher has a VCF or MAF file containing variants. They want to annotate each variant with gene name, consequence type, and amino acid change.

**Acceptance Scenarios**:

1. **Given** a VCF file with KRAS G12C variant, **When** the user runs `vibe-vep annotate input.vcf`, **Then** the output shows KRAS, missense_variant, G12C.

2. **Given** a MAF file, **When** the user runs `vibe-vep annotate input.maf`, **Then** all variants are annotated with consequences.

3. **Given** no local cache, **When** the user runs `vibe-vep annotate --rest input.vcf`, **Then** variants are annotated using Ensembl REST API.

---

### User Story 2 - Download Annotations (Priority: P1) ✅

A user wants to set up vibe-vep with minimal effort. They need a simple command to download required annotation data.

**Acceptance Scenarios**:

1. **Given** a fresh install, **When** the user runs `vibe-vep download --assembly GRCh38`, **Then** GENCODE files are downloaded to ~/.vibe-vep/ (~95MB).

2. **Given** downloaded annotations, **When** the user runs `vibe-vep annotate`, **Then** the tool auto-detects and uses the local cache.

---

### User Story 3 - Validate MAF Annotations (Priority: P2) ✅

A researcher has a MAF file with existing annotations (e.g., from cBioPortal). They want to compare those annotations against vibe-vep predictions to find discrepancies.

**Acceptance Scenarios**:

1. **Given** a MAF with HGVSp annotations, **When** the user runs `vibe-vep annotate --validate data_mutations.txt`, **Then** output shows matches and mismatches with summary statistics.

2. **Given** `--validate-all` flag, **When** validation runs, **Then** all variants are shown (not just mismatches).

---

### User Story 4 - Output Format Selection (Priority: P3) ⏳

**Status**: Partially implemented (tab output only, VCF output pending)

A researcher needs VCF output with CSQ INFO field for pipeline integration.

---

## Requirements

### Functional Requirements (Implemented)

- **FR-001**: ✅ Parse VCF files (VCF 4.1+, gzipped supported)
- **FR-002**: ✅ Parse MAF files (cBioPortal format)
- **FR-003**: ✅ Auto-detect input format from extension or content
- **FR-004**: ✅ Load transcript annotations from GENCODE GTF
- **FR-005**: ✅ Load CDS sequences from GENCODE FASTA
- **FR-006**: ✅ Download GENCODE files via `vibe-vep download`
- **FR-007**: ✅ Auto-detect GENCODE cache in ~/.vibe-vep/
- **FR-008**: ✅ Fall back to Ensembl REST API with `--rest` flag
- **FR-009**: ✅ Classify consequences using Sequence Ontology terms
- **FR-010**: ✅ Calculate amino acid changes for coding variants
- **FR-011**: ✅ Output tab-delimited format (VEP-compatible)
- **FR-012**: ✅ Validation mode comparing MAF annotations to predictions
- **FR-013**: ✅ Handle multi-allelic VCF sites
- **FR-014**: ✅ Normalize chromosome names (chr12 → 12)
- **FR-015**: ⏳ VCF output with CSQ INFO field (pending)

### Data Sources

| Source | Format | Size | Purpose |
|--------|--------|------|---------|
| GENCODE GTF | gzipped GTF | ~50MB | Transcript structures |
| GENCODE FASTA | gzipped FASTA | ~45MB | CDS sequences |
| Ensembl REST | JSON API | on-demand | Fallback when no cache |

### Key Entities

- **Variant**: Genomic change (chrom, pos, ref, alt) from VCF or MAF
- **Transcript**: Gene isoform with exon structure, CDS coordinates, canonical status
- **Annotation**: Predicted effect (consequence type, amino acid change, impact)
- **Cache**: GENCODE GTF + FASTA providing transcript annotations

## Success Criteria

| Criterion | Status | Result |
|-----------|--------|--------|
| SC-001: KRAS G12C annotates correctly | ✅ Pass | missense_variant, G12C |
| SC-002: Download < 5 minutes | ✅ Pass | ~95MB, ~2 min |
| SC-003: No VEP installation required | ✅ Pass | Standalone binary |
| SC-004: TCGA validation >60% match | ✅ Pass | 62.3% on PAAD |

## Architecture Decisions

### Why GENCODE instead of VEP Cache?

The original plan used VEP's Sereal cache format. We switched to GENCODE because:

| Aspect | VEP Sereal | GENCODE |
|--------|------------|---------|
| Download size | 17GB | 95MB |
| Format | Proprietary Perl | Open GTF/FASTA |
| Dependencies | Perl, Sereal | None |
| Data source | Derived from GENCODE | Primary source |

### Why MAF Support?

Cancer genomics workflows use MAF format (cBioPortal, GDC). Supporting MAF enables:
- Re-annotation of existing datasets
- Validation of existing annotations
- Integration with TCGA/GDC pipelines

### Why REST API Fallback?

Users can annotate without downloading anything:
```bash
vibe-vep annotate --rest input.vcf
```
Slower but useful for quick tests or CI environments.

## Validation Results

**TCGA PAAD (GRCh38)**: 24,849 variants tested
- Matches: 15,472 (62.3%)
- Mismatches: 9,377 (37.7%)

Primary mismatch causes:
- Different transcript versions (49% of mismatches)
- Different canonical transcript selection (14%)
- Splice site boundary definitions (4%)

See `docs/validation_report_paad.md` for detailed analysis.
