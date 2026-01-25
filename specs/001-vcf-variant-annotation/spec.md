# Feature Specification: Basic VCF Parsing and Variant Annotation

**Feature Branch**: `001-vcf-variant-annotation`
**Created**: 2026-01-24
**Status**: Draft
**Input**: User description: "basic vcf parsing and variant annotation. Use a common example for testing e.g. KRAS G12C"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Parse Single VCF File (Priority: P1)

A researcher has a VCF file containing variants from a tumor sample. They want to run vibe-vep
against the file and receive annotated output showing the predicted effect of each variant
on genes and proteins.

**Why this priority**: This is the core functionality - without VCF parsing and basic annotation,
no other features can work. This represents the minimal viable product.

**Independent Test**: Can be tested by running vibe-vep with a small VCF file containing the
KRAS G12C variant (chr12:25245350 G>C in GRCh38) and verifying the output correctly identifies
it as a missense variant affecting KRAS codon 12 (Gly→Cys).

**Acceptance Scenarios**:

1. **Given** a valid VCF file with one variant (KRAS G12C), **When** the user runs
   `vibe-vep annotate input.vcf`, **Then** the output shows the variant with gene name (KRAS),
   consequence type (missense_variant), and amino acid change (G12C).

2. **Given** a VCF file with multiple variants, **When** the user runs vibe-vep,
   **Then** all variants are annotated and output in the same order as input.

3. **Given** an empty VCF file (header only, no variants), **When** the user runs vibe-vep,
   **Then** the tool exits successfully with an informational message indicating zero variants
   processed.

---

### User Story 2 - Output Format Selection (Priority: P2)

A researcher needs to integrate vibe-vep output into their existing pipeline. They want to
choose between VEP-compatible tab-delimited output and VCF output with annotations in the
INFO field.

**Why this priority**: Output flexibility enables integration with existing workflows, but
basic annotation must work first.

**Independent Test**: Run vibe-vep with the same input file using different output format
flags and verify each format is correctly structured.

**Acceptance Scenarios**:

1. **Given** a VCF input file, **When** the user runs `vibe-vep annotate --output-format vcf input.vcf`,
   **Then** the output is a valid VCF file with annotations in the CSQ INFO field (VEP-compatible format).

2. **Given** a VCF input file, **When** the user runs `vibe-vep annotate --output-format tab input.vcf`,
   **Then** the output is tab-delimited text matching VEP's default output format.

3. **Given** no output format specified, **When** the user runs vibe-vep,
   **Then** the default output format is tab-delimited (matching VEP behavior).

---

### User Story 3 - Error Handling for Invalid Input (Priority: P3)

A researcher accidentally provides a malformed VCF file or specifies a non-existent file.
They want clear error messages that help them fix the problem.

**Why this priority**: Robust error handling improves user experience but is not core
functionality.

**Independent Test**: Run vibe-vep with various invalid inputs and verify error messages
are clear and actionable.

**Acceptance Scenarios**:

1. **Given** a file path that does not exist, **When** the user runs vibe-vep,
   **Then** the error message states the file was not found and includes the path attempted.

2. **Given** a file that is not valid VCF format, **When** the user runs vibe-vep,
   **Then** the error message indicates the file format problem and the line number where
   parsing failed.

3. **Given** a VCF with variants on chromosomes not in the reference data, **When** the user
   runs vibe-vep, **Then** those variants are marked as unannotated with a clear reason,
   and processing continues for valid variants.

---

### Edge Cases

- What happens when a variant position is outside gene boundaries (intergenic)?
  → Annotate as intergenic_variant with nearest gene information.
- What happens when multiple transcripts exist for a gene?
  → Report annotations for all transcripts (canonical transcript flagged).
- What happens when the VCF contains structural variants (SVs)?
  → For this initial version, skip SVs with an informational warning; SNVs and small indels only.
- What happens with multi-allelic sites (multiple ALT alleles)?
  → Split and annotate each alternate allele separately.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST parse VCF files conforming to VCF 4.1+ specification
- **FR-002**: System MUST identify the genomic location (chromosome, position, ref, alt) for each variant
- **FR-003**: System MUST determine which gene(s) and transcript(s) overlap each variant position
- **FR-004**: System MUST classify the consequence type using Sequence Ontology terms
  (e.g., missense_variant, synonymous_variant, stop_gained, frameshift_variant)
- **FR-005**: System MUST calculate amino acid change for coding variants (e.g., G12C)
- **FR-006**: System MUST output results in VEP-compatible tab-delimited format by default
- **FR-007**: System MUST support VCF output format with CSQ INFO field annotations
- **FR-008**: System MUST read gene/transcript annotations from VEP cache files
- **FR-009**: System MUST provide clear error messages for invalid input with line numbers
- **FR-010**: System MUST handle multi-allelic VCF sites by annotating each allele separately
- **FR-011**: System MUST indicate canonical transcript in output when multiple transcripts exist
- **FR-012**: System MUST support reading gzipped VCF files (.vcf.gz)

### Key Entities

- **Variant**: Represents a single genomic change with chromosome, position, reference allele,
  alternate allele, and optional quality/filter information from VCF
- **Gene**: Genomic region with name, chromosome, start/end positions, strand, and associated
  transcripts
- **Transcript**: Specific gene isoform with ID, exon structure, coding sequence coordinates,
  and canonical status
- **Annotation**: The predicted effect of a variant including consequence type(s), affected
  gene/transcript, and amino acid change if applicable
- **VEP Cache**: Pre-computed gene/transcript data from Ensembl in VEP's cache format

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Users can annotate a single-variant VCF file and receive correct output
  within 5 seconds (excluding cache loading time)
- **SC-002**: Annotation results for KRAS G12C match the output from Ensembl VEP
  for the same variant (consequence type, amino acid change, gene name)
- **SC-003**: 100% of variants in a test VCF are either successfully annotated or clearly
  marked as unannotated with a reason
- **SC-004**: Error messages for invalid input include specific remediation guidance
  (e.g., "Line 45: Expected 8 tab-separated columns, found 7")
- **SC-005**: Users with existing VEP cache files can use vibe-vep without downloading
  additional data

## Assumptions

- Users have already downloaded VEP cache files for GRCh38 human genome (setup documented separately)
- Input VCF files use GRCh38 coordinates (GRCh37 support is future scope)
- For this initial version, only SNVs and small indels (<50bp) are annotated; structural variants
  are explicitly out of scope
- The KRAS G12C variant (chr12:25245350 G>C, GRCh38) serves as the primary test case for validation
