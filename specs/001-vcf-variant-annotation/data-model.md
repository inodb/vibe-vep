# Data Model: VCF Parsing and Variant Annotation

**Feature**: 001-vcf-variant-annotation
**Date**: 2026-01-24

## Core Entities

### Variant

Represents a single genomic change from a VCF file.

| Field | Type | Description | Validation |
|-------|------|-------------|------------|
| Chrom | string | Chromosome name (e.g., "12", "chr12") | Required, must be valid chromosome |
| Pos | int64 | 1-based genomic position | Required, > 0 |
| ID | string | Variant identifier (e.g., rs ID) | Optional |
| Ref | string | Reference allele | Required, non-empty, valid nucleotides |
| Alt | string | Alternate allele | Required, non-empty, valid nucleotides |
| Qual | float64 | Quality score | Optional |
| Filter | string | Filter status (PASS or filter name) | Optional |
| Info | map[string]interface{} | INFO field key-value pairs | Optional |

**Derived Properties:**
- `IsSNV() bool`: len(Ref) == 1 && len(Alt) == 1
- `IsIndel() bool`: len(Ref) != len(Alt)
- `IsInsertion() bool`: len(Alt) > len(Ref)
- `IsDeletion() bool`: len(Ref) > len(Alt)

---

### Gene

Represents a genomic region with associated transcripts.

| Field | Type | Description | Validation |
|-------|------|-------------|------------|
| ID | string | Gene identifier (e.g., ENSG00000133703) | Required |
| Name | string | Gene symbol (e.g., KRAS) | Required |
| Chrom | string | Chromosome | Required |
| Start | int64 | Gene start position (1-based) | Required, > 0 |
| End | int64 | Gene end position (1-based, inclusive) | Required, >= Start |
| Strand | int8 | +1 (forward) or -1 (reverse) | Required, +1 or -1 |
| Biotype | string | Gene biotype (e.g., protein_coding) | Required |
| Transcripts | []Transcript | Associated transcripts | May be empty |

---

### Transcript

Represents a specific gene isoform.

| Field | Type | Description | Validation |
|-------|------|-------------|------------|
| ID | string | Transcript ID (e.g., ENST00000311936) | Required |
| GeneID | string | Parent gene ID | Required |
| GeneName | string | Parent gene symbol | Required |
| Chrom | string | Chromosome | Required |
| Start | int64 | Transcript start (1-based) | Required, > 0 |
| End | int64 | Transcript end (1-based, inclusive) | Required, >= Start |
| Strand | int8 | +1 or -1 | Required |
| Biotype | string | Transcript biotype | Required |
| IsCanonical | bool | Ensembl canonical flag | Required |
| IsMANESelect | bool | MANE Select transcript | Required (human only) |
| Exons | []Exon | Ordered exons | Required for protein_coding |
| CDSStart | int64 | CDS start (genomic, 1-based) | 0 if non-coding |
| CDSEnd | int64 | CDS end (genomic, 1-based) | 0 if non-coding |
| CDSSequence | string | Coding DNA sequence | Loaded on demand |
| ProteinSequence | string | Translated protein sequence | Loaded on demand |

**Relationships:**
- Many-to-one with Gene
- One-to-many with Exon

---

### Exon

Represents a single exon within a transcript.

| Field | Type | Description | Validation |
|-------|------|-------------|------------|
| Number | int | Exon number (1-based) | Required, > 0 |
| Start | int64 | Genomic start (1-based) | Required, > 0 |
| End | int64 | Genomic end (1-based, inclusive) | Required, >= Start |
| CDSStart | int64 | CDS portion start | 0 if entirely non-coding |
| CDSEnd | int64 | CDS portion end | 0 if entirely non-coding |
| Frame | int | Reading frame (0, 1, or 2) | -1 if non-coding |

**Notes:**
- Exons are ordered by Number, not genomic position
- For reverse strand, genomic Start < End but CDS positions count backward

---

### Annotation

Represents the predicted effect of a variant on a transcript.

| Field | Type | Description | Validation |
|-------|------|-------------|------------|
| VariantID | string | Source variant identifier | Required |
| TranscriptID | string | Affected transcript | Required |
| GeneName | string | Gene symbol | Required |
| GeneID | string | Gene identifier | Required |
| Consequence | string | SO consequence term | Required |
| ConsequenceID | string | SO accession (e.g., SO:0001583) | Required |
| Impact | string | HIGH, MODERATE, LOW, MODIFIER | Required |
| CDSPosition | int64 | Position in CDS | 0 if not in CDS |
| ProteinPosition | int64 | Amino acid position | 0 if not in CDS |
| AminoAcidChange | string | e.g., "G12C" | Empty if not missense |
| CodonChange | string | e.g., "GGT/TGT" | Empty if not coding |
| IsCanonical | bool | Annotation on canonical transcript | Required |

**Impact Levels:**
- HIGH: stop_gained, frameshift_variant, splice_acceptor/donor_variant, start/stop_lost
- MODERATE: missense_variant, inframe_insertion, inframe_deletion
- LOW: synonymous_variant, splice_region_variant
- MODIFIER: intron_variant, UTR variants, intergenic_variant

---

### VEPCache

Represents the VEP cache metadata and access.

| Field | Type | Description | Validation |
|-------|------|-------------|------------|
| Species | string | Species name (e.g., homo_sapiens) | Required |
| Version | int | Ensembl release version | Required, > 0 |
| Assembly | string | Genome assembly (e.g., GRCh38) | Required |
| CacheDir | string | Path to cache directory | Required, must exist |
| Format | string | "storable", "sereal", or "gtf" | Required |

---

## Entity Relationships

```
┌─────────────┐
│   VCF File  │
└──────┬──────┘
       │ parses to
       ▼
┌─────────────┐     annotates      ┌─────────────┐
│   Variant   │◄──────────────────►│  Annotation │
└─────────────┘                    └──────┬──────┘
                                          │ references
                                          ▼
┌─────────────┐     contains       ┌─────────────┐
│    Gene     │◄──────────────────►│ Transcript  │
└─────────────┘                    └──────┬──────┘
                                          │ contains
                                          ▼
                                   ┌─────────────┐
                                   │    Exon     │
                                   └─────────────┘
```

---

## State Transitions

### Variant Processing State

```
┌─────────┐     parse      ┌──────────┐     lookup      ┌───────────┐
│ Raw VCF │───────────────►│  Parsed  │────────────────►│ Annotated │
│  Line   │                │  Variant │                 │  Variant  │
└─────────┘                └──────────┘                 └───────────┘
     │                           │                            │
     │ invalid format            │ no overlapping             │ success
     ▼                           │ transcripts                ▼
┌─────────┐                      ▼                      ┌───────────┐
│  Error  │                ┌───────────┐                │  Output   │
│  State  │                │Intergenic │                │  Written  │
└─────────┘                └───────────┘                └───────────┘
```

---

## Validation Rules

### Variant Validation
1. Chromosome must be recognized (chr1-22, chrX, chrY, chrM, or without 'chr' prefix)
2. Position must be positive integer within chromosome bounds
3. Ref and Alt must contain only valid nucleotides (A, C, G, T, N)
4. Ref must not equal Alt

### Transcript Validation
1. CDSStart must be >= Start and <= End
2. CDSEnd must be >= CDSStart and <= End
3. Exons must not overlap
4. Exon numbers must be sequential starting from 1

### Annotation Validation
1. Consequence must be valid SO term
2. Impact must match consequence type
3. If consequence is missense_variant, AminoAcidChange must be non-empty
4. CDSPosition and ProteinPosition must be consistent: `ProteinPosition = (CDSPosition - 1) / 3 + 1`

---

## Indexing Strategy

### Transcript Lookup
- Primary index: Chromosome + Position interval tree
- Query: Find all transcripts overlapping variant position
- Performance: O(log n + k) where n = transcripts, k = matches

### Exon Lookup
- Per-transcript index: Position to exon mapping
- Precompute CDS position offsets for each exon
- Cache CDSSequence on first access per transcript
