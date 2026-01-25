# Research: VCF Parsing and Variant Annotation

**Feature**: 001-vcf-variant-annotation
**Date**: 2026-01-24

## 1. VCF Parsing Library

### Decision
Use **brentp/vcfgo** with **biogo/hts** for BGZF support.

### Rationale
- Battle-tested in production bioinformatics tools (vcfanno, mosdepth)
- True streaming with lazy sample parsing for performance
- MIT license (permissive)
- `SplitAlts()` function for multi-allelic decomposition
- Active maintenance (September 2025 updates)
- Handles malformed VCF files gracefully with error reporting

### Alternatives Considered
| Library | Pros | Cons | Decision |
|---------|------|------|----------|
| vertgenlab/gonomics/vcf | Channel-based streaming, strict VCF 4.3 | Heavier dependency, less adopted | Backup option |
| exascience/elprep/vcf | Exceptional performance | AGPL-3.0 license (copyleft) | License incompatible |
| mendelics/vcf | Simple API | Archived, no gzip support | Rejected |

### Implementation Notes
```go
import (
    "github.com/brentp/vcfgo"
    "github.com/biogo/hts/bgzf"
)

// For gzipped VCF with BGZF
f, _ := os.Open("variants.vcf.gz")
br, _ := bgzf.NewReader(f, 0)  // Parallel decompression
rdr, _ := vcfgo.NewReader(br, true)  // lazySamples=true
```

---

## 2. VEP Cache Format

### Decision
Use **Sereal format** (converted from Perl Storable) for transcript data, with option to
support GTF/GFF as alternative data source.

### Rationale
- VEP cache uses Perl Storable format which is Perl-specific
- Sereal has official Go implementation: `github.com/Sereal/Sereal/Go/sereal`
- One-time conversion: `perl convert_cache.pl --sereal -species homo_sapiens`
- Tabix-indexed variation files (`all_vars.gz`) can be read with `biogo/hts/tabix`

### VEP Cache Directory Structure
```
~/.vep/homo_sapiens/115_GRCh38/
├── info.txt                 # Cache metadata
├── chr_synonyms.txt         # Chromosome aliases
├── 12/                      # Chromosome directory
│   ├── 25000001-26000000.gz # Transcript data (Storable/Sereal)
│   ├── all_vars.gz          # Tabix-indexed variants
│   └── all_vars.gz.tbi      # Tabix index
```

### Transcript Data Fields Required
- `stable_id`: Transcript ID (e.g., ENST00000311936)
- `_gene_stable_id`, `_gene_symbol`: Gene identifiers
- `start`, `end`, `strand`: Genomic coordinates
- `biotype`: protein_coding, lncRNA, etc.
- `is_canonical`: Canonical transcript flag
- `_trans_exon_array`: Exon coordinates with phase information
- `cdna_coding_start`, `cdna_coding_end`: CDS boundaries
- `peptide`, `translateable_seq`: Protein and CDS sequences

### Alternatives Considered
| Approach | Pros | Cons | Decision |
|----------|------|------|----------|
| Sereal format | Official Go library, cross-language | Requires one-time conversion | Primary |
| GTF/GFF files | Standard format, Go parsers exist | Less data than cache | Fallback |
| Custom format | Optimized for Go | Development overhead | Future option |

---

## 3. Consequence Prediction Algorithm

### Decision
Implement hierarchical predicate-based consequence classification using Sequence Ontology
terms, matching VEP's approach.

### Consequence Types by Impact Level

**HIGH Impact:**
- `stop_gained` (SO:0001587): Premature stop codon
- `frameshift_variant` (SO:0001589): Indel length not divisible by 3
- `stop_lost` (SO:0001578): Stop codon removed
- `start_lost` (SO:0002012): Start codon removed
- `splice_donor_variant` (SO:0001575): GT dinucleotide at +1/+2
- `splice_acceptor_variant` (SO:0001574): AG dinucleotide at -1/-2

**MODERATE Impact:**
- `missense_variant` (SO:0001583): Different amino acid
- `inframe_insertion` (SO:0001821): Insertion length divisible by 3
- `inframe_deletion` (SO:0001822): Deletion length divisible by 3

**LOW Impact:**
- `synonymous_variant` (SO:0001819): Same amino acid
- `splice_region_variant` (SO:0001630): Within 3bp exon / 8bp intron of splice site

**MODIFIER Impact:**
- `5_prime_UTR_variant` (SO:0001623)
- `3_prime_UTR_variant` (SO:0001624)
- `intron_variant` (SO:0001627)
- `intergenic_variant` (SO:0001628)

### Decision Tree
1. Does variant overlap any transcript? NO → intergenic_variant
2. Within transcript bounds? NO → upstream/downstream_gene_variant (within 5kb)
3. In exon? NO → Check splice sites, then intron_variant
4. Protein-coding transcript? NO → non_coding_transcript_exon_variant
5. Within CDS? NO → 5'/3' UTR variant
6. Calculate codon change → Classify coding impact

### Coordinate Mapping
1. Genomic → Transcript position (accounting for exon structure)
2. Transcript → CDS position (subtract 5' UTR length)
3. CDS position → Codon number: `(cds_pos - 1) / 3 + 1`
4. Position in codon: `(cds_pos - 1) % 3` (0, 1, or 2)

### KRAS G12C Verification
- Position: chr12:25245350 G>C (GRCh38)
- KRAS is on **reverse strand**
- Maps to CDS position 34 → Codon 12, position 1
- Reference codon: GGT (Glycine)
- Alternate codon: TGT (Cysteine)
- Result: `missense_variant`, p.Gly12Cys (G12C)

---

## 4. Canonical Transcript Selection

### Decision
Follow Ensembl's priority order: MANE Select > MANE Plus Clinical > Ensembl Canonical.

### Rationale
- MANE Select provides one transcript per protein-coding gene for human
- Identical between RefSeq and GENCODE/Ensembl
- Perfectly aligned to GRCh38
- Standard for clinical annotation

### Fallback Algorithm (when no MANE)
1. APPRIS Principal isoform (highest score)
2. UniProt canonical isoform
3. Longest CDS
4. protein_coding biotype preference

---

## 5. Key Implementation Considerations

### Strand Handling
KRAS and many genes are on reverse strand. For reverse strand genes:
- Reference allele in VCF is forward strand
- Must reverse-complement to get coding strand sequence
- Exon order in transcript is reversed relative to genomic coordinates

### Multi-allelic Sites
VCF multi-allelic sites (multiple ALT alleles) should be split and annotated separately
using `vcfgo.SplitAlts()`.

### Performance Optimizations
1. Pre-predicate filtering: Skip expensive calculations for intronic variants
2. Cache transcript sequences per chromosome
3. Use interval trees for fast transcript lookup
4. Process variants in batches

---

## 6. Go Dependencies Summary

| Package | Purpose | License |
|---------|---------|---------|
| github.com/brentp/vcfgo | VCF parsing | MIT |
| github.com/biogo/hts | BGZF, tabix | BSD-3 |
| github.com/Sereal/Sereal/Go/sereal | VEP cache reading | Artistic-2.0 |

---

## Sources

- [Ensembl VEP Documentation](http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html)
- [VEP Cache Structure](https://www.ensembl.info/2020/10/23/cool-stuff-ensembl-vep-can-do-whats-in-the-cache-and-how-does-vep-use-it/)
- [Ensembl Calculated Consequences](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html)
- [The Ensembl VEP (Genome Biology)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4893825/)
- [Sequence Ontology](http://www.sequenceontology.org/)
- [Ensembl Canonical Transcripts](https://www.ensembl.org/info/genome/genebuild/canonical.html)
- [brentp/vcfgo](https://github.com/brentp/vcfgo)
- [biogo/hts](https://github.com/biogo/hts)
- [Sereal Go](https://github.com/Sereal/Sereal/tree/master/Go/sereal)
