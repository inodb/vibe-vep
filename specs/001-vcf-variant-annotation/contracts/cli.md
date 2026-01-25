# CLI Contract: vibe-vep

**Feature**: 001-vcf-variant-annotation
**Date**: 2026-01-24

## Command Structure

```
vibe-vep <command> [options] <input>
```

## Commands

### annotate

Annotate variants in a VCF file with consequence predictions.

```
vibe-vep annotate [options] <input.vcf>
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `<input.vcf>` | Input VCF file path (supports .vcf and .vcf.gz) |

#### Options

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--cache-dir` | `-c` | `~/.vep` | Path to VEP cache directory |
| `--species` | `-s` | `homo_sapiens` | Species name |
| `--assembly` | `-a` | `GRCh38` | Genome assembly |
| `--output` | `-o` | stdout | Output file path |
| `--output-format` | `-f` | `tab` | Output format: `tab` or `vcf` |
| `--canonical` | | false | Only report canonical transcript annotations |
| `--help` | `-h` | | Show help message |

#### Examples

```bash
# Basic annotation with defaults
vibe-vep annotate sample.vcf

# Annotate gzipped VCF, output to file
vibe-vep annotate -o annotated.txt sample.vcf.gz

# VCF output format with CSQ field
vibe-vep annotate -f vcf -o annotated.vcf sample.vcf

# Specify custom cache directory
vibe-vep annotate -c /data/vep_cache sample.vcf

# Only canonical transcripts
vibe-vep annotate --canonical sample.vcf
```

---

## Output Formats

### Tab-delimited Format (default)

VEP-compatible tab-delimited output.

**Header:**
```
## vibe-vep output
## Command: vibe-vep annotate sample.vcf
## VEP cache: homo_sapiens/115_GRCh38
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	Extra
```

**Columns:**

| Column | Description | Example |
|--------|-------------|---------|
| Uploaded_variation | Variant ID or generated ID | chr12_25245350_G/C |
| Location | Chromosome:position | 12:25245350 |
| Allele | Alternate allele | C |
| Gene | Gene symbol | KRAS |
| Feature | Transcript ID | ENST00000311936 |
| Feature_type | Feature type | Transcript |
| Consequence | SO consequence term | missense_variant |
| cDNA_position | Position in cDNA | 35 |
| CDS_position | Position in CDS | 34 |
| Protein_position | Amino acid position | 12 |
| Amino_acids | Reference/Alternate AA | G/C |
| Codons | Reference/Alternate codon | gGt/gTt |
| Existing_variation | Known variant IDs | - |
| Extra | Additional annotations | IMPACT=MODERATE;CANONICAL=YES |

**Example Output:**
```
chr12_25245350_G/C	12:25245350	C	KRAS	ENST00000311936	Transcript	missense_variant	35	34	12	G/C	gGt/gTt	-	IMPACT=MODERATE;CANONICAL=YES
```

### VCF Format

VCF output with CSQ INFO field containing annotations.

**CSQ Header:**
```
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from vibe-vep. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|CANONICAL">
```

**Example:**
```
12	25245350	.	G	C	.	.	CSQ=C|missense_variant|MODERATE|KRAS|ENSG00000133703|Transcript|ENST00000311936|protein_coding|2/5||35|34|12|G/C|gGt/gTt|YES
```

---

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | General error (file not found, parse error, etc.) |
| 2 | Invalid usage (bad arguments, missing required options) |

---

## Error Messages

### File Not Found
```
Error: input file not found: /path/to/missing.vcf
```

### Invalid VCF Format
```
Error: invalid VCF format at line 42: expected 8 columns, found 7
```

### Cache Not Found
```
Error: VEP cache not found at /home/user/.vep/homo_sapiens/115_GRCh38
Hint: Download cache from https://ftp.ensembl.org/pub/current_variation/vep/
```

### Unknown Chromosome
```
Warning: chromosome 'chrUn_gl000220' not found in cache, skipping 3 variants
```

---

## Stdin/Stdout Support

```bash
# Read from stdin
cat sample.vcf | vibe-vep annotate -

# Pipe output
vibe-vep annotate sample.vcf | grep "missense_variant"

# Full pipeline
zcat sample.vcf.gz | vibe-vep annotate - | gzip > annotated.txt.gz
```
