# vibe-vep

A fast, lightweight variant effect predictor for genomic annotation.

## Features

- **VCF/MAF Parsing**: Annotate variants from VCF or MAF files
- **GENCODE Annotations**: Uses GENCODE GTF/FASTA (~95MB download vs 17GB VEP cache)
- **Consequence Prediction**: Classifies variants using Sequence Ontology terms
- **Validation Mode**: Compare existing MAF annotations against predictions
- **REST API Fallback**: On-demand fetching from Ensembl REST API

## Installation

```bash
go install github.com/inodb/vibe-vep/cmd/vibe-vep@latest
```

Or build from source:

```bash
git clone https://github.com/inodb/vibe-vep.git
cd vibe-vep
go build -o vibe-vep ./cmd/vibe-vep
```

## Quick Start

### 1. Download GENCODE Annotations (one-time setup)

```bash
# For GRCh38 (default)
vibe-vep download --assembly GRCh38

# For GRCh37
vibe-vep download --assembly GRCh37
```

This downloads ~95MB of annotation data to `~/.vibe-vep/`.

### 2. Annotate Variants

```bash
# Annotate a VCF file
vibe-vep annotate input.vcf

# Annotate a MAF file
vibe-vep annotate input.maf

# Output to file
vibe-vep annotate -o output.txt input.vcf
```

### 3. Validate MAF Annotations

Compare existing MAF annotations against VEP predictions:

```bash
vibe-vep annotate --validate data_mutations.txt
```

## Usage

```
vibe-vep - Variant Effect Predictor

Commands:
  annotate    Annotate variants in a VCF or MAF file
  download    Download GENCODE annotation files
  convert     Convert VEP cache to DuckDB format
  help        Show help message

Annotate Options:
  --assembly      Genome assembly: GRCh37 or GRCh38 (default: GRCh38)
  --cache         Cache file (DuckDB) or directory
  --rest          Use Ensembl REST API (no local cache needed)
  --validate      Validate MAF annotations against predictions
  --validate-all  Show all variants in validation (default: mismatches only)
  -f, --output-format   Output format: tab, vcf (default: tab)
  -o, --output    Output file (default: stdout)
  --canonical     Only report canonical transcript annotations

Download Options:
  --assembly      Genome assembly: GRCh37 or GRCh38 (default: GRCh38)
  --output        Output directory (default: ~/.vibe-vep/)
```

## Examples

```bash
# Basic annotation
vibe-vep annotate sample.vcf

# Use REST API (slower, but no download needed)
vibe-vep annotate --rest sample.vcf

# Validate TCGA MAF file
vibe-vep annotate --validate data_mutations.txt

# Show all validation results (not just mismatches)
vibe-vep annotate --validate --validate-all data_mutations.txt

# Annotate with GRCh37
vibe-vep annotate --assembly GRCh37 sample.vcf
```

## Output Format

Default tab-delimited output includes:

| Column | Description |
|--------|-------------|
| #Uploaded_variation | Variant identifier (chr:pos ref>alt) |
| Location | Genomic location |
| Allele | Alternate allele |
| Gene | Gene symbol |
| Feature | Transcript ID |
| Consequence | SO term (e.g., missense_variant) |
| cDNA_position | Position in cDNA |
| CDS_position | Position in CDS |
| Protein_position | Position in protein |
| Amino_acids | Amino acid change (e.g., G/C) |
| Codons | Codon change |
| IMPACT | Impact level (HIGH, MODERATE, LOW, MODIFIER) |

## Data Sources

- **GENCODE**: Gene annotations from [GENCODE](https://www.gencodegenes.org/)
- **Ensembl REST API**: On-demand transcript data from [Ensembl](https://rest.ensembl.org/)

## Development

```bash
# Run tests
go test ./...

# Build
go build -o vibe-vep ./cmd/vibe-vep

# Run validation against TCGA data
vibe-vep annotate --validate /path/to/data_mutations.txt
```

## License

MIT License
