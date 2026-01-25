# Quickstart: vibe-vep

**Feature**: 001-vcf-variant-annotation
**Date**: 2026-01-24

## Prerequisites

1. **VEP Cache Files** (required)
   - Download from Ensembl FTP: https://ftp.ensembl.org/pub/current_variation/vep/
   - Human GRCh38 cache: `homo_sapiens_vep_115_GRCh38.tar.gz`

   ```bash
   # Download and extract (example)
   cd ~/.vep
   curl -O https://ftp.ensembl.org/pub/release-115/variation/vep/homo_sapiens_vep_115_GRCh38.tar.gz
   tar xzf homo_sapiens_vep_115_GRCh38.tar.gz
   ```

2. **Convert cache to Sereal format** (recommended for Go)
   ```bash
   # Requires Perl VEP installation for conversion
   perl convert_cache.pl --sereal -species homo_sapiens -version 115
   ```

## Installation

```bash
# Build from source
git clone https://github.com/your-org/vibe-vep.git
cd vibe-vep
go build -o vibe-vep ./cmd/vibe-vep

# Or install directly
go install github.com/your-org/vibe-vep/cmd/vibe-vep@latest
```

## Basic Usage

### Annotate a VCF file

```bash
# Default tab-delimited output to stdout
vibe-vep annotate sample.vcf

# Save to file
vibe-vep annotate -o annotated.txt sample.vcf

# Gzipped input
vibe-vep annotate sample.vcf.gz
```

### Output Formats

```bash
# Tab-delimited (default, VEP-compatible)
vibe-vep annotate -f tab sample.vcf

# VCF with CSQ INFO field
vibe-vep annotate -f vcf -o annotated.vcf sample.vcf
```

### Custom Cache Location

```bash
vibe-vep annotate -c /data/vep_cache sample.vcf
```

## Validation Test: KRAS G12C

Create a test VCF file with the KRAS G12C variant:

```bash
cat > test_kras.vcf << 'EOF'
##fileformat=VCFv4.2
##reference=GRCh38
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
12	25245350	.	G	C	.	.	.
EOF
```

Run annotation:

```bash
vibe-vep annotate test_kras.vcf
```

**Expected Output:**
```
chr12_25245350_G/C	12:25245350	C	KRAS	ENST00000311936	Transcript	missense_variant	35	34	12	G/C	gGt/gTt	-	IMPACT=MODERATE;CANONICAL=YES
```

**Key fields to verify:**
- Gene: `KRAS`
- Consequence: `missense_variant`
- Protein_position: `12`
- Amino_acids: `G/C` (Glycine to Cysteine)
- CANONICAL: `YES`

## Pipeline Integration

```bash
# Filter for high-impact variants
vibe-vep annotate sample.vcf | grep -E "(stop_gained|frameshift|splice_donor|splice_acceptor)"

# Count variants by consequence
vibe-vep annotate sample.vcf | cut -f7 | sort | uniq -c | sort -rn

# Process large files with gzip
zcat large_cohort.vcf.gz | vibe-vep annotate - | gzip > annotated.txt.gz
```

## Troubleshooting

### Cache not found
```
Error: VEP cache not found at /home/user/.vep/homo_sapiens/115_GRCh38
```
**Solution:** Verify cache path and download if missing:
```bash
ls ~/.vep/homo_sapiens/
# Should show: 115_GRCh38/
```

### Invalid VCF format
```
Error: invalid VCF format at line 42: expected 8 columns, found 7
```
**Solution:** Check VCF file format. Ensure tab-separated with 8+ columns.

### Chromosome not found
```
Warning: chromosome 'chr12' not found, trying '12'
```
**Note:** vibe-vep automatically handles chr prefix differences.

## Performance Tips

1. Use gzipped VCF files for large datasets (automatic decompression)
2. For batch processing, consider splitting by chromosome
3. Use `--canonical` flag to reduce output size if only interested in one transcript per gene
