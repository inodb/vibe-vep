# vibe-vep

[![CI](https://github.com/inodb/vibe-vep/actions/workflows/ci.yml/badge.svg)](https://github.com/inodb/vibe-vep/actions/workflows/ci.yml)

An experiment in vibe coding a lightweight variant effect predictor for use by the [cBioPortal](https://www.cbioportal.org/), [Genome Nexus](https://www.genomenexus.org/), and [OncoKB](https://www.oncokb.org/) teams. Unlike Ensembl VEP, vibe-vep is a single Go binary with no Perl dependencies and uses a smaller annotation cache (~95MB GENCODE download vs 17GB VEP cache). It incorporates transcript prioritization similar to Genome Nexus — selecting a single gene and protein change by prioritizing coding transcripts and highest-impact consequences — while also providing effect predictions for all overlapping transcripts. Achieves 99.8% concordance with GDC/VEP annotations across 1M+ TCGA variants.

## Features

- **VCF/MAF Parsing**: Annotate variants from VCF or MAF files
- **GENCODE Annotations**: Uses GENCODE GTF/FASTA (~95MB download vs 17GB VEP cache)
- **Consequence Prediction**: Classifies variants using Sequence Ontology terms
- **Validation Mode**: Compare existing MAF annotations against predictions
- **Fast**: ~4,000 variants/sec end-to-end (parse + annotate + write), single-threaded (see [validation report](testdata/tcga/validation_report.md))

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
  config      Manage configuration (show/set/get)
  download    Download GENCODE annotation files
  help        Show help message

Annotate Options:
  --assembly      Genome assembly: GRCh37 or GRCh38 (default: GRCh38)
  --validate      Validate MAF annotations against predictions
  --validate-all  Show all variants in validation (default: mismatches only)
  -f, --output-format   Output format: vcf, maf (default: auto-detect from input)
  -o, --output    Output file (default: stdout)
  --canonical     Only report canonical transcript annotations

Download Options:
  --assembly      Genome assembly: GRCh37 or GRCh38 (default: GRCh38)
  --output        Output directory (default: ~/.vibe-vep/)
```

## Examples

### Annotating a VCF file

The repository includes a small [example VCF](examples/example.vcf) containing the KRAS G12C (p.Gly12Cys) hotspot mutation:

```bash
vibe-vep annotate examples/example.vcf
```

VCF input produces VCF output by default — the original VCF lines are preserved and a `CSQ` INFO field is added with consequence annotations for all overlapping transcripts (canonical flagged with `YES`):

```
##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from vibe-vep. Format: Allele|Consequence|IMPACT|SYMBOL|...">
#CHROM  POS       ID           REF  ALT  QUAL  FILTER  INFO
12      25245351  rs121913529  C    A    100   PASS    DP=50;CSQ=A|missense_variant|MODERATE|KRAS|...|YES,...
```

### Annotating a MAF file

The repository also includes a small [example MAF](examples/example.maf) with the same KRAS G12C variant:

```bash
vibe-vep annotate examples/example.maf
```

When the input is MAF, the output defaults to MAF format — all original columns are preserved and annotation columns are updated with fresh VEP predictions:

```
Hugo_Symbol  Entrez_Gene_Id  Center  ...  Consequence       Variant_Classification  ...  HGVSc    HGVSp        HGVSp_Short  Transcript_ID
KRAS         3845            .       ...  missense_variant  Missense_Mutation       ...  c.34G>T  p.Gly12Cys   p.G12C       ENST00000311936
```

### Other examples

```bash
# Validate TCGA MAF file
vibe-vep annotate --validate data_mutations.txt

# Show all validation results (not just mismatches)
vibe-vep annotate --validate --validate-all data_mutations.txt

# Annotate with GRCh37
vibe-vep annotate --assembly GRCh37 sample.vcf
```

## Consequence Prediction

### Supported Consequences

The tool predicts variant consequences using Sequence Ontology (SO) terms, ordered by impact:

| Impact | Consequences |
|--------|-------------|
| HIGH | stop_gained, frameshift_variant, stop_lost, start_lost, splice_donor_variant, splice_acceptor_variant |
| MODERATE | missense_variant, inframe_insertion, inframe_deletion |
| LOW | synonymous_variant, splice_region_variant, stop_retained_variant |
| MODIFIER | intron_variant, 5_prime_UTR_variant, 3_prime_UTR_variant, upstream/downstream_gene_variant, non_coding_transcript_exon_variant, mature_miRNA_variant |

### How Consequence Prediction Works

For each variant, the tool determines which transcripts overlap the variant position, then classifies the effect on each transcript:

1. **Upstream/Downstream**: Variant outside transcript boundaries
2. **Intronic**: Variant between exons. Checks for splice site overlap:
   - **Splice donor/acceptor** (HIGH): within ±1-2bp of exon boundary on intron side (donor at 5' end of intron, acceptor at 3' end, strand-aware)
   - **Splice region** (LOW): within 3bp exon side or 3-8bp intron side of splice junction
3. **UTR**: Variant in 5' or 3' untranslated region
4. **Coding**: Variant in CDS — calculates codon/amino acid change to determine missense, synonymous, stop_gained, etc.

For **indels**, the tool checks the entire deletion span (not just the start position) for:
- Splice site overlap across the full `[pos, pos+len(ref)-1]` range
- Start codon deletion (for deletions spanning from UTR into CDS)
- Stop codon overlap (frameshift at stop codon → frameshift_variant,stop_lost)
- In-frame insertions that create stop codons

### Transcript Biotype Handling

The tool treats transcripts as protein-coding if they have defined CDS coordinates (`CDSStart > 0 && CDSEnd > 0`), which covers:
- `protein_coding`
- `nonsense_mediated_decay` (appends `NMD_transcript_variant` modifier)
- `IG_*_gene` / `TR_*_gene` (immunoglobulin/T-cell receptor segments)
- `protein_coding_LoF`, `non_stop_decay`

For `miRNA` biotype transcripts, exonic variants are classified as `mature_miRNA_variant` instead of `non_coding_transcript_exon_variant`.

## Validation

### How Validation Works

The `--validate` flag compares existing MAF annotations against fresh predictions. For each MAF variant, the tool:

1. Parses the variant and its existing annotation (consequence, gene, transcript ID)
2. Re-annotates the variant against all overlapping GENCODE transcripts
3. Selects the best matching annotation for comparison (see below)
4. Normalizes consequence terms and compares

### Transcript Selection for Validation

When comparing against a MAF entry, the tool selects the best VEP annotation using this priority:

1. **Exact transcript ID match** — if the MAF specifies a transcript ID (e.g., ENST00000311936), use the annotation for that transcript. However, if the transcript's biotype has changed between GENCODE versions (e.g., was `protein_coding` but is now `retained_intron`) and the MAF has a coding consequence, the match is skipped.

2. **Same gene, best transcript** — among annotations for the same gene (by Hugo symbol), prefer:
   - Canonical transcripts (marked by GENCODE/Ensembl)
   - Protein-coding biotypes over non-coding
   - Higher-impact consequences (HIGH > MODERATE > LOW > MODIFIER)

3. **Any transcript, best match** — if no gene match, use the same ranking across all annotations.

### Consequence Normalization

MAF files may use different consequence terms than SO standard. The validation normalizes both sides before comparison:

- **MAF to SO mapping**: `Missense_Mutation` → `missense_variant`, `Silent` → `synonymous_variant`, etc.
- **Modifier stripping**: Drops terms like `non_coding_transcript_variant`, `NMD_transcript_variant`, `coding_sequence_variant`, `start_retained_variant`, `stop_retained_variant` that are secondary modifiers when a higher-impact term is present
- **Splice normalization**: Maps `splice_donor_region_variant` and `splice_donor_5th_base_variant` to `splice_region_variant`; drops `splice_region_variant` when a primary consequence is present
- **Impact-based stripping**: Drops `intron_variant` when splice donor/acceptor is present; drops UTR terms when a HIGH-impact term is present; drops `stop_gained`/`stop_lost` when co-occurring with `frameshift_variant`
- **Inframe grouping**: Normalizes `protein_altering_variant`, `inframe_deletion`, and `inframe_insertion` to a common term
- **Upstream/downstream tolerance**: MAF upstream/downstream calls are always accepted as matching, since different canonical transcript sets produce different transcript boundaries
- **Sorting**: Comma-separated terms are sorted alphabetically for consistent comparison

### Validation Results

Tested against 7 TCGA GDC studies from [cBioPortal/datahub](https://github.com/cBioPortal/datahub) (1,052,366 total variants):

| Column | Match | Mismatch | Rate |
|--------|-------|----------|------|
| Consequence | 1,050,079 | 46 | 99.8% |
| HGVSp | 997,292 | 282 | 94.8% |
| HGVSc | 1,037,748 | 278 | 98.6% |

Mismatches that are due to GENCODE version differences (not algorithm bugs) are reclassified into separate categories:
- **transcript_model_change** (501): transcript biotype changed between versions (e.g. protein_coding → retained_intron)
- **gene_model_change** (4): gene boundary differences (coding ↔ intergenic)
- **position_shift** (695 consequence, 4,733 HGVSp, 5,753 HGVSc): CDS coordinate changes between GENCODE versions

The validation also tracks mismatches in [OncoKB cancer genes](https://www.oncokb.org/cancerGenes) — only 9 cancer genes have any mismatches (1 each) across 1M+ variants.

See the full [validation report](testdata/tcga/validation_report.md) for per-study breakdowns and category details.

To download the TCGA test data and regenerate the report:

```bash
make download-testdata
go test ./internal/output/ -run TestValidationBenchmark -v -count=1
```

Remaining 46 consequence mismatches are primarily:
- CDS sequence differences between GENCODE versions (synonymous vs missense)
- Transcript structure differences (exon/intron boundary changes)
- Complex multi-region indels with ambiguous consequence priority

## Output Format

### VCF output (default for VCF input)

When the input is a VCF file, vibe-vep outputs VCF format by default. Original VCF header and data lines are preserved, and a `CSQ` INFO field is appended with consequence annotations for all overlapping transcripts. The CSQ format follows the VEP convention — pipe-delimited fields per transcript, comma-separated between transcripts:

```
CSQ=Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|CANONICAL
```

Canonical transcripts are flagged with `CANONICAL=YES`.

### MAF output (default for MAF input)

When the input is a MAF file, vibe-vep outputs MAF format by default. All original columns are preserved and annotation columns (`Hugo_Symbol`, `Consequence`, `Variant_Classification`, `HGVSp`, `HGVSp_Short`, `HGVSc`, `Transcript_ID`) are updated with fresh VEP predictions. Original values are kept when VEP has no prediction (e.g., intergenic variants).

When an [OncoKB cancer gene list](https://www.oncokb.org/cancerGenes) is configured, a `Gene_Type` column is appended with oncogene/TSG classification (see [Configuration](#configuration)).

## Performance

- **End-to-end throughput**: ~4,000 variants/sec (MAF parse + annotate + compare + write, single-threaded)
- **Cache loading**: ~25 seconds to load 254k GENCODE v46 transcripts
- **Memory**: Proportional to transcript count (~254k for GENCODE v46)
- **Total benchmark**: 1,052,366 variants across 7 TCGA studies in ~4.5 minutes

## Configuration

vibe-vep reads configuration from `~/.vibe-vep.yaml` (or `--config` flag). Environment variables with `VIBE_VEP_` prefix also work (e.g. `VIBE_VEP_ONCOKB_CANCER_GENE_LIST`).

```yaml
# OncoKB cancer gene list — adds Gene_Type column to MAF output
oncokb:
  cancer-gene-list: /path/to/cancerGeneList.tsv

# Optional annotation sources
annotations:
  alphamissense: true  # AlphaMissense pathogenicity scores (CC BY 4.0)
```

The `cancerGeneList.tsv` file can be downloaded from [OncoKB](https://www.oncokb.org/cancerGenes) or is included in this repository.

Use `vibe-vep config set annotations.alphamissense true` to enable AlphaMissense. Then run `vibe-vep download` to fetch the data (~643MB) and `vibe-vep prepare` to load it into the annotation database.

## Data Sources

- **GENCODE**: Gene annotations from [GENCODE](https://www.gencodegenes.org/)
- **OncoKB**: Cancer gene classification from [OncoKB](https://www.oncokb.org/) (optional, for gene-level ONCOGENE/TSG annotations)
- **AlphaMissense**: Pathogenicity predictions from [AlphaMissense](https://github.com/google-deepmind/alphamissense) (Cheng et al., Science 2023). Data licensed under CC BY 4.0.

## Development

```bash
# Run tests
go test ./...

# Run benchmarks
go test ./internal/annotate/ -bench .

# Build
go build -o vibe-vep ./cmd/vibe-vep

# Download TCGA test data for validation (~1.6GB)
make download-testdata

# Run validation against TCGA data
vibe-vep annotate --validate testdata/tcga/chol_tcga_gdc_data_mutations.txt
```

## Roadmap

- [ ] **Feature parity for MAF annotation** — Match the annotation capabilities of the [genome-nexus-annotation-pipeline](https://github.com/genome-nexus/genome-nexus-annotation-pipeline), which currently relies on Genome Nexus Server + VEP
  - [x] Consequence prediction (99.8% concordance with GDC/VEP across 1M+ TCGA variants, see [Validation Results](#validation-results))
  - [x] HGVSp/HGVSc notation (94.8%/98.6% match rates)
  - [x] Full MAF output format (all required columns)
  - [x] Cancer gene annotations (OncoKB oncogene/TSG classification)
- [ ] **Additional annotation sources** — Pathogenicity predictions and external databases
  - [x] AlphaMissense pathogenicity scores (CC BY 4.0, Google DeepMind)
  - [ ] SIFT predictions (tolerated/deleterious)
  - [ ] PolyPhen-2 predictions (benign/possibly_damaging/probably_damaging)
  - [ ] gnomAD allele frequencies
  - [ ] ClinVar clinical significance
- [ ] **Re-annotate datahub GDC studies** — Validate by re-annotating [cBioPortal/datahub](https://github.com/cBioPortal/datahub) GDC studies with vibe-vep
- [ ] **Replace genome-nexus-annotation-pipeline for datahub** — Use vibe-vep as the annotation tool for datahub processing

## License

MIT License
