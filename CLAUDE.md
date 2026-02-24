# vibe-vep Development Guidelines

## Tech Stack

- Go 1.22+

## Project Structure

```
cmd/vibe-vep/       CLI entry point (annotate, download commands)
internal/
  annotate/         Consequence prediction (PredictConsequence, Annotator)
  cache/            Transcript cache (GENCODE GTF/FASTA loader)
  maf/              MAF file parser
  output/           Output formatting and validation comparison
  vcf/              VCF file parser
testdata/
  cache/            Test transcript data (JSON)
  tcga/             TCGA GDC MAF files for validation
```

## Commands

```bash
# Run tests
go test ./...

# Run benchmarks
go test ./internal/annotate/ -bench . -benchmem

# Build
go build -o vibe-vep ./cmd/vibe-vep

# Run validation against TCGA data
go run ./cmd/vibe-vep annotate --validate testdata/tcga/chol_tcga_gdc_data_mutations.txt
```

## Key Design Decisions

- GENCODE GTF/FASTA is the only cache backend (no DuckDB, REST API, or Sereal)
- Consequence prediction follows Sequence Ontology terms with VEP-compatible output
- Transcript prioritization: canonical > protein-coding biotype > highest impact
- Validation normalizes MAF and SO consequence terms before comparison
- Performance target: >100k variants/sec (currently ~720k/sec)
