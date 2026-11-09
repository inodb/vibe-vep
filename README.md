# vibe-vep

[![CI](https://github.com/inodb/vibe-vep/actions/workflows/ci.yml/badge.svg)](https://github.com/inodb/vibe-vep/actions/workflows/ci.yml)
[![Docs](https://github.com/inodb/vibe-vep/actions/workflows/docs.yml/badge.svg)](https://inodb.github.io/vibe-vep/)

A lightweight variant effect predictor for [cBioPortal](https://www.cbioportal.org/), [Genome Nexus](https://www.genomenexus.org/), and [OncoKB](https://www.oncokb.org/). Single Go binary, no Perl dependencies, ~95MB GENCODE download (vs 17GB VEP cache). Achieves 99.8% concordance with GDC/VEP annotations across 1M+ TCGA variants.

**[Documentation](https://inodb.github.io/vibe-vep/)**

## Install

```bash
go install github.com/inodb/vibe-vep/cmd/vibe-vep@latest
```

## Quick Start

```bash
# Download annotations (one-time)
vibe-vep download

# Annotate VCF
vibe-vep annotate vcf input.vcf

# Annotate MAF
vibe-vep annotate maf input.maf

# Convert VCF to MAF
vibe-vep convert vcf2maf input.vcf

# Validate MAF annotations
vibe-vep compare data_mutations.txt

# Serve the Genome Nexus API (no Genome Nexus service, no MongoDB)
vibe-vep serve --assembly GRCh37 --port 38080
```

## Features

- **VCF/MAF Parsing** — Annotate variants from VCF or MAF files
- **Consequence Prediction** — Sequence Ontology terms ([details](https://inodb.github.io/vibe-vep/docs/consequence-prediction/))
- **Annotation Sources** — AlphaMissense, ClinVar, Cancer Hotspots, SIGNAL, OncoKB ([details](https://inodb.github.io/vibe-vep/docs/annotation-sources/))
- **Validation** — 99.8% consequence match across 1M+ TCGA variants ([details](https://inodb.github.io/vibe-vep/docs/validation/))
- **Genome Nexus API** — `vibe-vep serve` answers it directly, so the annotation pipeline and other clients can be pointed at it by base URL alone ([details](https://inodb.github.io/vibe-vep/docs/annotation-server/))
- **Fast** — ~14,000 variants/sec parallel, ~5,000 single-threaded

## Documentation

Full documentation is available at **[inodb.github.io/vibe-vep](https://inodb.github.io/vibe-vep/)**:

- [Getting Started](https://inodb.github.io/vibe-vep/docs/getting-started/) — Installation, usage, configuration
- [Consequence Prediction](https://inodb.github.io/vibe-vep/docs/consequence-prediction/) — Supported consequences and how they work
- [Annotation Sources](https://inodb.github.io/vibe-vep/docs/annotation-sources/) — Data sources and match levels
- [Annotation Server](https://inodb.github.io/vibe-vep/docs/annotation-server/) — Serving the Genome Nexus API
- [Output Formats](https://inodb.github.io/vibe-vep/docs/output-formats/) — VCF and MAF output details
- [Validation](https://inodb.github.io/vibe-vep/docs/validation/) — Methodology and results
- [Development](https://inodb.github.io/vibe-vep/docs/development/) — Contributing and project structure

## License

MIT License
