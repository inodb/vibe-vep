---
title: Architecture
description: Technical architecture of vibe-vep.
weight: 2
aliases:
  - /docs/architecture/
---

## Overview

vibe-vep is a single Go binary that annotates cancer variants with consequence predictions, protein effects, and clinical annotations. It reads VCF/MAF files, predicts consequences against GENCODE transcripts, and enriches variants with data from multiple annotation sources.

```mermaid
graph LR
    subgraph Input
        VCF[VCF file]
        MAF[MAF file]
    end

    subgraph Core["Annotation Engine"]
        Parser["Parser<br/><small>internal/maf, internal/vcf</small>"]
        Cache["Transcript Cache<br/><small>254K transcripts</small>"]
        Annotator["Annotator<br/><small>internal/annotate</small>"]
        Sources["Annotation Sources<br/><small>internal/datasource/*</small>"]
    end

    subgraph Output
        MAFOUT["Annotated MAF"]
        VCFOUT["Annotated VCF"]
        Parquet["Parquet export"]
    end

    VCF --> Parser
    MAF --> Parser
    Parser --> Annotator
    Cache --> Annotator
    Annotator --> Sources
    Sources --> MAFOUT
    Sources --> VCFOUT
    Sources --> Parquet
```

## Annotation Pipeline

Each variant flows through a three-stage pipeline: parsing, consequence prediction, and source enrichment.

```mermaid
flowchart TB
    subgraph parse ["1. Parse"]
        direction LR
        P1["Read variant<br/>chr, pos, ref, alt"]
        P2["Normalize coordinates<br/>strip chr prefix,<br/>VCF→MAF indel style"]
        P1 --> P2
    end

    subgraph predict ["2. Predict Consequences"]
        direction TB
        T1["Find overlapping transcripts<br/><small>interval tree, O(log n)</small>"]
        T2["Classify per transcript:<br/>upstream/downstream,<br/>intronic/splice, UTR,<br/>coding (missense, fs, etc.)"]
        T3["Select best transcript<br/><small>canonical > protein-coding > highest impact</small>"]
        T1 --> T2 --> T3
    end

    subgraph enrich ["3. Enrich"]
        direction LR
        S1["Genomic sources<br/><small>AlphaMissense, ClinVar, SIGNAL,<br/>gnomAD, SIFT, PolyPhen-2, dbSNP</small><br/><small>SQLite point lookup, ~1-5 &mu;s</small>"]
        S2["Protein sources<br/><small>Cancer Hotspots</small><br/><small>transcript + AA position</small>"]
        S3["Gene sources<br/><small>OncoKB cancer gene list</small>"]
    end

    parse --> predict --> enrich
```

## Data Flow

### Startup

On first run, `vibe-vep download` fetches GENCODE annotations. Subsequent runs load from a gob-serialized transcript cache for fast startup (~10s for 254K transcripts vs ~50s from raw GTF/FASTA).

```mermaid
flowchart LR
    subgraph download ["vibe-vep download (one-time)"]
        DL1["GENCODE GTF<br/><small>49 MB gz</small>"]
        DL2["GENCODE FASTA<br/><small>47 MB gz</small>"]
        DL3["Canonical overrides<br/><small>26 MB</small>"]
    end

    subgraph prepare ["vibe-vep prepare (one-time)"]
        DL4["AlphaMissense TSV<br/><small>614 MB gz</small>"]
        DL5["ClinVar VCF<br/><small>182 MB gz</small>"]
        SQLite["genomic_annotations.sqlite<br/><small>3.8 GB, WITHOUT ROWID</small>"]
        DL4 --> SQLite
        DL5 --> SQLite
    end

    subgraph startup ["vibe-vep annotate (each run)"]
        GOB["transcripts.gob<br/><small>210 MB, ~10s load</small>"]
        IDX["Interval tree index<br/><small>in-memory, O(log n) lookup</small>"]
        GOB --> IDX
    end

    DL1 --> GOB
    DL2 --> GOB
    DL3 --> GOB
```

### Per-Variant

```mermaid
sequenceDiagram
    participant P as Parser
    participant A as Annotator
    participant T as Transcript Cache
    participant G as Genomic Index
    participant W as Writer

    P->>A: variant (chr, pos, ref, alt)
    A->>T: FindOverlapping(chr, pos, end)
    T-->>A: []Transcript
    loop Each transcript
        A->>A: PredictConsequence(variant, transcript)
    end
    A->>A: SelectBestAnnotation (canonical > coding > impact)
    A->>G: Lookup(chr, pos, ref, alt)
    G-->>A: AlphaMissense, ClinVar, Hotspots, ...
    A->>W: Write annotated row
```

## Storage Architecture

All data lives under `~/.vibe-vep/{assembly}/`:

```
~/.vibe-vep/grch38/
  gencode.v46.annotation.gtf.gz         # 49 MB  - GENCODE gene models (source)
  gencode.v46.pc_transcripts.fa.gz      # 47 MB  - protein-coding sequences (source)
  ensembl_biomart_canonical_*.txt        # 26 MB  - canonical transcript overrides (source)
  transcripts.gob                        # 210 MB - serialized transcript cache
  genomic_annotations.sqlite             # 3.8 GB+ - AlphaMissense + ClinVar + SIGNAL + gnomAD + SIFT/PP2 + dbSNP
  variant_cache.duckdb                   # var.    - annotation result cache
```

### Storage Engines

| Engine | File | Purpose | Access Pattern |
|--------|------|---------|---------------|
| **Gob** | `transcripts.gob` | 254K transcripts with exon/CDS/protein data | Full deserialize on startup |
| **SQLite** | `genomic_annotations.sqlite` | 75M+ annotation records (AM + ClinVar + SIGNAL + gnomAD + SIFT/PP2 + dbSNP) | mmap point lookups, ~1-5 &mu;s |
| **DuckDB** | `variant_cache.duckdb` | Previously annotated results + Parquet export | Columnar read/write, batch queries |

### SQLite vs DuckDB: Why Two Databases?

**SQLite** (`genomic_annotations.sqlite`) is a **read-only reference database** used *during* annotation. It holds pre-computed data from external sources (AlphaMissense scores, ClinVar significance, SIGNAL frequencies, gnomAD allele frequencies, SIFT/PolyPhen-2 predictions, dbSNP RS IDs) and is queried once per variant via a point lookup on `(chrom, pos, ref, alt)`. Built once by `vibe-vep prepare`, it uses a `WITHOUT ROWID` clustered primary key for single B-tree traversals served directly from mmap with near-zero Go heap allocation (~1-5 &mu;s per lookup).

**DuckDB** (`variant_cache.duckdb`) stores **annotation results** written *after* annotation. It serves two purposes:

1. **Variant result cache** (`--save-results`): Stores completed annotations so re-annotating the same variant skips prediction. Useful when re-running on overlapping datasets.
2. **Post-annotation analysis** (`--from-cache`, `export parquet`): Enables columnar queries over previously annotated results --- filter by gene, consequence, or clinical significance without re-processing the input file.

In short: SQLite provides input data for the annotation engine, DuckDB captures its output for downstream use.

## Parallelism

Annotation is CPU-bound (consequence prediction per transcript). The annotator supports parallel execution via Go channels:

```mermaid
flowchart LR
    subgraph input
        Parser["Parser<br/><small>sequential read</small>"]
    end

    subgraph workers ["Worker Pool (GOMAXPROCS)"]
        W1["Worker 1<br/><small>Annotate batch</small>"]
        W2["Worker 2<br/><small>Annotate batch</small>"]
        W3["Worker 3<br/><small>Annotate batch</small>"]
        W4["Worker N<br/><small>Annotate batch</small>"]
    end

    subgraph output
        Reorder["Ordered Collector<br/><small>restore input order</small>"]
        Writer["Writer<br/><small>sequential write</small>"]
    end

    Parser --> W1
    Parser --> W2
    Parser --> W3
    Parser --> W4
    W1 --> Reorder
    W2 --> Reorder
    W3 --> Reorder
    W4 --> Reorder
    Reorder --> Writer
```

Variants are dispatched to workers via a buffered channel. Each worker annotates independently (transcript cache and SQLite index are read-only, thread-safe). An ordered collector reassembles results in input order for deterministic output.

**Performance** (4 cores, 1.05M TCGA variants):
- Sequential: ~11K variants/sec
- Parallel: ~16-23K variants/sec (1.5-2x speedup)

## Package Dependencies

```mermaid
graph TB
    CLI["cmd/vibe-vep<br/><small>CLI entry point</small>"]
    WASM["cmd/vibe-vep-wasm<br/><small>Browser build</small>"]

    CLI --> output
    CLI --> annotate
    CLI --> cache
    CLI --> duckdb
    CLI --> genomicindex
    CLI --> datasource
    WASM --> annotate
    WASM --> cache

    annotate["internal/annotate<br/><small>Consequence prediction</small>"]
    cache["internal/cache<br/><small>Transcript loading</small>"]
    output["internal/output<br/><small>MAF/VCF writing,<br/>compare, categorize</small>"]
    duckdb["internal/duckdb<br/><small>Variant cache</small>"]
    genomicindex["internal/genomicindex<br/><small>SQLite annotation index</small>"]
    parquet["internal/parquet<br/><small>Parquet export</small>"]
    mafpkg["internal/maf<br/><small>MAF parser</small>"]
    vcfpkg["internal/vcf<br/><small>VCF parser</small>"]
    datasource["internal/datasource/*<br/><small>OncoKB, Hotspots,<br/>ClinVar, AlphaMissense,<br/>SIGNAL, gnomAD, dbNSFP, dbSNP</small>"]

    output --> annotate
    output --> mafpkg
    output --> vcfpkg
    output --> cache
    duckdb --> cache
    CLI --> parquet
    CLI --> mafpkg
    CLI --> vcfpkg
    annotate --> cache
    genomicindex --> datasource

    style CLI fill:#e1f5fe
    style WASM fill:#e1f5fe
    style annotate fill:#fff3e0
    style cache fill:#fff3e0
    style output fill:#e8f5e9
    style duckdb fill:#fce4ec
    style genomicindex fill:#fce4ec
    style datasource fill:#f3e5f5
```

## CLI Commands

```
vibe-vep
  annotate          Annotate variants with consequence predictions
    maf             Annotate MAF file
    vcf             Annotate VCF file
    variant         Annotate a single variant
  compare           Compare variant annotations
    maf             Compare two MAF files (with optional --categorize)
    vcf             Compare two VCF files
  convert           Format conversion
    vcf2maf         Convert VCF to MAF format
  download          Download GENCODE + annotation source data
  prepare           Build genomic annotation index (SQLite)
  export            Export to analytical formats
    parquet         Export annotated results to Parquet
  config            Manage configuration
  version           Show version and loaded annotation sources
```
