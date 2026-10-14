---
title: Annotation Server
weight: 6
aliases:
  - /docs/annotation-server/
---

# Annotation Server

`vibe-vep serve` answers the Genome Nexus API directly, so anything built
against Genome Nexus — the [annotation
pipeline](https://github.com/genome-nexus/genome-nexus-annotation-pipeline),
cBioPortal's mutation mapper, your own scripts — can be pointed at it by base
URL alone. There is no Genome Nexus service and no MongoDB behind it: one Go
process holds the transcript cache and the reference tables in memory.

```bash
vibe-vep serve --assembly GRCh37 --port 38080
```

A server with a single assembly answers at the root, the way a Genome Nexus
deployment does:

```bash
curl "http://localhost:38080/annotation/genomic/7,140453136,140453136,A,T?isoformOverrideSource=mskcc&fields=annotation_summary"
curl "http://localhost:38080/version"
```

Load more than one assembly and each is addressed by name instead, since there
is no longer one obvious answer:

```bash
vibe-vep serve --assembly GRCh37,GRCh38 --port 38080
curl "http://localhost:38080/genome-nexus/grch38/annotation/genomic/7,140753336,140753336,A,T"
```

## Annotating a MAF

Point the annotation pipeline at the server and run it unchanged:

```bash
GENOMENEXUS_BASE=http://localhost:38080 \
  java -jar annotationPipeline.jar \
    --filename input.maf --output-filename annotated.maf \
    --isoform-override mskcc --mode extended
```

## Reference data

`vibe-vep download` fetches everything the server needs. Three of those files
are what make its answers match a Genome Nexus deployment's rather than merely
resemble them, and the server logs at startup which of them it found:

| File | What it decides |
|---|---|
| `ensembl_biomart_canonical_transcripts_per_hgnc.txt` | Which transcript a variant is reported on, and the gene's approved symbol and Entrez id |
| `cancerGeneList.tsv` | Which gene a variant falling in two of them is reported against |
| `VUEs.json` | Where published evidence overrules the predicted protein effect |

Without the canonical transcript table the server falls back to Ensembl's own
canonical flags, which for a cancer panel is usually the wrong transcript.

`VUEs.json` is the one that moves: reVUE publishes new findings, and a Genome
Nexus deployment re-reads the catalogue every time it starts. This server
reads the downloaded copy instead, which keeps a run reproducible but means a
stale copy quietly reports the predicted effect for a variant the catalogue
has since revised. It logs the number of effects it loaded at startup, and
`vibe-vep download` refreshes it.
