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

## Prerequisites

| | |
|---|---|
| **Go** | 1.24+, CGO enabled (SQLite and DuckDB) |
| **Disk** | ~1.3 GB for the transcript cache alone; ~13 GB with every annotation source built. See the table below. |
| **Memory** | ~1 GB resident for GRCh37 with the transcript cache loaded. Annotation indexes are memory-mapped, so extra sources cost page cache rather than heap. |
| **Network** | Only `vibe-vep download`. The server itself makes outbound calls only for the `my_variant_info` and `oncokb` fields. |

```bash
vibe-vep download --assembly GRCh37    # fetch reference data
vibe-vep prepare  --assembly GRCh37    # build the transcript cache and indexes
vibe-vep serve    --assembly GRCh37 --port 38080
```

## Reference data

`vibe-vep download` fetches everything below into
`~/.vibe-vep/<assembly>/raw/` (override the root with `$VIBE_VEP_DATA_DIR`).
`vibe-vep prepare` turns it into the cache and indexes the server reads.

Three of these decide whether the server's answers match a Genome Nexus
deployment's rather than merely resembling them, and the server logs at
startup which of them it found:

| File | What it decides |
|---|---|
| `ensembl_biomart_canonical_transcripts_per_hgnc.txt` | Which transcript a variant is reported on, and the gene's approved symbol and Entrez id |
| `cancerGeneList.tsv` | Which gene a variant falling in two of them is reported against |
| `VUEs.json` | Where published evidence overrules the predicted protein effect |

Without the canonical transcript table the server falls back to Ensembl's own
canonical flags, which for a cancer panel is usually the wrong transcript.

### Where each file comes from

Sizes are for GRCh37 and were measured, not estimated.

| id | file | origin | size |
|---|---|---|---|
| gencode | `gencode.v19.annotation.gtf.gz`, `.pc_transcripts.fa.gz` | [GENCODE](https://www.gencodegenes.org/) via EBI FTP | 38 + 28 MB |
| gencode | `GRCh37.p13.genome.fa.gz` | GENCODE via EBI FTP | 805 MB |
| hgnc | `ensembl_biomart_canonical_transcripts_per_hgnc.txt` | [genome-nexus-importer](https://github.com/genome-nexus/genome-nexus-importer) `data/<bundle>/export/` | 23 MB |
| ensembl | `ensembl_biomart_transcripts.json.gz` | genome-nexus-importer `data/<bundle>/export/` | 27 MB |
| — | `annotation_version.txt` | genome-nexus-importer `data/<bundle>/export/` | 3 KB |
| pfam | `pfamA.txt`, `ensembl_biomart_pfam.txt` | genome-nexus-importer | 0.8 + 13 MB |
| ptm | `ptm.json.gz` | genome-nexus-importer `data/ptm/export/` | 18 MB |
| cancer_hotspots | `hotspots_v2_and_3d.txt` | genome-nexus-importer `data/<bundle>/export/` | 1.4 MB |
| signal | `mutations.json.gz` | genome-nexus-importer `data/signal/export/` | 33.6 MB |
| oncokb | `cancerGeneList.tsv` | [OncoKB](https://www.oncokb.org/cancerGenes); a copy ships in this repo | 0.2 MB |
| revue | `VUEs.json` | [reVUE-data](https://github.com/knowledgesystems/reVUE-data) | 0.4 MB |
| uniprot | `enst_to_uniprot_mapping_id.txt` | genome-nexus-importer | 4 MB |
| refseq | `ensembl_refseq_mrna.GRCh37.tsv.gz` | Ensembl | 0.2 MB |
| clinvar | `clinvar.vcf.gz` | [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) NCBI FTP | 192 MB |
| alphamissense | `AlphaMissense_hg19.tsv.gz` | [AlphaMissense](https://github.com/google-deepmind/alphamissense), CC BY 4.0 | 594 MB |
| sift, polyphen | `translation_md5.txt.gz`, `protein_function_predictions.txt.gz` | Ensembl variation FTP, release matched to the GENCODE release | 4.1 GB |

**Keep the downloaded files.** `vibe-vep prepare` builds its indexes from them,
and it rebuilds whenever an index looks out of date. Deleting a raw file after
its index exists means the next rebuild has nothing to read: it will now fail
loudly rather than leave an empty index behind, but you will still be waiting
on a re-download. Roughly 2.6 GB of raw files buys you a rebuild that takes
minutes instead of an hour.

### Why the Ensembl transcript records are a file and not the GENCODE cache

The `/ensembl` endpoints report a transcript's CCDS and RefSeq accessions, its
exon identifiers and versions, and its Pfam domain ranges. The GENCODE GTF
carries none of those, so they cannot be derived from the transcript cache
however much is read from it.

genome-nexus answers those endpoints from its `ensembl.transcript` collection,
which genome-nexus-importer publishes as
`ensembl_biomart_transcripts.json.gz`. vibe-vep reads the same file into
`ensembl_transcripts.sqlite` (215,170 transcripts, 336 MB, about 15 seconds to
build), so the two agree by construction rather than by effort.

Without it those endpoints still answer, from the GENCODE cache, with the
fields it can supply and no others.

`VUEs.json` is the one that moves: reVUE publishes new findings, and a Genome
Nexus deployment re-reads the catalogue every time it starts. This server
reads the downloaded copy instead, which keeps a run reproducible but means a
stale copy quietly reports the predicted effect for a variant the catalogue
has since revised. It logs the number of effects it loaded at startup, and
`vibe-vep download` refreshes it.

See [Annotation Sources](../annotation-sources/) for what each source
contributes to an annotation and how it is matched.
