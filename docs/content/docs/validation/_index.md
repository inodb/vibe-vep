---
title: Validation
description: Validation methodology and results against cBioPortal datahub datasets.
weight: 3
---

## How Validation Works

The `vibe-vep compare` command compares existing MAF annotations against fresh predictions. For each MAF variant, the tool:

1. Parses the variant and its existing annotation (consequence, gene, transcript ID)
2. Re-annotates the variant against all overlapping GENCODE transcripts
3. Selects the best matching annotation for comparison (see below)
4. Normalizes consequence terms and compares

## Transcript Selection for Validation

When comparing against a MAF entry, the tool selects the best VEP annotation using this priority:

1. **Exact transcript ID match** — if the MAF specifies a transcript ID (e.g., ENST00000311936), use the annotation for that transcript. However, if the transcript's biotype has changed between GENCODE versions (e.g., was `protein_coding` but is now `retained_intron`) and the MAF has a coding consequence, the match is skipped.

2. **Same gene, best transcript** — among annotations for the same gene (by Hugo symbol), prefer:
   - Canonical transcripts (marked by GENCODE/Ensembl)
   - Protein-coding biotypes over non-coding
   - Higher-impact consequences (HIGH > MODERATE > LOW > MODIFIER)

3. **Any transcript, best match** — if no gene match, use the same ranking across all annotations.

## Consequence Normalization

MAF files may use different consequence terms than SO standard. The validation normalizes both sides before comparison:

- **MAF to SO mapping**: `Missense_Mutation` to `missense_variant`, `Silent` to `synonymous_variant`, etc.
- **Modifier stripping**: Drops terms like `non_coding_transcript_variant`, `NMD_transcript_variant`, `coding_sequence_variant`, `start_retained_variant`, `stop_retained_variant` that are secondary modifiers when a higher-impact term is present
- **Splice normalization**: Maps `splice_donor_region_variant` and `splice_donor_5th_base_variant` to `splice_region_variant`; drops `splice_region_variant` when a primary consequence is present
- **Impact-based stripping**: Drops `intron_variant` when splice donor/acceptor is present; drops UTR terms when a HIGH-impact term is present; drops `stop_gained`/`stop_lost` when co-occurring with `frameshift_variant`
- **Inframe grouping**: Normalizes `protein_altering_variant`, `inframe_deletion`, and `inframe_insertion` to a common term
- **Upstream/downstream tolerance**: MAF upstream/downstream calls are always accepted as matching, since different canonical transcript sets produce different transcript boundaries
- **Sorting**: Comma-separated terms are sorted alphabetically for consistent comparison

## Datahub Benchmark

Comprehensive validation against [cBioPortal/datahub](https://github.com/cBioPortal/datahub) MAF annotations across both assemblies:

| Tier | Assembly | Studies | Variants | Conseq Match | Time |
|------|----------|---------|----------|-------------|------|
| `datahub_gdc` | GRCh38 | 32 | 2,489,269 | 99.8% | ~2 min |
| `datahub_all` | GRCh37 | 431 | 30,612,130 | 99.7% | ~19 min |
| **Total** | | **463** | **~33M** | **99.7%** | **~21 min** |

### Corner Case Tests

66 unit tests covering all mismatch patterns discovered during benchmark analysis, runs in **7ms**:

```bash
go test ./internal/annotate/ -run 'TestDatahub_|TestEdge_|TestCorner_' -v
```

| File | Tests | Source |
|------|-------|--------|
| `corner_case_test.go` | 13 | Manual edge cases + GRCh37 bug fixes |
| `vep_edge_cases_test.go` | 17 | Ensembl VEP test patterns |
| `datahub_mismatch_test.go` | 36 | Datahub mismatch mining |

### Mismatch Collection

Collect all mismatches for analysis:

```bash
# GRCh38 mismatches
go test ./internal/output/ -run TestMismatchCollection -v -count=1 -timeout 60m
# Output: testdata/datahub_all/GRCh38/mismatches.json

# GRCh37 mismatches
go test ./internal/output/ -run TestMismatchCollectionGRCh37 -v -count=1 -timeout 120m
# Output: testdata/datahub_all/GRCh37/mismatches.json
```

## GRCh38 Results

32-study GDC benchmark from [cBioPortal/datahub](https://github.com/cBioPortal/datahub):

{{< validation-report assembly="grch38" >}}

## GRCh37 Results

431-study benchmark from [cBioPortal/datahub](https://github.com/cBioPortal/datahub):

{{< validation-report assembly="grch37" >}}

## Reproducing Validation

```bash
# Download test data
make download-datahub-gdc     # 32 GDC studies (~5 GB)
make download-datahub-all     # 431+ studies (~50 GB, includes GRCh37)

# Download GENCODE reference
vibe-vep download --assembly GRCh38
vibe-vep download --assembly GRCh37   # only needed for GRCh37 benchmark

# Run fast benchmark (GDC only, ~2 min)
go test ./internal/output/ -run TestValidationBenchmarkGDC -v -count=1 -timeout 10m

# Run full benchmark (all assemblies, ~21 min)
go test ./internal/output/ -run 'TestValidationBenchmark(GDC|All|AllGRCh37)' -v -count=1 -timeout 120m
```

Full markdown reports are generated in each study directory:
- `testdata/datahub_gdc/validation_report.md` — GRCh38 (32 GDC studies)
- `testdata/datahub_all/GRCh37/validation_report.md` — GRCh37 (431 studies)

## Which VEP release the benchmark is measured against

The genome-nexus parity work is measured against `vep111.genomenexus.org`,
which runs VEP cache 111. That release carries a bug in HGVS c. notation on
reverse-strand transcripts, so part of the measured difference is the oracle
being wrong rather than vibe-vep.

Ensembl [issue #1609](https://github.com/Ensembl/ensembl-vep/issues/1609)
reports it: for the same variant VEP 109 gave `c.1988_1989dup` while VEP 111
gave `c.1989_1990insAG`, which is a duplication of CT reported as an insertion
of AG. The maintainer states the bug was **introduced in release 110**, says a
fix would come in e112, and closed it as a duplicate of
[#1488](https://github.com/Ensembl/ensembl-vep/issues/1488). The example,
`NM_017617.5` (NOTCH1), is on the reverse strand, which the issue text does
not mention.

What goes wrong is that VEP 110 and 111 resolve the 3' ambiguity once in
genomic forward orientation and then project the result into transcript
coordinates. On a forward-strand transcript that is the same direction HGVS
requires; on a reverse-strand transcript it is the opposite one, so the
duplication is never recognised and no `hgvs_offset` is reported.

Measured over 1,378 reverse-strand indel consequences:

| | agreement with vibe-vep |
|---|---|
| VEP 111 (`vep111.genomenexus.org`) | 83.4% |
| VEP 116 (`www.genomenexus.org`) | 91.9% |

VEP 111 and VEP 116 agree with each other on **100%** of forward-strand indels
and **78.7%** of reverse-strand ones, so the release changed reverse-strand
behaviour only, and changed it toward vibe-vep. One case:

```
7:g.2968322_2968323insG, ENST00000355508 (reverse strand)
  VEP 111    c.76_77insC
  VEP 116    c.76dup
  vibe-vep   c.76dup
```

So a reverse-strand `dup`/`ins` difference against VEP 111 is the oracle's
defect, not ours, and is not something to fix. Forward-strand differences are
ours: there VEP 111 and 116 agree, and vibe-vep is the one not shifting far
enough.

The harness takes the oracle as a parameter, so both can be measured:

```bash
ORACLE_BASE=https://vep111.genomenexus.org ...   # VEP 111, matches our bundle
ORACLE_BASE=https://www.genomenexus.org    ...   # VEP 116, the later behaviour
```
