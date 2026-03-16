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
