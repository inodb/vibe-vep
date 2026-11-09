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

So a reverse-strand `dup`/`ins` difference against VEP 111 is the oracle's
defect, not ours, and is not something to fix. Forward-strand differences are
ours: vibe-vep is the one not shifting far enough there, and that gap shows
against any oracle.

`7:g.2968322_2968323insG` on ENST00000355508, a reverse-strand transcript, is
the shape of it: VEP 111 writes `c.76_77insC` where vibe-vep writes `c.76dup`.
The inserted base duplicates its neighbour, so the duplication is the correct
name.

### A second VEP 111 defect: 5'-incomplete transcripts numbered two ways

Where a transcript's CDS is 5'-incomplete (GENCODE `cds_start_NF`, 11,145 of
196,520 transcripts), VEP numbers c.1 from the first base of the notional
complete first codon, so the annotated CDS starts at c.2 or c.3. vibe-vep does
the same -- `CDSStartOffset`, `(3 - phase) % 3` -- and that is what keeps its
own `hgvsc` and `protein_start` describing the same residue.

VEP applies the offset to an exonic substitution and to nothing else. Over
20,000 variants, 977 `hgvsc` differences across 168 transcripts are exactly
this, and the breakdown is clean:

```
496  substitution, intronic anchor
319  deletion, exonic
 96  deletion, intronic anchor
 54  duplication, exonic
  8  duplication, intronic anchor
  4  insertion, exonic
  0  substitution, exonic
```

Not one exonic substitution disagrees, and nothing else agrees. The gap is
always the offset: 2 where the first CDS phase is 1, 1 where it is 2.

RYBP ENST00000477973 (phase 1, offset 2) shows both answers on one
transcript:

```
exonic substitution   vep111  c.632C>A        vibe-vep  c.632C>A       agree
exonic deletion       vep111  c.506_508del    vibe-vep  c.508_510del
```

The oracle calls that base c.632 -- raw position 630 plus the offset -- and
then names a deletion three exons away from the raw count alone. ANKRD11
ENST00000562194 is the intronic form of the same thing: exon 2 starts at
c.154 on both sides, and the oracle anchors intron 1 to c.152.

So the padding is applied where the CDS coordinate is reached through the
codon, and skipped where it is reached by counting cDNA bases. vibe-vep is
consistent with itself and with the oracle's own exonic substitutions, and
matching the rest would mean reproducing the inconsistency: its `c.` names
would then contradict its own `protein_start` on these transcripts. These
differences are therefore not a gap to close.

An earlier version of this section called the class "intronic anchors" and
put it at roughly 900 differences. That was the visible half. The exonic
indels are the same defect and were being counted as a shift gap of ours,
which they are not.

### Genomic HGVS is renormalised; transcript-relative HGVS is not, yet

`hgvsg` is the one field that is not transcript-relative, and it is named
from the reference rather than from the query. HGVS puts a change at the
3'-most position the reference allows and calls a repeat of the preceding
bases a duplication, so `13:41134330 G>-` in a GGG run is `13:g.41134332del`
and a G inserted after a G is `dup`, not `insG`.

That needs reference bases at serve time, which the transcript cache cannot
supply: the shift runs along the chromosome and the variant need not be in
any transcript. `vibe-vep prepare` therefore builds
`genome_reference.bin`, the assembled chromosomes two bits to a base, about
740 MB, read a window at a time. `serve` builds it on first start if
`prepare` has not. Without it `hgvsg` falls back to naming the change where
the query put it, which is what every release before this one did.

Over 20,000 variants this closed 21,581 of 29,067 differing fields.

Two divergences remain, and in both of them VEP is the one not following
HGVS. A 57-base insertion that exactly repeats 18:48591886-48591942 is a
tandem duplication; VEP names it `48591936_48591937ins...`. This is the same
shape as the reverse-strand `dup`/`ins` defect above and is not something to
fix.

The transcript-relative names needed the same thing. `hgvsc` was shifted only
as far as the transcript cache held sequence -- the coding sequence, both UTR
tails and a hundred bases of each intron flank -- so a position deeper in an
intron than that stopped, and a non-coding transcript, which has none of the
three, never shifted at all. The shifts now fall back to the reference, and
the forward-strand non-coding gap, 964 transcript consequences, is gone.

Over 20,000 variants the whole series took differing fields from 29,067
(1.1371%) to 4,498 (0.1760%), and variants identical in every field from
16,503 to 18,319. On the full cohort it took the MAF from 98.4965% to
98.5138% against VEP 111, with HGVS_Offset going from 121 differing rows
to 2.

### What is left is the two oracle defects, and almost nothing else

Of the 2,512 `hgvsc` differences remaining over 20,000 variants:

```
1,530   reverse-strand dup/ins      Ensembl #1609
  977   cds_start_NF numbering      oracle applies its own offset to
                                    exonic substitutions and nothing else
    5   unparsed
```

Every one of the 977 differs from ours by exactly the 5'-incomplete offset,
1 or 2, and never in the other direction.

The 224 that looked like a shift gap of ours are the same dup/ins pattern
seen from the wrong angle: `c.269_270insCAG` against `c.267_269dup` is one
event named two ways, and the coordinates differ only because a duplication
names the duplicated bases while an insertion names the positions it falls
between. There is no HGVSc shift gap of ours left in this corpus.

### What the amino acids, the splice terms and the part-codons closed

Five gaps of ours were closed after the reference landed, each measured on
its own against the oracle and then against the full cohort.

**`amino_acids` came from the named protein change, not the codons.** That
change knows only the residue it renames, so a two-base substitution
straddling a codon boundary read `R/C` where VEP reads `VR/VC` -- both
codons, even though only the second residue moves. The codons field already
covered the right span and already matched the oracle, so it is the
authority now. One rule about part-codons had to be measured rather than
assumed: a frameshift leaves a part-codon named X, except where the change
translates to a stop and nothing else, where VEP writes `Y/*` and never
`Y/*X`. It really is only that case -- `*L` with a part-codon after it is
still `*LX`. Checked against 4,334 part-codon sides pulled from the oracle,
the rule and the oracle disagree nowhere. This took `amino_acids` from 382
differences to 0, and `aminoAcids`, `aminoAcidRef` and `aminoAcidAlt` in the
annotation summary from 285 to 0.

**A non-coding exon's splice site was never named.** A deletion running out
of a non-coding exon and over the donor dinucleotide destroys that site; the
non-coding path asked only whether an exon edge was touched and said
`splice_region_variant`. `spanSpliceTerms` already enumerates the right
terms for the coding path and contributes nothing for a non-coding
transcript's exonic positions, so joining it with the exon term is the whole
fix. An insertion is excluded: it sits between two bases and takes neither
out.

**A CDS that ends in a part-codon had no name for it.** VEP reports
`incomplete_terminal_codon_variant` there and still gives the codons it has,
`X` for the residue and the position as a terminator: MSH6 ENST00000455383
c.346G>A is `Gt/At`, `X`, `p.Ter116=`. Reporting the codons is what restores
the protein position as well, since a transcript consequence with no codons
withholds one.

**Coming the other way was the same gap.** A deletion starting inside an
intron and running over the acceptor into a non-coding exon was reported as
the bare site plus `non_coding_transcript_variant`, because `spanSpliceTerms`
names an exonic base from its CDS position or its UTR and a non-coding
transcript has neither. Saying so there also drops
`splice_polypyrimidine_tract_variant`, which was already encoded -- the tract
is suppressed whenever an exonic term is present, and the exonic term was the
one missing. A last pairing went with it: `mature_miRNA_variant` already says
the transcript codes for nothing, so VEP does not add the generic modifier
beside it, and neither does vibe-vep now.

**A 3'UTR insertion did not report how far it shifted.** The walk through
the cached 3'UTR found the right name and threw the distance away, so
`hgvs_offset` was absent where VEP reports it -- TCF7L2
10:114925316_114925317insA is `c.*30dup` with an offset of 9. 109 of the
110 missing offsets were this.

Together these took differing fields from 4,498 (0.1760%) to 3,474
(0.1359%), variants identical in every field from 18,319 to 18,627, and the
full MAF from 98.5138% to 98.5203% against VEP 111. Every cell that moved
in the MAF moved toward VEP: 35 better, none worse.

### Measuring the whole cohort, not a sample

The numbers above are over the first 20,000 variants of the list. Running
all **274,506 distinct variants** of the cohort takes 18 minutes and compares
**34,519,155 fields**, and it says two useful things.

The sample was honest: 0.1359% differing over 20,000 against **0.1449%** over
everything, with 92.8% of variants identical in every field and still **zero
missing transcripts** out of 1,973,100. And `hgvsc` keeps its shape at scale
-- of 36,358 differences, 21,560 are the reverse-strand dup/ins of Ensembl
\#1609 and 14,459 are `cds_start_NF` numbering, so **99.1% of them are the two
declared oracle defects**.

What the sample could not do is rank what is left. At 20,000 variants
`consequence_terms` is 71 differences across 27 unrelated shapes; at 274,506
it is 1,257, and more than a third of them turn out to be one thing.

### Where the annotation stops

A transcript GENCODE tags `cds_end_NF` ends at its last coding base: no stop
codon, no 3'UTR, nothing annotated past it. One whose start codon sits on its
first base has the same hole at the other end. A deletion running off such an
edge was being described by what the annotation lacks rather than by what VEP
reports there, and that single hole produced five different wrong answers:

| what we said | what VEP says |
|---|---|
| nothing, for the bases past the edge | the UTR that lies at that end |
| `splice_donor_variant`, two bases after the final exon | a splice site needs an intron to sit in |
| `stop_lost`, reaching the last three coding bases | there is no terminator in the annotation to lose |
| a protein consequence for the span | the regions it covers, the way it declines `start_lost` when the 5' end is missing |
| a symmetric protein range | lopsided: `10,10` off the 3' end, `?-11` off the 5' |

Two more came out of chasing those, both in the response rather than the
prediction. `GNTranscriptConsequence.MarshalJSON` re-derives the protein range
from the HGVSc and was reading "neither end is coding" out of an *empty*
string, discarding a position the annotator had got right. And
`GNIntegerRange` wrote an unknown start as `0` where genome-nexus writes
`null` -- which is what the annotation pipeline renders as the `?-` in
`?-15/136`.

That last one is the reason both measurements are run. The response diff saw
nothing wrong: our `transcript_consequences` entry already matched the oracle
field for field. The MAF goes through the annotation pipeline, which reads
`annotation_summary` instead, where it did not match. Either measurement alone
would have shipped it.

Together: **50,032 differing fields to 48,808** (0.1449% to 0.1414%), with
`consequence_terms` 1,257 to 560, `protein_end` -227, `protein_start` -159 and
`most_severe_consequence` 103 to 27. Variants identical in every field 254,806
to 254,975. On the full MAF, 98.5203% to **98.5215%**, `Consequence` 55 to 47
and `Variant_Classification` 9 to 8 -- 13 cells better, 1 neither, none worse.

### What is left in consequence_terms

560 differences over the full cohort, of which 220 are too long to classify
from the dump. Of the 340 that can be read, the largest group is 48 miRNA
transcripts where the variant sits in the transcript but outside the mature
product: VEP reports `non_coding_transcript_exon_variant` and vibe-vep reports
`mature_miRNA_variant`, and the mature coordinates come from Ensembl's miRNA
annotation, which the GENCODE GTF does not carry. After that the counts fall
to 43, 40, 35 and below, with no shared cause -- the next one costs its own
investigation.

### hgnc_id, and why the gene id is the wrong key

1,695 cells differ. GENCODE v19 carries no `hgnc_id` attribute at all, so on
GRCh37 the identifier is looked up in the biomart HGNC export **by gene
symbol**, while VEP reads it off the gene. Re-keying on the Ensembl gene id
looks obviously right and was measured **10x worse**: 1,695 to 17,508 cells,
whole-API 0.1449% to 0.1907%.

The reason is vintage, but the oracle's. The HGNC export is *newer* than VEP
111's gene-to-HGNC xrefs, so the gene id finds recent lncRNA assignments --
HGNC 55495, 56033 -- on 127 clone-named transcripts (`RP11-770J1.3`,
`AC002553.4`) that the oracle has nothing for. The symbol lookup fails on
exactly those genes, because a clone name is not an HGNC symbol, so **its
blind spot matches the oracle's**. Matching VEP here means not using the
better data.

The residue is four ways a symbol is not a key: one HGNC renamed (`MARCH9` to
`MARCHF9`, whose `previous_symbols` cell is Excel-mangled to `9-Mar`), one two
GENCODE loci share (`SNORA4`), one HGNC retired into another gene (`CASC7`
into `AGO2`, which would lend a protein-coding gene's identifier to a
lincRNA), and one real source disagreement (`CCND2-AS1` vs `CCND2-AS2`).

### Comparing against a later release

**`www.genomenexus.org` is not a clean VEP 116 measurement.** It answers from
genome-nexus's own annotation store, which was populated by an earlier VEP, so
what comes back is a mixture of releases rather than the behaviour of one.
Numbers taken from it were used here to ask whether the reverse-strand
difference was version-specific, and they pointed that way, but they do not
establish what any particular VEP release does.

The claim that the behaviour changed rests on Ensembl issue #1609 instead,
where the maintainer states the bug was introduced in release 110 and says a
fix was coming in e112. Settling it by measurement means running
`ensembl-vep` at a pinned release directly, which has not been done.

The harness takes the oracle as a parameter:

```bash
ORACLE_BASE=https://vep111.genomenexus.org ...   # VEP 111, matches our bundle
```
