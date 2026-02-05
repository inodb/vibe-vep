# Validation Mismatch Report: TCGA PAAD (GRCh38)

## Summary

| Metric | Value |
|--------|-------|
| Dataset | TCGA-PAAD (Pancreatic Adenocarcinoma) |
| Assembly | GRCh38 |
| Annotation Source | GENCODE v46 |
| Total Variants | 24,849 |
| **Matches** | **15,472 (62.3%)** |
| **Mismatches** | **9,377 (37.7%)** |

**Total mismatches analyzed:** 9,377

## 1. Consequence Type Transitions (Top 30)

| MAF Consequence | VEP Consequence | Count | % |
|-----------------|-----------------|-------|---|
| missense_variant | synonymous_variant | 4624 | 49.3% |
| synonymous_variant | missense_variant | 1302 | 13.9% |
| missense_variant | stop_gained | 651 | 6.9% |
| stop_gained | missense_variant | 539 | 5.7% |
| stop_gained | synonymous_variant | 352 | 3.8% |
| missense_variant | stop_lost | 322 | 3.4% |
| missense_variant | missense_variant | 263 | 2.8% |
| splice_acceptor_variant | intron_variant | 238 | 2.5% |
| synonymous_variant | stop_gained | 210 | 2.2% |
| splice_donor_variant | intron_variant | 157 | 1.7% |
| missense_variant | non_coding_transcript_exon_variant | 105 | 1.1% |
| synonymous_variant | stop_lost | 98 | 1.0% |
| splice_region_variant | synonymous_variant | 97 | 1.0% |
| synonymous_variant | non_coding_transcript_exon_variant | 46 | 0.5% |
| downstream_gene_variant | missense_variant | 39 | 0.4% |
| splice_region_variant | missense_variant | 37 | 0.4% |
| downstream_gene_variant | non_coding_transcript_exon_variant | 26 | 0.3% |
| intron_variant | intron_variant | 26 | 0.3% |
| upstream_gene_variant | missense_variant | 22 | 0.2% |
| 3_prime_UTR_variant | non_coding_transcript_exon_variant | 20 | 0.2% |
| downstream_gene_variant | synonymous_variant | 19 | 0.2% |
| upstream_gene_variant | synonymous_variant | 19 | 0.2% |
| mature_miRNA_variant | non_coding_transcript_exon_variant | 15 | 0.2% |
| frameshift_variant | frameshift_variant | 13 | 0.1% |
| upstream_gene_variant | non_coding_transcript_exon_variant | 10 | 0.1% |
| downstream_gene_variant | intron_variant | 10 | 0.1% |
| stop_gained | stop_gained | 9 | 0.1% |
| splice_donor_variant | frameshift_variant | 6 | 0.1% |
| missense_variant | intron_variant | 6 | 0.1% |
| upstream_gene_variant | intron_variant | 6 | 0.1% |

## 2. Genes with Most Mismatches (Top 20)

| Gene | Mismatches |
|------|------------|
| TP53 | 41 |
| TTN | 35 |
| SMAD4 | 19 |
| MUC16 | 18 |
| CDKN2A | 14 |
| LRP1B | 13 |
| SYNE1 | 12 |
| NEB | 11 |
| RELN | 9 |
| DIDO1 | 9 |
| PCDH15 | 9 |
| SACS | 9 |
| KMT2C | 9 |
| DST | 8 |
| ADAMTS12 | 8 |
| MYO9A | 8 |
| ASPM | 8 |
| MYO15A | 8 |
| SYNE2 | 7 |
| AKAP9 | 7 |

## 3. Mismatch Categories

| Category | Count | % |
|----------|-------|---|
| MAF=missense, VEP=synonymous | 4624 | 49.3% |
| MAF=synonymous, VEP=missense | 1339 | 14.3% |
| MAF=missense, VEP=stop | 973 | 10.4% |
| Other | 599 | 6.4% |
| MAF=stop_gained, VEP=missense | 539 | 5.7% |
| MAF=splice, VEP=intron | 400 | 4.3% |
| MAF=splice+coding, VEP=coding_only | 361 | 3.8% |
| MAF=stop_gained, VEP=synonymous | 352 | 3.8% |
| UTR/flanking region differences | 170 | 1.8% |
| Indel annotation differences | 20 | 0.2% |

## 4. Sample Mismatches by Category

### MAF=missense, VEP=synonymous (4624 cases)

| Variant | Gene | MAF_AA | VEP_AA |
|---------|------|--------|--------|
| 4:103589873 C>T | TACR3 | p.V403M |  |
| 15:42450147 G>A | ZNF106 | p.H686Y |  |
| 18:10472023 C>T | APCDD1 | p.R246W |  |
| X:139619651 G>C | MCF2 | p.T248R |  |
| 1:77565646 G>C | ZZZ3 | p.N902K |  |
| 2:42492838 C>G | KCNG3 | p.G222R |  |
| 3:142375918 A>C | XRN1 | p.V953G |  |
| 4:72339620 C>G | ADAMTS3 | p.M245I |  |
| 4:127895084 G>C | PLK4 | p.W898C |  |
| 5:124649066 A>T | ZNF608 | p.S440T |  |

### MAF=synonymous, VEP=missense (1339 cases)

| Variant | Gene | MAF_AA | VEP_AA |
|---------|------|--------|--------|
| 5:140877668 C>T | PCDHA12 | p.C732= | K732N |
| 14:70051102 C>T | SLC8A3 | p.T679= | D679E |
| 15:73236429 G>A | NEO1 | p.T458= | F458L |
| 4:41603899 G>A | LIMCH1 | p.L123= | S123R |
| 4:72548727 C>G | ADAMTS3 | p.T85= | K85N |
| 19:38100027 G>T | SIPA1L3 | p.G577= | Q577H |
| 22:21962499 C>T | TOP3B | p.E485= | F485L |
| 5:33984416 C>G | SLC45A2 | p.V56= | M56I |
| 7:21617648 C>A | DNAH11 | p.R1375= | F1375L |
| 7:75423073 C>A | POM121C | p.P393= | E393D |

### MAF=stop_gained, VEP=synonymous (352 cases)

| Variant | Gene | MAF_AA | VEP_AA |
|---------|------|--------|--------|
| 16:53231661 C>G | CHD9 | p.Y796* |  |
| 5:58454599 A>C | PLK2 | p.L681* |  |
| 3:74302716 G>T | CNTN3 | p.S587* |  |
| 8:14090449 C>T | SGCZ | p.W311* |  |
| 17:7676257 G>A | TP53 | p.Q38* |  |
| 3:391790 T>A | CHL1 | p.Y953* |  |
| 15:89606822 C>T | TICRR | p.Q907* |  |
| 6:7229036 C>T | RREB1 | p.Q313* |  |
| 7:123868773 C>G | HYAL4 | p.S167* |  |
| 18:58607414 C>T | ALPK2 | p.W45* |  |

### MAF=splice, VEP=intron (400 cases)

| Variant | Gene | MAF_AA | VEP_AA |
|---------|------|--------|--------|
| 5:114364003 T>C | KCNN2 | p.X406_splice |  |
| 1:22659644 G>A | C1QB | p.X61_splice |  |
| 8:28859564 C>A | INTS9 | p.X4_splice |  |
| 12:120660842 T>C | CABP1 | p.X313_splice |  |
| 21:44059102 G>A | TRAPPC10 | p.X227_splice |  |
| 16:89513039 T>C | SPG7 | p.X126_splice |  |
| 9:69336270 G>T | ENTREP1 | p.X157_splice |  |
| 14:61449131 ATAAATGTATAATTGCTGGATTTAATTTCCTAGTGTGCACCTGTGTC> | PRKCH | p.X205_splice |  |
| 14:70733343 C>G | MAP3K9 | p.X676_splice |  |
| 15:59067084 G>C | RNF111 | p.X562_splice |  |

### MAF=missense, VEP=stop (973 cases)

| Variant | Gene | MAF_AA | VEP_AA |
|---------|------|--------|--------|
| 3:39134218 G>A | TTC21A | p.V925M | *925R |
| 9:138121530 G>A | CACNA1B | p.R2184H | L2184* |
| 11:46873167 G>A | LRP4 | p.R1506W | R1506* |
| 8:28717972 C>A | EXTL3 | p.P638H | W638* |
| 15:40389884 C>T | KNSTRN | p.R214W | Q214* |
| 3:74302717 A>T | CNTN3 | p.S587T | *587K |
| 15:60455372 A>G | ICE2 | p.I246T | *246S |
| 16:23215206 C>T | SCNN1G | p.R563C | R563* |
| 17:7673776 G>A | TP53 | p.R282W | Q282* |
| 1:109017528 C>A | WDR47 | p.D78Y | E78* |

## 5. Likely Causes of Mismatches

| Cause | Evidence | Recommendation |
|-------|----------|----------------|
| **Different transcript versions** | missense<->synonymous swaps with different AA | Use same GENCODE version as original annotations |
| **Different canonical transcript** | Same position, completely different consequence | Filter to MANE Select transcripts |
| **Strand handling** | Stop gained vs synonymous patterns | Verify strand in GTF matches VEP |
| **Splice site definitions** | splice_donor -> intron_variant | VEP uses 2bp, may differ from original |
| **CDS boundary differences** | downstream_gene -> coding | Check transcript start/end coordinates |
