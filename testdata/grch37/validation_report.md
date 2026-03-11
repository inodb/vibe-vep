# GRCh37 Validation Report

Generated: 2026-03-11 16:10 UTC  
Assembly: GRCh37  
GENCODE transcripts: 196520 (loaded from gob cache in 1.75s)  
Workers: 4 (GOMAXPROCS)

## Match Rates

| Study | Variants | Conseq Match | Conseq Mismatch | Conseq Rate | HGVSp Match | HGVSp Mismatch | HGVSp Rate | HGVSc Match | HGVSc Mismatch | HGVSc Rate |
|-------|----------|-------------|-----------------|-------------|-------------|----------------|------------|-------------|----------------|------------|
| msk_impact_50k_2026 | 479096 | 476530 | 100 | 99.5% | 442364 | 105 | 92.3% | 459316 | 161 | 95.9% |
| **Total** | **479096** | **476530** | **100** | **99.5%** | **442364** | **105** | **92.3%** | **459316** | **161** | **95.9%** |

## Consequence Category Breakdown

| Study | delins_normalized | match | mismatch | position_shift | transcript_model_change | upstream_reclassified |
|-------|------|------|------|------|------|------|
| msk_impact_50k_2026 | 2343 | 469297 | 100 | 92 | 31 | 7233 |
| **Total** | **2343** | **469297** | **100** | **92** | **31** | **7233** |

## HGVSp Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | fuzzy_fs | maf_empty | maf_nonstandard | match | mismatch | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| msk_impact_50k_2026 | 6929 | 5191 | 41 | 15608 | 364 | 496 | 442364 | 105 | 5319 | 1496 | 726 | 30 | 26 | 401 |
| **Total** | **6929** | **5191** | **41** | **15608** | **364** | **496** | **442364** | **105** | **5319** | **1496** | **726** | **30** | **26** | **401** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | mismatch | position_shift | transcript_model_change |
|-------|------|------|------|------|------|------|------|------|
| msk_impact_50k_2026 | 6861 | 10125 | 407 | 414 | 459316 | 161 | 1811 | 1 |
| **Total** | **6861** | **10125** | **407** | **414** | **459316** | **161** | **1811** | **1** |

## Cancer Gene Mismatches

359/475 cancer genes have 100% match across all columns. Mismatches in 116 gene(s):

| Gene | Variants | Conseq Mismatches | HGVSp Mismatches | HGVSc Mismatches |
|------|----------|-------------------|------------------|------------------|
| TP53 | 26720 | 9 | 7 | 9 |
| KRAS | 9096 | 0 | 4 | 11 |
| EGFR | 4039 | 1 | 5 | 5 |
| DNMT3A | 1069 | 2 | 2 | 6 |
| KDM6A | 2144 | 3 | 3 | 3 |
| ZFHX3 | 3661 | 2 | 2 | 4 |
| ARID2 | 2282 | 2 | 2 | 3 |
| CDKN2A | 2805 | 3 | 2 | 2 |
| CHEK2 | 569 | 1 | 3 | 3 |
| PIK3CA | 8894 | 1 | 3 | 3 |
| INPPL1 | 1451 | 2 | 2 | 2 |
| ARID1B | 2276 | 2 | 1 | 2 |
| ATRX | 3286 | 1 | 1 | 3 |
| KMT2D | 6791 | 1 | 2 | 2 |
| NOTCH3 | 2573 | 1 | 1 | 3 |
| PDGFRA | 1521 | 1 | 2 | 2 |
| SPEN | 2688 | 1 | 1 | 3 |
| AR | 1486 | 0 | 2 | 2 |
| EP300 | 2280 | 1 | 1 | 2 |
| GATA3 | 1784 | 4 | 0 | 0 |
| KEAP1 | 1847 | 0 | 2 | 2 |
| MITF | 427 | 2 | 1 | 1 |
| PLK2 | 651 | 0 | 2 | 2 |
| WT1 | 705 | 1 | 1 | 2 |
| APC | 8976 | 1 | 1 | 1 |
| ATM | 3720 | 0 | 1 | 2 |
| AXIN2 | 1063 | 1 | 1 | 1 |
| DNAJB1 | 273 | 1 | 1 | 1 |
| ERCC2 | 846 | 1 | 1 | 1 |
| ESR1 | 1371 | 1 | 1 | 1 |
| FGFR2 | 1128 | 0 | 1 | 2 |
| FLCN | 583 | 1 | 1 | 1 |
| HLA-B | 589 | 2 | 0 | 1 |
| IKZF1 | 1094 | 1 | 1 | 1 |
| KMT2C | 4612 | 1 | 1 | 1 |
| MST1R | 914 | 1 | 1 | 1 |
| NF1 | 4189 | 1 | 1 | 1 |
| NKX3-1 | 216 | 1 | 1 | 1 |
| NUP93 | 626 | 1 | 1 | 1 |
| PIK3C3 | 635 | 1 | 1 | 1 |
| PPARG | 449 | 0 | 1 | 2 |
| PTPRS | 2136 | 0 | 1 | 2 |
| RECQL4 | 970 | 1 | 1 | 1 |
| RUNX1 | 940 | 1 | 1 | 1 |
| SLFN11 | 152 | 1 | 1 | 1 |
| SMAD3 | 652 | 3 | 0 | 0 |
| AKT3 | 548 | 0 | 1 | 1 |
| ATR | 1940 | 0 | 1 | 1 |
| B2M | 778 | 2 | 0 | 0 |
| BRCA1 | 1402 | 0 | 1 | 1 |
| CCND1 | 385 | 0 | 1 | 1 |
| CD79B | 230 | 0 | 1 | 1 |
| CDK12 | 1417 | 0 | 1 | 1 |
| EPHB1 | 1562 | 0 | 1 | 1 |
| ERBB2 | 2092 | 0 | 1 | 1 |
| EZH2 | 601 | 0 | 1 | 1 |
| FLT1 | 1530 | 0 | 1 | 1 |
| FOXA1 | 1334 | 0 | 1 | 1 |
| GRIN2A | 2745 | 0 | 1 | 1 |
| HGF | 1317 | 0 | 1 | 1 |
| JAK1 | 1526 | 0 | 1 | 1 |
| KDR | 1740 | 0 | 1 | 1 |
| MED12 | 2074 | 0 | 1 | 1 |
| MSH3 | 972 | 0 | 1 | 1 |
| MTOR | 1958 | 0 | 1 | 1 |
| NBN | 646 | 0 | 1 | 1 |
| NOTCH4 | 2020 | 1 | 0 | 1 |
| NRAS | 1238 | 0 | 1 | 1 |
| PDCD1 | 505 | 0 | 1 | 1 |
| PIK3C2G | 1816 | 0 | 1 | 1 |
| REST | 121 | 0 | 1 | 1 |
| SDHA | 503 | 0 | 1 | 1 |
| SMARCA4 | 2858 | 0 | 1 | 1 |
| STAG2 | 1576 | 1 | 0 | 1 |
| SYK | 569 | 0 | 1 | 1 |
| TBX3 | 1621 | 0 | 1 | 1 |
| TGFBR1 | 716 | 0 | 1 | 1 |
| TSHR | 708 | 0 | 1 | 1 |
| VEGFA | 180 | 0 | 1 | 1 |
| AURKB | 222 | 1 | 0 | 0 |
| BBC3 | 189 | 1 | 0 | 0 |
| BRCA2 | 2722 | 0 | 0 | 1 |
| CARD11 | 1718 | 0 | 0 | 1 |
| CARM1 | 362 | 1 | 0 | 0 |
| CREBBP | 2811 | 0 | 0 | 1 |
| CTLA4 | 261 | 0 | 0 | 1 |
| ELF3 | 939 | 1 | 0 | 0 |
| ERG | 628 | 0 | 0 | 1 |
| FLT3 | 1091 | 0 | 0 | 1 |
| FOXO1 | 497 | 0 | 0 | 1 |
| FUBP1 | 617 | 1 | 0 | 0 |
| IRS1 | 1132 | 0 | 0 | 1 |
| KNSTRN | 262 | 1 | 0 | 0 |
| LATS1 | 1115 | 0 | 0 | 1 |
| MAP2K2 | 369 | 1 | 0 | 0 |
| MSH2 | 968 | 1 | 0 | 0 |
| MSH6 | 1319 | 0 | 0 | 1 |
| MST1 | 274 | 0 | 0 | 1 |
| MYC | 385 | 1 | 0 | 0 |
| NTHL1 | 164 | 1 | 0 | 0 |
| NTRK1 | 1006 | 0 | 0 | 1 |
| PIK3CG | 1778 | 0 | 0 | 1 |
| POLE | 2041 | 0 | 0 | 1 |
| PTPRT | 3415 | 0 | 0 | 1 |
| RAD51C | 322 | 1 | 0 | 0 |
| RHOA | 501 | 0 | 0 | 1 |
| RNF43 | 1974 | 0 | 0 | 1 |
| ROS1 | 2445 | 0 | 0 | 1 |
| SESN1 | 231 | 0 | 0 | 1 |
| SETD2 | 2666 | 0 | 0 | 1 |
| SOS1 | 764 | 0 | 0 | 1 |
| STK19 | 258 | 0 | 0 | 1 |
| TERT | 8046 | 1 | 0 | 0 |
| TGFBR2 | 1102 | 0 | 0 | 1 |
| TP63 | 1155 | 1 | 0 | 0 |
| VHL | 781 | 0 | 0 | 1 |

## Performance

Transcript load: 1.75s from gob cache

| Study | Variants | Sequential | Seq v/s | Parallel | Par v/s | Speedup |
|-------|----------|-----------|---------|----------|---------|--------|
| msk_impact_50k_2026 | 479096 | 20.984s | 22832 | 11.948s | 40098 | 1.76x |
| **Total** | **479096** | **20.984s** | **22832** | **11.948s** | **40098** | **1.76x** |
