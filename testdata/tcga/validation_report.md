# TCGA Validation Report

Generated: 2026-02-27 21:19 UTC  
GENCODE transcripts loaded: 254070  
Cache load time: 43.679s  
Workers: 4 (GOMAXPROCS)

## Match Rates

| Study | Variants | Conseq Match | Conseq Mismatch | Conseq Rate | HGVSp Match | HGVSp Mismatch | HGVSp Rate | HGVSc Match | HGVSc Mismatch | HGVSc Rate |
|-------|----------|-------------|-----------------|-------------|-------------|----------------|------------|-------------|----------------|------------|
| blca_tcga_gdc | 115850 | 115650 | 8 | 99.8% | 109995 | 34 | 94.9% | 114177 | 40 | 98.6% |
| brca_tcga_gdc | 89012 | 88848 | 5 | 99.8% | 84169 | 19 | 94.6% | 87777 | 20 | 98.6% |
| chol_tcga_gdc | 3764 | 3757 | 0 | 99.8% | 3484 | 0 | 92.6% | 3695 | 0 | 98.2% |
| coad_tcga_gdc | 244552 | 244081 | 16 | 99.8% | 229190 | 60 | 93.7% | 240973 | 74 | 98.5% |
| gbm_tcga_gdc | 54870 | 54745 | 3 | 99.8% | 51914 | 12 | 94.6% | 54140 | 11 | 98.7% |
| luad_tcga_gdc | 190868 | 190442 | 9 | 99.8% | 181738 | 46 | 95.2% | 188231 | 55 | 98.6% |
| skcm_tcga_gdc | 353450 | 352553 | 9 | 99.7% | 336796 | 67 | 95.3% | 348780 | 78 | 98.7% |
| **Total** | **1052366** | **1050076** | **50** | **99.8%** | **997286** | **238** | **94.8%** | **1037773** | **278** | **98.6%** |

## Consequence Category Breakdown

| Study | delins_normalized | gene_model_change | match | mismatch | no_cds_data | position_shift | transcript_model_change | upstream_reclassified |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 1 | 0 | 114751 | 8 | 74 | 59 | 58 | 899 |
| brca_tcga_gdc | 2 | 0 | 88304 | 5 | 78 | 42 | 37 | 544 |
| chol_tcga_gdc | 0 | 0 | 3726 | 0 | 4 | 1 | 2 | 31 |
| coad_tcga_gdc | 8 | 2 | 242408 | 16 | 193 | 134 | 118 | 1673 |
| gbm_tcga_gdc | 0 | 0 | 54374 | 3 | 58 | 37 | 27 | 371 |
| luad_tcga_gdc | 1 | 0 | 189092 | 9 | 219 | 116 | 81 | 1350 |
| skcm_tcga_gdc | 5 | 2 | 349986 | 9 | 398 | 305 | 178 | 2567 |
| **Total** | **17** | **4** | **1042641** | **50** | **1024** | **694** | **501** | **7435** |

## HGVSp Category Breakdown

| Study | both_empty | delins_normalized | fuzzy_fs | maf_empty | maf_nonstandard | match | mismatch | no_cds_data | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 479 | 7 | 645 | 739 | 2620 | 109995 | 34 | 74 | 502 | 57 | 28 | 611 | 59 |
| brca_tcga_gdc | 336 | 19 | 1040 | 447 | 1567 | 84169 | 19 | 79 | 667 | 94 | 45 | 476 | 54 |
| chol_tcga_gdc | 24 | 2 | 77 | 23 | 85 | 3484 | 0 | 4 | 31 | 8 | 1 | 21 | 4 |
| coad_tcga_gdc | 2317 | 9 | 3024 | 1306 | 5448 | 229190 | 60 | 196 | 1266 | 43 | 34 | 1372 | 287 |
| gbm_tcga_gdc | 199 | 0 | 582 | 303 | 1117 | 51914 | 12 | 58 | 303 | 14 | 11 | 326 | 31 |
| luad_tcga_gdc | 728 | 10 | 1060 | 1073 | 4118 | 181738 | 46 | 223 | 690 | 53 | 21 | 1005 | 103 |
| skcm_tcga_gdc | 1414 | 17 | 493 | 2009 | 7930 | 336796 | 67 | 399 | 1333 | 52 | 21 | 2673 | 246 |
| **Total** | **5497** | **64** | **6921** | **5900** | **22885** | **997286** | **238** | **1033** | **4792** | **321** | **161** | **6484** | **784** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | mismatch | position_shift | vep_empty |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 81 | 33 | 892 | 114177 | 40 | 617 | 0 |
| brca_tcga_gdc | 8 | 106 | 40 | 538 | 87777 | 20 | 523 | 0 |
| chol_tcga_gdc | 0 | 7 | 3 | 31 | 3695 | 0 | 28 | 0 |
| coad_tcga_gdc | 25 | 90 | 217 | 1653 | 240973 | 74 | 1518 | 2 |
| gbm_tcga_gdc | 3 | 11 | 27 | 369 | 54140 | 11 | 309 | 0 |
| luad_tcga_gdc | 25 | 208 | 69 | 1329 | 188231 | 55 | 951 | 0 |
| skcm_tcga_gdc | 47 | 155 | 35 | 2525 | 348780 | 78 | 1827 | 3 |
| **Total** | **118** | **658** | **424** | **7337** | **1037773** | **278** | **5773** | **5** |

## Cancer Gene Mismatches

1205/1207 cancer genes have 100% match across all columns. Mismatches in 2 gene(s):

| Gene | Variants | Conseq Mismatches | HGVSp Mismatches | HGVSc Mismatches |
|------|----------|-------------------|------------------|------------------|
| CTCF | 87 | 1 | 0 | 0 |
| PPFIBP1 | 66 | 0 | 0 | 1 |

## Performance

Cache load time: 43.679s

| Study | Variants | Sequential | Seq v/s | Parallel | Par v/s | Speedup |
|-------|----------|-----------|---------|----------|---------|--------|
| blca_tcga_gdc | 115850 | 30.018s | 3859 | 12.162s | 9525 | 2.47x |
| brca_tcga_gdc | 89012 | 27.645s | 3220 | 9.694s | 9182 | 2.85x |
| chol_tcga_gdc | 3764 | 1.503s | 2505 | 396ms | 9501 | 3.79x |
| coad_tcga_gdc | 244552 | 59.163s | 4134 | 26.256s | 9314 | 2.25x |
| gbm_tcga_gdc | 54870 | 15.149s | 3622 | 5.515s | 9950 | 2.75x |
| luad_tcga_gdc | 190868 | 50.091s | 3810 | 19.909s | 9587 | 2.52x |
| skcm_tcga_gdc | 353450 | 1m24.082s | 4204 | 33.529s | 10542 | 2.51x |
| **Total** | **1052366** | **4m27.651s** | **3932** | **1m47.46s** | **9793** | **2.49x** |
