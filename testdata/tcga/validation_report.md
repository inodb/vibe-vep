# TCGA Validation Report

Generated: 2026-02-27 20:53 UTC  
GENCODE transcripts loaded: 254070  
Cache load time: 43.148s  
Workers: 4 (GOMAXPROCS)

## Match Rates

| Study | Variants | Conseq Match | Conseq Mismatch | Conseq Rate | HGVSp Match | HGVSp Mismatch | HGVSp Rate | HGVSc Match | HGVSc Mismatch | HGVSc Rate |
|-------|----------|-------------|-----------------|-------------|-------------|----------------|------------|-------------|----------------|------------|
| blca_tcga_gdc | 115850 | 115646 | 8 | 99.8% | 109990 | 34 | 94.9% | 114169 | 40 | 98.5% |
| brca_tcga_gdc | 89012 | 88844 | 6 | 99.8% | 84159 | 19 | 94.5% | 87764 | 20 | 98.6% |
| chol_tcga_gdc | 3764 | 3757 | 0 | 99.8% | 3484 | 0 | 92.6% | 3695 | 0 | 98.2% |
| coad_tcga_gdc | 244552 | 244080 | 16 | 99.8% | 229172 | 60 | 93.7% | 240955 | 74 | 98.5% |
| gbm_tcga_gdc | 54870 | 54746 | 3 | 99.8% | 51903 | 12 | 94.6% | 54130 | 11 | 98.7% |
| luad_tcga_gdc | 190868 | 190442 | 8 | 99.8% | 181720 | 46 | 95.2% | 188213 | 55 | 98.6% |
| skcm_tcga_gdc | 353450 | 352549 | 9 | 99.7% | 336742 | 67 | 95.3% | 348724 | 78 | 98.7% |
| **Total** | **1052366** | **1050064** | **50** | **99.8%** | **997170** | **238** | **94.8%** | **1037650** | **278** | **98.6%** |

## Consequence Category Breakdown

| Study | delins_normalized | gene_model_change | match | mismatch | no_cds_data | position_shift | transcript_model_change | upstream_reclassified |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 1 | 0 | 114747 | 8 | 74 | 63 | 58 | 899 |
| brca_tcga_gdc | 2 | 0 | 88299 | 6 | 78 | 45 | 37 | 545 |
| chol_tcga_gdc | 0 | 0 | 3726 | 0 | 4 | 1 | 2 | 31 |
| coad_tcga_gdc | 8 | 2 | 242407 | 16 | 193 | 135 | 118 | 1673 |
| gbm_tcga_gdc | 0 | 0 | 54375 | 3 | 58 | 36 | 27 | 371 |
| luad_tcga_gdc | 1 | 0 | 189092 | 8 | 219 | 117 | 81 | 1350 |
| skcm_tcga_gdc | 5 | 2 | 349982 | 9 | 398 | 309 | 178 | 2567 |
| **Total** | **17** | **4** | **1042628** | **50** | **1024** | **706** | **501** | **7436** |

## HGVSp Category Breakdown

| Study | both_empty | delins_normalized | fuzzy_fs | maf_empty | maf_nonstandard | match | mismatch | no_cds_data | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 479 | 7 | 645 | 739 | 2620 | 109990 | 34 | 74 | 507 | 57 | 28 | 611 | 59 |
| brca_tcga_gdc | 336 | 19 | 1040 | 447 | 1567 | 84159 | 19 | 79 | 677 | 94 | 45 | 476 | 54 |
| chol_tcga_gdc | 24 | 2 | 77 | 23 | 85 | 3484 | 0 | 4 | 31 | 8 | 1 | 21 | 4 |
| coad_tcga_gdc | 2317 | 9 | 3023 | 1306 | 5448 | 229172 | 60 | 196 | 1285 | 43 | 34 | 1372 | 287 |
| gbm_tcga_gdc | 199 | 0 | 583 | 303 | 1117 | 51903 | 12 | 58 | 313 | 14 | 11 | 326 | 31 |
| luad_tcga_gdc | 728 | 10 | 1060 | 1073 | 4118 | 181720 | 46 | 223 | 708 | 53 | 21 | 1005 | 103 |
| skcm_tcga_gdc | 1414 | 17 | 493 | 2009 | 7930 | 336742 | 67 | 399 | 1387 | 52 | 21 | 2673 | 246 |
| **Total** | **5497** | **64** | **6921** | **5900** | **22885** | **997170** | **238** | **1033** | **4908** | **321** | **161** | **6484** | **784** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | mismatch | position_shift | vep_empty |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 81 | 33 | 892 | 114169 | 40 | 625 | 0 |
| brca_tcga_gdc | 8 | 106 | 40 | 538 | 87764 | 20 | 536 | 0 |
| chol_tcga_gdc | 0 | 7 | 3 | 31 | 3695 | 0 | 28 | 0 |
| coad_tcga_gdc | 25 | 90 | 217 | 1653 | 240955 | 74 | 1536 | 2 |
| gbm_tcga_gdc | 3 | 11 | 27 | 369 | 54130 | 11 | 319 | 0 |
| luad_tcga_gdc | 25 | 208 | 69 | 1329 | 188213 | 55 | 969 | 0 |
| skcm_tcga_gdc | 47 | 155 | 35 | 2525 | 348724 | 78 | 1883 | 3 |
| **Total** | **118** | **658** | **424** | **7337** | **1037650** | **278** | **5896** | **5** |

## Cancer Gene Mismatches

| Gene | Variants | Conseq Mismatches | HGVSp Mismatches | HGVSc Mismatches |
|------|----------|-------------------|------------------|------------------|
| CTCF | 87 | 1 | 0 | 0 |
| PPFIBP1 | 66 | 0 | 0 | 1 |

## Performance

Cache load time: 43.148s

| Study | Variants | Sequential | Seq v/s | Parallel | Par v/s | Speedup |
|-------|----------|-----------|---------|----------|---------|--------|
| blca_tcga_gdc | 115850 | 33.267s | 3482 | 14.475s | 8004 | 2.30x |
| brca_tcga_gdc | 89012 | 29.24s | 3044 | 10.561s | 8428 | 2.77x |
| chol_tcga_gdc | 3764 | 1.507s | 2498 | 591ms | 6366 | 2.55x |
| coad_tcga_gdc | 244552 | 59.672s | 4098 | 25.257s | 9682 | 2.36x |
| gbm_tcga_gdc | 54870 | 16.716s | 3282 | 6.184s | 8874 | 2.70x |
| luad_tcga_gdc | 190868 | 58.533s | 3261 | 22.132s | 8624 | 2.64x |
| skcm_tcga_gdc | 353450 | 1m15.353s | 4691 | 34.813s | 10153 | 2.16x |
| **Total** | **1052366** | **4m34.288s** | **3837** | **1m54.013s** | **9230** | **2.41x** |
