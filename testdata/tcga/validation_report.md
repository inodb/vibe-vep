# TCGA Validation Report

Generated: 2026-03-02 02:18 UTC  
GENCODE transcripts loaded: 254070  
Cache load time: 29.881s  
Workers: 4 (GOMAXPROCS)

## Match Rates

| Study | Variants | Conseq Match | Conseq Mismatch | Conseq Rate | HGVSp Match | HGVSp Mismatch | HGVSp Rate | HGVSc Match | HGVSc Mismatch | HGVSc Rate |
|-------|----------|-------------|-----------------|-------------|-------------|----------------|------------|-------------|----------------|------------|
| blca_tcga_gdc | 115850 | 115648 | 1 | 99.8% | 109978 | 0 | 94.9% | 114159 | 0 | 98.5% |
| brca_tcga_gdc | 89012 | 88847 | 1 | 99.8% | 84153 | 0 | 94.5% | 87759 | 0 | 98.6% |
| chol_tcga_gdc | 3764 | 3757 | 0 | 99.8% | 3483 | 0 | 92.5% | 3694 | 0 | 98.1% |
| coad_tcga_gdc | 244552 | 244091 | 2 | 99.8% | 229147 | 0 | 93.7% | 240930 | 0 | 98.5% |
| gbm_tcga_gdc | 54870 | 54747 | 0 | 99.8% | 51902 | 0 | 94.6% | 54128 | 0 | 98.6% |
| luad_tcga_gdc | 190868 | 190442 | 1 | 99.8% | 181696 | 0 | 95.2% | 188188 | 0 | 98.6% |
| skcm_tcga_gdc | 353450 | 352552 | 3 | 99.7% | 336723 | 0 | 95.3% | 348704 | 0 | 98.7% |
| **Total** | **1052366** | **1050084** | **8** | **99.8%** | **997082** | **0** | **94.7%** | **1037562** | **0** | **98.6%** |

## Consequence Category Breakdown

| Study | delins_normalized | gene_model_change | match | mismatch | no_cds_data | position_shift | transcript_model_change | upstream_reclassified |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 1 | 0 | 114749 | 1 | 74 | 42 | 84 | 899 |
| brca_tcga_gdc | 2 | 0 | 88303 | 1 | 78 | 32 | 52 | 544 |
| chol_tcga_gdc | 0 | 0 | 3726 | 0 | 4 | 1 | 2 | 31 |
| coad_tcga_gdc | 8 | 2 | 242410 | 2 | 193 | 99 | 157 | 1681 |
| gbm_tcga_gdc | 0 | 0 | 54376 | 0 | 58 | 33 | 32 | 371 |
| luad_tcga_gdc | 1 | 0 | 189091 | 1 | 219 | 96 | 109 | 1351 |
| skcm_tcga_gdc | 5 | 3 | 349984 | 3 | 398 | 266 | 223 | 2568 |
| **Total** | **17** | **5** | **1042639** | **8** | **1024** | **569** | **659** | **7445** |

## HGVSp Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | fuzzy_fs | maf_empty | maf_nonstandard | match | mismatch | no_cds_data | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 479 | 7 | 1 | 645 | 739 | 2620 | 109978 | 0 | 74 | 505 | 57 | 28 | 611 | 47 | 59 |
| brca_tcga_gdc | 336 | 19 | 4 | 1040 | 447 | 1567 | 84153 | 0 | 79 | 673 | 94 | 45 | 476 | 25 | 54 |
| chol_tcga_gdc | 24 | 2 | 0 | 77 | 23 | 85 | 3483 | 0 | 4 | 31 | 8 | 1 | 21 | 1 | 4 |
| coad_tcga_gdc | 2317 | 8 | 2 | 3024 | 1306 | 5448 | 229147 | 0 | 196 | 1280 | 43 | 34 | 1372 | 88 | 287 |
| gbm_tcga_gdc | 199 | 0 | 2 | 582 | 303 | 1117 | 51902 | 0 | 58 | 303 | 14 | 11 | 326 | 22 | 31 |
| luad_tcga_gdc | 728 | 10 | 0 | 1060 | 1073 | 4118 | 181696 | 0 | 223 | 716 | 53 | 21 | 1005 | 62 | 103 |
| skcm_tcga_gdc | 1414 | 17 | 0 | 493 | 2009 | 7930 | 336723 | 0 | 399 | 1348 | 52 | 21 | 2673 | 125 | 246 |
| **Total** | **5497** | **63** | **9** | **6921** | **5900** | **22885** | **997082** | **0** | **1033** | **4856** | **321** | **161** | **6484** | **370** | **784** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | position_shift | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 81 | 25 | 892 | 114159 | 384 | 299 | 0 |
| brca_tcga_gdc | 8 | 106 | 34 | 538 | 87759 | 416 | 151 | 0 |
| chol_tcga_gdc | 0 | 7 | 3 | 31 | 3694 | 19 | 10 | 0 |
| coad_tcga_gdc | 25 | 90 | 175 | 1653 | 240930 | 1246 | 431 | 2 |
| gbm_tcga_gdc | 3 | 9 | 22 | 369 | 54128 | 242 | 97 | 0 |
| luad_tcga_gdc | 25 | 207 | 45 | 1329 | 188188 | 634 | 440 | 0 |
| skcm_tcga_gdc | 47 | 152 | 26 | 2525 | 348704 | 1362 | 631 | 3 |
| **Total** | **118** | **652** | **330** | **7337** | **1037562** | **4303** | **2059** | **5** |

## Cancer Gene Mismatches

No mismatches across all 1207 cancer genes tested.

## Performance

Cache load time: 29.881s

| Study | Variants | Sequential | Seq v/s | Parallel | Par v/s | Speedup |
|-------|----------|-----------|---------|----------|---------|--------|
| blca_tcga_gdc | 115850 | 23.079s | 5020 | 8.366s | 13848 | 2.76x |
| brca_tcga_gdc | 89012 | 18.8s | 4735 | 8.178s | 10884 | 2.30x |
| chol_tcga_gdc | 3764 | 787ms | 4784 | 310ms | 12132 | 2.54x |
| coad_tcga_gdc | 244552 | 38.123s | 6415 | 15.375s | 15906 | 2.48x |
| gbm_tcga_gdc | 54870 | 10.965s | 5004 | 5.406s | 10150 | 2.03x |
| luad_tcga_gdc | 190868 | 34.102s | 5597 | 14.582s | 13089 | 2.34x |
| skcm_tcga_gdc | 353450 | 54.629s | 6470 | 22.075s | 16011 | 2.47x |
| **Total** | **1052366** | **3m0.486s** | **5831** | **1m14.293s** | **14165** | **2.43x** |
