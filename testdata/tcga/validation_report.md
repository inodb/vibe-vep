# TCGA Validation Report

Generated: 2026-03-02 18:40 UTC  
GENCODE transcripts loaded: 254070  
Cache load time: 42.796s  
Workers: 4 (GOMAXPROCS)

## Match Rates

| Study | Variants | Conseq Match | Conseq Mismatch | Conseq Rate | HGVSp Match | HGVSp Mismatch | HGVSp Rate | HGVSc Match | HGVSc Mismatch | HGVSc Rate |
|-------|----------|-------------|-----------------|-------------|-------------|----------------|------------|-------------|----------------|------------|
| blca_tcga_gdc | 115850 | 115646 | 1 | 99.8% | 109978 | 0 | 94.9% | 114160 | 0 | 98.5% |
| brca_tcga_gdc | 89012 | 88845 | 1 | 99.8% | 84155 | 0 | 94.5% | 87762 | 0 | 98.6% |
| chol_tcga_gdc | 3764 | 3757 | 0 | 99.8% | 3481 | 0 | 92.5% | 3692 | 0 | 98.1% |
| coad_tcga_gdc | 244552 | 244091 | 2 | 99.8% | 229131 | 0 | 93.7% | 240913 | 0 | 98.5% |
| gbm_tcga_gdc | 54870 | 54748 | 0 | 99.8% | 51888 | 0 | 94.6% | 54114 | 0 | 98.6% |
| luad_tcga_gdc | 190868 | 190445 | 1 | 99.8% | 181684 | 0 | 95.2% | 188178 | 0 | 98.6% |
| skcm_tcga_gdc | 353450 | 352554 | 3 | 99.7% | 336681 | 0 | 95.3% | 348662 | 0 | 98.6% |
| **Total** | **1052366** | **1050086** | **8** | **99.8%** | **996998** | **0** | **94.7%** | **1037481** | **0** | **98.6%** |

## Consequence Category Breakdown

| Study | delins_normalized | gene_model_change | match | mismatch | no_cds_data | position_shift | transcript_model_change | upstream_reclassified |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 1 | 0 | 114747 | 1 | 74 | 42 | 86 | 899 |
| brca_tcga_gdc | 2 | 0 | 88300 | 1 | 78 | 32 | 54 | 545 |
| chol_tcga_gdc | 0 | 0 | 3726 | 0 | 4 | 1 | 2 | 31 |
| coad_tcga_gdc | 8 | 2 | 242410 | 2 | 193 | 99 | 157 | 1681 |
| gbm_tcga_gdc | 0 | 0 | 54377 | 0 | 58 | 33 | 31 | 371 |
| luad_tcga_gdc | 1 | 0 | 189094 | 1 | 219 | 96 | 106 | 1351 |
| skcm_tcga_gdc | 5 | 3 | 349986 | 3 | 398 | 266 | 221 | 2568 |
| **Total** | **17** | **5** | **1042640** | **8** | **1024** | **569** | **657** | **7446** |

## HGVSp Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | fuzzy_fs | maf_empty | maf_nonstandard | match | mismatch | no_cds_data | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 479 | 7 | 1 | 645 | 739 | 2620 | 109978 | 0 | 74 | 505 | 57 | 28 | 611 | 47 | 59 |
| brca_tcga_gdc | 336 | 19 | 4 | 1040 | 447 | 1567 | 84155 | 0 | 79 | 671 | 94 | 45 | 476 | 25 | 54 |
| chol_tcga_gdc | 24 | 2 | 0 | 77 | 23 | 85 | 3481 | 0 | 4 | 33 | 8 | 1 | 21 | 1 | 4 |
| coad_tcga_gdc | 2317 | 8 | 2 | 3023 | 1306 | 5448 | 229131 | 0 | 196 | 1297 | 43 | 34 | 1372 | 88 | 287 |
| gbm_tcga_gdc | 199 | 0 | 2 | 582 | 303 | 1117 | 51888 | 0 | 58 | 317 | 14 | 11 | 326 | 22 | 31 |
| luad_tcga_gdc | 728 | 10 | 0 | 1060 | 1073 | 4118 | 181684 | 0 | 223 | 728 | 53 | 21 | 1005 | 62 | 103 |
| skcm_tcga_gdc | 1414 | 17 | 0 | 493 | 2009 | 7930 | 336681 | 0 | 399 | 1390 | 52 | 21 | 2673 | 125 | 246 |
| **Total** | **5497** | **63** | **9** | **6920** | **5900** | **22885** | **996998** | **0** | **1033** | **4941** | **321** | **161** | **6484** | **370** | **784** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | position_shift | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 81 | 25 | 892 | 114160 | 384 | 298 | 0 |
| brca_tcga_gdc | 8 | 106 | 34 | 538 | 87762 | 415 | 149 | 0 |
| chol_tcga_gdc | 0 | 7 | 3 | 31 | 3692 | 21 | 10 | 0 |
| coad_tcga_gdc | 25 | 90 | 175 | 1653 | 240913 | 1263 | 431 | 2 |
| gbm_tcga_gdc | 3 | 9 | 22 | 369 | 54114 | 256 | 97 | 0 |
| luad_tcga_gdc | 25 | 207 | 45 | 1329 | 188178 | 646 | 438 | 0 |
| skcm_tcga_gdc | 47 | 152 | 26 | 2525 | 348662 | 1405 | 630 | 3 |
| **Total** | **118** | **652** | **330** | **7337** | **1037481** | **4390** | **2053** | **5** |

## Cancer Gene Mismatches

No mismatches across all 1207 cancer genes tested.

## Performance

Cache load time: 42.796s

| Study | Variants | Sequential | Seq v/s | Parallel | Par v/s | Speedup |
|-------|----------|-----------|---------|----------|---------|--------|
| blca_tcga_gdc | 115850 | 29.008s | 3994 | 10.353s | 11190 | 2.80x |
| brca_tcga_gdc | 89012 | 25.277s | 3521 | 14.458s | 6156 | 1.75x |
| chol_tcga_gdc | 3764 | 957ms | 3933 | 552ms | 6824 | 1.74x |
| coad_tcga_gdc | 244552 | 52.335s | 4673 | 24.265s | 10078 | 2.16x |
| gbm_tcga_gdc | 54870 | 15.22s | 3605 | 8.023s | 6839 | 1.90x |
| luad_tcga_gdc | 190868 | 47.233s | 4041 | 25.711s | 7423 | 1.84x |
| skcm_tcga_gdc | 353450 | 1m21.352s | 4345 | 37.616s | 9396 | 2.16x |
| **Total** | **1052366** | **4m11.382s** | **4186** | **2m0.979s** | **8699** | **2.08x** |
