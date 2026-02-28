# TCGA Validation Report

Generated: 2026-02-28 03:28 UTC  
GENCODE transcripts loaded: 254070  
Cache load time: 33.18s  
Workers: 4 (GOMAXPROCS)

## Match Rates

| Study | Variants | Conseq Match | Conseq Mismatch | Conseq Rate | HGVSp Match | HGVSp Mismatch | HGVSp Rate | HGVSc Match | HGVSc Mismatch | HGVSc Rate |
|-------|----------|-------------|-----------------|-------------|-------------|----------------|------------|-------------|----------------|------------|
| blca_tcga_gdc | 115850 | 115651 | 1 | 99.8% | 110000 | 0 | 95.0% | 114183 | 0 | 98.6% |
| brca_tcga_gdc | 89012 | 88849 | 1 | 99.8% | 84176 | 0 | 94.6% | 87784 | 0 | 98.6% |
| chol_tcga_gdc | 3764 | 3757 | 0 | 99.8% | 3483 | 0 | 92.5% | 3694 | 0 | 98.1% |
| coad_tcga_gdc | 244552 | 244089 | 2 | 99.8% | 229211 | 0 | 93.7% | 240995 | 0 | 98.5% |
| gbm_tcga_gdc | 54870 | 54748 | 0 | 99.8% | 51913 | 0 | 94.6% | 54140 | 0 | 98.7% |
| luad_tcga_gdc | 190868 | 190446 | 1 | 99.8% | 181741 | 0 | 95.2% | 188235 | 0 | 98.6% |
| skcm_tcga_gdc | 353450 | 352555 | 3 | 99.7% | 336835 | 0 | 95.3% | 348817 | 0 | 98.7% |
| **Total** | **1052366** | **1050095** | **8** | **99.8%** | **997359** | **0** | **94.8%** | **1037848** | **0** | **98.6%** |

## Consequence Category Breakdown

| Study | delins_normalized | gene_model_change | match | mismatch | no_cds_data | position_shift | transcript_model_change | upstream_reclassified |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 1 | 0 | 114752 | 1 | 74 | 42 | 81 | 899 |
| brca_tcga_gdc | 2 | 0 | 88305 | 1 | 78 | 32 | 50 | 544 |
| chol_tcga_gdc | 0 | 0 | 3726 | 0 | 4 | 1 | 2 | 31 |
| coad_tcga_gdc | 8 | 2 | 242408 | 2 | 193 | 99 | 159 | 1681 |
| gbm_tcga_gdc | 0 | 0 | 54377 | 0 | 58 | 33 | 31 | 371 |
| luad_tcga_gdc | 1 | 0 | 189095 | 1 | 219 | 96 | 105 | 1351 |
| skcm_tcga_gdc | 5 | 3 | 349987 | 3 | 398 | 266 | 220 | 2568 |
| **Total** | **17** | **5** | **1042650** | **8** | **1024** | **569** | **648** | **7445** |

## HGVSp Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | fuzzy_fs | maf_empty | maf_nonstandard | match | mismatch | no_cds_data | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 479 | 7 | 1 | 645 | 739 | 2620 | 110000 | 0 | 74 | 483 | 57 | 28 | 611 | 47 | 59 |
| brca_tcga_gdc | 336 | 19 | 4 | 1040 | 447 | 1567 | 84176 | 0 | 79 | 650 | 94 | 45 | 476 | 25 | 54 |
| chol_tcga_gdc | 24 | 2 | 0 | 77 | 23 | 85 | 3483 | 0 | 4 | 31 | 8 | 1 | 21 | 1 | 4 |
| coad_tcga_gdc | 2317 | 8 | 2 | 3023 | 1306 | 5448 | 229211 | 0 | 196 | 1217 | 43 | 34 | 1372 | 88 | 287 |
| gbm_tcga_gdc | 199 | 0 | 2 | 583 | 303 | 1117 | 51913 | 0 | 58 | 291 | 14 | 11 | 326 | 22 | 31 |
| luad_tcga_gdc | 728 | 10 | 0 | 1060 | 1073 | 4118 | 181741 | 0 | 223 | 671 | 53 | 21 | 1005 | 62 | 103 |
| skcm_tcga_gdc | 1414 | 17 | 0 | 493 | 2009 | 7930 | 336835 | 0 | 399 | 1236 | 52 | 21 | 2673 | 125 | 246 |
| **Total** | **5497** | **63** | **9** | **6921** | **5900** | **22885** | **997359** | **0** | **1033** | **4579** | **321** | **161** | **6484** | **370** | **784** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | position_shift | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 81 | 25 | 892 | 114183 | 361 | 298 | 0 |
| brca_tcga_gdc | 8 | 106 | 34 | 538 | 87784 | 392 | 150 | 0 |
| chol_tcga_gdc | 0 | 7 | 3 | 31 | 3694 | 19 | 10 | 0 |
| coad_tcga_gdc | 25 | 90 | 175 | 1653 | 240995 | 1181 | 431 | 2 |
| gbm_tcga_gdc | 3 | 9 | 22 | 369 | 54140 | 230 | 97 | 0 |
| luad_tcga_gdc | 25 | 207 | 45 | 1329 | 188235 | 588 | 439 | 0 |
| skcm_tcga_gdc | 47 | 152 | 26 | 2525 | 348817 | 1250 | 630 | 3 |
| **Total** | **118** | **652** | **330** | **7337** | **1037848** | **4021** | **2055** | **5** |

## Cancer Gene Mismatches

No mismatches across all 1207 cancer genes tested.

## Performance

Cache load time: 33.18s

| Study | Variants | Sequential | Seq v/s | Parallel | Par v/s | Speedup |
|-------|----------|-----------|---------|----------|---------|--------|
| blca_tcga_gdc | 115850 | 26.878s | 4310 | 10.878s | 10650 | 2.47x |
| brca_tcga_gdc | 89012 | 22.41s | 3972 | 8.574s | 10381 | 2.61x |
| chol_tcga_gdc | 3764 | 792ms | 4751 | 335ms | 11247 | 2.37x |
| coad_tcga_gdc | 244552 | 46.057s | 5310 | 15.774s | 15504 | 2.92x |
| gbm_tcga_gdc | 54870 | 12.969s | 4231 | 4.206s | 13046 | 3.08x |
| luad_tcga_gdc | 190868 | 39.5s | 4832 | 13.425s | 14218 | 2.94x |
| skcm_tcga_gdc | 353450 | 1m2.689s | 5638 | 23.899s | 14789 | 2.62x |
| **Total** | **1052366** | **3m31.295s** | **4981** | **1m17.09s** | **13651** | **2.74x** |
