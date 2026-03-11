# GRCh38 Validation Report

Generated: 2026-03-11 14:26 UTC  
Assembly: GRCh38  
GENCODE transcripts: 254070 (loaded from gob cache in 4.3s)  
Workers: 4 (GOMAXPROCS)

## Match Rates

| Study | Variants | Conseq Match | Conseq Mismatch | Conseq Rate | HGVSp Match | HGVSp Mismatch | HGVSp Rate | HGVSc Match | HGVSc Mismatch | HGVSc Rate |
|-------|----------|-------------|-----------------|-------------|-------------|----------------|------------|-------------|----------------|------------|
| blca_tcga_gdc | 115850 | 115638 | 1 | 99.8% | 109969 | 0 | 94.9% | 114146 | 0 | 98.5% |
| brca_tcga_gdc | 89012 | 88848 | 1 | 99.8% | 84155 | 0 | 94.5% | 87759 | 0 | 98.6% |
| chol_tcga_gdc | 3764 | 3758 | 0 | 99.8% | 3483 | 0 | 92.5% | 3694 | 0 | 98.1% |
| coad_tcga_gdc | 244552 | 244084 | 3 | 99.8% | 229134 | 0 | 93.7% | 240915 | 0 | 98.5% |
| gbm_tcga_gdc | 54870 | 54740 | 0 | 99.8% | 51886 | 0 | 94.6% | 54112 | 0 | 98.6% |
| luad_tcga_gdc | 190868 | 190435 | 2 | 99.8% | 181666 | 0 | 95.2% | 188156 | 0 | 98.6% |
| skcm_tcga_gdc | 353450 | 352553 | 3 | 99.7% | 336649 | 0 | 95.2% | 348622 | 0 | 98.6% |
| **Total** | **1052366** | **1050056** | **10** | **99.8%** | **996942** | **0** | **94.7%** | **1037404** | **0** | **98.6%** |

## Consequence Category Breakdown

| Study | delins_normalized | gene_model_change | match | mismatch | no_cds_data | position_shift | transcript_model_change | upstream_reclassified |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 1 | 2 | 114738 | 1 | 74 | 42 | 92 | 900 |
| brca_tcga_gdc | 2 | 0 | 88304 | 1 | 78 | 32 | 51 | 544 |
| chol_tcga_gdc | 0 | 0 | 3727 | 0 | 4 | 1 | 1 | 31 |
| coad_tcga_gdc | 8 | 3 | 242403 | 3 | 193 | 99 | 162 | 1681 |
| gbm_tcga_gdc | 0 | 0 | 54369 | 0 | 58 | 33 | 39 | 371 |
| luad_tcga_gdc | 1 | 2 | 189084 | 2 | 219 | 96 | 113 | 1351 |
| skcm_tcga_gdc | 5 | 5 | 349985 | 3 | 398 | 267 | 219 | 2568 |
| **Total** | **17** | **12** | **1042610** | **10** | **1024** | **570** | **677** | **7446** |

## HGVSp Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | fuzzy_fs | maf_empty | maf_nonstandard | match | no_cds_data | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 460 | 7 | 1 | 645 | 758 | 2607 | 109969 | 74 | 518 | 57 | 28 | 611 | 60 | 55 |
| brca_tcga_gdc | 322 | 19 | 4 | 1041 | 461 | 1565 | 84155 | 79 | 674 | 94 | 45 | 476 | 27 | 50 |
| chol_tcga_gdc | 24 | 2 | 0 | 77 | 23 | 84 | 3483 | 4 | 33 | 8 | 1 | 21 | 2 | 2 |
| coad_tcga_gdc | 2255 | 8 | 2 | 3024 | 1368 | 5441 | 229134 | 196 | 1300 | 43 | 34 | 1373 | 95 | 279 |
| gbm_tcga_gdc | 186 | 0 | 2 | 583 | 316 | 1114 | 51886 | 58 | 320 | 14 | 11 | 326 | 25 | 29 |
| luad_tcga_gdc | 698 | 10 | 0 | 1060 | 1103 | 4113 | 181666 | 223 | 750 | 53 | 21 | 1005 | 67 | 99 |
| skcm_tcga_gdc | 1348 | 17 | 0 | 493 | 2075 | 7923 | 336649 | 399 | 1431 | 52 | 21 | 2673 | 132 | 237 |
| **Total** | **5293** | **63** | **9** | **6923** | **6104** | **22847** | **996942** | **1033** | **5026** | **321** | **161** | **6485** | **408** | **751** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | position_shift | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 81 | 25 | 892 | 114146 | 398 | 296 | 2 |
| brca_tcga_gdc | 8 | 106 | 34 | 538 | 87759 | 420 | 147 | 0 |
| chol_tcga_gdc | 0 | 7 | 3 | 31 | 3694 | 21 | 8 | 0 |
| coad_tcga_gdc | 25 | 90 | 175 | 1653 | 240915 | 1268 | 422 | 4 |
| gbm_tcga_gdc | 3 | 9 | 22 | 369 | 54112 | 259 | 96 | 0 |
| luad_tcga_gdc | 25 | 207 | 45 | 1329 | 188156 | 670 | 433 | 3 |
| skcm_tcga_gdc | 47 | 152 | 26 | 2525 | 348622 | 1446 | 627 | 5 |
| **Total** | **118** | **652** | **330** | **7337** | **1037404** | **4482** | **2029** | **14** |

## Cancer Gene Mismatches

No mismatches across all 1207 cancer genes tested.

## Performance

Transcript load: 4.3s from gob cache

| Study | Variants | Sequential | Seq v/s | Parallel | Par v/s | Speedup |
|-------|----------|-----------|---------|----------|---------|--------|
| blca_tcga_gdc | 115850 | 7.877s | 14707 | 6.924s | 16730 | 1.14x |
| brca_tcga_gdc | 89012 | 6.457s | 13785 | 3.957s | 22497 | 1.63x |
| chol_tcga_gdc | 3764 | 246ms | 15272 | 164ms | 22972 | 1.50x |
| coad_tcga_gdc | 244552 | 15.331s | 15951 | 11.922s | 20513 | 1.29x |
| gbm_tcga_gdc | 54870 | 3.497s | 15691 | 1.562s | 35139 | 2.24x |
| luad_tcga_gdc | 190868 | 10.759s | 17740 | 8.331s | 22911 | 1.29x |
| skcm_tcga_gdc | 353450 | 24.809s | 14247 | 13.571s | 26045 | 1.83x |
| **Total** | **1052366** | **1m8.977s** | **15257** | **46.43s** | **22665** | **1.49x** |
