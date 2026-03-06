# TCGA Validation Report

Generated: 2026-03-06 17:08 UTC  
GENCODE transcripts: 254070 (loaded from gob cache in 9.979s)  
Workers: 4 (GOMAXPROCS)

## Match Rates

| Study | Variants | Conseq Match | Conseq Mismatch | Conseq Rate | HGVSp Match | HGVSp Mismatch | HGVSp Rate | HGVSc Match | HGVSc Mismatch | HGVSc Rate |
|-------|----------|-------------|-----------------|-------------|-------------|----------------|------------|-------------|----------------|------------|
| blca_tcga_gdc | 115850 | 115644 | 1 | 99.8% | 109968 | 0 | 94.9% | 114146 | 0 | 98.5% |
| brca_tcga_gdc | 89012 | 88845 | 1 | 99.8% | 84154 | 0 | 94.5% | 87758 | 0 | 98.6% |
| chol_tcga_gdc | 3764 | 3757 | 0 | 99.8% | 3483 | 0 | 92.5% | 3694 | 0 | 98.1% |
| coad_tcga_gdc | 244552 | 244086 | 3 | 99.8% | 229133 | 0 | 93.7% | 240914 | 0 | 98.5% |
| gbm_tcga_gdc | 54870 | 54746 | 0 | 99.8% | 51886 | 0 | 94.6% | 54113 | 0 | 98.6% |
| luad_tcga_gdc | 190868 | 190438 | 2 | 99.8% | 181662 | 0 | 95.2% | 188152 | 0 | 98.6% |
| skcm_tcga_gdc | 353450 | 352549 | 3 | 99.7% | 336647 | 0 | 95.2% | 348627 | 0 | 98.6% |
| **Total** | **1052366** | **1050065** | **10** | **99.8%** | **996933** | **0** | **94.7%** | **1037404** | **0** | **98.6%** |

## Consequence Category Breakdown

| Study | delins_normalized | gene_model_change | match | mismatch | no_cds_data | position_shift | transcript_model_change | upstream_reclassified |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 1 | 2 | 114744 | 1 | 74 | 42 | 86 | 900 |
| brca_tcga_gdc | 2 | 0 | 88301 | 1 | 78 | 32 | 54 | 544 |
| chol_tcga_gdc | 0 | 0 | 3726 | 0 | 4 | 1 | 2 | 31 |
| coad_tcga_gdc | 8 | 3 | 242405 | 3 | 193 | 99 | 160 | 1681 |
| gbm_tcga_gdc | 0 | 0 | 54375 | 0 | 58 | 33 | 33 | 371 |
| luad_tcga_gdc | 1 | 2 | 189087 | 2 | 219 | 96 | 110 | 1351 |
| skcm_tcga_gdc | 5 | 5 | 349981 | 3 | 398 | 266 | 224 | 2568 |
| **Total** | **17** | **12** | **1042619** | **10** | **1024** | **569** | **669** | **7446** |

## HGVSp Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | fuzzy_fs | maf_empty | maf_nonstandard | match | no_cds_data | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 479 | 7 | 1 | 645 | 739 | 2620 | 109968 | 74 | 513 | 57 | 28 | 611 | 47 | 61 |
| brca_tcga_gdc | 336 | 19 | 4 | 1040 | 447 | 1567 | 84154 | 79 | 672 | 94 | 45 | 476 | 25 | 54 |
| chol_tcga_gdc | 24 | 2 | 0 | 77 | 23 | 85 | 3483 | 4 | 31 | 8 | 1 | 21 | 1 | 4 |
| coad_tcga_gdc | 2317 | 8 | 2 | 3023 | 1306 | 5448 | 229133 | 196 | 1294 | 43 | 34 | 1372 | 88 | 288 |
| gbm_tcga_gdc | 199 | 0 | 2 | 583 | 303 | 1117 | 51886 | 58 | 318 | 14 | 11 | 326 | 22 | 31 |
| luad_tcga_gdc | 728 | 10 | 0 | 1060 | 1073 | 4118 | 181662 | 223 | 747 | 53 | 21 | 1005 | 62 | 106 |
| skcm_tcga_gdc | 1414 | 17 | 0 | 493 | 2009 | 7930 | 336647 | 399 | 1422 | 52 | 21 | 2673 | 125 | 248 |
| **Total** | **5497** | **63** | **9** | **6921** | **5900** | **22885** | **996933** | **1033** | **4997** | **321** | **161** | **6484** | **370** | **792** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | position_shift | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 81 | 25 | 892 | 114146 | 393 | 301 | 2 |
| brca_tcga_gdc | 8 | 106 | 34 | 538 | 87758 | 417 | 151 | 0 |
| chol_tcga_gdc | 0 | 7 | 3 | 31 | 3694 | 19 | 10 | 0 |
| coad_tcga_gdc | 25 | 90 | 175 | 1653 | 240914 | 1260 | 431 | 4 |
| gbm_tcga_gdc | 3 | 9 | 22 | 369 | 54113 | 257 | 97 | 0 |
| luad_tcga_gdc | 25 | 207 | 45 | 1329 | 188152 | 667 | 440 | 3 |
| skcm_tcga_gdc | 47 | 152 | 26 | 2525 | 348627 | 1437 | 631 | 5 |
| **Total** | **118** | **652** | **330** | **7337** | **1037404** | **4450** | **2061** | **14** |

## Cancer Gene Mismatches

No mismatches across all 1207 cancer genes tested.

## Performance

Transcript load: 9.979s from gob cache

| Study | Variants | Sequential | Seq v/s | Parallel | Par v/s | Speedup |
|-------|----------|-----------|---------|----------|---------|--------|
| blca_tcga_gdc | 115850 | 10.93s | 10599 | 7.562s | 15321 | 1.45x |
| brca_tcga_gdc | 89012 | 8.643s | 10298 | 4.298s | 20709 | 2.01x |
| chol_tcga_gdc | 3764 | 305ms | 12343 | 291ms | 12954 | 1.05x |
| coad_tcga_gdc | 244552 | 21.958s | 11137 | 12.208s | 20032 | 1.80x |
| gbm_tcga_gdc | 54870 | 5.128s | 10699 | 2.396s | 22897 | 2.14x |
| luad_tcga_gdc | 190868 | 15.938s | 11976 | 9.569s | 19947 | 1.67x |
| skcm_tcga_gdc | 353450 | 30.943s | 11423 | 22.407s | 15774 | 1.38x |
| **Total** | **1052366** | **1m33.846s** | **11214** | **58.73s** | **17919** | **1.60x** |
