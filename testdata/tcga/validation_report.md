# TCGA Validation Report

Generated: 2026-02-27 18:20 UTC  
GENCODE transcripts loaded: 254070  
Cache load time: 42.605s  
Workers: 4 (GOMAXPROCS)

## Match Rates

| Study | Variants | Conseq Match | Conseq Mismatch | Conseq Rate | HGVSp Match | HGVSp Mismatch | HGVSp Rate | HGVSc Match | HGVSc Mismatch | HGVSc Rate |
|-------|----------|-------------|-----------------|-------------|-------------|----------------|------------|-------------|----------------|------------|
| blca_tcga_gdc | 115850 | 115651 | 7 | 99.8% | 109996 | 37 | 94.9% | 114171 | 39 | 98.6% |
| brca_tcga_gdc | 89012 | 88848 | 5 | 99.8% | 84171 | 20 | 94.6% | 87776 | 20 | 98.6% |
| chol_tcga_gdc | 3764 | 3757 | 0 | 99.8% | 3484 | 0 | 92.6% | 3695 | 0 | 98.2% |
| coad_tcga_gdc | 244552 | 244082 | 15 | 99.8% | 229196 | 91 | 93.7% | 240970 | 74 | 98.5% |
| gbm_tcga_gdc | 54870 | 54745 | 4 | 99.8% | 51908 | 14 | 94.6% | 54133 | 11 | 98.7% |
| luad_tcga_gdc | 190868 | 190445 | 6 | 99.8% | 181750 | 50 | 95.2% | 188239 | 53 | 98.6% |
| skcm_tcga_gdc | 353450 | 352553 | 8 | 99.7% | 336796 | 70 | 95.3% | 348775 | 78 | 98.7% |
| **Total** | **1052366** | **1050081** | **45** | **99.8%** | **997301** | **282** | **94.8%** | **1037759** | **275** | **98.6%** |

## Consequence Category Breakdown

| Study | delins_normalized | gene_model_change | match | mismatch | no_cds_data | position_shift | transcript_model_change | upstream_reclassified |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 1 | 0 | 114752 | 7 | 74 | 59 | 58 | 899 |
| brca_tcga_gdc | 2 | 0 | 88304 | 5 | 78 | 42 | 37 | 544 |
| chol_tcga_gdc | 0 | 0 | 3726 | 0 | 4 | 1 | 2 | 31 |
| coad_tcga_gdc | 8 | 2 | 242409 | 15 | 193 | 134 | 118 | 1673 |
| gbm_tcga_gdc | 0 | 0 | 54374 | 4 | 58 | 36 | 27 | 371 |
| luad_tcga_gdc | 1 | 0 | 189095 | 6 | 219 | 116 | 81 | 1350 |
| skcm_tcga_gdc | 5 | 2 | 349986 | 8 | 398 | 306 | 178 | 2567 |
| **Total** | **17** | **4** | **1042646** | **45** | **1024** | **694** | **501** | **7435** |

## HGVSp Category Breakdown

| Study | both_empty | delins_normalized | fuzzy_fs | maf_empty | maf_nonstandard | match | mismatch | no_cds_data | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 479 | 7 | 645 | 739 | 2620 | 109996 | 37 | 74 | 498 | 57 | 28 | 611 | 59 |
| brca_tcga_gdc | 336 | 19 | 1040 | 447 | 1567 | 84171 | 20 | 79 | 664 | 94 | 45 | 476 | 54 |
| chol_tcga_gdc | 24 | 2 | 77 | 23 | 85 | 3484 | 0 | 4 | 31 | 8 | 1 | 21 | 4 |
| coad_tcga_gdc | 2317 | 9 | 3024 | 1306 | 5448 | 229196 | 91 | 196 | 1229 | 43 | 34 | 1372 | 287 |
| gbm_tcga_gdc | 199 | 0 | 582 | 303 | 1117 | 51908 | 14 | 58 | 307 | 14 | 11 | 326 | 31 |
| luad_tcga_gdc | 728 | 10 | 1059 | 1073 | 4118 | 181750 | 50 | 223 | 675 | 53 | 21 | 1005 | 103 |
| skcm_tcga_gdc | 1414 | 17 | 493 | 2009 | 7930 | 336796 | 70 | 399 | 1330 | 52 | 21 | 2673 | 246 |
| **Total** | **5497** | **64** | **6920** | **5900** | **22885** | **997301** | **282** | **1033** | **4734** | **321** | **161** | **6484** | **784** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | mismatch | position_shift | vep_empty |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 81 | 37 | 892 | 114171 | 39 | 620 | 0 |
| brca_tcga_gdc | 8 | 106 | 41 | 538 | 87776 | 20 | 523 | 0 |
| chol_tcga_gdc | 0 | 7 | 3 | 31 | 3695 | 0 | 28 | 0 |
| coad_tcga_gdc | 25 | 90 | 248 | 1653 | 240970 | 74 | 1490 | 2 |
| gbm_tcga_gdc | 3 | 11 | 29 | 369 | 54133 | 11 | 314 | 0 |
| luad_tcga_gdc | 25 | 208 | 73 | 1329 | 188239 | 53 | 941 | 0 |
| skcm_tcga_gdc | 47 | 155 | 38 | 2525 | 348775 | 78 | 1829 | 3 |
| **Total** | **118** | **658** | **469** | **7337** | **1037759** | **275** | **5745** | **5** |

## Cancer Gene Mismatches

| Gene | Variants | Conseq Mismatches | HGVSp Mismatches | HGVSc Mismatches |
|------|----------|-------------------|------------------|------------------|
| ASMTL | 50 | 1 | 0 | 0 |
| ATM | 295 | 0 | 1 | 0 |
| BRD2 | 73 | 0 | 1 | 0 |
| CREBBP | 247 | 0 | 1 | 0 |
| CTCF | 87 | 1 | 0 | 0 |
| ELF3 | 103 | 0 | 1 | 0 |
| PPFIBP1 | 66 | 0 | 0 | 1 |
| RB1 | 235 | 0 | 1 | 0 |
| SMYD3 | 53 | 0 | 1 | 0 |

## Performance

Cache load time: 42.605s

| Study | Variants | Sequential | Seq v/s | Parallel | Par v/s | Speedup |
|-------|----------|-----------|---------|----------|---------|--------|
| blca_tcga_gdc | 115850 | 31.503s | 3677 | 15.691s | 7383 | 2.01x |
| brca_tcga_gdc | 89012 | 28.326s | 3142 | 15.929s | 5588 | 1.78x |
| chol_tcga_gdc | 3764 | 1.238s | 3041 | 607ms | 6197 | 2.04x |
| coad_tcga_gdc | 244552 | 1m0.498s | 4042 | 37.507s | 6520 | 1.61x |
| gbm_tcga_gdc | 54870 | 15.203s | 3609 | 6.461s | 8493 | 2.35x |
| luad_tcga_gdc | 190868 | 52.478s | 3637 | 24.685s | 7732 | 2.13x |
| skcm_tcga_gdc | 353450 | 1m27.789s | 4026 | 47.924s | 7375 | 1.83x |
| **Total** | **1052366** | **4m37.035s** | **3799** | **2m28.804s** | **7072** | **1.86x** |
