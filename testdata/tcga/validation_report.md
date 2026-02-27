# TCGA Validation Report

Generated: 2026-02-27 17:39 UTC  
GENCODE transcripts loaded: 254070

## Match Rates

| Study | Variants | Conseq Match | Conseq Mismatch | Conseq Rate | HGVSp Match | HGVSp Mismatch | HGVSp Rate | HGVSc Match | HGVSc Mismatch | HGVSc Rate |
|-------|----------|-------------|-----------------|-------------|-------------|----------------|------------|-------------|----------------|------------|
| blca_tcga_gdc | 115850 | 115650 | 7 | 99.8% | 109981 | 37 | 94.9% | 114156 | 40 | 98.5% |
| brca_tcga_gdc | 89012 | 88848 | 5 | 99.8% | 84176 | 20 | 94.6% | 87781 | 20 | 98.6% |
| chol_tcga_gdc | 3764 | 3757 | 0 | 99.8% | 3483 | 0 | 92.5% | 3694 | 0 | 98.1% |
| coad_tcga_gdc | 244552 | 244082 | 15 | 99.8% | 229171 | 91 | 93.7% | 240945 | 74 | 98.5% |
| gbm_tcga_gdc | 54870 | 54745 | 4 | 99.8% | 51901 | 14 | 94.6% | 54126 | 11 | 98.6% |
| luad_tcga_gdc | 190868 | 190444 | 7 | 99.8% | 181762 | 50 | 95.2% | 188251 | 55 | 98.6% |
| skcm_tcga_gdc | 353450 | 352553 | 8 | 99.7% | 336818 | 70 | 95.3% | 348795 | 78 | 98.7% |
| **Total** | **1052366** | **1050079** | **46** | **99.8%** | **997292** | **282** | **94.8%** | **1037748** | **278** | **98.6%** |

## Consequence Category Breakdown

| Study | delins_normalized | gene_model_change | match | mismatch | no_cds_data | position_shift | transcript_model_change | upstream_reclassified |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 1 | 0 | 114751 | 7 | 74 | 60 | 58 | 899 |
| brca_tcga_gdc | 2 | 0 | 88304 | 5 | 78 | 42 | 37 | 544 |
| chol_tcga_gdc | 0 | 0 | 3726 | 0 | 4 | 1 | 2 | 31 |
| coad_tcga_gdc | 8 | 2 | 242409 | 15 | 193 | 134 | 118 | 1673 |
| gbm_tcga_gdc | 0 | 0 | 54374 | 4 | 58 | 36 | 27 | 371 |
| luad_tcga_gdc | 1 | 0 | 189094 | 7 | 219 | 116 | 81 | 1350 |
| skcm_tcga_gdc | 5 | 2 | 349986 | 8 | 398 | 306 | 178 | 2567 |
| **Total** | **17** | **4** | **1042644** | **46** | **1024** | **695** | **501** | **7435** |

## HGVSp Category Breakdown

| Study | both_empty | delins_normalized | fuzzy_fs | maf_empty | maf_nonstandard | match | mismatch | no_cds_data | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 479 | 7 | 646 | 739 | 2620 | 109981 | 37 | 74 | 512 | 57 | 28 | 611 | 59 |
| brca_tcga_gdc | 336 | 19 | 1040 | 447 | 1567 | 84176 | 20 | 79 | 659 | 94 | 45 | 476 | 54 |
| chol_tcga_gdc | 24 | 2 | 77 | 23 | 85 | 3483 | 0 | 4 | 32 | 8 | 1 | 21 | 4 |
| coad_tcga_gdc | 2317 | 9 | 3032 | 1306 | 5448 | 229171 | 91 | 196 | 1246 | 43 | 34 | 1372 | 287 |
| gbm_tcga_gdc | 199 | 0 | 583 | 303 | 1117 | 51901 | 14 | 58 | 313 | 14 | 11 | 326 | 31 |
| luad_tcga_gdc | 728 | 10 | 1059 | 1073 | 4118 | 181762 | 50 | 223 | 663 | 53 | 21 | 1005 | 103 |
| skcm_tcga_gdc | 1414 | 17 | 493 | 2009 | 7930 | 336818 | 70 | 399 | 1308 | 52 | 21 | 2673 | 246 |
| **Total** | **5497** | **64** | **6930** | **5900** | **22885** | **997292** | **282** | **1033** | **4733** | **321** | **161** | **6484** | **784** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | mismatch | position_shift | vep_empty |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 81 | 37 | 892 | 114156 | 40 | 634 | 0 |
| brca_tcga_gdc | 8 | 106 | 41 | 538 | 87781 | 20 | 518 | 0 |
| chol_tcga_gdc | 0 | 7 | 3 | 31 | 3694 | 0 | 29 | 0 |
| coad_tcga_gdc | 25 | 90 | 248 | 1653 | 240945 | 74 | 1515 | 2 |
| gbm_tcga_gdc | 3 | 11 | 29 | 369 | 54126 | 11 | 321 | 0 |
| luad_tcga_gdc | 25 | 208 | 73 | 1329 | 188251 | 55 | 927 | 0 |
| skcm_tcga_gdc | 47 | 155 | 38 | 2525 | 348795 | 78 | 1809 | 3 |
| **Total** | **118** | **658** | **469** | **7337** | **1037748** | **278** | **5753** | **5** |

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

| Study | Variants | Time | Variants/sec |
|-------|----------|------|-------------|
| blca_tcga_gdc | 115850 | 32.33s | 3583 |
| brca_tcga_gdc | 89012 | 23.729s | 3751 |
| chol_tcga_gdc | 3764 | 2.888s | 1303 |
| coad_tcga_gdc | 244552 | 50.272s | 4865 |
| gbm_tcga_gdc | 54870 | 22.314s | 2459 |
| luad_tcga_gdc | 190868 | 53.964s | 3537 |
| skcm_tcga_gdc | 353450 | 1m15.055s | 4709 |
| **Total** | **1052366** | **4m20.553s** | **4039** |
