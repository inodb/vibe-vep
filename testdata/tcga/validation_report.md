# TCGA Validation Report

Generated: 2026-02-27 16:37 UTC  
GENCODE transcripts loaded: 254070

## Consequence Match

| Study | Variants | Match | Mismatch | Match Rate |
|-------|----------|-------|----------|------------|
| blca_tcga_gdc | 115850 | 115647 | 41 | 99.8% |
| brca_tcga_gdc | 89012 | 88834 | 34 | 99.8% |
| chol_tcga_gdc | 3764 | 3757 | 0 | 99.8% |
| coad_tcga_gdc | 244552 | 244070 | 93 | 99.8% |
| gbm_tcga_gdc | 54870 | 54745 | 14 | 99.8% |
| luad_tcga_gdc | 190868 | 190439 | 57 | 99.8% |
| skcm_tcga_gdc | 353450 | 352549 | 85 | 99.7% |
| **Total** | **1052366** | **1050041** | **324** | **99.8%** |

## Consequence Category Breakdown

| Study | match | mismatch | no_cds_data | position_shift | upstream_reclassified |
|-------|------|------|------|------|------|
| blca_tcga_gdc | 114748 | 41 | 86 | 76 | 899 |
| brca_tcga_gdc | 88289 | 34 | 87 | 57 | 545 |
| chol_tcga_gdc | 3726 | 0 | 6 | 1 | 31 |
| coad_tcga_gdc | 242397 | 93 | 221 | 168 | 1673 |
| gbm_tcga_gdc | 54374 | 14 | 64 | 47 | 371 |
| luad_tcga_gdc | 189089 | 57 | 235 | 137 | 1350 |
| skcm_tcga_gdc | 349982 | 85 | 454 | 362 | 2567 |
| **Total** | **1042605** | **324** | **1153** | **848** | **7436** |

## HGVSp Category Breakdown

| Study | both_empty | delins_normalized | fuzzy_fs | maf_empty | maf_nonstandard | match | mismatch | no_cds_data | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 479 | 7 | 645 | 739 | 2620 | 109990 | 37 | 74 | 504 | 57 | 28 | 611 | 59 |
| brca_tcga_gdc | 336 | 19 | 1040 | 447 | 1567 | 84170 | 20 | 79 | 665 | 94 | 45 | 476 | 54 |
| chol_tcga_gdc | 24 | 2 | 77 | 23 | 85 | 3484 | 0 | 4 | 31 | 8 | 1 | 21 | 4 |
| coad_tcga_gdc | 2317 | 9 | 3025 | 1306 | 5448 | 229193 | 91 | 196 | 1231 | 43 | 34 | 1372 | 287 |
| gbm_tcga_gdc | 199 | 0 | 582 | 303 | 1117 | 51912 | 14 | 58 | 303 | 14 | 11 | 326 | 31 |
| luad_tcga_gdc | 728 | 10 | 1060 | 1073 | 4118 | 181749 | 50 | 223 | 675 | 53 | 21 | 1005 | 103 |
| skcm_tcga_gdc | 1414 | 17 | 493 | 2009 | 7930 | 336792 | 70 | 399 | 1334 | 52 | 21 | 2673 | 246 |
| **Total** | **5497** | **64** | **6922** | **5900** | **22885** | **997290** | **282** | **1033** | **4743** | **321** | **161** | **6484** | **784** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | mismatch | position_shift | vep_empty |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 81 | 37 | 892 | 114168 | 40 | 622 | 0 |
| brca_tcga_gdc | 8 | 106 | 41 | 538 | 87774 | 20 | 525 | 0 |
| chol_tcga_gdc | 0 | 7 | 3 | 31 | 3695 | 0 | 28 | 0 |
| coad_tcga_gdc | 25 | 90 | 248 | 1653 | 240967 | 74 | 1493 | 2 |
| gbm_tcga_gdc | 3 | 11 | 29 | 369 | 54137 | 11 | 310 | 0 |
| luad_tcga_gdc | 25 | 208 | 73 | 1329 | 188238 | 55 | 940 | 0 |
| skcm_tcga_gdc | 47 | 155 | 38 | 2525 | 348771 | 78 | 1833 | 3 |
| **Total** | **118** | **658** | **469** | **7337** | **1037750** | **278** | **5751** | **5** |

## Performance

| Study | Variants | Time | Variants/sec |
|-------|----------|------|-------------|
| blca_tcga_gdc | 115850 | 41.821s | 2770 |
| brca_tcga_gdc | 89012 | 31.208s | 2852 |
| chol_tcga_gdc | 3764 | 1.183s | 3182 |
| coad_tcga_gdc | 244552 | 57.093s | 4283 |
| gbm_tcga_gdc | 54870 | 16.524s | 3321 |
| luad_tcga_gdc | 190868 | 55.276s | 3453 |
| skcm_tcga_gdc | 353450 | 1m16.523s | 4619 |
| **Total** | **1052366** | **4m39.628s** | **3763** |
