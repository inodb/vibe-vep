# TCGA Validation Report

Generated: 2026-02-26 18:25 UTC  
GENCODE transcripts loaded: 254070

## Consequence Match

| Study | Variants | Match | Mismatch | Match Rate |
|-------|----------|-------|----------|------------|
| blca_tcga_gdc | 116684 | 116471 | 139 | 99.8% |
| brca_tcga_gdc | 89012 | 88825 | 109 | 99.8% |
| chol_tcga_gdc | 3764 | 3758 | 2 | 99.8% |
| coad_tcga_gdc | 244552 | 244018 | 341 | 99.8% |
| gbm_tcga_gdc | 54870 | 54743 | 69 | 99.8% |
| luad_tcga_gdc | 190868 | 190432 | 217 | 99.8% |
| skcm_tcga_gdc | 353450 | 352547 | 505 | 99.7% |
| **Total** | **1053200** | **1050794** | **1382** | **99.8%** |

## Consequence Category Breakdown

| Study | match | mismatch | no_cds_data | upstream_reclassified |
|-------|------|------|------|------|
| blca_tcga_gdc | 115566 | 139 | 74 | 905 |
| brca_tcga_gdc | 88281 | 109 | 78 | 544 |
| chol_tcga_gdc | 3728 | 2 | 4 | 30 |
| coad_tcga_gdc | 242346 | 341 | 193 | 1672 |
| gbm_tcga_gdc | 54372 | 69 | 58 | 371 |
| luad_tcga_gdc | 189082 | 217 | 219 | 1350 |
| skcm_tcga_gdc | 349981 | 505 | 398 | 2566 |
| **Total** | **1043356** | **1382** | **1024** | **7438** |

## HGVSp Category Breakdown

| Study | both_empty | fuzzy_fs | maf_empty | maf_nonstandard | match | mismatch | no_cds_data | position_shift | splice_no_protein | splice_vs_syn | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 481 | 605 | 743 | 2637 | 110685 | 157 | 74 | 576 | 57 | 616 | 53 |
| brca_tcga_gdc | 335 | 993 | 448 | 1567 | 84007 | 148 | 79 | 815 | 94 | 476 | 50 |
| chol_tcga_gdc | 24 | 73 | 23 | 85 | 3466 | 4 | 4 | 54 | 8 | 21 | 2 |
| coad_tcga_gdc | 2308 | 2894 | 1315 | 5448 | 228400 | 267 | 196 | 2029 | 44 | 1373 | 278 |
| gbm_tcga_gdc | 200 | 564 | 302 | 1117 | 51878 | 78 | 58 | 304 | 14 | 326 | 29 |
| luad_tcga_gdc | 723 | 1000 | 1078 | 4118 | 181551 | 260 | 223 | 761 | 53 | 1005 | 96 |
| skcm_tcga_gdc | 1413 | 454 | 2010 | 7931 | 336791 | 429 | 399 | 1063 | 52 | 2673 | 235 |
| **Total** | **5484** | **6583** | **5919** | **22903** | **996778** | **1343** | **1033** | **5602** | **322** | **6490** | **743** |

## HGVSc Category Breakdown

| Study | both_empty | dup_vs_ins | maf_empty | match | mismatch | position_shift | vep_empty |
|-------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 38 | 898 | 115008 | 122 | 608 | 0 |
| brca_tcga_gdc | 8 | 41 | 538 | 87793 | 126 | 506 | 0 |
| chol_tcga_gdc | 0 | 3 | 31 | 3695 | 7 | 28 | 0 |
| coad_tcga_gdc | 25 | 247 | 1653 | 241008 | 163 | 1454 | 2 |
| gbm_tcga_gdc | 3 | 29 | 369 | 54154 | 22 | 293 | 0 |
| luad_tcga_gdc | 25 | 73 | 1329 | 188249 | 263 | 929 | 0 |
| skcm_tcga_gdc | 47 | 38 | 2525 | 348870 | 232 | 1735 | 3 |
| **Total** | **118** | **469** | **7343** | **1038777** | **935** | **5553** | **5** |

## Performance

| Study | Variants | Time | Variants/sec |
|-------|----------|------|-------------|
| blca_tcga_gdc | 116684 | 44.704s | 2610 |
| brca_tcga_gdc | 89012 | 34.587s | 2574 |
| chol_tcga_gdc | 3764 | 1.803s | 2088 |
| coad_tcga_gdc | 244552 | 1m11.681s | 3412 |
| gbm_tcga_gdc | 54870 | 15.887s | 3454 |
| luad_tcga_gdc | 190868 | 54.363s | 3511 |
| skcm_tcga_gdc | 353450 | 1m24.44s | 4186 |
| **Total** | **1053200** | **5m7.466s** | **3425** |
