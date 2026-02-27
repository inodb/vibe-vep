# TCGA Validation Report

Generated: 2026-02-27 05:43 UTC  
GENCODE transcripts loaded: 254070

## Consequence Match

| Study | Variants | Match | Mismatch | Match Rate |
|-------|----------|-------|----------|------------|
| blca_tcga_gdc | 116684 | 116473 | 137 | 99.8% |
| brca_tcga_gdc | 89012 | 88830 | 104 | 99.8% |
| chol_tcga_gdc | 3764 | 3758 | 2 | 99.8% |
| coad_tcga_gdc | 244552 | 244021 | 338 | 99.8% |
| gbm_tcga_gdc | 54870 | 54741 | 71 | 99.8% |
| luad_tcga_gdc | 190868 | 190435 | 214 | 99.8% |
| skcm_tcga_gdc | 353450 | 352547 | 505 | 99.7% |
| **Total** | **1053200** | **1050805** | **1371** | **99.8%** |

## Consequence Category Breakdown

| Study | match | mismatch | no_cds_data | upstream_reclassified |
|-------|------|------|------|------|
| blca_tcga_gdc | 115568 | 137 | 74 | 905 |
| brca_tcga_gdc | 88286 | 104 | 78 | 544 |
| chol_tcga_gdc | 3728 | 2 | 4 | 30 |
| coad_tcga_gdc | 242349 | 338 | 193 | 1672 |
| gbm_tcga_gdc | 54370 | 71 | 58 | 371 |
| luad_tcga_gdc | 189085 | 214 | 219 | 1350 |
| skcm_tcga_gdc | 349981 | 505 | 398 | 2566 |
| **Total** | **1043367** | **1371** | **1024** | **7438** |

## HGVSp Category Breakdown

| Study | both_empty | fuzzy_fs | maf_empty | maf_nonstandard | match | mismatch | no_cds_data | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 481 | 574 | 743 | 2637 | 110757 | 125 | 74 | 539 | 57 | 28 | 616 | 53 |
| brca_tcga_gdc | 335 | 939 | 448 | 1567 | 84158 | 96 | 79 | 725 | 94 | 45 | 476 | 50 |
| chol_tcga_gdc | 24 | 69 | 23 | 85 | 3478 | 4 | 4 | 45 | 8 | 1 | 21 | 2 |
| coad_tcga_gdc | 2308 | 2785 | 1315 | 5448 | 229073 | 234 | 196 | 1465 | 44 | 33 | 1373 | 278 |
| gbm_tcga_gdc | 200 | 547 | 302 | 1117 | 51893 | 67 | 58 | 306 | 14 | 11 | 326 | 29 |
| luad_tcga_gdc | 723 | 962 | 1078 | 4118 | 181715 | 230 | 223 | 644 | 53 | 21 | 1005 | 96 |
| skcm_tcga_gdc | 1413 | 426 | 2010 | 7931 | 336744 | 409 | 399 | 1137 | 52 | 21 | 2673 | 235 |
| **Total** | **5484** | **6302** | **5919** | **22903** | **997818** | **1165** | **1033** | **4861** | **322** | **160** | **6490** | **743** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | mismatch | position_shift | vep_empty |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 82 | 38 | 898 | 114987 | 40 | 629 | 0 |
| brca_tcga_gdc | 8 | 106 | 41 | 538 | 87767 | 20 | 532 | 0 |
| chol_tcga_gdc | 0 | 7 | 3 | 31 | 3695 | 0 | 28 | 0 |
| coad_tcga_gdc | 25 | 90 | 247 | 1653 | 240972 | 73 | 1490 | 2 |
| gbm_tcga_gdc | 3 | 11 | 29 | 369 | 54129 | 11 | 318 | 0 |
| luad_tcga_gdc | 25 | 208 | 73 | 1329 | 188236 | 55 | 942 | 0 |
| skcm_tcga_gdc | 47 | 155 | 38 | 2525 | 348745 | 77 | 1860 | 3 |
| **Total** | **118** | **659** | **469** | **7343** | **1038531** | **276** | **5799** | **5** |

## Performance

| Study | Variants | Time | Variants/sec |
|-------|----------|------|-------------|
| blca_tcga_gdc | 116684 | 19.544s | 5970 |
| brca_tcga_gdc | 89012 | 17.441s | 5103 |
| chol_tcga_gdc | 3764 | 736ms | 5114 |
| coad_tcga_gdc | 244552 | 36.606s | 6681 |
| gbm_tcga_gdc | 54870 | 10.695s | 5130 |
| luad_tcga_gdc | 190868 | 31.3s | 6098 |
| skcm_tcga_gdc | 353450 | 47.883s | 7382 |
| **Total** | **1053200** | **2m44.206s** | **6414** |
