# TCGA Validation Report

Generated: 2026-02-27 12:58 UTC  
GENCODE transcripts loaded: 254070

## Consequence Match

| Study | Variants | Match | Mismatch | Match Rate |
|-------|----------|-------|----------|------------|
| blca_tcga_gdc | 116684 | 116474 | 136 | 99.8% |
| brca_tcga_gdc | 89012 | 88830 | 104 | 99.8% |
| chol_tcga_gdc | 3764 | 3758 | 2 | 99.8% |
| coad_tcga_gdc | 244552 | 244017 | 342 | 99.8% |
| gbm_tcga_gdc | 54870 | 54741 | 71 | 99.8% |
| luad_tcga_gdc | 190868 | 190433 | 216 | 99.8% |
| skcm_tcga_gdc | 353450 | 352547 | 505 | 99.7% |
| **Total** | **1053200** | **1050800** | **1376** | **99.8%** |

## Consequence Category Breakdown

| Study | match | mismatch | no_cds_data | upstream_reclassified |
|-------|------|------|------|------|
| blca_tcga_gdc | 115569 | 136 | 74 | 905 |
| brca_tcga_gdc | 88286 | 104 | 78 | 544 |
| chol_tcga_gdc | 3728 | 2 | 4 | 30 |
| coad_tcga_gdc | 242345 | 342 | 193 | 1672 |
| gbm_tcga_gdc | 54370 | 71 | 58 | 371 |
| luad_tcga_gdc | 189083 | 216 | 219 | 1350 |
| skcm_tcga_gdc | 349981 | 505 | 398 | 2566 |
| **Total** | **1043362** | **1376** | **1024** | **7438** |

## HGVSp Category Breakdown

| Study | both_empty | fuzzy_fs | maf_empty | maf_nonstandard | match | mismatch | no_cds_data | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 481 | 574 | 743 | 2637 | 110790 | 125 | 74 | 506 | 57 | 28 | 616 | 53 |
| brca_tcga_gdc | 335 | 939 | 448 | 1567 | 84171 | 91 | 79 | 717 | 94 | 45 | 476 | 50 |
| chol_tcga_gdc | 24 | 69 | 23 | 85 | 3484 | 4 | 4 | 39 | 8 | 1 | 21 | 2 |
| coad_tcga_gdc | 2308 | 2786 | 1315 | 5448 | 229185 | 231 | 196 | 1355 | 44 | 33 | 1373 | 278 |
| gbm_tcga_gdc | 200 | 547 | 302 | 1117 | 51902 | 65 | 58 | 299 | 14 | 11 | 326 | 29 |
| luad_tcga_gdc | 723 | 962 | 1078 | 4118 | 181750 | 225 | 223 | 614 | 53 | 21 | 1005 | 96 |
| skcm_tcga_gdc | 1413 | 426 | 2010 | 7931 | 336790 | 414 | 399 | 1086 | 52 | 21 | 2673 | 235 |
| **Total** | **5484** | **6303** | **5919** | **22903** | **998072** | **1155** | **1033** | **4616** | **322** | **160** | **6490** | **743** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | mismatch | position_shift | vep_empty |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 82 | 38 | 898 | 114990 | 40 | 626 | 0 |
| brca_tcga_gdc | 8 | 106 | 41 | 538 | 87775 | 20 | 524 | 0 |
| chol_tcga_gdc | 0 | 7 | 3 | 31 | 3695 | 0 | 28 | 0 |
| coad_tcga_gdc | 25 | 90 | 247 | 1653 | 240960 | 73 | 1502 | 2 |
| gbm_tcga_gdc | 3 | 11 | 29 | 369 | 54127 | 11 | 320 | 0 |
| luad_tcga_gdc | 25 | 208 | 73 | 1329 | 188238 | 55 | 940 | 0 |
| skcm_tcga_gdc | 47 | 155 | 38 | 2525 | 348768 | 77 | 1837 | 3 |
| **Total** | **118** | **659** | **469** | **7343** | **1038553** | **276** | **5777** | **5** |

## Performance

| Study | Variants | Time | Variants/sec |
|-------|----------|------|-------------|
| blca_tcga_gdc | 116684 | 24.997s | 4668 |
| brca_tcga_gdc | 89012 | 22.116s | 4025 |
| chol_tcga_gdc | 3764 | 1.06s | 3551 |
| coad_tcga_gdc | 244552 | 45.786s | 5341 |
| gbm_tcga_gdc | 54870 | 12.621s | 4348 |
| luad_tcga_gdc | 190868 | 40.967s | 4659 |
| skcm_tcga_gdc | 353450 | 1m2.411s | 5663 |
| **Total** | **1053200** | **3m29.957s** | **5016** |
