# TCGA Validation Report

Generated: 2026-02-26 13:26 UTC  
GENCODE transcripts loaded: 254070

## Consequence Match

| Study | Variants | Match | Mismatch | Match Rate |
|-------|----------|-------|----------|------------|
| blca_tcga_gdc | 116684 | 116518 | 166 | 99.9% |
| brca_tcga_gdc | 89012 | 88875 | 137 | 99.8% |
| chol_tcga_gdc | 3764 | 3760 | 4 | 99.9% |
| coad_tcga_gdc | 244552 | 244133 | 419 | 99.8% |
| gbm_tcga_gdc | 54870 | 54786 | 84 | 99.8% |
| luad_tcga_gdc | 190868 | 190582 | 286 | 99.9% |
| skcm_tcga_gdc | 353450 | 352783 | 667 | 99.8% |
| **Total** | **1053200** | **1051437** | **1763** | **99.8%** |

## Consequence Category Breakdown

| Study | match | mismatch | upstream_reclassified |
|-------|------|------|------|
| blca_tcga_gdc | 115613 | 166 | 905 |
| brca_tcga_gdc | 88331 | 137 | 544 |
| chol_tcga_gdc | 3730 | 4 | 30 |
| coad_tcga_gdc | 242461 | 419 | 1672 |
| gbm_tcga_gdc | 54415 | 84 | 371 |
| luad_tcga_gdc | 189232 | 286 | 1350 |
| skcm_tcga_gdc | 350217 | 667 | 2566 |
| **Total** | **1043999** | **1763** | **7438** |

## HGVSp Category Breakdown

| Study | both_empty | fuzzy_fs | maf_empty | maf_nonstandard | match | mismatch | position_shift | splice_no_protein | splice_vs_syn | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 481 | 606 | 743 | 2637 | 110671 | 162 | 584 | 57 | 616 | 127 |
| brca_tcga_gdc | 335 | 993 | 448 | 1567 | 83989 | 149 | 832 | 94 | 476 | 129 |
| chol_tcga_gdc | 24 | 73 | 23 | 85 | 3466 | 4 | 54 | 8 | 21 | 6 |
| coad_tcga_gdc | 2308 | 2902 | 1315 | 5448 | 228365 | 275 | 2048 | 44 | 1373 | 474 |
| gbm_tcga_gdc | 200 | 565 | 302 | 1117 | 51859 | 80 | 320 | 14 | 326 | 87 |
| luad_tcga_gdc | 723 | 999 | 1078 | 4118 | 181568 | 254 | 751 | 53 | 1005 | 319 |
| skcm_tcga_gdc | 1413 | 454 | 2010 | 7931 | 336747 | 427 | 1109 | 52 | 2673 | 634 |
| **Total** | **5484** | **6592** | **5919** | **22903** | **996665** | **1351** | **5698** | **322** | **6490** | **1776** |

## HGVSc Category Breakdown

| Study | both_empty | maf_empty | match | mismatch | position_shift | vep_empty |
|-------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 898 | 114990 | 160 | 626 | 0 |
| brca_tcga_gdc | 8 | 538 | 87773 | 167 | 526 | 0 |
| chol_tcga_gdc | 0 | 31 | 3695 | 10 | 28 | 0 |
| coad_tcga_gdc | 25 | 1653 | 240975 | 410 | 1487 | 2 |
| gbm_tcga_gdc | 3 | 369 | 54135 | 51 | 312 | 0 |
| luad_tcga_gdc | 25 | 1329 | 188265 | 336 | 913 | 0 |
| skcm_tcga_gdc | 47 | 2525 | 348823 | 270 | 1782 | 3 |
| **Total** | **118** | **7343** | **1038656** | **1404** | **5674** | **5** |

## Performance

| Study | Variants | Time | Variants/sec |
|-------|----------|------|-------------|
| blca_tcga_gdc | 116684 | 27.064s | 4311 |
| brca_tcga_gdc | 89012 | 25.47s | 3495 |
| chol_tcga_gdc | 3764 | 901ms | 4177 |
| coad_tcga_gdc | 244552 | 48.179s | 5076 |
| gbm_tcga_gdc | 54870 | 13.09s | 4192 |
| luad_tcga_gdc | 190868 | 42.279s | 4514 |
| skcm_tcga_gdc | 353450 | 1m9.323s | 5099 |
| **Total** | **1053200** | **3m46.305s** | **4654** |
