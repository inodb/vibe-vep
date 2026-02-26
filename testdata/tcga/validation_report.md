# TCGA Validation Report

Generated: 2026-02-26 12:49 UTC  
GENCODE transcripts loaded: 254070

## Consequence Match

| Study | Variants | Match | Mismatch | Match Rate |
|-------|----------|-------|----------|------------|
| blca_tcga_gdc | 116684 | 116513 | 171 | 99.9% |
| brca_tcga_gdc | 89012 | 88874 | 138 | 99.8% |
| chol_tcga_gdc | 3764 | 3760 | 4 | 99.9% |
| coad_tcga_gdc | 244552 | 244131 | 421 | 99.8% |
| gbm_tcga_gdc | 54870 | 54785 | 85 | 99.8% |
| luad_tcga_gdc | 190868 | 190581 | 287 | 99.8% |
| skcm_tcga_gdc | 353450 | 352780 | 670 | 99.8% |
| **Total** | **1053200** | **1051424** | **1776** | **99.8%** |

## Consequence Category Breakdown

| Study | match | mismatch | upstream_reclassified |
|-------|------|------|------|
| blca_tcga_gdc | 115608 | 171 | 905 |
| brca_tcga_gdc | 88330 | 138 | 544 |
| chol_tcga_gdc | 3730 | 4 | 30 |
| coad_tcga_gdc | 242459 | 421 | 1672 |
| gbm_tcga_gdc | 54414 | 85 | 371 |
| luad_tcga_gdc | 189231 | 287 | 1350 |
| skcm_tcga_gdc | 350214 | 670 | 2566 |
| **Total** | **1043986** | **1776** | **7438** |

## HGVSp Category Breakdown

| Study | both_empty | fuzzy_fs | maf_empty | maf_nonstandard | match | mismatch | position_shift | splice_no_protein | splice_vs_syn | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 483 | 605 | 741 | 2637 | 110668 | 158 | 591 | 57 | 616 | 128 |
| brca_tcga_gdc | 336 | 993 | 447 | 1567 | 83996 | 149 | 823 | 94 | 476 | 131 |
| chol_tcga_gdc | 24 | 73 | 23 | 85 | 3466 | 4 | 53 | 8 | 21 | 7 |
| coad_tcga_gdc | 2308 | 2894 | 1315 | 5448 | 228363 | 267 | 2065 | 44 | 1373 | 475 |
| gbm_tcga_gdc | 200 | 564 | 302 | 1117 | 51849 | 79 | 332 | 14 | 326 | 87 |
| luad_tcga_gdc | 725 | 999 | 1076 | 4118 | 181493 | 258 | 820 | 54 | 1005 | 320 |
| skcm_tcga_gdc | 1415 | 454 | 2008 | 7931 | 336671 | 429 | 1180 | 52 | 2673 | 637 |
| **Total** | **5491** | **6582** | **5912** | **22903** | **996506** | **1344** | **5864** | **323** | **6490** | **1785** |

## HGVSc Category Breakdown

| Study | both_empty | maf_empty | match | mismatch | position_shift | vep_empty |
|-------|------|------|------|------|------|------|
| blca_tcga_gdc | 94 | 814 | 113896 | 146 | 380 | 1354 |
| brca_tcga_gdc | 59 | 487 | 87054 | 157 | 406 | 849 |
| chol_tcga_gdc | 3 | 28 | 3654 | 10 | 20 | 49 |
| coad_tcga_gdc | 188 | 1490 | 238887 | 363 | 1197 | 2427 |
| gbm_tcga_gdc | 48 | 324 | 53574 | 43 | 248 | 633 |
| luad_tcga_gdc | 201 | 1153 | 186173 | 305 | 641 | 2395 |
| skcm_tcga_gdc | 349 | 2223 | 344896 | 251 | 1331 | 4400 |
| **Total** | **942** | **6519** | **1028134** | **1275** | **4223** | **12107** |

## Performance

| Study | Variants | Time | Variants/sec |
|-------|----------|------|-------------|
| blca_tcga_gdc | 116684 | 22.285s | 5236 |
| brca_tcga_gdc | 89012 | 21.689s | 4104 |
| chol_tcga_gdc | 3764 | 1.15s | 3274 |
| coad_tcga_gdc | 244552 | 36.835s | 6639 |
| gbm_tcga_gdc | 54870 | 11.274s | 4867 |
| luad_tcga_gdc | 190868 | 34.782s | 5488 |
| skcm_tcga_gdc | 353450 | 48.472s | 7292 |
| **Total** | **1053200** | **2m56.488s** | **5968** |
