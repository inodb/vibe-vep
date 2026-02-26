# TCGA Validation Report

Generated: 2026-02-26 04:07 UTC  
GENCODE transcripts loaded: 254070

## Consequence Match

| Study | Variants | Match | Mismatch | Match Rate |
|-------|----------|-------|----------|------------|
| blca_tcga_gdc | 116684 | 116516 | 168 | 99.9% |
| brca_tcga_gdc | 89012 | 88874 | 138 | 99.8% |
| chol_tcga_gdc | 3764 | 3760 | 4 | 99.9% |
| coad_tcga_gdc | 244552 | 244136 | 416 | 99.8% |
| gbm_tcga_gdc | 54870 | 54786 | 84 | 99.8% |
| luad_tcga_gdc | 190868 | 190582 | 286 | 99.9% |
| skcm_tcga_gdc | 353450 | 352782 | 668 | 99.8% |
| **Total** | **1053200** | **1051436** | **1764** | **99.8%** |

## HGVSp Match

| Study | Compared | Match | Mismatch | Match Rate |
|-------|----------|-------|----------|------------|
| blca_tcga_gdc | 113507 | 111277 | 2230 | 98.0% |
| brca_tcga_gdc | 87015 | 84984 | 2031 | 97.7% |
| chol_tcga_gdc | 3647 | 3539 | 108 | 97.0% |
| coad_tcga_gdc | 236750 | 231237 | 5513 | 97.7% |
| gbm_tcga_gdc | 53539 | 52418 | 1121 | 97.9% |
| luad_tcga_gdc | 185974 | 182542 | 3432 | 98.2% |
| skcm_tcga_gdc | 344052 | 337119 | 6933 | 98.0% |
| **Total** | **1024484** | **1003116** | **21368** | **97.9%** |

### HGVSp Match Breakdown

| Study | Exact | Fuzzy (fs) | Skipped: empty | Skipped: p.*N* | Skipped: splice |
|-------|-------|------------|----------------|----------------|------------------|
| blca_tcga_gdc | 110672 | 605 | 483 | 2637 | 57 |
| brca_tcga_gdc | 83991 | 993 | 336 | 1567 | 94 |
| chol_tcga_gdc | 3466 | 73 | 24 | 85 | 8 |
| coad_tcga_gdc | 228342 | 2895 | 2310 | 5448 | 44 |
| gbm_tcga_gdc | 51854 | 564 | 200 | 1117 | 14 |
| luad_tcga_gdc | 181542 | 1000 | 723 | 4118 | 53 |
| skcm_tcga_gdc | 336665 | 454 | 1415 | 7931 | 52 |
| **Total** | **996532** | **6584** | **5491** | **22903** | **322** |

## HGVSc Match

| Study | Variants | Match | Mismatch | Match Rate |
|-------|----------|-------|----------|------------|
| blca_tcga_gdc | 116590 | 113899 | 2691 | 97.7% |
| brca_tcga_gdc | 88953 | 87048 | 1905 | 97.9% |
| chol_tcga_gdc | 3761 | 3654 | 107 | 97.2% |
| coad_tcga_gdc | 244364 | 238864 | 5500 | 97.7% |
| gbm_tcga_gdc | 54822 | 53579 | 1243 | 97.7% |
| luad_tcga_gdc | 190667 | 186224 | 4443 | 97.7% |
| skcm_tcga_gdc | 353101 | 344891 | 8210 | 97.7% |
| **Total** | **1052258** | **1028159** | **24099** | **97.7%** |

## Performance

| Study | Variants | Time | Variants/sec |
|-------|----------|------|-------------|
| blca_tcga_gdc | 116684 | 19.875s | 5871 |
| brca_tcga_gdc | 89012 | 15.621s | 5698 |
| chol_tcga_gdc | 3764 | 797ms | 4723 |
| coad_tcga_gdc | 244552 | 35.573s | 6875 |
| gbm_tcga_gdc | 54870 | 9.994s | 5490 |
| luad_tcga_gdc | 190868 | 31.368s | 6085 |
| skcm_tcga_gdc | 353450 | 47.923s | 7375 |
| **Total** | **1053200** | **2m41.152s** | **6535** |
