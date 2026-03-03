# Annotation Sources Report

Generated: 2026-03-03 16:19 UTC  
GENCODE transcripts: 254070  
AlphaMissense variants in database: 71697556  
Cancer Hotspots: 679 genes, 4183 positions  
Workers: 4 (GOMAXPROCS)

## AlphaMissense Coverage

| Study | Variants | Missense | AM Hits | Coverage | likely_benign | ambiguous | likely_pathogenic |
|-------|----------|----------|---------|----------|---------------|-----------|-------------------|
| blca_tcga_gdc | 115850 | 75938 | 69937 | 92.1% | 41805 | 7662 | 20470 |
| brca_tcga_gdc | 89012 | 57703 | 53240 | 92.3% | 31416 | 5807 | 16017 |
| chol_tcga_gdc | 3764 | 2336 | 2151 | 92.1% | 1244 | 234 | 673 |
| coad_tcga_gdc | 244552 | 128655 | 118282 | 91.9% | 70342 | 12727 | 35213 |
| gbm_tcga_gdc | 54870 | 32684 | 30155 | 92.3% | 18286 | 3259 | 8610 |
| luad_tcga_gdc | 190868 | 120249 | 111597 | 92.8% | 65087 | 12682 | 33828 |
| skcm_tcga_gdc | 353450 | 204372 | 185786 | 90.9% | 115483 | 20014 | 50289 |
| **Total** | **1052366** | **621937** | **571148** | **91.8%** | **343663** | **62385** | **165100** |

## Cancer Hotspots Coverage

| Study | Variants | Checked | Hotspot Hits | Hit Rate | single residue | in-frame indel | 3d | splice |
|-------|----------|---------|--------------|----------|----------------|----------------|----|--------|
| blca_tcga_gdc | 115850 | 114037 | 572 | 0.50% | 301 | 1 | 248 | 0 |
| brca_tcga_gdc | 89012 | 87866 | 884 | 1.01% | 538 | 0 | 309 | 0 |
| chol_tcga_gdc | 3764 | 3695 | 41 | 1.11% | 30 | 0 | 11 | 0 |
| coad_tcga_gdc | 244552 | 239629 | 1099 | 0.46% | 608 | 1 | 468 | 0 |
| gbm_tcga_gdc | 54870 | 54144 | 397 | 0.73% | 153 | 0 | 213 | 0 |
| luad_tcga_gdc | 190868 | 187939 | 825 | 0.44% | 393 | 0 | 400 | 0 |
| skcm_tcga_gdc | 353450 | 348150 | 838 | 0.24% | 345 | 0 | 476 | 0 |
| **Total** | **1052366** | **1035460** | **4656** | **0.45%** | **2368** | **2** | **2125** | **0** |

## AlphaMissense Performance

| Study | Variants | Base Time | AM Lookup Time | AM Overhead | Lookups/sec |
|-------|----------|-----------|----------------|-------------|-------------|
| blca_tcga_gdc | 115850 | 1m11.254s | 9.158s | 12.9% | 8292 |
| brca_tcga_gdc | 89012 | 23.228s | 7.19s | 31.0% | 8025 |
| chol_tcga_gdc | 3764 | 1.06s | 4.012s | 378.6% | 582 |
| coad_tcga_gdc | 244552 | 42.547s | 10.916s | 25.7% | 11786 |
| gbm_tcga_gdc | 54870 | 11.985s | 5.924s | 49.4% | 5518 |
| luad_tcga_gdc | 190868 | 41.836s | 10.994s | 26.3% | 10938 |
| skcm_tcga_gdc | 353450 | 1m2.537s | 15.184s | 24.3% | 13460 |
| **Total** | **1052366** | **4m14.447s** | **1m3.378s** | **24.9%** | **9813** |
