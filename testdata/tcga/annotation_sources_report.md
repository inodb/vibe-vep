# Annotation Sources Report

Generated: 2026-02-28 21:04 UTC  
GENCODE transcripts: 254070  
AlphaMissense variants in database: 71697556  
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

## AlphaMissense Performance

| Study | Variants | Base Time | AM Lookup Time | AM Overhead | Lookups/sec |
|-------|----------|-----------|----------------|-------------|-------------|
| blca_tcga_gdc | 115850 | 51.2s | 6.976s | 13.6% | 10886 |
| brca_tcga_gdc | 89012 | 20.321s | 6.444s | 31.7% | 8954 |
| chol_tcga_gdc | 3764 | 851ms | 2.25s | 264.4% | 1038 |
| coad_tcga_gdc | 244552 | 34.441s | 8.975s | 26.1% | 14334 |
| gbm_tcga_gdc | 54870 | 10.215s | 4.881s | 47.8% | 6696 |
| luad_tcga_gdc | 190868 | 32.488s | 7.845s | 24.1% | 15329 |
| skcm_tcga_gdc | 353450 | 51.329s | 12.747s | 24.8% | 16033 |
| **Total** | **1052366** | **3m20.845s** | **50.118s** | **25.0%** | **12409** |
