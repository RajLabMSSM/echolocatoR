# echoverse modules

``` r
library(echolocatoR)
```

    ## 

    ## ── 🦇  🦇  🦇 e c h o l o c a t o R 🦇  🦇  🦇 ─────────────────────────────────

    ## 

    ## ── v2.0.5 ──────────────────────────────────────────────────────────────────────

    ## 

    ## ────────────────────────────────────────────────────────────────────────────────

    ## ⠊⠉⠡⣀⣀⠊⠉⠡⣀⣀⠊⠉⠢⣀⡠⠊⠉⠢⣀⡠⠊⠉⠢⣀⡠⠊⠉⠢⣀⡠⠊⠉⠢⣀⡠⠊⠉⠢⣀⡠                                    
    ## ⠌⢁⡐⠉⣀⠊⢂⡐⠑⣀⠊⢂⡐⠑⣀⠊⢂⡐⠑⣀⠊⢂⡐⠑⣀⠊⢂⡐⠑⣀⠉⢂⡈⠑⣀⠉⢄⡈⠡⣀                                    
    ## ⠌⡈⡐⢂⢁⠒⡈⡐⢂⢁⠒⡈⡐⢂⢁⠑⡈⡈⢄⢁⠡⠌⡈⠤⢁⠡⠌⡈⠤⢁⠡⠌⡈⡠⢁⢁⠊⡈⡐⢂                                    
    ## ⠌⡈⡐⢂⢁⠒⡈⡐⢂⢁⠒⡈⡐⢂⢁⠑⡈⡈⢄⢁⠡⠌⡈⠤⢁⠡⠌⡈⠤⢁⠡⠌⡈⡠⢁⢁⠊⡈⡐⢂                                    
    ## ⠌⢁⡐⠉⣀⠊⢂⡐⠑⣀⠊⢂⡐⠑⣀⠊⢂⡐⠑⣀⠊⢂⡐⠑⣀⠊⢂⡐⠑⣀⠉⢂⡈⠑⣀⠉⢄⡈⠡⣀                                    
    ## ⠊⠉⠡⣀⣀⠊⠉⠡⣀⣀⠊⠉⠢⣀⡠⠊⠉⠢⣀⡠⠊⠉⠢⣀⡠⠊⠉⠢⣀⡠⠊⠉⠢⣀⡠⠊⠉⠢⣀⡠                                    
    ## ⓞ If you use echolocatoR or any of the echoverse subpackages, please cite:      
    ##      ▶ Brian M Schilder, Jack Humphrey, Towfique                                
    ##      Raj (2021) echolocatoR: an automated                                       
    ##      end-to-end statistical and functional                                      
    ##      genomic fine-mapping pipeline,                                             
    ##      Bioinformatics; btab658,                                                   
    ##      https://doi.org/10.1093/bioinformatics/btab658                             
    ## ⓞ Please report any bugs/feature requests on GitHub:
    ##      ▶
    ##      https://github.com/RajLabMSSM/echolocatoR/issues
    ## ⓞ Contributions are welcome!:
    ##      ▶
    ##      https://github.com/RajLabMSSM/echolocatoR/pulls

## Overview

`echolocatoR` orchestrates a suite of modular R packages collectively
called the **echoverse**. While
[`echolocatoR::finemap_loci()`](https://rajlabmssm.github.io/echolocatoR/reference/finemap_loci.md)
runs the full pipeline automatically, each module can also be used
independently for custom workflows.

This vignette provides a quick reference for when to use each module.

## Module reference

### Level 1: Base packages

These provide foundational utilities used by all other modules.

#### `echodata`

**Purpose**: Data management, example datasets, column standardization.

**Use independently when you want to**:

- Load bundled fine-mapping results: `BST1`, `LRRK2`, `MEX3C`
- Download full GWAS summary stats:
  [`get_Nalls2019()`](https://rdrr.io/pkg/echodata/man/get_Nalls2019.html),
  [`get_Kunkle2019()`](https://rdrr.io/pkg/echodata/man/get_Kunkle2019.html)
- Standardize column names across different GWAS formats
- Access the Fine-Mapping Portal:
  [`portal_query()`](https://rdrr.io/pkg/echodata/man/portal_query.html)

``` r
## Bundled fine-mapping results (no download needed)
dat <- echodata::BST1
cat("BST1 locus:", nrow(dat), "SNPs,", ncol(dat), "columns\n")
```

    ## BST1 locus: 6216 SNPs, 26 columns

#### `echogithub`

**Purpose**: GitHub API utilities for downloading files and releases.

**Use independently when you want to**:

- Download files from GitHub repos programmatically
- Access GitHub release assets

#### `downloadR`

**Purpose**: File download helpers with caching.

**Use independently when you want to**:

- Download and cache large files with progress tracking
- Upload data to Zenodo

#### `devoptera`

**Purpose**: Developer utilities for the echoverse packages.

### Level 2: Mid-level utilities

#### `echotabix`

**Purpose**: Tabix-indexed file querying and liftover.

**Use independently when you want to**:

- Query tabix-indexed VCF or summary stats files by genomic region
- Perform genome build liftover (hg19/hg38)
- Convert between file formats

#### `echoLD`

**Purpose**: LD (linkage disequilibrium) matrix computation.

**Use independently when you want to**:

- Compute LD matrices from 1000 Genomes Phase 3, UK Biobank, or custom
  VCFs
- Filter variants by LD or MAF thresholds

``` r
## Bundled LD matrix
ld <- echodata::BST1_LD_matrix
cat("LD matrix:", nrow(ld), "x", ncol(ld), "SNPs\n")
```

    ## LD matrix: 95 x 95 SNPs

#### `echoconda`

**Purpose**: Conda environment management for external tools.

**Use independently when you want to**:

- Set up conda environments for FINEMAP, PolyFun, or other command-line
  tools
- Manage Python dependencies from R via basilisk

### Level 3: Analysis modules

#### `echofinemap`

**Purpose**: Statistical fine-mapping methods.

**Use independently when you want to**:

- Run individual fine-mapping methods (ABF, SUSIE, FINEMAP,
  PolyFun+SUSIE)
- Compare results across methods
- See available methods:
  [`echofinemap::lfm()`](https://rdrr.io/pkg/echofinemap/man/lfm.html)

``` r
## List available fine-mapping methods
echofinemap::lfm()
```

    ## Gathering method sources.

    ## Gathering method citations.

    ## [1] "ABF"              "COJO_conditional" "COJO_joint"       "COJO_stepwise"   
    ## [5] "FINEMAP"          "PAINTOR"          "POLYFUN_FINEMAP"  "POLYFUN_SUSIE"   
    ## [9] "SUSIE"

#### `echoannot`

**Purpose**: Functional annotation and enrichment.

**Use independently when you want to**:

- Overlap fine-mapped SNPs with ROADMAP epigenomic annotations
- Run GoShifter enrichment tests
- Annotate with NOTT2019/CORCES2020 single-cell data
- Create super summary plots

#### `echoplot`

**Purpose**: Locus visualization.

**Use independently when you want to**:

- Create multi-track locus zoom plots
- Visualize fine-mapping posterior probabilities alongside GWAS p-values
- Overlay LD, gene tracks, and annotations

#### `catalogueR`

**Purpose**: QTL colocalization via the eQTL Catalogue.

**Use independently when you want to**:

- Query the eQTL Catalogue for colocalization analysis
- Access tissue/cell-type-specific QTL datasets

### Level 4: Pipeline orchestrator

#### `echolocatoR`

**Purpose**: Unified fine-mapping pipeline.

This is the package you’re using now. It calls all modules above in the
correct order via
[`finemap_loci()`](https://rajlabmssm.github.io/echolocatoR/reference/finemap_loci.md)
/
[`finemap_locus()`](https://rajlabmssm.github.io/echolocatoR/reference/finemap_locus.md).

## Choosing your approach

| Scenario | Recommended approach |
|----|----|
| First-time fine-mapping | [`echolocatoR::finemap_loci()`](https://rajlabmssm.github.io/echolocatoR/reference/finemap_loci.md) — handles everything |
| Custom fine-mapping workflow | Use `echofinemap` directly with your own LD matrices |
| Just need LD matrices | Use [`echoLD::get_LD()`](https://rdrr.io/pkg/echoLD/man/get_LD.html) |
| Just need annotations | Use `echoannot` functions directly |
| Just need plots | Use [`echoplot::plot_locus()`](https://rdrr.io/pkg/echoplot/man/plot_locus.html) |
| Querying large summary stats | Use [`echotabix::query()`](https://rdrr.io/pkg/echotabix/man/query.html) |
| Post-hoc analysis of results | Use `echodata` to load + `echoannot` to annotate |

## Session info

``` r
utils::sessionInfo()
```

```
## R version 4.5.1 (2025-06-13)
## Platform: aarch64-apple-darwin20
## Running under: macOS Tahoe 26.3.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: America/New_York
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] echolocatoR_2.0.5 BiocStyle_2.38.0 
## 
## loaded via a namespace (and not attached):
##   [1] splines_4.5.1               aws.s3_0.3.22              
##   [3] BiocIO_1.20.0               bitops_1.0-9               
##   [5] filelock_1.0.3              tibble_3.3.1               
##   [7] R.oo_1.27.1                 cellranger_1.1.0           
##   [9] basilisk.utils_1.22.0       graph_1.88.1               
##  [11] rpart_4.1.24                XML_3.99-0.22              
##  [13] lifecycle_1.0.5             mixsqp_0.3-54              
##  [15] pals_1.10                   OrganismDbi_1.52.0         
##  [17] ensembldb_2.34.0            lattice_0.22-9             
##  [19] MASS_7.3-65                 backports_1.5.0            
##  [21] magrittr_2.0.4              Hmisc_5.2-5                
##  [23] openxlsx_4.2.8.1            sass_0.4.10                
##  [25] rmarkdown_2.30              jquerylib_0.1.4            
##  [27] yaml_2.3.12                 otel_0.2.0                 
##  [29] zip_2.3.3                   reticulate_1.45.0          
##  [31] ggbio_1.58.0                gld_2.6.8                  
##  [33] mapproj_1.2.12              DBI_1.3.0                  
##  [35] RColorBrewer_1.1-3          maps_3.4.3                 
##  [37] abind_1.4-8                 expm_1.0-0                 
##  [39] GenomicRanges_1.62.1        purrr_1.2.1                
##  [41] R.utils_2.13.0              AnnotationFilter_1.34.0    
##  [43] biovizBase_1.58.0           BiocGenerics_0.56.0        
##  [45] RCurl_1.98-1.17             nnet_7.3-20                
##  [47] VariantAnnotation_1.56.0    IRanges_2.44.0             
##  [49] S4Vectors_0.48.0            echofinemap_1.0.0          
##  [51] echoLD_0.99.12              catalogueR_2.0.1           
##  [53] irlba_2.3.7                 pkgdown_2.2.0              
##  [55] echodata_1.0.0              piggyback_0.1.5            
##  [57] codetools_0.2-20            DelayedArray_0.36.0        
##  [59] DT_0.34.0                   xml2_1.5.2                 
##  [61] tidyselect_1.2.1            UCSC.utils_1.6.1           
##  [63] farver_2.1.2                viridis_0.6.5              
##  [65] matrixStats_1.5.0           stats4_4.5.1               
##  [67] base64enc_0.1-6             Seqinfo_1.0.0              
##  [69] echotabix_1.0.0             GenomicAlignments_1.46.0   
##  [71] jsonlite_2.0.0              e1071_1.7-17               
##  [73] Formula_1.2-5               survival_3.8-6             
##  [75] systemfonts_1.3.2           ggnewscale_0.5.2           
##  [77] tools_4.5.1                 ragg_1.5.1                 
##  [79] DescTools_0.99.60           Rcpp_1.1.1                 
##  [81] glue_1.8.0                  gridExtra_2.3              
##  [83] SparseArray_1.10.9          xfun_0.56                  
##  [85] MatrixGenerics_1.22.0       GenomeInfoDb_1.46.2        
##  [87] dplyr_1.2.0                 withr_3.0.2                
##  [89] BiocManager_1.30.27         fastmap_1.2.0              
##  [91] basilisk_1.22.0             boot_1.3-32                
##  [93] digest_0.6.39               R6_2.6.1                   
##  [95] colorspace_2.1-2            textshaping_1.0.5          
##  [97] dichromat_2.0-0.1           RSQLite_2.4.6              
##  [99] cigarillo_1.0.0             R.methodsS3_1.8.2          
## [101] utf8_1.2.6                  tidyr_1.3.2                
## [103] generics_0.1.4              data.table_1.18.2.1        
## [105] rtracklayer_1.70.1          class_7.3-23               
## [107] httr_1.4.8                  htmlwidgets_1.6.4          
## [109] S4Arrays_1.10.1             pkgconfig_2.0.3            
## [111] gtable_0.3.6                Exact_3.3                  
## [113] blob_1.3.0                  S7_0.2.1                   
## [115] XVector_0.50.0              echoconda_1.0.0            
## [117] htmltools_0.5.9             susieR_0.14.2              
## [119] bookdown_0.46               RBGL_1.86.0                
## [121] ProtGenerics_1.42.0         scales_1.4.0               
## [123] Biobase_2.70.0              lmom_3.2                   
## [125] png_0.1-8                   knitr_1.51                 
## [127] rstudioapi_0.18.0           reshape2_1.4.5             
## [129] tzdb_0.5.0                  rjson_0.2.23               
## [131] checkmate_2.3.4             curl_7.0.0                 
## [133] proxy_0.4-29                cachem_1.1.0               
## [135] stringr_1.6.0               rootSolve_1.8.2.4          
## [137] parallel_4.5.1              foreign_0.8-91             
## [139] AnnotationDbi_1.72.0        restfulr_0.0.16            
## [141] desc_1.4.3                  pillar_1.11.1              
## [143] grid_4.5.1                  reshape_0.8.10             
## [145] vctrs_0.7.1                 cluster_2.1.8.2            
## [147] htmlTable_2.4.3             evaluate_1.0.5             
## [149] readr_2.2.0                 GenomicFeatures_1.62.0     
## [151] mvtnorm_1.3-3               cli_3.6.5                  
## [153] compiler_4.5.1              Rsamtools_2.26.0           
## [155] rlang_1.1.7                 crayon_1.5.3               
## [157] aws.signature_0.6.0         plyr_1.8.9                 
## [159] forcats_1.0.1               fs_1.6.7                   
## [161] stringi_1.8.7               coloc_5.2.3                
## [163] echoannot_1.0.1             viridisLite_0.4.3          
## [165] BiocParallel_1.44.0         Biostrings_2.78.0          
## [167] lazyeval_0.2.2              Matrix_1.7-4               
## [169] downloadR_1.0.0             echoplot_0.99.9            
## [171] dir.expiry_1.18.0           BSgenome_1.78.0            
## [173] hms_1.1.4                   patchwork_1.3.2            
## [175] bit64_4.6.0-1               ggplot2_4.0.2              
## [177] KEGGREST_1.50.0             SummarizedExperiment_1.40.0
## [179] haven_2.5.5                 memoise_2.0.1              
## [181] snpStats_1.60.0             bslib_0.10.0               
## [183] bit_4.6.0                   readxl_1.4.5
```
