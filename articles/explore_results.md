# Exploring fine-mapping results

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

This vignette demonstrates how to explore and interpret the fine-mapping
results produced by
[`finemap_loci()`](https://rajlabmssm.github.io/echolocatoR/reference/finemap_loci.md).
All examples use bundled data included in `echodata`, so **no internet
connection is required**.

## Load example results

`echodata` ships with pre-computed fine-mapping results for three
Parkinson’s disease loci from Nalls et al. (2019):

``` r
dat <- echodata::BST1
dim(dat)
```

    ## [1] 6216   26

Each row is a SNP within the locus window. Key columns include:

- **GWAS columns**: `SNP`, `CHR`, `POS`, `P`, `Effect`, `StdErr`, `MAF`
- **Fine-mapping posterior probabilities**: `ABF.PP`, `SUSIE.PP`,
  `FINEMAP.PP`
- **Credible set assignments**: `ABF.CS`, `SUSIE.CS`, `FINEMAP.CS`
- **Consensus columns**: `Support` (number of methods placing SNP in a
  CS), `Consensus_SNP`, `mean.PP`, `mean.CS`

``` r
## Fine-mapping specific columns
fm_cols <- c("SNP","P","ABF.PP","SUSIE.PP","FINEMAP.PP",
             "Support","Consensus_SNP","mean.PP")
knitr::kable(head(dat[order(dat$mean.PP, decreasing = TRUE), ..fm_cols], 10),
             digits = 4, row.names = FALSE,
             caption = "Top 10 SNPs by mean posterior probability")
```

| SNP        |   P | ABF.PP | SUSIE.PP | FINEMAP.PP | Support | Consensus_SNP | mean.PP |
|:-----------|----:|-------:|---------:|-----------:|--------:|:--------------|--------:|
| rs4541502  |   0 | 0.0000 |        1 |     1.0000 |       3 | TRUE          |  0.7500 |
| rs34559912 |   0 | 0.0621 |        1 |     0.0000 |       2 | TRUE          |  0.5155 |
| rs4389574  |   0 | 0.0000 |        1 |     0.0000 |       2 | TRUE          |  0.5000 |
| rs4698412  |   0 | 0.3871 |        0 |     1.0000 |       1 | FALSE         |  0.3468 |
| rs6852450  |   0 | 0.0000 |        0 |     1.0000 |       1 | FALSE         |  0.2500 |
| rs3756246  |   0 | 0.0000 |        0 |     1.0000 |       1 | FALSE         |  0.2500 |
| rs35519415 |   0 | 0.0000 |        0 |     0.9994 |       1 | FALSE         |  0.2498 |
| rs11724635 |   0 | 0.3064 |        0 |     0.0000 |       0 | FALSE         |  0.0766 |
| rs4698413  |   0 | 0.2427 |        0 |     0.0000 |       0 | FALSE         |  0.0607 |
| rs4613561  |   0 | 0.0014 |        0 |     0.0000 |       0 | FALSE         |  0.0004 |

Top 10 SNPs by mean posterior probability

## Consensus SNPs

Consensus SNPs are those placed in a credible set by multiple
fine-mapping methods. The `Support` column indicates how many methods
agree:

``` r
consensus <- subset(dat, Consensus_SNP == TRUE)
knitr::kable(
    consensus[, c("SNP","POS","P","ABF.PP","SUSIE.PP","FINEMAP.PP",
                   "Support","mean.PP")],
    digits = 4, row.names = FALSE,
    caption = "Consensus fine-mapped SNPs (BST1 locus)"
)
```

| SNP        |      POS |   P | ABF.PP | SUSIE.PP | FINEMAP.PP | Support | mean.PP |
|:-----------|---------:|----:|-------:|---------:|-----------:|--------:|--------:|
| rs4541502  | 15712787 |   0 | 0.0000 |        1 |          1 |       3 |  0.7500 |
| rs34559912 | 15730146 |   0 | 0.0621 |        1 |          0 |       2 |  0.5155 |
| rs4389574  | 15730398 |   0 | 0.0000 |        1 |          0 |       2 |  0.5000 |

Consensus fine-mapped SNPs (BST1 locus)

## SNP group counts

Summarise how many SNPs fall into each support level:

``` r
support_table <- table(Support = dat$Support)
knitr::kable(as.data.frame(support_table),
             caption = "SNP counts by support level")
```

| Support | Freq |
|:--------|-----:|
| 0       | 6209 |
| 1       |    4 |
| 2       |    2 |
| 3       |    1 |

SNP counts by support level

## Credible sets

Each fine-mapping method assigns SNPs to credible sets (CS). A value of
`0` means the SNP is not in any CS; `1` = first CS, `2` = second CS,
etc.

``` r
## Count SNPs in each method's credible sets
cs_summary <- data.frame(
    Method = c("ABF", "SUSIE", "FINEMAP"),
    SNPs_in_CS = c(
        sum(dat$ABF.CS > 0, na.rm = TRUE),
        sum(dat$SUSIE.CS > 0, na.rm = TRUE),
        sum(dat$FINEMAP.CS > 0, na.rm = TRUE)
    )
)
knitr::kable(cs_summary, caption = "Number of SNPs in credible sets per method")
```

| Method  | SNPs_in_CS |
|:--------|-----------:|
| ABF     |          0 |
| SUSIE   |          3 |
| FINEMAP |          5 |

Number of SNPs in credible sets per method

## Working with multiple loci

`echodata` also bundles LRRK2 and MEX3C loci. You can combine them to
simulate multi-locus analysis:

``` r
bst1 <- data.table::copy(echodata::BST1); bst1[, Locus := "BST1"]
lrrk2 <- data.table::copy(echodata::LRRK2); lrrk2[, Locus := "LRRK2"]
mex3c <- data.table::copy(echodata::MEX3C); mex3c[, Locus := "MEX3C"]
shared_cols <- Reduce(intersect, list(colnames(bst1), colnames(lrrk2), colnames(mex3c)))
all_loci <- rbind(bst1[, ..shared_cols], lrrk2[, ..shared_cols], mex3c[, ..shared_cols])
## Consensus SNPs per locus
consensus_per_locus <- do.call(rbind, lapply(
    split(all_loci, all_loci$Locus),
    function(x) data.frame(
        Locus = x$Locus[1],
        Total_SNPs = nrow(x),
        Consensus_SNPs = sum(x$Consensus_SNP, na.rm = TRUE),
        Min_P = min(x$P, na.rm = TRUE)
    )
))
knitr::kable(consensus_per_locus, row.names = FALSE,
             caption = "Summary across bundled loci")
```

| Locus | Total_SNPs | Consensus_SNPs | Min_P |
|:------|-----------:|---------------:|------:|
| BST1  |       6216 |              3 |     0 |
| LRRK2 |       3844 |              4 |     0 |
| MEX3C |       5372 |              2 |     0 |

Summary across bundled loci

## LD matrix

A pre-computed LD matrix for BST1 is also available:

``` r
ld <- echodata::BST1_LD_matrix
cat("LD matrix dimensions:", dim(ld), "\n")
```

    ## LD matrix dimensions: 95 95

``` r
cat("SNPs in LD matrix:", ncol(ld), "\n")
```

    ## SNPs in LD matrix: 95

``` r
## Show a small corner
knitr::kable(round(ld[1:5, 1:5], 3), caption = "LD matrix (top-left corner)")
```

|            | rs4698412 | rs6852450 | rs10000232 | rs10000250 | rs10000290 |
|:-----------|----------:|----------:|-----------:|-----------:|-----------:|
| rs4698412  |     1.000 |     0.885 |     -0.023 |     -0.045 |     -0.431 |
| rs6852450  |     0.885 |     1.000 |     -0.006 |     -0.044 |     -0.464 |
| rs10000232 |    -0.023 |    -0.006 |      1.000 |     -0.017 |      0.024 |
| rs10000250 |    -0.045 |    -0.044 |     -0.017 |      1.000 |      0.089 |
| rs10000290 |    -0.431 |    -0.464 |      0.024 |      0.089 |      1.000 |

LD matrix (top-left corner)

## Next steps

- Visualize loci:
  [`vignette("plotting")`](https://rajlabmssm.github.io/echolocatoR/articles/plotting.md)
- Run the full pipeline:
  [`vignette("echolocatoR")`](https://rajlabmssm.github.io/echolocatoR/articles/echolocatoR.md)
- Summarise across loci:
  [`vignette("summarise")`](https://rajlabmssm.github.io/echolocatoR/articles/summarise.md)
- Learn about sub-packages:
  [`vignette("echoverse_modules")`](https://rajlabmssm.github.io/echolocatoR/articles/echoverse_modules.md)

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
