# Summarise

``` r

library(echolocatoR)
```

    ## 

    ## ── 🦇  🦇  🦇 e c h o l o c a t o R 🦇  🦇  🦇 ─────────────────────────────────

    ## 

    ## ── v3.0.0 ──────────────────────────────────────────────────────────────────────

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
    ##      Raj (2021). echolocatoR: an automated                                      
    ##      end-to-end statistical and functional                                      
    ##      genomic fine-mapping pipeline.                                             
    ##      Bioinformatics, btab658.                                                   
    ##      https://doi.org/10.1093/bioinformatics/btab658                             
    ## ⓞ Please report any bugs/feature requests on GitHub:
    ##      ▶
    ##      https://github.com/RajLabMSSM/echolocatoR/issues
    ## ⓞ Contributions are welcome!:
    ##      ▶
    ##      https://github.com/RajLabMSSM/echolocatoR/pulls
    ## Tip: Run check_echoverse_setup() to diagnose any issues.

``` r

has_internet <- tryCatch(
    !is.null(curl::nslookup("github.com", error = FALSE)),
    error = function(e) FALSE
)
```

## Download data

## Summarise

Using pre-merged data for vignette speed.

``` r

merged_DT <- echodata::get_Nalls2019_merged()
```

### `get_SNPgroup_counts()`

Get the number of SNPs for each SNP group per locus. It also prints the
mean number of SNPs for each SNP group across all loci. **NOTE**: You
will need to make sure to set
`merge_finemapping_results(minimum_support=1)` in the above step to get
accurate counts for all SNP groups.

``` r

snp_groups <- echodata::get_SNPgroup_counts(merged_DT = merged_DT)
```

    ## All loci (75) :

    ##            Total.SNPs          nom.sig.GWAS              sig.GWAS 
    ##               4948.16                924.68                 82.36 
    ##                    CS             Consensus          topConsensus 
    ##                  7.88                  2.69                  1.47 
    ## topConsensus.leadGWAS 
    ##                  0.41

    ## Loci with at least one Consensus SNP (69) :

    ##            Total.SNPs          nom.sig.GWAS              sig.GWAS 
    ##               5019.07                911.41                 84.28 
    ##                    CS             Consensus          topConsensus 
    ##                  7.77                  2.93                  1.59 
    ## topConsensus.leadGWAS 
    ##                  0.45

### `get_CS_counts()`

Count the number of tool-specific and UCS Credible Set SNPs per locus.

``` r

UCS_counts <- echodata::get_CS_counts(merged_DT = merged_DT)
knitr::kable(head(UCS_counts, 10))
```

| Locus | ABF.CS_size | SUSIE.CS_size | POLYFUN_SUSIE.CS_size | FINEMAP.CS_size | mean.CS_size | UCS.CS_size |
|:---|---:|---:|---:|---:|---:|---:|
| GPNMB | 0 | 4 | 5 | 5 | 0 | 13 |
| MAP4K4 | 0 | 4 | 3 | 5 | 0 | 12 |
| MBNL2 | 0 | 4 | 4 | 5 | 0 | 12 |
| TMEM163 | 0 | 5 | 5 | 5 | 0 | 12 |
| CRLS1 | 0 | 3 | 4 | 5 | 0 | 11 |
| DNAH17 | 0 | 5 | 4 | 5 | 0 | 11 |
| GBF1 | 0 | 5 | 5 | 5 | 0 | 11 |
| GCH1 | 0 | 5 | 2 | 4 | 0 | 11 |
| MIPOL1 | 0 | 2 | 4 | 5 | 0 | 11 |
| FBRSL1 | 0 | 3 | 3 | 5 | 0 | 10 |

## Plot

- The following functions each return a list containing both the
  `...$plot` and the `...$data` used to make the plot.
- Where available, `snp_filter` allows user to use any filtering
  argument (supplied as a string) to subset the data they want to use in
  the plot/data.

### Colocalization results

If you ran colocalization tests with `echolocatoR` (via `catalogueR`)
you can use those results to come up with a top QTL nominated gene for
each locus (potentially implicating that gene in your phenotype).

``` r

coloc_res <- echodata::get_Nalls2019_coloc()
```

### Super summary plot

``` r

super_plot <- echoannot::super_summary_plot(merged_DT = merged_DT,
                                            coloc_results = coloc_res,
                                            plot_missense = FALSE)
```

    ## + SUMMARISE:: Nominating genes by top colocalized eQTL eGenes

    ## Warning in ggplot2::geom_bar(stat = "identity", color = "white", size = 0.05): Ignoring unknown parameters: `size`
    ## Ignoring unknown parameters: `size`

    ## Importing previously downloaded files: /github/home/.cache/R/echoannot/NOTT2019_epigenomic_peaks.rds

    ## ++ NOTT2019:: 634,540 ranges retrieved.

    ## Converting dat to GRanges object.

    ## 113 query SNP(s) detected with reference overlap.

    ## ++ NOTT2019:: Getting regulatory regions data.

    ## Importing Astrocyte enhancers ...

    ## Importing Astrocyte promoters ...

    ## Importing Neuronal enhancers ...

    ## Importing Neuronal promoters ...

    ## Importing Oligo enhancers ...

    ## Importing Oligo promoters ...

    ## Importing Microglia enhancers ...

    ## Importing Microglia promoters ...

    ## Converting dat to GRanges object.
    ## Converting dat to GRanges object.

    ## 48 query SNP(s) detected with reference overlap.

    ## ++ NOTT2019:: Getting interaction anchors data.

    ## Importing Microglia interactome ...

    ## Importing Neuronal interactome ...

    ## Importing Oligo interactome ...

    ## Converting dat to GRanges object.

    ## 52 query SNP(s) detected with reference overlap.

    ## Converting dat to GRanges object.

    ## 44 query SNP(s) detected with reference overlap.

    ## CORCES2020:: Extracting overlapping cell-type-specific scATAC-seq peaks

    ## Converting dat to GRanges object.

    ## 13 query SNP(s) detected with reference overlap.

    ## CORCES2020:: Annotating peaks by cell-type-specific target genes

    ## CORCES2020:: Extracting overlapping bulkATAC-seq peaks from brain tissue

    ## Converting dat to GRanges object.

    ## 4 query SNP(s) detected with reference overlap.

    ## CORCES2020:: Annotating peaks by bulk brain target genes

    ## Converting dat to GRanges object.

    ## 70 query SNP(s) detected with reference overlap.

    ## Converting dat to GRanges object.

    ## 72 query SNP(s) detected with reference overlap.

    ## + CORCES2020:: Found 142 hits with HiChIP_FitHiChIP coaccessibility loop anchors.

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## ℹ The deprecated feature was likely used in the echoannot package.
    ##   Please report the issue at <https://github.com/RajLabMSSM/echoannot/issues>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## ℹ The deprecated feature was likely used in the echoannot package.
    ##   Please report the issue at <https://github.com/RajLabMSSM/echoannot/issues>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation idioms with `aes()`.
    ## ℹ See also `vignette("ggplot2-in-packages")` for more information.
    ## ℹ The deprecated feature was likely used in the echoannot package.
    ##   Please report the issue at <https://github.com/RajLabMSSM/echoannot/issues>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: The dot-dot notation (`..count..`) was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `after_stat(count)` instead.
    ## ℹ The deprecated feature was likely used in the echoannot package.
    ##   Please report the issue at <https://github.com/RajLabMSSM/echoannot/issues>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: Removed 83 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](summarise_files/figure-html/super_summary_plot()-1.png)

## Session info

``` r

utils::sessionInfo()
```

```
## R Under development (unstable) (2026-03-12 r89607)
## Platform: x86_64-pc-linux-gnu
## Running under: Ubuntu 24.04.4 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: UTC
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] echolocatoR_3.0.0 BiocStyle_2.39.0 
## 
## loaded via a namespace (and not attached):
##   [1] splines_4.6.0               aws.s3_0.3.22              
##   [3] BiocIO_1.21.0               bitops_1.0-9               
##   [5] filelock_1.0.3              tibble_3.3.1               
##   [7] R.oo_1.27.1                 cellranger_1.1.0           
##   [9] basilisk.utils_1.23.1       graph_1.89.1               
##  [11] rpart_4.1.24                XML_3.99-0.22              
##  [13] lifecycle_1.0.5             mixsqp_0.3-54              
##  [15] pals_1.10                   OrganismDbi_1.53.2         
##  [17] ensembldb_2.35.0            lattice_0.22-9             
##  [19] MASS_7.3-65                 backports_1.5.0            
##  [21] magrittr_2.0.4              Hmisc_5.2-5                
##  [23] openxlsx_4.2.8.1            sass_0.4.10                
##  [25] rmarkdown_2.30              jquerylib_0.1.4            
##  [27] yaml_2.3.12                 otel_0.2.0                 
##  [29] zip_2.3.3                   reticulate_1.45.0          
##  [31] ggbio_1.59.0                gld_2.6.8                  
##  [33] mapproj_1.2.12              DBI_1.3.0                  
##  [35] RColorBrewer_1.1-3          maps_3.4.3                 
##  [37] abind_1.4-8                 expm_1.0-0                 
##  [39] GenomicRanges_1.63.1        purrr_1.2.1                
##  [41] R.utils_2.13.0              AnnotationFilter_1.35.0    
##  [43] biovizBase_1.59.0           BiocGenerics_0.57.0        
##  [45] RCurl_1.98-1.17             nnet_7.3-20                
##  [47] VariantAnnotation_1.57.1    IRanges_2.45.0             
##  [49] S4Vectors_0.49.0            echoLD_1.0.0               
##  [51] echofinemap_1.0.0           irlba_2.3.7                
##  [53] pkgdown_2.2.0               echodata_1.0.0             
##  [55] piggyback_0.1.5             codetools_0.2-20           
##  [57] DelayedArray_0.37.0         DT_0.34.0                  
##  [59] xml2_1.5.2                  tidyselect_1.2.1           
##  [61] UCSC.utils_1.7.1            farver_2.1.2               
##  [63] viridis_0.6.5               matrixStats_1.5.0          
##  [65] stats4_4.6.0                base64enc_0.1-6            
##  [67] Seqinfo_1.1.0               echotabix_1.0.1            
##  [69] GenomicAlignments_1.47.0    jsonlite_2.0.0             
##  [71] e1071_1.7-17                Formula_1.2-5              
##  [73] survival_3.8-6              systemfonts_1.3.2          
##  [75] ggnewscale_0.5.2            tools_4.6.0                
##  [77] ragg_1.5.1                  DescTools_0.99.60          
##  [79] Rcpp_1.1.1                  glue_1.8.0                 
##  [81] gridExtra_2.3               SparseArray_1.11.11        
##  [83] xfun_0.56                   MatrixGenerics_1.23.0      
##  [85] GenomeInfoDb_1.47.2         dplyr_1.2.0                
##  [87] withr_3.0.2                 BiocManager_1.30.27        
##  [89] fastmap_1.2.0               basilisk_1.23.0            
##  [91] boot_1.3-32                 digest_0.6.39              
##  [93] R6_2.6.1                    colorspace_2.1-2           
##  [95] textshaping_1.0.5           dichromat_2.0-0.1          
##  [97] RSQLite_2.4.6               cigarillo_1.1.0            
##  [99] R.methodsS3_1.8.2           utf8_1.2.6                 
## [101] tidyr_1.3.2                 generics_0.1.4             
## [103] data.table_1.18.2.1         rtracklayer_1.71.3         
## [105] class_7.3-23                httr_1.4.8                 
## [107] htmlwidgets_1.6.4           S4Arrays_1.11.1            
## [109] pkgconfig_2.0.3             gtable_0.3.6               
## [111] Exact_3.3                   blob_1.3.0                 
## [113] S7_0.2.1                    XVector_0.51.0             
## [115] echoconda_1.0.0             htmltools_0.5.9            
## [117] susieR_0.14.2               bookdown_0.46              
## [119] RBGL_1.87.0                 ProtGenerics_1.43.0        
## [121] scales_1.4.0                Biobase_2.71.0             
## [123] lmom_3.2                    png_0.1-8                  
## [125] knitr_1.51                  rstudioapi_0.18.0          
## [127] tzdb_0.5.0                  reshape2_1.4.5             
## [129] rjson_0.2.23                checkmate_2.3.4            
## [131] curl_7.0.0                  proxy_0.4-29               
## [133] cachem_1.1.0                stringr_1.6.0              
## [135] rootSolve_1.8.2.4           parallel_4.6.0             
## [137] foreign_0.8-91              AnnotationDbi_1.73.0       
## [139] restfulr_0.0.16             desc_1.4.3                 
## [141] pillar_1.11.1               grid_4.6.0                 
## [143] reshape_0.8.10              vctrs_0.7.1                
## [145] cluster_2.1.8.2             htmlTable_2.4.3            
## [147] evaluate_1.0.5              readr_2.2.0                
## [149] GenomicFeatures_1.63.1      mvtnorm_1.3-5              
## [151] cli_3.6.5                   compiler_4.6.0             
## [153] Rsamtools_2.27.1            rlang_1.1.7                
## [155] crayon_1.5.3                labeling_0.4.3             
## [157] aws.signature_0.6.0         plyr_1.8.9                 
## [159] forcats_1.0.1               fs_1.6.7                   
## [161] stringi_1.8.7               coloc_5.2.3                
## [163] echoannot_1.0.1             viridisLite_0.4.3          
## [165] BiocParallel_1.45.0         Biostrings_2.79.5          
## [167] lazyeval_0.2.2              Matrix_1.7-4               
## [169] downloadR_1.0.0             echoplot_1.0.0             
## [171] dir.expiry_1.19.0           BSgenome_1.79.1            
## [173] patchwork_1.3.2             hms_1.1.4                  
## [175] bit64_4.6.0-1               ggplot2_4.0.2              
## [177] KEGGREST_1.51.1             SummarizedExperiment_1.41.1
## [179] haven_2.5.5                 memoise_2.0.1              
## [181] snpStats_1.61.1             bslib_0.10.0               
## [183] bit_4.6.0                   readxl_1.4.5
```
