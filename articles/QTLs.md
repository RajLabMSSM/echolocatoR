# QTLs

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

## QTL pipeline

- Here, we will use GWAS-eQTL colocalization results provided via the
  [*echolocatoR Fine-mapping
  Portal*](https://github.com/RajLabMSSM/Fine_Mapping_Shiny), from the
  publication:

> Lopes, K.d.P., Snijders, G.J.L., Humphrey, J. et al. Genetic analysis
> of the human microglial transcriptome across brain regions, aging and
> disease pathologies. Nat Genet 54, 4–17 (2022).
> <https://doi.org/10.1038/s41588-021-00976-y>

### Import QTL data

This data is actually merged GWAS-QTL colocalization results, but it
contains all of the necessary columns from the original eQTL summary
stats that we need to perform eQTL fine-mapping.

``` r

coloc_res <- echodata::get_Kunkle2019_coloc(return_path = TRUE)
```

### Prepare `colmap`

Prepare a column mapping object for the summary statistics. We’ll reuse
this for both the `import_topSNPs` and `finemap_loci` steps.

``` r

colmap <- echodata::construct_colmap(
      CHR = "chr",
      POS = "pos",
      N = "qtl.N",
      SNP = "snp",
      P = "qtl.pvalues",
      Effect = "qtl.beta",
      StdErr = "qtl.varbeta",
      MAF = "qtl.MAF",
      Locus = "Locus",
      Gene = "gene")
```

### Prepare `top_SNPs` data.frame

- In this case, we don’t have a top SNPs file ready. So we’re just going
  to make one directly from the full summary stats file itself (*NOTE*:
  You can only use this approach if you can fit the entire file in
  memory).
- In this case, you’ll want to make sure to set
  `grouping_vars=c("Locus","Gene")` so that you get top SNPs for each
  eGene-locus pair (not just one SNP per locus).

``` r

topSNPs <- echodata::import_topSNPs(
  topSS = coloc_res$path,
  colmap = colmap,
  ## Important for QTLs: group by both Locus and Gene
  grouping_vars = c("Locus","Gene"))
```

    ## Loading required namespace: MungeSumstats

    ## Renaming column: chr ==> CHR

    ## Renaming column: pos ==> POS

    ## Renaming column: snp ==> SNP

    ## Renaming column: qtl.pvalues ==> P

    ## Renaming column: qtl.beta ==> Effect

    ## Renaming column: qtl.varbeta ==> StdErr

    ## Renaming column: qtl.MAF ==> MAF

    ## Renaming column: gene ==> Gene

    ## Renaming column: qtl.N ==> N

    ## [1] "+ Assigning Gene and Locus independently."

    ## Standardising column headers.

    ## First line of summary statistics file:

    ## Locus    SNP V.df1   z.df1   r.df1   lABF.df1    V.df2   z.df2   r.df2   lABF.df2    internal.sum.lABF   SNP.PP.H4   gwas.pvalues    gwas.beta   gwas.varbeta    gwas.MAF    CHR POS A1  A2  gwas.N  gwas.type   s   Gene    P   Effect  StdErr  MAF QTL_chr QTL_pos N   qtl.type    N_cases N_controls  

    ## Returning unmapped column names without making them uppercase.

    ## + Mapping colnames from MungeSumstats ==> echolocatoR

``` r

head(topSNPs)
```

    ## Key: <Locus>
    ##                     Locus        SNP      V.df1      z.df1     r.df1   lABF.df1
    ##                    <char>     <char>      <num>      <num>     <num>      <num>
    ## 1:  ABCA7_ENSG00000160953 rs76951864 0.00227529 -1.2138365 0.9461792 -0.7639979
    ## 2:   BIN1_ENSG00000136731  rs4663110 0.00021904  2.0405405 0.9945538 -0.5358563
    ## 3:    CR1_ENSG00000266094  rs1518110 0.00030276 -1.0977011 0.9924879 -1.8476694
    ## 4: INPP5D_ENSG00000168918  rs1881492 0.00036100  0.1631579 0.9910557 -2.3451794
    ## 5: MS4A6A_ENSG00000149476  rs6591611 0.00028900 -1.4176471 0.9928268 -1.4710500
    ## 6:  PILRA_ENSG00000106366  rs2227631 0.00020736 -0.1041667 0.9948427 -2.6282771
    ##          V.df2     z.df2     r.df2   lABF.df2 internal.sum.lABF    SNP.PP.H4
    ##          <num>     <num>     <num>      <num>             <num>        <num>
    ## 1: 0.007802613 -5.826752 0.2374864  3.8958871         3.1318892 6.231637e-08
    ## 2: 0.004149176  2.778507 0.4321435  1.3851521         0.8492958 1.013056e-25
    ## 3: 0.002913938 -3.976592 0.4592553  3.3237628         1.4760934 4.495411e-11
    ## 4: 0.012164546 -2.654756 0.1470877  0.4387679        -1.9064116 3.587352e-08
    ## 5: 0.015193274  2.675112 0.2112492  0.6372210        -0.8338290 8.201862e-13
    ## 6: 0.006858275  8.459513 0.5879401 20.5941901        17.9659130 9.907667e-01
    ##    gwas.pvalues gwas.beta gwas.varbeta gwas.MAF   CHR       POS     A1     A2
    ##           <num>     <num>        <num>    <num> <num>     <num> <char> <char>
    ## 1:      0.22500   -0.0579   0.00227529   0.1054    19    374604      A      G
    ## 2:      0.04153    0.0302   0.00021904   0.6262     2 127922463      T      C
    ## 3:      0.27270   -0.0191   0.00030276   0.7803     1 206944861      A      C
    ## 4:      0.87100    0.0031   0.00036100   0.7763     2 233406998      T      G
    ## 5:      0.15760   -0.0241   0.00028900   0.2207    11  60405172      A      C
    ## 6:      0.91860   -0.0015   0.00020736   0.4056     7 100769538      A      G
    ##    gwas.N gwas.type     s            Gene           P    Effect      StdErr
    ##     <int>    <char> <num>          <char>       <num>     <num>       <num>
    ## 1:  94437        cc  0.37 ENSG00000160953 8.31563e-09 -0.514691 0.007802613
    ## 2:  94437        cc  0.37 ENSG00000136731 7.93573e-05  0.178975 0.004149176
    ## 3:  94437        cc  0.37 ENSG00000266094 1.00516e-04 -0.214660 0.002913938
    ## 4:  94437        cc  0.37 ENSG00000168918 6.42575e-06 -0.292801 0.012164546
    ## 5:  94437        cc  0.37 ENSG00000149476 2.20256e-04  0.329737 0.015193274
    ## 6:  94437        cc  0.37 ENSG00000106366 5.45534e-17  0.700572 0.006858275
    ##       Freq QTL_chr   QTL_pos     N qtl.type N_cases N_controls
    ##      <num>  <char>     <int> <int>   <char>   <int>      <int>
    ## 1: 0.08770   chr19    374604    90    quant   34941      59496
    ## 2: 0.34385    chr2 127164887    90    quant   34941      59496
    ## 3: 0.25363    chr1 206771516    90    quant   34941      59496
    ## 4: 0.14447    chr2 232542288    90    quant   34941      59496
    ## 5: 0.18503   chr11  60637699    90    quant   34941      59496
    ## 6: 0.42796    chr7 101126257    90    quant   34941      59496

### Run fine-mapping pipeline

The following call runs the full fine-mapping pipeline on QTL data. This
is computationally expensive and may take **10–30 minutes per locus**.

``` r

res <- echolocatoR::finemap_loci(fullSS_path = coloc_res$path,
                                 topSNPs = topSNPs,
                                 ## Let's just fine-map 1 locus for demo purposes
                                 loci = topSNPs$Locus[1],
                                 dataset_name = "Kunkle_2019.microgliaQTL",
                                 dataset_type = "QTL",
                                 bp_distance = 1000,
                                 colmap = colmap,
                                 show_plot = TRUE,
                                 finemap_methods = c("ABF","FINEMAP","SUSIE") )
```

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
## [155] crayon_1.5.3                aws.signature_0.6.0        
## [157] ieugwasr_1.1.0              plyr_1.8.9                 
## [159] forcats_1.0.1               fs_1.6.7                   
## [161] stringi_1.8.7               coloc_5.2.3                
## [163] echoannot_1.0.1             viridisLite_0.4.3          
## [165] BiocParallel_1.45.0         Biostrings_2.79.5          
## [167] lazyeval_0.2.2              Matrix_1.7-4               
## [169] downloadR_1.0.0             echoplot_1.0.0             
## [171] dir.expiry_1.19.0           MungeSumstats_1.19.5       
## [173] BSgenome_1.79.1             patchwork_1.3.2            
## [175] hms_1.1.4                   bit64_4.6.0-1              
## [177] ggplot2_4.0.2               KEGGREST_1.51.1            
## [179] SummarizedExperiment_1.41.1 haven_2.5.5                
## [181] memoise_2.0.1               snpStats_1.61.1            
## [183] bslib_0.10.0                bit_4.6.0                  
## [185] readxl_1.4.5
```
