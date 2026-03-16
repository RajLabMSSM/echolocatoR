# QTLs

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

## Next steps

- Explore and interpret results:
  [`vignette("explore_results")`](https://rajlabmssm.github.io/echolocatoR/articles/explore_results.md)
- Visualize loci:
  [`vignette("plotting")`](https://rajlabmssm.github.io/echolocatoR/articles/plotting.md)
- Summarise across loci:
  [`vignette("summarise")`](https://rajlabmssm.github.io/echolocatoR/articles/summarise.md)
- GWAS fine-mapping (non-QTL):
  [`vignette("echolocatoR")`](https://rajlabmssm.github.io/echolocatoR/articles/echolocatoR.md)

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
## [157] aws.signature_0.6.0         ieugwasr_1.1.0             
## [159] plyr_1.8.9                  forcats_1.0.1              
## [161] fs_1.6.7                    stringi_1.8.7              
## [163] coloc_5.2.3                 echoannot_1.0.1            
## [165] viridisLite_0.4.3           BiocParallel_1.44.0        
## [167] Biostrings_2.78.0           lazyeval_0.2.2             
## [169] Matrix_1.7-4                downloadR_1.0.0            
## [171] echoplot_0.99.9             dir.expiry_1.18.0          
## [173] MungeSumstats_1.18.1        BSgenome_1.78.0            
## [175] hms_1.1.4                   patchwork_1.3.2            
## [177] bit64_4.6.0-1               ggplot2_4.0.2              
## [179] KEGGREST_1.50.0             SummarizedExperiment_1.40.0
## [181] haven_2.5.5                 memoise_2.0.1              
## [183] snpStats_1.60.0             bslib_0.10.0               
## [185] bit_4.6.0                   readxl_1.4.5
```
