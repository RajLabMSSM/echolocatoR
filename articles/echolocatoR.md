# Get started

``` r

library(echolocatoR)
```

    ## Registered S3 method overwritten by 'bit64':
    ##   method          from 
    ##   print.bitstring tools

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

## Full pipeline

All examples below use data from the Parkinson’s disease GWAS by Nalls
et al. (2019).

### Prepare `top_SNPs` data.frame

- To enable rapid fine-mapping of many loci, you can create a `top_SNPs`
  data.frame which contains the position of the lead/index SNP within
  each locus.
- [`finemap_loci()`](https://rajlabmssm.github.io/echolocatoR/reference/finemap_loci.md)
  (see next step) will then use this info to extract subsets of the full
  GWAS/QTL summary statistics using windows centered on each lead/index
  SNP.
- The `topSS` argument can either be a data.frame, or a path to a topSS
  file saved somewhere. Most common tabular data formats (e.g. .tsv,
  .csv, .xlsx) are accepted.

``` r

#### Load example top SNPs (pre-formatted) ####
topSS <- echodata::topSNPs_Nalls2019_raw
#### construct a column mapping object ####
colmap <- echodata::construct_colmap(P = "P, all studies",
                                     Effect = "Beta, all studies",
                                     Locus = "Nearest Gene",
                                     Gene = "QTL Nominated Gene (nearest QTL)")
#### Import top SNPs ####
topSNPs <- echodata::import_topSNPs(
    topSS = echodata::topSNPs_Nalls2019_raw,
    colmap = colmap,
    grouping_vars = "Locus Number")
```

    ## Loading required namespace: MungeSumstats

    ## Renaming column: P, all studies ==> P

    ## Renaming column: Beta, all studies ==> Effect

    ## Renaming column: Nearest Gene ==> Locus

    ## Renaming column: QTL Nominated Gene (nearest QTL) ==> Gene

    ## [1] "+ Assigning Gene and Locus independently."

    ## Standardising column headers.

    ## First line of summary statistics file:

    ## SNP  CHR BP  Locus   Gene    Effect allele   Other allele    Effect allele frequency Effect  SE, all studies P   P, COJO, all studies    P, random effects, all studies  P, Conditional 23AndMe only P, 23AndMe only I2, all studies Freq1, previous studies Beta, previous studies  StdErr, previous studies    P, previous studies I2, previous studies    Freq1, new studies  Beta, new studies   StdErr, new studies P, new studies  I2, new studies Passes pooled 23andMe QC    Known GWAS locus within 1MB Failed final filtering and QC   Locus within 250KB  Locus Number    

    ## Returning unmapped column names without making them uppercase.

    ## + Mapping colnames from MungeSumstats ==> echolocatoR

``` r

head(topSNPs)
```

    ## Key: <Locus>
    ##           SNP   CHR       POS   Locus    Gene     A2     A1   Freq  Effect
    ##        <char> <num>     <num>  <char>  <char> <char> <char>  <num>   <num>
    ## 1:  rs1941685    18  31304318   ASXL3    <NA>      t      g 0.4983  0.0531
    ## 2:  rs2280104     8  22525980    BIN3    BIN3      t      c 0.3604  0.0556
    ## 3: rs61169879    17  59917366   BRIP1   MED13      t      c 0.1641  0.0820
    ## 4:  rs4698412     4  15737348    BST1    BST1      a      g 0.5529  0.1035
    ## 5: rs11950533     5 134199105 C5orf24 TXNDC15      a      c 0.1020 -0.0916
    ## 6:  rs9568188    13  49927732  CAB39L  CAB39L      t      c 0.7397  0.0617
    ##    SE, all studies        P P, COJO, all studies P, random effects, all studies
    ##              <num>    <num>                <num>                          <num>
    ## 1:          0.0094 1.69e-08             1.61e-08                      1.690e-08
    ## 2:          0.0098 1.16e-08             1.40e-08                      1.860e-06
    ## 3:          0.0134 9.28e-10             9.40e-10                      6.210e-06
    ## 4:          0.0094 2.06e-28             9.73e-36                      1.680e-19
    ## 5:          0.0158 7.16e-09             6.73e-09                      2.680e-08
    ## 6:          0.0108 1.15e-08             1.11e-08                      2.458e-04
    ##    P, Conditional 23AndMe only P, 23AndMe only I2, all studies
    ##                          <num>           <num>           <num>
    ## 1:                    1.64e-08        1.60e-08             0.0
    ## 2:                    4.35e-02        4.94e-02             8.9
    ## 3:                    9.07e-07        2.05e-06            16.4
    ## 4:                    1.05e-07        1.15e-07            13.9
    ## 5:                    5.08e-04        4.22e-04             1.9
    ## 6:                    4.29e-06        4.41e-06            21.4
    ##    Freq1, previous studies Beta, previous studies StdErr, previous studies
    ##                      <num>                  <num>                    <num>
    ## 1:                  0.4978                 0.0507                   0.0113
    ## 2:                  0.3595                 0.0647                   0.0117
    ## 3:                  0.1641                 0.0849                   0.0163
    ## 4:                  0.5547                 0.1023                   0.0112
    ## 5:                  0.1008                -0.0986                   0.0190
    ## 6:                  0.7405                 0.0657                   0.0130
    ##    P, previous studies I2, previous studies Freq1, new studies
    ##                  <num>                <num>              <num>
    ## 1:            6.82e-06                  0.0             0.4995
    ## 2:            3.23e-08                 15.0             0.3626
    ## 3:            1.90e-07                  0.0             0.1641
    ## 4:            7.32e-20                 49.3             0.5485
    ## 5:            2.05e-07                 43.8             0.1048
    ## 6:            4.00e-07                 39.8             0.7380
    ##    Beta, new studies StdErr, new studies P, new studies I2, new studies
    ##                <num>               <num>          <num>           <num>
    ## 1:            0.0586              0.0171      6.147e-04            14.3
    ## 2:            0.0350              0.0177      4.724e-02             1.9
    ## 3:            0.0761              0.0236      1.242e-03            23.7
    ## 4:            0.1062              0.0170      4.160e-10            10.9
    ## 5:           -0.0756              0.0287      8.388e-03             0.0
    ## 6:            0.0526              0.0196      7.353e-03            22.3
    ##    Passes pooled 23andMe QC Known GWAS locus within 1MB
    ##                      <char>                       <num>
    ## 1:                        T                           0
    ## 2:                        T                           1
    ## 3:                        T                           0
    ## 4:                        T                           1
    ## 5:                        T                           0
    ## 6:                        T                           0
    ##    Failed final filtering and QC Locus within 250KB Locus Number
    ##                            <num>             <char>       <char>
    ## 1:                             0                  0           73
    ## 2:                             0                  0           39
    ## 3:                             0                  0           71
    ## 4:                             0                  0           20
    ## 5:                             0                  0           28
    ## 6:                             0                  0           53

### Path to full summary stats file

- Since a full GWAS summary stats file would be too large to include
  within *echolocatoR*, we instead provide an example subset of the full
  summary stats.

- To simulate how you’d actually use your own full summary stats file,
  we will save our example dataset to your computer (you can change the
  path to wherever you like).

- We highly recommend munging your full summary stats using the
  Bioconductor package
  [`MungeSumstats`](https://github.com/neurogenomics/MungeSumstats)
  first. It’s easy to use and very robust. It also means you don’t have
  to provide most column mapping arguments in `finemap_loci` when
  `munged=TRUE`.

Here’s an example of how to munge your full summary stats file:

    fullSS_path <- echodata::example_fullSS(munged = FALSE)
    fullSS_path <- MungeSumstats::format_sumstats(path = fullSS_path, ref_genome = "GRCH37")

We have already munged the following example summary stats for you.

``` r

fullSS_path <- echodata::example_fullSS(dataset = "Nalls2019")
```

    ## Writing file to ==> /tmp/RtmpIaFbFB/nalls2019.fullSS_subset.tsv

### Run fine-mapping pipeline

For a full description of all arguments, see
[`?finemap_loci`](https://rajlabmssm.github.io/echolocatoR/reference/finemap_loci.md).

Here are some key arguments:

- *results_dir*: Where you want to store all of your results.
- *finemap_methods*: Which fine-mapping methods you want to run. For a
  full list of currently supported methods, run the function
  [`echofinemap::lfm()`](https://rdrr.io/pkg/echofinemap/man/lfm.html).
- *bp_distance*: Controls window size. Specifically, `bp_distance` is
  the number of basepairs upstream/downstream you want to extract for
  each locus. For example, if you want a 2Mb window (+/- 1Mb from the
  lead/index SNP in `top_SNPs`), set `bp_distance=1e+06`.
- *plot_zoom*: Zoom in/out from the center of each locus when producing
  the multiview plot. You can adjust this separately from `bp_distance`
  so that you don’t have rerun the whole pipeline each time (locus
  subsets, LD matrices, and fine-mapping results are all automatically
  saved in locus-specific folders).

**Note**: Please use the full absolute paths (instead of relative paths)
wherever possible (e.g. `results_dir`). This is especially important for
the tool *FINEMAP*.

The following call fine-maps two loci (BST1 and MEX3C) using three
statistical fine-mapping methods. Depending on your machine and internet
speed, this may take **10-30 minutes per locus**.

``` r

results <- echolocatoR::finemap_loci(
 fullSS_path = fullSS_path,
 topSNPs = topSNPs,
 loci = c("BST1","MEX3C"),
 LD_reference = "1KGphase3",
 dataset_name = "Nalls23andMe_2019",
 fullSS_genome_build = "hg19",
 bp_distance = 1000,
 finemap_methods = c("ABF","SUSIE","FINEMAP"),
 munged = TRUE)
```

### Inspect results

The output of
[`finemap_loci()`](https://rajlabmssm.github.io/echolocatoR/reference/finemap_loci.md)
is a named list of data.tables, one per locus. Each data.table contains
the original GWAS summary statistics augmented with fine-mapping
posterior probabilities and credible set assignments.

Here we use the bundled BST1 example results to show what the output
looks like:

``` r

dat <- echodata::BST1
## Top fine-mapped SNPs (consensus across methods)
consensus <- subset(dat, Consensus_SNP == TRUE)
knitr::kable(
    consensus[, c("SNP","CHR","POS","P","mean.PP","Support","Consensus_SNP")],
    caption = "Consensus fine-mapped SNPs for the BST1 locus",
    row.names = FALSE
)
```

| SNP        | CHR |      POS |   P |   mean.PP | Support | Consensus_SNP |
|:-----------|----:|---------:|----:|----------:|--------:|:--------------|
| rs4541502  |   4 | 15712787 |   0 | 0.7500001 |       3 | TRUE          |
| rs34559912 |   4 | 15730146 |   0 | 0.5155283 |       2 | TRUE          |
| rs4389574  |   4 | 15730398 |   0 | 0.5000000 |       2 | TRUE          |

Consensus fine-mapped SNPs for the BST1 locus {.table}

## Next steps

- Explore and interpret results:
  [`vignette("explore_results")`](https://rajlabmssm.github.io/echolocatoR/articles/explore_results.md)
- Visualize loci:
  [`vignette("plotting")`](https://rajlabmssm.github.io/echolocatoR/articles/plotting.md)
- Summarise across loci:
  [`vignette("summarise")`](https://rajlabmssm.github.io/echolocatoR/articles/summarise.md)
- Fine-map QTL data:
  [`vignette("QTLs")`](https://rajlabmssm.github.io/echolocatoR/articles/QTLs.md)
- Learn about sub-packages:
  [`vignette("echoverse_modules")`](https://rajlabmssm.github.io/echolocatoR/articles/echoverse_modules.md)

------------------------------------------------------------------------

## Session info

``` r

utils::sessionInfo()
```

```
## R Under development (unstable) (2026-03-15 r89629)
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

\
