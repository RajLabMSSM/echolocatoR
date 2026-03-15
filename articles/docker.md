# Docker/Singularity Containers

## DockerHub

echolocator is now available via
[DockerHub](https://hub.docker.com/repository/docker/bschilder/echolocator)
as a containerised environment with Rstudio and all necessary
dependencies pre-installed.

### Installation

### Method 1: via Docker

First, [install Docker](https://docs.docker.com/get-docker/) if you have
not already.

Create an image of the [Docker](https://www.docker.com/) container in
command line:

    docker pull bschilder/echolocator

Once the image has been created, you can launch it with:

    docker run \
      -d \
      -e ROOT=true \
      -e PASSWORD="<your_password>" \
      -v ~/Desktop:/Desktop \
      -v /Volumes:/Volumes \
      -p 8787:8787 \
      bschilder/echolocator

#### NOTES

- Make sure to replace `<your_password>` above with whatever you want
  your password to be.\
- Change the paths supplied to the `-v` flags for your particular use
  case.
- The `-d` ensures the container will run in “detached” mode, which
  means it will persist even after you’ve closed your command line
  session.\
- The username will be *“rstudio”* by default.\
- Optionally, you can also install the [Docker
  Desktop](https://www.docker.com/products/docker-desktop) to easily
  manage your containers.

### Method 2: via Singularity

If you are using a system that does not allow Docker (as is the case for
many institutional computing clusters), you can instead [install Docker
images via
Singularity](https://sylabs.io/guides/2.6/user-guide/singularity_and_docker.html).

    singularity pull docker://bschilder/echolocator

### Usage

Finally, launch the containerised Rstudio by entering the following URL
in any web browser: *<http://localhost:8787/>*

Login using the credentials set during the Installation steps.

## Session Info

``` r

utils::sessionInfo()
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
    ## [157] plyr_1.8.9                  forcats_1.0.1              
    ## [159] fs_1.6.7                    stringi_1.8.7              
    ## [161] coloc_5.2.3                 echoannot_1.0.1            
    ## [163] viridisLite_0.4.3           BiocParallel_1.45.0        
    ## [165] Biostrings_2.79.5           lazyeval_0.2.2             
    ## [167] Matrix_1.7-4                downloadR_1.0.0            
    ## [169] echoplot_1.0.0              dir.expiry_1.19.0          
    ## [171] BSgenome_1.79.1             patchwork_1.3.2            
    ## [173] hms_1.1.4                   bit64_4.6.0-1              
    ## [175] ggplot2_4.0.2               KEGGREST_1.51.1            
    ## [177] SummarizedExperiment_1.41.1 haven_2.5.5                
    ## [179] memoise_2.0.1               snpStats_1.61.1            
    ## [181] bslib_0.10.0                bit_4.6.0                  
    ## [183] readxl_1.4.5

\
