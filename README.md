<img src='https://github.com/RajLabMSSM/echolocatoR/raw/echoverse/inst/hex/hex.png' height='300'><br><br>
[![](https://img.shields.io/badge/devel%20version-2.0.0-black.svg)](https://github.com/RajLabMSSM/echolocatoR)
[![R build
status](https://github.com/RajLabMSSM/echolocatoR/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RajLabMSSM/echolocatoR/actions)
[![](https://img.shields.io/github/last-commit/RajLabMSSM/echolocatoR.svg)](https://github.com/RajLabMSSM/echolocatoR/commits/master)
[![](https://codecov.io/gh/RajLabMSSM/echolocatoR/branch/master/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/echolocatoR)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
[![](https://img.shields.io/badge/doi-10.1093/bioinformatics/btab658-blue.svg)](https://doi.org/10.1093/bioinformatics/btab658)
<h5>
Author: <i>Brian M. Schilder</i>
</h5>
<h5>
README updated: <i>Mar-07-2022</i>
</h5>

## `echolocatoR`: Automated statistical and functional fine-mapping with extensive access to genome-wide datasets.

### The *echoverse*

`echolocatoR` is part of the
[***echoverse***](https://github.com/topics/echoverse), a suite of R
packages designed to facilitate different steps in genetic fine-mapping.

`echolocatoR` calls each of these other packages (i.e. “modules”)
internally to create a unified pipeline. However, you can also use each
module independently to create your own custom workflows.

#### ***echoverse*** dependency graph

<img src="./images/echoverse.png" height="400px" style="border-radius: 20px;">

> Made with [`echodeps`](https://github.com/RajLabMSSM/echodeps), yet
> another ***echoverse*** module. See [here for the interactive
> version](https://rajlabmssm.github.io/Fine_Mapping/echolocatoR.dep_graph.html)
> with package descriptions and links to each GitHub repo.

### Citation

If you use `echolocatoR`, or any of the **echoverse** subpackages,
please cite:

> Brian M Schilder, Jack Humphrey, Towfique Raj (2021) echolocatoR: an
> automated end-to-end statistical and functional genomic fine-mapping
> pipeline, *Bioinformatics*; btab658,
> <https://doi.org/10.1093/bioinformatics/btab658>

## Introduction

Fine-mapping methods are a powerful means of identifying causal variants
underlying a given phenotype, but are underutilized due to the technical
challenges of implementation. ***echolocatoR*** is an R package that
automates end-to-end genomics fine-mapping, annotation, and plotting in
order to identify the most probable causal variants associated with a
given phenotype.

It requires minimal input from users (a GWAS or QTL summary statistics
file), and includes a suite of statistical and functional fine-mapping
tools. It also includes extensive access to datasets (linkage
disequilibrium panels, epigenomic and genome-wide annotations, QTL).

The elimination of data gathering and preprocessing steps enables rapid
fine-mapping of many loci in any phenotype, complete with locus-specific
publication-ready figure generation. All results are merged into a
single per-SNP summary file for additional downstream analysis and
results sharing. Therefore ***echolocatoR*** drastically reduces the
barriers to identifying causal variants by making the entire
fine-mapping pipeline rapid, robust and scalable.

<img src="./images/echolocatoR_Fig1.png" style="border-radius: 10px;">

## Installation

``` r
if(!require("remotes")) install.packages("remotes")

remotes::install_github("RajLabMSSM/echolocatoR@echoverse")
library(echolocatoR)
```

## Documentation

### [Website](https://rajlabmssm.github.io/echolocatoR)

### [Getting started](https://rajlabmssm.github.io/echolocatoR/articles/echolocatoR)

### Bugs/requests

Please report any bugs/requests on [GitHub
Issues](https://github.com/RajLabMSSM/echolocatoR/issues).

[Contributions](https://github.com/RajLabMSSM/echolocatoR/pulls) are
welcome!

## Literature

### For applications of ***echolocatoR*** in the literature, please see:

> 1.  E Navarro, E Udine, K de Paiva Lopes, M Parks, G Riboldi, BM
>     Schilder…T Raj (2020) Dysregulation of mitochondrial and
>     proteo-lysosomal genes in Parkinson’s disease myeloid cells.
>     Nature Genetics. <https://doi.org/10.1101/2020.07.20.212407>
> 2.  BM Schilder, T Raj (2021) Fine-Mapping of Parkinson’s Disease
>     Susceptibility Loci Identifies Putative Causal Variants. Human
>     Molecular Genetics, ddab294,
>     <https://doi.org/10.1093/hmg/ddab294>  
> 3.  K de Paiva Lopes, G JL Snijders, J Humphrey, A Allan, M Sneeboer,
>     E Navarro, BM Schilder…T Raj (2022) Genetic analysis of the human
>     microglial transcriptome across brain regions, aging and disease
>     pathologies. Nature Genetics,
>     <https://doi.org/10.1038/s41588-021-00976-y>

## `echofinemap`: Fine-mapping tools

***echolocatoR*** will automatically check whether you have the
necessary columns to run each tool you selected in
`echolocatoR::finemap_loci(finemap_methods=...)`. It will remove any
tools that for which there are missing necessary columns, and produces a
message letting you know which columns are missing. Note that some
columns (e.g. `MAF`,`N`,`t-stat`) can be automatically inferred if
missing.  
For easy reference, we list the necessary columns here as well.  
See `?finemap_loci()` for descriptions of these columns.  
All methods require the columns: `SNP`,`CHR`,`POS`,`Effect`,`StdErr`

Additional required columns:

``` r
knitr::kable(echofinemap::required_cols())
```

| method        | required                                                 | suggested        |
|:--------------|:---------------------------------------------------------|:-----------------|
| ABF           | SNP , CHR , POS , Effect , StdErr , N , proportion_cases | MAF              |
| COJO          | SNP , CHR , POS , Effect, StdErr, A1 , A2                | Freq, P , N      |
| COLOC         | SNP , CHR , POS , Effect, StdErr                         | NULL             |
| FINEMAP       | SNP , CHR , POS , Effect, StdErr                         | A1 , A2 , MAF, N |
| PAINTOR       | SNP , CHR , POS , Effect, StdErr                         | MAF              |
| POLYFUN_SUSIE | SNP , CHR , POS , Effect, StdErr, P , A1 , A2            | MAF, N           |
| SUSIE         | SNP , CHR , POS , Effect, StdErr                         | N                |

### [ABF](https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html)

### [FINEMAP](http://www.christianbenner.com)

### [SuSiE](https://github.com/stephenslab/susieR)

### [PolyFun](https://github.com/omerwe/polyfun)

### [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0)

### [GCTA-COJO](https://cnsgenomics.com/software/gcta/#COJO)

### [coloc](https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html)

<br>

## Multi-finemap results files

The main output of ***echolocatoR*** are the multi-finemap files (for
example, `echodata::BST1`). They are stored in the locus-specific
*Multi-finemap* subfolders.

### Column descriptions

-   **Standardized GWAS/QTL summary statistics**: e.g.
    `SNP`,`CHR`,`POS`,`Effect`,`StdErr`. See `?finemap_loci()` for
    descriptions of each.  
-   **leadSNP**: The designated proxy SNP per locus, which is the SNP
    with the smallest p-value by default.
-   **\<tool\>.CS**: The 95% probability Credible Set (CS) to which a
    SNP belongs within a given fine-mapping tool’s results. If a SNP is
    not in any of the tool’s CS, it is assigned `NA` (or `0` for the
    purposes of plotting).  
-   **\<tool\>.PP**: The posterior probability that a SNP is causal for
    a given GWAS/QTL trait.  
-   **Support**: The total number of fine-mapping tools that include the
    SNP in its CS.
-   **Consensus_SNP**: By default, defined as a SNP that is included in
    the CS of more than `N` fine-mapping tool(s), i.e. `Support>1`
    (default: `N=1`).  
-   **mean.PP**: The mean SNP-wise PP across all fine-mapping tools
    used.
-   **mean.CS**: If mean PP is greater than the 95% probability
    threshold (`mean.PP>0.95`) then `mean.CS` is 1, else 0. This tends
    to be a very stringent threshold as it requires a high degree of
    agreement between fine-mapping tools.

### Notes

-   Separate multi-finemap files are generated for each LD reference
    panel used, which is included in the file name (e.g.
    *UKB_LD.Multi-finemap.tsv.gz*).

-   Each fine-mapping tool defines its CS and PP slightly differently,
    so please refer to the associated original publications for the
    exact details of how these are calculated (links provided above).

## Datasets

Datasets are now stored/retrieved via the following **echoverse**
sub-packages.

For more detailed information about each dataset, use `?`:

``` r
library(echoannot)   
?NOTT_2019.interactome # epigenomic annotations

library(echodata) 
?BST1 # fine-mapping results 
```

### `echoannot`: Epigenomic & genome-wide annotations

#### [Nott et al. (2019)](https://science.sciencemag.org/content/366/6469/1134.abstract)

-   Data from this publication contains results from cell type-specific
    (neurons, oligodendrocytes, astrocytes, microglia, & peripheral
    myeloid cells) epigenomic assays (H3K27ac, ATAC, H3K4me3) from human
    brain tissue.

#### [XGR](http://xgr.r-forge.r-project.org)

-   API access to a diverse library of cell type/line-specific
    epigenomic (e.g. ENCODE) and other genome-wide annotations.

#### [Roadmap](http://www.roadmapepigenomics.org)

-   API access to cell type-specific epigenomic data.

#### [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)

-   API access to various genome-wide SNP annotations (e.g. missense,
    nonsynonmous, intronic, enhancer).

#### [HaploR](https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html)

-   API access to known per-SNP QTL and epigenomic data hits.

### `catalogueR`: QTLs

#### [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/)

-   API access to full summary statistics from many standardized
    e/s/t-QTL datasets.  
-   Data access and colocalization tests facilitated through the
    [catalogueR](https://github.com/RajLabMSSM/catalogueR) R package.

<br>

## Enrichment tools

### [XGR](http://xgr.r-forge.r-project.org)

-   Binomial enrichment tests between customisable foreground and
    background SNPs.

### [GoShifter](https://github.com/immunogenomics/goshifter)

-   LD-informed iterative enrichment analysis.

### [S-LDSC](https://www.nature.com/articles/ng.3954)

-   Genome-wide stratified LD score regression.
-   Inlccles 187-annotation baseline model from [Gazal et al.
    2018](https://www.nature.com/articles/s41588-018-0231-8).  
-   You can alternatively supply a custom annotations matrix.

### [motifbreakR](https://github.com/Simon-Coetzee/motifBreakR)

-   Identification of transcript factor binding motifs (TFBM) and
    prediction of SNP disruption to said motifs.
-   Includes a comprehensive list of TFBM databases via
    [MotifDB](https://bioconductor.org/packages/release/bioc/html/MotifDb.html)
    (9,900+ annotated position frequency matrices from 14 public
    sources, for multiple organisms).

### [GARFIELD](https://www.bioconductor.org/packages/release/bioc/html/garfield.html) (**under construction**)

-   Genomic enrichment with LD-informed heuristics.

<br>

## `echoLD`: LD reference panels

### [UK Biobank](https://www.ukbiobank.ac.uk)

### [1000 Genomes Phase 1](https://www.internationalgenome.org)

### [1000 Genomes Phase 3](https://www.internationalgenome.org)

<hr>

## Contact

<a href="https://bschilder.github.io/BMSchilder/" target="_blank">Brian
M. Schilder, Bioinformatician II</a>  
<a href="https://rajlab.org" target="_blank">Raj Lab</a>  
<a href="https://icahn.mssm.edu/about/departments/neuroscience" target="_blank">Department
of Neuroscience, Icahn School of Medicine at Mount Sinai</a>

<br>

# Session info

<details>

``` r
utils::sessionInfo()
```

```
## R version 4.1.0 (2021-05-18)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## loaded via a namespace (and not attached):
##   [1] BiocFileCache_2.2.1         coloc_5.1.1                
##   [3] plyr_1.8.6                  splines_4.1.0              
##   [5] BiocParallel_1.28.3         usethis_2.1.5              
##   [7] GenomeInfoDb_1.30.1         ggplot2_3.3.5              
##   [9] digest_0.6.29               yulab.utils_0.0.4          
##  [11] htmltools_0.5.2             viridis_0.6.2              
##  [13] echofinemap_0.99.0          fansi_1.0.2                
##  [15] magrittr_2.0.2              memoise_2.0.1              
##  [17] BSgenome_1.62.0             gert_1.5.0                 
##  [19] openxlsx_4.2.5              credentials_1.3.2          
##  [21] Biostrings_2.62.0           matrixStats_0.61.0         
##  [23] R.utils_2.11.0              askpass_1.1                
##  [25] prettyunits_1.1.1           colorspace_2.0-3           
##  [27] blob_1.2.2                  rappdirs_0.3.3             
##  [29] gitcreds_0.1.1              xfun_0.29                  
##  [31] dplyr_1.0.8                 crayon_1.5.0               
##  [33] RCurl_1.98-1.6              dlstats_0.1.4              
##  [35] echodata_0.99.6             jsonlite_1.8.0             
##  [37] survival_3.3-0              VariantAnnotation_1.40.0   
##  [39] glue_1.6.2                  gtable_0.3.0               
##  [41] zlibbioc_1.40.0             XVector_0.34.0             
##  [43] DelayedArray_0.20.0         BiocGenerics_0.40.0        
##  [45] scales_1.1.1                DBI_1.1.2                  
##  [47] Rcpp_1.0.8                  viridisLite_0.4.0          
##  [49] progress_1.2.2              reticulate_1.24-9000       
##  [51] bit_4.0.4                   clisymbols_1.2.0           
##  [53] stats4_4.1.0                DT_0.21                    
##  [55] htmlwidgets_1.5.4           httr_1.4.2                 
##  [57] badger_0.1.0                RColorBrewer_1.1-2         
##  [59] ellipsis_0.3.2              pkgconfig_2.0.3            
##  [61] reshape_0.8.8               XML_3.99-0.9               
##  [63] R.methodsS3_1.8.1           seqminer_8.4               
##  [65] dbplyr_2.1.1                utf8_1.2.2                 
##  [67] tidyselect_1.1.2            rlang_1.0.2                
##  [69] AnnotationDbi_1.56.2        munsell_0.5.0              
##  [71] tools_4.1.0                 cachem_1.0.6               
##  [73] cli_3.2.0                   generics_0.1.2             
##  [75] RSQLite_2.2.10              evaluate_0.15              
##  [77] stringr_1.4.0               fastmap_1.1.0              
##  [79] yaml_2.3.5                  sys_3.4                    
##  [81] knitr_1.37                  bit64_4.0.5                
##  [83] fs_1.5.2                    zip_2.2.0                  
##  [85] purrr_0.3.4                 KEGGREST_1.34.0            
##  [87] gh_1.3.0                    R.oo_1.24.0                
##  [89] xml2_1.3.3                  biomaRt_2.50.3             
##  [91] compiler_4.1.0              rstudioapi_0.13            
##  [93] filelock_1.0.2              curl_4.3.2                 
##  [95] susieR_0.11.92              png_0.1-7                  
##  [97] tibble_3.1.6                stringi_1.7.6              
##  [99] highr_0.9                   GenomicFeatures_1.46.5     
## [101] desc_1.4.0                  lattice_0.20-45            
## [103] Matrix_1.4-0                vctrs_0.3.8                
## [105] pillar_1.7.0                lifecycle_1.0.1            
## [107] BiocManager_1.30.16         snpStats_1.44.0            
## [109] data.table_1.14.2           bitops_1.0-7               
## [111] irlba_2.3.5                 rtracklayer_1.54.0         
## [113] GenomicRanges_1.46.1        R6_2.5.1                   
## [115] BiocIO_1.4.0                gridExtra_2.3              
## [117] IRanges_2.28.0              assertthat_0.2.1           
## [119] SummarizedExperiment_1.24.0 openssl_1.4.6              
## [121] rprojroot_2.0.2             rjson_0.2.21               
## [123] GenomicAlignments_1.30.0    Rsamtools_2.10.0           
## [125] S4Vectors_0.32.3            GenomeInfoDbData_1.2.7     
## [127] parallel_4.1.0              hms_1.1.1                  
## [129] grid_4.1.0                  tidyr_1.2.0                
## [131] rmarkdown_2.11              rvcheck_0.2.1              
## [133] MatrixGenerics_1.6.0        echotabix_0.99.2           
## [135] echoLD_0.99.1               mixsqp_0.3-43              
## [137] piggyback_0.1.1             Biobase_2.54.0             
## [139] lubridate_1.8.0             restfulr_0.0.13
```

</details>
