<center><h1> )    )   )  ) ))) 🦇 echolocatoR  🦇 ((( (  (   (    ( </h1></center>  

### Automated statistical and functional fine-mapping with extensive access to genome-wide datasets.   

<hr>  

### When using __*echolocatoR*__, please cite:  
[Brian M. Schilder, Jack Humphrey, Towfique Raj (2020) *echolocatoR: an automated end-to-end statistical and functional genomic fine-mapping pipeline*
bioRxiv 2020.10.22.351221; doi: https://doi.org/10.1101/2020.10.22.351221](https://www.biorxiv.org/content/10.1101/2020.10.22.351221v1)  

<hr>

Fine-mapping methods are a powerful means of identifying causal variants underlying a given phenotype, but are underutilized due to the technical challenges of implementation. __*echolocatoR*__ is an R package that automates end-to-end genomics fine-mapping, annotation, and plotting in order to identify the most probable causal variants associated with a given phenotype.

It requires minimal input from users (a GWAS or QTL summary statistics file), and includes a suite of statistical and functional fine-mapping tools. It also includes extensive access to datasets (linkage disequilibrium panels, epigenomic and genome-wide annotations, QTL).

The elimination of data gathering and preprocessing steps enables rapid fine-mapping of many loci in any phenotype, complete with locus-specific publication-ready figure generation. All results are merged into a single per-SNP summary file for additional downstream analysis and results sharing. Therefore __*echolocatoR*__ drastically reduces the barriers to identifying causal variants by making the entire fine-mapping pipeline rapid, robust and scalable.  


## Documentation

### [Documentation website](https://rajlabmssm.github.io/echolocatoR/)  

### [Full pipeline vignette](https://rajlabmssm.github.io/echolocatoR/articles/full_pipeline_vignette.html)  

### [Plotting vignette](https://rajlabmssm.github.io/echolocatoR/articles/plotting_vignette.html)


<br>


## Workflow  

![echoFlow](./images/echolocatoR_Fig1.png)  


<br>


## Quick installation  

In R:  
```R
if(!"devtools" %in% installed.packages()){install.packages("devtools")}
devtools::install_github("RajLabMSSM/echolocatoR")
```   

## Robust installation (*conda*)

  As with most softwares, installation is half the battle.
The easiest way to install all of __*echolocatoR*__'s dependencies (which include R, Python, and command line tools) and make sure they play well together
is to create a [*conda*](https://docs.conda.io/en/latest/) environment.

1. If you haven't done so already, install [*conda*](https://docs.conda.io/en/latest/).  

2. Download the *echoR.yml* file found [here](https://github.com/RajLabMSSM/echolocatoR/blob/master/inst/conda/echoR.yml) (this file tells *conda* what to install).

3. In command line, create the env from the *.yml* file:  
```
conda env create -f <path_to_file>/echoR.yml
```

4. Activate the new env:  
```
conda activate echoR
```

5. Open Rstudio from the command line interface (not by clicking the Rstudio icon). 
This helps to ensure Rstudio can find the paths to the packages in the conda env.

6. In R, install __*echolocatoR*__:  
```R
if(!"devtools" %in% installed.packages()){install.packages("devtools")}
devtools::install_github("RajLabMSSM/echolocatoR")
```

To make sure __*echolocatoR*__ uses the packages in this env (esp. if using from RStudio), you can then supply the env name to the `finemap_loci()` function using `conda_env="echoR"`.


<br>


## Dependencies   

For a full list of suggested packages, see [DESCRIPTION](https://github.com/RajLabMSSM/echolocatoR/blob/master/DESCRIPTION).  

\* = _optional_  

### R  
```
- magrittr  
- R.utils  
- dplyr  
- BiocManager 
- tidyverse
- knitr
- rmarkdown  
- data.table  
- foreign  
- reticulate  
- ggplot2    
- ggrepel  
- coloc    
- RColorBrewer   
- patchwork   
- htmltools  
- stringr    
- openxlsx  
- EnsDb.Hsapiens.v75    
- ensembldb   
- ggbio    
- BSgenome  
- Ckmeans.1d.dp  
- refGenome   
```

### Python  
```
- python>=3.6.1  
- pandas>=0.25.0   
- pandas-plink  
- pyarrow  
- fastparquet  
- scipy  
- scikit-learn  
- tqdm  
- bitarray  
- networkx  
- rpy2  
- requests  
```

### Command line   

#### [Tabix](http://www.htslib.org/doc/tabix.html)  
  + Rapid querying of summary stats files.   
  + To use it, specify `query_by="tabix"` in `finemap_loci()`.   

#### [bcftools](http://samtools.github.io/bcftools/bcftools.html)  
  + Used here for filtering populations in vcf files.  
  
#### [Axel](https://github.com/axel-download-accelerator/axel)  *    
  + Rapid multi-core downloading of large files (e.g. LD matrices from UK Biobank).  
  + To use it, specify `download_method="axel"` in `finemap_loci()`.  
  + For more info on installing/using *axel* in general, see this [tutorial](https://www.tecmint.com/axel-commandline-download-accelerator-for-linux/). Depending on what kind of computer you're using, this process will look a bit different:
    - **Mac**: Install [brew](https://brew.sh/), then: `brew install axel`  
    - **CentOS/RHEL 7**: `yum install epel-release; yum install axel`  
    - **Fedora**: `yum install axel; dnf install axel`  
    - **Debian Jessie (e.g. Ubuntu, Linux Mint)**: `aptitude install axel`  
   
 
 
<br>


## Fine-mapping Tools  

__*echolocatoR*__ will automatically check whether you have the necessary columns 
to run each tool you selected in `finemap_loci(finemap_methods=...)`. 
It will remove any tools that for which there are missing necessary columns, 
and produces a message letting you know which columns are missing.
Note that some columns (e.g. `MAF`,`N`,`t-stat`) can be automatically inferred if missing.  
For easy reference, we list the necessary columns here as well.   
See `?finemap_loci()` for descriptions of these columns.  
All methods require the columns: `SNP`,`CHR`,`POS`,`Effect`,`StdErr`    

Additional required columns: 
### [ABF](https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html): `proportion_cases`,`MAF`  
### [FINEMAP](http://www.christianbenner.com):`A1`,`A2`,`MAF`,`N`  
### [SuSiE](https://github.com/stephenslab/susieR): `N`  
### [PolyFun](https://github.com/omerwe/polyfun): `A1`,`A2`,`P`,`N`   
### [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0): `A1`,`A2`,`t-stat`  
### [GCTA-COJO](https://cnsgenomics.com/software/gcta/#COJO): `A1`,`A2`,`Freq`,`P`,`N`  
### [coloc](https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html): `N`,`MAF`  


<br>


## Datasets

For more detailed information about each dataset, use `?`:   
  ```R
  library(echolocatoR)
  ?NOTT_2019.interactome # example dataset
  ```

### Epigenomic & Genome-wide Annotations

#### [Nott et al. (2019)](https://science.sciencemag.org/content/366/6469/1134.abstract)
- Data from this publication contains results from cell type-specific (neurons, oligodendrocytes, astrocytes, microglia, & peripheral myeloid cells) epigenomic assays (H3K27ac, ATAC, H3K4me3) from human brain tissue.  

- For detailed metadata, see:
  ```R
  data("NOTT_2019.bigwig_metadata")
  ```  
- Built-in datasets:  
  + Enhancer/promoter coordinates (as *GenomicRanges*)   
  ```R
  data("NOTT_2019.interactome")
  # Examples of the data nested in "NOTT_2019.interactome" object:
  NOTT_2019.interactome$`Neuronal promoters`
  NOTT_2019.interactome$`Neuronal enhancers`
  NOTT_2019.interactome$`Microglia promoters`
  NOTT_2019.interactome$`Microglia enhancers`
  ...
  ...
  ```
  + PLAC-seq enhancer-promoter interactome coordinates   
  ```R
  NOTT_2019.interactome$H3K4me3_around_TSS_annotated_pe
  NOTT_2019.interactome$`Microglia interactome`
  NOTT_2019.interactome$`Neuronal interactome`
  NOTT_2019.interactome$`Oligo interactome`
  ...
  ...
  ```
- API access to full bigWig files on UCSC Genome Browser, which includes  
  + Epigenomic reads (as *GenomicRanges*)  
  + Aggregate epigenomic *score* for each cell type - assay combination     
  
#### [Corces et al. (2020)](https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1)  
- Data from this preprint contains results from bulk and single-cell chromatin accessibility epigenomic assays in 39 human brains. 
  ```R
  data("CORCES_2020.bulkATACseq_peaks")
  data("CORCES_2020.cicero_coaccessibility")
  data("CORCES_2020.HiChIP_FitHiChIP_loop_calls")
  data("CORCES_2020.scATACseq_celltype_peaks")
  data("CORCES_2020.scATACseq_peaks")
  ```
  
#### [XGR](http://xgr.r-forge.r-project.org)    
- API access to a diverse library of cell type/line-specific epigenomic (e.g. ENCODE) and other genome-wide annotations.    

#### [Roadmap](http://www.roadmapepigenomics.org)  
- API access to cell type-specific epigenomic data.  

#### [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)  
- API access to various genome-wide SNP annotations (e.g. missense, nonsynonmous, intronic, enhancer).  

#### [HaploR](https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html)  
- API access to known per-SNP QTL and epigenomic data hits.  

### QTLs

#### [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/)  
- API access to full summary statistics from many standardized e/s/t-QTL datasets.  
- Data access and colocalization tests facilitated through the [catalogueR](https://github.com/RajLabMSSM/catalogueR) R package.  

<br>


## Enrichment Tools

### [XGR](http://xgr.r-forge.r-project.org)   
- Binomial enrichment tests between customisable foreground and background SNPs.  

### [GoShifter](https://github.com/immunogenomics/goshifter)  
- LD-informed iterative enrichment analysis.

### [S-LDSC](https://www.nature.com/articles/ng.3954)
- Genome-wide stratified LD score regression.
- Inlccles 187-annotation baseline model from [Gazal et al. 2018](https://www.nature.com/articles/s41588-018-0231-8).  
- You can alternatively supply a custom annotations matrix.

### [motifbreakR](https://github.com/Simon-Coetzee/motifBreakR)  
- Identification of transcript factor binding motifs (TFBM) and prediction of SNP disruption to said motifs. 
- Includes a comprehensive list of TFBM databases via [MotifDB](https://bioconductor.org/packages/release/bioc/html/MotifDb.html) (9,900+ annotated position frequency matrices from 14 public sources, for multiple organisms).  

### [GARFIELD](https://www.bioconductor.org/packages/release/bioc/html/garfield.html) (**under construction**)
- Genomic enrichment with LD-informed heuristics.   


<br>


## LD Reference Panels  

### [UK Biobank](https://www.ukbiobank.ac.uk)
### [1000 Genomes Phase 1](https://www.internationalgenome.org)  
### [1000 Genomes Phase 3](https://www.internationalgenome.org)  


<hr><hr>

## Author

<a href="https://bschilder.github.io/BMSchilder/" target="_blank">Brian M. Schilder, Bioinformatician II</a>  
<a href="https://rajlab.org" target="_blank">Raj Lab</a>  
<a href="https://icahn.mssm.edu/about/departments/neuroscience" target="_blank">Department of Neuroscience, Icahn School of Medicine at Mount Sinai</a>  
![Sinai](./images/sinai.png)
