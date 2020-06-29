<center><h1> )    )   )  ) ))) 🦇 echolocatoR  🦇 ((( (  (   (    ( </h1></center>  

__*echolocatoR*__ is an R package for automating statistical and functional fine-mapping with extensive API access to datasets.  


## Documentation

### [Full documentation website](https://rajlabmssm.github.io/echolocatoR/)  

### [Full pipeline vignette](https://rajlabmssm.github.io/echolocatoR/articles/full_pipeline_vignette.html)  

### [Plotting vignette](https://rajlabmssm.github.io/echolocatoR/articles/plotting_vignette.html)


## Workflow  

![echoFlow](./images/PD_Finemapping_Flowchart_plus.png)  


## Quick installation  

In R:  
```R
if(!"devtools" %in% installed.packages()){install.packages("devtools")}
devtools::install_github("RajLabMSSM/echolocatoR")
```  
*NOTE*: While this GitHub repo is still private, you need to use a token to install __*echolocatoR*__  using the `auth_token` argument (see [here for details](http://devtools.r-lib.org/reference/install_github)).  



## Robust installation (co*conda*
As with most softwares, installation is half the battle.
The easiest way to install all of __*echolocatoR*__'s dependencies (inwhich include, Python, and command line)  toolsand make sure they play well together
is to create a [*conda*](https://docs.conda.io/en/latest/) environment.

1. If you haven't done so already, install [*conda*](https://docs.conda.io/en/latest/).  

2. Download the echoR.yml file found [here](https://github.com/RajLabMSSM/echolocatoR/blob/master/inst/conda/echoR.yml) (this file tells *conda* what to install).

3. Create the env from the *.yml* file,   
by the entering the following into the command line:  
```
conda env create -f <path_to_file>/echoR.yml
```

4. Activate the new env:  
```
conda activate echoR
```

5. In R, install __*echolocatoR*__:  
```R
if(!"devtools" %in% installed.packages()){install.packages("devtools")}
devtools::install_github("RajLabMSSM/echolocatoR")
```

To make sure __*echolocatoR*__ uses the packages in this env, 
you can then supply the env name to the `finemap_loci()` function using `conda_env="echoR"`.


### Axel

*Axel* is a command line tool that drastically accelerates large file downloads 
by using multiple cores. This is very useful for the large UK Bibank pre-computed LD files.
Unfortunately, it does require a bit of extra setup to use, but you can find the [instructions here](https://github.com/axel-download-accelerator/axel).

Once you've successfully installed *Axel*, you can tell `finemap_loci()` to use it with `download_method = "axel"`.  

  

## Fine-mapping Tools  

### Currently implemented  
- [ABF](https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html)  
- [SuSiE](https://github.com/stephenslab/susieR)  
- [FINEMAP](http://www.christianbenner.com)  
- [PolyFun](https://github.com/omerwe/polyfun)
- [GCTA-COJO](https://cnsgenomics.com/software/gcta/#COJO)
- [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0)  
- [coloc](https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html)


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
  echolocatoR::NOTT_2019.bigwig_metadata
  # Or 
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
- Data from this preprint contains results from single-cell chromatin accessibility epigenomic assays in from 39 human brains. 
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
- Inlcudes 187-annotation baseline model from [Gazal et al. 2018](https://www.nature.com/articles/s41588-018-0231-8).  
- You can alternatively supply a custom annotations matrix.

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
![Sinai](./inst/images/sinai.png)
