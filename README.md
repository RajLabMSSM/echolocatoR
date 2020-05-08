# echolocatoR
Automated statistical and functional fine-mapping pipeline with extensive API access to datasets.


## Workflow  

![echoFlow](./inst/images/echolocatoR_flowchart.png)

## Fine-mapping Tools  

### Currently implemented  
- [ABF](https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html)  
- [susieR](https://github.com/stephenslab/susieR)  
- [FINEMAP](http://www.christianbenner.com)  
- [Polyfun+SusieR](https://github.com/omerwe/polyfun)
- [GCTA-COJO](https://cnsgenomics.com/software/gcta/#COJO)
- [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0)  
- [coloc](https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html)

### Planning to implement  
- [CAVIAR](http://genetics.cs.ucla.edu/caviar/)  
- [CAVIAR-BF](https://www.ncbi.nlm.nih.gov/pubmed/25948564)  
- [eCAVIAR](http://genetics.cs.ucla.edu/caviar/)  
- [DAP](https://github.com/xqwen/dap) 
- [fGWAS](https://github.com/joepickrell/fgwas)  


<br>


## Datasets

### [Nott et al. (2019)](https://science.sciencemag.org/content/366/6469/1134.abstract)
- Cell type-specific epigenomics from brain tissue, including:    
  + Epigenomic reads (as *GenomicRanges*)
  + Enhancer/promoter coordinates (as *GenomicRanges*)
  + PLAC-seq interactome  
  
### [XGR](http://xgr.r-forge.r-project.org)    
- API access to a diverse library of cell type/line-specific epigenomic (e.g. ENCODE) and other genome-wide annotations.  
### [Roadmap](http://www.roadmapepigenomics.org)
- API access to cell type-specific epigenomic data.

### [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
- API access to various genome-wide SNP annotations (e.g. missense, nonsynonmous, intronic, enhancer)

### [HaploR](https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html)  
- API access to known per-SNP QTL and epigenomic data hits.


<br>


## Enrichment Tools

### [XGR](http://xgr.r-forge.r-project.org)   
- Binomial enrichment tests between customisable foreground and background SNPs.  

### [GoShifter](https://github.com/immunogenomics/goshifter)  


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
![Sinai](./web/images/sinai.png)
