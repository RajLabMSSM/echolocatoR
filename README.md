<img src= 'https://github.com/RajLabMSSM/echolocatoR/raw/master/inst/hex/hex.png' height= '600' ><br><br><br><br>
[![](https://img.shields.io/badge/devel%20version-2.0.3-black.svg)](https://github.com/RajLabMSSM/echolocatoR)
[![R build
status](https://github.com/RajLabMSSM/echolocatoR/workflows/rworkflows/badge.svg)](https://github.com/RajLabMSSM/echolocatoR/actions)
[![](https://img.shields.io/github/last-commit/RajLabMSSM/echolocatoR.svg)](https://github.com/RajLabMSSM/echolocatoR/commits/master)
[![](https://img.shields.io/github/languages/code-size/RajLabMSSM/echolocatoR.svg)](https://github.com/RajLabMSSM/echolocatoR)
[![](https://app.codecov.io/gh/RajLabMSSM/echolocatoR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/RajLabMSSM/echolocatoR)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
[![](https://img.shields.io/badge/doi-10.1093/bioinformatics/btab658-blue.svg)](https://doi.org/10.1093/bioinformatics/btab658)
¶ <h4> ¶ Authors: <i>Brian Schilder, Jack Humphrey, Towfique Raj</i> ¶
</h4>
<h5> ¶ README updated: <i>Dec-22-2022</i> ¶ </h5>

## `echolocatoR`: Automated statistical and functional fine-mapping

with extensive access to genome-wide datasets.

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

If you use `echolocatoR`, or any of the **echoverse** modules, please
cite:

> Brian M Schilder, Jack Humphrey, Towfique Raj (2021) echolocatoR: an
> automated end-to-end statistical and functional genomic fine-mapping
> pipeline, *Bioinformatics*; btab658,
> <https://doi.org/10.1093/bioinformatics/btab658>

## Installation

``` r
if(!require("remotes")) install.packages("remotes")

remotes::install_github("RajLabMSSM/echolocatoR")
library(echolocatoR)
```

#### Installation troubleshooting

<details>

- Because `echolocatoR` now relies on many subpackages that rely on one
  another, sometimes errors can occur when R tries to update one R
  package before updating its *echoverse* dependencies (and thus is
  unable to find new functions). As *echoverse* stabilizes over time,
  this should happen less frequently. However, in the meantime the
  solution is to simply rerun
  `remotes::install_github("RajLabMSSM/echolocatoR")` until all
  subpackages are fully updates.
- `susieR`: Sometimes an older version of `susieR` is installed from
  CRAN (e.g. 0.11.92), but `echofinemap` requires version \>= 0.12.0. To
  get around this, you can install `susieR` directly from GitHub:
  `devtools::install_github("stephenslab/susieR")`
- System dependencies can sometimes cause issues when using different
  packages. I’ve tried to account for as many of these as possible
  automatically within the code, but using the **Docker/Singularity**
  provided below can further mitigate these issues.
- The R package `XML` (which some *echoverse* subpackages depend on) has
  some additional system dependencies that must be installed beforehand.
  If `XML` does not install automatically, try installing `lbxml` on
  your system using `brew install libxml2` (MacOS),
  `sudo apt-get install libxml2` (Linux) or `conda install r-xml` if you
  are running `echolocatoR` from within a conda environment.

</details>

### \[Optional\] [Docker/Singularity](https://rajlabmssm.github.io/echolocatoR/articles/docker)

`echolocatoR` now has its own dedicated Docker/Singularity container!
This greatly reduces issues related to system dependency conflicts and
provides a containerized interface for Rstudio through your web browser.
See [here for installation
instructions](https://rajlabmssm.github.io/echolocatoR/articles/docker).

## Documentation

### [Website](https://rajlabmssm.github.io/echolocatoR)

### [Get started](https://rajlabmssm.github.io/echolocatoR/articles/echolocatoR)

### [Bugs/requests](https://github.com/RajLabMSSM/echolocatoR/issues)

Please report any bugs/requests on [GitHub
Issues](https://github.com/RajLabMSSM/echolocatoR/issues).

[Contributions](https://github.com/RajLabMSSM/echolocatoR/pulls) are
welcome!

### All *echoverse* vignettes

<details>

``` r
echoverse <- c('echolocatoR','echodata','echotabix',
               'echoannot','echoconda','echoLD',
               'echoplot','catalogueR','downloadR',
               'echofinemap','echodeps', # under construction
               'echogithub')
toc <- echogithub::github_pages_vignettes(owner = "RajLabMSSM",
                                          repo = echoverse,
                                          as_toc = TRUE,
                                          verbose = FALSE)
```

<ul class="toc-list">
<li>
<h2>
🦇 <a href='https://rajlabmssm.github.io/echolocatoR/'>echolocatoR</a>
</h2>
<ul>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echolocatoR//articles/QTLs.html' target='blank'>QTLs</a>
</h3>
</li>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echolocatoR//articles/docker.html' target='blank'>docker</a>
</h3>
</li>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echolocatoR//articles/echolocatoR.html' target='blank'>echolocatoR</a>
</h3>
</li>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echolocatoR//articles/finemapping_portal.html' target='blank'>finemapping
portal</a>
</h3>
</li>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echolocatoR//articles/plot_locus.html' target='blank'>plot
locus</a>
</h3>
</li>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echolocatoR//articles/summarise.html' target='blank'>summarise</a>
</h3>
</li>
</ul>
</li>
<li>
<h2>
🦇 <a href='https://rajlabmssm.github.io/echodata/'>echodata</a>
</h2>
<ul>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echodata//articles/echodata.html' target='blank'>echodata</a>
</h3>
</li>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echodata//articles/echolocatoR_Finemapping_Portal.html' target='blank'>echolocatoR
Finemapping Portal</a>
</h3>
</li>
</ul>
</li>
<li>
<h2>
🦇 <a href='https://rajlabmssm.github.io/echotabix/'>echotabix</a>
</h2>
<ul>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echotabix//articles/echotabix.html' target='blank'>echotabix</a>
</h3>
</li>
</ul>
</li>
<li>
<h2>
🦇 <a href='https://rajlabmssm.github.io/echoannot/'>echoannot</a>
</h2>
<ul>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echoannot//articles/cell_type_specific_epigenomics.html' target='blank'>cell
type specific epigenomics</a>
</h3>
</li>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echoannot//articles/echoannot.html' target='blank'>echoannot</a>
</h3>
</li>
</ul>
</li>
<li>
<h2>
🦇 <a href='https://rajlabmssm.github.io/echoconda/'>echoconda</a>
</h2>
<ul>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echoconda//articles/echoconda.html' target='blank'>echoconda</a>
</h3>
</li>
</ul>
</li>
<li>
<h2>
🦇 <a href='https://rajlabmssm.github.io/echoLD/'>echoLD</a>
</h2>
<ul>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echoLD//articles/echoLD.html' target='blank'>echoLD</a>
</h3>
</li>
</ul>
</li>
<li>
<h2>
🦇 <a href='https://rajlabmssm.github.io/echoplot/'>echoplot</a>
</h2>
<ul>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echoplot//articles/echoplot.html' target='blank'>echoplot</a>
</h3>
</li>
</ul>
</li>
<li>
<h2>
🦇 <a href='https://rajlabmssm.github.io/echofinemap/'>echofinemap</a>
</h2>
<ul>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echofinemap//articles/echoverseTemplate.html' target='blank'>echoverseTemplate</a>
</h3>
</li>
</ul>
</li>
<li>
<h2>
🦇 <a href='https://rajlabmssm.github.io/echodeps/'>echodeps</a>
</h2>
<ul>
<li>
<h3>
<a href='https://rajlabmssm.github.io/echodeps//articles/echoverseTemplate.html' target='blank'>echoverseTemplate</a>
</h3>
</li>
</ul>
</li>
</ul>
</details>

## Introduction

Fine-mapping methods are a powerful means of identifying causal variants
underlying a given phenotype, but are underutilized due to the technical
challenges of implementation. `echolocatoR` is an R package that
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
results sharing. Therefore `echolocatoR` drastically reduces the
barriers to identifying causal variants by making the entire
fine-mapping pipeline rapid, robust and scalable.

<img src="./images/echolocatoR_Fig1.png" style="border-radius: 10px;">

## Literature

### For applications of `echolocatoR` in the literature, please see:

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

## `echolocatoR` v1.0 vs. v2.0

There have been a series of major updates between `echolocatoR` v1.0 and
v2.0. Here are some of the most notable ones (see **Details**):

<details>

- ***echoverse* subpackages**: `echolocatoR` has been broken into
  separate subpackages, making it much easier to edit/debug each step of
  the full `finemap_loci` pipeline, and improving robustness throughout.
  It also provides greater flexibility for users to construct their own
  custom pipelines from these modules.
- **`GITHUB_TOKEN`**: GitHub now requires users to create Personal
  Authentication Tokens (PAT) to avoid download limits. This is
  essential for installing `echolocatoR` as many resources from GitHub
  need to be downloaded. See
  [here](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)
  for further instructions. = `echodata::construct_colmap()`:
  Previously, users were required to input key column name mappings as
  separate arguments to `echolocatoR::finemap_loci`. This functionality
  has been deprecated and replaced with a single argument, `colmap=`.
  This allows users to save the `construct_colmap()` output as a single
  variable and reuse it later without having to write out each mapping
  argument again (and helps reduce an already crowded list of
  arguments).
- **`MungeSumstats`**: `finemap_loci` now accepts the output of
  [`MungeSumstats::format_sumstats`/`import_sumstats`](https://github.com/neurogenomics/MungeSumstats)
  as-is (without requiring `colmap=`, so long as `munged=TRUE`).
  Standardizing your GWAS/QTL summary stats this way greatly reduces (or
  eliminates) the time taken to do manual formatting.
- **`echolocatoR::finemap_loci` arguments**: Several arguments have been
  deprecated or had their names changed to be harmonized across all the
  subpackages and use a unified naming convention. See
  `?echolocatoR::finemap_loci` for details.
- **`echoconda`**: The *echoverse* subpackage `echoconda` now handles
  all conda environment creation/use internally and automatically,
  without the need for users to create the conda environment themselves
  as a separate step. Also, the default conda env `echoR` has been
  replaced by `echoR_mini`, which reduces the number of dependencies to
  just the bare minimum (thus greatly speeding up build time and
  reducing potential version conflicts).
- **`FINEMAP`**: More outputs from the tool `FINEMAP` are now recorded
  in the `echolocatoR` results (see `?echofinemap::FINEMAP` or [this
  Issue](https://github.com/RajLabMSSM/echofinemap/issues/7) for
  details). Also, a common dependency conflict between `FINEMAP`\>=1.4
  and MacOS has been resolved (see [this
  Issue](https://github.com/RajLabMSSM/echofinemap/issues/9) for
  details.
- **`echodata`**: All example data and data transformation functions
  have been moved to the *echoverse* subpackage
  [`echodata`](https://github.com/RajLabMSSM/echodata).
- **`LD_reference=`**: In addition to the *UKB*, *1KGphase1/3* LD
  reference panels, `finemap_loci()` can now take custom LD panels by
  supplying `finemap_loci(LD_reference=)` with a list of paths to VCF
  files (.vcf / vcf.gz / vcf.bgz) or pre-computed LD matrices with RSIDs
  as the row/col names (.rda / .rds / .csv / .tsv. / .txt / .csv.gz /
  tsv.gz / txt.gz).
- **Expanded fine-mapping methods**: “ABF”, “COJO_conditional”,
  “COJO_joint” “COJO_stepwise”,“FINEMAP”,“PAINTOR” (including multi-GWAS
  and multi-ancestry fine-mapping),“POLYFUN_FINEMAP”
  ,“POLYFUN_SUSIE”,“SUSIE”
- **`FINEMAP` fixed**: There were a number of issues with `FINEMAP` due
  to differing output formats across different versions, system
  dependency conflicts, and the fact that it can produce multiple
  Credible Sets. All of these have been fixed and the latest version of
  `FINEMAP` can be run on all OS platforms.  
- **Debug mode**: Within `finemap_loci()` I use a `tryCatch()` when
  iterating across loci so that if one locus fails, the rest can
  continue. However this prevents using traceback feature in R, making
  debugging hard. Thus I now enabled debugging mode via a new argument:
  `use_tryCatch=FALSE`.

</details>

## Output descriptions

By default, `echolocatoR::finemap_loci()` returns a nested list
containing grouped by locus names (e.g. `$BST1`, `$MEX3C`). The results
of each locus contain the following elements:

<details>

- `finemap_dat`: Fine-mapping results from all selected methods merged
  with the original summary statistics (i.e. **Multi-finemap results**).
- `locus_plot`: A nested list containing one or more zoomed views of
  locus plots.  
- `LD_matrix`: The post-processed LD matrix used for fine-mapping.
- `LD_plot`: An LD plot (if used).
- `locus_dir`: Locus directory results are saved in.
- `arguments`: A record of the arguments supplied to `finemap_loci`.

In addition, the following object summarizes the results from the
locus-specific elements:  
- `merged_dat`: A merged `data.table` with all fine-mapping results from
all loci.

### Multi-finemap results files

The main output of `echolocatoR` are the multi-finemap files (for
example, `echodata::BST1`). They are stored in the locus-specific
*Multi-finemap* subfolders.

#### Column descriptions

- **Standardized GWAS/QTL summary statistics**: e.g.
  `SNP`,`CHR`,`POS`,`Effect`,`StdErr`. See `?finemap_loci()` for
  descriptions of each.  
- **leadSNP**: The designated proxy SNP per locus, which is the SNP with
  the smallest p-value by default.
- **\<tool\>.CS**: The 95% probability Credible Set (CS) to which a SNP
  belongs within a given fine-mapping tool’s results. If a SNP is not in
  any of the tool’s CS, it is assigned `NA` (or `0` for the purposes of
  plotting).  
- **\<tool\>.PP**: The posterior probability that a SNP is causal for a
  given GWAS/QTL trait.  
- **Support**: The total number of fine-mapping tools that include the
  SNP in its CS.
- **Consensus_SNP**: By default, defined as a SNP that is included in
  the CS of more than `N` fine-mapping tool(s), i.e. `Support>1`
  (default: `N=1`).  
- **mean.PP**: The mean SNP-wise PP across all fine-mapping tools used.
- **mean.CS**: If mean PP is greater than the 95% probability threshold
  (`mean.PP>0.95`) then `mean.CS` is 1, else 0. This tends to be a very
  stringent threshold as it requires a high degree of agreement between
  fine-mapping tools.

### Notes

- Separate multi-finemap files are generated for each LD reference panel
  used, which is included in the file name (e.g.
  *UKB_LD.Multi-finemap.tsv.gz*).
- Each fine-mapping tool defines its CS and PP slightly differently, so
  please refer to the associated original publications for the exact
  details of how these are calculated (links provided above).

</details>

## Fine-mapping tools

Fine-mapping functions are now implemented via
[`echofinemap`](https://github.com/RajLabMSSM/echofinemap):

<details>

- `echolocatoR` will automatically check whether you have the necessary
  columns to run each tool you selected in
  `echolocatoR::finemap_loci(finemap_methods=...)`. It will remove any
  tools that for which there are missing necessary columns, and produces
  a message letting you know which columns are missing.
- Note that some columns (e.g. `MAF`,`N`,`t-stat`) will be automatically
  inferred if missing.  
- For easy reference, we list the necessary columns here as well.  
  See `?echodata::construct_colmap()` for descriptions of these
  columns.  
  All methods require the columns: `SNP`,`CHR`,`POS`,`Effect`,`StdErr`

</details>

<br>

``` r
fm_methods <- echofinemap::required_cols(add_versions = FALSE, 
                                         embed_links = TRUE,
                                         verbose = FALSE)
```

    ## Registered S3 method overwritten by 'GGally':
    ##   method from   
    ##   +.gg   ggplot2

``` r
knitr::kable(x = fm_methods)
```

| method           | required   | suggested  | source                                             | citation                                                  |
|:-----------------|:-----------|:-----------|:---------------------------------------------------|:----------------------------------------------------------|
| ABF              | SNP, CHR…. |            | [source](https://github.com/chr1swallace/coloc)    | [cite](https://doi.org/10.1086%2F519024)                  |
| COJO_conditional | SNP, CHR…. | Freq, P, N | [source](https://github.com/jianyangqt/gcta)       | [cite](https://doi.org/10.1038/ng.2213)                   |
| COJO_joint       | SNP, CHR…. | Freq, P, N | [source](https://github.com/jianyangqt/gcta)       | [cite](https://doi.org/10.1038/ng.2213)                   |
| COJO_stepwise    | SNP, CHR…. | Freq, P, N | [source](https://github.com/jianyangqt/gcta)       | [cite](https://doi.org/10.1038/ng.2213)                   |
| FINEMAP          | SNP, CHR…. | A1, A2, …. | [source](http://www.christianbenner.com/)          | [cite](https://doi.org/10.1093%2Fbioinformatics%2Fbtw018) |
| PAINTOR          | SNP, CHR…. | MAF        | [source](https://github.com/gkichaev/PAINTOR_V3.0) | [cite](https://doi.org/10.1093/bioinformatics/btw615)     |
| POLYFUN_FINEMAP  | SNP, CHR…. | MAF, N     | [source](https://github.com/omerwe/polyfun)        | [cite](https://doi.org/10.1038/s41588-022-01036-9)        |
| POLYFUN_SUSIE    | SNP, CHR…. | MAF, N     | [source](https://github.com/omerwe/polyfun)        | [cite](https://doi.org/10.1038/s41588-022-01036-9)        |
| SUSIE            | SNP, CHR…. | N          | [source](https://github.com/stephenslab/susieR)    | [cite](https://doi.org/10.1371/journal.pgen.1010299)      |

## Datasets

Datasets are now stored/retrieved via the following **echoverse**
subpackages:  
- [`echodata`](https://github.com/RajLabMSSM/echodata): Pre-computed
fine-mapping results. Also handles the semi-automated standardization of
summary statistics.  
- [`echoannot`](https://github.com/RajLabMSSM/echoannot): Annotates
GWAS/QTL summary statistics using epigenomics, pre-compiled annotation
matrices, and machine learning model predictions of variant-specific
functional impacts.  
- [`catalogueR`](https://github.com/RajLabMSSM/catalogueR): Large
compendium of fully standardized e/s/t-QTL summary statistics.

For more detailed information about each dataset, use `?`:

``` r
### Examples ###

library(echoannot)   
?NOTT_2019.interactome # epigenomic annotations
library(echodata) 
?BST1 # fine-mapping results 
```

<details>

### [**`MungeSumstats`**](https://github.com/neurogenomics/MungeSumstats):

- You can search, import, and standardize any GWAS in the [*Open
  GWAS*](https://gwas.mrcieu.ac.uk/) database via
  [`MungeSumstats`](https://github.com/neurogenomics/MungeSumstats),
  specifically the functions `find_sumstats` and `import_sumstats`.

### [`catalogueR`](https://github.com/RajLabMSSM/catalogueR): QTLs

#### [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/): `catalogueR::eQTL_Catalogue.query()`

- API access to full summary statistics from many standardized e/s/t-QTL
  datasets.  
- Data access and colocalization tests facilitated through the
  [`catalogueR`](https://github.com/RajLabMSSM/catalogueR) R package.

### [`echodata`](https://github.com/RajLabMSSM/catalogueR): fine-mapping results

#### [***echolocatoR Fine-mapping Portal***](https://rajlab.shinyapps.io/Fine_Mapping_Shiny): pre-computed fine-mapping results

- You can visit the *echolocatoR Fine-mapping Portal* to interactively
  visualize and download pre-computed fine-mapping results across a
  variety of phenotypes.
- This data can be searched and imported programmatically using
  `echodata::portal_query()`.

### [`echoannot`](https://github.com/RajLabMSSM/echoannot): Epigenomic & genome-wide annotations

#### [Nott et al. (2019)](https://science.sciencemag.org/content/366/6469/1134.abstract): `echoannot::NOTT2019_*()`

- Data from this publication contains results from cell type-specific
  (neurons, oligodendrocytes, astrocytes, microglia, & peripheral
  myeloid cells) epigenomic assays (H3K27ac, ATAC, H3K4me3) from *ex
  vivo* pediatric human brain tissue.

#### [Corces et al.2020](https://doi.org/10.1038/s41588-020-00721-x): `echoannot::CORCES2020_*()`

- Data from this publication contains results from single-cell and bulk
  chromatin accessibility assays (\[sc\]ATAC-seq) and chromatin
  interactions ( [`FitHiChIP`](https://ay-lab.github.io/FitHiChIP/))
  from *postmortem* adult human brain tissue.

#### [XGR](http://xgr.r-forge.r-project.org): `echoannot::XGR_download_and_standardize()`

- API access to a diverse library of cell type/line-specific epigenomic
  (e.g. **ENCODE**) and other genome-wide annotations.

#### [Roadmap](http://www.roadmapepigenomics.org): `echoannot::ROADMAP_query()`

- API access to cell type-specific epigenomic data.

#### [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html): `echoannot::annotate_snps()`

- API access to various genome-wide SNP annotations (e.g. missense,
  nonsynonmous, intronic, enhancer).

#### [HaploR](https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html): `echoannot::annotate_snps()`

- API access to known per-SNP QTL and epigenomic data hits.

</details>

## Enrichment tools

Annotation enrichment functions are now implemented via
[`echoannot`](https://github.com/RajLabMSSM/echoannot):

<details>

### Implemented

#### [XGR](http://xgr.r-forge.r-project.org): `echoannot::XGR_enrichment()`

- Binomial enrichment tests between customisable foreground and
  background SNPs.

#### [motifbreakR](https://github.com/Simon-Coetzee/motifBreakR): `echoannot::MOTIFBREAKR()`

- Identification of transcript factor binding motifs (TFBM) and
  prediction of SNP disruption to said motifs.
- Includes a comprehensive list of TFBM databases via
  [MotifDB](https://bioconductor.org/packages/release/bioc/html/MotifDb.html)
  (9,900+ annotated position frequency matrices from 14 public sources,
  for multiple organisms).

#### [regioneR](http://bioconductor.org/packages/release/bioc/html/regioneR.html): `echoannot::test_enrichment()`

- Iterative pairwise permutation testing of overlap between all
  combinations of two
  [`GRangesList`](https://biodatascience.github.io/compbio/bioc/GRL.html)
  objects.

### Under construction

#### [GARFIELD](https://www.bioconductor.org/packages/release/bioc/html/garfield.html)

- Genomic enrichment with LD-informed heuristics.

#### [GoShifter](https://github.com/immunogenomics/goshifter)

- LD-informed iterative enrichment analysis.

#### [S-LDSC](https://www.nature.com/articles/ng.3954)

- Genome-wide stratified LD score regression.
- Inlccles 187-annotation baseline model from [Gazal et al.
  2018](https://www.nature.com/articles/s41588-018-0231-8).  
- You can alternatively supply a custom annotations matrix.

</details>

## LD reference panels

LD reference panels are now queried/processed by
[`echoLD`](https://github.com/RajLabMSSM/echoLD), specifically the
function `get_LD()`:

<details>

### [UK Biobank](https://www.ukbiobank.ac.uk)

### [1000 Genomes Phase 1](https://www.internationalgenome.org)

### [1000 Genomes Phase 3](https://www.internationalgenome.org)

### Custom LD panel:

- From user-supplied VCFs

### Custom LD panel

- From user-supplied precomputed LD matrices

</details>

## Plotting

Plotting functions are now implemented via:  
- [`echoplot`](https://github.com/RajLabMSSM/echoplot): Multi-track
locus plots with GWAS, fine-mapping results, and functional annotations
(`plot_locus()`). Can also plot multi-GWAS/QTL and multi-ancestry
results (`plot_locus_multi()`).  
- [`echoannot`](https://github.com/RajLabMSSM/echoannot): Study-level
summary plots showing aggregted info across many loci at once
(`super_summary_plot()`).  
- [`echoLD`](https://github.com/RajLabMSSM/echoLD): Plot an LD matrix
using one of several differnt plotting methods (`plot_LD()`).

## Tabix queries

All queries of [`tabix`](http://www.htslib.org/doc/tabix.html)-indexed
files (for rapid data subset extraction) are implemented via
[`echotabix`](https://github.com/RajLabMSSM/echotabix).

<details>

- `echotabix::convert_and_query()` detects whether the GWAS summary
  statistics file you provided is already `tabix`-indexed, and it not,
  automatically performs all steps necessary to convert it (sorting,
  `bgzip`-compression, indexing) across a wide variety of scenarios.  
- `echotabix::query()` contains many different methods for making tabix
  queries
  (e.g. `Rtracklayer`,`echoconda`,`VariantAnnotation`,`seqminer`), each
  of which fail in certain circumstances. To avoid this, `query()`
  automatically selects the method that will work for the particular
  file being queried and your machine’s particular versions of
  R/Bioconductor/OS, taking the guesswork and troubleshooting out of
  `tabix` queries.

</details>

## Downloads

Single- and multi-threaded downloads are now implemented via
[`downloadR`](https://github.com/RajLabMSSM/downloadR).

<details>

- Multi-threaded downloading is performed using
  [`axel`](https://github.com/axel-download-accelerator/axel), and is
  particularly useful for speeding up downloads of large files.
- `axel` is installed via the official *echoverse*
  [conda](https://docs.conda.io/en/latest/) environment: “echoR_mini”.
  This environment is automatically created by the function
  `echoconda::yaml_to_env()` when needed.

</details>
<hr>

# Developer

<a href="https://bschilder.github.io/BMSchilder/" target="_blank">Brian
M. Schilder, Bioinformatician II</a>  
<a href="https://rajlab.org" target="_blank">Raj Lab</a>  
<a href="https://icahn.mssm.edu/about/departments/neuroscience" target="_blank">Department
of Neuroscience, Icahn School of Medicine at Mount Sinai</a>

<hr>

# Session info

<details>

``` r
utils::sessionInfo()
```

```
## R version 4.2.1 (2022-06-23)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur ... 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.2                  reticulate_1.26            
##   [3] R.utils_2.12.2              tidyselect_1.2.0           
##   [5] RSQLite_2.2.20              AnnotationDbi_1.60.0       
##   [7] htmlwidgets_1.6.0           grid_4.2.1                 
##   [9] BiocParallel_1.32.4         echogithub_0.99.1          
##  [11] XGR_1.1.8                   munsell_0.5.0              
##  [13] codetools_0.2-18            interp_1.1-3               
##  [15] DT_0.26                     colorspace_2.0-3           
##  [17] OrganismDbi_1.40.0          Biobase_2.58.0             
##  [19] filelock_1.0.2              highr_0.10                 
##  [21] knitr_1.41                  supraHex_1.36.0            
##  [23] rstudioapi_0.14             stats4_4.2.1               
##  [25] DescTools_0.99.47           gitcreds_0.1.2             
##  [27] MatrixGenerics_1.10.0       GenomeInfoDbData_1.2.9     
##  [29] mixsqp_0.3-48               bit64_4.0.5                
##  [31] echoconda_0.99.9            rprojroot_2.0.3            
##  [33] basilisk_1.10.2             vctrs_0.5.1                
##  [35] generics_0.1.3              xfun_0.36                  
##  [37] biovizBase_1.46.0           BiocFileCache_2.6.0        
##  [39] R6_2.5.1                    GenomeInfoDb_1.34.4        
##  [41] AnnotationFilter_1.22.0     bitops_1.0-7               
##  [43] cachem_1.0.6                reshape_0.8.9              
##  [45] DelayedArray_0.24.0         assertthat_0.2.1           
##  [47] BiocIO_1.8.0                scales_1.2.1               
##  [49] nnet_7.3-18                 rootSolve_1.8.2.3          
##  [51] gtable_0.3.1                ggbio_1.46.0               
##  [53] lmom_2.9                    ensembldb_2.22.0           
##  [55] rlang_1.0.6                 echodata_0.99.16           
##  [57] splines_4.2.1               lazyeval_0.2.2             
##  [59] rtracklayer_1.58.0          dichromat_2.0-0.1          
##  [61] hexbin_1.28.2               checkmate_2.1.0            
##  [63] reshape2_1.4.4              BiocManager_1.30.19        
##  [65] yaml_2.3.6                  backports_1.4.1            
##  [67] snpStats_1.48.0             GenomicFeatures_1.50.3     
##  [69] ggnetwork_0.5.10            Hmisc_4.7-2                
##  [71] RBGL_1.74.0                 tools_4.2.1                
##  [73] ggplot2_3.4.0               ellipsis_0.3.2             
##  [75] RColorBrewer_1.1-3          proxy_0.4-27               
##  [77] BiocGenerics_0.44.0         coloc_5.1.0.1              
##  [79] Rcpp_1.0.9                  plyr_1.8.8                 
##  [81] base64enc_0.1-3             progress_1.2.2             
##  [83] zlibbioc_1.44.0             purrr_1.0.0                
##  [85] RCurl_1.98-1.9              basilisk.utils_1.10.0      
##  [87] prettyunits_1.1.1           rpart_4.1.19               
##  [89] deldir_1.0-6                viridis_0.6.2              
##  [91] S4Vectors_0.36.1            cluster_2.1.4              
##  [93] SummarizedExperiment_1.28.0 ggrepel_0.9.2              
##  [95] fs_1.5.2                    here_1.0.1                 
##  [97] crul_1.3                    magrittr_2.0.3             
##  [99] data.table_1.14.6           echotabix_0.99.8           
## [101] dnet_1.1.7                  openxlsx_4.2.5.1           
## [103] gh_1.3.1                    mvtnorm_1.1-3              
## [105] ProtGenerics_1.30.0         matrixStats_0.63.0         
## [107] patchwork_1.1.2             hms_1.1.2                  
## [109] evaluate_0.19               rworkflows_0.99.3          
## [111] XML_3.99-0.13               jpeg_0.1-10                
## [113] readxl_1.4.1                IRanges_2.32.0             
## [115] gridExtra_2.3               testthat_3.1.6             
## [117] compiler_4.2.1              biomaRt_2.54.0             
## [119] tibble_3.1.8                crayon_1.5.2               
## [121] R.oo_1.25.0                 htmltools_0.5.4            
## [123] echoannot_0.99.10           tzdb_0.3.0                 
## [125] Formula_1.2-4               tidyr_1.2.1                
## [127] expm_0.999-6                Exact_3.2                  
## [129] DBI_1.1.3                   dbplyr_2.2.1               
## [131] MASS_7.3-58.1               rappdirs_0.3.3             
## [133] boot_1.3-28.1               dlstats_0.1.6              
## [135] Matrix_1.5-3                badger_0.2.2               
## [137] readr_2.1.3                 piggyback_0.1.4            
## [139] brio_1.1.3                  cli_3.5.0                  
## [141] R.methodsS3_1.8.2           parallel_4.2.1             
## [143] echofinemap_0.99.4          igraph_1.3.5               
## [145] GenomicRanges_1.51.5        pkgconfig_2.0.3            
## [147] rvcheck_0.2.1               GenomicAlignments_1.34.0   
## [149] dir.expiry_1.6.0            RCircos_1.2.2              
## [151] foreign_0.8-84              osfr_0.2.9                 
## [153] xml2_1.3.3                  XVector_0.38.0             
## [155] yulab.utils_0.0.6           echoLD_0.99.8              
## [157] stringr_1.5.0               VariantAnnotation_1.44.0   
## [159] digest_0.6.31               graph_1.76.0               
## [161] httpcode_0.3.0              Biostrings_2.66.0          
## [163] rmarkdown_2.19              cellranger_1.1.0           
## [165] htmlTable_2.4.1             gld_2.6.6                  
## [167] restfulr_0.0.15             curl_4.3.3                 
## [169] Rsamtools_2.14.0            rjson_0.2.21               
## [171] lifecycle_1.0.3             nlme_3.1-161               
## [173] jsonlite_1.8.4              desc_1.4.2                 
## [175] viridisLite_0.4.1           BSgenome_1.66.1            
## [177] fansi_1.0.3                 downloadR_0.99.5           
## [179] pillar_1.8.1                susieR_0.12.27             
## [181] GGally_2.1.2                lattice_0.20-45            
## [183] KEGGREST_1.38.0             fastmap_1.1.0              
## [185] httr_1.4.4                  survival_3.4-0             
## [187] glue_1.6.2                  zip_2.2.2                  
## [189] png_0.1-8                   bit_4.0.5                  
## [191] Rgraphviz_2.42.0            class_7.3-20               
## [193] stringi_1.7.8               blob_1.2.3                 
## [195] latticeExtra_0.6-30         memoise_2.0.1              
## [197] dplyr_1.0.10                irlba_2.3.5.1              
## [199] e1071_1.7-12                ape_5.6-2
```

</details>

<br>
