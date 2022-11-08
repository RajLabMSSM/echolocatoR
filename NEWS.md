# echolocatoR 2.0.3

## New features

* Implement new `rworkflows` version of GHA.
* Remove *Dockerfile*.
* Move "plot_locus.Rmd" vignette to `echoplot`.

## Bug fixes

* Push Docker.


# echolocatoR 2.0.2

## New features

* Add new internal function: `check_topSNPs`
* Impute MAF when.
* Fully documented all functions.
* Updated manual.
* `finemap_loci`:
  - Throw an error when no `loci` are available to fine-map after locus filtering.
* Passing all CRAN checks for the first time!!
* Allow more args to have >1 input 
  (and make these vars more consistently named):
  - `n_causal_i`
  - `bp_distance_i`
  - `LD_genome_build_i`
* Create dedicated `echolocatoR` *docker/singualarit* container:
  - Add *docker* vignette
  - Update README

## Bug fixes

* Improve handling of `conditioned_snps` in various scenarios.
* Fix GHA: @master --> @v2  
* Fix vignettes:
  - *plot_locus*
* Fix `Error in .new_IRanges_from_start_end(start, end): 'start' or 'end' cannot contain NAs`
  - Was due to `loci` not present in `topSNPs` being included. 
  Now handled by `echodata::gene_locus_list`.
* Add Suggests used in vignettes: `ggplot2`,`patchwork`
* `check_genome`:
  - Omit tests that require large optional databases to be installed.
  - Update dbSNP default to 155.
* README
  - Set `echofinemap::required_cols(add_versions = FALSE,...)` 
  to avoid issues with not being able to find PolyFun/PAINTOR 
  installations during rendering.
  
# echolocatoR 2.0.1

## New features

* New args in `finemap_loci`: 
  - `seed`
  - `priors_col`
* Pass `compute_n` to `echofinemap::multifinemap`.
* `echoplot` updates:
  - Passed up new arg from `echoplot::plot_locus(tx_biotypes=)`. 
  - Changed arg `qtl_prefixes` --> `qtl_suffixes`.
  - Passed up `show_plot` to `finemap_loci`.
  - Fixed 0 transcripts bug.
* Can now supply `LD_reference` with a list of vcf/csv/tsv/txt files
  (and their compressed versions), or rda/rds files. 

## Bug fixes

* Fix typo in `tryCatch` in `finemap_loci`. 
* Deprecate `PAINTOR_QTL_datasets` arg.
* `finemap_locus`: Recording args causing an error every time 
  (`object 'locus' not found`) even when wrapped in `tryCatch`. 
  Setting `arguments <- NULL` for now.
* Fix old arg in: `plot.types` --> `plot_types`

# echolocatoR 2.0.0

## New features

* Split echolocatoR into task-specific sub-packages ("modules"). See details [here](https://github.com/RajLabMSSM/echolocatoR/issues/62).  
* Collapse column-name mapping args into one argument: `construct_colmap`
* Added a `NEWS.md` file to track changes to the package.
* QTL studies now supported by `finemap_loci`.
* Remove *inst/tools* and moved to `echofinemap` as git submodules.
* Add startup message.
* Add messages with icons for each step. 
* Record arguments as a named list and include it in the results list. 
* Add (interim) Docker/Singularity instructions to README.

## Bug fixes

- FINEMAP:
  - File parsing issues, esp across FINEMAP versions. 
  - dylib errors
  - Multiple CS per SNP
