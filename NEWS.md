# echolocatoR 2.0.2

## New features

* Add new internal function: `check_topSNPs`
* Impute MAF when.
* Fully documented all functions.

## Bug fixes

* Improve handling of `conditioned_snps` in various scenarios.
* Fix GHA: @master --> @v2  

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
