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
