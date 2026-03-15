# Check echoverse setup

Diagnose the installation status of all echoverse packages, system
dependencies, external tools, and Python/conda environments. Prints an
actionable report with platform-specific fix commands.

## Usage

``` r
check_echoverse_setup(
  echoverse_pkgs = c("echolocatoR", "echodata", "echotabix", "echoannot", "echoconda",
    "echoLD", "echoplot", "echofinemap", "catalogueR", "downloadR", "echogithub",
    "devoptera", "echodeps", "echoAI", "echoverseTemplate"),
  key_deps = c("data.table", "BiocManager", "reticulate", "susieR", "coloc",
    "MungeSumstats", "VariantAnnotation", "rtracklayer", "GenomicRanges", "ggbio",
    "basilisk", "Rsamtools"),
  verbose = TRUE
)
```

## Arguments

- echoverse_pkgs:

  Character vector of echoverse package names to check.

- key_deps:

  Character vector of key non-echoverse dependency package names to
  check.

- verbose:

  Print detailed results. Default `TRUE`.

## Value

A list (invisibly) with elements:

- packages:

  data.frame of echoverse package status

- system:

  list of system dependency checks

- tools:

  list of external tool checks

- python:

  list of Python/conda checks

- pass:

  logical, TRUE if all critical checks pass

## Examples

``` r
if (FALSE) { # \dontrun{
## Full diagnostic report
results <- check_echoverse_setup()

## Quick check with a subset of packages
results <- check_echoverse_setup(
    echoverse_pkgs = c("echodata"),
    key_deps = c("data.table"),
    verbose = FALSE
)
} # }
```
