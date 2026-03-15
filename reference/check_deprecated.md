# Check deprecated arguments

Semi-automatically check all deprecated args in a given function.

## Usage

``` r
check_deprecated(
  fun = "finemap_loci",
  pkg = "echolocatoR",
  when = "2.0.0",
  args = match.call(),
  lifecycle_fun = lifecycle::deprecate_warn,
  reassign = FALSE,
  map = list(A1_col = "colmap", A2_col = "colmap", chrom_col = "colmap", position_col =
    "colmap", effect_col = "colmap", freq_col = "colmap", gene_col = "colmap", locus_col
    = "colmap", MAF_col = "colmap", N_cases = "colmap", N_controls = "colmap",
    N_cases_col = "colmap", N_controls_col = "colmap", sample_size = "colmap", MAF_col =
    "colmap", pval_col = "colmap", stderr_col = "colmap", tstat_col = "colmap", snp_col =
    "colmap", file_sep = NULL, probe_path = NULL, chrom_type = NULL, PAINTOR_QTL_datasets
    = NULL, 
     QTL_prefixes = "qtl_suffixes", proportion_cases = NULL, server = NULL,
    vcf_folder = NULL, top_SNPs = "topSNPs", PP_threshold = "credset_thresh",
    consensus_threshold = "consensus_thresh", plot.types = "plot_types", plot.Roadmap =
    "roadmap", plot.Roadmap_query = "roadmap_query", plot.XGR_libnames = "xgr_libnames",
    plot.zoom = "zoom", plot.zoom = "zoom", plot.Nott_epigenome = "nott_epigenome",
    plot.Nott_show_placseq = "nott_show_placseq")
)
```

## Arguments

- fun:

  Function to check.

- pkg:

  Package that the function is from.

- when:

  A string giving the version when the behaviour was deprecated.

- args:

  Argument calls to assess.

- lifecycle_fun:

  Which lifecycle function to use by default.

- reassign:

  Attempt to reassign deprecated variables to the corresponding new
  variable (if applicable).

- map:

  Mapping between old:new argument names. Use `NULL` if the argument is
  no longer used at all.

## Examples

``` r
if (FALSE) { # \dontrun{
topSNPs <- echodata::topSNPs_Nalls2019
fullSS_path <- echodata::example_fullSS()
testthat::expect_error(
  echolocatoR::finemap_loci(
    fullSS_path = fullSS_path,
    topSNPs = topSNPs,
    loci = c("BST1","MEX3C"),
    chrom_col = "CHR",
    position_col = "BP")
)
} # }
```
