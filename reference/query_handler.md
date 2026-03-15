# Query handler

Handles which query method to use

## Usage

``` r
query_handler(
  fullSS_path,
  colmap,
  locus_dir = NULL,
  locus = NULL,
  topSNPs,
  subset_path,
  bp_distance = 5e+05,
  query_by = c("tabix", "fullSS"),
  force_new_subset = FALSE,
  conda_env = "echoR_mini",
  nThread = 1,
  verbose = TRUE
)
```

## Arguments

- fullSS_path:

  Path to the full summary statistics file (GWAS or QTL) that you want
  to fine-map. It is usually best to provide the absolute path rather
  than the relative path.

- colmap:

  Column name mappings in in `fullSS_path`. Must be a named list. Can
  use
  [construct_colmap](https://rdrr.io/pkg/echodata/man/construct_colmap.html)
  to assist with this. This function can be used in two different ways:

  `munged=FALSE`

  :   When `munged=FALSE`, you will need to provide the necessary column
      names to the `colmap` argument (*default*).

  `munged=TRUE`

  :   Alternatively, instead of filling out each argument in
      [construct_colmap](https://rdrr.io/pkg/echodata/man/construct_colmap.html),
      you can simply set `munged=TRUE` if `fullSS_path` has already been
      munged with
      [format_sumstats](https://al-murphy.github.io/MungeSumstats/reference/format_sumstats.html).

- locus:

  Locus name to fine-map (e.g. `"BIN1"`). Can be named to indicate a
  specific gene within a QTL locus (e.g. `c(ENSG00000136731="BIN1")`).

- topSNPs:

  A data.frame with the genomic coordinates of the lead SNP for each
  locus. The lead SNP will be used as the center of the window when
  extracting subset from the full GWAS/QTL summary statistics file. Only
  one SNP per **Locus** should be included. At minimum, `topSNPs` should
  include the following columns:

  *Locus*

  :   A unique name for each locus. Often, loci are named after a
      relevant gene (e.g. LRRK2) or based on the name/coordinates of the
      lead SNP (e.g. locus_chr12_40734202)

  *CHR*

  :   The chromosome that the SNP is on. Can be "chr12" or "12" format.

  *POS*

  :   The genomic position of the SNP (in basepairs)

- subset_path:

  Path of the resulting locus subset file.

- bp_distance:

  Distance around the lead SNP to include.

- query_by:

  Choose which method you want to use to extract locus subsets from the
  full summary stats file. Methods include:

  "tabix"

  :   Convert the full summary stats file in an indexed tabix file.
      Makes querying lightning fast after the initial conversion is
      done. (*default*)

  "coordinates"

  :   Extract locus subsets using min/max genomic coordinates with
      *awk*.

- force_new_subset:

  By default, if a subset of the full summary stats file for a given
  locus is already present, then echolocatoR will just use the
  pre-existing file. Set `force_new_subset=T` to override this and
  extract a new subset. Subsets are saved in the following path
  structure:
  *Data/\\dataset_type\\/\\dataset_name\\/\\locus\\/Multi-finemap/
  \\locus\\\_\\dataset_name\\\_Multi-finemap.tsv.gz*

- conda_env:

  Conda environment to use.

- nThread:

  Number of threads to parallelise saving across.

- verbose:

  Print messages.

## See also

Other query functions:
[`extract_snp_subset()`](https://rajlabmssm.github.io/echolocatoR/reference/extract_SNP_subset.md)
