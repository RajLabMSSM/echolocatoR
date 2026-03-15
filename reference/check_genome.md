# Check genome build

If `fullSS_genome_build==NULL` and `munged=TRUE`, infers genome build
(hg19 vs. hg38) from summary statistics using
[get_genome_builds](https://al-murphy.github.io/MungeSumstats/reference/get_genome_builds.html).
This can only be done with summary statistics that have already been
munged by
[format_sumstats](https://al-murphy.github.io/MungeSumstats/reference/format_sumstats.html).
When `fullSS_genome_build` is a synonym of hg19 or hg38, this function
simply returns a standardized version of the user-provided genome build.

## Usage

``` r
check_genome(
  fullSS_genome_build = NULL,
  munged = FALSE,
  fullSS_path = NULL,
  sampled_snps = 10000,
  names_from_paths = TRUE,
  dbSNP = 155,
  nThread = 1,
  verbose = TRUE
)
```

## Arguments

- fullSS_genome_build:

  Genome build of the full summary statistics (`fullSS_path`). Can be
  "GRCH37" or "GRCH38" or one of their synonyms.. If
  `fullSS_genome_build==NULL` and `munged=TRUE`, infers genome build
  (hg19 vs. hg38) from summary statistics using
  [get_genome_builds](https://al-murphy.github.io/MungeSumstats/reference/get_genome_builds.html).

- munged:

  Whether `fullSS_path` have already been standardised/filtered full
  summary stats with
  [format_sumstats](https://al-murphy.github.io/MungeSumstats/reference/format_sumstats.html).
  If `munged=FALSE` you'll need to provide the necessary column names to
  the `colmap` argument.

- fullSS_path:

  Path to the full summary statistics file (GWAS or QTL) that you want
  to fine-map. It is usually best to provide the absolute path rather
  than the relative path.

- sampled_snps:

  Downsample the number of SNPs used when inferring genome build to save
  time.

- names_from_paths:

  Infer the name of each item in `sumstats_list` from its respective
  file path. Only works if `sumstats_list` is a list of paths.

- dbSNP:

  version of dbSNP to be used (144 or 155). Default is 155.

- nThread:

  Number of threads to parallelise saving across.

- verbose:

  Print messages.

## Value

Character string indicating genome build.

## Examples

``` r
## When the build is already known, simply standardizes the name
build <- check_genome(fullSS_genome_build = "hg19")
print(build)
#> [1] "GRCH37"
```
