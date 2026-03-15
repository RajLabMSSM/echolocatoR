# Construct the path of the locus subset

Construct the path of the locus subset

## Usage

``` r
construct_subset_path(
  subset_path = "auto",
  results_dir = tempdir(),
  dataset_type,
  dataset_name,
  locus = NULL,
  suffix = ".tsv.gz"
)
```

## See also

Other directory functions:
[`construct_locus_dir()`](https://rajlabmssm.github.io/echolocatoR/reference/construct_locus_dir.md),
[`get_locus_dir()`](https://rajlabmssm.github.io/echolocatoR/reference/get_locus_dir.md),
[`get_study_dir()`](https://rajlabmssm.github.io/echolocatoR/reference/get_study_dir.md)

## Examples

``` r
subset_path <- echolocatoR:::construct_subset_path(dataset_type="GWAS",
                                                   dataset_name="Nalls2019",
                                                   locus="BST1")
```
