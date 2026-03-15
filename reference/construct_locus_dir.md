# Make locus-specific results folder

Make locus-specific results folder

## Usage

``` r
construct_locus_dir(results_dir = tempdir(), dataset_type, dataset_name, locus)
```

## See also

Other directory functions:
[`construct_subset_path()`](https://rajlabmssm.github.io/echolocatoR/reference/construct_subset_path.md),
[`get_locus_dir()`](https://rajlabmssm.github.io/echolocatoR/reference/get_locus_dir.md),
[`get_study_dir()`](https://rajlabmssm.github.io/echolocatoR/reference/get_study_dir.md)

## Examples

``` r
locus_dir <- echolocatoR:::construct_locus_dir(dataset_type="GWAS",
                                               dataset_name="Nalls2019",
                                               locus="BST1")
```
