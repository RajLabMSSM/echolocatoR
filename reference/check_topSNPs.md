# Check `topSNPs`

Check `topSNPs`. Creates a data.frame from the full summary stats if
`topSNPs` is set to "auto".

## Usage

``` r
check_topSNPs(topSNPs, fullSS_path, verbose = TRUE, ...)
```

## Arguments

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

- fullSS_path:

  Path to the full summary statistics file (GWAS or QTL) that you want
  to fine-map. It is usually best to provide the absolute path rather
  than the relative path.

- verbose:

  Print messages.

- ...:

  Arguments passed on to
  [`echodata::import_topSNPs`](https://rdrr.io/pkg/echodata/man/import_topSNPs.html)

  `topSS`

  :   Can be a data.frame with the top summary stats per locus.
      Alternatively, you can provide a path to the stored top summary
      stats file. Can be in any tabular format (e.g. excel, .tsv, .csv,
      etc.). This file should have one lead GWAS/QTL hits per locus. If
      there is more than one SNP per locus, the one with the smallest
      p-value (then the largest effect size) is selected as the lead
      SNP. The lead SNP will be used as the center of the locus when
      constructing the locus subset files.

  `sheet`

  :   If the *topSS* file is an excel sheet, you can specify which tab
      to use. You can provide either a number to identify the tab by
      order, or a string to identify the tab by name.

  `startRow`

  :   first row to begin looking for data. Empty rows at the top of a
      file are always skipped, regardless of the value of startRow.

  `cols`

  :   A numeric vector specifying which columns in the Excel file to
      read. If NULL, all columns are read.

  `munge`

  :   Standardise column names.

  `colmap`

  :   Column mappings object. Uses
      [construct_colmap](https://rdrr.io/pkg/echodata/man/construct_colmap.html)
      by default.

  `min_POS`

  :   Column containing minimum genomic position (used instead of an
      arbitrary window size).

  `max_POS`

  :   Column containing maximum genomic position (used instead of an
      arbitrary window size).

  `grouping_vars`

  :   The variables that you want to group by such that each
      grouping_var combination has its own index SNP. For example, if
      you want one index SNP per QTL eGene - GWAS locus pair, you could
      supply: `grouping_vars=c("Locus","Gene")`.

  `remove_variants`

  :   SNPs to remove from `topSS`,

  `show_table`

  :   Create an interative data table.

## Value

data.frame
