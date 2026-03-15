# Identify SNPs to condition on

When running conditional analyses (e.g. *GCTA-COJO*), this functions
automatically identifies SNP to condition on.

## Usage

``` r
snps_to_condition(conditioned_snps, topSNPs, loci)
```

## Arguments

- conditioned_snps:

  Which SNPs to conditions on when fine-mapping with (e.g. *COJO*).

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

- loci:

  Character list of loci in **Locus** col of `topSNPs`.

## Value

Vector of SNPs
