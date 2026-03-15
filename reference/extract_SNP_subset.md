# Extract a subset of the summary stats

Use *tabix* to extract a locus subset from the full summary statistics
file.

## Usage

``` r
extract_snp_subset(
  subset_path,
  locus = NULL,
  colmap = echodata::construct_colmap(),
  fullSS_path,
  topSNPs,
  LD_reference,
  force_new_subset = FALSE,
  force_new_maf = FALSE,
  bp_distance = 5e+05,
  superpopulation = "EUR",
  compute_n = "ldsc",
  query_by = "tabix",
  download_method = "axel",
  nThread = 1,
  conda_env = "echoR_mini",
  verbose = TRUE
)
```

## Arguments

- subset_path:

  Path where the `query` should be saved after standardization.

- locus:

  Locus name to fine-map (e.g. `"BIN1"`). Can be named to indicate a
  specific gene within a QTL locus (e.g. `c(ENSG00000136731="BIN1")`).

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

- fullSS_path:

  Path to the full summary statistics file (GWAS or QTL) that you want
  to fine-map. It is usually best to provide the absolute path rather
  than the relative path.

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

- LD_reference:

  LD reference to use:

  1KGphase1

  :   1000 Genomes Project Phase 1 (genome build: hg19).

  1KGphase3

  :   1000 Genomes Project Phase 3 (genome build: hg19).

  UKB

  :   Pre-computed LD from a British European-decent subset of UK
      Biobank. *Genome build* : hg19

  \<vcf_path\>

  :   User-supplied path to a custom VCF file to compute LD matrix
      from.\
      *Accepted formats*: *.vcf* / *.vcf.gz* / *.vcf.bgz*\
      *Genome build* : defined by user with `target_genome`.

  \<matrix_path\>

  :   User-supplied path to a pre-computed LD matrix. *Accepted
      formats*: *.rds* / *.rda* / *.csv* / *.tsv* / *.txt*\
      *Genome build* : defined by user with `target_genome`.

- force_new_subset:

  By default, if a subset of the full summary stats file for a given
  locus is already present, then echolocatoR will just use the
  pre-existing file. Set `force_new_subset=T` to override this and
  extract a new subset. Subsets are saved in the following path
  structure:
  *Data/\\dataset_type\\/\\dataset_name\\/\\locus\\/Multi-finemap/
  \\locus\\\_\\dataset_name\\\_Multi-finemap.tsv.gz*

- force_new_maf:

  Download UKB_MAF file again.

- bp_distance:

  Distance around the lead SNP to include.

- superpopulation:

  Superpopulation to subset LD panel by (used only if `LD_reference` is
  "1KGphase1" or "1KGphase3"). See
  [popDat_1KGphase1](https://rdrr.io/pkg/echoLD/man/popDat_1KGphase1.html)
  and
  [popDat_1KGphase3](https://rdrr.io/pkg/echoLD/man/popDat_1KGphase3.html)
  for full tables of their respective samples.

- compute_n:

  How to compute per-SNP sample size (new column "N").\
  If the column "N" is already present in `dat`, this column will be
  used to extract per-SNP sample sizes and the argument `compute_n` will
  be ignored.\
  If the column "N" is *not* present in `dat`, one of the following
  options can be supplied to `compute_n`:

  `0`

  :   N will not be computed.

  `>0`

  :   If any number \>0 is provided, that value will be set as N for
      every row. \*\*Note\*\*: Computing N this way is incorrect and
      should be avoided if at all possible.

  `"sum"`

  :   N will be computed as: cases (N_CAS) + controls (N_CON), so long
      as both columns are present.

  `"ldsc"`

  :   N will be computed as effective sample size: Neff
      =(N_CAS+N_CON)\*(N_CAS/(N_CAS+N_CON)) /
      mean((N_CAS/(N_CAS+N_CON))(N_CAS+N_CON)==max(N_CAS+N_CON)).

  `"giant"`

  :   N will be computed as effective sample size: Neff = 2 / (1/N_CAS +
      1/N_CON).

  `"metal"`

  :   N will be computed as effective sample size: Neff = 4 / (1/N_CAS +
      1/N_CON).

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

- download_method:

  Download method to use:

  `"axel"`

  :   Multi-threaded

  `"wget"`

  :   Single-threaded

  `"download.file"`

  :   Single-threaded

  `"internal"`

  :   Single-threaded (passed to
      [download.file](https://rdrr.io/r/utils/download.file.html))

  `"wininet"`

  :   Single-threaded (passed to
      [download.file](https://rdrr.io/r/utils/download.file.html))

  `"libcurl"`

  :   Single-threaded (passed to
      [download.file](https://rdrr.io/r/utils/download.file.html))

  `"curl"`

  :   Single-threaded (passed to
      [download.file](https://rdrr.io/r/utils/download.file.html))

- nThread:

  Number of threads to parallelise saving across.

- conda_env:

  Conda environment to use.

- verbose:

  Print messages.

## See also

Other query functions:
[`query_handler()`](https://rajlabmssm.github.io/echolocatoR/reference/query_handler.md)
