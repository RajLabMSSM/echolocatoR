% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{merge_finemapping_results_each}
\alias{merge_finemapping_results_each}
\title{Create full cross-locus merged files for each dataset,
then return a subset of those files as one super-merged table.}
\usage{
merge_finemapping_results_each(
  study_dirs,
  LD_reference = "1KGphase3",
  minimum_support = 1,
  include_leadSNPs = T,
  return_filter = "!is.na(SNP)",
  merged_path = "merged_DT.csv.gz",
  force_new_merge = F,
  nThread = 4,
  verbose = T
)
}
\description{
Create full cross-locus merged files for each dataset,
then return a subset of those files as one super-merged table.
}
\seealso{
Other annotate: 
\code{\link{ANNOTATE.annotate_missense}()},
\code{\link{ANNOTATE.plot_missense}()},
\code{\link{SNPs_by_mutation_type}()},
\code{\link{biomart_geneInfo}()},
\code{\link{biomart_snp_info}()},
\code{\link{biomart_snps_to_geneInfo}()},
\code{\link{epigenetics_enrichment}()},
\code{\link{epigenetics_summary}()},
\code{\link{haploR.HaploReg}()},
\code{\link{haploR.regulomeDB}()}
}
\concept{annotate}
