% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CORCES_2020.scATACseq_celltype_peaks}
\alias{CORCES_2020.scATACseq_celltype_peaks}
\title{scATACseq cell type-specific peaks from human brain tissue}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 221062 rows and 13 columns.
}
\source{
\url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
}
\usage{
CORCES_2020.scATACseq_celltype_peaks
}
\description{
Each row represents an individual peak identified from the feature binarization analysis (see methods).
}
\details{
Data originally from \href{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}{Corces et al. (bioRxiv)}, as of May 2020.
Specifically: \emph{STable6_Features_scATAC-seq_celltype_Peaks}
}
\examples{
\dontrun{
dat <- readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Coceres_2020/STable6_Features_scATAC-seq_celltype_Peaks.xlsx", skip = 15)
CORCES_2020.scATACseq_celltype_peaks <- data.table::data.table(dat)
usethis::use_data(CORCES_2020.scATACseq_celltype_peaks)
}
}
\seealso{
Other CORCES_2020: 
\code{\link{CORCES_2020.HiChIP_FitHiChIP_loop_calls}},
\code{\link{CORCES_2020.bulkATACseq_peaks}},
\code{\link{CORCES_2020.cicero_coaccessibility}},
\code{\link{CORCES_2020.get_ATAC_peak_overlap}()},
\code{\link{CORCES_2020.get_HiChIP_FitHiChIP_overlap}()},
\code{\link{CORCES_2020.prepare_bulkATAC_peak_overlap}()},
\code{\link{CORCES_2020.prepare_scATAC_peak_overlap}()},
\code{\link{CORCES_2020.scATACseq_peaks}}
}
\concept{CORCES_2020}
\keyword{datasets}
