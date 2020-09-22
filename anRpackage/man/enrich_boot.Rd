\name{enrich_boot}
\alias{enrich_boot}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("enrich_boot")}
\format{
  A data frame with 16000 observations on the following 19 variables.
  \describe{
    \item{\code{name}}{a character vector}
    \item{\code{nAnno}}{a numeric vector}
    \item{\code{nOverlap}}{a numeric vector}
    \item{\code{fc}}{a numeric vector}
    \item{\code{zscore}}{a numeric vector}
    \item{\code{pvalue}}{a numeric vector}
    \item{\code{adjp}}{a numeric vector}
    \item{\code{or}}{a numeric vector}
    \item{\code{CIl}}{a numeric vector}
    \item{\code{CIu}}{a numeric vector}
    \item{\code{Study}}{a character vector}
    \item{\code{Cell_type}}{a character vector}
    \item{\code{Assay}}{a character vector}
    \item{\code{fg_filter}}{a character vector}
    \item{\code{bg_filter}}{a logical vector}
    \item{\code{fg_sample_size}}{a numeric vector}
    \item{\code{bg_sample_size}}{a numeric vector}
    \item{\code{SNP_Group}}{a character vector}
    \item{\code{FDR}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(enrich_boot)
## maybe str(enrich_boot) ; plot(enrich_boot) ...
}
\keyword{datasets}
