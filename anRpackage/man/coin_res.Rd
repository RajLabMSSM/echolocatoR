\name{coin_res}
\alias{coin_res}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("coin_res")}
\format{
  The format is:
Formal class 'MaxTypeIndependenceTest' [package "coin"] with 5 slots
  ..@ distribution:Formal class 'AsymptNullDistribution' [package "coin"] with 11 slots
  .. .. ..@ seed          : int [1:626] 10403 38 -1552939936 702769827 -120062527 1154308528 -53770034 -1732004295 -1353488309 -7734790 ...
  .. .. ..@ size          :function (alpha, type)  
  .. .. ..@ pvalueinterval:function (q)  
  .. .. ..@ midpvalue     :function (q)  
  .. .. ..@ q             :function (p, ...)  
  .. .. ..@ d             :function (x)  
  .. .. ..@ support       :function ()  
  .. .. ..@ parameters    :List of 1
  .. .. .. ..$ corr: num [1:4, 1:4] 1 -0.333 -0.333 -0.333 -0.333 ...
  .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. ..$ : chr [1:4] "Consensus" "GWAS lead" "Random" "UCS"
  .. .. .. .. .. ..$ : chr [1:4] "Consensus" "GWAS lead" "Random" "UCS"
  .. .. ..@ pvalue        :function (q, ...)  
  .. .. ..@ p             :function (q, ...)  
  .. .. ..@ name          : chr "Multivariate Normal Distribution"
  ..@ statistic   :Formal class 'MaxTypeIndependenceTestStatistic' [package "coin"] with 14 slots
  .. .. ..@ alternative                : chr "two.sided"
  .. .. ..@ teststatistic              : num 5.63
  .. .. ..@ standardizedlinearstatistic: Named num [1:4] -5.63 4.09 3.07 -1.54
  .. .. .. ..- attr(*, "names")= chr [1:4] "Consensus" "GWAS lead" "Random" "UCS"
  .. .. ..@ linearstatistic            : num [1:4] 26.3 52.8 50 37.5
  .. .. ..@ expectation                : Named num [1:4] 41.7 41.7 41.7 41.7
  .. .. .. ..- attr(*, "names")= chr [1:4] "Consensus" "GWAS lead" "Random" "UCS"
  .. .. ..@ covariance                 :Formal class 'Variance' [package "coin"] with 1 slot
  .. .. .. .. ..@ variance: Named num [1:4] 7.4 7.4 7.4 7.4
  .. .. .. .. .. ..- attr(*, "names")= chr [1:4] "Consensus" "GWAS lead" "Random" "UCS"
  .. .. ..@ xtrans                     : num [1:400, 1:4] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:400] "1" "2" "3" "4" ...
  .. .. .. .. ..$ : chr [1:4] "Consensus" "GWAS lead" "Random" "UCS"
  .. .. .. ..- attr(*, "assign")= int [1:4] 1 1 1 1
  .. .. ..@ ytrans                     : num [1:400, 1] 0.218 0.353 0.579 1 0.28 ...
  .. .. .. ..- attr(*, "assign")= int 1
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : chr ""
  .. .. ..@ xtrafo                     :function (data, numeric_trafo = id_trafo, factor_trafo = f_trafo, ordered_trafo = of_trafo, surv_trafo = logrank_trafo, 
    var_trafo = NULL, block = NULL)  
  .. .. ..@ ytrafo                     :function (data, numeric_trafo = id_trafo, factor_trafo = f_trafo, ordered_trafo = of_trafo, surv_trafo = logrank_trafo, 
    var_trafo = NULL, block = NULL)  
  .. .. ..@ x                          :'data.frame':	400 obs. of  1 variable:
  .. .. .. ..$ SNP_Group: Factor w/ 4 levels "Consensus","GWAS lead",..: 3 3 3 3 3 3 3 3 3 3 ...
  .. .. .. ..- attr(*, "terms")=Classes 'terms', 'formula'  language ~SNP_Group
  .. .. .. .. .. ..- attr(*, "variables")= language list(SNP_Group)
  .. .. .. .. .. ..- attr(*, "factors")= int [1, 1] 1
  .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. .. .. ..$ : chr "SNP_Group"
  .. .. .. .. .. .. .. ..$ : chr "SNP_Group"
  .. .. .. .. .. ..- attr(*, "term.labels")= chr "SNP_Group"
  .. .. .. .. .. ..- attr(*, "order")= int 1
  .. .. .. .. .. ..- attr(*, "intercept")= int 1
  .. .. .. .. .. ..- attr(*, "response")= int 0
  .. .. .. .. .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
  .. .. .. .. .. ..- attr(*, "predvars")= language list(SNP_Group)
  .. .. .. .. .. ..- attr(*, "dataClasses")= Named chr "factor"
  .. .. .. .. .. .. ..- attr(*, "names")= chr "SNP_Group"
  .. .. ..@ y                          :'data.frame':	400 obs. of  1 variable:
  .. .. .. ..$ p: num [1:400] 0.218 0.353 0.579 1 0.28 ...
  .. .. .. ..- attr(*, "terms")=Classes 'terms', 'formula'  language ~p
  .. .. .. .. .. ..- attr(*, "variables")= language list(p)
  .. .. .. .. .. ..- attr(*, "factors")= int [1, 1] 1
  .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. .. .. ..$ : chr "p"
  .. .. .. .. .. .. .. ..$ : chr "p"
  .. .. .. .. .. ..- attr(*, "term.labels")= chr "p"
  .. .. .. .. .. ..- attr(*, "order")= int 1
  .. .. .. .. .. ..- attr(*, "intercept")= int 1
  .. .. .. .. .. ..- attr(*, "response")= int 0
  .. .. .. .. .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
  .. .. .. .. .. ..- attr(*, "predvars")= language list(p)
  .. .. .. .. .. ..- attr(*, "dataClasses")= Named chr "numeric"
  .. .. .. .. .. .. ..- attr(*, "names")= chr "p"
  .. .. ..@ block                      : Factor w/ 1 level "0": 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..@ weights                    : num [1:400] 1 1 1 1 1 1 1 1 1 1 ...
  ..@ estimates   : list()
  ..@ method      : chr "General Independence Test"
  ..@ call        : language independence_test.IndependenceProblem(object = new("IndependenceProblem", x = list(SNP_Group = c(3L,  3L, 3L, 3L,| __truncated__ ...
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
data(coin_res)
## maybe str(coin_res) ; plot(coin_res) ...
}
\keyword{datasets}
