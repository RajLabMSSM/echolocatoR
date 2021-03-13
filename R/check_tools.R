



check_tools <- function(check_tabix=T,
                        stop_for_tabix=F,
                        check_bcftools=T,
                        stop_for_bcftools=F,
                        conda_env=NULL){
  if(check_tabix){
    printer("Checking for tabix installation...")
    tabix <- CONDA.find_package("tabix", conda_env=conda_env, verbose = F)
    if(tabix=="tabix"){
      tabix.out <- system("which tabix",intern = T)
      if(length(tabix.out)==0){
        if(stop_for_tabix){
          stop("No tabix installation detected. ",
               "This will cause downstream errors. ",
               "To avoid this, please install htslib via source, brew or conda.",
               "For details, see: https://github.com/RajLabMSSM/echolocatoR/#command-line ")
        }else {
          warning("No tabix installation detected. ",
                  "This may cause downstream errors in some cases. ",
                  "To avoid this, please install htslib via command line or conda. ",
                  "For details, see: https://github.com/RajLabMSSM/echolocatoR/#command-line ")
        }

      }
    }
  }
    if(check_bcftools){
      printer("Checking for bcftools installation...")
      bcftools <- CONDA.find_package("bcftools", conda_env=conda_env, verbose = F)
      if(bcftools=="bcftools"){
        bcftools.out <- system("which bcftools",intern = T)
        if(length(bcftools.out)==0){
          if(check_bcftools){
            stop("No bcftools installation detected. ",
                 "This will cause downstream errors. ",
                 "To avoid this, please install bcftools via command line or conda. ",
                 "For details, see: https://github.com/RajLabMSSM/echolocatoR/#command-line ")
          }else {
            warning("No bcftools installation detected. ",
                    "This may cause downstream errors in some cases. ",
                    "To avoid this, please install bcftools via command line or conda. ",
                    "For details, see: https://github.com/RajLabMSSM/echolocatoR/#command-line ")
          }
        }
      }
  }

}
