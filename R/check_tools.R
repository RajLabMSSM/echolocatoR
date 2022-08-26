check_tools <- function(check_tabix=FALSE,
                        stop_for_tabix=FALSE,
                        check_bcftools=FALSE,
                        stop_for_bcftools=FALSE,
                        conda_env=NULL){
  if(check_tabix){
    messager("Checking for tabix installation...")
    tabix <- echoconda::find_packages("tabix", conda_env=conda_env, verbose  = FALSE)
    if(tabix=="tabix"){
      tabix.out <- system("which tabix",intern = TRUE)
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
      messager("Checking for bcftools installation...")
      bcftools <- echoconda::find_packages("bcftools", conda_env=conda_env, verbose  = FALSE)
      if(bcftools=="bcftools"){
        bcftools.out <- system("which bcftools",intern = TRUE)
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
