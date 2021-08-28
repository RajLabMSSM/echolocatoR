readSparse <- function(LD_path,
                       convert_to_df=F){
  LD_sparse <- readRDS(LD_path)
  if(convert_to_df) LD_sparse <- data.frame(as.matrix(LD_sparse))
  return(LD_sparse)
}

