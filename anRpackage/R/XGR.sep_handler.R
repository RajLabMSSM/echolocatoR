XGR.sep_handler <-
function(lib.name){
  # "_(?=[^_]+$)" : Split by the last "_"
  sepDict <- list("ENCODE_TFBS_ClusteredV3_CellTypes"="[.]",
                  "ENCODE_DNaseI_ClusteredV3_CellTypes"="_(?=[^_]+$)",
                  "Broad_Histone"="_(?=[^_]+$)",
                  "FANTOM5_Enhancer"="_(?=[^_]+$)",
                  "TFBS_Conserved"="[$]")
  if(lib.name %in% names(sepDict)){
    sep <- sepDict[[lib.name]]
  }else{ sep <- "_(?=[^_]+$)"}
  return(sep)
}
