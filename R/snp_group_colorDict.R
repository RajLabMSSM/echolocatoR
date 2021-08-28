snp_group_colorDict <- function(invert=F){
  colorDict <- c(
    "Random"="grey",
    "All"="red4",
    "GWAS nom. sig."="red3",
    "GWAS sig."="red2",
    "GWAS lead"="red",

    "ABF CS"="springgreen4",
    "SUSIE CS"="springgreen3",
    "POLYFUN-SUSIE CS"="springgreen2",
    "FINEMAP CS"="springgreen",
    "UCS (-PolyFun)"="green4",
    "UCS"="green2",

    "Support==0"="darkorchid4",
    "Support==1"="darkorchid3",
    "Support==2"="darkorchid2",
    "Support==3"="darkorchid1",
    "Support==4"="darkorchid",

    "Consensus (-PolyFun)"="goldenrod4",
    "Consensus"="goldenrod2")
  if(invert) colorDict <- setNames(names(colorDict), unname(colorDict))
  return(colorDict)
}
