#' Remove old vcf files
#'
#' @examples
#' vcf_cleaning(root="Data/GWAS/Ripke_2014", LD_ref="1KGphase3", loci=top_SNPs$Locus)
vcf_cleaning <- function(root,
                         LD_ref="1KGphase3",
                         loci=NULL,
                         verbose=T){
  files <- list.files(root, paste0("*",LD_ref,".vcf.gz","*"),
                      recursive = T, full.names = T)
  locus_names <- unique(basename(dirname(dirname(files))))
  select_loci <- dplyr::intersect(loci, locus_names)
  if(length(select_loci)>0){
    files <- files[locus_names %in% select_loci]
  }
  printer("Removing",length(files),"files",v=verbose)
  out <- file.remove(files)
}
