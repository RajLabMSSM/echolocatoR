

#' Import top GWAS/QTL summary statistics
#'
#' The resulting top_SNPs data.frame can be used to guide
#' the \code{\link{finemap_loci}} in querying and fine-mapping loci.
#'
#' @param top_SNPs Can be a data.frame with the top sumary stats per locus.
#' Alternatively, you can provide a path to the stored top summary stats file.
#' Can be in any tabular format (e.g. excel, .tsv, .csv, etc.).
#' This file should have one lead GWAS/QTL hits per locus.
#' If there is more than one SNP per locus, the one with the smallest p-value (then the largest effect size) is selected as the lead SNP.
#' The lead SNP will be used as the center of the locus when constructing the locus subset files.
#'
#' @family query functions
#' @param show_table Create an interative data table.
#' @param sheet If the \emph{topSS} file is an excel sheet, you can specify which tab to use.
#' You can provide either a number to identify the tab by order,
#' or a string to identify the tab by name.
#' @param gene_col An optional column to keep track of the causal gene(s) in each locus (if known).
#' @param grouping_vars The variables that you want to group by
#' such that each grouping_var combination has its own index SNP.
#' For example, if you want one index SNP per QTL eGene - GWAS locus pair, you could supply:
#' \code{grouping_vars=c("Locus","Gene")}.
#' @inheritParams finemap_pipeline
#' @export
import_topSNPs <- function(topSS,
                           show_table=T,
                           sheet = 1,
                           munge=T,
                           ref_genome=NULL,
                           chrom_col=NULL,
                           position_col=NULL,
                           min_POS_col=NULL,
                           max_POS_col=NULL,
                           snp_col=NULL,
                           pval_col=NULL,
                           effect_col=NULL,
                           locus_col="Locus",
                           grouping_vars=c("Locus"),
                           gene_col="Gene",
                           remove_variants=NULL,
                           nThread=4,
                           verbose=T){
  # Reassign colnames bc read.xlsx won't let you prevent this in some versions....
  if(!is.data.frame(topSS)){
    if(endsWith(topSS, ".xlsx") | endsWith(topSS, ".xlsm")){
      locus_col=gsub(" ",".",trimws(locus_col))
      gene_col=gsub(" ",".",trimws(gene_col))
      chrom_col=gsub(" ",".",trimws(chrom_col))
      position_col=gsub(" ",".",trimws(position_col))
      snp_col=gsub(" ",".",trimws(snp_col))
      pval_col=gsub(" ",".",trimws(pval_col))
      effect_col=gsub(" ",".",trimws(effect_col))
    }
  }

  #### Import top SNPs ####
  top_SNPs <- topSNPs_reader(topSS, sheet)
  orig_top_SNPs <- top_SNPs

  top_SNPs <- standardise_gene_locus_cols(top_SNPs=top_SNPs,
                                          orig_top_SNPs=orig_top_SNPs,
                                          locus_col=locus_col,
                                          gene_col=gene_col)

  if(munge){
    #### Check for necessary cols ####
    top_SNPs <- MUNGESUMSTATS.check_syn(dat = top_SNPs, col_name=pval_col, col_type = "P")
    top_SNPs <- MUNGESUMSTATS.check_syn(top_SNPs, col_name=effect_col, col_type = "BETA")
    printer("+ Munging top_SNPs",v=verbose)
    top_SNPs <- suppressMessages( MungeSumstats:::standardise_sumstats_column_headers_crossplatform(sumstats_dt = top_SNPs)$sumstats_dt)
    # top_SNPs <- MUNGESUMSTATS.to_echolocatoR(top_SNPs)
    map <- column_map(package = "MungeSumstats")
    chrom_col <- map$chrom_col; position_col <- map$position_col; snp_col <- map$snp_col;
    pval_col <- map$pval_col; effect_col <- map$effect_col;
    top_SNPs <- dplyr::rename(top_SNPs,
                              Locus="LOCUS",
                              Gene="GENE")
  }



  # Fill in missing cols with nonsense
  if(is.null(chrom_col)) {top_SNPs$CHR <- NA; chrom_col<-"CHR";}
  if(is.null(position_col)) {top_SNPs$POS <- NA; position_col<-"POS";}
  if(!effect_col %in% colnames(top_SNPs)) {top_SNPs$Effect <- 1; effect_col<-"Effect";}

  top_SNPs <- dplyr::rename(top_SNPs,
                            CHR=all_of(chrom_col),
                            POS=all_of(position_col),
                            SNP=all_of(snp_col)) %>%
    dplyr::mutate(CHR=gsub("chr","",CHR))



  top_SNPs <- suppressMessages(top_SNPs %>%
                                 dplyr::select(Locus=Locus,
                                               Gene=Gene,
                                               CHR=CHR,
                                               POS=POS,
                                               SNP=SNP,
                                               P=pval_col,
                                               Effect=effect_col,
                                               min_POS=min_POS_col,
                                               max_POS=max_POS_col))
  if("min_POS" %in% colnames(top_SNPs) &
     "max_POS" %in% colnames(top_SNPs)){
    top_SNPs <- dplyr::mutate(top_SNPs, span_kb=(max_POS-min_POS)/1000)
  }

  # Remove specific variants
  if(!is.null(remove_variants)){
    top_SNPs <- subset(top_SNPs, !(SNP %in% remove_variants))
  }

  # Get the top representative SNP and Gene per locus (by lowest p-value and effect size)
  if(!is.null(grouping_vars)){
    top_SNPs <- suppressWarnings(top_SNPs %>%
                                   dplyr::arrange(P, dplyr::desc(Effect)) %>%
                                   dplyr::group_by(.dots=grouping_vars) %>%
                                   dplyr::slice(1) %>%
                                   replace(., .=="NA", NA) %>%
                                   subset(!is.na(Locus)) %>%
                                   dplyr::mutate(CHR=as.numeric(gsub("chr", "",CHR))) )
  }

  # Make sure cols are numeric
  top_SNPs <- top_SNPs %>% dplyr::mutate_at(.vars = vars(POS,P,Effect),
                                            .funs = function(x){as.numeric(gsub(",| ","",x))})
  if(show_table){
    createDT(top_SNPs, caption = "Top SNP per locus")
  }
  return(data.table::data.table(top_SNPs))
}



topSNPs_reader <- function(topSS,
                           sheet = 1){
  if(is.data.frame(topSS)){
    return(data.frame(topSS, check.names = FALSE))
  } else{
    if(endsWith(topSS, ".xlsx") | endsWith(topSS, ".xlsm")){
      topSS <- openxlsx::read.xlsx(topSS,
                                   sheet = sheet,
                                   # sep.names=" ", # Only in some versions?
                                   check.names = F)
    } else {
      topSS <- data.table::fread(file=topSS,
                                 header = T,
                                 stringsAsFactors = F,
                                 nThread=nThread)
    }
    return(topSS)
  }
}



standardise_gene_locus_cols <- function(top_SNPs,
                                        orig_top_SNPs,
                                        locus_col=NULL,
                                        gene_col=NULL,
                                        verbose=TRUE){
  # Add Locus/Gene columns
  bad_chars <- "/|[:]|[,]"
  if((is.null(gene_col) & is.null(locus_col))){
    printer("+ Constructing locus names from CHR and index SNP",v=verbose)
    # locus_chr1_rs10737496
    top_SNPs <- dplyr::mutate(top_SNPs,
                              Locus=paste0("locus_chr",CHR,"_",SNP),
                              Gene=paste0("locus_chr",CHR,"_",SNP))
    locus_col <- gene_col <- "Locus";
  } else if(!any(gene_col %in% colnames(top_SNPs)) & !any(locus_col %in% colnames(top_SNPs))){
    printer("+ Constructing locus names from CHR and index SNP",v=verbose)
    # locus_chr1_rs10737496
    top_SNPs <- dplyr::mutate(top_SNPs,
                              Locus=paste0("locus_chr",CHR,"_",SNP),
                              Gene=paste0("locus_chr",CHR,"_",SNP))
    locus_col <- gene_col <- "Locus";
  } else if(gene_col %in% colnames(top_SNPs) & locus_col %in% colnames(top_SNPs)){
    print("+ Assigning gene_col and locus_col independently",v=verbose)
    top_SNPs$Gene <- gsub(bad_chars,"_", orig_top_SNPs[[gene_col]])
    top_SNPs$Locus <- gsub(bad_chars,"_",orig_top_SNPs[[locus_col]]) # Get rid of problematic characters

  } else {
    if(gene_col %in% colnames(top_SNPs) & all(!is.na(top_SNPs[[gene_col]])) ){
      print("+ Assigning gene_col to locus_col",v=verbose)
      top_SNPs$Gene <- gsub(bad_chars,"_", orig_top_SNPs[[gene_col]])
      top_SNPs$Locus <- gsub(bad_chars,"_",top_SNPs$Gene) # Get rid of problematic characters
    } else {
      print("+ Assigning locus_col to gene_col",v=verbose)
      top_SNPs$Locus <- gsub(bad_chars,"_",orig_top_SNPs[[locus_col]])
      top_SNPs$Gene <- gsub(bad_chars,"_",top_SNPs$Locus) # Get rid of problematic characters
    }
  }
  return(top_SNPs)
}

