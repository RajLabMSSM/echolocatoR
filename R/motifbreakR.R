
# devtools::install_github("Simon-Coetzee/motifBreakR") **** USE THIS! BiocManager outdated
# BiocManager::install("motifbreakR")
# BiocManager::install(c("SNPlocs.Hsapiens.dbSNP142.GRCh37","BSgenome.Hsapiens.UCSC.hg19"))
# library(motifbreakR); library(BSgenome)



#' Filter by motif database metadata
#' @family motifbreakR
#' @keywords internal
MOTIFBREAKR.filter_by_metadata <- function(mb.results,
                                           Organism="Hsapiens"){
  meta <- subset(mcols(MotifDb::MotifDb), organism==Organism)
  mb.filtered <-subset(mb.results, providerId %in% meta$providerId)
  return(mb.filtered)
}




#' Run \code{\link{motifbreakR}}
#'
#' \code{\link{motifbreakR}} is a package to predict how much a SNP will disrupt
#' a transcription factor binding motif (if it falls witihn one).
#' @inheritParams motifbreakR::snps.from.rsid
#' @inheritParams motifbreakR::motifbreakR
#' @return motifbreakr results
#' @family motifbreakR
#' @source
#' \strong{Publication:}
#' \url{https://pubmed.ncbi.nlm.nih.gov/26272984/}
#' \strong{GitHub:}
#' \url{https://github.com/Simon-Coetzee/MotifBreakR}
#' \strong{Vignette:}
#' \url{http://simon-coetzee.github.io/motifBreakR}
#' @examples
#' \dontrun{
#' # data("merged_DT")
#' # snp_list <- unique(subset(merged_DT, Support>0 | leadSNP)$SNP)
#' snp_list <- "rs1006140"
#' mb.results <- MOTIFBREAKR(snp_list=snp_list)
#' }
MOTIFBREAKR <- function(snp_list,
                        save_rds=T,
                        dataset_dir,
                        pwmList=NULL,
                        organism="Hsapiens",
                        threshold=.85,
                        show.neutral=F,
                        method="default",
                        verbose=T,
                        calculate_all_pval=F){
  # library(echolocatoR); library(motifbreakR); library(BSgenome); save_rds=T; dataset_dir <- "./results/GWAS/Nalls23andMe_2019";
  # pwmList=NULL; organism="Hsapiens";  threshold=.85; method = "default"; verbose=T; show.neutral = F; calculate_all_pval=T; save_rds=T;
  # data("merged_DT"); snp_list <- unique(subset(merged_DT, Support>0 | leadSNP)$SNP);

  library(BSgenome)
  # Prepare input
  printer("+ MOTIFBREAKR:: Turning SNP list into motifbreakR input format.", v=verbose)
  variants <- motifbreakR::snps.from.rsid(rsid = snp_list,
                                          dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37::SNPlocs.Hsapiens.dbSNP142.GRCh37,
                                          search.genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19);
  # Subset motif databases
  if(is.null(pwmList)){pwmList <- MotifDb::MotifDb}
  select_dbs <- grep(paste0("^",organism),names(pwmList),value = T)
  pwmList <- pwmList[select_dbs]
  # Run motifbreakR
  printer("+ MOTIFBREAKR:: Identifying motifs and predicting disruptions.", v=verbose)
  mb.results <- motifbreakR::motifbreakR(snpList = variants,
                                         pwmList = pwmList,

                                         filterp = T,
                                         threshold = threshold,
                                         method = method,
                                         show.neutral = show.neutral,
                                         verbose = verbose);
  # Calculate p-values
  if(calculate_all_pval){
    printer("+ MOTIFBREAKR:: Calculating p-values for all SNPs...", v=verbose)
    mb.results <- motifbreakR::calculatePvalue(mb.results)
  }
  # Save results
  if(save_rds){
    rds_path <- file.path(dataset_dir,'_genome_wide','motifbreakR','motifbreakR_results.rds');
    printer("+ MOTIFBREAKR:: Saving reults ==>",rds_path, v=verbose)
    dir.create(dirname(rds_path),showWarnings = F, recursive = T);
    # printer("+ MOTIFBREAKR:: Saving results ==>", rds_path);
    saveRDS(mb.results, rds_path);
  }
  return(mb.results)
}




#' Plot \code{\link{motifbreakR}} results
#'
#' @source
#' \strong{Publication:}
#' \url{https://pubmed.ncbi.nlm.nih.gov/26272984/}
#'
#' \strong{GitHub:}
#' \url{https://github.com/Simon-Coetzee/MotifBreakR}
#'
#' @examples
#' \dontrun{
#' # mb.results <- readRDS("/sc/arion/projects/pd-omics/brian/motifbreakR/motifbreakR_results.rds")
#' mb.results <- readRDS("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/motifbreakR/motifbreakR_results.rds")
#' MOTIFBREAKR.plot(mb.results=mb.results, rsid="rs114528427")
#'
#' # Example from motifbreakR
#' library(motifbreakR)
#' data("example.results")
#' motifbreakR::plotMB(mb.results=example.results, "rs2661839", effect = "strong")
#' }
MOTIFBREAKR.plot <- function(mb.results,
                             rsid=NULL,
                             effect=c("strong","weak")){
  library(BSgenome); library(BSgenome.Hsapiens.UCSC.hg19); library(motifbreakR);
  printer("+ MOTIFBREAKR::",length(subset(mb.results, effect=="strong")),"strong effects detected.")
  printer("+ MOTIFBREAKR::",length(subset(mb.results, effect=="weak")),"weak effects detected.")

  # rsid<-"rs7294619"
  if(is.null(rsid)){rsid <- names(mb.results)[20]}
  motifbreakR::plotMB(results = mb.results[1:10],
                      rsid = names(mb.results)[5],
                      effect = c("strong")
                      )
}


#' Summarise \code{\link{motifbreakR}} + \code{\link{echolocatoR}} results
#'
#' For each SNP we have at least one allele achieving a p-value below 1e-4 threshold that we required.
#' The seqMatch column shows what the reference genome sequence is at that location,
#' with the variant position appearing in an uppercase letter.
#' pctRef and pctAlt display the the score for the motif in the sequence
#' as a percentage of the best score that motif could achieve on an ideal sequence.
#' In other words (scoreVariant−minscorePWM)/(maxscorePWM−minscorePWM).
#' We can also see the absolute scores for our method in scoreRef and scoreAlt
#' and thier respective p-values.
#' @examples
#' data("merged_DT")
MOTIFBREAKR.summarize <- function(merged_DT,
                                  mb.results,
                                  no_no_loci=NULL){
  # no_no_loci<- c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1",
  #                "RETREG3","UBTF","FAM171A2","MAP3K14","CRHR1","MAPT-AS1","KANSL1","NSF","WNT3")
  mb.results <- MOTIFBREAKR.filter_by_metadata(mb.results = mb.results,
                                               Organism = "Hsapiens")
  mb.results$SNP <- names(mb.results)
  pct_threshold <-.7
  mb.merge <- data.table::merge.data.table(merged_DT,
                               data.table::as.data.table(mb.results),
                               by="SNP") %>%
    subset(!(Locus %in% no_no_loci)) %>%
    dplyr::mutate(risk_allele=ifelse(A1==REF,"REF",ifelse(A1==ALT,"ALT",NA))) %>%
    dplyr::mutate(risk_pct=ifelse(risk_allele=="REF", pctRef, pctAlt),
                  nonrisk_pct=ifelse(risk_allele=="REF", pctAlt,pctRef)) %>%
    # Filter where pval is sig AND it's the effect allele in the GWAS
    subset((pctRef<pct_threshold & A1==REF) |
           (pctAlt<pct_threshold & A1==ALT))

  # Tally hits hits per
  db_tally <- mb.merge %>%
    dplyr::group_by(effect, dataSource) %>%
    dplyr::tally() %>%
    dplyr::arrange(effect, desc(n)) %>%
    data.frame()
  print(db_tally)

  top_rsids <- (mb.merge %>%
    dplyr::group_by(Locus) %>%
    dplyr::arrange(risk_pct) %>%
    dplyr::slice(1))$SNP %>%
    unique()

  locus_tally <- mb.merge %>%
    dplyr::group_by(Locus) %>%
    dplyr::summarise(n_lead=dplyr::n_distinct(SNP[leadSNP]),
                     n_UCS=dplyr::n_distinct(SNP[Support>0]),
                     n_consensus=dplyr::n_distinct(SNP[Consensus_SNP]),
                     lead_in_consensus=SNP[leadSNP] %in% SNP[Consensus_SNP],
                     # consensus_SNPs=paste(unique(SNP[Consensus_SNP]), collapse = ", "),
                     top_disrupting_SNP=paste(unique(SNP[SNP%in%top_rsids]), collapse = '; '),
                     top_TF=paste(unique(geneSymbol[SNP%in%top_rsids]), collapse='; '),
                     top_sequence=paste(unique(gsub(" ","",seqMatch[SNP%in%top_rsids])), collapse='; '),
                     # consensus_RefAlt=paste(unique(.[Consensus_SNP,c("SNP","risk_allele")])$risk_allele, collapse=",")
                     all_disrupting_SNPs=paste(unique(SNP), collapse = '; '),
                     all_TFs=paste(unique(geneSymbol), collapse='; '),
                     all_sequences=paste(unique(gsub(" ","",seqMatch)), collapse='; ')
                     ) %>%
    dplyr::mutate(top_disrupting_SNP_is_lead=top_disrupting_SNP %in% unique(subset(mb.merge, leadSNP)$SNP),
                  top_disrupting_SNP_in_UCS=top_disrupting_SNP %in% unique(subset(mb.merge, Support>0)$SNP),
                  top_disrupting_SNP_in_consensus=top_disrupting_SNP %in% unique(subset(mb.merge, Consensus_SNP)$SNP)
) %>%
    data.frame() %>% unique()
  data.table::fwrite(locus_tally, "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/motifbreakR/motifbreakR_locus_tally_ALL.csv", sep=",")

 return(locus_tally)
}




#' \code{\link{motifbreakR}} summary plot
#'
#'

