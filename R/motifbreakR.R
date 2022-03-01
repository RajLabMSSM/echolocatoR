
# devtools::install_github("Simon-Coetzee/motifBreakR") **** USE THIS! BiocManager outdated
# BiocManager::install("motifbreakR")
# BiocManager::install(c("SNPlocs.Hsapiens.dbSNP142.GRCh37","BSgenome.Hsapiens.UCSC.hg19"))
# library(motifbreakR); library(BSgenome)



#' Filter by motif database metadata
#' @family motifbreakR
#' @keywords internal
MOTIFBREAKR.filter_by_metadata <- function(mb.results,
                                           Organism="Hsapiens"){
  meta <- subset(GenomicRanges::mcols(MotifDb::MotifDb), organism==Organism)
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
#' # rsid_list <- unique(subset(merged_DT, Locus %in% c("LRRK2","MED12L","DYRK1A","FCGR2A") & (Consensus_SNP | leadSNP))$SNP)
#' # rsid_list <- unique(subset(merged_DT, Consensus_SNP | leadSNP)$SNP)
#' rsid_list <- c("rs11175620","rs7294619","rs74324737")
#' mb.results <- MOTIFBREAKR(rsid_list=rsid_list, calculate_all_pval=TRUE, force_new = TRUE)
#' }
#' @export
MOTIFBREAKR <- function(rsid_list,
                        save_rds=TRUE,
                        dataset_dir="./results",
                        pwmList=NULL,
                        organism=NULL,#"Hsapiens",
                        # If `filterp=T`, this is a p-value threshold. If `filterp=F` it is the pct threshold.
                        threshold=.85, # 1e-4 #  4e-8
                        show.neutral=FALSE,
                        method="default",
                        verbose=TRUE,
                        calculate_all_pval=TRUE,
                        force_new=FALSE){
  ##  data("merged_DT"); snp_list <- unique(subset(merged_DT, Support>0 | leadSNP)$SNP);
  # library(echolocatoR); library(motifbreakR); library(BSgenome); save_rds=T; dataset_dir <- "./results/motifbreakR"; force_new=F;
  # pwmList=NULL; organism=NULL;  method = "default"; verbose=T; show.neutral = F; calculate_all_pval=T; threshold=1e-4;

  library(BSgenome)
  rds_path <- file.path(dataset_dir,'_genome_wide','motifbreakR','motifbreakR_results.rds');
  # rds_path <- file.path(dataset_dir,'_genome_wide','motifbreakR','motifbreakR_results.p_values.rds');

  if(!file.exists(rds_path) | force_new){
    # Prepare input
    messager("+ MOTIFBREAKR:: Turning SNP list into motifbreakR input format.", v=verbose)
    gr.snps <- motifbreakR::snps.from.rsid(rsid = rsid_list,
                                           dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37::SNPlocs.Hsapiens.dbSNP142.GRCh37,
                                           search.genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19);
    # Subset motif databases
    if(is.null(pwmList)){pwmList <- MotifDb::MotifDb}
    if(!is.null(organism)){
      pwmList <- pwmList[grep(paste0("^",organism),names(pwmList),value = TRUE)]
    }

    # Run motifbreakR
    messager("+ MOTIFBREAKR:: Identifying motifs and predicting disruptions.", v=verbose)
    mb.results <- motifbreakR::motifbreakR(snpList = gr.snps,
                                           pwmList = pwmList,

                                           filterp = TRUE,
                                           threshold = threshold,
                                           method = method,
                                           show.neutral = show.neutral,
                                           verbose = verbose);
    # Save tmp results
    if(save_rds){
      dir.create(dirname(tmp_path),showWarnings = FALSE, recursive = TRUE);
      messager("+ MOTIFBREAKR:: Saving tmp results ==>", tmp_path);
      saveRDS(mb.results, tmp_path);
    }
  } else {
    messager( "+ MOTIFBREAKR:: Using pre-existing tmp file.", v=verbose)
    mb.results <- readRDS(tmp_path)
  }

  return(mb.results)
}



MOTIFBREAK.calc_pvals <- function(mb.results,
                                  remove_NA_TF=TRUE,
                                  effect_strengths=c("strong"),
                                  filter_by_locus=NULL,
                                  dataset_dir="./results",
                                  verbose=TRUE){
  # Filter
  if(remove_NA_TF){
    mb.results <- subset(mb.results, !is.na(geneSymbol))
  }
  if(!is.null(effect_strengths)){
   mb.results <- subset(mb.results, effect %in% effect_strengths)
   messager("+ MOTIFBREAKR",nrow(mb.results),"SNPs in results @",
           paste0("`effect_strengths=",paste(effect_strengths, collapse=", "),"`"),v=verbose)
  }
  if(!is.null(filter_by_locus)){
    data("merged_DT")
    mb.results = subset(mb.results, SNP_id %in% subset(merged_DT, Locus==filter_by_locus)$SNP)
  }

  # Calculate p-values
  if(calculate_all_pval){
    messager("+ MOTIFBREAKR:: Calculating p-values for all SNPs...", v=verbose)
    mb.results_p <- motifbreakR::calculatePvalue(mb.results)
  }
  if(save_rds){
    rds_path <- file.path(dataset_dir,'_genome_wide','motifbreakR','motifbreakR_results.rds');
    messager("+ MOTIFBREAKR:: Saving results ==>",rds_path, v=verbose)
    saveRDS(mb.results_p, rds_path);
  }
  return(mb.results_p)
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
#' \dontrun{
#' data("merged_DT")
#' microglia_TF <- read.csv("~/Desktop/Fine_Mapping/resources/microglia_TF.csv")
#'
#' root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/motifbreakR"
#' # mb.results <- readRDS(file.path(root, "motifbreakR_results.rds"))
#' mb.results <- readRDS("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/motifbreakR/motifbreakR_results.rds")
#'
#' mb.lrrk2 <- readRDS(file.path(root, "motifbreakR_results.p_values_LRRK2.rds"))
#' mb.encode <- readRDS(file.path(root, "motifbreakR_results.encode.lrrk2.rds"))
#' mb.DYRK1A_FCGR2A <- readRDS(file.path(root,"mb.results_p.DYRK1A_FCGR2A.RDS"))
#' mb.MED12L<- readRDS(file.path(root,"MED12L.pvalues.RDS"))
#'
#' mb.sub <- subset(mb.results, SNP_id %in% subset(merged_DT, Locus=="DNAH17" & (leadSNP | Consensus_SNP))$SNP)
#' mb.DNAH17 <- motifbreakR::calculatePvalue(results=mb.sub); saveRDS(mb.DNAH17, file.path(root,"DNAH17.pvalues.RDS"))
#' mb.DNAH17 <- readRDS(file.path(root,"DNAH17.pvalues.RDS"))
#'
#' mb.MBNL2 <- readRDS(file.path(root, "motifbreakR_results.p_values_MBNL2.rds"))
#'
#' MOTIFBREAKR.filter(merged_DT, mb.lrrk2, pct_threshold=NULL, effect_strengths=NULL)
#' }
#' @export
MOTIFBREAKR.filter <- function(merged_DT,
                                mb.results,
                                pct_threshold=NULL,
                                pvalue_threshold=1e-4,
                                qvalue_threshold=.05,
                                effect_strengths=c("strong"),
                                snp_filter="Consensus_SNP==T",
                                top_TF_hits=FALSE,
                                no_no_loci=NULL,
                                verbose=TRUE){
  # Quickstart
  # pct_threshold=.7; verbose=T;
  # pvalue_threshold=1e-4;
  # effect_strengths=c("strong");
  # snp_filter="Consensus_SNP==T";
  # no_no_loci<- c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1",
  #                "RETREG3","UBTF","FAM171A2","MAP3K14","CRHR1","MAPT-AS1","KANSL1","NSF","WNT3");


  # mb.results <- MOTIFBREAKR.filter_by_metadata(mb.results = mb.results,
  #                                              Organism = "Hsapiens")
  library(BSgenome); library(BSgenome.Hsapiens.UCSC.hg19); library(motifbreakR);
  mb.results$SNP <- names(mb.results)

  # if(!is.null(no_no_loci)) {
  #   merged_DT <- subset(merged_DT, !(Locus %in% no_no_loci))
  # }
  mb.results <- subset(mb.results, !is.na(geneSymbol) )

  mb.merge <- data.table::merge.data.table(x = merged_DT,
                                y = data.table::as.data.table(mb.results),
                                by="SNP", all = TRUE)  %>%
    dplyr::mutate(risk_allele=ifelse(A1==REF,"REF",ifelse(A1==ALT,"ALT",NA))) %>%
    subset(!is.na(risk_allele)) %>%
    dplyr::mutate(risk_pct=ifelse(risk_allele=="REF", pctRef, pctAlt),
                  # nonrisk_pct=ifelse(risk_allele=="ALT",  pctRef, pctAlt),

                  risk_score=ifelse(risk_allele=="REF", scoreRef, scoreAlt),
                  # nonrisk_score=ifelse(risk_allele=="ALT", scoreRef, scoreAlt),

                  risk_pvalue=ifelse(risk_allele=="REF", Refpvalue, Altpvalue),
                  # nonrisk_pvalue=ifelse(risk_allele=="ALT", scoreRef, scoreAlt)
                  ) %>%
    dplyr::mutate(risk_qvalue = stats::p.adjust(p = risk_pvalue, method = "bonferroni"))
  messager("+ MOTIFBREAKR",nrow(mb.merge),"SNPs in results.", v=verbose)

  if(!is.null(snp_filter)){
    mb.merge <- subset(mb.merge, eval(parse(text=snp_filter)))
    messager("+ MOTIFBREAKR",nrow(mb.merge),"SNPs in results @",
            paste0("`snp_filter='",snp_filter,"'`"),v=verbose)
  }

  if(!is.null(effect_strengths)){
    mb.merge <- subset(mb.merge, effect %in% effect_strengths)
    messager("+ MOTIFBREAKR",nrow(mb.merge),"SNPs in results @",
            paste0("`effect_strengths=",paste(effect_strengths, collapse=", "),"`"),v=verbose)
  }

  if(!is.null(pvalue_threshold)){
    mb.merge <- subset(mb.merge, risk_pvalue < pvalue_threshold)
    messager("+ MOTIFBREAKR",nrow(mb.merge),"SNPs in results @",
            paste0("`pvalue_threshold=",pvalue_threshold,"`"),v=verbose)
  }

  if(!is.null(qvalue_threshold)){
    mb.merge <- subset(mb.merge, risk_qvalue < qvalue_threshold)
    messager("+ MOTIFBREAKR",nrow(mb.merge),"SNPs in results @",
            paste0("`qvalue_threshold=",qvalue_threshold,"`"),v=verbose)
  }

  if(!is.null(pct_threshold)){
    # Filter where pval is sig AND it's the effect allele in the GWAS
    mb.merge <-  subset(mb.merge,
                        # (pctRef<pct_threshold & A1==REF) |
                        # (pctAlt<pct_threshold & A1==ALT)
                        risk_pct < pct_threshold
                        )
    messager("+ MOTIFBREAKR",nrow(mb.merge),"SNPs in results @",
            paste0("`pct_threshold=",pct_threshold,"`"),v=verbose)
  }
  if(top_TF_hits){
    mb.merge <- mb.merge %>%
      dplyr::group_by(Locus, SNP, geneSymbol) %>%
      dplyr::arrange(risk_pvalue, dplyr::desc(risk_score)) %>%
      dplyr::slice(1) %>% data.table::data.table()
  }
  # select_cols <- c("Locus","SNP","Consensus_SNP","geneSymbol","dataSource","seqMatch","risk_pct","risk_score","risk_pvalue")
  # if(!is.null(save_path)){
  #   data.table::fwrite(data.frame(mb.merge), "~/Downloads/motifbreakR_sig_rs6781790-rs759905.csv")
  # }
  return(mb.merge)
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
#' # Example from motifbreakR
#' library(motifbreakR)
#' data("example.results")
#' library(echolocatoR)
#' data("merged_DT")
#' root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/motifbreakR"
#'
#' mb.results<- readRDS(file.path(root,"DNAH17.pvalues.RDS"))
#'
#' # Get overlap with PEAKS
#' # PEAKS <- NOTT_2019.get_epigenomic_peaks(convert_to_GRanges = TRUE)
#' REGIONS <- NOTT_2019.get_regulatory_regions(as.granges = TRUE)
#' gr.hits <- GRanges_overlap(dat1 = subset(merged_DT, Consensus_SNP), chrom_col.1 = "CHR", start_col.1 = "POS", end_col.1 = "POS", dat2 = REGIONS)
#' mb.filter <- MOTIFBREAKR.filter(merged_DT = merged_DT, mb.results = mb.results, pct_threshold=NULL)
#'
#' subset(mb.filter, SNP_id %in% gr.hits$SNP)
#' plot_paths <- MOTIFBREAKR.plot(mb.results=mb.results, mb.filter=mb.filter, save_dir="~/Desktop")
#' }
#' @export
MOTIFBREAKR.plot <- function(mb.results,
                             mb.filter=NULL,
                             rsid=NULL,
                             effect=c("strong","weak"),
                             save_dir=NULL,
                             height=3,
                             width=7){
  library(BSgenome); library(BSgenome.Hsapiens.UCSC.hg19); library(motifbreakR);
  mb.results <- MOTIFBREAKR.make_id(mb.results = mb.results)
  if(!is.null(mb.filter)){
    mb.filter <- MOTIFBREAKR.make_id(mb.results = mb.filter)
    mb.results <- subset(mb.results, id %in% mb.filter$id)
  }
  if(is.null(rsid)) rsid <- unique(mb.results$SNP_id)
  plot_paths <- c()
  # par(mfrow=c(2,2))
  for(rs in rsid){
    save_path <- file.path(save_dir,paste0(rs,'.png'))
    plot_paths <- append(plot_paths, save_path)
    if(!is.null(save_dir)) png(filename = save_path,
                               height = height,
                               width = width)
    motifbreakR::plotMB(results = mb.results,
                        rsid = rs,
                        effect = effect)
    dev.off();
  }
  return(plot_paths)
}





MOTIFBREAKR.make_id <- function(mb.results){
  mb.results$id <- paste(mb.results$SNP_id, mb.results$dataSource, mb.results$providerId, sep='_')
  return(mb.results)
}





MOTIFBREAKR.summarize <- function(){

  # Tally hits hits per
  db_tally <- mb.merge %>%
    dplyr::group_by(effect, dataSource) %>%
    dplyr::tally() %>%
    dplyr::arrange(effect, desc(n)) %>%
    data.frame()
  print(db_tally)

  top_snps <- mb.merge %>%
    dplyr::group_by(Locus) %>%
    # dplyr::arrange(desc(risk_score), risk_pct) %>%
    dplyr::arrange(desc(alleleDiff)) %>%
    dplyr::slice(1)
  top_rsids <- unique(top_snps$SNP)

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
                     # all_disrupting_SNPs=paste(unique(SNP), collapse = '; '),
                     # all_TFs=paste(unique(geneSymbol), collapse='; '),
                     # all_sequences=paste(unique(gsub(" ","",seqMatch)), collapse='; ')
    ) %>%
    dplyr::mutate(top_disrupting_SNP_is_lead=top_disrupting_SNP %in% unique(subset(mb.merge, leadSNP)$SNP),
                  top_disrupting_SNP_in_UCS=top_disrupting_SNP %in% unique(subset(mb.merge, Support>0)$SNP),
                  top_disrupting_SNP_in_consensus=top_disrupting_SNP %in% unique(subset(mb.merge, Consensus_SNP)$SNP)
    ) %>%
    data.frame() %>% unique()

  # if(save_path!=FALSE){
  #   data.table::fwrite(locus_tally, "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/motifbreakR/motifbreakR_locus_tally_ALL.csv", sep=",")
  # }
  return(locus_tally)
}



MOTIFBREAKR.check_TF_overlap <- function(mb.merge,
                                         mb.encode,
                                         pvalue_threshold=1e-4){

  dat_merged <- base::merge(subset(mb.merge, risk_pvalue<pvalue_threshold),
                            data.frame(subset(mb.encode, Refpvalue<pvalue_threshold | Altpvalue<pvalue_threshold)),
                            by=c("SNP_id","geneSymbol"))
  if(all(dat_merged$risk_allele=="ALT")){
    dat_merged <- subset(dat_merged, Altpvalue.y<pvalue_threshold)
  } else {
    dat_merged <- subset(dat_merged, Refpvalue.y<pvalue_threshold)
  }

  TF_overlap <- unique(dat_merged$geneSymbol)
  return(TF_overlap)
}
