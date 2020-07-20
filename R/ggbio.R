# ----------------------- #
# ----- ggbio plots ------#
# ----------------------- #



#' Multi-fine-map summary dot plot
#'
#' @examples
#' data("BST1")
#' gp <- GGBIO.dot_summary(finemap_dat=BST1)
GGBIO.dot_summary <- function(finemap_dat,
                              PP_threshold=.95,
                              show_plot=T){
  library(data.table)
  snp.labs <- construct_SNPs_labels(subset_DT = finemap_dat,
                                    remove_duplicates = F) %>%
    dplyr::arrange(CHR,POS) %>%
    dplyr::mutate(SNP=factor(SNP, levels = unique(SNP), ordered = T))
  CS_cols <- grep(".CS$",colnames(snp.labs), value = T)
  PP_cols <- grep(".PP$",colnames(snp.labs), value = T)
  methodDict <- setNames(gsub("\\.PP","",PP_cols),1:length(PP_cols))

  snp.melt <- suppressWarnings(
    data.table::melt.data.table(data = data.table::as.data.table(snp.labs),
                                id.vars = c("SNP","CHR","POS","leadSNP","Consensus_SNP","color","size","shape"),
                                measure.vars =  patterns(CS_group=".CS$", PP=".PP$"),
                                variable.name = c("method")) %>%
      dplyr::mutate(method=factor(methodDict[method], levels = rev(unname(methodDict)), ordered = T),
                    CS=ifelse(CS_group>0,T,NA),
                    PP=ifelse(PP>PP_threshold,PP,NA))
  )
  gp <- ggplot() +
    # Lead SNP
    geom_point(data = subset(snp.melt, leadSNP), aes(x=SNP, y=method),
               color="red", shape=18, size=2, stroke=1) +
    # CS
    geom_point(data = subset(snp.melt, !is.na(CS)), aes(x=SNP, y=method),
               color="green3", shape=1, size=3, stroke=1) +
    # Consensus
    geom_point(data = subset(snp.melt, Consensus_SNP), aes(x=SNP, y=method),
               color="darkgoldenrod1", shape=2, size=5, stroke=1) +
    scale_y_discrete(position = "right") +
    theme_bw()
  if(show_plot)print(gp)
  return(gp)
}




#' Make legend invisible
#'
#' @family plot
#' @keywords internal
GGBIO.invisible_legend <- function(gg){
  gg <- gg + theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  ) +
    scale_color_continuous(
      guide = guide_legend(override.aes = list(color = "white"))
    ) +
    scale_color_continuous(
      guide = guide_legend(override.aes = list(color = "white", fill="white"))
    )
  return(gg)
}




#' Track plot for SNPs
#'
#' Uses \code{\link{ggbio}}.
#'
#' @param gr.snp \emph{echolocator} fine-mapping results
#' in \code{\link{GRanges}} format.
#' @param method Which method to plot.
#' @param labels_subset Which SNP groups to label.
#' @param r2 Whether color the points by LD r2.
#' @param show.legend Whether to show legend.
#' @inheritParams finemap_pipeline
#' @family plot
#' @keywords internal
#' @examples
#' data("BST1");
#' gr.snp <- DT_to_GRanges(dat)
#' snp.track <- GGBIO.SNP_track(gr.snp=gr.snp, method = "original", color_r2=F)
GGBIO.SNP_track <- function(gr.snp,
                            method = "original",
                            labels_subset = c("Lead SNP", "Credible Set", "Consensus SNP"),
                            color_r2=F,
                            show.legend=T,
                            PP_threshold=.95,
                            sig_cutoff=5e-8){
  if(color_r2==F){gr.snp$r2 <- NA}
  # Format data
  if(!(method %in% c("original"))){
    if(method=="COJO"){
      GenomicRanges::mcols(gr.snp)[,"PP"] <- GenomicRanges::mcols(gr.snp)[,paste0(method,".Conditioned_Effect")]
    } else{
      GenomicRanges::mcols(gr.snp)[,"PP"] <- GenomicRanges::mcols(gr.snp)[,paste0(method,".PP")]
    }
  }
  dat <- as.data.frame(gr.snp)
  ### Label set
  labelSNPs <- construct_SNPs_labels(subset_DT = dat,
                                     lead="Lead SNP" %in% labels_subset,
                                     consensus="Consensus SNP" %in% labels_subset,
                                     method=T,
                                     remove_duplicates = F)
  labelSNPs_labels <- construct_SNPs_labels(subset_DT = dat,
                                            lead="Lead SNP" %in% labels_subset,
                                            consensus="Consensus SNP" %in% labels_subset,
                                            method=T,
                                            remove_duplicates = T)
  leader_SNP <- subset(labelSNPs, type=="Lead SNP")
  CS_set <- subset(labelSNPs, type=="Credible Set")

  ## Make track
  if(method=="original"){
    sig_cutoff <- -log10(sig_cutoff)
    cutoff_lab <- paste("P <",sig_cutoff)
    ymax <- max(-log10(gr.snp$P))
    r2_multiply <- 150
    a1 <- ggbio::plotGrandLinear(gr.snp,
                   geom = "point",
                   coord = "genome",
                   size=2,
                   alpha=.8,
                   aes(y = -log10(P), x=POS, color=r2),
                   facets=SEQnames~.) +
      labs(y="-log10(P-value)") +
      ylim(c(0,ymax*1.1))
  } else {
    sig_cutoff <- PP_threshold
    cutoff_lab <- paste0(PP_threshold*100,"% probability")
    r2_multiply <- 5
    # Further filter label tags if plotting fine-mapping results
    labelSNPs_labels <- labelSNPs_labels[ (labelSNPs_labels[paste0(method,".CS")]>0),] # IMPORTANT!: >0 (not TRUE)
    a1 <- ggbio::plotGrandLinear(gr.snp,
                   geom = "point",
                   coord = "genome",
                   aes(y = PP, x=POS, color=r2),
                   legend = F,
                   size=2,
                   alpha=.8,
                   facets=SEQnames~.) +
      ylim(c(0,1.1)) # PP is always 0-1 scale
    if(method=="COJO"){a1 <- a1 + labs(y="Conditioned Effect")}
  }

  a1 <- a1 +
    scale_color_gradient(low="blue", high="red", limits = c(0,1)) +
    geom_hline(yintercept=0,alpha=.5, linetype=1, size=.5) +
    geom_hline(yintercept=sig_cutoff, alpha=.5, linetype=2, size=.5, color="black") +
    geom_text(aes(x=min(POS), y=sig_cutoff*1.1),
              label=cutoff_lab,
              size=3,color="grey", hjust = 0)
  # Only draw LD line on the first track
  if(color_r2 & method=="original"){
    a1 <- a1 + geom_line(data=dat, stat="smooth",
                         aes(x=POS,y=r2*r2_multiply), #y=scales::rescale(gr.snp$r2, to=c(0,20))),
                         se = F, formula = y ~ x,  method = 'loess',
                         span=.1, size=.5, color="firebrick1")
  }
  a1 <- a1 +
    # Add diamond overtop leadSNP
    geom_point(data=labelSNPs,
               pch=labelSNPs$shape,
               fill=NA,
               size=labelSNPs$size,
               color=labelSNPs$color) +
    ### Background color label
    ggrepel::geom_label_repel(data=labelSNPs_labels,
                     aes(label=SNP),
                     color=NA,
                     # nudge_x = .5,
                     fill="black",
                     box.padding = .25,
                     label.padding = .25,
                     label.size=NA,
                     alpha=.6,
                     seed = 1,
                     size = 3,
                     min.segment.length = 1) +
    ### Foreground color label
    ggrepel::geom_label_repel(data=labelSNPs_labels,
                     aes(label=SNP),
                     color=labelSNPs_labels$color,
                     segment.alpha = .5,
                     # nudge_x = .5,
                     box.padding = .25,
                     label.padding = .25,
                     segment.size = 1,
                     fill = NA,
                     alpha=1,
                     seed = 1,
                     size = 3,
                     min.segment.length = 1) +
    theme_classic() +
    theme(legend.title = element_text(size=8),
          legend.text = element_text(size=6),
          strip.text.y = element_text(angle = 0),
          strip.text = element_text(size=9),
          legend.key.width=grid::unit(.5,"line"),
          legend.key.height=grid::unit(.5,"line"))
  # if(show.legend==F){a1 <- GGBIO.invisible_legend(gg = a1)}
  return(a1)
}




guess_pvalue_col <- function(dat,
                             QTL_prefix=NULL){
  p_options <- c("p","p-value","p-values","pvalue","pvalues","pval")
  p_options <- c(p_options, stringr::str_to_sentence(p_options))
  pval_col <- paste0(QTL_prefix,p_options)[paste0(QTL_prefix,p_options) %in% colnames(dat)][1]
  return(pval_col)
}




#' Plot QTL data
#'
#' @inheritParams finemap_pipeline
#' @inheritParams GGBIO.SNP_track
 #' @family plot
#' @keywords internal
#' @examples
#' data("BST1")
#' gr.snp <- DT_to_GRanges(dat)
#' ## Create fake QTL P-values
#' gr.snp$fake_eQTL.P <- gr.snp$P  * c(1,.9,.7)
#' qtl.track <- GGBIO.QTL_track(gr.snp=gr.snp, QTL_prefix="fake_eQTL.", labels_subset=c("Lead SNP"), color_r2=F)
GGBIO.QTL_track <- function(gr.snp,
                            QTL_prefix="QTL",
                            labels_subset = c("Lead SNP", "Credible Set", "Consensus SNP"),
                            color_r2=T,
                            show.legend=T,
                            PP_threshold=.95,
                            sig_cutoff=5e-8){
  if(color_r2==F){gr.snp$r2 <- NA}
  dat <- as.data.frame(gr.snp)
  pval_col <- guess_pvalue_col(dat = dat,
                               QTL_prefix = QTL_prefix)
  ### Label set
  labelSNPs <- construct_SNPs_labels(subset_DT = dat,
                                     lead="Lead SNP" %in% labels_subset,
                                     consensus="Consensus SNP" %in% labels_subset,
                                     method=T,
                                     remove_duplicates = F)
  labelSNPs_labels <- construct_SNPs_labels(subset_DT = dat,
                                            lead="Lead SNP" %in% labels_subset,
                                            consensus="Consensus SNP" %in% labels_subset,
                                            method=T,
                                            remove_duplicates = T)
  leader_SNP <- subset(labelSNPs, type=="Lead SNP")
  CS_set <- subset(labelSNPs, type=="Credible Set")

  ## Make track
  sig_cutoff <- -log10(sig_cutoff)
  cutoff_lab <- paste("P <",sig_cutoff)
  ymax <- max(-log10(gr.snp$P))
  r2_multiply <- 150
  q1 <- ggbio::plotGrandLinear(gr.snp,
                               geom = "point",
                               coord = "genome",
                               size=2,
                               alpha=.8,
                               aes(y = -log10(eval(parse(text=pval_col))),
                                   x=POS, color=r2),
                               facets=SEQnames~.) +
    labs(y="-log10(P-value)") +
    ylim(c(0,ymax*1.1)) +
    geom_hline(yintercept=0,alpha=.5, linetype=1, size=.5) +
    geom_hline(yintercept=sig_cutoff, alpha=.5, linetype=2, size=.5, color="black") +
    geom_text(aes(x=min(POS), y=sig_cutoff*1.1),
              label=cutoff_lab,
              size=3,color="grey", hjust = 0)
  # Only draw LD line on the first track
  if(color_r2 & ("r2" %in% colnames(dat)) ){
    q1 <- q1 + geom_line(data=dat, stat="smooth",
                         aes(x=POS,y=r2*r2_multiply), #y=scales::rescale(gr.snp$r2, to=c(0,20))),
                         se = F, formula = y ~ x,  method = 'loess',
                         span=.1, size=.5, color="firebrick1") +
      scale_color_gradient(low="blue", high="red", limits = c(0,1))
  }
  q1 <- q1 +
    # Add diamond overtop leadSNP
    geom_point(data=labelSNPs,
               pch=labelSNPs$shape,
               fill=NA,
               size=labelSNPs$size,
               color=labelSNPs$color) +
    ### Background color label
    ggrepel::geom_label_repel(data=labelSNPs_labels,
                              aes(label=SNP),
                              color=NA,
                              # nudge_x = .5,
                              fill="black",
                              box.padding = .25,
                              label.padding = .25,
                              label.size=NA,
                              alpha=.6,
                              seed = 1,
                              size = 3,
                              min.segment.length = 1) +
    ### Foreground color label
    ggrepel::geom_label_repel(data=labelSNPs_labels,
                              aes(label=SNP),
                              color=labelSNPs_labels$color,
                              segment.alpha = .5,
                              # nudge_x = .5,
                              box.padding = .25,
                              label.padding = .25,
                              segment.size = 1,
                              fill = NA,
                              alpha=1,
                              seed = 1,
                              size = 3,
                              min.segment.length = 1) +
    theme_classic() +
    theme(legend.title = element_text(size=8),
          legend.text = element_text(size=6),
          strip.text.y = element_text(angle = 0),
          strip.text = element_text(size=9),
          legend.key.width=grid::unit(.5,"line"),
          legend.key.height=grid::unit(.5,"line"))
   return(q1)
}







#' Track plot for transcript models
#'
#' Uses \code{\link{ggbio}}.
#'
#' @param gr.snp_CHR \emph{echolocator} fine-mapping results
#' in \code{\link{GRanges}} format.
#' @param max_transcripts Max number of transcripts per gene to show.
#' @inheritParams GGBIO.SNP_track
#' @family plot
#' @keywords internal
#' @examples
#' data("BST1");
#' gr.snp_CHR <- biovizBase::transformDfToGr(BST1, seqnames = "CHR", start = "POS", end = "POS")
#' track.genes <- GGBIO.transcript_model_track(gr.snp_CHR, max_transcripts=1)
GGBIO.transcript_model_track <- function(gr.snp_CHR,
                                         max_transcripts=1,
                                         show.legend=T,
                                         show_plot=F){
  # library(ggbio)
  printer("++ GGBIO:: Gene Model Track")
  lead.index <- which(gr.snp_CHR$leadSNP==T)
  printer("+ Annotating at transcript-level.")
  db.gr <- ensembldb::transcripts(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75) %>%
    data.table::as.data.table() %>%
    dplyr::mutate(index=row.names(.)) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::slice(1:max_transcripts)
  if(!"symbol" %in% colnames(db.gr)){
    db.gr$symbol <- ensembl_to_hgnc(db.gr$gene_id)
  }
  # Subset to only the region encompassed by the sumstats
  db.gr <- subset(db.gr,
                    seqnames == unique(gr.snp_CHR$CHR) &
                    start >= min(gr.snp_CHR$POS) &
                    end <= max(gr.snp_CHR$POS)) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  db.gr$symbol <- factor(db.gr$symbol, levels = unique(db.gr$symbol), ordered = T)
  edb <-  ensembldb::addFilter( EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                AnnotationFilter::TxIdFilter(db.gr$tx_id))
  track.genes <- suppressWarnings(suppressMessages(
   ggbio::autoplot(edb,
                            # Have to limit (can only handle depth < 1000)
                            which = db.gr,
                            names.expr = "gene_name",
                            aes(fill=gene_name, color=gene_name),
                            show.legend=show.legend)  +
      theme_classic() +
      theme(strip.text.y = element_text(angle = 0),
            strip.text = element_text(size=9),
            legend.text = element_text(size=5),
            legend.key.width=grid::unit(.1,"line"),
            legend.key.height=grid::unit(.1,"line")) +
      guides(fill=guide_legend(override.aes = list(size=1), ncol=4),
             color=guide_legend(override.aes = list(size=1), ncol=4),
             size=.5)
  )) # end suppress
  if(show_plot){ print(track.genes)}
  return(track.genes)
}




get_window_limits <- function(finemap_dat,
                              plot.window){
  lead.pos <- subset(finemap_dat, leadSNP)$POS
  if(!is.null(plot.window)){
    xlims <- c(lead.pos-as.integer(plot.window/2),
               lead.pos+as.integer(plot.window/2))
  } else {
    xlims <- c(min(finemap_dat$POS),
               max(finemap_dat$POS))
  }
  return(xlims)
}




# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
############ PLOT ALL TRACKS TOGETHER ############
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


#' Multi-track plot of fine-mapping results
#'
#' Plots fine-mapping results from \emph{echolocatoR}.
#' Create separate tracks for the original results (-log10(P))
#' and each fine-mapping method (Posterior Probability).
#' Uses \code{\link{ggbio}}.
#' @param QTL_prefixes Prefixes to the columns that contain QTL data.
#' For exmaple, P-values for the \emph{Fairfax_2014} QTL study would be
#' stored in the column \emph{Fairfax_2014.P}.
#' @param color_r2 If an LD_matrix is supplied, color each SNP by its LD with the lead SNP from the GWAS ()
#' @inheritParams GGBIO.SNP_track
#' @inheritParams GGBIO.QTL_track
#' @inheritParams GGBIO.transcript_model_track
#' @inheritParams finemap_pipeline
#' @family plot
#' @keywords internal
#' @source
#' \url{http://bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf}
#' @examples
#' \dontrun{
#' data("BST1"); data("BST1_LD_matrix"); data("locus_dir");
#' locus_dir <- file.path("~/Desktop",locus_dir)
#'
#' # Using NO annotations
#' trk_plot <- GGBIO.plot(finemap_dat=BST1, LD_matrix=BST1_LD_matrix, locus_dir=locus_dir, XGR_libnames=NULL, save_plot=F, color_r2=T)
#'
#' Using NO annotations (dot plot summary instead of Manhattan plots for each fine-mapping tool)
#' ### WARNING: Currently doesn't align as well due to the summary plot having a differnt x-axis.
#' trk_plot <- GGBIO.plot(finemap_dat=BST1, LD_matrix=BST1_LD_matrix, locus_dir=locus_dir, method_list=NULL, mean.PP=F, dot_summary=T, XGR_libnames=NULL, save_plot=F)
#'
#'
#' # Using only XGR annotations
#' trk_plot <- GGBIO.plot(finemap_dat=BST1, LD_matrix=BST1_LD_matrix, locus_dir=locus_dir, XGR_libnames=c("ENCODE_TFBS_ClusteredV3_CellTypes"), save_plot=F)
#'
#' # Using only Roadmap annotations
#' trk_plot <- GGBIO.plot(finemap_dat=BST1, LD_matrix=BST1_LD_matrix, locus_dir=locus_dir, XGR_libnames=NULL, Roadmap=T, Roadmap_query="monocyte", save_plot=F)
#'
#' # Using only Nott_2019 annotations
#' trk_plot <- GGBIO.plot(finemap_dat=BST1, LD_matrix=BST1_LD_matrix, locus_dir=locus_dir, Nott_epigenome=T, XGR_libnames=NULL, plot.window=100000)
#' }
GGBIO.plot <- function(finemap_dat,
                       locus_dir,
                       LD_matrix=NULL,
                       color_r2=T,
                       method_list=c("ABF","FINEMAP","SUSIE","POLYFUN_SUSIE"),
                       dot_summary=F,
                       QTL_prefixes=NULL,
                       mean.PP=T,
                       PP_threshold=.95,
                       consensus_thresh=2,
                       sig_cutoff=5e-8,

                       XGR_libnames=c("ENCODE_TFBS_ClusteredV3_CellTypes",
                                      "ENCODE_DNaseI_ClusteredV3_CellTypes",
                                      "Broad_Histone"),
                       n_top_xgr=5,

                       Roadmap=F,
                       Roadmap_query=NULL,
                       n_top_roadmap=7,
                       annot_overlap_threshold=5,

                       Nott_epigenome=F,
                       Nott_regulatory_rects=T,
                       Nott_show_placseq=T,
                       Nott_binwidth=200, #2500
                       Nott_bigwig_dir=NULL,

                       save_plot=T,
                       show_plot=T,
                       max_transcripts=1,
                       plot.window=NULL,
                       dpi=300,
                       height=12,
                       width=10,
                       verbose=T){
  # consensus_thresh=2; XGR_libnames="ENCODE_TFBS_ClusteredV3_CellTypes";n_top_xgr=5; mean.PP=T; Roadmap=T; PP_threshold=.95;  Nott_epigenome=T;  save_plot=T; show_plot=T; method_list=c("ABF","SUSIE","POLYFUN_SUSIE","FINEMAP","mean"); full_data=T;  max_transcripts=3; plot.window=100000;
  # Nott_epigenome=T; Nott_regulatory_rects=T; Nott_show_placseq=T; Nott_binwidth=2500; max_transcripts=1; dpi=400; height=12; width=10; results_path=NULL;  n_top_roadmap=7; annot_overlap_threshold=5; Nott_bigwig_dir=NULL; locus="BST1"; Roadmap_query=NULL; sig_cutoff=5e-8; verbose=T; QTL_prefixes=NULL;

  # tracks <- ggbio::tracks
  tracks <- get("tracks", asNamespace("ggbio"))
  locus <- basename(locus_dir)
  dir.create(locus_dir, showWarnings = F, recursive = T)
  # Set up data
  finemap_dat <- find_consensus_SNPs(finemap_dat = finemap_dat,
                                    credset_thresh = PP_threshold,
                                    consensus_thresh = consensus_thresh,
                                    verbose = F)
  finemap_dat <- fillNA_CS_PP(finemap_dat = finemap_dat)

  available_methods <- gsub("\\.PP$","",grep("*\\.PP$",colnames(finemap_dat),
                                             value = T)) %>% unique()
  method_list <- unique(method_list[method_list %in% available_methods])
  if(mean.PP){method_list <- unique(c(method_list, "mean"))}

  # Set window limits
  xlims <- get_window_limits(finemap_dat = finemap_dat,
                             plot.window = plot.window)

  TRACKS_list <- NULL

  # Add LD into the dat
  LD_SNP <- subset(finemap_dat, leadSNP==T)$SNP
  if(is.null(LD_matrix)){
    print("GGBIO:: No LD_matrix detected. Setting color_r2=F");
    color_r2 <- F
    dat <- finemap_dat
  } else {
    print("GGBIO:: LD_matrix detected. Coloring SNPs by LD with lead SNP.")
    LD_sub <- LD_with_leadSNP(LD_matrix = LD_matrix,
                              LD_SNP = LD_SNP)
    dat <- data.table::merge.data.table(finemap_dat, LD_sub,
                                        by = "SNP",
                                        all.x = T)
  }

  # Convert to GRange object
  gr.snp <- DT_to_GRanges(dat)
  gr.snp_CHR <- biovizBase::transformDfToGr(dat, seqnames = "CHR",
                                            start = "POS", end = "POS")


  # Treack 0: Summary
  if(dot_summary){
    printer("++ GGBIO:: Creating dot plot summary of fine-mapping results.")
    gp <- GGBIO.dot_summary(finemap_dat = finemap_dat,
                            PP_threshold = PP_threshold,
                            show_plot = F)
    TRACKS_list <- append(TRACKS_list, list(gp))
    names(TRACKS_list)[ifelse(is.null(TRACKS_list),1,length(TRACKS_list))] <- "Fine-mapping\nSummary"
  }


  # Track 1: GWAS
  printer("++ GGBIO::","GWAS","track", v=verbose)
  track.gwas <- GGBIO.SNP_track(gr.snp = gr.snp,
                          method = "original",
                          sig_cutoff=sig_cutoff,
                          labels_subset = c("Lead SNP","Consensus SNP"),
                          color_r2 = color_r2)
  TRACKS_list <- append(TRACKS_list, track.gwas)
  names(TRACKS_list)[ifelse(is.null(TRACKS_list),1,length(TRACKS_list))] <- "GWAS"




  for (qtl in QTL_prefixes){
    printer("++ GGBIO::",qtl,"track", v=verbose)
    qtl_track <- GGBIO.QTL_track(gr.snp = gr.snp,
                                 QTL_prefix=qtl,
                                 labels_subset = c("Lead SNP", "Consensus SNP"),
                                 color_r2=color_r2,
                                 show.legend=F,
                                 PP_threshold=PP_threshold,
                                 sig_cutoff=sig_cutoff)
    TRACKS_list <- append(TRACKS_list, qtl_track)
    names(TRACKS_list)[length(TRACKS_list)] <- qtl
  }


  # Tracks 2n: Fine-mapping
  for(m in method_list){
    printer("++ GGBIO::",m,"track", v=verbose)
    track.finemapping <- GGBIO.SNP_track(gr.snp, method = m,
                                   labels_subset = c("Lead SNP", "Credible Set"),
                                   color_r2 = color_r2,
                                   show.legend = F)
    TRACKS_list <- append(TRACKS_list, track.finemapping)
    names(TRACKS_list)[length(TRACKS_list)] <- m
  }

  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Track 3: Gene Model Track
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # DB tutorial: https://rdrr.io/bioc/ensembldb/f/vignettes/ensembldb.Rmd
  track.genes <- GGBIO.transcript_model_track(gr.snp_CHR = gr.snp_CHR,
                                              max_transcripts = max_transcripts)
  TRACKS_list <- append(TRACKS_list, track.genes)
  names(TRACKS_list)[length(TRACKS_list)] <- "Gene Track"

  # Track 3: Annotation - XGR Annotations
  ## Download
  palettes <- c("Spectral","BrBG","PiYG", "PuOr")
  if(any(!is.null(XGR_libnames))){printer("++ GGBIO:: XGR Tracks")}
  i=1
  for(lib in XGR_libnames){
    annot_file <- annotation_file_name(locus_dir = locus_dir,
                                       lib_name = paste0("XGR.",lib))
    gr.lib <- XGR.download_and_standardize(lib.selections = lib,
                                           finemap_dat = finemap_dat,
                                           nCores=1)
    gr.filt <- XGR.filter_sources(gr.lib=gr.lib, n_top_sources=n_top_xgr)
    gr.filt <- XGR.filter_assays(gr.lib=gr.filt, n_top_assays=n_top_xgr)
    saveRDS(gr.filt, annot_file)
    xgr.track <- XGR.plot_peaks(gr.lib=gr.filt,
                                subset_DT=finemap_dat,
                                fill_var = "Assay",
                                facet_var = "Source",
                                geom = "density",
                                show_plot = F)
    colourCount <- length(unique(gr.filt$assay))
    xgr.track <- xgr.track +
      theme_classic() +
      theme(strip.text.y = element_text(angle = 0),
            strip.text = element_text(size=9)) +
      scale_fill_manual(values = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, palettes[i]))(colourCount) ) +
      scale_y_continuous(n.breaks = 3) +
      guides(fill = guide_legend(ncol = 2, keyheight = .5, keywidth = .5))
    TRACKS_list <- append(TRACKS_list, xgr.track)
    names(TRACKS_list)[length(TRACKS_list)] <- gsub("_|[.]","\n",lib)
    i = i+1
  }

  # Track 4: Roadmap Chromatin Marks API
  ## Download
  if(Roadmap){
    printer("+ GGBIO:: Creating Roadmap track")
    lib <- "Roadmap.ChromatinMarks_CellTypes"
    annot_file <- annotation_file_name(locus_dir = locus_dir,
                                       lib_name = lib)
    if(file.exists(annot_file)){
      printer("+ Saved annotation file detected. Loading...")
      grl.roadmap <- readRDS(annot_file)
    } else {
      grl.roadmap <- ROADMAP.query(results_path = dirname(annot_file),
                                   gr.snp = gr.snp,
                                   keyword_query = Roadmap_query,
                                   limit_files=F)
      saveRDS(grl.roadmap, annot_file)
    }
    grl.roadmap.filt <- ROADMAP.merge_and_process_grl(grl.roadmap = grl.roadmap,
                                                      gr.snp = gr.snp,
                                                      n_top_tissues=n_top_roadmap,
                                                      sep = "\n")
    track.roadmap <- ROADMAP.track_plot(grl.roadmap.filt=grl.roadmap.filt,
                                        gr.snp=gr.snp,
                                        show_plot = F)
    TRACKS_list <- append(TRACKS_list, track.roadmap)
    names(TRACKS_list)[length(TRACKS_list)] <- "Roadmap\nChromatinMarks\nCellTypes"
  }


  if(Nott_epigenome){
    # Epigenomic histograms
    track.Nott_histo <- NOTT_2019.epigenomic_histograms(finemap_dat = finemap_dat,
                                                        locus_dir = locus_dir,
                                                        show_plot=F,
                                                        save_plot=F,
                                                        full_data=T,
                                                        plot.window = plot.window,
                                                        return_assay_track=T,
                                                        binwidth=Nott_binwidth,
                                                        bigwig_dir=Nott_bigwig_dir,
                                                        save_annot=T)
    TRACKS_list <- append(TRACKS_list, track.Nott_histo)
    names(TRACKS_list)[length(TRACKS_list)] <- "Nott (2019)\nRead Densities"

    # PLAC-seq
    if(Nott_show_placseq){
      track.Nott_plac <- NOTT_2019.plac_seq_plot(finemap_dat = finemap_dat,
                                                 locus_dir=locus_dir,
                                                 title=locus,
                                                 show_regulatory_rects=Nott_regulatory_rects,
                                                 return_interaction_track=T,
                                                 show_arches=T,
                                                 save_annot=T)
      TRACKS_list <- append(TRACKS_list, track.Nott_plac)
      names(TRACKS_list)[length(TRACKS_list)] <- "Nott (2019)\nPLAC-seq"
    }

  }

  # ------- Fuse all tracks
  if(dot_summary){
    # Dot summary requires a special way of merging bc it doesn't share the same x-axis
    grob_list <- lapply(TRACKS_list, function(x){ggbio:::Grob(x)})
    TRKS_FINAL <- patchwork::wrap_plots(grob_list, ncol = 1)
  } else {
    heights <- GGBIO.track_heights_dict(TRACKS_list = TRACKS_list)
    params_list <- list(title = paste0(locus," locus [",length(GenomicRanges::seqnames(gr.snp))," SNPs]"),
                        track.bg.color = "transparent",
                        track.plot.color = "transparent",
                        label.text.cex = .7,
                        label.bg.fill = "gainsboro",
                        label.text.color = "grey12",
                        label.text.angle = 0,
                        label.width = grid::unit(5.5, "lines"),
                        xlim = xlims,
                        heights = heights)
    trks <- suppressMessages(do.call("tracks", append(TRACKS_list, params_list)))
    TRKS_FINAL <- GGBIO.add_lines(trks = trks,
                                  finemap_dat = finemap_dat)
  }

  if(save_plot){
    window_suffix <- ifelse(is.null(plot.window),"",paste0(plot.window/1000,"kb"))
    plot_path <- file.path(locus_dir,paste0("multiview_",locus,"_",window_suffix,".png"))
    printer("+ GGBIO:: Saving plot ==>",plot_path)
    ggbio::ggsave(filename = plot_path,
                  plot = TRKS_FINAL,
                  height = height,
                  width = width,
                  dpi = dpi,
                  bg = "transparent")
  }
  if(show_plot){print(TRKS_FINAL)}
  return(TRKS_FINAL)
}




#' Determine track heights
#' @family plot
#' @keywords internal
GGBIO.track_heights_dict <- function(TRACKS_list,
                                     Manhattan=.2,
                                     Gene=.2,
                                     XGR=.5,
                                     Roadmap=1,
                                     Nott_epigenome=1,
                                     Nott_placseq=.33){
  dict <- c("GWAS"=Manhattan,
            "QTL"=Manhattan,
            "ABF"=Manhattan,
            "SUSIE"=Manhattan,
            "POLYFUN_SUSIE"=Manhattan,
            "FINEMAP"=Manhattan,
            "mean"=Manhattan,
            "Gene Track"=Gene,
            "Roadmap\nChromatinMarks\nCellTypes"=Roadmap,
            "Nott (2019)\nRead Densities"=Nott_epigenome,
            "Nott (2019)\nPLAC-seq"=Nott_placseq
  )
  heights <- lapply(names(TRACKS_list),function(trk_name){
    if(endsWith(trk_name,"QTL")){trk_name <- "QTL" }
    if(trk_name %in% names(dict)){dict[[trk_name]]
    } else{.33}
  }) %>% unlist()
  return(heights)
}



#' Add vertical lines
#'
#' Adds vertical lines indicates key SNPs
#' @family plot
#' @keywords internal
GGBIO.add_lines <- function(trks,
                            finemap_dat,
                            alpha=.7,
                            size=.3){
  # Add lines
  lead.pos <- subset(finemap_dat, leadSNP)$POS
  consensus.pos <- subset(finemap_dat, Consensus_SNP==T)$POS

  TRKS_FINAL <- suppressWarnings(suppressMessages(
    trks +
     geom_vline(xintercept = consensus.pos, color="goldenrod2",
                 alpha=alpha, size=size, linetype='solid') +
     geom_vline(xintercept = lead.pos, color="red",
                alpha=alpha, size=size, linetype='solid') +
     theme(strip.text.y = element_text(angle = 0),
           strip.text = element_text(size = 7),
           panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
           plot.subtitle = element_text(color = "turquoise", size = 8)) +
     scale_x_continuous( labels=function(x)x/1000000)))
  return(TRKS_FINAL)
}

