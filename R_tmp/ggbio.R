# ----------------------- #
# ----- ggbio plots ------#
# ----------------------- #

invisible_legend <- function(gg){
  gg <- gg + theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  ) +
    scale_color_discrete(
      guide = guide_legend(override.aes = list(color = "white"))
    ) +
    scale_fill_discrete(
      guide = guide_legend(override.aes = list(color = "white", fill="white"))
    )
  return(gg)
}

### ALL LOCI

ggbio.all_loci <- function(FM = merge_finemapping_results()){
  # Convert to GRange object
  DT <- FM %>% dplyr::mutate(SEQnames = paste0("chr",CHR))
  gr.snp <- biovizBase::transformDfToGr(DT, seqnames = "SEQnames", start = "POS", end = "POS")
  gr.snp_CHR <- biovizBase::transformDfToGr(FM, seqnames = "CHR", start = "POS", end = "POS")


  snp.track <- SNP_track(gr.snp, r2 = NULL, labels_subset = c("Lead SNP", "Consensus SNP"))
}


### SINGLE LOCUS

GR.name_filter_convert <- function(GR.final, GR.names, min_hits=1){
  names(GR.final) <- GR.names
  grl <- GR.final[!as.logical(lapply(GR.final, is.null))]
  # Filter to those that had at least N hits
  grl <- grl[as.logical(lapply(grl, function(g, min_hits.=min_hits){length(seqnames(g)) >= min_hits.}))]
  # Convert to GRangesList (important)
  grl <- GRangesList(grl)
  return(grl)
}


SNP_track <- function(gr.snp,
                      method = "original",
                      labels_subset = c("Lead SNP", "Credible Set", "Consensus SNP"),
                      r2=NULL,
                      show.legend=T,
                      PP_threshold=.95){
  # Format data
  if(!(method %in% c("original"))){
    if(method=="COJO"){
      mcols(gr.snp)[,"PP"] <- mcols(gr.snp)[,paste0(method,".Conditioned_Effect")]
    } else{
      mcols(gr.snp)[,"PP"] <- mcols(gr.snp)[,paste0(method,".PP")]
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
    cutoff <- -log10(5e-8)
    cutoff_lab <- paste("P < 5e-8")
    ymax <- max(-log10(gr.snp$P))
    r2_multiply <- 150
    a1 <- plotGrandLinear(gr.snp,
                   geom = "point",
                   coord = "genome",
                   size=2,
                   alpha=.8,
                   aes(y = -log10(P), x=POS, color=r2),
                   facets=SEQnames~.) +
      labs(y="-log10(P-value)") +
      ylim(c(0,ymax*1.1))
  } else {
    cutoff <- PP_threshold
    cutoff_lab <- paste0(PP_threshold*100,"% probability")
    r2_multiply <- 5
    # Further filter label tags if plotting fine-mapping results
    labelSNPs_labels <- labelSNPs_labels[ (labelSNPs_labels[paste0(method,".Credible_Set")]>0),] # IMPORTANT!: >0 (not TRUE)
    a1 <- plotGrandLinear(gr.snp,
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
    geom_hline(yintercept=cutoff, alpha=.5, linetype=2, size=.5, color="black") +
    geom_text(aes(x=min(POS), y=cutoff*1.1),
              label=cutoff_lab,
              size=3,color="grey", hjust = 0)
  # Only draw LD line on the first track
  if(is.null(r2) & method=="original"){
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
          legend.key.width=unit(.5,"line"),
          legend.key.height=unit(.5,"line"))
  # if(show.legend==F){a1 <- a1 + scale_fill_discrete(guide=F) + guides(fill=F) + theme(legend.position="none")}
  return(a1)
}




####### XGR track


save_annotations <- function(gr, anno_path, libName){
  dir.create(dirname(anno_path), showWarnings = F, recursive = T)
  saveRDS(gr, file.path(anno_path))
}



transcript_model_track <- function(gr.snp_CHR,
                                   max_transcripts=1,
                                   show.legend=T){
  library(ggbio)
  printer("++ GGBIO:: Gene Model Track")
  lead.index <- which(gr.snp_CHR$leadSNP==T)
  print("+ Annotating at transcript-level.")
  db.gr <- ensembldb::transcripts(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75) %>%
    data.table::as.data.table() %>%
    dplyr::mutate(index=row.names(.)) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::slice(1:max_transcripts)
  if(!"symbol" %in% colnames(db.gr)){
    db.gr$symbol <- ensembl_to_hgnc(db.gr$gene_id)
  }
  # Subset to only the region encompassed by the sumstats
  db.gr <- subset(db.gr, seqnames == unique(gr.snp_CHR$CHR) &
                    start >= min(gr.snp_CHR$POS) &
                    end <= max(gr.snp_CHR$POS)) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  db.gr$symbol <- factor(db.gr$symbol, levels = unique(db.gr$symbol), ordered = T) 
  edb <-  ensembldb::addFilter( EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75, 
                                AnnotationFilter::TxIdFilter(db.gr$tx_id))
  track.genes <- autoplot(edb,
                          # Have to limit (can only handle depth < 1000)
                          which = db.gr,
                          names.expr = "gene_name",
                          aes(fill=gene_name, color=gene_name),
                          show.legend=show.legend)  +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 0),
          strip.text = element_text(size=9),
          legend.text = element_text(size=5),
          legend.key.width=unit(.1,"line"),
          legend.key.height=unit(.1,"line")) +
    guides(fill=guide_legend(override.aes = list(size=1), ncol=4),
           color=guide_legend(override.aes = list(size=1), ncol=4),
           size=.5)
  track.genes
  return(track.genes)
}



# gene_model_track <- function(gr.snp_CHR,
#                              max_transcripts=1,
#                              show.legend=T){
#   library(ggbio)
#   printer("++ GGBIO:: Gene Model Track")
#   lead.index <- which(gr.snp_CHR$leadSNP==T)
#   print("+ Annotating at transcript-level.")
#   db.gr <- ensembldb::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75) %>%
#     data.table::as.data.table() %>%
#     dplyr::mutate(index=row.names(.)) %>%
#     dplyr::group_by(gene_name) %>%
#     dplyr::slice(1:max_transcripts)
#   if(!"symbol" %in% colnames(db.gr)){
#     db.gr$symbol <- ensembl_to_hgnc(db.gr$gene_id)
#   }
#   # Subset to only the region encompassed by the sumstats
#   db.gr <- subset(db.gr, seqnames == unique(gr.snp_CHR$CHR) &
#                     start >= min(gr.snp_CHR$POS) &
#                     end <= max(gr.snp_CHR$POS)) %>%
#     GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
#   db.gr$symbol <- factor(db.gr$symbol, levels = unique(db.gr$symbol), ordered = T) 
#   edb <-  ensembldb::addFilter( EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75, 
#                                 AnnotationFilter::GeneNameFilter(db.gr$gene_name))
#   track.genes <- autoplot(edb,
#                           # Have to limit (can only handle depth < 1000)
#                           which = db.gr,
#                           names.expr = "gene_name",
#                           aes(fill=gene_name, color=gene_name),
#                           show.legend=show.legend)  +
#     theme_classic() +
#     theme(strip.text.y = element_text(angle = 0),
#           strip.text = element_text(size=9),
#           legend.text = element_text(size=5),
#           legend.key.width=unit(.1,"line"),
#           legend.key.height=unit(.1,"line")) +
#     guides(fill=guide_legend(override.aes = list(size=1), ncol=4),
#            color=guide_legend(override.aes = list(size=1), ncol=4),
#            size=.5)
#   track.genes
#   return(track.genes)
# }

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
############ PLOT ALL TRACKS TOGETHER ############
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


GGBIO.plot <- function(finemap_DT,
                       LD_matrix,
                       gene=NULL,
                       results_path,
                       method_list=c("SUSIE","FINEMAP","PAINTOR",
                                     "PAINTOR_Fairfax"),
                       mean.PP=T,
                       XGR_libnames=c("ENCODE_TFBS_ClusteredV3_CellTypes",
                                      "ENCODE_DNaseI_ClusteredV3_CellTypes",
                                      "Broad_Histone"),
                       n_top_xgr=5,
                       ROADMAP=F,
                       n_top_roadmap=7,
                       annot_overlap_threshold=5,
                       PP_threshold=.95,
                       consensus_thresh=2,
                       
                       Nott_sn_epigenome=F,
                       show_regulatory_rects=T,
                       show_placseq=T, 
                       plot_Nott_binwidth=2500,
                       Nott_bigwig_dir=NULL,
                       
                       save_plot=T,
                       show_plot=T,
                       max_transcripts=1,
                       plot_window=NULL,
                       dpi=200){
  # http://bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
  library(ggbio)
  require(GenomicRanges)
  require(biovizBase)
  # consensus_thresh=2; XGR_libnames=NULL; mean.PP=T; ROADMAP=F; PP_threshold=.95;  Nott_sn_epigenome=T;  save_plot=T; show_plot=T; method_list=c("ABF","SUSIE","POLYFUN_SUSIE","FINEMAP","mean"); full_data=T;  max_transcripts=3; plot_window=100000;
  # quick_finemap("HIP1R", consensus_thresh);
  # Nott_sn_epigenome=T; show_regulatory_rects=T; show_placseq=T; plot_Nott_binwidth=2500; max_transcripts=1; dpi=400
  
   
  
  # Set up data
  finemap_DT <- find_consensus_SNPs(finemap_DT = finemap_DT,
                                    credset_thresh = PP_threshold,
                                    consensus_thresh = consensus_thresh)

  available_methods <- gsub("\\.PP$","",grep("*\\.PP$",colnames(finemap_DT),
                                             value = T)) %>% unique()
  method_list <- unique(method_list[method_list %in% available_methods])
  if(mean.PP){method_list <- unique(c(method_list, "mean"))}
  
  # Set window limits
  lead.pos <- subset(finemap_DT, leadSNP)$POS
  if(!is.null(plot_window)){
    xlims <- c(lead.pos-as.integer(plot_window/2),
               lead.pos+as.integer(plot_window/2))
  } else {
    xlims <- c(min(finemap_DT$POS),
               max(finemap_DT$POS))
  }

  TRACKS_list <- NULL
  # Add LD into the DT
  LD_SNP <- subset(finemap_DT, leadSNP==T)$SNP
  LD_sub <- LD_with_leadSNP(LD_matrix, LD_SNP)
  DT <- data.table:::merge.data.table(finemap_DT, LD_sub,
                                      by = "SNP", all.x = T)
  # Convert to GRange object
  DT <- DT %>% dplyr::mutate(SEQnames = paste0("chr",CHR))
  gr.snp <- biovizBase::transformDfToGr(DT, seqnames = "SEQnames", start = "POS", end = "POS")
  gr.snp_CHR <- biovizBase::transformDfToGr(DT, seqnames = "CHR", start = "POS", end = "POS")

  # Track 1: GWAS
  track.gwas <- SNP_track(gr.snp = gr.snp,
                          method = "original",
                          labels_subset = c("Lead SNP", "Consensus SNP"))
  TRACKS_list <- append(TRACKS_list, track.gwas)
  names(TRACKS_list)[1] <- "GWAS"

  # Tracks 2n: Fine-mapping
  for(m in method_list){
    printer("++ GGBIO::",m,"Track")
    track.finemapping <- SNP_track(gr.snp, method = m,
                             labels_subset = c("Lead SNP", "Credible Set"),
                             show.legend = F)
    TRACKS_list <- append(TRACKS_list, track.finemapping)
    names(TRACKS_list)[length(TRACKS_list)] <- m
  }

  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Track 3: Gene Model Track
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # DB tutorial: https://rdrr.io/bioc/ensembldb/f/vignettes/ensembldb.Rmd
  
  # track.genes <- invisible_legend(track.genes)
  track.genes <- transcript_model_track(gr.snp_CHR, 
                                        max_transcripts = max_transcripts)
  TRACKS_list <- append(TRACKS_list, track.genes)
  names(TRACKS_list)[length(TRACKS_list)] <- "Gene Track"

  # Track 3: Annotation - XGR Annotations
  ## Download
  library(RColorBrewer)
  palettes <- c("Spectral","BrBG","PiYG", "PuOr")
  counter <- 1
  if(any(!is.null(XGR_libnames))){printer("++ GGBIO:: Gene Model Track")}
  for(lib in XGR_libnames){
    anno_data_path <- file.path("echolocatoR/tools/Annotations", paste0("XGR_",lib,".rds"))
    grl.xgr <- XGR.import_annotations(gr.snp = gr.snp,
                                      anno_data_path = anno_data_path,
                                      lib.name = lib,
                                      save_xgr=T,
                                      annot_overlap_threshold=1)
    grl.xgr.merged.filt <- XGR.merge_and_process(grl.xgr = grl.xgr,
                                                 lib = lib,
                                                 n_top_sources=n_top_xgr)
    colourCount <- length(unique(grl.xgr.merged.filt$Assay))
    # Plot
    facet <- lib!="ENCODE_DNaseI_ClusteredV3_CellTypes"
    adjust <- ifelse(lib=="ENCODE_TFBS_ClusteredV3_CellTypes", .5, .2)
    # facet <-T
    if(facet==F){
      track.xgr <- ggbio::autoplot(grl.xgr.merged.filt, 
                                   which = gr.snp,
                                   aes(fill=Assay),
                                   # fill = "magenta",
                                   color = "white",#NA
                                   geom = "density",
                                   adjust = adjust,
                                   position="stack",
                                   # bins=50,
                                   size=.1,
                                   alpha = 1) 
    } else {
      track.xgr <- ggbio::autoplot(grl.xgr.merged.filt, which = gr.snp,
                                   aes(fill=Assay),
                                   # fill = "magenta",
                                   color = "white",#NA
                                   geom = "density",
                                   adjust = adjust,
                                   position="stack",
                                   # bins=50,
                                   size=.1,
                                   alpha = 1,
                                   facets=Source~.)
    }
    track.xgr <- track.xgr+
      theme_classic() +
      theme(strip.text.y = element_text(angle = 0),
            strip.text = element_text(size=9)) +
      scale_fill_manual(values = colorRampPalette(brewer.pal(8, palettes[counter]))(colourCount) ) + #scale_fill_brewer(palette=palettes[1])
      scale_y_continuous(n.breaks = 3) +
      guides(fill = guide_legend(ncol = 2, keyheight = .5, keywidth = .5))
    
    TRACKS_list <- append(TRACKS_list, track.xgr)
    new_name <- paste(strsplit(lib,"_")[[1]], collapse="\n")
    names(TRACKS_list)[length(TRACKS_list)] <- new_name
    counter = counter+1
  }

  # Track 4: Roadmap Chromatin Marks API
  ## Download
  if(ROADMAP){
    printer("+ GGBIO:: Creating ROADMAP track")
    lib <- "Roadmap_ChromatinMarks_CellTypes"
    anno_path <- file.path(results_path, "Annotation",paste0("GRanges_",lib,".rds"))
    if(file.exists(anno_path)){
      printer("+ Saved annotation file detected. Loading...")
      grl.roadmap <- readRDS(anno_path)
    } else {
      grl.roadmap <- ROADMAP.track(results_path = results_path,
                                   gr.snp = gr.snp,
                                   limit_files=NA)
      save_annotations(gr = grl.roadmap, anno_path = anno_path, libName = lib)
    }
    grl.roadmap.filt <- ROADMAP.merge_and_process_grl(grl.roadmap,  n_top_tissues=n_top_roadmap)
    ## Make track
    track.roadmap <- autoplot(grl.roadmap.filt, which = gr.snp,
                              aes(fill=ChromState),
                              color="white",
                              size=.1,
                              geom = "density",# density
                              adjust = 1,
                              # bins=10,
                              position="stack",# stack, fill, dodge
                              facets=Source~.,
                              alpha=1) +
      theme_classic() +
      theme(strip.text.y = element_text(angle = 0),
            strip.text = element_text(size=9 )) +
      guides(fill = guide_legend(ncol = 2, keyheight = .5, keywidth = .5)) +
      scale_y_continuous(n.breaks = 3)
    
    TRACKS_list <- append(TRACKS_list, track.roadmap)
    names(TRACKS_list)[length(TRACKS_list)] <- "ROADMAP\nChromatinMarks\nCellTypes"
  }


  if(Nott_sn_epigenome){
    # Epigenomic histograms
    track.Nott_histo <- NOTT_2019.epigenomic_histograms(finemap_DT,
                                                        results_path,
                                                        show_plot=F,
                                                        save_plot=F,
                                                        full_data=T,
                                                        return_assay_track=T,
                                                        binwidth = plot_Nott_binwidth, 
                                                        bigwig_dir=Nott_bigwig_dir)
    TRACKS_list <- append(TRACKS_list, track.Nott_histo)
    names(TRACKS_list)[length(TRACKS_list)] <- "Nott (2019)\nRead Densities"
     
    # PLAC-seq 
    if(show_placseq){
      track.Nott_plac <- NOTT_2019.plac_seq_plot(finemap_DT,
                                                 title=gene,
                                                 return_interaction_track = T,
                                                 show_arches = T)
      TRACKS_list <- append(TRACKS_list, track.Nott_plac)
      names(TRACKS_list)[length(TRACKS_list)] <- "Nott (2019)\nPLAC-seq"
    }
    
  }


  # Fuse all tracks 
  heights <- c(.2, # GWAS track
               rep(.15,length(method_list)), # Fine-mapping tracks 
               .2, # transcript track
               if(!is.null(XGR_libnames))rep(.5,length(XGR_libnames)),
               ifelse(ROADMAP,1,0),
               if(Nott_sn_epigenome){1}else{NULL},
               if(Nott_sn_epigenome & show_placseq){.33}else{NULL})
  params_list <- list(title = paste0(gene," locus [",length(seqnames(gr.snp))," SNPs]"),
                      track.bg.color = "transparent",
                      track.plot.color = "transparent",
                      label.text.cex = .7,
                      label.bg.fill = "gainsboro",
                      label.text.color = "grey12",
                      label.text.angle = 0,
                      label.width = unit(5.5, "lines"),
                      xlim = xlims,
                      heights = heights)
  trks <- suppressWarnings(do.call("tracks", append(TRACKS_list, params_list)))

  # Add lines
  lead.pos <- subset(finemap_DT, leadSNP)$POS
  consensus.pos <- subset(finemap_DT, Consensus_SNP==T)$POS
  trks_plus_lines <- trks +
    geom_vline(xintercept = consensus.pos, color="goldenrod2",
               alpha=1, size=.3, linetype='solid') +
    geom_vline(xintercept = lead.pos, color="red",
               alpha=1, size=.3, linetype='solid') +
    theme(strip.text.y = element_text(angle = 0),
          strip.text = element_text(size = 7),
          panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
          plot.subtitle = element_text(color = "turquoise", size = 8)) +
    scale_x_continuous( labels=function(x)x/1000000)  
    
  # Save
  if(save_plot){
    window_suffix <- ifelse(is.null(plot_window),"",paste0(plot_window/1000,"kb"))
    plot.path <- file.path(results_path,"Multi-finemap",paste0(gene,"_ggbio",window_suffix,".png"))
    printer("+ GGBIO:: Saving plot ==>",plot.path)

    ggsave(filename = plot.path,
           plot = trks_plus_lines,
           height = 12,#+length(XGR_libnames)+n_roadmap+n_Nott,
           width = 10, dpi = dpi, bg = "transparent")
  }
  
  if(show_plot){print(trks_plus_lines)}
  
  return(trks_plus_lines)
}







