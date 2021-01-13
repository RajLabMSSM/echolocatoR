# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# -------- PLOT ALL TRACKS TOGETHER -------------#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


#' Generate a locus plot
#'
#' Generate a locus-specific plot with multiple selectable tracks.
#' Users can also generate multiple zoomed in views of the plot at multiple resolutions.
#'
#' #### TO-do
#' Add anchors to plac-seq arches
#'
#'
#' @examples
#' library(echolocatoR)
#' finemap_dat<-echolocatoR::BST1; LD_matrix <- echolocatoR::BST1_LD_matrix;
#' locus_dir <- file.path("~/Desktop","results/GWAS/Nalls23andMe_2019/BST1")
#'
#' locus_plot <- PLOT.locus(finemap_dat, locus_dir=locus_dir, LD_matrix=LD_matrix, Nott_epigenome=T, xtext=F, plot.zoom=c("5x"))
PLOT.locus <- function(finemap_dat,
                       locus_dir,
                       LD_matrix=NULL,
                       LD_reference=NULL,
                       dataset_type="GWAS",
                       color_r2=T,
                       method_list=c("ABF","FINEMAP","SUSIE","POLYFUN_SUSIE"),
                       track_order= NULL,
                       track_heights=NULL,
                       plot_full_window=T,
                       dot_summary=F,
                       QTL_prefixes=NULL,
                       mean.PP=T,
                       PP_threshold=.95,
                       consensus_threshold=2,
                       sig_cutoff=5e-8,
                       gene_track=T,
                       point_size=1,
                       point_alpha=.6,
                       snp_group_lines=c("Lead","UCS","Consensus"),
                       xtext=F,
                       show.legend_genes=T,

                       XGR_libnames=NULL,
                       #c("ENCODE_TFBS_ClusteredV3_CellTypes",
                       #   "ENCODE_DNaseI_ClusteredV3_CellTypes",
                       # "Broad_Histone"),
                       n_top_xgr=5,

                       Roadmap=F,
                       Roadmap_query=NULL,
                       n_top_roadmap=7,
                       annot_overlap_threshold=5,

                       Nott_epigenome=F,
                       Nott_regulatory_rects=T,
                       Nott_show_placseq=T,
                       Nott_binwidth=200,
                       Nott_bigwig_dir=NULL,

                       save_plot=T,
                       show_plot=T,
                       genomic_units="Mb",
                       strip.text.y.angle=0,
                       max_transcripts=1,
                       plot.zoom=c("1x"),
                       dpi=300,
                       height=12,
                       width=10,
                       plot_format="jpg",
                       nThread=4,
                       return_list=F,
                       verbose=T){
  # library(dplyr); library(ggplot2); LD_reference="UKB";
  # finemap_dat<-echolocatoR::BST1; LD_matrix <- echolocatoR::BST1_LD_matrix; locus="BST1";
  # consensus_threshold=2; XGR_libnames=NULL; n_top_xgr=5; mean.PP=T; Roadmap=F; PP_threshold=.95;  Nott_epigenome=T;  save_plot=T; show_plot=T; method_list=c("ABF","SUSIE","POLYFUN_SUSIE","FINEMAP","mean"); full_data=T;  max_transcripts=3; pz=plot.zoom="1x"; dataset_type="GWAS"; dot_summary=F; snp_group_lines=c("UCS","Consensus","Lead"); nThread=4;
  # Nott_epigenome=T; Nott_regulatory_rects=T; Nott_show_placseq=T; Nott_binwidth=200; max_transcripts=1; dpi=400; height=12; width=10; results_path=NULL;  n_top_roadmap=7; annot_overlap_threshold=5; Nott_bigwig_dir=NULL;
  #  Roadmap_query=NULL; sig_cutoff=5e-8; verbose=T; QTL_prefixes=NULL; gene_track=T; genomic_units="Mb";strip.text.y.angle=0; xtext=F; plot_format="jpg"; return_list=F;
  # track_order= c("Genes","GWAS full window","zoom_polygon","GWAS","Fine-mapping", "Roadmap\nchromatin marks\ncell-types", "Nott (2019)\nread densities", "Nott (2019)\nPLAC-seq"); track_heights=NULL; plot_full_window=T;
  message("+-------- Locus Plots --------+")
  locus <- basename(locus_dir)
  dir.create(locus_dir, showWarnings = F, recursive = T)
  # Set up data
  finemap_dat <- find_consensus_SNPs(finemap_dat = finemap_dat,
                                     credset_thresh = PP_threshold,
                                     consensus_thresh = consensus_threshold,
                                     verbose = F)
  finemap_dat <- fillNA_CS_PP(finemap_dat = finemap_dat)
  finemap_dat$Mb <- finemap_dat$POS/1000000

  available_methods <- gsub("\\.PP$","",grep("*\\.PP$",colnames(finemap_dat),
                                             value = T)) %>% unique()
  method_list <- unique(method_list[method_list %in% available_methods])
  if(mean.PP){method_list <- unique(c(method_list, "mean"))}
  # Add LD into the dat
  finemap_dat <- LD.get_lead_r2(finemap_dat = finemap_dat,
                                LD_matrix = LD_matrix,
                                LD_format = "matrix")
  # Begin constructing tracks
  TRKS <- NULL;
  # Treack: Summary
  if(dot_summary){
    printer("++ PLOT:: Creating dot plot summary of fine-mapping results.")
    TRKS[["Summary"]] <- PLOT.dot_summary(finemap_dat = finemap_dat,
                                           PP_threshold = PP_threshold,
                                           show_plot = F)
  }
  ####  Track: Main (GWAS) frozen ####
  if(plot_full_window){
    printer("++ PLOT::",dataset_type,"full window track", v=verbose)
    full_window_name <- paste(dataset_type,"full window")
    TRKS[[full_window_name]] <- PLOT.SNP_track_merged(finemap_dat = finemap_dat,
                                                        yvar = "-log10(P)",
                                                        sig_cutoff = sig_cutoff,
                                                        labels_subset = NULL,
                                                        xtext = T,
                                                        show.legend = F,
                                                        dataset_type = gsub(" ","\n",full_window_name),
                                                        strip.text.y.angle = strip.text.y.angle,
                                                        verbose = verbose) +
      geom_point(data = subset(finemap_dat, leadSNP),
                 aes(x=POS, y=-log10(P)),
                 color="red",pch=9, size=3, show.legend = F, alpha=1) +
      theme(axis.title.x = element_blank())
  }

  ####  Track: Main (GWAS) ####
  printer("++ PLOT::",dataset_type,"track", v=verbose)
  TRKS[[dataset_type]] <- PLOT.SNP_track_merged(finemap_dat = finemap_dat,
                                          yvar = "-log10(P)",
                                          sig_cutoff = sig_cutoff,
                                          labels_subset = NULL,
                                          xtext = xtext,
                                          dataset_type = "GWAS",
                                          strip.text.y.angle = strip.text.y.angle,
                                          verbose = verbose) +
    geom_point(data = subset(finemap_dat, leadSNP),
               aes(x=POS, y=-log10(P)),
               color="red",pch=9, size=3, show.legend = F, alpha=1)

  #### Track: QTL ####
  for (qtl in QTL_prefixes){
    printer("++ PLOT::",qtl,"track", v=verbose)
    TRKS[[qtl]]  <- PLOT.SNP_track_merged(finemap_dat = finemap_dat,
                                          yvar = paste0("-log10(",qtl,"P)"),
                                          sig_cutoff = sig_cutoff,
                                          labels_subset = NULL,
                                          xtext = xtext,
                                          dataset_type = qtl,
                                          strip.text.y.angle = strip.text.y.angle,
                                          verbose = verbose)
  }

  #### Tracks: Fine-mapping ####
  printer("++ PLOT:: Merged fine-mapping track", v=verbose)
  TRKS[["Fine-mapping"]] <- PLOT.SNP_track_merged(finemap_dat = finemap_dat,
                                                  yvar = "PP",
                                                  sig_cutoff = sig_cutoff,
                                                  absolute_labels = F,
                                                  label_type = "rsid_only",
                                                  labels_subset = c("Lead","CS"),
                                                  show.legend = F,
                                                  xtext = xtext,
                                                  strip.text.y.angle = strip.text.y.angle,
                                                  verbose = verbose)

  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##### Track: Gene Model Track ####
  # DB tutorial: https://rdrr.io/bioc/ensembldb/f/vignettes/ensembldb.Rmd
  if(gene_track){
    printer("++ PLOT:: Adding Gene model track.",v=verbose)
    try({
      TRKS[["Genes"]] <- PLOT.transcript_model_track(finemap_dat = finemap_dat,
                                                     show.legend = show.legend_genes,
                                                     xtext = xtext,
                                                     method = "ggplot",
                                                     max_transcripts = 1,
                                                     collapseTranscripts = F,
                                                     verbose=T)
    })
  }


  #### Track: Annotation - XGR ####
  palettes <- c("Spectral","BrBG","PiYG", "PuOr")
  if(any(!is.null(XGR_libnames))){printer("++ PLOT:: XGR Tracks")}
  i=1
  for(lib in XGR_libnames){
    annot_file <- annotation_file_name(locus_dir = locus_dir,
                                       lib_name = paste0("XGR.",lib))
    gr.lib <- XGR.download_and_standardize(lib.selections = lib,
                                           finemap_dat = finemap_dat,
                                           nCores = nThread)
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

    TRKS[[gsub("_|[.]","\n",lib)]] <- xgr.track
    i = i+1
  }

  #### Track: Roadmap Chromatin Marks ####
  ## Download
  if(Roadmap){
    printer("+ PLOT:: Creating Roadmap track")
    lib <- "Roadmap.ChromatinMarks_CellTypes"
    annot_file <- annotation_file_name(locus_dir = locus_dir,
                                       lib_name = lib)
    if(file.exists(annot_file)){
      printer("+ Saved annotation file detected. Loading...")
      grl.roadmap <- readRDS(annot_file)
    } else {
      grl.roadmap <- ROADMAP.query(results_path = dirname(annot_file),
                                   # Will convert data.table automatically
                                   gr.snp = finemap_dat,
                                   keyword_query = Roadmap_query,
                                   limit_files=F,
                                   nThread = nThread)
      saveRDS(grl.roadmap, annot_file)
    }
    grl.roadmap.filt <- ROADMAP.merge_and_process_grl(grl.roadmap = grl.roadmap,
                                                      gr.snp = gr.snp,
                                                      n_top_tissues=n_top_roadmap,
                                                      sep = "\n")
    track.roadmap <- ROADMAP.track_plot(grl.roadmap.filt=grl.roadmap.filt,
                                        gr.snp=gr.snp,
                                        show_plot = F)
    TRKS[["Roadmap\nchromatin marks\ncell-types"]] <- track.roadmap
  }


  if(Nott_epigenome){
    try({
      #### Track: NOTT_2019 histogram  ####
      track.Nott_histo <- NOTT_2019.epigenomic_histograms(finemap_dat = finemap_dat,
                                                          locus_dir = locus_dir,
                                                          geom = "histogram",
                                                          plot_formula = "Cell_type ~.",
                                                          show_plot=F,
                                                          save_plot=F,
                                                          full_data=T,
                                                          return_assay_track=T,
                                                          binwidth=Nott_binwidth,
                                                          bigwig_dir=Nott_bigwig_dir,
                                                          save_annot=T,
                                                          as_ggplot = T,
                                                          strip.text.y.angle = strip.text.y.angle,
                                                          xtext=xtext,
                                                          nThread=nThread,
                                                          verbose=verbose)
      TRKS[["Nott (2019)\nread densities"]] <- track.Nott_histo +
        labs(y="Nott (2019)\nread densities")
    })

    if(Nott_show_placseq){
      #### Track: NOTT_2019 PLAC-seq  ####
      try({
        track.Nott_plac <- NOTT_2019.plac_seq_plot(finemap_dat = finemap_dat,
                                                   locus_dir=locus_dir,
                                                   title=locus,
                                                   show_regulatory_rects=Nott_regulatory_rects,
                                                   return_interaction_track=T,
                                                   show_arches=T,
                                                   save_annot=T,
                                                   strip.text.y.angle = strip.text.y.angle,
                                                   xtext=xtext,
                                                   as_ggplot = T,
                                                   nThread=nThread,
                                                   verbose=verbose)
        TRKS[["Nott (2019)\nPLAC-seq"]] <- track.Nott_plac +
          labs(y="Nott (2019)\nPLAC-seq")
      })
    }
  }

  # WARNING:: The order of these adjustment functions matters!
  ## Some of them reset the parameters of others

  #### Remove plots margins to save space ####
  TRKS <- PLOT.remove_margins(TRKS = TRKS,
                              finemap_dat = finemap_dat,
                              verbose = verbose)
  #### Make sure last plot has xtext ####
  TRKS <- PLOT.add_back_xtext(TRKS = TRKS,
                              verbose = verbose)
  #### Add vertical lines  ####
  if(!is.null(snp_group_lines)){
    TRKS <- PLOT.add_multitrack_lines(TRKS = TRKS,
                                      finemap_dat = finemap_dat,
                                      snp_groups = snp_group_lines,
                                      line_alpha = .8,
                                      remove_duplicated_UCS_Consensus = T,
                                      verbose = verbose)
  }


  ##### Iterate over different window sizes #####
  plot_list <- list()
  for(pz in plot.zoom){
    # try() Allows (X11) errors to occur and still finish the loop
    try({
      message("+>+>+>+>+ plot.zoom = ",pz," +<+<+<+<+")
      TRKS_zoom <- TRKS
      window_suffix <- PLOT.get_window_suffix(finemap_dat=finemap_dat,
                                              plot.zoom=pz)
      if((plot_full_window) & (!window_suffix %in% c("1x","all"))){
        #### Add zoom polygon ####
        TRKS_zoom[["zoom_polygon"]] <- PLOT.zoom_polygon(finemap_dat = finemap_dat,
                                                         genomic_units = genomic_units,
                                                         plot.zoom = pz)
        TRKS_zoom[[full_window_name]] <- PLOT.zoom_highlight(gg = TRKS_zoom[[full_window_name]],
                                                             finemap_dat = finemap_dat,
                                                             plot.zoom = pz)
      }
      #### Reorder tracks ####
      TRKS_zoom <- PLOT.reorder_tracks(TRKS = TRKS_zoom,
                                       track_order = track_order,
                                       verbose = verbose)

      if(window_suffix=="1x"){
        # This track becomes redundant when you don't zoom in at all.
        printer("+ PLOT:: Removing",full_window_name,"track @ zoom=1x",v=verbose)
        TRKS_zoom[[full_window_name]] <- NULL
      }
      #### Check track heights ####
      track_heights <- PLOT.check_track_heights(TRKS = TRKS_zoom,
                                                track_heights = track_heights,
                                                default_height = 1,
                                                verbose = verbose)
      #### Define plot.zoom limits ####
      TRKS_zoom <- PLOT.set_window_limits(TRKS = TRKS_zoom,
                                          finemap_dat = finemap_dat,
                                          plot.zoom = pz,
                                          verbose = verbose)
      #### Construct title ####
      n_snps <- if(dataset_type %in% names(TRKS_zoom)){
        paste0("n SNPs: ", nrow(ggplot2::ggplot_build(TRKS_zoom[[dataset_type]])$data[[2]]),", ")
      } else {NULL}
      title_text <- paste0(basename(locus_dir),"   (",n_snps,"zoom: ",window_suffix,")")

      #### Fuse all tracks ####
      TRKS_FINAL <- patchwork::wrap_plots(TRKS_zoom, ncol = 1) +
        patchwork::plot_layout(heights = track_heights) +
        patchwork::plot_annotation(title = title_text)

      #### Add plot to list of zoomed plots ####
      plot_list[[pz]] <- if(return_list) TRKS_zoom else TRKS_FINAL

      #### Save plot ####
      if(save_plot){
        plot_path <- file.path(locus_dir,paste("multiview",locus,LD_reference,window_suffix,plot_format,sep="."))
        printer("+ PLOT:: Saving plot ==>",plot_path)
        ggplot2::ggsave(filename = plot_path,
                        plot = TRKS_FINAL,
                        height = height,
                        width = width,
                        dpi = dpi,
                        bg = "transparent")
      }
      if(show_plot){print(TRKS_FINAL)}
    })
  } # End plot.zoom loop
  return(plot_list)
}


#### + SUPPORT FUNCTIONS + ####
PLOT.remove_margins <- function(TRKS,
                                finemap_dat,
                                verbose=T){
  printer("+ Removing subplot margins...", v=verbose)
  TRKS_std <- list()
  for(x in names(TRKS)){
    # print(x)
    # Some plots only had the labels tranformed, not the actual values.
    ## Check which ones are which and set the limits accordingly.
    TRKS_std[[x]] <-  suppressMessages(suppressWarnings(
      TRKS[[x]] +
        theme(plot.margin = unit(rep(0,4),"cm") )
    ))
  }
  return(TRKS_std)
}





PLOT.SNP_track_merged <- function(finemap_dat,
                                  yvar="-log10(P)",
                                  labels_subset = c("Lead","UCS","Consensus"),
                                  absolute_labels=F,
                                  label_type="rsid_only",
                                  sig_cutoff=5e-8,
                                  cutoff_lab=paste("P <",sig_cutoff),
                                  point_alpha=.5,
                                  show.legend=T,
                                  xtext=T,
                                  facet_formula="Method~.",
                                  dataset_type=NULL,
                                  genomic_units="POS",
                                  strip.text.y.angle=0,
                                  show_plot=F,
                                  verbose=T){
  # labels_subset = c("Lead","UCS","Consensus"); yvar="-log10(P)"; absolute_labels=F; label_type="rsid_only";  sig_cutoff=5e-8;
  # cutoff_lab=paste("P <",sig_cutoff); point_alpha=.5; show.legend=T; xtext=T;  facet_formula="Method~."; track_heights=NULL;
  if(endsWith(yvar,"PP")) {
    finemap_melt <- melt_finemapping_results(finemap_dat = finemap_dat,
                                             verbose = verbose)
    yvar <- "PP"
    sig_cutoff <- .95
    cutoff_lab <- paste("PP ≥",sig_cutoff)
    remove_duplicates <- F
    melt_methods <- T
    grouping_vars <- c("SNP","Method")
  } else {
    finemap_melt <- finemap_dat
    finemap_melt$Method <- if(is.null(dataset_type)) yvar else dataset_type
    cutoff_lab <- paste("P <",sig_cutoff)
    sig_cutoff <- -log10(sig_cutoff)
    remove_duplicates <- T
    melt_methods <- F
    grouping_vars <- c("SNP")
  }
  finemap_melt$Mb <- finemap_melt$POS/1000000

  # Plot
  snp_plot <- ggplot(data = finemap_melt,
                     aes_string(x=genomic_units, y=yvar, color="r2")) +
    # Bottom plot delineator
    geom_hline(yintercept=0) +
    geom_point(alpha=point_alpha, show.legend = show.legend) +
    scale_color_gradient(low="blue",high ="red", breaks=c(0,.5,1), limits=c(0,1)) +
    ## Sig cutoff line
    geom_hline(yintercept=sig_cutoff, alpha=.5, linetype=2, size=.5, color="black") +
    geom_text(aes(x=min(Mb), y=sig_cutoff*1.1),
              label=cutoff_lab,
              size=3,color="grey", hjust = 0) +
    labs(color=bquote(r^2),
         y=if(startsWith(yvar,"-log10")) bquote("-log"[10]~"(p)") else yvar ) +
    theme_classic() +
    facet_grid(facets =if(is.null(facet_formula)) facet_formula else  as.formula(facet_formula)) +
    theme(strip.text.y = element_text(angle=strip.text.y.angle)
          ) +
    guides(fill = guide_colourbar(barwidth = 1, barheight = 3))

  #  Choose breaks and ylimss
  if(startsWith(yvar,"-log")){
    snp_plot <- suppressMessages(
      snp_plot + scale_y_continuous(n.breaks = 3,
                                    limits =  c(0,-log10(min(finemap_dat$P))*1.1))
    )
  }else {
    snp_plot <- suppressMessages(
      snp_plot + scale_y_continuous(breaks = c(0,.5,1),
                                    limits = c(0,1.15)) +
        labs(y="Fine-mapping PP")
    )
  }
  if(!is.null(labels_subset)){
    snp_plot <- PLOT.add_snp_labels(snp_plot = snp_plot,
                               finemap_dat = finemap_dat,
                               yvar = yvar,
                               genomic_units = genomic_units,
                               labels_subset = labels_subset,
                               remove_duplicates = remove_duplicates,
                               melt_methods = melt_methods,
                               grouping_vars = grouping_vars,
                               show.legend = show.legend)
  }
  if(xtext==F){
    snp_plot <- snp_plot +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank())
  }
  if(show_plot) print(snp_plot)
  return(snp_plot)
}




PLOT.add_snp_labels <- function(snp_plot,
                               finemap_dat,
                               labels_subset=c("CS"),
                               yvar="PP",
                               genomic_units="Mb",
                               grouping_vars=c("SNP","Method"),
                               remove_duplicates=F,
                               melt_methods=T,
                               show.legend=T){
  if(melt_methods){
    finemap_melt <- melt_finemapping_results(finemap_dat = finemap_dat,
                                             verbose = F)
  } else {finemap_melt <- finemap_dat}
  snp_points <- construct_SNPs_labels(subset_DT = finemap_melt,
                                      labels_subset = c("Lead","UCS","Consensus"),
                                      remove_duplicates = remove_duplicates,
                                      grouping_vars = grouping_vars)
  if(melt_methods){
    snp_labels <- construct_SNPs_labels_separate(finemap_melt = finemap_melt,
                                                 labels_subset = labels_subset)
  } else {snp_labels <- snp_points}
  # Add SNP labels to plot
  snp_plot_labeled <- snp_plot +
    geom_point(data=snp_labels,
               pch=snp_labels$shape,
               fill=NA,
               size=snp_labels$size,
               color=snp_labels$color) +
    ### Background color label
    ggrepel::geom_label_repel(data=snp_labels,
                              aes(label=SNP),
                              color=NA,
                              # nudge_x = .5,
                              fill="black",
                              box.padding = .5,
                              label.padding = .25,
                              point.padding = .5,
                              # min.segment.length = 1,
                              # nudge_x = if(absolute_labels).5 else 0,
                              # nudge_y= if(absolute_labels).5 else 0,
                              label.size=NA,
                              alpha=.75,
                              seed = 1,
                              size = 3) +
    ### Foreground color label
    ggrepel::geom_label_repel(data=snp_labels,
                              aes(label=SNP),
                              segment.colour = snp_labels$color,
                              color="white",#labelSNPs_noDupes$color,
                              segment.alpha = .5,
                              box.padding = .5,
                              label.padding = .25,
                              point.padding = .5,
                              # min.segment.length = 1,
                              # nudge_x = if(absolute_labels).5 else NULL,
                              # nudge_y= if(absolute_labels).5 else NULL,
                              segment.size = .75,
                              fill = NA,
                              alpha=1,
                              seed = 1,
                              size = 3) +
    # Enhance the colors of SNPs with labeled background (to see them better)
    geom_point(data =snp_points,
               aes_string(x=genomic_units, y=yvar,color="r2"),
               alpha=1, pch=snp_points$shape,
               show.legend = show.legend)
  return(snp_plot_labeled)
}



PLOT.guess_genomic_units <- function(gg,
                                     decimals_default="Mb"){
  gp <- ggplot2::ggplot_build(gg)
  data_dims <- lapply(gp$data, nrow) %>% unlist()
  i <- which(data_dims==max(data_dims))[1]
  plot_dat <- gp$data[[i]]
  x_var <- c("xend","x")[c("xend","x") %in% colnames(plot_dat)][1]
  ### Guess based on  the presence of decimals (which implies non-POS units)
  ### Currently cannot distinguish between different non-POS units (Gb,Mb,Kb) so we just assume Mb.
  genomic_units <- if(max(plot_dat[[x_var]], na.rm = T)%%1==0) "POS" else decimals_default
  return(genomic_units)
}


PLOT.heights_dict <- function(){
  c("Genes"=.33,
    "GWAS full window"=.33,
    "QTL full window"=.33,
    "zoom_polygon"=.15,
    "GWAS"=.33,
    "QTL"=.33,
    "Fine-mapping"=1,
    "Roadmap\nchromatin marks\ncell-types"=1,
    "Nott (2019)\nread densities"=.5,
    "Nott (2019)\nPLAC-seq"=1)
}

PLOT.check_track_heights <- function(TRKS,
                                     track_heights,
                                     default_height=1,
                                     verbose=T){
  printer("+ Checking track heights...",v=verbose)
  heights_dict <- PLOT.heights_dict()

  if(is.null(track_heights)){
    track_heights <- heights_dict[names(TRKS)]
  } else {
    # Ensures n TRKS == n track_heights,
    ## even if user-supplied track_heights is > or < than n TRKS.
    track_heights <- track_heights[1:length(TRKS)]
  }
  track_heights <- track_heights[!is.na(track_heights)]

  if(length(track_heights)<length(TRKS)){
    n_missing <- length(TRKS)-length(track_heights)
    track_heights <- c(track_heights, rep(default_height,n_missing))
  }
  return(track_heights)
}

PLOT.reorder_tracks <- function(TRKS,
                                track_order=NULL,
                                verbose=T){
  printer("+ Reordering tracks...", v=verbose)
  track_order <- if(is.null(track_order)) names(PLOT.heights_dict()) else track_order

  ## Find user-given tracks that are actually available
  actual_track_order <- track_order[track_order %in% names(TRKS)]
  ## Add back in any tracks user might have missed
  actual_track_order <- unique(c(actual_track_order, names(TRKS)))
  TRKS <- TRKS[actual_track_order]
  return(TRKS)
}


#' Get window size limits for plot
#'
#' @family plot
#' @keywords internal
#' @examples
#' data("BST1");
#' xlims <- PLOT.get_window_limits(finemap_dat=BST1, plot.zoom=50000)
#' xlims <- PLOT.get_window_limits(finemap_dat=BST1, plot.zoom="all")
#' xlims <- PLOT.get_window_limits(finemap_dat=BST1, plot.zoom="5x")
PLOT.get_window_limits <- function(finemap_dat,
                                   index_as_center=T,
                                   plot.zoom=NULL,
                                   genomic_units="Mb",
                                   verbose=T){
  # plot.zoom <- c("all","1x","4x",5000);
  # plot.zoom <- c("5000");
  # finemap_dat=echolocatoR::BST1;
  # index_as_center=T; genomic_units="Mb"; verbose=T;
  #
  plot.zoom <- if(is.null(plot.zoom)) "1x" else plot.zoom
  plot.zoom[!is.na(plot.zoom)] <- plot.zoom;
  plot.zoom[!is.null(plot.zoom)] <- plot.zoom;
  # Iterate over list of zooms
  xlims_list <- lapply(plot.zoom, function(pz,
                                           .finemap_dat=finemap_dat,
                                           .index_as_center=index_as_center,
                                           .genomic_units=genomic_units,
                                           .verbose=verbose){
    printer("+ Inferring genomic limits for window:",pz,v=.verbose)
    # Zoom #x as  input
    if(.index_as_center) {
      middle_pos <- subset(.finemap_dat, leadSNP)$POS[1]
    } else {
      # Lead Pos isn't always dead middle if manual xlims were used during querying
      middle_pos <- .finemap_dat[as.numeric(round(nrow(.finemap_dat))/2),]$POS
    }

    if(grepl("x$",tolower(pz))){
      if(tolower(pz)=="1x") {
        min_limit <- min(.finemap_dat$POS, na.rm = T)
        max_limit <- max(.finemap_dat$POS, na.rm = T)
      } else {
        total_bp_span <- (max(.finemap_dat$POS, na.rm = T) - min(.finemap_dat$POS, na.rm = T))
        new_window <- total_bp_span / as.numeric(gsub("x","",pz))
        # Prevent extending beyond the borders of the data (producing blank space)
        min_limit <- middle_pos - as.integer(new_window/2)
        max_limit <- middle_pos + as.integer(new_window/2)
      }
    } else {
      # Basepairs as input
      if(is.null(pz)){
        min_limit <- min(.finemap_dat$POS, na.rm = T)
        max_limit <- max(.finemap_dat$POS, na.rm = T)
      } else {
        # 'all' as input
        if(pz=="all"){
          min_limit <- min(.finemap_dat$POS, na.rm = T)
          max_limit <- max(.finemap_dat$POS, na.rm = T)
        } else {
          min_limit <- middle_pos - as.integer(as.numeric(pz)/2)
          max_limit <- middle_pos + as.integer(as.numeric(pz)/2)
        }
      }
    }
    xlims <- c(min_limit, max_limit)
    if(.genomic_units=="Mb"){
      xlims <- xlims/1000000
    }
    return(xlims)
  }) %>% `names<-`(plot.zoom)

  # For backwards compatibility
  if(length(xlims_list)==1) return(xlims_list[[1]]) else return(xlims_list)
}


#' Determine the plot suffix indicating its window size
#'
#' @family plot
#' @keywords internal
#' @examples
#' data("BST1")
#' window_suffix <- get_window_suffix(finemap_dat=BST1, plot.zoom=1000)
#' window_suffix <- get_window_suffix(finemap_dat=BST1, plot.zoom=NULL)
#' window_suffix <- get_window_suffix(finemap_dat=BST1, plot.zoom="all")
#' window_suffix <- get_window_suffix(finemap_dat=BST1, plot.zoom="2x")
PLOT.get_window_suffix <- function(finemap_dat,
                                   plot.zoom,
                                   verbose=T){
  printer("+ PLOT:: Get window suffix...",v=verbose)
  window_suffix <- if(is.null(plot.zoom)){
    return(paste0(DescTools::RoundTo((max(finemap_dat$POS) - min(finemap_dat$POS))/1000, 100),"kb"))
  } else {
    if(is.character(plot.zoom)){
      if(plot.zoom=="all"){
        return("1x")
        # return(paste0(DescTools::RoundTo((max(finemap_dat$POS) - min(finemap_dat$POS))/1000, 100),"kb"))
      } else {
        return(plot.zoom)
      }
    } else {
      return(paste0(plot.zoom/1000,"kb"))
    }
  }
}




PLOT.set_window_limits <- function(TRKS,
                                   finemap_dat,
                                   plot.zoom,
                                   exceptions_str="*full window$|zoom_polygon|^Genes$",
                                   verbose=T){
  printer("+ Aligning xlimits for each subplot...",v=verbose)
  ### Exclude the purposefully unzoomed tracks
  exceptions <- grep(exceptions_str, names(TRKS), value=T)

  for(x in names(TRKS)){
    # print(x)
    trk <- TRKS[[x]]
    genomic_units <- PLOT.guess_genomic_units(gg = trk)
    xlims <- PLOT.get_window_limits(finemap_dat=finemap_dat,
                                    plot.zoom = if(x %in% exceptions) "1x" else plot.zoom,
                                    genomic_units=genomic_units,
                                    verbose=F)
    unit_divisor <- if(genomic_units=="Mb") 1 else 1000000

    suppressMessages(suppressWarnings(
      TRKS[[x]] <- trk +
        scale_x_continuous(labels = function(x)x/unit_divisor,
                           limits = xlims)
    ))
  }
  return(TRKS)
}






#' Multi-fine-map summary dot plot
#'
#' @family plot
#' @examples
#' data("BST1")
#' gp <- PLOT.dot_summary(finemap_dat=BST1)
PLOT.dot_summary <- function(finemap_dat,
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



guess_pvalue_col <- function(dat,
                             QTL_prefix=NULL){
  p_options <- c("p","p-value","p-values","pvalue","pvalues","pval")
  p_options <- c(p_options, stringr::str_to_sentence(p_options))
  pval_col <- paste0(QTL_prefix,p_options)[paste0(QTL_prefix,p_options) %in% colnames(dat)][1]
  return(pval_col)
}



PLOT.zoom_polygon <- function(finemap_dat,
                              genomic_units="Mb",
                              plot.zoom="5x",
                              alpha=.15,
                              verbose=T){
  printer("+ Constructing zoom polygon...",v=verbose)
  xlims_orig <- PLOT.get_window_limits(finemap_dat = finemap_dat,
                                       genomic_units = genomic_units,
                                       verbose = F)
  xlims_zoom <- PLOT.get_window_limits(finemap_dat = finemap_dat,
                                       plot.zoom = plot.zoom,
                                       genomic_units = genomic_units,
                                       verbose = F)
  positions <- data.frame(x=c(xlims_zoom, rev(xlims_orig)),
                          y=c(c(1,1),c(0,0)))

  zp <- ggplot(positions, aes(x = x, y = y)) +
    geom_polygon(fill="blue",
                 color="transparent",
                 alpha=alpha, show.legend = F) +
    # scale_fill_gradient(limits=c(0.75, 4), low = "lightgrey", high = "red") +
    theme_void() +
    theme(plot.margin = unit(rep(0,4),"cm") )
  return(zp)
}


PLOT.zoom_highlight <- function(gg,
                                finemap_dat,
                                plot.zoom="5x",
                                alpha=.15,
                                verbose=T){
  printer("+ Highlighting zoom origin...",v=verbose)
  genomic_units <- PLOT.guess_genomic_units(gg = gg)
  xlims_zoom <- PLOT.get_window_limits(finemap_dat = finemap_dat,
                                       plot.zoom = plot.zoom,
                                       genomic_units = genomic_units,
                                       verbose = F)
  rect_dat <- data.frame(x=c(xlims_zoom, rev(xlims_zoom)), y=c(-Inf,-Inf,Inf,Inf))
  gg2 <- gg +
    ggplot2::geom_polygon(data = rect_dat,
                          aes(x=x, y=y),
                          color="transparent",
                          alpha=alpha,
                          fill="blue")
  return(gg2)
}



#' Plot gene/transcript models
#'
#' @source
#' \href{https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html}{ensembld tutorial}
#' \href{https://bioconductor.org/packages/devel/bioc/vignettes/Gviz/inst/doc/Gviz.html#45_GeneRegionTrack}{Gvix tutorial}
#' \href{http://bioconductor.org/packages/devel/bioc/vignettes/ggbio/inst/doc/ggbio.pdf}{ggbio tutorial}
#'
#' @examples
#' data("LRRK2")
#' gene_track <- PLOT.transcript_model_track(finemap_dat=LRRK2)
PLOT.transcript_model_track <- function(finemap_dat,
                                        max_transcripts=1,
                                        remove_RP11=T,
                                        show.legend=T,
                                        show_plot=F,
                                        fill="skyblue",
                                        shape=c( "arrow", "box", "ellipse","smallArrow"),
                                        transcriptAnnotation = c("symbol","transcript"),
                                        collapseTranscripts=c(F,T,"longest"),
                                        stacking=c("squish","hide", "dense", "pack","full"),
                                        method="ggplot",
                                        xtext=T,
                                        verbose=T){
  # method="ggbio"; verbose=T; stacking="squish";  fill="blue"; transcriptAnnotation = c("symbol","transcript");  collapseTranscripts=c(F,T,"longest"); shape="arrow";
  gr.snp_CHR <- biovizBase::transformDfToGr(data = finemap_dat,
                                            seqnames = "CHR",
                                            start = "POS", end = "POS")
  #### Gviz ####
  if(method=="gviz"){
    library(Gviz)
    library(EnsDb.Hsapiens.v75)
    suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr.snp_CHR) <- "NCBI")
    printer("+ PLOT:: Gene Model Track",v=verbose)
    edb <-  EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
    tx <- ensembldb::getGeneRegionTrackForGviz(edb,
                                               filter = GRangesFilter(value = gr.snp_CHR, feature = "tx"))
    tx_lengths <- transcriptLengths(edb,
                                    filter = GRangesFilter(value = gr.snp_CHR, feature = "tx"))
    tx$tx_len <- setNames(tx_lengths$tx_len,tx_lengths$tx_id)[tx$transcript]
    tx <- subset(tx, !startsWith(as.character(symbol),"RP11"))
    # gtrack <- GenomeAxisTrack()
    tx.models <- Gviz::GeneRegionTrack(tx,
                                       name = "Gene Model",
                                       fill=fill,
                                       stacking = stacking[1])
    # ggplotify::as.ggplot(tx.models)
    Gviz::plotTracks(list(tx.models),
                     transcriptAnnotation = transcriptAnnotation[1],
                     background.title = "grey20",
                     collapseTranscripts=collapseTranscripts[1],
                     shape = shape[1]  )
  } else {
    #### ggbio ####
    # Warning:: MUST load the full package bc
    # it loads other necessary packages into the namespace.
    library(Homo.sapiens)
    # columns(Homo.sapiens)
    # library(ggbio)
    suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr.snp_CHR) <- "UCSC")
    if("pals" %in% row.names(installed.packages())){
      palette <- unname(pals::alphabet())
    }
    if(collapseTranscripts[1]==T){
      fill_var <- NULL
      stat_opt <- "reduce"
      names.expr <- "SYMBOL"
    } else {
      fill_var <- "SYMBOL"
      stat_opt <- "identity"
      names.expr <- "SYMBOL (TXNAME)"
    }
    tx.models <- ggbio::autoplot( Homo.sapiens::Homo.sapiens,
                                 which=gr.snp_CHR,
                                 aes_string(fill=fill_var, color=fill_var),
                                 columns = c("SYMBOL","TXNAME"),
                                 names.expr = names.expr,
                                 stat = stat_opt) +
      theme_classic() +
      labs(y="Transcript")

    if(xtext==F){
      tx.models <- tx.models +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank())
    }
    # if(genomic_units=="Mb"){
    #   tx.models <- tx.models +
    #     scale_x_continuous( labels=function(x)x/1000000)
    # }
    if("pals" %in% row.names(installed.packages())){
      # Give it some nice fancy colors
      palette <- unname(pals::alphabet())
      tx.models <- tx.models +
        scale_fill_manual(values = palette) +
        scale_color_manual(values = palette)
    }
    if(method=="ggplot"){
      tx.models <- tx.models@ggplot
    }
  }
  if(show_plot) print(tx.models)
  return(tx.models)
}





#' Add vertical lines
#'
#' Adds vertical lines indicates key SNPs
#' @family plot
#' @keywords internal
PLOT.add_multitrack_lines <- function(TRKS,
                                      finemap_dat,
                                      snp_groups=c("Lead","UCS","Consensus"),
                                      line_alpha=1,
                                      line_size=.3,
                                      remove_duplicated_UCS_Consensus=T,
                                      verbose=F){
  printer("+ Adding vertical lines to highlight SNP groups...",v=verbose)
  snp_labels <- construct_SNPs_labels(subset_DT = finemap_dat,
                                      labels_subset = snp_groups,
                                      remove_duplicates = F,
                                      grouping_vars = "SNP")
  ## Remove any lines that are both UCS and Consensus.
  ### Otherwise, when line_alpha<1, they'll make a muddled green/gold color.
  if(remove_duplicated_UCS_Consensus){
    UCS_Cons <- subset(snp_labels, type%in%c("UCS","Consensus"))
    UCS_Cons <- UCS_Cons[base::duplicated(UCS_Cons$SNP),]
    snp_vlines <- rbind(UCS_Cons, subset(snp_labels, type=="Lead"))
  } else {
    snp_vlines <- snp_labels
  }

  TRKS_lines <- list()
  # iterate over each ggplot object in a named list
  for(x in names(TRKS)){
    # printer(x, v=verbose)
    trk <- TRKS[[x]]
    genomic_units <- PLOT.guess_genomic_units(gg = trk)
    if("Consensus" %in% snp_groups){
      consensus.pos <- subset(snp_vlines, type=="Consensus")[[genomic_units]]
      trk <- trk +  # Consensus
        geom_vline(xintercept = consensus.pos, color="goldenrod1",
                   alpha=line_alpha, size=line_size, linetype='solid')
    }
    if("UCS" %in% snp_groups){
      UCS.pos <- subset(snp_vlines, type=="UCS")[[genomic_units]]
      trk <- trk +  # Consensus
        geom_vline(xintercept = UCS.pos, color="green3",
                   alpha=line_alpha, size=line_size, linetype='dashed')
    }
    if("Lead" %in% snp_groups){
      lead.pos <- subset(snp_vlines, type=="Lead")[[genomic_units]]
      trk <- trk +  # Consensus
        geom_vline(xintercept = lead.pos, color="red",
                   alpha=line_alpha, size=line_size, linetype='dashed')
    }
    TRKS_lines[[x]] <- trk
  }
  return(TRKS_lines)
}




PLOT.add_back_xtext <- function(TRKS,
                                verbose=T){
  print("+ Ensuring last track shows genomic units...",v=verbose)
  i <- length(TRKS)
  unit_divisor <- if(PLOT.guess_genomic_units(gg = TRKS[[i]])=="Mb") 1 else 1000000
  TRKS[[i]] <- suppressMessages(
    TRKS[[i]] + theme(axis.text.x = element_text(),
                      axis.title.x = element_text()) +
      labs(x="Mb") +
      scale_x_continuous(labels = function(x)x/unit_divisor)
  )
  return(TRKS)
}




construct_SNPs_labels <- function(subset_DT,
                                  labels_subset=c("Lead","CS","UCS","Consensus"),
                                  remove_duplicates=T,
                                  grouping_vars=c("SNP"),
                                  merge_with_input=F,
                                  verbose=F){
  printer("+ PLOT:: Constructing SNP labels...", v=verbose)
  labelSNPs <- data.table::data.table()
  subset_DT <- data.table::as.data.table(subset_DT)
  subset_DT$Mb <- subset_DT$POS/1000000

  ## BEFORE fine-mapping
  if("lead" %in% tolower(labels_subset)){
    lead_snps <- subset(subset_DT %>% arrange(P), leadSNP == T)
    lead_snps$type <- "Lead"
    lead_snps$color <- "red"
    lead_snps$shape <- 9# 18
    lead_snps$size <- 3
    labelSNPs <- rbind(labelSNPs, lead_snps, fill=T)
  }
  if(("cs" %in% tolower(labels_subset)) & ("CS" %in% colnames(subset_DT) ) ){
    # AFTER fine-mapping
    CS_snps = subset(subset_DT, CS>0)
    if(dim(CS_snps)[1]>0){
      CS_snps$type <- "CS"
      CS_snps$color<- "green3"
      CS_snps$shape <- 16
      CS_snps$size=5
      labelSNPs <- rbind(labelSNPs, CS_snps, fill=T)
    }
  }
  if("ucs" %in% tolower(labels_subset)){
    # AFTER fine-mapping
    UCS_snps = subset(subset_DT, Support>0)
    if(dim(UCS_snps)[1]>0){
      UCS_snps$type <- "UCS"
      UCS_snps$color<- "green3"
      UCS_snps$shape <- 16
      UCS_snps$size=5
      labelSNPs <- rbind(labelSNPs, UCS_snps, fill=T)
    }
  }
  if(("consensus" %in% tolower(labels_subset)) & ("Consensus_SNP" %in% colnames(subset_DT) ) ){
    # Conensus across all fine-mapping tools
    consensus_SNPs <- subset(subset_DT, Consensus_SNP==T)
    if(dim(consensus_SNPs)[1]>0){
      consensus_SNPs$type <- "Consensus"
      consensus_SNPs$color <- "darkgoldenrod1"
      consensus_SNPs$shape <- 16
      consensus_SNPs$size=6
      labelSNPs <- rbind(labelSNPs, consensus_SNPs, fill=T)
    }
  }
  # If there's duplicates only show the last one
  labelSNPs$rowID <- 1:nrow(labelSNPs)
  if(remove_duplicates){
    labelSNPs <- labelSNPs %>%
      dplyr::group_by(.dots=grouping_vars) %>%
      dplyr::arrange(rowID) %>%
      dplyr::slice(n())
  }
  labelSNPs$type <- factor(labelSNPs$type, levels = c("UCS","CS","Consensus","Lead"), ordered = T)
  labelSNPs <- dplyr::arrange(labelSNPs, type)
  if("Method"%in%colnames(labelSNPs)){
    labelSNPs$Method <- factor(labelSNPs$Method, unique(labelSNPs$Method), ordered = T)
  }


  ## Merge with input df
  if(merge_with_input){
    plot_dat <- data.table::merge.data.table(subset_DT,
                                             subset(labelSNPs, select=c("SNP","Method","type","color","shape","size")),
                                             by = c("SNP","Method")[c("SNP","Method") %in% colnames(subset_DT)],
                                             all.x = T)
    plot_dat[is.na(plot_dat$color),"color"] <- "transparent";
    plot_dat[is.na(plot_dat$shape),"shape"] <- 16;
    plot_dat[is.na(plot_dat$size),"size"] <- 3;
    return(plot_dat)
  }else {return(as.data.frame(labelSNPs))}
}



construct_SNPs_labels_separate <- function(finemap_melt,
                                           labels_subset="CS"){
  # Process means and non-means separately so you can add extra labels to mean
  finemap_labels <- construct_SNPs_labels(subset_DT = finemap_melt,
                                          labels_subset = labels_subset,
                                          grouping_vars = c("SNP","Method")) %>%
    subset(Method!="mean")
  mean_labels <- construct_SNPs_labels(subset_DT = finemap_melt,
                                       labels_subset = c("Lead","UCS","Consensus"),
                                       grouping_vars = c("SNP","Method"))  %>%
    subset(Method=="mean")
  mod_labels <- rbind(finemap_labels, mean_labels)
  return(mod_labels)
}



PLOT.get_max_histogram_height <- function(gg,
                                          round_to=NULL,
                                          verbose=T){
  if(tolower(class(gg)[1])=="ggbio") gg <- gg@ggplot
  printer("+ PLOT:: Calculating max histogram height",v=verbose)
  dat <- ggplot2::ggplot_build(gg)$data[[1]]
  max_height <- max(dat$ymax)
  if(!is.null(round_to)){
    max_height <- DescTools::RoundTo(max_height, round_to)
  }
  return(max_height)
}





