#
# ^^^^^^^^^^^^^^ Nott et al. (2019) # ^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^ single-nucleus human brain epigenomic dataset # ^^^^^^^^^^^^^^
# Nott, Alexi, Inge R. Holtman, Nicole G. Coufal, Johannes C.M. Schlachetzki, Miao Yu, Rong Hu, Claudia Z. Han, et al. “Cell Type-Specific Enhancer-Promoter Connectivity Maps in the Human Brain and Disease Risk Association.” Science 0793, no. November (2019): 778183. https://doi.org/10.1101/778183.
#  ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^

# UCSC Nott tracks
# https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf


# UCSC Genome Browser: Command Line Utilities:
# http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/



#' Plot brain cell-specific epigenomic data
#'
#' Brain cell-specific epigenomic data from Nott et al. (2019).
#' @inheritParams finemap_pipeline
#' @keywords internal
#' @family NOTT_2019
#' @source
#' \href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
#' \url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
#' @examples
#' data("BST1"); data("locus_dir");
#' track.Nott_histo <- NOTT_2019.epigenomic_histograms(finemap_dat = BST1, locus_dir = locus_dir, save_plot=F, return_assay_track=T, save_annot=F)
NOTT_2019.epigenomic_histograms <- function(finemap_dat,
                                            locus_dir,
                                            show_plot=T,
                                            save_plot=T,
                                            full_data=T,
                                            return_assay_track=F,
                                            binwidth=200,
                                            density_adjust=.2,
                                            plot.zoom="1x",
                                            strip.text.y.angle=90,
                                            xtext=T,
                                            geom="density",
                                            plot_formula="Cell_type ~.",
                                            fill_var="Assay",
                                            bigwig_dir=NULL,
                                            genomic_units="Mb",
                                            as_ggplot=T,
                                            nThread=4,
                                            save_annot=F,
                                            verbose=T){
  printer("NOTT_2019:: Creating epigenomic histograms plot",v=verbose)
  # library(BiocGenerics)
  # library(GenomicRanges)
  # library(ggbio)
  # show_plot=T;save_plot=T;full_data=T;return_assay_track=F;binwidth=2500; geom="histogram"; plot_formula="Cell_type ~."; show_regulatory_rects=T;  bigwig_dir=NULL; verbose=T; nThread=4;
  # finemap_dat=echolocatoR::LRRK2; plot.zoom=500000; fill_var="Assay"; density_adjust=.2; strip.text.y.angle=0;

  # Import BigWig annotation files
  bigWigFiles <- echolocatoR::NOTT_2019.bigwig_metadata
  # Some bigWig files were initially loaded to UCSC GB, but then later taken down by the authors....
  # However I saved these files on Minerva beforehand.
  bigWigFiles <- subset(bigWigFiles, UCSC_available=="T")
  bigWigFiles <- dplyr::mutate(bigWigFiles, cell_type = gsub(" ",".",cell_type))
  # Convert finemap data to granges
  dat <- finemap_dat
  dat$seqnames <- dat$CHR
  dat$start.end <- dat$POS
  gr.dat <- GenomicRanges::makeGRangesFromDataFrame(df = dat,
                                                    seqnames.field = "seqnames",
                                                    start.field = "start.end",
                                                    end.field = "start.end",
                                                    keep.extra.columns = T)
  # ! IMPORTANT !: Needs to be in chr1 format in order to query!
  GenomeInfoDb::seqlevelsStyle(gr.dat) <- "UCSC"
  printer("NOTT_2019:: Importing bigWig subsets from UCSC...", v=verbose)
  bw.grlist <- parallel::mclapply(1:nrow(bigWigFiles), function(i){
    if(!is.null(bigwig_dir)){
      bw.file <- file.path(bigwig_dir,paste0(bigWigFiles$long_name[i],".ucsc.bigWig"))
    } else { bw.file <- bigWigFiles$data_link[i]}

    bw.name <- gsub("_pooled|pooled_","",bigWigFiles$name[i])
    printer("+ NOTT_2019:: Importing...",paste0("[",i,"]"),bw.name)
    bw.filt <- import.bw.filt(bw.file=bw.file,
                              gr.dat=gr.dat,
                              full_data=full_data)
    bw.filt$Cell_type <- bigWigFiles$cell_type[i]
    bw.filt$Assay <- bigWigFiles$assay[i]
    bw.filt$Experiment <- gsub("_"," ",bw.name)
    return(bw.filt)
  }, mc.cores = nThread)
  bw.cols <- bigWigFiles$name
  # names(bw.grlist) <- bw.cols
  bw.gr <- unlist(GenomicRanges::GRangesList(bw.grlist))
  bw.gr$Assay <- gsub("atac","ATAC",bw.gr$Assay)
  bw.gr$Cell_type <- gsub("oligodendrocytes","oligo",bw.gr$Cell_type)

  xlims <- PLOT.get_window_limits(finemap_dat = finemap_dat,
                                   plot.zoom = plot.zoom,
                                   genomic_units = "POS")
  bw.gr <- subset(bw.gr,
                  GenomicRanges::seqnames(bw.gr)==paste0("chr",gsub("chr","",finemap_dat$CHR[1])) &
                  GenomicRanges::start(bw.gr)>=xlims[1] &
                  GenomicRanges::end(bw.gr)<=xlims[2])
  # merge into a single granges object
  # gr.snp <- Reduce(function(x, y) GenomicRanges::merge(x, y, all.x=T),
  #                  append(bw.grlist, gr.dat))
  # GenomicRanges::findOverlaps(query = gr.dat,
  #                             subject = bw.gr)

  if(save_annot){
    annot_file <- annotation_file_name(locus_dir = locus_dir,
                                       lib_name = "Nott_2019.epigenomics")
    saveRDS(bw.gr, annot_file)
  }

  PEAKS <- NOTT_2019.get_epigenomic_peaks(nThread = nThread,
                                          verbose = verbose)
  gr.peaks <- GRanges_overlap(dat1 = bw.gr,
                              dat2 = PEAKS)
  GenomicRanges::mcols(gr.peaks)[,c("chr","start","end","score")] <- NULL
  gr.peaks <- unique(gr.peaks)
  # Adjust line width to make sure bars aren't just all white
  line_width <- 0#binwidth/1000 * .25

  #### Density/Histogram plot ####
  color_dict <- assay_color_dict()
  nott_tracks <-  suppressWarnings(
    ggbio::autoplot(object=bw.gr,
                    geom=geom,
                    binwidth=binwidth,
                    alpha=.7,
                    position="stack",
                    adjust=density_adjust,
                    color="white",
                    size=line_width,
                    aes_string(fill=fill_var), show.legend=T)
  )
  # Pause and calculate max histo height
  max_height <- PLOT.get_max_histogram_height(gg=nott_tracks)
  rect_height <- max_height / if(geom=="density") 8 else 10
  gr.peaks$y <- 0 - rect_height
  if(gsub(" ","",plot_formula)=="Cell_type~."){
    gr.peaks$Assay <- "peaks"
  }
  nott_tracks <- nott_tracks +
    ggbio::geom_rect(gr.peaks,
                     stat="identity",
                     # position="stack",
                     rect.height= rect_height,
                     aes_string(y="y", fill=fill_var),
                     alpha=.25,
                     hjust=1,
                     color="transparent") +
    facet_grid(facets = formula(plot_formula)) +
    theme_classic() +
    theme(legend.position="right",
          strip.text.y = element_text(angle = strip.text.y.angle)) +
    scale_y_continuous(n.breaks = 3) +
    scale_fill_manual(values = color_dict)

  #### Make adjustments ####
  if(genomic_units=="Mb"){
    printer("++ Converting label units to Mb...",v=verbose)
    nott_tracks <- suppressMessages(
      nott_tracks +
        ggbio::scale_x_sequnit(unit = "Mb")
    )
  }
  if(xtext==F){
    nott_tracks <- nott_tracks +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank())
  }
  if(return_assay_track){
    if(as_ggplot) return(nott_tracks@ggplot) else return(nott_tracks)
  }



  #### (optional) Fine-mapping tracks ####
  TRACKS_list <- list(
    "GWAS"=ggbio::plotGrandLinear(obj = gr.dat, aes(y=-log10(P), color=-log10(P))),
    "Fine_mapping"=ggbio::plotGrandLinear(obj = gr.dat, aes(y=mean.PP, color=mean.PP)) +
      scale_color_viridis_c(),
    "Nott_etal_2019" = nott_tracks
  )
  # Fuse all tracks
  params_list <- list(title = paste0(gene," [",length(GenomicRanges::seqnames(gr.dat))," SNPs]"),
                      track.bg.color = "transparent",
                      track.plot.color = "transparent",
                      label.text.cex = .7,
                      label.bg.fill = "grey12",
                      label.text.color = "white",
                      label.text.angle = 0,
                      label.width = unit(5.5, "lines"),
                      # xlim = c(min(start(gr.snp)), max(start(gr.snp))),
                      heights = c(.3,.3,1)
                      )
  tracks <- get("tracks", asNamespace("ggbio"))
  trks <- suppressWarnings(do.call("tracks", append(TRACKS_list, params_list)))
  # add lines
  lead.pos <- subset(finemap_dat,leadSNP)$POS
  consensus.pos <- subset(finemap_dat, Consensus_SNP==T)$POS

  trks_plus_lines <- trks +
    # theme_bw() +
    # ggbio::theme_genome() +
    theme(strip.text.y = element_text(angle = 0),
          strip.text = element_text(size = 9),
          panel.background = element_rect(fill = "white", colour = "black", linetype = "solid")) +
    geom_vline(xintercept = consensus.pos, color="goldenrod2", alpha=1, size=.1, linetype='solid') +
    geom_vline(xintercept = lead.pos, color="red", alpha=1, size=.1, linetype='solid')

  if(show_plot){ print(trks_plus_lines) }
  if(save_plot){
    save_path <- file.path(locus_dir,"annotations",
                           paste0(basename(locus_dir),"_Glass.snEpigenomics.png") )
    dir.create(dir(save_path), showWarnings = F, recursive = T)
    ggbio::ggsave(save_path,
                   plot = trks_plus_lines, dpi=400, height = 15, width = 8,
                   bg = "transparent")
  }
  if(as_ggplot) return(trks_plus_lines@ggplot) else return(trks_plus_lines)
}






assay_color_dict <- function(){
  color_dict <- c("ATAC"="magenta",
                  "H3K27ac"="blue",
                  "H3K4me3"="turquoise",
                  "peaks"="black")
  return(color_dict)
}



#' Get cell type-specific superenhancer data
#'
#' Brain cell-specific epigenomic data from Nott et al. (2019).
#' @keywords internal
#' @family NOTT_2019
#' @source
#' \href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
#' \url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
NOTT_2019.superenhancers <- function(){
  annot_sub <- subset(echolocatoR::NOTT_2019.superenhancer_interactome,
                      chr== paste0("chr",unique(finemap_dat$CHR)) &
                      start>=min(finemap_dat$POS) &
                      end<=max(finemap_dat$POS) )
  if(nrow(annot_sub)>0){
    merged_dat <- data.table:::merge.data.table(finemap_dat %>%
                                                 dplyr::mutate(chr=paste0("chr",CHR),
                                                               start=as.numeric(POS)) %>%
                                                 data.table::data.table(),
                                               data.table::data.table(s6),
                                               by = c("chr","start"))
  }
  return(merged_dat)
}




#' Get cell type-specific promoter/emhancer/interactome data
#'
#' Brain cell-specific epigenomic data from Nott et al. (2019).
#' @keywords internal
#' @family NOTT_2019
#' @source
#' \href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
#' \url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
NOTT_2019.get_promoter_interactome_data <- function(finemap_dat){
  # Subset to window
  annot_sub <- echolocatoR::NOTT_2019.interactome$H3K4me3_around_TSS_annotated_pe %>%
    dplyr::rename(chr=Chr, start=Start, end=End) %>%
    subset(chr== paste0("chr",unique(finemap_dat$CHR)) &
           start>=min(finemap_dat$POS) &
           end<=max(finemap_dat$POS) )
  return(annot_sub)
}




#' Get promoter cell types
#'
#' Brain cell-specific epigenomic data from Nott et al. (2019).
#' @keywords internal
#' @family NOTT_2019
NOTT_2019.get_promoter_celltypes <- function(annot_sub,
                                             marker_key){
  promoter.cols <- grep("*_active_promoter",  colnames(annot_sub), value = T)
  logical.list <- colSums(annot_sub[,promoter.cols])>0
  promoter_celltypes <- gsub("\\_.*","", promoter.cols[as.logical(logical.list)] )
  promoter_celltypes <- as.character(marker_key[promoter_celltypes])
  promoter_celltypes <- paste(promoter_celltypes,collapse="; ")
  return(promoter_celltypes)
}




#' Import cell type-specific interactomes
#'
#' Brain cell-specific epigenomic data from Nott et al. (2019).
#' @keywords internal
#' @family NOTT_2019
#' @source
#' \href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
#' \url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
NOTT_2019.get_interactome <- function(annot_sub,
                                      top.consensus.pos,
                                      marker_key,
                                      verbose=T){
  printer("+ NOTT_2019:: Getting interactome data.", v=verbose)
  interact.cols <- grep("*_interactions", colnames(annot_sub), value = T)
  interact.DT <- lapply(interact.cols, function(column){
    coords <- strsplit(annot_sub[,column][[1]], ",")
    coord.dt <- lapply(coords, function(coord, .column=column){
      data.table::data.table(Interaction=.column,
                             Cell_type=marker_key[gsub("\\_.*","",.column)],
                             Coordinates=coord)
    }) %>% data.table::rbindlist()
    return(coord.dt)
  } )  %>% data.table::rbindlist()
  interact.DT <- subset(interact.DT, !is.na(Coordinates) & Coordinates!="") %>%
    tidyr::separate(col = Coordinates,
                    into=c("chr","Start","End"), sep = ":|-") %>%
    tidyr::separate(col = Interaction, into=c("Marker","Element",NA), sep="_", remove = F )
  interact.DT <- interact.DT %>%
    # Standardize CHR (NCBI format)
    dplyr::mutate(chr=gsub("chr","",chr),
                  Cell_type_interaction=paste(Cell_type,"-",Element))
  interact.DT$Cell_type <- interact.DT$Cell_type %>% as.character()
  interact.DT$Start <- as.numeric(interact.DT$Start)
  interact.DT$End <- as.numeric(interact.DT$End)
  # Summarise distance from different celltype enhancer interactions
  summarise_top.consensus.dist <- interact.DT %>%
    dplyr::mutate(top.consensus.dist=End - top.consensus.pos) %>%
    dplyr::group_by(Cell_type) %>%
    dplyr::summarise(top.consensus.dist = mean(top.consensus.dist))
  # print(summarise_top.consensus.dist)
  return(interact.DT)
}




#' Import cell type-specific interactomes
#'
#' Brain cell-specific epigenomic data from Nott et al. (2019).
#' @keywords internal
#' @family NOTT_2019
#' @source
#' \href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
#' \url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
NOTT_2019.get_interactions <- function(finemap_dat,
                                       as.granges=F){
  NOTT_2019.interactome <- echolocatoR::NOTT_2019.interactome
  selected_sheets <- grep("interactome$",names(NOTT_2019.interactome), value = T)
  interactomes <- lapply(selected_sheets, function(s){
    printer("Importing",s,"...")
    # Read the sheet you want
    dat <- echolocatoR::NOTT_2019.interactome[[s]]
    dat$Name <- s
    return(dat)
  }) %>% data.table::rbindlist()
  # Anchor 1
  interactomes.anchor1 <- GRanges_overlap(dat1 = finemap_dat,
                                          chrom_col.1 = "CHR",
                                          start_col.1 = "POS",
                                          end_col.1 = "POS",
                                          dat2 = interactomes,
                                          chrom_col.2 = "chr1",
                                          start_col.2 = "start1",
                                          end_col.2 = "end1" )
  interactomes.anchor1$Anchor <-1
  # Anchor 2
  interactomes.anchor2 <- GRanges_overlap(dat1 = finemap_dat,
                                          chrom_col.1 = "CHR",
                                          start_col.1 = "POS",
                                          end_col.1 = "POS",
                                          dat2 = interactomes,
                                          chrom_col.2 = "chr2",
                                          start_col.2 = "start2",
                                          end_col.2 = "end2" )
  interactomes.anchor2$Anchor <-2
  # Merge
  interactomes.anchor <- c(interactomes.anchor1, interactomes.anchor2)
  # Modify
  interactomes.anchor$Assay <- "PLAC"
  interactomes.anchor <- cbind(interactomes.anchor,
                               tidyr::separate(data.frame(interactomes.anchor), Name,
                                               into=c("Cell_type","Element"))[,c("Cell_type","Element")])
  cell_dict <- c("Microglia"="microglia",
                 "Neuronal"="neurons",
                 "Oligo"="oligo")
  interactomes.anchor$Cell_type <- cell_dict[interactomes.anchor$Cell_type]

  if(as.granges){
    interactomes.anchor <- GenomicRanges::makeGRangesFromDataFrame(interactomes.anchor, keep.extra.columns = T)
  }
  return(interactomes.anchor)
}




# Import cell type-specific epigenomic peaks
#
# Brain cell-specific epigenomic data from Nott et al. (2019).
# @keywords internal
# @family NOTT_2019
# @source
# \href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
# \url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
# @examples
# \dontrun{
# PEAKS.merged <- NOTT_2019.get_epigenomic_peaks(peak.dir="/pd-omics/data/Nott_2019/peaks", narrow_peaks=T, broad_peaks=F)
# }
NOTT_2019.get_epigenomic_peaks_macs2 <- function(peak.dir="/pd-omics/data/Nott_2019/peaks",
                                                 narrow_peaks=T,
                                                 broad_peaks=T,
                                                 nThread=4){
  bigWigFiles <- echolocatoR::NOTT_2019.bigwig_metadata
  peak_types <- c(ifelse(narrow_peaks,".narrowPeak$", NA),
                  ifelse(broad_peaks,"_broad.bed12$", NA))
  peak_types <- peak_types[!is.na(peak_types)]
  peaks.paths <- list.files(peak.dir,
                            pattern =  paste(peak_types, collapse = "|"),
                            full.names = T,
                            recursive = T)
  PEAKS <- MACS2.import_peaks(peaks.paths = peaks.paths,
                              as_granges = T)
  PEAKS <- parallel::mclapply(PEAKS, function(peak){
    pk.name <- gsub(".ucsc_narrowPeak1|.ucsc_broadRegion1","",peak$name[1])
    meta <- subset(bigWigFiles, long_name==pk.name)
    peak$Cell_type <- meta$cell_type
    peak$Assay <- meta$assay
    peak$Fresh_frozen <- meta$fresh_frozen
    peak$Marker  <- meta$marker
    return(peak)
  },mc.cores = nThread) %>% GenomicRanges::GRangesList()
 PEAKS.merged <- unlist(PEAKS)
 PEAKS.merged$peak_type <- PEAKS.merged
 return(PEAKS.merged)
}




#' Download cell type-specific epigenomic peaks
#'
#' API access to brain cell type-specific epigenomic peaks (bed format)
#' from Nott et al. (2019).
#' @keywords internal
#' @family NOTT_2019
#' @source
#' \href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
#' \url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
#' @examples
#' PEAKS <- NOTT_2019.get_epigenomic_peaks(nThread=1)
NOTT_2019.get_epigenomic_peaks <- function(assays=c("ATAC","H3K27ac","H3K4me3"),
                                           cell_types=c("neurons","microglia","oligo","astrocytes"),
                                           convert_to_GRanges=T,
                                           nThread=4,
                                           verbose=T){
  baseURL <- "https://raw.githubusercontent.com/nottalexi/brain-cell-type-peak-files/master"
  cell_dict <- list(neurons="NeuN",
                    microglia="PU1",
                    oligo="Olig2",
                    astrocytes="LHX2",
                    periph="peripheral PU1+")
  cell_dict_invert <- as.list(setNames(names(cell_dict), cell_dict))
  assay_dict <- list(ATAC="_optimal_peak_IDR_ENCODE.ATAC.bed",
                     H3K27ac="_optimal_peak.H3K27.bed",
                     H3K4me3="_optimal_peak.H3K4me3.bed")
  file_names <- unlist(lapply(assays, function(assay){file.path(assay,paste0(cell_dict[cell_types],assay_dict[assay])) }))

  printer("++ NOTT_2019:: Downloading and merging",length(file_names),"peaks BED files.", v=verbose)
  PEAKS <- parallel::mclapply(file_names, function(f, .verbose=verbose){
    printer("++ NOTT_2019:: Downloading",f, v=.verbose)
    bed_dat <- data.table::fread(file.path(baseURL,f),
                                 col.names = c("chr","start","end"),
                                 nThread = 1) #IMPORTANT! must =1 if parallelizing
    bed_dat$Assay <- dirname(f)
    bed_dat$Marker <- strsplit(basename(f),"_")[[1]][1]
    bed_dat$Cell_type <- cell_dict_invert[[strsplit(basename(f),"_")[[1]][1]]]
    return(bed_dat)
  }, mc.cores = nThread) %>% data.table::rbindlist()

  if(convert_to_GRanges){
    printer("++ NOTT_2019:: Converting merged BED files to GRanges.", v=verbose)
    PEAKS <- biovizBase::transformDfToGr(PEAKS, seqnames = "chr", start = "start", end="end")
    suppressWarnings(GenomeInfoDb::seqlevelsStyle(PEAKS) <- "NCBI")
  }
  printer("++ NOTT_2019::",length(PEAKS),"ranges retrieved.", v=verbose)
  return(PEAKS)
}




NOTT_2019.prepare_placseq_overlap <- function(merged_dat,
                                              snp_filter="!is.na(SNP)",
                                              return_counts=T){
  finemap_dat <- subset(merged_dat, eval(parse(text=snp_filter)), .drop=F)

  if(return_counts){
    interactome <- NOTT_2019.get_interactions(finemap_dat = finemap_dat)
    dat_melt <- count_and_melt(merged_annot = interactome,
                               snp_filter = snp_filter)
    if(sum(dat_melt$Count==0 | is.na(dat_melt$Count), na.rm = T)>0){
      try({
        dat_melt[dat_melt$Count==0 | is.na(dat_melt$Count),"Count"] <- NA
      })
    }
    return(dat_melt)
  } else {
    interactome <- NOTT_2019.get_interactions(finemap_dat = finemap_dat, as.granges = T)
    return(interactome)
  }
}





NOTT_2019.prepare_peak_overlap <- function(merged_dat,
                                           snp_filter="!is.na(SNP)",
                                           return_counts=T){
  PEAKS <- NOTT_2019.get_epigenomic_peaks()
  # Get SNP groups
  finemap_dat <- subset(merged_dat, eval(parse(text=snp_filter)), .drop=F)
  # Get overlap with PEAKS
  gr.hits <- GRanges_overlap(dat1 = finemap_dat,
                             chrom_col.1 = "CHR",
                             start_col.1 = "POS",
                             end_col.1 = "POS",
                             dat2 = PEAKS)

  if(return_counts){
    merged_annot <- find_topConsensus(dat = data.frame(gr.hits),
                                      grouping_vars = c("Locus","Cell_type","Assay"))
    dat_melt <- count_and_melt(merged_annot = merged_annot,
                               snp_filter = snp_filter)
    if(sum(dat_melt$Count==0 | is.na(dat_melt$Count),na.rm = T)>0){
      try({
        dat_melt[dat_melt$Count==0 | is.na(dat_melt$Count),"Count"] <- NA
      })
    }
    return(dat_melt)
  } else {return(gr.hits)}
}




NOTT_2019.prepare_regulatory_overlap <- function(merged_dat,
                                                 snp_filter="!is.na(SNP)",
                                                 return_counts=T){
  gr.reg <- NOTT_2019.get_regulatory_regions(as.granges = T)
  finemap_sub <- subset(merged_dat, eval(parse(text=snp_filter)), .drop=F)
  gr.hits.reg <- GRanges_overlap(dat1 = finemap_sub,
                                 chrom_col.1 = "CHR",
                                 start_col.1 = "POS",
                                 end_col.1 = "POS",
                                 dat2 = gr.reg)
  if(return_counts){
    merged_annot.reg <- find_topConsensus(dat = data.frame(gr.hits.reg) %>% dplyr::rename(Assay=Element),
                                          grouping_vars = c("Locus","Cell_type","Assay"))
    dat_melt.reg <- count_and_melt(merged_annot = merged_annot.reg,
                                   snp_filter = snp_filter)
  return(dat_melt.reg)
  }else {
    gr.hits.reg$Assay <- gr.hits.reg$Element
    return(gr.hits.reg)
    }
}









#' Plot brain cell-specific epigenomic data
#'
#' @keywords internal
#' @family NOTT_2019
#' @source
#' \href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
#' \url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
NOTT_2019.get_regulatory_regions <- function(as.granges=F,
                                             nThread=1,
                                             verbose=T){
  selected_sheets <- grep("promoters$|enhancers$",names(echolocatoR::NOTT_2019.interactome), value = T)
  regions <- parallel::mclapply(selected_sheets, function(s){
    printer("Importing",s,"...")
    dat <- echolocatoR::NOTT_2019.interactome[[s]]
    dat$Name <- tolower(s)
    return(dat)
  }, mc.cores = nThread) %>% data.table::rbindlist(fill=T)

  cell_dict <- c("astrocyte"="astrocytes",
                 "neuronal"="neurons",
                 "oligo"="oligo",
                 "oligodendrocytes"="oligo",
                 "microglia"="microglia")
  regions_sub <- regions %>%
    tidyr::separate(Name, into=c("Cell_type","Element"), remove=F) %>%
    dplyr::mutate(middle = as.integer( end-abs(end-start)/2),
                  Cell_type = cell_dict[Cell_type])
  if(as.granges){
    printer("++ Converting to GRanges.",v=verbose)
    regions_sub <- GenomicRanges::makeGRangesFromDataFrame(df = regions_sub,#dplyr::mutate(regions_sub, chr = as.numeric(gsub("chr","",chr))),
                                                           seqnames.field = "chr",
                                                           start.field = "start",
                                                           end.field = "end",
                                                           keep.extra.columns = T)
    suppressWarnings(GenomeInfoDb::seqlevelsStyle(regions_sub) <- "NCBI")
  }
  return(regions_sub)
}



# ***************************** #




#' Plot brain cell-specific interactome data
#'
#' @keywords internal
#' @family NOTT_2019
#' @source
#' \href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
#' \url{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf}
#' @examples
#' \dontrun{
#' data("BST1"); data("locus_dir");
#'
#' trks_plus_lines <- NOTT_2019.plac_seq_plot(finemap_dat=BST1, locus_dir=file.path("~/Desktop",locus_dir), highlight_plac=T)
#' # Zoom in
#' trks_plus_lines <- NOTT_2019.plac_seq_plot(finemap_dat=BST1, locus_dir=file.path("~/Desktop",locus_dir), zoom_window=500000, highlight_plac=T)
#' }
NOTT_2019.plac_seq_plot <- function(finemap_dat=NULL,
                                    locus_dir=NULL,
                                    title=NULL,
                                    print_plot=T,
                                    save_plot=T,
                                    return_interaction_track=F,
                                    xlims=NULL,
                                    zoom_window=NULL,
                                    index_SNP=NULL,
                                    genomic_units="Mb",
                                    color_dict=c("enhancers"="springgreen2","promoters"="purple","anchors"="black"),
                                    return_consensus_overlap=T,
                                    show_arches=T,
                                    highlight_plac=F,
                                    show_regulatory_rects=T,
                                    show_anchors=T,
                                    strip.text.y.angle=0,
                                    xtext=T,
                                    save_annot=F,
                                    point_size=2,
                                    height=7,
                                    width=7,
                                    dpi=300,
                                    as_ggplot=T,
                                    nThread=4,
                                    verbose=T){
  # finemap_dat=echolocatoR::LRRK2; print_plot=T; save_plot=T; title=NULL; index_SNP=NULL; xlims=NULL; zoom_window=NULL; return_consensus_overlap =T; nThread=1; highlight_plac=F; point_size=2;   color_dict=c("enhancers"="springgreen","promoters"="purple","anchors"="black"); genomic_units="Mb"; verbose=T; save_annot=F;
  printer("NOTT_2019:: Creating PLAC-seq interactome plot",v=verbose)
  if(!"Mb" %in% colnames(finemap_dat)){
    finemap_dat$Mb <- finemap_dat$POS/1000000
  }
  xvar <- genomic_units
  if(!"Consensus_SNP" %in% colnames(finemap_dat)){finemap_dat <- find_consensus_SNPs(finemap_dat, verbose = F)}
  marker_key <- list(PU1 = "microglia", Olig2 = "oligo",
                     NeuN = "neurons", LHX2 = "astrocytes")
  if(is.null(index_SNP)) {
    lead.pos <- subset(finemap_dat, leadSNP)[[xvar]]
  } else {
    lead.pos <- subset(finemap_dat, SNP == index_SNP)[[xvar]]
  }
  consensus.pos <- subset(finemap_dat, Consensus_SNP == T)[[xvar]]

  if(length(consensus.pos) > 0) {
    top.consensus.pos <- (dplyr::top_n(subset(finemap_dat, Consensus_SNP == T), n = 1, wt = mean.PP) %>%
                            dplyr::top_n(1, wt = Effect))[[xvar]][1]
  } else {
    top.consensus.pos <- (dplyr::top_n(subset(finemap_dat, Support > 0), n = 1, wt = mean.PP) %>%
                            dplyr::top_n(1, wt = Effect))[[xvar]][1]
  }
  if (is.null(xlims)) {
    xlims = c(min(finemap_dat[[xvar]], na.rm = T),
              max(finemap_dat[[xvar]], na.rm = T))
  }
  if (!is.null(zoom_window)) {
    xlims = c(lead.pos - as.integer(zoom_window/2), lead.pos +
                as.integer(zoom_window/2))
  }
  annot_sub <- NOTT_2019.get_promoter_interactome_data(finemap_dat = finemap_dat)
  promoter_celltypes <- NOTT_2019.get_promoter_celltypes(annot_sub = annot_sub,
                                                         marker_key = marker_key)
  #### get PLAC-seq junctions ####
  interact.DT <- NOTT_2019.get_interactome(annot_sub = annot_sub,
                                           top.consensus.pos = top.consensus.pos,
                                           marker_key = marker_key)
  #### get promoter/enhancers ####
  regions <- NOTT_2019.get_regulatory_regions(nThread = nThread,
                                              as.granges = T,
                                              verbose = verbose)
  #### regions needs to be in ####
  regions <- subset(regions,
                    as.character(GenomicRanges::seqnames(regions)) == gsub("chr","",unique(finemap_dat$CHR)[1]) &
                      GenomicRanges::start(regions) >=
                      min(finemap_dat[["POS"]], na.rm = T) & GenomicRanges::end(regions) <=
                      max(finemap_dat[["POS"]], na.rm = T))
  #### Plot PLAC-seq anchors ####
  if(show_anchors){
    interact.anchors <- NOTT_2019.get_interactions(finemap_dat = finemap_dat,
                                                   as.granges = T)
    interact.anchors$Element <- "anchors"
    regions <- c(regions, interact.anchors)
  }


  if(highlight_plac){
    # JH - which PLAC-Seq junctions overlap (5kb) the consensus SNPs?
    #interact.DT
    consensus_snps <- filter(finemap_dat, Consensus_SNP == TRUE)
    consensus_snps$CHR <- paste0("chr", consensus_snps$CHR)
    # make GRanges
    consensus.gr <- GenomicRanges::GRanges( seqnames = consensus_snps$CHR,
                                            ranges = IRanges::IRanges(start = consensus_snps$POS-1, end = consensus_snps$POS) )
    plac_start.gr <- GenomicRanges::GRanges( seqnames = interact.DT$chr,
                                             ranges = IRanges::IRanges(start = interact.DT$Start - 5000, end = interact.DT$Start) )
    plac_end.gr <- GenomicRanges::GRanges( seqnames = interact.DT$chr,
                                           ranges = IRanges::IRanges(start = interact.DT$End, end = interact.DT$End + 5000) )
    # find overlaps
    end_overlaps <- GenomicRanges::findOverlaps(consensus.gr, plac_end.gr)
    start_overlaps <- GenomicRanges::findOverlaps(consensus.gr, plac_start.gr)
    all_overlaps <- unique(c( S4Vectors::subjectHits(end_overlaps), S4Vectors::subjectHits(start_overlaps)))

    interact.DT$consensus_snp_overlap <- FALSE
    interact.DT$consensus_snp_overlap[all_overlaps] <- TRUE
  }
  #### save annotations ####
  if(save_annot) {
    annot_file <- annotation_file_name(locus_dir = locus_dir,
                                       lib_name = "Nott_2019.interactome")
    saveRDS(interact.DT, annot_file)
    annot_file <- annotation_file_name(locus_dir = locus_dir,
                                       lib_name = "Nott_2019.enhancers_promoters")
    saveRDS(regions, annot_file)
  }

  max.height=10;
  interact_y=(1.25*3) # Start after rects


  if(highlight_plac){
    NOTT.interact_trk <-
      ggbio::ggbio() +
      ggbio::geom_arch(data = interact.DT, aes(x = Start, xend = End, alpha = consensus_snp_overlap),
                       max.height = max.height, colour = "black") +
      scale_alpha_manual("Consensus SNP overlaps", values = c(0.05, 1))
  }else{
    NOTT.interact_trk <-
      ggbio::ggbio() +
      ggbio::geom_arch(data = interact.DT, alpha = 0.25, color = "black",
                       max.height = max.height,
                       aes(x = Start, xend = End, y=interact_y)) +
      labs(y=NULL)
  }

  NOTT.interact_trk <- suppressMessages(
    NOTT.interact_trk +
      facet_grid(facets = Cell_type ~ .) +
      scale_y_reverse() +
      theme_classic() +
      theme(legend.key.width = unit(1.5, "line"),
            legend.key.height = unit(1.5, "line"),
            axis.text.y = element_blank(),
            strip.text.y = element_text(angle = strip.text.y.angle)
      )
  )



  # Show enhancers/promoters as rectangles
  if (show_regulatory_rects) {
    printer("++ NOTT_2019:: Adding enhancer/promoter rectangles",v=verbose)
    rect_height <- max.height/10
    # regions$y <- ifelse(regions$Element=="promoters",  0+(rect_height*2), 0)
    regions$y <- ifelse(regions$Element=="promoters", 0+(rect_height*2),
                        ifelse(regions$Element=="enhancers", 0,
                               interact_y+rect_height))
    regions$Element <- factor(regions$Element, levels=c("enhancers","promoters","anchors"), ordered = T)

    NOTT.interact_trk <- suppressMessages(
      NOTT.interact_trk +
      ggbio::geom_rect(data = regions,
                       stat="identity",
                       rect.height= rect_height,
                       aes(y=y, fill = Element),
                       alpha = .7, inherit.aes = F,
                       color="transparent",
                       hjust=0) +
      scale_fill_manual(values = color_dict) +
      # geom_point(data = data.frame(regions),
      #            aes(x = middle, y = 0, color = Element), size = 0.5,
      #            inherit.aes = F, alpha = 0.7) +
      scale_color_manual(values = color_dict) +
      geom_hline(yintercept = Inf, alpha = 0.2,
                 show.legend = F) +
      scale_y_reverse()
    )

    if(genomic_units=="Mb"){
      printer("++ Converting genomic units to Mb...",v=verbose)
      NOTT.interact_trk <- NOTT.interact_trk +
        ggbio::scale_x_sequnit(unit = "Mb")
    }
  }
  if(xtext==F){
    NOTT.interact_trk <- NOTT.interact_trk +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank())
  }

  if (return_interaction_track) {
    printer("++ NOTT_2019:: Returning PLAC-seq track.")
    if(as_ggplot)return(NOTT.interact_trk@ggplot) else return(NOTT.interact_trk)
  } else {
    # make the ggplots
    GWAS_trk <- ggplot(data = finemap_dat, aes(x = POS, y = -log10(P), color = -log10(P))) +
      geom_point(alpha=.5, shape=16, size=point_size)

    FM_trk <- ggplot(data = finemap_dat, aes(x = POS, y = mean.PP, color = mean.PP)) +
      geom_point(alpha=.5, shape=16, size=point_size) +
      scale_color_viridis_c(breaks = c(0,0.5, 1), limits = c(0, 1)) + ylim(0, 1)

    TRACKS_list <- list(GWAS = GWAS_trk + theme_classic(),
                        `Fine-mapping` = FM_trk + theme_classic(), `Nott (2019)\nInteractome` = NOTT.interact_trk)
    params_list <- list(title = paste0(title), track.bg.color = "transparent",
                        track.plot.color = "transparent", label.text.cex = 0.7,
                        label.bg.fill = "grey12", label.text.color = "white",
                        label.text.angle = 0, label.width = unit(5.5, "lines"),
                        xlim = xlims)
    TRACKS_list <- append(TRACKS_list, params_list)
    tracks <- get("tracks", asNamespace("ggbio"))
    trks <- suppressWarnings(do.call("tracks", TRACKS_list))
    trks_plus_lines <- trks + geom_vline(xintercept = lead.pos,
                                         color = "red", alpha = 1, size = 0.3, linetype = "solid") +
      geom_vline(xintercept = consensus.pos, color = "goldenrod2",
                 alpha = 1, size = 0.3, linetype = "solid")
    if (print_plot) {print(trks_plus_lines)}

    if (save_plot) {
      plot.path <- file.path(locus_dir, paste0("Nott.sn-epigenomics_ggbio.png"))
      dir.create(dirname(plot.path), showWarnings = F,
                 recursive = T)
      ggbio::ggsave(filename = plot.path, plot = trks_plus_lines,
             height = height, width = width, dpi = dpi, bg = "transparent")
    }
    # Return
    if(as_ggplot)return(trks_plus_lines@ggplot) else return(trks_plus_lines)
  }
}





