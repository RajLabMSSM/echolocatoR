#
# ^^^^^^^^^^^^^^ Nott et al. (2019) # ^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^ single-nucleus human brain epigenomic dataset # ^^^^^^^^^^^^^^
# Nott, Alexi, Inge R. Holtman, Nicole G. Coufal, Johannes C.M. Schlachetzki, Miao Yu, Rong Hu, Claudia Z. Han, et al. “Cell Type-Specific Enhancer-Promoter Connectivity Maps in the Human Brain and Disease Risk Association.” Science 0793, no. November (2019): 778183. https://doi.org/10.1101/778183.
#  ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^
# UCSC Genome Browser: Command Line Utilities:
# http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/


NOTT_2019.epigenomic_histograms <- function(finemap_DT,
                                            results_path,
                                            show_plot=T,
                                            save_plot=T,
                                            full_data=T,
                                            return_assay_track=F,
                                            binwidth=2500, 
                                            geom="histogram",
                                            plot_formula="Assay + Cell_type ~.",
                                            bigwig_dir=NULL){
  library(BiocGenerics)
  library(GenomicRanges)
  library(ggbio)
  # GLASS DATA: UCSC GB
  # https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf
  # show_plot=T;save_plot=T;full_data=T;return_assay_track=F;binwidth=2500; geom="histogram";plot_formula="Assay + Cell_type ~.";show_regulatory_rects=T
  
  # UCSC Tracks
  import.bw.filt <- function(bw.file, gr.dat, full_data=T){
    if(full_data){
      # Get all ranges within min/max
      gr.span <- gr.dat[1,]
      mcols(gr.span) <- NULL
      start(gr.span) <- min(gr.dat$POS)
      end(gr.span) <- max(gr.dat$POS)
    } else {
      # Otherwise, just use the score for the exact values
      gr.span <- gr.dat
      }
    # bw.dat <- rtracklayer::BigWigSelection(ranges = gr.dat,  colnames = "score")
    bw.filt <- rtracklayer::import.bw(con = bw.file, selection = gr.span)
    # plot(x = start(bw.filt), y=bw.filt$score)
    return(bw.filt)
  }

  # Import BigWig annotation files
  bigWigFiles <- readxl::read_excel("./echolocatoR/annotations/Nott_2019/Nott_2019.snEpigenomics.xlsx")
  # bigWigFiles <- subset(bigWigFiles, marker!="-" &  cell_type!="peripheral microglia")
  bigWigFiles <- dplyr::mutate(bigWigFiles, cell_type = gsub(" ",".",cell_type))
  # Authors have since removed this from UCSC, but I saved it on Minerva beforehand.
  subset(bigWigFiles, cell_type!="peripheral.PU1+") 


   
  
  # Convert finemap data to granges
  dat <- finemap_DT
  dat$seqnames <- paste0("chr",dat$CHR)
  dat$start.end <- dat$POS
  gr.dat <- GenomicRanges::makeGRangesFromDataFrame(df = dat,
                                                    seqnames.field = "seqnames",
                                                    start.field = "start.end",
                                                    end.field = "start.end",
                                                    keep.extra.columns = T)
  
  bw.grlist <- parallel::mclapply(1:nrow(bigWigFiles), function(i){
    if(!is.null(bigwig_dir)){
      bw.file <- file.path(bigwig_dir,paste0(bigWigFiles$long_name[i],".ucsc.bigWig"))
    } else { bw.file <- bigWigFiles$data_link[i]}
   
    bw.name <- gsub("_pooled|pooled_","",bigWigFiles$name[i])
    printer("GVIZ:: Importing...",bw.name)
    bw.filt <- import.bw.filt(bw.file=bw.file,
                              gr.dat=gr.dat,
                              full_data=full_data)
    bw.filt$Cell_type <- bigWigFiles$cell_type[i]
    bw.filt$Assay <- bigWigFiles$assay[i]
    bw.filt$Experiment <- gsub("_"," ",bw.name)
    # colnames(mcols(bw.filt))[1] <- bw.name
    # bw.filt$expt_name <- bw.name
    # bw.filt$cell_type <-strsplit(bw.name, "_")[[1]][[1]]
    # bw.filt$assay <- strsplit(bw.name, "_")[[1]][[2]]
    return(bw.filt)
  }, mc.cores = 4)
  bw.cols <- bigWigFiles$name
  # names(bw.grlist) <- bw.cols
  bw.gr <- unlist(GenomicRanges::GRangesList(bw.grlist))

  # merge into a single granges object
  # gr.snp <- Reduce(function(x, y) GenomicRanges::merge(x, y, all.x=T),
  #                  append(bw.grlist, gr.dat))
  
  # GenomicRanges::findOverlaps(query = subset(gr.dat, Support>1), 
  #                             subject = bw.gr)

  nott_tracks <-  ggbio::autoplot(object=bw.gr,
                                  geom=geom,
                                  # facets= Experiment~.,
                                  # facets=Cell_type ~ .,
                                  # scales="free_y",
                                  # bins=300,
                                  binwidth=binwidth, 
                                  alpha=.7,
                                  position="identity",
                                  aes(fill=Cell_type), show.legend=T) +
    facet_grid(facets = formula(plot_formula)) +  
    theme(legend.position="right", strip.text.y = element_text(angle = 0))
  
  if(return_assay_track){ return(nott_tracks)}


  # Fine-mapping tracks
  TRACKS_list <- list(
    "GWAS"=ggbio::plotGrandLinear(obj = gr.dat, aes(y=-log10(P), color=-log10(P))),
    "Fine_mapping"=ggbio::plotGrandLinear(obj = gr.dat, aes(y=mean.PP, color=mean.PP)) +
      scale_color_viridis_c(),
    "Nott_etal_2019" = nott_tracks
  )

  # gr.snp[start(gr.snp)!=Inf]
  # Fuse all tracks
  params_list <- list(title = paste0(gene," [",length(seqnames(gr.dat))," SNPs]"),
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
  trks <- suppressWarnings(do.call("tracks", append(TRACKS_list, params_list)))
  # add lines
  lead.pos <- subset(finemap_DT,leadSNP)$POS
  consensus.pos <- subset(finemap_DT, Consensus_SNP==T)$POS

  trks_plus_lines <- trks +
    # theme_bw() +
    # ggbio::theme_genome() +
    theme(strip.text.y = element_text(angle = 0),
          strip.text = element_text(size = 9),
          panel.background = element_rect(fill = "white", colour = "black", linetype = "solid")) +
    geom_vline(xintercept = consensus.pos, color="goldenrod2", alpha=1, size=.1, linetype='solid') +
    geom_vline(xintercept = lead.pos, color="red", alpha=1, size=.1, linetype='solid') +
    scale_x_continuous( labels=function(x)x/1000000)

  if(show_plot){ print(trks_plus_lines) }
  if(save_plot){
    ggsave(file.path(results_path,"Annotation",
                     paste0(basename(results_path),"_Glass.snEpigenomics.png") ),
           plot = trks_plus_lines, dpi=400, height = 15, width = 8,
           bg = "transparent")
  }
  return(trks_plus_lines)
}


NOTT_2019.superenhancers <- function(s6_path="./echolocatoR/annotations/Nott_2019/aay0793-Nott-Table-S6.xlsx"){
  s6 <- readxl::read_excel( , skip = 2)
  annot_sub <- subset(s6, chr== paste0("chr",unique(finemap_DT$CHR)) & start>=min(finemap_DT$POS) & end<=max(finemap_DT$POS) )
  if(nrow(annot_sub)>0){
    merged_DT <- data.table:::merge.data.table(finemap_DT %>%
                                                 dplyr::mutate(chr=paste0("chr",CHR),
                                                               start=as.numeric(POS)) %>%
                                                 data.table::data.table(),
                                               data.table::data.table(s6),
                                               by = c("chr","start"))
  }
}

NOTT_2019.get_promoter_interactome_data <- function(finemap_DT,
                                                     s5_path="./echolocatoR/annotations/Nott_2019/aay0793-Nott-Table-S5.xlsx"){
  sheets_s5 <- readxl::excel_sheets(s5_path)
  s5 <- readxl::read_excel(s5_path, sheet = sheets_s5[1], skip = 2)
  s5 <- s5 %>% dplyr::rename(chr=Chr, start=Start, end=End)
  # Subset to window
  annot_sub <- subset(s5, chr== paste0("chr",unique(finemap_DT$CHR)) &
                        start>=min(finemap_DT$POS) &
                        end<=max(finemap_DT$POS) )
  return(annot_sub)
}

NOTT_2019.get_promoter_celltypes <- function(annot_sub, marker_key){
  promoter.cols <- grep("*_active_promoter",  colnames(annot_sub), value = T)
  logical.list <- colSums(annot_sub[,promoter.cols])>0
  promoter_celltypes <- gsub("\\_.*","", promoter.cols[as.logical(logical.list)] )
  promoter_celltypes <- as.character(marker_key[promoter_celltypes])
  promoter_celltypes <- paste(promoter_celltypes,collapse="; ")
  return(promoter_celltypes)
}

NOTT_2019.get_interactome <- function(annot_sub, top.consensus.pos, marker_key){
  interact.cols <- grep("*_interactions", colnames(annot_sub), value = T)
  interact.DT <- lapply(interact.cols, function(column){
    coords <- strsplit(annot_sub[,column][[1]], ",")
    coord.dt <- lapply(coords, function(coord){
      data.table::data.table(Interaction=column,
                             Cell_type=marker_key[gsub("\\_.*","",column)],
                             Coordinates=coord) %>% return()
    }) %>% data.table::rbindlist()
    return(coord.dt)
  } )  %>% data.table::rbindlist()
  interact.DT <- subset(interact.DT, !is.na(Coordinates) & Coordinates!="") %>%
    tidyr::separate(col = Coordinates,
                    into=c("chr","Start","End"), sep = ":|-") %>%
    tidyr::separate(col = Interaction, into=c("Marker","Element",NA), sep="_", remove = F )
  interact.DT <- interact.DT %>%
    dplyr::mutate(Cell_type_interaction=paste(Cell_type,"-",Element))
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


NOTT_2019.get_interactions <- function(finemap_DT,
                                       s5_path="./echolocatoR/annotations/Nott_2019/aay0793-Nott-Table-S5.xlsx"){
  sheets_s5 <- readxl::excel_sheets(s5_path)
  selected_sheets <- grep("interactome$",sheets_s5, value = T)
  interactomes <- lapply(selected_sheets, function(s){
    printer("Importing",s,"...")
    # Read the sheet you want
    dat <- readxl::read_excel(s5_path, sheet = s, skip = 2)
    dat$Name <- s
    return(dat)
  }) %>% data.table::rbindlist()
  interactomes_sub <- subset(interactomes,
                         chr1== paste0("chr",unique(finemap_DT$CHR)) &
                         chr2== paste0("chr",unique(finemap_DT$CHR)) &
                         start1>=min(finemap_DT$POS) &
                         start2>=min(finemap_DT$POS) &
                         end1<=max(finemap_DT$POS) &
                         end2<=max(finemap_DT$POS)
                         ) %>% tidyr::separate(Name, into=c("Cell_type","Element"), remove=F)
  return(interactomes_sub)
}

NOTT_2019.get_epigenomic_peaks <- function(peak.dir="/Volumes/Scizor/Nott_2019/peaks",
                                           narrow_peaks=T,
                                           broad_peaks=T){
  bigWigFiles <- readxl::read_excel("./echolocatoR/annotations/Nott_2019/Nott_2019.snEpigenomics.xlsx")
  peak_types <- c(ifelse(narrow_peaks,".narrowPeak$", NA),
                  ifelse(broad_peaks,"_broad.bed12$", NA))
  peak_types <- peak_types[!is.na(peak_types)]

  peaks.paths <- list.files(peak.dir,
                            pattern =  paste(peak_types, collapse = "|"),
                            full.names = T,recursive = T)
  PEAKS <- MACS2.import_peaks(peaks.paths, as_granges=T)
  PEAKS <- lapply(PEAKS, function(peak){
    pk.name <- gsub(".ucsc_narrowPeak1|.ucsc_broadRegion1","",peak$name[1])
    meta <- subset(bigWigFiles, long_name==pk.name)
    peak$Cell_type <- meta$cell_type
    peak$Assay <- meta$assay
    peak$Fresh_frozen <- meta$fresh_frozen
    peak$marker  <- meta$marker
    return(peak)
  }) %>% GenomicRanges::GenomicRangesList()
 PEAKS.merged <- unlist(PEAKS)
 PEAKS.merged$peak_type <- PEAKS.merged
 return(PEAKS.merged)
}


NOTT_2019.get_regulatory_regions <- function(finemap_DT,
                                             s5_path="./echolocatoR/annotations/Nott_2019/aay0793-Nott-Table-S5.xlsx",
                                             as.granges=F){
  # Get sheet names
  sheets_s5 <- readxl::excel_sheets(s5_path)
  selected_sheets <- grep("promoters$|enhancers$",sheets_s5, value = T)
  regions <- parallel::mclapply(selected_sheets, function(s){
    printer("Importing",s,"...")
    # Read the sheet you want
    dat <- readxl::read_excel(s5_path, sheet = s, skip = 1, col_names = c("chr","start","end"))
    dat$Name <- s
    return(dat)
  }, mc.cores = 4) %>% data.table::rbindlist()
  regions_sub <- regions %>%
    # subset(chr== paste0("chr",unique(finemap_DT$CHR)) &
    #                     start>=min(finemap_DT$POS) &
    #                     end<=max(finemap_DT$POS) ) %>%
    tidyr::separate(Name, into=c("Cell_type","Element"), remove=F) %>% 
    dplyr::mutate(middle=as.integer( end-abs(end-start)/2) )
  if(as.granges){
    regions_sub <- GenomicRanges::makeGRangesFromDataFrame(df = regions_sub,#dplyr::mutate(regions_sub, chr = as.numeric(gsub("chr","",chr))),
                                                     seqnames.field = "chr",
                                                     start.field = "start",
                                                     end.field = "end",
                                                     keep.extra.columns = T)
    
  } 
  return(regions_sub)
}

NOTT_2019.report_regulatory_overlap <- function(finemap_DT, regions){
  library(GenomicRanges)
  library(BiocGenerics) 
  # consensus.snps <- subset(finemap_DT, Consensus_SNP==T) 
  gr.finemap <- GenomicRanges::makeGRangesFromDataFrame(finemap_DT, 
                                                          seqnames.field = "CHR", 
                                                          start.field = "POS", end.field = "POS",
                                                          ignore.strand = T,
                                                          keep.extra.columns = T)
  
  if(class(regions)[1]=="GRanges"){
    gr.regions <- regions 
    GenomeInfoDb::seqlevelsStyle(gr.regions) <- "NCBI"
  } else{
    gr.regions <- GenomicRanges::makeGRangesFromDataFrame(regions, 
                                                          seqnames.field = "CHR", 
                                                          start.field = "start", end.field = "end",
                                                          ignore.strand = T,
                                                          keep.extra.columns = T) 
  }
  
  hits <- GenomicRanges::findOverlaps(query = gr.finemap,
                                      subject = gr.regions) 
  
  gr.hits <- gr.regions[ S4Vectors::subjectHits(hits), ] 
  mcols(gr.hits) <- cbind(mcols(gr.hits),
                          mcols(gr.finemap[S4Vectors::queryHits(hits),]) )
  # gr.hits <- cbind(mcols(gr.regions[ S4Vectors::subjectHits(hits), ] ),
  #                         mcols(gr.consensus[S4Vectors::queryHits(hits),]) )
  message("",nrow(mcols(gr.hits))," query SNP(s) detected with reference overlap." )
  # print(data.frame(mcols(gr.hits[,c("Name","SNP")])) )
  return(gr.hits)
}

# ***************************** #
NOTT_2019.plac_seq_plot <- function(finemap_DT=NULL,
                                     title=NULL,
                                     print_plot=T,
                                     save_plot=T,
                                     return_interaction_track=F,
                                     xlims=NULL,
                                     zoom_window=NULL,
                                     index_SNP=NULL,
                                     return_consensus_overlap=T,
                                     show_arches=T, 
                                     show_regulatory_rects=T){
  #  quick_finemap()
  # print_plot=T; save_plot=T; title=NULL; index_SNP=NULL; xlims=NULL; zoom_window=NULL; return_consensus_overlap =T
  
  library(ggbio)
  marker_key <- list(PU1="Microglia",
                     Olig2="Oligo",
                     NeuN="Neuronal",
                     LHX2="Astrocyte")
  if(is.null(index_SNP)){
    lead.pos <- subset(finemap_DT, leadSNP)$POS #top_n(finemap_DT, n = 1, wt = -P)$POS
  } else {
    lead.pos <- subset(finemap_DT, SNP==index_SNP)$POS
  }

  #subset(finemap_DT,SNP=="rs76904798")$POS
  consensus.pos <- subset(finemap_DT, Consensus_SNP==T)$POS
  if(length(consensus.pos)>0){
    top.consensus.pos <- (top_n(subset(finemap_DT, Consensus_SNP==T),
                                n=1, wt = mean.PP) %>% top_n(1,wt=Effect))$POS[1]
  } else {
    top.consensus.pos <- (top_n(subset(finemap_DT, Support>0),
                                n=1, wt = mean.PP )%>% top_n(1,wt=Effect))$POS[1]
  }
  # Define window size
  if(is.null(xlims)){
    xlims = c(min(finemap_DT$POS), max(finemap_DT$POS))
  }
  if(!is.null(zoom_window)){
    xlims = c(lead.pos-as.integer(zoom_window/2),  lead.pos+as.integer(zoom_window/2))
  }


  # Subset interactome to relevant region
  annot_sub <- NOTT_2019.get_promoter_interactome_data(finemap_DT = finemap_DT)
  ## Extract active promoter celltypes
  promoter_celltypes <- NOTT_2019.get_promoter_celltypes(annot_sub = annot_sub,
                                                         marker_key = marker_key)
  ## Extract promoter interactions
  interact.DT <- NOTT_2019.get_interactome(annot_sub = annot_sub,
                                            top.consensus.pos =  top.consensus.pos,
                                            marker_key = marker_key)
  # interactions <- NOTT_2019.get_interactions(finemap_DT = finemap_DT)
  
  
  
  ## $$$$$ HERE $$$$$$ : use for fig 1, third col
  regions <- NOTT_2019.get_regulatory_regions(finemap_DT = finemap_DT, 
                                              as.granges = T)  
  regions <- subset(regions, seqnames(regions) %in% paste0("chr",unique(finemap_DT$CHR)) & 
                   start(regions)>=min(finemap_DT$POS) &
                   end(regions)<=max(finemap_DT$POS))
  # Report overlap 
  printer("+ NOTT_2019:: Evaluating Consensus SNP overlap with regulatory elements...") 
  gr.hits <- NOTT_2019.report_regulatory_overlap(finemap_DT=subset(finemap_DT, Consensus_SNP==T),
                                                regions=regions)
  
 

  ####### Assemble Tracks #######
  # GWAS track
  GWAS_trk <- ggplot(data=finemap_DT, aes(x=POS, y=-log10(P), color=-log10(P))) +
    geom_point()

  FM_trk <- ggplot(data=finemap_DT, aes(x=POS, y=mean.PP, color=mean.PP)) +
    geom_point() +
    scale_color_viridis_c(breaks=c(0,.5,1), limits=c(0,1)) +
    ylim(0,1)

  # Nott tracks
  # Nott_s5 <- ggbio::ggbio() +
  #   ggbio::geom_rect(data = annot_sub, aes(xmin=start, xmax=end, ymin=0, ymax=1),
  #                    fill="turquoise", alpha=.75) +
  #   facet_grid(facets = Annotation~.)
  # Nott_s5 <- invisible_legend(Nott_s5)

  # Nott:  interactions
  # Nott_interactions
  NOTT.interact_trk <- ggbio::ggbio() + 
    ggbio::geom_arch(data = interact.DT, alpha=.6, color="gray60", max.height = 10,
                                        aes(x=Start, xend=End)) + 
    facet_grid(facets = Cell_type~.) +
    # scale_y_reverse() +
    # scale_colour_brewer(palette = "Set2") +
    theme_classic() +
    # labs(subtitle = paste0(annot_sub$Annotation[[1]]," - ",promoter_celltypes) ) +
    theme(legend.key.width=unit(1.5,"line"),
          legend.key.height=unit(1.5,"line"),
          axis.text.y = element_blank() )  
    # xlim(c(40500000, 40700000))
  
  if(show_regulatory_rects){
    NOTT.interact_trk <- NOTT.interact_trk + 
      ggbio::geom_rect(data = regions,
                       aes(xmin=start, xmax=end, ymin=-1, ymax=1, fill=Element), alpha=.8, inherit.aes=F) +
      scale_fill_manual(values = c("turquoise2", "purple2")) +
      geom_point(data = data.frame(regions), aes(x=middle, y=0, color=Element), size=.5,
                 inherit.aes = F,  alpha=.8) +
      scale_color_manual(values = c("turquoise","purple")) +
      geom_hline(yintercept = Inf, alpha=.2, show.legend = F)
  }
  

  if(return_interaction_track){
    printer("++ Nott sn-epigenomics:: Returning PLAC-seq track.")
    return(NOTT.interact_trk
             # ggbio::geom_rect(data = annot_sub,
             #                  aes(xmin=start, xmax=end, ymin=-0, ymax=Inf),
             #                  fill="turquoise", alpha=.5, inherit.aes=F) +
            )
  } else{
    # Makes tracks list
    TRACKS_list <- list(
      "GWAS"=GWAS_trk + theme_classic(),
      "Fine-mapping"=FM_trk + theme_classic(),
      "Nott (2019)\nInteractome"=NOTT.interact_trk
      # "Nott et al. (2019)\nPromoter\ninteractome"=Nott_s5
    )
    # Parameters
    params_list <- list(title = paste0(title),
                        track.bg.color = "transparent",
                        track.plot.color = "transparent",
                        label.text.cex = .7,
                        label.bg.fill = "grey12",
                        label.text.color = "white",
                        label.text.angle = 0,
                        label.width = unit(5.5, "lines"),
                        xlim = xlims
                        # xlim = c(min(finemap_DT$POS), max(finemap_DT$POS))
                        # heights = c(rep(1,length(INIT_list)), rep(1,length(BW_tracks)) )
    )
    TRACKS_list <- append(TRACKS_list, params_list)
    trks <- suppressWarnings(do.call("tracks", TRACKS_list))

    trks_plus_lines <- trks +
      # Nott: promoter
      # ggbio::geom_rect(data = annot_sub, aes(xmin=start, xmax=end, ymin=-0, ymax=Inf),
      #                  fill="turquoise", alpha=.5, inherit.aes=F) +
      # Lead GWAS line
      geom_vline(xintercept = lead.pos, color="red", alpha=1, size=.3, linetype='solid') +
      # Consensus line
      geom_vline(xintercept = consensus.pos, color="goldenrod2", alpha=1, size=.3, linetype='solid')

      # theme_bw() +
      # theme(plot.subtitle = element_text(color = "turquoise", size = 8)) +
      # ggbio::geom_rect(data = regions,
      #                  aes(xmin=start, xmax=end, ymin=-1, ymax=1, fill=Element), alpha=.8, inherit.aes=F) +
      # scale_fill_manual(values = c("turquoise2", "purple2"))

    if(print_plot){print(trks_plus_lines)}
    # SAVE PLOT
    if(save_plot){
      plot.path <- file.path(results_path,"Multi-finemap",
                             paste0("Nott.sn-epigenomics_ggbio.png"))
      dir.create(dirname(plot.path), showWarnings = F, recursive = T)
      ggsave(filename = plot.path,
             plot = trks_plus_lines,
             height = 7, width = 7, dpi = 200, bg = "transparent")
    }
    return(trks_plus_lines)
  }
}




