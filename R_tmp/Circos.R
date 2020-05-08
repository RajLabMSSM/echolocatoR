
# --------------- Circos Plot Packages --------------- #

## 1. ggbio: http://girke.bioinformatics.ucr.edu/CSHL_RNAseq/mydoc/mydoc_Rgraphics_7/
# http://bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
###  - Limited number of plot types. Not very flexible.

## 2. BioCircos: https://cran.r-project.org/web/packages/BioCircos/vignettes/BioCircos.html
### - Specifying point size makes all points disappear.
### - Annotations via BioCircosArcTrack() function don't appear.
### - Benefit - Interactive, can add metadata to points.


## 3. OmicCircos: https://bioconductor.org/packages/release/bioc/vignettes/OmicCircos/inst/doc/OmicCircos_vignette.pdf
### - Intstructions opaque to the point of being useless. Do not adequately desribe input data format.

# ----------------------------------------------------#



# ======  BioCircos ======== #

BioCircos.add_SNP_track <- function(tracklist=list(),
                                    gr, 
                                    trackname="GWAS",
                                    minRadius=0.5,
                                    maxRadius=0.9,
                                    point.colors=c("limegreen"),
                                    background=T, 
                                    height.range=c(0,2),
                                    max.size=34){
  # Group By Locus and Support
  gr[[trackname]] <- gr[[trackname]] %>% 
    dplyr::group_by(Gene, Support) %>%  
    dplyr::select(Gene, Support, CHR, POS, Effect, SNP) %>% unique() %>%
    dplyr::summarise(CHR=mean(CHR),
                     POS=mean(POS),
                     Size=n(), 
                     Effect=mean(Effect), 
                     SNPs = paste0(unique(SNP), collapse = "<br>"))
  gr[[trackname]]$Label <- paste0("<strong>Locus</strong>: ",gr[[trackname]]$Gene,
                                  "<br><strong>Support Level</strong>: ",gr[[trackname]]$Support, 
                                  "<br><strong>n SNPs</strong>: ",gr[[trackname]]$Size,
                                  "<br><strong>rsids</strong>: <br>", gr[[trackname]]$SNPs)
  # gr[[trackname]]$Size <- scales::rescale(gr[[trackname]]$Size, to=c(0,1), from=c(0,max.size))
  
  printer("BioCircos:: Constructing SNP track for", trackname)
  printer("+ BioCircos:: n SNPs =", length(gr[[trackname]]$POS ))
  # printer("+ BioCircos:: n Sizes =", length(gr[[trackname]]$Size ))
  
  tracklist <- tracklist + BioCircos::BioCircosSNPTrack(trackname = trackname, 
                                                        chromosomes = gr[[trackname]]$CHR, 
                                                        # points_coordinates = gr[[trackname]]$POS, 
                                                        positions = gr[[trackname]]$POS, 
                                                        values = gr[[trackname]]$Support, 
                                                        labels = gr[[trackname]]$Label,
                                                        # size = gr[[trackname]]$Size,
                                                        colors = point.colors,
                                                        minRadius = minRadius,
                                                        maxRadius = maxRadius, 
                                                        opacities = .5,
                                                        range = height.range
  )
  # Background are always placed below other tracks
  if(background){ 
    tracklist <- tracklist + BioCircos::BioCircosBackgroundTrack("GWAS.background", 
                                                                 minRadius = minRadius,
                                                                 maxRadius = maxRadius,
                                                                 borderColors = "#AAAAAA", 
                                                                 borderSize = 0.6, 
                                                                 fillColors = "black")  
  } 
  return(tracklist)
}

BioCircos.SNP_Groups <- function(FM = merge_finemapping_results()
){
  library(BioCircos)
  tracklist =  BioCircos::BioCircosTracklist()
  increment <- .2
  minRadius <- .2
  maxRadius <- minRadius + increment
  color.list <- list(Consensus="goldenrod3",
                     CS="limegreen",
                     GWAS="red",
                     all="lightblue")
  
  gr <- list("Consensus"=subset(FM, Consensus_SNP==T),
             "CS"=subset(FM, Support>0),
             "GWAS"=subset(FM, P<0.5),
             "all"=FM
  )
  max.size <- max((FM %>% dplyr::group_by(Gene) %>% count())$n)
  height.range <- c(length(grep(".Credible_Set",colnames(FM))), 0)
  # Iterate tracks 
  for(group in c("Consensus","CS","GWAS")){
    minRadius <- maxRadius + .05
    maxRadius <- minRadius + increment
    tracklist <- tracklist + BioCircos.add_SNP_track(tracklist, 
                                                     gr = gr,
                                                     trackname = group, 
                                                     minRadius=minRadius,
                                                     maxRadius=maxRadius,
                                                     point.colors=color.list[[group]], 
                                                     height.range=height.range)
  }
  # Bar plots 
  # gr.list <- ggbio.prepare_SNPgroups(FM = FM)
  # # xgr.list <- XGR.gather_XGR_annotations(gr.all = gr.list$all, merge_all = F, overlap.cutoff = 1)
  # xgr <- xgr.list$ENCODE_TFBS_ClusteredV3_CellTypes 
  # chroms <- as.numeric( gsub("chr","",as.character(seqnames(xgr))))
  # tracklist = tracklist +  <-BioCircosBarTrack("XGR",
  #                         chromosomes = chroms,
  #                         starts = GenomicAlignments::start(xgr),
  #                         ends = GenomicAlignments::end(xgr), 
  #                         minRadius = 5,
  #                         maxRadius = 5.5
  #                     ) 
  # # Plot all tracks 
  # tracklist <- tracklist + BioCircos::BioCircosBackgroundTrack(trackname ="chrom" , 
  #                                                              fillColors = "BuPu",
  #                                                              minRadius=maxRadius+.5,
  #                                                              maxRadius=maxRadius+increment)
  BC <- BioCircos(tracklist, 
                  genomeFillColor = "BuPu",#"Spectral",#"PuOr",
                  chrPad = 0.005,
                  displayGenomeBorder = T,
                  yChr =  F,
                  genomeTicksDisplay = F,
                  genomeLabelTextSize = 12,
                  genomeLabelDy = 0,
                  minRadius=maxRadius+.5,
                  maxRadius=maxRadius+increment,
                  SNPMouseOverColor = "blue" 
  )  
  return(BC)
} 




# ======  ggbio ======== #
ggbio.prepare_SNPgroups <- function(FM = merge_finemapping_results(dataset = "./Data/GWAS/Nalls23andMe_2019")){ 
  DF.dat <- FM 
  DF.dat$Consensus <- FM$Consensus_SNP>0
  DF.dat$CS <- FM$Support>0
  DF.dat$GWAS <-FM$P<=0.05
  na.rm <- F
  
  gr.all <- DF.dat %>%  
    dplyr::mutate(SeqNames=paste0("chr",CHR)) %>%  data.frame() %>%
    biovizBase::transformDfToGr(seqnames = "SeqNames", start = "POS", end = "POS") 
  
  # Average by locus
  gr.loci <- DF.dat %>%  dplyr::group_by(Gene) %>% 
    dplyr::summarise(CHR=mean(CHR), 
                     mean.POS=mean(POS, na.rm = na.rm), 
                     min.POS=min(POS), 
                     max.POS=max(POS), 
                     Effect=mean(Effect, na.rm = na.rm), 
                     Support=mean(Support, na.rm = na.rm)) %>%
    dplyr::mutate(SeqNames=paste0("chr",CHR)) %>% data.frame() %>%
    data.frame() %>%
    biovizBase::transformDfToGr(seqnames = "SeqNames", start = "min.POS", end = "max.POS")
  gr.loci$Height <- rep(c(1:5), length(gr.loci$SeqNames))[1:length(gr.loci$SeqNames)]
  
  
  # Average by Chromosome
  gr.chr <- DF.dat %>% 
    dplyr::mutate(SeqNames=paste0("chr",CHR)) %>% 
    dplyr::group_by(SeqNames) %>% 
    dplyr::summarise(mean.POS=mean(POS, na.rm = na.rm)) %>%   
    data.frame() %>%
    biovizBase::transformDfToGr(seqnames = "SeqNames", start = "mean.POS", end = "mean.POS")  
  
  # GWAS
  gr.GWAS <- subset(FM, P<=0.05) %>% 
    dplyr::group_by(Gene) %>% 
    dplyr::summarise(CHR=mean(CHR, na.rm = na.rm), 
                     POS=mean(POS, na.rm = na.rm), 
                     Size=n(), 
                     Effect=mean(Effect, na.rm = na.rm), 
                     Support=mean(Support, na.rm = na.rm)) %>%
    dplyr::mutate(SeqNames=paste0("chr",CHR)) %>%  data.frame() %>%
    biovizBase::transformDfToGr(seqnames = "SeqNames", start = "POS", end = "POS")
  gr.GWAS$SNP.Group <- "GWAS"
  if(any(is.na(gr.GWAS$Support))){
    gr.GWAS[is.na(gr.GWAS$Support),]$Support <- 0 
  }
  
  
  #CredSet
  gr.CS <- subset(FM, Support>0) %>% 
    dplyr::group_by(Gene) %>% 
    dplyr::summarise(CHR=mean(CHR, na.rm = na.rm), 
                     POS=mean(POS, na.rm = na.rm),
                     Size=n(), 
                     Effect=mean(Effect, na.rm = na.rm), 
                     Support=mean(Support, na.rm = na.rm)) %>%
    dplyr::mutate(SeqNames=paste0("chr",CHR)) %>%  data.frame() %>%
    biovizBase::transformDfToGr(seqnames = "SeqNames", start = "POS", end = "POS")
  gr.CS$SNP.Group <- "Credible Set"
  
  # Consensus
  gr.Consensus <- subset(FM, Consensus_SNP==T) %>% 
    dplyr::group_by(Gene) %>% 
    dplyr::summarise(CHR=mean(CHR, na.rm = na.rm), 
                     POS=mean(POS, na.rm = na.rm), 
                     Size=n(), 
                     Effect=mean(Effect, na.rm = na.rm),
                     Support=mean(Support, na.rm = na.rm)) %>%
    dplyr::mutate(SeqNames=paste0("chr",CHR)) %>%  data.frame() %>%
    biovizBase::transformDfToGr(seqnames = "SeqNames", start = "POS", end = "POS")
  gr.Consensus$SNP.Group <- "Consensus" 
  
  # SNP-level
  gr.Cons.rsid <- subset(FM, Consensus_SNP==T) %>%
    # dplyr::group_by(Gene) %>% slice(1) %>%
    dplyr::mutate(SeqNames=paste0("chr",CHR)) %>% data.frame() %>%
    biovizBase::transformDfToGr(seqnames = "SeqNames", start = "POS", end = "POS") 
  gr.Cons.rsid$Height <- rep(c(1:5), length(gr.Cons.rsid$SeqNames))[1:length(gr.Cons.rsid$SeqNames)]
  return(list(
    "all"=gr.all,
    "loci"=gr.loci,
    "chr"=gr.chr,
    "GWAS"=gr.GWAS,
    "CS"=gr.CS,
    "Consensus"=gr.Consensus,
    "Cons.rsid"=gr.Cons.rsid)) 
}
# gr <- ggbio.prepare_SNPgroups(FM)



ggbio.add_fgwas_annnotations <- function(ggb, top_annot=10){
  printer("ggbio:: Creating annotation tracks using fGWAS processed files.")
  fgwas_annots <- data.table::fread("./Data/GWAS/Nalls23andMe_2019/_genome_wide/fGWAS/Input/annotations.Nalls23andMe_2019.txt") 
  fgwas_results <- data.table::fread("./Data/GWAS/Nalls23andMe_2019/_genome_wide/fGWAS/fGWAS_summary.Nalls23andMe_2019.txt")
  top_annot.names <- (subset(fgwas_results, SNP.Group=="Consensus") %>% 
                        dplyr::rename(AIC="AIC:") %>% 
                        arrange(desc(estimate), desc(AIC)))[1:top_annot,]$parameter %>%  gsub("_ln","", x=.)
  
  # col.sums <- colSums(fgwas_annots[,4:ncol(fgwas_annots)]) 
  annot.cols <- top_annot.names#colnames(fgwas_annots)[4:ncol(fgwas_annots)]
  color.list <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  for(i in 1:length(annot.cols)){ 
    printer(annot.cols[i])
    annot.sub <- subset(fgwas_annots, select = c("chr","pos",annot.cols[i])) %>%  
      dplyr::mutate("SeqNames" = chr) 
    
    # Fit line to data 
    lm.out <- loess(data = annot.sub, 
                    span = 0.25,
                    formula = "`HCPEpiC-DS12447.hotspot.twopass.fdr0.05.merge` ~ pos")
    
    
    annot.sub <- annot.sub[annot.sub[annot.cols[i]]==1,]
    annot.gr <- biovizBase::transformDfToGr(annot.sub, seqnames = "SeqNames", 
                                            start = "pos", end = "pos") 
    gg.new <- ggbio::circle(annot.gr, 
                            geom="rect", 
                            label = annot.cols[i],
                            color = color.list[i]) 
    ggb <- ggb + gg.new
    # ggbio::ggbio() +  gg.new
  } 
  return(ggb)
}






ggbio.circos <- function(gr){ 
  library(ggbio)
  gr <- ggbio.prepare_SNPgroups()
  xgr <- XGR.gather_XGR_annotations(gr.all = gr[["all"]], merge_all=F, overlap.cutoff = 1)
  
  # SNP grid
  alpha <- .75
  grid.background <-  "black"
  grid.line <-"darkgrey"
  grid.n <- 3
  # Annotation grid 
  grid.background.ann <- "black"
  grid.line.ann <-"black"
  grid.n.ann <- 3  
  # Karyogram
  # data(ideoCyto, package = "biovizBase")
  # biovizBase::isIdeogram(ideoCyto$hg19)
  # ideoCyto$hg19$SeqNames <- as.character(seqnames(ideoCyto$hg19))
  # data("CRC", package  = "biovizBase")
  # head(hg19sub) 
  # autoplot(ideoCyto$hg19, layout = "circle", cytobands = TRUE)
  #  
  
  
  ggb <- ggbio::ggbio() +  
    # SNP Groups
    ggbio::circle(gr[["Consensus"]], aes(fill=SNP.Group, size=Size, y=Support),
                  geom = "point", color="goldenrod2",colors="goldenrod2", alpha=alpha,
                  grid=T, grid.background=grid.background, grid.n=grid.n, grid.line=grid.line) + #, radius=20
    ggbio::circle(gr[["CS"]], aes(fill=SNP.Group, size=Size, y=Support),
                  geom = "point", color="limegreen", alpha=alpha,
                  grid=T, grid.background=grid.background, grid.n=grid.n, grid.line=grid.line) + # , radius=15 
    ggbio::circle(gr[["GWAS"]], aes(fill=SNP.Group, size=Size, y=Support),
                  geom = "point", color="red", alpha=alpha,
                  grid=T, grid.background=grid.background, grid.n=grid.n, grid.line=grid.line) #, radius=10
  
  # ggb <- ggbio.add_fgwas_annnotations(ggb)
  
  # ggbio::circle(GRangesList(gr[["GWAS"]], gr[["Consensus"]], gr[["CS"]]) %>% unlist(),
  #               aes(color=SNP.Group, size=Size, y=Support, fill=SNP.Group, group=SNP.Group),
  #               geom = "point", colors="goldenrod3", alpha=alpha,
  #               grid=T, grid.background=grid.color, grid.n=1) +  
  # RSID text
  # ggbio::circle(gr[["Cons.rsid"]], geom = "text", 
  # aes(label=SNP, size=5, y=Height), trackWidth=20) +  # angle=0
  
  # Annotations
  # ggbio::circle(xgr[["ENCODE_TFBS_ClusteredV3_CellTypes"]], geom="rect", color="purple",
  #               grid=T, grid.background=grid.background.ann, grid.n=grid.n.ann, grid.line=grid.line.ann) +
  # ggbio::circle(xgr[["ENCODE_DNaseI_ClusteredV3_CellTypes"]], geom="rect", color="magenta",
  #               grid=T, grid.background=grid.background.ann, grid.n=grid.n.ann, grid.line=grid.line.ann) +
  # ggbio::circle(xgr[["Broad_Histone"]], geom="rect", color="pink",
  #               grid=T, grid.background=grid.background.ann, grid.n=grid.n.ann, grid.line=grid.line.ann) +
  # ggbio::circle(xgr, geom="rect", aes(group=Dataset, fill=Dataset, color=Dataset)) +
  
  
  
  # Loci
  # circle(gr.loci, geom = "bar", color="black", aes(y=Height, label=Gene)) +
  ggb <- ggb +
    ggbio::circle(gr[["loci"]], geom = "ideogram", fill="turquoise3", color="turquoise4") +  
    ggbio::circle(gr[["loci"]], geom = "bar", color="turquoise3", trackWidth=10, 
                  aes(label=Gene, y=Height)) +
    ggbio::circle(gr[["loci"]], geom = "text", color="turquoise3", trackWidth=20, 
                  aes(label=Gene, size=5, y=Height)) + 
    
    # Chromosomes
    ggbio::circle(hg19sub, geom = "ideogram", cytobands=T, trackWidth=22) + 
    ggbio::circle(hg19sub, geom = "text", color="black", trackWidth=22,
                  aes(label=seqnames, size=8))
  print(ggb)  
  png(file.path("./Data/GWAS/Nalls23andMe_2019/_genome_wide/ggbio_circos.png"), 
      bg = "transparent", height = 1000, width=1000)
  ggb
  dev.off() 
  
  return(ggb)
}

# unit. = 360/length(unique(gr.all$CHR))
# lapply(1:length(unique(gr.all$CHR)), function(x){unit.+i*})
# setNames(unique(gr.all$CHR), rep(1,length(unique(gr.all$CHR))))

ggbio.mergeQTL <- function(FM, qtl_file="Fairfax_2014_IFN"){
  FM_merge <- mergeQTL.merge_handler(FM_all=FM, qtl_file) 
  FM_merge[is.na(FM_merge$QTL.FDR), "QTL.FDR"] <- 1 
  FM_merge[is.infinite(FM_merge$QTL.FDR), "QTL.FDR"] <- 1 
  FM_merge <- FM_merge %>% dplyr::group_by(SNP) %>% dplyr::slice(1) %>% data.table::data.table()
  # Convert to 
  gr.qtl <- FM_merge %>% dplyr::mutate(SeqNames=paste0("chr",CHR),
                                       QTL.neg.log= -log10(QTL.FDR)) %>%  
    data.frame() %>%
    biovizBase::transformDfToGr(seqnames = "SeqNames", start = "POS", end = "POS") 
  # gr.qtl <- subset(gr.qtl, QTL.FDR < 0.05) 
  if(any(is.na(gr.qtl$QTL.neg.log))){ gr.qtl[is.na(gr.qtl$QTL.neg.log),]$QTL.neg.log <- 0 }
  if(any(is.infinite(gr.qtl$QTL.neg.log))){ gr.qtl[is.infinite(gr.qtl$QTL.neg.log),]$QTL.neg.log <- 0 } 
  return(gr.qtl)
}

brew.colors <- function(group.list, palette="Blues"){
  suppressWarnings(
    colors <- RColorBrewer::brewer.pal(length(group.list), palette)[1:length(group.list)]
  )
  named.list <- setNames(colors, group.list) 
  return(named.list)
}


ggbio.line_circos <- function(){
  space.skip <- .005  #0.015
  rect.inter.n = 5
  trackWidth = 5
  # SNP grid
  alpha <- .75
  grid.background <-  "black"
  grid.line <-"darkgrey"
  grid.n <- 3
  # Annotation grid 
  grid.background.ann <- "black"
  grid.line.ann <-"black"
  grid.n.ann <- 3  
  # QTL 
  grid.n2 <- 1  
  grid.background2 <- "white"
  qtl.geom <- "rec"
  # seqlengths(gr.all) <- rep(1, length(gr.all))
  # "point", "line", "link", "ribbon", "rect", "bar", "segment", "hist", "scale", "heatmap", "ideogram", "text"
  FM[is.na(FM)] <-0
  gr.all <- FM %>% dplyr::mutate(SeqNames=paste0("chr",CHR)) %>%  
    data.frame() %>%
    biovizBase::transformDfToGr(seqnames = "SeqNames", start = "POS", end = "POS") 
  gr.all$GWAS.neg.log <- -log10(gr.all$P) 
  

  color.list <- c("Consensus SNP"="goldenrod2",
                  "Credible Set"="green",
                  "Lead GWAS SNP"="red",
                  "Nalls et al. (2019) GWAS -log10(P-value)"="red3",
                  "Fine-mapped Mean PP"="green4") 
  
  Fairfax <- c("Fairfax_2014_CD14",
               "Fairfax_2014_IFN",
               "Fairfax_2014_LPS2",
               "Fairfax_2014_LPS24")
  color.list <- append(color.list, 
                       brew.colors(Fairfax, palette="Oranges"))
  
  
  psychENCODE <- c("psychENCODE_eQTL")
  color.list <- append(color.list,  
                       brew.colors(psychENCODE, palette="Greys"))
  
  gtex <- c("GTEx_V7_Brain_Cerebellum",
              "GTEx_V7_Brain_Hippocampus",
              "GTEx_V7_Brain_Putamen_basal_ganglia",
              "GTEx_V7_Brain_Substantia_nigra")
  color.list <- append(color.list, 
                       brew.colors(gtex, palette="Purples"))
  
  xQTL <- c( "Brain_xQTL_Serve.eQTL",
             "Brain_xQTL_Serve.haQTL")
  color.list <- append(color.list, 
                       brew.colors(xQTL, palette="Blues"))
  
  loci <- unique(gr[["loci"]]$Gene)
  loci.subset <- loci[seq(2, length(loci), 2)] 
  # ggb <- ggbio::ggbio(radius=2)
  # for(qtl_file in c(xQTL, gtex, Fairfax, psychENCODE)){ 
  #   printer(qtl_file)
  #   gg.qtl <- ggbio.mergeQTL(FM, qtl_file=qtl_file)
  #   ggb <- ggb +  ggbio::circle(gg.qtl, 
  #                               aes(y=QTL.neg.log, fill=qtl_file), size=.5,
  #                               geom = qtl.geom, color=color.list[[qtl]],  alpha=alpha,
  #                               grid=T, grid.background=grid.background2, grid.n=grid.n2, grid.line=grid.line,
  #                               space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)
  # }
  ggb <- ggbio::ggbio(radius=2) +
    # Fairfax QTL
    ggbio::circle(ggbio.mergeQTL(FM, qtl_file="Fairfax_2014_CD14"), 
                  aes(y=QTL.neg.log, fill="Fairfax_2014_CD14"), size=.5,
                  geom = qtl.geom, alpha=alpha, color=color.list[["Fairfax_2014_CD14"]],
                  grid=T, grid.background=grid.background2, grid.n=grid.n2, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    ggbio::circle(ggbio.mergeQTL(FM, qtl_file="Fairfax_2014_IFN"), 
                  aes(y=QTL.neg.log, fill="Fairfax_2014_IFN"), size=.5,
                  geom = qtl.geom,  alpha=alpha, color=color.list[["Fairfax_2014_IFN"]],
                  grid=T, grid.background=grid.background2, grid.n=grid.n2, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    ggbio::circle(ggbio.mergeQTL(FM, qtl_file="Fairfax_2014_LPS2"), 
                  aes(y=QTL.neg.log, fill="Fairfax_2014_LPS2"), size=.5,
                  geom = qtl.geom,   alpha=alpha, color=color.list[["Fairfax_2014_LPS2"]],
                  grid=T, grid.background=grid.background2, grid.n=grid.n2, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    ggbio::circle(ggbio.mergeQTL(FM, qtl_file="Fairfax_2014_LPS24"), 
                  aes(y=QTL.neg.log, fill="Fairfax_2014_LPS24"), size=.5,
                  geom = qtl.geom,   alpha=alpha, color=color.list[["Fairfax_2014_LPS24"]],
                  grid=T, grid.background=grid.background2, grid.n=grid.n2, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +  
    # psychENCODE
    ggbio::circle(ggbio.mergeQTL(FM, qtl_file="psychENCODE_eQTL"), 
                  aes(y=QTL.neg.log, fill="psychENCODE_eQTL"), size=.5,
                  geom = qtl.geom,   alpha=alpha, color=color.list[["psychENCODE_eQTL"]],
                  grid=T, grid.background=grid.background2, grid.n=grid.n2, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    # ggbio::circle(ggbio.mergeQTL(FM, qtl_file="psychENCODE_cQTL"), aes(y=QTL.neg.log), size=.5,
    #               geom = "hist", color="purple2", alpha=alpha,
    #               grid=T, grid.background=grid.background, grid.n=grid.n, grid.line=grid.line,
    #               space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    # ggbio::circle(ggbio.mergeQTL(FM, qtl_file="psychENCODE_isoQTL"), aes(y=QTL.neg.log), size=.5,
    #               geom = "hist", color="purple3", alpha=alpha,
    #               grid=T, grid.background=grid.background, grid.n=grid.n, grid.line=grid.line,
    #               space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    # ggbio::circle(ggbio.mergeQTL(FM, qtl_file="psychENCODE_tQTL"), aes(y=QTL.neg.log), size=.5,
    #               geom = "hist", color="purple4", alpha=alpha,
    #               grid=T, grid.background=grid.background, grid.n=grid.n, grid.line=grid.line,
    #               space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    # Brain_xQTL_Serve
    ggbio::circle(ggbio.mergeQTL(FM, qtl_file="Brain_xQTL_Serve.eQTL"), 
                  aes(y=QTL.neg.log, fill="Brain_xQTL_Serve.eQTL"), size=.5,
                  geom = qtl.geom,   alpha=alpha, color=color.list[["Brain_xQTL_Serve.eQTL"]],
                  grid=T, grid.background=grid.background2, grid.n=grid.n2, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    ggbio::circle(ggbio.mergeQTL(FM, qtl_file="Brain_xQTL_Serve.haQTL"), 
                  aes(y=QTL.neg.log, fill="Brain_xQTL_Serve.haQTL"), size=.5,
                  geom = qtl.geom,   alpha=alpha, color=color.list[["Brain_xQTL_Serve.haQTL"]],
                  grid=T, grid.background=grid.background2, grid.n=grid.n2, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    # ggbio::circle(ggbio.mergeQTL(FM, qtl_file="Brain_xQTL_Serve.mQTL"), aes(y=QTL.neg.log), size=.5,
    #               geom = "point", color="cyan3", alpha=alpha,
    #               grid=T, grid.background=grid.background, grid.n=grid.n, grid.line=grid.line,
    #               space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    # ggbio::circle(ggbio.mergeQTL(FM, qtl_file="Brain_xQTL_Serve.cell-specificity-eQTL"), aes(y=QTL.neg.log), size=.5,
    #               geom = "point", color="cyan4", alpha=alpha,
    #               grid=T, grid.background=grid.background, grid.n=grid.n, grid.line=grid.line,
    #               space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    # GTEx
    ggbio::circle(ggbio.mergeQTL(FM, qtl_file="GTEx_V7_Brain_Substantia_nigra"), 
                  aes(y=QTL.neg.log, fill="GTEx_V7_Brain_Substantia_nigra"), size=.5,
                  geom = qtl.geom,  alpha=alpha, color=color.list[["GTEx_V7_Brain_Substantia_nigra"]],
                  grid=T, grid.background=grid.background2, grid.n=grid.n2, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    
    ggbio::circle(ggbio.mergeQTL(FM, qtl_file="GTEx_V7_Brain_Hippocampus"), 
                  aes(y=QTL.neg.log, fill="GTEx_V7_Brain_Hippocampus"), size=.5,
                  geom = qtl.geom,   alpha=alpha, color=color.list[["GTEx_V7_Brain_Hippocampus"]],
                  grid=T, grid.background=grid.background2, grid.n=grid.n2, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    
    ggbio::circle(ggbio.mergeQTL(FM, qtl_file="GTEx_V7_Brain_Putamen_basal_ganglia"), 
                  aes(y=QTL.neg.log, fill="GTEx_V7_Brain_Putamen_basal_ganglia"), size=.5,
                  geom = qtl.geom,   alpha=alpha, color=color.list[["GTEx_V7_Brain_Putamen_basal_ganglia"]],
                  grid=T, grid.background=grid.background2, grid.n=grid.n2, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    
    ggbio::circle(ggbio.mergeQTL(FM, qtl_file="GTEx_V7_Brain_Cerebellum"), 
                  aes(y=QTL.neg.log, fill="GTEx_V7_Brain_Cerebellum"), size=.5,
                  geom = qtl.geom,  alpha=alpha, color=color.list[["GTEx_V7_Brain_Cerebellum"]],
                  grid=T, grid.background=grid.background2, grid.n=grid.n2, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth) +

    # SNP Groups
    ggbio::circle(subset(gr.all, Consensus_SNP==T), aes(y=Support, fill="Consensus SNP"), size=.5,
                  geom = "bar",  alpha=alpha, color=color.list[["Consensus SNP"]],
                  grid=T, grid.background=grid.background, grid.n=grid.n, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth) +
    ggbio::circle(subset(gr.all, Support>0), aes(y=Support, fill="Credible Set"), size=.5,
                  geom = "bar", alpha=alpha, color=color.list[["Credible Set"]],
                  grid=T, grid.background=grid.background, grid.n=grid.n, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    ggbio::circle(subset(gr.all, leadSNP==T), aes(y=Support, fill="Lead GWAS SNP"), size=.5,
                  geom = "bar", alpha=alpha,  color=color.list[["Lead GWAS SNP"]],
                  grid=T, grid.background=grid.background, grid.n=grid.n, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth) + 
    # Fine-mapping 
    ggbio::circle(gr.all, aes(y=mean.PP, fill="Fine-mapped Mean PP"), size=.5,
                  geom = "point", alpha=alpha, color=color.list[["Fine-mapped Mean PP"]],
                  grid=T, grid.background=grid.background, grid.n=grid.n, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth) +
    # # GWAS
    ggbio::circle(gr.all, aes(y=GWAS.neg.log, fill="Nalls et al. (2019) GWAS -log10(P-value)"), size=.5,
                  geom = "point", alpha=alpha, color=color.list[["Nalls et al. (2019) GWAS -log10(P-value)"]],
                  grid=T, grid.background=grid.background, grid.n=grid.n, grid.line=grid.line,
                  space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth)  +
    # Loci
    # ggbio::circle(gr.all, geom = "ideogram", cytobands=T, size = .5) + 
    # ggbio::circle(gr.all, geom = "text", color="turquoise", aes(label=Gene), size=2, alpha=.7,
    #               space.skip = space.skip, rect.inter.n = rect.inter.n, trackWidth=trackWidth) + 
  
    ggbio::circle(gr[["loci"]], geom = "bar", color="red", aes(label=Gene, y=5)) +
    # ggbio::circle(subset(gr[["loci"]], Gene %in% loci.subset), geom = "text", color="red2", aes(label=Gene,  y=Height), size=5) +
    # Chromosome
    ggbio::circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
  
    ggb <- ggb + scale_fill_manual(values=color.list) + scale_color_manual(values=color.list)+
       theme(legend.key = element_rect(fill = "transparent"),
             rect = element_rect(fill = "transparent"),
             panel.background = element_rect(fill = "transparent"), # bg of the panel
             plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
             panel.border = element_blank(), 
             panel.grid.major = element_blank(), # get rid of major grid
             panel.grid.minor = element_blank(), # get rid of minor grid
             legend.background = element_rect(fill = "transparent"), # get rid of legend bg
             legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
             ) 
    
     
    # tiff(file.path("./Data/GWAS/Nalls23andMe_2019/_genome_wide/ggbio_circos_HD.tiff"), 
    #     bg = "transparent", res=1000, width = 11, height = 8, units = 'in')
    tiff(file.path("./Data/GWAS/Nalls23andMe_2019/_genome_wide/ggbio_circos.tiff"), 
             bg = "transparent", res=100, width = 11, height = 8, units = 'in')
    print(ggb)
    dev.off()
  return(ggb)
  }


XGR.gather_XGR_annotations <- function(gr.all, 
                                       merge_all=T, 
                                       overlap.cutoff=1){
  XGR.GR_length_filter <- function(gr.list, overlap.cutoff=1){
    printer("+ XGR:: Removing GRange objects with less than",overlap.cutoff,"overlapping region(s).")
    bools <- lapply(gr.list, function(e){length(e) >= overlap.cutoff}) %>% as.logical()
    return(gr.list[bools])
  }
  # Get XGR annotations
  xgr.encode_tfbs <- XGR_track(gr.all, lib.name = "ENCODE_TFBS_ClusteredV3_CellTypes") 
  xgr.encode_tfbs.merge <- unlist(XGR.GR_length_filter(xgr.encode_tfbs, overlap.cutoff=overlap.cutoff)) 
  xgr.encode_tfbs.merge$Dataset <- "ENCODE_TFBS_ClusteredV3_CellTypes"
  
  xgr.encode_DNase <- XGR_track(gr.all, lib.name = "ENCODE_DNaseI_ClusteredV3_CellTypes")
  xgr.encode_DNase.merge <- unlist(XGR.GR_length_filter(xgr.encode_DNase, overlap.cutoff=overlap.cutoff)) 
  xgr.encode_DNase.merge$Dataset <- "ENCODE_DNaseI_ClusteredV3_CellTypes"
  
  xgr.broad_histone <- XGR_track(gr.all, lib.name = "Broad_Histone")
  xgr.broad_histone.merge <- unlist(XGR.GR_length_filter(xgr.broad_histone, overlap.cutoff=overlap.cutoff)) 
  xgr.broad_histone.merge$Dataset <- "Broad_Histone"
  if(merge_all){
    xgr.all <- GRangesList(xgr.encode_tfbs.merge, xgr.encode_DNase.merge, xgr.broad_histone.merge) %>% unlist()
    return(xgr.all)
  } else{ 
    return(GRangesList("ENCODE_TFBS_ClusteredV3_CellTypes"=xgr.encode_tfbs.merge,
                       "ENCODE_DNaseI_ClusteredV3_CellTypes"=xgr.encode_DNase.merge,
                       "Broad_Histone"=xgr.broad_histone.merge))
  } 
}




####### OmicCircos ########

OmicCircos.finemapping_circos <- function(){
  library(OmicCircos)  # BiocManager::install("OmicCircos") 
  data(UCSC.hg19.chr) # Comes with OmicCircos
  # FM <- merge_finemapping_results()
  DF.list <- list("Consensus"=subset(FM, Consensus_SNP==T) %>% data.frame(),
                  "CS"=subset(FM, Support>0) %>% data.frame(),
                  "GWAS"=subset(FM, P<0.5) %>% data.frame(),
                  "all"=FM %>% data.frame()
  ) 
  gr <- ggbio.prepare_SNPgroups(FM)
  
  OmicCircos.prepare_segAnglePo <- function(DF){
    DF <- gr[["loci"]] %>% data.frame() %>%
      dplyr::select(chr=CHR, po=mean.POS, name=Gene,  min.POS=min.POS, max.POS=max.POS)
    db <- segAnglePo(seg.dat = DF, seg = colnames(DF))
    
  }
  
  plot(c(1,800) , c(1,800) , type="n", axes=FALSE, xlab="", ylab="", main="");
  OmicCircos::circos(cir= UCSC.hg19.chr) 
  
}




#  CIRCOS.get_locus_sizes <- function(dataset.path="./Data/GWAS/Nalls23andMe_2019/"){
#    multi.maps <- list.files(dataset.path, 
#                             pattern = "Multi-finemap_results.txt", 
#                             recursive = T, 
#                             full.names = T)
#    count.file.rows <- function(file.name){
#      out <- system(paste("wc -l",file.name), intern = T)
#      out.split <- strsplit(out, " ")[[1]]
#      rows <- as.integer(out.split[length(out.split)-1])
#      locus <- basename(dirname(dirname(file.name)))  
#      return( setNames(rows, locus) )
#    }
#    locus.list <- lapply(multi.maps, count.file.rows) %>% unlist()
#    locus.df <- data.frame(locus.list) %>% 
#      dplyr::mutate(Gene=row.names(.)) %>% 
#      `colnames<-`(c("Locus.size","Locus"))
#    return(locus.df)
#  }

# CIRCOS.SNPGroup_counts <- function(FM){
#   df0 <- FM %>% dplyr::group_by(Gene) %>% 
#     dplyr::summarise(CHR=mean(CHR), POS=mean(POS)) %>% 
#     dplyr::rename(Locus=Gene) # Mean position
#   
#   df1 <- FM %>% subset(Consensus_SNP==T) %>% 
#     dplyr::group_by(Gene) %>% 
#     count() %>% `colnames<-`(c("Locus","Consensus.size")) %>% data.frame()
#   df2 <- FM %>% subset(Support>0) %>% 
#     dplyr::group_by(Gene) %>%
#     count() %>% `colnames<-`(c("Locus","CS.size")) %>% data.frame()
#   df3 <- FM %>% subset(P<0.05) %>% 
#     dplyr::group_by(Gene) %>%
#     count() %>% `colnames<-`(c("Locus","GWAS.size")) %>% data.frame()
#   # df4 <- get.locus.sizes() 
#   merged.df <- base::merge(df0, df1,  all=T)
#   merged.df <- base::merge(merged.df,  df2,  all=T)
#   merged.df <- base::merge(merged.df, df3, all=T)
#   # merged.df <- base::merge(merged.df, df4, all=T)
#   melt.df <- reshape2::melt(merged.df,  
#                             id.vars=c("Locus","CHR","POS"),
#                             variable.name="SNP.Group",
#                             value.name="SNPs")
#   melt.df <- melt.df %>% arrange(CHR, POS)  
#   # melt.df[is.na(melt.df)] <- 0
#   return(melt.df)
# }



##################################################################
###################### circlize ######################
##################################################################
# library(circlize)
# 
# set.seed(999)
# n = 1000 
# FM_gene <- subset(FM, Gene=="LRRK2")
# df = data.frame(factors = sample(letters[1:8], n, replace = TRUE),
#                 x = rnorm(n), y = runif(n))
# circos.par("track.height" = 0.1)
# circos.initialize(factors = FM_gene$CHR, x = FM_gene$POS) 
# circos.track(factors = FM_gene$CHR, 
#              y = FM_gene$POS,
#              panel.fun = function(x, y) {
#                circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
#                            CELL_META$sector.index)
#                circos.axis(labels.cex = 0.6)
#              }) 
# col = rep(c("#FF0000", "#00FF00"), 4)
# circos.trackPoints(FM_gene$CHR, FM_gene$POS, FM_gene$P, col = col, pch = 16, cex = 0.5)
# circos.text(-1, 0.5, "text", sector.index = "a", track.index = 1)
# 
# 
# circos.clear()

