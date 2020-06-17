## ---- include = FALSE---------------------------------------------------------
root.dir <- "~/Desktop"
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  root.dir = root.dir
)  
knitr::opts_knit$set(root.dir = root.dir)
# knitr::opts_chunk$get("root.dir")

## ----setup, root.dir="~/Desktop"----------------------------------------------
library(echolocatoR) 

## -----------------------------------------------------------------------------
library(ggplot2)
data("BST1"); data("LD_matrix"); data("locus_dir");
locus_dir <- file.path("~/Desktop",locus_dir)
finemap_DT <- BST1 

## ----trk_plot-----------------------------------------------------------------
trk_plot <- GGBIO.plot(finemap_dat=finemap_DT, 
                       LD_matrix=LD_matrix, 
                       locus_dir=locus_dir, 
                       XGR_libnames=NULL, 
                       save_plot=F,
                       show_plot=F)
trk_plot <- trk_plot + theme(strip.background = element_rect(fill="gray"),
                             text = element_text(size = 3, color="blue"))

## ----trk_plot.xgr-------------------------------------------------------------
trk_plot.xgr <- GGBIO.plot(finemap_dat=finemap_DT, 
                           LD_matrix=LD_matrix, 
                           locus_dir=locus_dir, 
                           XGR_libnames=c("ENCODE_TFBS_ClusteredV3_CellTypes"), 
                           save_plot=F,
                           show_plot=F)

## ----trk_plot.roadmap---------------------------------------------------------
trk_plot.roadmap <- GGBIO.plot(finemap_dat=finemap_DT, 
                               LD_matrix=LD_matrix, 
                               locus_dir=locus_dir, 
                               XGR_libnames=NULL, 
                               Roadmap=T, 
                               Roadmap_query="monocyte", 
                               save_plot=F, 
                               show_plot=F)

## ----trk_plot.nott_2019-------------------------------------------------------
trk_plot.nott_2019 <- GGBIO.plot(finemap_dat=finemap_DT, 
                                 LD_matrix=LD_matrix, 
                                 locus_dir=locus_dir, 
                                 XGR_libnames=NULL, 
                                 Nott_epigenome=T, 
                                 Nott_binwidth = 100,
                                 Nott_regulatory_rects = T, 
                                 Nott_show_placseq = T,
                                 save_plot=F,
                                 show_plot=F)

