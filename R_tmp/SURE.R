
# ^^^^^^^^^^^^^^ SURE ^^^^^^^^^^^^^^
# Data download source:
# https://sure.nki.nl

# Reference:::
# Arensbergen, Joris van, Ludo Pagie, Vincent D. FitzPatrick, Marcel de Haas, Marijke P. Baltissen, Federico Comoglio, Robin H. van der Weide, et al. “High-Throughput Identification of Human SNPs Affecting Regulatory Element Activity.” Nature Genetics 51, no. July (2019). https://doi.org/10.1038/s41588-019-0455-2.
# ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^ 

SURE.preprocess_data <- function(results_path="./Data/GWAS/Nalls23andMe_2019/LRRK2", cor.plot=F){  
  FM_gene <- data.table::fread(file.path(results_path,"Multi-finemap","Multi-finemap_results.txt"), nThread = 4)
  FM_gene <- cbind(Locus=basename(results_path), FM_gene)
  sure.path <- file.path(results_path,"SURE")
  
  sure <- data.table::fread(list.files(sure.path, pattern = "snp_browser_selected",full.names = T), nThread = 4)
  sure <- data.table:::merge.data.table(FM_gene, sure, 
                                            by.x = "SNP", by.y = "id", 
                                            all.x = T)
  
  
  # sure.dat %>% arrange(desc(mean.PP))
  prob.cols <- grep(".PP",colnames(sure), value = T)
  sure.cols <- grep("_mean",colnames(sure), value = T)
  # Corrplot
  if(cor.plot){
    cor.dat <- (sure %>% dplyr::mutate(GWAS.neg.log.P=-log10(P)))[,c("GWAS.neg.log.P",prob.cols, sure.cols)]
    cor.dat[is.na(cor.dat)] <-0
    corrplot::corrplot(cor(cor.dat), 
                       tl.cex = .5, 
                       method = "ellipse", 
                       type = "upper",
                       diag = F,  
                       addCoef.col = "black", 
                       cl.cex = .8, )
  }
 
  # Melt
  suppressWarnings(sure.dat <- sure %>%
                       data.table::melt.data.table(id.vars = c("Locus","SNP","CHR","POS","leadSNP","Support","Consensus_SNP","mean.PP","P"), 
                                                   measure.vars = c(prob.cols, sure.cols), 
                                                   variable.name = "Metric",
                                                   value.name = "Value"))

  
  sure.dat$Allele <- as.character(NA)
  sure.dat[grep(pattern = "_ref_", sure.dat$Metric),"Allele"] <- "Ref"
  sure.dat[grep(pattern = "_alt_", sure.dat$Metric),"Allele"] <- "Alt"
  sure.dat$cell_line <- as.character(NA)
  sure.dat[grep(pattern = "k562_", sure.dat$Metric),"cell_line"] <- "k562"
  sure.dat[grep(pattern = "hepg2_", sure.dat$Metric),"cell_line"] <- "hepg2"
  
  sure.dat$Mb <- sure.dat$POS/1000000
  return(sure.dat)
}
 
SURE.report <- function(sure.dat, locus=NA){ 
  # CredSet
  CS <- subset(sure.dat, Support>0 & !is.na(Allele) ) %>% arrange(desc(mean.PP)) 
  notNA <- CS %>% dplyr::group_by(SNP) %>% dplyr::summarise(notNA =! any(is.na(Value)) )
  printer("SURE::",sum(notNA$notNA),"/",nrow(notNA),
          paste0("(",(sum(notNA$notNA)/nrow(notNA)*100),"%)"),
          "of Credible Set SNPs were tested by the SURE assay.")
  CS.df <- data.frame(n.CredSet=nrow(notNA), 
                      n.CredSet.SURE=sum(notNA$notNA),
                      perc.CredSet.SURE=(sum(notNA$notNA)/nrow(notNA)*100))
  # Consensus
  consensus <- subset(sure.dat, Consensus_SNP==T & !is.na(Allele) ) %>% arrange(desc(mean.PP)) 
  notNA <- consensus %>% dplyr::group_by(SNP) %>% dplyr::summarise(notNA =! any(is.na(Value)) )
  printer("SURE::",sum(notNA$notNA),"/",nrow(notNA),
          paste0("(",(sum(notNA$notNA)/nrow(notNA)*100),"%)"),
          "of Consensus SNPs were tested by the SURE assay.")
  Consensus.df <- data.frame(n.Consensus=nrow(notNA), 
                             n.Consensus.SURE=sum(notNA$notNA), 
                             perc.Consensus.SURE=(sum(notNA$notNA)/nrow(notNA)*100))
  report.df <- cbind(Locus=locus, CS.df, Consensus.df) %>% data.table::data.table()
  return(report.df)
}



SURE.track_plot <- function(sure.dat=NULL, 
                            results_path="./Data/GWAS/Nalls23andMe_2019/LRRK2"){
  if(is.null(sure.dat)){
    sure.dat <- SURE.preprocess_data(results_path=results_path)  
  }
  prob.cols <- grep(".PP",unique(sure.dat$Metric), value = T)
  snp.labels <- construct_SNPs_labels(sure.dat) %>%  
    dplyr::group_by(SNP) %>%  
    remove_missing(vars = c("Value")) %>% 
    slice(tail(row_number(), 1)) #%>% dplyr::select(SNP,CHR,POS,P,Value,Mb,type,color)
 
  library(patchwork)
  # GWAS
  sp <- ggplot(sure.dat) + 
    geom_point(aes(x=Mb, y=-log10(P), color=-log10(P))) +  
    geom_point(data=snp.labels, pch=21, fill=NA, size=4, color=snp.labels$color, stroke=1, alpha=0.8, 
               aes(x=Mb, y=-log10(P))) +
    scale_color_gradient(low="blue", high="red") +
    labs(y="-log(P-value)", color="GWAS\n-log(P-value)") + 
    theme_dark() + 
    theme(rect = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(angle = 0), 
          strip.text = element_text(size=9 )) + 
    facet_grid("Nalls et al. (2019)\nGWAS"~., scales = "free_y") + 
    geom_label_repel(data = snp.labels, aes(x=Mb, y=-log10(P), label=SNP),
                     segment.alpha = .5, 
                     color=snp.labels$color,
                     size=3,
                     # nudge_x = .5, 
                     # box.padding = .5,  
                     fill = "white", 
                     alpha=.8, 
                     seed = 1) + 
    # Fine-mapping  
    ggplot(subset(sure.dat, Metric %in% prob.cols)) +
    geom_point(aes(x=Mb, y=Value, color=Value)) +  
    scale_color_gradient(low="grey", high="green", breaks=c(0,.5,1)) + 
    scale_y_continuous(breaks = c(0,.5,1)) + 
    labs(y="Probability", color="Fine-mapping\nProbability") + 
    theme_dark() + 
    theme(rect = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(angle = 0), 
          strip.text = element_text(size=9 )) +
    facet_grid(gsub("\\.","\n",Metric)~.) + 
    # SuRE assay
    ggplot(subset(sure.dat, Metric %in% c("k562_ref_mean","k562_alt_mean","hepg2_ref_mean","hepg2_alt_mean"))) +
    geom_point(aes(x=Mb, y=Value, color=Value)) + 
    scale_color_viridis_c() +
    labs(y="Mean Effect", color="SuRE\nMean Effect") + 
    theme_dark() + 
    theme(rect = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(angle = 0), 
          strip.text = element_text(size=9 )) +
    facet_grid(gsub("\\.","\n",Metric)~.) +
    # Overall
    plot_layout(ncol = 1, heights = c(.3,1,1)) + 
    plot_annotation(title = paste(unique(sure.dat$Locus), collapse=","), 
                    theme = theme(plot.title = element_text(hjust = 0.5),
                                  panel.border = element_rect(fill = "transparent"),
                                  panel.background = element_rect(fill = "transparent"),
                                  rect = element_rect(fill = "transparent")))   
  print(sp)
  if(save_plot){
    ggsave(file.path(results_path,"SURE","SURE.png"), height = 13, width = 10, bg="transparent")
  }
  return(sp)
}
  
# need to first figure out how to programmatically gather all data...
# SURE.report_all_loci <- function(dataset.path="./Data/GWAS/Nalls23andMe_2019"){
#   locus.paths <- list.dirs(dataset.path, recursive = F) 
#   locus.paths <- locus.paths[locus.paths!=file.path(dataset.path,"_genome_wide")]
#   
#   report.DF <- lapply(locus.paths, function(locus.path){
#     locus <- basename(locus.path)
#     sure.dat <-  SURE.preprocess_data(results_path = locus.path)
#     report.df <- SURE.report(sure.dat, locus = locus)
#     return(report.df)
#   }) 
# }


