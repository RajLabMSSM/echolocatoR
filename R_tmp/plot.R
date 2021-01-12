# %%%%%%%%%%%%%%%%% #
####### PLOT #######
# %%%%%%%%%%%%%%%%% #


construct_SNPs_labels <- function(subset_DT,
                                  lead=T,
                                  method=T,
                                  consensus=T,
                                  remove_duplicates=T,
                                  verbose=F){
  printer("+ PLOT:: Constructing SNP labels...", v=verbose)
  labelSNPs <- data.table::data.table()
  subset_DT <- data.table::as.data.table(subset_DT)

  ## BEFORE fine-mapping
  if(lead){
    before <- subset( subset_DT %>% arrange(P), leadSNP == T)
    before$type <- "Lead"
    before$color <- "red"
    before$shape <- 18
    before$size <- 3
    labelSNPs <- rbind(labelSNPs, before, fill=T)
  }
  if(method){
    # AFTER fine-mapping
    after = subset(subset_DT, Support>0)
    if(dim(after)[1]>0){
      after$type <- "UCS"
      after$color<- "green3"
      after$shape <- 1
      after$size=3
      labelSNPs <- rbind(labelSNPs, after, fill=T)
    }
  }
  if(consensus & "Consensus_SNP" %in% colnames(subset_DT)){
    # Conensus across all fine-mapping tools
    cons_SNPs <- subset(subset_DT, Consensus_SNP==T)
    if(dim(cons_SNPs)[1]>0){
      cons_SNPs$type <- "Consensus"
      cons_SNPs$color <- "darkgoldenrod1"
      cons_SNPs$shape <- 1
      cons_SNPs$size=4
      labelSNPs <- rbind(labelSNPs, cons_SNPs, fill=T)
    }
  }
  # If there's duplicates only show the last one
  if(remove_duplicates){
    labelSNPs$rowID <- 1:nrow(labelSNPs)
    labelSNPs <- labelSNPs %>%
      dplyr::group_by(SNP) %>%
      dplyr::arrange(rowID) %>%
      dplyr::slice(n())
  } else {
    labelSNPs$rowID <- 1:nrow(labelSNPs)
    labelSNPs <- labelSNPs %>%
      dplyr::group_by(SNP, type) %>%
      dplyr::arrange(rowID) %>%
      dplyr::slice(n())
  }
  return(as.data.frame(labelSNPs))
}


# # cojo_DT = data.table::fread(file.path("Data/GWAS/Nalls23andMe_2019/LRRK2/COJO/COJO_results.txt"), sep="\t")
#
# # SNP_list = c("rs76904798","rs34637584","rs117073808")
# COJO_plot <- function(locus,
#                       cojo_DT,
#                       results_path,
#                       conditioned_snps,
#                       show_plot=T,
#                       save_plot=T){
#   locus <- basename(results_path)
#   # Label independent SNPs
#   ind_SNPs <- subset(cojo_DT, COJO.CS==T)
#   ind_SNPs <- cbind(ind_SNPs, type="Independent",color="purple")
#   labelSNPs <- ind_SNPs#rbind(labelSNPs, ind_SNPs)
#
#   effect_SNPs <- cojo_DT %>% arrange(desc(abs(COJO.Conditioned_Effect)))
#   effect_SNPs <- effect_SNPs[1:5,]
#   effect_SNPs$type <- "effect_SNPs"
#   effect_SNPs$color <- "blue3"
#   labelSNPs <- rbind(labelSNPs, effect_SNPs, fill=T)
#
#   # topEffect_snps <- cojo_DT %>% arrange(desc(COJO.Conditioned_Effect))
#   # topEffect_snps <- topEffect_snps$SNP[1:5]
#   # cojo_DT <- cojo_DT %>% dplyr::mutate(color = ifelse(COJO.CS | COJO.Conditioned_Effect %in% topEffect_snps,
#   #                                                     ifelse(COJO.CS & COJO.Conditioned_Effect %in% topEffect_snps, "purple",
#   #                                                            ifelse(COJO.CS & !(COJO.Conditioned_Effect %in% topEffect_snps), "lightpurple", "lightblue") ), NA)
#   # )
#
#   spacing <- if(length(cojo_DT$SNP)>1000){250000}else{50000}
#   # roundBreaks <- seq(plyr::round_any(min(cojo_DT$POS),10000), max(cojo_DT$POS),spacing)
#
#   cp <- ggplot(cojo_DT, aes(x=POS, y = -log10(P), label=SNP, color= -log10(P))) +
#     # ylim(yLimits1) +
#     geom_hline(yintercept=0, alpha=.5, linetype=1, size=.5) +
#     geom_point(alpha=.5) +
#     geom_segment(aes(xend=POS, yend=0, color= -log10(P) ), alpha=.5) +
#     geom_point(data=labelSNPs, pch=21, size=4, colour=labelSNPs$color, stroke=1) +
#     geom_label_repel(data=labelSNPs, aes(label=SNP), color=NA, nudge_x = .5, box.padding = .5,
#                      label.size=NA, alpha=.8, seed = 1) +
#     geom_label_repel(data=labelSNPs, aes(label=SNP), color=labelSNPs$color, segment.alpha = .5,
#                      nudge_x = .5, box.padding = .5, fill = NA, alpha=1, seed = 1, show.legend = T ) +
#     labs(title=paste(locus,": Conditional & Stepwise Results (COJO)"),
#          subtitle = paste("Purple = Independent signals from stepwise procedure\n",
#                           "Blue = Residual effects conditioned on:",conditioned_snps),
#          y="-log10(p-value)", x="Position", color="-log10(p-value)") +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5) )
#   # scale_x_continuous(breaks = roundBreaks)
#
#   if(save_plot){
#     png(filename=file.path(results_path,"COJO/COJO_plot.png"), width = 1000, 800)
#     print(cp)
#     dev.off()
#   }
#   if(show_plot){print(cp)}else{return(cp)}
# }
#
#
#
# snp_plot <- function(finemap_dat,
#                      LD_matrix,
#                      locus,
#                      method="original",
#                      show_plot=T,
#                      subtitle=NA,
#                      multi = T,
#                      LD_SNP = NA){
#   {
#   # finemap_dat <- data.table::fread("Data/GWAS/Nalls23andMe_2019/LRRK2/Multi-finemap/Multi-finemap_results.txt", sep="\t"); LD_matrix <- readRDS("Data/GWAS/Nalls23andMe_2019/LRRK2/plink/UKB_LD.RDS")
#   score_dict <- setNames(c("-log10(P-value)", "PIP","PIP","PP", "PIP", "Conditional Probability","-log10(P-value)"),
#                          c("original","SUSIE","POLYFUN_SUSIE", "ABF", "FINEMAP", "COJO","COLOC"))
#   # X-tick spacing
#   # spacing <- if(length(finemap_dat$SNP)>1000){250000}else{50000}
#   # roundBreaks <- seq(plyr::round_any(min(finemap_dat$POS),10000), max(finemap_dat$POS), spacing)
#   # yLimits2 <- c(0,1.1)
#
#   # Merge with LD info
#   # If not specified, identify a lead SNP by summing the PPs from each fine-mapping method
#
#   LD_SNP <- subset(finemap_dat, leadSNP==T)$SNP
#   LD_sub <- LD_with_leadSNP(LD_matrix, LD_SNP)
#   DT <- data.table:::merge.data.table(finemap_dat, LD_sub, by = "SNP") %>%
#     dplyr::mutate(Mb=POS/1000000)
#
#
#   if(method=="original"){
#     r2_multiply = 150
#     is.na(DT$Probability) <- 0
#     p <- ggplot(data = DT, aes(x=Mb, y= -log10(P), label=SNP, color= r2)) +
#       geom_hline(yintercept=0,alpha=.5, linetype=1, size=.5) +
#       stat_smooth(data=DT, aes(x=Mb, y=r2*r2_multiply, fill=r2),color="firebrick1",
#                   se = F, formula = y ~ x,
#                   method = 'loess', span=.1) +
#       geom_hline(yintercept=0,alpha=.5, linetype=1, size=.5) +
#       geom_hline(yintercept= -log10(5e-8), alpha=.5, linetype=2, size=.5, color="black")
#     title <- "GWAS"#paste0(locus," : Before fine-mapping")
#     tag_SNPs <- labelSNPs <- construct_SNPs_labels(DT, lead=T, method = F, consensus = T)
#     subtitle <- if(is.na(subtitle)){paste0(length(DT$SNP)," SNPs")}else{subtitle}
#   } else {
#
#     # COLOC
#     if(method=="COLOC"){
#       DT <- DT %>% dplyr::rename(CS = "Colocalized" )
#       title <- paste0(locus," : Colocalization (",method,")")
#       p <- ggplot(data = DT, aes(x=POS, y= -log10(P), label=SNP, color= r2 ))
#     # Fine-mapping methods
#     } else if (multi){
#       title <- method #paste0(locus," : After fine-mapping (",method,")")
#       DT <- DT %>% dplyr::rename(Probability = paste0(method,".PP"),
#                                  CS = paste0(method,".CS"))
#       is.na(DT$Probability) <- 0
#       p <- ggplot(data = DT, aes(x=POS, y=Probability, color= r2 )) +
#         ylim(c(0,1.2))
#       subtitle <- if(is.na(subtitle)){
#         paste0(length(subset(DT, CS>0)$SNP), " Candidate SNP(s)")
#         }else{subtitle}
#     } else {
#       title <- paste0(locus," : After fine-mapping (",method,")")
#       p <- ggplot(data = DT, aes(x=POS, y=Probability,  color= r2 )) +
#         ylim(c(0,1.2))
#     }
#     labelSNPs <- construct_SNPs_labels(DT, lead=T, method = T, consensus = F)
#     tag_SNPs <- subset(labelSNPs, CS>0)
#   }
#
#   p <- p +
#     geom_point() + # alpha=.5
#     scale_color_gradient(low="blue", high="red", limits = c(0,1)) +
#     # geom_segment(aes(xend=POS, yend=0, color= -log10(P)), alpha=.5) +
#     geom_point(data=labelSNPs, pch=21, fill=NA, size=4, color=labelSNPs$color, stroke=1) +
#     geom_point(data=subset(DT, SNP==LD_SNP), pch=18, fill=NA, size=4, color="red") +
#     # Labels (one for background, one for text)
#     geom_label_repel(data=tag_SNPs, aes(label=SNP),
#                      color=NA,
#                      nudge_x = .5,
#                      fill="black",
#                      box.padding = .5,
#                      label.size=NA,
#                      alpha=.6,
#                      seed = 1) +
#     geom_label_repel(data=tag_SNPs, aes(label=SNP),
#                      color=tag_SNPs$color,
#                      segment.alpha = .5,
#                      nudge_x = .5,
#                      box.padding = .5,
#                      fill = NA,
#                      alpha=1,
#                      seed = 1) +
#     labs(title = title,
#          subtitle = subtitle,
#          y = score_dict[[method]],
#          x = "Position",
#          color = bquote(paste(r^2," with ",.(LD_SNP) ) )) +
#     theme_bw()+
#     theme(plot.title = element_text(hjust = 0.5),
#           plot.subtitle = element_text(hjust = 0.5),
#           rect = element_rect(fill = "transparent") ) +
#     # scale_x_continuous(breaks = roundBreaks) +
#     scale_color_gradient(low="blue", high="red", limits = c(0,1))
#   if(show_plot){print(p)} else{ return(p) }
#   }
#
# }
#
#
#
# multi_finemap_plot <- function(finemap_dat,
#                                LD_matrix,
#                                results_path,
#                                finemap_method_list,
#                                locus,
#                                conditioned_snps,
#                                original=T,
#                                save_plot=T,
#                                ncols=1,
#                                width=500, #500,
#                                height=1000 #1000
#                                  ){
#   # locus <- "LRRK2"
#   # results_path <- file.path("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019",locus)
#   # finemap_dat <- data.table::fread(file.path(results_path,"Multi-finemap/Multi-finemap_results.txt"))
#   # LD_matrix <- readRDS(file.path(results_path,"plink/UKB_LD.RDS"))
#   # method_list=c("ABF","SUSIE","POLYFUN_SUSIE","FINEMAP")
#
#   method_list <- if(original){c("original", finemap_method_list)}else{finemap_method_list}
#
#   # Assemble plots in list
#    plot_list <- lapply(method_list, function(method,
#                                             finemap_dat.=finemap_dat,
#                                             LD_matrix.=LD_matrix){
#       printer("\n Plotting...",method)
#       if(method=="COJO"){
#         p <- COJO_plot(cojo_DT = finemap_dat.,
#                   results_path = results_path,
#                   show_plot = F,
#                   save_plot = T,
#                   conditioned_snps = conditioned_snps)
#       } else{
#         p <- snp_plot(finemap_dat = finemap_dat.,
#                       LD_matrix = LD_matrix.,
#                       locus = locus,
#                       method = method,
#                       show_plot = F,
#                       multi = T)
#       }
#       return(p)
#   })
#   # Multi-plot
#   # grDevices::graphics.off()
#   cp <- patchwork::wrap_plots(plot_list, ncol = ncols)
#   print(cp)
#   # plot_list <- lapply(1:3,function(e){ggplot(iris, aes(x=Sepal.Width,y=Sepal.Length))+geom_point()})
#
#   # library(gridExtra)
#   # cp <- do.call("grid.arrange", c(plot_list, ncol=ncols))
#   # plot(cp)
#
#   # Rmisc::multiplot(plotlist = plot_list, cols = ncols)
#   ## Adjust plotting parameters
#   if(length(finemap_method_list)>3){ncols <- 2; width <- width*2}
#
#   # Save plot
#   if(save_plot){
#     # png(filename=file.path(results_path,"Multi-finemap/multi_finemap_plot.png"),
#     #     width = width, height)
#     # print(cp)
#     # dev.off()
#     ggplot2::ggsave(file=file.path(results_path,"Multi-finemap/multi_finemap_plot.png"),
#            cp, width = 5, height = 12,
#            bg = "transparent", dpi=400)
#   }
#   return(cp)
# }
#
#
#
#
# eQTL_barplots <- function(subset_path_list, group_list, SNP_list,
#                           x_lab="Population", y_lab="Effect",
#                           snp_col="snps", effect_col="beta",
#                           title="", writeCSV=F){
#   library(ggplot2)
#   dt = data.table::data.table()
#   i=1
#   for (s in subset_path_list){
#     dt_sub <- subset(data.table::fread(s), eval(parse(text = snp_col)) %in% SNP_list)
#     dt_sub$Group = group_list[i] #strsplit(strsplit(s,"/")[[1]][4], "_")[[1]][2]
#     dt = rbind(dt, dt_sub)
#     i=i+1
#   }
#   g <- ggplot(data=dt, aes(x=Group, y=eval(parse(text = effect_col)), fill=Group)) +
#     geom_col() + facet_grid(~eval(parse(text = snp_col))) +
#     labs(title = title, x = x_lab, y = y_lab)
#   print(g)
#   if(writeCSV!=F){
#     data.table::fwrite(dt, writeCSV)
#   }
# }
#
#
#
#
#
#
# plot_mean.PP <- function(results_path, top_snps=20){
#   locus <- basename(results_path)
#   gene_DT <- data.table::fread(file.path(results_path,"Multi-finemap/Multi-finemap_results.txt"),
#                                nThread = 4)
#   gene_DT <- find_consensus_SNPs(gene_DT,)
#   lead.rsid <- subset(gene_DT, leadSNP==T)$SNP
#
#   LD_matrix <- readRDS(file.path(results_path,"plink/LD_matrix.RData"))
#   # corrplot.p3 <- data.table::fread(file.path(results_path,"LDlink_1KGphase3.txt")) %>%
#   #   data.frame(row.names = 1)  %>%
#   #   as.matrix()
#   # corrplot.p1 <- LD_matrix[c("rs11175620","rs7294619","rs76904798"),c("rs11175620","rs7294619","rs76904798")]
#   #
#   # corrplot::corrplot(corrplot.p3,
#   #                    tl.col = "black",
#   #                    col = colorRampPalette(c("blue","white","red"))(200),
#   #                    addCoef.col = "white")
#   ld.dat <- data.frame(r2=LD_matrix[lead.rsid,], SNP=names(LD_matrix[lead.rsid,]))
#
#   gene_DT <- data.table:::merge.data.table(gene_DT,
#                                 data.table::data.table(ld.dat),
#                                 by="SNP", all.x = T)
#
#   fulldat <- data.table::fread(file.path(results_path, "LRRK2_Nalls23andMe_2019_subset.txt"))
#   # Prepare data
#   pp.dat <- gene_DT %>% arrange(desc(mean.PP)) %>% unique() %>% head(top_snps)
#   pp.dat$SNP <- factor(pp.dat$SNP, levels = pp.dat$SNP)
#
#   # Plot
#   ppp <- ggplot(pp.dat, aes(x=SNP, y=mean.PP, fill=r2)) +
#     geom_col() +
#     labs(title=locus,
#          subtitle = "Mean Fine-mapping Posterior Probability",
#          fill=paste0("r2 with\n",lead.rsid)) +
#     geom_text(stat = 'identity',aes(label=SNP, hjust=-.2, vjust= -.5,   srt=45), color="green") +
#     geom_text(data = subset(pp.dat,SNP==lead.rsid),
#               stat = 'identity',aes(label=SNP, hjust=-.2, vjust= -.5,   srt=45),color="red") +
#     # geom_text(data = subset(pp.dat,SNP==lead.rsid),
#     #           stat = 'identity',aes(label=SNP, hjust=-.2, vjust= -.5,   srt=45),color="red") +
#     geom_text(data = subset(pp.dat, Support>=2),
#               stat = 'identity',aes(label=SNP, hjust=-.2, vjust= -.5,   srt=45),color="goldenrod2") +
#     theme_classic() +
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(),
#           plot.title = element_text(hjust = 0.5),
#           plot.subtitle = element_text(hjust = 0.5),
#           legend.position = c(0.95, 0.8),
#           rect = element_rect(fill = "transparent"),
#           panel.background = element_rect(fill = "transparent")) +
#     scale_fill_gradient(low="blue", high="red", limits = c(0,1)) +
#     ylim(c(0,min(max(pp.dat$mean.PP)*1.5,1.1) )) +
#     xlim(c(0,top_snps+2))
#   print(ppp)
#
#   ggsave(filename = file.path(results_path,"Multi-finemap","mean.PP.png"),
#          plot = ppp, dpi = 600, width=8, height=7, bg="transparent")
#   return(ppp)
# }
#
#

# remotes::install_github("tylermorganwall/rayshader")
# library(rayshader)
# library(ggplot2)
# gg = ggplot(diamonds, aes(x, depth)) +
#   stat_density_2d(aes(fill = stat(nlevel)),
#                   geom = "polygon",
#                   n = 100,bins = 10,contour = TRUE) +
#   facet_wrap(clarity~.) +
#   scale_fill_viridis_c(option = "A")
# plot_gg(gg,multicore=TRUE,width=5,height=5,scale=250)
#



