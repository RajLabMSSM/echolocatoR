REPLICATION.compare_LD_panels <- function(root="/sc/arion/projects/pd-omics/brian/Fine_Mapping",
                                           dataset="Data/GWAS/Nalls23andMe_2019",
                                           no_no_loci = c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1",
                                                          "RETREG3","UBTF","FAM171A2","MAP3K14","CRHR1","MAPT-AS1","KANSL1","NSF","WNT3")){
  requireNamespace("ggplot2")
  requireNamespace("patchwork")
  # merged_UKB <- data.table::fread(file.path(root,"Data/GWAS/Nalls23andMe_2019/_genome_wide/merged_UKB.csv.gz"))
  merged_UKB <- merge_finemapping_results(dataset = file.path(root,dataset),
                                          LD_reference = "UKB",
                                          minimum_support = 1, include_leadSNPs = TRUE)
  # merged_1KG <- data.table::fread(file.path(root,"Data/GWAS/Nalls23andMe_2019/_genome_wide/merged_1KGphase3.csv.gz"))
  merged_1KG <- merge_finemapping_results(dataset = file.path(root,dataset),
                                          LD_reference = "1KGphase3",
                                          minimum_support = 1, include_leadSNPs = TRUE)
  # Merge and remove non-overlapping loci
  common_loci <- dplyr::intersect(unique(merged_UKB$Locus),
                                  unique(merged_1KG$Locus))

  common_loci <- common_loci[!common_loci %in% no_no_loci]
  # merged_PD <- rbind(subset(merged_UKB, Locus %in% common_loci) %>%
  #                      dplyr::mutate(LD_panel="UKB"),
  #                    subset(merged_1KG, Locus %in% common_loci) %>%
  #                      dplyr::mutate(LD_panel="1KG")
  #                    )
  merged_PD <- data.table::merge.data.table(#UKB
    subset(merged_UKB, Locus %in% common_loci) %>%
      `colnames<-`(paste(colnames(.),"UKB",sep='.')) %>%
      dplyr::rename(Locus=Locus.UKB, SNP=SNP.UKB),
    # 1KG
    subset(merged_1KG, Locus %in% common_loci) %>%
      `colnames<-`(paste(colnames(.),"1KG",sep='.')) %>%
      dplyr::rename(Locus=Locus.1KG, SNP=SNP.1KG),
    by = c("Locus","SNP")
  ) %>%
    dplyr::rename(P=P.UKB)

  overlap <- merged_PD %>%
    dplyr::mutate(leadSNP_overlap=leadSNP.UKB>0 & leadSNP.1KG>0,
                  UCS_overlap=Support.UKB>0 & Support.1KG>0,
                  Consensus_overlap=Consensus_SNP.UKB>0 & Consensus_SNP.1KG>0)
  overlap_summary <- overlap %>%
    dplyr::group_by(Locus) %>%
    dplyr::summarise(leadSNP_overlap=sum(leadSNP_overlap, na.rm = TRUE),
                     leadSNP_prop=sum(leadSNP_overlap, na.rm = TRUE)/sum(leadSNP_overlap, na.rm = TRUE),

                     UCS_overlap=sum(UCS_overlap, na.rm = TRUE),
                     UCS_prop=sum(UCS_overlap, na.rm = TRUE)/sum(UCS_overlap, na.rm = TRUE),

                     Consensus_overlap=sum(Consensus_overlap, na.rm = TRUE),
                     Consensus_prop=sum(Consensus_overlap, na.rm = TRUE)/sum(Consensus_overlap, na.rm = TRUE)
    ) %>% data.frame()

  #### heatmap ####
  heat <- REPLICATION.compare_PP_heatmap(merged_PD = merged_PD)

  #### scatter plot ####
  methods <- gsub("\\.PP\\.UKB","",grep("*\\.PP\\.UKB", colnames(merged_PD), value = TRUE))
  pp_list <- lapply(methods, function(m){
    print(m)
    REPLICATION.compare_PP_scatterplot(merged_PD=merged_PD,
                                       col1=paste0(m,".PP.UKB"),
                                       col2=paste0(m,".PP.1KG"),
                                       title=paste(m,"PP"),
                                       max_labels = 2,
                                       show_plot  = FALSE)
  }) %>% `names<-`(methods)

  #### boxplot ####
  bp <- REPLICATION.compare_PP_pairedplot(merged_PD=merged_PD,
                                       col1 = "mean.PP.UKB",
                                       col2 = "mean.PP.1KG",
                                       title = "Paired plot: mean.PP")
  pp_list$paired_meanPP <- bp

  pp_wrap <- patchwork::wrap_plots(pp_list, ncol = 2) +
    patchwork::plot_annotation(tag_levels = "a", title = "UKB LD vs. 1KG LD")
  print(pp_wrap)

  if(save_plot!=FALSE){
    save_path <- file.path(results_dir, "GWAS/Nalls23andMe_2019/_genome_wide","LD_comparison/UKB_vs_1KG.PP.png")
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    ggsave(save_path, pp_wrap, dpi = 300, width = 10, height = 12)
  }


}




REPLICATION.compare_singleton_results <- function(singleton_url="/pd-omics/data/Singleton_FINEMAP_PD/pd_meta5_sum_stats_fm_results.csv.gz",
                                                  results_dir="/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data",
                                                  save_path=file.path(results_dir,"GWAS/Nalls23andMe_2019/_genome_wide/Singleton_comparison/Schilder_vs_Grenn.png")){

  requireNamespace("ggplot2")
  requireNamespace("patchwork")
  # Originally from: https://storage.googleapis.com/nihnialng-share-f23bef932184/pd_meta5_sum_stats_fm_results.csv.gz
  singleton <- data.table::fread(singleton_url, nThread=1)
  # Merging by SNP yields waaaaaay more hits than CHR/POS, even after liftover.
  merged_PD <- data.table::merge.data.table(merged_dat,
                                            singleton,
                                            by="SNP")
  overlap <- subset(merged_PD, FINEMAP.PP>=.95 & prob>=.95)
  # overlap <- subset(merged_PD,  Support>0 )
  overlap[is.na(overlap)] <- 0
  # overlap <- subset(overlap, (!is.na(FINEMAP.PP)) & (!is.na(prob)))

  pearson <- cor.test(merged_PD$FINEMAP.PP, merged_PD$prob, method="pearson")
  spearman <- cor.test(merged_PD$FINEMAP.PP, merged_PD$prob, method="spearman")
  spearman <- cor.test(overlap$FINEMAP.PP, overlap$prob, method="spearman")

  #### scatter plot ####
  methods <- gsub("\\.PP","",grep("*\\.PP$", colnames(merged_PD), value = TRUE))
  pp_list <- lapply(methods, function(m){
    print(m)
    REPLICATION.compare_PP_scatterplot(merged_PD=merged_PD,
                                       col1=paste0(m,".PP"),
                                       col2="prob",
                                       title=m,
                                       max_labels = 2,
                                       filter_str = "prob > 0",
                                       xlabel=paste0(m," PP (Schilder et al.)"),
                                       ylabel="FINEMAP PP (Grenn et al.)",
                                       show_count = FALSE,
                                       show_plot  = FALSE)
  }) %>% `names<-`(methods)

  bp <- REPLICATION.compare_PP_pairedplot(merged_PD=merged_PD,
                                          col1 = "FINEMAP.PP",
                                          col2 = "prob",
                                          title = "Paired plot: inter-study comparison",
                                          variable_key = c("FINEMAP.PP"="(Schilder et al.)",
                                                           "prob"="(Grenn et al.)"),
                                          filter_str="FINEMAP.PP>0 & prob>0",
                                          ylabel = "FINEMAP PP",
                                          show_plot  = FALSE)
  pp_list$paired_FINEMAP.PP <- bp


  pp_wrap <- patchwork::wrap_plots(pp_list, ncol = 2) +
    patchwork::plot_annotation(tag_levels = "a",
                               title = "Fine-mapping: Schilder et al. vs. Grenn et al.")
  if(show_plot) print(pp_wrap)
  if(save_path!=FALSE){
    dir.create(dirname(save_path),showWarnings = FALSE, recursive = TRUE)
    ggsave(save_path,
           plot = pp_wrap, dpi = 300,  width = 10, height = 12)
  }


}





REPLICATION.compare_PP_scatterplot <- function(merged_PD,
                                              col1="mean.PP.UKB",
                                              col2="mean.PP.1KG",
                                              title=NULL,
                                              filter_str="!is.na(SNP)",
                                              max_labels=5,
                                              xlabel=col1,
                                              ylabel=col2,
                                              show_count=TRUE,
                                              show_plot=TRUE){
  requireNamespace("ggplot2")
  #### Scatter plot ####
  label_snps <- subset(merged_PD, eval(parse(text=col1))>.75 &  eval(parse(text=col2))>.75) %>%
    dplyr::mutate(label=paste(Locus,SNP, sep = "\n"))
  label_snps <- label_snps[1:max_labels]
  pearson <- cor.test(merged_PD[[col1]], merged_PD[[col2]], method="pearson")
  spearman <- cor.test(merged_PD[[col1]], merged_PD[[col2]], method="spearman")

  pp <- ggplot(data = subset(merged_PD, eval(parse(text = filter_str))),
               aes_string(x=col1, y=col2)) +
    geom_smooth(method = "lm") +
    geom_point(alpha=.5, aes(color=P)) +
    ggrepel::geom_label_repel(data = label_snps, alpha=.75, #box.padding = 2,
                              #nudge_x = -.5,
                              nudge_y = .5,
                              color="black",  segment.alpha = .5,
                              aes_string(x=col1, y=col2, color="P", label="label")) +
    labs(title=title,
         subtitle = paste("Spearman: rho =",round(spearman$estimate,2),
                          if(spearman$p.value<2.2e-16) "p < 2.2e-16" else paste("p =",round(spearman$p.value,2)),
                          if(show_count)paste("(n =",dplyr::n_distinct(merged_PD$SNP),"SNPs)") else NULL
                          ),
         color="GWAS P",
         x=xlabel, y=ylabel) +
    scale_color_viridis_c() +
    xlim(c(0,1)) + ylim(c(0,1)) +
    theme_bw()
  if(show_plot) print(pp)
  return(pp)
}




REPLICATION.compare_PP_pairedplot <- function(merged_PD,
                                               col1="mean.PP.UKB",
                                               col2="mean.PP.1KG",
                                               title=NULL,
                                               filter_str="!is.na(SNP)",
                                               variable_key=NULL,
                                               line_alpha=.1,
                                               xlabel=NULL,
                                               ylabel=NULL,
                                               show_plot=TRUE){
  requireNamespace("ggplot2")
  melt_PD <- data.table::melt.data.table(data = subset(merged_PD, eval(parse(text = filter_str))),
                                         id.vars = c("SNP","Locus","P"),
                                         measure.vars = c(col1,col2),
                                         value.name = "PP")
  # ifelse(variable=="FINEMAP.PP","PP (Schilder et al.)","PP (Grenn et al.)")
  if(is.null(variable_key)) variable_key <- setNames(unique(melt_PD$variable),unique(melt_PD$variable))
  melt_PD <- melt_PD %>%
    dplyr::mutate(source = variable_key[variable],
                  paired=1)

  #### Paired boxplot ####
  gp <- ggplot(data =melt_PD,
               aes(x=source, y=PP, fill=source)) +
    # geom_boxplot() +
    geom_violin(alpha=.5) +
    geom_line(aes(group = SNP),
              alpha = line_alpha, colour = "darkgrey", position = position_dodge(0.2) ) +
    # geom_jitter(width=0.15, height = 0, alpha=.1) +
    geom_point(aes(fill=source, group=SNP), position = position_dodge(0.2),
               alpha=.5) +
    labs(title=title,
         subtitle = paste0("PP ≥ ",PP_thresh)) +
    ylim(c(0,1)) +
    theme_bw() +
    theme(legend.position = "none")
  if(!is.null(xlabel)) gp <- gp + labs(x=xlabel)
  if(!is.null(ylabel)) gp <- gp + labs(y=ylabel)
  if(show_plot) print(gp)
  return(gp)
}




REPLICATION.compare_PP_heatmap <- function(merged_PD,
                                           PP_cols= grep("\\.PP\\.UKB$|\\.PP\\.1KG$", colnames(merged_PD), value = TRUE)){
  requireNamespace("heatmaply")
  sources <- lapply(PP_cols, function(x) rev(strsplit(x,"\\.")[[1]])[1]) %>% unlist()
  corr_dat <- data.frame(merged_PD, row.names = merged_PD$SNP)[,PP_cols]
  corr_dat[is.na(corr_dat)] <- 0
  # corr_mat <- data.frame(corrplot::cor.mtest(corr_dat)$p,
  #                        row.names = PP_cols) %>%
  #   `colnames<-`(PP_cols)
  corr_mat <- cor(corr_dat,
                  use="na.or.complete",
                  method = "spearman")

  heat <- heatmaply::heatmaply(x = corr_mat,
                    col = viridis::magma(10), #rev(RColorBrewer::brewer.pal(10,"magma")),
                    # plot_method="ggplot",
                    # return_ppxpy = TRUE,
                    # file = file.path("Data/GWAS/Nalls23andMe_2019/_genome_wide","LD_comparison/UKB_vs_1KG.PP.heatmap.png"),
                    key.title = "Spearman's rho",
                    grid_gap = .75,

                    limits = c(-1,1),
                    row_side_colors = data.frame("LD panel"=sources, check.names  = FALSE),
                    col_side_colors = data.frame("LD panel"=sources, check.names  = FALSE)
  )
  print(heat)
  return(heat)
}




