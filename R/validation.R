



#' Merge validation assays into a single plot
#'
#' @family VALIDATION
#' @examples
#'
VALIDATION.super_plot <- function(base_url="/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide",
                                  height=10, width=12,
                                  show_plot=T,
                                  save_plot=F){
  snp.groups_list <- snp_group_filters()

  # S-LDSC h2
  res.h2 <- data.table::fread(file.path(base_url,"PolyFun/Nalls23andMe_2019.h2_enrich.snp_groups.csv.gz"))
  plt.h2 <- POLYFUN.h2_enrichment_SNPgroups_plot(RES = res.h2,
                                                 snp_groups = c("GWAS lead","UCS","Consensus"),
                                                 comparisons_filter = NULL,
                                                 # snp_groups = names(snp_group_filters()),
                                                 show_plot = F) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())
  # IMPACT
  res.IMPACT <- data.table::fread(file.path(base_url,"IMPACT/TOP_IMPACT_all.csv.gz"))
  plt.IMPACT <- IMPACT.snp_group_boxplot(TOP_IMPACT_all = res.IMPACT,
                                         snp_groups = c("GWAS lead","UCS","Consensus"),
                                         # snp_groups = names(snp_group_filters()),
                                         save_path = F,
                                         title = "IMPACT",
                                         ylabel = "mean IMPACT score",
                                         comparisons_filter = NULL,
                                         show_plot = F) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())
  # Deep learning
  res.DL <- data.table::fread(file.path(base_url,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.snp_groups.mean.csv.gz"))
  plt.DL <- DEEPLEARNING.plot(annot.melt = res.DL,
                              snp_groups = c("GWAS lead", "UCS","Consensus"),
                              comparisons_filter = NULL,
                              model.metric = "MEAN",
                              show_plot = F)
  # SURE
  res.SURE <- data.table::fread(file.path(base_url,"SURE/Nalls23andMe_2019.SURE.snp_groups.mean.csv.gz"))
  plt.SURE <- SURE.plot(sure.melt = sure.melt,
                        snp_groups = c("GWAS lead","UCS","Consensus"),
                        comparisons_filter = NULL,
                        show_plot = F
                        # snp_groups = names(snp_group_filters())
                        )
   # MERGE PLOTS
  library(patchwork)
  plt.ALL <-  ( ( (plt.h2 + plt.IMPACT) / plt.SURE) |
                (plt.DL) ) +
    patchwork::plot_layout(widths = c(.5,1)) +
    patchwork::plot_annotation(tag_levels = letters)

  if(show_plot) print(plt.ALL)
  if(save_plot!=F){
    ggsave(save_plot, plt.ALL, dpi=300,
           height=height, width=width)
  }
  return(plt.ALL)
}


