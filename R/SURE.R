
# ^^^^^^^^^^^^^^ SURE ^^^^^^^^^^^^^^
# Data download source:
# https://sure.nki.nl

# Reference:::
# Arensbergen, Joris van, Ludo Pagie, Vincent D. FitzPatrick, Marcel de Haas, Marijke P. Baltissen, Federico Comoglio, Robin H. van der Weide, et al. “High-Throughput Identification of Human SNPs Affecting Regulatory Element Activity.” Nature Genetics 51, no. July (2019). https://doi.org/10.1038/s41588-019-0455-2.
# ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^

#' Download \emph{SuRE} annotations
#'
#' @source
#' \href{https://www.nature.com/articles/s41588-019-0455-2}{Publication}
#' \href{https://github.com/vansteensellab/SuRE-SNV-code}{GitHub}
#' \href{https://osf.io/w5bzq/wiki/home/?view}{Full SuRE data}
#' \href{https://sure.nki.nl}{SNP-SuRE data browse}
SURE.download_annotations <- function(URL="https://osf.io/vxfk3/download",
                                      output_dir=".",
                                      nThread=4,
                                      v=verbose){
  echolocatoR::downloader(input_url = URL,
                          output_path = output_dir,
                          nThread = nThread)
  out_file <- file.path(output_dir, "SuRE_SNP_table_LP190708.txt.gz")
  printer("SURE:: Data downloaded ==>",out_file,v=verbose)
  return(output_path)
}




#' Merge with fine-mapping and \emph{SuRE} results
#'
#' @family SURE
#' @source
#' \href{https://www.nature.com/articles/s41588-019-0455-2}{Publication}
#' \href{https://github.com/vansteensellab/SuRE-SNV-code}{GitHub}
#' \href{https://osf.io/w5bzq/wiki/home/?view}{Full SuRE data}
#' \href{https://sure.nki.nl}{SNP-SuRE data browse}
#' @examples
#' sure <- data.table::fread("/sc/arion/projects/pd-omics/data/MPRA/SURE/SuRE_SNP_table_LP190708.txt.gz", nThread=4)
SURE.merge <- function(merged_DT,
                       sure){
  sure <- dplyr::mutate(sure,
                        chr = as.integer(gsub("chr","",chr)),
                        SNPabspos = as.integer(SNPabspos)
                        )
  sure_DT <- data.table::merge.data.table(merged_DT,
                                           sure,
                                           all.x = T,
                                           by.x = "SNP",
                                           by.y = "SNP_ID"
                                           # by.x = c("CHR","POS"),
                                           # by.y = c("chr","SNPabspos")
                                           )
  return(sure_DT)
}





#' Aggregate \emph{SuRE} p-values for each SNP group
#'
#' @family SURE
#' @source
#' \href{https://www.nature.com/articles/s41588-019-0455-2}{Publication}
#' \href{https://github.com/vansteensellab/SuRE-SNV-code}{GitHub}
#' \href{https://osf.io/w5bzq/wiki/home/?view}{Full SuRE data}
#' \href{https://sure.nki.nl}{SNP-SuRE data browse}
#' @examples
#' root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
#'
#' sure_DT <- data.table::fread("/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/SURE/Nalls23andMe_2019.SURE.csv.gz", nThread=4)
#' sure_DT <- find_consensus_SNPs_no_PolyFun(sure_DT)
#' sure.melt <- SURE.melt_snp_groups(sure_DT)
#' data.table::fwrite(sure.melt, file.path(root,"SURE/Nalls23andMe_2019.SURE.snp_groups.mean.csv.gz"))
SURE.melt_snp_groups <- function(sure_DT,
                                 grouping_vars=c("Locus"),
                                 metric_str="mean",
                                 replace_NA=NA){
  measure.vars <- c("k562.wilcox.p.value","hepg2.wilcox.p.value"
                    # "k562.wilcox.p.value.random","hepg2.wilcox.p.value.random",
                    # "ref.element.count","alt.element.count","k562.ref.mean", "k562.alt.mean",
                    #  "hepg2.ref.mean","hepg2.alt.mean"
                    )
  snp.groups_list <- snp_group_filters()
  column_substr <- paste(measure.vars, collapse="|")
  metric <- get(metric_str)
  sampling_df <- sure_DT
  sure.melt <- sure_DT %>%
    dplyr::group_by(.dots=grouping_vars) %>%
    dplyr::summarise_at(.vars = vars(measure.vars),
                        .funs = list("Random"= ~ metric(tidyr::replace_na(sample(.x, size=3, replace = T),replace_NA), na.rm = T),
                                     "All"= ~ metric(tidyr::replace_na(.x,replace_NA), na.rm = T),
                                     "GWAS nom. sig."= ~ metric(tidyr::replace_na(.x[P<.05],replace_NA), na.rm = T),
                                     "GWAS sig."= ~ metric(tidyr::replace_na(.x[P<5e-8],replace_NA), na.rm = T),
                                     "GWAS lead"= ~ metric(tidyr::replace_na(.x[leadSNP],replace_NA), na.rm = T),
                                     "ABF CS"= ~ metric(tidyr::replace_na(.x[ABF.CS>0],replace_NA), na.rm = T),
                                     "SUSIE CS"= ~ metric(tidyr::replace_na(.x[SUSIE.CS>0],replace_NA), na.rm = T),
                                     "POLYFUN-SUSIE CS"= ~ metric(tidyr::replace_na(.x[POLYFUN_SUSIE.CS>0],replace_NA), na.rm = T),
                                     "FINEMAP CS"= ~ metric(tidyr::replace_na(.x[FINEMAP.CS>0],replace_NA), na.rm = T),
                                     "UCS (-PolyFun)"= ~ metric(tidyr::replace_na(.x[Support_noPF>0],replace_NA), na.rm = T),
                                     "UCS"= ~ metric(tidyr::replace_na(.x[Support>0],replace_NA), na.rm = T),
                                     "Support==0"= ~ metric(tidyr::replace_na(.x[Support==0],replace_NA), na.rm = T),
                                     "Support==1"= ~ metric(tidyr::replace_na(.x[Support==1],replace_NA), na.rm = T),
                                     "Support==2"= ~ metric(tidyr::replace_na(.x[Support==2],replace_NA), na.rm = T),
                                     "Support==3"= ~ metric(tidyr::replace_na(.x[Support==3],replace_NA), na.rm = T),
                                     "Support==4"= ~ metric(tidyr::replace_na(.x[Support==4],replace_NA), na.rm = T),
                                     "Consensus (-PolyFun)"= ~ metric(tidyr::replace_na(.x[Consensus_SNP_noPF],replace_NA), na.rm = T),
                                     "Consensus"= ~ metric(tidyr::replace_na(.x[Consensus_SNP],replace_NA), na.rm = T)
                        ),
    ) %>%
    data.table::data.table() %>%
    data.table::melt.data.table(id.vars = grouping_vars, variable.name = "Annotation") %>%
    dplyr::mutate(Annotation =gsub("p\\.value","pval",Annotation)) %>%
    dplyr::mutate(Annotation =gsub("\\.","_",Annotation)) %>%
    tidyr::separate(col = "Annotation", sep="_", into=c("Cell_type","Test","Metric","SNP_group"), remove=F) %>%
    dplyr::mutate(Annotation = DescTools::StrTrim(Annotation, "_"),
                  SNP_group = factor(SNP_group, levels = names(snp.groups_list), ordered = T),
                  neg.log.value = log1p(value))
  return(sure.melt)
}





#' Plot \emph{SuRE} results across SNP groups
#'
#' @family SURE
#' @source
#' \href{https://www.nature.com/articles/s41588-019-0455-2}{Publication}
#' \href{https://github.com/vansteensellab/SuRE-SNV-code}{GitHub}
#' \href{https://osf.io/w5bzq/wiki/home/?view}{Full SuRE data}
#' \href{https://sure.nki.nl}{SNP-SuRE data browse}
#' @examples
#' sure.melt <- data.table::fread("/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/SURE/Nalls23andMe_2019.SURE.snp_groups.mean.csv.gz")
SURE.plot <- function(sure.melt,
                      snp_groups=c("GWAS lead","UCS","Consensus (-POLYFUN)","Consensus"),
                      comparisons_filter=function(x){if("Consensus" %in% x) return(x)},
                      title="SuRE MPRA",
                      xlabel="SNP Group",
                      facet_formula=". ~ Cell_type",
                      show_padj=T,
                      show_signif=T,
                      vjust_signif=0.5,
                      show_plot=T,
                      save_path=F,
                      height=5,
                      width=5){
  colorDict <- snp_group_colorDict()
  plot_dat <-  subset(sure.melt, SNP_group %in% snp_groups) %>%
    dplyr::mutate(SNP_group=factor(SNP_group, levels=names(colorDict), ordered = T))
  snp.groups <- unique(plot_dat$SNP_group)
  comparisons <- utils::combn(x = as.character(snp.groups),
                              m=2,
                              FUN = comparisons_filter,
                              simplify = F) %>% purrr::compact()
  method="wilcox.test"
  pb <-  ggplot(data = plot_dat, aes(x=SNP_group, y=-log1p(value), fill=SNP_group)) +
    geom_jitter(alpha=.1,width = .25, show.legend = F, shape=16, height=0) +
    geom_violin(alpha=.6, show.legend = F) +
    geom_boxplot(alpha=.6, color="black", show.legend = F) +
    geom_jitter(alpha=.1, width = .3, height=0) + #aes(color=SNP.Group)) +
    facet_grid(facets = as.formula(facet_formula) ,
               scales = "free") +
    labs(title=title, x=xlabel) +
    scale_fill_manual(values = colorDict) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          legend.position = "none",
          strip.background = element_rect(fill = "grey20"),
          strip.text= element_text(color = "white"))
  if(show_padj){
    pb <- pb + ggpubr::stat_compare_means(method = method,
                                          comparisons = comparisons,
                                          label = "p.adj", size=3,  vjust=2)
  }
  if(show_signif){
    pb <- pb + ggpubr::stat_compare_means(method = method,
                                          comparisons = comparisons,
                                          label = "p.signif", size=3, vjust = vjust_signif)
  }

  if(show_plot) print(pb)
  # save_path <- file.path("/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/SURE/",
  #                        paste("Nalls23andMe_2019","SURE","png",sep="."))
  if(save_path!=F){
    ggsave(save_path,
           pb, dpi = 400, height=height, width=width)
  }
  return(pb)
}





