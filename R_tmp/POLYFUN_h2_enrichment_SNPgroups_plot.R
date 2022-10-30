#' Plot heritability (h2) enrichment
#'
#' @family polyfun
#' @examples
#' \dontrun{
#' root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping"
#' merged_dat <- merge_finemapping_results(dataset = file.path("Data/GWAS/Nalls23andMe_2019"), LD_reference = "UKB", minimum_support = 0)
#' RES <- POLYFUN_h2_enrichment_SNPgroups(merged_dat=merged_dat, out.path="/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output")
#'
#' plot.h2 <- POLYFUN_h2_enrichment_SNPgroups_plot(RES = RES, show_plot = TRUE)
#' }
POLYFUN_h2_enrichment_SNPgroups_plot <- function(RES,
                                                 snp_groups=c("GWAS lead","UCS","Consensus (-PolyFun)","Consensus"),
                                                 comparisons_filter=function(x){if("Consensus" %in% x) return(x)},
                                                 method="wilcox.test",
                                                 remove_outliers=TRUE,
                                                 title= "S-LDSC heritability enrichment",
                                                 xlabel="SNP group",
                                                 ylabel=bquote(~'h'^2~'enrichment'),
                                                 show_xtext=TRUE,
                                                 show_padj=TRUE,
                                                 show_signif=TRUE,
                                                 vjust_signif=0.5,
                                                 show_plot=TRUE,
                                                 save_path=FALSE){
    requireNamespace("ggplot2")
    requireNamespace("ggpubr")
    requireNamespace("purrr")
    colorDict <- echodata::snp_group_colorDict()
    if(is.null(snp_groups)) snp_groups <- names(colorDict)
    if(remove_outliers){
        # https://www.r-bloggers.com/2020/01/how-to-remove-outliers-in-r/
        outliers <- boxplot(RES$h2.enrichment, plot=FALSE)$out
        RES <- RES[-which(RES$h2.enrichment %in% outliers),]
    }
    plot_dat <- subset(RES, SNP_group %in% snp_groups) |>
        dplyr::mutate(SNP_group=factor(SNP_group, levels=names(colorDict), ordered = TRUE))
    snp.groups <- unique(plot_dat$SNP_group)
    comparisons <- utils::combn(x = as.character(snp.groups),
                                m=2,
                                FUN = comparisons_filter,
                                simplify  = FALSE) |> purrr::compact()
    pb <- ggplot(data = plot_dat, aes(x=SNP_group, y=h2.enrichment, fill=SNP_group)) +
        geom_jitter(alpha=.1,width = .25, show.legend = FALSE, shape=16, height=0) +
        geom_violin(alpha=.6, show.legend  = FALSE) +
        geom_boxplot(alpha=.6, color="black", show.legend  = FALSE) +
        geom_hline(yintercept = log10(1), linetype=2, alpha=.5) +
        # scale_fill_manual(values =  c("red","green3","goldenrod3")) +
        labs(y=ylabel, x=xlabel,
             title=title) +
        theme(legend.position = "none")  +
        scale_fill_manual(values = colorDict) +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45, hjust=1))
    if(show_padj){
        pb <- pb + ggpubr::stat_compare_means(method = method,
                                              comparisons = comparisons,
                                              label = "p.adj", vjust=2, size=3)
    }
    if(show_signif){
        pb <- pb + ggpubr::stat_compare_means(method = method,
                                              comparisons = comparisons,
                                              label = "p.signif", size=3, vjust =  vjust_signif)
    }
    
    if(!show_xtext){
        pb <- pb + theme(axis.text.x = element_blank(),
                         axis.title.x = element_blank())
    }
    
    if(show_plot) print(pb)
    if(save_path!=FALSE){
        # save_path <- '~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/snp_group.h2_enrichment.png'
        ggsave(save_path,
               plot = pb, dpi = 300, height=9, width=11)
    }
    return(pb)
}

