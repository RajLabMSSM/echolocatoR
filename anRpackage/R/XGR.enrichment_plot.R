XGR.enrichment_plot <-
function(enrich_res,
                                title=NULL,
                                subtitle=NULL,
                                facet_formula=NULL,
                                line_formula = "y ~ exp(x)",
                                FDR_thresh=1,
                                plot_type="bar",
                                shape_var="Cell_type",
                                facet_scales="free",
                                show_plot=T,
                                save_plot=F,
                                height=5,
                                width=5){
  enrich_res <- dplyr::mutate(enrich_res,
                              SNP_Group=factor(SNP_Group, levels = unique(SNP_Group), ordered = T),
                              ## Make Random size smaller (otherwise will make everything else relatively tiny)
                              nOverlap=ifelse(SNP_Group=="Random",10,nOverlap)
                              )
 sum(enrich_res$fc==-Inf)
  colorDict <- c("Random"="grey",
                 "GWAS lead"="red",
                 "UCS"="green2",
                 "Consensus (-POLYFUN)"="goldenrod4",
                 "Consensus"="goldenrod2")
  if(plot_type=="bar"){
    gp <- ggplot(data=subset(enrich_res, FDR<=FDR_thresh),
                 aes(x=SNP_Group, y= fc, fill=SNP_Group)) +
      # geom_col(stat="identity", alpha=.5, show.legend = F) +
      geom_boxplot()+
      geom_jitter(height=0, width = 0, alpha=.1, show.legend = F) +
      scale_fill_manual(values = colorDict) +
      # ggpubr::stat_compare_means(method = method,
      #                            comparisons = comparisons,
      #                            label = "p.signif", size=3, vjust = 1.5) +
      facet_grid(facets = if(is.null(facet_formula)) facet_formula else as.formula(facet_formula),
                 scales="free_y") +
      labs(x="SNP Group", title=title, subtitle=subtitle) +
      theme_bw() +
      theme(strip.background = element_rect(fill="grey20"),
            strip.text = element_text(color="white"),
            axis.text.x = element_text(angle=45, hjust=1))
  }

  if(plot_type=="point"){
    gp <- ggplot(data=subset(enrich_res, FDR<=FDR_thresh),
                 aes(x= log1p(fc), y= -log10(FDR),
                                      size=nOverlap, color=SNP_Group, group=SNP_Group,
                                      fill=SNP_Group,
                                      shape=eval(parse(text=shape_var)))
                 ) +
      geom_smooth (alpha=0.1, size=0, span=1,
                   method = "lm", formula = line_formula) +
      stat_smooth (geom="line", alpha=0.3, size=1, #span=0.5,
                   method = "lm", formula = line_formula) +

      geom_point(alpha=.5) +
      scale_color_manual(values = colorDict) +
      scale_fill_manual(values = colorDict) +
      scale_shape_manual(values=12:(12+dplyr::n_distinct(enrich_res[[shape_var]]))) +
      geom_hline(yintercept = -log10(0.05), linetype=2, alpha=.5) +
      facet_grid(facets = if(is.null(facet_formula)) facet_formula else as.formula(facet_formula),
                 scales = facet_scales)  +
      labs(title=title, subtitle=subtitle, shape=shape_var) +
      theme_bw() +
      theme(strip.background = element_rect(fill="grey20"),
            strip.text = element_text(color="white"))
  }

  if(show_plot) print(gp)

  if(save_plot!=F){
    ggsave(save_plot, gp, dpi=400, height=height, width=width)
  }
  return(gp)
}
