



#' Merge validation assays into a single plot
#'
#' @family VALIDATION
#' @examples
#' plt.ALL <- VALIDATION.super_plot(base_url="/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide")
VALIDATION.super_plot <- function(base_url="/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide",
                                  height=10, width=12,
                                  show_plot=T,
                                  save_plot=F){
  snp.groups_list <- snp_group_filters()
  expanded_groups <-  grep("Support",names(snp_group_filters()),value = T, invert = T)

  #### S-LDSC h2 ####
  res.h2 <- data.table::fread(file.path(base_url,"PolyFun/Nalls23andMe_2019.h2_enrich.snp_groups.csv.gz"))
  plt.h2 <- POLYFUN.h2_enrichment_SNPgroups_plot(RES = res.h2,
                                                 # snp_groups = c("GWAS lead","UCS","Consensus"),
                                                 # comparisons_filter = NULL,
                                                 snp_groups = expanded_groups,
                                                 save_path = file.path(base_url,"PolyFun/Nalls23andMe_2019.h2_enrich.snp_groups_expanded.png"),
                                                 show_plot = T)# +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())

  # coin::independence_test(h2.enrichment ~ SNP.Group,
  #                         data = subset(res.h2, SNP.Group %in% c("GWAS lead","UCS"))%>%
  #                           dplyr::mutate(SNP.Group=factor(SNP.Group),
  #                                         Locus=factor(Locus)))


  #### IMPACT ####
  res.IMPACT <- data.table::fread(file.path(base_url,"IMPACT/TOP_IMPACT_all.csv.gz"))
  ## binarize
  # res.IMPACT <- dplyr::mutate(res.IMPACT, mean_IMPACT=ifelse(mean_IMPACT>=.95,1,0))
  plt.IMPACT <- IMPACT.snp_group_boxplot(TOP_IMPACT_all = res.IMPACT,
                                         # snp_groups = c("GWAS lead","UCS","Consensus",paste0("Support==",0:4)),
                                         snp_groups = expanded_groups,
                                         title = "IMPACT",
                                         ylabel = "mean IMPACT score",
                                         save_path = file.path(base_url,"IMPACT/Nalls23andMe_2019.IMPACT.snp_groups_expanded.png"),
                                         # comparisons_filter = NULL,
                                         show_plot = T) #+

    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())


  #### Deep learning ####
  res.DL <- data.table::fread(file.path(base_url,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.snp_groups.mean.csv.gz"))
  plt.DL <- DEEPLEARNING.plot(annot.melt = res.DL,
                              # snp_groups = c("GWAS lead", "UCS","Consensus"),
                              snp_groups = expanded_groups,
                              facet_formula = "Assay ~ Model + Tissue",
                              # comparisons_filter = NULL,
                              model.metric = "MEAN",
                              save_path = file.path(base_url,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.MEAN_expanded.png"),
                              show_plot = T)
  plt.DL_zoom <- plt.DL + ylim(c(0,0.006))
  ggsave(file.path(base_url,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.MEAN_expanded.png"),
         plt.DL_zoom, dpi=400)

  #### SURE ####
  res.SURE <- data.table::fread(file.path(base_url,"SURE/Nalls23andMe_2019.SURE.snp_groups.mean.csv.gz"))
  plt.SURE <- SURE.plot(sure.melt = res.SURE,
                        # snp_groups = c("GWAS lead","UCS","Consensus"),
                        snp_groups = expanded_groups,
                        # comparisons_filter = NULL,
                        save_path = file.path(base_url,"SURE/Nalls23andMe_2019.SURE.snp_groups_expanded.mean.png"),
                        show_plot = T,
                        width=8)
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




#' Check whether Support level correlates with some variable
#'
#' @examples
#' # S-LDSC h2
#' h2_stats <- data.table::fread("/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/Nalls23andMe_2019.grouped_snpvar_stats.csv.gz")
VALIDATION.compare_support <- function(){
  plot_dat <- h2_stats
  variable <- "SNPVAR_min"

  plot_dat$log <- log1p( plot_dat[[variable]])
  comparisons <- utils::combn(x = as.character(unique(plot_dat$Support)),
                              m=2,
                              simplify = F) %>% purrr::compact()
  method="wilcox.test"
  ggpubr::ggviolin(data = plot_dat,
                   x="Support", y="log",
                   fill="Support", alpha=.5,
                   add="boxplot", add.params = list(alpha=.1)) +
      ggpubr::stat_compare_means(method=method, comparisons = comparisons,
                                 label = "p.signif") +
    geom_jitter(alpha=.1, height=0)
}










VALIDATION.compare_singleton_results <- function(singleton_url="/pd-omics/data/Singleton_FINEMAP_PD/pd_meta5_sum_stats_fm_results.csv.gz"){
  # Originally from: https://storage.googleapis.com/nihnialng-share-f23bef932184/pd_meta5_sum_stats_fm_results.csv.gz
  singleton <- data.table::fread(singleton_url)
  # Merging by SNP yields waaaaaay more hits than CHR/POS, even after liftover.
  merged_PD <- data.table::merge.data.table(merged_dat,
                                            singleton,
                                            by="SNP" )
  # single <- singleton %>%
  #   tidyr::separate(col = "rsid", into=c("CHR","POS"), sep=":", remove=F) %>%
  #   dplyr::mutate(CHR=as.integer(gsub("chr","",CHR)),
  #                 POS=as.integer(POS)) %>%
  #   LIFTOVER(dat=singleton, return_as_granges=F) %>%
  #   data.table::data.table()
  # merged_PD <- data.table::merge.data.table(merged_dat,
  #                                           single,
  #                                           by.x=c("CHR","POS"),
  #                                           by.y =c("CHR","POS"))
  # merged_PD[is.na(merged_PD)] <- 0


  cor.test(merged_PD$FINEMAP.PP, merged_PD$prob, method="pearson")
  cor.test(merged_PD$FINEMAP.PP, merged_PD$prob, method="spearman")
}




#### MUST USE FULL DATASET (all loci merged) to get proper null distribution
# Alternaitvely, could use a uniform distribution (but this is an assumption)
## But for IMPACT, unif dist may atually be higher impact than reality

VALIDATION.permute <- function(){
  # permut
  coin::independence_test(Length ~ SNP_Group,
                          data = data)


}


VALIDATION.bootstrap <- function(metric_df,
                               metric,
                               validation_method=NULL,
                               synthesize_random=F,
                               snp_groups=c("Random","GWAS lead","UCS","Consensus"),
                               save_path=F,
                               nThread=4){
  # ANNOT_MELT <- data.table::fread("/sc/arion/projects/pd-omics/data/IMPACT/IMPACT707/Annotations/IMPACT_overlap.csv.gz", nThread=4)
  # res.IMPACT <- data.table::fread(file.path(base_url,"IMPACT/TOP_IMPACT_all.csv.gz"))

  # metric_df <- sampling_df <- ANNOT
  # metric <- grep("Basenji|DeepSEA", colnames(metric_df), value = T)[1] # "mean_IMPACT"
  sampling_df <- metric_df
  snp_filters <- snp_group_filters()
  # snp_filters <- setNames(paste0("SNP_group=='",snp_groups,"'"), snp_groups)

  # Bootstrap with replacement
  RES_GROUPS <- lapply(snp_groups, function(snp_group){
    message(snp_group)
    parallel::mclapply(1:1000, function(i,
                                       .snp_group=snp_group,
                                       .synthesize_random=synthesize_random){
        message(i)
        snp_filt <- snp_filters[[.snp_group]]
        bg_size = 3
        #### Random sample ####
        random_snps <-  metric_df %>%
          dplyr::sample_n(size=bg_size)%>%
          dplyr::mutate(SNP_Group="Random")
        #### optional: use random noise if you dont have the fulll sample
        if(.synthesize_random){
          # Synthesize a uniform distribution
          random_snps[[metric]] <- runif(n = bg_size,
                                         min=min(metric_df[[metric]], na.rm = T),
                                         max=max(metric_df[[metric]], na.rm = T))
        } else {
          # Sample from the real distribution
          random_snps[[metric]] <- sample(metric_df[[metric]],
                                          size = bg_size)
        }

        #### Target sample ####
        fg_size=3
        if(.snp_group=="Random"){
          .snp_group <- "Random_target"
          target_snps <-  metric_df %>%
            dplyr::sample_n(size=fg_size)%>%
            dplyr::mutate(SNP_Group=.snp_group)
          if(.synthesize_random){
            target_snps[[metric]] <- runif(n = fg_size,
                                           min=min(metric_df[[metric]], na.rm = T),
                                           max=max(metric_df[[metric]], na.rm = T))
          } else {
            target_snps[[metric]] <- sample(metric_df[[metric]],
                                            size = fg_size)
          }

        } else {
          target_snps <- subset(metric_df, eval(parse(text=snp_filt))) %>%
            dplyr::sample_n(size = fg_size) %>%
            dplyr::mutate(SNP_Group=.snp_group)
        }
        # Merge
        dat <- rbind(random_snps, target_snps) %>%
          dplyr::mutate(SNP_Group=as.factor(SNP_Group),
                        Locus=as.factor(Locus))
        # res <- coin::wilcox_test(data=dat,
        #                          IMPACT_score ~ SNP_Group|Locus)
        # data.frame(Z= coin::statistic(res),
        #            p= coin::pvalue(res))
        res <- ggpubr::compare_means(data = dat,
                                     formula = as.formula(paste(metric,"~ SNP_Group")),
                                     # group.by = "Locus",
                                     method="wilcox.test")
        return(res)
    }, mc.cores = nThread) %>%  data.table::rbindlist() %>%
        dplyr::mutate(SNP_Group=snp_group)
  })  %>% data.table::rbindlist()

  if(save_path!=F){
    # validation_method = "Dey_DeepLearning"
    data.table::fwrite(RES_GROUPS, file.path(save_path,
                                             validation_method,
                                             paste("Nalls23andMe_2019",validation_method,"permutations.csv.gz",sep=".")))
  }
  return(RES_GROUPS)
}




VALIDATION.bootstrap_multimetric <- function(metric_df,
                                             metric_names){
  # save_path <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
  # metric_df <- ANNOT; validation_method = "Dey_DeepLearning";
  # metric_names <- grep("Basenji.*_MEAN|DeepSEA.*_MEAN", colnames(metric_df), value = T)
  RES_GROUPS <- lapply(metric_names, function(metric){
    message(metric)
    VALIDATION.bootstrap(metric_df, metric=metric, save_path=F)
  }) %>% data.table::rbindlist(fill=T)

  if(save_path!=F){
    data.table::fwrite(RES_GROUPS, file.path(save_path,
                                             validation_method,
                                             paste("Nalls23andMe_2019",validation_method,"permutations.csv.gz",sep=".")))
  }
  return(RES_GROUPS)
}





#' Test whether permutation p-value distributions are different
#'
#' @family VALIDATION
#' @examples
#' root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"
#' permute.IMPACT <-  data.table::fread(file.path(root,"_genome_wide/IMPACT/Nalls23andMe_2019.IMPACT.permutations.csv.gz"))
#' res <- VALIDATION.permute_compare_results(permute_res=permute.IMPACT)
VALIDATION.compare_distributions <- function(permute_res){
  # GLM
  fit <- stats::glm(data = permute_res,
                    # family = stats::poisson # This would actually remove the signfiicance from the signal
                    formula = p ~ SNP_Group)
  print(summary(fit))
  ## Extract p-value
  coefs <-  data.frame(coef(summary(fit))) %>%
    `colnames<-`(c("Estimate","StdErr","z","p")) %>%
    tibble::rownames_to_column("SNP_Group") %>%
    dplyr::mutate(SNP_Group=gsub("SNP_Group","",SNP_Group),
                  signif = pvalues_to_symbols(p))
  ## Extract lambda w/ MASS
  lambda <- MASS::fitdistr(x = subset(permute_res, !is.na(p))$p,
                           densfun = "Poisson")
  coefs$lambda <- as.numeric(lambda[[1]])
  return(coefs)

  # Coin
  # RES_GROUPS$SNP_Group <- as.factor(RES_GROUPS$SNP_Group)
  # coin::independence_test(p ~ SNP_Group,
  #                   data=RES_GROUPS)
  # Now test the p-val distributions

  ## None of these account for the zero inflted nature of the p-value distributions
  # res <- ggpubr::compare_means(data = RES_GROUPS,
  #                              formula = p ~ SNP_Group,
  #                              method="wilcox.test")
  # res <- coin::independence_test(data= RES_GROUPS,
  #                                p ~ SNP_Group,
  #                                distribution="approximate")
  ## Negative-binomial distributions are intended for count data
  # fit <- MASS::glm.nb(data = RES_GROUPS,
  #              formula = p ~ SNP_Group)
  # summary(fit)

  ## Using a Poisson distribution, and estimating lambda with a Generalized Linear Model (GLM),
  ## will fit a model that much better appromxiates the p-value distributions.
}




#' Plot permutation test results
#'
#' @family VALIDATION
#' @examples
#' root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
#'
#'  ## IMPACT
#' permute_res <-  data.table::fread(file.path(root,"IMPACT/Nalls23andMe_2019.IMPACT.permutations.csv.gz"))
#' ## Dey_DeepLearning
#' permute_res <-  data.table::fread(file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.permutations.csv.gz"))
#' gp <- VALIDATION.permute_plot(permute_res=permute_res)
VALIDATION.permute_plot <- function(permute_res,
                                    validation_method="Dey_DeepLearning",
                                    save_plot=F,
                                    show_plot=T){
  library(ggplot2)
  plot_dat <- permute_res %>%
    dplyr::mutate(SNP_Group = factor(SNP_Group, levels = c("Random","GWAS lead","UCS","Consensus"), ordered = T),
                  Metric=.y.)
  colorDict <- setNames(c("grey","red","green2","goldenrod2"),levels(plot_dat$SNP_Group))

  # Conduct GLM on pval distributions
  permute_res$SNP_Group <- factor(permute_res$SNP_Group, levels = c("Random","GWAS lead","UCS","Consensus"))
  if(validation_method=="Dey_DeepLearning"){
    glm_res <- lapply(unique(permute_res$.y.), function(metric){
      print(metric)
      VALIDATION.compare_distributions(permute_res=subset(permute_res, .y.==metric)) %>%
        dplyr::mutate(Metric=metric)
    }) %>% data.table::rbindlist()
    # Postprocess
    glm_res <- glm_res %>%
      subset(SNP_Group!="(Intercept)") %>%
      tidyr::separate(Metric, sep="_", into=c("Model","Tissue","Assay"), extra="drop", remove=F) %>%
      dplyr::mutate(SNP_Group = factor(SNP_Group, levels = c("Random","GWAS lead","UCS","Consensus"), ordered = T))
    plot_dat <- plot_dat %>% tidyr::separate(Metric, sep="_", into=c("Model","Tissue","Assay"), extra="drop", remove=F)

    # ggplot(data = glm_res, aes(x=SNP_Group, y=-log10(p), fill=SNP_Group)) +
      # geom_bar(stat="identity",alpha=.5) +
      ggpubr::ggviolin(data = glm_res,
                       x="SNP_Group", y="p", fill="SNP_Group",
                       alpha=.5) +
      geom_boxplot(alpha=.5) +
      yscale("-log1p", .format = T) +
      geom_jitter(height = 0, alpha=.1) +
      scale_fill_manual(values = colorDict) +
      geom_hline(yintercept = -log10(0.05), alpha=.5, linetype=2) +
      facet_grid(facets = . ~   Tissue,
                 scales="free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle=45, hjust=1))

  } else {
    glm_res <- VALIDATION.compare_distributions(permute_res=permute_res)
  }

  # Plot
  gp <- ggplot(data = plot_dat, aes(x= p, fill=SNP_Group, color=SNP_Group,
                                    linetype=SNP_Group)) +
    # geom_histogram(alpha=.5,  position="identity", bins=100) +
    geom_density(position = "identity", adjust=1, alpha=.5,  ) +
    scale_fill_manual(values = colorDict) +
    scale_color_manual(values = colorDict) +
    labs(title=paste("Bootstrapped tests:",validation_method,"score"), x="bootstrapped p-values") +
    theme_bw()

  # Get density peaks
  b <- ggplot_build(gp)
  density_peaks <- b$data[[1]] %>% dplyr::group_by(fill) %>% dplyr::top_n(n = 1, wt = y) %>%
    dplyr::mutate(SNP_Group =  setNames(names(colorDict), unname(colorDict))[[fill]]) %>%
    data.table::data.table() %>%
    data.table::merge.data.table(glm_res, by="SNP_Group")
  # Annotate plot
  gp_lab <- gp  +
    ggrepel::geom_label_repel(data=density_peaks,
                              aes(x=x, y=y, label=paste("p =",format(as.numeric(p), digits = 3), signif)),
                              show.legend = F, point.padding = 1,
                              alpha=.75, segment.alpha = .5,
                              color="black")
  # facet_wrap(facets = Tissue ~ .)

  if(show_plot) print(gp_lab)
  if(save_plot!=F){
    ggsave(save_plot, gp_lab, dpi=400, width=6)
  }
  return(gp_lab)
}





