



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




#### MUST USE FULL DATASET (all loci merged) to get proper null distribution
# Alternaitvely, could use a uniform distribution (but this is an assumption)
## But for IMPACT, unif dist may atually be higher impact than reality


#' Perform validation bootstrap procedure
#'
#' @family VALIDATION
#' @source
#' \href{https://www.datanovia.com/en/lessons/wilcoxon-test-in-r/}{Wilcox test tutorial}
#' @examples
#' \dontrun{
#' save_path <- root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
#'
#' #### h2 ####
#' path <- file.path(root,"PolyFun/h2_merged.csv.gz")
#' metric <- "SNPVAR"
#'
#' #### IMPACT ####
#' ## mean_IMPACT
#' # res.IMPACT <- data.table::fread(file.path(root,"IMPACT/TOP_IMPACT_all.csv.gz"))
#' ## IMPACT_score (raw)
#' path <- "/sc/arion/projects/pd-omics/data/IMPACT/IMPACT707/Annotations/IMPACT_overlap.csv.gz"
#' metric <- "IMPACT_score"
#'
#' #### Dey_DeepLearning ####
#' ###  (see VALIDATION.bootstrap_multimetric()) ####
#'
#' ## Import and Process ##
#' metric_df <- data.table::fread(path, nThread=8)
#' # metric_df <- subset(metric_df, select=c(SNP,leadSNP,ABF.Credible_Set,ABF.PP,SUSIE.Credible_Set,SUSIE.PP,POLYFUN_SUSIE.Credible_Set,POLYFUN_SUSIE.PP,FINEMAP.Credible_Set,FINEMAP.PP,Consensus_SNP,Support,Locus,IMPACT_score))
#' metric_df$Consensus_SNP_noPF <- find_consensus_SNPs(metric_df, exclude_methods = "POLYFUN_SUSIE", sort_by_support = F)$Consensus_SNP
#'
#' #### run bootstrap ####
#' boot_res <- VALIDATION.bootstrap(metric_df=metric_df, metric=metric,   nThread=8)
#' data.table::fwrite(boot_res, file.path(root,"PolyFun/Nalls23andMe_2019.h2.bootstrap.coin_wilcox_test.csv.gz"))
#' ## data.table::fwrite(boot_res, file.path(root,"IMPACT/Nalls23andMe_2019.IMPACT_score.bootstrap.coin_wilcox_test.csv.gz"))
#' }
VALIDATION.bootstrap <- function(metric_df,
                                 metric,
                                 validation_method=NULL,
                                 synthesize_random=F,
                                 snp_groups=c("Random","GWAS lead","UCS","Consensus (-POLYFUN)","Consensus"),
                                 test_method="stats_wilcox.test",
                                 save_path=F,
                                 nThread=4,
                                 verbose=T){
  sampling_df <- metric_df
  if(!"Consensus_SNP_noPF" %in% colnames(metric_df)){
    try({
      metric_df$Consensus_SNP_noPF <- find_consensus_SNPs(metric_df,
                                                          exclude_methods = "POLYFUN_SUSIE",
                                                          sort_by_support = F)$Consensus_SNP
    })
  }

  snp_filters <- snp_group_filters()
  # snp_filters <- setNamecs(paste0("SNP_group=='",snp_groups,"'"), snp_groups)

  # Bootstrap without replacement (in a given iteration)
  printer("VALIDATION:: Using",test_method,v=verbose)
  RES_GROUPS <- lapply(snp_groups, function(snp_group,
                                            .test_method=test_method){
    message(snp_group)
    parallel::mclapply(1:1000, function(i,
                                       .snp_group=snp_group,
                                       .synthesize_random=synthesize_random,
                                       test_method=.test_method){
      res <- NULL
      try({
        message(i)
          snp_filt <- snp_filters[[.snp_group]]
          bg_size = 20
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
          }

          #### Target sample ####
          fg_size=20
          if(.snp_group=="Random"){
            .snp_group <- "Random_target"
            target_snps <-  metric_df %>%
              dplyr::sample_n(size=fg_size)%>%
              dplyr::mutate(SNP_Group=.snp_group)
            if(.synthesize_random){
              target_snps[[metric]] <- runif(n = fg_size,
                                             min=min(metric_df[[metric]], na.rm = T),
                                             max=max(metric_df[[metric]], na.rm = T))
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

          #### Permutation testing ####


          #### Technically you're only supposed to use this function after doing a regular independence test,
          ### but the way I'm using it here is different (only comparing two groups at once).
          ### The output is just convenient.
          #### Builds upon `coin::independence_test`

          # package <- "coin";
          #### coin ####
          if(test_method=="coin_independence_test"){
            coin_res <- coin::independence_test(data=dat,
                                     as.formula(paste(metric,"~ SNP_Group") ) )
            res <- data.frame(stat= coin::statistic(coin_res),
                              p= coin::pvalue(coin_res))
          }
          if(test_method=="coin_wilcox_test"){
            coin_res <- coin::wilcox_test(data=dat,
                                          as.formula(paste(metric,"~ SNP_Group") ),
                                          method="exact")
            res <- data.frame(stat= coin::statistic(coin_res),
                              p= coin::pvalue(coin_res))
          }

          #### rcompanion ####
          if(test_method=="rcompanion"){
            correction_method <- "fdr"
            res <- rcompanion::pairwisePermutationTest(data = dat,
                                                       formula=as.formula(paste(metric,"~ SNP_Group")),
                                                       method = correction_method) %>%
              dplyr::mutate(Comparison = gsub(" = 0","",Comparison),
                            adjust.method=correction_method) %>%
              tidyr::separate(col = "Comparison", sep=" - ", into=c("group1","group2"))
          }
          #### ggpubr ####
          if(test_method=="ggpubr_wilcox"){
            res <- ggpubr::compare_means(data = dat,
                                         formula = as.formula(paste(metric,"~ SNP_Group")),
                                         # group.by = "Locus",
                                         method="wilcox.test")
          }
          #### stats ####
          if(test_method=="stats_wilcox.test"){
            res_wilcox <- stats::wilcox.test(as.formula(paste(metric,"~ SNP_Group")),
                                      data=dat, conf.int=T )
            res <- data.frame(stat=res_wilcox$statistic,
                              p=res_wilcox$p.value,
                              confInt_lower=res_wilcox$conf.int[1],
                              confInt_upper=res_wilcox$conf.int[2],
                              estimate=res_wilcox$estimate)
          }
          if(test_method=="stats_pairwise.wilcox.test"){
            res <- stats::pairwise.wilcox.test(x=dat[[metric]],
                                               g = dat[["SNP_Group"]])
          }
          res$Metric <- metric
          res$test_method <- test_method
        }) ## End try
        return(res)
    }, mc.cores = nThread) %>%  data.table::rbindlist(fill=T) %>%
        dplyr::mutate(SNP_Group=snp_group)
  })  %>% data.table::rbindlist(fill=T)

  if(save_path!=F){
    # validation_method = "Dey_DeepLearning"
    data.table::fwrite(RES_GROUPS, file.path(save_path,
                                             validation_method,
                                             paste("Nalls23andMe_2019",metric,"permutations.csv.gz",sep=".")))
  }
  return(RES_GROUPS)
}



#' Conduct bootstrap procedure on multiple columns
#'
#' @family VALIDATION
#' @examples
#' save_path <- root <-  "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
#' path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.csv.gz")
#'
#' # Import and process
#' metric_df <- data.table::fread(path, nThread=8)
#' metric_df$Consensus_SNP_noPF <- find_consensus_SNPs(metric_df, exclude_methods = "POLYFUN_SUSIE", sort_by_support = F)$Consensus_SNP
#' metric_names <- grep("Basenji.*_MEAN|DeepSEA.*_MEAN", colnames(metric_df), value = T)
#' validation_method = "Dey_DeepLearning"
#'
#' enrich.boot <- VALIDATION.bootstrap_multimetric(metric_df=metric_df, metric_names=metric_names, validation_method=validation_method, save_path=gsub("\\.csv\\.gz",".bootstrap.stats_wilcox.test.csv.gz",path))
VALIDATION.bootstrap_multimetric <- function(metric_df,
                                             metric_names,
                                             validation_method=NULL,
                                             test_method="stats_wilcox.test",
                                             save_path=F,
                                             verbose=T){
  RES_GROUPS <- lapply(metric_names, function(metric,
                                              .validation_method=validation_method,
                                              .test_method=test_method){
    message(metric)
    VALIDATION.bootstrap(metric_df,
                         metric=metric,
                         validation_method = .validation_method,
                         synthesize_random=F,
                         test_method=.test_method,
                         save_path=F)
  }) %>% data.table::rbindlist(fill=T)

  if(save_path!=F){
    printer("VALIDATION:: Saving results ==>",save_path,v=verbose)
    data.table::fwrite(RES_GROUPS,save_path)
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
VALIDATION.compare_distributions <- function(permute_res,
                                             formula_str="stat ~ SNP_Group"){
  # GLM
  fit <- stats::glm(data = permute_res,
                    # family = stats::poisson # This would actually remove the signfiicance from the signal
                    formula = as.formula(formula_str))
  print(summary(fit))
  ## Extract p-value
  coefs <-  data.frame(coef(summary(fit))) %>%
    `colnames<-`(c("Estimate","StdErr","z","p")) %>%
    tibble::rownames_to_column("SNP_Group") %>%
    dplyr::mutate(SNP_Group=gsub("SNP_Group","",SNP_Group),
                  signif = pvalues_to_symbols(p))
  ## Extract lambda w/ MASS
  metric <- strsplit(formula_str," ~ ")[[1]][1]
  lambda <- MASS::fitdistr(x = subset(permute_res, !is.na(eval(parse(text=metric))))[[metric]],
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




VALIDATION.aggregate_permute_res <- function(permute_res){
  permute_agg <- permute_res %>%
    dplyr::mutate(SNP_Group = factor(SNP_Group, levels = unique(permute_res$SNP_Group), ordered = T),
                  Metric=.y.) %>%
    tidyr::separate(Metric, sep="_", into=c("Model","Tissue","Assay"), extra="drop", remove=F) %>%
    dplyr::mutate(.y.=gsub(paste(unique(Assay),collapse="|"),"",.y.),
                  SNP_Group=as.character(SNP_Group)) %>%
    dplyr::select(-c(Metric,Model,Tissue,Assay))
  return(permute_agg)
}



#' Plot permutation test results
#'
#' @family VALIDATION
#' @examples
#' root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
#'
#' ## h2
#' path <-  file.path(root,"PolyFun/Nalls23andMe_2019.h2.bootstrap.coin_wilcox_test.csv.gz")
#' validation_method <- "h2_score"
#'
#' ## IMPACT
#' path <- file.path(root,"IMPACT//Nalls23andMe_2019.IMPACT_score.bootstrap.coin_wilcox_test.csv.gz")
#' validation_method <- "IMPACT"
#'
#'  ## Dey_DeepLearning
#' path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.bootstrap.stats_wilcox.test.csv.gz")
#' validation_method <- "Dey_DeepLearning"
#'
#' ## SURE MPRA
#' path <- file.path(root,"SURE/Nalls23andMe_2019.SURE.permutations.csv.gz")
#' validation_method <- "SuRE MPRA"
#'
#' # ------  Import
#' permute_res <-  data.table::fread(path)
#' ## permute_res <- dplyr::rename(permute_res, stat=Stat, p=p.value)
#'
#' # ---Plot ---
#' library(patchwork)
#' permute_sub <- subset(permute_res, endsWith(Metric,"MAX"))
#' gp1 <- VALIDATION.permute_plot(permute_res=permute_sub, validation_method=validation_method, y_var="stat", save_plot=gsub("\\.csv\\.gz",".stat_values.png",path), width=9, facet_formula = "Assay ~ Model + Tissue")
#' gp2 <- VALIDATION.permute_plot(permute_res=subset(permute_sub, stat>0), validation_method=validation_method, y_var="p", save_plot=gsub("\\.csv\\.gz",".p_values.png",path), width=9,facet_formula = "Assay ~ Model + Tissue")
#' gp12 <- (gp1 / gp2) + patchwork::plot_annotation(tag_levels = "a")
#' ggsave(gsub("\\.csv\\.gz",".png",path),gp12, dpi=400, height=10, width=9)
#'
VALIDATION.permute_plot <- function(permute_res,
                                    validation_method=NULL,
                                    facet_formula=". ~ .",
                                    y_var="stat",
                                    show_plot=T,
                                    save_plot=F,
                                    height=5,
                                    width=7,
                                    verbose=T){
  library(ggplot2)
  plot_dat <- permute_res %>%
    dplyr::mutate(SNP_Group = factor(SNP_Group, levels = unique(permute_res$SNP_Group), ordered = T))
  colorDict <-  colorDict <- c("Random"="grey",
                               "GWAS lead"="red",
                               "UCS"="green2",
                               "Consensus (-POLYFUN)"="goldenrod4",
                               "Consensus"="goldenrod2")

  # Conduct GLM on pval distributions
  permute_res$SNP_Group <- factor(permute_res$SNP_Group, levels = unique(permute_res$SNP_Group))
  metric_count <- length(unique(permute_res$Metric))
  if(metric_count>1){
    printer("VALIDATION:: Facetting plots by metric.",v=verbose)
    glm_res <- lapply(unique(permute_res$Metric), function(metric){
      print(metric)
      VALIDATION.compare_distributions(permute_res=subset(permute_res, Metric==metric),
                                       formula_str = paste(y_var,"~ SNP_Group")) %>%
        dplyr::mutate(Metric=metric)
    }) %>% data.table::rbindlist()
    # Postprocess
    glm_res <- glm_res %>%
      subset(SNP_Group!="(Intercept)") %>%
      dplyr::mutate(SNP_Group = factor(SNP_Group, levels = unique(SNP_Group), ordered = T))

    if(validation_method=="Dey_DeepLearning"){
      glm_res <- glm_res %>%
        tidyr::separate(Metric, sep="_", into=c("Model","Tissue","Assay"), extra="drop", remove=F)
        plot_dat <- plot_dat %>%
          tidyr::separate(Metric, sep="_", into=c("Model","Tissue","Assay"), extra="drop", remove=F)
    }
      glm_violin <-
        ggplot(data=glm_res, aes(x=SNP_Group,y=p, fill=SNP_Group)) +
      geom_violin(alpha=.5) +
      geom_boxplot(alpha=.5) +
      yscale("-log1p", .format = T) +
      geom_jitter(height = 0, alpha=.1) +
      scale_fill_manual(values = colorDict) +
      geom_hline(yintercept = -log10(0.05), alpha=.5, linetype=2) +
      facet_grid(facets = Tissue ~ Model,
                 scales="free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle=45, hjust=1))
      print(glm_violin)
  } else {
    glm_res <- VALIDATION.compare_distributions(permute_res=permute_res,
                                                formula_str = paste(y_var,"~ SNP_Group"))
  }


  # Plot
  # plot_dat$SNP_Group <- forcats::fct_rev( plot_dat$SNP_Group )
  gp <- ggplot(data = plot_dat, aes(x=eval(parse(text=y_var)), fill=SNP_Group, color=SNP_Group,
                                    linetype=SNP_Group)) +
    # geom_histogram(alpha=.5,  position="identity", bins=100) +
    geom_density(position = "identity", adjust=1, alpha=.1) +
    scale_fill_manual(values = colorDict) +
    scale_color_manual(values = colorDict) +
    labs(title=paste("Bootstrapped tests:",validation_method,"tests"),
         x=paste0("bootstrapped ",y_var,"-values")) +
    theme_bw() +
    facet_grid(facets = as.formula(facet_formula))
  # gp + scale_fill_manual(values = rep("transparent",dplyr::n_distinct(plot_dat$SNP_Group)))

  # Get density peaks
  if(metric_count==1){
    b <- ggplot_build(gp)
    density_peaks <- b$data[[1]] %>% dplyr::group_by(fill) %>% dplyr::top_n(n = 1, wt = y) %>%
      dplyr::mutate(SNP_Group =  setNames(names(colorDict), unname(colorDict))[[fill]]) %>%
      data.table::data.table() %>%
      data.table::merge.data.table(glm_res, by="SNP_Group") %>%
      dplyr::mutate(p=ifelse(p*dplyr::n_distinct(PANEL)>1, 1,p*dplyr::n_distinct(PANEL))) %>%
      dplyr::mutate(signif = pvalues_to_symbols(p))
    # Annotate plot
    gp_lab <- gp  +
      ggrepel::geom_label_repel(data=density_peaks,
                                aes(x=x, y=y, label=paste("p =",format(as.numeric(p), digits = 3), signif)),
                                show.legend = F, point.padding = 1,
                                alpha=.75, segment.alpha = .5,
                                color="black") +

      theme(strip.background = element_rect(fill="grey20"),
            strip.text = element_text(color='white'))
  } else {gp_lab <- gp}


  if(show_plot) print(gp_lab)

  if(save_plot!=F){
    ggsave(save_plot, gp_lab, dpi=400, height=height, width=width)
  }
  return(gp_lab)
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


