
#
# DEEPLEARNING.index_annotations <- function(base_url="/sc/arion/projects/pd-omics/data/Dey_DeepLearning"){
#
# }

# *******
# Original Data Source:
# https://alkesgroup.broadinstitute.org/LDSCORE/DeepLearning/Dey_DeepLearning.tgz
# *******


#' Query deep learning annotations and LDscores
#'
#' Query deep learning annotations and LDscores, and
#' then merge with \code{subset_DT} by \emph{SNP}.
#'
#' @family DEEPLEARNING
#' @examples
#' \dontrun{
#' BST1 <- echodata::BST1
#' annot.dat <- DEEPLEARNING.query_one_chr(subset_DT=BST1, tissue="NTS", model="Basenji", type="annot")
#' }
DEEPLEARNING.query_one_chr <- function(subset_DT,
                                       base_url="/sc/arion/projects/pd-omics/data/Dey_DeepLearning",
                                       level=c("Variant_Level","Allelic_Effect"),
                                       tissue=c("NTS","Blood","Brain"),
                                       model=c("Basenji","BiClassCNN","DeepSEA","ChromHMM","Roadmap","Others"),
                                       # assay=c("DNASE","H3K27ac","H3K4me1","H3K4me3"),
                                       mean_max=c("MEAN","MAX"),
                                       type=c("annot","ldscore"),
                                       nThread=1,
                                       verbose=TRUE){
  level=level[1]; tissue=tissue[1]; model=model[1]; mean_max=mean_max[1]; type=type[1];

  if(level=="Allelic_Effect"){
    messager("DEEPLEARNING:: Only the following models are available when 'level=Variant_Level': Basenji and DeepSEA",v=verbose)
    model <- model[model %in% c("Basenji","DeepSEA")]
  }
  subset_DT$CHR <- gsub("chr","",subset_DT$CHR)
  chrom <- subset_DT$CHR[1]

  messager("DEEPLEARNING:: Searching",base_url,"for matching files",v=verbose)
  # Naming and file formats are inconsistent across folders,
  ## so you have to do a mix of file searching and filtering
  search_dir = file.path(base_url, level, tissue, model)
  search_files = list.files(search_dir,
                            paste("*",chrom,paste0(type,".gz"),sep="\\."),
                            recursive=TRUE, full.names=TRUE)
  # Filter assay
  # search_files <- search_files[grepl(assay, basename(search_files))]
  # Filter mean/max
  mean_max_pattern <- if(mean_max=="MEAN") paste(c("MEAN","AVG"), collapse="|") else mean_max
  search_files <- search_files[grepl(mean_max_pattern, toupper(basename(search_files)))]

  if(length(search_files)==0)stop("DEEPLEARNING:: No files meeting criterion were found.");
  messager("DEEPLEARNING:: Importing file <<=",search_files,v=verbose)
  annot.dat  <- data.frame(subset_DT)
  for(x in search_files){
    # Get assay
    ## Assign here bc diff tissue-model combinations have somewhat different assays
    assay <- toupper(strsplit(basename(x),"_")[[1]][1])
    new_col_name <- paste(model,tissue,assay,type,mean_max,sep="_")
    # print(x)
    # print(new_col_name)
    annot <- data.table::fread(x, nThread = 1)
    if("L2" %in% colnames(annot)) annot <- dplyr::rename(AN=L2)
    merged <- data.table::merge.data.table(subset_DT,
                                              subset(annot, select=-c(SNP,AN)),
                                              all.x = TRUE,
                                              by="SNP") %>% data.frame()
    annot.dat[ new_col_name ] <-  merged$AN
  }
  return(data.table::data.table(annot.dat))
}





#' Query deep learning annotations and LDscores (iterate)
#'
#' Query deep learning annotations and LDscores, and
#' then merge with \code{subset_DT} by \emph{SNP}.
#'  Repeat for each locus,
#'
#' @family DEEPLEARNING
#' @examples
#' \dontrun{
#' data("merged_DT")
#' ANNOT.DAT <- DEEPLEARNING.query_multi_chr(merged_dat=merged_DT, tissue="NTS", model="Basenji", type="annot")
#' }
DEEPLEARNING.query_multi_chr <- function(merged_dat,
                                        base_url="/sc/arion/projects/pd-omics/data/Dey_DeepLearning",
                                        level=c("Variant_Level","Allelic_Effect"),
                                        tissue=c("NTS","Blood","Brain"),
                                        model=c("Basenji","BiClassCNN","DeepSEA","ChromHMM","Roadmap","Others"),
                                        # assay=c("DNASE","H3K27ac","H3K4me1","H3K4me3"),
                                        mean_max=c("MEAN","MAX"),
                                        type=c("annot","ldscore"),
                                        nThread=1,
                                        verbose=TRUE){
  ANNOT.DAT <- parallel::mclapply(unique(merged_dat$CHR), function(chrom){
    messager("Chromosome:",chrom)
    DEEPLEARNING.query_one_chr(subset_DT=subset(merged_dat, CHR==chrom),
                                base_url=base_url,
                                level=level,
                                tissue=tissue,
                                model=model,
                                # assay=assay,
                                mean_max=mean_max,
                                type=type,
                                nThread = 1)
  }, mc.cores = nThread) %>% data.table::rbindlist()
  return(ANNOT.DAT)
}




#' Iteratively collect deep learning annotations
#'
#' @family DEEPLEARNING
#' @examples
#' \dontrun{
#' root = "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
#' merged_dat <- merge_finemapping_results(dataset = "Data/GWAS/Nalls23andMe_2019", LD_reference = "UKB", minimum_support = 0)
#' merged_dat <- echodata::find_consensus_snps_no_polyfun(merged_dat)
#'
#' #### Allelic_Effect ####
#' ANNOT.ae <- DEEPLEARNING.query (merged_dat=merged_dat, level="Allelic_Effect", type="annot")
#' data.table::fwrite(ANNOT.ae, file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Allelic_Effect.csv.gz"))
#'
#'
#' #### Variant_Level ####
#' ANNOT.vl <- DEEPLEARNING.query (merged_dat=merged_dat, level="Variant_Level", type="annot")
#' data.table::fwrite(ANNOT.vl, file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Variant_Level.csv.gz"))
#'
#' }
DEEPLEARNING.query <- function(merged_dat,
                               base_url="/sc/arion/projects/pd-omics/data/Dey_DeepLearning",
                               level=c("Variant_Level","Allelic_Effect"),
                               tissue=c("NTS","Blood","Brain"),
                               model=c("Basenji","BiClassCNN","DeepSEA","ChromHMM","Roadmap","Others"),
                               # assay=c("DNASE","H3K27ac","H3K4me1","H3K4me3"),
                               mean_max=c("MEAN","MAX"),
                               type=c("annot","ldscore"),
                               nThread=1,
                               verbose=TRUE){
  # Get every combination of each argument
  param_combs <- expand.grid(level=level, tissue=tissue, model=model, mean_max=mean_max, type=type)

  col.list <- lapply(1:nrow(param_combs), function(i){
    message("\nParameter combination: ",i)
    print(param_combs[i,])
    message("\n")
    annot_cols <- NULL
    try({
      annot <- DEEPLEARNING.query_multi_chr(merged_dat = merged_dat,
                                           base_url = base_url,
                                           level = param_combs[i,"level"],
                                           tissue = param_combs[i,"tissue"],
                                           model = param_combs[i,"model"],
                                           # assay = param_combs[i,"assay"],
                                           mean_max = param_combs[i,"mean_max"],
                                           type = param_combs[i,"type"],
                                           nThread = 1)
      annot_cols <- dplyr::select(annot, !dplyr::any_of(c("BP","CM",colnames(merged_dat))) )
    })
    return(annot_cols)
  })
  ANNOT <- cbind(merged_dat, do.call("cbind",col.list))
  return(ANNOT)
}




#' Melt deep learning annotations into long-format
#'
#' @examples
#' \dontrun{
#'  root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
#' ## merged_dat <- merge_finemapping_results(dataset = "Data/GWAS/Nalls23andMe_2019", minimum_support = 0, LD_reference = "UKB")
#' ## ANNOT <- DEEPLEARNING.query(merged_dat=merged_dat, level="Allelic_Effect", type="annot")
#'
#' #### Allelic_Effect ####
#' path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Allelic_Effect.csv.gz")
#'
#' #### Variant_Level ####
#' path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Variant_Level.csv.gz")
#'
#' ANNOT <- data.table::fread(path, nThread=8)
#' ANNOT <- echodata::find_consensus_snps_no_polyfun(ANNOT)
#' annot.melt <- DEEPLEARNING.melt(ANNOT=ANNOT, aggregate_func="mean", save_path=gsub("\\.csv\\.gz",".snp_groups_mean.csv.gz",path))
#' }
DEEPLEARNING.melt <- function(ANNOT,
                              model=c("Basenji","BiClassCNN","DeepSEA","ChromHMM","Roadmap","Others"),
                              aggregate_func="mean",
                              replace_NA=NA,
                              replace_negInf=NA,
                              save_path=FALSE){
  snp.groups_list <- echodata::snp_group_filters()
  agg_func <- get(aggregate_func)
  sampling_df <- ANNOT
  annot.melt <- ANNOT %>%
    dplyr::group_by(Locus) %>%
    dplyr::summarise_at(.vars =dplyr::vars(grep(paste(model, collapse="|"),colnames(.), value = TRUE)),
                        .funs = list("Random"= ~ agg_func(tidyr::replace_na(sample(.x, size=3, replace = TRUE),replace_NA), na.rm = TRUE),
                                     "All"= ~ agg_func(tidyr::replace_na(.x,replace_NA), na.rm = TRUE),
                                     "GWAS nom. sig."= ~ agg_func(tidyr::replace_na(.x[P<.05],replace_NA), na.rm = TRUE),
                                     "GWAS sig."= ~ agg_func(tidyr::replace_na(.x[P<5e-8],replace_NA), na.rm = TRUE),
                                     "GWAS lead"= ~ agg_func(tidyr::replace_na(.x[leadSNP],replace_NA), na.rm = TRUE),
                                     "ABF CS"= ~ agg_func(tidyr::replace_na(.x[ABF.CS>0],replace_NA), na.rm = TRUE),
                                     "SUSIE CS"= ~ agg_func(tidyr::replace_na(.x[SUSIE.CS>0],replace_NA), na.rm = TRUE),
                                     "POLYFUN-SUSIE CS"= ~ agg_func(tidyr::replace_na(.x[POLYFUN_SUSIE.CS>0],replace_NA), na.rm = TRUE),
                                     "FINEMAP CS"= ~ agg_func(tidyr::replace_na(.x[FINEMAP.CS>0],replace_NA), na.rm = TRUE),
                                     "UCS (-PolyFun)"= ~ agg_func(tidyr::replace_na(.x[Support_noPF>0],replace_NA), na.rm = TRUE),
                                     "UCS"= ~ agg_func(tidyr::replace_na(.x[Support>0],replace_NA), na.rm = TRUE),
                                     "Support==0"= ~ agg_func(tidyr::replace_na(.x[Support==0],replace_NA), na.rm = TRUE),
                                     "Support==1"= ~ agg_func(tidyr::replace_na(.x[Support==1],replace_NA), na.rm = TRUE),
                                     "Support==2"= ~ agg_func(tidyr::replace_na(.x[Support==2],replace_NA), na.rm = TRUE),
                                     "Support==3"= ~ agg_func(tidyr::replace_na(.x[Support==3],replace_NA), na.rm = TRUE),
                                     "Support==4"= ~ agg_func(tidyr::replace_na(.x[Support==4],replace_NA), na.rm = TRUE),
                                     "Consensus (-PolyFun)"= ~ agg_func(tidyr::replace_na(.x[Consensus_SNP_noPF],replace_NA), na.rm = TRUE),
                                     "Consensus"= ~ agg_func(tidyr::replace_na(.x[Consensus_SNP],replace_NA), na.rm = TRUE)
                        ),
    ) %>%
    data.table::data.table() %>%
    data.table::melt.data.table(id.vars = "Locus", variable.name = "Annotation") %>%
    # dplyr::mutate(value = as.numeric(gsub(-Inf,replace_negInf,value))) %>%
    tidyr::separate(col = "Annotation", sep="_", into=c("Model","Tissue","Assay","Type","Metric","SNP_group"), remove=FALSE) %>%
    dplyr::mutate(Annotation = DescTools::StrTrim(Annotation, "_"),
                  SNP_group = factor(SNP_group,
                                     levels = names(snp.groups_list), ordered = TRUE))
  if(save_path!=FALSE){
    messager("DEEPLEARNING:: Saving aggregated SNP_group values",aggregate_func,"==>",save_path)
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    data.table::fwrite(annot.melt, save_path)
  }
  return(annot.melt)
}





#' Plot deep learning predictions
#'
#' @examples
#' \dontrun{
#' root = "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
#'
#' #### Allelic_Effect ####
#' path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Allelic_Effect.snp_groups_mean.csv.gz")
#'
#' #### Variant_Level ####
#' path <- file.path(root,"Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.annot.Variant_Level.snp_groups_mean.csv.gz")
#'
#' annot.melt <-  data.table::fread(path, nThread=8)
#'
#' gp <- DEEPLEARNING.plot(annot.melt=annot.melt, facet_formula="Tissue ~ Model", comparisons_filter=NULL, save_path=gsub("\\.csv\\.gz",".png",path))
#' }
DEEPLEARNING.plot <- function(annot.melt,
                              snp_groups=c("GWAS lead","UCS","Consensus (-POLYFUN)","Consensus"),
                              comparisons_filter=function(x){if("Consensus" %in% x) return(x)},
                              model_metric=c("MAX"),
                              facet_formula=". ~ Model",
                              remove_outliers=TRUE,
                              show_padj=TRUE,
                              show_signif=TRUE,
                              vjust_signif=0,
                              show_plot=TRUE,
                              save_path=FALSE,
                              height=6,
                              width=8){

  if(remove_outliers){
    # https://www.r-bloggers.com/2020/01/how-to-remove-outliers-in-r/
    outliers <- boxplot(annot.melt$value, plot=FALSE)$out
    annot.melt <- annot.melt[-which(annot.melt$value %in% outliers),]
  }
  colorDict <- echodata::snp_group_colorDict()
  dat_plot <-  subset(annot.melt,
                      Metric %in% model_metric &
                        SNP_group %in%  snp_groups) %>%
    dplyr::mutate(SNP_group=factor(SNP_group, levels = names(colorDict), ordered = TRUE))
    # dplyr::group_by(SNP_group, Locus, Model, Tissue) %>%
    # dplyr::summarise(value=mean(value,na.rm = TRUE))
  snp.groups <- unique(dat_plot$SNP_group)
  comparisons <- utils::combn(x = as.character(snp.groups),
                              m=2,
                              FUN = comparisons_filter,
                              simplify  = FALSE) %>% purrr::compact()
  method="wilcox.test"
  gp <- ggplot(data = dat_plot, aes(x=SNP_group, y=value, fill=SNP_group)) +
    geom_jitter(alpha=.1, width = .3, height=0) +
    # geom_bar(stat = "identity",alpha=.6) +
    geom_violin(alpha = .6) +
    geom_boxplot(alpha=.6, color="black") +
    facet_grid(facets = as.formula(facet_formula),
               scales = "free") +
    labs(title=paste0("Deep learning annotations (",tolower(model_metric),")"),
         y='value', x="SNP Group") +
    scale_fill_manual(values = colorDict) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          legend.position = "none",
          strip.background = element_rect(fill = "grey20"),
          strip.text= element_text(color = "white"))
    # coord_cartesian(ylim = quantile(dat_plot$value, c(0.1, 0.9), na.rm = TRUE))
  if(show_padj){
    gp <- gp + ggpubr::stat_compare_means(method = method,
                                          comparisons = comparisons,
                                          label = "p.adj", size=3,  vjust=2)
  }
  if(show_signif){
    gp <- gp + ggpubr::stat_compare_means(method = method,
                                          comparisons = comparisons,
                                          label = "p.signif", size=3, vjust=vjust_signif)
  }
  if(show_plot) print(gp)
  # ggplot(data =dat_plot, aes(x=value, fill=SNP_group)) +
  #   geom_density(position="identity", alpha=.5) +
  #   theme_bw() +
  #   scale_fill_manual(values = colorDict) +
  #   xlim(c(0,.05))
  if(save_path!=FALSE){
    # save_path <- file.path("/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/Dey_DeepLearning",
    #                        paste("Nalls23andMe_2019","Dey_DeepLearning",model.metric,"png",sep="."))
    ggsave(save_path,
           gp, dpi = 300, height=height, width=width)
  }
  return(gp)
}




