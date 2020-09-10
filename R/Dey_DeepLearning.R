
#
# DEEPLEARNING.index_annotations <- function(base_url="/sc/arion/projects/pd-omics/data/Dey_DeepLearning"){
#
# }


#' Query deep learning annotations and LDscores
#'
#' Query deep learning annotations and LDscores, and
#' then merge with \code{subset_DT} by \emph{SNP}.
#'
#' @family DEEPLEARNING
#' @example
#' data("BST1")
#' annot.dat <- DEEPLEARNING.query_one_chr(subset_DT=BST1, tissue="NTS", model="Basenji", type="annot")
DEEPLEARNING.query_one_chr <- function(subset_DT,
                                       base_url="/sc/arion/projects/pd-omics/data/Dey_DeepLearning",
                                       level=c("Variant_Level","Allelic_Effect"),
                                       tissue=c("NTS","Blood","Brain"),
                                       model=c("Basenji","BiClassCNN","DeepSEA","ChromHMM","Roadmap","Others"),
                                       # assay=c("DNASE","H3K27ac","H3K4me1","H3K4me3"),
                                       mean_max=c("MEAN","MAX"),
                                       type=c("annot","ldscore"),
                                       nThread=4,
                                       verbose=T){
  level=level[1]; tissue=tissue[1]; model=model[1]; mean_max=mean_max[1]; type=type[1];

  if(level=="Allelic_Effect"){
    printer("DEEPLEARNING:: Only the following models are available when 'level=Variant_Level': Basenji and DeepSEA",v=verbose)
    model <- model[model %in% c("Basenji","DeepSEA")]
  }
  subset_DT$CHR <- gsub("chr","",subset_DT$CHR)
  chrom <- subset_DT$CHR[1]

  printer("DEEPLEARNING:: Searching",base_url,"for matching files",v=verbose)
  # Naming and file formats are inconsistent across folders,
  ## so you have to do a mix of file searching and filtering
  search_dir = file.path(base_url, level, tissue, model)
  search_files = list.files(search_dir,
                            paste("*",chrom,paste0(type,".gz"),sep="\\."),
                            recursive=T, full.names=T)
  # Filter assay
  # search_files <- search_files[grepl(assay, basename(search_files))]
  # Filter mean/max
  mean_max_pattern <- if(mean_max=="MEAN") paste(c("MEAN","AVG"), collapse="|") else mean_max
  search_files <- search_files[grepl(mean_max_pattern, basename(search_files))]

  if(length(search_files)==0)stop("DEEPLEARNING:: No files meeting criterion were found.");
  printer("DEEPLEARNING:: Importing file <<=",search_files,v=verbose)
  annot.dat  <- data.frame(subset_DT)
  for(x in search_files){
    # Get assay
    ## Assign here bc diff tissue-model combinations have somewhat different assays
    assay <- toupper(strsplit(basename(x),"_")[[1]][1])
    new_col_name <- paste(model,tissue,assay,type,mean_max,sep="_")
    # print(x)
    # print(new_col_name)
    annot <- data.table::fread(x, nThread = 1)
    merged <- data.table::merge.data.table(subset_DT,
                                              subset(annot, select=c(SNP,AN)),
                                              all.x = T,
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
#' @example
#' data("merged_DT")
#' ANNOT.DAT <- DEEPLEARNING.query_multi_chr(merged_dat=merged_DT, tissue="NTS", model="Basenji", type="annot")
DEEPLEARNING.query_multi_chr <- function(merged_dat,
                                        base_url="/sc/arion/projects/pd-omics/data/Dey_DeepLearning",
                                        level=c("Variant_Level","Allelic_Effect"),
                                        tissue=c("NTS","Blood","Brain"),
                                        model=c("Basenji","BiClassCNN","DeepSEA","ChromHMM","Roadmap","Others"),
                                        # assay=c("DNASE","H3K27ac","H3K4me1","H3K4me3"),
                                        mean_max=c("MEAN","MAX"),
                                        type=c("annot","ldscore"),
                                        nThread=4,
                                        verbose=T){
  ANNOT.DAT <- parallel::mclapply(unique(merged_dat$CHR), function(chrom){
    printer("Chromosome:",chrom)
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
#' @example
#' data("merged_DT")
#'
#' ANNOT <- DEEPLEARNING.query (merged_dat=merged_DT, level="Allelic_Effect", type="annot")
#'
#' base_url = "/sc/arion/projects/pd-omics/brian/Fine_Mapping"
#' data.table::fwrite(ANNOT , file.path(base_url,"Data/GWAS/Nalls23andMe_2019/_genome_wide/Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.csv.gz"))
DEEPLEARNING.query <- function(merged_dat,
                               base_url="/sc/arion/projects/pd-omics/data/Dey_DeepLearning",
                               level=c("Variant_Level","Allelic_Effect"),
                               tissue=c("NTS","Blood","Brain"),
                               model=c("Basenji","BiClassCNN","DeepSEA","ChromHMM","Roadmap","Others"),
                               # assay=c("DNASE","H3K27ac","H3K4me1","H3K4me3"),
                               mean_max=c("MEAN","MAX"),
                               type=c("annot","ldscore"),
                               nThread=4,
                               verbose=T){
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



#' Melt deep learning anntations into long-format
#'
#' @examples
#' ANNOT <- DEEPLEARNING.query (merged_dat=merged_DT, level="Allelic_Effect", type="annot")
#'
#' annot.melt <- DEEPLEARNING.melt(ANNOT=ANNOT, metric_str="mean")
#' base_url = "/sc/arion/projects/pd-omics/brian/Fine_Mapping"
#' data.table::fwrite(annot.melt, file.path(base_url,"Data/GWAS/Nalls23andMe_2019/_genome_wide/Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.snp_groups.mean.csv.gz"))
DEEPLEARNING.melt <- function(ANNOT,
                              model=c("Basenji","BiClassCNN","DeepSEA","ChromHMM","Roadmap","Others"),
                              metric_str="max",
                              replace_NA=NA,
                              replace_negInf=NA){
  metric <- get(metric_str)
  annot.melt <- ANNOT %>%
    dplyr::group_by(Locus) %>%
    dplyr::summarise_at(.vars = vars(grep(paste(model, collapse="|"),colnames(.), value = T)),
                        .funs = list("Random"= ~ metric(tidyr::replace_na(sample(.x, size=3),replace_NA), na.rm = T),
                                     "All"= ~ metric(tidyr::replace_na(.x,replace_NA), na.rm = T),
                                     "GWAS nom. sig."= ~ metric(tidyr::replace_na(.x[P<.05],replace_NA), na.rm = T),
                                     "GWAS sig."= ~ metric(tidyr::replace_na(.x[P<5e-8],replace_NA), na.rm = T),
                                     "GWAS lead"= ~ metric(tidyr::replace_na(.x[leadSNP],replace_NA), na.rm = T),
                                     "ABF CS"= ~ metric(tidyr::replace_na(.x[ABF.CS>0],replace_NA), na.rm = T),
                                     "SUSIE CS"= ~ metric(tidyr::replace_na(.x[SUSIE.CS>0],replace_NA), na.rm = T),
                                     "POLYFUN-SUSIE CS"= ~ metric(tidyr::replace_na(.x[POLYFUN_SUSIE.CS>0],replace_NA), na.rm = T),
                                     "FINEMAP CS"= ~ metric(tidyr::replace_na(.x[FINEMAP.CS>0],replace_NA), na.rm = T),
                                     "UCS"= ~ metric(tidyr::replace_na(.x[Support>0],replace_NA), na.rm = T),
                                     "Consensus"= ~ metric(tidyr::replace_na(.x[Consensus_SNP],replace_NA), na.rm = T)
                        ),
    ) %>%
    data.table::data.table() %>%
    data.table::melt.data.table(id.vars = "Locus", variable.name = "Annotation") %>%
    # dplyr::mutate(value = as.numeric(gsub(-Inf,replace_negInf,value))) %>%
    tidyr::separate(col = "Annotation", sep="_", into=c("Model","Tissue","Assay","Type","Metric","SNP.Group"), remove=F) %>%
    dplyr::mutate(Annotation = DescTools::StrTrim(Annotation, "_"),
                  SNP.Group = factor(`SNP.Group`,
                                     levels = c("Random","All","GWAS nom. sig.","GWAS sig.","GWAS lead",
                                                "ABF CS","SUSIE CS","POLYFUN-SUSIE CS","FINEMAP CS",
                                                "UCS","Consensus"), ordered = T),
                  log.value = log1p(value))
  return(annot.melt)
}





#' Plot deep learning predictions
#'
#' @examples
#' base_url = "/sc/arion/projects/pd-omics/brian/Fine_Mapping"
#' annot.melt <- data.table::fread(file.path(base_url,"Data/GWAS/Nalls23andMe_2019/_genome_wide/Dey_DeepLearning/Nalls23andMe_2019.Dey_DeepLearning.snp_groups.mean.csv.gz"))
DEEPLEARNING.plot <- function(annot.melt,
                              snp_groups=c("GWAS lead","UCS","Consensus"),
                              comparisons_filter=function(x){if("Consensus" %in% x) return(x)},
                              model.metric=c("MEAN"),
                              show_plot=T,
                              save_path=F,
                              height=9,
                              width=8){
  dat_plot <-  subset(annot.melt,
                      Metric %in% model.metric &
                        SNP.Group %in%  snp_groups)
  snp.groups <- unique( dat_plot$SNP.Group)
  comparisons <- utils::combn(x = as.character(snp.groups),
                              m=2,
                              FUN = comparisons_filter,
                              simplify = F) %>% purrr::compact()
  method="wilcox.test"
  pb <- ggpubr::ggviolin(data = dat_plot,
                   x = "SNP.Group", y="log.value", fill = "SNP.Group",
                   add = "boxplot", alpha = .6, trim = T,
                   add.params = list(alpha=.1) )  +
    ggpubr::stat_compare_means(method = method,
                               comparisons = comparisons,
                               label = "p.signif", size=3, vjust = 1.5) +
    # ggpubr::stat_compare_means(method = method,
    #                            comparisons = comparisons,
    #                            label = "p.adj", hjust=5, size=3) +
    geom_jitter(alpha=.1, width = .3) + #aes(color=SNP.Group)) +
    facet_grid(facets =   Assay ~ Model + Tissue ,
               scales = "free") +
    labs(title=paste0("Deep learning annotations (",tolower(model.metric),")"),
         y='log10(value)', x="SNP Group") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          legend.position = "none",
          strip.background = element_rect(fill = "grey20"),
          strip.text= element_text(color = "white"))
  if(length(dplyr::union(snp.groups, c("GWAS lead","UCS","Consensus")))==3){
    pb <- pb + scale_fill_manual(values =  c("red","green2","goldenrod2"))
  }
  if(show_plot) print(pb)
  if(save_path!=F){
    # save_path <- file.path("/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/Dey_DeepLearning",
    #                        paste("Nalls23andMe_2019","Dey_DeepLearning",model.metric,"png",sep="."))
    ggsave(save_path,
           pb, dpi = 300, height=height, width=width)
  }
  return(pb)
}





DEEPLEARNING.permut_test <- function(ANNOT){
  metric_str = "max"
  annot.melt <- DEEPLEARNING.melt(ANNOT = ANNOT,
                                  metric_str = metric_str)
  # coin::independence_test(data = as.data.frame(annot.melt),
  #                         formula = value ~ SNP.Group)
  # oneway_test(data= as.data.frame(annot.melt),
  #             formula = value ~ SNP.Group + Tissue + Model)
  #
  # independence_test(colonies~place,
  #             data=ants, )

}
