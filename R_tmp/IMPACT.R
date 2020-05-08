
# IMPACT #
# https://github.com/immunogenomics/IMPACT



IMPACT.get_annotation_key <- function(URL="https://github.com/immunogenomics/IMPACT/raw/master/IMPACT707/IMPACT_annotation_key.txt",
                                      save_path="./echolocatoR/annotations/IMPACT/IMPACT_annotation_key.txt.gz",
                                      force_new_download=F){
  if(file.exists(save_path) & force_new_download==F){
    print("+ IMPACT:: Importing local anotation key...")
    annot.key <- data.table::fread(save_path)
  } else {
    print("+ IMPACT:: Downloading annotation key from GitHub...")
    annot.key <- data.table::fread(URL)
    data.table::fwrite(annot.key, save_path, sep="\t")
    R.utils::gzip(save_path)
  }
  annot.key$Annot <- as.factor(paste0("Annot",annot.key$IMPACT))
  return(annot.key)
}


IMPACT.get_annotations <- function(baseURL="https://github.com/immunogenomics/IMPACT/raw/master/IMPACT707/Annotations",
                                   chrom=NULL, 
                                   subset_DT=NULL,
                                   nThread=4,
                                   all_snps_in_range=F){ 
  # These are large files stored via GitHub's large file storage (lfs) 
  # Install LFS: https://git-lfs.github.com
  # Getting started with LFS: https://www.atlassian.com/git/tutorials/git-lfs
  # Download metadata / raw data for specific files with curl:  https://docs.gitlab.com/ee/api/repository_files.html#get-file-from-repository
  
  if(!is.null(subset_DT)){
    chrom <- subset_DT$CHR[1]
  }
  # https://github.com/immunogenomics/IMPACT.git
  # "curl --header 'https://git-lfs.github.com/spec/v1/projects/13083/repository/files/app%2Fmodels%2Fkey%2Erb/raw?ref=master'"
  
  
  URL <- file.path(baseURL, paste0("IMPACT707_EAS_chr",chrom,".annot.gz"))
  annot <- data.table::fread(URL, nThread = nThread)
  
  if(!is.null(subset_DT)){ 
    annot_merge <- data.table:::merge.data.table(data.table::data.table(subset_DT),
                                                 annot,
                                                 by.x = c("SNP","CHR","POS"), 
                                                 by.y = c("SNP","CHR","BP"), 
                                                 all.x = T,
                                                 all.y = all_snps_in_range)
  } else {annot_merge <- annot}
  annot_merge <- subset(annot_merge, POS>=min(subset_DT$POS) & POS<=max(subset_DT$POS))
  # Merge with metadata
  annot.key <- IMPACT.get_annotation_key()
  annot_cols <- grep("^Annot*",colnames(annot_merge), value = T)
  annot_melt <- data.table:::melt.data.table(annot_merge, measure.vars = annot_cols,
                                             variable.name = "Annot",
                                             value.name = "IMPACT_score", 
                                             na.rm=F) %>%
    data.table:::merge.data.table(annot.key,
                                  by="Annot", 
                                  all = T,
                                  allow.cartesian = T)
  return(annot_melt)
} 



# quick_finemap()
# annot_melt <- IMPACT.get_annotations(baseURL = "/Volumes/Steelix/IMPACT/IMPACT707/Annotations", subset_DT = subset_DT)

IMPACT.iterate_get_annotations <- function(merged_DT,
                                           IMPACT_score_thresh=.1,
                                           baseURL="/Volumes/Steelix/IMPACT/IMPACT707/Annotations",
                                           all_snps_in_range=T,
                                           top_annotations=F,
                                           force_one_annot_per_locus=F,
                                           snp.filter="!is.na(SNP)"){  
  ANNOT_MELT <- lapply(unique(merged_DT$Locus), function(locus){
    message("+ IMPACT:: Gathering annotations for Locus = ",locus)
    try({
      subset_DT <- subset(merged_DT, Locus==locus)
      annot_melt <- IMPACT.get_annotations(
                                           baseURL = baseURL,
                                           # baseURL = "../../data/IMPACT/IMPACT707/Annotations",
                                           subset_DT = subset_DT, 
                                           all_snps_in_range = all_snps_in_range,
                                           nThread = 4)
      if(top_annotations!=F){
        top_impact <- IMPACT.get_top_annotations(ANNOT_MELT = annot_melt,
                                                 snp.filter = snp.filter,
                                                 top_annotations = top_annotations, 
                                                 force_one_annot_per_locus = force_one_annot_per_locus)
        annot_melt <- subset(annot_melt, Tissue %in% top_impact$Tissue & 
                                         CellDeriv %in% top_impact$CellDeriv & 
                                         Cell %in% top_impact$Cell &
                                         TF %in% top_impact$TF)
      }
      annot_melt <- subset(annot_melt, IMPACT_score >= IMPACT_score_thresh)
      printer("+ IMPACT::",nrow(annot_melt),"annotations found at IMPACT_score ≥", IMPACT_score_thresh)
      return(annot_melt)
    })
  }) %>% data.table::rbindlist()
 return(ANNOT_MELT)
}




find_topConsensus <- function(dat){
  # CGet top consensus SNPs
  top.consensus <- (dat %>% 
                      dplyr::group_by(Locus) %>% 
                      subset(Consensus_SNP) %>% 
                      dplyr::arrange(-mean.PP) %>%
                      # dplyr::arrange(-IMPACT_score) %>%
                      dplyr::slice(1))$SNP %>% unique()
  dat$topConsensus <- dat$SNP %in% top.consensus
  return(dat)
}


IMPACT.postprocess_annotations <- function(ANNOT_MELT,
                                           order_loci = T,
                                           no_no_loci = NULL){ 
  # ANNOT_MELT <- data.table::fread("Data/GWAS/Nalls23andMe_2019/_genome_wide/IMPACT/IMPACT_overlap.csv.gz") 
  sort_loci_by_topConsensus_impact <- function(ANNOT_MELT){
    locus_sort <- ANNOT_MELT %>%  
      dplyr::group_by(Locus, .drop=F) %>%  
      dplyr::summarise(meanIMPACT=mean(IMPACT_score[topConsensus], na.rm = T)) %>%  
      dplyr::arrange(meanIMPACT) %>% 
      dplyr::mutate(Locus=factor(Locus, ordered = T)) 
    locus_sort[is.na(locus_sort$meanIMPACT),"meanIMPACT"] <- 0
    ANNOT_MELT$Locus <- factor(ANNOT_MELT$Locus, levels = locus_sort$Locus, ordered = T)
    
    ANNOT_MELT <- ANNOT_MELT%>%
      dplyr::mutate(#Tissue = ifelse(Tissue=="GI", Tissue, tolower(Tissue)),
                    Cell_Type = paste0(ifelse(CellDeriv=="mesendoderm","hESC derived mesendodermal cells",CellDeriv),
                                     ifelse(Cell %in% c("hESC derived mesendodermal cells","."),"",paste0(" (",Cell,")") ))  )
    return( ANNOT_MELT)
  }
  
  # sort_loci_by_annot <- function(ANNOT_MELT){
  #   locus_sort <- ANNOT_MELT %>%
  #     dplyr::group_by(Locus, .drop=F) %>%
  #     dplyr::arrange(TF, Tissue, CellDeriv, Cell) %>%  
  #     dplyr::mutate(Locus=factor(Locus, ordered = T)) 
  #   ANNOT_MELT$Locus <- factor(ANNOT_MELT$Locus, levels = locus_sort$Locus, ordered = T)
  #   return(ANNOT_MELT)
  # }
  # 
  # Fill NA IMPACT_score (<.1) with 0
  ANNOT_MELT <- find_topConsensus(ANNOT_MELT)

  # Check that there's actually 1/locus
  # ANNOT_MELT %>% dplyr::group_by(Locus) %>%  
  #   dplyr::summarise(count = n_distinct(SNP[topConsensus])) %>% data.frame()
  if(order_loci){
    ANNOT_MELT <- sort_loci_by_topConsensus_impact(ANNOT_MELT)
  }
 
  if(!is.null(no_no_loci)){
    ANNOT_MELT <- subset(ANNOT_MELT, !Locus %in% no_no_loci)
  }
  
  
  ANNOT_MELT$uniqueID <- 1:nrow(ANNOT_MELT) 
  return(ANNOT_MELT)
}




# get_mean_IMPACT
IMPACT.get_top_annotations <- function(ANNOT_MELT,
                                       snp.filter="!is.na(IMPACT_score)",
                                       top_annotations=1,
                                       force_one_annot_per_locus=F){  
  top_impact <- ANNOT_MELT %>%  
    # Group all SNPs within the snp group (snp filter)
    dplyr::group_by(Locus, Tissue, CellDeriv, Cell, TF, .drop=F) %>%
    dplyr::summarise(SNPs = n_distinct(SNP[eval(parse(text = snp.filter))], na.rm = T),
                     # These metrics are often very similar or the same 
                     ## when you only have 1 snp (in the snp group)p er locus
                     mean_IMPACT = mean(IMPACT_score[eval(parse(text = snp.filter))], na.rm=T),
                     max_IMPACT = max(IMPACT_score[eval(parse(text = snp.filter))], na.rm=T),
                     median_IMPACT = median(IMPACT_score[eval(parse(text = snp.filter))], na.rm = T)
                     ) 
  if(any(is.infinite(top_impact$max_IMPACT))){
    # Max function produces -Inf values when there's only NA values and na.rm=T
    top_impact[is.infinite(top_impact$max_IMPACT),"max_IMPACT"] <- NA
  } 
  
  # Get the annotation with the top mean/max/mode 
  if(top_annotations!=F){
    top_impact <- top_impact %>% 
      # dplyr::ungroup() %>%
      dplyr::group_by(Locus) %>%
      dplyr::top_n(n=top_annotations, wt=max_IMPACT)
  } 
  # Sometimes top_n gives multiple rows/group when they have the same value. 
  ## Use slice to ensure only one row/group.
  if(force_one_annot_per_locus){
    top_impact <- top_impact %>% 
      dplyr::group_by(Locus) %>% 
      dplyr::arrange(-mean_IMPACT, -max_IMPACT, -median_IMPACT) %>%
      dplyr::slice(1)
  } 
  top_impact$snp.filter <- snp.filter
  return(top_impact)
}
  
  
prepare_mat_meta <- function(TOP_IMPACT,
                             TOP_IMPACT_all,
                             snp.group="Consensus",
                             value.var="mean_IMPACT",
                             fill_na=0){
  # Reorder loci
  TOP_IMPACT <- TOP_IMPACT %>% dplyr::arrange(Tissue, CellDeriv, Cell, TF) 
  TOP_IMPACT$Locus <- factor(as.character(TOP_IMPACT$Locus), 
                             levels=unique(as.character(TOP_IMPACT$Locus)), ordered = T)
  # Merge the meanIMPACT scores with the top annotations so that you only have one row/locus
  mat_meta <- TOP_IMPACT %>% 
    dplyr::group_by(Locus, SNP_group) %>%
    dplyr::summarise(mean_IMPACT=mean(mean_IMPACT, na.rm=T),
                     max_IMPACT=gsub(-Inf,NA, max(max_IMPACT, na.rm=T)),
                     median_IMPACT=mean(median_IMPACT, na.rm=T)) %>%
    data.table::data.table() %>%
    data.table::dcast(formula ="Locus ~ SNP_group",
                      value.var=value.var) %>% 
    # Get the annotations from here
    merge(subset(TOP_IMPACT, SNP_group==snp.group), sort = F) %>%
    # REPLACE NA with 0!: in actuality, these are just snps with IMPACT score <.01 during query
    dplyr::mutate_at(.vars = names(snp.groups),
                     .funs = function(.){as.numeric(ifelse(is.na(.),fill_na,.))}) %>%
    dplyr::mutate(Tissue=gsub("STEMCELL","STEM CELL",Tissue)) %>%
    dplyr::mutate(Tissue = ifelse(Tissue=="GI","GI", 
                                  stringr::str_to_sentence(Tissue))) %>%
    arrange(Tissue, CellDeriv, Cell, TF)  %>% 
    `rownames<-`(.[["Locus"]]) %>%
    subset(select=-Locus)
  return(mat_meta)
}


  
IMPACT_heatmap <- function(ANNOT_MELT){
  library(dplyr)
  library(ComplexHeatmap);# devtools::install_github("jokergoo/ComplexHeatmap")
  # merged_DT <- merge_finemapping_results(minimum_support = 1, include_leadSNPs = T,xlsx_path = F,dataset = "Data/GWAS/Nalls23andMe_2019")
  # merged_DT$Locus <- merged_DT$Gene
  # ANNOT_MELT <- IMPACT.iterate_get_annotations(merged_DT = merged_DT,
  #                                              IMPACT_score_thresh=0,
  #                                              baseURL="../../data/IMPACT/IMPACT707/Annotations",
  #                                              all_snps_in_range=T,
  #                                              top_annotations_only=F)
  # data.table::fwrite(ANNOT_MELT,"../../data/IMPACT/IMPACT707/Annotations/IMPACT_overlap.csv.gz")
  ANNOT_MELT <- data.table::fread("Data/GWAS/Nalls23andMe_2019/_genome_wide/IMPACT/IMPACT_overlap.csv.gz")
  no_no_loci =  c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1",  
                  "RETREG3","UBTF","FAM171A2","MAP3K14","CRHR1","MAPT-AS1","KANSL1","NSF","WNT3")
  ANNOT_MELT <- IMPACT.postprocess_annotations(ANNOT_MELT, no_no_loci = no_no_loci)
  
  
  
  
  snp.groups <- c("Lead_GWAS" = "leadSNP==T",
                  "UCS"="Support>0",
                  "Consensus"="Consensus_SNP"
                  # "Top_Consensus"="topConsensus"
                  )
  TOP_IMPACT <- lapply(names(snp.groups), function(x){
    print(x)   
    top_impact <- IMPACT.get_top_annotations(ANNOT_MELT = ANNOT_MELT, 
                                             snp.filter = snp.groups[[x]], 
                                             top_annotations = 1,
                                             force_one_annot_per_locus = T)  
    top_impact$SNP_group <- x
    return(top_impact)
  }) %>% data.table::rbindlist()
  TOP_IMPACT <- subset(TOP_IMPACT, !Locus %in% no_no_loci)
  mat_meta <- prepare_mat_meta(TOP_IMPACT = TOP_IMPACT, 
                               # TOP_IMPACT_all = TOP_IMPACT_all,
                               snp.group = "Consensus",
                               value.var = "mean_IMPACT", 
                               fill_na = 0)
  # Prepare data for boxplot and t-tests
  TOP_IMPACT_all <- lapply(names(snp.groups), function(x){
    print(x)   
    top_impact <- IMPACT.get_top_annotations(ANNOT_MELT = ANNOT_MELT, 
                                             snp.filter = snp.groups[[x]], 
                                             # The number of top annots affects the significance of the t-test
                                             top_annotations = 1,
                                             force_one_annot_per_locus = F)  
    top_impact$SNP_group <- x
    return(top_impact)
  }) %>% data.table::rbindlist()
  TOP_IMPACT_all <- subset(TOP_IMPACT_all, !Locus %in% no_no_loci) 
  # Use a less averaged version of the data to gain power  
  TOP_IMPACT_all$rowID <- 1:nrow(TOP_IMPACT_all)
  # Depending on if and at what stage you fill na with 0, you get very different boxplot and tstats
  TOP_IMPACT_all[is.na(TOP_IMPACT_all)] <- 0
  boxplot_mat <- data.table::dcast(TOP_IMPACT_all,
                                  formula ="rowID ~ SNP_group",
                                  value.var="mean_IMPACT") %>%
   dplyr::select(Lead_GWAS, UCS, Consensus)
 
  IMPACT.snp_group_boxplot(TOP_IMPACT_all)
  res <- pairwise.t.test(x =  TOP_IMPACT_all$mean_IMPACT, 
                         g =  TOP_IMPACT_all$SNP_group)
  res
  
  
  # ComplexHeatmap
  # https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#colors
  main_heatmap <- function(mat, cex=.5, width=4.5, row_split=NULL){
    colnames(mat) <- gsub("_"," ",colnames(mat))
    ha = Heatmap(mat, name = "IMPACT score",
                 width = unit(width, "cm"),
                 # na_col = "gainsborough", 
                       heatmap_legend_param = list(
                         title = "mean\nIMPACT\nscore",  
                         at = c(1, .5, 0), 
                         labels = c("1", "0.5", "< 0.1"),
                         legend_height = unit(3, "cm")
                       ),
                      
                      # Columns
                      col = rev(RColorBrewer::brewer.pal(n=11, "Spectral")),
                      column_title="SNP group", 
                      # column_gap = 1,
                      # column_split = colnames(mat),
                      # column_order = colnames(mat),
                      cluster_columns = F, 
                      # column_dend_side = "top", 
                      column_names_side = "top",
                      column_names_rot = 0,
                      column_names_gp = gpar(col = c("red","green3","goldenrod3","goldenrod2"), cex=.75),
                      column_names_centered = T,
                      column_gap = unit(.5,"cm"),
                      
                      # Rows
                      show_row_names = T,
                      row_title = "Locus", 
                      row_names_side = "left",
                      row_names_gp = gpar(cex=cex), 
                      row_dend_reorder = F,  
                      row_split=row_split,
                 
                      cluster_rows = F, 
                      # row_km = 4, 
                      # row_dend_width =  unit(2, "cm"),
                      top_annotation = HeatmapAnnotation("mean IMPACT scores" = anno_boxplot(boxplot_mat, gp = gpar(fill = c("red","green3","goldenrod3","goldenrod2"), col = "black", border = "gray", alpha=.75)),
                                                         annotation_name_side = "left", 
                                                         height=unit(3, "cm"),
                                                         annotation_name_gp = gpar(cex=.5))
                      # top_annotation = columnAnnotation("} Mean IMPACT score (all loci)"=colMeans(mat,na.rm = T ) ) 
    )
    return(ha)
  }
  
  
  annotation_table <- function(mat_meta, cex=.6, text_color="white", rot=0, width_factor=.6){ 
    meta <- mat_meta[,c("Tissue","CellDeriv","Cell")] 
    meta <- meta %>% dplyr::rename(`Cell type/origin`=CellDeriv, `Cell subtype`=Cell)
    # widths <- lapply(colnames(mat_meta),function(variable){max_text_width(mat_meta[[variable]])}) %>% unlist()
    # widths <- widths*width_factor
     
    master_dict <- hierarchical_colors(mat_meta)
    ha = Heatmap(meta, name = "Metadata_table",  
                 # width =widths,
                 show_heatmap_legend = F,
                 col = master_dict,
                 
                 row_names_side = "right",
                 row_split = meta$Tissue, 
                 # show_row_names = F,
                 
                 column_title = "Top IMPACT annotations",
                 column_names_rot = rot,
                 column_names_side = "top", 
                 column_names_centered = T,
                 column_names_gp = gpar(cex=.75),
                 cluster_rows = F ,
                 row_km = 4,
                 # cell_fun = function(j, i, x, y,){grid.text(meta[[variable]], gp=gpar(fontsize=6))},
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%s", meta[i, j]), x, y, 
                             gp = gpar(cex=cex, col=text_color))
                 }
    )
    draw(ha)
    return(ha)
  }
  
  
  annotation_col <- function(meta, variable, width_factor=.6, rot=0, na_fill=NA, cex=.5, text_color="white"){ 
    meta <- data.frame(mat_meta) 
    meta[is.na(mat_meta[[variable]]), variable] <-  na_fill   
    
    ha = Heatmap(meta[variable], name = variable,  
                 width = max_text_width(meta[[variable]])*width_factor, 
                 show_heatmap_legend = F, 
                 
                 row_names_side = "right",
                 row_names_gp = gpar(cex=.5), 
                 row_split = meta[variable],  
                 # show_row_names = F,
                  
                 column_names_rot = rot,
                 column_names_side = "top", 
                 column_names_centered = T,
                 column_names_gp = gpar(cex=.75),
                 cluster_rows = F ,
                 row_km = 4,
                 # cell_fun = function(j, i, x, y,){grid.text(meta[[variable]], gp=gpar(fontsize=6))},
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%s", as.matrix(meta[[variable]])[i, j]), x, y, 
                             gp = gpar(cex=cex, col=text_color))
                 }
                 )
    # dev.off(); 
    # draw(ha)
    return(ha)
  }
  
  
  # png("Data/GWAS/Nalls23andMe_2019/_genome_wide/IMPACT/IMPACT_summary_heatmap.png")
  {
    set.seed(2019)
    ht_list = main_heatmap(mat = mat_meta[,names(snp.groups)], 
                           row_split = mat_meta$Tissue)  +  
      annotation_table(mat_meta) + 
      # annotation_col(mat_meta, "Tissue") +  
      # annotation_col(mat_meta, "CellDeriv") + 
      # annotation_col(mat_meta, "Cell" ) +  
      annotation_col(mat_meta, variable = "TF")
    # ht_list = ht_list +  annotation_col(meta, "TF") 
    draw(ht_list, ht_gap = unit(.5, "mm"), 
         heatmap_legend_side="left")
         # height = unit(30, "cm"))
  }
  # dev.off();
  
  
  
  
  # text_col <- function(mat_meta, variable, width_factor=1.1, rot=45, palette=F, na_fill=NA){
  #   meta <- data.frame(mat_meta)
  #   meta[is.na(mat_meta[[variable]]), variable] <-  na_fill 
  #   if(palette!=F){
  #     color_list =  RColorBrewer::brewer.pal(n = n_distinct(meta[[variable]]), name = palette)
  #   } else{ color_list <- NULL } 
  #   ha <- rowAnnotation(TFBS = anno_text(meta[[variable]],
  #                                        location = 0.5, just = "center", 
  #                                        gp = gpar(cex=1, fill = color_list, col = "black", border = "gray"),
  #                                        width = max_text_width(meta[[variable]])*width_factor))
  #   draw(ha)
  #   return(ha)
  # }
  # library(data.tree)
  # library(treemap)
  # hier_dat <- subset(mat_meta,SNP_group=="Top_Consensus", 
  #                    select=c("Tissue","CellDeriv","Cell","TF")) %>%
  #   dplyr::mutate(pathString=file.path(Tissue, CellDeriv,Cell,TF))
  # # hier_dat <- hier_dat[complete.cases(hier_dat),]
  # 
  # pop <- data.tree::as.Node(hier_dat)
  # plot(population)
  # treemap::treemap(hier_dat%>% dplyr::mutate(size=1),
  #                  index=c("Tissue","CellDeriv","Cell","TF"),
  #                  vSize ="size")
  
  
  # Superheat
  # https://rlbarter.github.io/superheat-examples/Organ/
  # Vignette: https://rlbarter.github.io/superheat/saving-superheatmaps.html
  # superheat(X = mat,
  #           row.dendrogram = T, 
  #           
  #           # heat.lim = c(0,1),
  #           heat.pal = rev( RColorBrewer::brewer.pal(5, "Spectral")),
  #           heat.na.col = "white",
  #           grid.vline.col = "white",
  #           
  #           # left labels
  #           left.label.size = 0.5,
  #           left.label.text.size = 3, 
  #           
  #           # Top plot
  #           yt = colMeans(mat, na.rm = T),
  #           yt.plot.type = "boxplot",
  #           yt.axis.name="Mean IMPACT\nacross loci",
  #           n.clusters.cols = 3,
  #           
  #           yr = meta$x,
  #           )
  
  
  
  # heatmaply::ggheatmap(mat, ) 
  
  # meta_colors <-list(Tissue=setNames(RColorBrewer::brewer.pal(n_distinct(meta$Tissue), "Set3"),unique(meta$Tissue)),
  #                    CellDeriv=setNames(RColorBrewer::brewer.pal(n_distinct(meta$CellDeriv), "Set3"),unique(meta$CellDeriv)),
  #                    Cell=setNames(RColorBrewer::brewer.pal(n_distinct(meta$Cell), "Set3"),unique(meta$Cell)),
  #                    TF=setNames(RColorBrewer::brewer.pal(n_distinct(meta$TF), "Blues"),unique(meta$TF)))
  # 
  # pheatmap::pheatmap(mat = mat, angle_col = 45,
  #                    annotation_row = meta, 
  #                    annotation_legend = F,
  #                    annotation_colors = )
  return(mat_meta)
}



IMPACT.snp_group_boxplot <- function(TOP_IMPACT_all){
  bp <- ggplot(data=TOP_IMPACT_all, aes(x=SNP_group, y=max_IMPACT, fill=SNP_group)) + 
    geom_boxplot(alpha=.5,notch = T, outlier.alpha = 0) + 
    geom_violin(alpha=.5) +
    geom_jitter(alpha=.1,width = .25)
  print(bp)
}
 
  
  
hierarchical_colors <- function(mat_meta){
  # http://hughjonesd.github.io/tweaking-colours-with-the-shades-package.html
  # Start by assigning group colors to each Tissue
  meta <- data.frame(mat_meta)
  Tissue_dict <- setNames(pals::cols25(n = n_distinct(meta$Tissue)), #RColorBrewer::brewer.pal(n=n_distinct(meta$Tissue), "Dark2"),
                          unique(meta$Tissue))
  # color_list <- Tissue_dict[meta$Tissue]
  # names(color_list) <- meta[[variable]]
  get_new_colors <- function(mat_meta, 
                             group="Tissue",
                             variable="CellDeriv",
                             dict){ 
    counts <- mat_meta %>% 
      dplyr::group_by(eval(parse(text=group))) %>% 
      dplyr::summarise_at(.vars=variable,
                          .funs=c(count=n_distinct) ) %>%
      `colnames<-`(c(group,"count"))
    
    
    new_color_dict <- lapply(1:nrow(counts), function(i){
      row <- counts[i,]
      # print(row[[group]])
      start_color <- dict[row[[group]] ]
      
      incrementer <- rev(seq(.1, .8, by=.1) )#seq(.8,.1,length.out = row$count)
      colurz <- lapply(1:row$count, function(j){
        if(is.na(j)){ return(setNames("gray",NA))
        } else {shades::brightness(unname(start_color), incrementer[j])[[1]]   } 
      }) %>% unlist()
      colur_namez <- mat_meta[mat_meta[[group]]==row[[group]] & !is.na(mat_meta[[group]]),][[variable]] %>% unique()
      names(colurz) <- colur_namez
      return(colurz)
    }) %>% unlist()
    return(new_color_dict)
  } 
  CellDeriv_dict <- get_new_colors(mat_meta,group = "Tissue", variable="CellDeriv", dict=Tissue_dict)
  Cell_dict <- get_new_colors(mat_meta,group = "Tissue", variable="Cell", dict=Tissue_dict)
  TF_dict <- get_new_colors(mat_meta,group = "Tissue", variable="TF", dict=Tissue_dict)
  
  master_dict <- c(Tissue_dict, CellDeriv_dict, Cell_dict, TF_dict)
  return(master_dict)
}

  
  
  
  
IMPACT.plot_top_annotations <- function(){
  # merged_DT <- quick_merged_DT()
  # consensus.snps <- subset(merged_DT, Consensus_SNP)$SNP
  
  
  
  # TILE PLOT 
  # ggplot(top_impact, aes(x=Locus, y=TF, fill=IMPACT_score)) +
  #   geom_tile() + 
  #   facet_grid(facets = Tissue ~ ., 
  #              scales = "free_y", 
  #              space = "free_y") +
  #   theme_bw() +
  #   scale_x_discrete(position = "top") +
  #   scale_y_discrete(position = "right") + 
  #   theme(axis.text.x = element_text(angle=45, hjust=0), 
  #         strip.placement = "outside", 
  #         strip.text.y = element_text(angle=0))
  
  # library(ggsci)
  # barplot_layer <- function(ANNOT_MELT, 
  #                           ytext=F, 
  #                           title=NULL,
  #                           xlabel="Mean IMPACT score", 
  #                           palette=NULL, 
  #                           snp.filter="!is.na(IMPACT_score)"){
  #   
  #   top_impact <- get_mean_IMPACT(ANNOT_MELT, snp.filter=snp.filter)
     
  #   bp <- ggplot(top_impact, aes(x=meanIMPACT, y=Locus, fill=meanIMPACT)) + 
  #     geom_col() + 
  #     theme_bw() + 
  #     labs(fill="Mean IMPACT score", x=xlabel, title=title) +
  #     scale_fill_continuous(limits=c(0,1), breaks=c(0,.5,1)) +
  #     scale_x_continuous(limits = c(0,1), breaks = c(0,.5,1) ,position = "top") +  
  #     scale_y_discrete(drop=FALSE) +  
  #     geom_text(aes(x=meanIMPACT, label = SNPs, color=meanIMPACT),
  #               size=3, show.legend = F, nudge_x = .01, hjust=0) +
  #     theme(legend.position = "left") 
  #   if(ytext==F){
  #     bp <- bp + 
  #       ylab(NULL) + 
  #       theme(axis.text.y = element_blank(),
  #             plot.title = element_text(hjust = .5),
  #             legend.position = "None")
  #   }
  #   return(bp)
  # }
  # 
  # top_impact <- get_mean_IMPACT(ANNOT_MELT, snp.filter="topConsensus")
  # meta_table <- top_impact[,c("Locus","Tissue","CellDeriv","Cell","TF")]
  # # meta_table[, c("Tissue","CellDeriv","Cell")][is.na(meta_table[, c("Tissue","CellDeriv","Cell")])] <- "."
  # meta_table <- meta_table %>% 
  #   dplyr::mutate_at(.vars = c("Tissue"),#,"CellDeriv","Cell"), 
  #                    .funs = tolower) %>% 
  #   dplyr::mutate(Tissue = gsub("^gi$","GI",Tissue),
  #                 Cell_Type = paste0(ifelse(CellDeriv=="mesendoderm","hESC derived mesendodermal cells",CellDeriv), 
  #                                    ifelse(Cell %in% c("hESC derived mesendodermal cells","."),"",paste0(" (",Cell,")") ))  )
 
  # meta_table <- meta_table[,c("Locus","Tissue","Cell_Type","TF")]
  # meta_table <- data.table::melt.data.table(data.table::data.table(meta_table),
  #                                           id.vars = "Locus")
  # #
  # g_table <- ggplot(data = meta_table, aes(x="1", y=Locus)) +
  #   geom_tile(aes(fill=value)) +
  #   geom_text(aes(label=value), size=3) +
  #   facet_grid(facets = . ~ variable, space = "free_x") +
  #   labs(title="Top IMPACT annotation per top Consensus SNP", x=NULL, y=NULL) +
  #   theme_bw() +
  #   theme(axis.text = element_blank(),
  #         plot.title = element_text(hjust=.5)) 
  #   
  # metadata_tiles <- function(variable="Tissue", 
  #                            palette="Set3",
  #                            show.legend=F){ 
  #   tiles <- ggplot(data=meta_table[,c("Locus",variable)], aes(x="1", y=Locus)) +
  #     geom_tile(aes(fill=eval(parse(text=variable))), show.legend = show.legend, alpha=.5) +
  #     # scale_fill_discrete(na.value = "transparent") +
  #     scale_fill_brewer(palette = palette, na.value="transparent") +
  #     geom_text(aes(label=eval(parse(text=variable))), color="black") +
  #     labs(x=NULL, y=NULL, title=variable, fill=variable) + 
  #     theme_bw() +
  #     theme(axis.text = element_blank(),
  #           plot.title = element_text(hjust=.5), 
  #           legend.position = "bottom")
  #   return(tiles)
  # }
  # mt_Tissue <- metadata_tiles(variable="Tissue", palette="Set3")
  # mt_CellDeriv <- metadata_tiles(variable="CellDeriv", palette="Spectral")
  # mt_Cell <- metadata_tiles(variable="Cell", palette="Paired")
  # mt_Tissue + mt_CellDeriv + mt_Cell
  
 
  
  library(patchwork)
  bpl <- barplot_layer(ANNOT_MELT, snp.filter = "leadSNP==T", ytext=T, title = "Lead GWAS SNPs") + 
    barplot_layer(ANNOT_MELT, snp.filter = "Support>0", title = "UCS SNPs") + 
    barplot_layer(ANNOT_MELT, snp.filter = "Consensus_SNP", title = "Consensus SNPs") +  
    barplot_layer(ANNOT_MELT, snp.filter = "topConsensus", title = "Top Consensus SNP") +
    g_table +
    patchwork::plot_layout(nrow = 1, widths = c(rep(.3,4),1))  
  print(bpl)
  
  if(save_path!=F){
    save_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide/IMPACT"
    dir.create(save_path,showWarnings = F, recursive = T)
    ggsave(file.path(save_path,"IMPACT_SNPgroups_summary.png"), dpi = 400, height = 12, width = 17)
  }
  
}
 


 
IMPACT.compute_enrichment <- function(annot_melt, locus=NULL){  
  sum.IMPACT <- sum(annot_melt$IMPACT_score, na.rm=T)
  len.SNPs <- n_distinct(annot_melt$SNP, na.rm = T)
  # SNP.groups <- list("leadGWAS"=subset(annot_melt, leadSNP),
  #                    "UCS"=subset(annot_melt, Consensus_SNP),
  #                    "ABF_CS"=subset(annot_melt, ABF.Credible_Set>0),
  #                    "FINEMAP_CS"=subset(annot_melt, FINEMAP.Credible_Set>0),
  #                    "SUSIE_CS"=subset(annot_melt, SUSIE.Credible_Set>0),
  #                    "POLYFUN_CS"=subset(annot_melt, POLYFUN_SUSIE.Credible_Set>0),
  #                    "Consensus"=subset(annot_melt, Consensus_SNP))
  # enrich <- lapply(names(SNP.groups), function(snp.group){
  #   print(snp.group)
  #   e <- SNP.groups[[snp.group]] %>% 
  #     dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
  #     dplyr::summarise(enrichment = (sum(IMPACT_score, na.rm = T) / sum.IMPACT) /
  #                        (n_distinct(SNP, na.rm = T) / len.SNPs) ) %>% 
  #     dplyr::arrange(-enrichment) %>% 
  #     data.table::data.table()
  #   e <- cbind(SNP.group=snp.group, e)
  #   return(e)
  # }) %>% data.table::rbindlist()
  annot_melt[is.na(annot_melt$IMPACT_score),"IMPACT_score"] <- 0

  SNP.groups <- list(
    "leadGWAS" = annot_melt %>% 
      dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
      dplyr::summarise(enrichment = (sum(IMPACT_score[leadSNP], na.rm = T) / sum(IMPACT_score, na.rm = T)) /
                         (n_distinct(SNP[leadSNP], na.rm = T) / n_distinct(SNP, na.rm = T)) ),
    "UCS" = annot_melt %>% 
      dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
      dplyr::summarise(enrichment = (sum(IMPACT_score[Support>0], na.rm = T) / sum(IMPACT_score, na.rm = T)) /
                         (n_distinct(SNP[Support>0], na.rm = T) / n_distinct(SNP, na.rm = T)) ),
    "ABF_CS" = annot_melt %>% 
      dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
      dplyr::summarise(enrichment = (sum(IMPACT_score[ABF.Credible_Set>0], na.rm = T) / sum(IMPACT_score, na.rm = T)) /
                         (n_distinct(SNP[ABF.Credible_Set>0], na.rm = T) / n_distinct(SNP, na.rm = T)) ), 
    "FINEMAP_CS" = annot_melt %>% 
      dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
      dplyr::summarise(enrichment = (sum(IMPACT_score[FINEMAP.Credible_Set>0], na.rm = T) / sum(IMPACT_score, na.rm = T)) /
                         (n_distinct(SNP[FINEMAP.Credible_Set>0], na.rm = T) / n_distinct(SNP, na.rm = T)) ), 
    "SUSIE_CS" = annot_melt %>% 
      dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
      dplyr::summarise(enrichment = (sum(IMPACT_score[SUSIE.Credible_Set>0], na.rm = T) / sum(IMPACT_score, na.rm = T)) /
                         (n_distinct(SNP[SUSIE.Credible_Set>0], na.rm = T) / n_distinct(SNP, na.rm = T)) ),
    "POLYFUN_CS" = annot_melt %>% 
      dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
      dplyr::summarise(enrichment = (sum(IMPACT_score[POLYFUN_SUSIE.Credible_Set>0], na.rm = T) / sum(IMPACT_score, na.rm = T)) /
                         (n_distinct(SNP[POLYFUN_SUSIE.Credible_Set>0], na.rm = T) / n_distinct(SNP, na.rm = T)) ), 
    "Consensus" = annot_melt %>% 
      dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
      dplyr::summarise(enrichment = (sum(IMPACT_score[Consensus_SNP], na.rm = T) / sum(IMPACT_score, na.rm = T)) /
                         (n_distinct(SNP[Consensus_SNP], na.rm = T) / n_distinct(SNP, na.rm = T)) ) 
  )
  enrich <- data.table::rbindlist(SNP.groups, idcol = "SNP.group") %>% dplyr::arrange(-enrichment)
  enrich <- cbind(Locus=locus, enrich)
  enrich$TF <- factor(enrich$TF, ordered = T)
  enrich$SNP.group <- factor(enrich$SNP.group, levels=names(SNP.groups), ordered = T)
  return(enrich)
}



IMPACT.iterate_enrichment <- function(gwas_paths,
                                      annot_baseURL="../../data/IMPACT/IMPACT707/Annotations"){
  # gwas_paths <- list.files(path = "./Data/GWAS/Nalls23andMe_2019", pattern = "Multi-finemap_results.txt", recursive = T, full.names = T)
  # no_no_loci <- c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1")
  # gwas_paths <- gwas_paths[!basename(dirname(dirname(gwas_paths))) %in% no_no_loci]
 
  # ENRICH <- lapply(gwas_paths, function(x){
  #   locus <- basename(dirname(dirname(x)))
  #   message(locus)
  #   enrich <- NULL
  #   try({
  #     subset_DT <- data.table::fread(x, nThread = 4)
  #     if(!"Locus" %in% colnames(subset_DT)){
  #       subset_DT <- cbind(Locus=locus, subset_DT) 
  #     }
  #     subset_DT <- find_consensus_SNPs(finemap_DT = subset_DT)
  #     annot_melt <- IMPACT.get_annotations(baseURL = annot_baseURL, 
  #                                          subset_DT = subset_DT, 
  #                                          nThread = 4) 
  #     enrich <- IMPACT.compute_enrichment(annot_melt = annot_melt,
  #                                         locus = locus)
  #   }) 
  #   return(enrich)
  # }) %>% data.table::rbindlist(fill=T)
  
  ENRICH <- lapply(unique(ANNOT_MELT$Locus), function(locus){
    message("+ IMPACT:: Locus = ",locus)
    annot_melt <- subset(ANNOT_MELT, Locus==locus)
    enrich <- IMPACT.compute_enrichment(annot_melt = annot_melt,
                                        locus = locus)
    return(enrich)
  }) %>% data.table::rbindlist()
  
  
  # ENRICH
  return(ENRICH)
}

 

IMPACT.plot_enrichment <- function(ENRICH){
  enrich_dat <- ENRICH %>%
    # subset(!is.na(enrichment) & enrichment>=1) %>% 
    dplyr::group_by(SNP.group, Tissue, CellDeriv) %>%
    dplyr::top_n(n=1, wt=enrichment)
  mean_dat <- enrich_dat %>% 
    dplyr::group_by(SNP.group) %>% 
    dplyr::summarise(mean.enrichment= mean(enrichment))
  
  ep <- ggplot() +
    # geom_boxplot(data=enrich_dat, aes(x=SNP.group, y=enrichment)) +
    geom_col(data=mean_dat, aes(x=SNP.group, y=mean.enrichment, fill=SNP.group), position = "dodge") +
    geom_jitter(data=enrich_dat, aes(x=SNP.group, y=enrichment),
                size=1, alpha=.25, color="cyan") +
    geom_hline(yintercept = 1, linetype="dashed", alpha=.8) + 
    theme_bw() +
    theme(strip.text = element_text(angle=0),
          axis.text.x = element_text(angle=45, hjust=1))
  print(ep)
  
  
  
  ep <- ggplot(ENRICH, aes(x=Tissue, y=enrichment, fill=SNP.group)) +
    geom_violin(position = "dodge") + 
    geom_jitter(size=1, alpha=.25, color="cyan") +
    geom_hline(yintercept = 1, linetype="dashed", alpha=.8) +
    facet_grid(facets = SNP.group ~ .,#Tissue + CellDeriv,
               switch = "y", space = "free_x",
               scales = "free_x") +
    theme_bw() +
    theme(strip.text = element_text(angle=0),
          axis.text.x = element_text(angle=45, hjust=1))
  print(ep)
  
  ep <- ggplot(enrich_dat, aes(x=TF, y=SNP.group, fill=enrichment)) +
    geom_col(position = "dodge") +
    geom_hline(yintercept = 1, linetype="dashed", alpha=.8) +
    facet_grid(facets = . ~ Tissue,# + CellDeriv, 
               switch = "y", space = "free_x", scales = "free_x") +
    theme_bw() + 
    theme(
      # strip.text = element_text(angle=0), 
          axis.text.x = element_text(angle=45, hjust=1))
  print(ep)
}


IMPACT.plot_impact_score <- function(annot_melt, 
                                     save_path=F,
                                     show_plot=T){
  
  # subset_DT <- data.table::fread("Data/GWAS/Nalls23andMe_2019/CD19/Multi-finemap/Multi-finemap_results.txt")
  # subset_DT <- find_consensus_SNPs(subset_DT)
  # subset_DT$Locus <- "CD19"
  # annot_melt <- IMPACT.get_annotations(baseURL = "/Volumes/Steelix/IMPACT/IMPACT707/Annotations", subset_DT = subset_DT)
  
  library(patchwork)
  library(ggridges) 
  
  annot_melt$Mb <- annot_melt$POS/1000000
  # Get the SNP w/ the top impact score for each annotation
  annot_top <- annot_melt %>% 
    dplyr::group_by(Tissue, Cell, CellDeriv, TF) %>% 
    top_n(n=1, wt=IMPACT_score) %>% 
    data.table::data.table()
  # subset(annot_top,SNP %in% unique(subset(annot_melt, Consensus_SNP)$SNP))
  # annot_top
  
  
  # Reduce to smaller df to make plotting faster
  finemap_cols <- grep("*.PP$|*.Credible_Set$",colnames(annot_melt),value=T)
  annot_snp <- subset(annot_melt, select=c("SNP","CHR","POS","Mb","P","Consensus_SNP","leadSNP","Support",finemap_cols)) %>% unique()
  annot_snp <- dplyr::mutate(annot_snp, SNP.Group = ifelse(Consensus_SNP,"Consensus SNP",ifelse(leadSNP,"Lead GWAS SNP",ifelse(Support>0,"Credible Set SNP",NA))))
  labelSNPs <- construct_SNPs_labels(DT = annot_snp, lead=T, method=T, consensus=T)
  leader_SNP <- subset(labelSNPs, type=="Lead SNP")
  CS_set <- subset(labelSNPs, type=="Credible Set")
  # ggb <- GGBIO.plot(finemap_DT = annot_snp, LD_matrix = LD_matrix, 
  #            XGR_libnames = NULL, 
  #            save_plot=F,
  #            Nott_sn_epigenome=F) 
  # GWAS row
  gwas <- ggplot(annot_snp, aes(x=Mb, y=-log10(P), color=-log10(P))) +
    geom_point(size=1) +
    geom_point(data=leader_SNP, pch=18, fill=NA, size=2.5, color=leader_SNP$color) +
    # Green rings aronud Credible Set SNPs
    geom_point(data=CS_set, pch=21, fill=NA, size=2.5, color=CS_set$color, stroke=1, alpha=0.8) +
    ### Background color label
    ggrepel::geom_label_repel(data=labelSNPs,
                              aes(label=SNP),
                              color=NA,
                              # nudge_x = .5,
                              fill="black",
                              box.padding = .25,
                              label.padding = .25,
                              label.size=NA,
                              alpha=.6,
                              seed = 1,
                              size = 3,
                              min.segment.length = 1) +
    ### Foreground color label
    ggrepel::geom_label_repel(data=labelSNPs,
                              aes(label=SNP),
                              color=labelSNPs$color,
                              segment.alpha = .5,
                              # nudge_x = .5,
                              box.padding = .25,
                              label.padding = .25,
                              segment.size = 1,
                              fill = NA,
                              alpha=1,
                              seed = 1,
                              size = 3) + 
    ylim(c(0,max(-log10(annot_snp$P)))*1.1) +
    geom_vline(xintercept = unique(subset(labelSNPs, Consensus_SNP)$Mb), color="goldenrod2",
               alpha=1, size=.3, linetype='solid') +
    geom_vline(xintercept = unique(subset(labelSNPs, leadSNP)$Mb), color="red",
               alpha=1, size=.3, linetype='solid') +
    theme_bw()
  
  # finemaping rows
  finemap <- ggplot(annot_snp, aes(Mb, y=POLYFUN_SUSIE.PP, color=POLYFUN_SUSIE.PP)) +
    geom_point(size=1) +  
    scale_color_viridis_c( breaks=c(0,.5,1)) +
    geom_point(data=subset(annot_snp,POLYFUN_SUSIE.Credible_Set>0), pch=21, fill=NA, size=2.5, 
               color="green3", stroke=1, alpha=0.8) +
    ylim(c(0,1.1)) +
    geom_vline(xintercept = unique(subset(labelSNPs, Consensus_SNP)$Mb), color="goldenrod2",
               alpha=1, size=.3, linetype='solid') +
    geom_vline(xintercept = unique(subset(labelSNPs, leadSNP)$Mb), color="red",
               alpha=1, size=.3, linetype='solid') + 
    theme_bw()
  # print(finemap)
  
  # IMPACT rows
  impact <- ggplot(subset(annot_melt, IMPACT_score>0.5), 
                   aes(x=Mb, y=IMPACT_score, color=TF)) + 
    geom_point(show.legend = F) +
    # geom_col(position = "identity", show.legend = T) +
    facet_grid(facets = Tissue ~ ., switch = "y") +
    # ggridges::geom_ridgeline(aes(height = IMPACT_score), na.rm = T, size=.1, show.legend = F) +
    # ggridges::theme_ridges() +
    theme_bw() + 
    labs(y="IMPACT score per tissue") + 
    theme(strip.text.y = element_text(angle = 0), 
          panel.grid = element_blank(), 
          axis.title.y = element_text(vjust = .5))
  # impact 
  
  impact_plot <- gwas + finemap + impact + patchwork::plot_layout(ncol = 1, heights = c(.2,.2,1)) +
    geom_vline(xintercept = unique(subset(labelSNPs, Consensus_SNP)$Mb), color="goldenrod2",
               alpha=1, size=.3, linetype='solid') +
    geom_vline(xintercept = unique(subset(labelSNPs, leadSNP)$Mb), color="red",
               alpha=1, size=.3, linetype='solid')
  # impact_plot
  
  if(show_plot){print(impact_plot)}
  
  if(save_path!=F){
    # save_path="./Data/GWAS/Nalls23andMe_2019/LRRK2/IMPACT/LRRK2_IMPACT_plot.png"
    dir.create(dirname(save_path), showWarnings = F, recursive = T)
    printer("IMPACT:: Saving plot ==>",save_path)
    ggsave(save_path, plot=impact_plot, height=10, width=10)
  }
  return(impact_plot)
}


##### xxxxxxxxxxxxxxxxxxxxxxxxx TRACK PLOTS xxxxxxxxxxxxxxxxxxxxxxx ####
IMPACT.plot_impact_score_compare <- function(){
  
  subset_DT <- lapply(c("CD19","TRIM40","NUCKS1","LRRK2","MED12L","MEX3C"), function(locus){
     dat <- data.table::fread(file.path("Data/GWAS/Nalls23andMe_2019", locus,
                                 "Multi-finemap/Multi-finemap_results.txt")) 
     dat$Locus <- locus 
     dat <- assign_lead_SNP(dat) 
     return(dat)
  }) %>% data.table::rbindlist(fill = T)
  subset_DT <- find_consensus_SNPs(subset_DT)
  subset_DT <- find_topConsensus(subset_DT)
  
  
  annot_melt <- IMPACT.iterate_get_annotations(subset_DT, 
                                               IMPACT_score_thresh=0,
                                               baseURL = "../../data/IMPACT/IMPACT707/Annotations",
                                               # baseURL="/Volumes/Steelix/IMPACT/IMPACT707/Annotations",
                                               all_snps_in_range=T, 
                                               top_annotations_only=T,
                                               snp.filter = "Consensus_SNP")
  # data.table::fwrite(annot_melt,"Data/GWAS/Nalls23andMe_2019/_genome_wide/IMPACT/annot_melt_topAnnot_subset.csv")
  annot_melt <- data.table::fread("/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/IMPACT/annot_melt_topAnnot_subset.csv.gz", nThread = 4)
  annot_melt <- IMPACT.postprocess_annotations(annot_melt) 
  top_impact <- IMPACT.get_top_annotations(ANNOT_MELT = annot_melt, 
                                           snp.filter = "Consensus_SNP",
                                           top_annotations = 1,
                                           force_one_annot_per_locus = T)
  # top_impact <- subset(top_impact, Locus %in% c("CD19","NUCKS1","MEX3C"))
  
  # When data was originally merged, kept all subset_DT rows
  annot_sub <- subset(annot_melt,
                       Locus %in% top_impact$Locus & # somehow this happens sometimes...
                       Tissue %in% top_impact$Tissue & 
                       CellDeriv %in% top_impact$CellDeriv & 
                       Cell %in% top_impact$Cell & 
                       TF %in% top_impact$TF) %>% 
    dplyr::mutate(Mb=POS/1000000)
  
  
  
  # Convert to GRange object
  # gr.snp_CHR <- biovizBase::transformDfToGr(subset_DT, seqnames = "CHR", start = "POS", end = "POS")
  # gene_model <- transcript_model_track(gr.snp_CHR,
  #                                      max_transcripts=1,
  #                                      show.legend=T)
  # # Remove any pseudogenes
  # db.gr <- db.gr[grep("*pseudogene*",db.gr$tx_biotype, invert = T, fixed = F),]
   
  # autoplot(edb,
  #          # Have to limit (can only handle depth < 1000)
  #          which = db.gr,
  #          names.expr = "gene_name",
  #          aes(fill=gene_name, color=gene_name),
  #          show.legend=show.legend)  +
  #   theme_classic() +
  #   theme(strip.text.y = element_text(angle = 0),
  #         strip.text = element_text(size=9),
  #         legend.text = element_text(size=5),
  #         legend.key.width=unit(.1,"line"),
  #         legend.key.height=unit(.1,"line")) +
  #   guides(fill=guide_legend(override.aes = list(size=1), ncol=4),
  #          color=guide_legend(override.aes = list(size=1), ncol=4),
  #          size=.5)
  add_snp_labels <- function(p, annot_sub, y_var="-log10(P)"){
    label_tags <- construct_SNPs_labels(subset_DT = annot_sub, lead=T, method=T, consensus=T,
                                        remove_duplicates = F) 
    label_tags_unique <- construct_SNPs_labels(subset_DT = annot_sub, lead=T, method=T, consensus=T,
                                        remove_duplicates = T) 
    p <- p +
      # Circles
      geom_point(data = label_tags, aes(x=Mb, y=eval(parse(text=y_var)) ),
                 color=label_tags$color,
                 shape=label_tags$shape, 
                 size=label_tags$size, 
                 stroke=1) + 
      # Labels
      ggrepel::geom_label_repel(data=label_tags_unique,
                                aes(label=SNP),
                                color=NA,
                                # nudge_x = .5,
                                fill="black",
                                box.padding = .25,
                                label.padding = .25,
                                label.size=NA,
                                alpha=.6,
                                seed = 1,
                                size = 3,
                                min.segment.length = 1) +
      ### Foreground color label
      ggrepel::geom_label_repel(data=label_tags_unique,
                                aes(label=SNP),
                                color=label_tags_unique$color,
                                segment.alpha = .5,
                                # nudge_x = .5,
                                box.padding = .25,
                                label.padding = .25,
                                segment.size = 1,
                                fill = NA,
                                alpha=1,
                                seed = 1,
                                size = 3,
                                min.segment.length = 1)  
    return(p)
  }
   
  library(patchwork)
  gp1 <- ggplot(data = annot_sub, aes(x=Mb, y=-log10(P), color=-log10(P))) + 
    geom_point(size=.5, alpha=.25) +
    theme_bw() +
    facet_grid(facets = . ~ Locus,
               scales = "free_x", switch="y") +
    labs(title="Fine-mapping results overlaid on GWAS and IMPACT results", x=NULL, y="GWAS -log10(P)") + 
    theme(axis.text.x = element_blank(),
          plot.title = element_text(hjust=.5),
          strip.background = element_rect(fill="transparent"))   
  
  # Get the density of SNPs with a high IMPACT_score
  density_data <- annot_sub
  density_data[density_data$IMPACT_score<.1, "IMPACT_score"] <- NA
  density_data <- subset(density_data, !is.na(IMPACT_score))
  
  gp2 <- ggplot(data = annot_sub, aes(x=Mb, y=IMPACT_score, color=IMPACT_score)) +   
    geom_density(data = density_data, 
                 aes(x=Mb, y=..scaled..), 
                 color="transparent", alpha=.25, fill="green",
                 show.legend = F, adjust=.5) + 
      geom_point(size=.5, alpha=.5) + 
      scale_color_viridis_c() + 
      # geom_smooth(data=density_data, 
      #             stat="smooth", method = 'loess', 
      #             aes(x=Mb),
      #             span=.2,
      #             se=F) +
      labs(y="IMPACT score") +
      facet_grid(facets =  . ~ Locus,
                 scales = "free_x", switch = "y") +  
      theme_bw() + 
      theme(strip.background = element_blank(),
            strip.text.x = element_blank())   
  # gp2
  
  # Combine plots
  cp <- add_snp_labels(gp1, annot_sub, y_var = "-log10(P)") + 
    add_snp_labels(gp2, annot_sub, y_var = "IMPACT_score") + 
    patchwork::plot_layout(ncol = 1)
  cp
  
  ggsave("Data/GWAS/Nalls23andMe_2019/_genome_wide/IMPACT/IMPACT_example_tracks.png", 
         plot=cp, dpi=400, height=6, width=12)
}



IMPACT.get_ldscores <- function(chrom=NULL, 
                                subset_DT=NULL,
                                nThread=4){ 
  warning("LDSCores do not include any SNPs with MAF<0.5%, as they are restricted to HapMap3 SNPs. \
This may affect subsequent analyss (e.g. fine-mapping).")
  if(!is.null(subset_DT)){
    chrom <- subset_DT$CHR[1]
  }
  baseURL <- "https://github.com/immunogenomics/IMPACT/raw/master/IMPACT707/LDscores"
  URL <- file.path(baseURL, paste0("IMPACT707_EAS_chr",chrom,".l2.ldscore.gz"))
  ldscore <- data.table::fread(URL, nThread = nThread)
  
  if(!is.null(subset_DT)){ 
    ldscore_merge <- data.table:::merge.data.table(data.table::data.table(subset_DT),
                                                   ldscore,
                                                   by.x = c("SNP","CHR","POS"), 
                                                   by.y = c("SNP","CHR","BP"),
                                                   all = F)  
    return(ldscore_merge)
  } else {
    return(ldscore)
  } 
}

