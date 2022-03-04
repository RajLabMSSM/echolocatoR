
### Summarise fine-mapping results from all loci ###
# merged_dat <- merge_finemapping_results(consensus_thresh = 3)
# (merged_dat %>% subset(P<5e-8) %>%dplyr::group_by(Gene) %>% count())$n %>% mean()


lead.SNP.coords <- function(consensus_thresh=2){
  annot <- readxl::read_excel("./Data/annotated_finemapping_results.xlsx")
  annot <- find_consensus_SNPs(annot,
                               consensus_thresh = consensus_thresh,
                               verbose = F)
  annot[is.na(annot$mean.PP),"mean.PP"] <-0

  annot.sub <- subset(annot, leadSNP==T, select=c(Gene, SNP, CHR, POS)) %>%
    dplyr::rename(lead.SNP=SNP) %>%
    dplyr::mutate(min.POS=POS - 1e+06, max.POS=POS + 1e+06)
  data.table::fwrite(annot.sub, "./Data/lead.SNP.coords.csv", sep=",")
  annot[is.na(annot)] <-0

  SNPgroup.summary <- function(DF, group.name=""){
    n.total = nrow(DF)
    n.per.locus <- (DF %>% group_by(Gene) %>% tally())$n %>% mean()
    means <- DF[,c("P","Effect","mean.PP","MAF")] %>% dplyr::summarise_all(function(x){abs(mean(x, na.rm = T))})
    means <- cbind(`SNP Group`=group.name, means, `SNPs / locus` = n.per.locus, `Total SNPs` = n.total)
    # Count distribution
    mean.size.statCS <- table((statCS %>% tally())$n )
    print(mean.size.statCS)
    return(means)
  }

  # All SNPs
  allSNPs.means <- SNPgroup.summary(annot, group.name="All GWAS")
  # Nom sig GWAS SNPs
  nomSigSNPs <- subset(annot, P<0.05)
  nomSigSNPs.means <- SNPgroup.summary(nomSigSNPs, group.name="nom. sig. GWAS")
  # FDR sig GWAS SNPs
  fdrSigSNPs <- subset(annot, P<5e-8)
  fdrSigSNPs.means <- SNPgroup.summary(fdrSigSNPs, group.name="FDR sig. GWAS")

  #CS stats

  # Lead GWAS SNPs
  leadSNPs <- subset(annot, leadSNP)
  leadSNP.means <- SNPgroup.summary(leadSNPs, group.name="Lead GWAS SNP")
  #CS stats
  ## Statistical FM
  statCS <- annot %>%
    group_by(Gene) %>%
    subset(SUSIE.CS>0 | FINEMAP.CS > 0| PAINTOR.CS >0)
  statCS.means <- SNPgroup.summary(statCS, group.name="Statistical.CS")
  # Functional FM
  funcCS <- annot %>% group_by(Gene) %>% subset(PAINTOR.CS>0)
  funcCS.means <- SNPgroup.summary(funcCS, group.name="Functional.CS")
  # Consensus stats
  consensus <- annot %>% group_by(Gene) %>% subset(Consensus_SNP)
  consensus.means <- SNPgroup.summary(consensus, group.name="Consensus")
  percent.loci.w.consensus <- length(unique(consensus$Gene)) / length(unique(annot$Gene))

  # Merged
  merged.means <- rbind(allSNPs.means, nomSigSNPs.means, fdrSigSNPs.means, leadSNP.means,statCS.means, funcCS.means, consensus.means)
  merged.means[,-c(1:2)] <-round(merged.means[,-c(1:2)],3)
  merged.means$P <- formatC(merged.means$P , format = "e", digits = 3)
  data.table::fwrite(merged.means, "./Data/GWAS/Nalls23andMe_2019/_genome_wide/SNP.summary.csv")

  library(patchwork)
  bins=150
  x.var="Effect"
  alpha=1
  ggplot() +
    geom_histogram(data=leadSNPs, aes(x=eval(parse(text=x.var)), fill="Lead GWAS SNPs"), alpha=alpha,   fill="red", bins = bins) +
    labs(title = "Lead GWAS SNPs", x=x.var) +
    theme_classic() +
    xlim(min(annot[[x.var]]), max(annot[[x.var]])) +

    ggplot() +
    geom_histogram(data = statCS, aes(x=eval(parse(text=x.var)), fill="Statistical Credible Set"), alpha=alpha,  fill="green", bins = bins) +
    labs(title = "Statistical Credible Set", x=x.var) +
    theme_classic() +
    xlim(min(annot[[x.var]]), max(annot[[x.var]])) +

    ggplot() +
    geom_histogram(data = funcCS, aes(x=eval(parse(text=x.var)), fill="Functional Credible Set"), alpha=alpha,   fill="green4", bins = bins) +
    labs(title = "Functional Credible Set", x=x.var) +
    theme_classic() +
    xlim(min(annot[[x.var]]), max(annot[[x.var]])) +

    ggplot() +
    geom_histogram(data = consensus, aes(x=eval(parse(text=x.var)), fill="Consensus SNPs"), alpha=alpha,  fill="goldenrod2", bins = bins) +
    labs(title = "Consensus SNPs", x=x.var) +
    theme_classic() +
    xlim(min(annot[[x.var]]), max(annot[[x.var]])) +

    patchwork::plot_layout(ncol = 1)

  return(merged.means)
}





compare_finemapping_methods <- function(dataset="./Data/GWAS/Nalls23andMe_2019"){
  FM_orig <- merge_finemapping_results(minimum_support = 0,
                                       top_CS_only=F,
                                       dataset = dataset,
                                       exclude_methods = NULL)
  # counts <- (FM_orig %>% group_by(Gene) %>% count())
  # plotrix::std.error(counts$n)
  # min(counts$n)
  # max(counts$n)

  # Remove loci that were manually added
  ## Also remove TRIM40 for now bc I'm having issues getting its LD
  FM_all <- subset(FM_orig, !(Gene %in%c("ATG14","SP1","LMNB1","ATP6V0A1","TRIM40")) )
  # Identify any un-finemapped loci
  FM_all %>% group_by(Gene) %>% summarise(leadGWAS=sum(leadSNP==T)) %>% arrange(leadGWAS)

  # Proportion of CS SNPs that are the leadSNP (by tool)
  if(top_CS_only){
    dat <- FM_all %>%
      summarise_at(.vars = vars(ends_with(".CS")),
                   .funs = list(CS.sum=sum(.==1, na.rm = T),
                                leadGWAS.sum=sum(leadSNP==T, na.rm = T),
                                leadGWAS.CS.sum=sum(leadSNP==T & .==1, na.rm = T),
                                PROP.SNPS=sum(leadSNP==T & .==1, na.rm = T)/sum(.==1, na.rm = T)) )
    colnames(dat) <- gsub("\\.CS","", colnames(dat))

  } else {
    dat <- FM_all %>%
      summarise_at(.vars = vars(ends_with(".CS")),
                   .funs = list(CS.sum=sum(.>0, na.rm = T),
                                leadGWAS.sum=sum(leadSNP==T, na.rm = T),
                                leadGWAS.CS.sum=sum(leadSNP==T & .>0, na.rm = T),
                                PROP.SNPS=sum(leadSNP==T & .>0, na.rm = T)/sum(.>0, na.rm = T)) )
    colnames(dat) <- gsub("\\.CS","", colnames(dat))
  }

  # Proportion of loci in which the CS contains the leadSNP (by tool)
  # dat3.1 <- FM_all %>% group_by(Gene) %>%
  #   summarise_at(vars(ends_with(".CS")),
  #                funs(leadGWAS.CS.size=sum(leadSNP==T & .>0, na.rm = T)))
  # dat3.2 <- FM_all %>% group_by(Gene) %>%
  #   summarise_at(vars(ends_with(".CS")),
  #                funs(leadGWAS.CS.size=sum(.>0, na.rm = T)))
  # colSums(dat3.1[,-1]) / colSums(dat3.2[,-1])

  # Number of CS SNPs (by tool)
  dat2 <- FM_all %>% group_by(Gene) %>% summarise_at(.vars = vars(ends_with(".CS")),
                                                     .funs = list(CS=sum(.,na.rm=T)))
  colnames(dat2) <- gsub("\\.CS","", colnames(dat2))
  # Proportion of loci with at least one CS SNP (by tool)
  colSums(dat2[,-1]>0) / nrow(dat2)




  # subset(dat, ABF._CS.size>0) %>%

  # dat %>% dplyr::select(ends_with("_SUM"))
  #
  # dat %>%  summarise_at(vars(ends_with("_SUM")),
  #                       funs(PROP_LOCI=sum(.,na.rm = T)/ sum(.>0) ) )
  #
}


plot_snpGroupPP_by_tool <- function(FM_all){

  dat <- FM_all %>% summarise_at(.vars = vars(ends_with(".PP")),
                                .funs = list(Overall=mean(.,na.rm=T),
                                      GWAS.nom.sig=mean(subset(.,P<.05),na.rm=T),
                                      GWAS.sig=mean(subset(.,P<5e-8),na.rm=T),
                                      GWAS.lead=mean(subset(.,leadSNP==T),na.rm=T),
                                      Credible.Set=mean(subset(.,Support>0),na.rm=T),
                                      Consensus=mean(subset(.,Consensus_SNP==T),na.rm=T)
                                 ))
  pdat <- data.frame(PP=t(dat))
  pdat$Method <- gsub("\\.PP_.*", "", rownames(pdat))
  pdat$SNP.Group <- gsub(".*\\.PP_", "", rownames(pdat))
  pdat$SNP.Group <- factor(pdat$SNP.Group, levels = unique(pdat$SNP.Group), ordered = T)

  ggplot(pdat, aes(x=SNP.Group, y=PP, fill=Method)) +
    geom_col(position = "dodge") +
    labs(y="mean PP") +
    theme_bw()
}


top_finemapped_loci <- function(dataset="./Data/GWAS/Nalls23andMe_2019",
                                save_results=T,
                                biomart=T,
                                check_if_known_locus=F,
                                force_new_subset=F){

  FM_orig <- merge_finemapping_results(minimum_support = 0,
                                       dataset = dataset,
                                       exclude_methods = NULL)
  FM_all <- FM_orig
  FM_all[is.na(FM_all)] <- 0
  # Remove manually added loci
  FM_all <- subset(FM_all, !(Gene %in%c("ATG14","SP1","LMNB1","ATP6V0A1")) )
  # Identify which loci had the smallest Union Credile Sets between SUSIE and POLYFUN+SUSIE.
  # subset(FM_all, SUSIE.CS > 0 | POLYFUN_SUSIE.CS > 0) %>%
  #   dplyr::group_by(Gene) %>% count(name = "N") %>% arrange(N)



  # List all available QTL datasets
  list_Data_dirs() %>% dplyr::filter(grepl("QTL",type)) %>% dplyr::select(Dataset, type)
  qlt.list <- c("psychENCODE_eQTL",
                "Fairfax_2014_CD14",
                "Fairfax_2014_IFN",
                "Fairfax_2014_LPS2",
                "Fairfax_2014_LPS24",
                "Cardiogenics_macrophages",
                "Cardiogenics_monocytes",
                "MESA_CAU"
                # "Brain.xQTL.Serve_eQTL",
                # "Brain.xQTL.Serve_haQTL",
                # "Brain.xQTL.Serve_mQTL"
  )
  FM_tmp <- FM_all
  for(qtl in qlt.list){
    FM_tmp <- mergeQTL.merge_handler(FM_all = FM_tmp,
                                     qtl_file = qtl,
                                     force_new_subset = force_new_subset)
    FM_tmp$QTL.sig <- ifelse(data.frame(FM_tmp)[,"QTL.P"]<=5e-8,"Y","N")
    colnames(FM_tmp) <- gsub("^QTL\\.",paste0(qtl,"."), colnames(FM_tmp))
  }
  FM_all <- FM_tmp
  qtl.sig.cols <- grep("\\.sig$",colnames(FM_all),value = T)
  FM_all$QTL.count <- rowSums(data.frame(FM_all)[,qtl.sig.cols]=="Y", na.rm = T)

  # FM_all %>% dplyr::mutate()
  # Biomart Annotations
  if(biomart){
    query_snps <- unique(subset(FM_all, Consensus_SNP)$SNP)
    printer("+ BIOMART:: Gathering annotations for",length(query_snps),"SNPs...")
    SNP.info <- biomart_snp_info(snp_list = query_snps)
    SNP.info.collapse <- SNP.info %>%
      dplyr::rename(SNP=refsnp_id) %>%
      dplyr::select(SNP, consequence_type_tv, reg_consequence_types) %>%
      dplyr::group_by(SNP) %>%
      dplyr::summarise(consequence_type_tv= paste0(unique(consequence_type_tv), collapse = "/"),
                       reg_consequence_types= paste0(unique(reg_consequence_types), collapse = "/") ) %>%
      dplyr::mutate(consequence_type_tv=gsub(", ,|NA, ", "",consequence_type_tv),
                    reg_consequence_types=gsub(", ,|NA, ", "",reg_consequence_types))
    FM_all <- data.table:::merge.data.table(FM_all,
                                            SNP.info.collapse,
                                            by = "SNP")
  }

  # Create summary data.frame
  library(tidyverse)
  FM_all$GWAS.lead <- ifelse(FM_all$leadSNP==T,"Y","N")
  grouped.dat <- FM_all %>% group_by(Gene)
  cols <- list(
    ## SNP Group counts
    grouped.dat %>%
      summarise(Consensus.RSID=paste(SNP[Consensus_SNP==T], collapse=", "),
                Consensus.ID=paste(SNP_id[Consensus_SNP==T], collapse=", "),
                CredSet.RSID=paste(SNP[Support>0], collapse=", "),
                Total.size=n(),
                GWAS.nom.sig.size=sum(P<0.05),
                GWAS.sig.size=sum(P<5e-8),
                CredSet.size=sum(Support>0),
                Consensus.size=sum(Consensus_SNP==T)) %>%
      dplyr::select(-Gene),
    ## Text cols
    ### Is the Consensus SNP the GWAS lead?
    grouped.dat %>%
      summarise_at(.vars = vars(GWAS.lead, ends_with(".sig")),
                   .funs = list(Consensus=paste(replace_na(subset(., Consensus_SNP), "N"),collapse=", "),
                                CredSet=paste(replace_na(subset(., Support>0), "N"),collapse=", ")) ),

    # Numeric cols
    ## As separated text
    grouped.dat %>%
      summarise_at(.vars = vars(ends_with(".PP"), ends_with("QTL.count"),ends_with("Effect"), MAF),
                   .funs = list(paste(round(subset(., Consensus_SNP),3),collapse=", ")) ) %>%
      dplyr::select(-Gene),
    # As means
    grouped.dat %>%
      summarise_at(.vars = vars(mean.PP, ends_with("QTL.count")),
                   .funs = list(avg=mean(subset(., Consensus_SNP), na.rm=T)) ) %>%
      dplyr::select(-Gene)
  )
  top_loci <- do.call(cbind, cols)

  if(check_if_known_locus){
    # Check whether the locus is novel according to the most recent Nalls et al (2019) PD GWAS
    ## `Known GWAS locus within 1MB` (locus-level)***
    ## `Locus within 250KB` (SNP-level?)
    Nalls <- readxl::read_excel("./Data/GWAS/Nalls23andMe_2019/Nalls2019_TableS2.xlsx")
    Nalls.novel <- (Nalls %>% dplyr::mutate(Novel.Locus=ifelse(`Known GWAS locus within 1MB`==1,"N","Y")))[c("Nearest Gene","Novel.Locus") ]%>% dplyr::rename(Gene=`Nearest Gene`) %>% unique()
    # Nalls.novel <-Nalls.novel %>%
    #   group_by(Gene) %>%
    #   summarise_each(funs(paste(., collapse = ", ")))
    top_loci <- data.table:::merge.data.table(data.table::data.table(top_loci),
                                              data.table::data.table(Nalls.novel),
                                              by="Gene") %>% dplyr::rename(Locus=Gene)
  }
  # Sort
  top_loci <- top_loci %>%
    dplyr::mutate(GWAS.lead_Consensus.any = grepl("Y",GWAS.lead_Consensus) ,
                  GWAS.lead_CredSet.any = grepl("Y",GWAS.lead_CredSet)) %>%
    arrange(Consensus.size,
            GWAS.lead_Consensus.any,
            CredSet.size,
            desc(QTL.count_avg))
  # Move loci with no consensus SNPs to the bottom
  top_loci_sort <- rbind(subset(top_loci, Consensus.size>0), subset(top_loci,Consensus.size==0))


  # Sort/filter by criterion
  # top_loci_filt <- top_loci %>%
  #   ## [0] It's one of the PD GWAS loci
  #   subset(GWAS.sig.size>0) %>%
  #   ## [1] There is at least one consensus SNP
  #   subset(Consensus.size>=1) %>%
  #   ## [2] None of the consensus SNPs are the lead GWAS SNP
  #   # dplyr::filter(grepl("N",GWAS.lead) & !grepl("Y",GWAS.lead)) %>%
  #   subset(GWAS.lead_Consensus.any==F) %>%
  #   ## [3] There's at least one QTL
  #   subset(QTL.count_avg>0) %>%
  #   ## [3] Just one consensus SNP and small Credible Set
  #   arrange(Consensus.size, CredSet.size, desc(QTL.count))

  if(save_results){
    topLoci.path <- file.path(dataset,"_genome_wide/top_loci.csv")
    printer("+ Saving top loci ==>",topLoci.path)
    data.table::fwrite(top_loci_sort, topLoci.path)
  }

  return(top_loci)
}




# multi_finemap_results_table <- function(multi_finemap_DT,
#                                         finemap_method_list,
#                                         fancy_table=F,
#                                         minimum_support=0,
#                                         include_leadSNPs=T){
#   # finemap_dat <- data.table::fread(file.path(results_path,"Multi-finemap/Multi-finemap_results.txt"),stringsAsFactors = F)
#   CS_cols <- colnames(multi_finemap_DT)[endsWith(colnames(multi_finemap_DT), ".CS")]
#   if(include_leadSNPs){
#     support_DT <- subset(multi_finemap_DT, Support >= minimum_support | leadSNP==T)
#   } else {
#     support_DT <- subset(multi_finemap_DT, Support >= minimum_support)
#   }
#   # support_DT <- subset(support_DT, select=c("Gene","SNP","CHR","POS","P","leadSNP",CS_cols,"Support"))  %>%
#   #   arrange(desc(Support))
#   support_DT <- dplyr::select(support_DT, -dplyr::one_of(c("Dataset"))) %>%
#       arrange(desc(Support))
#   # Plot table
#   if(fancy_table){
#     customGreen0 = "#DeF7E9"
#     customGreen = "#71CA97"
#     customRed = "#ff7f7f"
#     CS_formatter <-
#       formattable::formatter("span",
#                              style = x ~ style(
#                                color = ifelse(x > 0, customGreen, ifelse(x == 0, "black", "black"))))
#     formattable::formattable(support_DT,
#                              align =c("l","c","c","c","c", "c", "c", "c", "r"),
#                              list( P = formattable::color_tile(customGreen, customGreen0),
#                                    SUSIE.CS = CS_formatter,
#                                    ABF.CS = CS_formatter,
#                                    FINEMAP.CS = CS_formatter,
#                                    COJO.CS = CS_formatter,
#                                    Support = formattable::color_tile("white", "green"))
#     )
#
#   }
#   return(support_DT)
# }




leadSNP_comparison <- function(top_SNPs, merged_results){
  leadSNP_summary_table <- data.table:::merge.data.table(
    top_SNPs %>% dplyr::select(leadSNP=SNP, Gene),
    merged_results %>% dplyr::select(finemappedSNP=SNP, Gene),
    by="Gene", all=T) %>%
    dplyr::group_by(Gene, leadSNP) %>%
    dplyr::mutate(Overlap = leadSNP %in% finemappedSNP) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(Overlap = sum(Overlap)) %>% dplyr::mutate(leadSNP_in_CredSet = Overlap > 0 )

  percent_leadSNPs <- round( sum(leadSNP_summary_table$leadSNP_in_CredSet) /
                               length(leadSNP_summary_table$leadSNP_in_CredSet) * 100,2)
  return(leadSNP_summary_table[,c("Gene","leadSNP_in_CredSet")])
}




#' Tally tool-specific and union CS sizes
#'
#' @family summarise
#' @export
#' @examples
#' data("merged_DT");
#' locus_order <- SUMMARISE.get_CS_counts(merged_dat=merged_DT)
SUMMARISE.get_CS_counts <- function(merged_dat,
                                    top_CS_only=F){
  UCS_count <- suppressMessages(merged_dat %>%
    dplyr::group_by(Locus, .drop=F)  %>%
    dplyr::summarise(UCS.CS_size=dplyr::n_distinct(SNP[Support>0])))

  if(top_CS_only){
    tmp <- suppressWarnings(merged_dat %>%
                       dplyr::group_by(Locus, .drop=F) %>%
                       dplyr::summarise_at(.vars = vars(dplyr::ends_with("CS")),
                                           .funs=funs(size=dplyr::n_distinct(SNP[.==1], na.rm = T)
                                           ) ))
  }else {
    tmp <- suppressWarnings(merged_dat %>%
                              dplyr::group_by(Locus, .drop=F) %>%
                              dplyr::summarise_at(.vars = vars(dplyr::ends_with("CS")),
                                                  .funs=funs(size=dplyr::n_distinct(SNP[.>0], na.rm = T)
                                                  ) ))
  }
  locus_order <- tmp %>%
    base::merge(UCS_count, by="Locus") %>%
    dplyr::arrange(-UCS.CS_size)
  # locus_order <- locus_order[!endsWith(colnames(locus_order), ".CS_size")]
  locus_order$Locus <- factor(locus_order$Locus,  levels = locus_order$Locus, ordered = T)
  return(data.frame(locus_order))
}




#' Count bins of tool-specific and union CS sizes
#'
#' @family summarise
#' @examples
#' data("merged_DT");
#' bin_counts <- SUMMARISE.get_CS_bins(merged_dat=merged_DT)
SUMMARISE.get_CS_bins <- function(merged_dat){
  locus_order <- SUMMARISE.get_CS_counts(merged_dat = merged_dat)
  max_CS_size <- sapply(locus_order[,-1], max, na.rm=T) %>% max()
  labels = c("0","1","2-4","5-7","8-10","11-15","16+")
  bin_counts <-
    locus_order %>%
    reshape2:::melt.data.frame(measure.vars = grep("*_size$", colnames(locus_order), value = T),
                               value.name = "CS_size") %>%
    dplyr::mutate(Method=gsub("\\.CS_size$|_size$","",variable))  %>%
    dplyr::group_by(Method, .drop=F) %>%
    dplyr::mutate(bin = case_when(
      CS_size == 0 ~ labels[1],
      CS_size == 1 ~ labels[2],
      CS_size > 1 & CS_size <= 4 ~ labels[3],
      CS_size > 4 & CS_size <= 7 ~ labels[4],
      CS_size > 7 & CS_size <= 10 ~ labels[5],
      CS_size > 10 & CS_size <= 15 ~ labels[6],
      CS_size >= 16  ~ labels[7]
    ))
  bin_counts$bin <- factor(bin_counts$bin, levels = rev(labels), ordered = T)
  return(data.frame(bin_counts))
}




#' Plot CS bin counts
#'
#' @family summarise
#' @export
#' @examples
#' data("merged_DT");
#' bin_plot <- SUMMARISE.CS_bin_plot(merged_dat=merged_DT)
SUMMARISE.CS_bin_plot <- function(merged_dat,
                                  show_plot=T){
  bin_counts <- SUMMARISE.get_CS_bins(merged_dat = merged_dat)
  # Assign bin colors
  used_bins <- levels(bin_counts$bin)[levels(bin_counts$bin) %in% unique(bin_counts$bin)]
  custom_colors <- RColorBrewer::brewer.pal(n=length(levels(bin_counts$bin)), "GnBu")
  custom_colors_dict <- setNames(custom_colors[1:length(used_bins)], rev(used_bins))
  custom_colors_dict[names(custom_colors_dict)=="0"] <- "lightgray"

  bin_plot <- ggplot(subset(bin_counts, Method!="mean"), aes(x=Method, fill=bin)) +
    geom_bar(stat="count",show.legend = T, position = position_stack(reverse = F), color="white") +
    # scale_fill_brewer(palette = "Spectral", direction = -1) +
    scale_fill_manual(values = custom_colors_dict) +
    # geom_text(aes(label = paste(bin,"SNPs")), position =  position_stack(vjust = .5), vjust=-1, stat = "count") +
    geom_text(aes(label = ..count..),  position =  position_stack(vjust = .5), vjust=.5, stat = "count") +
    theme_bw() +
    labs(x=NULL, y="Loci", fill="CS size") +
    coord_flip() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          rect = element_blank(),
          axis.text.x =element_blank(),
          axis.ticks = element_blank(),
          legend.position = "top") +
    guides(fill = guide_legend(nrow = 1, reverse = T))
  if(show_plot)print(bin_plot)
  return(list(plot=bin_plot,
              data=bin_counts))
}




#' Tally locus-specific SNP group sizes
#'
#' @family summarise
#' @export
#' @examples
#' data("merged_DT");
#' snp_groups <- SUMMARISE.get_SNPgroup_counts(merged_dat=merged_DT)
SUMMARISE.get_SNPgroup_counts <- function(merged_dat,
                                          grouping_vars="Locus"){
  snp_groups <- suppressMessages(merged_dat %>%
    dplyr::group_by(.dots=grouping_vars) %>%
    dplyr::summarise(Total.SNPs=n_distinct(SNP, na.rm = T),
                     nom.sig.GWAS=n_distinct(SNP[P<.05], na.rm = T),
                     sig.GWAS=n_distinct(SNP[P<5e-8], na.rm = T),
                     CS=n_distinct(SNP[Support>0], na.rm = T),
                     Consensus=n_distinct(SNP[Consensus_SNP], na.rm = T),
                     topConsensus=n_distinct(SNP[Consensus_SNP & mean.PP==max(mean.PP)], na.rm = T ),
                     topConsensus.leadGWAS=n_distinct(SNP[Consensus_SNP & leadSNP], na.rm = T )) )
  message("Report:: all loci:")
  print( snp_groups[,!colnames(snp_groups) %in% grouping_vars] %>% colSums() / n_distinct(snp_groups$Locus))
  message("Report:: loci with at least one Consensus SNP:")
  consensus_present <- subset(snp_groups, Consensus > 0)
  print(consensus_present[,!colnames(consensus_present) %in% grouping_vars] %>% colSums() / n_distinct(consensus_present$Locus))
  return(data.frame(snp_groups))
}





#' Bar plot of tool-specific CS sizes
#'
#' Loci ordered by UCS size (smallest to largest).
#' @export
#' @family summarise
#' @examples
#' data("merged_DT")
#' gg_CS <- SUMMARISE.CS_counts_plot(merged_dat=merged_DT)
SUMMARISE.CS_counts_plot <- function(merged_dat,
                                     show_numbers=T,
                                     ylabel="Locus",
                                     legend_nrow=3,
                                     label_yaxis=T,
                                     top_CS_only=F,
                                     show_plot=T){
  locus_order <- SUMMARISE.get_CS_counts(merged_dat,
                                         top_CS_only = top_CS_only)
  melt.dat <-
    locus_order %>%
    dplyr::mutate(Locus_UCS=paste0(Locus,"  (",UCS.CS_size,")")) %>%
    reshape2:::melt.data.frame(measure.vars = grep(".CS_size$", colnames(locus_order), value = T),
                               variable.name = "CS",
                               value.name = "Credible Set size") %>%
  dplyr::mutate(Method=gsub(".CS_size$","", CS)) %>%
  dplyr::arrange(Locus, Method) %>%
  dplyr::mutate(Method=factor(Method)) %>%
  subset(Method!="mean")
  melt.dat <- order_loci(dat = melt.dat, merged_dat = merged_dat)
  melt.dat[melt.dat$`Credible Set size`==0 | is.na(melt.dat$`Credible Set size`),"Credible Set size"] <- NA


ggplot(data=melt.dat, aes(y=Locus, x=`Credible Set size`, fill=Method)) +
  geom_bar(stat = "identity", color="white", size=.05) +
  geom_text(aes(label = `Credible Set size`), color="grey20",
            size=3, show.legend = F, position = position_stack(vjust = .5)) +
  geom_text(aes(x=sum(`Credible Set size`), label = Locus_UCS),
            size=3, show.legend = F, position = position_stack(vjust = 1))

  ## Method-specific CS
  gg_CS <- ggplot(data = melt.dat,
                  aes(y=Locus, x=`Credible Set size`, fill=Method)) +
    geom_bar(stat = "identity", color="white", size=.05) +
    labs(x=NULL, y=ylabel) +
    theme_bw() +
    theme(legend.position = "top",
          # axis.text.y = element_text(),
          # axis.title.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.text = element_text(size = 8),
          legend.key.size = unit(.5, units = "cm" )) +
    guides(fill = guide_legend(nrow = legend_nrow,
                               title.position = "top",
                               title.hjust = .5))
  if(show_numbers){
    gg_CS <- gg_CS + geom_text(aes(label = `Credible Set size`), color="grey20",
              size=3, show.legend = F, position = position_stack(vjust = .5)) +
      geom_text(aes(x=sum(`Credible Set size`), label = Locus_UCS),
                size=3, show.legend = F, position = position_stack(vjust = 1))
  }

  if(label_yaxis==F){
    gg_CS <- gg_CS + theme(axis.text.y = element_blank())
  }
  if(show_plot)print(gg_CS)
  return(list(plot=gg_CS,
              data=melt.dat))
}



clean_granges <- function(gr){
  no_no_cols <- c("seqnames", "ranges", "strand", "seqlevels", "seqlengths",
                  "isCircular", "start", "end", "width", "element")
  metadat <- GenomicRanges::elementMetadata(gr)
  GenomicRanges::elementMetadata(gr) <- metadat[,!colnames(metadat) %in% no_no_cols]
  return(gr)
}


SUMMARISE.peak_overlap <- function(merged_dat,
                                   snp_filter="!is.na(SNP)",
                                   include.NOTT_2019_peaks=T,
                                   include.NOTT_2019_enhancers_promoters=T,
                                   include.NOTT_2019_PLACseq=T,
                                   include.CORCES_2020_scATACpeaks=T,
                                   include.CORCES_2020_Cicero_coaccess=T,
                                   include.CORCES_2020_bulkATACpeaks=T,
                                   include.CORCES_2020_HiChIP_FitHiChIP_coaccess=T,
                                   include.CORCES_2020_gene_annotations=T,
                                   verbose=T){
  gr.hits <- GenomicRanges::GRanges()
  ######## NOTT et al. 2019 #########
  if(include.NOTT_2019_peaks){
    try({
      NOTTpeaks <- NOTT_2019.prepare_peak_overlap(merged_dat = merged_dat,
                                                 snp_filter = snp_filter,
                                                 return_counts = F)
      NOTTpeaks <- clean_granges(NOTTpeaks)
      NOTTpeaks$Study <- "Nott et al. (2019)"
      gr.hits <- c(gr.hits, NOTTpeaks)
    })
  }

  if(include.NOTT_2019_enhancers_promoters){
    try({
      NOTTreg <- NOTT_2019.prepare_regulatory_overlap(merged_dat = merged_dat,
                                                     snp_filter = snp_filter,
                                                     return_counts = F)
      NOTTreg <- clean_granges(NOTTreg)
      NOTTreg$background <- 1
      NOTTreg$Study <- "Nott et al. (2019)"
      gr.hits <- c(gr.hits, NOTTreg)
    })
  }

  if(include.NOTT_2019_PLACseq){
    try({
      NOTTplac <- NOTT_2019.prepare_placseq_overlap(merged_dat = merged_dat,
                                                     snp_filter = snp_filter,
                                                     return_counts = F)
      NOTTplac <- clean_granges(NOTTplac)
      NOTTplac$background <- NA
      NOTTplac$Study <- "Nott et al. (2019)"
      gr.hits <- c(gr.hits, NOTTplac)
    })
  }

  ######## CORCES et al. 2020 #########
  if(include.CORCES_2020_scATACpeaks){
    try({
      CORCES_scPeaks <- CORCES_2020.prepare_scATAC_peak_overlap(merged_dat = merged_dat,
                                                                 snp_filter = snp_filter,
                                                                 add_cicero = include.CORCES_2020_Cicero_coaccess,
                                                                 annotate_genes = include.CORCES_2020_gene_annotations,
                                                                 verbose = verbose,
                                                                 return_counts = F)
      CORCES_scPeaks <- clean_granges(CORCES_scPeaks)
      CORCES_scPeaks$background <- NA
      CORCES_scPeaks$Study <- "Corces et al. (2020)"
      gr.hits <- c(gr.hits, CORCES_scPeaks)
    })
  }
  if(include.CORCES_2020_bulkATACpeaks){
    try({
      CORCES_bulkPeaks <- CORCES_2020.prepare_bulkATAC_peak_overlap(merged_dat = merged_dat,
                                                                     snp_filter = snp_filter,
                                                                     add_HiChIP_FitHiChIP = include.CORCES_2020_HiChIP_FitHiChIP_coaccess,
                                                                     annotate_genes = include.CORCES_2020_gene_annotations,
                                                                     verbose = verbose,
                                                                     return_counts = F)
      CORCES_bulkPeaks <- clean_granges(CORCES_bulkPeaks)
      CORCES_bulkPeaks$background <- NA
      CORCES_bulkPeaks$Study <- "Corces et al. (2020)"
      gr.hits <- c(gr.hits, CORCES_bulkPeaks)
    })
  }
  printer(length(gr.hits),"hits across",length(unique(gr.hits$Assay)),"assays in",
          length(unique(gr.hits$Study)),"studies found.",v=verbose)
  return(gr.hits)
}



#' Plot overlap between some SNP group and various epigenomic data
#'
#' @param include.NOTT_2019_peaks Plot SNP subset overlap with
#'  peaks from cell-type-specific bulk ATAC, H3K27ac, and H3K4me3 assays.
#' @param include.NOTT_2019_enhancers_promoters Plot SNP subset overlap with
#' cell enhancers and promoters.
#' @param include.CORCES_2020_scATACpeaks Plot SNP subset overlap with
#' cell-type-specific scATAC-seq peaks.
#' @param include.CORCES_2020_Cicero_coaccess Plot SNP subset overlap with
#' Cicero coaccessibility peaks (derived from scATACseq).
#' @export
#' @family summarise
#' @source
#' \href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}
#' \href{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}{Corces et al. (2020/bioRxiv)}
#' @examples
#' data("merged_DT");
#'
#' ... Consensus SNPs ...
#' gg_peaks <- SUMMARISE.peak_overlap_plot(merged_dat=merged_DT, snp_filter="Consensus_SNP==T", fill_title="Consensus SNPs in epigenomic peaks")
#' ... UCS SNPs ...
#' gg_peaks <- SUMMARISE.peak_overlap_plot(merged_dat=merged_DT, snp_filter="Support>0", fill_title="UCS SNPs in epigenomic peaks")
SUMMARISE.peak_overlap_plot <- function(merged_DT,
                                        snp_filter="Consensus_SNP==T",
                                        include.NOTT_2019_peaks=T,
                                        include.NOTT_2019_enhancers_promoters=T,
                                        include.NOTT_2019_PLACseq=T,
                                        include.CORCES_2020_scATACpeaks=T,
                                        include.CORCES_2020_Cicero_coaccess=T,
                                        include.CORCES_2020_bulkATACpeaks=T,
                                        include.CORCES_2020_HiChIP_FitHiChIP_coaccess=T,
                                        include.CORCES_2020_gene_annotations=T,
                                        plot_celltype_specificity=T,
                                        plot_celltype_specificity_genes=F,
                                        facets_formula=". ~ Cell_type",
                                        show_plot=T,
                                        label_yaxis=T,
                                        x_strip_angle=90,
                                        x_tick_angle=40,
                                        drop_empty_cols=F,
                                        fill_title=paste(snp_filter,"\nin epigenomic peaks"),
                                        save_path=F,
                                        height=11,
                                        width=12,
                                        subplot_widths = c(1,.5),
                                        verbose=T){
  # verbose=T;include.NOTT_2019_peaks=T; include.NOTT_2019_enhancers_promoters=T; include.NOTT_2019_PLACseq=T; include.CORCES_2020_scATACpeaks=T;
  # include.CORCES_2020_Cicero_coaccess=T;include.CORCES_2020_bulkATACpeaks=T;include.CORCES_2020_HiChIP_FitHiChIP_coaccess=T;include.CORCES_2020_gene_annotations=T;
  # no_no_loci<- c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1",
  #                "RETREG3","UBTF","FAM171A2","MAP3K14","CRHR1","MAPT-AS1","KANSL1","NSF","WNT3")
  # merged_dat <- subset(merged_dat, !Locus %in% no_no_loci)
  # x_strip_angle=90;  facets_formula=". ~ Cell_type";  drop_empty_cols=F; x_tick_angle=40;  snp_filter="Consensus_SNP==T";
  # fill_title=paste(snp_filter,"\nin epigenomic peaks"); subplot_widths = c(1,.5);plot_celltype_specificity_genes=F;

  dat_melt <- data.frame()
  ######## NOTT et al. 2019 #########
  if(include.NOTT_2019_peaks){
   try({
     dat_melt.NOTTpeaks <- NOTT_2019.prepare_peak_overlap(merged_dat = merged_dat,
                                                          snp_filter = snp_filter)
     dat_melt.NOTTpeaks$background <- NA
     dat_melt.NOTTpeaks$Study <- "Nott et al. (2019)"
     dat_melt <- base::rbind(dat_melt, dat_melt.NOTTpeaks)
   })
  }

  if(include.NOTT_2019_enhancers_promoters){
    try({
      dat_melt.NOTTreg <- NOTT_2019.prepare_regulatory_overlap(merged_dat = merged_dat,
                                                               snp_filter = snp_filter)
      dat_melt.NOTTreg$background <- 1
      dat_melt.NOTTreg$Study <- "Nott et al. (2019)"
      dat_melt <- base::rbind(dat_melt, dat_melt.NOTTreg)
    })
  }

  if(include.NOTT_2019_PLACseq){
    try({
      dat_melt.NOTTplac <- NOTT_2019.prepare_placseq_overlap(merged_dat = merged_dat,
                                                             snp_filter = snp_filter)
      dat_melt.NOTTplac$background <- NA
      dat_melt.NOTTplac$Study <- "Nott et al. (2019)"
      dat_melt <- base::rbind(dat_melt, dat_melt.NOTTplac)
    })
  }

  ######## CORCES et al. 2020 #########
  if(include.CORCES_2020_scATACpeaks){
    try({
      dat_melt.CORCES_scPeaks <- CORCES_2020.prepare_scATAC_peak_overlap(merged_dat = merged_dat,
                                                                         snp_filter = snp_filter,
                                                                         add_cicero = include.CORCES_2020_Cicero_coaccess,
                                                                         annotate_genes = include.CORCES_2020_gene_annotations,
                                                                         verbose = verbose)
      dat_melt.CORCES_scPeaks$background <- NA
      dat_melt.CORCES_scPeaks$Study <- "Corces et al. (2020)"
      dat_melt <- base::rbind(dat_melt, dat_melt.CORCES_scPeaks, fill=T)
    })
  }
  if(include.CORCES_2020_bulkATACpeaks){
    try({
      dat_melt.CORCES_bulkPeaks <- CORCES_2020.prepare_bulkATAC_peak_overlap(merged_dat = merged_dat,
                                                                             snp_filter = snp_filter,
                                                                             add_HiChIP_FitHiChIP = include.CORCES_2020_HiChIP_FitHiChIP_coaccess,
                                                                             annotate_genes = include.CORCES_2020_gene_annotations,
                                                                             verbose = verbose)
      dat_melt.CORCES_bulkPeaks$background <- NA
      dat_melt.CORCES_bulkPeaks$Study <- "Corces et al. (2020)"
      dat_melt <- base::rbind(dat_melt, dat_melt.CORCES_bulkPeaks, fill=T)
    })
  }
  ## Account for situations where include.CORCES_2020_bulkATACpeaks=F or no overlap was found
  if(!"Gene_Symbol" %in% colnames(dat_melt)) dat_melt$Gene_Symbol <- NA
  if(!"Annotation" %in% colnames(dat_melt)) dat_melt$Annotation <- NA



  plot_dat <- order_loci(dat = dat_melt,
                         merged_dat = merged_dat)
  plot_dat$Assay <- factor(plot_dat$Assay,
                           levels = c("H3K27ac","H3K4me3","ATAC","bulkATAC","scATAC","PLAC",
                                      "Cicero","HiChIP_FitHiChIP","enhancers","promoters"),
                           ordered = T)
  if(x_strip_angle!=90) plot_dat$Cell_type <- gsub(" ","\n",plot_dat$Cell_type);
  neuronal_cols <- grep("neuron",unique(plot_dat$Cell_type), value = T)
  plot_dat$Cell_type <- factor(plot_dat$Cell_type, levels = c("astrocytes","microglia","oligo","OPCs",neuronal_cols,"brain"), ordered = T)
  plot_dat$background <- as.numeric(plot_dat$background)
  ### Double check there's no errors
  plot_dat <- subset(plot_dat, !is.na(Locus))
  ### Make sure there's no missing loci
  # unique(plot_dat$Locus)==unique(merged_dat$Locus)

  # Plot
  gg_pks <- ggplot(data=plot_dat, aes(x=Assay, y=Locus, fill=Count)) +
    geom_tile(color="white") +
    # scale_fill_manual(values = consensus_colors) +
    # scale_fill_discrete(na.value = "transparent") +
    scale_fill_gradient(na.value = "transparent",
                        low = scales::alpha("blue",.7),
                        high = scales::alpha("red",.7)) +
    # geom_point(aes(size=ifelse(Count>0, "dot", "no_dot")), show.legend = F, alpha=.8, color="white") +

    # geom_rect( aes(xmin = Assay, xmax = dplyr::lead(Assay), ymin = -0.5, ymax = Inf, fill = background),
    #           alpha = 0.5, color="grey") +
    geom_tile(data = subset(plot_dat, !is.na(background)), aes( width=0.9, height=0.9), color="cyan", size=.7) +

    facet_grid(facets = formula(facets_formula),
               scales = if(drop_empty_cols) "free_x" else "fixed",
               space = "free_x") +
    scale_size_manual(values=c(dot=.5, no_dot=NA), guide="none") +
    labs(fill = fill_title) +
    theme_bw() +
    theme(legend.position = "top",
          legend.title.align = .5,
          axis.text.x = element_text(angle = x_tick_angle,
                                     hjust = if(x_tick_angle>0) 1 else NULL),
          # legend.background =  element_rect(fill = "lightgray"),
          legend.key = element_rect(colour = "gray60"),
          axis.title.y = element_blank(),
          strip.background = element_rect(fill="grey90"),
          strip.text.x = element_text(angle = x_strip_angle),
          legend.text = element_text(size = 8),
          legend.text.align = .5,
          # legend.key.size = unit(.5, units = "cm" ),
          legend.box="horizontal",
          panel.background = element_rect(fill = 'transparent'),
          # panel.grid = element_line(color="gray", size=5),
          # panel.grid.major.x = element_line(color="gray", size=.5),
          # panel.grid.major.x = element_line(color="gray", size=.5),

          panel.grid.minor = element_line(color="white", size=.5),
          plot.margin = unit(rep(.1,4), "cm")) +
    guides(color = guide_legend(nrow = 1, reverse = F,
                                title.position = "top",
                                # label.position = "top",
                                title.hjust = .5,
                                label.hjust = -1)) +
    # Keep unused levels/Loci
    scale_y_discrete(drop=FALSE)
  if(label_yaxis==F){
    gg_pks <- gg_pks + theme(axis.text.y = element_blank())
  }

  if(plot_celltype_specificity){
    try({
      library(patchwork)
      gg_cells <- SUMMARISE.cell_type_specificity(plot_dat = plot_dat,
                                                  merged_dat = merged_dat,
                                                  label_yaxis = F,
                                                  show_genes = plot_celltype_specificity_genes,
                                                  y_lab = NULL,
                                                  x_strip_angle = x_strip_angle,
                                                  show_plot = F)
      gg_pks <- gg_pks + gg_cells$plot + patchwork::plot_layout(nrow = 1, widths = subplot_widths)
    })
  }

  if(show_plot) print(gg_pks)
  if(save_path!=F){
    printer("+ Saving plot ==>",save_path,v=verbose)
    ggplot2::ggsave(save_path, gg_pks, height=height, width=width)
  }
  return(list(data=dat_melt,
              plot=gg_pks))
}






#' Get cell-type-specifity score for each cell type
#'
#' Aggregate SNP overlap across various epigenomic datasets
#' and then identify the number of SNPs overlapping by each cell type
#'
SUMMARISE.cell_type_specificity <- function(plot_dat,
                                            merged_dat,
                                            min_count=NULL,
                                            top_celltype_only=F,
                                            label_yaxis=T,
                                            y_lab=NULL,
                                            show_genes=F,
                                            x_strip_angle=40,
                                            show_plot=T){
  Cell_group_dict <- c("astrocytes"="astrocytes",
                       "microglia"="microglia",
                       "oligo"="oligo",
                       "OPCs"="oligo",
                       "neurons"="neurons",
                       "neurons (+)"="neurons",
                       "neurons (-)"="neurons",
                       "neurons (nigral)"="neurons",
                       "brain"="brain")

  cell_tally <- plot_dat %>%
    dplyr::mutate(Assay_count=ifelse(Count>0,1,0), # Set any overlap ==1
                  Cell_group=factor(Cell_group_dict[Cell_type],
                                    levels = unique(unname(Cell_group_dict)), ordered = T )) %>%
    dplyr::group_by(Locus, Cell_group, .drop=F) %>%
    dplyr::summarise(Assay_count=sum(Assay_count,na.rm=T),
                     SNP_Count=sum(Count,na.rm = T),
                     Gene_Symbol=gsub("^NA$",NA,Gene_Symbol),
                     Annotation=Annotation )

  if(top_celltype_only){
    cell_tally <- dplyr::top_n(cell_tally, n = 1, wt="Cell_group")
  }
  if(!is.null(min_count)){
    cell_tally[cell_tally$Count < min_count,"Count"] <- NA
  }

  cell_tally <- order_loci(dat = cell_tally,
                           merged_dat = merged_dat)
  gg_tally <- ggplot(data = cell_tally, aes(x=Cell_group, y=Locus, fill=Assay_count)) +
    geom_tile(color="white") +
    facet_grid(facets = . ~ Cell_group,
               scales = "free_x") +
    scale_fill_viridis_c(na.value = "transparent") +
    labs(y=y_lab) +
    theme_bw() +
    # scale_x_discrete(position = "top") +
    theme(axis.text.x = element_blank(),#element_text(angle = x_text_angle, hjust = 0),
          legend.box="horizontal",
          legend.position = "top",
          legend.text = element_text(size = 8),
          legend.text.align = .5,
          strip.text.x = element_text(angle=x_strip_angle, color="white"),
          strip.background.x = element_rect(fill="black"),
          panel.spacing = unit(.1, "lines"),
          plot.margin = unit(c(.1,2,.1,.1), "cm")) +
    scale_y_discrete(drop=FALSE)
  if(show_genes){
    gg_tally <- gg_tally +
      geom_text(aes(label=eval(parse(text="Gene_Symbol"))), size=3, color="cyan")
  }
  if(label_yaxis==F){
    gg_tally <- gg_tally + theme(axis.text.y = element_blank())
  }
  if(label_yaxis=="right"){
    gg_tally <- gg_tally + scale_y_discrete(position = "right")
  }
  if(show_plot) print(gg_tally)
  return(list(plot=gg_tally,
              data=cell_tally))
}


#' Nominate target genes within each locus
#'
#' Across all GWAS-QTL colocalization tests across all studies,
#' take the eGene with the highest colocalziation probability (PP.H4)
#' and assign it as the most likely causal gene in that locus.
#'
#' eQTL queries and colocalization test done with \pkg{catalogueR}.
#' @export
#' @examples
#' \dontrun{
#' data("merged_DT")
#' base_url <- "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"
#' coloc_results_path <- file.path(base_url,"_genome_wide/COLOC/coloc.eQTL_Catalogue_ALL.csv.gz")
#' gg_egene <- SUMMARISE.coloc_nominated_eGenes(coloc_results, merged_dat=merged_DT, fill_var=NULL)
#'
#' # QTL
#' base_url <- "/sc/hydra/projects/ad-omics/microglia_omics/Fine_Mapping"
#' coloc_results_path <- file.path(base_url, "Kunkle_Microglia_all_regions/QTL_merged_coloc_results.snp.tsv.gz")
#' merged_dat <- data.table::fread("/pd-omics/brian/Fine_Mapping/Data/QTL/Microglia_all_regions/multiGWAS.microgliaQTL_finemapping.csv.gz")
#' gg_egene <- SUMMARISE.coloc_nominated_eGenes(coloc_results, merged_dat=merged_dat, fill_var=NULL)
#' }
SUMMARISE.coloc_nominated_eGenes <- function(coloc_results,
                                             merged_dat,
                                             label_yaxis=T,
                                             y_lab="Locus",
                                             x_lab=NULL,
                                             fill_var='PP.H4',
                                             text_size=2,
                                             PP_threshold=NULL,
                                             nThread=4,
                                             show_plot=T,
                                             verbose=T){
  # Check Corces gene annotations against eQTL/coloc eGenes
  printer("+ SUMMARISE:: Nominating genes by top colocalized eQTL eGenes",v=verbose)
  if(is.data.frame(coloc_results)){
    dat <- coloc_results
  } else {
    dat <- data.table::fread(coloc_results, nThread = nThread)
  }

  # for(column %in% c("gene","snp","chr"))
  # if('gene' %in%)

  top_eGenes <- dat %>%
    subset(PP.H4>if(is.null(PP_threshold)) 0 else PP_threshold) %>%
    # Remove RP11 and other spurious RP genes
    subset(!(startsWith(eGene,"RP") | eGene=="NA" | is.na(eGene)) ) %>%
    dplyr::group_by(Locus.GWAS) %>%
    dplyr::top_n(n = 1, wt = PP.H4) %>%
    # Ensure only 1 eGene per Locus
    dplyr::arrange(desc(PP.H4),desc(SNP.PP.H4)) %>%
    dplyr::group_by(Locus.GWAS) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(Locus=as.character(Locus.GWAS))

  top_eGenes <- order_loci(dat = top_eGenes,
                           merged_dat = merged_dat)
  top_eGenes <- subset(top_eGenes, !is.na(Locus))
  top_eGenes$dummy <- "Top\ncolocalized\neGene"

  if(is.null(fill_var)){
    text_color="grey20"
  }else {text_color="white"}

  gg_egene <- ggplot(top_eGenes, aes(x=dummy, y=Locus)) +
    labs(x=x_lab, y=y_lab) +
    geom_tile(fill='transparent') +
    geom_text(aes(label=eGene), color=text_color, size=text_size) +
    # scale_fill_viridis_c(end = .8, na.value = "transparent") +
    # scale_fill_gradient(low = "blue", high = "red", na.value = "transparent") +
    theme_bw() +
    # scale_x_discrete(position = "top") +
    theme(#axis.text.x = element_blank(),
          legend.box="vertical",
          legend.position = "top",
          legend.text = element_text(size = 8),
          legend.text.align = .5,
          plot.margin = unit(rep(.1,4), "cm")) +
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
           size = guide_legend(title.position="top", title.hjust = 0.5)) +
    scale_y_discrete(drop=F)
  if(label_yaxis==F){
    gg_egene <- gg_egene + theme(axis.text.y = element_blank())
  }
  if(show_plot) print(gg_egene)
  return(list(data=data.table::data.table(top_eGenes),
              plot=gg_egene))
}





count_and_melt <- function(merged_annot,
                           snp_filter="Consensus_SNP==T",
                           grouping_vars=c("Locus","Cell_type","Assay")){
  consensus_melt <-
    data.table::setDT(merged_annot)[, .(Count = dplyr::n_distinct(SNP[eval(parse(text=snp_filter))], na.rm = T)),
                                    by=grouping_vars]
  if(length(grouping_vars)>=3){
    consensus_melt <- subset(consensus_melt,
                             !is.na(dplyr::vars(grouping_vars[2])) &
                               !is.na(dplyr::vars(grouping_vars[3]) ), .drop=F)%>%
      dplyr::mutate(Celltype_Assay = paste0(eval(parse(text=grouping_vars[2])),"_",eval(parse(text=grouping_vars[3]))))
  }
  return(consensus_melt)
}







#' Give a quick summary report of the fine-mapping results
#' @family summarise
#' @examples
#' data("merged_DT");
#' results_report(merged_DT)
results_report <- function(merged_dat){
  message("echolocatoR results report (all loci):")
  printer("+ Overall report:")
  printer("++",length(unique(merged_dat$Locus)),"Loci.")
  printer("++",nrow(merged_dat),"SNPs.")
  cat('\n')
  printer("+ Lead SNP report:")
  printer("++", length(subset(merged_dat, leadSNP)$SNP),"lead SNPs.")
  printer("++ Lead SNP mean PP =", round(mean(subset(merged_dat, leadSNP)$mean.PP, na.rm = T),2))
  cat('\n')
  printer("+ Union Credible Set report:")
  printer("++", length(subset(merged_dat, Support>0)$SNP),"UCS SNPs.")
  printer("++ UCS mean PP =", round(mean(subset(merged_dat, Support>0)$mean.PP, na.rm = T),2))
  printer("++",nrow(subset(merged_dat, Support>0 & leadSNP)),"UCS SNPs that are also lead SNPs")
  cat('\n')
  printer("+ Consensus SNP report:")
  printer("++", length(subset(merged_dat, Consensus_SNP)$SNP),"Consensus SNPs.")
  printer("++ Consensus SNP mean PP =", round(mean(subset(merged_dat, Consensus_SNP)$mean.PP, na.rm = T),2))
  printer("++",nrow(subset(merged_dat, Consensus_SNP & leadSNP)),"Consensus SNPs that are also lead SNPs")
}




#' Merge all summary plots into one super plot
#'
#' @family summarise
#' @export
super_summary_plot <- function(merged_dat,
                               snp_filter="Consensus_SNP==T",
                               coloc_results=NULL,
                               plot_missense=T,
                               show_plot=T,
                               save_plot=F,
                               height=15,
                               width=13,
                               dpi=500){
  library(patchwork)
  bin_plot <- SUMMARISE.CS_bin_plot(merged_dat = merged_dat,
                                    show_plot = F)

  if(!is.null(coloc_results)){
    gg_egene <- SUMMARISE.coloc_nominated_eGenes(coloc_results = coloc_results,
                                                 merged_dat=merged_dat,
                                                 PP_threshold = .8,
                                                 fill_var = NULL,
                                                 text_size = 2.5,
                                                 y_lab = "Locus",
                                                 x_lab = NULL,
                                                 label_yaxis = T,
                                                 show_plot = F)
    gg_egene_width <- .125
  } else {
    gg_egene<-c();
    gg_egene$plot <- patchwork::plot_spacer();
    gg_egene_width <- 0
  }

  gg_CS <- SUMMARISE.CS_counts_plot(merged_dat = merged_dat,
                                    show_numbers=F,
                                    label_yaxis=F,
                                    ylabel=NULL,
                                    show_plot = F)
  #### plot_missense ####
  if(plot_missense){
    try({
      gg_missense <- ANNOTATE.plot_missense(merged_dat = merged_dat,
                                            snp_filter=snp_filter,
                                            show.legend = F,
                                            show_plot = F)
    })
  }
  # In case biomart times out
  if(!exists("gg_missense")){
    gg_missense<-c();
    gg_missense$plot <- patchwork::plot_spacer();
    gg_missense_width <- 0
  } else {gg_missense_width <- .01}

  gg_peaks <- SUMMARISE.peak_overlap_plot(merged_dat = merged_dat,
                                          snp_filter=snp_filter,
                                          include.NOTT_2019_peaks=T,
                                          include.NOTT_2019_enhancers_promoters=T,
                                          include.NOTT_2019_PLACseq=T,
                                          include.CORCES_2020_scATACpeaks=T,
                                          include.CORCES_2020_Cicero_coaccess=F,
                                          include.CORCES_2020_bulkATACpeaks=T,
                                          include.CORCES_2020_HiChIP_FitHiChIP_coaccess=T,
                                          include.CORCES_2020_gene_annotations=T,
                                          plot_celltype_specificity=T,
                                          facets_formula=". ~ Cell_type",
                                          show_plot=F,
                                          label_yaxis=T,
                                          subplot_widths = c(1,.2),
                                          x_strip_angle = 90,
                                          drop_empty_cols = T,
                                          fill_title=paste(snp_filter,"SNPs\nin epigenomic peaks"),
                                          # save_path="~/Desktop/super_peak_plot.png",
                                          verbose=T)
  # Merge
  gg_merged <- (patchwork::plot_spacer() + bin_plot$plot + patchwork::plot_layout(widths = c(.4,.6)))  /
    (gg_egene$plot + gg_CS$plot + gg_missense$plot + gg_peaks$plot +  patchwork::plot_layout(widths = c(gg_egene_width,.3,gg_missense_width,1))) +
    patchwork::plot_layout(heights = c(.15,1),ncol = 1)

  if(show_plot) print(gg_merged)
  if(save_plot!=F){
    ggplot2::ggsave(save_plot,
                    gg_merged,
                    dpi = dpi,
                    height=height, width=width)
  }
 return(list(data=gg_peaks$data,
             plot=gg_merged))
}








#' Plot inter-study SNP overlap
#'
#' Cross-tabulate SNP overlap (after applying filter)
#' between each pair of studies.
#' @family summarise
SUMMARISE.plot_dataset_overlap <- function(merged_dat,
                                           snp_filter="!is.na(SNP)",
                                           filename=NA,
                                           formula_str="~ SNP + Dataset",
                                           triangle=F,
                                           proxies=NULL){
  snp_xtab <- subset(merged_dat,  eval(parse(text = snp_filter)), .drop=F) %>%
    stats::xtabs(formula = stats::as.formula(formula_str),
                 sparse = F,
                 drop.unused.levels = F)
  snp_xprod <- crossprod(snp_xtab)
  diag(snp_xprod) <- NA
  mode(snp_xprod) <- "integer"

  if(triangle){
    max_count <- max(snp_xprod, na.rm = T)
    printer("max_count =",max_count)
    cl.length <- if(max_count<=10) max_count else median(DescTools::Divisors(max_count)[[1]]) + 1
    printer("cl.length =",cl.length)
    png(filename,height=500,width = 500,type = "cairo")
    dat <- corrplot::corrplot(corr = snp_xprod,
                              method = "color",
                              type = "lower",
                              addgrid.col = "grey",
                              tl.col = "black",
                              hclust.method = "ward.D2",
                              title = paste("SNP overlap:",gsub("[|]","\nOR",snp_filter)),
                              order = "hclust",
                              cl.length = cl.length,
                              mar = c(0,0,4,4),
                              # tl.pos = "lt",
                              diag = F,
                              is.corr = F)
    dev.off();
  } else {
    pheatmap::pheatmap(snp_xprod,
                       display_numbers=T,
                       filename = filename,
                       # number_color = "white",
                       main = paste("SNP overlap:",snp_filter),
                       angle_col = 45,
                       cluster_cols = T,
                       cluster_rows = T,
                       drop_levels = F,
                       na_col = "white")
  }
  return(snp_xprod)
}



