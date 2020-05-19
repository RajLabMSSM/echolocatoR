
### Summarise fine-mapping results from all loci ###
# merged_DT <- merge_finemapping_results(consensus_thresh = 3)
# (merged_DT %>% subset(P<5e-8) %>%dplyr::group_by(Gene) %>% count())$n %>% mean()


lead.SNP.coords <- function(){
  annot <- readxl::read_excel("./Data/annotated_finemapping_results.xlsx")
  annot <- find_consensus_SNPs(annot, support_thresh = 2)
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
    subset(SUSIE.Credible_Set>0 | FINEMAP.Credible_Set > 0| PAINTOR.Credible_Set >0)
  statCS.means <- SNPgroup.summary(statCS, group.name="Statistical.CS")
  # Functional FM
  funcCS <- annot %>% group_by(Gene) %>% subset(PAINTOR.Credible_Set>0)
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
  dat <- FM_all %>%
    summarise_at(vars(ends_with(".Credible_Set")),
                 funs(CS.sum=sum(.>0, na.rm = T),
                      leadGWAS.sum=sum(leadSNP==T, na.rm = T),
                      leadGWAS.CS.sum=sum(leadSNP==T & .>0, na.rm = T),
                      PROP.SNPS=sum(leadSNP==T & .>0, na.rm = T)/sum(.>0, na.rm = T)) )
  colnames(dat) <- gsub("\\.Credible_Set","", colnames(dat))
  dat

  # Proportion of loci in which the CS contains the leadSNP (by tool)
  # dat3.1 <- FM_all %>% group_by(Gene) %>%
  #   summarise_at(vars(ends_with(".Credible_Set")),
  #                funs(leadGWAS.CS.size=sum(leadSNP==T & .>0, na.rm = T)))
  # dat3.2 <- FM_all %>% group_by(Gene) %>%
  #   summarise_at(vars(ends_with(".Credible_Set")),
  #                funs(leadGWAS.CS.size=sum(.>0, na.rm = T)))
  # colSums(dat3.1[,-1]) / colSums(dat3.2[,-1])

  # Number of CS SNPs (by tool)
  dat2 <- FM_all %>% group_by(Gene) %>% summarise_at(vars(ends_with(".Credible_Set")),
                                                     funs(CS=sum(.,na.rm=T)))
  colnames(dat2) <- gsub("\\.Credible_Set","", colnames(dat2))
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

  dat <- FM_all %>% summarise_at(vars(ends_with(".PP")),
                                 funs(Overall=mean(.,na.rm=T),
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
  # subset(FM_all, SUSIE.Credible_Set > 0 | POLYFUN_SUSIE.Credible_Set > 0) %>%
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
      summarise_at(vars(GWAS.lead, ends_with(".sig")),
                   funs(Consensus=paste(replace_na(subset(., Consensus_SNP), "N"),collapse=", "),
                        CredSet=paste(replace_na(subset(., Support>0), "N"),collapse=", ")) ),

    # Numeric cols
    ## As separated text
    grouped.dat %>%
      summarise_at(vars(ends_with(".PP"), ends_with("QTL.count"),ends_with("Effect"), MAF),
                   funs(paste(round(subset(., Consensus_SNP),3),collapse=", ")) ) %>%
      dplyr::select(-Gene),
    # As means
    grouped.dat %>%
      summarise_at(vars(mean.PP, ends_with("QTL.count")),
                   funs(avg=mean(subset(., Consensus_SNP), na.rm=T)) ) %>%
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
#   # finemap_DT <- data.table::fread(file.path(results_path,"Multi-finemap/Multi-finemap_results.txt"),stringsAsFactors = F)
#   CS_cols <- colnames(multi_finemap_DT)[endsWith(colnames(multi_finemap_DT), ".Credible_Set")]
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
#                                    SUSIE.Credible_Set = CS_formatter,
#                                    ABF.Credible_Set = CS_formatter,
#                                    FINEMAP.Credible_Set = CS_formatter,
#                                    COJO.Credible_Set = CS_formatter,
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


