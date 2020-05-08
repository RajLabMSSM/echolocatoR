
##################################################################
#################       psychENCODE Data         #################
##################################################################

# Data Repository:
# http://resource.psychencode.org/?_ga=2.202944787.381602893.1563834178-1901423550.1563307711

# Reference:
# https://science.sciencemag.org/content/362/6420/eaat8464

# library(dplyr)
# library(ggplot2)

psychENCODE.download_summary_stats <- function(
  output_dir = "./Data/QTL/psychENCODE",
  url.list = file.path("http://resource.psychencode.org/Datasets/",
                       c("Derived/QTLs/DER-08a_hg19_eQTL.significant.txt",
                         "Derived/QTLs/DER-09_hg19_cQTL.significant.txt",
                         "Derived/QTLs/DER-10a_hg19_isoQTL.significant.txt",
                         "Integrative/INT-16_HiC_EP_linkages_cross_assembly.csv"))
                                               ){
  printer("psychENCODE:: Downloading summary stats.")
  # QTLs
  for(URL in url.list){
    output.path <- file.path(output_dir, basename(URL))
    if(file.exists(output.path) | file.exists(paste0(output.path,".gz")) ){
      printer("psychENCODE:: File already exists.")
    } else{
      printer("psychENCODE:: Downloading:")
      printer("              ",basename(URL))
      dat <- data.table::fread(URL)
      data.table::fwrite(dat,  output.path)
      R.utils::gzip(output.path)
    }
  }
}

summarise_SNPgroup_overlap <- function(FM_sub, ASSAY_sub, assay_type, gene){
  # % CS SNPs that are sig HITS
  CS <- subset(FM_sub, Support > 0)$SNP_id %>% unique()
  # % Consensus SNPs that are sig HITS
  Consensus <- subset(FM_sub, Consensus_SNP==T)$SNP_id %>% unique()
  suppressWarnings(
    # % Lead SNPs that are sig HITS
    if(assay_type=="HiC.ENH"){
      CS <- gsub(":*", "",CS)
      Consensus <- gsub(":*", "",Consensus)
      top.hit <- NA
      GWAS_lead <- subset(FM_sub, leadSNP==T)$POS %>% unique()
      GWAS_lead.HIT <- NA
      HIT.total <- dim(ASSAY_sub)[1]
      CS.HIT <- sum(CS >= min(ASSAY_sub$Enhancer_Start_hg19) &
                      CS <= max(ASSAY_sub$Enhancer_End_hg19) )
      Consensus.HIT <- sum(Consensus >= min(ASSAY_sub$Enhancer_Start_hg19) &
                             Consensus <= max(ASSAY_sub$Enhancer_End_hg19) )
    } else{
      top.hit <- ASSAY_sub[1,]$SNP_id
      GWAS_lead <- subset(FM_sub, leadSNP==T)$SNP_id %>% unique()
      GWAS_lead.HIT <- sum(GWAS_lead %in% unique(ASSAY_sub["SNP_id"]))
      HIT.total <- length(unique(ASSAY_sub$SNP_id))
      CS.HIT = sum(CS %in% unique(ASSAY_sub$SNP_id))
      Consensus.HIT <- sum(Consensus %in% unique(ASSAY_sub$SNP_id))
    }
  )

  # Create report
  dat <- data.frame(HIT.total = HIT.total,
                    GWAS_lead.total = length(GWAS_lead),
                    GWAS_lead.HIT = GWAS_lead.HIT,
                    GWAS_lead.has_top_HIT = top.hit %in% GWAS_lead,
                    CS.total = length(CS),
                    CS.HIT = CS.HIT,
                    CS.has_top_HIT = top.hit %in% CS,
                    Consensus.total = length(Consensus),
                    Consensus.HIT = Consensus.HIT,
                    Consensus.has_top_HIT = top.hit %in% Consensus
  ) %>% dplyr::mutate(GWAS_lead.HIT_percent = round(GWAS_lead.HIT /GWAS_lead.total*100, 1),
                      CS.HIT_percent = round(CS.HIT /CS.total*100, 2),
                      Consensus.HIT_percent = round(Consensus.HIT /Consensus.total*100, 2)
  )
  dat[is.na(dat)] <- 0
  dat$Feature <- gene
  dat$ASSAY_type <- assay_type
  return(dat)
}

psychENCODE.assay_summary <- function(){
  root <- "./Data/QTL/psychENCODE"
  ASSAY_files <- file.path(root,
                         c("DER-08a_hg19_eQTL.significant.txt.gz",
                           "DER-09_hg19_cQTL.significant.txt.gz",
                           "DER-10a_hg19_isoQTL.significant.txt.gz",
                           "INT-16_HiC_EP_linkages_cross_assembly.csv.gz"))
  ASSAY_files = setNames(ASSAY_files, nm = c("eQTL","cQTL","isoQTL","HiC"))

  # Import annotated fine-mapping results across all loci
  # FM_results <- data.table::as.data.table(readxl::read_excel("Data/annotated_results_table.xlsx"))
  FM_results <- merge_finemapping_results(minimum_support = 0)
  FM_results$SNP_id <- paste0(FM_results$CHR,":",FM_results$POS)
  # assay_type = names(ASSAY_files)[1]
  percent_df <- lapply(names(ASSAY_files), function(assay_type){
    print(assay_type)
    # Import ASSAY data
    print(paste("+ Importing",assay_type,"data..."))
    ASSAY <- data.table::fread(ASSAY_files[assay_type], check.names = T)

    # Iterate over each locus
    feature_df <- lapply(unique(FM_results$Gene), function(gene, ASSAY.=ASSAY){
      print(gene)
      FM_sub <- subset(FM_results, Gene==gene)
      # HiC
      if(assay_type=="HiC"){
        ENH_sub <- subset(ASSAY., Enhancer_Chromosome_hg19 == paste0("chr",unique(FM_sub$CHR)) &
                            Enhancer_Start_hg19>=min(FM_sub$POS) & Enhancer_End_hg19<= max(FM_sub$POS)) %>%
          mutate(SNP_chr=Enhancer_Chromosome_hg19,
                 SNP_id = paste0("chr",Enhancer_Chromosome_hg19,":",Transcription_Start_Site_hg19))
        TSS_sub <- subset(ASSAY., Enhancer_Chromosome_hg19 == paste0("chr",unique(FM_sub$CHR)) &
                            Transcription_Start_Site_hg19 %in% FM_sub$POS) %>%
          mutate(SNP_chr=Enhancer_Chromosome_hg19,
                 SNP_id = paste0("chr",Enhancer_Chromosome_hg19,":",Transcription_Start_Site_hg19))

        dat.ENH <- summarise_SNPgroup_overlap(FM_sub, ENH_sub, assay_type="HiC.ENH", gene=gene)
        dat.TSS <- summarise_SNPgroup_overlap(FM_sub, TSS_sub, assay_type="HiC.TSS", gene=gene)
        dat <- rbind(dat.ENH, dat.TSS)
      }
      # QTL
      if(assay_type %like% "QTL"){
        ## Subset by gene to capture whole QTL region
        # if(n=="eQTL"){
        #   ASSAY_sub <- subset(QTL, ID == g)
        # }

        # Subset by max-max locus coordinates to make sure we're capturing
        ## all SNPs with a given region tested by the QTL
        ASSAY_sub <- subset(ASSAY., SNP_chr == paste0("chr",unique(FM_sub$CHR)) &
                            SNP_start >= min(FM_sub$POS) & SNP_end <= max(FM_sub$POS))
        # Subset by FDR
        ASSAY_sub <- subset(ASSAY_sub, FDR<=0.05) %>% arrange(FDR)
        # ASSAY_sub <- subset(QTL, SNP_id %in% FM_sub$SNP_id)
        dat <- summarise_SNPgroup_overlap(FM_sub, ASSAY_sub, assay_type=assay_type, gene=gene)
      }
      return(data.table::as.data.table(dat))
    }) %>% data.table::rbindlist()
    # # Merge with FM data
    # merged_dat <- QTL_merge(FM_results[,c("SNP","CHR","POS","leadSNP","Consensus_SNP","SNP_id")],
    #                         QTL_path = ASSAY_files[n])
    # percent_summ <- percent.HIT.summary(FM_results, merged_dat, QTL_type = n)
    # return(data.table::as.data.table(percent_summ))
    return(feature_df)
  }) %>% data.table::rbindlist()
}


get_group_sums <- function(per_df){
  group_sums <- per_df[,c("ASSAY_type","GWAS_lead.HIT","Consensus.HIT","CS.HIT")] %>%
    dplyr::group_by(ASSAY_type) %>% summarise_each(funs(sum(., na.rm=T) )) %>%
    data.table::melt(id.vars="ASSAY_type",variable.name="Variable", value.name="Total")
  group_sums <- tidyr::separate(group_sums, "Variable", c("SNP.Group",NA), "\\.")
  # group_sums$Total <- paste0("n=",group_sums$Total)
  return(group_sums)
}

plot_percent_HITs <- function(per_df){
  # Consider the fact that only some loci have Consensus SNP (count those as NA, not zeros)
  Consensus.cols <- grep("Consensus",colnames(per_df), value = T)
  per_df[per_df$Consensus.total==0, Consensus.cols] <- NA

  percent_cols <- grep("percent", colnames(per_df), value = T)
  mean_df <- data.table::melt(per_df,
                              id.vars = c("Feature","ASSAY_type"),
                              measure.vars = percent_cols,
                              variable.name = c("Variable"),
                              value.name = c("Percent"))
  mean_df$SNP.Group <- gsub("\\..*","",mean_df$Variable)
  mean_DF <- mean_df %>% dplyr::group_by(SNP.Group, ASSAY_type) %>%
    dplyr::summarise(Percent = mean(Percent, na.rm = T))
  # Get sums
  group_sums <- get_group_sums(per_df)
  # Merge mean percent and total sum DTs
  merged_DF <- data.table:::merge.data.table(data.table::data.table(mean_DF),
                                data.table::data.table(group_sums),
                                by = c("ASSAY_type","SNP.Group"))
  # Plot
  g <- ggplot(data=merged_DF, aes(x=ASSAY_type, y=Percent, fill=SNP.Group ), color=SNP.Group ) +
    geom_col(show.legend = T, alpha=.9, position = "dodge") +
    labs(title=paste0("psychENCODE Frontal Cortex"),
         subtitle = paste0("% SNP group containing any significant hit"),
         x="Assay Type", y="% SNPs") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.text = element_text(color ="gray14", size=8),
          legend.title = element_text(color = "gray14", size=10)) +
    geom_text(aes(label=Total), vjust=-0.3, size=3.5, position = position_dodge(1)) +
    scale_fill_brewer(palette="Paired")
  # print(g)
  return(g)
}


plot_has_top_HIT <- function(per_df){
  # per_df <- percent_df
  per_df <- subset(per_df, !(ASSAY_type %like% 'HiC'))
  top.hit_cols <- grep(".has_top_HIT", colnames(per_df), value = T)
  # Consider the fact that only some loci have Consensus SNP (count those as NA, not zeros)
  Consensus.cols <- grep("Consensus",colnames(per_df), value = T)
  per_df[per_df$Consensus.total==0, Consensus.cols] <- NA

  mean_df <- data.table::melt(per_df,
                              id.vars = c("Feature","ASSAY_type"),
                              measure.vars = top.hit_cols,
                              variable.name = c("Variable"),
                              value.name = c("Bool"))
  mean_df$SNP.Group <- gsub("\\..*","",mean_df$Variable)

  bool_df <- mean_df %>% dplyr::group_by(SNP.Group, ASSAY_type, Feature) %>%
    dplyr::summarise(Max_Bool = max(Bool)) %>%
    dplyr::group_by(SNP.Group, ASSAY_type) %>%
    dplyr::summarise(Percent = mean(Max_Bool, na.rm=T)*100)
  # Get sums
  group_sums <- get_group_sums(per_df)
  # Merge mean percent and total sum DTs
  merged_DF <- data.table:::merge.data.table(data.table::data.table(bool_df),
                                             data.table::data.table(group_sums),
                                             by = c("ASSAY_type","SNP.Group"))
  # Plot
  g2 <- ggplot(data=merged_DF, aes(x=ASSAY_type, y=Percent, fill=SNP.Group ), color=SNP.Group ) +
    geom_col(show.legend = T, alpha=.9, position = "dodge") +
    labs(title=paste0("psychENCODE Frontal Cortex"),
         subtitle = paste0("% SNP group containing the top hit"),
         x="Assay Type", y="% SNPs") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.text = element_text(color ="gray14", size=8),
          legend.title = element_text(color = "gray14", size=10)) +
    geom_text(aes(label=Total), vjust=-0.3, size=3.5, position = position_dodge(1)) +
    scale_fill_brewer(palette="Paired")
  return(g2)
}



psychENCODE.overlap_plots <- function(percent_df){
  # Plot
  p1 <- plot_percent_HITs(percent_df) + ylim(0,50)
  p2 <- plot_has_top_HIT(percent_df) + ylim(0,50)
  cp <- cowplot::plot_grid(p1,p2,labels = LETTERS)
  print(cp)
}

 

psychENCODE.ENS_to_HGNC <- function(QTL, reference_genome="grch37"){
  print("+ Converting Ensembl IDs to HGNC Gene Symbols")
  # Identify and rename feature column
  feature_cols <- c("gene_id","Peak_id","transcript_id")
  feat_col <- feature_cols[feature_cols %in% colnames(QTL)]
  QTL <- QTL %>% dplyr::rename("feature_ID" = feat_col)
  # Remove info after dot
  QTL$feature_ID <- gsub("\\..*","",QTL$feature_ID)
  # Specify host (grch38 is the default)
  if(reference_genome=="grch38"){
    host <- "www.ensembl.org"
  } else {host <- paste0(reference_genome,".ensembl.org")}
  # BIOMART
  mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                          host=host,
                          dataset = "hsapiens_gene_ensembl")
  # View(listFilters(mart))
  # View(listAttributes(mart))
  res <- biomaRt::getBM(filters= "ensembl_gene_id",
                        attributes= c("ensembl_gene_id","ensembl_transcript_id","hgnc_symbol"),
                        values = unique(QTL$feature_ID),
                        mart = mart)
  # Translate from feature to gene (or transcript)
  gene_dict <- setNames(res$hgnc_symbol, nm = res$ensembl_gene_id)
  QTL$ID <- gene_dict[QTL$feature_ID]
  return(QTL)
}


 

#
# SNPs_that_are_QTLs <- function(QTL_type="eQTL"){
#   root <- "./echolocatoR/tools/Annotations/psychENCODE"
#   ASSAY_files <- file.path(root,
#                          c("DER-08a_hg19_eQTL.significant.txt",
#                            "DER-09_hg19_cQTL.significant.txt",
#                            "DER-10a_hg19_isoQTL.significant.txt"))
#   ASSAY_files = setNames(ASSAY_files, nm = c("eQTL","cQTL","isoQTL"))
#
#
#   # Gather all merged data
#   merged.dat <- merge_finemapping_results(minimum_support = 1)
#   merged.dat["SNP_id"] <- paste0(merged.dat$CHR,":",merged.dat$POS)
#
#   # QTL_type <- names(ASSAY_files)[1]
#   QTL.results <- lapply(names(ASSAY_files)[1], function(QTL_type){
#     printer("Summarizing fine-mapping vs.",QTL_type,"results")
#     # Import QTL data
#     eQTL <- data.table::fread(ASSAY_files[QTL_type])
#     merge.HIT.all <- data.table:::merge.data.table(eQTL,
#                                                data.table::as.data.table(merged.dat),
#                                                by = "SNP_id",
#                                                all.y = T)
#     all.fractions <- lapply(unique(merged.dat$Gene), function(gene){
#       printer("+ Gene",gene)
#       merge.HIT <- subset(merge.HIT.all, Gene==gene)
#
#
#       N.leadSNP = subset(merge.HIT, leadSNP)$SNP %>% unique() %>% length()
#       N.leadSNPs.x.HIT = subset(merge.HIT, leadSNP==T & FDR <= 0.05)$SNP %>% unique() %>% length()
#
#       N.CS = subset(merge.HIT, Support>0)$SNP %>% unique() %>% length()
#       N.CS.x.HIT = subset(merge.HIT, Support>0 & FDR <= 0.05)$SNP %>% unique() %>% length()
#
#       N.Consensus = subset(merge.HIT, Consensus_SNP==T)$SNP %>% unique() %>% length()
#       N.Consensus.x.HIT = subset(merge.HIT, Consensus_SNP==T & FDR <= 0.05)$SNP %>% unique() %>% length()
#
#
#       fractions <- data.frame(N.leadSNP=N.leadSNP,
#                               percent.leadSNPs.x.HIT = N.leadSNPs.x.HIT/N.leadSNP,
#                               N.CS=N.CS,
#                               percent.CS.x.HIT = N.CS.x.HIT/N.CS,
#                               N.Consensus=N.Consensus,
#                               percent.Consensus.x.HIT = N.Consensus.x.HIT/N.Consensus)
#       fractions <- cbind(fractions, Gene=gene)
#       return(fractions)
#     }) %>%  data.table::rbindlist()
#
#     all.fractions <- cbind(all.fractions, QTL_type=QTL_type)
#     return(all.fractions)
#   }) %>% data.table::rbindlist()
#   return(QTL.results)
# }
#

#
# plot_QTL.results <- function(QTL.results){
#   QTL.results[is.na(QTL.results)] <- 0
#   per_cols <- grep("percent.",colnames(QTL.results), value = T)
#   mean_df <- data.table::melt(QTL.results,
#                               id.vars = c("Gene","QTL_type"),
#                               measure.vars = per_cols,
#                               variable.name = c("Measure"),
#                               value.name = c("Fraction"))
#   mean_df <- mean_df %>% dplyr::group_by(Measure, QTL_type) %>%
#     dplyr::summarise(Fraction = mean(Fraction) )
#
#   ggplot(mean_df,aes( x=QTL_type, y=Fraction,
#                       fill=Measure, color=Measure)) +
#     geom_col(position = "dodge")
# }
#

#
#
# percent.HIT.summary <- function(locus,
#                                 merged_dat,
#                                 QTL_type="eQTL"){
#   loc_path <- paste0("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/",locus,"/Multi-finemap/Multi-finemap_results.txt")
#   FM_results <- data.table::fread(loc_path)
#
#   dat <- data.frame(
#     GWAS_leadSNPs.total = length(subset(FM_results, leadSNP==T)$SNP %>% unique()),
#     GWAS_leadSNPs.HIT = length(subset(merged_dat, leadSNP==T)$SNP %>% unique()),
#     CS.total =  length(FM_results$SNP %>% unique()),
#     CS.HIT = length(merged_dat$SNP %>% unique()),
#     Consensus.total =  length(subset(FM_results, Consensus_SNP==T)$SNP %>% unique()),
#     Consensus.HIT = length(subset(merged_dat, Consensus_SNP==T)$SNP %>% unique())
#   ) %>% `row.names<-`("n SNPs")
#
#   percent_df <- dat %>% dplyr::mutate("% GWAS Lead SNPs" = round(GWAS_leadSNPs.HIT / GWAS_leadSNPs.total*100, 2),
#                                        "% Credible Set SNPs" = round(CS.HIT / CS.total*100, 2),
#                                        "% Consensus SNPs" = round(Consensus.HIT / Consensus.total*100, 2)) %>%
#     dplyr::select(c("% GWAS Lead SNPs","% Credible Set SNPs","% Consensus SNPs")) %>% t() %>% `colnames<-`("Percent SNPs")
#   percent_df <- data.frame(percent_df)
#   percent_df['SNP.Group'] <- row.names(percent_df)
#   percent_df["QTL_type"] <- QTL_type
#   percent_df["Locus"] <- locus
#   return(percent_df)
# }









# Attempts to download and manipulate files on server
# file_name = paste0("/sc/orga/projects/ad-omics/rajt/psych_encode/",files[1])
# output_dir = "./echolocatoR/tools/Annotations/psychENCODE/"
# dir.create(output_dir, recursive = T, showWarnings = F)
# process <- paste0(" | grep ", position," ")
#
# cmd <- paste0("scp -r schilb03@data4.hpc.mssm.edu:",file_name," ",output_dir)
#
# # cmd <- paste0("ssh -tt schilb03@data4.hpc.mssm.edu:/sc/orga/projects/ad-omics/rajt/psych_encode/",process, file_name)
# # library(RCurl)
# RCurl::scp(host = "data4.hpc.mssm.edu", user="schilb03", path = file_name)
# download.file(url = "schilb03@data4.hpc.mssm.edu:/sc/orga/projects/ad-omics/rajt/psych_encode/DER-08a_hg19_eQTL.significant.txt",
#               destfile = paste0(output_dir,files[2]))
#
# QTL <- system(cmd, intern = F)


