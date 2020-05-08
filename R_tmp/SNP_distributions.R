
library(dplyr) 
library(ggplot2)


# 
# merged_results <- readxl::read_excel("~/Desktop/Fine_Mapping/Data/annotated_results_table.xlsx") 
# finemap_DT <- data.table::fread("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/LRRK2/Multi-finemap/Multi-finemap_results.txt")
# merged_results <- merge_finemapping_results(minimum_support=1,
#                                             include_leadSNPs=T,
#                                             xlsx_path="./Data/annotated_finemapping_results.xlsx",
#                                             from_storage=T,
#                                             haploreg_annotation=T,
#                                             biomart_annotation=T, 
#                                             verbose = F)

gather_locus_report <- function(){
  multi_files <- list.files(path = "./Data", pattern = "Multi-finemap_results.txt", 
                            recursive = T, full.names = T)
  # For each locus... 
  locus_report <- lapply(multi_files, function(f){
    locus <- basename(dirname(dirname(f)))
    printer("+ Extracting SNP count data for locus:", locus)
    dat <- data.table::fread(f)
    # Count signficiant SNPs from GWAS
    snp_count <- subset(dat, P<=5e-8) %>% count()
    report <- data.frame(Locus = locus,
                         Sig_SNPs.count = snp_count$n)
    # Count Credible Set SNPs
    CS_cols <- grep(pattern = "Credible_Set", x = colnames(dat), value = T)
    for(s in 1:length(CS_cols)){ 
      report[paste0("CS_support==",s)] <- (subset(dat, Support == s) %>% count())$n
    } 
    report["CS_SNP.count"] <-  (subset(dat, Support > 0) %>% count())$n
    
    # Count Consensus SNPs
    report$Consensus_SNP.count <- (subset(dat, Consensus_SNP==T) %>% count())$n 
    return(report)
  }) %>% data.table::rbindlist()
  return(locus_report)
}


 
# PLOTS
SNPgroups_plot <- function(locus_report){
  # GWAS Significant SNPs
  p1 <- ggplot(data=locus_report, aes(x=Sig_SNPs.count)) +
    geom_histogram(color="black", fill="grey") +
    labs(title = "Significant GWAS SNPs / Locus",
         x = "n SNPs", y = "Counts / Locus") +
    theme(plot.title = element_text(hjust = 0.5))  
  
  # Statistical Credible Sets
  support_cols <- grep("support==", colnames(locus_report), value = T)
  locus_melt <- data.table::melt(locus_report, id.vars=c("Locus","CS_SNP.count"),
                                 measure.vars = support_cols,
                                 variable.name = "Support_level",
                                 value.name = "Count" )
  p2 <- ggplot(data=subset(locus_melt, CS_SNP.count>0), aes(x=CS_SNP.count, fill=Support_level )) +
    geom_histogram(color="darkgreen") +
    labs(title = "Statistical Credible Set SNPs / Locus",
         x = "n SNPs", y = "Counts / Locus") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
    scale_fill_brewer(palette="Greens")
  
  
  # Functional Credible Sets
  p3 <- ggplot(data=subset(locus_melt, CS_SNP.count>0), aes(x=CS_SNP.count, fill=Support_level)) +
    geom_histogram(color="darkgreen") +
    labs(title = "Functional Credible Set SNPs / Locus",
         x = "n SNPs", y = "Counts / Locus") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
    scale_fill_brewer(palette="Greens")
  
  # Consensus SNPs
  p4 <- ggplot(data=locus_report, aes(x=Consensus_SNP.count )) +
    geom_histogram(fill="goldenrod2", color="goldenrod4") +
    labs(title = "Consensus SNPs / Locus",
         x = "n SNPs", y = "Counts / Locus") +
    theme(plot.title = element_text(hjust = 0.5))  
  
  
  # Annotations
  ## Credible Set Annotations
  p5 <- ggplot(data=subset(merged_results, Support > 0),
               aes(x=consequence_type_tv, fill=consequence_type_tv)) +  
    geom_histogram(stat="count",  color='darkblue') + 
    labs(title="Genic Annotations:\nCredible Set SNPs",
         x="Genic Annotation (Biomart)",
         y="n SNPs") + coord_flip() +  
    scale_fill_brewer(palette="Blues") +
    theme_classic() + 
    theme(legend.position = "none")
  
  ## Consensus SNP Annotations
  CS_cols <- grep("Credible_Set",colnames(merged_results), value = T)
  # merged_results$consequence_type_tv <- as.character(merged_results$consequence_type_tv)
  p6 <- ggplot(data=subset(merged_results, Support == length(CS_cols), drop = F),
               aes(x=consequence_type_tv, fill=consequence_type_tv))+  
    geom_histogram(stat="count",  color='darkblue') + 
    labs(title="Genic Annotations:\nConsensus SNPs",
         x="Genic Annotation (Biomart)",
         y="n SNPs")+ coord_flip() +  
    scale_fill_brewer(palette="Blues") + 
    theme_classic() + 
    theme(legend.position = "none")
  
  # Merge plots
  cowplot::plot_grid(p1, p4, p2, p3, p5, p6, ncol = 2, labels = c("A","B", "C","D", "E","F"))
  
}


# Plot 1
## Break counts into groups
# break_cut <- function(locus_report, 
#                       column="Sig_SNPs.count", 
#                       column.new="Sig_SNPs.groups", 
#                       cutsize=100, 
#                       max_value=NA){
#   if(is.na(max_value)){
#     max_value <- round( max(locus_report[,column]), digits = -3)
#   }
#   locus_report <- as.data.frame(locus_report)
#   cuts <- cut_interval(locus_report[,column], length=cutsize)
#   # breaks <- seq(0, max_value, by = cutsize)   
#   # cuts <- Hmisc::cut2(locus_report[,column], cuts = breaks, include.lowest=T, right=F)
#   # e=cuts[3]
#   intervals <- lapply(cuts, function(e){  
#       splitter <- suppressWarnings( strsplit(toString(e), ",|)|\\(|\\]|\\[|] ")[[1]] %>% as.integer())
#       # print(splitter)
#       numbers <- splitter[!is.na(splitter)]
#       interval <- paste0(numbers[1] ,"-",numbers[2]) 
#       # print(interval)
#       return(interval)
#   }) %>% unlist()
#   locus_report[column.new] <- as.factor(intervals) 
#   return(locus_report)
# } 
# locus_report  <- break_cut(locus_report, 
#                            column="Sig_SNPs.count", 
#                            column.new="Sig_SNPs.groups", 
#                            cutsize=100)
# locus_report  <- break_cut(locus_report, 
#                            column="CS_SNP.count", 
#                            column.new="CS_SNP.groups", 
#                            cutsize=3,
#                            max_value = max(locus_report$Sig_SNPs.count)) 

# p2 <- ggplot(data=subset(locus_report, CS_SNP.count>0), aes(x=0, fill=CS_SNP.groups)) +
#   geom_histogram(color="white") +
#   labs(title = "Credible Set SNPs / Locus",
#        x = "Count", y = "") +
#   theme(plot.title = element_text(hjust = 0.5)) + coord_flip()

# p1 <- ggplot() + 
#   geom_histogram(data=locus_report, aes(x=0, fill=Sig_SNPs.groups), 
#                  color="white", binwidth = 50) + 
#   labs(title = "Significant GWAS SNPs / Locus",
#        x = "Count", y = "") + 
#   theme(plot.title = element_text(hjust = 0.5)) + coord_flip()

