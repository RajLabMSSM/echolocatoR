
check.allel_direction <- function(){
  full.dat <- data.table::fread("./Data/GWAS/Nalls23andMe_2019/nallsEtAl2019_allSamples_allVariants.mod.txt",
                                nThread = 4)
  top.dat <- readxl::read_excel("./Data/GWAS/Nalls23andMe_2019/Nalls2019_TableS2.xlsx") %>% 
    data.table::data.table()
  full.dat <- subset(full.dat, RSID %in% top.dat$SNP) 
  
  merged.dat <- data.table:::merge.data.table(top.dat, full.dat, by.x = "SNP", by.y = "RSID")
  sum(merged.dat$A1==toupper(merged.dat$`Effect allele`)) / nrow(merged.dat)*100
  sum(merged.dat$A2==toupper(merged.dat$`Other allele`)) /  nrow(merged.dat)*100
  
  # A1 is always the Effect allele
  # A2 is always the Other allele
  
}



createDT <- function(DF, caption="", scrollY=400){
  data <- DT::datatable(DF, caption=caption,
                        extensions = 'Buttons',
                        options = list( dom = 'Bfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        scrollY = scrollY, scrollX=T, scrollCollapse = T, paging = F,
                                        columnDefs = list(list(className = 'dt-center', targets = "_all"))
                        )
  )
  return(data)
}

subset_genotype_data <- function(snp_list, 
                                 genotype_path = "volunteer_421.impute2.dosage"){  
  # NOTE: "volunteer_421.impute2.dosage" is the file you want to subset from! 
   cmd <- paste0("grep -E '",paste0(snp_list %>% unique(), collapse="|"),"' ", genotype_path)
   geno <- data.table::fread(text = system(cmd, intern = T),
                               sep = " ", header = F) 
  return(geno)
}


# Expression
## Get probe IDs for gene 
probes_mapping <- function(probe_path., gene_list){
  printer("++ Extracting probe info")
  cmd <- paste0("grep -E '", paste0(gene_list, collapse="|"),"' ",probe_path.)
  col_names <- data.table::fread(probe_path., sep="\t", nrows = 0) %>% colnames()
  probe_map <- data.table::fread(text = system(cmd, intern = T), 
                                 sep = "\t", header = F, col.names = col_names)
  if(dim(probe_map)[1]==0){
    stop("Could not identify gene in the probe mapping file: ",paste(gene_list, collapse=", "))
  }
  # Subset just to make sure you didn't accidentally grep other genes
  probe_map <- subset(probe_map, GENE %in% gene_list)
  return(probe_map)
}
## Subset expression data
get_expression_data <- function(expression_paths, gene., probe_path){
  printer("")
  printer("+ Processsing Expression data")
  # Find which probes for search for 
  probe_map <- probes_mapping(probe_path. = probe_path, gene_list = gene.)
  probe_list <- probe_map$PROBE_ID %>% unique()
  # Find the subjects that exist in all datasets
  # common_subjects <- lapply(expression_paths, function(x){
  #   col_names <- data.table::fread(x, sep="\t", header = T, nrows = 0) %>% colnames() 
  #   col_names <- col_names[-1]
  # }) 
  # Reduce(intersect, common_subjects)
  
  exp_dat <- lapply(expression_paths, function(x, gene.= gene, probe_list.=probe_list){
    condition = basename(dirname(x))
    printer("++",condition) 
    # grep rows with the probe names
    cmd <- paste0("gzcat ",x," | grep -E '", paste0(probe_list., collapse="|"),"'")
    col_names <- data.table::fread(x, sep="\t", header = T, nrows = 0) %>% colnames()
    exprs <- data.table::fread(text = system(cmd, intern = T), 
                               sep = "\t", header = F, col.names = col_names)
    ## Subset just to be sure contains probes
    exp_sub <- subset(exprs, PROBE_ID %in% probe_list.)
    # Add condition col
    exp_sub <- cbind( data.table::data.table(Condition = condition), exp_sub)
    # Melt subject IDs into single column
    subjects_IDs <- colnames(exp_sub)[!colnames(exp_sub) %in% c("Condition", "PROBE_ID")]
    exp_melt <- reshape2::melt(exp_sub, id.vars = c("Condition", "PROBE_ID"), 
                               measure.vars = subjects_IDs,
                               variable.name = "Subject_ID", 
                               value.name = "Expression")
  }) %>% data.table::rbindlist() 
  # Only use the probe with the highest average expression across individuals
  probe_means <- exp_dat %>% dplyr::group_by(PROBE_ID) %>% summarise(meanExp = mean(Expression))
  winning_probe <- subset(probe_means,  meanExp == max(meanExp))$PROBE_ID
  exp_dat <- subset(exp_dat, PROBE_ID == winning_probe)
  return(exp_dat)
} 

# eQTL SUMMARY STATS
get_SumStats_data <- function(eQTL_SS_paths, gene){
  printer("")
  printer("+ Processing Summary Stats data")
  SS <- lapply(eQTL_SS_paths, function(qx, gene.=gene){
    dat <- data.table::fread(qx, sep="\t", header=T) 
    condition <- basename(dirname(dirname(qx))) 
    dat <- dat %>% dplyr::rename(PROBE_ID="gene", Effect="beta", P="p-value")
    dat <- cbind(data.table::data.table(Condition = condition, Gene = gene.),
                 dat)
  }) %>% data.table::rbindlist()
  return(SS)
}

# GENOTYPE 
get_genotype_data <- function(genotype_path, .fam_path, subset_genotype_file = F, probe_ID. = NA){
  printer("")
  printer("+ Processing Genotype data") 
  if(subset_genotype_file){
    geno_subset <- subset_genotype_data(probe_ID. = probe_ID)
  } else { 
    geno_subset <- data.table::fread(genotype_path, sep=" ", header=F, stringsAsFactors = F)
  }
  # Add column names
  first_cols <- c("CHR","SNP","POS","A1","A2") 
  subject_IDs <- data.table::fread(.fam_path, stringsAsFactors = F)$V1
  colnames(geno_subset) <- c(first_cols,  subject_IDs) 
  ## Melt subject IDs into single column
  geno_melt <- geno_subset %>% reshape2::melt(geno_subset, id.vars = c("CHR","SNP","POS","A1","A2"),
                                              measure.vars = as.character(subject_IDs),
                                              variable.name = "Subject_ID", 
                                              value.name = "Genotype" ) %>% 
    dplyr::mutate(Genotype = round(as.numeric(Genotype),0) %>% as.factor())
  # Translate numeric genotypes to letters 
  # geno_melt <- geno_melt %>%  dplyr::rename(Genotype_int = Genotype) %>% 
  #   dplyr::mutate(Genotype = ifelse(Genotype_int == 0, paste0(A1,"/",A1),
  #                                   ifelse(Genotype_int == 1, paste0(A1,"/",A2),
  #                                          ifelse(Genotype_int == 2, paste0(A2,"/",A2), NA))) )
  return(geno_melt)
}


# MERGE: SS + EXP + GENO
merge_SS.EXP.GENO <- function(SS_data, geno_data, exp_data){
  printer("")
  printer("+ Merging Summary Stats, Genotype, and Expression data")
  ## SS + Genotype
  unique_SS <- unique(SS_data[,c("Gene","CHR","POS","Condition","Effect","t-stat","P","FDR")])
  SS_geno <- data.table:::merge.data.table(unique_SS,
                                           geno_data,
                                           by = c("CHR","POS"), 
                                           allow.cartesian = T)
  ## SS/Genotype + Expression
  SS_geno_exp <- data.table:::merge.data.table(SS_geno, 
                                               exp_data,
                                               by = c("Condition","Subject_ID")) 
  # SS_geno_exp <- dplyr::mutate(SS_geno_exp, Expression = Expression %>% as.numeric()) 
  return(SS_geno_exp)
}



##### Merge all eQTL files together into one file (for the selected SNPs) #####
# snp_list <- c("rs7294619","rs76904798","rs11175620")
merge_QTL_data <- function(snp_list,
                            eQTL_SS_paths = file.path("Data/QTL/Fairfax_2014",
                                                      c("CD14/LRRK2/LRRK2_Fairfax_CD14.txt",
                                                        "IFN/LRRK2/LRRK2_Fairfax_IFN.txt",
                                                        "LPS2/LRRK2/LRRK2_Fairfax_LPS2.txt",
                                                        "LPS24/LRRK2/LRRK2_Fairfax_LPS24.txt")),
                            expression_paths = file.path("Data/QTL/Fairfax_2014",
                                                         c("CD14/CD14.47231.414.b.txt.gz",
                                                           "IFN/IFN.47231.367.b.txt.gz",
                                                           "LPS2/LPS2.47231.261.b.txt.gz",
                                                           "LPS24/LPS24.47231.322.b.txt.gz")),
                            genotype_path = "Data/QTL/Fairfax_2014/geno.subset.txt", 
                            subset_genotype_file = F,
                            probe_path = "Data/QTL/Fairfax_2014/gene.ILMN.map",
                            .fam_path = "Data/QTL/Fairfax_2014/volunteers_421.fam",
                            gene = "LRRK2",
                            save_merged=T){ 
  # Expression 
  exp_data <- get_expression_data(expression_paths, gene, probe_path)

  # Summary Stats  
  SS_data <- get_SumStats_data(eQTL_SS_paths, gene)
  ## IMPORTANT! Subset according to the probes in the expression data.
  SS_data <- subset(SS_data, PROBE_ID == exp_data$PROBE_ID %>% unique() )
  
  # Genotype
  geno_data <- suppressWarnings(get_genotype_data(genotype_path, 
                                                   .fam_path, 
                                                   probe_ID. = unique(SS_data$PROBE_ID)) )
  # Merge x3
  SS_geno_exp <- merge_SS.EXP.GENO(SS_data, geno_data, exp_data)
  
  # Save separate standardized files in their respective folders
  if(save_merged){
    printer("++ Splitting and writing merged files to storage...")
    for(i in 1:length(unique(SS_geno_exp$Condition))){ 
      results_path <- dirname(eQTL_SS_paths[i])
      QTL.condition <- basename(dirname(results_path))
      subset_path <- get_subset_path(results_path, gene = gene, subset_path="auto")
      c_sub <- subset(SS_geno_exp, Condition==QTL.condition)
      printer("+++",subset_path)
      data.table::fwrite(c_sub, subset_path, quote = F, sep="\t")
    } 
  }
  return(SS_geno_exp)
}



eQTL_boxplots <- function(snp_list,
                         eQTL_SS_paths = file.path("Data/QTL/Fairfax_2014",
                                                   c("CD14/LRRK2/LRRK2_Fairfax_CD14.txt",
                                                     "IFN/LRRK2/LRRK2_Fairfax_IFN.txt",
                                                     "LPS2/LRRK2/LRRK2_Fairfax_LPS2.txt",
                                                     "LPS24/LRRK2/LRRK2_Fairfax_LPS24.txt")),
                         expression_paths = file.path("Data/QTL/Fairfax_2014",
                                                      c("CD14/CD14.47231.414.b.txt.gz",
                                                        "IFN/IFN.47231.367.b.txt.gz",
                                                        "LPS2/LPS2.47231.261.b.txt.gz",
                                                        "LPS24/LPS24.47231.322.b.txt.gz")),
                         genotype_path = "Data/QTL/Fairfax_2014/geno.subset.txt", 
                         subset_genotype_file = F,
                         probe_path = "Data/QTL/Fairfax_2014/gene.ILMN.map",
                         .fam_path = "Data/QTL/Fairfax_2014/volunteers_421.fam",
                         gene = "LRRK2",
                         show_plot = T,
                         SS_annotations = T,
                         interact = F,
                         save_merged = T){
  # Helper function
  printer <- function(..., v=T){if(v){print(paste(...))}}
  
 
  SS_geno_exp <- merge_QTL_data(snp_list=snp_list,
                                 eQTL_SS_paths=eQTL_SS_paths,
                                 expression_paths=expression_paths,
                                 genotype_path=genotype_path,
                                 subset_genotype_file=subset_genotype_file,
                                 probe_path=probe_path,
                                 .fam_path=.fam_path,
                                 gene=gene,
                                 save_merged=save_merged)
  
  if(show_plot){
    printer("")
    printer("+ Plotting eQTLs") 
    # encode genotypes
    DAT <- SS_geno_exp %>%
      mutate(genotype = case_when(Genotype == 0 ~ paste(A1,A1,sep="/"),
                                  Genotype == 1 ~ paste(A1,A2,sep="/"),
                                  Genotype == 2 ~ paste(A2,A2,sep="/")))
    
    expression.genotypes <- DAT %>% 
      dplyr::group_by(Gene, Condition, SNP, genotype, Genotype) %>%
      dplyr::summarise_each(c(Effect, P, FDR, Expression), funs = mean) 
    
    # Import PD GWAS data
    results_path <- "./Data/GWAS/Nalls23andMe_2019/LRRK2/"
    finemap_DT <- data.table::fread(file.path(results_path, "Multi-finemap/Multi-finemap_results.txt"))
    finemap_DT <- subset(finemap_DT, SNP %in% unique(SS_geno_exp$SNP)) 
    ## If the PD GWAS effect size is negative, flip the alleles and make the effect size positive.
    finemap_flip <- finemap_DT %>% 
      dplyr::mutate(Effect.flip = case_when(Effect != abs(Effect) ~ -Effect,
                                               Effect == abs(Effect) ~ Effect),
                        Risk.allele = case_when(Effect != abs(Effect) ~ A2,
                                               Effect == abs(Effect) ~ A1),
                        Non_risk.allele = case_when(Effect != abs(Effect) ~ A1,
                                               Effect == abs(Effect) ~ A2))
    finemap_flip <- finemap_flip[,c("SNP","MAF","Effect.flip","Risk.allele","Non_risk.allele")]
    
    # Merge QTL and GWAS
    DAT <- merge(DAT, finemap_flip, by="SNP")
    DAT <- DAT %>% dplyr::mutate(Risk.level = case_when(genotype == paste(Risk.allele,Risk.allele,sep="/") ~ "High",
                                                   genotype == paste(Non_risk.allele,Risk.allele,sep="/") |
                                                     genotype == paste(Risk.allele,Non_risk.allele,sep="/") ~ "Mid", 
                                                   genotype == paste(Non_risk.allele,Non_risk.allele,sep="/") ~ "Low"))
    
    DAT$Risk.level <- factor(DAT$Risk.level, levels=c("High","Mid","Low"), ordered = T)
    # DAT$MAF <- paste("MAF =",DAT$MAF)
     
      
    # Get 1 effect size per Condition x SNP combination
    # d <- 4
    
    labels <- DAT %>% 
      dplyr::group_by(Condition, SNP) %>% 
      dplyr::summarise(Effect=mean(Effect), FDR=mean(FDR), MAF=unique(MAF)) %>% 
      dplyr::mutate( Effect=round(Effect,4), FDR=formatC(FDR, format = "e", digits = 2)) %>% 
      dplyr::mutate(sig = ifelse(FDR<5e-9,"***",""))
    # labels <- subset(SS_geno_exp , select = c("CHR","POS","PROBE_ID","Condition","SNP","Effect","P","FDR")) %>%
    #   dplyr::mutate(FDR_sig = ifelse(as.numeric(FDR) < 1.34e-05, "FDR < 1.34e-05**",
    #                                  paste0("FDR = ", formatC(FDR, format = "e", digits = 2))),
    #                 FDR_scient = paste0("FDR = ", formatC(FDR, format = "e", digits = 2), 
    #                                     ifelse(as.numeric(FDR) < 1.34e-05, "**","") ),
    #                 P_sig = ifelse(as.numeric(P) < 0.05, "P < 0.05",paste0("P = ", formatC(P, format = "e", digits = 2))),
    #                 Beta = format(round(Effect,d), nsmall=d),
    #                 P = paste("P =",formatC(P, format = "e", digits = 2)),
    #                 FDR = paste("FDR =",formatC(FDR, format = "e", digits = 2))
    #                             ) %>%
    #   arrange(SNP, Condition) %>%
    #   unique()
    
    bp <- ggplot(data = DAT, aes(x = genotype, y = Expression, fill=Risk.level)) + 
      geom_boxplot(show.legend = T) + 
      geom_jitter(alpha=.5, width =.2,  show.legend = F) +  
      facet_grid(facets = Condition~SNP+MAF, scales = "free_x", drop = T) +
      theme_bw() + 
      theme(strip.text.y = element_text(angle = 0), 
            strip.text = element_text(colour = "white"),
            strip.background = element_rect(fill="grey10"),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5), 
            plot.background = element_rect(fill = "transparent"), 
            panel.background = element_rect(fill = "transparent")) + 
      # scale_fill_manual(values = c("red","yellow","green"))
      scale_fill_brewer(palette = "Spectral") + 
      labs(title=gene,subtitle = "PD Risk SNPs in Fairfax eQTL", x="Genotype") + 
      geom_text(data = labels, inherit.aes = F, size = 3, color="firebrick",
                aes(x = 2, y = 8.25,  label = paste0("Effect = ",Effect,"\nFDR = ",FDR," ",sig))) 
      # ylim(c(NA,max(SS_geno_exp$Expression)*1.1))
    print(bp)
    
      
    ggsave("./Data/QTL/Fairfax_2014/PD.Risk_Fairfax.eQTL.png", bp, dpi = 400, height = 12, width = 9)
    
    if(interact){
     print(plotly::ggplotly(bp)) 
    } else {
      print(bp) 
    }
   
  } 
  return(SS_geno_exp)
}
