
# ********************************
# ************ fGWAS# ************ 
# ********************************
# """
# The R version of fGWAS appears to be extremely limited. 
# https://github.com/wzhy2000/fGWAS
# 
# Using the command line version instead.
# https://github.com/joepickrell/fgwas
## User Guide:
## https://github.com/joepickrell/fgwas/blob/master/man/fgwas_manual.pdf

# Pre-formatted annotations available here:
# https://github.com/joepickrell/1000-genomes

# """

fGWAS.install <- function(echo_path="./echolocatoR/tools", 
                          download.annotations=F){
  printer("fGWAS:: Installing and compiling fGWAS.")
  # Install R wrapper
  # devtools::install_github("wzhy2000/fGWAS/pkg")
  # # Install and compile tool
  # system(paste0("git submodule add  https://github.com/wzhy2000/fGWAS.git ", echo_path))
  # system( paste0("cd ",file.path(echo_path,"fGWAS")," & R CMD INSTALL pkg") )w
  tar.gz <- "0.3.6.tar.gz"
  tar <- file.path(echo_path,gsub(".gz","",tar.gz))
  # Download from github
  system( paste0("wget ",
                file.path("https://github.com/joepickrell/fgwas/archive/",tar.gz),
                " ",echo_path) )
  # Decompress
  system(paste0("tar -xf ",tar))
  # Compile
  fgwas.folder <- file.path(echo_path,paste0("fgwas-",gsub(".tar.gz","",tar.gz)) )
  cmd <- paste0("cd ",fgwas.folder," & ./configure & make")
  system(cmd)
  # Create symlink?...
  if(download.annotations){
    fGWAS.download_annotations()
  }
}


fGWAS.download_annotations <- function(FM_all,
                                       results_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide",
                                       force_new_annot = F,
                                       dataset = ""
                                       ){
  # Path/file names
  Input_path <- file.path(results_path,"fGWAS","Input")
  dir.create(Input_path, showWarnings = F, recursive = T)
  annot.file.name <- file.path(Input_path, paste0("annotations.",basename(dataset),".txt")) 
  
  # Takes a long time: Redo only if it doesn't exist (or its being forced to)
  if(!file.exists(annot.file.name) | force_new_annot==T){
    printer("fGWAS: Downloading annotations from GitHub repo:")
    printer("       https://github.com/joepickrell/1000-genomes") 
    # system(paste("git submodule add https://github.com/joepickrell/1000-genomes.git",annot.folder))  
    base.url <- "https://github.com/joepickrell/1000-genomes/raw/master"
    FM.snps <- unique(FM_all$SNP) #paste(unique(FM_all$SNP), collapse="|")
    # Download each annotation and merge with data
    annots.filt <- lapply(unique(FM_all$CHR), function(chr){
      printer("fGWAS: Downloading annotations for Chrom",chr)
      file.url <- file.path(base.url, paste0("chr",chr,".annot.wdist.wcoding.gz"))
      chr.dat <- data.table::fread(file.url)
      CHR.dat <- subset(chr.dat, rs %in% FM.snps) 
      return(CHR.dat)
    }) %>% data.table::rbindlist()
    # Save filtered annotations 
    data.table::fwrite(annots.filt, annot.file.name, sep="\t")
  } else {
    printer("+ fGWAS: Existing annotations found.")
    printer("    Importing:",annot.file.name)
    annots.filt <- data.table::fread(annot.file.name, nThread = 4)
    } 
  # Merge annotations afterwards
  FM_annot <- data.table:::merge.data.table(FM_all, annots.filt, 
                                            by.x = "SNP",
                                            by.y = "rs", 
                                            all.x = T)
  # Fill NAs with 0s in certain columns
  for (j in names(FM_annot)[24:ncol(FM_annot)]){
    set(FM_annot,which(is.na(FM_annot[[j]])),j,0)
  }  
  return(FM_annot)
}



fGWAS.annotation_names <- function(fgwas="./echolocatoR/tools/fgwas-0.3.6/src/fgwas"){ 
  ## Get annotation names from summary file
  annot_files <- data.table::fread(file.path( dirname(dirname(fgwas)),
                                              "annot", "annotation_list_winfo.txt"),
                                   col.names = c('bed.name',"type","description"))
  annot_files$Name <- gsub(".bed.gz","",basename(annot_files$bed.name)) 
  annot_files$Source <- lapply(annot_files$bed.name, function(s){strsplit(s,"/")[[1]][1]}) %>% unlist() 
  annot_files <- tidyr::separate(annot_files, col = "Name", 
                                 into=c("Annot"), 
                                 sep="\\.",remove=F, extra="drop")
  annot_files$Annot <- make.unique(annot_files$Annot)
  return(annot_files)
}

fGWAS.top_annotations <- function(dat.fgwas, annot_files, SNP.Group=""){
  ## Count how many SNPs have hits per annotation
  overlapping.snps <- subset(dat.fgwas, select=as.character(annot_files$Name) ) %>% 
    colSums() 
  printer("fGWAS: Annotations with the most overlapping SNPs:",SNP.Group)
  print(head(overlapping.snps %>% sort(decreasing = T), 5))
  overlap.df <- data.frame(overlapping.snps) %>% `colnames<-`(SNP.Group) 
  # data.table::data.table(overlap.df, keep.rownames = "Annot", key="Annot")
  return(overlap.df)
}
 


fGWAS.prepare_input <- function(FM_annot,
                                results_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide",
                                SNP.Groups = c("Consensus", "CredibleSet", "GWAS.lead", "Selected", "Random"), 
                                dataset = "./Data/GWAS/Nalls23andMe_2019", 
                                selected_SNPs = F,
                                random.seed = F,
                                Locus=NA
                                ){
  printer("fGWAS: Preparing input files.")
  # """
  # 1. SNPID: a SNP identifier
  # 2. CHR: the chromosome for the SNP
  # 3. POS: the position of the SNP
  # 4. F: the allele frequency of one of the alleles of the SNP
  # 5. Z: the Z-score for from the association study
  # 6. N: the sample size used in the association study at this SNP (this can vary from SNP to SNP due to, e.g., missing data).
  # """ 
  dir.create(file.path(results_path, "fGWAS"), showWarnings = F, recursive = F) 
  # List annotation names
  annot_files <- fGWAS.annotation_names()
  
  Locus_counts <- dat.fgwas %>% dplyr::group_by(Gene) %>% tally() 
  block_size <- min(subset(Locus_counts, n>10)$n)
  
  # Construct input for each SNP.Group
  fgwas.inputs <- data.frame()
  overlap.df.list <- list() 
  for(group in SNP.Groups){
    printer("+ fGWAS: Creating", group, "file") 
    if(group=="GWAS"){
      dat.fgwas <- FM_annot
    }
    if(group=="Multi-finemap"){
      dat.fgwas <- FM_annot
      dat.fgwas$Effect <- rowMeans(subset(FM_annot, select=grep(".PP",colnames(FM_annot))))
      # FM_merge$Adjusted.Effect <- FM_merge$Effect * FM_merge$mean.PP
    }
    # Subset according to which group we're looking at
    if(group=="Consensus"){
      dat.fgwas <- subset(FM_annot, Consensus_SNP==T)
    }
    if(group=="CredibleSet"){
      dat.fgwas <- subset(FM_annot, Support>0)
    }
    if(group=="GWAS.lead"){
      dat.fgwas <- subset(FM_annot, leadSNP==T)
    }
    if(group=="Selected"){
      if(any(selected_SNPs!=F)){
        dat.fgwas <- subset(FM_annot, SNP %in% selected_SNPs)
      } else { 
        stop("+ fGWAS:: Please provide a list of SNPs to the argument 'selected_SNPs',",
             "or remove 'Selected SNPs' from the SNP.Groups argument.")
      } 
    }
    if(group=="Random"){
      CS.size <- nrow(subset(FM_annot, Consensus_SNP==T))
      if(random.seed!=F){set.seed(random.seed)}
      dat.fgwas <- FM_annot[sample(nrow(FM_annot), CS.size), ] 
    }
  
    # Construct data in fGWAS format
    dat.fgwas <- dat.fgwas  %>%
      dplyr::mutate(N=N_cases+N_controls) %>%
      dplyr::rename(SNPID=SNP, 
                    CHR=CHR, 
                    POS=POS, 
                    "F"=Freq, 
                    Z=Effect, 
                    N=N, 
                    NCASE=N_cases, 
                    NCONTROL=N_controls, 
                    Gene=Gene) %>% 
      arrange(POS)  
    # Assign each locus a segment ID
    if(!is.na(Locus)){ dat.fgwas <- subset(dat.fgwas, Gene==Locus)}
    seg.table <- data.table::data.table(Gene = unique(FM_annot$Gene), 
                                        SEGNUMBER = (1:length(unique(FM_annot$Gene)))+1 )
    dat.fgwas <- data.table:::merge.data.table(dat.fgwas,
                                               seg.table,
                                               by = "Gene")
    dat.fgwas <- unique(dat.fgwas) %>% arrange(SEGNUMBER, CHR, POS) 
    dat.fgwas$SEGNUMBER <- as.numeric(dat.fgwas$SEGNUMBER)
    
    
    # Save file
    dat.path <- file.path(results_path,"fGWAS","Input",paste0("dat.fgwas.",group,".txt"))
    dir.create(dirname(dat.path), showWarnings = F, recursive = T) 
    data.table::fwrite(dat.fgwas, dat.path, sep = " ", quote = F, nThread = 4)
    gzip(dat.path, overwrite=T) 
    # Add to summary dataframe
    fgwas.inputs <- rbind(fgwas.inputs, data.frame(SNP.Group=group,
                                                   File=paste0(dat.path,".gz"),
                                                   N.SNPs=block_size))
   
    # Report annotations with the most overlap
    overlap.df <- fGWAS.top_annotations(dat.fgwas, annot_files, SNP.Group=group)
    overlap.df.list <- append(overlap.df.list, overlap.df)
  } 
  overlap.DF <- cbind.data.frame(overlap.df.list) %>% `row.names<-`(rownames(overlap.df)) 
  return(list(fgwas.inputs=fgwas.inputs, 
              overlap.DF=overlap.DF))
}

 

fGWAS.gather_results <- function(results_path, 
                                 output_dir=file.path(results_path,"fGWAS/Output")
                                 ){
  printer("+ fGWAS: Gathering all results...")
  #### Default output:
  # 1. fgwas.params. The maximum likelihood parameter estimates for each parameter in the model. The columns are the name of the parameter (“pi” is the parameter for the prior probability that any given genomic region contains an association), the lower bound of the 95% confidence interval on the parameter, the maximum likelihood estimate of the parameter, and the upper bound of the 95% confidence interval on the parameter.
  # 2. fgwas.llk. The likelihood of the model. The lines are the log-likelihood of the model under the maximum likelihood parameter estimates, the number of parameters, and the AIC.
  
#### With "-print" flag on:
  # 1.fgwas.ridgeparams. The estimates of the parameters under the penalized likelihood. The first line of this file is the penalty used, then the penalized parameters estimates follow. If you have used the -xv flag, the last line will be the cross-validation likelihood.
  # 2. fgwas.segbfs.gz. The association statistics in each region of the genome defined in the model. The columns of this file are the block number, chromosome, start position, end posi- tion, maximum absolute value of Z-score, log regional Bayes factor, regional prior probability of association, log regional posterior odds for association, and the regional posterior proba- bility of association. The annotations of the region (if any) are in the remaining columns.
  # 3. fgwas.bfs.gz. The association statistics for each SNP in the genome as estimated by the model. The columns of this file are the SNP ID, the chromosome, genomic position, log Bayes factor, Z-score, estimated variance in the effect size, prior probability of association, two columns (pseudologPO and pseudoPPA) for internal use only, the posterior probability of association (conditional on there being an association in the region), and the region number. The annotations in the model (if any) then follow.
  
  # Gather file names
  llk.files <- list.files(output_dir, pattern = ".llk")
  params.files <-  list.files(output_dir, pattern = ".params")
  # Set up preliminary df
  df <- data.frame(llk.file=llk.files, 
                   params.file=params.files)
  df <- tidyr::separate(df, col = "llk.file", into=c("SNP.Group","Annot",NA), 
                        sep="__|\\.",remove=F, extra="drop")
  # Gather results stats for each file
  res_df <- lapply(df$llk.file, function(llk){
    # llk file
    res <- data.table::fread(file.path(output_dir,llk), nThread = 4) %>% data.table::transpose(fill=0)
    colnames(res) <- as.character(res[1,])
    res <- res[-1,]
    res$llk.file <- llk
    # params file
    param <- gsub(".llk",".params",llk)
    res.p <- data.table::fread(file.path(output_dir, param), nThread = 4)
    res$parameter <- res.p$parameter[-1]
    res$CI_lo <- res.p$CI_lo[-1]
    res$CI_hi <- res.p$CI_hi[-1]
    res$estimate <- res.p$estimate[-1]  
    return(res)
  }) %>% data.table::rbindlist(fill=T)
  # Merge stats with preliminary df
  DF <- data.table:::merge.data.table(df, res_df, by = "llk.file")
  printer("+ Data for",nrow(DF),"results gathered.")
  return(DF)
}


fGWAS.run <- function(results_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide",
                      fgwas = "./echolocatoR/tools/fgwas-0.3.6/src/fgwas",
                      overlap.DF,
                      fgwas.inputs){
  # Iterate over each annotation file,
  ## so you can figure out where are significant.
  fgwas.out = file.path(results_path,"fGWAS","Output")
  dir.create(file.path(fgwas.out), showWarnings = F, recursive = T) 
  fg.start <- Sys.time()
  for(annot in rownames(overlap.DF)){
    printer("")
    printer("------------------")
    printer("------------------")
    printer("+ fGWAS: Using annotation =",annot)
    # Iterate over each SNP.group file, 
    ## so you can compare their enrichment results to one another. 
    for(i in 1:nrow(fgwas.inputs)){ 
      # Gather info
      group <- fgwas.inputs[i,"SNP.Group"]
      f <- fgwas.inputs[i,"File"]
      N.SNPs <- fgwas.inputs[i,"N.SNPs"]
      message("++ fGWAS: Running enrichment on SNP.Group:  ",group)
      # Construct command
      cmd <- paste0(
        fgwas,
        " -i ", f,
        " -cc ", # Suppposed to have this on but makes everything super slow and doesn't seem to affect results.
        " -o ",file.path(fgwas.out, paste0(group,"__", annot)), 
        " -fine",
        # " -xv", # Get cross-validated likelihood
        " -w ",annot,
        # " -print", # Produces extra files
        " -k ",N.SNPs
        # " -onlyp -print" # Turn this on to avoid "WARNING: failed to converge" message.
      ) 
      system(cmd) 
    } 
  }
  fg.end <- Sys.time()
  printer("fGWAS:",nrow(overlap.DF)*ncol(overlap.DF),"enrichment tests run in",
          round(fg.end-fg.start,2))
}

 


fGWAS <- function(results_path = "./Data/GWAS/Nalls23andMe_2019/_genome_wide", 
                  fgwas = "./echolocatoR/tools/fgwas-0.3.6/src/fgwas",
                  dataset = "./Data/GWAS/Nalls23andMe_2019",
                  SNP.Groups = c("GWAS","Multi-finemap"),#c("Consensus", "CredibleSet", "GWAS.lead", "Selected", "Random"),
                  selected_SNPs = F,
                  remove_tmps = T,
                  force_new_fgwas = F,
                  force_new_annot = F,
                  random.seed = F,
                  random.iterations = 100){
  fgwas.out = file.path(results_path,"fGWAS","Output")
  # [0] Check if results already exist for this dataset
  output.summary <- file.path(results_path,"fGWAS",paste0("fGWAS_summary.",basename(dataset),".txt"))
  if(file.exists(output.summary) & force_new_fgwas==F){
    printer("fGWAS:: Results file already exists.")
    printer("  Importing: ",output.summary)
    RESULTS.DF <- data.table::fread(output.summary)
  } else { 
    # [1] Create FM annot file
    print("fGWAS:: Gathering fine-mapping results from all loci and merging into one data.table.")
    FM_all <- merge_finemapping_results(minimum_support = 0, 
                                        include_leadSNPs = T, 
                                        verbose = F) 
    FM_all <- subset(FM_all, Dataset==dataset)
    ## Gather annotations data and merge with FM_all
    FM_annot <- fGWAS.download_annotations(FM_all,
                                           force_new_annot = force_new_annot, 
                                           dataset = dataset)
    
    # [2] Prepare inputs
    if("Random" == paste(SNP.Groups, collapse="")){
      printer("+ fGWAS:: Conducting",random.iterations,"for Random SNPs.")
    } else if("Random" %in% SNP.Groups){
        printer("+ fGWAS:: More than one SNP Group specified. Only conducting one Random iteration.")
        random.iterations <- 1
    }
    RESULTS.DF <- lapply(1:random.iterations, function(iter){
      printer("+ fGWAS:: Iteration =",iter)
      prepare_input.list <- fGWAS.prepare_input(FM_annot = FM_annot,
                                                results_path = results_path,
                                                SNP.Groups = SNP.Groups,
                                                selected_SNPs = selected_SNPs,
                                                dataset = dataset,
                                                random.seed = random.seed,
                                                Locus = NA) 
      fgwas.inputs <- prepare_input.list$fgwas.inputs
      overlap.DF <- prepare_input.list$overlap.DF   
      
      # [3] Run fGWAS
      fGWAS.run(results_path=results_path,
                overlap.DF=overlap.DF,
                fgwas.inputs=fgwas.inputs) 
      
      # [4] Gather results
      printer("fGWAS:: Saving all fGWAS results to ==>",output.summary)
      results.DF <- fGWAS.gather_results(results_path = results_path) 
      data.table::fwrite(results.DF, output.summary, sep="\t", nThread = 4)
      # wut <- data.table::fread( "./Data/GWAS/Nalls23andMe_2019/_genome_wide/fGWAS/fGWAS_summary.Nalls23andMe_2019.txt")
      
      # [5] Delete tmp files (input/output)
      if(remove_tmps){
        printer("fGWAS: Removing tmp input files...")
        # Remove Input  
        suppressWarnings(file.remove(list.files(file.path(dirname(fgwas.out),"Input"), pattern = ".gz")))
        # Remove Output
        removeDirectory(fgwas.out, recursive = T, mustExist = F) 
      } 
      return(results.DF)
    }) %>% data.table::rbindlist(fill=T) 
  } 
  
  # annot.info <- fGWAS.annotation_names()
  # RESULTS.DF <- data.table:::merge.data.table(RESULTS.DF,
  #                               annot.info[,c("Annot","type","Source","description")], by="Annot")
  return(RESULTS.DF)
}






##### PLOTS #######

# Boxplot
fGWAS.boxplot <- function(results.DF, 
                          title="fGWAS Enrichment Results", 
                          subtitle = "451 Annotations", 
                          show_plot = T, 
                          interact = T){
  DF <- results.DF
  DF$Enrichment <- ifelse(results.DF$estimate==0, "None", 
                     ifelse(results.DF$estimate>0, "Enriched","Depleted"))
  DF$Enrichment <- factor(DF$Enrichment, levels = c("Enriched","Depleted","None"),
                     labels = c("Enriched","Depleted","None"), 
                     ordered = T) 
  # counts <- DF %>% dplyr::group_by(SNP.Group, Enrichment, .drop=F) %>% count() %>% 
  #   arrange(Enrichment, SNP.Group) %>% 
  #   subset(Enrichment!="None")
  # x.ticks <- paste0(counts$SNP.Group,"\n(n=",counts$n,")")
  
  stat_box_data <- function(y, lower_limit = min(DF$estimate) * 1.15) {
    return( 
      data.frame(
        y = 0.95 * lower_limit,
        label = paste0('n=', length(y))
      )
    )
  }
  
  bp <- ggplot(DF, aes(x=SNP.Group, y=estimate, fill=SNP.Group, 
                       text = paste("Annotation:",Annot))) +
    geom_point(show.legend = F, aes(alpha=0.7)) +
    geom_jitter(width = .2, show.legend = F) +
    geom_boxplot(show.legend = F, aes(alpha=0.7)) + 
    facet_grid("~Enrichment") +
    # stat_summary(
    #   fun.data = stat_box_data, 
    #   geom = "text", 
    #   hjust = 0.5,
    #   vjust = 0.9
    # ) + 
    theme_classic() + 
    scale_fill_brewer(palette = "Dark2") + 
    labs(title=title, subtitle = subtitle) + 
    theme(plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5))  
  if(interact){
    bp <- plotly::ggplotly(bp)
    # bp <- htmltools::tagList(list(bp))
  }
  if(show_plot){print(bp)}
  return(bp)
}


# fGWAS.line_plot_consensus <- function(RESULTS.DF, top_n=10){
#   RESULTS.DF %>% dplyr::group_by(SNP.Group, type)
#   DF <- (subset(RESULTS.DF, SNP.Group=="Consensus") %>% dplyr::arrange(desc(estimate)))[1:top_n,]
#   fGWAS.line_plot(DF)
#   
# }



fGWAS.line_plot <- function(RESULTS.DF, 
                            show_plot=T, 
                            remove_random=T, 
                            results_path, 
                            save_plot=T,
                            top_annots=F){
  DF <- RESULTS.DF
  if(remove_random){ DF <- subset(DF, !(SNP.Group %in% c("Random", "Selected"))) %>% droplevels() }
  annot_files <- fGWAS.annotation_names()
  annot_files[annot_files$Source %in% c("syn","nonsyn", "ensembl_genes", "segmentation"),"Source"] <- "Ensembl"
  annot_files[annot_files$Source == 'stam_dnase',"Source"] <- "Stam Lab"
  annot_files$type <- paste0(annot_files$Source,"-", annot_files$type)
  annot.df <- data.table:::merge.data.table(DF, annot_files, 
                                            by = "Annot",
                                            all.x = T) 
  annot.df$CI_low <-  as.numeric( gsub("<|>","",annot.df$CI_lo))
  annot.df$CI_high <-  as.numeric( gsub("<|>","",annot.df$CI_hi)) # "fail" automatically converted to NA
  annot.df$Short_Description <- gsub(" *\\(.*?\\)", "", annot.df$description) 
  annot.df$Short_Description <- gsub(" *\\(.*?\\)", "", annot.df$description) 
  
  annot.df$Direction <- ifelse(abs(annot.df$estimate)==annot.df$estimate, "+","-") 
  annot.df$Type <- factor(annot.df$type,  
                          levels = unique(annot.df$type))
  annot.df$dummy <- ""
 
  # Order by the top Consensus SNP annotations
  annot.df <- dplyr::rename(annot.df, AIC = "AIC:")
  annot.SORT <- annot.df %>%
    subset(SNP.Group == "Multi-finemap") %>% 
    dplyr::group_by(Type) %>%
    arrange(estimate, desc(AIC))
  
  
  # Make variable ordered factors
  annot.df$Short_Description <- factor(annot.df$Short_Description,  
                                       levels = unique(annot.SORT$Short_Description)) 
  annot.df$description <- factor(annot.df$description,  
                                       levels = unique(annot.SORT$description)) 
  annot.df <- arrange(annot.df, estimate, desc(AIC))
  
  # if(top_annots != F){
  #   annot.df <- annot.df %>% dplyr::group_by(factor(SNP.Group)) %>% 
  #                             top_n(n=5, wt=estimate) %>% 
  #                             data.table::data.table()
  #   ( subset(annot.df, SNP.Group %in% c("Consensus")) %>%  
  #       arrange(desc(estimate), desc(AIC)) )[1:top_annots,]
  # }
  
  # annot.df <- annot.df %>% dplyr::group_by(Type, Short_Description, SNP.Group, Direction) %>%
  #   dplyr::summarise_at(.vars=c("estimate","CI_low","CI_high","ln(lk):","AIC:"), .funs = mean)
  # cor(annot.means$`ln(lk):`, annot.means$estimate)
  # Plot
  ## Point w/ CIs
  point.plot <- ggplot(annot.df, aes(x=log2(estimate),  y=description, fill=SNP.Group, color=SNP.Group)) + 
    # xlim(c(-max(abs(annot.df$estimate)), max(abs(annot.df$estimate)) )) +
    geom_point(show.legend = F, alpha = 0.5) +
    geom_errorbarh(aes(xmin=CI_low, xmax=CI_high), show.legend = F, alpha = 0.5 ) + 
    theme_minimal() + 
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + 
    facet_grid(~SNP.Group, drop = T) +
    labs(title = "fGWAS Results by SNP Group", x="log2(enrichment)") + 
    scale_color_manual(values = c("Consensus"="goldenrod3",
                                   "CredibleSet"="green3",
                                   "GWAS"="red2",
                                   "Selected"="purple2"))
  
  ## Tile labels
  tile.plot <- ggplot(annot.df, aes(x=" ", 
                                    y=description,#Short_Description, 
                                    fill=Type, width=0.9, height=0.9), size=1.5) + 
    geom_tile() + 
    theme_minimal() + 
    theme(legend.position="left")+#, strip.text.x = element_text(angle=90)) +
    scale_color_brewer(palette="Dark2") +
    scale_fill_brewer(palette="Dark2") + 
    labs(title="", x="") + 
    facet_grid(~dummy)
  # Merge plots
  cp <- cowplot::plot_grid(tile.plot, point.plot, nrow = 1, rel_widths = c(1/2, 1)) 
  # gridExtra::grid.arrange(tile.plot, point.plot, nrow = 1)
  if(show_plot){print(cp)}
  
  if(save_plot){
    cowplot::save_plot(file.path(results_path,"fGWAS","lineplot.png"), cp, 
                       base_height = 15, base_width = 13)
  }
  
  return(cp)
}


fGWAS.heatmap <- function(RESULTS.DF, 
                          annotation_or_tissue="tissue",
                          show_plot=T){
  DF <- RESULTS.DF
  # library(heatmaply)
  annot_files <- fGWAS.annotation_names()
  colors <- RColorBrewer::brewer.pal(11,"Spectral")  
  if(annotation_or_tissue=="annotation"){
   
    ## By each annotation
    mat <- reshape2::acast(DF, Annot~SNP.Group, value.var="estimate", 
                           fun.aggregate = mean, drop = F, fill = 0)  
    # Add annotation metadata
    mat.annot <- data.table:::merge.data.table(
      data.table(mat, keep.rownames = "Annot", key="Annot"),
      annot_files, 
      by = "Annot",
      all.x = T) 
    mat.annot <- data.frame(mat.annot[,-c("bed.name","Name")], 
                            row.names = mat.annot$Annot) 
    hm <- heatmaply::heatmaply(mat.annot, 
              height = 10, width = 5,
              cexRow = .7,
              # column_text_angle = 0,
              main = "fGWAS Enrichment: by Annotation",
              ylab = "Annotation", 
              xlab = "SNP Group", 
              key.title = "Estimate", 
              dendrogram = "row", 
              k_row = 5,
              colors = colors) 
  } else {
    ## By Tissue 
    DF.annot <- data.table:::merge.data.table(DF, annot_files, 
                                              by = "Annot",
                                              all.x = T)
    mat <- reshape2::acast(DF.annot, description~SNP.Group, value.var="estimate", 
                           fun.aggregate = mean, drop = F) 
    hm <- heatmaply::heatmaply(mat,  
              height = 10, width = 5,
              cexRow = .7,
              # column_text_angle = 0,
              main = "fGWAS Enrichment: by Tissue",
              ylab = "Tissue", 
              xlab = "SNP Group", 
              key.title = "Estimate", 
              dendrogram = "row",
              column_text_angle = 0,
              k_row = 5,
              colors = colors)
  }
  # Print and return plot
  # hm <- htmltools::tagList(list(hm))
  if(show_plot){print(hm)} 
  return(hm)
}

fGWAS.plots <- function(results.DF){
  DF <- results.DF 
  # Boxplot
  bp <- fGWAS.boxplot(DF,
                      title="fGWAS Enrichment Results", 
                      subtitle = "451 Annotations")
  print(bp)
  
  # Heatmap
  hm <- fGWAS.heatmap(DF, annotation_or_tissue="tissue")
  print(hm)
}

fGWAS.estimate_summary <- function(RESULTS.DF){
  thresh <- 1
  DF <- RESULTS.DF
  sig.subset <- subset(DF, (estimate>thresh) & (SNP.Group!="Random"))
  consensus.sig.subset <- sig.subset %>% group_by(Annot) %>% 
    subset((abs(estimate) == max(abs(estimate)) ) & (SNP.Group=="GWAS"))
  
  percent <- round(nrow(consensus.sig.subset)/nrow(sig.subset)*100,2)
  printer("fGWAS:: Consensus SNPs had the largest absolute estimate in",
          nrow(consensus.sig.subset),"/",nrow(sig.subset),"(",percent,"%)",
          "tests where the estimate was >",thresh)
  
  
  consensus <- subset(DF, SNP.Group =="Consensus")
  consensus_per <- round(nrow(subset(consensus, abs(estimate) >thresh)) / nrow(consensus)*100, 2) 
  leadGWAS <- subset(DF, SNP.Group =="GWAS")
  leadGWAS_per <- round(nrow(subset(leadGWAS, abs(estimate) >thresh)) / nrow(leadGWAS)*100, 2) 
  
  
}

# dat.path <- file.path(results_path,"fGWAS/fgwas.data.txt")
# data.table::fwrite(dat.fgwas, dat.path, sep="\t", 
#                    row.names = F, col.names = T)
# 
# r <- fGWAS::fg.load.simple(file.simple.snp = dat.path )





