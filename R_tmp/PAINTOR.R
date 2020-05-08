# [GitHub](https://github.com/gkichaev/PAINTOR_V3.0)
# [LD Tutorial](https://github.com/gkichaev/PAINTOR_V3.0/wiki/2a.-Computing-1000-genomes-LD)

PAINTOR.install <- function(paintor_path="./echolocatoR/tools/PAINTOR_V3.0/"){
  cmd <- paste("cd", paintor_path,"&& bash install")
  system(cmd)
}

PAINTOR.help_menu <- function(paintor_path="./echolocatoR/tools/PAINTOR_V3.0/"){
  system(file.path(paintor_path,"PAINTOR"))
}

zscore <- function(vec){
  z <- scale(vec, center = T, scale = T)
  return(z)
}


Zscore.get_mean_and_sd <- function(fullSS="./Data/GWAS/Nalls23andMe_2019/nallsEtAl2019_allSamples_allVariants.mod.txt",
                                   target_col="statistic",
                                   effect_col="beta",
                                   stderr_col="se",
                                   use_saved=T,
                                   output_path="./Data/GWAS/Nalls23andMe_2019/z.info.RDS"){
  if(use_saved & file.exists(output_path)){
    printer("Reading in:",output_path,"...")
    z.info <- readRDS(output_path) 
  } else { 
    printer("Extracting mean and standard deviation from",fullSS,"...")
    if(target_col=="calculate"){
      target_col <- "t_stat"
      sample_x <- data.table::fread(fullSS, nThread = 4, 
                                    select=c(effect_col, stderr_col), 
                                    col.names = c("Effect","StdErr"))
      sample_x <- subset(calculate.tstat(sample_x), select = target_col) 
    }else {
      sample_x <- data.table::fread(fullSS, nThread = 4, select=c(target_col))
    } 
    sample.mean <- mean(sample_x[[1]], na.rm = T)
    sample.stdv <- sd(sample_x[[1]])
    z.info <- list(file.name=fullSS,
                   colname=target_col,
                   sample.mean=sample.mean, 
                   sample.stdv=sample.stdv,
                   sample.min=min(sample_x[[1]]),
                   sample.max=max(sample_x[[1]])) 
    saveRDS(z.info, file = output_path )
  } 
  return(z.info)
}


Zscore <- function(x, z.info){ 
  # Need to use the mean and standard deviation of the FULL dataset (i.e. all beta fomr the full summary stats file)
  sample.stdv <- z.info$sample.stdv
  sample.mean <- z.info$sample.mean 
  z <- (x - sample.mean) / sample.stdv 
  return(z)
}

PAINTOR.create_locusFile <- function(results_path,
                                     PT_results_path,
                                     GWAS_datasets="Nalls23andMe_2019",
                                     QTL_datasets=c("Fairfax_2014_CD14","Fairfax_2014_IFN",
                                                    "Fairfax_2014_LPS2","Fairfax_2014_LPS24"), #NULL
                                     gene,
                                     locus_name){
  printer("++ Creating Locus File...")
  # The locus file should at the very minimum contain the Z-scores for all the populations,
  # though metadata on each SNP such as chromosome, position, and rsid are recommended.
  # The top line of the locus file must contain a header with the names of the fields.
  # The Z-score of a SNP is the Wald statistic, obtained from standard regression of
  # the phenotype onto the SNP.
  
  # Import GWAS data 
  GWAS_path <- file.path(results_path,"Multi-finemap/Multi-finemap_results.txt") 
  z.info.gwas <- Zscore.get_mean_and_sd(fullSS="./Data/GWAS/Nalls23andMe_2019/nallsEtAl2019_allSamples_allVariants.mod.txt",
                                        Effect_col="beta", 
                                        use_saved=T,
                                        output_path="./Data/GWAS/Nalls23andMe_2019/z.info.RDS")
  # Get the party started using the GWAS file
  finemap_DT <- data.table::fread(GWAS_path, nThread = 4) %>%
    dplyr::mutate(CHR=paste0("chr",CHR), 
                  RSID=SNP,
                  ZSCORE.P1=Zscore(x = Effect, z.info = z.info.gwas))
  merged_DT <- finemap_DT 
  
  # Merge QTL data (for loop won't run if QTL_datasets=NULL)
  i=1
  for(qtl in QTL_datasets){
    printer("PAINTOR:: Merging QTL data - ",qtl)
    fulSS.path <- Directory_info(qtl, "fullSS.local") 
    z.info.qtl <- Zscore.get_mean_and_sd(fullSS=fulSS.path,
                                         target_col="beta", 
                                         use_saved=T,
                                         output_path=file.path(dirname(fulSS.path),paste0("z.info.",qtl,".RDS")) )
    # Merge QTL data together
    merged_DT <- mergeQTL.merge_handler(FM_all = merged_DT, qtl_file = qtl)
    merged_DT <- dplyr::mutate(merged_DT, 
                               QTL.ZSCORE = Zscore(x=QTL.Effect, z.info = z.info.qtl ))
    names(merged_DT)[names(merged_DT) == "QTL.ZSCORE"] <- paste0("ZSCORE.P",i+1)
    QTL.cols <- grep("QTL.",colnames(merged_DT), value = T)
    merged_DT <- dplyr::select(merged_DT, -QTL.cols)
    i = i +1
  } 
  # Post-processing
  ## Subset to only necessary columns
  z.cols <- grep("ZSCORE.P",colnames(merged_DT), value=T)
  merged_DT <- subset(merged_DT, select=c("RSID","CHR","POS",z.cols))
  ## Remove NAs (PAINTOR doesn't tolerate missing values)
  merged_DT <- merged_DT[complete.cases(merged_DT),]
  ## Save 
  data.table::fwrite(merged_DT, 
                     file.path(PT_results_path, locus_name),
                     sep=" ", quote = F, na = NA, nThread = 4) 
  return(merged_DT)
}





PAINTOR.create_locusFile.QTL <- function(finemap_DT, 
                                         GWAS_datasets,
                                         QTL_datasets=NULL,
                                         PT_results_path,
                                         locus_name,
                                         metric_suffix=".t_stat",
                                         NA_method=c("fill","drop"),
                                         force_new_zscore=F){
  locus_DT <- finemap_DT
  i=1
  # Iterate GWAS
  for(gwas in GWAS_datasets){
    fullSS <- Directory_info(gwas,"fullSS.local")
    z.info.gwas <- Zscore.get_mean_and_sd(fullSS=fullSS,
                                          target_col = "calculate", 
                                          effect_col = "beta",
                                          stderr_col = "se",
                                          use_saved=(!force_new_zscore),
                                          output_path=file.path(dirname(fullSS),"z.info.RDS"))
    t_stat_name <- paste0(gwas,metric_suffix)
    if(!t_stat_name %in% colnames(locus_DT)){t_stat_name<-"t_stat"}
    locus_DT[,paste0("ZSCORE.P",i)] <- Zscore(x = locus_DT[[t_stat_name]], 
                                              z.info = z.info.gwas)
    i = i+1
  } 
  # Iterate QTL
  for(qtl in QTL_datasets){
    fullSS <- Directory_info(qtl,"fullSS.local")
    z.info.gwas <- Zscore.get_mean_and_sd(fullSS=fullSS,
                                          target_col = "statistic", 
                                          use_saved=(!force_new_zscore),
                                          output_path=file.path(dirname(fullSS),"z.info.RDS")) 
    locus_DT[,paste0("ZSCORE.P",i)] <- Zscore(x = locus_DT[[paste0(qtl,metric_suffix)]], 
                                              z.info = z.info.gwas)
    i = i+1
  }
  # Remove the un-zscored cols
  z.cols <- grep("ZSCORE.P",colnames(locus_DT), value=T)
  locus_DT <- locus_DT %>% dplyr::select(CHR, POS, RSID=SNP, z.cols) %>% 
    dplyr::mutate(CHR=paste0("chr",CHR))  
  ## Remove NAs (PAINTOR doesn't tolerate missing values) 
  if(NA_method[1]=="fill"){
    locus_DT[is.na(locus_DT)] <- 0
  } else if (NA_method[1]=="drop"){
    locus_DT <- locus_DT[complete.cases(locus_DT),]
  }
  ## Save 
  data.table::fwrite(locus_DT, 
                     file.path(PT_results_path, locus_name),
                     sep=" ", quote = F, na = NA, nThread = 4) 
  
  return(data.table::data.table(locus_DT))
}



PAINTOR.prepare_LD.transethnic <- function(subset_DT,
                                           PT_results_path, 
                                           gene,
                                           GWAS_populations=NULL,
                                           QTL_datasets,
                                           QTL_populations=c("AFA","CAU","HIS"), 
                                           LD_reference="1KG_Phase1",
                                           force_new_LD=F, 
                                           shared_snps_only=T,
                                           fillNA = 0){ 
  printer("PAINTOR:: Preparing multi-ethnic LD files...")
  subset_DT$CHR <- paste0("chr",subset_DT$CHR)
  # Download ld matrices
  ld.mat.list <- lapply(1:length(QTL_datasets), function(i){
    qtl <- QTL_datasets[i]
    pop <- QTL_populations[i] 
    printer("+ PAINTOR::",pop)
    LD_matrix <- LD.load_or_create(results_path=file.path(dirname(Directory_info(qtl, "fullSS.local")), gene),
                                   subset_DT=subset_DT,
                                   gene=gene,
                                   download_reference = T,
                                   LD_reference=LD_reference,
                                   superpopulation=translate_population(superpopulation = pop), 
                                   force_new_LD=force_new_LD, 
                                   fillNA = fillNA)   
    return(LD_matrix) 
  })   
  
  # Filter to only SNPs shared between all 3 LD matrices and the input GWAS/QTL data
  if(shared_snps_only){
    shared.snps <- Reduce(intersect, lapply(ld.mat.list,  colnames)) 
    shared.snps <- intersect(shared.snps, subset_DT$SNP)
    printer("+ PAINTOR::",length(shared.snps), "shared SNPs identified.")
  } else {
    shared.snps <- Reduce(intersect, lapply(ld.mat.list,  colnames))
    shared.snps <- intersect(shared.snps,  subset_DT$SNP)
  }
  
  
  # Filter to only SNPs in the dataset 
  .LD_file.paths <- lapply(1:length(ld.mat.list), function(i){ 
    qtl <- QTL_datasets[i]
    pop <- populations[i] 
    .LD_file <- file.path(PT_results_path, paste0(locus_name,".ld",i))
    printer("+ PAINTOR::",pop)
    LD_matrix <- ld.mat.list[[i]]
    ld.mat <- LD_matrix[shared.snps, shared.snps] %>% data.table::data.table()
    printer("++ PAINTOR::",paste(dim(ld.mat),collapse=" x "),"LD matrix.")
    # sum(is.na(ld.mat)); sum(abs(ld.mat)==Inf); sum(is.null(abs(ld.mat)))
    ## Write
    printer("++ PAINTOR:: Writing LD file to ==> ",.LD_file)
    data.table::fwrite(ld.mat,
                       .LD_file,
                       sep = " ", quote = F, 
                       col.names = F, row.names = F, 
                       na = 0.0,
                       nThread = 4)
    
    return(.LD_file)
  }) %>% unlist()
  
  if(!is.null(GWAS_populations)){
    printer("PAINTOR:: Identifying refrence pop for GWAS LD.")
    GWAS.LD.paths <- .LD_file.paths[match(GWAS_populations, 
                                          translate_population(superpopulation = QTL_populations))]
    .LD_file.paths <- append(GWAS.LD.paths, .LD_file.paths)
  } 
  return(list(.LD_file.paths=.LD_file.paths,
              shared.snps=shared.snps))
}


PAINTOR.prepare_LD <- function(results_path,
                               PT_results_path,
                               locus_name,
                               locus_DT, 
                               gene,
                               LD_matrix=NULL){
  ### "VERY IMPORTANT!
  # Particular care must be taken when computing LD from a reference panel such the 1000 genomes.
  ## It is imperative that all the reference and alternate alleles for SNPs from which the Z-scores
  ## were computed match the reference and alternate alleles of the reference panel.
  ## The output of PAINTOR will not be correct if there are mismatches of this type in the data.
  ##Please see wiki section 2a for instructions on how to use the LD utility provided with the software."
  printer("++ PAINTOR:: Creating LD Matrix File...")
  finemap_DT <- data.table::fread(file.path(results_path,"Multi-finemap/Multi-finemap_results.txt"), nThread = 4)
  if(is.null(LD_matrix)){
    LD_matrix <- LD.load_or_create(results_path=results_path,
                                   subset_DT=finemap_DT,
                                   gene=gene,
                                   download_reference = T,
                                   force_new_LD = F,
                                   LD_reference="1KG_Phase1",
                                   superpopulation="EUR") 
  }
  
  ## Make sure SNPs are in the same order as the Locus File
  .LD_file <- file.path(PT_results_path, paste0(locus_name,".ld1"))
  ld.mat <- LD_matrix[locus_DT$RSID, locus_DT$RSID] %>%
    data.table::as.data.table()
  ## Write
  printer("++ PAINTOR:: Writing LD file to ==> ",.LD_file)
  data.table::fwrite(ld.mat,
                     .LD_file,
                     sep = " ", quote = F,
                     col.names = F, row.names = F)
  return(.LD_file)
}






PAINTOR.locusName_handler <- function(locus_name=NULL, 
                                      gene, 
                                      GWAS_datasets=NULL, 
                                      QTL_datasets=NULL){
  if(is.null(locus_name)){ 
    locus_name <- paste(gene,GWAS_datasets,paste(QTL_datasets,collapse = "--"), sep = ".") 
    return(locus_name)
  } 
  return(locus_name)
}

PAINTOR.datatype_handler <- function(GWAS_datasets=NULL, 
                                     QTL_datasets=NULL, 
                                     gene){ 
  if(!all(is.null(GWAS_datasets)) & !all(is.null(QTL_datasets))){ 
    printer("++ PAINTOR::",length(GWAS_dataset_name),"GWAS and",
            length(QTL_datasets),"QTL input datasets detected. Feeding both types into PAINTOR...")
    results_path <- file.path("./Data/QTL",paste0(c(GWAS_datasets,QTL_datasets),collapse="--"), gene)
    
  } else if(!all(is.null(GWAS_datasets)) & all(is.null(QTL_datasets))){
    printer("++ PAINTOR:: Only GWAS input detected.",
            "Feeding",length(GWAS_datasets),"GWAS dataset(s) into PAINTOR...")
    results_path <- file.path("./Data/GWAS",GWAS_datasets, gene)
    
  } else if(all(is.null(GWAS_datasets)) & !all(is.null(QTL_datasets))){
    printer("++ PAINTOR:: Only QTL input detected.",
            "Feeding",length(QTL_datasets),"QTL dataset(s) into PAINTOR...")
    results_path <- file.path("./Data/QTL",paste0(QTL_datasets,collapse="--"), gene)
    
  } else {
    stop("++ PAINTOR:: Neither GWAS nor QTL datasets detected. Please enter at least one valid dataset.")
  } 
  PT_results_path <- file.path(results_path,"PAINTOR")
  dir.create(PT_results_path, showWarnings = F, recursive = T)
  printer("++ PAINTOR:: Results will be stored in  ==> ", PT_results_path) 
  return(PT_results_path)
} 



PAINTOR.download_annotations <- function(PT_results_path, 
                                         locus_name,
                                         locus_DT,
                                         XGR_dataset=NA,
                                         ROADMAP_search=NA,
                                         chromatin_state="TssA",
                                         no_annotations=F){
  printer("++ PAINTOR:: Creating Annotation Matrix File...") 
  .annotations_file <- file.path(PT_results_path, paste0(locus_name,".annotations.txt"))
  if(no_annotations){
    printer("+++ PAINTOR:: no_annotations=T; Prior Probability set to 1 for all SNPs.")
    data.table::fwrite(data.frame(Default_Priors = rep(1,nrow(locus_DT))),
                       .annotations_file,
                       sep=" ", quote = F)
    BED_paths.gz <- NA
  } else{
    ## Download GRs and convert to BED
    if(!is.na(XGR_dataset)){
      printer("+++ PAINTOR:: Gathering annotations via XGR for:",paste(XGR_dataset,collapse=", "))
      GR.annotations <- XGR::xRDataLoader(XGR_dataset)
      printer("+++ PAINTOR:: Writing BED file(s) to  ==> ",file.path(PT_results_path,"annotation_files"))
      BED_paths.gz <- GRs.to.BED(GR.annotations = GR.annotations,
                                 output_path = file.path(PT_results_path,"annotation_files"),
                                 sep=" ")
    }
    if(!is.na(ROADMAP_search)){
      printer("+++ PAINTOR:: Gathering annotations via ROADMAP server for:",paste(ROADMAP_search,collapse=", "))
      # Gather roadmap annotations
      if(ROADMAP_search=="all_cell_types"){ROADMAP_search <- ""}
      # Version 1: Subset locally 
      RoadMap_ref <- GoShifter.search_ROADMAP(fuzzy_search = ROADMAP_search)
      bed_names <- GoShifter.bed_names(RoadMap_ref)
      BED_paths.gz <- GoShifter.get_roadmap_annotations(bed.list = bed_names, 
                                                        chromatin_state = chromatin_state, 
                                                        verbose = T)  
    }
  }
  return(BED_paths.gz)
}


PAINTOR.list_paintor_annotations <- function(annotations_dir=file.path("/sc/orga/projects/pd-omics/brian/PAINTOR_V3.0",
                                                                       "Annotation_directory/Functional_Annotations")){
  types <- basename(list.dirs(annotations_dir, recursive = F))
  printer("PAINTOR:: Available annotations provided by PAINTOR_V3.0...")
  for(atype in types){printer(atype)} 
}


PAINTOR.select_paintor_annotations <- function(annotations_dir=file.path("/sc/orga/projects/pd-omics/brian/PAINTOR_V3.0",
                                                                         "Annotation_directory/Functional_Annotations"),
                                               annotation_types=c("Roadmap_ChromeHMM_15state",
                                                                  "RoadMap_Enhancers",
                                                                  "RoadMap_Promoter",
                                                                  "TFBS"),
                                               locally=T){
  if(locally){
    scp_paths <- file.path(paste0("schilb03@data4.hpc.mssm.edu:", annotations_dir), annotation_types)
    for(path in scp_paths){
      cmd <- paste0("scp -r ",path," ../PAINTOR_annotation/",basename(path))
      printer(cmd)
      # system(cmd)
    }
    BED_paths <- list.files("../PAINTOR_annotation", full.names = T, recursive = T)
  } else { BED_paths <- list.files(file.path(annotations_dir, annotation_types), full.names = T) }
  printer("PAINTOR::",length(BED_paths),"annotations pulled.")
  return(BED_paths)
}



PAINTOR.prepare_annotations <- function(paintor_path="./echolocatoR/tools/PAINTOR_V3.0",
                                        BED_paths,
                                        PT_results_path,
                                        locus_name,
                                        remove_BED=F){
  .locus_file <- file.path(PT_results_path, locus_name)
  .annotations_file <- file.path(PT_results_path, paste0(locus_name,".annotations.txt"))
  
  printer("++ PAINTOR:: Merging and formatting BED files using PAINTOR utility.")
  ## Use PAINTOR utility to format merge BED files
  printer("+++ PAINTOR:: Decompressing BED files.")
  gz.files <- grep("*.gz",BED_paths,value = T)
  if(length(gz.files)>0){
    for (gz in gz.files){
      # Unzip but keep original files
      try({gunzip(gz, overwrite=T, remove=F)})
    }
  }
  
  BED_paths <- gsub("*.gz","", BED_paths) 
  annotation_paths <- file.path(PT_results_path,"annotation_paths.txt")
  # Wait until the BED files are decompressed before trying to write them to a file
  while(any(!file.exists(BED_paths))){
    Sys.sleep(.001)
  }
  printer("+++ PAINTOR:: Writing annotations paths file to  ==> ",annotation_paths)
  data.table::fwrite(list(BED_paths),
                     annotation_paths,
                     sep="\n")
  
  
  cmd <-  paste("python2.7",file.path(paintor_path,"PAINTOR_Utilities","AnnotateLocus.py"),
                "--input", file.path(PT_results_path,"annotation_paths.txt"),
                "--locus", .locus_file,
                "--out", .annotations_file,
                "--chr","CHR",
                "--pos","POS")
  chimera=F 
  if(chimera){cmd <- paste("ml python/2.7.10 &&",cmd)}
  system(cmd)
  printer("+++ PAINTOR:: Annotation--SNP overlap summaries:")
  colSums( data.table::fread(.annotations_file)) 
  if(remove_BED){
    printer("+++ PAINTOR:: Removing temporary decompressed BED files.")
    file.remove(BED_paths)
  } 
}



PAINTOR.survey_annotation <- function(PT_results_path, 
                                      locus_name="Locus1"){
  anno <- data.table::fread(file.path(PT_results_path, paste0(locus_name,".annotations.txt")))
  column_sums <- sort(colSums(anno), decreasing = T)
  signals <- column_sums[column_sums > 0]
  print(signals)
}

PAINTOR.process_results <- function(PT_results_path, 
                                    locus_name="Locus1"){ 
  paintor.results <- data.table::fread(file.path(PT_results_path, paste0(locus_name,".results.txt")), 
                                       nThread = 4) %>% 
    arrange(desc(Posterior_Prob)) 
  return(data.table::data.table(paintor.results))
}

PAINTOR.merge_results <- function(finemap_DT, 
                                  paintor.results,
                                  PP_threshold=.5,
                                  multi_finemap_col_name="PAINTOR"){
  printer("PAINTOR:: Merging PAINTOR results with multi-finemap file.") 
  merged_DT <- data.table:::merge.data.table(finemap_DT, paintor.results[,c("RSID","Posterior_Prob")],
                                             by.x="SNP", by.y="RSID", all.x = T)  
  PP.col.name <- paste0(multi_finemap_col_name,".PP")
  names(merged_DT)[names(merged_DT) == "Posterior_Prob"] <- PP.col.name
  merged_DT[,paste0(multi_finemap_col_name,".Credible_Set")] <- ifelse(subset(merged_DT,select=PP.col.name) > PP_threshold, 1, 0) 
  printer("PAINTOR:: Credible Set size =",sum(subset(merged_DT, select=paste0(multi_finemap_col_name,".Credible_Set")),na.rm=T))
  return(merged_DT)
}


PAINTOR.run <- function(paintor_path,
                        PT_results_path,
                        .LD_file.paths,
                        n_datasets,
                        method="enumerate",#"mcmc"
                        n_causal=2 ){
  printer("+ PAINTOR:: Running PAINTOR...") 
  PT.start <- Sys.time()
  # Enumerate or mcmc sampling
  method_command <- ifelse(method=="enumerate", paste("-enumerate",n_causal),"-mcmc")
  .LD_suffixes <- lapply(.LD_file.paths, function(x){tail(strsplit(x, "[.]")[[1]], 1) }) %>% unlist()
  ## RUN
  # https://github.com/gkichaev/PAINTOR_V3.0/wiki/3.-Running-Software-and-Suggested-Pipeline
  cmd <- paste(
    file.path(paintor_path,"PAINTOR"),
    
    #### REQUIRED ####
    # (required) Filename of the input file containing the
    ## list of the fine-mapping loci [default: N/A]
    "-input",file.path(PT_results_path,"input.files"),
    
    #  (required) The name(s) of the Zscore column
    ## in the header of the locus file (comma separated) [default: N/A]
    "-Zhead", paste(paste0("ZSCORE.P",1:n_datasets), collapse=","),
    
    # (required) Suffix(es) for LD files. Must match the order of
    ## Z-scores in which the -Zhead flag is specified (comma separated) [Default:N/A]
    "-LDname", paste(.LD_suffixes, collapse=","), 
    
    # specify this flag if you want to enumerate all possible configurations
    ## followed by the max number of causal SNPs (eg. -enumerate 3 considers
    ## up to 3 causals at each locus) [Default: not specified]
    method_command, 
    #  should the algorithm be run with MCMC? [Default: not specified]
    # "-mcmc",
    
    #	specify the number of causals to pre-compute enrichments with [default: 2]
    "-max_causal",n_causal,
    
    #### OPTIONAL ####
    # The names of the annotations to include in model (comma separated)
    ## [default: N/A]
    # "-annotations",paste(basename(bed),collapse=","),
    
    # Input directory with all run files [default: ./ ]
    "-in", paste0(PT_results_path),
    
    #  Output directory where output will be written [default: ./ ]
    "-out",file.path(PT_results_path),
    
    # Output Filename for enrichment estimates [default: Enrichment.Estimate]
    "-Gname","Enrichment.Estimates.txt",
    
    # Suffix for ouput files of results [Default: results]
    "-RESname","results.txt",
    
    # Suffix for annotation files [Default: annotations]
    "-ANname","annotations.txt",
    
    "-set_seed","2019"
  )
  printer(cmd)
  system(cmd)
  # }
  # EXAMPLE
  # cmd <- paste("cd",paintor_path,"&& ./PAINTOR -input SampleData/input.files -in SampleData/ -out SampleData/ -Zhead Zscore -LDname ld -enumerate 2 -annotations DHS")
  # system(cmd) 
  PT.end <- Sys.time()
  printer("PAINTOR:: Completed fine-mapping in",round((PT.end-PT.start)/60, 2),"minutes.")
  
  # Check the time it took to see if it didn't actually run.
  ## Re-enter command until it does.
  if(PT.end-PT.start<1){
    PAINTOR.run(paintor_path,
                PT_results_path,
                .LD_file.paths,
                n_datasets,
                method,
                n_causal)
  }
}

# ------------------------------#
# ----------- PAINTOR ----------#
# ------------------------------#


PAINTOR.import_QTL_DT <- function(QTL_datasets, 
                                  gene, 
                                  trim_gene_limits=F, 
                                  force_new_subset=F,
                                  metric="t_stat"){ 
  QTL.dat <- lapply(QTL_datasets, function(qtl){
    printer("++ PAINTOR:: Importing",qtl,"...") 
    fullSS_path <- Directory_info(qtl)
    qtl.dir <- basename(dirname(dirname(fullSS_path)))
    qtl.subdir <- basename(dirname(fullSS_path))
    subset_path <- file.path(dirname(fullSS_path),paste0(gene,"_",qtl,".txt"))
    
    # Save subset
    # force_new_subset=T
    if(file.exists(subset_path) & force_new_subset==F){
      subset_DT <- data.table::fread(subset_path, nThread = 4)
    } else {
      printer("PAINTOR:: Creating subset file for",gene)
      qtl.dat <- data.table::fread(fullSS_path, nThread = 4) 
      ## Remove the "gene" column bc it confuses subsetting functions
      # qtl.dat <- dplyr::select(qtl.dat, select = -gene) %>% 
      #   subset(gene_name==gene)  
      data.table::fwrite(qtl.dat, subset_path, sep="\t")
      subset_DT <- preprocess_subset(gene = gene,  
                                     subset_path = subset_path,
                                     chrom_col = "chr",
                                     position_col = "pos_snps",
                                     snp_col = "snps", 
                                     effect_col = "beta", 
                                     pval_col = "pvalue", 
                                     tstat_col = "statistic",
                                     stderr_col = "calculate",
                                     A1_col = "ref", 
                                     A2_col = "alt")  %>% data.table::data.table()   
    }  
    subset_DT <- cbind(Dataset=qtl, subset_DT) 
    return(subset_DT)
  }) %>% data.table::rbindlist()
  # Trim by coordinates
  if(trim_gene_limits){
    QTL.dat <- gene_trimmer(subset_DT = QTL.dat, gene = gene)
  }  
  # Spread data
  lead.snps <- (QTL.dat %>% dplyr::group_by(Dataset) %>% top_n(n=1, wt=-P))
  QTL.dat$Dataset <- paste(QTL.dat$Dataset,metric,sep=".")
  QTL.spread <-  tidyr::spread(data = data.frame(QTL.dat)[,c("Dataset","CHR","POS","SNP","leadSNP","A1","A2",metric)], 
                               key="Dataset", value = metric) %>% 
    data.table::data.table()
  
  # Find lead SNP for each QTL condition 
  
  for (d in unique(lead.snps$Dataset)){
    lead <- subset(lead.snps, Dataset==d)$SNP
    QTL.spread[,paste0(d,".leadSNP")] <- ifelse(QTL.spread$SNP==lead, T, F)
  }  
  return(QTL.spread)
}



PAINTOR <- function(finemap_DT=NULL,
                    paintor_path = "./echolocatoR/tools/PAINTOR_V3.0",
                    GWAS_datasets=NULL,
                    QTL_datasets=NULL,
                    gene,
                    locus_name=NULL,
                    n_causal=5,
                    XGR_dataset=NA,#"FANTOM5_Enhancer",
                    ROADMAP_search=NA,
                    chromatin_state="TssA",
                    use_annotations=F,
                    PP_threshold=.95,
                    multi_finemap_col_name="PAINTOR",
                    trim_gene_limits=F, 
                    force_new_subset=F,
                    GWAS_populations="EUR",
                    QTL_populations="EUR",
                    LD_reference="1KG_Phase1",
                    force_new_LD=F,
                    LD_matrix=NULL){
  #@@@@@@@@@ Multi-condition, Multi-GWAS, multi-QTL, Multi-ethnic 
  #    and/or Annotated Fine-mapping @@@@@@@@@
  # Note: All file formats are assumed to be single space delimited.
  
  ## Quick setup
  # chromatin_state="TssA"; paintor_path = "./echolocatoR/tools/PAINTOR_V3.0";  QTL_datasets = c("Fairfax_2014_IFN","Fairfax_2014_LPS24");locus_name=NA; gene="LRRK2";  no_annotations=F;  XGR_dataset=NA; ROADMAP_search="monocyte"; n_causal=5; finemap_DT=NA; multi_finemap_col_name="PAINTOR_Fairfax";  PP_threshold=.95; QTL_datasets=c("MESA_CAU","MESA_HIS"); trim_gene_limits=T; multi_finemap_col_name <- "PAINTOR_MESA_transethnic"; QTL_populations = c("CAU","HIS"); GWAS_populations="EUR"; use_annotations=F; force_new_LD=F; force_new_subset=F; LD_reference="1KG_Phase1"; GWAS_datasets = "Nalls23andMe_2019";
  
  # if(is.na(GWAS_dataset_name)){
  #   GWAS_dataset_name <- basename(dirname(results_path)); 
  #   printer("PAINTOR:: Inferring GWAS dataset name:", GWAS_dataset_name)
  #   } 
  
  locus_name <- PAINTOR.locusName_handler(locus_name, 
                                          gene, 
                                          GWAS_datasets, 
                                          QTL_datasets)
  printer("++ locus_name =",locus_name)
  PT_results_path <- PAINTOR.datatype_handler(GWAS_datasets=GWAS_datasets, 
                                              QTL_datasets=QTL_datasets, 
                                              gene=gene)
  printer("****** Double checking PT_results_path",PT_results_path)
  results_path <- dirname(PT_results_path)
  
  if(!is.null(GWAS_datasets)){ 
      fullSS <- Directory_info(GWAS_datasets, "fullSS.local")
      results_path <- file.path(dirname(fullSS),gene)
      mfm_path <- file.path(results_path,"Multi-finemap/Multi-finemap_results.txt")
      if(is.null(finemap_DT)){
        printer("PAINTOR:: No finemap_DT supplied. Retrieving from storage:",mfm_path)
        finemap_DT <- data.table::fread(mfm_path, nThread = 4) 
      } 
      finemap_DT <- calculate.tstat(finemap_DT=finemap_DT)
  } 
  if(!is.null(QTL_datasets)){
    qtl_DT <- PAINTOR.import_QTL_DT(QTL_datasets = QTL_datasets,
                                        gene = gene,
                                        metric="t_stat",
                                        trim_gene_limits = trim_gene_limits, 
                                        force_new_subset = force_new_subset)
    if(!is.null(GWAS_datasets) & length(GWAS_datasets)==1){
      finemap_DT <- finemap_DT[,c("SNP","t_stat","A1","A2","P","leadSNP")]
      colnames(finemap_DT)[-1] <- paste0("GWAS.",colnames(finemap_DT)[-1])
      finemap_DT <- data.table:::merge.data.table(finemap_DT, qtl_DT, by="SNP")
      # Flip alleles
      finemap_DT <- finemap_DT %>% 
        dplyr::mutate(GWAS.t_stat = ifelse((GWAS.A1!=A1 & GWAS.A2!=A2),
                                           -GWAS.t_stat,GWAS.t_stat)) %>%
        data.table::data.table()
      # Rename GWAS cols
      colnames(finemap_DT) <- gsub("GWAS.",paste0(GWAS_datasets[1],"."),colnames(finemap_DT))
    } else{
      finemap_DT <- qtl_DT
    } 
  } 
  
 
  locus_DT <- PAINTOR.create_locusFile.QTL(finemap_DT=finemap_DT,
                                           GWAS_datasets=GWAS_datasets,
                                           QTL_datasets=QTL_datasets, 
                                           locus_name = locus_name,
                                           PT_results_path = PT_results_path,
                                           metric_suffix = ".t_stat",
                                           force_new_zscore = F,
                                           NA_method = "drop") 
  n_datasets <- sum(grepl(colnames(locus_DT), pattern = "ZSCORE")) 
  
  # 2. LD Matrix File
  if(length(populations)==1 & !is.null(GWAS_datasets)){
    .LD_file.paths <- PAINTOR.prepare_LD(results_path = results_path,
                                         PT_results_path = PT_results_path,
                                         locus_name = locus_name,
                                         locus_DT = locus_DT,
                                         LD_matrix=LD_matrix)
  } else{
    LD.list <- PAINTOR.prepare_LD.transethnic(subset_DT=subset(finemap_DT, SNP %in% locus_DT$RSID),
                                              PT_results_path=PT_results_path, 
                                              gene=gene,
                                              QTL_datasets=QTL_datasets,
                                              GWAS_populations=GWAS_populations,
                                              QTL_populations=QTL_populations, 
                                              LD_reference=LD_reference,
                                              force_new_LD=force_new_LD, 
                                              fillNA = NA)
    .LD_file.paths <-  LD.list$.LD_file.paths
    # .LD_file.paths[2] <- gsub(".ld2",".ld1",.LD_file.paths[2])
    # Make sure SNP order matches that of LD matrices  
    locus_DT <- subset(locus_DT, RSID %in% shared.snps)  
    locus_DT <- locus_DT[match(shared.snps, locus_DT$RSID),]
    locus_DT  <- locus_DT[complete.cases(locus_DT),]
    
    ## Save 
    data.table::fwrite(locus_DT, 
                       file.path(PT_results_path, locus_name),
                       sep=" ", quote = F, na = NA, nThread = 4)
  }
  
  
  # 3. Annotation Matrix File 
  if(use_annotations){
    PAINTOR.list_paintor_annotations(annotations_dir="../PAINTOR_annotation/") 
    BED_paths <- PAINTOR.select_paintor_annotations(annotations_dir=file.path("/sc/orga/projects/pd-omics/brian/PAINTOR_V3.0",
                                                                              "Annotation_directory/Functional_Annotations"),
                                                    annotation_types=c("Roadmap_ChromeHMM_15state",
                                                                       "RoadMap_Enhancers",
                                                                       "RoadMap_Promoter",
                                                                       "TFBS"), locally = T)
    PAINTOR.prepare_annotations(paintor_path = paintor_path,
                                BED_paths = BED_paths, 
                                PT_results_path = PT_results_path, 
                                locus_name = locus_name)
  } else {
    # Run with all priors set to 1
    BED_paths <- PAINTOR.download_annotations(PT_results_path,
                                              locus_name = locus_name,
                                              locus_DT = locus_DT, 
                                              no_annotations=T)
  }
  
  # 4. Input File
  printer("+ PAINTOR:: Preparing input.files")
  data.table::fwrite(list(locus_name), 
                     file.path(PT_results_path, "input.files"), 
                     sep="\n")
  
  
  # 5. Run PAINTOR!
  
  PAINTOR.run(paintor_path = paintor_path, 
              PT_results_path = PT_results_path, 
              n_datasets = n_datasets, 
              #enumerate is actually faster when n_causal is small (<3)
              # but far larger n_causal use mcmc
              method = "mcmc", 
              n_causal = 2, 
              .LD_file.paths = .LD_file.paths)
  
  
  
  # Summarise the results
  # locus_name="LRRK2.Nalls23andMe_2019."; multi_finemap_col_name="PAINTOR"
  # locus_name="LRRK2.Nalls23andMe_2019.Fairfax_2014_CD14--Fairfax_2014_IFN--Fairfax_2014_LPS2--Fairfax_2014_LPS24"; multi_finemap_col_name="PAINTOR_Fairfax"
  
  paintor.results <- PAINTOR.process_results(PT_results_path=PT_results_path, 
                                             locus_name=locus_name)
  # Remove columns with the same name
  ff.cols <- grep(paste0(multi_finemap_col_name,"."), colnames(finemap_DT), value = T)
  if(length(ff.cols)>0){finemap_DT[,c(ff.cols):=NULL]} 
  # Merge results 
  # LD_matrix <- readRDS("./Data/QTL/MESA/CAU/LRRK2/plink/LD_matrix.RData")
  finemap_DT.P <- PAINTOR.import_QTL_DT(QTL_datasets = QTL_datasets,
                                      gene = gene,
                                      trim_gene_limits = trim_gene_limits,
                                      force_new_subset = force_new_subset,
                                      metric = "P")
  finemap_DT.P <- data.table:::merge.data.table(finemap_DT.P, 
                                subset(finemap_DT, select=c("SNP", 
                                                            paste0(GWAS_datasets,".P"),
                                                            paste0(GWAS_datasets,".leadSNP")))
                                )
  
  merged_DT <- PAINTOR.merge_results(finemap_DT = finemap_DT.P,
                                     paintor.results = paintor.results,
                                     PP_threshold = PP_threshold,
                                     multi_finemap_col_name = multi_finemap_col_name)   
  merged_DT <- find_consensus_SNPs(merged_DT, support_thresh = 1)
   
  # Update Consensus SNP col and Summarise 
  # merged_DT <- find_consensus_SNPs(merged_DT, support_thresh = 2)
  # top_snps <- (finemap_DT %>% arrange(desc(Support)))[,c("SNP","Support","Consensus_SNP")] %>% head(10)
  # data.table::fwrite(merged_DT, mfm_path, nThread = 4, sep="\t")
  
  # PLOT 
  transethnic_plot(merged_DT = merged_DT,
                   save_path = file.path(PT_results_path,"track_plot.enumerate2.png"),
                   PAINTOR.label="PAINTOR\nTrans-ethnic",
                   conditions = c("Nalls23andMe_2019","MESA_CAU","MESA_HIS"))
  
  return(merged_DT)
}





### ### ### PLOT ### ### ### 

gather.transethnic.LD <- function(merged_DT, 
                                  gene,
                                  conditions=c("Nalls23andMe_2019", 
                                               "MESA_AFA","MESA_CAU","MESA_HIS")){
  for(cond in conditions){
    printer("+ PAINTOR:: Retrieving LD with lead SNP for",cond)
    fullSS <- dirname(Directory_info(cond, variable = "fullSS.local"))
    ld.path <- file.path(fullSS, gene,"plink","LD_matrix.RData")
    LD_matrix <- readRDS(ld.path) 
    lead.snp <- merged_DT[merged_DT[,paste0(cond,".leadSNP")]==T,]$SNP
    dat <- data.table(SNP=names(LD_matrix[lead.snp,]),
                      r2=LD_matrix[lead.snp,]^2)
    merged_DT <- data.table:::merge.data.table(merged_DT, dat, by = "SNP")
    colnames(merged_DT)[colnames(merged_DT)=="r2"] <- paste0(cond,".r2")
  }
  merged_DT <- dplyr::mutate(merged_DT,  
                             Mb=round(POS/1000000,3)) 
  
  return(data.table::data.table(merged_DT))
}



transethnic_plot <- function(merged_DT,
                             save_path,
                             title=gene,
                             subtitle="Trans-ethnic Fine-mapping",
                             PAINTOR.label="PAINTOR\nTrans-ethnic",
                             conditions=c("MESA_AFA","MESA_CAU","MESA_HIS")){ 
  # cons.snp <-  "rs7294619"; subset(plot_DT, SNP==cons.snp);
  plot_DT <- gather.transethnic.LD(merged_DT, conditions=conditions)
  plot_DT$PAINTOR.label <- PAINTOR.label
  # Melt P
  P.vars <- paste0(conditions,".P")
  r2.vars <- paste0(conditions,".r2") 
  leadSNP.vars <-  paste0(conditions,".leadSNP") 
  id.vars <- grep(paste(c(P.vars, r2.vars), collapse="|"), colnames(plot_DT), value = T, invert = T)
  dat <- data.table::melt.data.table(plot_DT, 
                                     id.vars = id.vars, 
                                     measure.vars = list(P.vars, r2.vars, leadSNP.vars), 
                                     variable.name = c("Condition_number"),
                                     value.name = c("P","r2","lead.snp")) 
  cond.dict <- setNames(conditions, 1:length(conditions))
  dat$Condition <- cond.dict[dat$Condition_number]
  
  
  
  library(patchwork)
  gg <- ggplot(data=dat, aes(x=Mb, y=-log10(P), color=r2)) + 
    geom_point() + 
    facet_grid(Condition~., scales="free_y") + 
    theme_classic() + 
    labs(y="-log10(P-value)", color=paste0("r2 with\n", subset(dat, lead.snp==T)$SNP)) + 
    # Lead SNP
    geom_point(dat=subset(dat, lead.snp==T), shape=1, size=6, color="red") +
    ggrepel::geom_label_repel(dat=subset(dat, lead.snp==T) , aes(label=SNP), 
                              alpha=0.8, point.padding = 1) + 
    # Credible Set
    geom_point(dat=subset(dat, Support>0), shape=1, size=6, color="green") +  
    #
    scale_colour_gradient(low = "blue", high = "red", limits=c(0,1), breaks=c(0,.5,1)) +
    theme(strip.background = element_rect(fill = "grey10"),
          strip.text.y = element_text(color = "white", angle = 0)) +  
    
    # Trans-ethnic PP layer 
    ggplot(data=plot_DT, aes(x=Mb, y=PAINTOR_MESA_transethnic.PP, 
                             color=PAINTOR_MESA_transethnic.PP)) + 
    geom_point() + 
    scale_color_viridis_c(limits=c(0,1), breaks=c(0,.5,1)) +
    # Credible Set
    geom_point(dat=subset(dat, Support>0), shape=1, size=6, color="green") +
    # Trans-ethnic fine-mapped CS SNPs
    geom_point(dat=subset(dat, Support>0), 
               shape=1, size=6, color="green") +
    ggrepel::geom_label_repel(dat=subset(dat, PAINTOR_MESA_transethnic.Credible_Set>0,
                                         select=c("SNP","CHR","POS","Mb","leadSNP","PAINTOR_MESA_transethnic.PP")) %>% unique(), aes(label=SNP), 
                              alpha=0.8, point.padding = 1, color="green") +  
    theme_classic() + 
    facet_grid(PAINTOR.label~.) + 
    labs(y="PP", color=paste0("PP")) +
    theme(strip.background = element_rect(fill = "grey10"),
          strip.text.y = element_text(color = "white", angle = 0)) +   
    # Overall layers
    patchwork::plot_layout(ncol = 1, heights = c(1,.5)) + 
    patchwork::plot_annotation(title = title, 
                               subtitle = subtitle, 
                               theme =  theme(plot.title = element_text(hjust = 0.5),
                                              plot.subtitle = element_text(hjust = 0.5)))
  # Display
  print(gg)  
  
  if(!is.null(save_path)){
    printer("PAINTOR:: Saving plot to =>",save_path)
    ggsave(save_path, plot = gg, 
           dpi = 400, height = 9, width = 9)
  }
}






