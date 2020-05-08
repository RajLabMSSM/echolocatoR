
  ##------------------------------------------------------##
 #########  ------------ GoShifter  ------------ #########
##----------------------------------------------------##

 # GitHub Repo:
 # https://github.com/immunogenomics/goshifter
 
 #




GoShifter.create_snpmap <- function(snp_df, GS_results, verbose=T){
  printer("++ GoShifter: Creating snpmap file...", v=verbose)  
  snp_df %>%
    dplyr::rename(Chrom = CHR, BP = POS) %>% 
    dplyr::mutate(Chrom = paste0("chr",Chrom)) %>%
    dplyr::select(SNP, Chrom, BP) %>%
    data.table::fwrite(file.path(GS_results,"snpmap.txt"),sep="\t") 
  printer("+++ snpmap written to :",file.path(GS_results,"snpmap.txt"))
}



GoShifter.create_LD <- function(results_path, verbose=T){
  printer("++ GoShifter: Creating LD file(s)...", v=verbose)
  # (chrA\tposA\trsIdA\tposB\trsIdB\tRsquared\tDPrime)
  GS_results_path <- file.path(results_path,"GoShifter")
  plink.ld <- data.table::fread(file.path(results_path,"plink","saved","plink.ld"))
  ld_file <- plink.ld %>% dplyr::rename(chrA = CHR_A, 
                                        posA = BP_A,
                                        rsIdA = SNP_A,
                                        posB = BP_B,
                                        rsIdB = SNP_B, 
                                        Rsquared = R,
                                        DPrime = DP) %>% 
    dplyr::mutate(Rsquared = Rsquared^2, chrA = paste0("chr",chrA)) %>% 
    dplyr::select(chrA,posA,rsIdA,posB,rsIdB,Rsquared,DPrime)
  # Create tabix file(s) 
  LD_folder <- file.path(GS_results_path,"LD")
  dir.create(LD_folder, showWarnings = F, recursive = T)
  for(chr in unique(ld_file$chrA)){
    ld_path <- file.path(LD_folder, paste0(chr,".EUR.tsv"))
    gz_path <- paste0(ld_path,".gz") 
    file.remove(gz_path) 
    data.table::fwrite(subset(ld_file, chrA==chr), ld_path,
                       sep="\t", quote = F, col.names = F) 
    # gzip(ld_path, destname = gz_path)
    system(paste("bgzip",ld_path))
    cmd <- paste("tabix",
                 "--begin 2",
                 "--end 2",
                 "--force",
                 gz_path)
    system(cmd) 
  }
  printer("+++ LD file(s) written to :",LD_folder)
}


GoShifter.search_ROADMAP <- function(Roadmap_reference = "./echolocatoR/tools/Annotations/ROADMAP/ROADMAP_Epigenomic.js",
                                   EID_filter = NA,
                                   GROUP_filter = NA,
                                   ANATOMY_filter = NA,
                                   GR_filter = NA,
                                   fuzzy_search = NA){ 
  printer("++ GoShifter: Searching for Roadmap annotation BED files...", v=verbose)
  fuzzy_search <- paste(fuzzy_search, collapse="|")
  RM_ref <- suppressWarnings(data.table::fread(Roadmap_reference, skip = 1 , header = F, 
                                               col.names = c("EID",
                                                             "GROUP",
                                                             "ANATOMY",
                                                             "GR",
                                                             "Epigenome.Mnemonic",
                                                             "Standardized.Epigenome.name",
                                                             "Epigenome.name",
                                                             "TYPE")) ) 
  if(!is.na(EID_filter)){RM_ref <- subset(RM_ref, EID %in% EID_filter)}
  if(!is.na(GROUP_filter)){RM_ref <- subset(RM_ref, GROUP %in% GROUP_filter)}
  if(!is.na(ANATOMY_filter)){RM_ref <- subset(RM_ref, ANATOMY %in% ANATOMY_filter)}
  if(!is.na(GR_filter)){RM_ref <- subset(RM_ref, GR %in% GR_filter)}
  if(!is.na(fuzzy_search)){
    RM_ref <- RM_ref[(tolower(GROUP) %like% fuzzy_search) | (tolower(ANATOMY)  %like% fuzzy_search) | (tolower(GR) %like% fuzzy_search) | (tolower(Standardized.Epigenome.name) %like% fuzzy_search) | (tolower(Epigenome.name) %like% fuzzy_search) | (tolower(TYPE) %like% fuzzy_search )]
  } 
  total_found <- dim(RM_ref)[1]
  if(total_found > 0){
    printer("++ GoShifter: Found",dim(RM_ref)[1],"annotation BED files matching search query.", v=verbose)
  } else {stop("No annotation BED files found :(")} 
  return(RM_ref)
}


GoShifter.bed_names <- function(RM_ref, suffix="_15_coreMarks_mnemonics.bed.gz"){ 
  # All files found here:
  # https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/
  
  # Alternative suffixes 
  # _15_coreMarks_dense.bed.gz
  # _15_coreMarks_dense.bed.gz
  
  # Create URLs
  bed_names <- lapply(unique(RM_ref$EID), function(eid, 
                                                   v=verbose, 
                                                   suffix.=suffix){ 
    return(paste0(eid,suffix)) 
  }) %>% unlist() 
  return(bed_names)
}
 
GoShifter.list_chromatin_states <- function(annotations_path = "./echolocatoR/tools/Annotations"){
  # Term key for chromatin states
  chromState_key <- data.table::fread(file.path(annotations_path, "ROADMAP", 
                                                "ROADMAP_chromatinState_HMM.tsv"))
  return(chromState_key)
}



GoShifter.get_roadmap_annotations <- function(annotations_path = "./echolocatoR/tools/Annotations", 
                                              bed.list,
                                              chromatin_state = "TssA",
                                              verbose = T){
  output_paths <- lapply(bed.list, function(bed, v=verbose){
    eid <- strsplit(bed, "_")[[1]][1]
    roadmap.annot.path <- file.path(annotations_path,"ROADMAP/Chromatin_Marks")
    dir.create(roadmap.annot.path, showWarnings = F, recursive = T)
    output_path <- file.path(roadmap.annot.path, bed)
    bed_subset <- file.path(dirname(output_path), paste0(eid,"_",chromatin_state,"_subset.bed"))
    bed_subset.gz <- paste0(bed_subset,".gz") 
    bed_url <- file.path("https://egg2.wustl.edu/roadmap/data/byFileType",
                         "chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final",
                         bed)
   
    # If you have the exact file subset you need, just load it
    if(file.exists(bed_subset.gz)){
      printer("+++ GoShifter:: Importing Previously downloaded BED subset:",bed_subset.gz, v=v);
      dat <- data.table::fread(bed_subset.gz, col.names = c("Chrom", "Start","End","State","NA") ) 
    } else {
      # If you have the right file but it needs to be subset
      if (file.exists(output_path)){
        printer("+++ GoShifter:: Importing Previously downloaded full BED:",bed_subset.gz, v=v);
        dat <- data.table::fread(output_path, col.names = c("Chrom", "Start","End","State") )
      } else {
        printer("+++ GoShifter:: Downloading annotation BED file from Roadmap server...", v=v)
        dat <- data.table::fread(bed_url, col.names = c("Chrom", "Start","End","State") )  
        data.table::fwrite(dat, file.path(roadmap.annot.path, bed), col.names = F, sep = " ")
      }
      # Subset data to just the chromatin state(s) you want to check enrichment for 
      dat[, c("Num", "State") := tstrsplit(State, "_", fixed=TRUE)] # inplace transform
      if(!is.na(chromatin_state)){
        printer("+++ GoShifter:: Subsetting BED file by chromatin state: ",chromatin_state)
        dat <- dat[State == chromatin_state,]
        data.table::fwrite(dat, bed_subset, col.names = F,sep = "\t") 
        gzip(bed_subset, destname=bed_subset.gz, overwrite = T)
      } else {
        printer("+++ GoShifter:: Including all",length(unique(dat$State)),"chromatin states.")
        bed_subset.gz <- file.path(roadmap.annot.path, bed)
      }
      printer("+++ GoShifter:: BED file saved to  ==> ", bed_subset.gz)
      
    } 
    return(bed_subset.gz)
  }) %>% unlist() 
  return(output_paths)
}

GoShifter.process_results <- function(RoadMap_ref, 
                                      results_path, 
                                      output_bed,
                                      out.prefix){
  # From GitHub Repo: https://github.com/immunogenomics/goshifter
  # ------- .LOCUSSCORE -------
  ## *.locusscore – the likelihood of a locus to overlap an annotation under the null. 
  ##The smaller the value the more likely a locus overlaps an annotation not by chance.
  ##Loci not overlapping any annotation are denoted as “N/A”.
  ## ** NOTE **: This file is terribly named. Each SNP is not necessarily a "locus".
  .locusscore<- data.table::fread(file.path(results_path,"GoShifter",
                                            paste0(out.prefix,".nperm1000.locusscore")),sep="\t")
  .locusscore$Annotation <- basename(output_bed)
  eid <- strsplit(basename(output_bed),"_")[[1]][1]
  .locusscore$EID <- eid 
  data <- subset(RoadMap_ref, EID == eid)
  GS.stats <- data.table:::merge.data.table(.locusscore, data.table::data.table(data), by="EID")
  
  # ------- .ENRICH -------
  ## *.enrich – output file with observed and permuted overlap values
  # nperm = 0 is the observed overlap nSnpOverlap – 
  ## number of loci where at least one SNP overlaps an annotation allSnps 
  ## – total number of tested loci enrichment – nSnpOverlap/allSnps
  # ***Note***, P-value is the number of times the “enrichment” is greater or equal to the observed overlap divided by total number of permutations. 
  .enrich <- data.table::fread(file.path(results_path,"GoShifter",
                                                    paste0(out.prefix,".nperm1000.enrich")) ,sep="\t") 
  # Manually calculate p-value...annoying that they don't just provide this in output files
  GS.stats$pval <- sum(.enrich$enrichment >= .enrich$nSnpOverlap) / max(.enrich$nperm)
  
  mean.enrich <- .enrich %>% summarise(mean.nSnpOverlap = mean(nSnpOverlap), 
                          allSnps=mean(allSnps), 
                          mean.enrichment = mean(enrichment))
  GS.stats$mean.nSnpOverlap <- mean.enrich$mean.nSnpOverlap
  GS.stats$allSnps <- mean.enrich$allSnps
  GS.stats$mean.enrichment <- mean.enrich$mean.enrichment
  GS.stats <- dt.replace(GS.stats, "N/A", NA)
  return(GS.stats)
}
 

GoShifter.check_overlap <- function(output_bed, 
                                    GS_results_path, 
                                    overlap_threshold = 1){
  bed.file <- data.table::fread(output_bed, col.names = c("Chrom", "Start","End","State","NA")) 
  snpmap <- data.table::fread(file.path(GS_results_path,"snpmap.txt")) 
  # Find overlappying bed
  bed.overlap <- lapply(1:nrow(snpmap), function(i){
    row <- snpmap[i,]
    bed.sub <- subset(bed.file, (Chrom==row$Chrom) & (Start<=row$BP & End>=row$BP))
    return(bed.sub)
  }) %>% data.table::rbindlist(fill=T)  
  return(bed.overlap)
}

GoShifter.run <- function(RoadMap_ref,
                           results_path,
                           permutations = 1000, 
                           goshifter_path = "./echolocatoR/tools/goshifter",
                           chromatin_state = "TssA",
                           R2_filter = 0.8,
                           remove_tmps = T,
                           verbose = T,
                           overlap_threshold = 1){
  bed_names <- GoShifter.bed_names(RoadMap_ref) 
  GS_results_path <- file.path(results_path,"GoShifter")
  # Iterate over each bed file
  GS_results <- lapply(bed_names, function(bed,
                                           remove_tmps. = remove_tmps,
                                           results_path. = results_path,
                                           RoadMap_ref. = RoadMap_ref,
                                           R2_filter. = R2_filter,
                                           overlap_threshold. = overlap_threshold){ 
      output_bed <- GoShifter.get_roadmap_annotations(bed.list = bed, 
                                                      verbose = verbose, 
                                                      chromatin_state = chromatin_state)
      # Check if there's any overlap first
      bed.overlap <- GoShifter.check_overlap(output_bed, 
                                             GS_results_path,
                                             overlap_threshold.)
      if(nrow(bed.overlap) >= overlap_threshold){
        printer("GoShifter:: SNPs overlapping between snpmap and bed files =",nrow(bed.overlap) )
        printer("GoShifter: Testing for enrichment of SNPs within annotation:",
                basename(output_bed), v=verbose) 
        out.prefix <- gsub(".bed.gz","",basename(output_bed)) 
        cmd <- paste(
          file.path(goshifter_path,"goshifter.py"),
          "--snpmap",file.path(GS_results_path,"snpmap.txt"),
          "--annotation",output_bed,
          "--permute",permutations,
          "--ld",file.path(GS_results_path,"LD"),
          "--out",file.path(GS_results_path, out.prefix),
          "--rsquared",R2_filter. # Optional
          # "--window",window, # Optional
          # "--min-shift",min_shift, # Optional
          # "--max-shift",max_shift, # Optional
          # "--ld-extend",ld_extend, # Optional
          # "--nold" # Optional
        ) 
        cmd.out <- system(cmd, intern = T) 
        cat(paste(cmd.out,collapse="\n"))
        # Gather results from written output files
        GS.stats <- GoShifter.process_results(RoadMap_ref = RoadMap_ref.,
                                              results_path = results_path., 
                                              output_bed = output_bed,
                                              out.prefix = out.prefix)
        GS.stats$chromatin_state <- chromatin_state
        # Extract p-value from console output 
        GS.pval <- gsub("p-value = ","",cmd.out[startsWith(cmd.out, "p-value")]) %>% as.numeric()
        GS.stats$GS.pval <- GS.pval
        GS.stats$true.overlap <- nrow(bed.overlap)
        if(remove_tmps.){ file.remove(output_bed)}
      } else{
        printer("GoShifter:: There is no overlap between the snpmap and bed files ( at overlap_threshold =",
                overlap_threshold,")." ) 
        GS.stats <- data.table::data.table(Ennotation=output_bed, 
                                           pval=NA,
                                           GS.pval=NA, 
                                           true.overlap=0)
      }
    return(GS.stats)
  }) %>% data.table::rbindlist(fill=T) 
  printer(".")
  printer(".")
  printer(".")
  return(GS_results)
}



GoShifter.summarise_results <- function(GS_results){
  # Summarise results
  sig_results <- GS_results %>%  
    # Apply bonferonni correction
    # subset(enrichment <= 0.05/nrow(GS_results) ) %>%
    subset(enrichment <= 0.05) %>%
    arrange(desc(enrichment), desc(nSnpOverlap)) %>%
    dplyr::mutate(enrichment = formatC(enrichment,format="e", digits=7) )  
  createDT(sig_results) %>% print()
  
  # Get sig enrichments per tissue per chomatin_state
  sig_count <- sig_results %>% 
    dplyr::group_by(EID) %>% 
    slice(1) %>%
    dplyr::group_by(chromatin_state) %>% 
    count() 
  printer("+++",sig_count$n,"/",nrow(RoadMap_ref),
          "tissues had significant enrichment for the chomatin state(s):",
          paste(chromatin_state, collapse = ", "), v = verbose)
  return(sig_results)
}



####----- Go Shifter: Main Function -----####

GoShifter <- function(results_path,
                      snp_df,
                      SNP.Group = "",
                      goshifter_path = "./echolocatoR/tools/goshifter", 
                      permutations = 1000,
                      ROADMAP_search = "",
                      ROADMAP_type = NA, # PrimaryTissue
                      chromatin_states = c("TssA"),
                      R2_filter = 0.8, # Include LD SNPs at rsquared >= NUM [default: 0.8]
                      overlap_threshold = 1,
                      force_new_goshifter = F,
                      remove_tmps = T,
                      verbose = T,
                      save_results = T){ 
  # results_path <-  "./Data/GWAS/Nalls23andMe_2019/LRRK2"; Roadmap_chromatinMarks_search <- "monocytes"; CredibleSet_only = T; permutations = 1000; chromatin_state = NA;  R2_filter = 0.8
  # Cleanup pyc tmp files
  suppressWarnings(file.remove(
    file.path(goshifter_path, 
              c("chromtree.pyc",
                "data.pyc",
                "docopt.pyc",
                "functions.pyc")) ) )
  # Term key for chromatin states
  chromState_key <- data.table::fread(file.path("./echolocatoR/tools/Annotations/ROADMAP/ROADMAP_chromatinState_HMM.tsv"))

  # Create GoShifter folder
  # results_path="Data/GWAS/Nalls23andMe_2019/LRRK2"
  GS_results_path <- file.path(results_path,"GoShifter")
  dir.create(GS_results_path, showWarnings = F, recursive = T) 
  GS.RESULTS.path <- file.path(GS_results_path,
                               paste0("GoShifter.",paste(chromatin_states, collapse = "-"),".",SNP.Group,".txt"))
  
  # snpmap file  
  GoShifter.create_snpmap(snp_df, GS_results_path, verbose = verbose)
  # OR, copy file directly from example
  # file.copy(file.path(goshifter_path,"test_data","bc.snpmappings.hg19.txt"), 
  #           file.path(GS_results_path,"snpmap.txt"),overwrite = T)
  
  # LD file 
  GoShifter.create_LD(results_path, verbose = verbose) 
  # Gather Annotation BED files 
  RoadMap_ref <- GoShifter.search_ROADMAP(fuzzy_search = ROADMAP_search)
  if(any(!is.na(ROADMAP_type))){RoadMap_ref <- subset(RoadMap_ref, TYPE %in% ROADMAP_type)} 
  
  # Use existing results
  if(file.exists(GS.RESULTS.path) & force_new_goshifter==F){
    GS.RESULTS <- data.table::fread(GS.RESULTS.path)
  } else{
    # Run GS over each chromatin mark subset
    GS.RESULTS <- lapply(chromatin_states, function(chromatin_state){ 
      try({ 
        printer("+ GoShifter:: Running on chromatin state subset:",chromatin_state)
        GS_results <- GoShifter.run(goshifter_path = goshifter_path,
                                    RoadMap_ref = RoadMap_ref, 
                                    results_path = results_path, 
                                    chromatin_state = chromatin_state,
                                    R2_filter = R2_filter,
                                    permutations = permutations, 
                                    remove_tmps = remove_tmps,  
                                    verbose = verbose,
                                    overlap_threshold = overlap_threshold)
        message("GoShifter results data.frame: ",nrow(GS_results), " rows x ", ncol(GS_results)," cols" )
      })  
      return(GS_results)
    }) %>% data.table::rbindlist(fill=T) 
  }
  
  if(save_results){ 
    data.table::fwrite(GS.RESULTS, 
                       GS.RESULTS.path,
                       sep="\t", quote = F)
  } 
  if(remove_tmps){
    # Remove GS output files to avoid reading in old files in future runs
    suppressWarnings(
      file.remove(file.path(GS_results_path,
                            paste0("GS_results",c(".nperm1000.enrich",
                                                  ".nperm1000.locusscore",
                                                  ".nperm1000.snpoverlap"))) ) )
  }
  return(GS.RESULTS)
}






# ######## PLOTS  # ######## 

GoShifter.histograms_pvals <- function(GS_results){
  GS <- GS_results[,c("Annotation","pval","GS.pval")] %>% unique() %>% 
    data.table::melt.data.table(id.vars = "Annotation", variable.name = "P.value Source", value.name = "P.value")
  hst <- ggplot(data = GS, aes(x=P.value, fill=`P.value Source`)) + 
          geom_histogram(alpha=.5) +  
          theme_bw()
  print(hst)
}

GoShifter.histograms_SNPgroups <- function(GS.groups, show_plot = T){
  GS <- GS.groups[,c("SNP.Group","Annotation","pval","GS.pval")] %>% unique()
  hst <- ggplot(GS, aes(x=pval, fill=SNP.Group)) + 
    geom_histogram(alpha=.7)
  if(show_plot){print(hst)}
  return(hst)
}

GoShifter.heatmap <- function(GS.groups, show_plot = T){
  # GS_results <- GS_results_CS
  sig.tests <- subset(GS.groups, pval<=0.05)
  mat <- reshape2::acast(sig.tests, 
                         Epigenome.name ~ chromatin_state + SNP.Group, 
                         value.var="mean.enrichment", 
                         fun.aggregate = mean, drop = F, fill = 0)  
  # library(heatmaply)
  hm <- heatmaply::heatmaply(mat, dendrogram = "row", k_row=3)
  # hm <- htmltools::tagList(list(hm))
  if(show_plot){print(hm)}
  return(hm)
}



 