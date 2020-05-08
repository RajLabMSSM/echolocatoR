###### STANDARDIZE ######




check_if_empty <- function(file_path, file_sep="\t"){
  rowCheck <- dim(data.table::fread(file_path, nrows = 2, sep = file_sep))[1]
  if(rowCheck==0){
    stop("No SNPs identified within the summary stats file that met your criterion. :o")
  } else {printer("+ Subset file looks good! :)")}
}



column_dictionary <- function(file_path){
  # Get the index of each column name
  # f <- data.table::fread(file_path, nrows = 0, header = T)

  cmd <- paste(ifelse(endsWith(file_path,".gz"),"gunzip -c","cat"),
               file_path,"| head -1")
  # system(cmd)
  cNames <- colnames(data.table::fread(cmd=cmd))
  colDict <- setNames(1:length(cNames), cNames  )
  return(colDict)
}


auto_topSNPs_sub <- function(top_SNPs, query, gene){
  # If no top_SNPs dataframe is supplied,
  ## this function will sort by p-value and then effect size,
  ## and use the SNP in the first row.
  if(toString(top_SNPs)=="auto"){
    top_SNPs <- query %>% mutate(Gene=gene) %>%
      arrange(P, desc(Effect)) %>% group_by(Gene) %>% dplyr::slice(1)
  }
  topSNP_sub <- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),][1,]
  return(topSNP_sub)
}


calculate.tstat <- function(finemap_DT, tstat_col="t_stat"){
  if(tstat_col %in% colnames(finemap_DT)){
    finemap_DT <- finemap_DT %>% dplyr::rename(t_stat = tstat_col)
  } else if(("Effect" %in% colnames(finemap_DT)) & ("StdErr" %in% colnames(finemap_DT))){
    printer("+ Calculating t-statistic from Effect and StdErr...")
    finemap_DT <- finemap_DT %>% dplyr::mutate(t_stat =  Effect/StdErr)
  } else {
    printer("+ Could not calculate t-stat due to missing Effect and/or StdErr columns. Returning input data.")
  }
  return(data.table::data.table(finemap_DT))
}



UKB.get_MAF <- function(subset_DT,
                        output.path = "./Data/Reference/UKB_MAF",
                        force_new_maf = F){
  printer("UKB MAF:: Extracting MAF from UKB reference.")
  # Documentation: http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=22801
  # subset_DT = data.table::fread("Data/GWAS/Kunkle_2019/PTK2B/PTK2B_Kunkle_2019_subset.tsv.gz")
  chrom <- unique(subset_DT$CHR)
  input.url <- paste0("biobank.ctsu.ox.ac.uk/showcase/showcase/auxdata/ukb_mfi_chr",chrom,"_v3.txt")
  out.file <- file.path(output.path, basename(input.url))
  if(file.exists(out.file) & force_new_maf==F){
    printer("+ UKB MAF:: Importing pre-existing file")
  } else{
    out.file <- axel(input.url = input.url,
                     output.path = output.path,
                     background = F)
  }
  maf <- data.table::fread(out.file, nThread = 4,
                           select = c(3,6),
                           col.names = c("POS","MAF"))
  maf <- subset(maf, POS %in% subset_DT$POS)
  merged_DT <- data.table:::merge.data.table(subset_DT, maf,
                                             by = "POS")
  return(merged_DT)
}




standardize_subset <- function(gene,
                              top_SNPs=NULL,
                              subset_path="./Data",
                              file_sep="\t",
                              chrom_col="CHR",
                              position_col="POS",
                              snp_col="SNP",
                              pval_col="P",
                              effect_col="Effect",
                              stderr_col="StdErr",
                              tstat_col="t_stat",
                              MAF_col="MAF",
                              freq_col="Freq",
                              N_cases_col="N_cases",
                              N_controls_col="N_controls",
                              N_cases=NULL,
                              N_controls=NULL,
                              proportion_cases="calculate",
                              A1_col="A1",
                              A2_col="A2",
                              return_dt=T,
                              verbose=T){
  printer("",v=verbose)
  message("---------------- Step 1.5: Standardize ----------")
  query_check <- data.table::fread(subset_path, nrows = 2)

  if(dim(query_check)[1]==0){
    file.remove(subset_path)
    stop("Could not find any rows in full data that matched query :(")
  } else{
    query <- data.table::fread(subset_path, header=T, stringsAsFactors = F)
    ## Calculate StdErr
    if(stderr_col=="calculate"){
      printer("Calculating Standard Error...")
      query$StdErr <- subset(query, select=effect_col) / subset(query, select=tstat_col)
      stderr_col <- "StdErr"
    }

    ## Rename subset DF
    query_mod <- query %>% subset(select=c(chrom_col, position_col, snp_col, pval_col, effect_col, stderr_col)) %>%
      dplyr::rename(CHR=chrom_col,POS=position_col, SNP=snp_col, P=pval_col,
                    Effect=effect_col, StdErr=stderr_col)

    # Add ref/alt alleles if available
    if(A1_col %in% colnames(query) & A2_col %in% colnames(query)){
      query2 <- query %>% dplyr::rename(A1=A1_col, A2=A2_col)
      query_mod$A1 <- query2$A1
      query_mod$A2 <- query2$A2
    }

    # ------ Optional columns ------ #
    ## Infer MAF from freq (assuming MAF is alway less than 0.5)
    if(MAF_col %in% colnames(query)){
      query <- query %>% dplyr::rename(MAF=MAF_col)
      query_mod$MAF <- as.numeric(query$MAF)
    } else if(freq_col %in% colnames(query)){
      query <- query %>% dplyr::rename(Freq=freq_col)
      query_mod$Freq <- as.numeric(query$Freq)
      if(MAF_col=="calculate" | !(MAF_col %in% colnames(query)) ){
        printer("+Inferring MAF from frequency column...")
        query_mod$MAF <- ifelse(query$Freq<0.5, query$Freq, abs(1-query$Freq))
      }
    } else {
      # As a last resort download UKB MAF
      query_mod <- UKB.get_MAF(subset_DT = query_mod,
                                output.path = "./Data/Reference/UKB_MAF",
                                force_new_maf = F)
    }
    query_mod$MAF <- abs(query_mod$MAF)



    ## Add proportion of cases if available
    if(N_cases_col %in% colnames(query) & N_controls_col %in% colnames(query)){
      query <- query %>% dplyr::rename(N_cases=N_cases_col, N_controls=N_controls_col)
      query_mod$N_cases <- query$N_cases
      query_mod$N_controls <- query$N_controls
    } else {
      query_mod$N_cases <- N_cases
      query$N_cases <- N_cases
      N_cases_col <- "N_cases"
      query_mod$N_controls <- N_controls
      query$N_controls <- N_controls
      N_controls_col <- "N_controls"
    }

    if(proportion_cases !="calculate"){
      query_mod$proportion_cases <- query[proportion_cases]
    } else if(proportion_cases=="calculate" &
              "N_cases" %in% colnames(query_mod) &
              "N_controls" %in% colnames(query_mod)){
      printer("+ Standardize:: Calculating proportion of cases.")
      ### Calculate proportion of cases if N_cases and N_controls available
      query_mod$proportion_cases <- query_mod$N_cases / (query_mod$N_controls + query_mod$N_cases)
    } else {
      ### Otherwise don't include this col
      printer("'proportion of cases' not included in data subset.")
    }



    query_mod$t_stat <- calculate.tstat(finemap_DT = query_mod,
                                        tstat_col = tstat_col)$t_stat


    if(is.null(top_SNPs)){top_SNPs <- cbind(Gene=gene,(query_mod %>% arrange(P))[1,])}
    topSNP_sub <- auto_topSNPs_sub(top_SNPs, query_mod, gene)

    ## Remove SNPs with NAs in stats
    query_mod[(query_mod$P<=0)|(query_mod$P>1),"P"] <- 1

    # Add leadSNP col
    ## Get just one SNP per location (just pick the first one)
    query_mod <- query_mod %>% group_by(CHR, POS) %>% dplyr::slice(1)
    ## Mark lead SNP
    query_mod$leadSNP <- ifelse(query_mod$SNP==topSNP_sub$SNP, T, F)
    if(sum(query_mod$leadSNP)==0){
      printer("+ leadSNP missing. Assigning new one by min p-value.")
      top.snp <- head(arrange(query_mod, P, desc(Effect)))[1,]$SNP
      query_mod$leadSNP <- ifelse(query_mod$SNP==top.snp,T,F)
    }

    # Only convert to numeric AFTER removing NAs (otherwise as.numeric will turn them into 0s)
    query_mod <- query_mod  %>%
      mutate(Effect=as.numeric(Effect), StdErr=as.numeric(StdErr), P=as.numeric(P))

    
    # Trim whitespaces
    ## Extra whitespace causes problems when you try to make space-delimited files
    cols_to_be_rectified <- names(query_mod)[vapply(query_mod, is.character, logical(1))]
    query_mod <- query_mod %>% mutate_at(vars(cols_to_be_rectified),
                                   funs(trimws) )
    
    data.table::fwrite(query_mod, subset_path, sep = "\t", nThread = 4)
    if(return_dt==T){return(query_mod)}
  }
}




extract_SNP_subset <- function(gene,
                               fullSS_path,
                               results_path,
                               force_new_subset=F,
                               top_SNPs="auto",
                               bp_distance=500000,
                               chrom_col="CHR",
                               position_col="POS",
                               snp_col="SNP",
                               gene_col="Gene",
                               pval_col="P",
                               effect_col="Effect",
                               stderr_col="StdErr",
                               MAF_col="MAF",
                               freq_col = "Freq",
                               tstat_col="t-stat",
                               A1_col = "A1",
                               A2_col = "A2",

                               N_cases_col="N_cases",
                               N_controls_col="N_controls",
                               N_cases=NULL,
                               N_controls=NULL,
                               proportion_cases="calculate",

                               superpopulation="",
                               min_POS=NA,
                               max_POS=NA,
                               file_sep="\t",
                               query_by="coordinates",
                               probe_path = "./Data/eQTL/gene.ILMN.map",
                               remove_tmps=T,
                               verbose=T){
  printer("",v=verbose)
  message("------------------ Step 1: Query ---------------")
  subset_path <- get_subset_path(results_path=results_path, gene=gene, subset_path="auto")
  # topSNP_sub <- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),][1,]
  # if(is.na(min_POS)){min_POS <- topSNP_sub$POS - bp_distance}
  # if(is.na(max_POS)){max_POS <- topSNP_sub$POS + bp_distance}
  # printer("---Min snp position:",min_POS, "---")
  # printer("---Max snp position:",max_POS, "---")
  multi_path <- file.path(dirname(subset_path),"Multi-finemap/Multi-finemap_results.txt")
  if(file.exists(subset_path) & force_new_subset==F){
    printer("+ Subset file already exists. Importing",subset_path,"...")
    check_if_empty(subset_path, file_sep=file_sep)
    query <- data.table::fread(subset_path, header=T, stringsAsFactors=F, sep="\t")
  } else if (file.exists(multi_path) & force_new_subset==F){
    printer("+ Importing previous Multi-finemap file...Importing summary stats.")
    check_if_empty(multi_path, file_sep=file_sep)
    query <- data.table::fread(multi_path, header=T, stringsAsFactors=F, sep="\t")
  } else {
    # Extract subset with awk
    printer("+ Extracting relevant variants from fullSS...")
    start_query <- Sys.time()
    # Function selects different methods of querying your SNPs
    query_handler(gene=gene,
                  fullSS_path=fullSS_path,
                  subset_path=subset_path,
                  top_SNPs=top_SNPs,
                  gene_col=gene_col,
                  chrom_col=chrom_col,
                  position_col=position_col,
                  file_sep=file_sep,
                  min_POS=min_POS,
                  max_POS=max_POS,
                  bp_distance=bp_distance,
                  query_by=query_by,
                  probe_path=probe_path)
    # Clean file
    query <- standardize_subset(gene=gene,
                               top_SNPs=top_SNPs,
                               subset_path=subset_path,
                               file_sep=file_sep,
                               chrom_col=chrom_col,
                               position_col=position_col,
                               snp_col=snp_col,
                               pval_col=pval_col,
                               effect_col=effect_col,
                               stderr_col=stderr_col,
                               tstat_col=tstat_col,
                               MAF_col=MAF_col,
                               freq_col=freq_col,
                               A1_col = A1_col,
                               A2_col = A2_col,

                               N_cases_col=N_cases_col,
                               N_controls_col=N_controls_col,
                               N_cases=N_cases,
                               N_controls=N_controls,
                               proportion_cases = proportion_cases)
    end_query <- Sys.time()
    printer("+ Extraction completed in", round(end_query-start_query, 2),"seconds")
    printer("+", dim(query)[1], "SNPs x ",dim(query)[2],"columns")
  }
  if(remove_tmps){
    printer("+ Removing subset tmp...")
    suppressWarnings(file.remove(subset_path))
  }
  return(query)
}


