###### STANDARDIZE ######


#' Automatically identify top SNP per locus
#'
#'  If no \code{top_SNPs} dataframe is supplied,
#'  this function will sort by p-value and then effect size,
#'  and use the SNP in the first row.
#'
#' @family standardization functions
auto_topSNPs_sub <- function(top_SNPs,
                             query,
                             locus){
  if(toString(top_SNPs)=="auto"){
    top_SNPs <- query %>% dplyr::mutate(Locus=locus) %>%
      dplyr::arrange(P, dplyr::desc(Effect)) %>% dplyr::group_by(Locus) %>% dplyr::slice(1)
  }
  topSNP_sub <- top_SNPs[top_SNPs$Locus==locus & !is.na(top_SNPs$Locus),][1,]
  return(topSNP_sub)
}




#' Compute t-stat
#'
#' If \strong{tstat} column is missing,
#' compute t-statistic from: \code{Effect / StdErr}.
#' @family standardization functions
calculate_tstat <- function(finemap_dat,
                            tstat_col="t_stat"){
  if(tstat_col %in% colnames(finemap_dat)){
    finemap_dat <- finemap_dat %>% dplyr::rename(t_stat = tstat_col)
  } else {
    if(("Effect" %in% colnames(finemap_dat)) & ("StdErr" %in% colnames(finemap_dat))){
      printer("+ Calculating t-statistic from Effect and StdErr...")
      finemap_dat <- finemap_dat %>% dplyr::mutate(t_stat =  Effect/StdErr)
    } else {
      printer("+ Could not calculate t-stat due to missing Effect and/or StdErr columns. Returning input data.")
    }
  }
  return(data.table::data.table(finemap_dat))
}


#' Get MAF from UK Biobank.
#'
#' If \strong{MAF} column is missing,
#' download MAF from UK Biobank and use that instead.
#'
#' @family standardizing functions
#' @examples
#' data("BST1");
#' subset_DT <- data.frame(BST1)[,colnames(BST1)!="MAF"]
#' BST1 <- get_UKB_MAF(subset_DT=subset_DT )
get_UKB_MAF <- function(subset_DT,
                        output_path = "./Data/Reference/UKB_MAF",
                        force_new_maf = F,
                        download_method="axel",
                        nThread=4){
  printer("UKB MAF:: Extracting MAF from UKB reference.")
  # Documentation: http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=22801
  # subset_DT = data.table::fread("Data/GWAS/Kunkle_2019/PTK2B/PTK2B_Kunkle_2019_subset.tsv.gz")
  chrom <- unique(subset_DT$CHR)
  input_url <- paste0("biobank.ctsu.ox.ac.uk/showcase/showcase/auxdata/ukb_mfi_chr",chrom,"_v3.txt")
  out_file <- file.path(output_path, basename(input_url))
  if(file.exists(out_file) & force_new_maf==F){
    printer("+ UKB MAF:: Importing pre-existing file")
  } else{
    out_file <- downloader(input_url = input_url,
                           output_path = output_path,
                           background = F,
                           download_method = download_method,
                           nThread = nThread)
  }
  maf <- data.table::fread(out_file, nThread = nThread,
                           select = c(3,6),
                           col.names = c("POS","MAF"))
  maf <- subset(maf, POS %in% subset_DT$POS)
  merged_dat <- data.table::merge.data.table(subset_DT, maf,
                                             by = "POS") %>%
    # Make sure each SNP just appears once
    dplyr::group_by(merged_dat, SNP) %>%
    dplyr::slice(1)
  return(merged_dat)
}



#' Stanardize the locus subset
#'
#' After querying a subset of the full summary statistics,
#' this function converts it into a standardized format
#' that the rest of \emph{echolocatoR} can work with.
#'
#' @family standardizing functions
#' @inheritParams finemap_pipeline
#' @examples
#' data("BST1")
#' ... Screw up Freq to see if function can fix it and infer MAF ...
#' BST1$rsid <- BST1$SNP
#' BST1 <- data.frame(BST1)[,!colnames(BST1) %in% c("MAF","SNP")]
#' BST1[c(10,30,55),"Freq"] <- 0
#' BST1[c(12,22),"Freq"] <- NA
#' data.table::fwrite(BST1, "~/Desktop/results/GWAS/Nalls23andMe_2019/BST1/BST1.tsv")
#' query_mod <- standardize_subset(locus="BST1", subset_path="~/Desktop/results/GWAS/Nalls23andMe_2019/BST1/BST1.tsv", MAF_col="calculate", snp_col="rsid")
standardize_subset <- function(locus,
                              top_SNPs=NULL,
                              subset_path="./Data",
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
                              sample_size=NULL,
                              A1_col="A1",
                              A2_col="A2",
                              gene_col="Gene",
                              QTL_prefixes=NULL,
                              return_dt=T,
                              nThread=4,
                              download_method="axel",
                              verbose=T){
  printer("",v=verbose)
  message("---------------- Step 1.5: Standardize ----------")
  query_check <- data.table::fread(subset_path, nrows = 2)

  if(dim(query_check)[1]==0){
    file.remove(subset_path)
    stop("Could not find any rows in full data that matched query :(")
  } else{
    query <- data.table::fread(subset_path,
                               header=T, stringsAsFactors = F,
                               nThread = nThread)
    ## Calculate StdErr
    if(stderr_col=="calculate"){
      printer("Calculating Standard Error...",v=verbose)
      query$StdErr <- subset(query, select=effect_col) / subset(query, select=tstat_col)
      stderr_col <- "StdErr"
    }

    ## Rename subset DF
    query_mod <- query %>% subset(select=c(chrom_col, position_col, snp_col, pval_col, effect_col, stderr_col)) %>%
      dplyr::rename(CHR=chrom_col,POS=position_col, SNP=snp_col, P=pval_col,
                    Effect=effect_col, StdErr=stderr_col)

    printer("++ Preparing Gene col", v=verbose)
    if(gene_col %in% colnames(query)){
      query <- dplyr::rename(query, Gene=gene_col)
      query_mod$Gene <- query$Gene
      if(detect_genes(loci = locus, verbose = F)){
        printer("+ Subsetting to gene =",names(locus), v=verbose)
        query_mod <- subset(query_mod, Gene==names(locus))
        query <- subset(query, Gene==names(locus))
        if(dplyr::n_distinct(query_mod$SNP)!=nrow(query_mod)) stop("N rows must be equal to N unique SNPS.")
      }
    }

    # Add ref/alt alleles if available
    printer("++ Preparing A1,A1 cols", v=verbose)
    if(A1_col %in% colnames(query) & A2_col %in% colnames(query)){
      query2 <- query %>% dplyr::rename(A1=A1_col, A2=A2_col)
      query_mod$A1 <- query2$A1
      query_mod$A2 <- query2$A2
    }

    # ------ Optional columns ------ #
    ## Infer MAF from freq (assuming MAF is alway less than 0.5)
    printer("++ Preparing MAF,Freq cols", v=verbose)
    if(MAF_col %in% colnames(query)){
      query <- query %>% dplyr::rename(MAF=MAF_col)
      query_mod$MAF <- as.numeric(query$MAF)
    } else {
      if(freq_col %in% colnames(query)){
        query <- query %>% dplyr::rename(Freq=freq_col)
        query_mod$Freq <- as.numeric(query$Freq)
        if(MAF_col=="calculate" | !(MAF_col %in% colnames(query)) ){
          printer("++ Inferring MAF from frequency column...")
          query_mod$MAF <- ifelse(abs(query$Freq<0.5), abs(query$Freq), abs(1-query$Freq))
        }
      } else {
        # As a last resort download UKB MAF
        query_mod <- get_UKB_MAF(subset_DT = query_mod,
                                 output_path = file.path(dirname(dirname(dirname(subset_path))),
                                                         "Reference/UKB_MAF"),
                                 force_new_maf = F,
                                 nThread = nThread,
                                 download_method = download_method)
      }
    }
    query_mod$MAF <- abs(query_mod$MAF)
    printer("++ Removing SNPs with MAF== 0 | NULL | NA", v=verbose)
    query_mod <- subset(query_mod, !(is.na(MAF) | is.null(MAF) | MAF==0))
    query <- subset(dplyr::rename(query, SNP=snp_col), SNP %in% unique(query_mod$SNP))

    ## Add proportion of cases if available
    printer("++ Preparing N_cases,N_controls cols", v=verbose)
    if(N_cases_col %in% colnames(query) & N_controls_col %in% colnames(query)){
      query <- query %>% dplyr::rename(N_cases=N_cases_col, N_controls=N_controls_col)
      query_mod$N_cases <- query$N_cases
      query_mod$N_controls <- query$N_controls
    } else {
      if(!is.null(N_cases)){
        query_mod$N_cases <- N_cases
        query$N_cases <- N_cases
        N_cases_col <- "N_cases"
      }
      if(!is.null(N_controls)){
        query_mod$N_controls <- N_controls
        query$N_controls <- N_controls
        N_controls_col <- "N_controls"
      }
    }

    printer("++ Preparing `proportion_cases` col", v=verbose)
    if(proportion_cases !="calculate"){
      if(length(proportion_cases)==nrow(query)){
        query_mod$proportion_cases <- query[[proportion_cases]]
      } else {
        query_mod$proportion_cases <- proportion_cases
      }
    } else if(proportion_cases=="calculate" &
              "N_cases" %in% colnames(query_mod) &
              "N_controls" %in% colnames(query_mod)){
      printer("++ Calculating `proportion_cases`.")
      ### Calculate proportion of cases if N_cases and N_controls available
      query_mod$proportion_cases <- query_mod$N_cases / (query_mod$N_controls + query_mod$N_cases)
    } else {
      ### Otherwise don't include this col
      printer("++ 'proportion_cases' not included in data subset.")
    }

    # Calculate sample size
    printer("++ Preparing N col", v=verbose)
    query_mod$N <- get_sample_size(subset_DT = query_mod,
                                   sample_size = sample_size,
                                   verbose=verbose)[["N"]]

    printer("++ Preparing t-stat col", v=verbose)
    query_mod$t_stat <- calculate_tstat(finemap_dat = query_mod,
                                        tstat_col = tstat_col)[["t_stat"]]

    if(any(query_mod$P<=0)){
      printer("++ Replacing P-values==0 with", .Machine$double.xmin, v=verbose)
      query_mod[(query_mod$P<=0),"P"] <- .Machine$double.xmin
    }
    if(any(is.na(query_mod$P))){
      printer("++ Removing SNPs with P-values==NA", v=verbose)
      query_mod <- subset(query_mod, !is.na(P))
    }

    printer("++ Assigning lead SNP", v=verbose)
    # Add leadSNP col
    if(is.null(top_SNPs)){
      top_SNPs <- cbind(Locus=locus,(query_mod %>% arrange(P))[1,])
    }
    topSNP_sub <- auto_topSNPs_sub(top_SNPs = top_SNPs,
                                   query = query_mod,
                                   locus = locus)
    query_mod$leadSNP <- query_mod$SNP==topSNP_sub$SNP

    printer("++ Ensuring Effect, StdErr, P are numeric", v=verbose)
    # Only convert to numeric AFTER removing NAs (otherwise as.numeric will turn them into 0s)
    query_mod <- query_mod  %>%
      dplyr::mutate(Effect=as.numeric(Effect),
                    StdErr=as.numeric(StdErr),
                    P=as.numeric(P))

    # Add QTL cols
    if(!is.null(QTL_prefixes)){
      printer("++ Adding back cols starting with:",paste(QTL_prefixes,collapse = ","), v=verbose)
      qtl_cols <- grep(paste(QTL_prefixes,collapse = "|"),colnames(query), value = T)
      query_mod <- cbind(data.frame(query_mod), data.frame(query)[,qtl_cols])
    }

    printer("++ Ensuring 1 SNP per row", v=verbose)
    ## Get just one SNP per location (just pick the first one)
    query_mod <- query_mod %>%
      dplyr::group_by(CHR, POS) %>%
      dplyr::slice(1) %>%
      dplyr::group_by(SNP) %>%
      dplyr::slice(1)
    # dplyr doesn't like working with grouped tables downstream.
    query_mod <- dplyr::ungroup(query_mod)
    # Trim whitespaces
    ## Extra whitespace causes problems when you try to make space-delimited files
    printer("++ Removing extra whitespace", v=verbose)
    cols_to_be_rectified <- names(query_mod)[vapply(query_mod, is.character, logical(1))]
    query_mod <- query_mod %>% dplyr::mutate_at(.vars = vars(cols_to_be_rectified),
                                                .funs = trimws )
    printer("++ Saving subset ==>",subset_path, v=verbose)
    dir.create(dirname(subset_path), showWarnings = F, recursive = T)
    data.table::fwrite(query_mod, subset_path,
                       sep = "\t", nThread = nThread)
    if(return_dt==T){return(query_mod)}
  }
}



