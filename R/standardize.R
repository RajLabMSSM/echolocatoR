#' Standardize the locus subset
#'
#' After querying a subset of the full summary statistics,
#' this function converts it into a standardized format
#' that the rest of \emph{echolocatoR} can work with.
#'
#' @family standardizing functions
#' @inheritParams finemap_locus
#' @importFrom echodata get_sample_size
#' @importFrom echotabix liftover
#' @examples
#' BST1 <- echodata::BST1
#' #### Screw up Freq to see if function can fix it and infer MAF ####
#' BST1$rsid <- BST1$SNP
#' BST1 <- data.frame(BST1)[,!colnames(BST1) %in% c("MAF","SNP")]
#' BST1[c(10,30,55),"Freq"] <- 0
#' BST1[c(12,22),"Freq"] <- NA
#'
#' subset_path <- file.path(tempdir(),"BST1.tsv")
#' data.table::fwrite(BST1, subset_path)
#' query_mod <- echolocatoR:::standardize_subset(subset_path=subset_path,
#'                                               locus="BST1",
#'                                               snp_col = "rsid",
#'                                               MAF_col="calculate")
standardize_subset <- function(subset_path,
                               locus,
                               top_SNPs=NULL,
                               fullSS_genome_build="hg19",
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
                               qtl_prefixes=NULL,
                               return_dt=TRUE,
                               nThread=1,
                               download_method="axel",
                               verbose=TRUE){
  messager("LD:: Standardizing summary statistics subset.",v=verbose)
  query_check <- data.table::fread(subset_path, nrows = 2)

  if(dim(query_check)[1]==0){
    file.remove(subset_path)
    stop("Could not find any rows in full data that matched query :(")
  } else{
    query <- data.table::fread(subset_path,
                               header=TRUE, stringsAsFactors = FALSE,
                               nThread = nThread)
    #### Impute StdErr ####
    if(stderr_col=="calculate"){
      messager("Calculating Standard Error...",v=verbose)
      query$StdErr <- subset(query, select=effect_col) /
        subset(query, select=tstat_col)
      stderr_col <- "StdErr"
    }
    #### Rename subset DF ####
    query_mod <- query %>%
      subset(select=c(chrom_col, position_col, snp_col,
                      pval_col, effect_col, stderr_col)) %>%
      dplyr::rename(CHR=chrom_col,POS=position_col, SNP=snp_col, P=pval_col,
                    Effect=effect_col, StdErr=stderr_col)
    if(!"SNP" %in% colnames(query)){
      query <- dplyr::rename(query, SNP=snp_col)
    }
    #### Gene col ####
    messager("++ Preparing Gene col", v=verbose)
    if(gene_col %in% colnames(query)){
      query <- dplyr::rename(query, Gene=gene_col)
      query_mod$Gene <- query$Gene
      # Subset by eGene
      if(detect_genes(loci = locus, verbose  = FALSE)){
        messager("+ Subsetting to gene =",names(locus), v=verbose)
        query_mod <- subset(query_mod, Gene==names(locus))
        query <- subset(query, Gene==names(locus))
        if(dplyr::n_distinct(query_mod$SNP)!=nrow(query_mod)) {
          stop("N rows must be equal to N unique SNPS.")
        }
      }
    }
    #### Liftover if needed ####
    ## Do this step BEFORE inferring MAF from external source
    if(!toupper(fullSS_genome_build) %in% c("HG19","HG37","GRCH37")){
      query_mod <- echotabix::liftover(
        sumstats_dt = query_mod,
        ref_genome = "HG38",
        convert_ref_genome = "HG37",
        chrom_col = "CHR",
        start_col = "POS",
        as_granges = FALSE,
        verbose = verbose)
    }
    # Have to subset both ways bc some SNPs only available in certain builds
    query <- subset(query, SNP %in% unique(query_mod$SNP))
    query_mod <- subset(query_mod, SNP %in% unique(query$SNP))
    #### Add ref/alt alleles if available ####
    messager("++ Preparing A1,A1 cols", v=verbose)
    if(any(A1_col %in% colnames(query)) & any(A2_col %in% colnames(query))){
      query2 <- query %>% dplyr::rename(A1=A1_col, A2=A2_col)
      query_mod$A1 <- toupper(query2$A1)
      query_mod$A2 <- toupper(query2$A2)
    }

    # ------ Optional columns ------ #
    #### Infer MAF from freq ####
    ## Assumes MAF is always less than 0.5 and that all SNPs are biallelic
    messager("++ Preparing MAF,Freq cols", v=verbose)
    if(any(MAF_col %in% colnames(query))){
      query <- query %>% dplyr::rename(MAF=MAF_col)
      query_mod$MAF <- as.numeric(query$MAF)
    } else {
      if(any(freq_col %in% colnames(query))){
        query <- query %>% dplyr::rename(Freq=freq_col)
        query_mod$Freq <- as.numeric(query$Freq)
        if(MAF_col=="calculate" | !(MAF_col %in% colnames(query)) ){
          messager("++ Inferring MAF from frequency column...")
          query_mod$MAF <- ifelse(abs(query$Freq<0.5),
                                  abs(query$Freq),
                                  abs(1-query$Freq))
        }
      } else {messager("++ Could not infer MAF",v=verbose)}
    }
    if(any("MAF" %in% colnames(query))){
      query_mod$MAF <- abs(query_mod$MAF)
      messager("++ Removing SNPs with MAF== 0 | NULL | NA", v=verbose)
      query_mod <- subset(query_mod, !(is.na(MAF) | is.null(MAF) | MAF==0))
      query <- subset(query, SNP %in% unique(query_mod$SNP))
    }

    ## Add proportion of cases if available
    messager("++ Preparing N_cases,N_controls cols", v=verbose)
    if(any(N_cases_col %in% colnames(query)) &
       any(N_controls_col %in% colnames(query))){
      query <- query %>%
        dplyr::rename(N_cases=N_cases_col,
                      N_controls=N_controls_col)
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

    messager("++ Preparing `proportion_cases` col", v=verbose)
    if(proportion_cases !="calculate"){
      if(length(proportion_cases)==nrow(query)){
        query_mod$proportion_cases <- query[[proportion_cases]]
      } else {
        query_mod$proportion_cases <- proportion_cases
      }
    } else if(proportion_cases=="calculate" &
              "N_cases" %in% colnames(query_mod) &
              "N_controls" %in% colnames(query_mod)){
      messager("++ Calculating `proportion_cases`.")
      ### Calculate proportion of cases if N_cases and N_controls available
      query_mod$proportion_cases <- query_mod$N_cases /
        (query_mod$N_controls + query_mod$N_cases)
    } else {
      ### Otherwise don't include this col
      messager("++ 'proportion_cases' not included in data subset.")
    }
    #### Impute sample size ####
    messager("++ Preparing N col", v=verbose)
    query_mod <- echodata::get_sample_size(dat = query_mod,
                                   method = sample_size,
                                   force_new = FALSE,
                                   verbose=verbose)
    #### Impute t-stat ####
    messager("++ Preparing t-stat col", v=verbose)
    query_mod$t_stat <- calculate_tstat(finemap_dat = query_mod,
                                        tstat_col = tstat_col)[["t_stat"]]

    if(any(query_mod$P<=0)){
      messager("++ Replacing P-values==0 with", .Machine$double.xmin, v=verbose)
      query_mod[(query_mod$P<=0),"P"] <- .Machine$double.xmin
    }
    if(any(is.na(query_mod$P))){
      messager("++ Removing SNPs with P-values==NA", v=verbose)
      query_mod <- subset(query_mod, !is.na(P))
    }

    messager("++ Assigning lead SNP", v=verbose)
    # Add leadSNP col
    if(is.null(top_SNPs)){
      top_SNPs <- cbind(Locus=locus,(query_mod %>% dplyr::arrange(P))[1,])
    }
    topSNP_sub <- auto_topSNPs_sub(top_SNPs = top_SNPs,
                                   query = query_mod,
                                   locus = locus)
    query_mod$leadSNP <- query_mod$SNP==topSNP_sub$SNP


    # Remove any duplicate columns
    query_mod <- data.frame(query_mod)[,!duplicated(colnames(query_mod))]

    messager("++ Ensuring Effect, StdErr, P are numeric", v=verbose)
    ## Only convert to numeric AFTER removing NAs
    ## (otherwise as.numeric will turn them into 0s)
    query_mod <- suppressWarnings(
      query_mod  %>%
      dplyr::mutate(CHR=as.integer(gsub("chr","",CHR)),
                    Effect=as.numeric(Effect),
                    StdErr=as.numeric(StdErr),
                    P=as.numeric(P))
      )

    ### Add QTL cols ####
    if(!is.null(qtl_prefixes)){
      messager("++ Adding back cols starting with:",
               paste(qtl_prefixes,collapse = ","), v=verbose)
      qtl_cols <- grep(paste(qtl_prefixes,collapse = "|"),
                       colnames(query), value = TRUE)
      query_mod <- cbind(data.frame(query_mod),
                         data.frame(query)[,qtl_cols])
    }

    messager("++ Ensuring 1 SNP per row", v=verbose)
    ## Get just one SNP per location (just pick the first one)
    query_mod <- query_mod %>%
      dplyr::group_by(CHR, POS) %>%
      dplyr::slice(1) %>%
      dplyr::group_by(SNP) %>%
      dplyr::slice(1)
    # dplyr doesn't like working with grouped tables downstream.
    query_mod <- dplyr::ungroup(query_mod)
    # Trim whitespaces
    ## Extra whitespace causes problems when you try
    ## to make space-delimited files.
    messager("++ Removing extra whitespace", v=verbose)
    cols_to_be_rectified <- names(query_mod)[
      vapply(query_mod, is.character, logical(1))
    ]
    query_mod <- query_mod %>%
      dplyr::mutate_at(.vars =dplyr::vars(cols_to_be_rectified),
                       .funs = trimws )
    messager("++ Saving subset ==>",subset_path, v=verbose)
    dir.create(dirname(subset_path), showWarnings = FALSE, recursive = TRUE)
    data.table::fwrite(query_mod,
                       subset_path,
                       sep = "\t", nThread = nThread)
    if(return_dt==TRUE){
      return(data.table::data.table(query_mod))
    } else {
        return(subset_path)
      }
  }
}
