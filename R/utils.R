
#### #### #### #### #### #### #### ####
#### GENERAL AND UTILITY FUNCTIONS ####
#### #### #### #### #### #### #### ####


#' @section general functions:
#' General use functions for various points in the \emph{echolocatoR} pipeline.



#' printer
#'
#' Concatenate and print any number of items.
#'
#' @family general functions
#' @examples
#' n.snps <- 50
#' printer("echolocatoR::","Processing",n.snps,"SNPs...")
#' @keywords internal
printer <- function(..., v=T){if(v){print(paste(...))}}

#' Identify current operating system (OS).
#'
#' @family general functions
#' @examples
#' get_os()
#' @keywords internal
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}


#' Ensembl IDs -> Gene symbols
#'
#' Print the echolocatoR image to the R console.
#'
#' @family general functions
#' @examples
#' startup_image()
#' @keywords internal
startup_image  <- function(){
  library(dplyr)
  try({
    col.text <- function(txt){
      c(txt,"\n") %>%
        crayon::blurred() %>%
        crayon::bgBlack() %>%
        # crayon::col_align(align = "left") %>%
        crayon::cyan() %>%
        cat()
    }
    col.text("))))))))))>>))))))))))>  E c h o l o c a t o R  <((((((((((<<((((((((((")
    col.text("")
    col.text("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ V1.0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    col.text("~~~~~~~~~~~~~~~~~~~~~~ © 2019 - Brian M. Schilder ~~~~~~~~~~~~~~~~~~~~~")
    col.text("~Department of Neuroscience, Department of Genetics & Genomic Sciences~")
    col.text("~~~~~~~~~~Icahn School of Medicine at Mount Sinai, NY, NYC, USA~~~~~~~~")
    col.text("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    grid::grid.newpage()
    img <- png::readPNG("./echolocatoR/images/echo_logo.png")
    grid::grid.raster(img)
  })
}





#' @family general functions
#' @keywords internal
.arg_list_handler <- function(arg, i){
  output <- if(length(arg)>1){arg[i]}else{arg}
  return(output)
}



#' Bind stored files
#'
#' Rapidly read a list of files from storage and concatenate them by rows.
#'
#' @family general functions
#' @examples
#' file.list <- c("data/file1.tsv", "data/file2.tsv","new_data/file3/tsv")
#' dat <- .rbind.file.list(file.list = file.list)
#' @keywords internal
.rbind.file.list <- function(file.list,
                            verbose=T,
                            nCores=4){
  merged.dat <- parallel::mclapply(file.list, function(x){
    printer(x, v = verbose)
    dat <- data.table::fread(x)
    return(dat)
  }, mc.cores = nCores) %>% data.table::rbindlist(fill=T)
  return(merged.dat)
}




#' tryCatch extension
#'
#' Extension of tryCatch function.
#'
#' @family general functions
#' @param input Function input.
#' @param func Function.
#' @keywords internal
tryFunc <- function(input, func, ...) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one
      # R expression in the "try" part then you'll have to
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression
      # in case the "try" part was completed successfully

      func(input, ...)
      # The return value of `readLines()` is the actual value
      # that will be returned in case there is no condition
      # (e.g. warning or error).
      # You don't need to state the return value via `return()` as code
      # in the "try" part is not wrapped insided a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message(paste("URL does not seem to exist:", input))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("URL caused a warning:", input))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>'
      message(paste("Processed URL:", input))
      message("Some other message at the end")
    }
  )
  return(out)
}




#' Store info necessary for Z-score
#'
#' Store the info about the full vector (length, mean, standard deviation),
#' so that you can later accurately compute Z-scores
#' when you're only working with a subset of the data.
#'
#' These functions are necessary for \code{\link{PAINTOR}}.
#' @keywords internal
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


#' Compute Z-score
#'
#' Computes Z-score when you don't have the full vector,
#' but you have the necessary info about the full vector stored in \code{z.info}.
#'
#'These functions are necessary for \code{\link{PAINTOR}}.
#' @keywords internal
Zscore <- function(x, z.info){
  # Need to use the mean and standard deviation of the FULL dataset (i.e. all beta fomr the full summary stats file)
  sample.stdv <- z.info$sample.stdv
  sample.mean <- z.info$sample.mean
  z <- (x - sample.mean) / sample.stdv
  return(z)
}


#' Compute Z-score
#'
#' Computes Z-score when you have the full vector of values (not just a subset).
#'
#' These functions are necessary for \code{\link{PAINTOR}}.
#' @keywords internal
zscore <- function(vec){
  z <- scale(vec, center = T, scale = T)
  return(z)
}

#' Interactive DT
#'
#' Generate an interactive data table with download buttons.
#'
#' @family general functions
#' @keywords internal
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




#' Interactive DT (html)
#'
#' Generate an interactive data table with download buttons.
#' Use this function when manually constructing rmarkdown chunks using cat() in a for loop.
#'
#' @family general functions
#' @keywords internal
createDT_html <- function(DF, caption="", scrollY=400){
  htmltools::tagList( createDT(DF, caption, scrollY))
}




#' Replace items within DT object
#'
#' Annoyingly, there is no native function to do simple find-and-replace in the `DT` library.
#'
#' @family general functions
#' @return data.frame
#' @keywords internal
dt.replace <- function(DT, target, replacement){
  for(col in names(DT)) set(DT, i=which(DT[[col]]==target), j=col, value=replacement)
  return(DT)
}



#' @family general functions
#' @keywords internal
view_gz_head <- function(gz_path, nrow=10){
  system(paste("zcat",gz_path,"| head",nrows))
}



#' @family general functions
#' @keywords internal
get_nrows <- function(large_file){
  printer("+ Calculating the number of rows in",basename(large_file),"...")
  if(endsWith(large_file,".gz")){
    out <- system(paste("zcat",large_file,"| wc -l"), intern=T)
  } else {
    out <- system(paste("wc -l",large_file), intern=T)
  }
  file_nrows <- as.numeric(strsplit(out," ")[[1]][1])
  printer("++ File contains",file_nrows,"rows.")
  return(file_nrows)
}



#' @family general functions
#' @keywords internal
get_header <- function(large_file){
  if(endsWith(large_file,".gz")){
    header <- system(paste("zcat",large_file,"| head -1 -"), intern=T)
  } else {
    header <- system(paste("head -1",large_file), intern=T)
  }
  return(header)
}

#' @family general functions
#' @keywords internal
check_if_empty <- function(file_path){
  rowCheck <- dim(data.table::fread(file_path, nrows = 2))[1]
  if(rowCheck==0){
    stop("No SNPs identified within the summary stats file that met your criterion. :o")
  } else {printer("+ Subset file looks good! :)")}
}




#' Map column names to positions.
#'
#' Useful in situations where you need to specify columns by index instead of name (e.g. awk queries).
#'
#' @param file_path Path to full summary stats file
#' (or any really file you want to make a column dictionary for).
#' @return Named list of column positions.
#' @keywords internal
column_dictionary <- function(file_path){
  # Get the index of each column name
  f <- data.table::fread(file_path, nrows = 0, header = T)
  cNames <- colnames(f)
  colDict <- setNames(1:length(cNames), cNames  )
  return(colDict)
}


# column_dictionary <- function(file_path){
#   # Get the index of each column name
#   # f <- data.table::fread(file_path, nrows = 0, header = T)
#
#   cmd <- paste(ifelse(endsWith(file_path,".gz"),"gunzip -c","cat"),
#                file_path,"| head -1")
#   # system(cmd)
#   cNames <- colnames(data.table::fread(cmd=cmd))
#   colDict <- setNames(1:length(cNames), cNames  )
#   return(colDict)
# }




#' Ensembl IDs -> Gene symbols
#'
#' Convert HGNC gene symbols to Ensembl IDs.
#'
#' @param gene_symbols List of HGNC gene symbols.
#' @return List of ensembl IDs.
#' @examples
#' gene_symbols <- c("FOXP2","BDNF","DCX","GFAP")
#' ensembl_IDs <- hgnc_to_ensembl(gene_symbols)
#' @keywords internal
hgnc_to_ensembl <- function(gene_symbols){
  gene_symbols[is.na(gene_symbols)] <- "NA"
  conversion <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                      keys = gene_symbols,
                                      keytype = "SYMBOL",
                                      column = "GENEID")
  return(conversion)
}

#' Gene symbols -> Ensembl IDs
#'
#' Convert Ensembl IDs to HGNC gene symbols.
#'
#' @param ensembl_ids List of ensembl IDs.
#' @return List of HGNC gene symbols.
#' @examples
#' ensembl_IDs <- c("ENSG00000128573","ENSG00000176697","ENSG00000077279","ENSG00000131095" )
#' gene_symbols <- ensembl_to_hgnc(ensembl_IDs)
#' @keywords internal
ensembl_to_hgnc <- function(ensembl_ids){
  ensembl_ids[is.na(ensembl_ids)] <- "NA"
  conversion <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                      keys = ensembl_ids,
                                      keytype = "GENEID",
                                      column = "SYMBOL")
  return(conversion)
}




#' Calculate effective sample size
#'
#' This is an often-cited approximation for the effective sample size in a case-control study.
#' i.e., this is the sample size required to identify an association with a
#'  quantitative trait with the same power as in your present study.
#'  (from email correpsondence with Omer Weissbrod)
#'
#' @param finemap_DT Preprocessed \emph{echolocatoR} locus subset file.
#' Requires the columns \strong{N_cases} and \strong{N_controls}.
#' @examples
#' finemap_DT$effective_sample_size <- effective_sample_size(finemap_DT)
#' @keywords internal
effective_sample_size <- function(finemap_DT){
  finemap_DT$N <- (4.0 / (1.0/finemap_DT$N_cases + 1.0/finemap_DT$N_controls) )
  sample_size <- as.integer(median(finemap_DT$N))
  return(sample_size)
}



#' Find the top Consensus SNP
#'
#' Identify the \code{top_N} Consensus SNP(s) per Locus,
#' defined as the Consensus SNPs with the highest mean PP across all fine-mapping tools used.
#' @keywords internal
find_topConsensus <- function(dat, top_N=1){
  # CGet top consensus SNPs
  top.consensus <- (dat %>%
                      dplyr::group_by(Locus) %>%
                      subset(Consensus_SNP) %>%
                      dplyr::arrange(-mean.PP) %>%
                      # dplyr::arrange(-IMPACT_score) %>%
                      dplyr::slice(top_N))$SNP %>% unique()
  dat$topConsensus <- dat$SNP %in% top.consensus
  return(dat)
}




#' Assign a lead GWAS SNP to a locus
#'
#' If none of the SNPs in the data.frame have \code{leadSNP==T},
#' then sort by lowest p-value (and then highest Effect size) and assign the top SNP as the lead SNP.
#'
#' @param data.frame Fine-mapping results data.frame.
#' @return Fine-mapping results data.frame with new boolean \strong{leadSNP} column,
#'  indicating whether each SNPs is the lead GWAS SNP in that locus or not.
#' @keywords internal
assign_lead_SNP <- function(new_DT, verbose=T){
  if(sum(new_DT$leadSNP)==0){
    printer("+ leadSNP missing. Assigning new one by min p-value.", v=verbose)
    top.snp <- head(arrange(new_DT, P, desc(Effect)))[1,]$SNP
    new_DT$leadSNP <- ifelse(new_DT$SNP==top.snp,T,F)
  }
  return(new_DT)
}



# --------------------------

#' @section SNP filters:
#' Functions to filter subset_DT/finemap_DT with.



#' Subset LD matrix and dataframe to only their shared SNPs
#'
#' Find the SNPs that are shared between an LD matrix and another data.frame with a `SNP` column.
#' Then remove any non-shared SNPs from both objects.
#'
#' @family SNP filters
#' @return data.frame
#' @keywords internal
subset_common_snps <- function(LD_matrix, finemap_DT, verbose=F){
  printer("+ Subsetting LD matrix and finemap_DT to common SNPs...", v=verbose)
  # Remove duplicate SNPs
  dups <- which(!duplicated(LD_matrix))
  LD_matrix <- LD_matrix[dups,dups]
  ld.snps <- row.names(LD_matrix)

  # Remove duplicate SNPs
  finemap_DT <- finemap_DT[which(!duplicated(finemap_DT$SNP)),]
  fm.snps <- finemap_DT$SNP
  common.snps <- base::intersect(ld.snps, fm.snps)
  printer("+ LD_matrix =",length(ld.snps),"SNPs.", v=verbose)
  printer("+ finemap_DT =",length(fm.snps),"SNPs.", v=verbose)
  printer("+",length(common.snps),"SNPs in common.", v=verbose)
  # Subset/order LD matrix
  new_LD <- LD_matrix[common.snps, common.snps]
  new_LD[is.na(new_LD)] <- 0
  # Subset/order finemap_DT
  finemap_DT <- data.frame(finemap_DT)
  row.names(finemap_DT) <- finemap_DT$SNP
  new_DT <- data.table::as.data.table(finemap_DT[common.snps, ])
  new_DT <- unique(new_DT)
  # Reassign the lead SNP if it's missing
  new_DT <- assign_lead_SNP(new_DT)
  # Check dimensions are correct
  if(nrow(new_DT)!=nrow(new_LD)){
    warning("+ LD_matrix and finemap_DT do NOT have the same number of SNPs.",v=verbose)
    warning("+ LD_matrix SNPs = ",nrow(new_LD),"; finemap_DT = ",nrow(finemap_DT), v=verbose)
  }
  printer("++ Subsetting complete.",v=verbose)
  return(list(LD=new_LD,
              DT=new_DT))
}




#' Remove all SNPs outside of of a given gene.
#'
#' Get the min/max coordinates of a given gene (including known regulatory regions, introns, and exons).
#' Remove any SNPs from the data.frame that fall outside these coordinates.
#'
#' @family SNP filters
#' @keywords internal
gene_trimmer <- function(subset_DT,
                         gene,
                         min_POS=NULL,
                         max_POS=NULL){
  printer("BIOMART:: Trimming data to only include SNPs within gene coordinates.")
  gene_info <- biomart_geneInfo(gene)
  gene_info_sub <- subset(gene_info, hgnc_symbol==gene)
  # Take most limiting min position
  min_POS <- max(min_POS, gene_info_sub$start_position, na.rm = T)
  # Take most limiting max position
  max_POS <- min(max_POS, gene_info_sub$end_position, na.rm = T)
  subset_DT <- subset(subset_DT, CHR==gene_info$chromosome_name[1] & POS>=min_POS & POS<=max_POS)
  printer("BIOMART::",nrow(subset_DT),"SNPs left after trimming.")
  return(subset_DT)
}


#' Limit the number of SNPs per locus.
#'
#' Start with the lead SNP and keep expanding the window until you reach the desired number of snps.
#' \code{subset_DT} should only contain one locus from one chromosome.
#'
#' @family SNP filters
#' @param max_snps The maximum number of SNPs to keep in the resulting data.frame.
#' @param subset_DT A data.frame that contains at least the following columns:
#' \describe{
#'   \item{SNP}{RSID for each SNP.}
#'   \item{POS}{Each SNP's genomic position (in basepairs).}
#' }
#' @keywords internal
limit_SNPs <- function(max_snps=500, subset_DT){
  printer("echolocator:: Limiting to only",max_snps,"SNPs.")
  if(nrow(subset_DT)>max_snps){
    orig_n <- nrow(subset_DT)
    lead.index <- which(subset_DT$leadSNP==T)
    i=1
    tmp.sub<-data.frame()
    while(nrow(tmp.sub)<max_snps-1){
      # print(i)
      snp.start <- max(1, lead.index-i)
      snp.end <- min(nrow(subset_DT), lead.index+i)
      tmp.sub <- subset_DT[snp.start:snp.end]
      i=i+1
    }
    # Need to add that last row on only one end
    snp.start <- max(1, lead.index-i+1) # +1 keep it the same index as before
    snp.end <- min(nrow(subset_DT), lead.index+i)
    tmp.sub <- subset_DT[snp.start:snp.end]

    printer("+ Reduced number of SNPs:",orig_n,"==>",nrow(tmp.sub))
    return(tmp.sub)
  } else {
    printer("+ Data already contains less SNPs than limit (",nrow(subset_DT),"<",max_snps,")")
    return(subset_DT)
  }
}




#' Filter SNPs
#'
#' Filter SNps by MAF, window size, min/max position, maxmimum number of SNPs, or gene coordinates.
#' You can also explicitly remove certain variants.
#'
#' @family SNP filters
#' @keywords internal
filter_snps <- function(subset_DT,
                        bp_distance,
                        remove_variants,
                        locus,
                        verbose=T,
                        min_POS=NULL,
                        max_POS=NULL,
                        max_snps=NULL,
                        min_MAF=NULL,
                        trim_gene_limits=F){
  if(remove_variants!=F){
    printer("Removing specified variants:",paste(remove_variants, collapse=','), v=verbose)
    try({subset_DT <- subset(subset_DT, !(SNP %in% remove_variants) )})
  }
  # Trim subset according to annotations of where the gene's limit are
  if(trim_gene_limits!=F){
    subset_DT <- gene_trimmer(subset_DT=subset_DT,
                              gene=trim_gene_limits,
                              min_POS=min_POS,
                              max_POS=min_POS)
  }
  if(!is.null(max_snps)){
    subset_DT <- limit_SNPs(max_snps = max_snps, subset_DT = subset_DT)
  }
  if(!is.null(min_MAF) & length(min_MAF>0)>0){
    printer("echolocatoR:: Removing SNPs w/ MAF <",min_MAF)
    subset_DT <- subset(subset_DT, MAF>=min_MAF)
  }
  # Limit range
  if(!is.null(bp_distance)){
    lead.snp <- subset(subset_DT, leadSNP)
    subset_DT <- subset(subset_DT,
                        POS >= lead.snp$POS - bp_distance &
                          POS <= lead.snp$POS + bp_distance)
  }
  if(!is.na(min_POS)){subset_DT <- subset(subset_DT, POS>=min_POS)}
  if(!is.na(max_POS)){subset_DT <- subset(subset_DT, POS<=max_POS)}
  printer("++ Post-filtered data:",paste(dim(subset_DT), collapse=" x "))
  return(subset_DT)
}

#' Identify SNPs to condition on.
#'
#' When running conditional analyses (e.g. \emph{GCTA-COJO}),
#' this functions automatically identifies SNP to condition on.
#'
#' @family SNP filters
snps_to_condition <- function(conditioned_snps, top_SNPs, loci){
  if(conditioned_snps=="auto"){
    lead_SNPs_DT <- subset(top_SNPs, Locus %in% loci)
    # Reorder
    lead_SNPs_DT[order(factor(lead_SNPs_DT$Locus,levels= loci)),]
    return(lead_SNPs_DT$SNP)
  } else {return(conditioned_snps)}
}


# ---------------


#' @section directory functions:
#' Functions for storing and extracting info about summary statistic dataset paths and metadata attributes.


#' @family directory functions
directory_info <- function(info_path=NULL,
                           dataset_name,
                           variable="fullSS.local"){
  Data_dirs <- data.table::fread(info_path)
  directory = subset(Data_dirs, Dataset==dataset_name, select=variable) %>% as.character()
  return(directory)
}




#' @family directory functions
#' @keywords internal
get_dataset_name <- function(file_path){
  dataset_name <- tail(strsplit(dirname(file_path), "/")[[1]],n = 1)
  return(dataset_name)
}




#' @family directory functions
#' @keywords internal
delete_subset <- function (force_new_subset, subset_path){
  # Force new file to be made
  if(force_new_subset==T){
    printer("\n + Removing existing summary stats subset...\n")
    # dataset_name <- get_dataset_name(subset_path)
    # subset_path <- paste(dirname(subset_path),"/",gene,"_",superpopulation,"_",dataset_name,"_subset.txt",sep="")
    suppressWarnings(file.remove(subset_path))
  }
}



#' Make paths for results and subsets
#'
#' @family directory functions
#' @keywords internal
make_results_path <- function(dataset_name,
                              dataset_type,
                              locus){
  results_path <- file.path("Data",dataset_type, dataset_name, locus)
  dir.create(results_path, showWarnings = F, recursive = T)
  return(results_path)
}

#' @family directory functions
#' @keywords internal
.get_subset_path <- function(results_path,
                            locus,
                            subset_path="auto",
                            suffix=".tsv.gz"){
  # Specify subset file name
  if(subset_path=="auto"){
    dataset_name <- basename(dirname(results_path))
    created_sub_path <- file.path(results_path, paste0(locus,"_",dataset_name,"_subset",suffix) )
    return(created_sub_path)
  } else{return(subset_path)}
}


