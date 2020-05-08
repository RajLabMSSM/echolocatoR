#  ((((((((((((((((()))))))))))))))))  #
 # ((((((((((((((((())))))))))))))))) #
 #    ------- echolocatoR -------    #
 # ((((((((((((((((())))))))))))))))) #
#              ((((()))))              #

# Author: Brian M. Schilder
# Bioinformatician II
# Icahn School of Medicine at Mount Sinai
# New York City, New York, USA
# https://bschilder.github.io/BMSchilder

     # =/\                  /\=
    #  / \'._   (\_/)   _.'/ \
   #  / .''._'--(o.o)--'_.''. \
  #  /.' _/ |`'=/ " \='`| \_ `.\
 #  /` .' `\;-,'\___/',-;/` '. '\
#  /.-'       `\(-V-)/`       `-.\
# `              v  v           `
#
#   message("      =/\                  /\=    ")
#   message("      / \'._   (\_/)   _.'/ \     ")
#   message("     / .''._'--(o.o)--'_.''. \    ")
#   message("    /.' _/ |`'=/ { \='`| \_ `.\   ")
#   message("   /` .' `\;-,'\___/',-;/` '. '\  ")
#   message("  /.-'       `\(-V-)/`       `-.\ ")
#   message(" `              v  v           `  ")
#





## susieR Function
#
# * Notes on L parameter
# + L is the expected number of causal variants
# + Increasing L increases computational time
# + L=1: Gives a good amount of variation in PIP.
# + L=2: Warns "IBSS algorithm did not converge in 100 iterations!", but gives good variation in PIP.
# + L=3: Warns "IBSS algorithm did not converge in 100 iterations!". All PIPs 1s and 0s.
# + These results seem to be at least partially dependent on whether the ethnic composition of the LD matrix.
# * Notes on variance:
#   + If 'estimate_residual_variance' = TRUE _without_ providing 'var_y' _and_ L>1, susieR will throw error:
#   __"Estimating residual variance failed: the estimated value is negative"__
# + Running susieR with 'var_y = var(b)' provides _exactly_ the same results.
# * Statistical Terms:
#   + posterior inclusion probability (PIP)
# + coefficient estimate (Beta)
# + Effect allele frequency (EAF)
# + The I^2 statistic describes the percentage of variation across studies that seems not to be due to chance.






# R package creation:
## http://r-pkgs.had.co.nz
# R package remote dependencies:
## https://cran.r-project.org/web/packages/devtools/vignettes/dependencies.html


# You can learn more about package authoring with RStudio at:
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# regex commands: http://www.endmemo.com/program/R/gsub.php

# Load libraries
# .libPaths()

# Add libraries to Imports secion of DESCRIPTION file:

# CRAN Imports:
# for(p in c("dplyr","tidyr","data.table","ggplot2","patchwork",
#            "pbmcapply","ggrepel","BiocManager",
#            "coloc","haploR","XGR","reticulate")){
#   usethis::use_package(package = p, type = "Imports")
# }
# # CRAN Suggests:
# for(p in c("crayon")){
#   usethis::use_package(package = p, type = "Suggests")
# }
# # Bioconductor Imports:
# for(p in c("ggbio","XGR",
#            "rtracklayer")){
#   usethis::use_package(package = p, type = "Imports")
# }
# # GiHub Imports:
# for(p in c("stephenslab/susieR")){
#   remotes::
#   usethis::use_package(package = p, type = "Imports")
# }


# library(R.utils)
# library(devtools)
# library(readxl)
# library(DT)
# library(data.table)
# library(dplyr)
# library(ggplot2)
# library(plotly)
# library(patchwork) #devtools::install_github("thomasp85/patchwork")

# library(pbmcapply); #devtools::install_github("kvnkuang/pbmcapply", ref = "dev")

# library(cowplot)
# library(ggrepel)
# library(curl)
# library(gaston)
# library(tidyr)
# library(BiocManager)
# library(biomaRt) # BiocManager::install("biomaRt")
# library(refGenome)
# library(crayon)
# library(ChIPQC); BiocManager::install('ChIPQC')

# library(coloc)
# install.packages("haploR", dependencies = TRUE)
# library(haploR)
# library(GeneOverlap) #BiocManager::install("GeneOverlap")
# library(rtracklayer) #BiocManager::install("rtracklayer")
# BiocManager::install(c("supraHex","graph","Rgraphviz","dnet"))
# library(XGR)# install.packages("XGR")



# library(Rsamtools) # BiocManager::install("Rsamtools")

# install_fGWAS <- function(){
#   devtools::install_github("wzhy2000/fGWAS/pkg")
#   system("git clone https://github.com/wzhy2000/fGWAS.git & cd fGWAS & R CMD INSTALL pkg")
# }
# library(fGWAS)
# library(snpStats) #BiocManager::install("snpStats")
# library(bigsnpr) # BiocManager::install("bigsnpr")

# *** SUSIE ****
# library(knitrBootstrap) #install_github('jimhester/knitrBootstrap')
# library(susieR) # devtools::install_github("stephenslab/susieR")

# *** finemapr ****
## finemapr contains: finemap, CAVIAR, and PAINTOR
# library(finemapr) # devtools::install_github("variani/finemapr")
# library(roxygen2) #roxygenize()

# *** locuscomparer ****
# https://github.com/boxiangliu/locuscomparer
# library(locuscomparer); #devtools::install_github("boxiangliu/locuscomparer")

# thm <- knitr::knit_theme$get("bipolar")
# knitr::knit_theme$set(thm)



#' @section \emph{echolocatoR} helper functions:
#' General use functions for various points in the \emph{echolocatoR} pipeline.



#' Ensembl IDs -> Gene symbols
#'
#' Convert HGNC gene symbols to Ensembl IDs.
#'
#' @param gene_symbols List of HGNC gene symbols.
#' @return List of ensembl IDs.
#' @examples
#' gene_symbols <- c("FOXP2","BDNF","DCX","GFAP")
#' ensembl_IDs <- hgnc_to_ensembl(gene_symbols)
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
effective_sample_size <- function(finemap_DT){
  finemap_DT$N <- (4.0 / (1.0/finemap_DT$N_cases + 1.0/finemap_DT$N_controls) )
  sample_size <- as.integer(median(finemap_DT$N))
  return(sample_size)
}




#' Assign a lead GWAS SNP to a locus
#'
#' If none of the SNPs in the data.frame have \code{leadSNP==T},
#' then sort by lowest p-value (and then highest Effect size) and assign the top SNP as the lead SNP.
#'
#' @param data.frame Fine-mapping results data.frame.
#' @return Fine-mapping results data.frame with new boolean \strong{leadSNP} column,
#'  indicating whether each SNPs is the lead GWAS SNP in that locus or not.
#' @examples
#' finemap_DT <- assign_lead_SNP(finemap_DT)
assign_lead_SNP <- function(new_DT, verbose=T){
  if(sum(new_DT$leadSNP)==0){
    printer("+ leadSNP missing. Assigning new one by min p-value.", v=verbose)
    top.snp <- head(arrange(new_DT, P, desc(Effect)))[1,]$SNP
    new_DT$leadSNP <- ifelse(new_DT$SNP==top.snp,T,F)
  }
  return(new_DT)
}




#' Subset LD matrix and data.frame to only their shared SNPs
#'
#' Find the SNPs that are shared between an LD matrix and another data.frame with a `SNP` column.
#' Then remove any non-shared SNPs from both objects.
#'
#' @param data.frame
#' @return data.frame
#' @example
#' finemap_DT <- assign_lead_SNP(finemap_DT)
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
#' @examples
#' subset_DT <- gene_trimmer(subset_DT, gene, min_POS=NULL, max_POS=NULL)
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
#' @param max_snps The maximum number of SNPs to keep in the resulting data.frame.
#' @param subset_DT A data.frame that contains at least the following columns:
#' \describe{
#'   \item{SNP}{RSID for each SNP.}
#'   \item{POS}{Each SNP's genomic position (in basepairs).}
#' }
#' @examples
#' subset_DT <- limit_SNPs(max_snps=500, subset_DT)
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
#' @examples
#' echoR.filter_snps <- function(subset_DT, bp_distance, remove_variants="res, gene, min_POS=10000060, max_POS=1000260, max_snps=NULL, min_MAF=0.01, trim_gene_limits=F)
echoR.filter_snps <- function(subset_DT,
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




arg_list_handler <- function(arg, i){
  output <- if(length(arg)>1){arg[i]}else{arg}
  return(output)
}




snps_to_condition <- function(conditioned_snps, top_SNPs, loci){
  if(conditioned_snps=="auto"){
    lead_SNPs_DT <- subset(top_SNPs, Locus %in% loci)
    # Reorder
    lead_SNPs_DT[order(factor(lead_SNPs_DT$Locus,levels= loci)),]
    return(lead_SNPs_DT$SNP)
  } else {return(conditioned_snps)}
}






#' @section \emph{echolocatoR} main functions:
#' The primary functions of \emph{echolocatoR} that expedite fine-mapping
#'  by wrapping many other \emph{echolocatoR} functions into one.
#'  Encompasses steps including:
#'  \describe{
#'  \item{Subset & standardize}{Extract subsets of the full summmary stats GWAS or QTL file and reformat them to be compatible with \emph{echolocatoR}'s various functions }
#'  \item{Calculate linkage disequilibrium}{Download and prepare the necessary LD matrix.}
#'  \item{Fine-map}{Run various fine-mapping tools and merge the results into a single multi-finemap data.frame.}
#'  \item{Plot}{Summarise the results in a multi-track plot for each locus.}
#'  }



#' Run echolocatoR pipeline on a single locus.
#'
#' Unlike \code{finemap_loci}, you don't need to provide a \code{top_SNPs} data.frame.
#' Instead, just manually provide the coordinates of the locus you want to fine-map.
#'
#' @param loci Character list of loci in \strong{Locus} col of \code{top_SNPs}.
#' @param fullSS_path Path to the full summary statistics file (GWAS or QTL) that you want to fine-map.
#' @param dataset_name The name you want to assign to the dataset being fine-mapped,
#' This will be used to name the subdirectory where your results will be stored (e.g. \code{Data/GWAS/<dataset_name>}).
#' Don't use special characters (e.g.".", "/").
#' @param dataset_type The kind dataset you're fine-mapping (e.g. GWAS, eQTL, tQTL).
#' This will also be used when creating the subdirectory where your results will be stored (e.g. \code{Data/<dataset_type>/Kunkle_2019}).
#' @param top_SNPs A data.frame with the genomic coordinates of the lead SNP for each locus.
#' The lead SNP will be used as the center of the window when extracting subset from the full GWAS/QTL summary statistics file.
#' Only one SNP per \strong{Locus} should be included.
#' At minimum, \code{top_SNPs} should include the following columns:
#' \describe{
#' \item{Locus}{A unique name for each locus. Often, loci are named after a relevant gene (e.g. LRRK2) or based on the name/coordinates of the lead SNP (e.g. locus_chr12_40734202) }
#' \item{CHR}{The chromosome that the SNP is on. Can be "chr12" or "12" format.}
#' \item{POS}{The genomic position of the SNP (in basepairs)}
#' }
#' @param force_new_subset By default, if a subset of the full summary stats file for a given locus is already present,
#' then \emph{echolocatoR} will just use the preexisting file.
#' Set \code{force_new_subset=T} to override this and extract a new subset.
#' Subsets are saved in the following path structure: \url{Data/<dataset_type>/<dataset_name>/<locus>/Multi-finemap/<locus>_<dataset_name>_Multi-finemap.tsv.gz}
#' @param force_new_LD  By default, if an LD matrix file for a given locus is already present,
#' then \emph{echolocatoR} will just use the preexisting file.
#' Set \code{force_new_LD=T} to override this and extract a new subset.
#' @param force_new_finemap By default, if an fine-mapping results file for a given locus is already present,
#' then \emph{echolocatoR} will just use the preexisting file.
#' Set \code{force_new_finemap=T} to override this and re-run fine-mapping.
#' @param finemap_methods Which fine-mapping methods you want to use.
#' @param bp_distance The width of the window size you want each locus to be.
#' For example, if \code{bp_distance=500000} then the locus will span 500kb from the lead SNP in either direction,
#' resulting in a locus that is ~1Mb long (depending on the dataset).
#' @param n_causal The maximum number of potential causal SNPs per locus.
#' This parameter is used somewhat differntly by different fine-mapping tools.
#' See tool-specific functions for details.
#' @param sample_size The overall sample size of the study.
#' If none is given, and \strong{N_cases} and \sontrg{N_controls} columns are present,
#' then sample_size is inferred to be \code{max(N_cases) + max(N_controls)}.
#' @param chrom_col Name of the chromosome column in the full summary stats file.
#' @param position_col Name of the genomic position column in the full summary stats file.
#' @param snp_col Name of the SNP RSID column in the full summary stats file.
#' @param pval_col Name of the p-value column in the full summary stats file.
#' Raw p-values are preferred, but if not available corrected p-values (e.g. FDR) can be used instead.
#' @param effect_col Name of the effect size column in the full summary stats file.
#' Effect size is preferred, but if not available other metrics like Beta for Odds Ratio can be used instead.
#' @param stderr_col Name of the standard error  column in the full summary stats file.
#' You can also set \code{stderr_col="calculate"} to infer standard error using: \code{effect / tstat}.
#' @param tstat_col Name of the t-statistic column in the full summary stats file.
#' This column is not necessary unless \code{stderr_col="calculate"} or the standard error column is missing.
#' @param locus_col Name of the locus column in the full summary stats file.
#' @param freq_col Name of the allele frequency column in the full summary stats file.
#' Effect allele frequency is preferred, but the non-effect allele can be provided instead (though this may be less accurate).
#' This column is not necessary unless \code{MAF_col="calculate"} or the MAF column is missing.
#' @param MAF_col Name of the minor allele frequency column in the full summary stats file.
#' Can be inferred from \strong{freq_col} if missing from the dataset.
#' @param A1_col Name of the effect/risk allele column in the full summary stats.
#' Unfortunately, different studies report different kinds of allele information in a non-standardized way.
#' Meaning that A1/A2 can refer to any number of things:
#'  \desribe{
#'  \item{risk/non-risk alleles}{in the case of diseases}
#'  \item{ref/alt alleles}{where ref is the reference genome being used}
#'  \item{major/minor alleles}}
#'  }
#'  This makes comparing summary stats across GWAS/QTL datasets very confusing for several reasons:
#'  \desribe{
#'  \item{Multi-allelic SNPs}{SNPs can have more than just 2 possible alleles (multi-allelic SNPs). Even if you compare the same SNP between two studies, you may accidentally be comparing totally different alleles.}
#'  \item{Valence}{The valence (+/-) of per-SNP GWAS effect sizes/beta can be relative to different allele types between studies.
#'  For example, let's say in one GWAS study your effect size for SNP A is 1.5 relative to the major allele in one study,
#'   and the minor allele happens to be the one found in the reference genome.
#'   You then try to compare that effect size to that of the same SNP in another GWAS.
#'   But, the valence of the effect sizes in the 2nd GWAS study are all relative to the reference genome (instead of the minor allele),
#'   giving the same SNP a value of -1.2. If you took the effect sizes at face value you'd say the signals are in opposite directions.
#'   But once you take into account how the valences were determined in each study you realize that they're actually both positive relative to the major allele.}
#'  }
#'  This process of reversing per-SNP valences based on aligning the alleles is known as allele flipping. This is important when comparing individual SNPs, but can also have an impact on colocalization results.
#'  @param N_cases_col Name of the column in the full summary stats that has the number of case subjects in the study.
#'  This can either be per SNP sample sizes, or one number repeated across all rows.
#'  Proxy cases (e.g. relatives of people with the disease being investigated) should be included in this estimate if any were used in the study.
#'  This column is not necesssary if \code{N_cases} parameter is provided.
#'  @param N_controls_col Name of the column in the full summary stats that has the number of control subjects in the study.
#'  This can either be per SNP sample sizes, or one number repeated across all rows.
#'  This column is not necesssary if \code{N_controls} parameter is provided.
#'   @param N_cases The number of case subjects in the study.
#'  Instead of providing a reudundant \strong{N_cases_col} column, you can simply enter one value here.
#'  @param N_controls The number of control subjects in the study.
#'  Instead of providing a reudundant \strong{N_controls_col} column, you can simply enter one value here.
#'  @param proportion_cases The proportion of total subjects in the study that were cases.
#'  if \code{proportion_cases="calculate"} then this is inferred:  \code{N_controls / N_controls}.
#'  @param LD_reference Which linkage disequilibrium reference panel do you want to use.
#'  Options include:
#'  \describe{
#'  \item{"UKB"}{A pre-caclulated LD reference matrix from a subset of caucasian British individuals from the UK Biobank. See \link[https://www.biorxiv.org/content/10.1101/807792v2]{Wiessbrod et al. (2019)} for more details.}
#'  \item{"1KG_Phase1"}{Download a subset of the 1000 Genomes Project Phase 1 vcf and calculate LD on the fly with plink.}
#'  \item{"1KG_Phase3"}{Download a subset of the 1000 Genomes Project Phase 3 vcf and calculate LD on the fly with plink.}
#'  }
#'  @param superpopulation Subset your LD reference panel by superopulation.
#'  Setting the superpopulation is not currently possible when \code{LD_reference="UKB"}.
#'  \link[https://www.internationalgenome.org/faq/which-populations-are-part-your-study/]{1KG_Phase1 options} include:
#'  \describe{
#'  \item{"AFR"}{African [descent]}
#'  \item{"AMR"}{Ad-mixed American}
#'  \item{"EAS"}{East Asian}
#'  \item{"EUR"}{European}
#'  \item{"SAS"}{South Asian}
#'  }
#' @param download_reference When acquiring LD matrixes,
#'  the default is to delete the full vcf or npz files after \emph{echolocator} has extracted the necssary subset.
#'  However, if you wish to keep these full files (which can be quite large) set \code{download_reference=T}.
#' @param min_POS Manually set the minimum genomic position for your locus subset.
#' \code{min_POS} can clip the window size set by \code{bp_distance}.
#' Only use this parameter when fine-mapping one locus at a time.
#' @param max_POS Manually set the maximum genomic position for your locus subset.
#' \code{max_POS} can clip the window size set by \code{bp_distance}.
#' Only use this parameter when fine-mapping one locus at a time.
#' @param min_MAF Remove any SNPs with \strong{MAF} < \code{min_MAF}.
#' @param trim_gene_limits If a valid gene symbol is provided to \code{trim_gene_limits},
#' the gene's canonical coordinates are pulled from \code{biomaRt}.
#' This includes introns, exons, and proximal regulatory regions (e.g. promoters).
#' Any SNPs that fall outside these coordinates are remove from downstream fine-mapping.
#' Set \code{trim_gene_limits=F} to not limit by gene coordinates (\emph{default}).
#' @param max_snps The maximum number of SNPs to include in the locus.
#' If the current window size yields > \code{max_snps},
#'  then the outer edges of the of the locus are trimmed until the number of SNPs ≤ \code{max_snps}.
#' @param file_sep The separator in the full summary stats file.
#' This parameter is only necessary if \code{query_by!="tabix"}.
#' @param min_r2 Remove any SNPs are below the LD r2 threshold with the lead SNP within their respective locus.
#' @param LD_block Calculate LD blocks with \emph{plink} and only include the block to which the lead SNP belongs.
#' @param block_size Adjust the granularity of block sizes when \code{LD_block=T}.
#' @param min_Dprime Remove any SNPs are below the LD D' threshold with the lead SNP within their respective locus.
#' This is paramter currently only works when \code{LD_reference!="UKB"}.
#' @param query_by Choose which method you want to use to extract locus subsets from the full summary stats file.
#' Methods include:
#' \describe{
#' \item{"tabix"}{Convert the full summary stats file in an indexed tabix file. Makes querying lightning fast after the initial conversion is done. (\emph{default})}
#' \item{"coordinates"}{Extract locus subsets using min/max genomic coordinates with \emph{awk}.}
#' }
#' @param remove_variants A list of variants to remove from the locus subset file.
#' @param remove_correlates If \code{remove_correlates} is set to a value between 0-1,
#' removes any SNPs that are in LD with any of the \code{remove_variants} above the threshold provided by \code{remove_correlates}.
#' @param probe_path The location of the file containing translations between probe IDs and gene symbols.
#' Only used for certain eQTL datasets.
#' @param conditioned_snps Which SNPs to conditions on when fine-mapping with \emph{COJO}.
#' @param plot_LD Whether to plot a subset of the LD matix.
#' @param verbose Whether \emph{echolocatoR} should be verbose or silent.
#' @param remove_tmps Whether to remove any temporary files (e.g. FINEMAP output files) after the pipeline is done running.
#' @param plot_types Which kinds of plots to include.
#' Options:
#' \describe{
#' \item{"simple"}{Just plot the following tracks: GWAS, fine-mapping, gene models, and brain cell type-specific epigenomic data from Nott et al. (2019).}
#' \item{"fancy"}{Additionally plot XGR annotation tracks.}
#' }
#' @param PAINTOR_QTL_datasets A list of QTL datasets to be used when conducting joint functional fine-mapping with \emph{PAINTOR}.
#' @param server
#' Whether \emph{echolocatoR} is being run on a computing cluster/server or on a local machine.
#' @param PP_threshold The minimum fine-mapped posterior probability for a SNP to be considered part of a Credible Set.
#' For example, \code{PP_threshold=.95} means that all Credible Set SNPs will be 95\% Credible Set SNPs.
#' @param plot_window Zoom into the center of the locus when plotting (without editing the fine-mapping results file).
#' @param plot_Nott_binwidth When including Nott et al. (2019) epigenomic data in the track plots,
#' adjust the bin width of the histograms.
#' @param Nott_bigwig_dir Instead of pulling Nott et al. (2019) epigenomic data from the UCSC Genome Browsers, use a set of local bigwig files.
finemap_pipeline <- function(locus,
                             fullSS_path,
                             dataset_name,
                             dataset_type="GWAS",
                             top_SNPs="auto",
                             force_new_subset=F,
                             force_new_LD=F,
                             force_new_finemap=T,
                             finemap_methods=c("ABF","SUSIE","FINEMAP"),
                             bp_distance=500000,
                             n_causal=5,
                             sample_size=NA,
                             chrom_col="CHR",
                             position_col="POS",
                             snp_col="SNP",
                             pval_col="P",
                             effect_col="Effect",
                             stderr_col="StdErr",
                             tstat_col="t-stat",
                             locus_col="Locus",
                             freq_col="Freq",
                             MAF_col="MAF",
                             A1_col = "A1",
                             A2_col = "A2",
                             N_cases_col="N_cases",
                             N_controls_col="N_controls",
                             N_cases=NULL,
                             N_controls=NULL,
                             proportion_cases="calculate",

                             LD_reference="1KG_Phase1",
                             superpopulation="EUR",
                             download_reference=T,
                             min_POS=NA,
                             max_POS=NA,
                             min_MAF=NA,
                             trim_gene_limits=F,
                             max_snps=NULL,

                             file_sep="\t",
                             min_r2=0,
                             LD_block=F,
                             block_size=.7,
                             min_Dprime=F,
                             query_by="coordinates",
                             remove_variants=F,
                             remove_correlates=F,
                             probe_path = "./Data/eQTL/gene.ILMN.map",
                             conditioned_snps,
                             plot_LD = F,
                             verbose=T,
                             remove_tmps=T,
                             plot_types=c("simple"),
                             PAINTOR_QTL_datasets=NULL,
                             server=F,
                             PP_threshold=.95,
                             plot_window=NULL,
                             plot_Nott_binwidth=2500,
                             Nott_bigwig_dir=NULL){
   # Create paths
   results_path <- make_results_path(dataset_name, dataset_type, locus)
   # Extract subset
   subset_DT <- extract_SNP_subset(gene = locus,
                                    top_SNPs = top_SNPs,
                                    fullSS_path = fullSS_path,
                                    results_path  =  results_path,
                                    force_new_subset = force_new_subset,

                                    chrom_col = chrom_col,
                                    position_col = position_col,
                                    snp_col = snp_col,
                                    pval_col = pval_col,
                                    effect_col = effect_col,
                                    stderr_col = stderr_col,
                                    locus_col = locus_col,
                                    tstat_col = tstat_col,
                                    MAF_col = MAF_col,
                                    freq_col = freq_col,
                                    A1_col = A1_col,
                                    A2_col = A2_col,

                                    N_cases_col = N_cases_col,
                                    N_controls_col = N_controls_col,
                                    N_cases = N_cases,
                                    N_controls = N_controls,
                                    proportion_cases = proportion_cases,

                                    bp_distance = bp_distance,
                                    superpopulation = superpopulation,
                                    min_POS = min_POS,
                                    max_POS = max_POS,
                                    file_sep = file_sep,
                                    query_by = query_by,
                                    probe_path = probe_path,
                                    remove_tmps = remove_tmps)

  ### Compute LD matrix
  message("--- Step 2: Calculate Linkage Disequilibrium ---")
  LD_matrix <- LD.load_or_create(results_path=results_path,
                                 subset_DT=subset_DT,
                                 gene=locus,
                                 force_new_LD=force_new_LD,
                                 LD_reference=LD_reference,
                                 superpopulation=superpopulation,
                                 download_reference=download_reference,
                                 min_r2=min_r2,
                                 LD_block=LD_block,
                                 block_size=block_size,
                                 min_Dprime=min_Dprime,
                                 remove_correlates=remove_correlates,
                                 verbose=verbose,
                                 server=server,
                                 remove_tmps=remove_tmps)

  #### ***** SNP Filters ***** ###
  # Remove pre-specified SNPs
  ## Do this step AFTER saving the LD to disk so that it's easier to re-subset in different ways later without having to redownload LD.
  message("-------------- Step 3: Filter SNPs -------------")
  subset_DT <- echoR.filter_snps(subset_DT=subset_DT,
                                  bp_distance=bp_distance,
                                  remove_variants=remove_variants,
                                  locus=locus,
                                  verbose=verbose,
                                  min_POS=min_POS,
                                  max_POS=max_POS,
                                  max_snps=max_snps,
                                  trim_gene_limits=trim_gene_limits)
  # Subset LD and df to only overlapping SNPs
  sub.out <- subset_common_snps(LD_matrix, subset_DT)
  LD_matrix <- sub.out$LD
  subset_DT <- sub.out$DT
  # finemap
  finemap_DT <- finemap_handler(results_path = results_path,
                                fullSS_path = fullSS_path,
                                finemap_methods = finemap_methods,
                                force_new_finemap = force_new_finemap,
                                dataset_type = dataset_type,
                                subset_DT = subset_DT,
                                LD_matrix = LD_matrix,
                                n_causal = n_causal,
                                sample_size = sample_size,
                                conditioned_snps = conditioned_snps,

                                snp_col = snp_col,
                                freq_col = freq_col,
                                effect_col = effect_col,
                                stderr_col = stderr_col,
                                pval_col = pval_col,
                                N_cases_col = N_cases_col,
                                N_controls_col = N_controls_col,
                                A1_col = A1_col,
                                A2_col = A2_col,
                                PAINTOR_QTL_datasets = PAINTOR_QTL_datasets,
                                PP_threshold = PP_threshold)
  finemap_DT <- find_consensus_SNPs(finemap_DT, credset_thresh = PP_threshold)
  # Step 6: COLOCALIZE
  # Step 7: Functionally Fine-map

  # Plot
  message("--------------- Step 7: Visualize --------------")
  if("simple" %in% plot_types){
    try({
      mf_plot <- GGBIO.plot(finemap_DT=finemap_DT,
                            LD_matrix=LD_matrix,
                            gene=locus,
                            results_path=results_path,
                            method_list=finemap_methods,
                            Nott_sn_epigenome = T,
                            XGR_libnames = NULL,
                            max_transcripts = 1,
                            plot_window = plot_window,
                            save_plot = T,
                            show_plot = T,
                            plot_Nott_binwidth = plot_Nott_binwidth,
                            Nott_bigwig_dir = Nott_bigwig_dir)
    })
  }
  if("fancy" %in% plot_types){
    try({
      trx <- GGBIO.plot(finemap_DT=finemap_DT,
                        LD_matrix=LD_matrix,
                        gene=locus,
                        results_path=results_path,
                        method_list=finemap_methods,
                        Nott_sn_epigenome = T,
                        XGR_libnames = c("ENCODE_TFBS_ClusteredV3_CellTypes",
                                         "ENCODE_DNaseI_ClusteredV3_CellTypes",
                                         "Broad_Histone"),
                        ROADMAP = T,
                        max_transcripts = 1,
                        plot_window = plot_window,
                        save_plot = T,
                        show_plot = T,
                        plot_Nott_binwidth = plot_Nott_binwidth,
                        Nott_bigwig_dir = Nott_bigwig_dir)
    })
  }

  # Plot LD
  if(plot_LD){
    try({
      LD_plot(LD_matrix=LD_matrix, subset_DT=subset_DT, span=10)
    })
  }




  # Show results table
  printer("+",locus,"Credible Set SNPs", v=verbose)
  print(createDT_html( subset(finemap_DT, Support >0) ))


  # Cleanup:
  if(remove_tmps){
    tmp_files <- file.path(results_path,"plink",
                           c("plink.bed",
                             "plink.bim",
                             "plink.fam",
                             "plink.ld",
                             "plink.ld.bin",
                             "plink.log",
                             "plink.nosex",
                             "SNPs.txt") )
    suppressWarnings(file.remove(tmp_files))
  }
  return(finemap_DT)
}




#' Fine-map multiple loci.
#'
#' \emph{echolocatoR} will separately fine-map each locus included in the \code{top_SNPs} data.frame.
#'
#' @param loci The list of loci you want to fine-map
#' @param subset_path The file you want your locus subset saved as.
#' Only use when fine-mapping one locus at a time.
#' If \code{subset_path="auto"} (\emph{default}), a locus subset file name is automatically constructed as:
#' \url{Data/<dataset_type>/<dataset_name>/<locus>/Multi-finemap/<locus>_<dataset_name>_Multi-finemap.tsv.gz}
#' @inheritParams finemap_pipeline
#' @return A merged data.frame with all fine-mapping results from all loci.
finemap_loci <- function(loci,
                         fullSS_path,
                         dataset_name,
                         dataset_type="GWAS",
                         force_new_subset=F,
                         force_new_LD=F,
                         force_new_finemap=T,
                         subset_path="auto",
                         top_SNPs="auto",
                         finemap_methods=c("ABF","SUSIE","FINEMAP"),
                         bp_distance=500000,
                         n_causal=5,
                         sample_size=NA,
                         chrom_col="CHR",
                         position_col="POS",
                         snp_col="SNP",
                         pval_col="P",
                         effect_col="Effect",
                         stderr_col="StdErr",
                         tstat_col = "t-stat",
                         MAF_col="MAF",
                         locus_col="Locus",
                         freq_col="Freq",
                         A1_col = "A1",
                         A2_col = "A2",
                         N_cases_col="N_cases",
                         N_controls_col="N_controls",
                         N_cases=NULL,
                         N_controls=NULL,
                         proportion_cases="calculate",

                         LD_reference="1KG_Phase1",
                         superpopulation="EUR",
                         topVariants=3,
                         download_reference=T,
                         min_POS=NA,
                         max_POS=NA,
                         min_MAF=NA,
                         trim_gene_limits=F,
                         max_snps=NULL,
                         file_sep="\t",
                         min_r2=0, LD_block=F, block_size=.7, min_Dprime=F,
                         query_by="coordinates",
                         remove_variants=F,
                         remove_correlates=F,
                         probe_path = "./Data/eQTL/gene.ILMN.map",
                         conditioned_snps="auto",
                         plot_LD=F,
                         verbose=T,
                         remove_tmps=T,
                         plot_types = c("simple"),
                         PAINTOR_QTL_datasets=NULL,
                         server=F,
                         PP_threshold=.95,
                         plot_window=NULL,
                         plot_Nott_binwidth=2500,
                         Nott_bigwig_dir=NULL){
  fineMapped_topSNPs <- data.table()
  fineMapped_allResults <- data.table()
  lead_SNPs <- snps_to_condition(conditioned_snps, top_SNPs, loci)

  for (i in 1:length(loci)){
    start_gene <- Sys.time()
    try({
      locus <- loci[i]
      message("🦇 🦇 🦇 ",locus," 🦇 🦇 🦇 ")
      lead_SNP <- arg_list_handler(lead_SNPs, i)
      gene_limits <- arg_list_handler(trim_gene_limits, i)
      conditioned_snp <- arg_list_handler(conditioned_snps, i)
      # message("^^^^^^^^^ Running echolocatoR on: ",locus," ^^^^^^^^^")
      # cat('  \n###', locus, '  \n')
      # Delete the old subset if force_new_subset == T
      finemap_DT <- finemap_pipeline(locus=locus,
                                     top_SNPs=top_SNPs,
                                     fullSS_path=fullSS_path,
                                     finemap_methods=finemap_methods,
                                     force_new_subset=force_new_subset,
                                     force_new_LD=force_new_LD,
                                     force_new_finemap=force_new_finemap,
                                     dataset_name=dataset_name,
                                     dataset_type=dataset_type,
                                     n_causal=n_causal,
                                     bp_distance=bp_distance,

                                     chrom_col=chrom_col,
                                     position_col=position_col,
                                     snp_col=snp_col,
                                     pval_col=pval_col,
                                     effect_col=effect_col,
                                     stderr_col=stderr_col,
                                     tstat_col=tstat_col,
                                     locus_col=locus_col,
                                     MAF_col=MAF_col,
                                     freq_col=freq_col,
                                     A1_col = A1_col,
                                     A2_col = A2_col,
                                     N_cases_col = N_cases_col,
                                     N_controls_col = N_controls_col,
                                     N_cases = N_cases,
                                     N_controls = N_controls,
                                     proportion_cases = proportion_cases,

                                     LD_reference=LD_reference,
                                     superpopulation=superpopulation,
                                     min_POS=min_POS,
                                     max_POS=max_POS,
                                     min_MAF=min_MAF,

                                     trim_gene_limits=gene_limits,
                                     max_snps=max_snps,
                                     file_sep=file_sep,
                                     min_r2=min_r2,
                                     LD_block=LD_block,
                                     block_size=block_size,
                                     min_Dprime=min_Dprime,
                                     query_by=query_by,
                                     remove_variants=remove_variants,
                                     remove_correlates=remove_correlates,
                                     probe_path=probe_path,
                                     conditioned_snps=lead_SNP,
                                     plot_LD=plot_LD,
                                     remove_tmps=remove_tmps,
                                     plot_types=plot_types,
                                     PAINTOR_QTL_datasets=PAINTOR_QTL_datasets,
                                     server=server,
                                     PP_threshold=PP_threshold,
                                     plot_window=plot_window,
                                     plot_Nott_binwidth=plot_Nott_binwidth,
                                     Nott_bigwig_dir=Nott_bigwig_dir)

      # Create summary table for all loci
      printer("Generating summary table...", v=verbose)
      newEntry <- cbind(data.table(Locus=locus), finemap_DT) %>% as.data.table()
      fineMapped_allResults <- rbind(fineMapped_allResults, newEntry, fill=T)
      cat('  \n')
    })
    end_gene <- Sys.time()
    message(gene,"fine-mapped in", round(end_gene-start_gene, 2),"seconds", v=verbose)

  }
  return(fineMapped_allResults)
}

