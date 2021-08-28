# %%%%%%%%%%%%%%%%% #
####### LD #######
# %%%%%%%%%%%%%%%%% #


# * Nalls et al. 2018: Imputation Panel Notes
# + _"One of the limitations of this study is the use of multiple imputation panels, due to logistic constraints.
# Adding datasets from non-European populations would be helpful to further improve our granularity in association
# testing and ability to fine-map loci through integration of more variable LD signatures."_
# + _"Post-Chang 23andMe samples were imputed using a combination of Finch for phasing (an in-house developed fork of Beagle)
# and miniMac2 for imputation with all-ethnicity samples from the __September 2013 release of 1000 Genomes Phase1__
# as reference haplotypes."_
# + _"The Nalls et al . 2014 and Chang et al . 2017 samples were imputed with Minimac2 using
# __1000 Genomes phase 1 haplotypes__.
# All additional sample series except for the post-Chang et al . 2017 samples from 23andMe were imputed using the
# __Haplotype Reference Consortium (HRC)__  on the University of Michigan imputation server under default settings
# with Eagle v2.3 phasing based on reference panel HRC r1.1 2016"_




#' Procure an LD matrix for fine-mapping
#'
#' Calculate and/or query linkage disequilibrium (LD) from reference panels (UK Biobank, 1000 Genomes),
#' or user-supplied datasets.
#'
#' Options:
#' \itemize{
#' \item Download pre-computed LD matrix from UK Biobank.
#' \item Download raw vcf file from 1KG and compute LD on the fly.
#' \item Compute LD on the fly from a user-supplied vcf file.
#' \item Use a user-supplied pre-computed LD-matrix.
#' }
#'
#' @param subset_DT The locus subset of the full summary stats file.
#' @inheritParams finemap_pipeline
#' @return A symmetric LD matrix of pairwise \emph{r} values.
#' @family LD
#' @keywords internal
#' @examples
#' \dontrun{
#' data("BST1"); data("locus_dir");
#' locus_dir <- file.path("~/Desktop",locus_dir)
#'  # BST1 <- limit_SNPs(500, BST1)
#'
#' # UK Biobank LD
#' LD_matrix <- LD.load_or_create(locus_dir=locus_dir, subset_DT=BST1, LD_reference="UKB")
#'
#' # 1000 Genomes
#' LD_matrix <- LD.load_or_create(locus_dir=locus_dir, subset_DT=BST1, LD_reference="1KGphase1", force_new_LD=T)
#'
#' # Local vcf file
#' LD_reference="~/Desktop/results/GWAS/Nalls23andMe_2019/BST1/LD/BST1.1KGphase1.vcf.gz"
#' LD_matrix <- LD.load_or_create(locus_dir=locus_dir, subset_DT=BST1, LD_reference=LD_reference, force_new_LD=T)
#' }
LD.load_or_create <- function(locus_dir,
                              subset_DT,
                              force_new_LD=F,
                              LD_reference="1KGphase1",
                              LD_genome_build="hg19",
                              superpopulation="EUR",
                              remote_LD=T,
                              download_method="direct",
                              vcf_folder=NULL,
                              # min_r2=0,
                              LD_block=F,
                              LD_block_size=.7,
                              # min_Dprime=F,
                              remove_correlates=F,
                              fillNA=0,
                              verbose=T,
                              server=F,
                              remove_tmps=T,
                              conda_env="echoR",
                              nThread=4){
  RDS_path <- LD.get_rds_path(locus_dir=locus_dir,
                              LD_reference=LD_reference)

  if(file.exists(RDS_path) & force_new_LD==F){
    #### Import existing LD ####
    printer("+  LD:: Previously computed LD_matrix detected. Importing...", RDS_path, v=verbose)
    LD_matrix <- readSparse(LD_path = RDS_path,
                            convert_to_df = F)
    LD_list <- list(DT=subset_DT,
                    LD=LD_matrix,
                    RDS_path=RDS_path)
  } else if(LD_reference=="UKB"){
    #### UK Biobank ####
    LD_list <- LD.UKBiobank(subset_DT = subset_DT,
                            locus_dir = locus_dir,
                            force_new_LD = force_new_LD,
                            chimera = server,
                            download_method = download_method,
                            fillNA = fillNA,
                            nThread = nThread,
                            return_matrix = T,
                            conda_env = conda_env,
                            remove_tmps = remove_tmps)
  } else if (LD_reference == "1KGphase1" |
             LD_reference == "1KGphase3") {
    #### 1000 Genomes ####
    LD_list <- LD.1KG(locus_dir = locus_dir,
                      subset_DT = subset_DT,
                      vcf_folder = vcf_folder,
                      LD_reference = LD_reference,
                      superpopulation = superpopulation,
                      remote_LD = remote_LD,

                      LD_block = LD_block,
                      LD_block_size = LD_block_size,
                      # min_Dprime = min_Dprime,
                      remove_correlates = remove_correlates,
                      fillNA = fillNA,
                      nThread = nThread,
                      conda_env = conda_env,
                      download_method = download_method)
  } else if (endsWith(tolower(LD_reference),".vcf") |
             endsWith(tolower(LD_reference),".vcf.gz")){
    #### Custom vcf ####
    LD_list <- LD.custom_panel(LD_reference=LD_reference,
                               LD_genome_build=LD_genome_build,
                               subset_DT=subset_DT,
                               locus_dir=locus_dir,
                               force_new_LD=force_new_LD,
                               # min_r2=min_r2,
                               # min_Dprime=min_Dprime,
                               # remove_correlates=remove_correlates,
                               fillNA=fillNA,
                               LD_block=LD_block,
                               LD_block_size=LD_block_size,
                               remove_tmps=remove_tmps,
                               nThread=nThread,
                               conda_env=conda_env,
                               verbose=verbose)
  } else {
    stop("LD:: LD_reference input not recognized. Please supply: '1KGphase1', '1KGphase3', 'UKB', or the path to a .vcf[.gz] file.")
  }
  return(LD_list)
}




LD.get_rds_path <- function(locus_dir,
                            LD_reference){
  RDS_path <- file.path(locus_dir,"LD",paste0(basename(locus_dir),".",basename(LD_reference),"_LD.RDS"))
  return(RDS_path)
}



#' Filter LD
#'
#' @family LD
#' @keywords internal
#' @examples
#' data("BST1"); data("LD_matrix");
#' LD_list <- list(LD=LD_matrix, DT=BST1)
#' LD_list <- LD.filter_LD(LD_list, min_r2=.2)
LD.filter_LD <- function(LD_list,
                         remove_correlates=F,
                         min_r2=0,
                         verbose=F){
  printer("+ FILTER:: Filtering by LD features.", v=verbose)
  subset_DT <- LD_list$DT
  LD_matrix <- LD_list$LD
  if(any(remove_correlates!=F)){
    # remove_correlates <- c("rs76904798"=.2, "rs10000737"=.8)
    for(snp in names(remove_correlates)){
      thresh <- remove_correlates[[snp]]
      printer("+ FILTER:: Removing correlates of",snp,"at r2 ≥",thresh,v=verbose)
      if(snp %in% row.names(LD_matrix)){
        correlates <- colnames(LD_matrix[snp,])[LD_matrix[snp,]>=sqrt(thresh)]
        LD_matrix <- LD_matrix[(!row.names(LD_matrix) %in% correlates),
                               (!colnames(LD_matrix) %in% correlates)]
      }
    }
  }
  if(min_r2!=0 & min_r2!=F){
    printer("+ FILTER:: Removing SNPs that don't correlate with lead SNP at r2 ≤",min_r2,v=verbose)
    LD_list <- subset_common_snps(LD_matrix = LD_matrix,
                                  finemap_dat = subset_DT,
                                  verbose = F)
    subset_DT <- LD_list$DT
    LD_matrix <- LD_list$LD
    lead.snp <- subset(subset_DT, leadSNP)$SNP[1]
    correlates <- colnames(LD_matrix[lead.snp,])[LD_matrix[lead.snp,]>=sqrt(min_r2)]
    LD_matrix <- LD_matrix[(row.names(LD_matrix) %in% correlates),
                           (colnames(LD_matrix) %in% correlates)]
  }
  LD_list <- subset_common_snps(LD_matrix = LD_matrix,
                                finemap_dat = subset_DT,
                                verbose = F)
  return(LD_list)
}




#' Compute LD from user-supplied vcf file
#'
#' @family LD
#' @keywords internal
#' @examples
#' \dontrun{
#' if(!"gaston" %in% row.names(installed.packages())){install.packages("gaston")}
#' data("BST1"); data("locus_dir")
#' LD_reference="~/Desktop/results/Reference/custom_panel_chr4.vcf"
#' LD_matrix <- LD.custom_panel(LD_reference=LD_reference, subset_DT=BST1, locus_dir=locus_dir)
#'
#' locus_dir <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/QTL/Microglia_all_regions/BIN1"
#' subset_DT <- data.table::fread(file.path(locus_dir,"Multi-finemap/BIN1.Microglia_all_regions.1KGphase3_LD.Multi-finemap.tsv.gz"))
#' LD_reference = "/sc/hydra/projects/pd-omics/glia_omics/eQTL/post_imputation_filtering/eur/filtered_variants/AllChr.hg38.sort.filt.dbsnp.snpeff.vcf.gz"
#' LD_matrix <- LD.custom_panel(LD_reference=LD_reference, subset_DT=BST1, locus_dir=locus_dir, LD_genome_build="hg38")
#' }
LD.custom_panel <- function(LD_reference,
                            fullSS_genome_build="hg19",
                            LD_genome_build="hg19",
                            subset_DT,
                            locus_dir,
                            force_new_LD=F,
                            min_r2=F,
                            # min_Dprime=F,
                            remove_correlates=F,
                            fillNA=0,
                            LD_block=F,
                            LD_block_size=.7,
                            remove_tmps=T,
                            nThread=4,
                            conda_env="echoR",
                            verbose=T){
  printer("LD:: Computing LD from local vcf file:",LD_reference)

  if(!toupper(LD_genome_build) %in% c("HG19","HG37","GRCH37")){
    printer("LD:: LD panel in hg38. Handling accordingly.",v=verbose)
    if(!toupper(fullSS_genome_build) %in% c("HG19","HG37","GRCH37")){
      ## If the query was originally in hg38,
      # that means it's already been lifted over to hg19.
      # So you can use the old stored POS.hg38 when the
      subset_DT <- subset_DT %>%
        dplyr::rename(POS.hg19=POS) %>%
        dplyr::rename(POS=POS.hg38)

    } else {
      ## If the query was originally in hg19,
      # that means no liftover was done.
      # So you need to lift it over now.
      subset_DT <- LIFTOVER(dat = subset_DT,
                            build.conversion = "hg19.to.hg38",
                            return_as_granges = F,
                            verbose = verbose)
    }
  }
  vcf_file <- LD.index_vcf(vcf_file=LD_reference,
                           force_new_index=F,
                           conda_env=conda_env,
                           verbose=verbose)
  # Make sure your query's chr format is the same as the vcf's chr format
  has_chr <- LD.determine_chrom_type_vcf(vcf_file = vcf_file)
  subset_DT <- dplyr::mutate(subset_DT,
                             CHR=if(has_chr) paste0("chr",gsub("chr","",CHR)) else gsub("chr","",CHR))

  vcf_subset <- LD.query_vcf(subset_DT=subset_DT,
                             locus_dir=locus_dir,
                             LD_reference=LD_reference,
                             vcf_URL=LD_reference,
                             whole_vcf=F,
                             remove_original_vcf=F,
                             force_new_vcf=force_new_LD,
                             query_by_regions=F,
                             nThread=nThread,
                             conda_env=conda_env,
                             verbose=verbose)

  bed_bim_fam <- LD.vcf_to_bed(vcf.gz.subset = vcf_subset,
                               locus_dir = locus_dir,
                               plink_prefix = "plink",
                               verbose =  verbose)
  # Calculate LD
  LD_matrix <- LD.snpstats_get_LD(LD_folder=file.path(locus_dir,"LD"),
                                  plink_prefix="plink",
                                  select.snps=unique(subset_DT$SNP),
                                  stats=c("R"),
                                  symmetric=T,
                                  depth="max",
                                  verbose=verbose)
  # Get MAF (if needed)
  subset_DT <- LD.snpstats_get_MAF(subset_DT=subset_DT,
                                   LD_folder=file.path(locus_dir,"LD"),
                                   plink_prefix="plink",
                                   force_new_MAF=F,
                                   verbose=verbose)
  # Filter out SNPs not in the same LD block as the lead SNP
  # Get lead SNP rsid
  leadSNP = subset(subset_DT, leadSNP==T)$SNP
  if(LD_block){
    block_snps <- LD.leadSNP_block(leadSNP = leadSNP,
                                   LD_folder = "./plink_tmp",
                                   LD_block_size = LD_block_size)
    LD_matrix <- LD_matrix[row.names(LD_matrix) %in% block_snps, colnames(LD_matrix) %in% block_snps]
    LD_matrix <- LD_matrix[block_snps, block_snps]
  }
  # IMPORTANT! Remove large data.ld file after you're done with it
  if(remove_tmps){
    suppressWarnings(file.remove(vcf_subset))
  }
  # Save LD matrix
  LD_list <- LD.save_LD_matrix(LD_matrix=LD_matrix,
                                subset_DT=subset_DT,
                                locus_dir=locus_dir,
                                fillNA = fillNA,
                                LD_reference=gsub(".vcf|.gz","",LD_reference),
                                sparse = T,
                                verbose=verbose)
  return(LD_list)
}




LD.determine_chrom_type_vcf <- function(vcf_file,
                                        conda_env="echoR",
                                        verbose=T){

  vcf <- gaston::read.vcf(vcf_file, max.snps = 1, convert.chr = F)
  has_chr <- grepl("chr",vcf@snps$chr[1])
  # bcftools <- CONDA.find_package(package = "bcftools",
  #                                conda_env = conda_env,
  #                                verbose = verbose)
  # # bcf_cmd <- paste("bcftools view -f '%CHROM' -H",vcf_file,"|head -1")
  # header <- data.table::fread(cmd=bcf_cmd)
  return(has_chr)
}




#' Save LD_matrix
#'
#' @family LD
#' @keywords internal
#' @examples
#' data("BST1"); data("LD_matrix"); data("locus_dir");
#' LD_list <- LD.save_LD_matrix(LD_matrix=LD_matrix, subset_DT=BST1, locus_dir=file.path("~/Desktop",locus_dir), LD_reference="UKB")
#' LD_list <- LD.save_LD_matrix(LD_matrix=LD_matrix, subset_DT=BST1, locus_dir=file.path("~/Desktop",locus_dir), LD_reference="custom_vcf")
LD.save_LD_matrix <- function(LD_matrix,
                              subset_DT,
                              locus_dir,
                              fillNA=0,
                              LD_reference,
                              subset_common=T,
                              sparse=T,
                              verbose=T){
  RDS_path <- LD.get_rds_path(locus_dir = locus_dir,
                              LD_reference = basename(LD_reference))
  printer("+ LD:: Saving LD matrix ==>",RDS_path,v=verbose)
  if(subset_common){
    sub.out <- subset_common_snps(LD_matrix = LD_matrix,
                                  fillNA = fillNA,
                                  finemap_dat = subset_DT,
                                  verbose = F)
    LD_matrix <- sub.out$LD
    subset_DT <- sub.out$DT
  }
  printer(dim(LD_matrix)[1],"x",dim(LD_matrix)[2],"LD_matrix",if(sparse) "(sparse)"else NULL, v=verbose)

  dir.create(dirname(RDS_path), showWarnings = F, recursive = T)
  if(sparse){
    saveSparse(LD_matrix = LD_matrix,
               LD_path = RDS_path,
               verbose = F)
  } else {
    saveRDS(LD_matrix, file = RDS_path)
  }
  return(list(LD=LD_matrix,
              DT=subset_DT,
              RDS_path=RDS_path))
}




#' Translate superopulation acronyms
#'
#' Ensures a common ontology for synonynmous superpopulation names.
#' @family LD
#' @keywords internal
#' @param superpopulation Three-letter superpopulation name.
LD.translate_population <- function(superpopulation){
  pop_dict <- list("AFA"="AFR", "CAU"="EUR", "HIS"="AMR",
                   "AFR"="AFR","EUR"="EUR", "AMR"="AMR")
  translated.list <- as.character(pop_dict[superpopulation])
  return(translated.list)
}




#' Plot a subset of the LD matrix
#'
#' Uses \code{gaston} to plot a SNP-annotated LD matrix.
#' @inheritParams finemap_pipeline
#' @family LD
#' @keywords internal
#' @param span This is very computationally intensive,
#' so you need to limit the number of SNPs with span.
#' If \code{span=10}, only 10 SNPs upstream and 10 SNPs downstream of the lead SNP will be plotted.
#' @examples
#' \dontrun{
#' data("BST1");
#' LD_matrix <- readRDS("/Volumes/Steelix/fine_mapping_files/GWAS/Nalls23andMe_2019/BST1/plink/UKB_LD.RDS")
#' LD.plot_LD(LD_matrix=LD_matrix, subset_DT=BST1)
#' }
LD.plot_LD <- function(LD_matrix,
                       subset_DT,
                       span=10,
                       method=c("gaston","heatmap","image")){
  leadSNP = subset(subset_DT, leadSNP==T)$SNP
  lead_index = match(leadSNP, row.names(LD_matrix))

  start_pos = lead_index - min(span, dim(LD_matrix)[1],na.rm = T)
  end_pos = lead_index +  min(span, dim(LD_matrix)[1],na.rm = T)
  sub_DT <- subset(subset_DT, SNP %in% rownames(LD_matrix))

  if(method[1]=="gaston"){
    gaston::LD.plot(LD = LD_matrix[start_pos:end_pos,
                                   start_pos:end_pos],
                     snp.positions = sub_DT$POS[start_pos:end_pos] )
  }
  if(method[1]=="heatmap"){
    heatmap(as.matrix(LD_sparse)[start_pos:end_pos,
                                 start_pos:end_pos])
  }
  if(method[1]=="image"){
    image(as.matrix(LD_sparse)[start_pos:end_pos,
                               start_pos:end_pos])
  }

}





#' Download vcf subset from 1000 Genomes
#'
#' @family LD
#' @keywords internal
#' @param query_by_regions You can make queries with \code{tabix} in two different ways:
#' \describe{
#' \item{\code{query_by_regions=F} \emph{(default)}}{Return a vcf with all positions between the min/max in \code{subset_DT} Takes up more storage but is MUCH faster}
#' \item{\code{query_by_regions=T}}{Return a vcf with only the exact positions present in \code{subset_DT}. Takes up less storage but is MUCH slower}
#' }
#' @inheritParams finemap_pipeline
#' @examples
#' \dontrun{
#' data("BST1");
#' subset_DT <- BST1
#' vcf_subset.popDat <- LD.1KG_download_vcf(subset_DT=BST1, LD_reference="1KGphase1", locus_dir=file.path("~/Desktop",locus_dir))
#' }
LD.1KG_download_vcf <- function(subset_DT,
                                LD_reference="1KGphase1",
                                remote_LD=T,
                                vcf_folder=NULL,
                                locus_dir,
                                locus=NULL,
                                whole_vcf=F,
                                download_method="wget",
                                force_new_vcf=F,
                                query_by_regions=F,
                                nThread=4,
                                conda_env="echoR",
                                verbose=T){
  # throw error if anything but phase 1 or phase 3 are specified
  if( ! LD_reference %in% c("1KGphase1", "1KGphase3" )){
    stop("LD_reference must be one of \"1KGphase1\" or \"1KGphase3\" ")
  }

  # Old FTP (deprecated?)
  ## http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
  # New FTP
  ## ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/
    # Download portion of vcf from 1KG website
  # vcf_folder <- LD.get_locus_vcf_folder(locus_dir=locus_dir)
  # Don't use the 'chr' prefix for 1KG queries:
  ## https://www.internationalgenome.org/faq/how-do-i-get-sub-section-vcf-file/
  subset_DT$CHR <- gsub("chr","",subset_DT$CHR)
  chrom <- unique(subset_DT$CHR)

  # PHASE 3 DATA
  if(LD_reference=="1KGphase3"){
    FTP <- "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/"
    # FTP <- "/sc/arion/projects/ad-omics/data/references/1KGPp3v5/"
    popDat <- echolocatoR::popDat_1KGphase3
    printer("LD Reference Panel = 1KGphase3", v=verbose)
    if(remote_LD){## With internet
      printer("+ LD:: Querying 1KG remote server.",v=verbose)
      vcf_URL <- paste0(FTP,"/ALL.chr",chrom,
                       ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
    }else{## WithOUT internet
      # vcf_folder <-  "/sc/arion/projects/ad-omics/data/references/1KGPp3v5/"
      printer("+ LD:: Querying 1KG local vcf files.",v=verbose)
      vcf_URL <- paste(vcf_folder, "/ALL.chr",chrom,
                       ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",sep="")
    }

    # PHASE 1 DATA
  } else if (LD_reference=="1KGphase1") {
    FTP <- "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/"
    popDat <- echolocatoR::popDat_1KGphase1
    printer("LD Reference Panel = 1KGphase1", v=verbose)
    if(remote_LD){## With internet
      printer("+ LD:: Querying 1KG remote server.",v=verbose)
      vcf_URL <- paste(FTP,"/ALL.chr",chrom,
                       ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")
    }else{## WithOUT internet
      printer("+ LD:: Querying 1KG local vcf files.",v=verbose)
      vcf_URL <- paste(vcf_folder,"/ALL.chr",chrom,
                       ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")
    }
  }
  phase <- gsub("1KG","",LD_reference)
  # phase 1 has no header whereas phase 3 does
  if( LD_reference == "1KGphase1" ){ use_header <- FALSE }
  if( LD_reference == "1KGphase3" ){ use_header <- TRUE }
  # Download and subset vcf if the subset doesn't exist already
  vcf_subset <- LD.query_vcf(subset_DT=subset_DT,
                             vcf_URL=vcf_URL,
                             locus_dir=locus_dir,
                             LD_reference=LD_reference,
                             whole_vcf=whole_vcf,
                             force_new_vcf=force_new_vcf,
                             download_method=download_method,
                             query_by_regions=query_by_regions,
                             nThread=nThread,
                             conda_env=conda_env,
                             verbose=verbose)
  return(list(vcf_subset = vcf_subset,
              popDat = popDat))
}




#' Get VCF storage folder
#' @family LD
#' @keywords internal
#' @examples
#' data("locus_dir")
#' vcf_folder <- LD.get_locus_vcf_folder(locus_dir=locus_dir)
LD.get_locus_vcf_folder <- function(locus_dir=NULL){
    vcf_folder <- file.path(locus_dir,"LD")
    out <- dir.create(vcf_folder, showWarnings = F, recursive = T)
  return(vcf_folder)
}




#' Construct the path to vcf subset
#'
#' @family LD
#' @keywords internal
#' @examples
#' data("locus_dir"); data("BST1");
#' vcf_subset <- LD.construct_subset_vcf_name(subset_DT=BST1, locus_dir=locus_dir, LD_reference="1KGlocal")
LD.construct_subset_vcf_name <- function(subset_DT,
                                         LD_reference=NULL,
                                         locus_dir,
                                         whole_vcf=F){
  vcf_folder <- LD.get_locus_vcf_folder(locus_dir=locus_dir)
  # Don't use the chr prefix: https://www.internationalgenome.org/faq/how-do-i-get-sub-section-vcf-file/
  subset_DT$CHR <- gsub("chr","",subset_DT$CHR)
  chrom <- unique(subset_DT$CHR)
  vcf_subset <- file.path(vcf_folder,
                          if(whole_vcf){
                            paste(basename(LD_reference), paste0("chr",chrom),sep=".")
                          } else {paste(basename(locus_dir),basename(LD_reference),sep=".")} )
  dir.create(path = dirname(vcf_subset), recursive = T, showWarnings = F)
  if(!(endsWith(vcf_subset,".vcf.gz")|
       endsWith(vcf_subset,".vcf"))){
    vcf_subset <- paste0(vcf_subset,".vcf")
  }
  return(vcf_subset)
}




#' Index vcf file if it hasn't been already
#'
#' @family LD
#' @keywords internal
#' @examples
#' \dontrun{
#' LD_reference <- "~/Desktop/results/Reference/custom_panel_chr4.vcf.gz"
#' vcf_file <- LD.index_vcf(vcf_file=LD_reference)
#' }
LD.index_vcf <- function(vcf_file,
                         force_new_index=F,
                         conda_env="echoR",
                         verbose=T){
  if(!endsWith(vcf_file,".gz")){
    printer("+ LD:: Compressing vcf with bgzip",v=verbose)
    bgzip <- CONDA.find_package(package="bgzip",
                                conda_env=conda_env,
                                verbose = verbose)
    cmd1 <- paste(bgzip,
                  vcf_file)
    printer("++ LD::",cmd1,v=verbose)
    system(cmd1)
    vcf_file <- paste0(vcf_file,".gz")
  }else { printer("+ LD::",vcf_file,"already compressed",v=verbose)}

  if(!file.exists(paste0(vcf_file,".tbi")) | force_new_index){
    printer("+ LD:: Indexing",vcf_file,v=verbose)
    tabix <- CONDA.find_package(package="tabix",
                                conda_env=conda_env,
                                verbose = verbose)
    cmd <- paste(tabix,
                 "-fp vcf",
                 vcf_file)
    printer("++ LD::",cmd,v=verbose)
    system(cmd)
  } else { printer("+ LD::",vcf_file,"already indexed.",v=verbose) }
  return(vcf_file)
}




#' Query vcf file
#'
#' @family LD
#' @keywords internal
#' @examples
#' \dontrun{
#' data("locus_dir"); data("BST1");
#' # Custom
#' LD_reference <- "~/Desktop/results/Reference/custom_panel_chr4.vcf"
#' vcf_file <- LD.index_vcf(vcf_file=LD_reference)
#' vcf_subset <- LD.query_vcf(subset_DT=BST1, locus_dir=locus_dir, vcf_URL=vcf_file, LD_reference=LD_reference, force_new_vcf=T)
#' }
LD.query_vcf <- function(subset_DT,
                         vcf_URL,
                         locus_dir,
                         LD_reference,
                         whole_vcf=F,
                         force_new_vcf=F,
                         remove_original_vcf=F,
                         download_method="wget",
                         query_by_regions=F,
                         nThread=4,
                         conda_env="echoR",
                         verbose=T){
  # vcf_subset <- "/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/BRIP1/LD/BRIP1.1KGphase3.vcf.gz"
  vcf_subset <- LD.construct_subset_vcf_name(subset_DT=subset_DT,
                                             locus_dir=locus_dir,
                                             LD_reference=LD_reference,
                                             whole_vcf=whole_vcf)
  # CHECK FOR EMPTY VCF FILES!
  ## These can be created if you stop the query early, or if the url fails.
  if(file.exists(vcf_subset)){
    if(file.size(vcf_subset) < 100){ # Less than 100 bytes
      printer("+ LD:: Removing empty vcf file and its index", v=verbose)
      file.remove(paste0(vcf_subset,"*")) # Remove both
    }
  }
  tabix <- CONDA.find_package(package = "tabix",
                              conda_env = conda_env,
                              verbose = verbose)
  if((!file.exists(vcf_subset)) | force_new_vcf){
    printer("LD:: Querying VCF subset", v=verbose)
    if(whole_vcf){
      region <- ""
      locus <- ""
      out.file <- downloader(input_url = vcf_URL,
                             output_path = dirname(vcf_subset),
                             download_method = download_method,
                             nThread = nThread,
                             conda_env=conda_env)
    } else {
      # Download tabix subset
      if(query_by_regions){
        ### Using region file (-R flag)
        regions.bed <- file.path(locus_dir,"LD","regions.tsv")
        data.table::fwrite(list(paste0(subset_DT$CHR), sort(subset_DT$POS)),
                           file=regions.bed, sep="\t")
        regions <- paste("-R",regions.bed)
        tabix_cmd <- paste(tabix,
                           "-fh",
                           "-p vcf",
                           vcf_URL,
                           gsub("\\./","",regions),
                           ">",
                           gsub("\\./","",vcf_subset) )
        printer(tabix_cmd)
        system(tabix_cmd)
      } else {
        ### Using coordinates range (MUCH faster!)
        coord_range <- paste0(unique(subset_DT$CHR)[1],":",
                              min(subset_DT$POS),"-",max(subset_DT$POS))
        tabix_cmd <- paste(tabix,
                           "-fh",
                           "-p vcf",
                           vcf_URL,
                           coord_range,
                           ">",
                           gsub("\\./","",vcf_subset) )
        printer(tabix_cmd)
        system(tabix_cmd)
      }
    }

    if(remove_original_vcf){
      vcf_name <- paste0(basename(vcf_URL), ".tbi")
      out <- suppressWarnings(file.remove(vcf_name))
    }
  } else {printer("+ Identified existing VCF subset file. Importing...", vcf_subset, v=verbose)}
  return(vcf_subset)
}


LD.query_vcf_variantannotation <- function(vcf_url,
                                           subset_DT){
  # vcf_url <- "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521//ALL.chr18.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
  # subset_DT <- echolocatoR::BST1
  # https://github.com/Bioconductor/VariantAnnotation/issues/29
  gr <-  DT_to_GRanges(subset_DT)

  #### Rsamtools ####
  # tab <- Rsamtools::scanTabix(file = vcf_url, gr)

  #### seqminer ####
  # tab <- seqminer::tabix.read(tabixFile = vcf_url, tabixRange = coord_range )

  #### VariantAnnotation::scanVcf ####
  vcf <- VariantAnnotation::VcfFile(file = vcf_url,
                                    index = paste0(vcf_url,".tbi"),
                                    yieldSize = 5000)
  line <- VariantAnnotation::scanVcf(file=vcf)
  #### VariantAnnotation::readVcf ####
  vcf <- VariantAnnotation::VcfFile(file = vcf_url,
                                    index = paste0(vcf_url,".tbi"))
  line2 <- VariantAnnotation::readVcf(file=vcf,
                                      genome = "hg19", param = gr)

  open(vcf)
  repeat {
    result <- VariantAnnotation::readVcf(vcf)
    if (length(result) == 0)
      break
    ## work on this chunk of the file
    message(nrow(result))
  }
  close(vcf)
  return(vcf)
}




#' Filter a vcf by min/max coordinates
#'
#' Uses \emph{bcftools} to filter a vcf by min/max genomic coordinates
#' (in basepairs).
#' @param vcf_subset Path to the locus subset vcf.
#' @param popDat The metadata file listing the superpopulation
#' to which each sample belongs.
#' @inheritParams finemap_pipeline
#' @family LD
#' @keywords internal
LD.filter_vcf <- function(vcf_subset,
                          popDat,
                          superpopulation,
                          remove_tmp=T,
                          verbose=T){
  vcf.gz <- paste0(vcf_subset,".gz")
  vcf.gz.subset <- gsub("_subset","_samples_subset",vcf.gz)
  # Compress vcf
  if(!file.exists(vcf.gz)){
    printer("LD:BCFTOOLS:: Compressing vcf file...", v=verbose)
    system(paste("bgzip -f",vcf_subset))
  }
  # Re-index vcf
  printer("LD:TABIX:: Re-indexing vcf.gz...", v=verbose)
  system(paste("tabix -f -p vcf",vcf.gz))
  # Subset samples
  selectedInds <- subset(popDat, superpop == superpopulation)$sample %>% unique()
  printer("LD:BCFTOOLS:: Subsetting vcf to only include",superpopulation,"individuals (",length(selectedInds), "/",length(popDat$sample%>%unique()),").", v=verbose)
  cmd <- paste("bcftools view -s",paste(selectedInds, collapse=","), vcf.gz, "| bgzip > tmp && mv tmp",vcf.gz.subset)
  system(cmd)
  # Remove old vcf
  if(remove_tmp){out <- suppressWarnings(file.remove(vcf_subset))}
  return(vcf.gz.subset)
}




#' Subset a vcf by superpopulation
#'
#' @inheritParams LD.filter_vcf
#' @family LD
#' @keywords internal
LD.filter_vcf_gaston <- function(vcf_subset,
                                 subset_DT,
                                 locus_dir,
                                 superpopulation,
                                 popDat,
                                 verbose=T){
  # Import w/ gaston and further subset
  printer("+ Importing VCF as bed file...", v=verbose)
  bed.file <- gaston::read.vcf(vcf_subset, verbose = F)
  ## Subset rsIDs
  bed <- gaston::select.snps(bed.file, id %in% subset_DT$SNP & id !=".")
  # Create plink sub-dir
  dir.create(file.path(locus_dir, "LD"), recursive = T, showWarnings = F)
  gaston::write.bed.matrix(bed, file.path(locus_dir, "LD/plink"), rds = NULL)
  # Subset Individuals
  selectedInds <- subset(popDat, superpop == superpopulation)
  bed <- gaston::select.inds(bed, id %in% selectedInds$sample)
  # Cleanup extra files
  remove(bed.file)
  # file.remove("vcf_subset")
  return(bed)
}




#' Convert vcf file to BED file
#'
#' Uses plink to convert vcf to BED.
#' @param vcf.gz.subset Path to the gzipped locus subset vcf.
#' @param locus_dir Locus-specific results directory.
#' @family LD
#' @keywords internal
LD.vcf_to_bed <- function(vcf.gz.subset,
                          locus_dir,
                          plink_prefix="plink",
                          verbose=T){
  plink <- LD.plink_file()
  printer("LD:PLINK:: Converting vcf.gz to .bed/.bim/.fam", v=verbose)
  LD_dir <- file.path(locus_dir,"LD")
  dir.create(LD_dir, recursive = T, showWarnings = F)
  cmd <- paste(plink,
               "--vcf",vcf.gz.subset,
               "--out", file.path(LD_dir,plink_prefix))
  system(cmd)

  return(
    list(bed=file.path(LD_dir,paste0(plink_prefix,".bed")),
         bim=file.path(LD_dir,paste0(plink_prefix,".bim")),
         fam=file.path(LD_dir,paste0(plink_prefix,".fam")))
  )
}




#' Calculate LD
#'
#' Calculate a pairwise LD matrix from a vcf file using \emph{plink}.
#' @param locus_dir Locus-specific results directory.
#' @param ld_window Set --r/--r2 max variant ct pairwise distance (usu. 10).
#' @param ld_format Whether to produce an LD matrix with
#' r (\code{ld_format="r"}) or D' (\code{ld_format="D"}) as the pairwise SNP correlation metric.
#' @family LD
#' @keywords internal
LD.calculate_LD <- function(locus_dir,
                            ld_window=1000, # 10000000
                            ld_format="r",
                            plink_prefix="plink",
                            verbose=T){
  plink <- LD.plink_file()
  printer("LD:PLINK:: Calculating LD ( r & D'-signed; LD-window =",ld_window,")", v=verbose)
  plink_path_prefix <- file.path(locus_dir,"LD",plink_prefix)
  dir.create(file.path(locus_dir,"LD"), recursive = T, showWarnings = F)
  out_prefix <- paste0(plink_path_prefix,".r_dprimeSigned")
  if(ld_format=="r"){
    cmd <- paste(plink,
                 "--bfile",plink_path_prefix,
                 "--r square bin",
                 "--out",out_prefix)
    ld.path <- paste0(out_prefix,".ld.bin")
  } else {
    cmd <- paste(plink,
                 "--bfile",plink_path_prefix,
                 "--r dprime-signed",
                 "--ld-window",ld_window,
                 "--ld-window-kb",ld_window,
                 "--out",out_prefix)
    ld.path <- paste0(out_prefix,".ld")
  }
  system(cmd)
  return(ld.path)
}





#' Create LD matrix from plink output
#'
#' Depending on which parameters you give \emph{plink} when calculating LD, you get different file outputs.
#' When it produces bin and bim files, use this function to create a proper LD matrix.
#' For example, this happens when you try to calculate D' with the \code{--r dprime-signed} flag (instead of just r).
#' @param LD_dir Directory that contains the bin/bim files.
#' @family LD
#' @keywords internal
#' @examples
#' \dontrun{
#' locus_dir <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/QTL/Microglia_all_regions/BIN1"
#' ld.matrix <- LD.read_bin(LD_dir=file.path(locus_dir, "LD"))
#' }

LD.read_bin <- function(LD_dir){
  bim <- data.table::fread(file.path(LD_dir, "plink.bim"), col.names = c("CHR","SNP","V3","POS","A1","A2"))
  bin.vector <- readBin(file.path(LD_dir, "plink.ld.bin"), what = "numeric", n=length(bim$SNP)^2)
  ld.matrix <- matrix(bin.vector, nrow = length(bim$SNP), dimnames = list(bim$SNP, bim$SNP))
  return(ld.matrix)
}





#' Create LD matrix from plink output.
#'
#' Depending on which parameters you give \emph{plink} when calculating LD, you get different file outputs.
#' When it produces an LD table, use this function to create a proper LD matrix.
#' @family LD
#' @keywords internal
LD.read_ld_table <- function(ld.path,
                             snp.subset=F,
                             fillNA=0,
                             verbose=T){
  # subset_DT <- data.table::fread(file.path(locus_dir,"Multi-finemap/Multi-finemap_results.txt")); snp.subset <- subset_DT$SNP
  ld.table <- data.table::fread(ld.path, nThread = 4)
  if(any(snp.subset!=F)){
    printer("LD:PLINK:: Subsetting LD data...", v=verbose)
    ld.table <- subset(ld.table, SNP_A %in% snp.subset | SNP_B %in% snp.subset)
  }
  printer("LD:PLINK:: Casting data.matrix...", v=verbose)
  ld.cast <- data.table::dcast.data.table(ld.table,
                                          formula = SNP_B ~ SNP_A,
                                          value.var="R",
                                          fill=0,
                                          drop=T,
                                          fun.agg = function(x){mean(x,na.rm = T)})
  ld.cast <- subset(ld.cast, SNP_B !=".", select = -`.`)
  ld.mat <- data.frame(ld.cast, row.names = ld.cast$SNP_B) %>% data.table() %>% as.matrix()
  # ld.mat[1:10,1:10]
  ld.mat <- LD.fill_NA(LD_matrix = ld.mat,
                       fillNA = fillNA,
                       verbose = verbose)
  return(ld.mat)
}




#' Find correct plink file
#'
#' @family LD
#' @keywords internal
#' @examples
#' plink <- LD.plink_file()
LD.plink_file <- function(plink="plink",
                          conda_env=NULL){
  if(!is.null(conda_env)){
    plink <- CONDA.find_package("plink",conda_env = conda_env)
  }
  if(plink=="plink"){
    base_url=system.file("tools/plink",package = "echolocatoR")
    os <- get_os()
    if (os=="osx") {
      plink <- file.path(base_url, "plink1.9_mac");
    } else if (os=="linux") {
      plink <- file.path(base_url, "plink1.9_linux");
    } else {
      plink <- file.path(base_url, "plink1.9_windows.exe");
    }
  }
  return(plink)
}




#' Compute LD from 1000 Genomes
#'
#' Downloads a subset vcf of the 1KG database that matches your locus coordinates.
#' Then uses \emph{plink} to calculate LD on the fly.
#'
#' This approach is taken, because other API query tools have limitations with the window size being queried.
#' This approach does not have this limitations, allowing you to fine-map loci more completely.
#'
#' @param fillNA When pairwise LD (r) between two SNPs is \code{NA}, replace with 0.
#' @inheritParams finemap_pipeline
#' @family LD
#' @keywords internal
#' @examples
#' \dontrun{
#' data("BST1"); data("locus_dir");
#' BST1 <- limit_SNPs(max_snps = 500, subset_DT = BST1)
#' LD_matrix <- LD.1KG(locus_dir=file.path("~/Desktop",locus_dir), subset_DT=BST1, LD_reference="1KGphase1")
#'
#' ## Kunkle et al 2019
#' locus_dir <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Kunkle_2019/ACE"
#'
#' }
LD.1KG <- function(locus_dir,
                   subset_DT,
                   LD_reference="1KGphase1",
                   superpopulation="EUR",
                   vcf_folder=NULL,
                   remote_LD=T,
                   # min_r2=F,
                   LD_block=F,
                   LD_block_size=.7,
                   # min_Dprime=F,
                   remove_correlates=F,
                   remove_tmps=T,
                   fillNA=0,
                   download_method="wget",
                   nThread=4,
                   conda_env="echoR",
                   verbose=T){
  # data("BST1"); data("locus_dir"); subset_DT=BST1; LD_reference="1KGphase1"; vcf_folder=NULL; superpopulation="EUR";
  # min_r2=F; LD_block=F; LD_block_size=.7; min_Dprime=F;  remove_correlates=F; remote_LD=T; verbose=T; nThread=4; conda_env="echoR";
  printer("LD:: Using 1000Genomes as LD reference panel.", v=verbose)
  locus <- basename(locus_dir)
  vcf_info <- LD.1KG_download_vcf(subset_DT=subset_DT,
                                  locus_dir=locus_dir,
                                  LD_reference=LD_reference,
                                  vcf_folder=vcf_folder,
                                  locus=locus,
                                  remote_LD=remote_LD,
                                  download_method=download_method,
                                  nThread=nThread,
                                  conda_env=conda_env,
                                  verbose=verbose)
  vcf_subset <- vcf_info$vcf_subset
  popDat <- vcf_info$popDat
  vcf.gz.path <- LD.filter_vcf(vcf_subset = vcf_subset,
                               popDat = popDat,
                               superpopulation = superpopulation,
                               remove_tmp = T,
                               verbose = verbose)
  bed_bim_fam <- LD.vcf_to_bed(vcf.gz.subset = vcf.gz.path,
                               locus_dir = locus_dir,
                               verbose = verbose)
  # Calculate LD
  LD_matrix <- LD.snpstats_get_LD(LD_folder=file.path(locus_dir,"LD"),
                                  plink_prefix="plink",
                                  select.snps=unique(subset_DT$SNP),
                                  stats=c("R"),
                                  symmetric=T,
                                  depth="max",
                                  verbose=verbose)
  # Get MAF (if needed)
  subset_DT <- LD.snpstats_get_MAF(subset_DT=subset_DT,
                                   LD_folder=file.path(locus_dir,"LD"),
                                   plink_prefix="plink",
                                   force_new_MAF=T,
                                   nThread=nThread,
                                   verbose=verbose)
  # Get lead SNP rsid
  leadSNP = subset(subset_DT, leadSNP==T)$SNP
  # Filter out SNPs not in the same LD block as the lead SNP
  if(LD_block){
    block_snps <- LD.leadSNP_block(leadSNP = leadSNP,
                                   LD_folder = file.path(locus_dir,"LD","plink_tmp"),
                                   LD_block_size = LD_block_size)
    LD_matrix <- LD_matrix[row.names(LD_matrix) %in% block_snps, colnames(LD_matrix) %in% block_snps]
    LD_matrix <- LD_matrix[block_snps, block_snps]
  }
  # IMPORTANT! Remove large data.ld file after you're done with it
  if(remove_tmps){
    suppressWarnings(file.remove(vcf_subset))
  }
  # Save LD matrix
  LD_list <- LD.save_LD_matrix(LD_matrix=LD_matrix,
                               subset_DT=subset_DT,
                               locus_dir=locus_dir,
                               subset_common = T,
                               sparse = T,
                               fillNA=fillNA,
                               LD_reference=LD_reference,
                               verbose=verbose)
  return(LD_list)
}




#' Calculate LD (D')
#'
#' This appriach computes an LD matrix of D' (instead of r or r2) from a vcf.
#' See \code{\link{LD.run_plink_LD}} for a faster (but less flexible) alternative to computing LD.
#' @family LD
#' @keywords internal
LD.dprime_table <- function(SNP_list,
                            LD_folder,
                            conda_env){
  plink <- LD.plink_file(conda_env)
  printer("+ Creating DPrime table")
  system( paste(plink, "--bfile",file.path(LD_folder,"plink"),
                "--ld-snps", paste(SNP_list, collapse=" "),
                "--r dprime-signed",
                "--ld-window 10000000", # max out window size
                "--ld-window-kb 10000000",
                "--out",file.path(LD_folder,"plink")) )
  #--ld-window-r2 0

  # # Awk method: theoretically faster?
  # if(min_Dprime==F){Dprime = -1}else{Dprime=min_Dprime}
  # if(min_r2==F){r = -1}else{r = round(sqrt(min_r2),2) }
  # columns <- data.table::fread(file.path(LD_folder, "plink.ld"), nrows = 0) %>% colnames()
  # col_dict <- setNames(1:length(columns), columns)
  # awk_cmd <- paste("awk -F \"\t\" 'NR==1{print $0}{ if(($",col_dict["DP"]," >= ",Dprime,")",
  #                  " && ($",col_dict["R"]," >= ",r,")) { print } }' ",file.path(LD_folder, "plink.ld"),
  #                  " > ",file.path(LD_folder, "plink.ld_filtered.txt"),  sep="")
  # system(awk_cmd)
  plink.ld <- data.table::fread(file.path(LD_folder, "plink.ld"), select = c("SNP_A", "SNP_B","DP","R"), )
  plink.ld <- plink.ld[complete.cases(plink.ld) ]
  return(plink.ld)
}




#' Get LD using \pkg{snpStats} package
#'
#' @param LD_folder Locus-specific LD output folder.
#' @inheritParams snpStats::ld
#' @family LD
#' @keywords internal
#' @source
#' \href{https://www.bioconductor.org/packages/release/bioc/html/snpStats.html}{snpStats Bioconductor page}
#' \href{https://www.bioconductor.org/packages/release/bioc/vignettes/snpStats/inst/doc/ld-vignette.pdf}{LD tutorial}
#' @examples
#' subset_DT <- data.table::fread("/pd-omics/brian/Fine_Mapping/Data/GWAS/Kunkle_2019/ABCA7/Multi-finemap/ABCA7.Kunkle_2019.1KGphase3_LD.Multi-finemap.tsv.gz")
#' LD_folder <- "/pd-omics/brian/Fine_Mapping/Data/GWAS/Kunkle_2019/ABCA7/LD"
#' LD_matrix <- LD.snpstats_get_LD(LD_folder=LD_folder, select.snps=subset_DT$SNP)
LD.snpstats_get_LD <- function(LD_folder,
                               plink_prefix="plink",
                               select.snps=NULL,
                               stats=c("R"),
                               symmetric=T,
                               depth="max",
                               nThread=4,
                               verbose=T){
  printer("LD:snpStats:: Computing LD",paste0("(stats = ",paste(stats,collapse=', '),")"),v=verbose)
  # select.snps= arg needed bc otherwise read.plink() sometimes complains of
  ## duplicate RSID rownames. Also need to check whether these SNPs exist in the plink files.
  ## (snpStats doesn't have very good error handling for these cases).
  select.snps <- LD.snpstats_ensure_nonduplicates(select.snps=select.snps,
                                                  LD_folder=LD_folder,
                                                  plink_prefix=plink_prefix,
                                                  nThread=nThread,
                                                  verbose=verbose)
  # Only need to give bed path (infers bin/fam paths)
  ss <- snpStats::read.plink(bed = file.path(LD_folder,plink_prefix),
                             select.snps = select.snps)
  # Compute LD from snpMatrix
  ld_list <- snpStats::ld(x = ss$genotypes,
                          y = ss$genotypes,
                          depth = if(depth=="max") ncol(ss$genotypes) else depth,
                          stats = stats,
                          symmetric = symmetric)
  if(length(stats)==1) return(ld_list) else return(ld_list$R)
}




LD.snpstats_ensure_nonduplicates <- function(select.snps=NULL,
                                             LD_folder,
                                             plink_prefix="plink",
                                             nThread=4,
                                             verbose=T){
  if(!is.null(select.snps)){
    bim_path <- file.path(LD_folder,paste0(plink_prefix,".bim"))
    bim <- data.table::fread(bim_path,
                             col.names = c("CHR","SNP","V3","POS","A1","A2"),
                             stringsAsFactors = F,
                             nThread=nThread)
    printer("+ LD:snpStats::",nrow(bim),"rows in bim file.",v=verbose)
    bim <- bim[!duplicated(bim$SNP),]
    select.snps <- select.snps[select.snps %in% unique(bim$SNP)]
    printer("+ LD:snpStats::",length(select.snps),"SNPs in select.snps.",v=verbose)
    select.snps <- if(length(select.snps)==0) NULL else unique(select.snps);
  }
  return(select.snps)
}




#' Get MAF using \pkg{snpStats} package
#'
#' @param LD_folder Locus-specific LD output folder.
#' @inheritParams snpStats::ld
#' @family LD
#' @keywords internal
#' @source
#' \href{https://www.bioconductor.org/packages/release/bioc/html/snpStats.html}{snpStats Bioconductor page}
LD.snpstats_get_MAF <- function(subset_DT,
                                 LD_folder,
                                 plink_prefix="plink",
                                 force_new_MAF=F,
                                 nThread=4,
                                 verbose=T){
  if(!"MAF" %in% colnames(subset_DT) | force_new_MAF){
    printer("LD::snpStats:: Filling `MAF` column with MAF from LD panel.",v=verbose)
    select.snps <- LD.snpstats_ensure_nonduplicates(select.snps=subset_DT$SNP,
                                                    LD_folder=LD_folder,
                                                    plink_prefix=plink_prefix,
                                                    nThread=nThread,
                                                    verbose=F)
    ss <- snpStats::read.plink(bed = file.path(LD_folder,plink_prefix),
                               select.snps = select.snps)

    MAF_df <- data.frame(SNP=row.names(snpStats::col.summary(ss$genotypes)),
                         MAF=snpStats::col.summary(ss$genotypes)$MAF)
    if("MAF" %in% colnames(subset_DT)) subset_DT <- subset(subset_DT,select=-MAF)
    subset_merge <- data.table::merge.data.table(data.table::data.table(subset_DT),
                                                 data.table::data.table(MAF_df),
                                                 by="SNP")
    return(subset_merge)
  } else {
      printer("LD::snpStats:: `MAF` column already present.",v=verbose);
      return(subset_DT)
  }
}


#' Calculate LD (r or r2)
#'
#' This appriach computes and LD matrix of r or r2 (instead of D') from a vcf.
#' See \code{\link{LD.dprime_table}} for a slower (but more flexible) alternative to computing LD.
#' @param bim A bim file produced by \emph{plink}
#' @param LD_folder Locus-specific LD output folder.
#' @param r_format Whether to fill the matrix with \code{r} or \code{r2}.
#' @family LD
#' @keywords internal
#' @examples
#' \dontrun{
#' data("LRRK2")
#' LD_folder <- "/Users/schilder/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/LRRK2/plink/saved"
#' bim_path <- file.path(LD_folder, "plink.bim");
#' bim <- data.table::fread(bim_path, col.names = c("CHR","SNP","V3","POS","A1","A2"), stringsAsFactors = F)
#' bim <- subset(bim, SNP %in% LRRK2$SNP)
#' ld.bin <- file.path(LD_folder, paste0("plink",".ld.bin"))
#' SNPs <- data.table::fread(file.path(LD_folder,"SNPs.txt"), col.names = 'RSID')
#' bin.vector <- readBin(ld.bin, what = "numeric", n=length(SNPs$RSID)^2)
#' }
LD.run_plink_LD <- function(bim,
                            LD_folder,
                            plink_prefix="plink",
                            r_format="r",
                            extract_file=NULL){
  plink <- LD.plink_file()
  # METHOD 2 (faster, but less control over parameters. Most importantly, can't get Dprime)
  system( paste(plink,
                "--bfile",file.path(LD_folder,plink_prefix),
                if(is.null(extract_file)) NULL else"--extract",extract_file,
                paste0("--",r_format," square bin"),
                "--out", file.path(LD_folder,plink_prefix)) )
  ld.bin <- file.path(LD_folder, paste0(plink_prefix,".ld.bin"))
  bin.vector <- readBin(ld.bin, what = "numeric", n=length(bim$SNP)^2)
  ld.matrix <- matrix(bin.vector, nrow = length(bim$SNP), dimnames = list(bim$SNP, bim$SNP))
  return(ld.matrix)
}




#' Calculate LD
#'
#' Use \emph{plink} to calculate LD from a vcf.
#' @family LD
#' @keywords internal
#' @examples
#' locus_dir <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Kunkle_2019/ACE"
#' LD_folder <- file.path(locus_dir,"LD")
#' ld.matrix <- LD.plink_LD(subset_DT=BST1, LD_folder=LD_folder)
LD.plink_LD <-function(leadSNP=NULL,
                       subset_DT,
                       bim_path=NULL,
                       remove_excess_snps=T,
                       # IMPORTANT! keep this F
                       merge_by_RSID=F,
                       LD_folder,
                       min_r2=F,
                       min_Dprime=F,
                       remove_correlates=F,
                       fillNA=0,
                       plink_prefix="plink",
                       verbose=T,
                       conda_env=NULL){
  # Dprime ranges from -1 to 1
  start <- Sys.time()
  if(is.null(leadSNP))leadSNP <- subset(subset_DT, leadSNP)$SNP[1]
  # Calculate LD
  printer("++ Reading in BIM file...", v=verbose)
  if(is.null(bim_path)) bim_path <- file.path(LD_folder, "plink.bim");
  bim <- data.table::fread(bim_path,
                           col.names = c("CHR","SNP","V3","POS","A1","A2"),
                           stringsAsFactors = F)
  if(remove_excess_snps){
    orig_n <- nrow(bim)
    if(merge_by_RSID){
      bim.merged <- data.table::merge.data.table(bim,
                                                 subset_DT,
                                                 by=c("SNP"))
    } else {
      # Standardize format adn merge
      bim.merged <- data.table::merge.data.table(dplyr::mutate(bim,
                                                               CHR=as.integer(gsub("chr","",CHR)),
                                                               POS=as.integer(POS)),
                                                 dplyr::mutate(subset_DT,
                                                               CHR=as.integer(gsub("chr","",CHR)),
                                                               POS=as.integer(POS)),
                                                 by=c("CHR","POS"))
    }
    bim <- subset(bim, SNP %in% bim.merged$SNP.x)
    printer("LD:PLINK:: Removing RSIDs that don't appear in locus subset:",orig_n,"==>",nrow(bim),"SNPs",v=verbose)
  }
  extract_file <- file.path(LD_folder,"SNPs.txt")
  data.table::fwrite(subset(bim, select="SNP"),
                     extract_file, col.names = F)

  printer("++ Calculating LD", v=verbose)
  ld.matrix <- LD.run_plink_LD(bim = bim,
                               LD_folder = LD_folder,
                               plink_prefix = plink_prefix,
                               extract_file = file.path(LD_folder,"SNPs.txt"))

  if((min_Dprime != F) | (min_r2 != F) | (remove_correlates != F)){
    plink.ld <- LD.dprime_table(SNP_list = row.names(ld.matrix),
                                LD_folder,
                                conda_env=conda_env)

    # DPrime filter
    if(min_Dprime != F){
      printer("+++ Filtering LD Matrix (min_Dprime): Removing SNPs with D' <=",min_Dprime,"for",leadSNP,"(lead SNP).", v=verbose)
      plink.ld <- subset(plink.ld, (SNP_A==leadSNP & DP>=min_Dprime) | (SNP_B==leadSNP & DP>=min_Dprime))
    } else{printer("+ min_Dprime == FALSE", v=verbose)}

    # R2 filter
    if(min_r2 != F ){
      printer("+++ Filtering LD Matrix (min_r2): Removing SNPs with r <=",min_r2,"for",leadSNP,"(lead SNP).", v=verbose)
      r = sqrt(min_r2) # PROBLEM: this doesn't give you r, which you need for SUSIE
      plink.ld <- subset(plink.ld, (SNP_A==leadSNP & R>=r) | (SNP_B==leadSNP & R>=r))
    } else{printer("+ min_r2 == FALSE", v=verbose)}

    # Correlates filter
    if(remove_correlates != F){
      r2_threshold <- remove_correlates# 0.2
      r <- sqrt(r2_threshold)
      printer("+++ Filtering LD Matrix (remove_correlates): Removing SNPs with R2 >=",r2_threshold,"for",paste(remove_correlates,collapse=", "),".", v=verbose)
      plink.ld <- subset(plink.ld, !(SNP_A %in% remove_correlates & R>=r) | (SNP_B %in% remove_correlates & R>=r))
    } else{printer("+ remove_correlates == FALSE", v=verbose)}

    # Apply filters
    A_list <- unique(plink.ld$SNP_A)
    B_list <- unique(plink.ld$SNP_B)
    snp_list <-   unique(c(A_list, B_list))
    ld.matrix <- ld.matrix[row.names(ld.matrix) %in% snp_list, colnames(ld.matrix) %in% snp_list]
    ## Manually remove rare variant
    # ld.matrix <- ld.matrix[rownames(ld.matrix)!="rs34637584", colnames(ld.matrix)!="rs34637584"]
  }
  # !IMPORTANT!: Fill NAs (otherwise susieR will break)
  ld.matrix <- LD.fill_NA(LD_matrix = ld.matrix,
                          fillNA = fillNA,
                          verbose = verbose)
  end <- Sys.time()
  printer("+ LD matrix calculated in",round(as.numeric(end-start),2),"seconds.", v=verbose)
  return(ld.matrix)
}


#' Fill NAs in an LD matrix
#'
#' Trickier than it looks.
#' @examples
#' \dontrun{
#' data("LD_matrix");
#' LD_matrix <- LD.fill_NA(LD_matrix)
#' }
LD.fill_NA <- function(LD_matrix,
                       fillNA=0,
                       verbose=F){
  printer("+ LD:: Removing unnamed rows/cols", v=verbose)
  # First, filter any rows/cols without names
  LD_matrix <- data.frame(LD_matrix)
  LD_matrix <- LD_matrix[rownames(LD_matrix)!=".", colnames(LD_matrix)!="."]
  LD_matrix_orig <- LD_matrix

  if(!is.null(fillNA)){
    printer("+ LD:: Replacing NAs with",fillNA, v=verbose)
    if(sum(is.na(LD_matrix))>0){
      LD_matrix[is.na(LD_matrix)] <- 0
    }
  }
  # Check for duplicate SNPs
  LD_matrix <- LD_matrix[row.names(LD_matrix)[!duplicated(row.names(LD_matrix))],
                         colnames(LD_matrix)[!duplicated(colnames(LD_matrix))]]
  return(LD_matrix)
}




#' Calculate LD blocks.
#'
#' Uses \emph{plink} to group highly correlated SNPs into LD blocks.
#'
#' @family LD
#' @keywords internal
LD.LD_blocks <- function(LD_folder,
                         LD_block_size=.7){
  printer("++ Calculating LD blocks...")
  # PLINK 1.07 LD: http://zzz.bwh.harvard.edu/plink/ld.shtml
  # PLINK 1.9 LD: https://www.cog-genomics.org/plink/1.9/ld
  # system("plink", "-h")
  # Identify duplicate snps
  # system("plink", "--vcf subset.vcf --list-duplicate-vars")
  # Convert vcf to plink format
  # system("plink", "--vcf subset.vcf --exclude ./plink_tmp/plink.dupvar --make-bed --out PTK2B")

  # Estimate LD blocks
  # Defaults: --blocks-strong-lowci = 0.70, --blocks-strong-highci .90

  # Reducing "--blocks-inform-frac" is the only parameter that seems to make the block sizes larger
  plink <- LD.plink_file()
  system( paste(plink, "--bfile",file.path(LD_folder,"plink"),
                "--blocks no-pheno-req no-small-max-span --blocks-max-kb 100000",
                # "--blocks-strong-lowci .52 --blocks-strong-highci 1",
                "--blocks-inform-frac",LD_block_size," --blocks-min-maf 0 --out",file.path(LD_folder,"plink")) )
  # system( paste("plink", "--bfile plink --ld-snp-list snp_list.txt --r") )
  blocks <- data.table::fread("./plink_tmp/plink.blocks.det")
  return(blocks)
}


#' Identify the LD block in which the lead SNP resides
#' @family LD
#' @keywords internal
LD.leadSNP_block <- function(leadSNP, LD_folder, LD_block_size=.7){
  printer("Returning lead SNP's block...")
  blocks <- LD.LD_blocks(LD_folder, LD_block_size)
  splitLists <- strsplit(blocks$SNPS,split = "[|]")
  block_snps <- lapply(splitLists, function(l, leadSNP){if(leadSNP %in% l){return(l)} }, leadSNP=leadSNP) %>% unlist()
  printer("Number of SNPs in LD block =", length(block_snps))
  return(block_snps)
}


# LD_clumping <- function(vcf_subset, subset_SS){
#   # PLINK clumping: http://zzz.bwh.harvard.edu/plink/clump.shtml
#   # Convert vcf to .map (beagle)
#   ## https://www.cog-genomics.org/plink/1.9/data
#   system(paste("plink", "--vcf",vcf_subset,"--recode beagle --out ./plink_tmp/plink"))
#   # Clumping
#   system("plink", "--file ./plink_tmp/plink.chr-8 --clump",subset_SS,"--out ./plink_tmp")
# }



# LD_with_leadSNP <- function(LD_matrix,
#                             LD_SNP){
#   printer("LD Matrix dimensions:", paste(dim(LD_matrix), collapse=" x "))
#   printer("Extracting LD subset for lead SNP:",LD_SNP)
#   LD_sub <- subset(LD_matrix, select=LD_SNP) %>%
#     data.table::as.data.table(keep.rownames = T) %>%
#     `colnames<-`(c("SNP","r")) %>%
#     dplyr::mutate(r2 = r^2)
#   return(LD_sub)
# }
#
#


#' Find correlates of the lead GWAS/QTL SNP
LD.get_lead_r2 <- function(finemap_dat,
                           LD_matrix=NULL,
                           fillNA=0,
                           LD_format="matrix",
                           verbose=T){
  if(any(c("r","r2") %in% colnames(finemap_dat)) ){
    finemap_dat <- dplyr::select(finemap_dat, -c(r,r2))
  }
  LD_SNP <- unique(subset(finemap_dat, leadSNP==T)$SNP)
  if(length(LD_SNP)>1){
    LD_SNP <- LD_SNP[1]
    warning("More than one lead SNP found. Using only the first one:",LD_SNP,v=verbose)
  }
  # Infer LD data format
  if(LD_format=="guess"){
    LD_format <- if(nrow(LD_matrix)==ncol(LD_matrix) | class(LD_matrix)[1]=="dsCMatrix") "matrix" else "df"
  }

  if(LD_format=="matrix"){
    if(is.null(LD_matrix)){
      printer("+ LD:: No LD_matrix detected. Setting r2=NA",v=verbose);
      dat <- finemap_dat
      dat$r2 <- NA
    } else {
      printer("+ LD:: LD_matrix detected. Coloring SNPs by LD with lead SNP.", v=verbose)
      LD_sub <- LD_matrix[,LD_SNP] %>%
        # subset(LD_matrix, select=LD_SNP) %>%
        # subset(select = -c(r,r2)) %>%
        data.table::as.data.table(keep.rownames = T) %>%
        `colnames<-`(c("SNP","r")) %>%
        dplyr::mutate(r2 = r^2) %>%
        data.table::as.data.table()
      dat <- data.table::merge.data.table(finemap_dat, LD_sub,
                                          by = "SNP",
                                          all.x = T)
    }
  }
  if(LD_format=="df"){
    LD_sub <- subset(LD_matrix, select=c("SNP",LD_SNP)) %>%
      `colnames<-`(c("SNP","r")) %>%
      dplyr::mutate(r2 = r^2) %>%
      data.table::as.data.table()
    dat <- data.table::merge.data.table(finemap_dat, LD_sub,
                                        by = "SNP",
                                        all.x = T)
  }

  if(fillNA!=F){
    dat$r <- tidyr::replace_na(dat$r, fillNA)
    dat$r2 <- tidyr::replace_na(dat$r2, fillNA)
  }
  return(dat)
}




#' Extract LD proxies from 1KGphase3
#'
#' Wrapper for \code{LDlinkR::LDproxy_batch}.
#' Eeasy to use but doesn't scale up well to many SNPs (takes way too long).
#' @family LD
#' @source
#' \href{https://www.rdocumentation.org/packages/LDlinkR/versions/1.0.2}{website}
#' @examples
#' data("merged_DT")
#' lead.snps <- setNames(subset(merged_DT, leadSNP)$Locus, subset(merged_DT, leadSNP)$SNP)
#' proxies <- LDlinkR.LDproxy_batch(snp=lead.snps)
LDlinkR.LDproxy_batch <- function(snp,
                                  pop="CEU",
                                  r2d = "r2",
                                  min_corr=F,
                                  save_dir=NULL,
                                  verbose=T){
  printer("LD:LDlinkR:: Retrieving proxies of",length(snp),"SNPs",v=verbose)
  res <- LDlinkR::LDproxy_batch(snp = snp,
                                pop = pop,
                                r2d = r2d,
                                append = T,
                                token = "df4298d58dc4")
  printer("+ LD:LDlinkR::",length(unique(res$RS_Number)),"unique proxies returned.",v=verbose)
  if(min_corr!=F){
    res <- subset(res, eval(parse(text = toupper(r2d)))>=min_corr)
    printer("+ LD:LDlinkR::",length(unique(res$RS_Number)),"remaining at",r2d," ≥",min_corr,v=verbose)
  }
  proxy_files <- 'combined_query_snp_list.txt' # Automatically named
  if(!is.null(save_dir)){
    # LDproxy_batch() saves all results  as individual .txt files in the cwd by default.
    ## It's pretty dumb that they don't let you control if and where these are saved,
    ## so we have to do this manually afterwards.
    # local_dir <- "~/Desktop"
    # proxy_files <- list.files(path= "./",
    #                           pattern = "^rs.*\\.txt$", full.names = T)
    new_path <- file.path(save_dir,basename(proxy_files))
    out <- file.rename(proxy_files, new_path)
    return(new_path)
  }else{return(proxy_files)}
}

