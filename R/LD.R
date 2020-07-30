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
                              superpopulation="EUR",
                              download_reference=T,
                              download_method="direct",
                              vcf_folder=NULL,
                              min_r2=0,
                              LD_block=F,
                              LD_block_size=.7,
                              min_Dprime=F,
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
    printer("+  LD:: Previously computed LD_matrix detected. Importing...", RDS_path, v=verbose)
    LD_matrix <- readRDS(RDS_path)
  } else if(LD_reference=="UKB"){
    LD_matrix <- LD.UKBiobank(subset_DT = subset_DT,
                              locus_dir = locus_dir,
                              force_new_LD = force_new_LD,
                              chimera = server,
                              download_method = download_method,
                              nThread = nThread,
                              return_matrix = T,
                              conda_env = conda_env,
                              remove_tmps = remove_tmps)
  } else if (LD_reference == "1KGphase1" |
             LD_reference == "1KGphase3") {
      LD_matrix <- LD.1KG(locus_dir = locus_dir,
                          subset_DT = subset_DT,
                          vcf_folder = vcf_folder,
                          LD_reference = LD_reference,
                          superpopulation = superpopulation,
                          download_reference = download_reference,

                          min_r2 = min_r2,
                          LD_block = LD_block,
                          LD_block_size = LD_block_size,
                          min_Dprime = min_Dprime,
                          remove_correlates = remove_correlates,
                          fillNA = fillNA,
                          nThread = nThread,
                          conda_env = conda_env,
                          download_method = download_method)
  } else if (endsWith(tolower(LD_reference),".vcf") |
             endsWith(tolower(LD_reference),".vcf.gz")){
    LD_matrix <- LD.custom_panel(LD_reference=LD_reference,
                                 subset_DT=subset_DT,
                                 locus_dir=locus_dir,
                                 min_r2=min_r2,
                                 min_Dprime=min_Dprime,
                                 remove_correlates=remove_correlates,
                                 fillNA=fillNA,
                                 LD_block=LD_block,
                                 LD_block_size=LD_block_size,
                                 remove_tmps=remove_tmps,
                                 nThread=nThread,
                                 conda_env=conda_env,
                                 verbose=verbose)
  } else {
    stop("LD:: LD_reference input not recognized. Please supply: '1KGphase1', '1KGphase3', 'UKB', or a custom .vcf[.gz] file.")
  }
  return(LD_matrix)
}




LD.get_rds_path <- function(locus_dir,
                            LD_reference){
  RDS_path <- file.path(locus_dir,"LD",paste0(basename(locus_dir),".",basename(LD_reference),"_LD.RDS"))
  return(RDS_path)
}




#' Compute LD from user-supplied vcf file
#'
#' @family LD
#' @keywords internal
#' @examples
#' data("BST1"); data("locus_dir")
#' LD_reference="~/Desktop/results/Reference/custom_panel_chr4.vcf"
#' LD_matrix <- LD.custom_panel(LD_reference=LD_reference, subset_DT=BST1, locus_dir=locus_dir)
LD.custom_panel <- function(LD_reference,
                            subset_DT,
                            locus_dir,
                            min_r2=F,
                            min_Dprime=F,
                            remove_correlates=F,
                            fillNA=0,
                            LD_block=F,
                            LD_block_size=.7,
                            remove_tmps=T,
                            nThread=4,
                            conda_env="echoR",
                            verbose=T){
  printer("LD:: Computing LD from local vcf file:",LD_reference)
  vcf_file <- LD.index_vcf(vcf_file=LD_reference)

  vcf_subset <- LD.query_vcf(subset_DT=subset_DT,
                             locus_dir=locus_dir,
                             LD_reference=LD_reference,
                             vcf_URL=LD_reference,
                             whole_vcf=F,
                             remove_original_vcf=F,
                             force_new_vcf=F,
                             query_by_regions=F,
                             nThread=nThread,
                             conda_env=conda_env,
                             verbose=verbose)
  LD.vcf_to_bed(vcf.gz.subset = vcf_subset,
                locus_dir = locus_dir)
  # Get lead SNP rsid
  leadSNP = subset(subset_DT, leadSNP==T)$SNP
  # Plink LD method
  LD_matrix <- LD.plink_LD(LD_folder = file.path(locus_dir,"LD"),
                           subset_DT = subset_DT,
                           remove_excess_snps=T,
                           leadSNP = leadSNP,
                           min_r2 = min_r2,
                           min_Dprime = min_Dprime,
                           remove_correlates = remove_correlates,
                           fillNA = fillNA)
  # Filter out SNPs not in the same LD block as the lead SNP
  if(LD_block){
    block_snps <- LD.leadSNP_block(leadSNP, "./plink_tmp", LD_block_size)
    LD_matrix <- LD_matrix[row.names(LD_matrix) %in% block_snps, colnames(LD_matrix) %in% block_snps]
    LD_matrix <- LD_matrix[block_snps, block_snps]
  }
  # IMPORTANT! Remove large data.ld file after you're done with it
  if(remove_tmps){
    suppressWarnings(file.remove(vcf_subset))
  }
  printer("Saving LD matrix of size:", dim(LD_matrix)[1],"rows x",dim(LD_matrix)[2],"columns.", v=verbose)
  return(LD_matrix)
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
#' LD.plot(LD_matrix=LD_matrix, subset_DT=BST1)
#' }
LD.plot <- function(LD_matrix,
                    subset_DT,
                    span=10){
  leadSNP = subset(subset_DT, leadSNP==T)$SNP
  lead_index = match(leadSNP, row.names(LD_matrix))

  if(dim(LD_matrix)[1]<span){
    start_pos = lead_index - dim(LD_matrix)[1]
    end_pos = lead_index + dim(LD_matrix)[1]
  } else{
    start_pos = lead_index - span
    end_pos = lead_index + span
  }
  sub_DT <- subset(subset_DT, SNP %in% rownames(LD_matrix))
  gaston::LD.plot( LD_matrix[start_pos:end_pos, start_pos:end_pos],
                   snp.positions = sub_DT$POS[start_pos:end_pos] )
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
                                vcf_folder=NULL,
                                locus_dir,
                                locus=NULL,
                                download_reference=T,
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
  vcf_folder <- LD.get_locus_vcf_folder(locus_dir=locus_dir)
  # Don't use the chr prefix: https://www.internationalgenome.org/faq/how-do-i-get-sub-section-vcf-file/
  subset_DT$CHR <- gsub("chr","",subset_DT$CHR)
  chrom <- unique(subset_DT$CHR)

  # PHASE 3 DATA
  if(LD_reference=="1KGphase3"){
    FTP <- "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/"
    printer("LD Reference Panel = 1KGphase3", v=verbose)
    if(download_reference){## With internet
      vcf_URL <- paste0(FTP,"/ALL.chr",chrom,
                       ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
      popDat_URL = paste0(FTP, "integrated_call_samples_v3.20130502.ALL.panel")
    }else{## WithOUT internet
      vcf_URL <- paste(vcf_folder, "/ALL.chr",chrom,
                       ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",sep="")
      popDat_URL = file.path(vcf_folder,"integrated_call_samples_v3.20130502.ALL.panel")
    }

    # PHASE 1 DATA
  } else if (LD_reference=="1KGphase1") {
    FTP <- "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/"

    printer("LD Reference Panel = 1KGphase1", v=verbose)

    if(download_reference){## With internet
      vcf_URL <- paste(FTP,"/ALL.chr",chrom,
                       ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")
      # popDat_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel"
      popDat_URL <- file.path(FTP,"phase1_integrated_calls.20101123.ALL.panel")
    }else{## WithOUT internet
      vcf_URL <- paste(vcf_folder,"/ALL.chr",chrom,
                       ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")
      popDat_URL = file.path(vcf_folder, "phase1_integrated_calls.20101123.ALL.panel")
    }
  }
  printer(paste("Reading population data from", popDat_URL), v=verbose)
  phase <- gsub("1KG","",LD_reference)
  # phase 1 has no header whereas phase 3 does
  if( LD_reference == "1KGphase1" ){ use_header <- FALSE }
  if( LD_reference == "1KGphase3" ){ use_header <- TRUE }
  popDat <-  data.table::fread(text=trimws(gsub(",\t",",",readLines(popDat_URL))),
                               header = use_header, sep="\t",  fill=T, stringsAsFactors = FALSE,
                               col.names = c("sample","population","superpop","platform"),
                               nThread = nThread)
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
#' LD_reference <- "~/Desktop/results/Reference/custom_panel_chr4.vcf.gz"
#' vcf_file <- LD.index_vcf(vcf_file=LD_reference)
LD.index_vcf <- function(vcf_file,
                         force_new_index=T,
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
#' data("locus_dir"); data("BST1");
#'
#' # Custom
#' LD_reference <- "~/Desktop/results/Reference/custom_panel_chr4.vcf"
#' vcf_file <- LD.index_vcf(vcf_file=LD_reference)
#' vcf_subset <- LD.query_vcf(subset_DT=BST1, locus_dir=locus_dir, vcf_URL=vcf_file, LD_reference=LD_reference, force_new_vcf=T)
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
  vcf_subset <- LD.construct_subset_vcf_name(subset_DT=subset_DT,
                                             locus_dir=locus_dir,
                                             LD_reference=LD_reference,
                                             whole_vcf=whole_vcf)
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
                             nThread = nThread)
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




#' Filter a vcf by min/max coordinates
#'
#' Uses \emph{bcftools} to filter a vcf by min/max genomic coordinates (in basepairs).
#' @param vcf_subset Path to the locus subset vcf.
#' @param popDat The metadata file listing the superpopulation to which each sample belongs.
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
                          verbose=T){
  plink <- LD.plink_file()
  printer("LD:PLINK:: Converting vcf.gz to .bed/.bim/.fam", v=verbose)
  dir.create(file.path(locus_dir,"LD"), recursive = T, showWarnings = F)
  cmd <- paste(plink,
               "--vcf",vcf.gz.subset,
               "--out", file.path(locus_dir,"LD","plink"))
  system(cmd)
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
                            verbose=T){
  plink <- LD.plink_file()
  printer("LD:PLINK:: Calculating LD ( r & D'-signed; LD-window =",ld_window,")", v=verbose)
  plink_path_prefix <- file.path(locus_dir,"LD","plink")
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
#' @param ld.path Directory that contains the bin/bim files.
#' @family LD
#' @keywords internal
LD.read_bin <- function(ld.path){
  bim <- data.table::fread(file.path(dirname(ld.path), "plink.bim"), col.names = c("CHR","SNP","V3","POS","A1","A2"))
  bin.vector <- readBin(ld.path, what = "numeric", n=length(bim$SNP)^2)
  ld.matrix <- matrix(bin.vector, nrow = length(bim$SNP), dimnames = list(bim$SNP, bim$SNP))
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
LD.plink_file <- function(base_url=system.file("tools/plink",package = "echolocatoR")){
  os <- get_os()
  if (os=="osx") {
    plink <- file.path(base_url, "plink1.9_mac");
  } else if (os=="linux") {
    plink <- file.path(base_url, "plink1.9_linux");
  } else {
    plink <- file.path(base_url, "plink1.9_windows.exe");
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
#' }
LD.1KG <- function(locus_dir,
                   subset_DT,
                   LD_reference="1KGphase1",
                   superpopulation="EUR",
                   vcf_folder=NULL,
                   download_reference=T,
                   min_r2=F,
                   LD_block=F,
                   LD_block_size=.7,
                   min_Dprime=F,
                   remove_correlates=F,
                   remove_tmps=T,
                   fillNA=0,
                   download_method="wget",
                   nThread=4,
                   conda_env="echoR",
                   verbose=T){
  # data("BST1"); data("locus_dir"); subset_DT=BST1; LD_reference="1KGphase1"; vcf_folder=NULL; superpopulation="EUR";
  # min_r2=F; LD_block=F; LD_block_size=.7; min_Dprime=F;  remove_correlates=F; download_reference=T; verbose=T; nThread=4; conda_env="echoR";
  printer("LD:: Using 1000Genomes as LD reference panel.", v=verbose)
  locus <- basename(locus_dir)
  vcf_folder <- LD.get_locus_vcf_folder(locus_dir = locus_dir)
  vcf_info <- LD.1KG_download_vcf(subset_DT=subset_DT,
                                  locus_dir=locus_dir,
                                  LD_reference=LD_reference,
                                  vcf_folder=vcf_folder,
                                  locus=locus,
                                  download_reference=download_reference,
                                  download_method=download_method,
                                  nThread=nThread,
                                  conda_env=conda_env,
                                  verbose=verbose)
  vcf_subset <- vcf_info$vcf_subset
  popDat <- vcf_info$popDat
  vcf.gz.path <- LD.filter_vcf(vcf_subset = vcf_subset,
                               popDat = popDat,
                               superpopulation = superpopulation,
                               remove_tmp = T)
  LD.vcf_to_bed(vcf.gz.subset = vcf.gz.path,
                locus_dir = locus_dir)
  # Calculate pairwise LD for all SNP combinations
  #### "Caution that the LD matrix has to be correlation matrix" -SuSiER documentation
  ### https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html
  # Gaston LD method
  # LD_matrix <- gaston::LD(bed, lim = c(1,ncol(bed)), measure ="r") #"D"
  # LD_matrix[!is.finite(LD_matrix)] <- 0

  # Get lead SNP rsid
  leadSNP = subset(subset_DT, leadSNP==T)$SNP #rs76904798
  # Plink LD method
  LD_matrix <- LD.plink_LD(LD_folder = file.path(locus_dir,"LD"),
                           subset_DT = subset_DT,
                           remove_excess_snps=T,
                           leadSNP = leadSNP,
                           min_r2 = min_r2,
                           min_Dprime = min_Dprime,
                           remove_correlates = remove_correlates,
                           fillNA = fillNA)
  # Filter out SNPs not in the same LD block as the lead SNP
  if(LD_block){
    block_snps <- LD.leadSNP_block(leadSNP, "./plink_tmp", LD_block_size)
    LD_matrix <- LD_matrix[row.names(LD_matrix) %in% block_snps, colnames(LD_matrix) %in% block_snps]
    LD_matrix <- LD_matrix[block_snps, block_snps]
  }
  # IMPORTANT! Remove large data.ld file after you're done with it
  if(remove_tmps){
    suppressWarnings(file.remove(vcf_subset))
  }
  printer("Saving LD matrix of size:", dim(LD_matrix)[1],"rows x",dim(LD_matrix)[2],"columns.", v=verbose)

  # Save LD matrix
  RDS_path <- LD.get_rds_path(locus_dir = locus_dir,
                              LD_reference = LD_reference)
  printer("+ LD:: Saving LD_matrix to:",RDS_path, v=verbose)
  saveRDS(LD_matrix, file = RDS_path)
  return(LD_matrix)
}




#' Calculate LD (D')
#'
#' This appriach computes an LD matrix of D' (instead of r or r2) from a vcf.
#' See \code{\link{LD.run_plink_LD}} for a faster (but less flexible) alternative to computing LD.
#' @family LD
#' @keywords internal
LD.dprime_table <- function(SNP_list, LD_folder){
  plink <- LD.plink_file()
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




#' Calculate LD (r or r2)
#'
#' This appriach computes and LD matrix of r or r2 (instead of D') from a vcf.
#' See \code{\link{LD.dprime_table}} for a slower (but more flexible) alternative to computing LD.
#' @param bim A bim file produced by \emph{plink}
#' @param LD_folder Locus-specific LD output folder.
#' @param r_format Whether to fill the matrix with \code{r} or \code{r2}.
#' @family LD
#' @keywords internal
LD.run_plink_LD <- function(bim,
                            LD_folder,
                            r_format="r"){
  plink <- LD.plink_file()
  # METHOD 2 (faster, but less control over parameters. Most importantly, can't get Dprime)
  system( paste(plink,
                "--bfile",file.path(LD_folder,"plink"),
                "--extract",file.path(LD_folder,"SNPs.txt"),
                paste0("--",r_format," square bin"), "--out", file.path(LD_folder,"plink")) )
  bin.vector <- readBin(file.path(LD_folder, "plink.ld.bin"), what = "numeric", n=length(bim$SNP)^2)
  ld.matrix <- matrix(bin.vector, nrow = length(bim$SNP), dimnames = list(bim$SNP, bim$SNP))
  return(ld.matrix)
}




#' Calculate LD
#'
#' Use \emph{plink} to calculate LD from a vcf.
#' @family LD
#' @keywords internal
LD.plink_LD <-function(leadSNP,
                       subset_DT,
                       remove_excess_snps=T,
                       LD_folder,
                       min_r2=F,
                       min_Dprime=F,
                       remove_correlates=F,
                       fillNA=0,
                       verbose=T){
  # Dprime ranges from -1 to 1
  start <- Sys.time()
  # Calculate LD
  printer("++ Reading in BIM file...", v=verbose)
  bim <- data.table::fread(file.path(LD_folder, "plink.bim"), col.names = c("CHR","SNP","V3","POS","A1","A2"))
  if(remove_excess_snps){
    orig_n <- nrow(bim)
    bim.merged <- data.table::merge.data.table(bim,
                                 subset_DT,
                                 by=c("CHR","POS"))
    bim <- subset(bim, SNP %in% bim.merged$SNP.x)
    printer("LD:PLINK:: Removing RSIDs that don't appear in locus subset:",orig_n,"==>",nrow(bim),"SNPs",v=verbose)
  }
  data.table::fwrite(subset(bim, select="SNP"), file.path(LD_folder,"SNPs.txt"), col.names = F)

  printer("++ Calculating LD", v=verbose)
  ld.matrix <- LD.run_plink_LD(bim, LD_folder)

  if((min_Dprime != F) | (min_r2 != F) | (remove_correlates != F)){
    plink.ld <- LD.dprime_table(SNP_list = row.names(ld.matrix), LD_folder)

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
  printer("LD:: Removing unnamed rows/cols", v=verbose)
  # First, filter any rows/cols without names
  LD_matrix <- data.frame(LD_matrix)
  LD_matrix <- LD_matrix[rownames(LD_matrix)!=".", colnames(LD_matrix)!="."]
  LD_matrix_orig <- LD_matrix
  LD_matrix <- data.table::data.table(LD_matrix)

  if(!is.null(fillNA)){
    # Replace NAs
    printer("LD:: Replacing NAs with",fillNA, v=verbose)
    replace_NAs = function(DT) {
      # Does inplace
      for (i in names(DT))
        DT[is.na(get(i)), (i):=0]
    }
    replace_NAs(LD_matrix)
  }
  # Check for duplicate SNPs
  LD_fill <- data.frame(LD_matrix,
                        row.names = rownames(LD_matrix_orig))
  if(dplyr::n_distinct(rownames(LD_fill))!=nrow(LD_fill) |
     dplyr::n_distinct(colnames(LD_fill))!=ncol(LD_fill)){
    stop("There's more rows/cols than unique SNPs.")
  }
  return(LD_fill)
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





