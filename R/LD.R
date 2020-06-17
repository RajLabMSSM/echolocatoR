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
#' \item Compute LD fon the fly from a user-supplied vcf file.
#' \item Use a user-supplied pre-computed LD-matrix.
#' }
#'
#' @param subset_DT The locus subset of the full summary stats file.
#' @inheritParams finemap_pipeline
#' @return A symmetric LD matrix of pairwise \emph{r} values.
LD.load_or_create <- function(subset_path,
                              subset_DT,
                              locus,
                              force_new_LD=F,
                              LD_reference="1KG_Phase1",
                              superpopulation="EUR",
                              download_reference=T,
                              download_method="direct",
                              min_r2=0,
                              LD_block=F,
                              LD_block_size=.7,
                              min_Dprime=F,
                              remove_correlates=F,
                              fillNA=0,
                              verbose=T,
                              server=F,
                              remove_tmps=T,
                              nThreads=4){
  locus_dir <- get_locus_dir(subset_path=subset_path)
  if(LD_reference=="UKB"){
    printer("LD:: Using UK Biobank LD reference panel...")
    LD_matrix <- LD.UKBiobank(subset_DT = subset_DT,
                              locus_dir = locus_dir,
                               force_new_LD = force_new_LD,
                               chimera = server,
                               download_method = download_method,
                               nThreads = nThreads,
                               return_matrix = T,
                               remove_tmps = remove_tmps)
  } else if (LD_reference == "1KG_Phase1" | LD_reference == "1KG_Phase3") {
    LD_path <- file.path(locus_dir,"LD",paste0(LD_reference,"_LD.RDS"))
    if(!file.exists(LD_path) | force_new_LD==T){
      printer("+ Computing LD matrix...", verbose)
      LD_matrix <- LD.1KG(locus_dir = locus_dir,
                           subset_DT = subset_DT,
                           locus = locus,
                           LD_reference = LD_reference,
                           superpopulation = superpopulation,
                           download_reference = download_reference,

                           min_r2 = min_r2,
                           LD_block = LD_block,
                           LD_block_size = LD_block_size,
                           min_Dprime = min_Dprime,
                           remove_correlates = remove_correlates,
                           fillNA = fillNA)
      # Save LD matrix
      printer("+ Saving LD matrix to:",LD_path, v=verbose)
      saveRDS(LD_matrix, file = LD_path)
    } else {
      printer("+ Previously computed LD matrix detected. Importing...", LD_path, v=verbose)
      LD_matrix <- readRDS(LD_path)
    }
  }
  return(LD_matrix)
}




#' Translate superopulation acronyms
#'
#' Ensures a common ontology for synonynmous superpopulation names.
#'
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
#' @param span This is very computationally intensive,
#' so you need to limit the number of SNPs with span.
#' If \code{span=10}, only 10 SNPs upstream and 10 SNPs downstream of the lead SNP will be plotted.
LD.plot <- function(LD_matrix,
                    subset_DT,
                    span=10){
  leadSNP = subset(subset_DT, leadSNP==T)$SNP
  lead_index = match(leadSNP, row.names(LD_matrix))
  if(dim(LD_matrix)[1]<span){
    start = lead_index - dim(LD_matrix)[1]
    end = lead_index + dim(LD_matrix)[1]
  } else{
    start = lead_index - span
    end = lead_index + span
  }
  sub_DT <- subset(subset_DT, SNP %in% rownames(LD_matrix))
  gaston::LD.plot( LD_matrix[start:end, start:end], snp.positions = sub_DT$POS[start:end] )
}




#' Download vcf subset from 1000 Genomes
#'
#' @inheritParams finemap_pipeline
LD.1KG_download_vcf <- function(subset_DT,
                                LD_reference="1KG_Phase1",
                                vcf_folder="./Data/Reference/1000_Genomes",
                                locus_dir,
                                locus,
                                download_reference=T,
                                whole_vcf=F){
  ## http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
  # Download portion of vcf from 1KG website

  # Don't use the chr prefix: https://www.internationalgenome.org/faq/how-do-i-get-sub-section-vcf-file/
  subset_DT$CHR <- gsub("chr","",subset_DT$CHR)
  chrom <- unique(subset_DT$CHR)

  # PHASE 3 DATA
  if(LD_reference=="1KG_Phase3"){
    printer("LD Reference Panel = 1KG_Phase3")
    if(download_reference){## With internet
      vcf_URL <- paste("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr",chrom,
                       ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",sep="")
      popDat_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
    }else{## WithOUT internet
      vcf_URL <- paste(vcf_folder, "/ALL.chr",chrom,
                       ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",sep="")
      popDat_URL = file.path(vcf_folder,"integrated_call_samples_v3.20130502.ALL.panel")
    }

    # PHASE 1 DATA
  } else if (LD_reference=="1KG_Phase1") {
    printer("LD Reference Panel = 1KG_Phase1")
    if(download_reference){## With internet
      vcf_URL <- paste("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr",chrom,
                       ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")
      popDat_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel"
    }else{## WithOUT internet
      vcf_URL <- paste(vcf_folder,"/ALL.chr",chrom,
                       ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")
      popDat_URL = file.path(vcf_folder, "phase1_integrated_calls.20101123.ALL.panel")
    }
  }
  phase <- gsub("1KG_","",LD_reference)
  popDat <-  data.table::fread(text=gsub(",\t",",",readLines(popDat_URL)),
                               header = F, sep="\t",  fill=T, stringsAsFactors = F,
                               col.names = c("sample","population","superpop","platform"))
  # library(Rsamtools); #BiocManager::install("Rsamtools")
  subset_vcf <- file.path(vcf_folder, phase,
                          ifelse(whole_vcf, paste(LD_reference, paste0("chr",chrom),"vcf", sep="."),
                                                 paste(locus,"subset.vcf",sep="_")))
  # Create directory if it doesn't exist
  dir.create(path = dirname(subset_vcf), recursive = T, showWarnings = F)
  # Download and subset vcf if the subset doesn't exist already
  if(!file.exists(subset_vcf)){
    if(whole_vcf){
      region <= ""
      locus=""
    } else {
      # Specify exact positions
      regions.bed <- file.path(locus_dir,"LD","regions.tsv")
      data.table::fwrite(list(paste0(subset_DT$CHR), sort(subset_DT$POS)),
                         file=regions.bed, sep="\t")
      # data.table::fwrite(list(paste0(subset_DT$CHR,":", subset_DT$POS,"-",subset_DT$POS)),
      #                    file=regions.bed, sep="\t")
      region <- paste("-R",regions.bed)
    }
    # Download tabix subset
    tabix_cmd <- paste("tabix -fh",vcf_URL,
                       # Specifying the SNP subset drastically reduces compute time
                       gsub("\\./","",region),
                       ">",
                       gsub("\\./","",subset_vcf) )
    printer(tabix_cmd)
    # system("ml tabix")
    system(tabix_cmd)
    vcf_name <- paste(basename(vcf_URL), ".tbi", sep="")
    file.remove(vcf_name)
  } else {printer("+ Identified matching VCF subset file. Importing...", subset_vcf)}
  return(list(subset_vcf = subset_vcf,
              popDat = popDat))
}






#' Filter a vcf by min/max coordinates
#'
#' Uses the \code{Rsamtools} package to filter a vcf fiile.
#' @keywords internal
LD.filter_vcf_rsamtools <- function(){
  gr <- GenomicRanges::makeGRangesFromDataFrame(df = subset_DT[,c("CHR","POS")],
                                                seqnames.field = "CHR",
                                                start.field = "POS",
                                                end.field = "POS")
  # Rsamtools::scanTabix(file = vcf_URL, param = gr)
  # library(Rsamtools)
  tbx <- Rsamtools::TabixFile(file = vcf_URL,  index = paste0(vcf_URL,".tbi"))
  countTabix(tbx, param=gr)
}


#' Filter a vcf by min/max coordinates
#'
#' Uses \emph{bcftools} to filter a vcf by min/max genomic coordinates (in basepairs).
#' @param subset_vcf Path to the locus subset vcf.
#' @param popDat The metadata file listing the superpopulation to which each sample belongs.
#' @inheritParams finemap_pipeline
LD.filter_vcf <- function(subset_vcf,
                          popDat,
                          superpopulation,
                          remove_tmp=F){
  # vcf.bgz <- Rsamtools::bgzip(subset_vcf, overwrite = T)
  # Rsamtools::indexTabix(file = vcf.bgz, seq = 12, begin=subset_DT$POS[1])
  # tbx <- Rsamtools::TabixFile(file = subset_vcf)

  vcf.gz <- paste0(subset_vcf,".gz")
  vcf.gz.subset <- gsub("_subset","_samples_subset",vcf.gz)
  # Compress vcf
  if(!file.exists(vcf.gz)){
    printer("LD:BCFTOOLS:: Compressing vcf file...")
    system(paste("bgzip -f",subset_vcf))
  }
  # Re-index vcf
  printer("LD:TABIX:: Re-indexing vcf.gz...")
  system(paste("tabix -f -p vcf",vcf.gz))

  # Subset samples
  selectedInds <- subset(popDat, superpop == superpopulation)$sample %>% unique()
  printer("LD:BCFTOOLS:: Subsetting vcf to only include",superpopulation,"individuals (",length(selectedInds), "/",length(popDat$sample%>%unique()),").")
  cmd <- paste("bcftools view -s",paste(selectedInds, collapse=","), vcf.gz, "| bgzip > tmp && mv tmp",vcf.gz.subset)
  system(cmd)
  # mega_cmd <- paste(subset_vcf,
  #                   "| bgzip",
  #                   "| tabix -f -p vcf",
  #                   "| bcftools view -s",paste(selectedInds, collapse=","), vcf.gz,
  #                   "| bgzip > tmp && mv tmp",vcf.gz.subset)
  return(vcf.gz.subset)
}




#' Subset a vcf by superpopulation
#'
#' @inheritParams LD.filter_vcf
LD.filter_vcf_population <- function(subset_vcf,
                                     subset_DT,
                                     locus_dir,
                                     superpopulation,
                                     popDat){
  # Import w/ gaston and further subset
  printer("+ Importing VCF as bed file...")
  bed.file <- gaston::read.vcf(subset_vcf, verbose = F)
  ## Subset rsIDs
  bed <- gaston::select.snps(bed.file, id %in% subset_DT$SNP & id !=".")
  # Create plink sub-dir
  dir.create(file.path(locus_dir, "LD"), recursive = T, showWarnings = F)
  gaston::write.bed.matrix(bed, file.path(locus_dir, "plink/plink"), rds = NULL)
  # Subset Individuals
  selectedInds <- subset(popDat, superpop == superpopulation)
  bed <- gaston::select.inds(bed, id %in% selectedInds$sample)
  # Cleanup extra files
  remove(bed.file)
  # file.remove("subset_vcf")
  return(bed)
}




#' Convert vcf file to BED file
#'
#' Uses plink to convert vcf to BED.
#' @param vcf.gz.subset Path to the gzipped locus subset vcf.
#' @param locus_dir Locus-specific results directory.
LD.vcf_to_bed <- function(vcf.gz.subset,
                          locus_dir){
  printer("LD:PLINK:: Converting vcf.gz to .bed/.bim/.fam")
  dir.create(file.path(locus_dir,"LD"), recursive = T, showWarnings = F)
  cmd <- paste("plink","--vcf",vcf.gz.subset, "--out", file.path(locus_dir,"LD","LD"))
  system(cmd)
}



#' Calculate LD
#'
#' Calculate a pairwise LD matrix from a vcf file using \emph{plink}.
#' @param locus_dir Locus-specific results directory.
#' @param ld_window Set --r/--r2 max variant ct pairwise distance (usu. 10).
#' @param ld_format Whether to produce an LD matrix with
#' r (\code{ld_format="r"}) or D' (\code{ld_format="D"}) as the pairwise SNP correlation metric.
LD.calculate_LD <- function(locus_dir,
                            ld_window=1000, # 10000000
                            ld_format="r"){
  printer("LD:PLINK:: Calculating LD ( r & D'-signed; LD-window =",ld_window,")")
  plink_path_prefix <- file.path(locus_dir,"LD","LD")
  dir.create(file.path(locus_dir,"LD"), recursive = T, showWarnings = F)
  out_prefix <- paste0(plink_path_prefix,".r_dprimeSigned")
  if(ld_format=="r"){
    cmd <- paste("plink",
                 "--bfile",plink_path_prefix,
                 "--r square bin",
                 "--out",out_prefix)
    ld.path <- paste0(out_prefix,".ld.bin")
  } else {
    cmd <- paste("plink",
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
LD.read_bin <- function(ld.path){
  bim <- data.table::fread(file.path(dirname(ld.path), "plink.bim"), col.names = c("CHR","SNP","V3","POS","A1","A2"))
  bin.vector <- readBin(ld.path, what = "numeric", n=length(bim$SNP)^2)
  ld.matrix <- matrix(bin.vector, nrow = length(bim$SNP), dimnames = list(bim$SNP, bim$SNP))
}

#' Create LD matrix from plink output.
#'
#' Depending on which parameters you give \emph{plink} when calculating LD, you get different file outputs.
#' When it produces an LD table, use this function to create a proper LD matrix.
LD.read_ld_table <- function(ld.path, snp.subset=F){
  # subset_DT <- data.table::fread(file.path(locus_dir,"Multi-finemap/Multi-finemap_results.txt")); snp.subset <- subset_DT$SNP
  ld.table <- data.table::fread(ld.path, nThread = 4)
  if(any(snp.subset!=F)){
    printer("LD:PLINK:: Subsetting LD data...")
    ld.table <- subset(ld.table, SNP_A %in% snp.subset | SNP_B %in% snp.subset)
  }
  printer("LD:PLINK:: Casting data.matrix...")
  ld.cast <- data.table::dcast.data.table(ld.table,
                                          formula = SNP_B ~ SNP_A,
                                          value.var="R",
                                          fill=0,
                                          drop=T,
                                          fun.agg = function(x){mean(x,na.rm = T)})
  ld.cast <- subset(ld.cast, SNP_B !=".", select = -`.`)
  ld.mat <- data.frame(ld.cast, row.names = ld.cast$SNP_B) %>% data.table() %>% as.matrix()
  # ld.mat[1:10,1:10]
  ld.mat[is.na(ld.mat)] <- 0
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
LD.1KG <- function(locus_dir,
                    subset_DT,
                    locus,
                    LD_reference="1KG_Phase1",
                    superpopulation="EUR",
                    vcf_folder="./Data/Reference/1000_Genomes",
                    download_reference=T,
                    min_r2=F,
                    LD_block=F,
                    LD_block_size=.7,
                    min_Dprime=F,
                    remove_correlates=F,
                    remove_tmps=T,
                    fillNA=0){
          # Quickstart
  # locus <- "LRRK2"; locus_dir <- file.path("./Data/GWAS/Nalls23andMe_2019",locus); LD_reference="1KG_Phase1"; vcf_folder="./Data/Reference/1000_Genomes"; superpopulation="EUR"; vcf_folder="./Data/Reference/1000_Genomes"; min_r2=F; LD_block=F; LD_block_size=.7; min_Dprime=F;  remove_correlates=F; download_reference=T; subset_DT <- data.table::fread(file.path(locus_dir,"Multi-finemap/Multi-finemap_results.txt"))
    printer("LD:: Using 1000Genomes LD reference panel...")
    vcf_info <- LD.1KG_download_vcf(subset_DT=subset_DT,
                                     locus_dir,
                                     LD_reference=LD_reference,
                                     vcf_folder=vcf_folder,
                                     locus=locus,
                                     download_reference=download_reference)
    subset_vcf <- vcf_info$subset_vcf
    popDat <- vcf_info$popDat
    vcf.gz.path <- BCFTOOLS.filter_vcf(subset_vcf = subset_vcf,
                                       popDat = popDat,
                                       superpopulation = superpopulation,
                                       remove_tmp = F)
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
    LD_matrix <- LD.plink_LD(plink_folder = file.path(locus_dir,"LD"),
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
      suppressWarnings(file.remove(subset_vcf))
    }
  return(LD_matrix)
  printer("Saving LD matrix of size:", dim(LD_matrix)[1],"rows x",dim(LD_matrix)[2],"columns.")
}




#' Calculate LD (D')
#'
#' This appriach computes an LD matrix of D' (instead of r or r2) from a vcf.
#' See \code{\link{LD.run_plink_LD}} for a faster (but less flexible) alternative to computing LD.
LD.dprime_table <- function(SNP_list, plink_folder){
  printer("+ Creating DPrime table")
  system( paste("plink", "--bfile",file.path(plink_folder,"LD"),
                "--ld-snps", paste(SNP_list, collapse=" "),
                "--r dprime-signed",
                "--ld-window 10000000", # max out window size
                "--ld-window-kb 10000000",
                "--out",file.path(plink_folder,"LD")) )
  #--ld-window-r2 0

  # # Awk method: theoretically faster?
  # if(min_Dprime==F){Dprime = -1}else{Dprime=min_Dprime}
  # if(min_r2==F){r = -1}else{r = round(sqrt(min_r2),2) }
  # columns <- data.table::fread(file.path(plink_folder, "plink.ld"), nrows = 0) %>% colnames()
  # col_dict <- setNames(1:length(columns), columns)
  # awk_cmd <- paste("awk -F \"\t\" 'NR==1{print $0}{ if(($",col_dict["DP"]," >= ",Dprime,")",
  #                  " && ($",col_dict["R"]," >= ",r,")) { print } }' ",file.path(plink_folder, "plink.ld"),
  #                  " > ",file.path(plink_folder, "plink.ld_filtered.txt"),  sep="")
  # system(awk_cmd)
  plink.ld <- data.table::fread(file.path(plink_folder, "plink.ld"), select = c("SNP_A", "SNP_B","DP","R"), )
  plink.ld <- plink.ld[complete.cases(plink.ld) ]
  return(plink.ld)
}




#' Calculate LD (r or r2)
#'
#' This appriach computes and LD matrix of r or r2 (instead of D') from a vcf.
#' See \code{\link{LD.dprime_table}} for a slower (but more flexible) alternative to computing LD.
#' @param bim A bim file produced by \emph{plink}
#' @param plink_folder Locus-specific LD output folder.
#' @param r_format Whether to fill the matrix with \code{r} or \code{r2}.
LD.run_plink_LD <- function(bim,
                            plink_folder,
                            r_format="r"){
  # METHOD 2 (faster, but less control over parameters. Most importantly, can't get Dprime)
  system( paste("plink", "--bfile",file.path(plink_folder,"LD"),
                "--extract",file.path(plink_folder,"SNPs.txt"),
                paste0("--",r_format," square bin"), "--out", file.path(plink_folder,"LD")) )
  bin.vector <- readBin(file.path(plink_folder, "plink.ld.bin"), what = "numeric", n=length(bim$SNP)^2)
  ld.matrix <- matrix(bin.vector, nrow = length(bim$SNP), dimnames = list(bim$SNP, bim$SNP))
  return(ld.matrix)
}




#' Calculate LD
#'
#' Use \emph{plink} to calculate LD from a vcf.
LD.plink_LD <-function(leadSNP,
                       plink_folder,
                       min_r2=F,
                       min_Dprime=F,
                       remove_correlates=F,
                       fillNA=0){
  # Dprime ranges from -1 to 1
  start <- Sys.time()

  # Calculate LD
  printer("++ Reading in BIM file...")
  bim <- data.table::fread(file.path(plink_folder, "plink.bim"), col.names = c("CHR","SNP","V3","POS","A1","A2"))
  data.table::fwrite(subset(bim, select="SNP"), file.path(plink_folder,"SNPs.txt"), col.names = F)

  printer("++ Calculating LD")
  ld.matrix <- run_plink_LD(bim, plink_folder)

  if((min_Dprime != F) | (min_r2 != F) | (remove_correlates != F)){
    plink.ld <- LD.dprime_table(SNP_list = row.names(ld.matrix), plink_folder)

    # DPrime filter
    if(min_Dprime != F){
      printer("+++ Filtering LD Matrix (min_Dprime): Removing SNPs with D' <=",min_Dprime,"for",leadSNP,"(lead SNP).")
      plink.ld <- subset(plink.ld, (SNP_A==leadSNP & DP>=min_Dprime) | (SNP_B==leadSNP & DP>=min_Dprime))
    } else{printer("+ min_Dprime == FALSE")}

    # R2 filter
    if(min_r2 != F ){
      printer("+++ Filtering LD Matrix (min_r2): Removing SNPs with r <=",min_r2,"for",leadSNP,"(lead SNP).")
      r = sqrt(min_r2) # PROBLEM: this doesn't give you r, which you need for SUSIE
      plink.ld <- subset(plink.ld, (SNP_A==leadSNP & R>=r) | (SNP_B==leadSNP & R>=r))
    } else{printer("+ min_r2 == FALSE")}

    # Correlates filter
    if(remove_correlates != F){
      r2_threshold <- remove_correlates# 0.2
      r <- sqrt(r2_threshold)
      printer("+++ Filtering LD Matrix (remove_correlates): Removing SNPs with R2 >=",r2_threshold,"for",paste(remove_correlates,collapse=", "),".")
      plink.ld <- subset(plink.ld, !(SNP_A %in% remove_correlates & R>=r) | (SNP_B %in% remove_correlates & R>=r))
    } else{printer("+ remove_correlates == FALSE")}

    # Apply filters
    A_list <- unique(plink.ld$SNP_A)
    B_list <- unique(plink.ld$SNP_B)
    snp_list <-   unique(c(A_list, B_list))
    ld.matrix <- ld.matrix[row.names(ld.matrix) %in% snp_list, colnames(ld.matrix) %in% snp_list]
    ## Manually remove rare variant
    # ld.matrix <- ld.matrix[rownames(ld.matrix)!="rs34637584", colnames(ld.matrix)!="rs34637584"]
    }
    # !IMPORTANT!: Fill NAs (otherwise susieR will break)
    ld.matrix[is.na(ld.matrix)] <- fillNA
    end <- Sys.time()
    printer("+ LD matrix calculated in",round(as.numeric(end-start),2),"seconds.")
    return(ld.matrix)
}




#' Calculate LD blocks.
#'
#' Uses \emph{plink} to group highly correlated SNPs into LD blocks.
#'
LD.LD_blocks <- function(plink_folder,
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

  # Reucing "--blocks-inform-frac" is the only parameter that seems to make the block sizes larger
  system( paste("plink", "--bfile",file.path(plink_folder,"LD"),
                "--blocks no-pheno-req no-small-max-span --blocks-max-kb 100000",
                # "--blocks-strong-lowci .52 --blocks-strong-highci 1",
                "--blocks-inform-frac",LD_block_size," --blocks-min-maf 0 --out",file.path(plink_folder,"LD")) )
  # system( paste("plink", "--bfile plink --ld-snp-list snp_list.txt --r") )
  blocks <- data.table::fread("./plink_tmp/plink.blocks.det")
  return(blocks)
}


LD.leadSNP_block <- function(leadSNP, plink_folder, LD_block_size=.7){
  printer("Returning lead SNP's block...")
  blocks <- LD.LD_blocks(plink_folder, LD_block_size)
  splitLists <- strsplit(blocks$SNPS,split = "[|]")
  block_snps <- lapply(splitLists, function(l, leadSNP){if(leadSNP %in% l){return(l)} }, leadSNP=leadSNP) %>% unlist()
  printer("Number of SNPs in LD block =", length(block_snps))
  return(block_snps)
}


# LD_clumping <- function(subset_vcf, subset_SS){
#   # PLINK clumping: http://zzz.bwh.harvard.edu/plink/clump.shtml
#   # Convert vcf to .map (beagle)
#   ## https://www.cog-genomics.org/plink/1.9/data
#   system(paste("plink", "--vcf",subset_vcf,"--recode beagle --out ./plink_tmp/plink"))
#   # Clumping
#   system("plink", "--file ./plink_tmp/plink.chr-8 --clump",subset_SS,"--out ./plink_tmp")
# }





