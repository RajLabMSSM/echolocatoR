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


#' @section Linkage disequilibrium
#'
#' Calculate and/or query LD from reference panels, or use-supplied datasets.
#'
#' @param subset_DT The locus subset of the full summary stats file.
#' @inheritParams finemap_pipeline
#' @return A symmetric LD matrix of pairwise \emph{r} values.
LD.load_or_create <- function(results_path,
                              subset_DT,
                              locus,
                              force_new_LD=F,
                              LD_reference="1KG_Phase1",
                              superpopulation="EUR",
                              download_reference=T,
                              min_r2=0,
                              LD_block=F,
                              block_size=.7,
                              min_Dprime=F,
                              remove_correlates=F,
                              fillNA=0,
                              verbose=T,
                              server=F,
                              remove_tmps=T){
  if(LD_reference=="UKB"){
    printer("LD:: Using UK Biobank LD reference panel...")
    LD_matrix <- LD.UKBiobank(subset_DT = subset_DT,
                              results_path = results_path,
                               force_new_LD = force_new_LD,
                               chimera = server,
                               download_method = "direct",
                               nThreads = 4,
                               return_matrix = T,
                               remove_tmps = remove_tmps)
  } else if (LD_reference == "1KG_Phase1" | LD_reference == "1KG_Phase3") {
    LD_path <- file.path(results_path,"plink",paste0(LD_reference,"_LD.RDS"))
    if(!file.exists(LD_path) | force_new_LD==T){
      printer("+ Computing LD matrix...", verbose)
      LD_matrix <- LD.1KG(results_path = results_path,
                           subset_DT = subset_DT,
                           gene = locus,
                           reference = LD_reference,
                           superpopulation = superpopulation,
                           download_reference = download_reference,

                           min_r2 = min_r2,
                           LD_block = LD_block,
                           block_size = block_size,
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

plink_file <- function(base_url="./echolocatoR/tools/plink"){
  os <- get_os()
  if (os=="osx") {
    plink_version <- file.path(base_url, "plink1.9_mac");
  } else if  (os=="linux") {
    plink_version <- file.path(base_url, "plink1.9_linux");
  } else {
    plink_version <- file.path(base_url, "plink1.9_windows.exe");
  }
  return(plink_version)
}
# plink_file()


translate_population <- function(superpopulation){
  pop_dict <- list("AFA"="AFR", "CAU"="EUR", "HIS"="AMR",
                   "AFR"="AFR","EUR"="EUR", "AMR"="AMR")
  translated.list <- as.character(pop_dict[superpopulation])
  return(translated.list)
}

list_all_vcfs <-function(){
  all_vcfs <- paste0("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr",c(1:22),
         ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", collapse="\n")
  data.table::fwrite(list(all_vcfs), "Data/Reference/1000_Genomes/1KG-P3_vcfs.txt", quote = F)
}

LD_plot <- function(LD_matrix, subset_DT, span=10){
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


LD.1KG_download_vcf <- function(subset_DT,
                                 reference="1KG_Phase1",
                                 vcf_folder="./Data/Reference/1000_Genomes",
                                 results_path,
                                 gene,
                                 download_reference=T,
                                 whole_vcf=F){
  ## http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
  # Download portion of vcf from 1KG website

  # Don't use the chr prefix: https://www.internationalgenome.org/faq/how-do-i-get-sub-section-vcf-file/
  subset_DT$CHR <- gsub("chr","",subset_DT$CHR)
  chrom <- unique(subset_DT$CHR)

  # PHASE 3 DATA
  if(reference=="1KG_Phase3"){
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
  } else if (reference=="1KG_Phase1") {
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
  phase <- gsub("1KG_","",reference)
  popDat <-  data.table::fread(text=gsub(",\t",",",readLines(popDat_URL)),
                               header = F, sep="\t",  fill=T, stringsAsFactors = F,
                               col.names = c("sample","population","superpop","platform"))
  # library(Rsamtools); #BiocManager::install("Rsamtools")
  subset_vcf <- file.path(vcf_folder, phase,
                          ifelse(whole_vcf, paste(reference, paste0("chr",chrom),"vcf", sep="."),
                                                 paste(gene,"subset.vcf",sep="_")))
  # Create directory if it doesn't exist
  dir.create(path = dirname(subset_vcf), recursive = T, showWarnings = F)
  # Download and subset vcf if the subset doesn't exist already
  if(!file.exists(subset_vcf)){
    if(whole_vcf){
      region <= ""
      gene=""
    } else {
      # Specify exact positions
      regions.bed <- file.path(results_path,"plink","regions.tsv")
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


RSAMTOOLS.download_filter_vcf <- function(){
  gr <- GenomicRanges::makeGRangesFromDataFrame(df = subset_DT[,c("CHR","POS")],
                                                seqnames.field = "CHR",
                                                start.field = "POS",
                                                end.field = "POS")
  # Rsamtools::scanTabix(file = vcf_URL, param = gr)
  library(Rsamtools)
  tbx <- Rsamtools::TabixFile(file = vcf_URL,  index = paste0(vcf_URL,".tbi"))
  countTabix(tbx, param=gr)
}


BCFTOOLS.filter_vcf <- function(subset_vcf,
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

PLINK.vcf_to_bed <- function(vcf.gz.subset,
                             results_path){
  printer("LD:PLINK:: Converting vcf.gz to .bed/.bim/.fam")
  dir.create(file.path(results_path,"plink"), recursive = T, showWarnings = F)
  cmd <- paste(plink_file(),"--vcf",vcf.gz.subset, "--out", file.path(results_path,"plink","plink"))
  system(cmd)
}

PLINK.calculate_LD <- function(results_path,
                               ld_window=1000, # 10000000
                               ld_format="bin"
                               ){
  printer("LD:PLINK:: Calculating LD ( r & D'-signed; LD-window =",ld_window,")")
  plink_path_prefix <- file.path(results_path,"plink","plink")
  dir.create(file.path(results_path,"plink"), recursive = T, showWarnings = F)
  out_prefix <- paste0(plink_path_prefix,".r_dprimeSigned")
  if(ld_format=="bin"){
    cmd <- paste(plink_file(),
                 "--bfile",plink_path_prefix,
                 "--r square bin",
                 "--out",out_prefix)
    ld.path <- paste0(out_prefix,".ld.bin")
  } else {
    cmd <- paste(plink_file(),
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

PLINK.read_bin <- function(ld.path){
  bim <- data.table::fread(file.path(dirname(ld.path), "plink.bim"), col.names = c("CHR","SNP","V3","POS","A1","A2"))
  bin.vector <- readBin(ld.path, what = "numeric", n=length(bim$SNP)^2)
  ld.matrix <- matrix(bin.vector, nrow = length(bim$SNP), dimnames = list(bim$SNP, bim$SNP))
}


PLINK.read_ld_table <- function(ld.path, snp.subset=F){
  # subset_DT <- data.table::fread(file.path(results_path,"Multi-finemap/Multi-finemap_results.txt")); snp.subset <- subset_DT$SNP
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







GASTON.filter_vcf <- function(subset_vcf,
                               subset_DT,
                               results_path,
                               superpopulation,
                               popDat){
  # Import w/ gaston and further subset
  printer("+ Importing VCF as bed file...")
  bed.file <- gaston::read.vcf(subset_vcf, verbose = F)
  ## Subset rsIDs
  bed <- gaston::select.snps(bed.file, id %in% subset_DT$SNP & id !=".")
  # Create plink sub-dir
  dir.create(file.path(results_path, "./plink"), recursive = T, showWarnings = F)
  gaston::write.bed.matrix(bed, file.path(results_path, "./plink/plink"), rds = NULL)
  # Subset Individuals
  selectedInds <- subset(popDat, superpop == superpopulation)
  bed <- gaston::select.inds(bed, id %in% selectedInds$sample)
  # Cleanup extra files
  remove(bed.file)
  # file.remove("subset_vcf")
  return(bed)
}

LD.1KG <- function(results_path,
                    subset_DT,
                    gene,
                    reference="1KG_Phase1",
                    superpopulation="EUR",
                    vcf_folder="./Data/Reference/1000_Genomes",
                    download_reference=T,
                    min_r2=F,
                    LD_block=F,
                    block_size=.7,
                    min_Dprime=F,
                    remove_correlates=F,
                    remove_tmps=T,
                    fillNA=0){
          # Quickstart
  # gene <- "LRRK2"; results_path <- file.path("./Data/GWAS/Nalls23andMe_2019",gene); reference="1KG_Phase1"; vcf_folder="./Data/Reference/1000_Genomes"; superpopulation="EUR"; vcf_folder="./Data/Reference/1000_Genomes"; min_r2=F; LD_block=F; block_size=.7; min_Dprime=F;  remove_correlates=F; download_reference=T; subset_DT <- data.table::fread(file.path(results_path,"Multi-finemap/Multi-finemap_results.txt"))
    printer("LD:: Using 1000Genomes LD reference panel...")
    vcf_info <- LD.1KG_download_vcf(subset_DT=subset_DT,
                                     results_path,
                                     reference=reference,
                                     vcf_folder=vcf_folder,
                                     gene=gene,
                                     download_reference=download_reference)
    subset_vcf <- vcf_info$subset_vcf
    popDat <- vcf_info$popDat
    vcf.gz.path <- BCFTOOLS.filter_vcf(subset_vcf = subset_vcf,
                                       popDat = popDat,
                                       superpopulation = superpopulation,
                                       remove_tmp = F)
    PLINK.vcf_to_bed(vcf.gz.subset = vcf.gz.path,
                     results_path = results_path)
    # Calculate pairwise LD for all SNP combinations
    #### "Caution that the LD matrix has to be correlation matrix" -SuSiER documentation
    ### https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html
    # Gaston LD method
    # LD_matrix <- gaston::LD(bed, lim = c(1,ncol(bed)), measure ="r") #"D"
    # LD_matrix[!is.finite(LD_matrix)] <- 0

    # Get lead SNP rsid
    leadSNP = subset(subset_DT, leadSNP==T)$SNP #rs76904798
    # Plink LD method
    LD_matrix <- plink_LD(plink_folder = file.path(results_path,"plink"),
                          leadSNP = leadSNP,
                          min_r2 = min_r2,
                          min_Dprime = min_Dprime,
                          remove_correlates = remove_correlates,
                          fillNA = fillNA)
    # Filter out SNPs not in the same LD block as the lead SNP
    if(LD_block){
      block_snps <- leadSNP_block(leadSNP, "./plink_tmp", block_size)
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



############# ############# #############
#############    PLINK      #############
############# ############# #############

Dprime_table <- function(SNP_list, plink_folder){
  printer("+ Creating DPrime table")
  system( paste(plink_file(), "--bfile",file.path(plink_folder,"plink"),
                "--ld-snps", paste(SNP_list, collapse=" "),
                "--r dprime-signed",
                "--ld-window 10000000",
                "--ld-window-kb 10000000",
                "--out",file.path(plink_folder,"plink")) )
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

# complex_LD <- function(bim, plink_folder, min_Dprime=F){
#   # METHOD 1
#   plink.ld <- Dprime_table(bim)
#   # Filter NaNs
#   plink.ld <- plink.ld[!is.nan(plink.ld$DP),]
#   plink.ld <- plink.ld[!is.nan(plink.ld$R),]
#
#   if(min_Dprime!=F){
#     printer("\n++++++++++ Filtering by DPrime ++++++++++\n")
#     plink.ld <- subset(plink.ld, DP>=min_Dprime)
#   }
#   ld.matrix <- data.table::dcast.data.table(plink.ld, formula = SNP_B ~ SNP_A, value.var="R",
#                                             fill=0, drop=F, fun.aggregate = mean)
#   ld.matrix <-  data.frame(ld.matrix, row.names = ld.matrix$SNP_B) %>% subset(select = -SNP_B) %>%
#     data.table() %>% as.matrix()
#   return(ld.matrix)
# }

run_plink_LD <-function(bim, plink_folder, r_format="r"){
  # METHOD 2 (faster, but less control over parameters. Most importantly, can't get Dprime)
  system( paste(plink_file(), "--bfile",file.path(plink_folder,"plink"),
                "--extract",file.path(plink_folder,"SNPs.txt"),
                paste0("--",r_format," square bin"), "--out", file.path(plink_folder,"plink")) )
  bin.vector <- readBin(file.path(plink_folder, "plink.ld.bin"), what = "numeric", n=length(bim$SNP)^2)
  ld.matrix <- matrix(bin.vector, nrow = length(bim$SNP), dimnames = list(bim$SNP, bim$SNP))
  return(ld.matrix)
}


plink_LD <-function(leadSNP,
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
    plink.ld <- Dprime_table(SNP_list = row.names(ld.matrix), plink_folder)

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

LD_blocks <- function(plink_folder, block_size=.7){
  printer("++ Calculating LD blocks...")
  # PLINK 1.07 LD: http://zzz.bwh.harvard.edu/plink/ld.shtml
  # PLINK 1.9 LD: https://www.cog-genomics.org/plink/1.9/ld
  # system(plink_file(), "-h")
  # Identify duplicate snps
  # system(plink_file(), "--vcf subset.vcf --list-duplicate-vars")
  # Convert vcf to plink format
  # system(plink_file(), "--vcf subset.vcf --exclude ./plink_tmp/plink.dupvar --make-bed --out PTK2B")

  # Estimate LD blocks
  # Defaults: --blocks-strong-lowci = 0.70, --blocks-strong-highci .90

  # Reucing "--blocks-inform-frac" is the only parameter that seems to make the block sizes larger
  system( paste(plink_file(), "--bfile",file.path(plink_folder,"plink"),
                "--blocks no-pheno-req no-small-max-span --blocks-max-kb 100000",
                # "--blocks-strong-lowci .52 --blocks-strong-highci 1",
                "--blocks-inform-frac",block_size," --blocks-min-maf 0 --out",file.path(plink_folder,"plink")) )
  # system( paste(plink_file(), "--bfile plink --ld-snp-list snp_list.txt --r") )
  blocks <- data.table::fread("./plink_tmp/plink.blocks.det")
  return(blocks)
}


leadSNP_block <- function(leadSNP, plink_folder, block_size=.7){
  printer("Returning lead SNP's block...")
  blocks <- LD_blocks(plink_folder, block_size)
  splitLists <- strsplit(blocks$SNPS,split = "[|]")
  block_snps <- lapply(splitLists, function(l, leadSNP){if(leadSNP %in% l){return(l)} }, leadSNP=leadSNP) %>% unlist()
  printer("Number of SNPs in LD block =", length(block_snps))
  return(block_snps)
}


# LD_clumping <- function(subset_vcf, subset_SS){
#   # PLINK clumping: http://zzz.bwh.harvard.edu/plink/clump.shtml
#   # Convert vcf to .map (beagle)
#   ## https://www.cog-genomics.org/plink/1.9/data
#   system(paste(plink_file(), "--vcf",subset_vcf,"--recode beagle --out ./plink_tmp/plink"))
#   # Clumping
#   system(plink_file(), "--file ./plink_tmp/plink.chr-8 --clump",subset_SS,"--out ./plink_tmp")
# }





snpStats_LDblocks <- function(){
  # https://cran.r-project.org/web/packages/adjclust/vignettes/snpClust.html
  library(coloc)
  library(snpStats)
  library(matrixStats)
  library(adjclust)
  data("ld.example", package = "snpStats")

  # snpStats::read.plink(bed = "Data/eQTL/Fairfax_CD14/plink/eQTL_283_hg19_040111.bed",
  #                      bim = "Data/eQTL/Fairfax_CD14/plink/eQTL_283_hg19_040111.bed",
  #                      fam = "Data/eQTL/Fairfax_CD14/plink/eQTL_283_hg19_040111.fam")
  geno <- ceph.1mb[, -316]  ## drop one SNP leading to one missing LD value
  p <- ncol(geno)
  nSamples <- nrow(geno)
  geno
  ld.ceph <- snpStats::ld(geno, stats = "R.squared", depth = p-1)
  image(ld.ceph, lwd = 0)

  fitH <- adjclust::snpClust(geno, h = 100, stats = "R.squared")
  plot(fitH, type = "rectangle", leaflab = "perpendicular")

  head(cbind(fitH$merge, fitH$gains))
  fitH$merge

  str(fitH)
}





GSA.download_vcf <- function(){
  lead.pos <- 40614434
  one.Mb <- 500000*2
  cmd <- paste0("tabix chr12.dose.vcf.gz chr12:", lead.pos-one.Mb,"-",lead.pos+one.Mb)
}

# Turn LD matrix into positive semi-definite matrix
# LD_matrix_check <- function(LD_matrix){
#   mat <- LD_matrix %>% data.table() %>% as.matrix()
#   mat <- Matrix::forceSymmetric(mat)
#   if(matrixcalc::is.positive.semi.definite(mat)){
#     return(LD_matrix)
#   } else{
#     # Makes matrix symmetric and positive definite
#     ## (but super slow...)
#     mat <- Matrix::nearPD(mat)$mat %>% as.matrix()
#   }
# }



