


################## QTL Data ################## 
# V- psychENCODE
# V- Fairfax
# V- MESA
# V- Cardiogenics
# V- Brain_xQTL_Serve
# - ROSMAP + CommmonMind (larget than Brain_xQTL_Serve)
# - STARNET
# V- GTEx V7 & V8 (49 tissues)


mergeQTL.add_SNP_id <- function(FM_all){
  FM_all <- FM_all %>% 
    dplyr::mutate(SNP_id=paste0(gsub("chr","",CHR),":",POS)) %>% 
    data.table::data.table() 
  return(FM_all)
}


mergeQTL.psychENCODE_SNPinfo <- function(QTL.sub, 
                                         SNPinfo_path="/sc/orga/projects/ad-omics/data/psychENCODE/SNP_Information_Table_with_Alleles.txt.gz"){
  printer("+ mergeQTL: Merging external SNP info...")
  SNPinfo <- data.table::fread(SNPinfo_path, nThread = 4)
  SNPinfo.sub <- subset(SNPinfo, PEC_id %in% QTL.sub$SNP_id, select=c("PEC_id","Rsid","REF","ALT"))
  snp_merge <- data.table:::merge.data.table(QTL.sub, SNPinfo.sub, 
                                             by.x = "SNP_id", by.y = "PEC_id", all.x = T)
  return(snp_merge)
}

# psychENCODE
## Need: StdErr, MAF
mergeQTL.psychENCODE <- function(FM_all=merge_finemapping_results(minimum_support = 0), 
                                 local_files=T,
                                 force_new_subset=F,
                                 ASSAYS=c("eQTL","cQTL","isoQTL","tQTL") # 
                                 ){ 
  FM_all <- mergeQTL.add_SNP_id(FM_all)
  for(assay in ASSAYS){
    printer("psychENCODE:: Checking for overlap with",assay) 
    output_path <- file.path("./Data/QTL/psychENCODE", assay, paste0("psychENCODE_",assay,".finemap.txt"))
    output_path_gz <- paste0(output_path,".gz")
    if(file.exists(output_path_gz) & force_new_subset==F){
      printer("psychENCODE:: Importing pre-existing file.")
      QTL.sub <- data.table::fread(output_path_gz, nThread = 4)
      # Drop duplicate columns (Peak_center)
      QTL.sub <- QTL.sub[,unique(names(QTL.sub)),with=F] 
    } else {
      server_file <- Directory_info(dataset_name = paste0("psychENCODE_",assay), "fullSS.local")
      QTL <- data.table::fread(server_file, nThread = 4)
      fullSS_nrow <- get_nrows(server_file)
      # Subset QTL data
      if(assay %in% c("fQTL")){
        QTL <- QTL %>% dplyr::rename(gene_chr=Chromosome_of_variant, 
                                      gene_start=Locus_of_variant,
                                      regression_slope=Regression_slope) %>% 
        dplyr::mutate(FDR=p.adjust(Nominal_p_val_of_association,method = "fdr", n = n_snps),
                      SNP_id=paste0(gsub("chr","",gene_chr),":",gene_start))
        }
      QTL.sub <- QTL %>% subset(SNP_id %in% FM_all$SNP_id)
      QTL.sub <- mergeQTL.psychENCODE_SNPinfo(QTL.sub) 
                                               
      # Write
      dir.create(dirname(output_path), showWarnings = F, recursive = T)
      data.table::fwrite(QTL.sub, output_path, sep="\t", nThread = 4)
      R.utils::gzip(output_path, overwrite=T, remove=T)
    }
    
    if(!("gene_id" %in% colnames(QTL.sub))){QTL.sub$gene_id <- NA}
    # QTLs have can multiple probes per SNP location. 
    ## Pick only the best one (lowest FDR and highest Effect) keep each row as a unique genomic position:
    QTL.sub <- QTL.sub %>%
      dplyr::group_by(SNP_id) %>%
      arrange(FDR, desc(regression_slope)) %>% 
      dplyr::slice(1) %>% 
      data.table::data.table()
    # Select and rename columns
    QTL.sub <-  (QTL.sub %>%  
      dplyr::mutate(SNP=Rsid,
                    QTL.Effect=regression_slope,
                    QTL.StdErr=NA,
                    QTL.P=nominal_pval,
                    QTL.FDR=FDR,
                    QTL.MAF=NA, 
                    QTL.A1=REF,
                    QTL.A2=ALT,
                    QTL.SampleSize=1866, #Assuming same size for each (not sure if true)
                    QTL.Gene=gene_id) %>%
      dplyr::select(SNP, QTL.Effect, QTL.StdErr, QTL.P, QTL.FDR, QTL.MAF, QTL.A1, QTL.A2, QTL.SampleSize, QTL.Gene) %>%
      arrange(QTL.FDR, desc(QTL.Effect)) ) %>% 
      data.table::data.table()  
    FM_all <- data.table:::merge.data.table(FM_all,
                                            QTL.sub,
                                            by = "SNP",
                                            all.x = T, allow.cartesian = T)
  }
  return(FM_all)
}

  

mergeQTL.Fairfax_probemap <- function(gene_list, probe_path.="./Data/QTL/Fairfax_2014/gene.ILMN.map"){
  printer("++ Extracting probe info")
  cmd <- paste0("grep -E '", paste0(gene_list, collapse="|"),"' ",probe_path.)
  col_names <- data.table::fread(probe_path., sep="\t", nrows = 0) %>% colnames()
  probe_map <- data.table::fread(text = system(cmd, intern = T), 
                                 sep = "\t", header = F, col.names = col_names)
  if(dim(probe_map)[1]==0){
    stop("Could not identify gene in the probe mapping file: ",paste(gene_list, collapse=", "))
  }
  # Subset just to make sure you didn't accidentally grep other genes
  probe_map <- subset(probe_map, GENE %in% gene_list)
  probe_map <- setNames(probe_map$GENE, probe_map$PROBE_ID) 
  return(probe_map)
}

mergeQTL.Fairfax_extract_SNPinfo <- function(genotype_path="/sc/orga/projects/ad-omics/data/fairfax/volunteer_421.impute2.dosage",
                                             rsid.key_path="/sc/orga/projects/ad-omics/data/fairfax/rsid.key.tsv.gz"){
  ### Extract SNP, CHR, and POS from first several cols of Fairfax genotype file
  printer("mergeQTL::Fairfax: Extracting RSID, CHR, POS, REF and ALT from:",basename(genotype_path))
  if(file.exists(rsid.key_path)){
    printer("+ mergeQTL::Fairfax: Pre-existing file detected. Importing...") 
  } else {
    system(paste0("awk -F ' ' 'BEGIN { print \"RSID\tCHR\tPOS\tREF\tALT\" }; { print $2\"\t\"$1\"\t\"$3\"\t\"$4\"\t\"$5 }' ", genotype_path," | gzip > ",rsid.key_path))
    printer("+ mergeQTL::Fairfax: File saved to ==>", rsid.key_path)
  } 
  rsid.key <- data.table::fread(rsid.key_path, nThread = 4)
  return(rsid.key)
}

mergeQTL.Fairfax_extract_MAF <- function(dat, 
                                         MAF_dir=file.path("/sc/orga/projects/ad-omics/data/fairfax/",
                                                           "imputation/imputation-hrc","MAF")){
  #NOTE!: There still a lot of mismatched alleles. Need to figure out how to fix.
  dat <- dat %>% tidyr::separate(SNP_id, c('CHR', 'POS'), sep=":", extra="drop", remove=F)
  printer("mergeQTL::Fairfax: Merging with imputation files to get REF(0), ALT(1) and MAF...")
  merged.DAT <- lapply(unique(dat$CHR), function(chr){
    printer("mergeQTL::Fairfax: Chrom",chr)
    dat.chr = subset(dat, CHR==chr)
    MAF_path <- file.path(MAF_dir,paste0("chr",chr,".info"))
    # MAF_path <- file.path(MAF_dir,paste0("chr",chr,".filtered.tab"))
    maf.info <- data.table::fread(MAF_path, nThread = 4)
    maf.info <- maf.info %>% dplyr::rename(SNP_id=SNP)
    merge.dat <- data.table:::merge.data.table(dat.chr, 
                                  maf.info[,c("SNP_id","REF(0)","ALT(1)","MAF")], 
                                  by="SNP_id", all.x = T)
    # Check if these are the same ref/alt in the other Fairfax files
    ref_bool <- merge.dat$REF!=merge.dat$`REF(0)`
    alt_bool <- merge.dat$ALT!=merge.dat$`ALT(1)`
    printer(" +", sum(ref_bool, na.rm = T), "REF allleles did not match.")
    printer(" +", sum(alt_bool, na.rm = T), "ALT allleles did not match.")
    return(merge.dat)
  }) %>% data.table::rbindlist(fill=T)
  return(merged.DAT)
}


# Fairfax: eQTL
## All fields complete for coloc!
mergeQTL.Fairfax <- function(FM_all, CONDITIONS=c("CD14","IFN","LPS2","LPS24"), force_new_subset=F){
  FM_all <- mergeQTL.add_SNP_id(FM_all)
  for(condition in CONDITIONS){
    dataset <- paste0("Fairfax_2014_",condition)
    printer("mergeQTL::Fairfax:: Processing:",dataset)
    server_path <- Directory_info(paste0("Fairfax_2014_", condition),variable = "fullSS.local") 
    output_path <- file.path("./Data/QTL/Fairfax_2014",condition,paste0(dataset,".finemap.txt"))
    output_path_gz <- paste0(output_path,".gz")
    dir.create(dirname(output_path), showWarnings = F, recursive = T)
    
    if(file.exists(output_path_gz)& force_new_subset==F){
      printer("+ mergeQTL::Fairfax: Pre-existing file detected. Importing",output_path_gz,"...")
      dat <- data.table::fread(output_path_gz, nThread = 4)
    } else { 
      probe_map <- mergeQTL.Fairfax_probemap(gene_list = unique(FM_all$Gene), probe_path.="./Data/QTL/Fairfax_2014/gene.ILMN.map")
      header <- get_header(server_path)
      printer("+ mergeQTL::Fairfax: Subsetting full summary stats file...")
      system(paste0("grep -E '", paste(c(header, unique(names(probe_map))), collapse="|"),"' " ,server_path," > ",output_path))
      
      # Import data
      dat <- data.table::fread(output_path, nThread = 4)
      
      ## Rename vars
      dat <- dplyr::rename(dat, SNP_id=SNP, probe=gene)
      dat$QTL.Gene <- probe_map[dat$probe] 
      ## Add SNP col
      rsid.key <- mergeQTL.Fairfax_extract_SNPinfo()
      rsid.key <- rsid.key %>% dplyr::mutate(SNP_id=paste0(CHR,":",POS))
      ### Extract SNPs from FM_all
      dat <- data.table:::merge.data.table(dat, rsid.key[,c("RSID","REF","ALT","SNP_id")], 
                                                  by="SNP_id", all.x = T)
      ### Extract MAF from imputation files
      merged.dat <- mergeQTL.Fairfax_extract_MAF(dat)
      dat <- merged.dat
      ## Overwrite merged subset
      data.table::fwrite(dat, output_path, sep = "\t", nThread = 4) 
      R.utils::gzip(output_path, overwrite=T)
    }  
    dat$StdErr <- dat$beta / dat$`t-stat` 
    dat.sub <- subset(dat, SNP_id %in% FM_all$SNP_id) 
    # Take only the top QTL per SNP location
    ## Rename columns
    dat.sub <- (dplyr::mutate(dat.sub, 
                             SNP_id,
                             QTL.Effect=beta, 
                             QTL.StdErr=StdErr,
                             QTL.P=`p-value`,
                             QTL.FDR=FDR, 
                             QTL.MAF=MAF, 
                             QTL.A1=`REF(0)`,
                             QTL.A2=`ALT(1)`,
                             QTL.SampleSize=432,
                             QTL.Gene=QTL.Gene) %>% 
      dplyr::select(SNP_id, QTL.Effect, QTL.StdErr, QTL.P, QTL.FDR, QTL.MAF, QTL.A1, QTL.A2, QTL.SampleSize, QTL.Gene) %>%
      dplyr::group_by(SNP_id) %>%
      arrange(QTL.FDR, desc(QTL.Effect)) )%>% 
      dplyr::slice(1) %>% 
      data.table::data.table()
    FM_all <- data.table:::merge.data.table(FM_all, dat.sub, 
                                              by="SNP_id", 
                                              all.x = T, allow.cartesian = T)
    
  } 
  return(FM_all) 
}


# MESA
## Need: MAF
mergeQTL.MESA <- function(FM_all, force_new_subset=F, POPULATIONS=c("AFA","CAU","HIS")){ 
  FM_all <- mergeQTL.add_SNP_id(FM_all)
  sample_size.dict <- list(AFA=233, CAU=578, HIS=352)
  for(pop in POPULATIONS){
    printer("MESA:: Extracting",pop,"data.")
    dataset <- paste0("MESA_",pop)
    output_path <- file.path("./Data/QTL/MESA",pop,paste0(dataset,".finemap.txt.gz"))
    
    if(file.exists(output_path) & force_new_subset==F){
      printer("MESA:: Pre-existing file detect. Importing...")
      # for some reason, saved files lose their header. Have to explicitly provide them instead...
      dat <- data.table::fread(output_path, nThread = 4, 
                               col.names = c("snps","gene","statistic","pvalue","FDR","beta",
                                             "chr","gene_name","start","end","gene_type","pos_snps","ref","alt"))
    } else {
      server_path <- Directory_info(paste0("MESA_",pop), "fullSS.local")
      output_path_txt <- gsub(".gz","",output_path)
      system(paste0("grep -E '", paste(unique(FM_all$Gene), collapse="|"),"' " ,server_path," > ",output_path_txt))
      # Import data 
      dat <- data.table::fread(output_path_txt, 
                               col.names = colnames(data.table::fread(server_path, nrow=0)), 
                               nThread = 4)#file_path
      dat <- subset(dat, snps %in% FM_all$SNP)
      ## Overwrite subset w/ fixed header (and smaller size)
      dir.create(dirname(output_path), recursive = T, showWarnings = F)
      data.table::fwrite(dat, output_path, sep = "\t", col.names = T)
      R.utils::gzip(output_path_txt, overwrite=T, remove=T)
    }  
    
    # Calculare StdErr
    dat$QTL.StdErr <- dat$beta / dat$statistic
    # Take only the top QTL per SNP location
    ## Rename columns
    dat.sub <- (dplyr::mutate(dat, 
                             SNP=snps, 
                             QTL.Effect=beta, 
                             QTL.StdErr=QTL.StdErr,
                             QTL.P=pvalue,
                             QTL.FDR=FDR, 
                             QTL.MAF=NA,
                             QTL.A1=ref,
                             QTL.A2=alt,
                             QTL.SampleSize=sample_size.dict[[pop]],
                             QTL.Gene=gene_name) %>% 
      dplyr::select(SNP, QTL.Effect, QTL.StdErr, QTL.P, QTL.FDR, QTL.MAF, QTL.A1, QTL.A2, QTL.SampleSize, QTL.Gene) %>%
      dplyr::group_by(SNP) %>%
      arrange(QTL.FDR, desc(QTL.Effect)) ) %>% 
      dplyr::slice(1) %>% 
      data.table::data.table()
    FM_all <- data.table:::merge.data.table(FM_all, dat.sub, 
                                            by="SNP", 
                                            all.x = T)
  } 
  return(FM_all)
}


# Cardiogenics
## Need: MAF, A1, A2
mergeQTL.Cardiogenics <- function(FM_all, 
                                  force_new_subset=F, 
                                  cis_only=T, 
                                  CELLTYPES=c("macrophages","monocytes")){
  FM_all <- mergeQTL.add_SNP_id(FM_all)
  for(celltype in CELLTYPES){
    printer("Cardiogenics:: Processing",celltype,"data...")
    dataset <- paste0("Cardiogenics_",celltype)
    output_path <- file.path("./Data/QTL/Cardiogenics",celltype,paste0(dataset,".finemap.txt"))
    output_path_gz <- paste0(output_path, ".gz")
    if(file.exists(output_path_gz) & force_new_subset==F){
      print("Cardiogenics:: Pre-existing file detected. Importing...")
      dat.sub <- data.table::fread(output_path_gz, nThread = 4)
    } else{ 
      server_path <- Directory_info(paste0("Cardiogenics_",celltype), "fullSS.local")
      dir.create(dirname(output_path), recursive = T, showWarnings = F)
      DAT <- data.table::fread(server_path,  nThread = 4)
      dat.sub <- subset(DAT, SNPID %in% unique(FM_all$SNP))
      data.table::fwrite(dat.sub, output_path, sep="\t", nThread = 4)
      R.utils::gzip(output_path)
    }
    
    dat.sub$StdErr <- dat.sub$beta / dat.sub$t.test
    
    if(cis_only){dat.sub <- subset(dat.sub, relativePosition=="cis")}
    dat.sub <- (dplyr::mutate(dat.sub, 
                             SNP=SNPID,
                             QTL.Effect=beta,
                             QTL.StdErr=StdErr,
                             QTL.P=NA,
                             QTL.FDR=FDR, 
                             QTL.MAF=NA, 
                             QTL.A1=NA,
                             QTL.A2=NA,
                             QTL.SampleSize=n.samples, # 758 individuals
                             QTL.Gene=reporterID) %>% 
      dplyr::select(SNP, QTL.Effect, QTL.StdErr, QTL.P, QTL.FDR, QTL.MAF, QTL.A1, QTL.A2, QTL.SampleSize, QTL.Gene) %>%
      dplyr::group_by(SNP) %>%
      arrange(QTL.FDR, desc(QTL.Effect)) ) %>% 
      dplyr::slice(1) %>%
      data.table::data.table()
    FM_all <- data.table:::merge.data.table(FM_all, dat.sub, 
                                            by="SNP", 
                                            all.x = T)
  } 
  return(FM_all)
}


mergeQTL.GTEx_list_files <- function(output_dir, 
                                     server_path, 
                                     GTEx_version="GTEx_V7", 
                                     fuzzy_search=F, 
                                     local_files=T){ 
  printer("GTEx:: Constructing reference file of available single-tissue eQTL files...")
  if(local_files){
    SS_files <- list.files(output_dir, pattern="*.finemap.txt.gz")
    
  } else {
    SS_files <- list.files(server_path, 
                           pattern=ifelse(GTEx_version=="GTEx_V7",".allpairs.txt.gz","*.egenes.*"), 
                           full.names = T) 
  }
 
  all_tissues <- lapply(SS_files, function(e){strsplit(basename(e),"[.]")[[1]][1] }) %>% unlist() %>% unique()
  all_tissues <- gsub(paste0(GTEx_version,"_"),"", all_tissues)
  dir.create(output_dir, showWarnings = F, recursive = T)
  tissues_df <- data.frame(tissue=all_tissues, sum_stats=SS_files)
  if(fuzzy_search!=F){
    tissues_df = tissues_df[base::grep(fuzzy_search, tissues_df$tissue),]
  }
  return(tissues_df)
}

mergeQTL.HGNC_to_ENS <- function(geneList){ 
  print("mergeQTL:: Gathering Ensembl_gene_id mappings for HGNC symbols...")
  info <- biomart_geneInfo(geneList = geneList)
  geneDict <- setNames(info$ensembl_gene_id, info$hgnc_symbol)
  converted_genes <- geneDict[geneList]
  na_genes <-  info$hgnc_symbol[as.logical(is.na(converted_genes))]
  printer("Could not find mappings for",length(na_genes),"genes:",paste(na_genes, collapse=", "))
  return(converted_genes) 
}

mergeQTL.GTEx_snp_dict <- function(snp_list){
  ref <- data.table::fread(file.path(server_path,"GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz"), nThread = 4)
  snp_dict <- setNames(ref$rs_id_dbSNP147_GRCh37p13, ref$variant_id)
  return(snp_dict)
}
 

# GTEx
## All fields complete for coloc!
mergeQTL.GTEx <- function(FM_all, fuzzy_search="Brain", GTEx_version="GTEx_V7", force_new_subset=F, local_files=T){
  FM_all <- mergeQTL.add_SNP_id(FM_all)
  server_path <- Directory_info(GTEx_version, "fullSS.local")
  output_dir <- file.path("./Data/QTL",GTEx_version)
  tissues_df <- mergeQTL.GTEx_list_files(output_dir, server_path, GTEx_version, fuzzy_search, local_files = local_files)
 
  # Gather QTL data
  for(tiss in tissues_df$tissue){
    printer(GTEx_version,":: Processsing eQTL for",tiss)
    dataset <- paste0(GTEx_version,"_",tiss) 
    output_path <- file.path(output_dir, paste0(dataset,".finemap.txt"))
    output_path_gz <- paste0(output_path,".gz")
    if(file.exists(output_path_gz) & force_new_subset==F){
      dat <- data.table::fread(output_path_gz, nThread = 4)
    } else { 
      geneDict <- mergeQTL.HGNC_to_ENS(geneList = unique(FM_all$Gene))
      ensembl_genes <- as.character(geneDict)[!is.na(geneDict)]
      server_file <- subset(tissues_df, tissue==tiss)$sum_stats %>% as.character() 
      header <- get_header(server_file)
      fullSS_nrows <- get_nrows(server_file)#Differs for each SS file, so need to calculate for each.
      system(paste0("zcat ",server_file," | grep -E '",paste(c(header, ensembl_genes), collapse="|"),"' > ",output_path))
      
      dat <- data.table::fread(output_path, nThread = 4)
      # Rewrite with FDR included
      dat <-  dat %>% tidyr::separate(variant_id, c('CHR', 'POS', 'A1','A2','build'), sep="_", extra="drop") %>%
                      dplyr::mutate(FDR=p.adjust(pval_nominal, method = "fdr", n=fullSS_nrows),
                                    CHR=as.numeric(CHR), 
                                    POS=as.numeric(POS))
      data.table::fwrite(dat, output_path)
      R.utils::gzip(output_path, overwrite=T, remove=T)
    }
    # Subset
    dat.sub <- (dat %>% 
      dplyr::select(CHR, POS,  
                    QTL.Effect=slope,
                    QTL.StdErr=slope_se, 
                    QTL.P=pval_nominal,
                    QTL.FDR=FDR, 
                    QTL.MAF=maf,
                    QTL.A1=A1,
                    QTL.A2=A2,
                    QTL.SampleSize=ma_samples, # 1000?
                    QTL.Gene=gene_id) %>%
        arrange(QTL.FDR, desc(QTL.Effect))) %>%
      data.table::data.table() 
    FM_all <- data.table:::merge.data.table(FM_all, dat.sub, 
                                            by=c("CHR","POS"), 
                                            all.x = T, allow.cartesian = T)  
  }
  return(FM_all)
}


# Brain_xQTL_Serve
## Need: MAF, A1, A2, StdErr, Sample Size

mergeQTL.Brain_xQTL_Serve <- function(FM_all, 
                                      force_new_subset=F, 
                                      ASSAYS=c("eQTL","haQTL","mQTL","cell-specificity-eQTL")){
  FM_all <- mergeQTL.add_SNP_id(FM_all)
  for(assay in ASSAYS){
    printer("Brain_xQTL_Serve:: Processing",assay,"...")
    dataset <- paste0("Brain_xQTL_Serve_",assay) 
    output_path <- file.path("./Data/QTL/Brain_xQTL_Serve",assay,paste0("Brain_xQTL_Serve.",assay,".finemap.txt"))
    output_path_gz <- paste0(output_path,".gz")
    if(file.exists(output_path_gz) & force_new_subset==F){
      printer("+ Brain_xQTL_Serve:: Pre-existing file detected. Importing...")
      QTL <- data.table::fread(output_path_gz, nThread = 4)
    } else {
      printer("+ Brain_xQTL_Serve:: Importing and subsettting full",assay,"dataset...")
      dir.create(dirname(output_path), showWarnings = F, recursive = T)
      server_file <- Directory_info(paste0("Brain.xQTL.Serve_",assay), "fullSS.local")
      header <- get_header(server_file)
      col1 <- strsplit(header, "\t")[[1]][1]
      
      if(assay %in% c("haQTL","mQTL")){
        printer("+ Brain_xQTL_Serve:: Importing large",assay,"file and subsetting via data.table...")
        QTL <- data.table::fread(server_file, nThread = 4)
        fullSS_nrows <- nrow(QTL)
        QTL <- subset(QTL, SNPid %in% unique(FM_all$SNP))
        QTL <- (QTL %>% dplyr::group_by(SNPid) %>% arrange(pValue, desc(SpearmanRho))) %>% 
          dplyr::slice(1) %>% data.table()
        data.table::fwrite(QTL, output_path, nThread=4, sep="\t")
      } else {
        fullSS_nrows <- get_nrows(server_file)
        system(paste0("zcat ",server_file," | grep -E '",paste(c(col1, unique(FM_all$Gene)), collapse="|"),
                    "' > ",output_path))
      }
      
      # Calculate FDR
      QTL <- data.table::fread(output_path, nThread = 4) 
      if(assay=="cell-specificity-eQTL"){
        QTL <- QTL %>% dplyr::rename(SNP=`SNP name`,
                                     FDR=`FDR (specific to each cell type)`,
                                     cell_type=`Cell type`)
      } else { 
        QTL$FDR <- p.adjust(QTL$pValue, method="fdr", n=fullSS_nrows) 
        QTL <- QTL %>% dplyr::rename(SNP=SNPid)
        } 
      # Rewrite after adding FDR
      data.table::fwrite(QTL, output_path, nThread = 4, sep="\t")
      R.utils::gzip(output_path, overwrite=T, remove=T)
    }
    
    # Subset
    dat.sub <- subset(QTL, SNP %in% unique(FM_all$SNP))
    if(assay=="cell-specificity-eQTL"){  
      dat.sub$Effect <- NA
      dat.sub <- (dat.sub %>% dplyr::mutate(SNP, 
                                            QTL.Effect=NA,  
                                            QTL.StdErr=NA,
                                            QTL.P=`Interaction pvalue`,
                                            QTL.FDR=FDR,  
                                            QTL.MAF=NA, 
                                            QTL.A1=NA,
                                            QTL.A2=NA,
                                            QTL.SampleSize=NA,
                                            # QTL.CellType=cell_type,
                                            QTL.Gene=Gene) %>% 
                    dplyr::select(SNP, QTL.Effect, QTL.StdErr, QTL.P, QTL.FDR, QTL.MAF, QTL.A1, QTL.A2, QTL.SampleSize, QTL.Gene) %>% 
        arrange(QTL.FDR, desc(QTL.Effect))) %>% 
        data.table::data.table()
    } else {
      dat.sub <- (dat.sub %>%
        dplyr::mutate(SNP, 
                      QTL.Effect=SpearmanRho, 
                      QTL.StdErr=NA,
                      QTL.P=pValue,
                      QTL.FDR=FDR,  
                      QTL.MAF=NA, 
                      QTL.A1=NA,
                      QTL.A2=NA,
                      QTL.SampleSize=NA,
                      QTL.Gene=featureName) %>% 
          dplyr::select(SNP, QTL.Effect, QTL.StdErr, QTL.P, QTL.FDR, QTL.MAF, QTL.A1, QTL.A2, QTL.SampleSize, QTL.Gene) %>%
        arrange(QTL.FDR, desc(QTL.Effect))) %>% 
        data.table::data.table()
    }
    FM_all <- data.table:::merge.data.table(FM_all, dat.sub, 
                                            by="SNP", 
                                            all.x = T, allow.cartesian = T)
  }
  return(FM_all)
}




############### Post-processing ################

mergeQTL.melt_FDR <- function(FM_merge){ 
  dim(FM_merge) 
  Effect.cols <- base::grep(".Effect", colnames(FM_merge), value = T)
  FDR.cols <- base::grep(".FDR", colnames(FM_merge), value = T)
  id.vars <- c("Gene","SNP","CHR","POS","P","Consensus_SNP","Support","leadSNP")
  FM_sub <- subset(FM_merge, select=c(id.vars, FDR.cols))
  colnames(FM_sub)[colnames(FM_sub) %in% FDR.cols] <- gsub(".FDR","",FDR.cols)
  FM_melt <- data.table::melt.data.table(FM_sub, 
                                         id.vars = id.vars, 
                                         variable.name = "QTL.Source",
                                         value.name = "FDR")
  FM_melt[is.infinite(FM_melt$FDR),"FDR"] <- NA 
  FM_melt[FM_melt$FDR %in% c(Inf,-Inf),"FDR"] <- NA 
  FM_melt[FM_melt$FDR==0,"FDR"] <- .Machine$double.xmin 
  return(FM_melt)
}

####### PLOTS #######


mergeQTL.count_overlap <- function(){
  library(dplyr)
  # FM_merge <- FM_all
  FM_merge <- data.table::fread(file.path("./Data/GWAS/Nalls23andMe_2019/_genome_wide",
                                          "Nalls23andMe_2019.QTL_overlaps.txt.gz"),
                                nThread = 4)
  FM_melt <- mergeQTL.melt_FDR(FM_merge) 
  
  mergedQTL.get_count <- function(SNP.subset){
    count.df <- subset(SNP.subset, (FDR<0.05), drop=F) %>% 
      dplyr::group_by(QTL.Source, .drop=F) %>% 
      tally(sort = F) %>% 
      arrange(QTL.Source)
    count.df$all.SNPs <- nrow(SNP.subset)
    count.df$Proportion <- count.df$n / count.df$all.SNPs
    return(count.df)
  }
  
  consensus <- mergedQTL.get_count(subset(FM_melt, (Consensus_SNP==T)))
  cred.set <- mergedQTL.get_count(subset(FM_melt, (Support>0)))
  GWAS.lead.SNP <- mergedQTL.get_count(subset(FM_melt, (leadSNP==T)))
  GWAS.sig.SNP <- mergedQTL.get_count(subset(FM_melt, (P<=0.05)))
  
  
  merged.data <- data.table::data.table(QTL.Source = consensus$QTL.Source,
                                        Consensus = consensus$Proportion,
                                        Credible.Set = cred.set$Proportion,
                                        GWAS.lead.SNP = GWAS.lead.SNP$Proportion,
                                        GWAS.sig.SNP = GWAS.sig.SNP$Proportion) %>% 
    data.table::melt.data.table(id.vars = "QTL.Source", 
                                variable.name = "SNP.Group", 
                                value.name = "Proportion.Overlap")
  # Plot all groups
  library(ggplot2)
  ggplot(merged.data, aes(x=QTL.Source, y=Proportion.Overlap*100, fill=SNP.Group)) + 
    geom_col(position = position_dodge(preserve = "single")) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 9),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs(title="Proportion of Overlapping QTL", y="% Overlap", x="QTL Source") 
}
 

mergeQTL.stacked_manhattan_plot <- function(FM_melt, SNP.Group=""){  
  GENE_df <- FM_melt
  op <- ggplot(GENE_df, aes(x=POS, y=-log10(FDR), color=-log10(FDR))) + 
    geom_point() + 
    facet_grid(QTL.Source~., drop = F) + 
    theme(strip.text.y = element_text(angle=0), 
          strip.background = element_rect(fill="slategray1"),
          strip.text = element_text(colour = 'black'),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +  
    labs(subtitle=paste0(SNP.Group,": Significant QTL Overlap")) + 
    geom_hline(yintercept = -log10(0.05), linetype="dashed", alpha=.9, color="purple", size=.2) + 
    ylim(c(0, -log10(min(GENE_df$FDR))*1.2  )) + 
    geom_point(data = subset(GENE_df, FDR>=0.05), aes(x=POS, y=-log10(FDR)), color="gray")
  return(op)
}

mergeQTL.QTL_distributions_plot <- function(){
  FM_merge <- data.table::fread(file.path("./Data/GWAS/Nalls23andMe_2019/_genome_wide",
                                          "Nalls23andMe_2019.QTL_overlaps.txt.gz"), 
                                nThread = 4) 
  FM_melt <- mergeQTL.melt_FDR(FM_merge)
  
  
  plot.QTL_distributions.subset <- function(SNP.Group="Consensus_SNP", cardiogenics=F){
    FM_melt <-  FM_melt %>% arrange(QTL.Source, FDR)
    if(cardiogenics==F){
      FM_melt <- subset(FM_melt, !(QTL.Source %in% c("Cardiogenics_macrophages","Cardiogenics_monocytes")))
    }
    
    if(SNP.Group=="Consensus_SNP"){
      GENE_df <- subset(FM_melt, (Consensus_SNP==T) & (FDR<=0.05), drop=F) 
    }
    if(SNP.Group=="Credible_Set"){
      GENE_df <- subset(FM_melt, (Support>0) & (FDR<=0.05), drop=F)  
    }
    if(SNP.Group=="GWAS_lead_SNP"){
      GENE_df <- subset(FM_melt, (leadSNP==T) & (FDR<=0.05), drop=F)  
    } 
    if(SNP.Group=="GWAS_sig_SNP"){
      GENE_df <- subset(FM_melt, (P<=0.05) & (FDR<=0.05), drop=F)  
    } 
    GENE_df$SNP <- factor(GENE_df$SNP, levels = unique(GENE_df$SNP), ordered = T) 
    
    # op <- plot.QTL_overlap(GENE_df, SNP.Group = SNP.Group)
    # print(op)
    
    gp <- ggplot(GENE_df , aes(x=QTL.Source, y=-log10(FDR), fill=QTL.Source, group=SNP)) + 
      geom_col(position = position_dodge(preserve = "single"), show.legend = F, aes(group=SNP), alpha=.7) + 
      geom_hline(yintercept = -log10(0.05), linetype="dashed", alpha=.9, color="purple", size=.2) + 
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) + 
      labs(title=NULL, subtitle=paste0(SNP.Group, " : Significant QTL Overlap"), x=NULL) + 
      scale_fill_discrete(drop=F) + 
      scale_x_discrete(drop=F)  
    # ylim(c(0,-log10( min(FM_melt$FDR, na.rm = T)) ))
    # ylim(c(0,150))
    # print(gp)
    return(gp) 
  } 
  consensus <- plot.QTL_distributions.subset(SNP.Group = "Consensus_SNP")
  cred.set <- plot.QTL_distributions.subset(SNP.Group = "Credible_Set")
  lead.snp <- plot.QTL_distributions.subset(SNP.Group = "GWAS_lead_SNP")
  sig.snp <- plot.QTL_distributions.subset(SNP.Group = "GWAS_sig_SNP")
  
  cowplot::plot_grid(consensus, cred.set, lead.snp, sig.snp, ncol = 1, align = "h")
}




##################################################

####----------- Gather QTL Overlap -----------####

##################################################

mergeQTL.flip_alleles <- function(FM_merge){
  printer("mergeQTL:: Checking ref/alt allele matchup between GWAS and QTL data...") 
  flip_check <- toupper(FM_merge$A1) == toupper(FM_merge$QTL.A1)
  printer("mergeQTL:: Flipping alleles at", sum(!flip_check, na.rm = T),"SNPs")
  FM_merge[!flip_check,"QTL.Effect"] <- -FM_merge[!flip_check,"QTL.Effect"]
  return(FM_merge)
}


mergeQTL.merge_handler <- function(FM_all, qtl_file, force_new_subset=F){
    qtl_name <- gsub(".finemap.txt.gz","",basename(qtl_file))
    message("mergeQTL::",qtl_name) 
    
    if(grepl("GTEx_V7",qtl_name)){
      tissue <- gsub("GTEx_V7_","",qtl_name)
      FM_merge <- mergeQTL.GTEx(FM_all, 
                                GTEx_version="GTEx_V7", 
                                fuzzy_search=tissue, 
                                force_new_subset=force_new_subset) 
      
    } else if(grepl("GTEx_V8",qtl_name)){
      tissue <- gsub("GTEx_V8_","",qtl_name)
      FM_merge <- mergeQTL.GTEx(FM_all, 
                                GTEx_version="GTEx_V8", 
                                fuzzy_search=tissue, 
                                force_new_subset=force_new_subset) 
    
    } else if(grepl("psychENCODE",qtl_name)){
      assay = strsplit(qtl_name,"_")[[1]][2]
      FM_merge <- mergeQTL.psychENCODE(FM_all, 
                                       ASSAYS=assay, 
                                       force_new_subset=force_new_subset) 
      
    } else if(grepl("Brain_xQTL_Serve",qtl_name)){
      assay = strsplit(qtl_name, "[.]")[[1]][2]
      FM_merge <- mergeQTL.Brain_xQTL_Serve(FM_all, 
                                            ASSAYS=assay, 
                                            force_new_subset=force_new_subset) 
      
    } else if(grepl("MESA",qtl_name)){
      pop <- strsplit(qtl_name,"_")[[1]][2]
      FM_merge <- mergeQTL.MESA(FM_all, 
                                POPULATIONS=pop, 
                                force_new_subset=force_new_subset) 
      
    } else if(grepl("Fairfax_2014",qtl_name)){
      condition <- strsplit(qtl_name,"_")[[1]][3]
      FM_merge <- mergeQTL.Fairfax(FM_all, 
                                   CONDITIONS=condition, 
                                   force_new_subset=force_new_subset) 
      
    } else if(grepl("Cardiogenics",qtl_name)){
      celltype <- strsplit(qtl_name,"_")[[1]][2]
      FM_merge <- mergeQTL.Cardiogenics(FM_all, 
                                        CELLTYPES=celltype, 
                                        force_new_subset=force_new_subset) 
    } 
    return(FM_merge)
}

mergeQTL <- function(dataset = "./Data/GWAS/Nalls23andMe_2019"){
  # Gather all Fine-mapping results
  FM_all <- merge_finemapping_results(minimum_support = 0, 
                                      include_leadSNPs = T, 
                                      dataset = "./Data/GWAS/Nalls23andMe_2019",
                                      xlsx_path = F) 
  FM_orig <- FM_all 
  # List all datasets
  printer("+ mergeQTL:: Listing available QTL datasets:")
  list_Data_dirs() %>% dplyr::filter(grepl("QTL",type)) %>% dplyr::select(Dataset, type)
  
  # psychENCODE eQTL, cQTL, isoQTL, tQTL, fQTL, HiC: DLPFC
  FM_merge <- mergeQTL.merge_handler(FM_all = FM_all, )
  # FM_merge <- mergeQTL.psychENCODE(FM_all=FM_all, local_files = T,  force_new_subset = T)
  # Fairfax eQTL: monocytes
  FM_merge <- mergeQTL.Fairfax(FM_all, force_new_subset = F)
  # MESA eQTL: monocytes in AFA, CAU, HIS
  FM_merge <- mergeQTL.MESA(FM_merge, force_new_subset = F) 
  # Cardiogenics eQTL: macrophages, monocytes
  FM_merge <- mergeQTL.Cardiogenics(FM_merge, force_new_subset = F, cis_only = T)
  # Brain_xQTL_Serve
  FM_merge <- mergeQTL.Brain_xQTL_Serve(FM_merge, force_new_subset = F)
  ## Write file and compress
  QTL_merged_path <- file.path(dataset,"_genome_wide",
                               "Nalls23andMe_2019.QTL_overlaps.txt")
  data.table::fwrite(FM_merge, QTL_merged_path, sep="\t", nThread = 4)
  R.utils::gzip(QTL_merged_path, overwrite=T, remove=T)
  
  
  #-------------------------------#
  # GTEx eQTL: 49 different tissues 
  FM_GTEx <- mergeQTL.GTEx_merge_subset(FM_all, GTEx_version="GTEx_V7", fuzzy_search="nigra")
  
  # FM_GTEx <- mergeQTL.GTEx(FM_all, fuzzy_search="Brain", GTEx_version = "GTEx_V7", force_new_subset = F)
  QTL_merged_path.GTEx <- file.path("./Data/GWAS/Nalls23andMe_2019/_genome_wide",
                               "Nalls23andMe_2019.GTEx_QTL_overlaps.txt")
  data.table::fwrite(FM_GTEx, QTL_merged_path.GTEx, sep="\t", nThread = 4)
  R.utils::gzip(QTL_merged_path.GTEx, overwrite=T, remove=T)
  
  return(FM_merge)
}

