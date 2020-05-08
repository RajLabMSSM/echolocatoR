

#' Quickstart (global)
#'
#' Developer function to import some fine-mapping results for testing code.
#' Defaults to Nalls et al. (2019) results. Imports files as global variables.
#'
#' @examples
#' quick_finemap(locus="LRRK2", consensus_thresh = 2)
quick_finemap <- function(locus="LRRK2", consensus_thresh = 2){
  gene <<- locus
  # locus <<- locus
  results_path <<- file.path("./Data/GWAS/Nalls23andMe_2019",gene)
  finemap_DT <<- data.table::fread(file.path(results_path, "Multi-finemap/Multi-finemap_results.txt"))
  finemap_DT <<- cbind(Locus=locus, Gene=locus, finemap_DT)
  finemap_DT <<- find_consensus_SNPs(finemap_DT, consensus_thresh = consensus_thresh)

  if(file.exists(file.path(results_path,"plink","UKB_LD.RDS"))){
    LD_matrix <<- readRDS(file.path(results_path,"plink","UKB_LD.RDS"))
  }
  if(file.exists(file.path(results_path,"plink","LD_matrix.RData")) ){
    file.path(results_path,"plink","LD_matrix.RData")
  }
  sub.out <- subset_common_snps(LD_matrix, finemap_DT)
  LD_matrix <<- sub.out$LD
  finemap_DT <<- sub.out$DT
  subset_DT <<- finemap_DT
}

#' Quickstart (local)
#'
#' Developer function to import some fine-mapping results for testing code.
#' Defaults to Nalls et al. (2019) results.
#' Imports `gene` and `results_path` variables as global variables.
#' @examples
#' finemap_DT <- quick_finemap(locus="LRRK2", consensus_thresh = 2)
quick_finemap_soft <- function(locus="LRRK2", consensus_thresh=2){
  gene <<- locus
  # locus <<- locus
  results_path <<- file.path("./Data/GWAS/Nalls23andMe_2019",gene)
  finemap_DT <- data.table::fread(file.path(results_path, "Multi-finemap/Multi-finemap_results.txt"))
  finemap_DT <- find_consensus_SNPs(finemap_DT, consensus_thresh=consensus_thresh)
  finemap_DT <- data.table::data.table(Gene=locus, finemap_DT)
  return(finemap_DT)
}


#' Quickstart merged
#'
#' Developer function to import and merge all fine-mapping results.
#' Defaults to Nalls et al. (2019) results.
#'
#' @return data.frame
#' @examples
#' finemap_DT <- quick_finemap <- (locus="LRRK2", consensus_thresh = 2)
#' dat <- rbind.file.list(file.list = file.list)
#' @export
quick_merged_DT <- function(minimum_support = 1,
                            dataset = "./Data/GWAS/Nalls23andMe_2019",
                            no_no_loci = c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1",
                                           # Tau region
                                           "RETREG3","UBTF","FAM171A2","MAP3K14","CRHR1","MAPT-AS1","KANSL1","NSF","WNT3")){
  merged_DT <- merge_finemapping_results(minimum_support=minimum_support,
                                         include_leadSNPs=T,
                                         dataset = dataset,
                                         xlsx_path=F,
                                         from_storage=T,
                                         consensus_thresh = 2,
                                         haploreg_annotation=F,
                                         biomart_annotation=F,
                                         verbose = F) %>%
    subset(!(Locus %in% no_no_loci))
  return(merged_DT)
}




#' Quickstart (global, many vars)
#'
#' Assign global variables for rapid testing.
#' @family developer functions
quickstart <- function(){
  # reload()
  allResults <<- list()
  gene <<- "LRRK2"
  leadSNP <<- "rs76904798"
  chrom_col <<- "CHR"
  position_col <<- "POS"
  snp_col <<- "RSID"
  pval_col <<- "p"
  effect_col <<- "beta"
  stderr_col <<- "se"
  freq_col <<- "freq"
  MAF_col<<-"calculate"
  A1_col <<- "A1"
  A2_col <<- "A2"
  finemap_method_list <<- c("SUSIE","ABF","FINEMAP","COJO")
  finemap_methods <<- c("SUSIE","FINEMAP")
  method <<- finemap_methods[1]
  method_list <<- c("SUSIE","FINEMAP")
  force_new_LD <<- F
  before_var <<- "P"

  download_reference <<- T#"../1000_Genomes_VCFs/Phase1"
  superpopulation <<- "EUR"
  force_new_subset <<- F
  min_POS <<- NA
  max_POS <<- NA
  min_MAF <<- NA
  file_sep <<- "\t"
  min_r2 <<- F
  LD_block <<- F
  block_size <<- .7
  min_Dprime <<- F
  plink_folder <<- "./Data/GWAS/Nalls23andMe_2019/LRRK2/plink"
  reference <<- "1KG_Phase1"
  LD_reference <<- "1KG_Phase1"
  bp_distance <<- 100000
  n_causal <<- 10
  vcf_URL <<- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
  popDat_URL <<- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
  # chr <<- 12
  vcf_name <<- "ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"
  vcf_folder <<- "./Data/Reference/1000_Genomes"
  query_by <<- "coordinates"
  remove_variants <<- "rs34637584"
  remove_correlates <<- "rs34637584"

  location_sep <<- ":"
  query_by <<- ""
  dataset_name <<- "Nalls23andMe_2019"
  dataset_type <<- "GWAS"
  N_cases_col <<- "N_cases"
  N_controls_col <<- "N_controls"
  proportion_cases <<- "calculate"


  top_SNPs <- Nalls_top_SNPs <- import_topSNPs(
    topSS_path = Directory_info(dataset_name, "topSS"),
    chrom_col = "CHR", position_col = "BP", snp_col="SNP",
    pval_col="P, all studies", effect_col="Beta, all studies", gene_col="Nearest Gene",
    caption= "Nalls et al. (2018) w/ 23andMe PD GWAS Summary Stats",
    group_by_locus = T,
    locus_col = "Locus Number",
    remove_variants = "rs34637584")
  topSNP_sub <<- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),]


  results_path <<- make_results_path(dataset_name, dataset_type, gene)
  subset_path <<- get_subset_path(results_path, gene)
  subset_DT <<- finemap_DT <<- data.table::fread("Data/GWAS/Nalls23andMe_2019/LRRK2/Multi-finemap/Multi-finemap_results.txt", sep="\t")
  multi <<- T
  subtitle=""
  # subset_path <<- 'Data/GWAS/Nalls23andMe_2019/LRRK2/LRRK2_combined_meta_subset.txt'
  # subset_DT <<- data.table::fread(subset_path, sep="\t")
  LD_path <<- file.path(results_path, "plink/LD_matrix.RData")

  fullSS_path <<- Directory_info(dataset_name, "fullSS.local")
  show_plot<<-T
  subtitle<<-NA
  multi<<-T
  LD_SNP<<- NA


  colDict <<- column_dictionary(fullSS_path)
  if(is.na(min_POS)){min_POS <<- topSNP_sub$POS - bp_distance}
  if(is.na(max_POS)){max_POS <<- topSNP_sub$POS + bp_distance}

  # COJO
  conditioned_snps <<- "rs76904798"#"rs34637584"
  excluded_snps <<- "rs34637584"
  min_MAF <<- 0
  GCTA_path <<- "echolocatoR/tools/gcta_1.92.1beta5_mac/bin/gcta64"
  bfiles <<- "plink_tmp/plink"
  # load("Data/GWAS/Nalls23andMe_2019/LRRK2/plink/LD_matrix.RData")
  LD_matrix <<- readRDS(file.path(results_path,"plink/UKB_LD.RDS"))
  # subset_DT <- subset(subset_DT, SNP %in% unique(row.names(LD_matrix), colnames(LD_matrix) ) )
  # subset_DT <- subset_DT[complete.cases(subset_DT),] # Remove any NAs
  # LD_matrix <- LD_matrix[row.names(LD_matrix) %in% subset_DT$SNP,  colnames(LD_matrix) %in% subset_DT$SNP]

  scaled_prior_variance <<- 0.1
  sample_size <<- NA
  freq_cutoff <<- 0.1

  finemap_method_list <<- c("SUSIE","FINEMAP")#  "ABF", "COJO"

  consensus <<- T

  dataset1_path <<- "./Data/GWAS/Nalls23andMe_2019/LRRK2/LRRK2_Nalls23andMe_2019_subset_500kb.txt"
  dataset2_path <<- "./Data/eQTL/MESA_CAU/LRRK2/LRRK2_MESA_CAU_subset.txt"
  # shared_MAF <<- data.table::fread("Data/GWAS/Nalls23andMe_2019/LRRK2/LRRK2_Nalls23andMe_2019_subset_500kb.txt.gz", sep="\t")$MAF
  plot_subtitle <<- "Fairfax (2014) + CD14 eQTL"
  dataset2_proportion_cases <<- 5e-324
  PP_threshold <<- 0.8
  force_new_finemap <<- F
  dataset1_type <<- "cc"
  dataset2_type <<- "quant"
  shared_MAF <<- 1
  which_merge <<-1
  show_plot <<- T

  goshifter_path <<- "./echolocatoR/tools/goshifter"
  permutations <<- 1000
  remove_tmps <<- T
  chromatin_states <<- c("TssA","EnhA1","EnhA2")
  diff_freq <<- 0.1

  paintor_path <<- "./echolocatoR/tools/PAINTOR_V3.0"
  locus_name<<- NULL
  GWAS_datasets <<- dataset_name
  QTL_datasets <<- NULL
  populations <<- "EUR"
  use_annotations <<- F
}

#' @title echolocatoR
#' @family developer functions
#' @description Assign global variables for rapid testing (Alzheimer's Disease version).
quickstart_AD <- function(locus="PTK2B", dataset_name="Kunkle_2019"){
  loci <<- "PTK2B"
  trim_gene_limits <<- F
  dataset_name <<- dataset_name
  dataset_type <<- "GWAS"
  query_by <<-"tabix"
  finemap_method <<- c("ABF","SUSIE","POLYFUN_SUSIE","FINEMAP")
  force_new_subset <<- T
  force_new_LD <<- F
  force_new_finemap <<- T
  # file_sep <<- " "
  fullSS_path <<- Directory_info(dataset_name, "fullSS.local")
  chrom_col <<- "Chromosome"
  position_col <<- "Position"
  snp_col <<- "MarkerName"
  pval_col <<- "Pvalue"
  effect_col <<- "Beta"
  stderr_col <<- "SE"
  A1_col <<- "Effect_allele"
  A2_col <<- "Non_Effect_allele"
  N_cases <<- 21982
  N_controls <<- 41944
  proportion_cases <<- "calculate"
  MAF_col <<- "MAF"
  freq_col <<- "freq"
  N_cases_col<<-"N_cases"
  N_controls_col<<-"N_controls"
  tstat_col <<-"calculate"


  bp_distance <<- 500000
  download_reference <<- T
  LD_reference <<- "UKB"
  superpopulation <<- "EUR"
  LD_block <<- F
  min_MAF <<- 0.001
  PP_threshold <<- .95
  n_causal <<- 5
  remove_tmps <<- F
  # server <<- F
  gene <<- locus
  results_path <<- file.path("Data/GWAS",dataset_name,gene)
  subset_path <<- file.path(results_path,"Multi-finemap" ,paste0(gene,"_",dataset_name,"_Multi-finemap.tsv.gz"))
  # finemap_DT <<- data.table::fread(subset_path)
  finemap_DT <<- data.table::data.table(standardize_subset(subset_path = subset_path, gene=gene))
  subset_DT <<- finemap_DT
  LD_matrix <<- readRDS(file.path(results_path,"plink/UKB_LD.RDS"))
  sub.out <- subset_common_snps(LD_matrix=LD_matrix,
                                finemap_DT=subset_DT)
  LD_matrix <- sub.out$LD
  subset_DT <- sub.out$DT

  polyfun<<-"./echolocatoR/tools/polyfun"
  force_new_priors <<- F
}
