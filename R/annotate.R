# %%%%%%%%%%%%%%%%% #
####### Annotate #######
# %%%%%%%%%%%%%%%%% #


#' Merge fine-mapping results from all loci
#'
#' Gather fine-mapping results from \emph{echolocatoR} across all loci
#' and merge into a single data.frame.
#'
#' @param dataset Path to the folder you want to recursively search for results files within
#'  (e.g. \url{"Data/GWAS/Nalls23andMe_2019"}).
#' Set this to a path that includes multiple subfolders if you want to gather results
#' from multiple studies at once
#' (e.g. \url{"Data/GWAS"}).
#' @param minimum_support Filter SNPs by the minimum number
#'  of fine-mapping tools that contained the SNP in their Credible Set.
#' @param include_leadSNPs Include lead GWAS/QTL SNPs per locus
#'  (regardless of other filtering criterion).
#' @param xlsx_path Save merged data.frame as excel file.
#' @param from_storage Search for stored results files.
#' @param haploreg_annotation Annotate SNPs with HaploReg (using \code{HaploR}).
#' @param regulomeDB_annotation Annotate SNPs with regulaomeDB (using \code{HaploR}).
#' @param biomart_annotation Annotate SNPs with \code{biomart}.
#' @param PP_threshold Mean posterior probability threshold to include SNPs in mean PP Credible Set
#'  (averaged across all fine-mapping tools).
#' @param consensus_thresh The minimum number of tools that have the SNPs in their Credible Set
#' to classify it as a \strong{Consensus_SNP}.
#' @param exclude_methods Exclude certain fine-mapping methods when estimating
#' \strong{mean.CS} and \strong{Consensus_SNP}.
#' @param verbose Print messages.
#' @family annotatate
#' @examples
#' dataset_dir <- "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"
#' # UCS and lead SNPs: No annotation
#' merged_DT <- merge_finemapping_results(dataset=dataset_dir, minimum_support=1, include_leadSNPs=T)
#'
#' # UCS and lead SNPs: With annotations
#' merged_DT <- merge_finemapping_results(dataset=dataset_dir, minimum_support=1, include_leadSNPs=T, haploreg_annotation=T, biomart_annotation=T)
merge_finemapping_results <- function(dataset="./Data/GWAS",
                                      minimum_support=1,
                                      include_leadSNPs=T,
                                      LD_reference=NULL,
                                      #="./Data/annotated_finemapping_results.xlsx",
                                      xlsx_path=F,
                                      from_storage=T,
                                      haploreg_annotation=F,
                                      regulomeDB_annotation=F,
                                      biomart_annotation=F,
                                      PP_threshold=.95,
                                      consensus_threshold=2,
                                      exclude_methods=NULL,
                                      top_CS_only=F,
                                      verbose=T,
                                      nThread=4){
  if(from_storage){
    printer("+ Gathering all fine-mapping results from storage...", v=verbose)
    # Find all multi-finemap_results files
    multifinemap_pattern <- create_method_path(locus_dir = file.path(dataset,"*"),
                                               finemap_method = "Multi-finemap",
                                               LD_reference = LD_reference,
                                               compress = F,
                                               create_dir = F)
    multifinemap_pattern_gz <- create_method_path(locus_dir = file.path(dataset,"*"),
                                                  finemap_method = "Multi-finemap",
                                                  LD_reference = LD_reference,
                                                  compress = T,
                                                  create_dir = F)
    multi_dirs <- list.files(path = dataset,
                             pattern = paste0( c(basename(multifinemap_pattern),
                                                 basename(multifinemap_pattern_gz)), collapse="|"),
                             recursive = T,
                             full.names = T)
    if(length(multi_dirs)==0){
      stop("+ No multi-finemap files found.")
    } else {
      printer("+",length(multi_dirs),"multi-finemap files found.",v=verbose)
    }
    loci <- basename(dirname(dirname(multi_dirs)))
    if(length(loci)>length(unique(loci))){
      printer("+ Removing duplicate Multi-finemap files per locus.",v=verbose)
      loci <- loci[!duplicated(loci)]
    }
    dataset_names <- dirname(dirname(dirname(multi_dirs))) %>% unique()
    # Loop through each GENE
    finemap_results <- lapply(dataset_names, function(dn, multi_dirs.=multi_dirs){
      # gene_dirs <- dirname(dirname(multi_dirs.))
      # Loop through each gene folder
      all_results <- lapply(multi_dirs., function(md){
        gene <- basename(dirname(dirname(md)))
        printer("+ Importing results...",gene, v=verbose)
        multi_data <- data.table::fread(md, nThread = nThread)
        multi_data <- update_cols(multi_data)
        multi_data <- cbind(data.table::data.table(Dataset=dn, Gene=gene), multi_data)
        return(multi_data)
      }) %>% data.table::rbindlist(fill=TRUE) # Bind genes
    }) %>% data.table::rbindlist(fill=TRUE) # Bind datasets
  }


  # Add/Update Support/Consensus cols
  merged_results <- find_consensus_SNPs(finemap_dat = finemap_results,
                                        credset_thresh = PP_threshold,
                                        consensus_thresh = consensus_threshold,
                                        exclude_methods = exclude_methods,
                                        top_CS_only = top_CS_only,
                                        verbose = verbose)
  merged_results <- subset(merged_results, Support>=minimum_support)
  if(!"Locus" %in% colnames(merged_results)){
    merged_results <- merged_results %>% dplyr::rename(Locus=Gene) %>% data.table::data.table()
  }

  # Loop through each DATASET
  # merged_results <- lapply(unique(finemap_results$Dataset), function(dname, include_leadSNPs.=include_leadSNPs){
  #   multi_finemap_dat <- subset(finemap_results, Dataset==dname)
  #   CS_cols <- colnames(multi_finemap_dat)[endsWith(colnames(multi_finemap_dat), ".CS")]
  #   finemap_method_list <- gsub(".CS","",CS_cols)
  #   # Create support table
  #   support_DT <- multi_finemap_results_table(multi_finemap_dat,
  #                                             finemap_method_list,
  #                                             fancy_table = F,
  #                                             minimum_support = minimum_support,
  #                                             include_leadSNPs = include_leadSNPs.)
  #   support_DT <- cbind(Dataset = dname, support_DT) %>% arrange(Gene, desc(Support))
  # }) %>% data.table::rbindlist()




  # Annotate with haplorR
  if(haploreg_annotation){
    HR_query <- haploR.HaploReg(snp_list = unique(merged_results$SNP), verbose = verbose)
    merged_results <- data.table:::merge.data.table(merged_results, HR_query,
                                                    by.x = "SNP",
                                                    by.y = "rsID",
                                                    all = T,
                                                    allow.cartesian=TRUE)
  }
  if(regulomeDB_annotation){
    regDB_query <- haploR.regulomeDB(snp_list = unique(merged_results$SNP), verbose = verbose)
    merged_results <- data.table:::merge.data.table(merged_results, regDB_query,
                                                    by.x = "SNP",
                                                    by.y = "rsID",
                                                    all = T,
                                                    allow.cartesian=TRUE)
  }
  if(biomart_annotation){
    biomart_query <- biomart_snp_info(snp_list = merged_results$SNP, verbose = verbose)
    merged_results <- data.table:::merge.data.table(merged_results, biomart_query,
                                                      by.x = "SNP",
                                                      by.y = "refsnp_id",
                                                      all = T,
                                                      allow.cartesian=TRUE)
  }
  if(xlsx_path!=F){
    # data.table::fwrite(merged_results, file = csv_path, quote = F, sep = ",")
    openxlsx::write.xlsx(merged_results, xlsx_path)
  }
  # createDT_html(merged_results) %>% print()
  return(merged_results)
}






counts_summary <- function(top_SNPs,
                           merged_results,
                           verbose=T){
  # Get total # of SNPs per gene per dataset
  candidate_counts <- merged_results %>% dplyr::group_by(Dataset, Locus) %>%
    count(name = "Total_SNPs")
  max_consensus <- sum(endsWith(colnames(merged_results),".CS"))
  candidate_counts <- merge.data.frame(candidate_counts,
                                       dplyr::group_by(merged_results, Locus) %>%
                                         count(name = "CS"),
                                       by = "Locus", all = T)
  candidate_counts <- merge.data.frame(candidate_counts,
                                       subset(merged_results, Support == max_consensus) %>%
                                         dplyr::group_by(Locus) %>%
                                         count(name = "Consensus_SNPs"),
                                       all = T)
  # Add lead rsid column
  candidate_counts <-  merge.data.frame(candidate_counts,
                                        top_SNPs[,c("Locus","SNP")] %>% dplyr::rename(leadSNP=SNP),
                                        on = "Locus", all = T)
  # Add "is lead SNP in Credible Set" column
  candidate_counts <-  merge.data.frame(candidate_counts,
                                        leadSNP_comparison(top_SNPs, merged_results),
                                        on = "Locus", all = T)
  # Gather credible set rsids
  CredSet_rsids <- merged_results %>% dplyr::group_by(Dataset, Locus) %>%
    subset(Support==max_consensus) %>%
    dplyr::summarise(CredSet_rsids = paste(SNP,collapse="; ")) %>%
    data.table::data.table()
  candidate_counts <- merge.data.frame(candidate_counts,
                                       CredSet_rsids,
                                       on = "Locus", all = T)
  # Gather consensus rsids
  consensus_rsids <- merged_results %>% dplyr::group_by(Dataset, Locus) %>%
    subset(Support==T) %>%
    dplyr::summarise(Consensus_rsids = paste(SNP,collapse="; ")) %>%
    data.table::data.table()
  candidate_counts <- merge.data.frame(candidate_counts,
                                       consensus_rsids,
                                       on = "Locus", all = T)
  # Fill 0s
  candidate_counts$Consensus_SNPs[is.na(candidate_counts$Consensus_SNPs)] <- 0
  means <- c(Locus=" ",
             Dataset="[Average]",
             candidate_counts[,c("Total_SNPs","CS","Consensus_SNPs")] %>% colMeans() %>% round(1),
             leadSNP_in_CredSet = paste0(round(sum(candidate_counts$leadSNP_in_CredSet) / dim(candidate_counts)[1]*100,2),"%"),
             CredSet_rsids = "",
             Conensus_rsids = ""
  )
  # Add averages
  candidate_counts <- suppressWarnings(rbind(candidate_counts,  means))
  percent_model_convergence <- round(sum(candidate_counts$Consensus_SNPs>0)  / length(candidate_counts$Consensus_SNPs) *100, 2)
  max_consensus_set <- max(candidate_counts$Consensus_SNPs)
  # Check if lead SNP is in the credible sets for each locus
  printer("\n + All",max_consensus,"models converged upon 1 to",
          max_consensus_set,"SNPs in",percent_model_convergence,"% of loci.",
          v=verbose)
  # createDT_html(candidate_counts) %>% print()
  return(candidate_counts)
}









#' Download SNP-wise annotations from HaploReg
#'
#' @keywords internal
#' @family annotate
#' @source
#' \href{https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html}{HaploR}
haploR.HaploReg <- function(snp_list, verbose=T, chunk_size=NA){
  printer("+ Gathering annotation data from HaploReg...", v=verbose)
  # Break into smaller chunks
  snp_list <- unique(snp_list)
  if(is.na(chunk_size)){chunk_size <- length(snp_list)}
  chunked_list <- split(snp_list, ceiling(seq_along(snp_list)/chunk_size))

  HR_query <- lapply(names(chunked_list), function(i){
    printer("++ Submitting chunk",i,"/",length(chunked_list))
    chunk <- chunked_list[[i]]
    HR_query <-  haploR::queryHaploreg(query = chunk, file = NULL, study = NULL, ldThresh = NA,
                                       ldPop = "EUR", epi = "vanilla", cons = "siphy", genetypes = "gencode",
                                       url = "https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php",
                                       timeout = 500, encoding = "UTF-8", verbose = FALSE)

    return(data.table::as.data.table(HR_query))
  }) %>% data.table::rbindlist()

  return(HR_query)
}




#' Download SNP-wise annotations from RegulomeDB
#'
#' @keywords internal
#' @family annotate
#' @source
#' \href{https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html}{HaploR}
haploR.regulomeDB <- function(snp_list, verbose=T, chunk_size=NA){
  printer("+ Gathering annotation data from HaploReg...", v=verbose)
  # Break into smaller chunks
  snp_list <- unique(snp_list)
  if(is.na(chunk_size)){chunk_size <- length(snp_list)}
  chunked_list <- split(snp_list, ceiling(seq_along(snp_list)/chunk_size))

  rDB_query <- lapply(names(chunked_list), function(i){
    printer("++ Submitting chunk",i,"/",length(chunked_list))
    chunk <- chunked_list[[i]]
    rdb_query <-  haploR::queryRegulome(query = chunk,
                                       timeout = 500,
                                       verbose = F)
    return(data.table::as.data.table(rdb_query))
  }) %>% data.table::rbindlist()

  return(rDB_query)
}




#' Download SNP-wise annotations from Biomart
#'
#' @keywords internal
#' @family annotate
#' @source
#' \href{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}{biomaRt}
biomart_snp_info <- function(snp_list,
                             reference_genome="grch37",
                             verbose=T){
  printer("+ Gathering annotation data from Biomart...", v=verbose)
  mart = biomaRt::useMart(biomart="ENSEMBL_MART_SNP",
                           host=paste0(reference_genome,".ensembl.org"),
                           path="/biomart/martservice",
                           dataset="hsapiens_snp",)
  # View(biomaRt::listFilters(mart))
  # View(biomaRt::listAttributes(mart))
  biomart_query = biomaRt::getBM(attributes =  c('refsnp_id',
                                 'allele',
                                 'chr_name',
                                 'chrom_start',
                                 'chrom_end',
                                 'chrom_strand',
                                 'ensembl_gene_stable_id',
                                 'consequence_type_tv',
                                 'polyphen_prediction',
                                 'polyphen_score',
                                 'sift_prediction',
                                 'sift_score',
                                 'reg_consequence_types',
                                 'validated'
                                 ),
                  filters = c('snp_filter'),
                  values = unique(snp_list),
                  mart = mart)
  biomart_query <- data.table::as.data.table(biomart_query)
  biomart_query[biomart_query$consequence_type_tv=="",]$consequence_type_tv <- NA
  # Only take the first annotation per variant
  # annotated_results %>% dplyr::group_by(Dataset, Gene, SNP) %>% slice(1)
  return(biomart_query)
}





#' Identify which genes SNPs belong to using Biomart
#'
#' @keywords internal
#' @family annotate
#' @source
#' \href{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}{biomaRt}
biomart_snps_to_geneInfo <- function(snp_list, reference_genome="grch37"){
  # listMarts()
  snp_mart = useMart("ENSEMBL_MART_SNP",
                     dataset="hsapiens_snp",
                     host =  paste0(reference_genome,".ensembl.org"))
  # View(listFilters(snp_mart))
  # View(listAttributes(snp_mart))
  snp_results <- biomaRt::getBM(snp_mart, filters="snp_filter",
                                values=snp_list,
                                attributes=c("refsnp_id","snp","chr_name", "chrom_start","chrom_end",
                                             "associated_gene","ensembl_gene_stable_id" ) )
  # # Split ensembl IDs
  gene_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  gene_results <- biomaRt::getBM(mart = gene_mart,
                                 filters = "ensembl_gene_id",
                                 # values = unlist(strsplit(snp_results$ensembl, ";")),
                                 values = snp_results$ensembl_gene_stable_id,
                                 attributes = c("hgnc_symbol","external_gene_name","ensembl_gene_id",
                                                "chromosome_name", "start_position", "end_position") )
  snp_results <-snp_results %>%
    mutate(ensembl = strsplit(as.character(ensembl_gene_stable_id), ";")) %>%
    tidyr::unnest(ensembl)
  merged_df <- data.table(gene_results, key = "ensembl_gene_id")[data.table(snp_results, key = "ensembl")]
  return(merged_df)
}
# biomart_snps_to_geneInfo(c("rs114360492"))




#' Get gene info using Biomart
#'
#' @keywords internal
#' @family annotate
#' @source
#' \href{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}{biomaRt}
#' @examples
#' gene_info <- biomart_geneInfo(c("PTK2B","CLU","APOE"))
biomart_geneInfo <- function(geneList,
                             reference_genome="grch37"){
  # listDatasets(useMart("ENSEMBL_MART_ENSEMBL") )
  gene_mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                               dataset="hsapiens_gene_ensembl",
                               host = paste0(reference_genome,".ensembl.org"))
  # View(listFilters(gene_mart))
  # View(listAttributes(gene_mart))
  gene_results <- biomaRt::getBM(mart = gene_mart,
                                 filters = "hgnc_symbol",
                                 # values = unlist(strsplit(snp_results$ensembl, ";")),
                                 values = geneList,
                                 attributes = c("hgnc_symbol","external_gene_name","ensembl_gene_id",
                                                "chromosome_name", "start_position", "end_position") )
  return(gene_results)
}



#' Return only the missense SNPs
#'
#' @keywords internal
#' @family annotate
#' @source
#' \href{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}{biomaRt}
SNPs_by_mutation_type <- function(merged_results,
                                  mutation_type="missense_variant"){
  potential_missense <- subset(merged_results, consequence_type_tv == mutation_type) %>%
    dplyr::group_by(Dataset, Gene, SNP) %>%
    dplyr::select(Dataset, Gene, SNP, consequence_type_tv) %>%
    unique()
  potential_missense_full <- subset(merged_results, SNP %in% potential_missense$SNP) %>%
    dplyr::group_by(Dataset, Gene, SNP) %>%
    dplyr::select(Dataset, Gene, SNP, consequence_type_tv) %>%
    unique()
  return(potential_missense_full)
}



#' Summarise \code{HaploR} annotations
#'
#' @keywords internal
#' @family annotate
#' @source
#' \href{https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html}{HaploR}
epigenetics_summary <- function(merged_results,
                                tissue_list = c("BRN","BLD"),
                                # Chromatin_Marks
                                epigenetic_variables = c("Promoter_histone_marks","Enhancer_histone_marks")
                                ){
  merged_results <- data.table(merged_results)
  summary_func <- function(ev){
    boolean <- lapply(ev, function(x){ intersect(tissue_list, strsplit(as.character(x), ", ")[[1]]) %>% length() > 0 }) %>% unlist()
    n_hits <- dim(merged_results[boolean,])[1]
    Total <- dim(merged_results)[1]
    Percent_Total <- round(n_hits / Total*100,2)
    return(list(Hits = n_hits,
                Total = Total,
                Percent_Total = Percent_Total))
  }
  epi_summary <- merged_results[, lapply(.SD, summary_func), .SDcols = epigenetic_variables] %>% t()
  colnames(epi_summary) <- c("Hits","Total_SNPs","Percent_Total")
  print(epi_summary)
}



#' Test for enrichment of \code{HaploR} annotations
#'
#' @keywords internal
#' @family annotate
#' @source
#' \href{https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html}{HaploR}
epigenetics_enrichment <- function(snp_list1,
                                   snp_list2,
                                   chisq=T,
                                   fisher=T,
                                   epigenetic_variables = c("Promoter_histone_marks",
                                                            "Enhancer_histone_marks"),
                                   tissue_list = c("BRN","BLD")){
  printer("Conducting SNP epigenomic annotation enrichment tests...")
  # Foreground
  printer("+++ SNP list 1 :")
  HR1 <- haploR.HaploReg(snp_list1)
  summ1 <- epigenetics_summary(HR1, tissue_list = tissue_list)
  # Background
  printer("+++ SNP list 2 :")
  HR2 <- haploR.HaploReg(snp_list2)
  summ2 <- epigenetics_summary(HR2, tissue_list = tissue_list)

  for(epi_var in epigenetic_variables){
    printer("++ Testing for enrichment of '",epi_var,
            "' in the tissues '",paste(tissue_list, collapse=" & "),"'")
    # Create contingency table
    cont_tab <- rbind(list1 = summ1[epi_var, c("Hits","Total_SNPs")] %>% unlist,
                      list2 = summ2[epi_var, c("Hits","Total_SNPs")] %>% unlist() ) %>% as.table()
    # Conduct tests
    if(chisq){
      chisq.results <- chisq.test(cont_tab, simulate.p.value = TRUE)
      print(chisq.results)
    }
    if(fisher){
      fisher.results <- fisher.test(cont_tab)
      print(fisher.results)
    }
  }
}




#' Annotate any missense variants
#'
#' @family annotate
#' @examples
#' data("merged_DT");
#' annotated_DT <- ANNOTATE.annotate_missense(merged_DT=merged_DT, snp_filter="Support>0")
ANNOTATE.annotate_missense <- function(merged_DT,
                                       snp_filter="Support>0"){
  snp_info <- biomart_snp_info(unique(subset(merged_DT, eval(parse(text=snp_filter)))$SNP))
  # unique(snp_info$consequence_type_tv)
  missense <- suppressMessages(snp_info %>%
                                 dplyr::group_by(refsnp_id) %>%
                                 dplyr::summarise(Missense = ifelse(any(consequence_type_tv=="missense_variant",na.rm = T),T,F)) %>%
                                 data.table::data.table())

  merged_DT <- data.table::merge.data.table(merged_DT, missense,
                                            all.x=T,
                                            by.x="SNP",
                                            by.y = "refsnp_id")
  # missense_counts <- suppressMessages(merged_DT %>% dplyr::group_by(Locus) %>%
  #   dplyr::summarise(Missense=sum(Missense, na.rm=T)))
  printer(sum(subset(merged_DT, Support>0)$Missense, na.rm = T),"missense mutations detected in UCS.")
  printer(sum(subset(merged_DT, Consensus_SNP)$Missense, na.rm = T),"missense mutations detected in Consensus SNPs")
  return(merged_DT)
}


#' Plot any missense variants
#'
#' @family annotate
#' @examples
#' data("merged_DT");
#' gg_missense <- ANNOTATE.plot_missense(merged_DT=merged_DT, snp_filter="Support>0")
#' gg_missense <- ANNOTATE.plot_missense(merged_DT=merged_DT, snp_filter="Consensus_SNP==T")
ANNOTATE.plot_missense <- function(merged_DT,
                                   snp_filter="Support>0",
                                   label_yaxis=F,
                                   x_label="UCS missense\nmutations",
                                   show.legend=T,
                                   show_numbers=F,
                                   show_plot=T){
  locus_order <- SUMMARISE.get_CS_counts(merged_DT = merged_DT)
  annotated_DT <- ANNOTATE.annotate_missense(merged_DT=merged_DT,
                                             snp_filter=snp_filter)
  dat_melt <-
    data.table::setDT(annotated_DT)[, .(Missense = n_distinct(SNP[Missense==T], na.rm = T)),
                        by=c("Locus")]  %>%
    dplyr::mutate(Locus=factor(Locus, levels = unique(locus_order$Locus), ordered = T),
                  Missense = as.integer(Missense))
  dat_melt$dummy <- x_label
  dat_melt[dat_melt$Missense==0,"Missense"] <- NA

  gg_missense <- ggplot(data=dat_melt, aes(x=dummy, y=Locus, fill=Missense)) +
    geom_tile(show.legend = show.legend, alpha=.7) +
    # scale_fill_continuous(na.value = "transparent") +
    # scale_fill_gradient(low = scales::alpha("blue",.7),
    #                     high = scales::alpha("red",.7),
    #                     na.value = "transparent",
    #                     n.breaks=max(dat_melt$Missense, na.rm = T)) +
    # scale_fill_fermenter(palette = "Spectral", na.value = "transparent") +
    scale_fill_viridis_c(na.value = "transparent") +
    theme_bw() +
    labs(x=NULL, fill="Missense\nmutations") +
    theme(axis.text.x = element_text(angle=40, hjust=1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "top",
          legend.title.align = .5,
          legend.text.align = .5,
          legend.box="horizontal",
          legend.key = element_rect(colour = "gray60"))
  if(show_numbers){
    gg_missense <- gg_missense +
      geom_text(aes(label = Missense), color="grey70")
  }
  if(label_yaxis==F){
    gg_missense <- gg_missense + theme(axis.text.y = element_blank()) +
      labs(y=NULL)
  }
  if(show_plot)print(gg_missense)
  return(list(plot=gg_missense,
              data=dat_melt))
}

