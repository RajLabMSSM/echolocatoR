# %%%%%%%%%%%%%%%%% #
####### QUERY #######
# %%%%%%%%%%%%%%%%% #


#' Import top GWAS/QTL summary statistics
#'
#' The resulting top_SNPs data.frame can be used to guide
#' the \code{\link{finemap_loci}} in querying and fine-mapping loci.
#'
#' @param top_SNPs Can be a data.frame with the top sumary stats per locus.
#' Alternatively, you can provide a path to the stored top summary stats file.
#' Can be in any tabular format (e.g. excel, .tsv, .csv, etc.).
#' This file should have one lead GWAS/QTL hits per locus.
#' If there is more than one SNP per locus, the one with the smallest p-value (then the largest effect size) is selected as the lead SNP.
#' The lead SNP will be used as the center of the locus when constructing the locus subset files.
#'
#' @family query functions
#' @param show_table Create an interative data table.
#' @param sheet If the \emph{topSS} file is an excel sheet, you can specify which tab to use.
#' You can provide either a number to identify the tab by order,
#' or a string to identify the tab by name.
#' @param gene_col An optional column to keep track of the causal gene(s) in each locus (if known).
#' @param grouping_vars The variables that you want to group by
#' such that each grouping_var combination has its own index SNP.
#' For example, if you want one index SNP per QTL eGene - GWAS locus pair, you could supply:
#' \code{grouping_vars=c("Locus","Gene")}.
#' @inheritParams finemap_pipeline
import_topSNPs <- function(topSS,
                           show_table=T,
                           sheet = 1,
                           chrom_col="CHR",
                           position_col="POS",
                           min_POS_col=NULL,
                           max_POS_col=NULL,
                           snp_col="SNP",
                           pval_col="P",
                           effect_col="Effect",
                           locus_col="Locus",
                           grouping_vars=c("Locus"),
                           gene_col="Gene",
                           remove_variants=NULL,
                           nThread=4,
                           verbose=T){
  # Import top SNPs
  topSNPs_reader <- function(topSS,
                             sheet = 1){
    if(is.data.frame(topSS)){
      return(topSS)
    } else{
      if(endsWith(topSS, ".xlsx") | endsWith(topSS, ".xlsm")){
        topSS <- openxlsx::read.xlsx(topSS,
                                     sheet = sheet,
                                     # sep.names=" ", # Only in some versions?
                                     check.names = F)
      } else {
        topSS <- data.table::fread(file=topSS,
                                   header = T,
                                   stringsAsFactors = F,
                                   nThread=nThread)
      }
      return(topSS)
    }
  }

  top_SNPs <- topSNPs_reader(topSS, sheet)
  orig_top_SNPs <- top_SNPs
  # Reassign colnames bc read.xlsx won't let you prevent this in some versions....
  if(!is.data.frame(topSS)){
    if(endsWith(topSS, ".xlsx") | endsWith(topSS, ".xlsm")){
      locus_col=gsub(" ",".",trimws(locus_col))
      gene_col=gsub(" ",".",trimws(gene_col))
      chrom_col=gsub(" ",".",trimws(chrom_col))
      position_col=gsub(" ",".",trimws(position_col))
      snp_col=gsub(" ",".",trimws(snp_col))
      pval_col=gsub(" ",".",trimws(pval_col))
      effect_col=gsub(" ",".",trimws(effect_col))
    }
  }

  # Fill in missing cols with nonsense
  if(is.na(chrom_col)) {top_SNPs$CHR <- NA; chrom_col<-"CHR";}
  if(is.na(position_col)) {top_SNPs$POS <- NA; position_col<-"POS";}
  if(!effect_col %in% colnames(top_SNPs)) {top_SNPs$Effect <- 1; effect_col<-"Effect";}

  top_SNPs <- dplyr::rename(top_SNPs,
                            CHR=chrom_col,
                            POS=position_col,
                            SNP=snp_col) %>%
    dplyr::mutate(CHR=gsub("chr","",CHR))

  # Add Locus/Gene columns
  bad_chars <- "/|[:]|[,]"
  if((is.null(gene_col) & is.null(locus_col))){
    printer("+ Constructing locus names from CHR and index SNP")
    # locus_chr1_rs10737496
    top_SNPs <- dplyr::mutate(top_SNPs,
                              Locus=paste0("locus_chr",CHR,"_",SNP),
                              Gene=paste0("locus_chr",CHR,"_",SNP))
    locus_col <- gene_col <- "Locus";
  } else if(!any(gene_col %in% colnames(top_SNPs)) & !any(locus_col %in% colnames(top_SNPs))){
    printer("+ Constructing locus names from CHR and index SNP")
    # locus_chr1_rs10737496
    top_SNPs <- dplyr::mutate(top_SNPs,
                              Locus=paste0("locus_chr",CHR,"_",SNP),
                              Gene=paste0("locus_chr",CHR,"_",SNP))
    locus_col <- gene_col <- "Locus";
  } else if(gene_col %in% colnames(top_SNPs) & locus_col %in% colnames(top_SNPs)){
    print("+ Assigning gene_col and locus_col independently",v=verbose)
    top_SNPs$Gene <- gsub(bad_chars,"_", orig_top_SNPs[[gene_col]])
    top_SNPs$Locus <- gsub(bad_chars,"_",orig_top_SNPs[[locus_col]]) # Get rid of problematic characters

  } else {
    if(gene_col %in% colnames(top_SNPs) & all(!is.na(top_SNPs[[gene_col]])) ){
      print("+ Assigning gene_col to locus_col",v=verbose)
      top_SNPs$Gene <- gsub(bad_chars,"_", orig_top_SNPs[[gene_col]])
      top_SNPs$Locus <- gsub(bad_chars,"_",top_SNPs$Gene) # Get rid of problematic characters
    } else {
      print("+ Assigning locus_col to gene_col",v=verbose)
      top_SNPs$Locus <- gsub(bad_chars,"_",orig_top_SNPs[[locus_col]])
      top_SNPs$Gene <- gsub(bad_chars,"_",top_SNPs$Locus) # Get rid of problematic characters
    }
  }


  top_SNPs <- suppressMessages(top_SNPs %>%
    dplyr::select(Locus=Locus,
                  Gene=Gene,
                  CHR=CHR,
                  POS=POS,
                  SNP=SNP,
                  P=pval_col,
                  Effect=effect_col,
                  min_POS=min_POS_col,
                  max_POS=max_POS_col))
    if("min_POS" %in% colnames(top_SNPs) &
       "max_POS" %in% colnames(top_SNPs)){
      top_SNPs <- dplyr::mutate(top_SNPs, span_kb=(max_POS-min_POS)/1000)
    }

    # Remove specific variants
    if(!is.null(remove_variants)){
      top_SNPs <- subset(top_SNPs, !(SNP %in% remove_variants))
    }

    # Get the top representative SNP and Gene per locus (by lowest p-value and effect size)
   if(!is.null(grouping_vars)){
      top_SNPs <- suppressWarnings(top_SNPs %>%
        dplyr::arrange(P, dplyr::desc(Effect)) %>%
        dplyr::group_by(.dots=grouping_vars) %>%
        dplyr::slice(1) %>%
        replace(., .=="NA", NA) %>%
        subset(!is.na(Locus)) %>%
        dplyr::mutate(CHR=as.numeric(gsub("chr", "",CHR))) )
    }

  # Make sure cols are numeric
  top_SNPs <- top_SNPs %>% dplyr::mutate_at(.vars = vars(POS,P,Effect),
                                            .funs = function(x){as.numeric(gsub(",| ","",x))})
  if(show_table){
    createDT(top_SNPs, caption = "Top SNP per locus")
  }
  return(data.table::data.table(top_SNPs))
}




#' Generate a named list of [e]gene-locus pairs
#'
gene_locus_list <- function(top_SNPs){
  setNames(top_SNPs$Locus, top_SNPs$Gene)
}


#' Detect QTL genes in full summary stats file
#'
#' Allows summary stats from different genes to be
#' fine-mapped separately.
#' @examples
#' loci <- c("BST1","LRKR2","MEX3C")
#' detect_genes(loci)
#' loci <- c(BST1="BST1", LRRK2="LRRK2", MEX3C="MEX3C")
#' detect_genes(loci)
detect_genes <- function(loci,
                          verbose=T){
  if(!is.null(names(loci))){
    printer("Fine-mapping",dplyr::n_distinct(loci),"gene:Locus pairs.", v=verbose)
    return(T)
  } else {printer("Fine-mapping",dplyr::n_distinct(loci),"loci.", v=verbose) }
  return(F)
}




#' Extract a subset of the summary stats
#'
#' Use either \emph{tabix} or \emph{awk} to extract a locus subset
#'  from the full summary statistics file.
#'
#' @family query functions
#' @keywords internal
extract_SNP_subset <- function(locus=NULL,
                               locus_dir,
                               results_dir=NULL,
                               fullSS_path,
                               fullSS_genome_build="hg19",
                               subset_path,
                               LD_reference,
                               force_new_subset=F,
                               top_SNPs="auto",
                               bp_distance=500000,
                               chrom_col="CHR",
                               chrom_type=NULL,
                               position_col="POS",
                               snp_col="SNP",
                               locus_col="Locus",
                               pval_col="P",
                               effect_col="Effect",
                               stderr_col="StdErr",
                               MAF_col="MAF",
                               freq_col = "Freq",
                               tstat_col="t-stat",
                               A1_col = "A1",
                               A2_col = "A2",
                               gene_col="gene",

                               N_cases_col="N_cases",
                               N_controls_col="N_controls",
                               N_cases=NULL,
                               N_controls=NULL,
                               proportion_cases="calculate",
                               sample_size=NULL,

                               superpopulation="",
                               min_POS=NA,
                               max_POS=NA,
                               genes_detected=F,

                               file_sep="\t",
                               query_by="coordinates",
                               probe_path = "./Data/eQTL/gene.ILMN.map",
                               QTL_prefixes=NULL,
                               remove_tmps=T,
                               conda_env = "echoR",
                               verbose=T){
  if(is.null(locus)) locus <- basename(locus_dir)
  multi_path <- create_method_path(locus_dir = locus_dir,
                                   LD_reference = LD_reference,
                                   finemap_method = "Multi-finemap",
                                   compress = T)

  if(file.exists(subset_path) & force_new_subset==F){
    printer("+ Importing pre-existing file:",subset_path, v=verbose)
    check_if_empty(subset_path)
    query <- data.table::fread(subset_path, header=T, stringsAsFactors=F)
  } else if (file.exists(multi_path) & force_new_subset==F){
    printer("+ Importing  pre-existing file:",multi_path, v=verbose)
    check_if_empty(multi_path)
    query <- data.table::fread(multi_path, header=T, stringsAsFactors=F)
  } else {
    # Extract subset with awk
    printer("+ Extracting relevant variants from fullSS...", v=verbose)
    start_query <- Sys.time()
    # Function selects different methods of querying your SNPs
    query_handler(locus_dir=locus_dir,
                  fullSS_path=fullSS_path,
                  subset_path=subset_path,
                  top_SNPs=top_SNPs,
                  locus_col=locus_col,
                  chrom_col=chrom_col,
                  chrom_type=chrom_type,
                  position_col=position_col,
                  file_sep=file_sep,
                  min_POS=min_POS,
                  max_POS=max_POS,
                  bp_distance=bp_distance,
                  query_by=query_by,
                  probe_path=probe_path,
                  conda_env = conda_env,
                  verbose=verbose)
    # Clean file
    query <- standardize_subset(locus=locus,
                                top_SNPs=top_SNPs,
                                fullSS_genome_build=fullSS_genome_build,
                                subset_path=subset_path,
                                chrom_col=chrom_col,
                                position_col=position_col,
                                snp_col=snp_col,
                                pval_col=pval_col,
                                effect_col=effect_col,
                                stderr_col=stderr_col,
                                tstat_col=tstat_col,
                                MAF_col=MAF_col,
                                freq_col=freq_col,
                                A1_col=A1_col,
                                A2_col=A2_col,
                                gene_col=gene_col,

                                N_cases_col=N_cases_col,
                                N_controls_col=N_controls_col,
                                N_cases=N_cases,
                                N_controls=N_controls,
                                proportion_cases=proportion_cases,
                                sample_size=sample_size,
                                QTL_prefixes=QTL_prefixes)

    end_query <- Sys.time()
    printer("+ Extraction completed in", round(end_query-start_query, 2),"seconds", v=verbose)
    printer("+", dim(query)[1], "SNPs x ",dim(query)[2],"columns", v=verbose)
  }
  if(remove_tmps){
    printer("+ Removing subset tmp...", v=verbose)
    suppressWarnings(file.remove(subset_path))
  }
  return(query)
}



#' Use \emph{awk} to query locus subsets.
#'
#' Search full summary stats file by genomic coordinates.
#'
#' @family query functions
#' @inheritParams finemap_pipeline
#' @keywords internal
query_by_coordinates <- function(top_SNPs,
                                 locus,
                                 subset_path,
                                 fullSS_path,
                                 file_sep,
                                 chrom_col,
                                 position_col,
                                 min_POS,
                                 max_POS,
                                 bp_distance){
  gz.reader <- ifelse(endsWith(fullSS_path,".gz"), " zless "," cat ")
  topSNP_sub <- top_SNPs[top_SNPs$Locus==locus & !is.na(top_SNPs$Locus),]
  if(detect_genes(loci = locus, verbose = F)){
    topSNP_sub <- subset(topSNP_sub, Gene==names(locus))
  }
  if(is.na(min_POS)){min_POS <- topSNP_sub$POS - bp_distance}
  if(is.na(max_POS)){max_POS <- topSNP_sub$POS + bp_distance}
  printer("---Min snp position:",min_POS, "---")
  printer("---Max snp position:",max_POS, "---")
  colDict <- column_dictionary(fullSS_path)
  awk_cmd <- paste0(gz.reader,fullSS_path," | awk -F '",file_sep,"' 'NR==1 {print $0} NR>1 { if($",colDict[[chrom_col]]," == ",topSNP_sub$CHR,
                   " && ($", colDict[[position_col]]," >= ",min_POS," && $",colDict[[position_col]]," <= ",max_POS,")) {print $0} }'"," > ",subset_path)
  printer(awk_cmd)
  system(awk_cmd)
}




#' Use \emph{awk} to query locus subsets.
#'
#' Search full summary stats file by genomics coordinates.
#' To be used when \strong{CHR} and \strong{POS} are merged into one column (e.g. chr12:12209944).
#'
#' @family query functions
#' @param location_sep The separator character when \strong{CHR} and \strong{POS} are merged into one column (e.g. ":" when formatted like chr12:12209944).
#' @inheritParams finemap_pipeline
#' @keywords internal
query_by_coordinates_merged <- function(top_SNPs,
                                        fullSS_path,
                                        subset_path,
                                        locus,
                                        chrom_col,
                                        file_sep=" ",
                                        location_sep=":",
                                        min_POS,
                                        max_POS,
                                        bp_distance){
  topSNP_sub <- top_SNPs[top_SNPs$Locus==locus & !is.na(top_SNPs$Locus),][1,]
  if(is.na(min_POS)){min_POS <- topSNP_sub$POS - bp_distance}
  if(is.na(max_POS)){max_POS <- topSNP_sub$POS + bp_distance}
  printer("---Min snp position:",min_POS, "---")
  printer("---Max snp position:",max_POS, "---")
  colDict <- column_dictionary(fullSS_path)
  awk_cmd <- paste("cat ",fullSS_path," | tr -s '",location_sep,"' '",file_sep,"'",
                   " | awk -F '",file_sep,"' 'NR==1 {print \"CHR POS \" $2\" \" $3\" \" $4\" \" $5\" \" $6 }",
                   " NR>1 {if($",colDict[[chrom_col]]," == ",topSNP_sub$CHR," && ",
                   "($",colDict[[chrom_col]]+1," >=",min_POS,"&& $",colDict[[chrom_col]]+1," <=",max_POS,")) {print $0}}'",
                   " | tr -s '",file_sep,"' '\t' > ",subset_path, sep="")
  # awk_cmd <- paste("cat ",fullSS_path," | tr -s '",location_sep,"' '",file_sep,"'",
  #                  " | awk -F '",file_sep,"' 'NR==1 {print \"CHR POS \" $2\" \" $3\" \" $4\" \" $5\" \" $6 }",
  #                  " NR>1 {if($1 == 1 && ($2 >=",min_POS,"&& $2 <=",max_POS,")) {print $0}}'",
  #                  " | tr -s '",location_sep,"' '\t' > ",subset_path, sep="")
  printer(awk_cmd)
  system(awk_cmd)
}




#' Use \emph{awk} to query locus subsets.
#'
#' To be used when gene name is a column in the full summary stats file.
#' More commonly useful for QTL full summary stats files.
#'
#' @family query functions
#' @param gene Gene symbol (e.g. FOXP2) of the gene you want to query.
#' @param gene_col Name of the column that contains the gene name.
#' @param subset_path Path of the resulting locus subset file.
#' @inheritParams finemap_pipeline
#' @keywords internal
query_by_gene <- function(fullSS_path,
                           subset_path,
                           gene,
                           gene_col,
                           file_sep){
  colDict <- column_dictionary(fullSS_path)
  # if(endsWith(fullSS_path,".gz")){fullSS_path <- paste("<(gzcat",fullSS_path,")")}
  awk_cmd <- paste("awk -F '",file_sep,"' 'NR==1{print $0} NR>1{if($",colDict[[gene_col]]," == \"",gene,"\"){print $0}}' ",fullSS_path,
                   "| tr -s '",file_sep,"' '\t'  > ",subset_path, sep="")
  printer(awk_cmd)
  system(awk_cmd)
}

#' Use \emph{awk} to query locus subsets.
#'
#' To be used when probe name (but not gene name) is a column in the full summary stats file.
#' Requires a gene-probe key file (\strong{probe_path}).
#' More commonly useful for QTL full summary stats files.
#'
#' @family query functions
#' @param coordinates_merged Whether \strong{CHR} and \strong{POS} are merged into one column (e.g. chr12:12209944).
#' @param location_sep The separator character when \strong{CHR} and \strong{POS} are merged into one column (e.g. ":" when formatted like chr12:12209944).
#' @inheritParams finemap_pipeline
#' @keywords internal
query_by_probe <- function(fullSS_path,
                           subset_path,
                           gene,
                           gene_col,
                           chrom_col,
                           file_sep,
                           probe_path,
                           coordinates_merged=T,
                           location_sep=":"){
  ## EXAMPLE COMMAND LINE
  # awk -F ' ' 'NR==1{print $0} NR>1{if($2 == "ILMN_2226015" || $2 == "ILMN_1776649") {print $0}}' cis.eqtls.fairfax.all.chr.IFN.47231.367.b.qced.f.txt > LRRK2_Fairfax_IFN.txt
  # probe_path="Data/eQTL/Fairfax/locus.ILMN.map"
  # subset_path="Data/eQTL/Fairfax/LRRK2_Fairfax_CD14.txt"
  # file_sep="\t"
  # locus="LRRK2"
  find_probes <- function(map_file, genes){
    df  = data.table::fread(map_file)
    return(subset(df, GENE %in% genes)$PROBE_ID)
  }
  colDict <- column_dictionary(fullSS_path)
  probes <- find_probes(map_file = probe_path, genes = gene)
  probe_string <- paste(paste(paste("$",colDict[[gene_col]], sep="")," == \"" ,probes,"\"", sep=""), collapse=" || ")

  awk_cmd <- paste("awk -F '",file_sep,"' 'NR==1 {print $0} NR>1 if(",probe_string,") {print}' ",fullSS_path,
                   " > ",subset_path, sep="")
  printer(awk_cmd)
  system(awk_cmd)

  if(coordinates_merged){
    ##  EXAMPLE COMMAND LINE

    # awk -F ' ' 'NR==1{print "Coord","CHR","POS",$2,$3,$4,$5,$6 } NR>1{split($1,a,":"); print $1, a[1], a[2], $2, $3, $4, $5, $6}' LRRK2_Fairfax_CD14.txt | tr -s " " "\t" > tmp.txt && mv tmp.txt LRRK2_Fairfax_CD14.txt
    awk_cmd <- paste("awk -F '",file_sep,"' 'NR==1{print \"Coord\",\"CHR\",\"POS\",$2,$3,$4,$5,$6 }",
                     "NR>1{split($",colDict[[chrom_col]],",a,\":\"); print $1, a[1], a[2], $2, $3, $4, $5, $6}' ",
                     subset_path, " | tr -s ' ' '\t' > tmp.txt && mv tmp.txt ",subset_path, sep="")
    printer(awk_cmd)
    system(awk_cmd)
  }

}

#' Munge the full sumary stats file.
#'
#' Read in a standardize the entire full summary stats file at once.
#'
#' @family query functions
#' @inheritParams finemap_pipeline
#' @keywords internal
query_fullSS <- function(fullSS_path,
                         subset_path){
  file.copy(fullSS_path, subset_path)
}

#' Handles which query method to use
#'
#' @family query functions
#' @param subset_path Path of the resulting locus subset file.
#' @inheritParams finemap_pipeline
#' @inheritParams finemap_loci
#' @keywords internal
query_handler <- function(fullSS_path,
                          locus_dir=NULL,
                          top_SNPs=NULL,
                          subset_path,
                          min_POS=NA,
                          max_POS=NA,
                          bp_distance=500000,
                          locus_col="Gene",
                          chrom_col="CHR",
                          chrom_type=NULL,
                          position_col="POS",
                          file_sep="\t",
                          query_by="coordinates",
                          probe_path = "./Data/eQTL/gene.ILMN.map",
                          conda_env="echoR",
                          verbose=T){
  printer("+ Query Method:",query_by,v=verbose)
  locus <- basename(locus_dir)

  if(query_by=="tabix"){
    topSNP_sub <- top_SNPs[top_SNPs$Locus==locus & !is.na(top_SNPs$Locus),]
    if(detect_genes(loci = locus, verbose = F)){
      topSNP_sub <- subset(topSNP_sub, Gene==names(locus))
    }
    if(is.na(min_POS)){min_POS <- topSNP_sub$POS - bp_distance}
    if(is.na(max_POS)){max_POS <- topSNP_sub$POS + bp_distance}
    printer("+ QUERY: Chromosome =",topSNP_sub$CHR[1],
            "; Min position =", min_POS,
            "; Max position =", max_POS, v=verbose)
    study_dir <- get_study_dir(locus_dir)

    query <- TABIX(fullSS_path=fullSS_path,
                   subset_path=subset_path,
                   study_dir=study_dir,
                   chrom_col=chrom_col,
                   chrom_type=chrom_type,
                   position_col=position_col,
                   min_POS=min_POS,
                   max_POS=max_POS,
                   chrom=topSNP_sub$CHR[1],
                   conda_env=conda_env,
                   verbose=verbose
          )
  }
  if(query_by=="coordinates"){
    query_by_coordinates(top_SNPs=top_SNPs, locus=locus,
                         subset_path=subset_path, fullSS_path=fullSS_path,
                         chrom_col=chrom_col, position_col=position_col,
                         min_POS=min_POS, max_POS=max_POS, bp_distance=bp_distance,
                         file_sep=file_sep)
  }
  if(query_by=="coordinates_merged"){
    query_by_coordinates_merged(top_SNPs=top_SNPs, fullSS_path=fullSS_path, subset_path=subset_path,
                                file_sep=file_sep,  locus=locus, location_sep= ":",
                                min_POS=min_POS, max_POS=max_POS, bp_distance=bp_distance,
                                chrom_col=chrom_col)
  }
  if(query_by=="locus"){
    query_by_locus(fullSS_path=fullSS_path, subset_path=subset_path,
                  locus=locus, locus_col=locus_col, file_sep=file_sep)
  }
  if(query_by=="probes"){
    query_by_probe(fullSS_path=fullSS_path, subset_path=subset_path,
                   gene=gene, gene_col=gene_col, chrom_col=chrom_col,
                   file_sep=file_sep, probe_path=probe_path)
  }
  if(query_by=="fullSS"){
    query_fullSS(fullSS_path=fullSS_path,
                 subset_path = subset_path)
  }
}







