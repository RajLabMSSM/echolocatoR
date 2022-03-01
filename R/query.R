#' Query handler
#'
#' Handles which query method to use
#'
#' @family query functions
#' @param subset_path Path of the resulting locus subset file.
#' @inheritParams finemap_locus
#' @inheritParams finemap_loci
#'
#' @keywords internal
#' @importFrom echotabix convert_and_query
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
                          conda_env="echoR",
                          verbose=TRUE){
  messager("+ Query Method:",query_by,v=verbose)
  locus <- basename(locus_dir)

  if(query_by=="tabix"){
    topSNP_sub <- top_SNPs[top_SNPs$Locus==locus & !is.na(top_SNPs$Locus),]
    if(detect_genes(loci = locus, verbose  = FALSE)){
      topSNP_sub <- subset(topSNP_sub, Gene==names(locus))
    }
    if(is.na(min_POS)){min_POS <- topSNP_sub$POS - bp_distance}
    if(is.na(max_POS)){max_POS <- topSNP_sub$POS + bp_distance}
    messager("+ QUERY: Chromosome =",topSNP_sub$CHR[1],
            "; Min position =", min_POS,
            "; Max position =", max_POS, v=verbose)
    study_dir <- get_study_dir(locus_dir)

    query <- echotabix::convert_and_query(
       fullSS_path=fullSS_path,
       subset_path=subset_path,
       study_dir=study_dir,
       chrom_col=chrom_col,
       start_col=position_col,
       min_POS=min_POS,
       max_POS=max_POS,
       chrom=topSNP_sub$CHR[1],
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
                   file_sep=file_sep)
  }
  if(query_by=="fullSS"){
    query_fullSS(fullSS_path=fullSS_path,
                 subset_path = subset_path)
  }
}


#' Munge the full sumary stats file.
#'
#' Read in a standardize the entire full summary stats file at once.
#'
#' @family query functions
#' @inheritParams finemap_locus
#' @keywords internal
query_fullSS <- function(fullSS_path,
                         subset_path){
  file.copy(fullSS_path, subset_path)
}


#' Use \emph{awk} to query locus subsets.
#'
#' Search full summary stats file by genomic coordinates.
#'
#' @family query functions
#' @inheritParams finemap_locus
#' @keywords internal
#' @importFrom echodata column_dictionary
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
  if(detect_genes(loci = locus, verbose  = FALSE)){
    topSNP_sub <- subset(topSNP_sub, Gene==names(locus))
  }
  if(is.na(min_POS)){min_POS <- topSNP_sub$POS - bp_distance}
  if(is.na(max_POS)){max_POS <- topSNP_sub$POS + bp_distance}
  messager("---Min snp position:",min_POS, "---")
  messager("---Max snp position:",max_POS, "---")
  colDict <- echodata::column_dictionary(fullSS_path)
  awk_cmd <- paste0(gz.reader,fullSS_path," | awk -F '",file_sep,"' 'NR==1 {print $0} NR>1 { if($",colDict[[chrom_col]]," == ",topSNP_sub$CHR,
                   " && ($", colDict[[position_col]]," >= ",min_POS," && $",colDict[[position_col]]," <= ",max_POS,")) {print $0} }'"," > ",subset_path)
  messager(awk_cmd)
  system(awk_cmd)
}




#' Use \emph{awk} to query locus subsets.
#'
#' Search full summary stats file by genomics coordinates.
#' To be used when \strong{CHR} and \strong{POS} are merged into one column (e.g. chr12:12209944).
#'
#' @family query functions
#' @param location_sep The separator character when \strong{CHR} and \strong{POS} are merged into one column (e.g. ":" when formatted like chr12:12209944).
#' @inheritParams finemap_locus
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
  messager("---Min snp position:",min_POS, "---")
  messager("---Max snp position:",max_POS, "---")
  colDict <- echodata::column_dictionary(fullSS_path)
  awk_cmd <- paste("cat ",fullSS_path," | tr -s '",location_sep,"' '",file_sep,"'",
                   " | awk -F '",file_sep,"' 'NR==1 {print \"CHR POS \" $2\" \" $3\" \" $4\" \" $5\" \" $6 }",
                   " NR>1 {if($",colDict[[chrom_col]]," == ",topSNP_sub$CHR," && ",
                   "($",colDict[[chrom_col]]+1," >=",min_POS,"&& $",colDict[[chrom_col]]+1," <=",max_POS,")) {print $0}}'",
                   " | tr -s '",file_sep,"' '\t' > ",subset_path, sep="")
  # awk_cmd <- paste("cat ",fullSS_path," | tr -s '",location_sep,"' '",file_sep,"'",
  #                  " | awk -F '",file_sep,"' 'NR==1 {print \"CHR POS \" $2\" \" $3\" \" $4\" \" $5\" \" $6 }",
  #                  " NR>1 {if($1 == 1 && ($2 >=",min_POS,"&& $2 <=",max_POS,")) {print $0}}'",
  #                  " | tr -s '",location_sep,"' '\t' > ",subset_path, sep="")
  messager(awk_cmd)
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
#' @inheritParams finemap_locus
#' @keywords internal
query_by_gene <- function(fullSS_path,
                           subset_path,
                           gene,
                           gene_col,
                           file_sep){
  colDict <- echodata::column_dictionary(fullSS_path)
  # if(endsWith(fullSS_path,".gz")){fullSS_path <- paste("<(gzcat",fullSS_path,")")}
  awk_cmd <- paste("awk -F '",file_sep,"' 'NR==1{print $0} NR>1{if($",colDict[[gene_col]]," == \"",gene,"\"){print $0}}' ",fullSS_path,
                   "| tr -s '",file_sep,"' '\t'  > ",subset_path, sep="")
  messager(awk_cmd)
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
#' @inheritParams finemap_locus
#' @keywords internal
query_by_probe <- function(fullSS_path,
                           subset_path,
                           gene,
                           gene_col,
                           chrom_col,
                           file_sep,
                           probe_path=NULL,
                           coordinates_merged=TRUE,
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
  colDict <- echodata::column_dictionary(fullSS_path)
  probes <- find_probes(map_file = probe_path, genes = gene)
  probe_string <- paste(paste(paste("$",colDict[[gene_col]], sep="")," == \"" ,probes,"\"", sep=""), collapse=" || ")

  awk_cmd <- paste("awk -F '",file_sep,"' 'NR==1 {print $0} NR>1 if(",probe_string,") {print}' ",fullSS_path,
                   " > ",subset_path, sep="")
  messager(awk_cmd)
  system(awk_cmd)

  if(coordinates_merged){
    ##  EXAMPLE COMMAND LINE

    # awk -F ' ' 'NR==1{print "Coord","CHR","POS",$2,$3,$4,$5,$6 } NR>1{split($1,a,":"); print $1, a[1], a[2], $2, $3, $4, $5, $6}' LRRK2_Fairfax_CD14.txt | tr -s " " "\t" > tmp.txt && mv tmp.txt LRRK2_Fairfax_CD14.txt
    awk_cmd <- paste("awk -F '",file_sep,"' 'NR==1{print \"Coord\",\"CHR\",\"POS\",$2,$3,$4,$5,$6 }",
                     "NR>1{split($",colDict[[chrom_col]],",a,\":\"); print $1, a[1], a[2], $2, $3, $4, $5, $6}' ",
                     subset_path, " | tr -s ' ' '\t' > tmp.txt && mv tmp.txt ",subset_path, sep="")
    messager(awk_cmd)
    system(awk_cmd)
  }

}
