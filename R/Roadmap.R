
# ----------- ROADMAP ------------- #


#' Gather Roadmap annotation metadata
#'
#' @param keyword_query Search all columns in the Roadmap annotations metadata
#' and only query annotations that contain your keywords.
#' Can provide multiple keywords in list form:
#' \code{c("placenta","liver","monocytes")}
#' @family ROADMAP
ROADMAP.construct_reference <- function(ref_path = system.file("extdata/ROADMAP","ROADMAP_Epigenomic.js", package = "echolocatoR"),
                                        keyword_query=NULL){
  # %like% is from data.table
  ref <- suppressWarnings(data.table::fread(ref_path))
  colnames(ref)[1] <- "EID"
  if(!is.null(keyword_query)){
    rows <- grep(paste(keyword_query,collapse = "|"), data.table::transpose(ref), ignore.case = T)
    ref <- ref[rows,]
    printer("+ ROADMAP::",nrow(ref),"annotation(s) identified that match `keyword_query`.")
  }
  return(ref)
}



#' Query Roadmap API
#'
#' Query Roadmap epigenomic annotations (chromatin marks)
#' using a range of genomic coordinates.
#'
#' @param results_path Where to store query results.
#' @param chrom Chromosome to query
#' @param min_pos Minimum genomic position
#' @param max_pos Maximum genomic position
#' @param eid Roadmap annotation ID
#' @param convert_to_GRanges Whether to return query
#' as a \code{data.frame} or \code{\link[GenomicRanges]{GRanges}}.
#' @family ROADMAP
ROADMAP.tabix <- function(results_path,
                          chrom,
                          min_pos,
                          max_pos,
                          eid,
                          convert_to_GRanges=T){
  dir.create(results_path, showWarnings = F, recursive = T)
  chrom <- paste0("chr",gsub("chr","",base::tolower(chrom)))
  tbx_start = Sys.time()
  printer("++ Downloading Roadmap Chromatin Marks:",eid)
  fname <- paste0(eid,"_15_coreMarks_dense.bed.bgz")
  URL <- file.path("https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final",
                   fname) # _15_coreMarks_stateno.bed.gz
  cmd <- paste0("cd ",results_path," && tabix -p bed --begin 2 --end 3 ", URL," ",
                chrom,":",min_pos,"-",max_pos)
  # out <- system(cmd, intern = T)
  dat <- data.table::fread(cmd = cmd, select = 1:4, col.names = c("Chrom","Start","End","State"))
  dat$EID <- eid
  dat$File <- fname
  if(convert_to_GRanges){
    dat <- biovizBase::transformDfToGr(dat, seqnames = "Chrom", start = "Start", end="End")
  }
  tbx_end =  Sys.time()
  printer("BED subset downloaded in",round(tbx_end-tbx_start,3),"seconds")
  return(dat)
}





#' Query Roadmap by genomic coordinates
#'
#' @param gr.snp \code{\link[GenomicRanges]{GRanges}} object of SNPs to query Roadmap with.
#' @param limit_files Limit the number of annotation files queried (for faster testing).
#' @inheritParams  ROADMAP.tabix
#' @family ROADMAP
#' @examples
#' data("finemap_DT")
#' gr.snp <- DT_to_GRanges(subset_DT = finemap_DT)
#' grl.roadmap <- ROADMAP.query(results_path="./Roadmap", gr.snp=gr.snp, keyword_query="placenta")
ROADMAP.query <- function(results_path,
                          gr.snp,
                          keyword_query=NULL,
                          limit_files=NULL){
  rm_start = Sys.time()
  roadmap_ref <- ROADMAP.construct_reference(keyword_query=keyword_query)
  if(!is.null(limit_files)){
    roadmap_ref <- roadmap_ref[1:limit_files,]
  }
  # Download via tabix (fast)
  counter <- 1
  gr.roadmap <- lapply(unique(roadmap_ref$EID), function(eid,
                                                         gr.snp.=gr.snp,
                                                         results_path.=results_path){
    printer("+ Querying subset from Roadmap API:",
            eid," - ",counter,"/",length(unique(roadmap_ref$EID)))
    counter <<- counter+1
    dat <- GenomicRanges::GRanges()
    try({
      dat <- ROADMAP.tabix(results_path=results_path.,
                           chrom = gsub("chr","",GenomicRanges::seqnames(gr.snp.)[1]),
                           min_pos = min(gr.snp.$POS),
                           max_pos = max(gr.snp.$POS),
                           eid=eid,
                           convert_to_GRanges=T)
    })
    if(length(GenomicRanges::seqnames(dat))>0){
      return(dat)
    } else{return(NULL)}
  })
  remove(counter)
  grl.roadmap <- GR.name_filter_convert(GR.final = gr.roadmap,
                                        GR.names =  roadmap_ref$`Epigenome name (from EDACC Release 9 directory)`,
                                        min_hits=1)
  rm_end = Sys.time()
  printer("All downloads complete in",round(rm_end-rm_start,1),"minutes")
  return(grl.roadmap)
}


#' Standardize Roadmap query
#'
#' @param grl.roadmap Roadmap query results
#' @param n_top_tissues The number of top tissues to include,
#' sorted by greatest number of rows
#' (i.e. the number of genomic ranges within the window).
#' @family ROADMAP
ROADMAP.merge_and_process_grl <- function(grl.roadmap,
                                          gr.snp,
                                          n_top_tissues=5){
  grl.roadmap.merged <- unlist(grl.roadmap)
  grl.roadmap.merged$Source <- names(grl.roadmap.merged)
  grl.roadmap.merged$Source <- gsub("_"," ", grl.roadmap.merged$Source)
  grl.roadmap.merged$ChromState <- lapply(grl.roadmap.merged$State, function(ROW){base::strsplit(ROW, "_")[[1]][2]})%>% unlist()
  # Tally ChromStates
  # chromState_key <- data.table::fread(file.path("./echolocatoR/tools/Annotations/ROADMAP/ROADMAP_chromatinState_HMM.tsv"))
  # snp.pos <- subset(gr.snp, SNP %in% c("rs7294619"))$POS
  # snp.sub <- subset(grl.roadmap.merged, Start<=snp.pos & End>=snp.pos) %>%  data.frame()
  # chrom_tally <- snp.sub %>%
  #   dplyr::group_by(ChromState) %>%
  #   tally(sort = T) %>%
  #   merge(y=chromState_key[,c("MNEMONIC","DESCRIPTION")],
  #         by.x="ChromState", by.y="MNEMONIC", all.x=T, sort=F) %>% arrange(desc(n))
  # createDT(chrom_tally)

  grl.roadmap.filt <- grl.roadmap.merged[unlist( lapply(grl.roadmap, function(e){IRanges::overlapsAny(e, gr.snp, minoverlap = 1)}) )]
  if(!is.null(n_top_tissues)){
    top_tissues <-  data.frame(grl.roadmap.filt) %>% dplyr::group_by(Source) %>% dplyr::tally(sort = T)
    grl.roadmap.filt <- subset(grl.roadmap.filt, Source %in% unique(top_tissues$Source[1:n_top_tissues]))
  }
  return(grl.roadmap.filt)
}




#' Plot Roadmap query
#'
#' @param grl.roadmap.filt Roadmap query results.
#' @param gr.snp Optionally, can include an extra \code{\link[GenomicRanges]{GRanges}} object
#'  to ensure the plot does not extend beyond certain coordinates.
#' @param geom The type of plot to create.
#' Options include "density" and "histogram".
#' @param adjust The granularity of the peaks.
#' @param show_plot Whether to print the plot.
#' @examples
#' data("finemap_DT")
#' gr.snp <- DT_to_GRanges(finemap_DT)
#' grl.roadmap <- ROADMAP.query(results_path="./Roadmap", gr.snp=gr.snp, keyword_query="monocyte")
#' grl.roadmap.filt <- ROADMAP.merge_and_process_grl(grl.roadmap=grl.roadmap, gr.snp=gr.snp, n_top_tissues=5)
#' track.roadmap <- ROADMAP.track_plot(grl.roadmap.filt, gr.snp=gr.snp)
ROADMAP.track_plot <- function(grl.roadmap.filt,
                               gr.snp=NULL,
                               geom="density",
                               adjust=.2,
                               show_plot=T){
  track.roadmap <-  ggbio::autoplot(grl.roadmap.filt,
                                    which = gr.snp,
                            aes(fill=ChromState),
                            color="white",
                            size=.1,
                            geom = geom,
                            adjust = adjust,
                            # bins=10,
                            position="stack",# stack, fill, dodge
                            facets=Source~.,
                            alpha=1) +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 0),
          strip.text = element_text(size=9 )) +
    guides(fill = guide_legend(ncol = 2, keyheight = .5, keywidth = .5)) +
    scale_y_continuous(n.breaks = 3)
  if(show_plot){print(track.roadmap)}
  return(track.roadmap)
}




#' Query and plot Roadmap epigenomic annotations
#'
#' @param subset_DT Data.frame with at least the following columns:
#' \describe{
#' \item{SNP}{SNP RSID}
#' \item{CHR}{chromosome}
#' \item{POS}{position}
#' }
#' @param force_new_query Force a new query from the XGR database.
#' @inheritParams ROADMAP.construct_reference
#' @inheritParams ROADMAP.tabix
#' @inheritParams ROADMAP.merge_and_process_grl
#' @return A named list containing:
#' \itemize{
#' \item{Roadmap_plot}
#' \item{Roadmap_query}
#' }
#' @family ROADMAP
#' @return List containing:
#' \itemize{
#' \item{\code{ggbio} plot}
#' \item{\code{GRanges} object within the queried coordinates}
#' }
#' @examples
#' data("finemap_DT")
#' roadmap_plot_query <- ROADMAP.query_and_plot(subset_DT=finemap_DT, keyword_query="monocytes")
ROADMAP.query_and_plot <- function(subset_DT,
                                   results_path="./ROADMAP",
                                   n_top_tissues=NULL,
                                   keyword_query=NULL,
                                   adjust=.2,
                                   force_new_query=F,
                                   remove_tmps=T){
  # Convert subset to GRanges
  if(all(class(subset_DT)!="GRanges")){
    printer("ROADMAP:: Converting data to GRanges...")
    gr.snp <- biovizBase::transformDfToGr(data = dplyr::mutate(subset_DT, SEQnames = paste0("chr",CHR)),
                                          seqnames = "SEQnames",
                                          start = "POS",
                                          end = "POS")
  } else {
    gr.snp <- subset_DT
  }
  # Roadmap query
  lib <- "Roadmap_ChromatinMarks_CellTypes"
  anno_path <- file.path(results_path, "Annotation",paste0("GRanges_",lib,".rds"))
  if(file.exists(anno_path) & force_new_query==F){
    printer("+ Saved annotation file detected. Loading...")
    grl.roadmap <- readRDS(anno_path)
  } else {
    dir.create(dirname(anno_path), showWarnings = F, recursive = T)
    grl.roadmap <- ROADMAP.query(results_path = results_path,
                                 gr.snp = gr.snp,
                                 keyword_query = keyword_query,
                                 limit_files=NULL)
    save_annotations(gr = grl.roadmap,
                     anno_path = anno_path,
                     libName = lib)
  }
  grl.roadmap.filt <- ROADMAP.merge_and_process_grl(grl.roadmap = grl.roadmap,
                                                    gr.snp = gr.snp,
                                                    n_top_tissues=n_top_tissues)
  # Plot
  track.roadmap <- ROADMAP.track_plot(grl.roadmap.filt=grl.roadmap.filt,
                                      gr.snp=gr.snp,
                                      adjust=adjust)
  if(remove_tmps){
    tbi <- list.files(path = results_path, pattern = ".tbi$", full.names = T)
    dummy <- suppressWarnings(file.remove(tbi))
  }
  return(list(Roadmap_plot=track.roadmap,
              Roadmap_query=grl.roadmap.filt))
}


