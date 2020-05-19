






ROADMAP.tabix <- function(results_path, 
                          chrom, 
                          min_pos, 
                          max_pos, 
                          eid, 
                          convert_to_GRanges=T){
  tbx_start = Sys.time()
  printer("++ Downloading Roadmap Chromatin Marks:",eid)
  fname <- paste0(eid,"_15_coreMarks_dense.bed.bgz")
  URL <- file.path("https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final",
                   fname) # _15_coreMarks_stateno.bed.gz
  cmd <- paste0("cd ",results_path,"&& tabix -p bed --begin 2 --end 3 ", URL," ",
                chrom,":",min_pos,"-",max_pos)
  # out <- system(cmd, intern = T)
  dat <- data.table::fread(cmd = cmd , sep = "\t", select = 1:4, col.names = c("Chrom","Start","End","State"))
  dat$EID <- eid
  dat$File <- fname
  if(convert_to_GRanges){
    dat <- biovizBase::transformDfToGr(dat, seqnames = "Chrom", start = "Start", end="End")
  }
  tbx_end =  Sys.time()
  printer("BED subset downloaded in",round(tbx_end-tbx_start,3),"seconds")
  return(dat)
}


ROADMAP.track <- function(results_path, gr.snp, limit_files=NA){
  rm_start = Sys.time()
  RoadMap_ref <- GS_construct_reference()
  if(!is.na(limit_files)){
    RoadMap_ref <- RoadMap_ref[1:limit_files,]
  }
  # Download via tabix (fast)
  counter <- 1
  gr.roadmap <- lapply(unique(RoadMap_ref$EID), function(eid,
                                                         gr.snp.=gr.snp,
                                                         results_path.=results_path){
    printer("+ Querying subset from Roadmap API:",
            eid," - ",counter,"/",length(unique(RoadMap_ref$EID)))
    counter <<- counter+1
    dat <- GRanges()
    try({
      dat <- ROADMAP.tabix(results_path=results_path,
                           chrom = unique(seqnames(gr.snp.)),
                           min_pos = min(gr.snp.$POS),
                           max_pos = max(gr.snp.$POS),
                           eid=eid,
                           convert_to_GRanges=T)
    })
    if(length(seqnames(dat))>0){
      return(dat)
    } else{return(NULL)}
  })
  remove(counter)
  grl.roadmap <- GR.name_filter_convert(gr.roadmap, RoadMap_ref$Epigenome.name, min_hits=1)
  rm_end = Sys.time()
  printer("All downloads complete in",round(rm_end-rm_start,1),"minutes")
  return(grl.roadmap)
}



ROADMAP.merge_and_process_grl <- function(grl.roadmap,
                                          n_top_tissues=5){
  library(IRanges)
  grl.roadmap.merged <- unlist(grl.roadmap)
  grl.roadmap.merged$Source <- names(grl.roadmap.merged)
  grl.roadmap.merged$Source <- gsub("_"," ", grl.roadmap.merged$Source)
  grl.roadmap.merged$ChromState <- lapply(grl.roadmap.merged$State, function(ROW){strsplit(ROW, "_")[[1]][2]})%>% unlist()
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
  
  grl.roadmap.filt <- grl.roadmap.merged[unlist( lapply(grl.roadmap, function(e){overlapsAny(e, gr.snp, minoverlap = 1)}) )]
  top_tissues <- grl.roadmap.filt %>% data.frame() %>% dplyr::group_by(Source) %>% dplyr::tally(sort = T)
  grl.roadmap.filt <- subset(grl.roadmap.filt, Source %in% unique(top_tissues$Source[1:n_top_tissues]))
  return(grl.roadmap.filt)
}
