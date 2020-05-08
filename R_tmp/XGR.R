
# DeepBlueR: Alternate source of annotation files
## https://bioconductor.org/packages/release/bioc/vignettes/DeepBlueR/inst/doc/DeepBlueR.html#listing-experiments


# TUTORIALS
# http://xgr.r-forge.r-project.org/#xenrichergenes
# http://galahad.well.ox.ac.uk:3020/R/ds2




## Download BED files via XGR
GRs.to.BED <- function(GR.annotations, output_path, sep="\t"){
  BED_paths <- lapply(names(GR.annotations), function(name){ 
    GR <- GR.annotations[[name]] 
    BED <- GR %>% as.data.table() %>% 
      dplyr::select(chrom=seqnames,
                    chromStart=start,
                    chromEnd=end,
                    strand) 
    BED_path <- file.path(output_path,paste0(gsub(":","-",name),".bed.txt"))
    dir.create(dirname(BED_path), recursive = T, showWarnings = F)
    data.table::fwrite(BED, BED_path, sep=sep, col.names = F, quote = F) 
    return(BED_path)
  }) %>% unlist()
  return(BED_paths)
}

XGR.prepare_foreground_background  <- function(subset_DT,
                                               foreground_filter="Support>0",
                                               background_filter=NULL){
  # Foreground
  fg <- subset(subset_DT, eval(parse(text=foreground_filter))) %>% 
    dplyr::mutate(chrom = paste0("chr",CHR), 
                  chromStart = POS, 
                  chromEnd = POS, 
                  name = SNP) %>% 
    dplyr::select(chrom, chromStart, chromEnd, name)
  # Background
  if(!is.null(background_filter)){
    bg_DT <- subset(subset_DT, eval(parse(text=background_filter)))
  } else { bg_DT <- subset_DT  }
  
  bg <- bg_DT %>% dplyr::mutate(chrom = paste0("chr",CHR), 
                                chromStart = POS, 
                                chromEnd = POS, 
                                name = SNP) %>% 
    dplyr::select(chrom, chromStart, chromEnd, name)
  printer("XGR::",nrow(fg),"SNPs in foreground.")
  printer("XGR::",nrow(bg),"SNPs in background")
  return(list("foreground"=fg,
              "background"=bg))
}

XGR.iterate_overlap <- function(lib.selections=c("ENCODE_TFBS_ClusteredV3_CellTypes",
                                                 "TFBS_Conserved",
                                                 "ReMap_PublicAndEncode_TFBS",
                                                 "Uniform_TFBS"),
                                subset_DT){
  OVERLAP <- lapply(lib.selections, function(lib.name){
    GR.annotations <- xRDataLoader(RData.customised = lib.name) 
    all_GRL <- parallel::mclapply(names(GR.annotations), function(n1){
      grl <- GR.annotations[[n1]]
      # Handle both nested and unnested entries
      if(class(grl)=="list"){
          GRL <- lapply(names(grl), function(n2){
            gr <- grl[[n2]]  
            gr$source <- n1
            gr$assay <- n2
            return(gr)
          }) %>% unlist()   
      } else {return(grl)} 
    }, mc.cores = 4) %>% unlist() # return all_GRL
    if(!is.null(all_GRL)){
      ALL_GRL <- unlist(GenomicRanges::GRangesList(all_GRL)) 
      ALL_GRL$library <- lib.name 
      return(ALL_GRL)
    } else{return(NULL)} 
  }) # return OVERLAP
  ALL_OVERLAP <- unlist(GenomicRanges::GRangesList(unlist(OVERLAP))) 
  
  no_no_loci <- c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1",  
                  "RETREG3","UBTF","FAM171A2","MAP3K14","CRHR1","MAPT-AS1","KANSL1","NSF","WNT3")
  subset_DT <- subset(subset_DT, !Locus %in% no_no_loci)
  subset_DT <- assign_lead_SNP(subset_DT)
  gr.hits <- NOTT_2019.report_regulatory_overlap(finemap_DT = subset_DT,
                                                 regions = ALL_OVERLAP) 
  
  ucs.hits <-subset(gr.hits, Consensus_SNP) 
  length(unique(subset(subset_DT, Consensus_SNP)$SNP))
  length(unique(subset(subset_DT, Consensus_SNP)$Locus))
  
  length(unique(ucs.hits$SNP))
  length(unique(ucs.hits$Locus))
  
  if(save_path!=F){
    save_path  <-  "Data/GWAS/Nalls23andMe_2019/_genome_wide/XGR/XGR_TFBS_UCS_overlap.csv.gz"
    dir.create(dirname(save_path), showWarnings = F, recursive = T)
    data.table::fwrite(data.frame(ucs.hits),save_path)
  }
  return(gr.hits)
}



XGR.iterate_enrichment <- function(foreground_snps,
                                   background_snps=NULL,
                                   lib.selections = c("ENCODE_TFBS_ClusteredV3_CellTypes",
                                                      "ENCODE_DNaseI_ClusteredV3_CellTypes", 
                                                      "Broad_Histone",
                                                      "FANTOM5_Enhancer",
                                                      "Segment_Combined_Gm12878",
                                                      "TFBS_Conserved",
                                                      "ReMap_PublicAndEncode_TFBS",
                                                      "Blueprint_VenousBlood_Histone",
                                                      "Blueprint_DNaseI",
                                                      # "Blueprint_Methylation_hyper",
                                                      # "Blueprint_Methylation_hypo",
                                                      # "Genic_anno",
                                                      "FANTOM5_CAT_Cell",
                                                      "FANTOM5_CAT_MESH",
                                                      "GWAScatalog_alltraits"),  
                                   save_path=F){ 
  # Description of all datasets
  # https://www.rdocumentation.org/packages/XGR/versions/1.1.5/topics/xDefineGenomicAnno 
  fg_bg <- XGR.prepare_foreground_background(subset_DT, 
                                             foreground_filter = "Consensus_SNP",
                                             background_filter = "leadSNP")
  
  # lib.selections = c("ENCODE_TFBS_ClusteredV3_CellTypes",
  #                    "ENCODE_DNaseI_ClusteredV3_CellTypes", 
  #                    "Broad_Histone",
  #                    "UW_Histone",
  #                    "SYDH_Histone",
  #                    "FANTOM5_Enhancer", 
  #                    "TFBS_Conserved",
  #                    "Uniform_TFBS",
  #                    "Uniform_DNaseI_HS")
  # roadmap_grl <- lapply(unique(subset_DT$Locus), function(locus){
  #       locus_DT <- subset(subset_DT, Locus==locus) 
  #       dat <- ROADMAP.tabix(results_path=results_path,
  #                            chrom = locus_DT$CHR[1],
  #                            min_pos = min(locus_DT$POS),
  #                            max_pos = max(locus_DT$POS),
  #                            eid=eid,
  #                            convert_to_GRanges=T)
  #       return(dat)
  # })
  
  
 
  database_results <- parallel::mclapply(lib.selections, function(lib.name){  
    printer("XGR:: Testing enrichment: ",lib.name)
    eTerm <- NULL
    try({
      GR.annotations <- xRDataLoader(RData.customised = lib.name)
      eTerm <- lapply(GR.annotations, function(grl){
         et <-xGRviaGenomicAnno(data.file = fg_bg$foreground,
                                   background.file = fg_bg$background,
                                   format.file="data.frame",
                                   GR.annotation = grl)  
        return(et)
      }) %>% data.table::rbindlist() 
      eTerm$lib <- lib.name
      eTerm$fullname <- names(unlist(GR.annotations))
      eTerm$source <-  lapply(eTerm$fullname, function(e){strsplit(e, "[.]")[[1]][1]}) 
      eTerm$assay <-  lapply(eTerm$fullname, function(e){strsplit(e, "[.]")[[1]][2]})   
    }) 
    return(eTerm)
  }, mc.cores = 4) 
  
  # Re-calculate corrected p-val to account for multiple dbs tested
  DB_res <- data.table::rbindlist(database_results) %>% 
    dplyr::mutate(FDR = p.adjust(p=pvalue, method = "fdr"),
                  Bonf = p.adjust(p=pvalue, method = "bonferroni")) %>%
    dplyr::arrange(FDR, -nOverlap, -fc) %>%
    subset(adjp<0.05) %>%
    data.table::data.table()
  
  # DB_res.conensusVSall <- DB_res
  
  
  if(save_path!=F){
    data.table::fwrite(DB_res, save_path, quote = F)
  }
  
  
}

XGR.plot_enrichment <- function(DB_res){
  DB_res$annotation <- factor(gsub(".*[$]","",DB_res$name), levels = rev(unique(gsub(".*[$]","",DB_res$name))), ordered = T)
  
  ggplot(data = DB_res, aes(x=fc, y=annotation, fill=fc)) + 
    geom_col() +
    labs(title="Epigenomic annotation enrichment", 
         subtitle="Foreground = Consensus SNPs\nBackground = Lead GWAS SNPs") +
    facet_grid(facets = lib~., 
               scales = "free_y",
               space = "free_y", 
               switch = "y") + 
    theme_bw() +
    theme(strip.placement = "outside",
          strip.text.y.left  = element_text(angle = 0), 
          strip.background = element_rect(color="black",fill="white"),  
          plot.title = element_text(hjust = .5),
          plot.subtitle = element_text(hjust = .5))
}









# SNP-based enrichment analysis using GWAS Catalog traits (mapped to EF) 
# data <- merged_results$SNP %>% unique()
# merged_DT <- merge_finemapping_results(minimum_support=0,
#                                        include_leadSNPs=T,
#                                        dataset = "./Data/GWAS/Nalls23andMe_2019/",
#                                        xlsx_path=F,
#                                        from_storage=T,
#                                        consensus_thresh = 2,
#                                        haploreg_annotation=F,
#                                        biomart_annotation=F, 
#                                        verbose = F) 
# subset_DT <- merged_DT
# results_path <- "Data/GWAS/Nalls23andMe_2019/_genome_wide"
# 
# xgr_snpEnrich <- function(subset_DT){
#   library(XGR)  
#   snps <- subset(subset_DT, Support>0)$SNP
#   
#   eTerm <- xEnricherSNPs(data=snps, 
#                          ontology="EF",
#                          path.mode=c("all_paths"), 
#                          min.overlap = 1)
#   xEnrichViewer(eTerm)
#   # e) barplot of significant enrichment results
#   bp <- xEnrichBarplot(eTerm, top_num="auto", displayBy="adjp")
#   print(bp)
#   # f) visualise the top 10 significant terms in the ontology hierarchy
#   # color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
#   xEnrichDAGplot(eTerm, top_num=10, 
#                  displayBy="adjp",
#                  node.info=c("full_term_name"))
#   # color-code terms according to the z-scores
#   xEnrichDAGplot(eTerm, top_num=10, 
#                  displayBy="zscore",
#                  node.info=c("full_term_name"))
#   # Circos plot
#   library(RCircos)
#   SNP.g <- xSocialiserSNPs(data = snps, 
#                            include.LD=NA, 
#                            LD.r2 = 1) 
#   xCircos(g=SNP.g, 
#           entity="SNP")
#   # b') optionally, enrichment analysis for input SNPs plus their LD SNPs
#   ## LD based on European population (EUR) with r2>=0.8
#   # eTerm_LD <- xEnricherSNPs(data=data, 
#   #                        include.LD="EUR", 
#   #                        LD.r2=0.8,  
#   #                        min.overlap = 1)
#   # xEnrichViewer(eTerm_LD)
# }

# xgr_geneEnrich <- function(top_SNPs){
#   library(XGR) 
#   # GENE-LEVEL ENRICHMENT 
#   data <- top_SNPs$Gene %>% unique()
#   eTerm_gene <- xEnricherGenes(data = data, 
#                                ontology = "MsigdbC2REACTOME", 
#                                min.overlap = 1) 
#   bp <- xEnrichBarplot(eTerm_gene, top_num="auto", displayBy="adjp")
#   print(bp)
#   
#   xEnrichDAGplot(eTerm, top_num=10, displayBy="zscore",
#                  node.info=c("full_term_name"), graph.node.attrs=list(fontsize=200)) 
#   } 

# 
# xgr_annoEnrich <- function(subset_DT){
#   library(XGR) 
#   RData.location <- "http://galahad.well.ox.ac.uk/bigdata" 
#   CS <- subset(subset_DT, Support>0)
#   
#   ## a) provide input data
#   # data.file <- "http://galahad.well.ox.ac.uk/bigdata/ImmunoBase_GWAS.bed" 
#   data.file_full <- subset_DT %>% dplyr::mutate(chrom = paste0("chr",CHR), 
#                                          chromStart = POS, 
#                                          chromEnd = POS, 
#                                          name = SNP) %>% 
#     dplyr::select(chrom, chromStart, chromEnd, name)
#   data.file_CS <- CS %>% dplyr::mutate(chrom = paste0("chr",CHR), 
#                                             chromStart = POS, 
#                                             chromEnd = POS, 
#                                             name = SNP) %>% 
#     dplyr::select(chrom, chromStart, chromEnd, name)
#   
#   
#   
#   # Built-in database
#   ## b) perform enrichment analysis using FANTOM expressed enhancers
#   ### one-tail p-value calculation (by default) 
#   eTerm_CS <- xGRviaGenomicAnno(data.file = data.file_CS, 
#                                 background.file = data.file_full,
#                                 format.file = "data.frame",
#                                 GR.annotation = "FANTOM5_Enhancer_Cell")
#   xEnrichViewer(eTerm_CS, 10)
#   # Downloaded database
#   ## Download
#   gr_broad <- xRDataLoader(RData.customised="Broad_Histone")
#   # "ENCODE_TFBS_ClusteredV3_CellTypes"
#   # "ENCODE_DNaseI_ClusteredV3_CellTypes"
#   gr_reactome <- xRDataLoader(RData.customised ="org.Hs.egMsigdbC2REACTOME")
#   
#   # Run
#   eTerm_broad <- xGRviaGenomicAnno(data.file = data.file_CS, 
#                                      background.file = data.file_full,
#                                      format.file="data.frame", 
#                                      GR.annotation = gr_broad)
#   subset(xEnrichViewer(eTerm_broad), adjp<0.05) 
#   
#   
#   # Gtex
#   ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase')
#   grl <- xDefineGenomicAnno("Uniform_TFBS")
#   
#   gr_gtex <- xRDataLoader(RData.customised="org.Hs.egGTExV6")
#   eTerm_ImmunoBase <- xGRviaGenomicAnno(data.file = data.file_CS, 
#                                    background.file = data.file_full,
#                                    format.file="data.frame",  
#                                    GR.annotation = ImmunoBase)
#   subset(xEnrichViewer(eTerm_ImmunoBase), adjp<=0.05) 
# }  


# xgr_annoEnrich_mass <- function(foreground_snps,
#                                 background_snps=NULL,
#                                 lib.selections = c("ENCODE_TFBS_ClusteredV3_CellTypes",
#                                                    "ENCODE_DNaseI_ClusteredV3_CellTypes",
#                                                     # "Broad_Histone",
#                                                     "FANTOM5_Enhancer",
#                                                     "Segment_Combined_Gm12878",
#                                                     "TFBS_Conserved",
#                                                     "ReMap_PublicAndEncode_TFBS",
#                                                     "Blueprint_VenousBlood_Histone",
#                                                     "Blueprint_DNaseI",
#                                                     # "Blueprint_Methylation_hyper",
#                                                     # "Blueprint_Methylation_hypo",
#                                                     # "Genic_anno",
#                                                     "FANTOM5_CAT_Cell",
#                                                     "FANTOM5_CAT_MESH",
#                                                     "GWAScatalog_alltraits"),
#                                 save_path=F){
#   # Description of all datasets
#   # https://www.rdocumentation.org/packages/XGR/versions/1.1.5/topics/xDefineGenomicAnno
# 
#   data.file <- foreground_snps
#   data.file_full <- background_snps
# 
#   database_results <- lapply(lib.selections, function(lib.name){
#     tryCatch({
#       # lib.name <- "ENCODE_TFBS_ClusteredV3_CellTypes"
#       message("+ Testing enrichment: ",lib.name)
#       # Import annotation files from server
#       GR.annotations <- xRDataLoader(lib.name,
#                                      RData.location=RData.location)
#       ## Un-nest files list (e.g. some are organized per tissue per TF)
#       GR.annotations <- unlist(GR.annotations)
#       # Iterate over each file
#       ls_df <- lapply(1:length(GR.annotations), function(i,
#                                                          library = lib.name,
#                                                          data.file. = data.file){
#         GR.annotation <- GR.annotations[i]
#         message(sprintf("Analysing '%s' (%s) ...", names(GR.annotation),
#                         as.character(Sys.time())), appendLF=T)
#         df <- xGRviaGenomicAnno(data.file=data.file.,
#                                 background.file = data.file_full,
#                                 format.file="data.frame",
#                                 GR.annotation=GR.annotation,
#                                 RData.location=RData.location,
#                                 verbose=F)
#         df$library <- library
#         df$total_snps <- nrows(data.file.)
#         return(df)
#       })
#       df <- do.call(rbind, ls_df) %>% data.table::as.data.table()
#       return(df)
#     },
#     error=function(e){data.table::data.table(name=NA, nAnno=NA,
#                                              nOverlap=NA,fc=NA,
#                                              zscore=NA,pvalue=NA,
#                                              adjp=NA,or=NA,
#                                              CIl=NA,CIu=NA,
#                                              total_snps=NA,
#                                              library=lib.name)})
#   })
#   # Process results
#   DT <- database_results %>%
#     data.table::rbindlist()
#   # %>% subset(!is.na(name))
#   # bonf <- 0.05 / nrow(DT)
#   # subset(DT, adjp<=bonf & nOverlap > 1) %>% arrange(desc(nOverlap))
# 
#   if(save_path!=F){
#     data.table::fwrite(DT, save_path, quote = F)
#   }
#   # Heatmap
#   # gp <- xEnrichHeatmap(database_results, fdr.cutoff=0.05, displayBy="fdr",
#   #                      reorder="both")
#   # gp
#   # # Barplot
#   # bp <- xEnrichBarplot(eTerm, top_num='auto', displayBy="fc")
#   # bp
#   # # Forest plot
#   # gp <- xEnrichForest(eTerm)
#   # gp
#   return(DT)
# }

  
# gather_databases <- function(lib.selections){
#   databases <- lapply(lib.selections, function(lib.name){ 
#     tryCatch({
#       xRDataLoader(lib.name)
#     }, 
#     error=function(e) NULL) 
#   })
#   names(databases) <- paste0(lib.selections,"..")
#   dbs <- unlist(databases)
#   names(dbs)
#   return(dbs)
# }
# GR.annotations <- gather_databases(lib.selections)

XGR.merge_and_process <- function(grl.xgr,
                                  lib,
                                  n_top_sources=10){
  # grl.xgr <- check_saved_XGR(results_path, lib)
  ## Make track
  ## Add and modify columns
  grl.xgr.merged <- unlist(grl.xgr)
  names(grl.xgr.merged) <- gsub("Broad_Histone_","",names(grl.xgr.merged))
  sep <- list("ENCODE_TFBS_ClusteredV3_CellTypes"="[.]",
              "ENCODE_DNaseI_ClusteredV3_CellTypes"="_",
              "Broad_Histone"="_",
              "FANTOM5_Enhancer"="_")[[lib]] 
  grl.xgr.merged$Source <- lapply(names(grl.xgr.merged), function(e){strsplit(e, sep)[[1]][1]}) %>% unlist()
  # grl.xgr.merged$Source <- gsub("_","\n", grl.xgr.merged$Source)
  grl.xgr.merged$Assay <- lapply(names(grl.xgr.merged), function(e){strsplit(e, sep)[[1]][2]}) %>% unlist()
  grl.xgr.merged$Start <- GenomicRanges::start(grl.xgr.merged)
  grl.xgr.merged$End <- GenomicRanges::end(grl.xgr.merged)
  # Filter
  top_sources <- grl.xgr.merged %>% data.frame() %>% dplyr::group_by(Source) %>% dplyr::tally(sort = T)
  grl.xgr.merged.filt <- subset(grl.xgr.merged, Source %in% unique(top_sources$Source[1:n_top_sources]))
  # Count
  # snp.pos <- subset(gr.snp, SNP %in% c("rs7294619"))$POS
  # snp.sub <- subset(grl.xgr.merged, Start<=snp.pos & End>=snp.pos) %>% data.frame() 
  grl.xgr.merged.filt$Source_Assay <- paste0(grl.xgr.merged.filt$Source,"_",grl.xgr.merged.filt$Assay)
  return(grl.xgr.merged.filt)
}


XGR.import_annotations <- function(gr.snp,
                                    anno_data_path=file.path("echolocatoR/tools/Annotations", paste0("XGR_",lib.name,".rds")),
                                    lib.name,
                                    save_xgr=T,
                                    annot_overlap_threshold=5){
  library(GenomicRanges)
  if(file.exists(anno_data_path)){
    printer("")
    printer("+ Saved annotation file detected. Loading...")
    GR.annotations <- readRDS(anno_data_path)
  } else {
    printer("")
    printer("+ XGR: Downloading...",lib.name)
    GR.annotations <- XGR::xRDataLoader(RData.customised=lib.name)
    if(save_xgr & !is.null(GR.annotations)){
      dir.create( dirname(anno_data_path), showWarnings = F, recursive = T)
      saveRDS(GR.annotations, file = anno_data_path)
    }
  }
  GR.orig <- unlist(GR.annotations)
  
  gr.xgr <- lapply(names(GR.orig), function(g, gr.snp. = gr.snp){
    # printer("Finding overlap for:", g)
    GR.overlap <- subsetByOverlaps(GR.orig[[g]], gr.snp.)
    len <- length(seqnames(GR.overlap) )
    # printer("   - Overlapping annotations = ",len)
    if(len>0){
      return(GR.overlap)
    } else{return(NULL)}
  })
  grl.xgr <- GR.name_filter_convert(gr.xgr, names(GR.orig), min_hits=annot_overlap_threshold)
  return(grl.xgr)
}


