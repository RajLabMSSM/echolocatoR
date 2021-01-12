# DeepBlueR: Alternate source of annotation files
## https://bioconductor.org/packages/release/bioc/vignettes/DeepBlueR/inst/doc/DeepBlueR.html#listing-experiments

# INSTALLATION
## CRAN version no longer works. Follow these instructions instead:
## https://github.com/hfang-bristol/XGR/issues/10
# > BiocManager::install("remotes")
# > BiocManager::install("hfang-bristol/XGR", dependencies=T)


# TUTORIALS
# http://xgr.r-forge.r-project.org/#xenrichergenes
# http://galahad.well.ox.ac.uk:3020/R/ds2



#' Name annotation file
#'
#' @keywords internal
annotation_file_name <- function(locus_dir,
                                 lib_name){
  annot_dir <- file.path(locus_dir,"annotations")
  dir.create(annot_dir, showWarnings = F, recursive = T)
  annot_file <- file.path(annot_dir, paste0(lib_name,".RDS"))
  return(annot_file)
}




#' Convert GRanges object to BED format and save
#'
#' @family XGR
#' @keywords internal
GRanges_to_BED <- function(GR.annotations,
                           output_path,
                           sep="\t",
                           nThread=4,
                           gzip=F){
  BED_paths <- parallel::mclapply(names(GR.annotations), function(name,
                                                                  .gzip=gzip){
    GR <- GR.annotations[[name]]
    BED <- data.table::as.data.table(GR) %>%
      dplyr::select(chrom=seqnames,
                    chromStart=start,
                    chromEnd=end,
                    strand)
    BED_path <- file.path(output_path,paste0(gsub(":","-",name),".bed.txt"))
    dir.create(dirname(BED_path), recursive = T, showWarnings = F)
    data.table::fwrite(BED, BED_path, sep=sep, col.names = F, quote = F)
    if(.gzip){
      R.utils::gzip(BED_path, overwrite=T)
      BED_path <- paste0(BED_path,".gz")
    }
    return(BED_path)
  }, mc.cores = nThread) %>% unlist()
  return(BED_paths)
}




#' Convert data.table to GRanges object
#'
#' @family XGR
#' @keywords internal
DT_to_GRanges <- function(subset_DT){
  dat <-  dplyr::mutate(subset_DT, SEQnames = paste0("chr",CHR))
  gr.snp <- biovizBase::transformDfToGr(dat,
                                        seqnames = "SEQnames",
                                        start = "POS",
                                        end = "POS")
  GenomeInfoDb::seqlevelsStyle(gr.snp) <- "NCBI"
  return(gr.snp)
}

#' Prepare SNP sets for enrichment
#'
#' Prepare custom foreground and background SNPs sets for enrichment tests with XGR annotations.
#'
#' @param subset_DT Data.frame with at least the following columns:
#' \describe{
#' \item{SNP}{SNP RSID}
#' \item{CHR}{chromosome}
#' \item{POS}{position}
#' }
#' @param foreground_filter Specify foreground by filtering SNPs in \code{subset_DT}.
#' Write filter as a string (or \code{NULL} to include all SNPs).
#' @param background_filter Specify background by filtering SNPs in \code{subset_DT}.
#' Write filter as a string (or \code{NULL} to include all SNPs).
#' @family XGR
#' @keywords internal
#' @examples
#' data("merged_DT")
#' merged_DT <- data.table::fread("/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/merged_UKB.csv.gz", nThread=4)
#' fg_bg <- XGR.prepare_foreground_background(subset_DT=merged_DT, foreground_filter="Consensus_SNP==T", background_filter="leadSNP==T")
XGR.prepare_foreground_background  <- function(subset_DT,
                                               foreground_filter="Support>0",
                                               background_filter=NULL,
                                               fg_sample_size=NULL,
                                               bg_sample_size=NULL,
                                               verbose=T){
  if(!exists("sampling_df")) sampling_df <- subset_DT
  printer("XGR:: Preparing foreground/background for enrichment test",v=verbose)
  #### Foreground ####
  fg <- subset(subset_DT, eval(parse(text=foreground_filter))) %>%
    dplyr::mutate(chrom = paste0(gsub("chr","",CHR)),
                  chromStart = POS,
                  chromEnd = POS,
                  name = SNP) %>%
    dplyr::select(chrom, chromStart, chromEnd, name)

  #### Background ####
  if(any(is.na(background_filter))){
    ## Optionally, can supply no background at all to XGR
    bg <- NULL
  }else {
    if(!is.null(background_filter)){
      bg_DT <- subset(subset_DT, eval(parse(text=background_filter)))
    } else { bg_DT <- subset_DT  }
    bg <- bg_DT %>% dplyr::mutate(chrom = paste0(gsub("chr","",CHR)),
                                  chromStart = POS,
                                  chromEnd = POS,
                                  name = SNP) %>%
      dplyr::select(chrom, chromStart, chromEnd, name)
  }


  #### Sample fg/bg (for bootstrapping) ####
  if(!is.null(fg_sample_size)){
    printer("XGR:: Sampling",fg_sample_size,"foreground SNPs",v=verbose)
    fg <- fg %>% dplyr::sample_n(size=fg_sample_size)
  }
  if(!is.null(bg_sample_size)){
    printer("XGR:: Sampling",bg_sample_size,"background SNPs",v=verbose)
    bg <- bg %>% dplyr::sample_n(size=bg_sample_size)
  }

  printer("XGR::",nrow(fg),"SNPs in foreground.")
  printer("XGR::",nrow(bg),"SNPs in background")
  return(list("foreground"=fg,
              "background"=bg))
}




#' XGR enrichment
#'
#' @family XGR
#' @examples
#' root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
#'
#' merged_dat <- merge_finemapping_results(dataset = dirname(root), LD_reference = "UKB", minimum_support = 0)
#' merged_dat <- find_consensus_SNPs_no_PolyFun(merged_dat)
#'
#' gr.merged <- merge_celltype_specific_epigenomics()
#' grouping_vars <- c("Study","Cell_type","Assay")
#'
#' enrich.lead <- XGR.enrichment(gr=gr, merged_dat=merged_dat, foreground_filter="leadSNP==T",  grouping_vars=grouping_vars)
XGR.enrichment <- function(gr,
                           merged_dat,
                           foreground_filter="Consensus_SNP==T",
                           background_filter=NULL,
                           grouping_vars=c("Study","Assay","Cell_type"),
                           fg_sample_size=NULL,
                           bg_sample_size=NULL,
                           background.annotatable.only=F,
                           verbose=T){
  fg_bg <- XGR.prepare_foreground_background(subset_DT=merged_dat,
                                             foreground_filter=foreground_filter,
                                             background_filter=background_filter,
                                             fg_sample_size = fg_sample_size,
                                             bg_sample_size = bg_sample_size)
  # Create all combinations
  # if(!is.null(grouping_vars)){
  #   combos <- expand.grid(sapply( subset(data.frame(gr), select=grouping_vars), unique)) %>%
  #     `colnames<-`(grouping_vars)
  #   if(length(grouping_vars)<2) {combos$dummy1 <- 1; gr$dummy1 <- 1; }
  # }else {
  #   combos <- data.frame(dummy1=1, dummy2=2);
  #   gr$dummy1 <- 1;  gr$dummy2 <- 2;
  # }

  combos <- unique(data.frame(gr)[,grouping_vars])
  combos[is.na(combos)] <- "NA"

  printer("+ XGR:: Conducting enrichment tests for",nrow(combos),"combinations of `grouping_vars`.",v=verbose)
  RES <- lapply(1:nrow(combos), function(i,
                                         .background.annotatable.only=background.annotatable.only){
    ROW <- combos[i,]
    # printer("+ XGR::",ROW)
    gr.sub <- gr
    for(column in colnames(combos)){
      gr.sub <- subset(gr.sub, eval(parse(text=column))==ROW[[column]])
    }

    res <- suppressMessages(XGR::xGRviaGenomicAnno(data.file = fg_bg$foreground,
                                background.file = fg_bg$background,
                                format.file="data.frame",
                                GR.annotation = gr.sub,
                                background.annotatable.only = .background.annotatable.only,
                                verbose=F))
    for(column in colnames(combos)){
      res[[column]] <- ROW[[column]]
    }
    return(res)
  }) %>% data.table::rbindlist()
  RES$fg_filter <- if(is.null(foreground_filter)) NA else foreground_filter
  RES$bg_filter <- if(is.null(background_filter)) NA else background_filter
  RES$fg_sample_size <- if(is.null(fg_sample_size)) nrow(fg_bg$foreground) else fg_sample_size
  RES$bg_sample_size <- if(is.null(bg_sample_size)) nrow(merged_dat) else bg_sample_size
  return(RES)
}





#' XGR enrichment (bootstrapped)
#'
#' @family XGR
#' @examples
#' root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
#'
#' merged_dat <- merge_finemapping_results(dataset = dirname(root), LD_reference = "UKB", minimum_support = 0)
#' merged_dat <- find_consensus_SNPs_no_PolyFun(merged_dat)
#' no_no_loci<- c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1", "RETREG3","UBTF","FAM171A2","MAP3K14","CRHR1","MAPT-AS1","KANSL1","NSF","WNT3");
#' merged_dat <- subset(merged_dat, !Locus %in% no_no_loci)
#'
#' gr.merged <- merge_celltype_specific_epigenomics()
#' grouping_vars <- c("Study","Cell_type","Assay")
#'
#' enrich_res <- XGR.enrichment_bootstrap(gr=gr.merged, merged_dat=merged_dat, grouping_vars=grouping_vars,  bootstrap=F, snp_groups=c("Random","GWAS lead","UCS (-PolyFun)","UCS","Consensus (-PolyFun)","Consensus"), nThread=12, save_path=file.path(root,"XGR/celltypespecific_epigenomics.SNP_groups.csv.gz"))
#'
#' enrich_boot <- XGR.enrichment_bootstrap(gr=gr.merged, merged_dat=merged_dat, grouping_vars=grouping_vars,  bootstrap=T, fg_sample_size=NULL, bg_sample_size=100, iterations=100)
XGR.enrichment_bootstrap <- function(gr,
                                     merged_dat,
                                     snp_groups=c("Random","GWAS lead","UCS (-PolyFun)","UCS","Consensus (-PolyFun)","Consensus"),
                                     background_filter=NULL,
                                     grouping_vars=c("Study","Assay","Cell_type"),
                                     iterations=1000,
                                     fg_sample_size=20,
                                     bg_sample_size=NULL,
                                     bootstrap=T,
                                     save_path=F,
                                     nThread=4,
                                     verbose=T){
  if(bootstrap){
    printer("XGR:: Initiating bootstrap enrichment procedure",v=verbose)
  } else {
    iterations <-1; bg_sample_size <- NULL; fg_sample_size <- NULL
  }
  sampling_df <- merged_dat;

  RES_GROUPS <- lapply(snp_groups, function(snp_group,
                                            .merged_dat=merged_dat,
                                            .grouping_vars=grouping_vars,
                                            .background_filter=background_filter,
                                            .fg_sample_size=fg_sample_size,
                                            .bg_sample_size=bg_sample_size){
    snp_filters <- snp_group_filters(random_sample_size=.fg_sample_size)
    .foreground_filter <- snp_filters[snp_group]
    message(snp_group," :: ",.foreground_filter)
    RES <- parallel::mclapply(1:iterations, function(i,
                                                     merged_dat=.merged_dat,
                                                     grouping_vars=.grouping_vars,
                                                    foreground_filter=.foreground_filter,
                                                    background_filter=.background_filter,
                                                    fg_sample_size=.fg_sample_size,
                                                    bg_sample_size=.bg_sample_size){
      try({
        XGR.enrichment(gr=gr,
                       merged_dat=merged_dat,
                       foreground_filter=foreground_filter,
                       background_filter=background_filter,
                       grouping_vars=grouping_vars,
                       fg_sample_size = fg_sample_size,
                       bg_sample_size = bg_sample_size)
      })
    }, mc.cores = nThread) %>% data.table::rbindlist(fill=T)
    RES$SNP_group <- snp_group
    return(RES)
  }) %>% data.table::rbindlist(fill=T)

  # Post-process
  RES_GROUPS <- RES_GROUPS %>%
    dplyr::mutate( SNP_group=factor(SNP_group, levels=unique(SNP_group), ordered = T),
                   FDR=stats::p.adjust(p = pvalue, method = "fdr"))
  if(save_path!=F){
    printer("XGR:: Saving enrichment results ==>",save_path,v=verbose)
    data.table::fwrite(RES_GROUPS,save_path)
  }
  return( RES_GROUPS  )
}





#' Plot enrichment results
#'
#' @family XGR
#' @examples
#' root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
#' ### merged enrichment results
#' enrich_res <- data.table::fread( file.path(root,"XGR/celltypespecific_epigenomics.SNP_groups.csv.gz"))
#' enrich_res <- data.table::fread( file.path(root,"XGR/celltypespecific_epigenomics.snp_groups.csv.gz"))
#' enrich_boot <- data.table::fread(file.path(root,"XGR/celltypespecific_epigenomics.snp_groups.permute.csv.gz"))
#' enrich_assay <- data.table::fread(file.path(root,"XGR/celltypespecific_epigenomics.snp_groups.assay.csv.gz"))
#'
#' # Merged volcano plot
#' enrich_res <- subset(enrich_res, SNP_Group != "Consensus (-PolyFun)") %>% dplyr::rename(SNP_group=SNP_Group)
#' gp <- XGR.enrichment_plot(enrich_res=subset(enrich_res , !Assay %in%  c("HiChIP_FitHiChIP","PLAC")) , title="Enrichment: Cell-type-specific epigenomics", plot_type="point",save_plot=file.path(root,"XGR/celltypespecific_epigenomics.enrich_volcano.png"), height=6, width=8, shape_var="Assay")
#' ## Merged bar plot
#' gp <- XGR.enrichment_plot(enrich_res=enrich_res,  plot_type="bar", facet_formula=".~Assay",FDR_thresh=.05)
#'
#'
#' # Merged volcano plot (permuted)
#' gp <- XGR.enrichment_plot(enrich_res=enrich.scATAC.permute, title="Permuted enrichment: Cell-type-specific peaks and elements", plot_type="point")
XGR.enrichment_plot <- function(enrich_res,
                                title=NULL,
                                subtitle=NULL,
                                facet_formula=NULL,
                                line_formula = "y ~ x",
                                line_method="lm",
                                line_span=1,
                                FDR_thresh=1,
                                plot_type="bar",
                                shape_var="Cell_type",
                                facet_scales="free",
                                show_plot=T,
                                save_plot=F,
                                height=5,
                                width=5){
  enrich_res <- dplyr::mutate(enrich_res,
                              SNP_group=factor(SNP_group, levels = unique(SNP_group), ordered = T),
                              ## Make Random size smaller (otherwise will make everything else relatively tiny)
                              nOverlap=ifelse(SNP_group=="Random",10,nOverlap)
                              )
 sum(enrich_res$fc==-Inf)
  colorDict <- snp_group_colorDict()
  if(plot_type=="bar"){
    gp <- ggplot(data=subset(enrich_res, FDR<=FDR_thresh),
                 aes(x=SNP_group, y= fc, fill=SNP_group)) +
      # geom_col(stat="identity", alpha=.5, show.legend = F) +
      geom_boxplot()+
      geom_jitter(height=0, width = 0, alpha=.1, show.legend = F) +
      scale_fill_manual(values = colorDict) +
      # ggpubr::stat_compare_means(method = method,
      #                            comparisons = comparisons,
      #                            label = "p.signif", size=3, vjust = 1.5) +
      facet_grid(facets = if(is.null(facet_formula)) facet_formula else as.formula(facet_formula),
                 scales="free_y") +
      labs(x="SNP Group", title=title, subtitle=subtitle) +
      theme_bw() +
      theme(strip.background = element_rect(fill="grey20"),
            strip.text = element_text(color="white"),
            axis.text.x = element_text(angle=45, hjust=1))
  }

  if(plot_type=="point"){
    gp <- ggplot(data=subset(enrich_res, FDR<=FDR_thresh),
                 aes(x= log1p(fc), y= -log10(pvalue),
                                      size=nOverlap, color=SNP_group, group=SNP_group,
                                      fill=SNP_group,
                                      shape=eval(parse(text=shape_var)))
                 ) +
      geom_smooth (alpha=0.1, size=0, span=line_span,
                   method = line_method, formula = line_formula) +
      stat_smooth (geom="line", alpha=0.3, size=1, span=line_span,
                   method = line_method, formula = line_formula) +

      geom_point(alpha=.5) +
      scale_color_manual(values = colorDict) +
      scale_fill_manual(values = colorDict) +
      scale_shape_manual(values=12:(12+dplyr::n_distinct(enrich_res[[shape_var]]))) +
      geom_hline(yintercept = -log10(0.05), linetype=2, alpha=.5) +
      facet_grid(facets = if(is.null(facet_formula)) facet_formula else as.formula(facet_formula),
                 scales = facet_scales)  +
      labs(title=title, subtitle=subtitle, shape=shape_var) +
      theme_bw() +
      theme(strip.background = element_rect(fill="grey20"),
            strip.text = element_text(color="white"))
  }

  if(show_plot) print(gp)

  if(save_plot!=F){
    ggsave(save_plot, gp, dpi=400, height=height, width=width)
  }
  return(gp)
}






#' @keywords internal
#' @family XGR
XGR.sep_handler <- function(lib.name){
  # "_(?=[^_]+$)" : Split by the last "_"
  sepDict <- list("ENCODE_TFBS_ClusteredV3_CellTypes"="[.]",
                  "ENCODE_DNaseI_ClusteredV3_CellTypes"="_(?=[^_]+$)",
                  "Broad_Histone"="_(?=[^_]+$)",
                  "FANTOM5_Enhancer"="_(?=[^_]+$)",
                  "TFBS_Conserved"="[$]")
  if(lib.name %in% names(sepDict)){
    sep <- sepDict[[lib.name]]
  }else{ sep <- "_(?=[^_]+$)"}
  return(sep)
}




#' @keywords internal
#' @family XGR
XGR.parse_metadata <- function(gr.lib,
                               lib.name=NA){
  # https://stackoverflow.com/questions/50518137/separate-a-column-into-2-columns-at-the-last-underscore-in-r
  sep <- XGR.sep_handler(lib.name = lib.name)
  GenomicRanges::mcols(gr.lib) <- tidyr::separate(data.frame(GenomicRanges::mcols(gr.lib)),
                                                  sep=sep,
                                                  col="fullname",
                                                  into=c("Source","Assay"),
                                                  extra="merge")
  return(gr.lib)
}




#' Download, standardize, and merge XGR annotations
#'
#' Merges a list of XGR annotations into a single GRanges object
#'
#' @param lib.selections Which XGR annotations to check overlap with.
#' For full list of libraries see:
#'  \url{http://xgr.r-forge.r-project.org/#annotations-at-the-genomic-region-level}
#' @param as_GRangesList Return as a \code{GRangesList}, instead of a single merged \code{GRanges} object.
#' @return GRangesList
#' @family XGR
#' @examples
#' data("BST1")
#' gr.lib <- XGR.download_and_standardize(lib.selections=c("ENCODE_DNaseI_ClusteredV3_CellTypes"), finemap_dat=BST1, nCores=1)
#' @keywords internal
XGR.download_and_standardize <- function(lib.selections=c("ENCODE_TFBS_ClusteredV3_CellTypes",
                                                          "TFBS_Conserved",
                                                          "Uniform_TFBS"),
                                         as_GRangesList=F,
                                         finemap_dat,
                                         nCores=4){
  # Iterate over XGR libraries
  gr.lib <- lapply(lib.selections, function(lib.name){
    GR.annotations <- XGR::xRDataLoader(RData.customised = lib.name)
    # Iterate over lists within each library
    all_GRL <- parallel::mclapply(names(GR.annotations), function(n1){
      grl <- GR.annotations[[n1]]

      # Handle both nested and unnested entries
      if(class(grl)=="list"){
        # grl$name <- names(unlist(grl))
        GRL <- lapply(names(grl), function(n2){
          gr <- grl[[n2]]
          gr$source <- n1
          gr$assay <- n2
          return(gr)
        }) %>% unlist()
      } else {
        grl$name <- names(grl)
        return(grl)
        }
    }, mc.cores = nCores) %>% unlist() # return all_GRL
    # Rename GRanges after they've been unnested
    names(all_GRL) <- names(unlist(GR.annotations))

    # Merge lists together
    if(!is.null(all_GRL)){
      ALL_GRL <- unlist(GenomicRanges::GRangesList(all_GRL))
      ALL_GRL <- GRanges_overlap(dat1 = finemap_dat,
                                 chrom_col.1 = "CHR",
                                 start_col.1 = "POS",
                                 end_col.1 = "POS",
                                 dat2 = ALL_GRL)
      # Parse metadata
      ALL_GRL$library <- lib.name
      ALL_GRL$fullname <- names(ALL_GRL)
      ALL_GRL <- XGR.parse_metadata(gr.lib = ALL_GRL, lib.name = lib.name)
      return(ALL_GRL)
    } else{return(NULL)}
  }) # return OVERLAP
  # Merge
  gr.lib <- GenomicRanges::GRangesList(unlist(gr.lib))
  if(as_GRangesList==F){ gr.lib <- unlist(gr.lib) }
  return(gr.lib)
}



#' Check overlap with XGR annotations
#'
#'
#' Automatically handles different file formats provided by XGR
#'  (e.g. varying kinds of nested/unnested \code{GRanges}).
#' Then returns a \code{Granges} object with only the XGR annotation ranges
#' that overlap with the SNPs in \code{subset_DT}.
#' The \code{GRanges} merges hits from \code{subset_DT}.
#'
#' @param nCores Multi-thread across libraries
#' @param save_path Save the results as a data.frame
#' @inheritParams XGR.prepare_foreground_background
#' @inheritParams XGR.download_and_standardize
#' @family XGR
#' @keywords internal
#' @examples
#' data("BST1")
#' finemap_DT <- BST1
#' gr.hits <- XGR.iterate_overlap(lib.selections=c("ENCODE_TFBS_ClusteredV3_CellTypes"), subset_DT=finemap_DT, nCores=1)
XGR.iterate_overlap <- function(lib.selections=c("ENCODE_TFBS_ClusteredV3_CellTypes",
                                                 "TFBS_Conserved",
                                                 "ReMap_PublicAndEncode_TFBS",
                                                 "Uniform_TFBS"),
                                subset_DT,
                                save_path=F,
                                nCores=4){
  gr.lib <- XGR.download_and_standardize(lib.selections = lib.selections,
                                         finemap_dat = subset_DT,
                                         nCores = nCores)
  # no_no_loci <- c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1",
  #                 "RETREG3","UBTF","FAM171A2","MAP3K14","CRHR1","MAPT-AS1","KANSL1","NSF","WNT3")
  # subset_DT <- subset(subset_DT, !Locus %in% no_no_loci)
  # subset_DT <- assign_lead_SNP(subset_DT)
  gr.hits <- GRanges_overlap(dat1 = subset_DT,
                             chrom_col.1 = "CHR",
                             start_col.1 = "POS",
                             end_col.1 = "POS",
                             dat2 = gr.lib)

  ucs.hits <-subset(gr.hits, Consensus_SNP)
  length(unique(subset(subset_DT, Consensus_SNP)$SNP))
  length(unique(subset(subset_DT, Consensus_SNP)$Locus))

  length(unique(ucs.hits$SNP))
  length(unique(ucs.hits$Locus))

  if(save_path!=F){
    # save_path  <-  "Data/GWAS/Nalls23andMe_2019/_genome_wide/XGR/XGR_TFBS_UCS_overlap.csv.gz"
    dir.create(dirname(save_path), showWarnings = F, recursive = T)
    data.table::fwrite(data.frame(ucs.hits),save_path)
  }
  return(gr.hits)
}


#' Conduct enrichment tests for each annotation
#'
#' XGR uses a binomial enrichment tests for each annotation.
#'
#' @examples
#' data("merged_DT")
#' enrich_res <- XGR.iterate_enrichment(subset_DT=merged_DT, foreground_filter = "Consensus_SNP", background_filter = "leadSNP", lib.selections = c("ENCODE_TFBS_ClusteredV3_CellTypes"), nCores=1)
#' @inheritParams XGR.prepare_foreground_background
#' @family XGR
#' @keywords internal
#' @export
XGR.iterate_enrichment <- function(subset_DT,
                                   foreground_filter = "Consensus_SNP",
                                   background_filter = "leadSNP",
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
                                   save_path=F,
                                   nCores=4){
  # Description of all datasets
  # https://www.rdocumentation.org/packages/XGR/versions/1.1.5/topics/xDefineGenomicAnno
  fg_bg <- XGR.prepare_foreground_background(subset_DT,
                                             foreground_filter = foreground_filter,
                                             background_filter = background_filter)

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
      GR.annotations <- XGR::xRDataLoader(RData.customised = lib.name)
      eTerm <- lapply(GR.annotations, function(grl){
        et <-XGR::xGRviaGenomicAnno(data.file = fg_bg$foreground,
                                    background.file = fg_bg$background,
                                    format.file="data.frame",
                                    GR.annotation = grl)
        return(et)
      }) %>% data.table::rbindlist()
      eTerm$lib <- lib.name
      eTerm$fullname <- names(unlist(GR.annotations))
      eTerm$source <-  lapply(eTerm$fullname, function(e){strsplit(e, "[.]")[[1]][1]}) %>% as.character()
      eTerm$assay <-  lapply(eTerm$fullname, function(e){strsplit(e, "[.]")[[1]][2]}) %>% as.character()
    })
    return(eTerm)
  }, mc.cores = nCores)

  # Re-calculate corrected p-val to account for multiple dbs tested
  enrich_res <- data.table::rbindlist(database_results) %>%
    dplyr::mutate(FDR = p.adjust(p=pvalue, method = "fdr"),
                  Bonf = p.adjust(p=pvalue, method = "bonferroni")) %>%
    dplyr::arrange(FDR, -nOverlap, -fc) %>%
    subset(adjp<0.05) %>%
    data.table::data.table()
  if(save_path!=F){
    data.table::fwrite(enrich_res, save_path, quote = F)
  }
  return(enrich_res)
}



#' Plot XGR enrichment
#'
#' @family XGR
#' @keywords internal
#' @examples
#' data("merged_DT")
#' enrich_res <- XGR.iterate_enrichment(subset_DT=merged_DT, foreground_filter = "Consensus_SNP", background_filter = "leadSNP", lib.selections = c("ENCODE_TFBS_ClusteredV3_CellTypes"), nCores=1)
#' XGR.plot_enrichment(enrich_res)
XGR.plot_enrichment <- function(enrich_res,
                                adjp_thresh=0.05,
                                top_annotations=NULL,
                                show_plot=T){
  # enrich_res$annotation <- factor(gsub(".*[$]","",enrich_res$name), levels = rev(unique(gsub(".*[$]","",enrich_res$name))), ordered = T)
  enrich_res <- dplyr::arrange(enrich_res, desc(fc))
  enrich_res$source <- factor(enrich_res$source, unique(enrich_res$source), ordered = T)
  enrich_res$assay <- factor(enrich_res$assay, unique(enrich_res$assay), ordered = T)
  if(is.null(top_annotations)){top_annotations <- nrow(enrich_res)}

  gp <- ggplot(data = subset(enrich_res, adjp<adjp_thresh)[1:top_annotations,],
               aes(y=fc, x=assay, fill=fc)) +
    geom_col() +
    labs(title=paste0("Epigenomic annotation enrichment (FDR < ",FDR_thresh,")"),
         subtitle="Foreground = Consensus SNPs\nBackground = Lead GWAS SNPs",
         y="Fold-change") +
    facet_grid(facets = Cell_type ~  Assay,
               scales = "free",
               space = "free") +
    theme_bw() +
    theme(strip.placement = "outside",
          strip.text.y.left  = element_text(angle = 0),
          strip.background = element_rect(color="black",fill="white"),
          plot.title = element_text(hjust = .5),
          plot.subtitle = element_text(hjust = .5))
  if(show_plot){print(gp)}
  return(gp)
}



#' Filter sources
#'
#' Identify the sources with the most annotations in the locus.
#' Then only keep these sources.
#' @keywords internal
#' @family XGR
XGR.filter_sources <- function(gr.lib,
                               n_top_sources=5){
  top_sources <-  data.frame(gr.lib) %>%
    dplyr::group_by(library, Source) %>%
    dplyr::tally(sort = T)
  if(!is.null(n_top_sources)){
    gr.filt <- subset(gr.lib, Source %in% unique(top_sources$Source[1:min(n_top_sources, dplyr::n_distinct(top_sources$Source))]))
  }
  return(gr.filt)
}

#' Filter assays
#'
#' Identify the assays with the most annotations in the locus.
#' Then only keep these assays
#' @keywords internal
#' @family XGR
XGR.filter_assays <- function(gr.lib,
                               n_top_assays=5){
  top_assays <-  data.frame(gr.lib) %>%
    dplyr::group_by(library, Assay) %>%
    dplyr::tally(sort = T)
  if(!is.null(n_top_assays)){
    gr.filt <- subset(gr.lib, Assay %in% unique(top_assays$Assay[1:min(n_top_assays, dplyr::n_distinct(top_assays$Assay))]))
  }
  return(gr.filt)
}



#' Plot XGR peaks
#'
#' Plots the distribution of annotations across a genomic region (x-axis).
#'
#' @family XGR
#' @keywords internal
#' @param gr.lib \code{GRanges} object of annotations.
#' @param geom Plot type ("density", or "histogram").
#' @param locus Locus name (\emph{optional}).
#' @param adjust The granularity of the peaks.
#' @param show_plot Print the plot.
#' @return \code{ggbio} track plot.
#' @inheritParams XGR.prepare_foreground_background
#' @examples
#' data("BST1")
#' finemap_DT <- BST1
#' gr.lib <- XGR.download_and_standardize(c("ENCODE_DNaseI_ClusteredV3_CellTypes"), finemap_dat=finemap_DT, nCores=1)
#' gr.filt <- XGR.filter_sources(gr.lib=gr.lib, n_top_sources=5)
#' gr.filt <- XGR.filter_assays(gr.lib=gr.filt, n_top_assays=5)
#' xgr.track <- XGR.plot_peaks(gr.lib=gr.filt, subset_DT=finemap_DT, fill_var="Assay", facet_var="Source")
XGR.plot_peaks <- function(gr.lib,
                           subset_DT,
                           fill_var="Assay",
                           facet_var="Source",
                           geom="density",
                           locus=NULL,
                           adjust=.2,
                           show_plot=T,
                           show.legend=T){
  # data("BST1"); subset_DT <- BST1; show.legend=T;  fill_var="Assay"; facet_var="Source"; geom="density"; adjust=.2;
  gr.lib$facet_label <- gsub("_","\n",GenomicRanges::mcols(gr.lib)[,facet_var])
  xgr.track <- ggbio::autoplot(gr.lib,
                               # which = gr.snp,
                               ggplot2::aes(fill=eval(parse(text=fill_var))),
                               facets = formula("facet_label ~ ."), #formula(paste0(facet_var," ~ .")),
                               # fill = "magenta",
                               color = "white",#NA
                               geom = geom,
                               adjust = adjust,
                               position="stack",
                               # bins=50,
                               size=.1,
                               alpha = .7,
                               show.legend=show.legend)  +
    ggplot2::theme_bw() +
    labs(fill=fill_var) +
    xlim(min(subset_DT$POS), max(subset_DT$POS))
  # ggbio::tracks(list("XGR"=xgr.track))

  if(show_plot) print(xgr.track)
  return(xgr.track)
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




#' Standardize XGR annotations
#'
#' Parses the metadata and adds it as columns,
#' and then merges the results into a single
#' \code{\link[GenomicRanges]{GenomicRangesList}}
#'
#' @param grl.xgr \code{\link[GenomicRanges]{GenomicRangesList}} of XGR queries.
#' @family XGR
#' @keywords internal
XGR.merge_and_process <- function(grl.xgr,
                                  lib,
                                  n_top_sources=10){
  # grl.xgr <- check_saved_XGR(results_path, lib)
  ## Make track
  ## Add and modify columns
  grl.xgr.merged <- unlist(grl.xgr)
  names(grl.xgr.merged) <- gsub("Broad_Histone_","",names(grl.xgr.merged))
  sep <- XGR.sep_handler(lib.name = lib)
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



#' Download XGR annotations
#'
#' @family XGR
#' @keywords internal
XGR.import_annotations <- function(gr.snp,
                                   anno_data_path=file.path("annotations", paste0("XGR_",lib.name,".rds")),
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
