



#' Get overlap between datatable of SNPs and scATAC peaks
#'
#' Can optionally add \code{Cicero} coaccessibility scores,
#' which are also derived from scATAC-seq data.
#'
#' @family CORCES_2020
#' @source \url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
CORCES_2020.get_ATAC_peak_overlap <- function(finemap_dat,
                                              FDR_filter=NULL,
                                              add_cicero=T,
                                              cell_type_specific=T,
                                              verbose=T){
  if(cell_type_specific){
    printer("CORCES_2020:: Extracting overlapping cell-type-specific scATAC-seq peaks",v=verbose)
    dat <- echolocatoR::CORCES_2020.scATACseq_celltype_peaks
    Assay <- 'scATAC'
  } else{
    printer("CORCES_2020:: Extracting overlapping bulkATAC-seq peaks from brain tissue",v=verbose)
    dat <- echolocatoR::CORCES_2020.bulkATACseq_peaks
    Assay <- 'bulkATAC'
  }

  gr.peaks_lifted <- LIFTOVER(dat = dat,
                              build.conversion = "hg38.to.hg19",
                              chrom_col = "hg38_Chromosome",
                              start_col = "hg38_Start",
                              end_col = "hg38_Stop",
                              verbose=F)
  # Get overlap with PEAKS
  gr.hits <- GRanges_overlap(dat1 = finemap_dat,
                             chrom_col.1 = "CHR",
                             start_col.1 = "POS",
                             end_col.1 = "POS",
                             dat2 = gr.peaks_lifted)
  gr.hits$Assay <- Assay
  if(!is.null(FDR_filter)){
    gr.hits <- subset(gr.hits, FDR < FDR_filter)
  }

  if(add_cicero & cell_type_specific){
    try({
      # Pretty sure the Peak_IDs are shared between the sc-ATACseq data and cicero,
      # because Cicero derives coaccess from sc-ATAC-seq data:
      # http://www.cell.com/molecular-cell/retrieve/pii/S1097276518305471?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1097276518305471%3Fshowall%3Dtrue
      ## Also pretty sure that checking for cicero overlap only in the scATACseq gr.hits object is ok
      # bc you can only test for coaccessibility if there's a peak to begin with.
      cicero <- echolocatoR::CORCES_2020.cicero_coaccessibility
      cicero_dict <- c(setNames(cicero$Coaccessibility, cicero$Peak_ID_Peak1),
                      setNames(cicero$Coaccessibility, cicero$Peak_ID_Peak2))
      gr.hits$Cicero <- cicero_dict[gr.hits$Peak_ID]
      gr.cicero <- subset(gr.hits, !is.na(Cicero))
      gr.cicero$Assay <- "Cicero"
      printer("+ CORCES_2020:: Cicero coaccessibility scores identified for",
              length(gr.cicero),"/",length(gr.hits),"peak hits.",v=verbose)
      gr.hits <- rbind_GRanges(gr1 = gr.hits, gr2 = gr.cicero)
    })
  }
  if(cell_type_specific==F){
    gr.hits$brain <- 1
  }
  return(gr.hits)
}






#' Get overlap between data table of SNPs and HiChIP_FitHiChIP coaccessibility anchors
#'
#' Anchors are the genomic regions that have evidence of being
#' functionally connected to one another (coaccessible),
#'  e.g. enhancer-promoter interactions.
#'
#' @family CORCES_2020
#' @source \url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
CORCES_2020.get_HiChIP_FitHiChIP_overlap <- function(finemap_dat,
                                                     verbose=T){
    loops <- echolocatoR::CORCES_2020.HiChIP_FitHiChIP_loop_calls
    # Anchor 1
    gr.anchor1 <- LIFTOVER(dat = loops,
                           build.conversion = "hg38.to.hg19",
                           chrom_col = "hg38_Chromosome_Anchor1",
                           start_col = "hg38_Start_Anchor1",
                           end_col = "hg38_Stop_Anchor1",
                           verbose=F)
    gr.anchor1_hits <- GRanges_overlap(dat1 = finemap_dat,
                                      chrom_col.1 = "CHR",
                                      start_col.1 = "POS",
                                      end_col.1 = "POS",
                                      dat2 = gr.anchor1)
    gr.anchor1_hits$Anchor <- 1

    # Anchor 2
    gr.anchor2 <- LIFTOVER(dat = loops,
                           build.conversion = "hg38.to.hg19",
                           chrom_col = "hg38_Chromosome_Anchor2",
                           start_col = "hg38_Start_Anchor2",
                           end_col = "hg38_Stop_Anchor2",
                           verbose=F)
    gr.anchor2_hits <- GRanges_overlap(dat1 = finemap_dat,
                                      chrom_col.1 = "CHR",
                                      start_col.1 = "POS",
                                      end_col.1 = "POS",
                                      dat2 = gr.anchor2)
    gr.anchor2_hits$Anchor <- 2
    # Merge and report
    gr.anchor <- rbind_GRanges(gr.anchor1_hits, gr.anchor2_hits)
    gr.anchor$Assay <- "HiChIP_FitHiChIP"
    # Have to make a pseudo cell-type col bc (i think) this analysis was done on bulk data
    gr.anchor$brain <- 1
    printer("+ CORCES_2020:: Found",length(gr.anchor),
            "hits with HiChIP_FitHiChIP coaccessibility loop anchors.",v=verbose)
    return(gr.anchor)
}









#' Prepare data to plot overlap between datatable of SNPs and
#' cell-type-specific epigenomic peaks and coaccessibility data.
#'
#' @family CORCES_2020
#' @source \url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
#' @examples
#' data("merged_DT");
#' finemap_dat <- subset(merged_DT, Consensus_SNP)
#' dat_melt <- CORCES_2020.prepare_scATAC_peak_overlap(merged_DT=merged_DT)
CORCES_2020.prepare_scATAC_peak_overlap <- function(merged_DT,
                                                    FDR_filter=NULL,
                                                    snp_filter="Consensus_SNP==T",
                                                    add_cicero=T,
                                                    annotate_genes=T,
                                                    verbose=T){
  cell_dict <- c(ExcitatoryNeurons="neurons (+)",
                 InhibitoryNeurons="neurons (-)",
                 NigralNeurons="neurons (nigral)",
                 Microglia="microglia",
                 Oligodendrocytes="oligo",
                 Astrocytes="astrocytes",
                 OPCs="OPCs")
  # Get SNP groups
  finemap_dat <- subset(merged_DT, eval(parse(text=snp_filter)), .drop=F) #finemap_dat[eval(parse(text=snp_filter))]
  # Get overlap with PEAKS and merge
  gr.hits <- CORCES_2020.get_ATAC_peak_overlap(finemap_dat = finemap_dat,
                                               add_cicero = add_cicero,
                                               cell_type_specific = T,
                                               verbose = verbose)

  annot_cols=NULL
  if(annotate_genes){
    annot_cols <- c("Gene_Symbol","CTCF","Distance_To_TSS","Annotation")
      printer("CORCES_2020:: Annotating peaks by cell-type-specific target genes",v=verbose)
      peak_overlap <- subset(echolocatoR::CORCES_2020.scATACseq_peaks, Peak_ID %in% unique(gr.hits$Peak_ID))
      prefix=""
    for(column in annot_cols){
      dict <- setNames(peak_overlap[[column]], peak_overlap$Peak_ID)
      GenomicRanges::mcols(gr.hits)[[paste0(prefix,column)]] <-  dict[gr.hits$Peak_ID]
    }
    annot_cols <- paste0(prefix,annot_cols)
  }

  # Melt cell type into one col
  cell_melt <- data.table::melt.data.table(data.table::data.table(data.frame(gr.hits)),
                              measure.vars = names(cell_dict),
                              variable.name = "cell_name",
                              value.name = "cell_value")
  cell_melt[cell_melt$cell_value==1, "Cell_type"] <- cell_melt[cell_melt$cell_value==1, "cell_name"]
  cell_melt <- subset(cell_melt, !is.na(Cell_type), .drop=F)

  dat_melt <- count_and_melt(merged_annot = cell_melt,
                             grouping_vars = c("Locus","Cell_type","Assay",annot_cols),
                             snp_filter = snp_filter)
  dat_melt$Cell_type <- cell_dict[dat_melt$Cell_type]
  if(sum(dat_melt$Count==0 | is.na(dat_melt$Count), na.rm = T)>0){
    try({
      dat_melt[dat_melt$Count==0 | is.na(dat_melt$Count),"Count"] <- NA
    })
  }
  dat_melt <- subset(dat_melt, !is.na(Count))
  dat_melt$background <- NA

  # Make sure locus order kept
  locus_order <- SUMMARISE.get_CS_counts(merged_DT)
  dat_melt$Locus <- factor(dat_melt$Locus,  levels = locus_order$Locus, ordered = T)
  return(dat_melt)
}




#' Prepare data to plot overlap between datatable of SNPs and
#' cell-type-specific epigenomic peaks and coaccessibility data.
#'
#' @family CORCES_2020
#' @source \url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
#' @examples
#' data("merged_DT");
#' finemap_dat <- subset(merged_DT, Consensus_SNP)
#' dat_melt <- CORCES_2020.prepare_bulkATAC_peak_overlap(merged_DT=merged_DT)
CORCES_2020.prepare_bulkATAC_peak_overlap <- function(merged_DT,
                                                      FDR_filter=NULL,
                                                      snp_filter="Consensus_SNP==T",
                                                      add_HiChIP_FitHiChIP=T,
                                                      annotate_genes=F,
                                                      verbose=T){
  # Get SNP groups
  finemap_dat <- subset(merged_DT, eval(parse(text=snp_filter)), .drop=F) #finemap_dat[eval(parse(text=snp_filter))]
  # Get overlap with PEAKS and merge
  gr.hits <- CORCES_2020.get_ATAC_peak_overlap(finemap_dat = finemap_dat,
                                               add_cicero = F,
                                               cell_type_specific = F,
                                               verbose = verbose)
  annot_cols=NULL
  if(annotate_genes){
    annot_cols <- c("Gene_Symbol","CTCF","Distance_To_TSS","Annotation")
    printer("CORCES_2020:: Annotating peaks by bulk brain target genes",v=verbose)
    peak_overlap <- subset(echolocatoR::CORCES_2020.bulkATACseq_peaks, Peak_ID %in% unique(gr.hits$Peak_ID))
    prefix=""
    for(column in annot_cols){
      dict <- setNames(peak_overlap[[column]], peak_overlap$Peak_ID)
      GenomicRanges::mcols(gr.hits)[[paste0(prefix,column)]] <-  dict[gr.hits$Peak_ID]
    }
    annot_cols <- paste0(prefix,annot_cols)
  }

  if(add_HiChIP_FitHiChIP){
    gr.anchor_hits <- CORCES_2020.get_HiChIP_FitHiChIP_overlap(finemap_dat = finemap_dat,
                                                               verbose=T)
    gr.hits <- c(gr.hits, gr.anchor_hits)
  }

  # Melt cell type into one col
  cell_melt <- data.table::melt.data.table(data.table::data.table(data.frame(gr.hits)),
                                           measure.vars = c("brain"),
                                           variable.name = "cell_name",
                                           value.name = "cell_value")
  cell_melt[cell_melt$cell_value==1, "Cell_type"] <- cell_melt[cell_melt$cell_value==1, "cell_name"]
  cell_melt <- subset(cell_melt, !is.na(Cell_type), .drop=F)

  dat_melt <- count_and_melt(merged_annot = cell_melt,
                             grouping_vars = c("Locus","Cell_type","Assay",annot_cols),
                             snp_filter = snp_filter)
  if(sum(dat_melt$Count==0 | is.na(dat_melt$Count),na.rm = T)>0){
    try({
      dat_melt[dat_melt$Count==0 | is.na(dat_melt$Count),"Count"] <- NA
    })
  }

  dat_melt <- subset(dat_melt, !is.na(Count))
  dat_melt$background <- NA

  # Make sure locus order kept
  locus_order <- SUMMARISE.get_CS_counts(merged_DT)
  dat_melt$Locus <- factor(dat_melt$Locus,  levels = locus_order$Locus, ordered = T)
  return(dat_melt)
}



CORCES_2020.scATAC_to_GRanges <- function(standardize_cellTypes=F){
  scATAC <- data.table::melt.data.table(echolocatoR::CORCES_2020.scATACseq_celltype_peaks,
                                        measure.vars = c("ExcitatoryNeurons","InhibitoryNeurons","NigralNeurons","Microglia","Oligodendrocytes","Astrocytes","OPCs"),
                                        variable.name = "Cell_type") %>% subset(value==1)
  scATAC$Assay <- "scATAC"
  scATAC$Study <- "Corces2020.peaks"
  gr.Corces2020.peaks <- LIFTOVER(dat = scATAC,
                                  build.conversion = "hg38.to.hg19",
                                  chrom_col = "hg38_Chromosome",
                                  start_col = "hg38_Start",
                                  end_col = "hg38_Stop")
  if(standardize_cellTypes){
    gr.Corces2020.peaks$Cell_type <- standardize_celltypes(gr.Corces2020.peaks$Cell_type)
  }
  return(gr.Corces2020.peaks)
}




