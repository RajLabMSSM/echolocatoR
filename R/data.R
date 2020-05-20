

#' echolocatoR output example
#'
#' An example results file after running \code{\link{finemap_loci}} on the \emph{BST1} locus.
#' Data originally comes from the Parkinson's disease GWAS
#' by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
#'
#' @format data.table
#' \describe{
#'   \item{SNP}{SNP RSID}
#'   \item{CHR}{Chromosome}
#'   \item{POS}{Genomic positiion (in basepairs)}
#'   ...
#' }
#' @source \url{https://www.biorxiv.org/content/10.1101/388165v3}
"finemap_DT"




#' TopSS example file
#'
#' Summary stats of the top SNP(s) per locus.
#' Used to query locus subsets.for fine-mapping.
#'
#' Data from \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}, Table S2.
#' @source  \url{https://www.biorxiv.org/content/10.1101/388165v3}
"Nalls_top_SNPs"


#-------Corces et al. (bioRxiv) data --------

# dat <- readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Coceres_2020/STable2_Features_bulkATAC-seq_Peaks.xlsx", skip = 18)
# CORCES_2020.bulkATACseq_peaks <- data.table::data.table(dat)
# usethis::use_data(CORCES_2020.bulkATACseq_peaks)
#' bulkATACseq peaks from Alzheimer's disease brain tissue
#'
#' Each row represents an individual peak identified in the bulk ATAC-seq data.
#'
#' Data originally from \href{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}{Corces et al. (bioRxiv)}, as of May 2020.
#' Specifically: \emph{STable2_Features_bulkATAC-seq_Peaks}
#'
#' @family CORCES_2020
#' @source \url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
"CORCES_2020.bulkATACseq_peaks"




# dat <- readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Coceres_2020/STable5_Features_scATAC-seq_Peaks_all.xlsx", skip = 18)
# CORCES_2020.scATACseq_peaks <- data.table::data.table(dat)
# usethis::use_data(CORCES_2020.scATACseq_peaks, overwrite = T)
#' scATACseq peaks from Alzheimer's disease brain tissue
#'
#' Each row represents an individual peak identified in the single-cell ATAC-seq data.
#'
#' Data originally from \href{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}{Corces et al. (bioRxiv)}, as of May 2020.
#' Specifically: \emph{STable5_Features_scATAC-seq_Peaks_all}
#' @family CORCES_2020
#' @source \url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
"CORCES_2020.scATACseq_peaks"




# dat <- readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Coceres_2020/STable6_Features_scATAC-seq_celltype_Peaks.xlsx", skip = 15)
# CORCES_2020.scATACseq_celltype_peaks <- data.table::data.table(dat)
# usethis::use_data(CORCES_2020.scATACseq_celltype_peaks)
#' scATACseq cell type-specific peaks from Alzheimer's disease brain tissue
#'
#' Each row represents an individual peak identified from the feature binarization analysis (see methods).
#'
#' Data originally from \href{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}{Corces et al. (bioRxiv)}, as of May 2020.
#' Specifically: \emph{STable6_Features_scATAC-seq_celltype_Peaks}
#'
#' @family CORCES_2020
#' @source \url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
"CORCES_2020.scATACseq_celltype_peaks"




# dat <- readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Coceres_2020/STable10_Coacessibility_Peak_loop_connection.xlsx", skip = 19, sheet=1)
# CORCES_2020.HiChIP_FitHiChIP_loop_calls <- data.table::data.table(dat)
# usethis::use_data(CORCES_2020.HiChIP_FitHiChIP_loop_calls)
#' FitHiChIP loop calls from Alzheimer's disease brain tissue
#'
#' FitHiChIP loop calls that overlap SNPs derived from analysis of H3K27ac HiChIP data.
#' Each row represents an individual peak identified from the feature binarization analysis (see methods).
#'
#' Data originally from \href{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}{Corces et al. (bioRxiv)}, as of May 2020.
#' Specifically: \emph{STable10_Coacessibility_Peak_loop_connection}, \emph{HiChIP FitHiChIP Loop Calls} sheet.
#' \strong{Column dictionary:}
#' hg38_Chromosome_Anchor1 - The hg38 chromosome of the first loop Anchor.
#' hg38_Start_Anchor1 - The hg38 start position of the first loop Anchor.
#' hg38_Stop_Anchor1 - The hg38 stop position of the first loop Anchor.
#' Width_Anchor1 - The width of the first loop Anchor.
#' hg38_Chromosome_Anchor2 - The hg38 chromosome of the second loop Anchor.
#' hg38_Start_Anchor2 - The hg38 start position of the second loop Anchor.
#' hg38_Stop_Anchor2 - The hg38 stop position of the second loop Anchor.
#' Width_Anchor2 - The width of the second loop Anchor.
#' Score - The -log10(q-value) of the loop call from FitHiChIP.
#' Anchor1_hasSNP - A boolean variable determining whether the first anchor overlaps a SNP from our AD/PD GWAS analyses.
#' Anchor2_hasSNP - A boolean variable determining whether the second anchor overlaps a SNP from our AD/PD GWAS analyses.
#'
#' @family CORCES_2020
#' @source \url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
"CORCES_2020.scATACseq_celltype_peaks"




# dat <- readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Coceres_2020/STable10_Coacessibility_Peak_loop_connection.xlsx", skip = 21, sheet=2)
# CORCES_2020.cicero_coaccessibility <- data.table::data.table(dat)
# usethis::use_data(CORCES_2020.cicero_coaccessibility)
#' Cicero_coaccessibility from Alzheimer's disease brain tissue
#'
#' Cicero coaccessibility analysis for peaks that overlap SNPs derived from analysis of scATAC-seq data.
#' Each row represents an individual peak identified from the feature binarization analysis (see methods).
#'
#' Data originally from \href{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}{Corces et al. (bioRxiv)}, as of May 2020.
#' Specifically: \emph{STable10_Coacessibility_Peak_loop_connection}, \emph{Cicero Coaccessibility} sheet.
#' Peak_ID_Peak1 - A unique number that identifies the peak across supplementary tables.
#' \strong{Columns dictionary}:
#' hg38_Chromosome_Peak1 - The hg38 chromosome of the first loop Peak.
#' hg38_Start_Peak1 - The hg38 start position of the first loop Peak.
#' hg38_Stop_Peak1 - The hg38 stop position of the first loop Peak.
#' Width_Peak1 - The width of the first loop Peak.
#' Peak_ID_Peak2 - A unique number that identifies the peak across supplementary tables.
#' hg38_Chromosome_Peak2 - The hg38 chromosome of the second loop Peak.
#' hg38_Start_Peak2 - The hg38 start position of the second loop Peak.
#' hg38_Stop_Peak2 - The hg38 stop position of the second loop Peak.
#' Width_Peak2 - The width of the second loop Peak.
#' Coaccessibility - The coaccessibility correlation for the given peak pair.
#' Peak1_hasSNP - A boolean variable determining whether the first peak overlaps a SNP from our AD/PD GWAS analyses.
#' Peak2_hasSNP - A boolean variable determining whether the second peak overlaps a SNP from our AD/PD GWAS analyses.
#'
#' @family CORCES_2020
#' @source \url{https://www.biorxiv.org/content/10.1101/2020.01.06.896159v1}
"CORCES_2020.scATACseq_celltype_peaks"




# -------------IMPACT ------------ #

# IMPACT_annotation_key <- data.table::fread("~/Desktop/Fine_Mapping/echolocatoR/annotations/IMPACT/IMPACT_annotation_key.txt.gz")
# usethis::use_data(IMPACT_annotation_key, overwrite = T)
#' IMPACT annotation key
#'
#' Metadata for each of the IMPACT (Inference and Modeling of Phenotype-related ACtive Transcription) annotation files.
#' Originally from \href{https://www.cell.com/ajhg/supplemental/S0002-9297%2819%2930108-9}{Amariuta et al. (2019)}.
#'
#' @family IMPACT
#' @source \url{https://github.com/immunogenomics/IMPACT}
"IMPACT_annotation_key"




# --------------Nott et al. (2019) -------------

# file <- "~/Desktop/Fine_Mapping/echolocatoR/annotations/Nott_2019/aay0793-Nott-Table-S5.xlsx"
# sheets <- readxl::excel_sheets(file)
# NOTT_2019.interactome <- lapply(sheets, function(s){readxl::read_excel(file, sheet=s, skip=2)})
# names(NOTT_2019.interactome) <- sheets
# usethis::use_data(NOTT_2019.interactome)
#' Brain cell type-specific enhancers, promoters, and interactomes
#'
#' Originally from \href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}.
#' Specifically: \emph{aay0793-Nott-Table-S5.xlsx}.
#'
#' @family NOTT_2019
#' @source \url{https://science.sciencemag.org/content/366/6469/1134}
"NOTT_2019.interactome"




# NOTT_2019.superenhancer_interactome <- data.table::data.table(readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Nott_2019/aay0793-Nott-Table-S6.xlsx", skip=2)  )
# usethis::use_data(NOTT_2019.superenhancer_interactome)
#' Brain cell type-specific interactomes with superenhancers
#'
#' Originally from \href{https://science.sciencemag.org/content/366/6469/1134}{Nott et al. (2019)}.
#' Specifically: \emph{aay0793-Nott-Table-S6.xlsx}.
#'
#' @family NOTT_2019
#' @source \url{https://science.sciencemag.org/content/366/6469/1134}
"NOTT_2019.superenhancer_interactome"




# NOTT_2019.bigwig_metadata <- data.table::data.table(readxl::read_excel("~/Desktop/Fine_Mapping/echolocatoR/annotations/Nott_2019/Nott_2019.snEpigenomics.xlsx"))
# usethis::use_data(NOTT_2019.bigwig_metadata)
#' Metadata and links to data
#'
#' Metadata for cell type-specific epigenomic bigWig files hosted on UCSC Genome Browser.
#' bigWig files contain the genomic ranges from each epigenomic assay,
#' as well as a Score column which describes the peaks of the aggregate reads.
#' @family NOTT_2019
#' @source \url{https://science.sciencemag.org/content/366/6469/1134}
"NOTT_2019.bigwig_metadata"

