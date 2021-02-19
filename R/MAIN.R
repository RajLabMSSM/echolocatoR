#  ((((((((((((((((()))))))))))))))))  #
# ((((((((((((((((())))))))))))))))) #
#    ------- echolocatoR -------    #
# ((((((((((((((((())))))))))))))))) #
#              ((((()))))              #

# Author: Brian M. Schilder
# Bioinformatician II
# Icahn School of Medicine at Mount Sinai
# New York City, New York, USA
# https://bschilder.github.io/BMSchilder

# =/\                  /\=
#  / \'._   (\_/)   _.'/ \
#  / .''._'--(o.o)--'_.''. \
#  /.' _/ |`'=/ " \='`| \_ `.\
#  /` .' `\;-,'\___/',-;/` '. '\
#  /.-'       `\(-V-)/`       `-.\
# `              v  v           `
#
#   message("      =/\                  /\=    ")
#   message("      / \'._   (\_/)   _.'/ \     ")
#   message("     / .''._'--(o.o)--'_.''. \    ")
#   message("    /.' _/ |`'=/ { \='`| \_ `.\   ")
#   message("   /` .' `\;-,'\___/',-;/` '. '\  ")
#   message("  /.-'       `\(-V-)/`       `-.\ ")
#   message(" `              v  v           `  ")
#


# R package creation:
## http://r-pkgs.had.co.nz

# R package remote dependencies:
## https://cran.r-project.org/web/packages/devtools/vignettes/dependencies.html

# Turn on rmarkdown syntax within Roxygen
## https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html

# regex commands:
## http://www.endmemo.com/program/R/gsub.php



# library(R.utils)
# library(devtools)
# library(DT)
# library(data.table)
# library(dplyr)
# library(ggplot2)
# library(plotly)
# library(patchwork) #devtools::install_github("thomasp85/patchwork")
# library(pbmcapply); #devtools::install_github("kvnkuang/pbmcapply", ref = "dev")


# library(ggrepel)
# library(curl)
# library(gaston)
# library(tidyr)
# library(BiocManager)
# library(biomaRt) # BiocManager::install("biomaRt")
# library(refGenome)
# library(crayon)
# library(ChIPQC); BiocManager::install('ChIPQC')

# library(coloc)
# install.packages("haploR", dependencies = TRUE)
# library(haploR)
# library(GeneOverlap) #BiocManager::install("GeneOverlap")
# library(rtracklayer) #BiocManager::install("rtracklayer")
# BiocManager::install(c("supraHex","graph","Rgraphviz","dnet"))
# library(XGR)# install.packages("XGR")



# library(Rsamtools) # BiocManager::install("Rsamtools")

# install_fGWAS <- function(){
#   devtools::install_github("wzhy2000/fGWAS/pkg")
#   system("git clone https://github.com/wzhy2000/fGWAS.git & cd fGWAS & R CMD INSTALL pkg")
# }
# library(fGWAS)
# library(snpStats) #BiocManager::install("snpStats")
# library(bigsnpr) # BiocManager::install("bigsnpr")

# *** SUSIE ****
# library(knitrBootstrap) #install_github('jimhester/knitrBootstrap')
# library(susieR) # devtools::install_github("stephenslab/susieR")

# *** finemapr ****
## finemapr contains: finemap, CAVIAR, and PAINTOR
# library(finemapr) # devtools::install_github("variani/finemapr")
# library(roxygen2) #roxygenize()

# *** locuscomparer ****
# https://github.com/boxiangliu/locuscomparer
# library(locuscomparer); #devtools::install_github("boxiangliu/locuscomparer")

# thm <- knitr::knit_theme$get("bipolar")
# knitr::knit_theme$set(thm)



#' Run \pkg{echolocatoR} pipeline on a single locus
#'
#' Unlike \code{finemap_loci}, you don't need to provide a \code{top_SNPs} data.frame.
#' Instead, just manually provide the coordinates of the locus you want to fine-map.
#'
#' The primary functions of \pkg{echolocatoR} that expedite fine-mapping
#'  by wrapping many other \pkg{echolocatoR} functions into one.
#'  Encompasses steps including:
#'  \describe{
#'  \item{Subset & standardize}{Extract subsets of the full summmary stats GWAS or QTL file and reformat them to be compatible with \pkg{echolocatoR}'s various functions }
#'  \item{Calculate linkage disequilibrium}{Download and prepare the necessary LD matrix.}
#'  \item{Fine-map}{Run various fine-mapping tools and merge the results into a single multi-finemap data.frame.}
#'  \item{Plot}{Summarise the results in a multi-track plot for each locus.}
#'  }
#'
#'
#' @section input file parameters:
#'
#' @param loci Character list of loci in \strong{Locus} col of \code{top_SNPs}.
#' @param fullSS_path Path to the full summary statistics file (GWAS or QTL) that you want to fine-map.
#' It is usually best to provide the absolute path rather than the relative path.
#' @param results_dir Where to store all results.
#' \strong{IMPORTANT!:} It is usually best to provide the absolute path rather than the relative path.
#' This is especially important for \emph{FINEMAP}.
#' @param file_sep The separator in the full summary stats file.
#' This parameter is only necessary if \code{query_by!="tabix"}.
#' @param query_by Choose which method you want to use to extract locus subsets from the full summary stats file.
#' Methods include:
#' \describe{
#' \item{"tabix"}{Convert the full summary stats file in an indexed tabix file. Makes querying lightning fast after the initial conversion is done. (\emph{default})}
#' \item{"coordinates"}{Extract locus subsets using min/max genomic coordinates with \emph{awk}.}
#' }
#' @param dataset_name The name you want to assign to the dataset being fine-mapped,
#' This will be used to name the subdirectory where your results will be stored
#' (e.g. \emph{Data/GWAS/<dataset_name>}).
#' Don't use special characters (e.g.".", "/").
#' @param dataset_type The kind dataset you're fine-mapping (e.g. GWAS, eQTL, tQTL).
#' This will also be used when creating the subdirectory where your results will be stored
#' (e.g. \emph{Data/<dataset_type>/Kunkle_2019}).
#' @param top_SNPs A data.frame with the genomic coordinates of the lead SNP for each locus.
#' The lead SNP will be used as the center of the window when extracting subset from the full GWAS/QTL summary statistics file.
#' Only one SNP per \strong{Locus} should be included.
#' At minimum, \code{top_SNPs} should include the following columns:
#' \describe{
#' \item{\emph{Locus}}{A unique name for each locus. Often, loci are named after a relevant gene (e.g. LRRK2) or based on the name/coordinates of the lead SNP (e.g. locus_chr12_40734202) }
#' \item{\emph{CHR}}{The chromosome that the SNP is on. Can be "chr12" or "12" format.}
#' \item{\emph{POS}}{The genomic position of the SNP (in basepairs)}
#' }
#'
#' @section input file column names:
#'
#' @param chrom_col Name of the chromosome column in the full summary stats file.
#' Can be "chr1" or "1" format.
#' (\emph{default: ="CHR"})
#' @param position_col Name of the genomic position column in the full summary stats file.
#' Must be in units of basepairs.
#' (\emph{default: ="POS"})
#' @param snp_col Name of the SNP RSID column in the full summary stats file.
#' (\emph{default: ="SNP"})
#' @param pval_col Name of the p-value column in the full summary stats file.
#' Raw p-values are preferred, but if not available corrected p-values (e.g. FDR) can be used instead.
#' (\emph{default: ="P"})
#' @param effect_col Name of the effect size column in the full summary stats file.
#' Effect size is preferred, but if not available other metrics like Beta for Odds Ratio can be used instead.
#' (\emph{default: ="Effect"})
#' @param stderr_col Name of the standard error  column in the full summary stats file.
#' You can also set \code{stderr_col="calculate"} to infer standard error using: \code{effect / tstat}.
#' (\emph{default: ="StdErr"})
#' @param tstat_col Name of the t-statistic column in the full summary stats file.
#' This column is not necessary unless \code{stderr_col="calculate"} or the standard error column is missing.
#' (\emph{default: ="t-stat"})
#' @param locus_col Name of the locus column in the full summary stats file.
#' (\emph{default: ="Locus"})
#' @param freq_col Name of the allele frequency column in the full summary stats file.
#' Effect allele frequency is preferred, but the non-effect allele can be provided instead (though this may be less accurate).
#' This column is not necessary unless \code{MAF_col="calculate"} or the MAF column is missing.
#' (\emph{default: ="Freq"})
#' @param MAF_col Name of the minor allele frequency column in the full summary stats file.
#' Can be inferred from \strong{freq_col} if missing from the dataset.
#' (\emph{default: ="MAF"})
#' @param A1_col Name of the effect/risk allele column in the full summary stats.
#'  \strong{\emph{IMPORTANT}}: Make sure this actually the case for your full summary stats file.
#' Unfortunately, different studies report different kinds of allele information in a non-standardized way.
#' Meaning that A1/A2 can refer to any number of things:
#'  \describe{
#'  \item{effect/other alleles}{in the case of diseases}
#'  \item{ref/alt alleles}{where ref is the reference genome being used}
#'  \item{major/minor alleles}{This dichotomy holds true for bi-allelic SNPs but not necessary multi-allelic SNPs}
#'  }
#'  This makes comparing summary stats across GWAS/QTL datasets very confusing for several reasons:
#'  \describe{
#'  \item{Multi-allelic SNPs}{SNPs can have more than just 2 possible alleles (multi-allelic SNPs). Even if you compare the same SNP between two studies, you may accidentally be comparing totally different alleles.}
#'  \item{Valence}{The valence (+/-) of per-SNP GWAS effect sizes/beta can be relative to different allele types between studies.
#'  For example, let's say in one GWAS study your effect size for SNP A is 1.5 relative to the major allele in one study,
#'   and the minor allele happens to be the one found in the reference genome.
#'   You then try to compare that effect size to that of the same SNP in another GWAS.
#'   But, the valence of the effect sizes in the 2nd GWAS study are all relative to the reference genome (instead of the minor allele),
#'   giving the same SNP a value of -1.2. If you took the effect sizes at face value you'd say the signals are in opposite directions.
#'   But once you take into account how the valences were determined in each study you realize that they're actually both positive relative to the major allele.}
#'  }
#' This process of reversing per-SNP valences based on aligning the alleles is known as allele flipping.
#' This is important when comparing individual SNPs, but can also have an impact on colocalization results.
#' @param gene_col For QTL studies, the name of the [e]gene column in the full summary stats file (\emph{default: "gene"}).
#' This column will be used for filtering summary stats if supplying a named list of gene:Locus pairs to \code{loci}.
#' @param N_cases_col Name of the column in the full summary stats that has the number of case subjects in the study.
#' This can either be per SNP sample sizes, or one number repeated across all rows.
#' Proxy cases (e.g. relatives of people with the disease being investigated) should be included in this estimate if any were used in the study.
#' This column is not necesssary if \code{N_cases} parameter is provided.
#' (\emph{default: ="N_cases"})
#' @param N_controls_col Name of the column in the full summary stats that has the number of control subjects in the study.
#'  This can either be per SNP sample sizes, or one number repeated across all rows.
#'  This column is not necesssary if \code{N_controls} parameter is provided.
#' (\emph{default: ="N_controls"})
#' @param N_cases The number of case subjects in the study.
#'  Instead of providing a redundant \strong{N_cases_col} column, you can simply enter one value here.
#' @param N_controls The number of control subjects in the study.
#'  Instead of providing a redundant \strong{N_controls_col} column, you can simply enter one value here.
#' @param proportion_cases The proportion of total subjects in the study that were cases.
#'  if \code{proportion_cases="calculate"} then this is inferred:  \code{N_controls / N_controls}.
#' @param sample_size The overall sample size of the study.
#' If none is given, and \strong{N_cases} and \strong{N_controls} columns are present,
#' then sample_size is inferred to be:  \code{max(N_cases) + max(N_controls)}.
#'
#' @section overwrite existing files:
#'
#' @param force_new_subset By default, if a subset of the full summary stats file for a given locus is already present,
#' then \pkg{echolocatoR} will just use the pre-existing file.
#' Set \code{force_new_subset=T} to override this and extract a new subset.
#' Subsets are saved in the following path structure:
#' \emph{Data/<dataset_type>/<dataset_name>/<locus>/Multi-finemap/<locus>_<dataset_name>_Multi-finemap.tsv.gz}
#' @param force_new_LD  By default, if an LD matrix file for a given locus is already present,
#' then \pkg{echolocatoR} will just use the preexisting file.
#' Set \code{force_new_LD=T} to override this and extract a new subset.
#' @param force_new_finemap By default, if an fine-mapping results file for a given locus is already present,
#' then \pkg{echolocatoR} will just use the preexisting file.
#' Set \code{force_new_finemap=T} to override this and re-run fine-mapping.
#'
#' @section fine-mapping parameters:
#'
#' @param finemap_methods Which fine-mapping methods you want to use.
#' @param bp_distance The width of the window size you want each locus to be.
#' For example, if \code{bp_distance=500000} then the locus will span 500kb from the lead SNP in either direction,
#' resulting in a locus that is ~1Mb long (depending on the dataset).
#' @param n_causal The maximum number of potential causal SNPs per locus.
#' This parameter is used somewhat differntly by different fine-mapping tools.
#' See tool-specific functions for details.
#' @param probe_path The location of the file containing translations between probe IDs and gene symbols.
#' Only used for certain eQTL datasets.
#' @param conditioned_snps Which SNPs to conditions on when fine-mapping with \emph{COJO}.
#' @param PAINTOR_QTL_datasets A list of QTL datasets to be used when conducting joint functional fine-mapping with \emph{PAINTOR}.
#' @param PP_threshold The minimum fine-mapped posterior probability for a SNP to be considered part of a Credible Set.
#' For example, \code{PP_threshold=.95} means that all Credible Set SNPs will be 95\% Credible Set SNPs.
#' @param consensus_threshold The minimum number of fine-mapping tools that include a SNP
#'  in their 95\% Credible Sets to consider that it a "Consensus SNP" (\emph{default=2}).
#'
#' @param min_POS Manually set the minimum genomic position for your locus subset.
#' \code{min_POS} can clip the window size set by \code{bp_distance}.
#' Can also be a list of positions (one for each locus) (e.g. \code{min_POS=top_SNPs$min_POS}).
#' @param max_POS Manually set the maximum genomic position for your locus subset.
#' \code{max_POS} can clip the window size set by \code{bp_distance}.
#' Can also be a list of positions (one for each locus) (e.g. \code{max_POS=top_SNPs$max_POS}).
#' @param min_MAF Remove any SNPs with \strong{MAF} < \code{min_MAF}.
#' @param trim_gene_limits If a valid gene symbol is provided to \code{trim_gene_limits},
#' the gene's canonical coordinates are pulled from \code{biomaRt}.
#' This includes introns, exons, and proximal regulatory regions (e.g. promoters).
#' Any SNPs that fall outside these coordinates are remove from downstream fine-mapping.
#' Set \code{trim_gene_limits=F} to not limit by gene coordinates (\emph{default}).
#' @param max_snps The maximum number of SNPs to include in the locus.
#' If the current window size yields > \code{max_snps},
#'  then the outer edges of the of the locus are trimmed until the number of SNPs ≤ \code{max_snps}.
#' @param min_r2 Remove any SNPs are below the LD r2 threshold with the lead SNP within their respective locus.
#' @param LD_block Calculate LD blocks with \emph{plink} and only include the block to which the lead SNP belongs.
#' @param LD_block_size Adjust the granularity of block sizes when \code{LD_block=T}.
#' @param min_Dprime Remove any SNPs are below the LD D' threshold with the lead SNP within their respective locus.
#' This is paramter currently only works when \code{LD_reference!="UKB"}.
#' @param remove_variants A list of variants to remove from the locus subset file.
#' @param remove_correlates A named list, where the names are the RSIDs of SNPs
#' whose LD correlates you wish to remove,
#' and the value is the absolute r2 threshold you wish to filter at for each RSID respectively
#' (e.g. \code{ remove_correlates = c("rs76904798"=.2, "rs10000737"=.8)}).
#' This will also remove the SNPs in \code{remove_correlates} themselves.
#'
#' @param LD_reference Which linkage disequilibrium reference panel do you want to use.
#'  Options include:
#'  \describe{
#'  \item{"UKB"}{A pre-caclulated LD reference matrix from a subset of caucasian British individuals from the UK Biobank. See \href{https://www.biorxiv.org/content/10.1101/807792v2}{Wiessbrod et al. (2019)} for more details.}
#'  \item{"1KGphase1"}{Download a subset of the 1000 Genomes Project Phase 1 vcf and calculate LD on the fly with plink.}
#'  \item{"1KGphase3"}{Download a subset of the 1000 Genomes Project Phase 3 vcf and calculate LD on the fly with plink.}
#'  \item{"<path>/*.vcf" or "<path>/*.vcf.gz"}{Alternatively, users can provide their own custom panel by supplying a list of \emph{.vcf} file path (one per locus) which \pkg{echolocatoR} will use to compute LD (using \emph{plink}).}
#'  }
#' @param superpopulation Subset your LD reference panel by superopulation.
#'  Setting the superpopulation is not currently possible when \code{LD_reference="UKB"}.
#'  \href{https://www.internationalgenome.org/faq/which-populations-are-part-your-study/}{1KGphase1 options} include:
#'  \describe{
#'  \item{"AFR"}{African [descent]}
#'  \item{"AMR"}{Ad-mixed American}
#'  \item{"EAS"}{East Asian}
#'  \item{"EUR"}{European}
#'  \item{"SAS"}{South Asian}
#'  }
#' @param remote_LD When acquiring LD matrixes,
#'  the default is to delete the full vcf or npz files after \pkg{echolocatoR} has extracted the necssary subset.
#'  However, if you wish to keep these full files (which can be quite large) set \code{remote_LD=T}.
#' @param plot_LD Whether to plot a subset of the LD matix.
#'
#'
#' @param plot.types Which kinds of plots to include.
#' Options:
#' \describe{
#' \item{"simple"}{Just plot the following tracks: GWAS, fine-mapping, gene models}
#' \item{"fancy"}{Additionally plot XGR annotation tracks (XGR, Roadmap, Nott).}
#' }
#' @param plot.zoom Zoom into the center of the locus when plotting (without editing the fine-mapping results file).
#' You can provide either:
#' \itemize{
#' \item{The size of your plot window in terms of basepairs (e.g. \code{plot.zoom=50000} for a 50kb window)}.
#' \item{How much you want to zoom in (e.g. \code{plot.zoom="1x"} for the full locus, \code{plot.zoom="2x"} for 2x zoom into the center of the locus, etc.)}.
#' }
#' You can pass a list of window sizes (e.g. \code{c(50000,100000,500000)}) to automatically generate
#' multiple views of each locus.
#' This can even be a mix of different style inputs: e.g. \code{c("1x","4.5x",25000)}.
#' @param plot.Nott_binwidth When including Nott et al. (2019) epigenomic data in the track plots,
#' adjust the bin width of the histograms.
#' @param plot.Nott_bigwig_dir Instead of pulling Nott et al. (2019) epigenomic data
#' from the \emph{UCSC Genome Browser}, use a set of local bigwig files.
#' @param plot.Roadmap Find and plot annotations from Roadmap.
#' @param plot.Roadmap_query Only plot annotations from Roadmap whose metadata contains a string or any items from  a list of strings
#' (e.g. \code{"brain"} or \code{c("brain","liver","monocytes")}).
#'
#' @param verbose Whether \pkg{echolocatoR} should be verbose or silent.
#' @param remove_tmps Whether to remove any temporary files (e.g. FINEMAP output files) after the pipeline is done running.
#' @param server
#' Whether \pkg{echolocatoR} is being run on a computing cluster/server or on a local machine.
#' @param conda_env The name of a conda environment to use.
#'
#' @family MAIN
#' @export
finemap_pipeline <- function(locus,
                             fullSS_path,
                             fullSS_genome_build="hg19",
                             LD_genome_build="hg19",
                             results_dir,
                             dataset_name="dataset_name",
                             dataset_type="GWAS",
                             top_SNPs="auto",
                             force_new_subset=F,
                             force_new_LD=F,
                             force_new_finemap=T,
                             finemap_methods=c("ABF","FINEMAP","SUSIE","POLYFUN_SUSIE"),
                             finemap_args=NULL,
                             bp_distance=500000,
                             n_causal=5,
                             chrom_col="CHR",
                             chrom_type=NULL,
                             position_col="POS",
                             snp_col="SNP",
                             pval_col="P",
                             effect_col="Effect",
                             stderr_col="StdErr",
                             tstat_col="t-stat",
                             locus_col="Locus",
                             freq_col="Freq",
                             MAF_col="MAF",
                             A1_col = "A1",
                             A2_col = "A2",
                             gene_col="Gene",
                             N_cases_col="N_cases",
                             N_controls_col="N_controls",
                             N_cases=NULL,
                             N_controls=NULL,
                             proportion_cases="calculate",
                             sample_size=NULL,

                             LD_reference="1KGphase1",
                             superpopulation="EUR",
                             remote_LD=T,
                             download_method="direct",
                             min_POS=NA,
                             max_POS=NA,
                             min_MAF=NA,
                             trim_gene_limits=F,
                             max_snps=NULL,

                             file_sep="\t",
                             min_r2=0,
                             LD_block=F,
                             LD_block_size=.7,
                             vcf_folder=NULL,
                             # min_Dprime=F,
                             query_by="coordinates",
                             remove_variants=F,
                             remove_correlates=F,
                             probe_path = "./Data/eQTL/gene.ILMN.map",
                             conditioned_snps,
                             plot_LD = F,
                             remove_tmps=T,
                             plot.types=c("simple"),
                             PAINTOR_QTL_datasets=NULL,
                             server=F,
                             PP_threshold=.95,
                             consensus_threshold=2,
                             case_control=T,
                             QTL_prefixes=NULL,
                             fillNA=0,

                             plot.zoom="1x",
                             plot.Nott_epigenome=F,
                             plot.Nott_show_placseq=F,
                             plot.Nott_binwidth=200,
                             plot.Nott_bigwig_dir=NULL,
                             plot.XGR_libnames=NULL,
                             plot.Roadmap=F,
                             plot.Roadmap_query=NULL,

                             conda_env="echoR",
                             nThread=4,
                             verbose=T){

  #### Create paths ####
  subset_path <- get_subset_path(results_dir = results_dir,
                                 dataset_type = dataset_type,
                                 dataset_name = dataset_name,
                                 locus = locus)
  locus_dir <- get_locus_dir(subset_path = subset_path)

  ####  Query ####
  subset_DT <- extract_SNP_subset(locus = locus,
                                  results_dir = results_dir,
                                  locus_dir = locus_dir,
                                  top_SNPs = top_SNPs,
                                  fullSS_path = fullSS_path,
                                  fullSS_genome_build = fullSS_genome_build,
                                  subset_path  =  subset_path,
                                  LD_reference = LD_reference,
                                  force_new_subset = force_new_subset,

                                  chrom_col = chrom_col,
                                  chrom_type = chrom_type,
                                  position_col = position_col,
                                  snp_col = snp_col,
                                  pval_col = pval_col,
                                  effect_col = effect_col,
                                  stderr_col = stderr_col,
                                  locus_col = locus_col,
                                  tstat_col = tstat_col,
                                  MAF_col = MAF_col,
                                  freq_col = freq_col,
                                  A1_col = A1_col,
                                  A2_col = A2_col,
                                  gene_col = gene_col,

                                  N_cases_col = N_cases_col,
                                  N_controls_col = N_controls_col,
                                  N_cases = N_cases,
                                  N_controls = N_controls,
                                  proportion_cases = proportion_cases,
                                  sample_size = sample_size,

                                  bp_distance = bp_distance,
                                  superpopulation = superpopulation,
                                  min_POS = min_POS,
                                  max_POS = max_POS,

                                  file_sep = file_sep,
                                  query_by = query_by,
                                  probe_path = probe_path,
                                  QTL_prefixes=QTL_prefixes,

                                  remove_tmps = remove_tmps,
                                  conda_env = conda_env,
                                  verbose = verbose)
  #### Extract LD ####
  LD_list <- LD.load_or_create(locus_dir=locus_dir,
                               subset_DT=subset_DT,
                               LD_reference=LD_reference,
                               # Optional args (with defaults)
                               force_new_LD=force_new_LD,
                               LD_genome_build=LD_genome_build,
                               superpopulation=superpopulation,
                               remote_LD=remote_LD,
                               download_method=download_method,
                               LD_block=LD_block,
                               LD_block_size=LD_block_size,
                               remove_correlates=remove_correlates,
                               server=server,
                               remove_tmps=remove_tmps,
                               vcf_folder=vcf_folder,
                               conda_env=conda_env,
                               nThread=nThread,
                               verbose=verbose)

  #### Filter SNPs####
  # Remove pre-specified SNPs
  ## Do this step AFTER saving the LD to disk so that it's easier to re-subset
  ## in different ways later without having to redownload LD.


  LD_list <- LD.filter_LD(LD_list=LD_list,
                          remove_correlates=remove_correlates,
                          min_r2=min_r2,
                          verbose=verbose)
  LD_matrix <- LD_list$LD
  subset_DT <- LD_list$DT

  subset_DT <- filter_snps(subset_DT=subset_DT,
                           min_MAF=min_MAF,
                           bp_distance=bp_distance,
                           remove_variants=remove_variants,
                           min_POS=min_POS,
                           max_POS=max_POS,
                           max_snps=max_snps,
                           trim_gene_limits=trim_gene_limits,

                           verbose=verbose)
  # Subset LD and df to only overlapping SNPs
  LD_list <- subset_common_snps(LD_matrix = LD_matrix,
                                finemap_dat = subset_DT,
                                fillNA = fillNA)
  LD_matrix <- LD_list$LD
  subset_DT <- LD_list$DT

  #### Fine-map ####
  finemap_dat <- finemap_handler(locus_dir = locus_dir,
                                 fullSS_path = fullSS_path,
                                 LD_reference = LD_reference,
                                 LD_matrix = LD_matrix,
                                 subset_DT = subset_DT,
                                 sample_size = sample_size,
                                 # General args (with defaults)
                                 finemap_methods = finemap_methods,
                                 finemap_args = finemap_args,
                                 force_new_finemap = force_new_finemap,
                                 dataset_type = dataset_type,
                                 PP_threshold = PP_threshold,
                                 consensus_threshold = consensus_threshold,
                                 case_control = case_control,
                                 n_causal = n_causal,
                                 # Tool-specific args
                                 conditioned_snps = conditioned_snps,
                                 PAINTOR_QTL_datasets = PAINTOR_QTL_datasets,
                                 # Optional args
                                 conda_env = conda_env,
                                 verbose = verbose)
  #### Visualize ####
  locus_plots <- list()
  if("simple" %in% plot.types){
    try({
      TRKS_simple <- PLOT.locus(finemap_dat = finemap_dat,
                                LD_matrix = LD_matrix,
                                LD_reference = LD_reference,
                                locus_dir = locus_dir,
                                dataset_type = dataset_type,
                                method_list = finemap_methods,
                                PP_threshold = PP_threshold,
                                consensus_threshold = consensus_threshold,
                                QTL_prefixes = QTL_prefixes,
                                Nott_epigenome = F,
                                plot_full_window = F,
                                mean.PP = T,
                                XGR_libnames = NULL,
                                plot.zoom = plot.zoom,
                                save_plot = T,
                                show_plot = T,

                                conda_env = conda_env,
                                nThread = nThread,
                                verbose = verbose)
      locus_plots[["simple"]] <- TRKS_simple
    })
  };

  if("fancy" %in% plot.types){
    try({
      TRKS_fancy <- PLOT.locus(finemap_dat = finemap_dat,
                        LD_matrix = LD_matrix,
                        LD_reference = LD_reference,
                        locus_dir = locus_dir,
                        dataset_type = dataset_type,
                        method_list = finemap_methods,
                        PP_threshold = PP_threshold,
                        consensus_threshold = consensus_threshold,
                        QTL_prefixes = QTL_prefixes,
                        plot.zoom = plot.zoom,
                        save_plot = T,
                        show_plot = T,

                        XGR_libnames = plot.XGR_libnames,

                        Roadmap = plot.Roadmap,
                        Roadmap_query = plot.Roadmap_query,

                        Nott_epigenome = plot.Nott_epigenome,
                        Nott_show_placseq = plot.Nott_show_placseq,
                        Nott_binwidth = plot.Nott_binwidth,
                        Nott_bigwig_dir = plot.Nott_bigwig_dir,

                        conda_env = conda_env,
                        nThread = nThread,
                        verbose = verbose)
      locus_plots[["fancy"]] <- TRKS_fancy
    })
  };

  # Plot LD
  ld_plot <- if(plot_LD){
    try({LD.plot_LD(LD_matrix=LD_matrix,
                            subset_DT=finemap_dat,
                            span=10)
    })
  } else {NULL};

  # Cleanup:
  if(remove_tmps){
    roadmap_tbi <- list.files(locus_dir,pattern=".bgz.tbi", recursive = T, full.names = T)
    tmp_files <- file.path(locus_dir,"LD",
                           c("plink.bed",
                             "plink.bim",
                             "plink.fam",
                             "plink.ld",
                             "plink.ld.bin",
                             "plink.log",
                             "plink.nosex",
                             "SNPs.txt") )
    out <- suppressWarnings(file.remove(tmp_files, roadmap_tbi))
  }
  return(list(finemap_dat=finemap_dat,
              locus_plot=locus_plots,
              LD_matrix=LD_matrix,
              LD_plot=ld_plot,
              locus_dir=locus_dir,
              arguments=as.list(match.call(expand.dots=F))
              ))
}






#' Fine-map multiple loci
#'
#' \pkg{echolocatoR} will automatically fine-map each locus.
#' Uses the \code{top_SNPs} data.frame to define locus coordinates.
#'
#' @family MAIN
#' @param loci The list of loci you want to fine-map.
#' If \code{subset_path="auto"} (\emph{default}), a locus subset file name is automatically constructed as:
#' \emph{Data/<dataset_type>/<dataset_name>/<locus>/Multi-finemap/<locus>_<dataset_name>_Multi-finemap.tsv.gz}
#' @inheritParams finemap_pipeline
#' @return A merged data.frame with all fine-mapping results from all loci.
#' @export
finemap_loci <- function(loci,
                         fullSS_path,
                         fullSS_genome_build="hg19",
                         LD_genome_build="hg19",
                         dataset_name="dataset_name",
                         dataset_type="GWAS",
                         force_new_subset=F,
                         force_new_LD=F,
                         force_new_finemap=T,
                         results_dir="./results",
                         top_SNPs="auto",
                         finemap_methods=c("ABF","FINEMAP","SUSIE","POLYFUN_SUSIE"),
                         finemap_args=NULL,
                         bp_distance=500000,
                         n_causal=5,
                         chrom_col="CHR",
                         chrom_type=NULL,
                         position_col="POS",
                         snp_col="SNP",
                         pval_col="P",
                         effect_col="Effect",
                         stderr_col="StdErr",
                         tstat_col = "t-stat",
                         MAF_col="MAF",
                         locus_col="Locus",
                         freq_col="Freq",
                         A1_col = "A1",
                         A2_col = "A2",
                         gene_col="Gene",

                         N_cases_col="N_cases",
                         N_controls_col="N_controls",
                         N_cases=NULL,
                         N_controls=NULL,
                         proportion_cases="calculate",
                         sample_size=NULL,

                         LD_reference="1KGphase1",
                         superpopulation="EUR",
                         download_method="direct",
                         vcf_folder=NULL,
                         remote_LD=T,

                         topVariants=3,
                         min_POS=NA,
                         max_POS=NA,
                         min_MAF=NA,
                         trim_gene_limits=F,
                         max_snps=NULL,
                         file_sep="\t",
                         min_r2=0,
                         LD_block=F, LD_block_size=.7,
                         # min_Dprime=F,
                         query_by="coordinates",
                         remove_variants=F,
                         remove_correlates=F,
                         probe_path = "./Data/eQTL/gene.ILMN.map",
                         conditioned_snps="auto",
                         plot_LD=F,
                         remove_tmps=T,
                         PAINTOR_QTL_datasets=NULL,
                         server=F,
                         PP_threshold=.95,
                         consensus_threshold=2,
                         case_control=T,
                         QTL_prefixes=NULL,

                         plot.types = c("simple"),
                         plot.zoom="1x",
                         plot.Nott_epigenome=F,
                         plot.Nott_show_placseq=F,
                         plot.Nott_binwidth=200,
                         plot.Nott_bigwig_dir=NULL,
                         plot.XGR_libnames=NULL,
                         plot.Roadmap=F,
                         plot.Roadmap_query=NULL,

                         conda_env="echoR",
                         nThread=4,
                         verbose=T){
  CONDA.activate_env(conda_env = conda_env)
  data.table::setDTthreads(threads = nThread);
  conditioned_snps <- snps_to_condition(conditioned_snps, top_SNPs, loci);

  if(detect_genes(loci = loci, verbose = F)){
    printer("Reassigning gene-specific locus names",v=verbose)
    loci <- setNames(paste(unname(loci),names(loci),sep="_"), names(loci))
  }

  FINEMAP_DAT <- lapply(1:length(unique(loci)), function(i){
    start_gene <- Sys.time()
    finemap_dat <- NULL
    locus <- loci[i]
    message("\n)   )  ) ))))))}}}}}}}} {{{{{{{{{(((((( (  (   (");
    message(locus," (",i ," / ",length(loci),")");
    message(")   )  ) ))))))}}}}}}}} {{{{{{{{{(((((( (  (   (");
    try({
      # lead_SNP <- .arg_list_handler(conditioned_snps, i)
      gene_limits <- .arg_list_handler(trim_gene_limits, i)
      conditioned_snp <- .arg_list_handler(conditioned_snps, i)
      min_pos <- .arg_list_handler(min_POS, i)
      max_pos <- .arg_list_handler(max_POS, i)
      LD_ref <- .arg_list_handler(LD_reference, i)

      out_list <- finemap_pipeline(locus=locus,
                                      top_SNPs=top_SNPs,
                                      fullSS_path=fullSS_path,
                                      fullSS_genome_build=fullSS_genome_build,
                                      LD_genome_build=LD_genome_build,
                                      results_dir=results_dir,
                                      finemap_methods=finemap_methods,
                                      finemap_args=finemap_args,
                                      force_new_subset=force_new_subset,
                                      force_new_LD=force_new_LD,
                                      force_new_finemap=force_new_finemap,
                                      dataset_name=dataset_name,
                                      dataset_type=dataset_type,
                                      n_causal=n_causal,
                                      bp_distance=bp_distance,

                                      chrom_col=chrom_col,
                                      chrom_type=chrom_type,
                                      position_col=position_col,
                                      snp_col=snp_col,
                                      pval_col=pval_col,
                                      effect_col=effect_col,
                                      stderr_col=stderr_col,
                                      tstat_col=tstat_col,
                                      locus_col=locus_col,
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

                                      LD_reference=LD_ref,
                                      superpopulation=superpopulation,
                                      download_method=download_method,
                                      min_POS=min_pos,
                                      max_POS=max_pos,
                                      min_MAF=min_MAF,
                                      vcf_folder=vcf_folder,
                                      remote_LD=remote_LD,

                                      trim_gene_limits=gene_limits,
                                      max_snps=max_snps,
                                      file_sep=file_sep,
                                      min_r2=min_r2,
                                      LD_block=LD_block,
                                      LD_block_size=LD_block_size,
                                      # min_Dprime=min_Dprime,
                                      query_by=query_by,
                                      remove_variants=remove_variants,
                                      remove_correlates=remove_correlates,
                                      probe_path=probe_path,
                                      conditioned_snps=conditioned_snps,
                                      plot_LD=plot_LD,
                                      remove_tmps=remove_tmps,
                                      plot.types=plot.types,
                                      PAINTOR_QTL_datasets=PAINTOR_QTL_datasets,
                                      server=server,
                                      PP_threshold=PP_threshold,
                                      consensus_threshold=consensus_threshold,
                                      case_control=case_control,
                                      QTL_prefixes=QTL_prefixes,

                                      plot.zoom=plot.zoom,
                                      plot.Nott_epigenome=plot.Nott_epigenome,
                                      plot.Nott_show_placseq=plot.Nott_show_placseq,
                                      plot.Nott_binwidth=plot.Nott_binwidth,
                                      plot.Nott_bigwig_dir=plot.Nott_bigwig_dir,
                                      plot.XGR_libnames=plot.XGR_libnames,
                                      plot.Roadmap=plot.Roadmap,
                                      plot.Roadmap_query=plot.Roadmap_query,

                                      conda_env=conda_env,
                                      nThread=nThread,
                                      verbose=verbose)
      # Output reminder
      # list(finemap_dat
      #     locus_plot
      #     LD_matrix
      #     LD_plot
      #     locus_dir
      #     arguments)
      finemap_dat <- out_list$finemap_dat
      if(!"Locus" %in% colnames(finemap_dat)){
        finemap_dat <- data.table::data.table(Locus=locus,
                                              finemap_dat)
      }
      cat('  \n')
    }) ## end try()
    end_gene <- Sys.time()
    message("Fine-mapping complete in:")
    print(round(end_gene-start_gene,1))
    return(finemap_dat)
  }) # end for loop
  FINEMAP_DAT <- data.table::rbindlist(FINEMAP_DAT, fill = T)
  try({
    FINEMAP_DAT <- find_consensus_SNPs(finemap_dat = FINEMAP_DAT,
                                       credset_thresh = PP_threshold,
                                       consensus_thresh = consensus_threshold,
                                       verbose = verbose)
  })
  # print(createDT_html( subset(FINEMAP_DAT, Support >0) ))
  return(FINEMAP_DAT)
}

