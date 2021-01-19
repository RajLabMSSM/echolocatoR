


# &&&&&&&&&&&&&&&&&& SPLICEAI &&&&&&&&&&&&&&&&&&&&
# GitHub: https://github.com/Illumina/SpliceAI

#####------ baseSpaceCLI instructions ------- #####
# To download large files from basepace, you can either:
# 1) use the BaseSpace Download on your local computer (extremely slow and unreliable).
# 2) Use the basespace command line interface (CLI):
# Full Documentation: https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview
# basespaceCLI is now installed on Minerva. To use it:
# 2.0) Request a bunch of cores (to speed up download)

# 2.1) load:
# ml BaseSpaceCLI
# 2.2) authenticate (use bs alias)
# bs auth
# 2.3) List your projects (and their IDs)
# bs list projects
# 2.4) Download entire project (automatically multi-threads)
# bs download project --id 66029966  -o spliceai_download/

# Prepare 'genome_scores_v1.3/spliceai_scores.raw.snv.hg19.vcf.gz' file:
# 1) convert vcf to tsv
# vk vcf2tsv wide --print-header spliceai_scores.raw.snv.hg19.vcf.gz > spliceai_scores.raw.snv.hg19.tsv
# 2) split by |
# awk -F'\t' '{gsub(/[|]/, "\t", $8)} 1' OFS='\t' spliceai_scores.raw.snv.hg19.tsv > splice_split.tsv
# 3) Remove empty cols
#  awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$4,$5,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' splice_split.tsv | tail -n +2 > tmp
# 4) Replace header names
# nano....
# cat tmp_header.txt tmp > spliceai_scores.raw.snv.hg19_mod.tsv
# rm tmp





#' Run pre-trained \emph{SpliceAI} model to get predictions
#'
#' @family SpliceAI
#' @source
#' \href{https://github.com/Illumina/SpliceAI}{GitHub}
#' \href{https://www.sciencedirect.com/science/article/pii/S0092867418316295}{Publication}
SPLICEAI.run <- function(vcf_path="./GWAS_converted.vcf",
                         output_path="spliceai_predictions.vcf",
                         # Ref fasta MUST be unzipped currently
                         reference_fasta="/pd-omics/tools/polyfun/reference_fasta/hg19.fa",
                         gene_annotation="./echolocatoR/tools/spliceAI/grch37.txt",
                         distance=50,
                         mask=0){
  # optional arguments:
  #   -h, --help     show this help message and exit
  #   -I [input]     path to the input VCF file, defaults to standard in
  #   -O [output]    path to the output VCF file, defaults to standard out
  #   -R reference   path to the reference genome fasta file
  #   -A annotation  "grch37" (GENCODE V24lift37 canonical annotation file in
  #                            package), "grch38" (GENCODE V24 canonical annotation file in
  #                                                package), or path to a similar custom gene annotation file
  #   -D [distance]  maximum distance between the variant and gained/lost splice
  #   site, defaults to 50
  #   -M [mask]      mask scores representing annotated acceptor/donor gain and
  #                  unannotated acceptor/donor loss, defaults to 0

  cmd <- paste("spliceai",
               "-I",vcf_path,
               "-O",output_path,
               "-R",reference_fasta,
               "-A",gene_annotation,
               "-D", distance,
               "-M",mask)
  print(cmd)
  system(cmd)
}




#' Query genome-wide \emph{SpliceAI} results file (vcf format)
#'
#' @family SpliceAI
#' @source
#' \href{https://github.com/Illumina/SpliceAI}{GitHub}
#' \href{https://www.sciencedirect.com/science/article/pii/S0092867418316295}{Publication}
SPLICEAI.subset_precomputed_vcf <- function(subset_DT,
                                            precomputed_path="/pd-omics/data/spliceAI/whole_genome_filtered_spliceai_scores.vcf.gz",
                                            subset_vcf="subset.vcf"){
  chrom = gsub("chr","",subset_DT$CHR[1] )
  min_POS = min(subset_DT$POS)
  max_POS = max(subset_DT$POS)

  cmd <- paste("bcftools query",precomputed_path,
               "-i",paste0("'CHROM=\"",chrom,"\" & (POS>=",min_POS," & ", "POS<=",max_POS,")'"),
               # "-i'CHROM="10" & (POS>=93821 & POS<=93825)'",
               "-f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO\\n'",
               "-H",
               ">",subset_vcf) # Include header
  cat(cmd)
  system(cat(cmd))
}




#' Query genome-wide \emph{SpliceAI} results file (tsv format)
#'
#' @family SpliceAI
#' @source
#' \href{https://github.com/Illumina/SpliceAI}{GitHub}
#' \href{https://www.sciencedirect.com/science/article/pii/S0092867418316295}{Publication}
SPLICEAI.subset_precomputed_tsv <- function(subset_DT,
                                            # precomputed_path="/pd-omics/data/spliceAI/spliceai_scores.raw.snv.hg19.tsv.gz",
                                            precomputed_path="/pd-omics/data/spliceAI/whole_genome_filtered_spliceai_scores.tsv.gz",
                                            merge_data=T,
                                            drop_na=T,
                                            filtered=T){
  dat <- TABIX.query(fullSS.gz = precomputed_path,
                     chrom = subset_DT$CHR[1],
                     start_pos = min(subset_DT$POS),
                     end_pos = max(subset_DT$POS))
  if(filtered){
    colnames(dat) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","SYMBOL","STRAND","TYPE","DIST","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL")
    dat <- subset(dat, select=-c(ID,QUAL,FILTER))
    by.x = c("CHR","POS")
    by.y = c("CHROM","POS")

  } else {
    colnames(dat) <- c("CHROM","POS","REF","ALT","MUT","SYMBOL","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL")
    # Have to merge on CHR/POS (RSID not given)
    by.x = c("CHR","POS")# "A1"
    by.y = c("CHROM","POS") # "MUT"
  }
  # summary(dat, na.rm=T)
  # hist(dat$DS_DG, breaks = 100)
  if(merge_data){
    dat_merged <- data.table:::merge.data.table(x = subset_DT,
                                                y = dat,
                                                by.x = by.x,
                                                by.y = by.y,
                                                all.x = !drop_na)
    return(dat_merged)
  } else {return(dat)}
}




#' Make multiple queries to genome-wide \emph{SpliceAI} results file (tsv format)
#'
#' @family SpliceAI
#' @source
#' \href{https://github.com/Illumina/SpliceAI}{GitHub}
#' \href{https://www.sciencedirect.com/science/article/pii/S0092867418316295}{Publication}
#' @examples
#' \dontrun{
#' root.pd <- "/sc/arion/projects/pd-omics"
#' precomputed_path <- file.path(root.pd,"data/spliceAI/spliceai_scores.raw.snv.hg19.tsv.gz")
#' sumstats_paths <- list.files(file.path(root.pd, "/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"), pattern = "*.UKB_LD.Multi-finemap.tsv.gz", recursive = T, full.names = T)
#' DAT <- SPLICEAI.subset_precomputed_tsv_iterate(sumstats_paths=sumstats_paths)
#' }
SPLICEAI.subset_precomputed_tsv_iterate <- function(sumstats_paths,
                                                    # precomputed_path="/pd-omics/data/spliceAI/whole_genome_filtered_spliceai_scores.tsv.gz",
                                                    precomputed_path="/pd-omics/data/spliceAI/spliceai_scores.raw.snv.hg19.tsv.gz",
                                                    nThread=4,
                                                    merge_data=T,
                                                    drop_na=T,
                                                    filtered=F,
                                                    save_path="./spliceAI_subset.tsv.gz"){
  # no_no_loci <- c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1")
  # sumstats_paths <- sumstats_paths[!basename(dirname(dirname(sumstats_paths))) %in% no_no_loci]
  DAT <- parallel::mclapply(sumstats_paths, function(x){
    subset_DT <- data.table::fread(x)
    dat_merged <- SPLICEAI.subset_precomputed_tsv(subset_DT,
                                                  precomputed_path=precomputed_path,
                                                  merge_data=merge_data,
                                                  drop_na=drop_na,
                                                  filtered = filtered)
    if(!"Locus" %in% colnames(dat_merged)){
      locus <- basename(dirname(dirname(x)))
      printer("Adding Locus column:", locus)
      dat_merged <- cbind(Locus=locus, dat_merged)
    }
    return(dat_merged)
  }, mc.cores = nThread) %>% data.table::rbindlist(fill=T)

  if(save_path!=F){
    printer("Saving SpliceAI subset ==>",save_path)
    dir.create(dirname(save_path),showWarnings = F, recursive = T)
    data.table::fwrite(DAT, save_path, nThread = nThread, sep="\t")
  }
  return(DAT)
}




#' Postprocess \emph{SpliceAI} results after querying
#'
#' @family SpliceAI
#' @source
#' \href{https://github.com/Illumina/SpliceAI}{GitHub}
#' \href{https://www.sciencedirect.com/science/article/pii/S0092867418316295}{Publication}
#' @examples
#' \dontrun{
#' root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping"
#' DAT <- data.table::fread(file.path(root, "Data/GWAS/Nalls23andMe_2019/_genome_wide/SpliceAI/spliceAI_Nalls23andMe_2019.hits.csv.gz"))
#' }
SPLICEAI.snp_probs <- function(DAT,
                               save_path=F){
  # merged_DT <- merge_finemapping_results(dataset = "./Data/GWAS/Nalls23andMe_2019",minimum_support = 0)
  # DAT <- data.table::fread("/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/spliceAI_Nalls2019.overlap.tsv.gz")
  DF <- DAT[,c("DS_AG","DS_AL","DS_DG","DS_DL")]
  DAT$max_spliceAI_group <- colnames(DF)[max.col(DF,ties.method="first")]
  DAT$max_spliceAI_prob <- apply(DF, 1, max)
  matchDAT <- DAT %>% dplyr::rename(Risk_allele=A1, Nonrisk_allele=A2,
                                    REF.spliceAI=REF, ALT.spliceAI=ALT) %>%
    subset(Risk_allele==ALT.spliceAI, select = c("Locus","SNP","CHR","POS","Effect","P",
                                                 "leadSNP","Consensus_SNP","Support","mean.PP","Risk_allele","Nonrisk_allele",
                                                 "SYMBOL","REF.spliceAI","ALT.spliceAI","DS_AG","DS_AL","DS_DG","DS_DL","max_spliceAI_group","max_spliceAI_prob")) %>%
    dplyr::mutate(GWAS.sig=P<5e-8) %>%
    subset(max_spliceAI_prob>.1 & ( leadSNP))
  matchDAT
  # data.table::fwrite(matchDAT, "./Data/GWAS/Nalls23andMe_2019/_genome_wide/SpliceAI/spliceAI_raw_subset_matched.tsv", sep="\t")
  if(save_path!=F){
    # save_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide/spliceAI_Nalls2019.matches.tsv"
    dir.create(dirname(save_path), showWarnings = F, recursive = T)
    data.table::fwrite(matchDAT, save_path, sep = "\t")
  }

  return(matchDAT)
}




#' Plot \emph{SpliceAI} predictions
#'
#' @family SpliceAI
#' @source
#' \href{https://github.com/Illumina/SpliceAI}{GitHub}
#' \href{https://www.sciencedirect.com/science/article/pii/S0092867418316295}{Publication}
SPLICEAI.plot <- function(dat_merged){
  library(patchwork)
  plt <-
    ggplot(dat_merged, aes(x=POS, y=-log10(P), color=-log10(P))) +
    geom_point() +
    theme_classic() +
    ggplot(dat_merged, aes(x=POS, y=DS_AG, color=DS_AG)) +
    geom_point() +
    scale_color_viridis_c() +
    theme_classic() +
    ggplot(dat_merged, aes(x=POS, y=DS_AL, color=DS_AL))+
    geom_point() +
    scale_color_viridis_c() +
    theme_classic() +
    ggplot(dat_merged, aes(x=POS, y=DS_DG, color=DS_DG))+
    geom_point() +
    scale_color_viridis_c() +
    theme_classic() +
    ggplot(dat_merged, aes(x=POS, y=DS_DL, color=DS_DL))+
    geom_point() +
    scale_color_viridis_c() +
    theme_classic() +
    patchwork::plot_layout(ncol = 1)
  print(plt)
  return(plt)

  dat_melt <- data.table::melt.data.table(data = dat_merged,
                                          measure.vars = c("DS_AG","DS_AL","DS_DG","DS_DL"),
                                          variable.name = "spliceAI_variable",
                                          value.name = "spliceAI_score",
                                          na.rm = T)
  ggplot(dat_melt, aes(x=POS, y=spliceAI_variable,   height=spliceAI_score)) +
    ggridges::geom_ridgeline()

  ggplot(dat_melt, aes(x=POS, y=-log10(P), color=-log10(P))) +
    geom_point() +
    theme_classic() +
    ggplot(dat_melt, aes(x=POS, y=spliceAI_score, color=spliceAI_score)) +
    geom_point() +
    scale_color_viridis_c() +
    theme_classic() +
    facet_grid(facets = spliceAI_variable~. ) +
    patchwork::plot_layout(ncol = 1)
}

