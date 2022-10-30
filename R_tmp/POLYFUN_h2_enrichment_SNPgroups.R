#' Run heritability enrichment tests across SNP groups
#' @source
#' https://www.nature.com/articles/s41588-020-00735-5
#' @keywords internal
#' @family polyfun
#' @examples
#' \dontrun{
#' root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
#' # IMPORTANT! For this to make sense, you need to merge the full data ("merged_DT" only includes Support>0 and leadSNPs)
#' merged_dat <- merge_finemapping_results(dataset = dirname(root), LD_reference = "UKB", minimum_support = 0)
#' merged_dat <- echodata::find_consensus_snps_no_polyfun(merged_dat)
#'
#' RES <- POLYFUN_h2_enrichment_SNPgroups(merged_dat=merged_dat, ldsc_dir=file.path(root,"PolyFun/output"),  save_enrich=file.path(root,"PolyFun/Nalls23andMe_2019.h2_enrich.snp_groups.csv.gz"))
#' }
POLYFUN_h2_enrichment_SNPgroups <- function(merged_dat,
                                            chrom="*",
                                            ldsc_dir="/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output",
                                            ldsc_suffix="*.snpvar_constrained.gz",
                                            save_enrich=FALSE,
                                            nThread=1){
    # Gather your heritability
    ldsc.files <- list.files(ldsc_dir, pattern = ldsc_suffix, full.names = TRUE) |>
        grep(pattern = paste0(".",chrom,"."), value = TRUE)
    h2_DF <- rbind_filelist(ldsc.files)
    
    h2_merged <- data.table::merge.data.table(merged_dat,
                                              subset(h2_DF, select=c(SNP,SNPVAR)),
                                              all.x = TRUE,
                                              by=c("SNP"))
    if(!"Support_noPF" %in% colnames(merged_dat)){
        merged_dat <- echodata::find_consensus_snps_no_polyfun(merged_dat)
    }
    
    # Iterate over loci
    RES <-parallel::mclapply(unique(merged_dat$Locus), function(locus){
        print(locus)
        dat <- subset(merged_dat, Locus==locus)
        h2_df <- subset(h2_DF, SNP %in% unique(dat$SNP))
        
        # Random subset  of the same size as the Consenus SNPs
        size <- dplyr::n_distinct(subset(dat, Consensus_SNP)$SNP)
        random <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                        target_SNPs= sample(dat$SNP, size = size) )
        # All SNPs
        GWAS.all <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                          target_SNPs=dat$SNP )
        # GWAS nominally sig hits
        GWAS.nom.sig <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                              target_SNPs=subset(dat, P<.05)$SNP )
        # GWAS sig hits
        GWAS.sig <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                          target_SNPs=subset(dat, P<5e-8)$SNP)
        # Lead GWAS
        GWAS.lead <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                           target_SNPs=subset(dat, leadSNP)$SNP)
        # Credible Set
        UCS <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                     target_SNPs = subset(dat, Support>0)$SNP)
        # # PAINTOR CS
        # PAINTOR_credset <- POLYFUN_h2_enrichment(h2_df=h2_df,
        #                                          target_SNPs = subset(dat, PAINTOR_CS>0)$SNP)
        ABF.credset <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                             target_SNPs = subset(dat, ABF.CS>0)$SNP)
        
        FINEMAP.credset <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                                 target_SNPs = subset(dat, FINEMAP.CS>0)$SNP)
        
        SUSIE.credset <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                               target_SNPs = subset(dat, SUSIE.CS>0)$SNP)
        
        POLYFUN_SUSIE.credset <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                                       target_SNPs = subset(dat, POLYFUN_SUSIE.CS>0)$SNP)
        # Support levels
        ## Support==0
        support0 <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                          target_SNPs = subset(dat, Support==0)$SNP)
        support1 <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                          target_SNPs = subset(dat, Support==1)$SNP)
        support2 <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                          target_SNPs = subset(dat, Support==2)$SNP)
        support3 <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                          target_SNPs = subset(dat, Support==3)$SNP)
        support4 <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                          target_SNPs = subset(dat, Support==4)$SNP)
        
        # Consenus SNPs
        Finemap.consensus <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                                   target_SNPs = subset(dat, Consensus_SNP)$SNP)
        
        # Consensus SNPs (no PolyFun)
        Finemap.consensus_noPF <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                                        target_SNPs = subset(dat, Consensus_SNP_noPF)$SNP)
        # UCS SNPs (no PolyFun)
        UCS_noPF <- POLYFUN_h2_enrichment(h2_df=h2_df,
                                          target_SNPs = subset(dat, Support_noPF>0)$SNP)
        
        res <- data.frame(SNP_group=c("Random",
                                      "All",
                                      "GWAS nom. sig.",
                                      "GWAS sig.",
                                      "GWAS lead",
                                      "ABF CS",
                                      "SUSIE CS",
                                      "POLYFUN-SUSIE CS",
                                      "FINEMAP CS",
                                      "UCS (-PolyFun)",
                                      "UCS",
                                      "Support==0",
                                      "Support==1",
                                      "Support==2",
                                      "Support==3",
                                      "Support==4",
                                      "Consensus (-PolyFun)",
                                      "Consensus"),
                          h2.enrichment=c(random,
                                          GWAS.all,
                                          GWAS.nom.sig,
                                          GWAS.sig,
                                          GWAS.lead,
                                          ABF.credset,
                                          SUSIE.credset,
                                          POLYFUN_SUSIE.credset,
                                          FINEMAP.credset,
                                          UCS_noPF,
                                          UCS,
                                          support0,
                                          support1,
                                          support2,
                                          support3,
                                          support4,
                                          Finemap.consensus_noPF,
                                          Finemap.consensus))
        res <- cbind(Locus=locus, res)
        return(res)
    }, mc.cores = nThread) |> data.table::rbindlist(fill = TRUE)
    
    if(save_enrich!=FALSE){
        messager("POLFUN:: Saving enrichment results ==>",save_enrich)
        dir.create(dirname(save_enrich), showWarnings = FALSE, recursive = TRUE)
        data.table::fwrite(RES, save_enrich, nThread=nThread)
    }
    return(RES)
}
