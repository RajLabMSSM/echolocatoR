
Data_dirs_to_table <- function(Data_dirs, writeCSV=FALSE){
  df <- data.table::rbindlist(Data_dirs, fill=TRUE)
  df <- cbind(Dataset=names(Data_dirs), df)
  echodata::createDT(df)
  if(writeCSV!=FALSE){
    data.table::fwrite(df, writeCSV, quote = FALSE, sep = ",", row.names  = FALSE)
  }
  return(df)
}


list_Data_dirs <- function(writeCSV = "~/Desktop/Fine_Mapping/Data/directories_table.csv"){
  root <- "/sc/arion/projects"
  Data_dirs = list(
    # ++++++++ GWAS SUMMARY STATS ++++++++ #

    # Unpublished results from Stahl/de Witte/Daner?
    "Daner_2020" = list(type="Bipiolar Disorder",
                        topSS="~/Desktop/Fine_Mapping/Data/GWAS/Daner_2020/daner_bip_pgc3.64indexsnps.tsv",
                        fullSS=NA,
                        fullSS.local="~/Desktop/Fine_Mapping/Data/GWAS/Daner_2020/daner_bip_pgc3_nm.tsv.gz",
                        reference=NA),

    # ++++++++ GWAS SUMMARY STATS ++++++++ #
    # Nall et al. (2019) w/ 23andMe
    "Nalls23andMe_2019" = list(type="Parkinson's GWAS",
                               topSS="~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/Nalls2019_TableS2.xlsx",
                               # fullSS=file.path(root,"pd-omics/data/nallsEtAl2019/23andme/PD_all_post30APRIL2015_5.2_extended.txt")),
                               fullSS=file.path(root,"pd-omics/data/nallsEtAl2019/combined_meta/nallsEtAl2019_allSamples_allVariants.mod.txt"),
                               fullSS.local="~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/nallsEtAl2019_allSamples_allVariants.mod.txt.gz",
                               reference="https://www.biorxiv.org/content/10.1101/388165v3"),

    ## IGAP
    "Lambert_2013" = list(type="Alzheimer's GWAS",
                          topSS="~/Desktop/Fine_Mapping/Data/GWAS/Lambert_2013/Lambert_2019_AD_GWAS.xlsx",
                          fullSS=file.path(root,"ad-omics/data/AD_GWAS/Lambert_2013/IGAP_stage_1.1000G.phase3.20130502.tsv"),
                          reference="https://www.nature.com/articles/ng.2802"),

    ## Marioni et al. (2018)
    "Marioni_2018" = list(type="Alzheimer's GWAS",
                          topSS="~/Desktop/Fine_Mapping/Data/GWAS/Marioni_2018/Marioni2018_supplementary_tables.xlsm",
                          fullSS=file.path(root,"ad-omics/data/AD_GWAS/Marioni_2018/Marioni2018.4_UK_Biobank_IGAP_17May2018.1000G.phase3.20130502.tsv"),
                          fullSS.local="~/Desktop/Fine_Mapping/Data/GWAS/Marioni_2018/Marioni2018.4_UK_Biobank_IGAP_17May2018.1000G.phase3.20130502.tsv.gz",
                          reference="https://www.nature.com/articles/s41398-018-0150-6"),

    ## Jansen et al. (2018)
    "Posthuma_2018" = list(type="Alzheimer's GWAS",
                           topSS="~/Desktop/Fine_Mapping/Data/GWAS/Posthuma_2018/Posthuma_2018_Table1.xlsx",
                           fullSS=file.path(root,"ad-omics/data/AD_GWAS/Posthuma_2018/AD_sumstats_Jansenetal_2019sept.txt.gz"),
                           fullSS.local="~/Desktop/Fine_Mapping/Data/GWAS/Posthuma_2018/AD_sumstats_Jansenetal_2019sept.txt.gz",
                           reference="https://www.nature.com/articles/s41588-018-0311-9"),

    ## Kunkle et al. (2018) Alzheimer's GWAS
    "Kunkle_2019" = list(type="Alzheimer's GWAS",
                         topSS="~/Desktop/Fine_Mapping/Data/GWAS/Kunkle_2019/Kunkle2019_supplementary_tables.xlsx",
                         fullSS=file.path(root,"ad-omics/data/AD_GWAS/Kunkle_2019/Kunkle_etal_Stage1_results.txt.gz"),
                         # tr " " "\t" <  Kunkle_etal_Stage1_results.txt > Kunkle_etal_Stage1_results.ts
                         fullSS.local="~/Desktop/Fine_Mapping/Data/GWAS/Kunkle_2019/Kunkle_etal_Stage1_results.tsv",
                         # fullSS_stage2=file.path(root,"ad-omics/data/AD_GWAS/Kunkle_etal_Stage2_results.txt"),
                         reference="https://www.nature.com/articles/s41588-019-0358-2"),

    # ++++++++ eQTL SUMMARY STATS ++++++++ #
    ## MESA eQTLs: African Americans
    "MESA_AFA" = list(type="eQTL",
                      topSS="~/Desktop/Fine_Mapping/Data/eQTL/MESA/AFA/AFA_eQTL_PTK2B.txt",
                      fullSS=file.path(root,"ad-omics/data/MESA/AFA_cis_eqtl_summary_statistics.txt"),
                      fullSS.local="~/Desktop/Fine_Mapping/Data/QTL/MESA/AFA/MESA_AFA.finemap.txt.gz",
                      reference="https://www.nhlbi.nih.gov/science/multi-ethnic-study-atherosclerosis-mesa"),
    ## MESA eQTLs: Caucasians
    "MESA_CAU" = list(type="eQTL",
                      topSS="~/Desktop/Fine_Mapping/Data/eQTL/MESA/CAU/CAU_eQTL_PTK2B.txt",
                      fullSS=file.path(root,"ad-omics/data/MESA/CAU_cis_eqtl_summary_statistics.txt"),
                      fullSS.local="~/Desktop/Fine_Mapping/Data/QTL/MESA/CAU/MESA_CAU.finemap.txt.gz",
                      reference="https://www.nhlbi.nih.gov/science/multi-ethnic-study-atherosclerosis-mesa"),
    ## MESA eQTLs: Hispanics
    "MESA_HIS" = list(type="eQTL",
                      topSS="~/Desktop/Fine_Mapping/Data/eQTL/MESA/HIS/HIS_eQTL_PTK2B.txt",
                      fullSS=file.path(root,"ad-omics/data/MESA/HIS_cis_eqtl_summary_statistics.txt"),
                      fullSS.local="~/Desktop/Fine_Mapping/Data/QTL/MESA/HIS/MESA_HIS.finemap.txt.gz",
                      reference="https://www.nhlbi.nih.gov/science/multi-ethnic-study-atherosclerosis-mesa"),

    ## Fairfax eQTLs: CD14
    "Fairfax_2014_CD14" = list(type="eQTL",
                               topSS=NA,
                               fullSS=file.path(root,"ad-omics/data/fairfax/sumstats/cis.eqtls.fairfax.all.chr.CD14.47231.414.b.qced.f.txt"),
                               fullSS.local="~/Desktop/Fine_Mapping/Data/QTL/Fairfax_2014/CD14/Fairfax_2014_CD14.finemap.txt.gz",
                               reference="https://science.sciencemag.org/content/343/6175/1246949"),
    ## Fairfax eQTLs: IFN
    "Fairfax_2014_IFN" = list(type="eQTL",
                              topSS=NA,
                              fullSS=file.path(root,"ad-omics/data/fairfax/sumstats/cis.eqtls.fairfax.all.chr.IFN.47231.367.b.qced.f.txt"),
                              fullSS.local="~/Desktop/Fine_Mapping/Data/QTL/Fairfax_2014/IFN/Fairfax_2014_IFN.finemap.txt.gz",
                              reference="https://science.sciencemag.org/content/343/6175/1246949"),

    ## Fairfax eQTLs: IFN
    "Fairfax_2014_LPS2" = list(type="eQTL",
                               topSS=NA,
                               fullSS=file.path(root,"ad-omics/data/fairfax/sumstats/cis.eqtls.fairfax.all.chr.LPS2.47231.261.b.qced.f.txt"),
                               fullSS.local="~/Desktop/Fine_Mapping/Data/QTL/Fairfax_2014/LPS2/Fairfax_2014_LPS2.finemap.txt.gz",
                               reference="https://science.sciencemag.org/content/343/6175/1246949"),
    ## Fairfax eQTLs: IFN
    "Fairfax_2014_LPS24" = list(type="eQTL",
                                topSS=NA,
                                fullSS=file.path(root,"ad-omics/data/fairfax/sumstats/cis.eqtls.fairfax.all.chr.LPS24.47231.322.b.qced.f.txt"),
                                fullSS.local="~/Desktop/Fine_Mapping/Data/QTL/Fairfax_2014/LPS24/Fairfax_2014_LPS24.finemap.txt.gz",
                                reference="https://science.sciencemag.org/content/343/6175/1246949"),

    ## Cardiogenics: Macrophages
    "Cardiogenics_macrophages" = list(type="eQTL",
                                      topSS=NA,
                                      fullSS=file.path(root,"ad-omics/data/cardiogenics/Cardiogenics/Macrophages.REPORT.fdr-0.5.tab"),
                                      reference="https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003240"),
    ## Cardiogenics: Monocytes
    "Cardiogenics_monocytes" = list(type="eQTL",
                                    topSS=NA,
                                    fullSS=file.path(root,"ad-omics/data/cardiogenics/Cardiogenics/Monocites.REPORT.fdr-0.5.tab"),
                                    reference="https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003240"),

    ## eqtlpsych (psychgen??? neeed to check)
    "eqtlpsych" = list(type="eQTL",
                       topSS=NA,
                       fullSS=file.path(root,"ad-omics/data/eqtlgen/trans-eQTLs_full_20180905.txt.gz"),
                       reference=NA),

    ## psychENCODE
    "psychENCODE_eQTL" = list(type="eQTL",
                              topSS=NA,
                              fullSS=file.path(root,"ad-omics/data/psychENCODE/eQTL/DER-08a_hg19_eQTL.significant.txt.gz"),
                              reference="http://resource.psychencode.org"),
    ## psychENCODE
    "psychENCODE_cQTL" = list(type="cQTL",
                              topSS=NA,
                              fullSS=file.path(root,"ad-omics/data/psychENCODE/cQTL/DER-09_hg19_cQTL.significant.txt.gz"),
                              reference="http://resource.psychencode.org"),
    ## psychENCODE
    "psychENCODE_isoQTL" = list(type="isoQTL",
                                topSS=NA,
                                fullSS=file.path(root,"ad-omics/data/psychENCODE/isoQTL/DER-10a_hg19_isoQTL.significant.txt.gz"),
                                reference="http://resource.psychencode.org"),
    ## psychENCODE
    "psychENCODE_tQTL" = list(type="tQTL",
                              topSS=NA,
                              fullSS=file.path(root,"ad-omics/data/psychENCODE/tQTL/DER-10c_hg19_tQTL.all.txt.gz"),
                              reference="http://resource.psychencode.org"),
    ## psychENCODE
    "psychENCODE_fQTL" = list(type="fQTL",
                              topSS=NA,
                              fullSS=file.path(root,"ad-omics/data/psychENCODE/fQTL/DER-11_hg19_fQTL.significant.txt.gz"),
                              reference="http://resource.psychencode.org"),
    ## psychENCODE
    "psychENCODE_HiC" = list(type="HiC",
                             topSS=NA,
                             fullSS=file.path(root,"ad-omics/data/psychENCODE/HiC/INT-16_HiC_EP_linkages_cross_assembly.csv.gz"),
                             reference="http://resource.psychencode.org"),

    # (eQTL, sQTL) x (CD4, CD14) x (AFR, ASN, EUR)
    ## Just selecting some
    # ImmVar
    "ImmVar_eQTL_CD4_EUR" = list(type="eQTL",
                                 topSS=NA,
                                 fullSS=file.path(root,"ad-omics/data/immvar/qtls/eur_cd4_chrALL_cis1mb_adj_spearman_pALL.out"),
                                 reference=NA),
    # ImmVar
    "ImmVar_eQTL_CD14_EUR" = list(type="eQTL",
                                  topSS=NA,
                                  fullSS=file.path(root,"ad-omics/data/immvar/qtls/eur_cd14_chrALL_cis1mb_adj_spearman_pALL.out"),
                                  reference=NA),

    # Brain.xQTL.Serve
    "Brain.xQTL.Serve_eQTL" = list(type="eQTL",
                                   topSS=NA,
                                   fullSS=file.path(root,"ad-omics/data/Brain_xQTL_Serve/eQTL/eQTLs_all.txt.gz"),
                                   reference="https://www.ncbi.nlm.nih.gov/pubmed/28869584"),
    # Brain.xQTL.Serve
    "Brain.xQTL.Serve_haQTL" = list(type="haQTL",
                                    topSS=NA,
                                    fullSS=file.path(root,"ad-omics/data/Brain_xQTL_Serve/haQTL/haQTLs_all.txt.gz"),
                                    reference="https://www.ncbi.nlm.nih.gov/pubmed/28869584"),
    # Brain.xQTL.Serve
    "Brain.xQTL.Serve_mQTL" = list(type="mQTL",
                                   topSS=NA,
                                   fullSS=file.path(root,"ad-omics/data/Brain_xQTL_Serve/mQTL/mQTLs_all.txt.gz"),
                                   reference="https://www.ncbi.nlm.nih.gov/pubmed/28869584"),
    # Brain.xQTL.Serve
    "Brain.xQTL.Serve_cell-specificity-eQTL" = list(type="mQTL",
                                                    topSS=NA,
                                                    fullSS=file.path(root,"ad-omics/data/Brain_xQTL_Serve/cell-specificity-eQTL/cell-specificity-eQTLs.tsv.gz"),
                                                    reference="https://www.ncbi.nlm.nih.gov/pubmed/28869584"),

    # GTEx: many different single-tissue eQTLs
    ## NOTE: this is a folder (not the actual file)
    "GTEx_V7" = list(type="eQTL",
                     topSS=NA,
                     fullSS=file.path(root,"ad-omics/data/GTEx_QTL/GTEx_Analysis_v7_eQTL_all_associations"),
                     reference="https://www.nature.com/articles/nature24277"),
    # GTEx: many different single-tissue eQTLs
    ## NOTE: this is a folder (not the actual file)
    "GTEx_V8" = list(type="eQTL",
                     topSS=NA,
                     fullSS=file.path(root,"ad-omics/data/GTEx_QTL/GTEx_Analysis_v8_eQTL_all_associations"),
                     reference="https://www.nature.com/articles/nature24277")
  )
  Data_dirs_table <- Data_dirs_to_table(Data_dirs, writeCSV)
  return(Data_dirs_table)
}


#' Retrieve the location of summary stats files
#'
#' @family directory
Directory_info <- function(dataset_name, variable="fullSS.local"){
  Data_dirs <- list_Data_dirs(writeCSV  = FALSE)
  directory = subset(Data_dirs, Dataset==dataset_name, select=variable) %>% as.character()
  return(directory)
}
