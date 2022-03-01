
library(echolocatoR)

subset_DT <- finemap_DT <- echolocatoR::BST1
LD_matrix <- echolocatoR::LD_matrix
results_dir="~/Desktop/results";
fullSS_genome_build="hg19";
LD_genome_build="hg19";
dataset_name="dataset_name";
dataset_type="GWAS";
top_SNPs="auto";
force_new_subset=F;
force_new_LD=F;
force_new_finemap=T;
bp_distance=500000;
n_causal=5;
chrom_col="CHR";
chrom_type=NULL;
position_col="POS";
snp_col="SNP";
pval_col="P";
effect_col="Effect";
stderr_col="StdErr";
tstat_col="t-stat";
locus_col="Locus";
freq_col="Freq";
MAF_col="MAF";
A1_col = "A1";
A2_col = "A2";
gene_col="Gene";
N_cases_col="N_cases";
N_controls_col="N_controls";
N_cases=NULL;
N_controls=NULL;
proportion_cases="calculate";
sample_size=NULL;

LD_reference="1KGphase1";
superpopulation="EUR";
remote_LD=T;
download_method="direct";
min_POS=NA;
max_POS=NA;
min_MAF=NA;
trim_gene_limits=F;
max_snps=NULL;

file_sep="\t";
min_r2=0;
leadSNP_LD_block=F;
# min_Dprime=F;
query_by="coordinates";
remove_variants=F;
remove_correlates=F;
probe_path = "./Data/eQTL/gene.ILMN.map";
# conditioned_snps;
remove_tmps=T;
plot_types=c("simple");
PAINTOR_QTL_datasets=NULL;
server=F;
PP_threshold=.95;
consensus_threshold=2;
case_control=T;
qtl_prefixes=NULL;
fillNA=0;

zoom="1x";
nott_epigenome=F;
nott_show_placseq=F;
nott_binwidth=200;
nott_bigwig_dir=NULL;
xgr_libnames=NULL;
roadmap=F;
roadmap_query=NULL;

conda_env="echoR";
nThread=1;
verbose=T;
chrom=NULL;






root.dir <- "~/Desktop"
data("Nalls_top_SNPs");
top_SNPs <- import_topSNPs(
  # topSS = "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/Nalls2019_TableS2.xlsx";
  topSS = Nalls_top_SNPs,
  chrom_col = "CHR",
  position_col = "BP",
  snp_col="SNP",
  pval_col="P, all studies",
  effect_col="Beta, all studies",
  gene_col="Nearest Gene",
  locus_col = "Nearest Gene",
  grouping_vars = c("Locus"),
  remove_variants = "rs34637584")
head(top_SNPs)

results_dir = "~/Desktop/results";
loci = c("BST1","MEX3C");# top_SNPs$Locus;
dataset_name = "Nalls23andMe_2019";
dataset_type = "GWAS";
force_new_subset = F;
force_new_LD = F;
force_new_finemap = F;
remove_tmps = T;

# SUMMARY STATS ARGUMENTS
fullSS_path = fullSS_path;
query_by ="tabix";
chrom_col = "CHR"; position_col = "POS"; snp_col = "RSID";
pval_col = "p"; effect_col = "beta"; stderr_col = "se";
freq_col = "freq"; MAF_col = "calculate";
A1_col = "A1";
A2_col = "A2";

# FILTERING ARGUMENTS
bp_distance = 500000*2;
min_MAF = 0.001;
trim_gene_limits = F;

# FINE-MAPPING ARGUMENTS
finemap_methods = c("ABF","FINEMAP","SUSIE","POLYFUN_SUSIE");
n_causal = 5;
PP_threshold = .95;

# LD ARGUMENTS
LD_reference = "UKB";#"1KGphase1";
superpopulation = "EUR";
download_method = "axel";

# PLOT ARGUMENTS
## general
plot_types=c("fancy");
## Generate multiple plots of different window sizes;
### all SNPs; 4x zoomed-in; and a 50000bp window
zoom = c("all","4x",50000);
## XGR
# xgr_libnames=c("ENCODE_TFBS_ClusteredV3_CellTypes");
## Roadmap
roadmap = F;
roadmap_query = NULL;
# Nott et al. (2019)
nott_epigenome = T;
nott_binwidth = 100;
verbose=T;


start_pos=min_POS;
end_pos=max_POS;
chrom="4";

subset_path <- get_subset_path(results_dir = results_dir,
                               dataset_type = dataset_type,
                               dataset_name = dataset_name,
                               locus = locus)
locus_dir <- get_locus_dir(subset_path = subset_path)



