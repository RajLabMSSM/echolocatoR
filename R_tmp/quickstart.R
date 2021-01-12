
library(echolocatoR)


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
LD_block=F;
LD_block_size=.7;
vcf_folder=NULL;
# min_Dprime=F;
query_by="coordinates";
remove_variants=F;
remove_correlates=F;
probe_path = "./Data/eQTL/gene.ILMN.map";
# conditioned_snps;
plot_LD = F;
remove_tmps=T;
plot.types=c("simple");
PAINTOR_QTL_datasets=NULL;
server=F;
PP_threshold=.95;
consensus_threshold=2;
case_control=T;
QTL_prefixes=NULL;
fillNA=0;

plot.zoom="1x";
plot.Nott_epigenome=F;
plot.Nott_show_placseq=F;
plot.Nott_binwidth=200;
plot.Nott_bigwig_dir=NULL;
plot.XGR_libnames=NULL;
plot.Roadmap=F;
plot.Roadmap_query=NULL;

conda_env="echoR";
nThread=4;
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
plot.types=c("fancy");
## Generate multiple plots of different window sizes;
### all SNPs; 4x zoomed-in; and a 50000bp window
plot.zoom = c("all","4x",50000);
## XGR
# plot.XGR_libnames=c("ENCODE_TFBS_ClusteredV3_CellTypes");
## Roadmap
plot.Roadmap = F;
plot.Roadmap_query = NULL;
# Nott et al. (2019)
plot.Nott_epigenome = T;
plot.Nott_binwidth = 100;
verbose=T;


start_pos=min_POS;
end_pos=max_POS;
chrom=topSNP_sub$CHR[1];


