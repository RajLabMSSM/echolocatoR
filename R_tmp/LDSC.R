# devtools::install_github("MichelNivard/GenomicSEM")
# https://github.com/MichelNivard/GenomicSEM/wiki
require(GenomicSEM)

setwd("~/Desktop/LDSC_test")
# root <- "~/Desktop/LDSC_test"

# MUNGE
n_cases <- 37.7 + 18.6
n_controls <- 1400
n_total <- n_cases + n_controls

## DOWNLOAD GWAS 
# gwas <- data.table::fread(,
#                           data.table=F) 
# write.table(gwas, file = "LRRK2_finemap.txt", 
#             sep = "\t", quote = F, row.names = F, col.names = T)

# Download H3 snp reference files from:
# https://utexas.app.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v/file/289808334208
files = file.path("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/",
                   c("LRRK2/Multi-finemap/Multi-finemap_results.txt",
                     "SNCA/Multi-finemap/Multi-finemap_results.txt"))

# MUNGE
genes <- basename(dirname(dirname(files)))
gz_files <- paste0(genes, ".sumstats.gz")
suppressWarnings(file.remove(gz_files))
GenomicSEM::munge(files = files, 
                  hm3 = "w_hm3.noMHC.snplist",
                  trait.names = genes,
                  N = rep(n_total, length(files)), 
                  info.filter = 0, 
                  maf.filter = 0 )
SS <- data.table::fread("LRRK2.sumstats.gz") 
 

# DOWNLOAD "LD score package"
# Your LD folder should be in a subdirectory of 
## the folder where your munged summary statistics are stored.
ref_url <- "https://data.broadinstitute.org/alkesgroup/LDSCORE"
ref_file <- "eur_w_ld_chr.tar.bz2" 
# download.file(url = file.path(ref_url, ref_file),
#               destfile = ref_file)
# untar(ref_dir, list = T, exdir = dirname(ref_dir))

##### LDSC
# 1. A vector of file names/paths to files which point to the munged sumstats.
traits <- c("LRRK2.sumstats.gz")
# 2. A vector of samples prevalences. Again, if the trait is continuous 
# the values should equal NA. The sample prevalence is equal to number of 
# cases/total sample size (possible range: 0-1).
sample.prev <- c(n_cases/n_total)
# A vectors of population prevalences of equal length. 
# If the trait is continuous, the values should equal NA.
population.prev <-c(0.005) # From Nalls et al. 2019 (bioRxiv)
# LD file paths 
ld <- gsub(".tar.bz2","",ref_file)
wld <- ld
trait.names <- genes
# run LDSC
LDSCoutput <- GenomicSEM::ldsc(traits, 
                               sample.prev, 
                               population.prev, 
                               ld, wld, 
                               trait.names)

##optional command to save the ldsc output in case you want to use it in a later R session. 
save(LDSCoutput, file="Pfactor.RData")
 