

# Pre-computed LD (European samples - UK10K sequence data), MAF, TSS distance, p-value files for two example traits (Crohn’s Disease from the IBD Consortium and Height from the GIANT consortium) and annotation files for 1005 GENCODE, ENCODE and Roadmap Epigenomics an- notations can be downloaded from http://www.ebi.ac.uk/birney-srv/GARFIELD/package/garfield- data.tar.gz. Note the data is 5.9Gb in compressed format and needs to be uncompressed prior to analysis (83Gb). Variant genomic position (build 37) is used as an identifier in all data files.
#
# Documentation:
# https://bioconductor.org/packages/release/bioc/manuals/garfield/man/garfield.pdf
#
# R function:
# https://rdrr.io/bioc/garfield/man/garfield.run.html
#
# Tutorial:
# https://www.ebi.ac.uk/birney-srv/GARFIELD/
#
# Chimera path:*******
# /sc/orga/projects/ad-omics/wongg05/garfields


# Load in Chimera
# system("ml garfield")
# Load in R
library(garfield) # BiocManager::install("garfield")

# garfield.run()



# INPUT FILES
dir.create("./Data/GWAS/Nalls23andMe_2019/_genome_wide/GARFIELD", showWarnings = FALSE, recursive = T)

garfield.run(out.file = "GARFIELD.",
             data.dir = "./Data/GWAS/Nalls23andMe_2019/_genome_wide/GARFIELD",
             trait="trait",
             run.option = "complete", # prep => perm
             chrs = c(1:22,'X'),
             exclude = c(895, 975, 976, 977, 978, 979, 98),
             nperm = 100000,
             thresh = c(0.1, 0.01, 0.001, 1e-04, 1e-05, 1e-06, 1e-07, 1e-08),
             pt_thresh = c(1e-05, 1e-06, 1e-07, 1e-08),
             maf.bins = 5,
             tags.bins = 5,
             tss.bins = 5,
             prep.file = "",
             optim_mode = 1,
             minit = 100,
             thresh_perm = 1e-04
             )

garfield.run("tmp", data.dir=system.file("extdata",package = "garfield"),
             run.option = "perm", nperm = 1000, thresh = c(0.001, 1e-04, 1e-05),
             pt_thresh = c(1e-04, 1e-05), maf.bins = 2, tags.bins = 3, tss.bins = 3,
             prep.file = "tmp.prep", optim_mode = TRUE, minit = 100, thresh_perm = 0.05)

if (file.exists("tmp.perm")){
  perm = read.table("tmp.perm", header=TRUE)
  head(perm)
} else { print("Error: tmp.perm does not exist!") }
