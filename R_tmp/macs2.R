############ ############ ############
############ MACS2 ############
############ ############ ############

# must have macs2 installed first
## GitHub repo: https://github.com/taoliu/MACS/
## Conda source: https://anaconda.org/bioconda/macs2

MACS2.download_full_bigwig <- function(bigWig.paths='Nott_2019',
                                       output.path = "/Volumes/Scizor/Nott_2019/",
                                       nThreads=4,
                                       force_overwrite = F){
  dir.create(output.path, showWarnings = F, recursive = T)
  if(all(bigWig.paths=='Nott_2019')){
    bigWig.meta <- readxl::read_excel("./echolocatoR/annotations/Nott_2019/Nott_2019.snEpigenomics.xlsx")
    biWig.links <- bigWig.meta$data_link
  }

  bigWig.paths <- lapply(biWig.links, function(URL){
    printer("MACS2:: Downloading", basename(URL))
    out.path <- axel(input.url = URL, output.path = output.path,
                     background = F, force_overwrite = force_overwrite, 
                     nThreads = nThreads)
    return(out.path)
  }) %>% unlist()
    return(bigWig.paths)
}
# bigWig.paths <- MACS2.download_full_bigwig(output.path = "../Nott_2019")
# bigWig.paths <- list.files("/pd-omics/brian/Nott_2019", full.names = T)

MACS2.bigwig_to_bedgraph <- function(bigWig.paths,
                                     converter_path="./echolocatoR/tools/UCSC_utilities",
                                     force_new_download=F,
                                     OS="MAC"){
  # Download executable
  destfile <- file.path(converter_path,OS,"bigWigToBedGraph.dms")
  
  if(!file.exists(destfile) | force_new_download){ 
      download.file("http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigWigToBedGraph",
                    destfile = destfile) 
    # Change permissions
    system(paste("chmod u+x",destfile))
  }

  # Convert to bedGraph
  bedGraph.paths <- lapply(bigWig.paths, function(bw){
    printer("MACS2:: Processing",basename(bw),"...")
    bg <- gsub(".bigWig$",".bedGraph",bw)
    if(file.exists(bg)){
      printer("MACS2:: Existing bedGraph file detected. Skipping...")
    } else {
      cmd <- paste(destfile, bw, bg)
      print(cmd)
      system(cmd)
    }
    print("--------")
    return(bg)
  }) %>% unlist()
  printer("MACS2:: Returning",length(bedGraph.paths),"bedGraph paths.")
  return(bedGraph.paths)
}
# bedGraph.paths <- MACS2.bigwig_to_bedgraph(bigWig.paths)


MACS2.call_peaks <- function(bedGraph.paths,
                             broad_peaks=F,
                             out.dir=NULL){
  # Activate conda env
  # reticulate::conda_install(envname = "echolocatoR",
  #                           conda =  "/opt/anaconda3/bin/conda",
  #                           packages = c("MACS2"))
  # reticulate::use_condaenv(condaenv = "echolocatoR",
  #                          conda =  "/opt/anaconda3/bin/conda")

  # Set up output dir
  if(is.null(out.dir)){
    out.dir <- file.path(dirname(bedGraph.paths[1]),"peaks")
  }
  dir.create(out.dir, showWarnings = F, recursive = T)

  # Call dem peaks
  peak.paths <- lapply(bedGraph.paths, function(bg){
    broad_narrow <- ifelse(broad_peaks, "_c2.0_C1.00_l200_g30_G800_broad.bed12","_c5.0_l200_g30_peaks.narrowPeak")
    pk.path <- file.path(out.dir, gsub(".bedGraph", broad_narrow, basename(bg)))

    if(file.exists(pk.path)){
      printer("MACS2:: Pre-existing peak file detected", basename(pk.path))
      printer("Skipping...")
    } else {
      printer("MACS2:: Calling peaks for",basename(bg))
      cmd <- paste("macs2",
                   ifelse(broad_peaks, "bdgbroadcall","bdgpeakcall"),
                   "--ifile",bg,
                   "--o-prefix",gsub(".bedGraph","",basename(bg)),
                   "--outdir",out.dir)
      print(cmd)
      system(cmd)
    }
    return(pk.path)
  }) %>% unlist()
  return(peak.paths)
}



MACS2.import_peaks <- function(peaks.paths,
                               as_granges=T){
  # Column format found here: https://genome.ucsc.edu/FAQ/FAQformat.html#format12
        #   chrom - Name of the chromosome (or contig, scaffold, etc.).
        #   chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
        #   chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
        #   name - Name given to a region (preferably unique). Use "." if no name is assigned.
        #   score - Indicates how dark the peak will be displayed in the browser (0-1000). If all scores were "'0"' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
        # strand - +/- to denote strand or orientation (whenever applicable). Use "." if no orientation is assigned.
        # signalValue - Measurement of overall (usually, average) enrichment for the region.
        # pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
        # qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
        # peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.


  # chrom - Name of the chromosome (or contig, scaffold, etc.).
  # chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
  # chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. If all scores were "0" when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
  # name - Name given to a region (preferably unique). Use "." if no name is assigned.
  # score - Indicates how dark the peak will be displayed in the browser (0-1000).
  # strand - +/- to denote strand or orientation (whenever applicable). Use "." if no orientation is assigned.
  # signalValue - Measurement of overall (usually, average) enrichment for the region.
  # pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
  # qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
  col.names = c("chrom","chromStart","chromEnd",
                "name","score","Strand",
                "signalValue","pValue","qValue","peak")
  PEAKS <- parallel::mclapply(peaks.paths, function(pk){
    # Select import type
    if(endsWith(pk,suffix = ".narrowPeak")){
      # NarrowPeaks
      peaks <- data.table::fread(pk, skip = 1,  col.names =  col.names)
    } else {
      # BroadPeaks are in a different format (not provided in macs2 documentation for some reason...)
      # https://www.biostars.org/p/174551/
      peaks <- data.table::fread(pk, skip = 1,
                                       col.names = c("chrom","chromStart","chromEnd",
                                                     "name","score","Strand",
                                                     "thickStart","thickEnd","itemRgb",
                                                     "blockCount","blockSizes","blockStarts",
                                                     "signalValue","pValue","qValue") ) %>%
        dplyr::mutate(peak=NA) %>%
        dplyr::select(col.names)
    }
    peaks$peak_type  <- gsub("1|Region|Peak","",strsplit(peaks$name[1],"\\.ucsc_")[[1]][2])
    printer("MACS2:: Importing peaks - ",basename(pk))
    if(as_granges){
      peaks <- GenomicRanges::makeGRangesFromDataFrame(df=peaks, keep.extra.columns = T, ignore.strand = T,
                                              seqnames.field = "chrom",
                                              start.field = "chromStart",end.field = "chromEnd")
    }
    return(peaks)
  }, mc.cores = 4)
  return(PEAKS)
}
