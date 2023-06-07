## ---- install,  eval = FALSE--------------------------------------------------
#  # install.packages("remotes")
#  remotes::install_github("shuye2009/GenomicPlot",
#                          build_manual = TRUE,
#                          build_vignettes = TRUE)

## ---- global code, eval = TRUE, fig.show = 'hold', fig.keep = 'all', fig.align = 'center', fig.dim = c(7,7), fig.ncol=1, fig.sep="\n\n"----
suppressPackageStartupMessages(library(GenomicPlot, quietly = TRUE))
gtffile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf", 
                       package = "GenomicPlot")
gff <- suppressMessages(RCAS::importGtf(saveObjectAsRds = TRUE, filePath = gtffile))
txdb <- suppressWarnings(makeTxDbFromGRanges(gff))

## ---- metagene code, eval = FALSE, fig.show = 'hold', fig.keep = 'all', fig.align = 'center', fig.dim = c(7,7), fig.ncol=1, fig.sep="\n\n"----
#  
#  gf <- prepare_5parts_genomic_features(txdb,
#    meta = TRUE, nbins = 100, fiveP = -2000, threeP = 1000,
#    longest = TRUE
#  )
#  
#  queryfiles <- system.file("extdata", "treat_chr19.bam", package = "GenomicPlot")
#  names(queryfiles) <- "clip_bam"
#  inputfiles <- system.file("extdata", "input_chr19.bam", package = "GenomicPlot")
#  names(inputfiles) <- "clip_input"
#  
#  bamimportParams <- list(
#    offset = -1, fix_width = 0, fix_point = "start", norm = TRUE,
#    useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#  )
#  
#  plot_5parts_metagene(
#    queryFiles = queryfiles,
#    gFeatures_list = list("metagene" = gf),
#    inputFiles = inputfiles,
#    scale = FALSE,
#    verbose = FALSE,
#    transform = NA,
#    smooth = TRUE,
#    stranded = TRUE,
#    outPrefix = "test_plot_5parts_metagene2",
#    importParams = bamimportParams,
#    heatmap = TRUE,
#    rmOutlier = 0,
#    nc = 2
#  )

## ----metagene2, eval=TRUE, echo = FALSE, fig.cap="Metagene plot with UTR", include = TRUE, out.width="50%"----
knitr::include_graphics("../longtests/test_output/test_plot_5parts_metagene2.pdf")

## ---- region code, eval = FALSE, fig.show = 'hold', fig.keep = 'all', fig.align = 'center', fig.dim = c(7,7)----
#  centerfiles <- system.file("extdata", "test_chip_peak_chr19.narrowPeak",
#                             package = "GenomicPlot")
#  names(centerfiles) <- c("NarrowPeak")
#  queryfiles <- c(
#    system.file("extdata", "chip_treat_chr19.bam", package = "GenomicPlot")
#  )
#  names(queryfiles) <- c("chip_bam")
#  inputfiles <- c(
#    system.file("extdata", "chip_input_chr19.bam", package = "GenomicPlot")
#  )
#  names(inputfiles) <- c("chip_input")
#  
#  chipimportParams <- list(
#    offset = 0, fix_width = 150, fix_point = "start", norm = TRUE,
#    useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#  )
#  
#  plot_region(
#    queryFiles = queryfiles,
#    centerFiles = centerfiles,
#    inputFiles = inputfiles,
#    nbins = 100,
#    heatmap = TRUE,
#    scale = FALSE,
#    regionName = "narrowPeak",
#    importParams = chipimportParams,
#    verbose = FALSE,
#    fiveP = -500,
#    threeP = 500,
#    smooth = TRUE,
#    transform = NA,
#    stranded = TRUE,
#    outPrefix = "test_plot_region",
#    rmOutlier = 0,
#    nc = 2
#  )

## ----region, eval=TRUE, echo = FALSE, fig.cap = "Coverage of custom genomic regions", include = TRUE, out.width="50%"----
knitr::include_graphics("../longtests/test_output/test_plot_region.pdf")

## ---- locus code, eval = FALSE, fig.show = 'hold', fig.keep = 'all', fig.align = 'center', fig.dim = c(7,7)----
#  centerfiles <- c(
#    system.file("extdata", "test_clip_peak_chr19.bed", package = "GenomicPlot"),
#    system.file("extdata", "test_chip_peak_chr19.bed", package = "GenomicPlot")
#  )
#  names(centerfiles) <- c("iCLIPPeak", "SummitPeak")
#  queryfiles <- c(
#    system.file("extdata", "chip_treat_chr19.bam", package = "GenomicPlot")
#  )
#  names(queryfiles) <- c("chip_bam")
#  inputfiles <- c(
#    system.file("extdata", "chip_input_chr19.bam", package = "GenomicPlot")
#  )
#  names(inputfiles) <- c("chip_input")
#  
#  plot_locus(
#    queryFiles = queryfiles,
#    centerFiles = centerfiles,
#    ext = c(-500, 500),
#    hl = c(-100, 100),
#    shade = TRUE,
#    smooth = TRUE,
#    importParams = chipimportParams,
#    binSize = 10,
#    refPoint = "center",
#    Xlab = "Center",
#    inputFiles = inputfiles,
#    stranded = TRUE,
#    scale = FALSE,
#    outPrefix = "test_plot_locus2",
#    verbose = FALSE,
#    transform = NA,
#    rmOutlier = 0,
#    statsMethod = "wilcox.test",
#    heatmap = TRUE,
#    nc = 2
#  )

## ----locus2, eval=TRUE, echo = FALSE, fig.cap = "Coverage around center of features", include = TRUE, out.width="50%"----
knitr::include_graphics("../longtests/test_output/test_plot_locus2.pdf")

## ---- annotation code, eval = FALSE, fig.show = 'hold', fig.keep = 'all', fig.align = 'center', fig.dim = c(7,7)----
#  gtffile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#                         package = "GenomicPlot")
#  
#  centerfile <- system.file("extdata", "test_chip_peak_chr19.bed",
#                            package = "GenomicPlot")
#  names(centerfile) <- c("SummitPeak")
#  
#  bedimportParams <- list(
#    offset = 0, fix_width = 100, fix_point = "center", norm = FALSE,
#    useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19"
#  )
#  
#  pa <- plot_peak_annotation(
#    peakFile = centerfile,
#    gtfFile = gtffile,
#    importParams = bedimportParams,
#    fiveP = -2000,
#    dsTSS = 300,
#    threeP = 1000,
#    simple = FALSE,
#    verbose = TRUE,
#    outPrefix = "test_plot_peak_annotation1"
#  )

## ----annotation1, eval=TRUE, echo = FALSE, include = TRUE, fig.cap = "Distribution of ChIPseq peaks in various types of gene", out.width="50%"----
knitr::include_graphics("../longtests/test_output/test_plot_peak_annotation1.pdf")

## ---- bam correlation, eval = TRUE, fig.show = 'hold', fig.keep = 'last', fig.align = 'center', fig.dim = c(7,7)----
bamQueryFiles <- c(
      system.file("extdata", "chip_input_chr19.bam", package = "GenomicPlot"),
      system.file("extdata", "chip_treat_chr19.bam", package = "GenomicPlot")
    )
    names(bamQueryFiles) <- c("chip_input", "chip_treat")
    
    bamImportParams <- list(
       offset = 0, fix_width = 150, fix_point = "start", norm = FALSE,
       useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19"
    )
   
    plot_bam_correlation(
      bamfiles = bamQueryFiles, binSize = 100000, outPrefix = NULL,
      importParams = bamImportParams, nc = 2, verbose = FALSE
    )

## ---- bed overlap, eval = TRUE, fig.show = 'hold', fig.keep = 'first', fig.align = 'center', fig.dim = c(7,7)----
  queryFiles <- c(
    system.file("extdata", "test_chip_peak_chr19.narrowPeak", package = "GenomicPlot"),
    system.file("extdata", "test_chip_peak_chr19.bed", package = "GenomicPlot"),
    system.file("extdata", "test_clip_peak_chr19.bed", package = "GenomicPlot")
  )
  names(queryFiles) <- c("narrowPeak", "summitPeak", "clipPeak")
 
  bedimportParams <- list(
    offset = 0, fix_width = 100, fix_point = "center", norm = FALSE,
    useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19"
  )
 
  plot_overlap_bed(
    bedList = queryFiles, importParams = bedimportParams, pairOnly = FALSE,
    stranded = FALSE, outPrefix = NULL
  )

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

