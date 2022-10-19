## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE, fig.show='hold', fig.keep='all', fig.align='center', fig.dim=c(7,7), fig.ncol=1, fig.sep="\n\n"----
#  library(GenomicPlot, quietly=TRUE)
#  txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
#  gf <- prepare_5parts_genomic_features(txdb, meta=TRUE, nbins=100, fiveP=500, threeP=500,
#                                        longest=TRUE)
#  
#  queryfiles <- system.file("data", "test_clip_chr19.bam", package="GenomicPlotData")
#  names(queryfiles) <- "query"
#  inputfiles <- system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData")
#  names(inputfiles) <- "input"
#  
#  handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE,
#                            useScore=FALSE, outRle=TRUE, useSizeFactor=TRUE, genome="hg19")
#  
#  plot_5parts_metagene(queryfiles, gFeatures=gf, inputfiles=inputfiles, scale=FALSE, verbose=FALSE,
#                       transform=FALSE, smooth=TRUE, stranded=TRUE, outPrefix=NULL,
#                       handleInputParams=handleInputParams, heatmap=TRUE, rmOutlier=TRUE, nc=5)
#  

## ---- eval=FALSE, fig.show='hold', fig.keep='all', fig.align='center', fig.dim=c(7,7)----
#  gf <- prepare_3parts_genomic_features(txdb, meta=FALSE, nbins=100, fiveP=3000, threeP=2000,
#                                        longest=TRUE)
#  
#  queryfiles <- system.file("data", "test_chip_chr19.bam", package="GenomicPlotData")
#  names(queryfiles) <- "query"
#  inputfiles <- system.file("data", "test_chip_input_chr19.bam", package="GenomicPlotData")
#  names(inputfiles) <- "input"
#  
#  handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE,
#                            useScore=FALSE, outRle=TRUE, useSizeFactor=TRUE, genome="hg19")
#  
#  plot_3parts_metagene(queryfiles, gFeatures=gf, inputfiles=inputfiles, scale=FALSE, verbose=FALSE,
#                       transform=FALSE, smooth=TRUE, stranded=TRUE, outPrefix=NULL,
#                       handleInputParams=handleInputParams, heatmap=TRUE, rmOutlier=TRUE, nc=5)
#  

## ---- eval=FALSE, fig.show='hold', fig.keep='last', fig.align='center', fig.dim=c(7,7)----
#  centerfiles <- system.file("data", "test_chip_peak_chr19.narrowPeak", package="GenomicPlotData")
#   names(centerfiles) <- c("narrowpeak")
#   queryfiles <- c(system.file("data", "test_clip_chr19.bam", package="GenomicPlotData"),
#                   system.file("data", "test_chip_chr19.bam", package="GenomicPlotData"))
#   names(queryfiles) <- c("clip_bam", "chip_bam")
#   inputfiles <- c(system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData"),
#                   system.file("data", "test_chip_input_chr19.bam", package="GenomicPlotData"))
#   names(inputfiles) <- c("clip_input", "chip_input")
#  
#   handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, useScore=FALSE,
#    outRle=TRUE, useSizeFactor=TRUE, genome="hg19")
#  
#   plot_reference_region(queryfiles, centerfiles, inputfiles, nbins=100, handleInputParams=handleInputParams,
#    verbose=FALSE, scale=FALSE, heatmap=TRUE, regionName="narrowPeak", fiveP=500, threeP=500, smooth=TRUE,
#    transform=FALSE, stranded=TRUE, outPrefix=NULL, rmOutlier=FALSE, nc=5)
#  
#  

