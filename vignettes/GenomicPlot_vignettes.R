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

## ---- eval=FALSE, fig.show='hold', fig.keep='all', fig.align='center', fig.dim=c(7,7)----
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
#                             outRle=TRUE, useSizeFactor=TRUE, genome="hg19")
#  
#   plot_reference_region(queryfiles, centerfiles, inputfiles, nbins=100, heatmap=TRUE, scale=FALSE,
#                         regionName="narrowPeak", handleInputParams=handleInputParams, verbose=FALSE,
#                         fiveP=500, threeP=500, smooth=TRUE, transform=FALSE, stranded=TRUE, outPrefix=NULL,
#                         rmOutlier=FALSE, nc=5)

## ---- eval=FALSE, fig.show='hold', fig.keep='all', fig.align='center', fig.dim=c(7,7)----
#  txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
#  queryfiles <- system.file("data", "test_clip_chr19.bam", package="GenomicPlotData")
#  names(queryfiles) <- "query"
#  inputfiles <- system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData")
#  names(inputfiles) <- "input"
#  ext <- c(-500, 200, -200, 500)
#  hl <- c(-50, 50, -50, 50)
#  
#  handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, useScore=FALSE,
#                            outRle=TRUE, useSizeFactor=TRUE, genome="hg19")
#  plot_start_end_feature(queryfiles=queryfiles, inputfiles=inputfiles, txdb=txdb, featureName="intron",
#                         binsize=10,handleInputParams=handleInputParams, longest=TRUE, ext=ext, hl=hl,
#                         randomize=TRUE, insert=100, stranded=TRUE, scale=FALSE, smooth=TRUE,
#                         outPrefix=NULL, nc=5)

## ---- eval=FALSE, fig.show='hold', fig.keep='all', fig.align='center', fig.dim=c(7,7)----
#   centerfiles <- c(system.file("data", "test_clip_peak_chr19.bed", package="GenomicPlotData"),
#                    system.file("data", "test_chip_peak_chr19.bed", package="GenomicPlotData"))
#   names(centerfiles) <- c("clip_peak", "chip_peak")
#   queryfiles <- c(system.file("data", "test_clip_chr19.bam", package="GenomicPlotData"),
#                   system.file("data", "test_chip_chr19.bam", package="GenomicPlotData"))
#   names(queryfiles) <- c("clip_bam", "chip_bam")
#   inputfiles <- c(system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData"),
#                   system.file("data", "test_chip_input_chr19.bam", package="GenomicPlotData"))
#   names(inputfiles) <- c("clip_input", "chip_input")
#  
#   handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, useScore=FALSE,
#                             outRle=TRUE, useSizeFactor=TRUE, genome="hg19")
#  
#   plot_reference_locus(queryfiles, centerfiles, ext=c(-500,500), hl=c(-100,100), shade=TRUE, smooth=TRUE,
#                       handleInputParams=handleInputParams, binsize=10, refPoint="center", Xlab="Center",
#                       inputfiles=inputfiles, stranded=TRUE, scale=FALSE, outPrefix=NULL, verbose=FALSE,
#                       transform=FALSE, rmOutlier=FALSE, stats.method="wilcox.test", heatmap=TRUE, nc=5)

## ---- eval=FALSE, fig.show='hold', fig.keep='all', fig.align='center', fig.dim=c(7,7)----
#  txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
#   centerfiles <- c(system.file("data", "test_clip_peak_chr19.bed", package="GenomicPlotData"))
#   names(centerfiles) <- c("clip_peak")
#   queryfiles <- c(system.file("data", "test_clip_chr19.bam", package="GenomicPlotData"))
#   names(queryfiles) <- c("clip_bam")
#   inputfiles <- c(system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData"))
#   names(inputfiles) <- c("clip_input")
#  
#  handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, useScore=FALSE,
#                            outRle=TRUE, useSizeFactor=TRUE, genome="hg19")
#  
#  plot_reference_locus_with_random(queryfiles, centerfiles, txdb, ext=c(-500,500), hl=c(-100,100),
#                                   shade=TRUE, handleInputParams=handleInputParams, binsize=10,
#                                   refPoint="center", Xlab="Center", smooth=TRUE,inputfiles=inputfiles,
#                                   stranded=TRUE, scale=FALSE, outPrefix=NULL, verbose=FALSE,
#                                   transform=FALSE, rmOutlier=FALSE, n_random=1, stats.method="wilcox.test",
#                                   nc=5)

## ---- eval=FALSE, fig.show='hold', fig.keep='all', fig.align='center', fig.dim=c(7,7)----
#  suppressPackageStartupMessages(library(GenomicPlot, quietly=TRUE))
#  gtffile <- system.file("data", "gencode.v19.annotation_chr19.gtf", package="GenomicPlotData")
#  
#  centerfile <- system.file("data", "test_chip_peak.narrowPeak", package="GenomicPlotData")
#  names(centerfile) <- c("narrowPeak")
#  
#  handleBedparams <- list(fix_width=0, fix_point="center", useScore=FALSE, outRle=FALSE, CLIP_reads=FALSE,
#                          norm=FALSE, useSizeFactor=FALSE, genome="hg19")
#  
#  pa <- plot_peak_annotation(peakfile=centerfile, gtfFile=gtffile, handleInputParams=handleBedparams,
#                             fiveP=2000, threeP=1000, simple=FALSE, RNA=FALSE, verbose=TRUE,
#                             outPrefix=NULL)

