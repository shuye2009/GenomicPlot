## start up ####
library(GenomicPlot)
rm(list=ls())

currentdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currentdir)
options(error=recover)

## testing parallel clusters
cl <- start_parallel(5L)
stop_parallel(cl)

## extract_longest_tx ####
#gtffile <- "C:/GREENBLATT/Rscripts/GenomicPlotData/data/gencode.v19.annotation_chr19.gtf"
#txdb_chr19 <-  makeTxDbFromGFF(gtffile)
#AnnotationDbi::saveDb(txdb_chr19, file="C:/GREENBLATT/Rscripts/GenomicPlotData/data/txdb_chr19.sql")
#txdb_chr19 <- AnnotationDbi::loadDb("C:/GREENBLATT/Rscripts/GenomicPlotData/data/txdb_chr19.sql")

txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))

longestTx <- extract_longest_tx(txdb, plot=FALSE)

output <- get_genomic_feature_coordinates(txdb, featureName="cds", featureSource="gencode", export=FALSE,
                                          longest=TRUE, protein_coding=TRUE)

gf <- prepare_5parts_genomic_features(txdb, meta=TRUE, nbins=100, fiveP=1000, threeP=1000, longest=TRUE)

txdb <- GenomicFeatures::makeTxDbFromEnsembl("Saccharomyces cerevisiae", server="useastdb.ensembl.org")
longestTx <- extract_longest_tx(txdb, plot=TRUE)

output <- get_genomic_feature_coordinates(txdb, featureName="intron", featureSource="Ensembl", export=TRUE,
                                          longest=TRUE, protein_coding=TRUE)
output

# handle_input ####

centerfiles <- c(system.file("data", "test_B.bed", package="GenomicPlotData"),
system.file("data", "test_C.bed", package="GenomicPlotData"))
names(centerfiles) <- c("TestB", "TestC")

queryfiles <- system.file("data", "test_clip_chr19.bam", package="GenomicPlotData")
names(queryfiles) <- "query"

handleInputParams <- list(fix_width=0, fix_point="center", useScore=FALSE, outRle=FALSE,
                          CLIP_reads=FALSE, norm=FALSE, useSizeFactor=FALSE, genome="hg19")

output <- handle_input(centerfiles, handleInputParams)
lapply(output, function(x){print(names(x)); print(length(x$query)); print(x$size);
   print(x$type); print(x$weight)})

inputfiles <- system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData")
names(inputfiles) <- "input"

handleInputParams <- list(fix_width=0, fix_point="center", useScore=FALSE, outRle=TRUE,
CLIP_reads=FALSE, norm=TRUE, useSizeFactor=TRUE, genome="hg19")

output <- handle_input(c(inputfiles, queryfiles), handleInputParams)

lapply(output, function(x){print(names(x)); print(length(x$query)); print(x$size);
print(x$type); print(x$weight)})

## plot_start_end_feature ####
txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
queryfiles <- system.file("data", "test_clip_chr19.bam", package="GenomicPlotData")
names(queryfiles) <- "query"
inputfiles <- system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData")
names(inputfiles) <- "input"
ext <- c(-500, 200, -200, 500)
hl <- c(-50, 50, -50, 50)
op <- "plot_start_end_feature1"
handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, useScore=FALSE, outRle=TRUE,
                          useSizeFactor=TRUE, genome="hg19")
plot_start_end_feature(queryfiles=queryfiles, inputfiles=inputfiles, txdb=txdb, featureName="intron", binsize=10,
                       handleInputParams=handleInputParams, longest=TRUE, ext=ext, hl=hl, randomize=TRUE, insert=100,
                       stranded=TRUE, scale=FALSE, smooth=TRUE, outPrefix=op, nc=5)

txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
queryfiles <- system.file("data", "test_chip_chr19.bam", package="GenomicPlotData")
names(queryfiles) <- "query"
inputfiles <- system.file("data", "test_chip_input_chr19.bam", package="GenomicPlotData")
names(inputfiles) <- "input"
ext <- c(-500, 200, -200, 500)
hl <- c(-50, 50, -50, 50)
op <- "plot_start_end_feature2"
handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, useScore=FALSE, outRle=TRUE,
                          useSizeFactor=TRUE, genome="hg19")
plot_start_end_feature(queryfiles=queryfiles, inputfiles=inputfiles, txdb=txdb, featureName="intron", binsize=10,
                       handleInputParams=handleInputParams, longest=TRUE, ext=ext, hl=hl, randomize=TRUE, insert=100,
                       stranded=TRUE, scale=FALSE, smooth=TRUE, outPrefix=op, nc=5)

txdb <- system.file("data", "test_chip_peak_chr19.narrowPeak", package="GenomicPlotData")
names(txdb) <- "narrowPeak"
queryfiles <- system.file("data", "test_chip_chr19.bam", package="GenomicPlotData")
names(queryfiles) <- "query"
inputfiles <- system.file("data", "test_chip_input_chr19.bam", package="GenomicPlotData")
names(inputfiles) <- "input"
ext <- c(-500, 200, -200, 500)
hl <- c(-50, 50, -50, 50)
op <- "plot_start_end_feature3"
plot_start_end_feature(queryfiles=queryfiles, inputfiles=inputfiles, txdb=txdb, featureName="custom",
                       insert=100, binsize=10, longest=TRUE, ext=ext, hl=hl, randomize=TRUE, stranded=TRUE,
                       scale=FALSE, smooth=TRUE, handleInputParams=handleInputParams,
                       outPrefix=op, nc=5)

txdb <- system.file("data", "test_chip_peak_chr19.narrowPeak", package="GenomicPlotData")
names(txdb) <- "narrowPeak"
queryfiles <- system.file("data", "test_clip_chr19.bam", package="GenomicPlotData")
names(queryfiles) <- "query"
inputfiles <- system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData")
names(inputfiles) <- "input"
ext <- c(-500, 200, -200, 500)
hl <- c(-50, 50, -50, 50)
op <- "plot_start_end_feature4"
plot_start_end_feature(queryfiles=queryfiles, inputfiles=inputfiles, txdb=txdb, featureName="custom",
                       insert=100, binsize=10, longest=TRUE, ext=ext, hl=hl, randomize=TRUE, stranded=TRUE,
                       scale=FALSE, smooth=TRUE, handleInputParams=handleInputParams,
                       outPrefix=op, nc=5)


# plot_start_end_reference_region ####
centerfiles <- system.file("data", "test_chip_peak_chr19.narrowPeak", package="GenomicPlotData")
 names(centerfiles) <- c("narrowPeak")
 queryfiles <- system.file("data", "test_chip_chr19.bam", package="GenomicPlotData")
 names(queryfiles) <- "query"
 inputfiles <- system.file("data", "test_chip_input_chr19.bam", package="GenomicPlotData")
 names(inputfiles) <- "input"
 ext <- c(-500, 200, -200, 500)
 hl <- c(-50, 50, -50, 50)
 op <- "plot_start_end_reference_region"

 handleInputParams <- list(CLIP_reads=FALSE, fix_width=100, fix_point="start", norm=TRUE,
  useScore=FALSE, outRle=TRUE, useSizeFactor=TRUE, genome="hg19")

 plot_start_end_reference_region(queryfiles=queryfiles, inputfiles=inputfiles, transform=FALSE,
  centerfiles=centerfiles, handleInputParams=handleInputParams, binsize=10, insert=100,
  ext=ext, hl=hl, stranded=TRUE, scale=FALSE, smooth=TRUE, rmOutlier=FALSE, outPrefix=op, nc=5)

# plot_3parts_metagene ####

txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
gf <- prepare_3parts_genomic_features(txdb, meta=TRUE, nbins=100, fiveP=500, threeP=500, longest=TRUE,
protein_coding=TRUE)
gf2 <- prepare_3parts_genomic_features(txdb, meta=FALSE, nbins=100, fiveP=2000, threeP=2000, longest=TRUE,
                                      protein_coding=TRUE)


queryfiles <- system.file("data", "test_clip_chr19.bam", package="GenomicPlotData")
names(queryfiles) <- "query"
inputfiles <- system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData")
names(inputfiles) <- "input"

handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, useScore=FALSE,
outRle=TRUE, useSizeFactor=TRUE, genome="hg19")

op <- "test_plot_3parts_metagene1"
plot_3parts_metagene(queryfiles, gFeatures=gf, inputfiles=inputfiles, scale=FALSE, verbose=FALSE, transform=FALSE,
                     smooth=TRUE, stranded=TRUE, outPrefix=op, handleInputParams=handleInputParams, heatmap=TRUE,
                     rmOutlier=TRUE, nc=5)
op <- "test_plot_3parts_metagene2"
plot_3parts_metagene(queryfiles, gFeatures=gf2, inputfiles=inputfiles, scale=FALSE, verbose=FALSE, transform=FALSE,
                     smooth=TRUE, stranded=TRUE, outPrefix=op, handleInputParams=handleInputParams, heatmap=TRUE,
                     rmOutlier=TRUE, nc=5)

queryfiles <- c(system.file("data", "test_chip_peak_chr19.narrowPeak", package="GenomicPlotData"),
                system.file("data", "test_chip_peak_chr19.bed", package="GenomicPlotData"),
                system.file("data", "test_clip_peak_chr19.bed", package="GenomicPlotData"))
names(queryfiles) <- c("narrowPeak", "summitPeak", "iCLIPPeak")
op <- "test_plot_3parts_metagene3"

handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="center", norm=FALSE, useScore=FALSE,
outRle=TRUE, useSizeFactor=FALSE, genome="hg19")

plot_3parts_metagene(queryfiles, gFeatures=gf2, inputfiles=NULL, handleInputParams=handleInputParams,
verbose=FALSE, smooth=TRUE, scale=FALSE, stranded=TRUE, outPrefix=op, transform=FALSE, heatmap=TRUE,
rmOutlier=TRUE, nc=5)

queryfiles <- system.file("data", "test_chip_chr19.bam", package="GenomicPlotData")
names(queryfiles) <- "chip"
inputfiles <- system.file("data", "test_chip_input_chr19.bam", package="GenomicPlotData")
names(inputfiles) <- "chip_input"
op <- "test_plot_3parts_metagene4"

handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, useScore=FALSE,
outRle=TRUE, useSizeFactor=TRUE, genome="hg19")

plot_3parts_metagene(queryfiles, gFeatures=gf2, inputfiles=inputfiles, handleInputParams=handleInputParams,
verbose=FALSE, smooth=TRUE, scale=FALSE, stranded=TRUE, outPrefix=op, transform=FALSE, heatmap=TRUE,
rmOutlier=TRUE, nc=5)


# plot_5parts_metagene ####
library(AnnotationHub)
ah = AnnotationHub()
hubInfo <- mcols(ah)
gencode_query <- query(ah, c("GENCODE", "GRCh37", "TxDb", "v31"))
hubInfo[names(gencode_query),]
txdb_v31 <- ah[[names(gencode_query)]]

txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
gf <- prepare_5parts_genomic_features(txdb, meta=TRUE, nbins=100, fiveP=500, threeP=500, longest=TRUE)
gf2 <- prepare_5parts_genomic_features(txdb, meta=FALSE, nbins=100, fiveP=2000, threeP=2000, longest=TRUE)


queryfiles <- system.file("data", "test_clip_chr19.bam", package="GenomicPlotData")
names(queryfiles) <- "query"
inputfiles <- system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData")
names(inputfiles) <- "input"

handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, useScore=FALSE,
                          outRle=TRUE, useSizeFactor=TRUE, genome="hg19")

op <- "test_plot_5parts_metagene1"
plot_5parts_metagene(queryfiles, gFeatures=gf, inputfiles=inputfiles, scale=FALSE, verbose=FALSE, transform=FALSE,
                     smooth=TRUE, stranded=TRUE, outPrefix=op, handleInputParams=handleInputParams, heatmap=TRUE,
                     rmOutlier=TRUE, nc=5)
op <- "test_plot_5parts_metagene2"
plot_5parts_metagene(queryfiles, gFeatures=gf2, inputfiles=inputfiles, scale=FALSE, verbose=FALSE, transform=FALSE,
                     smooth=TRUE, stranded=TRUE, outPrefix=op, handleInputParams=handleInputParams, heatmap=TRUE,
                     rmOutlier=TRUE, nc=5)

queryfiles <- c(system.file("data", "test_chip_peak_chr19.narrowPeak", package="GenomicPlotData"),
                system.file("data", "test_chip_peak_chr19.bed", package="GenomicPlotData"),
                system.file("data", "test_clip_peak_chr19.bed", package="GenomicPlotData"))
names(queryfiles) <- c("narrowPeak", "summitPeak", "iCLIPPeak")
op <- "test_plot_5parts_metagene3"

handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="center", norm=FALSE, useScore=FALSE,
                          outRle=TRUE, useSizeFactor=FALSE, genome="hg19")

plot_5parts_metagene(queryfiles, gFeatures=gf2, inputfiles=NULL, handleInputParams=handleInputParams,
                     verbose=FALSE, smooth=TRUE, scale=FALSE, stranded=TRUE, outPrefix=op, transform=FALSE, heatmap=TRUE,
                     rmOutlier=TRUE, nc=5)

queryfiles <- system.file("data", "test_chip_chr19.bam", package="GenomicPlotData")
names(queryfiles) <- "query"
inputfiles <- system.file("data", "test_chip_input_chr19.bam", package="GenomicPlotData")
names(inputfiles) <- "input"

handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, useScore=FALSE,
                          outRle=TRUE, useSizeFactor=TRUE, genome="hg19")

op <- "test_plot_5parts_metagene4"
plot_5parts_metagene(queryfiles, gFeatures=gf2, inputfiles=inputfiles, scale=FALSE, verbose=FALSE, transform=FALSE,
                     smooth=TRUE, stranded=TRUE, outPrefix=op, handleInputParams=handleInputParams, heatmap=TRUE,
                     rmOutlier=TRUE, nc=5)

# plot_reference_region ####
centerfiles <- system.file("data", "test_chip_peak_chr19.narrowPeak", package="GenomicPlotData")
 names(centerfiles) <- c("narrowpeak")
 queryfiles <- c(system.file("data", "test_clip_chr19.bam", package="GenomicPlotData"),
                 system.file("data", "test_chip_chr19.bam", package="GenomicPlotData"))
 names(queryfiles) <- c("clip_bam", "chip_bam")
 inputfiles <- c(system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData"),
                 system.file("data", "test_chip_input_chr19.bam", package="GenomicPlotData"))
 names(inputfiles) <- c("clip_input", "chip_input")
 op <- "test_plot_reference_region"

 handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, useScore=FALSE,
  outRle=TRUE, useSizeFactor=TRUE, genome="hg19")

 plot_reference_region(queryfiles, centerfiles, inputfiles, nbins=100, handleInputParams=handleInputParams,
  verbose=FALSE, scale=FALSE, heatmap=TRUE, regionName="narrowPeak", fiveP=500, threeP=500, smooth=TRUE,
  transform=FALSE, stranded=TRUE, outPrefix=op, rmOutlier=FALSE, nc=5)

 # plot_reference_locus ####
 centerfiles <- c(system.file("data", "test_clip_peak_chr19.bed", package="GenomicPlotData"),
                  system.file("data", "test_chip_peak_chr19.bed", package="GenomicPlotData"))
 names(centerfiles) <- c("clip_peak", "chip_peak")
 queryfiles <- c(system.file("data", "test_clip_chr19.bam", package="GenomicPlotData"),
                 system.file("data", "test_chip_chr19.bam", package="GenomicPlotData"))
 names(queryfiles) <- c("clip_bam", "chip_bam")
 inputfiles <- c(system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData"),
                 system.file("data", "test_chip_input_chr19.bam", package="GenomicPlotData"))
 names(inputfiles) <- c("clip_input", "chip_input")

 op <- "test_plot_reference_locus1"

 handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, useScore=FALSE,
                           outRle=TRUE, useSizeFactor=TRUE, genome="hg19")

 plot_reference_locus(queryfiles, centerfiles, ext=c(-500,500), hl=c(-100,100), shade=TRUE, smooth=TRUE,
                        handleInputParams=handleInputParams, binsize=10, refPoint="center", Xlab="Center",
                        inputfiles=inputfiles, stranded=TRUE, scale=FALSE, outPrefix=op, verbose=FALSE,
                        transform=FALSE, rmOutlier=FALSE, stats.method="wilcox.test", heatmap=TRUE, nc=5)

 centerfiles <- c(system.file("data", "test_chip_peak_chr19.bed", package="GenomicPlotData"))
 names(centerfiles) <- c("chip_peak")
 queryfiles <- c(system.file("data", "test_clip_peak_chr19.bed", package="GenomicPlotData"),
                 system.file("data", "test_chip_peak_chr19.narrowPeak", package="GenomicPlotData"))
 names(queryfiles) <- c("clip_peak", "narrowPeak")

 op <- "test_plot_reference_locus2"

 handleInputParams <- list(CLIP_reads=FALSE, fix_width=100, fix_point="center", norm=FALSE, useScore=FALSE,
                           outRle=TRUE, useSizeFactor=FALSE, genome="hg19")

 plot_reference_locus(queryfiles, centerfiles, ext=c(-3000,3000), hl=c(-100,100), shade=TRUE, smooth=TRUE,
                      handleInputParams=handleInputParams, binsize=10, refPoint="center", Xlab="Center",
                      inputfiles=NULL, stranded=TRUE, scale=FALSE, outPrefix=op, verbose=FALSE,
                      transform=FALSE, rmOutlier=FALSE, stats.method="wilcox.test", heatmap=TRUE, nc=5)

# plot_reference_locus_with_random ####
txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
 centerfiles <- c(system.file("data", "test_clip_peak_chr19.bed", package="GenomicPlotData"))
 names(centerfiles) <- c("clip_peak")
 queryfiles <- c(system.file("data", "test_clip_chr19.bam", package="GenomicPlotData"))
 names(queryfiles) <- c("clip_bam")
 inputfiles <- c(system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData"))
 names(inputfiles) <- c("clip_input")

op <- "test_plot_reference_locus_with_random1"
# note: this example may take a long time to finish (in hours)

handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, useScore=FALSE,
 outRle=TRUE, useSizeFactor=TRUE, genome="hg19")

plot_reference_locus_with_random(queryfiles, centerfiles, txdb, ext=c(-500,500), hl=c(-100,100), shade=TRUE,
handleInputParams=handleInputParams, binsize=10, refPoint="center", Xlab="Center", smooth=TRUE,
inputfiles=inputfiles, stranded=TRUE, scale=FALSE, outPrefix=op, verbose=FALSE, transform=FALSE,
rmOutlier=FALSE, n_random=1, stats.method="wilcox.test", nc=5)

# plot_bam_correlation ####
queryfiles <- c(system.file("data", "test_chip_chr19.bam", package="GenomicPlotData"),
                system.file("data", "test_clip_chr19.bam", package="GenomicPlotData"))
names(queryfiles) <- c("chip_query", "clip_query")
inputfiles <- c(system.file("data", "test_chip_input_chr19.bam", package="GenomicPlotData"),
                system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData"))
names(inputfiles) <- c("chip_input", "clip_input")

op="test_bam_correlation"
handleInputParams <- list(CLIP_reads=FALSE, fix_width=0, fix_point="start", norm=FALSE, useScore=FALSE,
outRle=FALSE, useSizeFactor=FALSE, genome="hg19")
plot_bam_correlation(bamfiles=c(queryfiles, inputfiles), binsize=10000, outPrefix=op,
                     handleInputParams=handleInputParams, nc=4)

# plot_peak_annotation ####
library(BiocFileCache)
bfc <- BiocFileCache()
url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
gz <- bfcrpath(bfc, url)
gtffile <- gsub(".gz", "", gz)
if(!file.exists(gtffile)) gtffile <- R.utils::gunzip(gz, remove=FALSE)

gtffile <- system.file("data", "gencode.v19.annotation_chr19.gtf", package="GenomicPlotData")

centerfile <- system.file("data", "test_chip_peak.narrowPeak", package="GenomicPlotData")
names(centerfile) <- c("narrowPeak")
op <- "test_plot_peak_annotation"
handleBedparams <- list(fix_width=0, fix_point="center", useScore=FALSE, outRle=FALSE, CLIP_reads=FALSE,
norm=FALSE, useSizeFactor=FALSE, genome="hg19")

pa <- plot_peak_annotation(peakfile=centerfile, gtfFile=gtffile, handleInputParams=handleBedparams, fiveP=2000,
threeP=1000, simple=FALSE, RNA=FALSE, verbose=TRUE, outPrefix=op)
pb <- plot_peak_annotation(peakfile=centerfile, gtfFile=gtffile, handleInputParams=handleBedparams, fiveP=2000,
                           threeP=1000, simple=TRUE, RNA=FALSE, verbose=TRUE, outPrefix=op)

# plot_overlap_bed ####
queryfiles <- c(system.file("data", "test_chip_peak_chr19.narrowPeak", package="GenomicPlotData"),
system.file("data", "test_chip_peak_chr19.bed", package="GenomicPlotData"),
system.file("data", "test_clip_peak_chr19.bed", package="GenomicPlotData"))
names(queryfiles) <- c("narrowPeak", "summitPeak", "iCLIPPeak")
op <- "test_plot_overlap_bed"
handleBedParams <- list(fix_width=100, fix_point="center", useScore=FALSE, outRle=FALSE,
CLIP_reads=FALSE, norm=FALSE, useSizeFactor=FALSE, genome="hg19")

plot_overlap_bed(bedList=queryfiles, outPrefix=op, handleInputParams=handleBedParams, pairOnly=FALSE, stranded=FALSE)

# rm_outlier ####
source("https://www.rickwash.com/papers/cscw08-appendix/powerlaw.R")
fullmatrix <- matrix(rpowerlaw(100, alpha=1.25, xmin=1), ncol=10)
normmatrix <- matrix(rnorm(1000), ncol=10)
fullmatrix <- rbind(fullmatrix, normmatrix)

rm_outlier(fullmatrix, verbose=TRUE)

