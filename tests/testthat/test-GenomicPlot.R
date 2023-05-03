
library(GenomicPlot)
library(magick)

Sys.setenv("R_TESTS" = "")

pdf_to_png <- function(op){
   img <- image_read_pdf(paste0(op,".pdf"))
   for(i in seq_along(img)){image_write(img[i], path=paste0(op,"_",i, ".png"), format="png")}
}

outdir <- "../test_output"
setwd(outdir)

gtffile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf", package="GenomicPlot")
txdb <- AnnotationDbi::loadDb(system.file("extdata", "txdb_chr19.sql", package="GenomicPlot"))

bedQueryFiles <- c(system.file("extdata", "test_chip_peak_chr19.narrowPeak", package="GenomicPlot"),
                   system.file("extdata", "test_chip_peak_chr19.bed", package="GenomicPlot"),
                   system.file("extdata", "test_clip_peak_chr19.bed", package="GenomicPlot"))
names(bedQueryFiles) <- c("Narrow", "Summit", "iCLIP")

bedHandleInputParams <- list(offset=0, fix_width=100, fix_point="center", norm=FALSE, 
                             useScore=FALSE, outRle=TRUE, useSizeFactor=FALSE, genome="hg19")

bamQueryFiles <- system.file("extdata", "treat_chr19.bam", package="GenomicPlot")
names(bamQueryFiles) <- "query"
bamInputFiles <- system.file("extdata", "input_chr19.bam", package="GenomicPlot")
names(bamInputFiles) <- "input"

bamHandleInputParams <- list(offset=-1, fix_width=0, fix_point="start", norm=TRUE, 
                             useScore=FALSE, outRle=TRUE, useSizeFactor=FALSE, genome="hg19")

chipQureyFiles <- system.file("extdata", "chip_treat_chr19.bam", package="GenomicPlot")
names(chipQureyFiles) <- "chip_query"
chipInputFiles <- system.file("extdata", "chip_input_chr19.bam", package="GenomicPlot")
names(chipInputFiles) <- "chip_input"

chipHandleInputParams <- list(offset=0, fix_width=150, fix_point="start", norm=TRUE, 
                              useScore=FALSE, outRle=TRUE, useSizeFactor=TRUE, genome="hg19")


test_that("testing plot_5parts_metagene", {
   gf <- prepare_5parts_genomic_features(txdb, meta=TRUE, nbins=100, fiveP=-2000, threeP=1000, longest=TRUE)
   op <- "test_plot_5parts_metagene1"
   plot_5parts_metagene(queryFiles=bedQueryFiles,
                        gFeatures_list=list("metagene"=gf),
                        inputFiles=NULL,
                        handleInputParams=bedHandleInputParams,
                        verbose=FALSE,
                        smooth=TRUE,
                        scale=FALSE,
                        stranded=TRUE,
                        outPrefix=op,
                        transform=NA,
                        heatmap=TRUE,
                        rmOutlier=0,
                        heatRange=c(0,1),
                        nc=2)
   pdf_to_png(op)
   
   op <- "test_plot_5parts_metagene2"
   plot_5parts_metagene(queryFiles=bamQueryFiles, 
                        gFeatures_list=list("metagene"=gf), 
                        inputFiles=bamInputFiles, 
                        scale=FALSE, 
                        verbose=FALSE, 
                        transform=NA, 
                        smooth=TRUE, 
                        stranded=TRUE, 
                        outPrefix=op, 
                        handleInputParams=bamHandleInputParams, 
                        heatmap=TRUE, 
                        rmOutlier=0, 
                        nc=2)
   pdf_to_png(op)
})
   
test_that("testing plot_3parts_metagene", {
   gf_gene <- prepare_3parts_genomic_features(txdb, meta=FALSE, nbins=100, fiveP=-3000, threeP=2000, 
                                              longest=TRUE)
   op <- "test_plot_3parts_metagene"
   plot_3parts_metagene(queryFiles=chipQureyFiles, 
                        gFeatures=gf_gene, 
                        inputFiles=chipInputFiles, 
                        scale=FALSE, 
                        verbose=FALSE, 
                        transform=NA, 
                        smooth=TRUE, 
                        stranded=TRUE, 
                        outPrefix=op, 
                        handleInputParams=chipHandleInputParams, 
                        heatmap=TRUE, 
                        rmOutlier=0, 
                        nc=2)
   pdf_to_png(op)
})

test_that("testing plot_locus", {
   op <- "test_plot_locus1"
   plot_locus(queryFiles=bedQueryFiles[c(1,3)],
              centerFiles=bedQueryFiles[2],
              ext=c(-1000, 1000),
              hl=c(-100, 100),
              inputFiles=NULL,
              handleInputParams=bedHandleInputParams,
              shade=TRUE,
              binSize=10,
              refPoint="center",
              Xlab="Summit",
              verbose=FALSE,
              smooth=TRUE,
              scale=FALSE,
              stranded=TRUE,
              outPrefix=op,
              transform=NA,
              heatmap=TRUE,
              heatRange=c(0,1),
              rmOutlier=0,
              nc=2)
   pdf_to_png(op)
   
   queryfiles <- c(bamQueryFiles, chipQureyFiles)
   inputfiles <- c(bamInputFiles, chipInputFiles)
   
   handleInputParams <- list(offset=0, fix_width=150, fix_point="start", norm=TRUE, 
                             useScore=FALSE, outRle=TRUE, useSizeFactor=TRUE, genome="hg19")
   
   op <- "test_plot_locus2"
   plot_locus(queryFiles=queryfiles, 
              centerFiles=bedQueryFiles[2:3], 
              ext=c(-500,500), 
              hl=c(-100,100), 
              shade=TRUE, 
              smooth=TRUE, 
              handleInputParams=handleInputParams, 
              binSize=10, 
              refPoint="center", 
              Xlab="Center", 
              inputFiles=inputfiles, 
              stranded=TRUE, 
              scale=FALSE, 
              outPrefix=op, 
              verbose=FALSE, 
              transform=NA, 
              rmOutlier=0, 
              statsMethod="wilcox.test", 
              heatmap=TRUE, 
              nc=2)
   pdf_to_png(op)
})
test_that("testing plot_peak_annotation", {
   op <- "test_plot_peak_annotation1"
   peakHandleInputParams <- list(offset=0, fix_width=100, fix_point="center", norm=FALSE, 
                                 useScore=FALSE, outRle=FALSE, useSizeFactor=FALSE, genome="hg19")
   
   plot_peak_annotation(peakFile=bedQueryFiles[2],
                        gtfFile=gtffile,
                        handleInputParams=peakHandleInputParams,
                        fiveP=-2000,
                        dsTSS=200,
                        threeP=1000,
                        outPrefix=op,
                        verbose=FALSE)
   pdf_to_png(op)
   
   op <- "test_plot_peak_annotation2"
   peakHandleInputParams <- list(offset=0, fix_width=21, fix_point="center", norm=FALSE, 
                                 useScore=FALSE, outRle=FALSE, useSizeFactor=FALSE, genome="hg19")
   
   plot_peak_annotation(peakFile=bedQueryFiles[3],
                        gtfFile=gtffile,
                        handleInputParams=peakHandleInputParams,
                        fiveP=-1000,
                        dsTSS=0,
                        threeP=2000,
                        outPrefix=op,
                        verbose=FALSE)
   pdf_to_png(op)
})

test_that("testing plot_region", {
   queryfiles <- c(bamQueryFiles, chipQureyFiles)
   inputfiles <- c(bamInputFiles, chipInputFiles)
                   
   handleInputParams <- list(offset=0, fix_width=150, fix_point="start", norm=TRUE, 
                             useScore=FALSE, outRle=TRUE, useSizeFactor=FALSE, genome="hg19")
   
   op <- "test_plot_region"
   plot_region(queryFiles=queryfiles, 
               centerFiles=bedQueryFiles[1], 
               inputFiles=inputfiles, 
               nbins=100, 
               heatmap=TRUE, 
               scale=FALSE,  
               regionName="narrowPeak", 
               handleInputParams=handleInputParams, 
               verbose=FALSE, 
               fiveP=-500, 
               threeP=500, 
               smooth=TRUE, 
               transform=NA, 
               stranded=TRUE, 
               outPrefix=op, 
               rmOutlier=0, 
               nc=2)
   pdf_to_png(op)
})

test_that("testing plot_start_end", {
   op <- "test_plot_start_end"
   plot_start_end(queryFiles=bamQueryFiles, 
                  inputFiles=bamInputFiles, 
                  txdb=txdb, 
                  centerFiles="intron", 
                  binSize=10,
                  handleInputParams=bamHandleInputParams, 
                  ext=c(-500, 200, -200, 500), 
                  hl=c(-100, 100, -100, 100), 
                  insert=100, 
                  stranded=TRUE, 
                  scale=FALSE, 
                  smooth=TRUE, 
                  outPrefix=op, 
                  nc=2)
   pdf_to_png(op)
})

test_that("testing plot_start_end_with_random", {
   op <- "test_plot_start_end_with_random"
   plot_start_end_with_random(queryFiles=bamQueryFiles, 
                  inputFiles=bamInputFiles, 
                  txdb=txdb, 
                  centerFile="intron", 
                  binSize=10,
                  handleInputParams=bamHandleInputParams, 
                  ext=c(-500, 200, -200, 500), 
                  hl=c(-100, 100, -100, 100), 
                  insert=100, 
                  stranded=TRUE, 
                  scale=FALSE, 
                  smooth=TRUE, 
                  outPrefix=op, 
                  nc=2)
   pdf_to_png(op)
})

test_that("testing plot_locus_with_random", {

   op <- "test_plot_locus_with_random"
   plot_locus_with_random(queryFiles=bamQueryFiles, 
                          centerFiles=bedQueryFiles[3], 
                          txdb, 
                          ext=c(-500,500), 
                          hl=c(-100,100), 
                          shade=TRUE, 
                          handleInputParams=bamHandleInputParams, 
                          binSize=10, 
                          refPoint="center", 
                          Xlab="Center", 
                          smooth=TRUE,
                          inputFiles=bamInputFiles, 
                          stranded=TRUE, 
                          scale=FALSE, 
                          outPrefix=op, 
                          verbose=FALSE, 
                          transform=NA, 
                          rmOutlier=0, 
                          n_random=1, 
                          statsMethod="wilcox.test", 
                          nc=2)
   pdf_to_png(op)
})
