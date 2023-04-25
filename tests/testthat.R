library(GenomicPlot)
library(magick)

Sys.setenv("R_TESTS" = "")

pdf_to_png <- function(op){
   img <- image_read_pdf(paste0(op,".pdf"))
   for(i in seq_along(img)){image_write(img[i], path=paste0(op,"_",i, ".png"), format="png")}
}

outdir <- "./test_output"
setwd(outdir)

gtffile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf", package="GenomicPlot")
txdb <- AnnotationDbi::loadDb(system.file("extdata", "txdb_chr19.sql", package="GenomicPlot"))

## plot_5parts_metagene ####
gf <- prepare_5parts_genomic_features(txdb, meta=TRUE, nbins=100, fiveP=-2000, threeP=1000, longest=TRUE)

bedQueryFiles <- c(system.file("extdata", "test_chip_peak_chr19.narrowPeak", package="GenomicPlot"),
                system.file("extdata", "test_chip_peak_chr19.bed", package="GenomicPlot"),
                system.file("extdata", "test_clip_peak_chr19.bed", package="GenomicPlot"))
names(bedQueryFiles) <- c("Narrow", "Summit", "iCLIP")

bedHandleInputParams <- list(CLIP_reads=FALSE, fix_width=100, fix_point="center", norm=FALSE, 
                             useScore=FALSE, outRle=TRUE, useSizeFactor=FALSE, genome="hg19")

op <- "test_plot_5parts_metagene1"
plot_5parts_metagene(queryFiles=bedQueryFiles,
                     gFeatures_list=list(metaF=gf),
                     inputFiles=NULL,
                     handleInputParams=bedHandleInputParams,
                     verbose=FALSE,
                     smooth=TRUE,
                     scale=FALSE,
                     stranded=TRUE,
                     outPrefix=op,
                     transform=NA,
                     heatmap=TRUE,
                     rmOutlier=TRUE,
                     heatRange=c(0,1),
                     nc=2)
pdf_to_png(op)


bamQueryFiles <- system.file("extdata", "treat_chr19.bam", package="GenomicPlot")
names(bamQueryFiles) <- "query"
bamInputFiles <- system.file("extdata", "input_chr19.bam", package="GenomicPlot")
names(bamInputFiles) <- "input"

bamHandleInputParams <- list(CLIP_reads=TRUE, fix_width=0, fix_point="start", norm=TRUE, 
                          useScore=FALSE, outRle=TRUE, useSizeFactor=FALSE, genome="hg19")

op <- "test_plot_5parts_metagene2"
plot_5parts_metagene(queryFiles=bamQueryFiles, 
                     gFeatures=list("metagene"=gf), 
                     inputFiles=bamInputFiles, 
                     scale=FALSE, 
                     verbose=FALSE, 
                     transform=NA, 
                     smooth=TRUE, 
                     stranded=TRUE, 
                     outPrefix=op, 
                     handleInputParams=bamHandleInputParams, 
                     heatmap=TRUE, 
                     rmOutlier=TRUE, 
                     nc=2)
pdf_to_png(op)

## plot_3parts_metagene ####
gf_gene <- prepare_3parts_genomic_features(txdb, meta=FALSE, nbins=100, fiveP=-3000, threeP=2000, 
                                      longest=TRUE)

chipQueryFiles <- system.file("extdata", "chip_treat_chr19.bam", package="GenomicPlot")
names(chipQueryFiles) <- "query"
chipInputFiles <- system.file("extdata", "chip_input_chr19.bam", package="GenomicPlot")
names(chipInputFiles) <- "input"

chipHandleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, 
                          useScore=FALSE, outRle=TRUE, useSizeFactor=TRUE, genome="hg19")

op <- "test_plot_3parts_metagene"
plot_3parts_metagene(queryFiles=chipQueryFiles, 
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
                     rmOutlier=TRUE, 
                     nc=2)
pdf_to_png(op)

## plot_locus ####
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
           rmOutlier=TRUE,
           nc=2)
pdf_to_png(op)

centerfiles <- c(system.file("extdata", "test_clip_peak_chr19.bed", package="GenomicPlot"),
                 system.file("extdata", "test_chip_peak_chr19.bed", package="GenomicPlot"))
names(centerfiles) <- c("clip_peak", "chip_peak")
queryfiles <- c(system.file("extdata", "treat_chr19.bam", package="GenomicPlot"),
                system.file("extdata", "chip_treat_chr19.bam", package="GenomicPlot"))
names(queryfiles) <- c("clip_bam", "chip_bam")
inputfiles <- c(system.file("extdata", "input_chr19.bam", package="GenomicPlot"),
                system.file("extdata", "chip_input_chr19.bam", package="GenomicPlot"))
names(inputfiles) <- c("clip_input", "chip_input")

handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, 
                          useScore=FALSE, outRle=TRUE, useSizeFactor=TRUE, genome="hg19")

op <- "test_plot_locus2"
plot_locus(queryFiles=queryfiles, 
           centerFiles=centerfiles, 
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
           rmOutlier=FALSE, 
           statsMethod="wilcox.test", 
           heatmap=TRUE, 
           nc=2)
pdf_to_png(op)

## plot_peak_annotation ####
op <- "test_plot_peak_annotation1"
peakHandleInputParams <- list(CLIP_reads=FALSE, fix_width=100, fix_point="center", norm=FALSE, 
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
peakHandleInputParams <- list(CLIP_reads=FALSE, fix_width=21, fix_point="center", norm=FALSE, 
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

## plot_region ####
centerfiles <- system.file("extdata", "test_chip_peak_chr19.narrowPeak", package="GenomicPlot")
names(centerfiles) <- c("narrowpeak")
queryfiles <- c(system.file("extdata", "treat_chr19.bam", package="GenomicPlot"),
                system.file("extdata", "chip_treat_chr19.bam", package="GenomicPlot"))
names(queryfiles) <- c("clip_bam", "chip_bam")
inputfiles <- c(system.file("extdata", "input_chr19.bam", package="GenomicPlot"),
                system.file("extdata", "chip_input_chr19.bam", package="GenomicPlot"))
names(inputfiles) <- c("clip_input", "chip_input")

handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, 
                          useScore=FALSE, outRle=TRUE, useSizeFactor=FALSE, genome="hg19")


op <- "test_plot_region"
plot_region(queryFiles=queryfiles, 
            centerFiles=centerfiles, 
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
            rmOutlier=FALSE, 
            nc=2)
pdf_to_png(op)

## plot_start_end ####
op <- "test_plot_stat_end"
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

## plot_locus_with_random ####
centerfiles <- c(system.file("extdata", "test_clip_peak_chr19.bed", package="GenomicPlot"))
names(centerfiles) <- c("clip_peak")
queryfiles <- c(system.file("extdata", "treat_chr19.bam", package="GenomicPlot"))
names(queryfiles) <- c("clip_bam")
inputfiles <- c(system.file("extdata", "input_chr19.bam", package="GenomicPlot"))
names(inputfiles) <- c("clip_input")

handleInputParams <- list(CLIP_reads=TRUE, fix_width=150, fix_point="start", norm=TRUE, 
                          useScore=FALSE, outRle=TRUE, useSizeFactor=FALSE, genome="hg19")

op <- "test_plot_locus_with_random"
plot_locus_with_random(queryFiles=queryfiles, 
                       centerFiles=centerfiles, 
                       txdb, 
                       ext=c(-500,500), 
                       hl=c(-100,100), 
                       shade=TRUE, 
                       handleInputParams=handleInputParams, 
                       binSize=10, 
                       refPoint="center", 
                       Xlab="Center", 
                       smooth=TRUE,
                       inputFiles=inputfiles, 
                       stranded=TRUE, 
                       scale=FALSE, 
                       outPrefix=op, 
                       verbose=FALSE, 
                       transform=NA, 
                       rmOutlier=FALSE, 
                       n_random=1, 
                       statsMethod="wilcox.test", 
                       nc=2)
pdf_to_png(op)
