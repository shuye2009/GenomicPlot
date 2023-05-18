## Generate RDS files with consistent import parameters

library(GenomicPlot)

setwd("C:/GREENBLATT/Rscripts/GenomicPlotData")

bedQueryFiles <- c(
   "test_chip_peak_chr19.narrowPeak",
   "test_chip_peak_chr19.bed",
   "test_clip_peak_chr19.bed"
)
names(bedQueryFiles) <- c("NarrowPeak", "SummitPeak", "iCLIPPeak")

bedimportParams <- list(
   offset = 0, fix_width = 100, fix_point = "center", norm = FALSE,
   useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
)

out_list <- handle_input(
   inputFiles = bedQueryFiles,
   importParams = bedimportParams, verbose = TRUE, nc = 2
)

queryFiles <- "treat_chr19.bam"
names(queryFiles) <- "query"

inputFiles <- "input_chr19.bam"
names(inputFiles) <- "input"

bamimportParams <- list(
  offset = -1, fix_width = 0, fix_point = "start", norm = TRUE,
  useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
)

out_list <- handle_input(
  inputFiles = c(queryFiles, inputFiles),
  importParams = bamimportParams, verbose = TRUE, nc = 2
)

queryFiles <- "chip_treat_chr19.bam"
names(queryFiles) <- "query"

inputFiles <- "chip_input_chr19.bam"
names(inputFiles) <- "input"

chipimportParams <- list(
  offset = 0, fix_width = 150, fix_point = "start", norm = TRUE,
  useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
)
 
out_list <- handle_input(
   inputFiles = c(queryFiles, inputFiles),
   importParams = chipimportParams, verbose = TRUE, nc = 2
)

queryFiles <- "test_wig_chr19_+.wig"
names(queryFiles) <- "test_wig"

wigimportParams <- list(
  offset = 0, fix_width = 0, fix_point = "start", norm = FALSE,
  useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
)

out <- handle_input(queryFiles, wigimportParams, verbose = TRUE) 

queryFiles <- "test_wig_chr19_+.bw"
names(queryFiles) <- "test_bw" 
 
out <- handle_input(queryFiles, wigimportParams, verbose = TRUE) 
