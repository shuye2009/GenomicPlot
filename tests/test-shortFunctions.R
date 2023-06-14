library(GenomicPlot)
library(testthat)

Sys.setenv("R_TESTS" = "")

test_that("testing plot_bam_correlation", {
   bamQueryFiles <- c(
      system.file("extdata", "chip_input_chr19.bam", package = "GenomicPlot"),
      system.file("extdata", "chip_treat_chr19.bam", package = "GenomicPlot")
    )
    names(bamQueryFiles) <- c("chip_input", "chip_treat")
    
    bamImportParams <- list(
       offset = 0, fix_width = 150, fix_point = "start", norm = FALSE,
       useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19"
    )
   
    op <- "test_plot_bam_correlation"
    plot_bam_correlation(
      bamFiles = bamQueryFiles, binSize = 100000, outPrefix = op,
      importParams = bamImportParams, nc = 2, verbose = FALSE
    )
})
 
test_that("testing plot_overlap_bed", {
  op <- "test_plot_overlap_bed"
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
    stranded = FALSE, outPrefix = op
  )
})

test_that("testing plot_argument_list", {
   
    gtffile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf", package = "GenomicPlot")
    gff <- RCAS::importGtf(saveObjectAsRds = TRUE, filePath = gtffile)
    txdb <- makeTxDbFromGRanges(gff)
   
    queryfiles <- system.file("extdata", "treat_chr19.bam", package = "GenomicPlot")
    names(queryfiles) <- "query"
   
    inputfiles <- system.file("extdata", "input_chr19.bam", package = "GenomicPlot")
    names(inputfiles) <- "input"
   
    gfeatures <- prepare_5parts_genomic_features(txdb,
      meta = TRUE, nbins = 100, fiveP = -1000,
      threeP = 1000, longest = TRUE, verbose = FALSE
    )
   
    bamimportParams <- list(
      offset = -1, fix_width = 0, fix_point = "start", norm = TRUE,
      useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
    )
   
    alist <- list(
      "txdb" = txdb, "treat" = queryfiles, "control" = inputfiles, 
      "feature" = gfeatures, "param" = bamimportParams
    )
    
    p <- GenomicPlot:::plot_named_list(alist)
    pdf("test_plot_argument_list.pdf")
    print(p)
    dev.off()
})

unlink("*.pdf")
 
