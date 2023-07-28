library(GenomicPlot)
library(testthat)

Sys.setenv("R_TESTS" = "")

data(gf5_meta)
data(gf5_genomic)

gtffile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf", 
                       package = "GenomicPlot")
txdb <- AnnotationDbi::loadDb(system.file("extdata", "txdb.sql", 
                                          package = "GenomicPlot"))

bedQueryFiles <- c(
   system.file("extdata", "test_chip_peak_chr19.narrowPeak", 
               package = "GenomicPlot"),
   system.file("extdata", "test_chip_peak_chr19.bed", package = "GenomicPlot"),
   system.file("extdata", "test_clip_peak_chr19.bed", package = "GenomicPlot")
)
names(bedQueryFiles) <- c("NarrowPeak", "SummitPeak", "iCLIPPeak")

bedImportParams <- setImportParams(
   offset = 0, fix_width = 100, fix_point = "center", norm = FALSE,
   useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
)

bamQueryFiles <- system.file("extdata", "treat_chr19.bam", 
                             package = "GenomicPlot")
names(bamQueryFiles) <- "clip_bam"
bamInputFiles <- system.file("extdata", "input_chr19.bam", 
                             package = "GenomicPlot")
names(bamInputFiles) <- "clip_input"

bamImportParams <- setImportParams(
   offset = -1, fix_width = 0, fix_point = "start", norm = TRUE,
   useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
)

chipQueryFiles <- system.file("extdata", "chip_treat_chr19.bam",
                              package = "GenomicPlot")
names(chipQueryFiles) <- "chip_bam"
chipInputFiles <- system.file("extdata", "chip_input_chr19.bam",
                              package = "GenomicPlot")
names(chipInputFiles) <- "chip_input"

chipImportParams <- setImportParams(
   offset = 0, fix_width = 150, fix_point = "start", norm = TRUE,
   useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
)


test_that("testing parallel_countOverlaps", {
   importParams <- setImportParams(fix_width = 100, outRle = FALSE)
   out_list <- handle_input(
      inputFiles = bedQueryFiles,
      importParams = importParams, verbose = FALSE, nc = 2
   )
   
   seqi <- GenomeInfoDb::Seqinfo(genome = "hg19") 
   grange_list <- lapply(out_list, function(x) x$query) 
   tilewidth <- 100000 
   tileBins <- tileGenome(seqi, tilewidth = tilewidth, 
                          cut.last.tile.in.chrom = TRUE)
   
   score_list1 <- parallel_countOverlaps(grange_list, tileBins, nc = 2)
})

test_that("testing parallel_scoreMatrixBin", {
   
   queryRegion <- handle_input(chipQueryFiles, chipImportParams, 
                               verbose = TRUE)[[1]]$query
   
   importParams <- setImportParams(outRle = FALSE)
   
   windowRegion <- handle_bed(bedQueryFiles[1], importParams, verbose = TRUE)$query
   
   out <- parallel_scoreMatrixBin(
      queryRegions = queryRegion,
      windowRs = windowRegion, 
      bin_num = 50, 
      bin_op = "mean",
      weight_col = "score",
      stranded = TRUE,
      nc = 2
   )
})

test_that("testing handle_bed", {
   out <- handle_bed(bedQueryFiles[1], bedImportParams, verbose = TRUE)
})

test_that("testing effective_size", {
   importParams <- setImportParams(outRle = FALSE)
   out_list <- handle_input(
      inputFiles = c(chipQueryFiles, chipInputFiles),
      importParams = importParams, verbose = TRUE, nc = 2
   )
   
   out <- effective_size(out_list, outRle = TRUE)
})

test_that("testing handle_input", {
   
   queryFiles2 <- system.file("extdata", "test_wig_chr19_+.wig", 
                              package = "GenomicPlot")
   names(queryFiles2) <- "test_wig"
   
   queryFiles3 <- system.file("extdata", "test_wig_chr19_+.bw", 
                              package = "GenomicPlot")
   names(queryFiles3) <- "test_bw" 
   
   importParams <- setImportParams()
   
   out <- handle_input(c(bamQueryFiles, queryFiles2, queryFiles3), 
                       importParams, verbose = TRUE)
})

test_that("testing plot_bam_correlation", {
   
   importParams <- setImportParams(fix_width = 150, outRle = FALSE)
   
   plot_bam_correlation(
      bamFiles = c(chipQueryFiles, chipInputFiles), binSize = 100000, 
      outPrefix = NULL, importParams = importParams, nc = 2, verbose = FALSE
   )
})

test_that("testing plot_overlap_bed", {
   importParams <- setImportParams(fix_width = 100, outRle = FALSE)
   plot_overlap_bed(
      bedList = bedQueryFiles, importParams = importParams, pairOnly = FALSE,
      stranded = FALSE, outPrefix = NULL
   )
})

test_that("testing plot_argument_list", {
   
   alist <- list(
      "txdb" = txdb, "treat" = bamQueryFiles, "control" = bamInputFiles, 
      "feature" = gf5_meta, "param" = bamImportParams
   )
   
   p <- GenomicPlot:::plot_named_list(alist)
})


test_that("testing plot_peak_annotation", {
   plot_peak_annotation(
      peakFile = bedQueryFiles[2], gtfFile = gtffile, importParams = bedImportParams,
      fiveP = -2000, dsTSS = 200, threeP = 2000, simple = FALSE
   )  
})


test_that("testing plot_overlap_genes", {
   testfile1 <- system.file("extdata", "test_file1.txt",  
                            package = "GenomicPlot")
   testfile2 <- system.file("extdata", "test_file2.txt",  
                            package = "GenomicPlot")
   testfile3 <- system.file("extdata", "test_file3.txt",  
                            package = "GenomicPlot")
   testfile4 <- system.file("extdata", "test_file4.txt",  
                            package = "GenomicPlot")
   testfiles <- c(testfile1, testfile2, testfile3, testfile4)
   names(testfiles) <- c("test1", "test2", "test3", "test4") 
   
   plot_overlap_genes(testfiles, c(3,2,1,1), pairOnly = FALSE)
})

test_that("testing plot_5parts_metagene", {
   plot_5parts_metagene(
      queryFiles = bedQueryFiles,
      gFeatures_list = list("metagene" = gf5_meta),
      inputFiles = NULL,
      importParams = bedImportParams,
      verbose = FALSE,
      smooth = TRUE,
      scale = FALSE,
      stranded = TRUE,
      outPrefix = NULL,
      transform = NA,
      heatmap = TRUE,
      rmOutlier = 0,
      heatRange = NULL,
      nc = 2
   )
})

test_that("testing plot_locus", {
   plot_locus(
      queryFiles = bedQueryFiles[c(1, 3)],
      centerFiles = bedQueryFiles[2],
      ext = c(-500, 500),
      hl = c(-100, 100),
      inputFiles = NULL,
      importParams = bedImportParams,
      shade = TRUE,
      binSize = 10,
      refPoint = "center",
      Xlab = "Summit",
      verbose = FALSE,
      smooth = TRUE,
      scale = FALSE,
      stranded = TRUE,
      outPrefix = NULL,
      transform = NA,
      heatmap = TRUE,
      heatRange = NULL,
      rmOutlier = 0,
      Ylab = "Coverage/base/peak",
      nc = 2
   )
})

test_that("testing plot_region", {
   plot_region(
      queryFiles = chipQueryFiles,
      centerFiles = bedQueryFiles[1],
      inputFiles = chipInputFiles,
      nbins = 100,
      heatmap = TRUE,
      scale = FALSE,
      regionName = "narrowPeak",
      importParams = chipImportParams,
      verbose = FALSE,
      fiveP = -200,
      threeP = 200,
      smooth = TRUE,
      transform = "log2",
      stranded = TRUE,
      Ylab = "Coverage/base/peak",
      outPrefix = NULL,
      rmOutlier = 0,
      nc = 2
   )
})

test_that("testing plot_start_end", {
   plot_start_end(
      queryFiles = bamQueryFiles,
      inputFiles = bamInputFiles,
      txdb = txdb,
      centerFiles = "intron",
      binSize = 10,
      importParams = bamImportParams,
      ext = c(-100, 100, -100, 100),
      hl = c(-50, 50, -50, 50),
      insert = 100,
      stranded = TRUE,
      scale = FALSE,
      smooth = TRUE,
      transform = "log2",
      outPrefix = NULL,
      nc = 2
   )
})

test_that("testing plot_start_end_with_random", {
   
   plot_start_end_with_random(
      queryFiles = bamQueryFiles,
      inputFiles = bamInputFiles,
      txdb = txdb,
      centerFile = "intron",
      binSize = 10,
      importParams = bamImportParams,
      ext = c(-100, 100, -100, 100),
      hl = c(-20, 20, -20, 20),
      insert = 100,
      stranded = TRUE,
      scale = FALSE,
      smooth = TRUE,
      transform = "log2",
      outPrefix = NULL,
      randomize = TRUE,
      nc = 2
   )
})

test_that("testing plot_locus_with_random", {
   plot_locus_with_random(
      queryFiles = bamQueryFiles,
      centerFiles = bedQueryFiles[3],
      txdb = txdb,
      ext = c(-200, 200),
      hl = c(-100, 100),
      shade = FALSE,
      importParams = bamImportParams,
      verbose = FALSE,
      smooth = FALSE,
      transform = NA,
      binSize = 10,
      refPoint = "center",
      Xlab = "Center",
      Ylab = "Coverage/base/peak",
      inputFiles = bamInputFiles,
      stranded = TRUE,
      scale = FALSE,
      outPrefix = NULL,
      rmOutlier = 0,
      n_random = 1,
      hw = c(8, 8),
      detailed = FALSE,
      statsMethod = "wilcox.test",
      nc = 2)
})
