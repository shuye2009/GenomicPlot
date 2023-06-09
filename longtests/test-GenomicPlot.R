library(GenomicPlot)
library(magick)
library(testthat)

Sys.setenv("R_TESTS" = "")

pdf_to_png <- function(op) {
  img <- image_read_pdf(paste0(op, ".pdf"))
  for (i in seq_along(img)) {
    image_write(img[i], path = paste0(op, "_", i, ".pdf"), format = "pdf", quality = 50)
  }
}

outdir <- "./test_output"
setwd(outdir)

gtffile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf", package = "GenomicPlot")
gff <- RCAS::importGtf(readFromRds = TRUE, filePath = gtffile)
txdb <- makeTxDbFromGRanges(gff)

bedQueryFiles <- c(
  system.file("extdata", "test_chip_peak_chr19.narrowPeak", package = "GenomicPlot"),
  system.file("extdata", "test_chip_peak_chr19.bed", package = "GenomicPlot"),
  system.file("extdata", "test_clip_peak_chr19.bed", package = "GenomicPlot")
)
names(bedQueryFiles) <- c("NarrowPeak", "SummitPeak", "iCLIPPeak")

bedimportParams <- list(
  offset = 0, fix_width = 100, fix_point = "center", norm = FALSE,
  useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
)

bamQueryFiles <- system.file("extdata", "treat_chr19.bam", package = "GenomicPlot")
names(bamQueryFiles) <- "clip_bam"
bamInputFiles <- system.file("extdata", "input_chr19.bam", package = "GenomicPlot")
names(bamInputFiles) <- "clip_input"

bamimportParams <- list(
  offset = -1, fix_width = 0, fix_point = "start", norm = TRUE,
  useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
)

chipQureyFiles <- system.file("extdata", "chip_treat_chr19.bam", package = "GenomicPlot")
names(chipQureyFiles) <- "chip_bam"
chipInputFiles <- system.file("extdata", "chip_input_chr19.bam", package = "GenomicPlot")
names(chipInputFiles) <- "chip_input"

chipimportParams <- list(
  offset = 0, fix_width = 150, fix_point = "start", norm = TRUE,
  useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
)


test_that("testing plot_5parts_metagene", {
  gf <- prepare_5parts_genomic_features(txdb, meta = TRUE, nbins = 100, fiveP = -2000, threeP = 1000, longest = TRUE)
  op <- "test_plot_5parts_metagene1"
  plot_5parts_metagene(
    queryFiles = bedQueryFiles,
    gFeatures_list = list("metagene" = gf),
    inputFiles = NULL,
    importParams = bedimportParams,
    verbose = FALSE,
    smooth = TRUE,
    scale = FALSE,
    stranded = TRUE,
    outPrefix = op,
    transform = NA,
    heatmap = TRUE,
    rmOutlier = 0,
    heatRange = NULL,
    nc = 2
  )
 # pdf_to_png(op)

  op <- "test_plot_5parts_metagene2"
  plot_5parts_metagene(
    queryFiles = bamQueryFiles,
    gFeatures_list = list("metagene" = gf),
    inputFiles = bamInputFiles,
    scale = FALSE,
    verbose = FALSE,
    transform = "log2",
    smooth = TRUE,
    stranded = TRUE,
    outPrefix = op,
    importParams = bamimportParams,
    heatmap = TRUE,
    rmOutlier = 0,
    nc = 2
  )
  #pdf_to_png(op)
})

test_that("testing plot_locus", {
  op <- "test_plot_locus1"
  plot_locus(
    queryFiles = bedQueryFiles[c(1, 3)],
    centerFiles = bedQueryFiles[2],
    ext = c(-1000, 1000),
    hl = c(-100, 100),
    inputFiles = NULL,
    importParams = bedimportParams,
    shade = TRUE,
    binSize = 10,
    refPoint = "center",
    Xlab = "Summit",
    verbose = FALSE,
    smooth = TRUE,
    scale = FALSE,
    stranded = TRUE,
    outPrefix = op,
    transform = NA,
    heatmap = TRUE,
    heatRange = NULL,
    rmOutlier = 0,
    nc = 2
  )
  #pdf_to_png(op)
  
  op <- "test_plot_locus2"
  plot_locus(
    queryFiles = chipQureyFiles,
    centerFiles = bedQueryFiles[2:3],
    ext = c(-500, 500),
    hl = c(-100, 100),
    shade = TRUE,
    smooth = TRUE,
    importParams = chipimportParams,
    binSize = 10,
    refPoint = "center",
    Xlab = "Center",
    inputFiles = chipInputFiles,
    stranded = TRUE,
    scale = FALSE,
    outPrefix = op,
    verbose = FALSE,
    transform = "log2",
    rmOutlier = 0,
    statsMethod = "wilcox.test",
    heatmap = TRUE,
    nc = 2
  )
 # pdf_to_png(op)
})
test_that("testing plot_peak_annotation", {
  op <- "test_plot_peak_annotation1"
  plot_peak_annotation(
    peakFile = bedQueryFiles[2],
    gtfFile = gtffile,
    importParams = bedimportParams,
    fiveP = -2000,
    dsTSS = 200,
    threeP = 1000,
    outPrefix = op,
    verbose = FALSE
  )
  #pdf_to_png(op)

  op <- "test_plot_peak_annotation2"
  plot_peak_annotation(
    peakFile = bedQueryFiles[3],
    gtfFile = gtffile,
    importParams = bedimportParams,
    fiveP = -1000,
    dsTSS = 0,
    threeP = 2000,
    outPrefix = op,
    verbose = FALSE
  )
  #pdf_to_png(op)
})

test_that("testing plot_region", {
  op <- "test_plot_region"
  plot_region(
    queryFiles = chipQureyFiles,
    centerFiles = bedQueryFiles[1],
    inputFiles = chipInputFiles,
    nbins = 100,
    heatmap = TRUE,
    scale = FALSE,
    regionName = "narrowPeak",
    importParams = chipimportParams,
    verbose = FALSE,
    fiveP = -500,
    threeP = 500,
    smooth = TRUE,
    transform = "log2",
    stranded = TRUE,
    outPrefix = op,
    rmOutlier = 0,
    nc = 2
  )
  #pdf_to_png(op)
})

test_that("testing plot_start_end", {
  op <- "test_plot_start_end"
  plot_start_end(
    queryFiles = bamQueryFiles,
    inputFiles = bamInputFiles,
    txdb = txdb,
    centerFiles = "intron",
    binSize = 10,
    importParams = bamimportParams,
    ext = c(-500, 200, -200, 500),
    hl = c(-100, 100, -100, 100),
    insert = 100,
    stranded = TRUE,
    scale = FALSE,
    smooth = TRUE,
    transform = "log2",
    outPrefix = op,
    nc = 2
  )
  #pdf_to_png(op)
})

test_that("testing plot_start_end_with_random", {
  op <- "test_plot_start_end_with_random"
  plot_start_end_with_random(
    queryFiles = bamQueryFiles,
    inputFiles = bamInputFiles,
    txdb = txdb,
    centerFile = "intron",
    binSize = 10,
    importParams = bamimportParams,
    ext = c(-500, 200, -200, 500),
    hl = c(-100, 100, -100, 100),
    insert = 100,
    stranded = TRUE,
    scale = FALSE,
    smooth = TRUE,
    transform = "log2",
    outPrefix = op,
    nc = 2
  )
  #pdf_to_png(op)
})

if(0){
test_that("testing plot_locus_with_random", {
  op <- "test_plot_locus_with_random"
  plot_locus_with_random(
    queryFiles = bamQueryFiles,
    centerFiles = bedQueryFiles[3],
    txdb,
    ext = c(-500, 500),
    hl = c(-100, 100),
    shade = TRUE,
    importParams = bamimportParams,
    binSize = 10,
    refPoint = "center",
    Xlab = "Center",
    smooth = TRUE,
    inputFiles = bamInputFiles,
    stranded = TRUE,
    scale = FALSE,
    outPrefix = op,
    verbose = FALSE,
    transform = "log2",
    rmOutlier = 0,
    n_random = 1,
    statsMethod = "wilcox.test",
    nc = 2
  )
  #pdf_to_png(op)
})
}
# remove files not used by README.md and GenomicPlot_vignettes.rmd

pdfs <- list.files(pattern = "*.pdf")
used_pdfs <- c("test_plot_5parts_metagene1.pdf",
               "test_plot_5parts_metagene2.pdf",
               "test_plot_region.pdf",
               "test_plot_locus2.pdf",
               "test_plot_peak_annotation1.pdf"
              )
unused_pdfs <- pdfs[!pdfs %in% used_pdfs]
unlink(unused_pdfs)

