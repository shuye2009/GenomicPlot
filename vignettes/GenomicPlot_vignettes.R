## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE, fig.show = 'hold', fig.keep = 'all', fig.align = 'center', fig.dim = c(7,7), fig.ncol=1, fig.sep="\n\n"----
#  library(GenomicPlot, quietly = TRUE)
#  txdb <- AnnotationDbi::loadDb(system.file("extdata", "txdb_chr19.sql", package = "GenomicPlot"))
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
#  importParams <- list(
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
#    outPrefix = NULL,
#    importParams = importParams,
#    heatmap = TRUE,
#    rmOutlier = 0,
#    nc = 5
#  )

## ----metagene1, echo = FALSE, fig.cap = "Individual sample profile and heatmap", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_5parts_metagene2_1.png")

## ----metagene2, echo = FALSE, fig.cap = "Profile overlay", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_5parts_metagene2_2.png")

## ----metagene3, echo = FALSE, fig.cap = "Ratio-over-input profile and heatmap", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_5parts_metagene2_3.png")

## ---- eval = FALSE, fig.show = 'hold', fig.keep = 'all', fig.align = 'center', fig.dim = c(7,7)----
#  txdb <- AnnotationDbi::loadDb(system.file("extdata", "txdb_chr19.sql", package = "GenomicPlot"))
#  gf <- prepare_3parts_genomic_features(txdb,
#    meta = FALSE, nbins = 100, fiveP = -3000, threeP = 2000,
#    longest = TRUE
#  )
#  
#  queryfiles <- system.file("extdata", "chip_treat_chr19.bam", package = "GenomicPlot")
#  names(queryfiles) <- "chip_bam"
#  inputfiles <- system.file("extdata", "chip_input_chr19.bam", package = "GenomicPlot")
#  names(inputfiles) <- "chip_input"
#  
#  importParams <- list(
#    offset = 0, fix_width = 150, fix_point = "start", norm = TRUE,
#    useScore = FALSE, outRle = TRUE, useSizeFactor = TRUE, genome = "hg19"
#  )
#  
#  plot_3parts_metagene(
#    queryFiles = queryfiles,
#    gFeatures = gf,
#    inputFiles = inputfiles,
#    scale = FALSE,
#    verbose = FALSE,
#    transform = NA,
#    smooth = TRUE,
#    stranded = TRUE,
#    outPrefix = NULL,
#    importParams = importParams,
#    heatmap = TRUE,
#    rmOutlier = 0,
#    nc = 5
#  )

## ----gene1, echo = FALSE, fig.cap = "Individual sample profile and heatmap", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_3parts_metagene_1.png")

## ----gene2, echo = FALSE, fig.cap = "Profile overlay", out.width = '75%'------
knitr::include_graphics("../inst/tests/test_output/test_plot_3parts_metagene_2.png")

## ----gene3, echo = FALSE, fig.cap = "Ratio-over-input profile and heatmap", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_3parts_metagene_3.png")

## ---- eval = FALSE, fig.show = 'hold', fig.keep = 'all', fig.align = 'center', fig.dim = c(7,7)----
#  centerfiles <- system.file("extdata", "test_chip_peak_chr19.narrowPeak", package = "GenomicPlot")
#  names(centerfiles) <- c("NarrowPeak")
#  queryfiles <- c(
#    system.file("extdata", "treat_chr19.bam", package = "GenomicPlot"),
#    system.file("extdata", "chip_treat_chr19.bam", package = "GenomicPlot")
#  )
#  names(queryfiles) <- c("clip_bam", "chip_bam")
#  inputfiles <- c(
#    system.file("extdata", "input_chr19.bam", package = "GenomicPlot"),
#    system.file("extdata", "chip_input_chr19.bam", package = "GenomicPlot")
#  )
#  names(inputfiles) <- c("clip_input", "chip_input")
#  
#  importParams <- list(
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
#    importParams = importParams,
#    verbose = FALSE,
#    fiveP = -500,
#    threeP = 500,
#    smooth = TRUE,
#    transform = NA,
#    stranded = TRUE,
#    outPrefix = NULL,
#    rmOutlier = 0,
#    nc = 5
#  )

## ----region1, echo = FALSE, fig.cap = "Individual sample profile", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_region_21.png")

## ----region2, echo = FALSE, fig.cap = "Individual sample coverage stats in narrowPeak", out.width = '100%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_region_22.png")

## ----region3, echo = FALSE, fig.cap = "Individual sample coverage profiles and heatmaps in narrowPeak", out.width = '100%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_region_23.png")

## ----region4, echo = FALSE, fig.cap = "Ratio-over-input profile", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_region_24.png")

## ----region5, echo = FALSE, fig.cap = "Ratio-over-input stats", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_region_25.png")

## ----region6, echo = FALSE, fig.cap = "Ratio-over-input profiles and heatmaps", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_region_26.png")

## ----region7, echo = FALSE, fig.cap = "Program names and arguments", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_region_27.png")

## ---- eval = FALSE, fig.show = 'hold', fig.keep = 'all', fig.align = 'center', fig.dim = c(7,7)----
#  txdb <- AnnotationDbi::loadDb(system.file("extdata", "txdb_chr19.sql", package = "GenomicPlot"))
#  queryfiles <- system.file("extdata", "treat_chr19.bam", package = "GenomicPlot")
#  names(queryfiles) <- "clip_bam"
#  inputfiles <- system.file("extdata", "input_chr19.bam", package = "GenomicPlot")
#  names(inputfiles) <- "clip_input"
#  ext <- c(-500, 200, -200, 500)
#  hl <- c(-50, 50, -50, 50)
#  
#  importParams <- list(
#    offset = -1, fix_width = 0, fix_point = "start", norm = TRUE,
#    useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#  )
#  plot_start_end(
#    queryFiles = queryfiles,
#    inputFiles = inputfiles,
#    txdb = txdb,
#    centerFiles = "intron",
#    binSize = 10,
#    importParams = importParams,
#    ext = ext,
#    hl = hl,
#    insert = 100,
#    stranded = TRUE,
#    scale = FALSE,
#    smooth = TRUE,
#    outPrefix = NULL,
#    nc = 5
#  )

## ----intron1, echo = FALSE, fig.cap = "Query and input sample profiles in start, center and end of intron", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_stat_end_3.png")

## ----intron2, echo = FALSE, fig.cap = "Ratio-over profile in start, center and end of intron", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_stat_end_5.png")

## ---- eval = FALSE, fig.show = 'hold', fig.keep = 'all', fig.align = 'center', fig.dim = c(7,7)----
#  centerfiles <- c(
#    system.file("extdata", "test_clip_peak_chr19.bed", package = "GenomicPlot"),
#    system.file("extdata", "test_chip_peak_chr19.bed", package = "GenomicPlot")
#  )
#  names(centerfiles) <- c("iCLIPPeak", "SummitPeak")
#  queryfiles <- c(
#    system.file("extdata", "treat_chr19.bam", package = "GenomicPlot"),
#    system.file("extdata", "chip_treat_chr19.bam", package = "GenomicPlot")
#  )
#  names(queryfiles) <- c("clip_bam", "chip_bam")
#  inputfiles <- c(
#    system.file("extdata", "input_chr19.bam", package = "GenomicPlot"),
#    system.file("extdata", "chip_input_chr19.bam", package = "GenomicPlot")
#  )
#  names(inputfiles) <- c("clip_input", "chip_input")
#  
#  importParams <- list(
#    offset = 0, fix_width = 150, fix_point = "start", norm = TRUE,
#    useScore = FALSE, outRle = TRUE, useSizeFactor = TRUE, genome = "hg19"
#  )
#  
#  plot_locus(
#    queryFiles = queryfiles,
#    centerFiles = centerfiles,
#    ext = c(-500, 500),
#    hl = c(-100, 100),
#    shade = TRUE,
#    smooth = TRUE,
#    importParams = importParams,
#    binSize = 10,
#    refPoint = "center",
#    Xlab = "Center",
#    inputFiles = inputfiles,
#    stranded = TRUE,
#    scale = FALSE,
#    outPrefix = NULL,
#    verbose = FALSE,
#    transform = NA,
#    rmOutlier = 0,
#    statsMethod = "wilcox.test",
#    heatmap = TRUE,
#    nc = 5
#  )

## ----locus11, echo = FALSE, fig.cap = "Ratio-over-input for clip signal around clip and chip peaks", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_locus2_21.png")

## ----locus12, echo = FALSE, fig.cap = "Ratio-over-input for clip signal around clip and chip peaks", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_locus2_22.png")

## ----locus21, echo = FALSE, fig.cap = "Ratio-over-input for chip signal around clip and chip peaks", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_locus2_23.png")

## ----locus22, echo = FALSE, fig.cap = "Ratio-over-input for chip signal around clip and chip peaks", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_locus2_24.png")

## ----locus31, echo = FALSE, fig.cap = "Ratio-over-input for chip and clip signal around clip peaks", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_locus2_25.png")

## ----locus32, echo = FALSE, fig.cap = "Ratio-over-input for chip and clip signal around clip peaks", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_locus2_26.png")

## ----locus41, echo = FALSE, fig.cap = "Ratio-over-input for chip and clip signal around chip peaks", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_locus2_27.png")

## ----locus42, echo = FALSE, fig.cap = "Ratio-over-input for chip and clip signal around chip peaks", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_locus2_28.png")

## ----locus5, echo = FALSE, fig.cap = "Ratio-over-input profiles for chip and clip signal around clip and chip peaks", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_locus2_29.png")

## ----locus6, echo = FALSE, fig.cap = "Ratio-over-input profiles and heatmap for chip and clip signal", out.width = '100%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_locus2_30.png")

## ---- eval = FALSE, fig.show = 'hold', fig.keep = 'all', fig.align = 'center', fig.dim = c(7,7)----
#  txdb <- AnnotationDbi::loadDb(system.file("extdata", "txdb_chr19.sql", package = "GenomicPlot"))
#  centerfiles <- c(system.file("extdata", "test_clip_peak_chr19.bed", package = "GenomicPlot"))
#  names(centerfiles) <- c("iCLIPPeak")
#  queryfiles <- c(system.file("extdata", "treat_chr19.bam", package = "GenomicPlot"))
#  names(queryfiles) <- c("clip_bam")
#  inputfiles <- c(system.file("extdata", "input_chr19.bam", package = "GenomicPlot"))
#  names(inputfiles) <- c("clip_input")
#  
#  importParams <- list(
#    offset = -1, fix_width = 150, fix_point = "start", norm = TRUE,
#    useScore = FALSE, outRle = TRUE, useSizeFactor = TRUE, genome = "hg19"
#  )
#  
#  plot_locus_with_random(
#    queryFiles = queryfiles,
#    centerFiles = centerfiles,
#    txdb,
#    ext = c(-500, 500),
#    hl = c(-100, 100),
#    shade = TRUE,
#    importParams = importParams,
#    binSize = 10,
#    refPoint = "center",
#    Xlab = "Center",
#    smooth = TRUE,
#    inputFiles = inputfiles,
#    stranded = TRUE,
#    scale = FALSE,
#    outPrefix = NULL,
#    verbose = FALSE,
#    transform = NA,
#    rmOutlier = 0,
#    n_random = 1,
#    statsMethod = "wilcox.test",
#    nc = 5
#  )

## ----locusr11, echo = FALSE, fig.cap = "Ratio-over-input for clip signal around clip peaks in 5'UTR", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_locus_with_random_27.png")

## ----locusr12, echo = FALSE, fig.cap = "Ratio-over-input for clip signal around clip peaks in 5'UTR", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_locus_with_random_28.png")

## ----locusr21, echo = FALSE, fig.cap = "Ratio-over-input for clip signal around clip peaks", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_locus_with_random_35.png")

## ----locusr22, echo = FALSE, fig.cap = "Ratio-over-input for clip signal around clip peaks", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_locus_with_random_36.png")

## ---- eval = FALSE, fig.show = 'hold', fig.keep = 'all', fig.align = 'center', fig.dim = c(7,7)----
#  gtffile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf", package = "GenomicPlot")
#  
#  centerfile <- system.file("extdata", "test_chip_peak_chr19.bed", package = "GenomicPlot")
#  names(centerfile) <- c("SummitPeak")
#  
#  handleBedparams <- list(
#    fix_width = 0, fix_point = "center", useScore = FALSE, outRle = FALSE,
#    offset = 0, norm = FALSE, useSizeFactor = FALSE, genome = "hg19"
#  )
#  
#  pa <- plot_peak_annotation(
#    peakFile = centerfile,
#    gtfFile = gtffile,
#    importParams = handleBedparams,
#    fiveP = -2000,
#    dsTSS = 300,
#    threeP = 1000,
#    simple = FALSE,
#    verbose = TRUE,
#    outPrefix = NULL
#  )

## ----annotation1, echo = FALSE, fig.cap = "Distribution of chip peak in various types of gene", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_peak_annotation1_1.png")

## ----annotation2, echo = FALSE, fig.cap = "Distribution of chip peak in various parts of gene", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_peak_annotation1_2.png")

## ----annotation3, echo = FALSE, fig.cap = "Density of chip peak in various parts of gene", out.width = '75%'----
knitr::include_graphics("../inst/tests/test_output/test_plot_peak_annotation1_3.png")

## ----last---------------------------------------------------------------------
sessionInfo()

