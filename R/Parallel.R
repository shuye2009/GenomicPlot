#' @title Prepare for parallel processing
#'
#' @description Creating a virtual cluster for parallel processing
#'
#' @param nc a positive integer greater than 1, denoting number of cores requested
#' @param verbose logical, whether to output additional information
#'
#' @return an object of class c("SOCKcluster", "cluster"), depending on platform
#' @author Shuye Pu
#'
#' @examples
#' 
#' cl <- start_parallel(2L)
#' stop_parallel(cl)
#'
#' @export start_parallel
#'
start_parallel <- function(nc = 2,
                           verbose = FALSE) {
  n.cores <- detectCores()
  fnc <- min(nc, as.integer(n.cores - 1))

  switch(Sys.info()[["sysname"]],
    Windows = {
      my.cluster <- makeCluster(fnc, type = "PSOCK")
    },
    Linux = {
      my.cluster <- makeCluster(fnc, type = "FORK")
    }
  )

  if (verbose) message("Using ", fnc, " out of ", n.cores, " cores!\n")

  invisible(my.cluster)
}

#' @title Stop parallel processing
#
#' @description Stopping a virtual cluster after parallel processing is finished
#'
#' @param cl a cluster or SOCKcluster object depending on platform
#' @return 0 if the cluster is stopped successfully, 1 otherwise.
#' @author Shuye Pu
#'
#' @examples
#' cl <- start_parallel(2L)
#' stop_parallel(cl)
#'
#' @export stop_parallel
#'

stop_parallel <- function(cl) {
  r <- stopCluster(cl = cl)
  if (is.null(r)) {
     return(0)
  } else {
     return(1)
  }
}

#' @title Parallel execution of scoreMatrixBin on a huge target windows object split into chunks
#'
#' @description Function for parallel computation of scoreMatrixBin. The 'windows' parameter of the scoreMatrixBin method is split into nc chunks, and scoreMatrixBin is called on each chunk simultaneously to speed up the computation.
#'
#' @param windowRs a single GRangesList object.
#' @param queryRegions a RleList object or Granges object providing input for the 'target' parameter of the scoreMatrixBin method
#' @param bin_num number of bins the windows should be divided into
#' @param bin_op operation on the signals in a bin, a string in c("mean", "max", "min", "median", "sum") is accepted.
#' @param weight_col if the queryRegions is a GRanges object, a numeric column in meta data part can be used as weights.
#' @param stranded logical, indicating if the strand of the windows should be considered to determine upstream and downstream
#' @param nc an integer denoting the number of cores requested, 2 is the default number that is allowed by CRAN but 5 gives best trade-off between speed and space
#'
#' @return a numeric matrix
#' @author Shuye Pu
#' 
#' @examples 
#' queryFiles <- system.file("extdata", "chip_treat_chr19.bam", package = "GenomicPlot")
#' names(queryFiles) <- "query"
#'
#' importParams <- list(
#'   offset = 0, fix_width = 140, fix_point = "start", norm = FALSE,
#'   useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' queryRegion <- handle_bam(queryFiles, importParams, verbose = TRUE)$query
#' 
#' windowFiles <- system.file("extdata", "test_chip_peak_chr19.narrowPeak", 
#'   package = "GenomicPlot"
#' )
#' names(windowFiles) <- "narrowPeak"
#'
#' importParams <- list(
#'   offset = 0, fix_width = 0, fix_point = "start", norm = FALSE,
#'   useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' windowRegion <- handle_bed(windowFiles, importParams, verbose = TRUE)$query
#'
#' out <- parallel_scoreMatrixBin(
#'   queryRegions = queryRegion,
#'   windowRs = windowRegion, 
#'   bin_num = 50, 
#'   bin_op = "mean",
#'   weight_col = "score",
#'   stranded = TRUE,
#'   nc = 2
#' )
#' 
#' \donttest{
#' draw_matrix_heatmap(out, dataName = "chipData", labels_col = NULL, 
#'   levels_col = NULL, ranking = "Sum", ranges = NULL, verbose = FALSE)
#' } 
#'  
#' @export parallel_scoreMatrixBin
#'
#'
parallel_scoreMatrixBin <- function(queryRegions,
                                    windowRs,
                                    bin_num,
                                    bin_op,
                                    weight_col,
                                    stranded,
                                    nc = 2) {
  call_scoreMatrixBin <- function(windowR) {
    ScoreMatrixBin(target = queryRegions, windows = windowR, bin.num = bin_num, bin.op = bin_op, weight.col = weight_col, strand.aware = stranded)
  }

  cl <- start_parallel(nc)
  wRs <- split(windowRs, factor(cut(seq_along(windowRs), breaks = nc)))

  clusterExport(cl, varlist = c("ScoreMatrixBin"), envir = environment())
  clusterExport(cl, varlist = c("queryRegions", "bin_num", "bin_op", "weight_col", "stranded"), envir = environment())
  smc <- parLapply(cl, wRs, call_scoreMatrixBin)
  stop_parallel(cl)
  sm <- lapply(smc, function(x) {
    y <- as(x, "matrix")
    colnames(y) <- seq(ncol(y))
    y
  })
  smdt <- Reduce(rbind, sm)

  invisible(smdt)
}

#' @title Parallel execution of countOverlaps
#'
#' @description Function for parallel computation of countOverlaps function in the GenomicRanges package
#'
#' @param grange_list a list of GRanges objects.
#' @param tileBins a GRanges object of tiled genome
#' @param switch logical, switch the order of query and feature
#' @param nc integer, number of cores for parallel processing
#'
#' @return a list of numeric vectors
#' @author Shuye Pu
#' 
#' @examples 
#' queryFiles <- system.file("extdata", "treat_chr19.bam", package = "GenomicPlot")
#' names(queryFiles) <- "query"
#'
#' inputFiles <- system.file("extdata", "input_chr19.bam", package = "GenomicPlot")
#' names(inputFiles) <- "input"
#'
#' importParams <- list(
#'   offset = -1, fix_width = 0, fix_point = "start", norm = FALSE,
#'   useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' out_list <- handle_input(
#'   inputFiles = c(queryFiles, inputFiles),
#'   importParams = importParams, verbose = TRUE, nc = 2
#' )
#' 
#' seqi <- GenomeInfoDb::Seqinfo(genome = "hg19") 
#' grange_list <- lapply(out_list, function(x) x$query) 
#' tilewidth <- 100000 
#' tileBins <- tileGenome(seqi, tilewidth = tilewidth, cut.last.tile.in.chrom = TRUE)
#'  
#' score_list1 <- parallel_countOverlaps(grange_list, tileBins, nc = 2)
#' dplyr::glimpse(score_list1)
#' score_list2 <- parallel_countOverlaps(grange_list, tileBins, nc = 2, switch = TRUE)
#' dplyr::glimpse(score_list2)
#'
#' @export parallel_countOverlaps
#'
#'
parallel_countOverlaps <- function(grange_list,
                                   tileBins,
                                   nc = 2,
                                   switch = FALSE) {
  
  cl <- start_parallel(min(nc, length(grange_list)))

  clusterExport(cl, c("countOverlaps", "tileBins", "switch"), envir = environment())
  score_list <- parLapply(cl, grange_list, function(x) {
    if (switch) {
      binned_count <- GenomicRanges::countOverlaps(x, tileBins)
    } else {
      binned_count <- GenomicRanges::countOverlaps(tileBins, x)
    }

    binned_count
  })
  stop_parallel(cl)

  invisible(score_list)
}
