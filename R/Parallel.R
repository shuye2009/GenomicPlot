
#' @title Prepare for parallel processing
#'
#' @description Method for starting a virtual cluster needed for parallel processing
#'
#' @param nc a positive integer greater than 1, denoting number of cores requested
#' @param verbose logical, whether to output additional information
#' 
#' @return an object of class c("SOCKcluster", "cluster"), depending on platform
#' @author Shuye Pu
#'
#' @examples
#' cl <- start_parallel(2L)
#'
#' @export start_parallel
#'
start_parallel <- function(nc=2, 
                           verbose=FALSE){

   n.cores <- detectCores()
   fnc <- min(nc, as.integer(n.cores-1))

   switch(Sys.info()[['sysname']],
          Windows = {my.cluster <- makeCluster(fnc, type = "PSOCK")},
          Linux = {my.cluster <- makeCluster(fnc, type = "FORK")}
   )

   if(verbose) print(paste("Using", fnc, "out of", n.cores, "cores!"))

   invisible(my.cluster)
}

#' @title Stop parallel processing
#
#' @description Method for stopping a virtual cluster needed for parallel processing
#'
#' @param cl a cluster or SOCKcluster object depending on platform
#' @return NULL
#' @author Shuye Pu
#'
#' @examples
#' cl <- start_parallel(2L)
#' stop_parallel(cl)
#'
#' @export stop_parallel
#'

stop_parallel <- function(cl){
   stopCluster(cl = cl)
}

#' @title Parallel execution of get_genomic_feature_coordinates
#' @description Get genomic coordinates for multiple features at once
#'
#' @param txdb a TxDb object
#' @param featureNames  a vector of gene features in c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene")
#' @param longest logical, indicating whether the output should be limited to the longest transcript of each gene
#' @param protein_coding logical, indicating whether to limit to protein_coding genes
#' @param nc number of cores to use
#'
#' @return a list of lists(Granges, GRangesList)
#' @author Shuye Pu
#' @keywords internal
#'
#' @note not working because txdb cannot be exported to nodes, need to find a solution
#'

parallel_feature_coordinates <- function(txdb, 
                                         featureNames, 
                                         longest=TRUE, 
                                         protein_coding=TRUE, 
                                         nc=2){
   tic()
   force(txdb)

   cl <- start_parallel(min(length(featureNames),nc))

   clusterExport(cl, varlist=c("get_genomic_feature_coordinates"), envir=environment())
   clusterExport(cl, varlist=c("txdb", "longest", "protein_coding"), envir=environment())
   alist <- clusterApply(cl, featureNames, get_genomic_feature_coordinates, txdb=txdb, longest=longest, protein_coding=protein_coding)
   stop_parallel(cl)
   toc()
   invisible(alist)
}

#' @title Parallel execution of scoreMatrixBin on a huge target windows object split into chunks
#'
#' @description Function for parallel computation of scoreMatrixBin. The 'windows' parameter of the scoreMatrixBin method is split into 5 chunks, and scoreMatrixBin is called on each chunk simultaneously to speed up the computation.
#'
#' @param windowRs, a single GRangesList object.
#' @param queryRegions, a RleList object or Granges object providing input for the 'target' parameter of the scoreMatrixBin method
#' @param bin_num, number of bins the windows should be divided into
#' @param bin_op, operation on the signals in a bin, a string in c("mean", "max", "min", "median", "sum") is accepted.
#' @param weight_col, if the queryRegions is a GRanges object, a numeric column in meta data part can be used as weights.
#' @param stranded, logical, indicating if the strand of the windows should be considered to determine upstream and downstream
#' @param nc, an integer denoting the number of cores requested, 2 is the default number that is allowed by CRAN but 5 gives best trade-off between speed and space
#'
#' @return a numeric matrix
#' @author Shuye Pu
#'
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
                                    nc=2){

   call_scoreMatrixBin <- function(windowR){
      ScoreMatrixBin(target = queryRegions, windows = windowR, bin.num=bin_num, bin.op=bin_op, weight.col=weight_col, strand.aware = stranded)
   }

   #print(system.time({
      cl <- start_parallel(nc)
      wRs <- split(windowRs, factor(cut(seq_along(windowRs), breaks=nc)))

      clusterExport(cl, varlist=c("ScoreMatrixBin"), envir=environment())
      clusterExport(cl, varlist=c("queryRegions", "bin_num", "bin_op", "weight_col", "stranded"), envir=environment())
      smc <- parLapply(cl, wRs, call_scoreMatrixBin)
      stop_parallel(cl)
      sm <- lapply(smc, function(x){
         y <- as(x, "matrix")
         colnames(y) <- seq(ncol(y))
         y
         })
      smdt <- Reduce(rbind, sm)
   #}))

   invisible(smdt)
}

#' @title Parallel execution of binnedAverage
#'
#' @description Function for parallel computation of binnedAverage function in the GenomicRanges package
#'
#' @param Rle_list a list of RleList objects.
#' @param tileBins, a GRanges object of tiled genome
#' @param nc integer, number of cores for parallel processing
#'
#' @return a list of numeric vectors
#' @author Shuye Pu
#'
#'
#' @export parallel_binnedAverage
#'
#'
parallel_binnedAverage <- function(Rle_list, 
                                   tileBins, 
                                   nc=2){

   #print(system.time({
      cl <- start_parallel(min(length(Rle_list), nc))

      clusterExport(cl, varlist=c("binnedAverage", "tileBins"), envir=environment())
      score_list <- parLapply(cl, Rle_list, function(x){
         binAverage <- binnedAverage(tileBins, x, varname="binned_score", na.rm=FALSE)
         binAverage$binned_score
      })

      stop_parallel(cl)
   #}))

   invisible(score_list)
}

#' @title Parallel execution of countOverlaps
#'
#' @description Function for parallel computation of countOverlaps function in the GenomicRanges package
#'
#' @param grange_list a list of GRanges objects.
#' @param tileBins, a GRanges object of tiled genome
#' @param switch, logical, switch the order of query and feature
#' @param nc integer, number of cores for parallel processing
#'
#' @return a list of numeric vectors
#' @author Shuye Pu
#'
#'
#' @export parallel_countOverlaps
#'
#'
parallel_countOverlaps <- function(grange_list, 
                                   tileBins, 
                                   nc=2, 
                                   switch=FALSE){

   #print(system.time({
      cl <- start_parallel(min(nc,length(grange_list)))

      clusterExport(cl, c("countOverlaps", "tileBins", "switch"), envir=environment())
      score_list <- parLapply(cl, grange_list, function(x){
         if(switch){
            binned_count <- GenomicRanges::countOverlaps(x, tileBins)
         }else{
            binned_count <- GenomicRanges::countOverlaps(tileBins, x)
         }
        
         binned_count
      })
      stop_parallel(cl)
   #}))

   invisible(score_list)
}
