
#' @title Rank rows of a matrix based on user input
#' @description The rows of a input numeric matrix is ordered based row sum, row maximum, or hierarchical clustering of the rows with euclidean
#' distannce and centroid linkage.
#'
#' @param fullmatrix a numeric matrix
#' @param ranking a string in c("Sum", "Max", "Hierarchical", "None")
#'
#' @return a numerical matrix
#' @author Shuye Pu
#'
#' @export rank_rows
#'
#'

rank_rows <- function(fullmatrix, ranking="Hierarchical"){
   fullmatrix <- data.matrix(fullmatrix)
   if(ranking == "None"){
      invisible(fullmatrix)
   }else if(ranking == "Sum"){
      fullmatrix <- arrange(as.data.frame(fullmatrix), desc(rowSums(fullmatrix)))
      invisible(data.matrix(fullmatrix))
   }else if(ranking == "Max"){
      fullmatrix <- arrange(as.data.frame(fullmatrix), desc(rowMax(fullmatrix)))
      invisible(data.matrix(fullmatrix))
   }else{
      clust <- hclust(dist(fullmatrix, method="euclidean"), method="centroid")
      invisible(data.matrix(fullmatrix[clust$order,]))
   }
}


#' @title Inspect a numeric matrix
#' @description Check the matrix for NA, NaN, INF, -INF and 0 values
#'
#' @param fullmatrix a numeric matrix
#' @param verbose logical, indicating whether to print out the stats in the console
#' @return NULL
#' @keywords internal

inspect_matrix <- function(fullmatrix, verbose=FALSE){

   if(verbose) print("Inspecting matrix")
   size <- nrow(fullmatrix)*ncol(fullmatrix)
   n_infinite <- sum(is.infinite(fullmatrix))
   n_NA <- sum(is.na(fullmatrix))
   n_NaN <- sum(is.nan(fullmatrix))
   n_zero <- sum(fullmatrix == 0)

   n_invalid <- c(n_infinite, n_NA, n_NaN, n_zero)
   fraction_invalid <- n_invalid/size

   stat_df <- data.frame(n_invalid, fraction_invalid)
   rownames(stat_df) <- c("infinite", "NA", "NaN", "zero")

   if(verbose) print(stat_df)

   return(NULL)

}

#' @title Impute missing values
#' @description Replace 0 and missing values in a matrix with half of minimum, to avoid use of arbitrary pseudo numbers
#' and to allow log transformation
#' @param fullmatrix a numeric matrix
#' @param verbose logical, whether to output additional information
#' 
#' @return a numeric matrix
#' @keywords internal
#'

impute_hm <- function(fullmatrix, verbose=FALSE){

   if(verbose){
      message("Imputing missing values. Matrix quartiles:")
      print(quantile(fullmatrix))
   }

   #minv <- median(fullmatrix[fullmatrix != 0])
   minv <- median(apply(fullmatrix, 2, mean))
   halfmin <- ifelse(minv>0, minv/2, minv*2)
   fullmatrix[fullmatrix < halfmin] <- halfmin ##  to avoid take log of zero and use of pseudo numbers

   if(verbose) print(paste("Imputed value", halfmin, "\n\n"))

   return(fullmatrix)
}

#' @title Preprocess scoreMatrix before plotting
#'
#' @description  This is a helper function for manipulate the score matrix produced by ScoreMatrix or ScoreMatrinBin functions defined in the 'genomation' package.
#'
#' @param fullmatrix a numeric matrix, with bins in columns and genomic windows in rows
#' @param scale logical, indicating whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param rmOutlier logical, indicating whether a row with abnormally high values in the score matrix should be removed
#' @param verbose logical, indicating whether to output additional information (data used for plotting or statistical test results)
#' @param transform a string in c("log", "log2", "log10"), default = NA indicating no transformation of data matrix
#' @param pc pseudo-count added to the data matrix before log transformation to avoid taking log of zero
#'
#' @return a numeric matrix with the same dimension as the fullmatrix
#' @author Shuye Pu
#'
#'
#' @export process_scoreMatrix
#'
#'
process_scoreMatrix <- function(fullmatrix, scale=FALSE, rmOutlier=FALSE, transform=NA, pc=0, verbose=FALSE){

   #rn <- rownames(fullmatrix)
   inspect_matrix(fullmatrix, verbose=verbose)

   fullmatrix[is.infinite(fullmatrix)] <- 0
   fullmatrix[is.na(fullmatrix)] <- 0

   fullmatrix <- impute_hm(fullmatrix, verbose)
   
   inspect_matrix(fullmatrix, verbose)

   if(!is.na(transform)) {
      if(min(fullmatrix) < 0){
         message("Negative values are found in the matrix, log transformation cannot be applied!")
      }else if(tranform == "log"){
         fullmatrix <- log(fullmatrix + pc)
      }else if(tranform == "log2"){
         fullmatrix <- log2(fullmatrix + pc)
      }else if(tranform == "log10"){
         fullmatrix <- log10(fullmatrix + pc)
      }
   }
   ## remove outliers from reference regions, using Hampel filter with 1000mad instead of 3mad.
   ## if outliers are detected, replace the outliers with up bound
   if(rmOutlier){
      fullmatrix <- rm_outlier(fullmatrix, verbose=verbose)
   }

   if(scale){
      #smc <- t(apply(fullmatrix, 1, scales::rescale)) ## rescale to 0:1 range
      smc <- t(base::scale(t(fullmatrix))) ## rescale to zscore by row
      fullmatrix <- as.matrix(smc)
   }
   fullmatrix[is.na(fullmatrix)] <- 0

   #rownames(fullmatrix) <- rn
   invisible(fullmatrix)
}

#' @title Remove outliers from scoreMatrix
#'
#' @description This is a helper function for dealing with excessively high values using Hampel filter. If outliers are detected, replace the outliers with the up bound = median(rowmax) + multiplier*mad(rowmax). This function is experimental. For data with normal distribution, the multiplier is usually set at 3. As the read counts data distribution is highly skewed, it is difficult to define a boundary for outliers, try the multiplier values between 10 to 1000.
#'
#' @param fullmatrix a numeric matrix, with bins in columns and genomic windows in rows
#' @param verbose logical, whether to output the outlier information to a log file
#' @param multiplier a numeric value to multiple the 'mad', default 1000, maybe adjusted based on data
#'
#' @return a numeric matrix
#' @author Shuye Pu
#'
#' @examples
#' fullmatrix <- matrix(rnorm(100), ncol=10)
#' fullmatrix[3,9] <- max(fullmatrix) + 1000000
#' rm_outlier(fullmatrix)
#'
#'
#' @export rm_outlier
#'

rm_outlier <- function(fullmatrix, verbose=FALSE, multiplier=1000){

   fullmatrix[is.na(fullmatrix)] <- 0
   rowmax <- apply(fullmatrix, 1, max)
   M <- mad(rowmax)
   if(M > 0){
      up_bound <- median(rowmax) + multiplier*M
   }else{
      up_bound <- rowmax
   }

   fullmatrix <- as.matrix(fullmatrix)
   fn <- ecdf(fullmatrix)
   percentile <- fn(up_bound)

   if(length(which(rowmax > up_bound)) > 0){
      if(verbose) print("Outlier detected:")
      outliers <- fullmatrix[fullmatrix>up_bound]

      if(verbose){
         print(paste("maximum of the matrix", max(rowmax)))
         print(paste("median of row max", median(rowmax)))
         print(paste("median absolute deviation (mad) of row max", M))
         print(paste("mulitplier of mad", multiplier))
         print(paste("up_bound and replace value", up_bound))
         print(paste("percentile of up_bound", percentile))
         print(paste("number of outlier rows", length(which(rowmax > up_bound))))
         print(paste("number of outliers", length(outliers)))
         print(paste("fraction of outliers", length(outliers)/(nrow(fullmatrix)*ncol(fullmatrix))))
         print("values of outliers")
         print(outliers)
      }

      fullmatrix[fullmatrix>up_bound] <- up_bound
   }

   invisible(fullmatrix)
}

#' @title Perform one-way ANOVA and post hoc TukeyHSD tests
#'
#' @description This is a helper function for performing one-way ANOVA analysis and post hoc Tukey's Honest Significant Differences tests
#'
#' @param df a dataframe with c("Intensity", "Group") in column names
#' @param xc a string denoting column name for grouping
#' @param yc a string denoting column name for numeric data to be plotted
#' @param op output prefix for statistical analysis results
#' @param verbose logical, to indicate whether a file should be produced to save the test results
#'
#' @return a list of two elements, the first is the p-value of ANOVA test and the second is a matrix of the output of TukeyHSD tests
#' @author Shuye Pu
#'
#' @note used in \code{plot_reference_locus}
#' @export aov_TukeyHSD
#'

aov_TukeyHSD <- function(df, xc="Group", yc="Intensity", op=NULL, verbose=FALSE){
   if(verbose){
      sink(file=paste0(op,"_TukeyHSD.txt"), append=TRUE, split=TRUE)
      cat(paste("Performing one-way ANOVA analysis for", op, "\n"))
   }
   formu <- as.formula(paste(yc, "~", xc))
   res.aov <- aov(formu, data = df)
   s <- summary(res.aov)
   p <-  s[[1]][1, "Pr(>F)"]
   if(verbose){
      print(s) ## anova summary
      cat("\nPost hoc Tukey Honest Significant Differences test\n")
   } 
   v <- TukeyHSD(res.aov, which = xc)
   if(verbose){
      print(v)
      sink()
   }
   invisible(list("ANOVA"=p, "HSD"=v[[xc]]))
}

#' @title Convert GRanges to dataframe
#' @description Convert GRanges object with metacolumns to dataframe
#'
#' @param gr a GRanges object
#'
#' @return a dataframe
#'
#' @author Shuye Pu
#'
#' @examples
#' gr2 <- GenomicRanges::GRanges(c("chr1", "chr1"), IRanges::IRanges(c(7,13), width=3),
#'  strand=c("+", "-"))
#' GenomicRanges::mcols(gr2) <- data.frame(score=c(0.3, 0.9), cat=c(TRUE, FALSE))
#' df2 <- gr2df(gr2)
#'
#' @export gr2df

gr2df <- function(gr){

   chr <- as.vector(seqnames(gr))
   start <- start(gr)
   end <- end(gr)
   width <- width(gr)
   strand <- as.vector(strand(gr))
   meta <- mcols(gr, us.names=TRUE)
   
   if("name" %in% colnames(meta)){
      name <- gr$name
      meta <- meta[, !colnames(meta) %in% c("name"), drop=FALSE]
   }else if("id" %in% colnames(meta)){
         name <- gr$id
         meta <- meta[, !colnames(meta) %in% c("id"), drop=FALSE]
   }else if(!is.null(names(gr))){
      name <- names(gr)
   }else{
      name <- seq_along(gr)
   }
   
   if("score" %in% colnames(meta)){
      score <- gr$score
      meta <- meta[, !colnames(meta) %in% c("score"), drop=FALSE]
   }else{
      score <- width
   }

   df <- as.data.frame(cbind(chr, start, end, name, score, strand, meta))
   
   return(df)
}


