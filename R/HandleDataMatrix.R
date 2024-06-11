#' @title Rank rows of a matrix based on user input
#' @description The rows of a input numeric matrix is ordered based row sum,
#' row maximum, or hierarchical clustering of the rows with euclidean distance
#' and centroid linkage. This a helper function for drawing matrix heatmaps.
#'
#' @param fullmatrix a numeric matrix
#' @param ranking a string in c("Sum", "Max", "Hierarchical", "None")
#'
#' @return a numeric matrix
#' @author Shuye Pu
#'
#' @examples
#'
#' fullMatrix <- matrix(rnorm(100), ncol = 10)
#' for (i in 5:8) {
#'     fullMatrix[i, 4:7] <- runif(4) + i
#' }
#' apply(fullMatrix, 1, sum)
#' ranked <- rank_rows(fullMatrix, ranking = "Sum")
#' apply(ranked, 1, sum)
#'
#' @export rank_rows

rank_rows <- function(fullmatrix,
                      ranking = "Hierarchical") {
    stopifnot(ranking %in% c("Sum", "Max", "Hierarchical", "None"))
    fullmatrix <- data.matrix(fullmatrix)
    if (ranking == "None") {
        invisible(fullmatrix)
    } else if (ranking == "Sum") {
        fullmatrix <- arrange(as.data.frame(fullmatrix),
                              desc(rowSums(fullmatrix)))
        invisible(data.matrix(fullmatrix))
    } else if (ranking == "Max") {
        fullmatrix <- arrange(as.data.frame(fullmatrix),
                              desc(Biobase::rowMax(fullmatrix)))
        invisible(data.matrix(fullmatrix))
    } else {
        clust <- hclust(dist(fullmatrix, method = "euclidean"),
                        method = "centroid")
        invisible(data.matrix(fullmatrix[clust$order, ]))
    }
}


#' @title Inspect a numeric matrix
#' @description Check the matrix for NA, NaN, INF, -INF and 0 values
#'
#' @param fullmatrix a numeric matrix
#' @param verbose logical, indicating whether to print out the stats in the
#'  console
#' @return a numerical matrix summarizing the unusual values
#'
#' @author Shuye Pu
#'
#' @examples
#' fullMatrix <- matrix(rnorm(100), ncol = 10)
#' for (i in 5:6) {
#'     fullMatrix[i, 4:7] <- NaN
#'     fullMatrix[i + 1, 4:7] <- NA
#'     fullMatrix[i + 2, 4:7] <- -Inf
#'     fullMatrix[i - 1, 4:7] <- 0
#'     fullMatrix[i - 2, 1:3] <- Inf
#' }
#'
#' GenomicPlot:::inspect_matrix(fullMatrix, verbose = TRUE)
#'
#' @keywords internal

inspect_matrix <- function(fullmatrix,
                           verbose = FALSE) {
    if (verbose) message("[inspect_matrix] started\n")
    size <- nrow(fullmatrix) * ncol(fullmatrix)
    n_infinite <- sum(is.infinite(fullmatrix))
    n_NA <- sum(is.na(fullmatrix))
    n_NaN <- sum(is.nan(fullmatrix))
    n_zero <- sum(fullmatrix == 0.0, na.rm = TRUE)
    n_negative <- sum(fullmatrix < 0, na.rm = TRUE)

    n_invalid <- c(n_infinite, n_NA, n_NaN, n_zero, n_negative)
    fraction_invalid <- n_invalid / size

    stat_df <- data.frame(n_invalid, fraction_invalid)
    rownames(stat_df) <- c("infinite", "NA", "NaN", "zero", "negative")

    if (verbose){
        print(stat_df)
        message("[inspect_matrix] finished\n")
    }

    invisible(stat_df)
}

#' @title Impute missing values
#' @description Replace 0 and missing values in a sparse non-negative matrix
#' with half of minimum of non-zero values, to avoid use of arbitrary pseudo
#' numbers, and to allow computing ratios and log transformation of matrices.
#' When a matrix is sparse (assuming it has many all-zero rows and few all-zero
#' columns), the half of minimum of non-zero values is a number that is small
#' enough so that is will not distort the data too much (comparing to a pseudo
#' count = 1), but large enough to avoid huge ratios when used as a denominator.
#'
#' @param fullmatrix a numeric matrix
#' @param verbose logical, whether to output additional information
#'
#' @return a numeric matrix
#'
#' @author Shuye Pu
#'
#' @examples
#' fullMatrix <- matrix(rlnorm(100), ncol = 10)
#' for (i in 5:6) {
#'     fullMatrix[i - 1, 4:7] <- 0
#' }
#'
#' imp <- GenomicPlot:::impute_hm(fullMatrix, verbose = TRUE)
#'
#' @keywords internal
#'

impute_hm <- function(fullmatrix,
                      verbose = FALSE) {
    if (min(fullmatrix) < 0) {
        message("Cannot impute because the matrix has negative values!\n")
        return(fullmatrix)
    }
    if (verbose) {
        message("[impute_hm] Imputing missing values...\nMatrix quartiles:\n")
        message(paste(quantile(fullmatrix), collapse = " "), "\n")
    }

    minv <- min(fullmatrix[fullmatrix != 0])

    halfmin <- minv / 2
    fullmatrix[fullmatrix < halfmin] <- halfmin
    ##  to avoid taking log of zero and to avoid using of pseudo numbers

    if (verbose) {
        message("\nMatrix quantiles after imputing:\n")
        message(paste(quantile(fullmatrix), collapse = " "), "\n")
        message("The imputed value is: ", halfmin, "\n")
        message("[impute_hm] finished!\n")
    }

    return(fullmatrix)
}

#' @title Preprocess scoreMatrix before plotting
#'
#' @description  This is a helper function for manipulate the score matrix
#' produced by ScoreMatrix or ScoreMatrinBin functions defined in the
#' `genomation` package. To facilitate downstream analysis, imputation of
#' missing values is performed implicitly when log transformation is required,
#' otherwise missing values are replaced with 0.
#'
#' @param fullmatrix a numeric matrix, with bins in columns and genomic windows
#'  in rows
#' @param scale logical, indicating whether the score matrix should be scaled
#'  to the range 0:1, so that samples with different baseline can be compared
#' @param rmOutlier a numeric value to multiple the 'mad' when detecting
#'  outliers, can be adjusted based on data.  Default 0, indicating not to
#'  remove outliers.
#' @param verbose logical, indicating whether to output additional information
#'  (data used for plotting or statistical test results)
#' @param transform a string in c("log", "log2", "log10"), default = NA
#'  indicating no transformation of data matrix
#'
#' @details If inputFiles for the plotting function is null, all operations
#'  (scale, rmOutlier and transform) can be applied to the score matrix, in
#'  the order of rmOutlier -> transform -> scale. When inputFiles are provided,
#'  only rmOutlier can be applied to the score matrix, as transform and scale
#'  will affect ratio calculation, especially when log2 transformation of the
#'  ratio is intended. However, all these operations can be applied to the
#'  resulting ratio matrix. In order to avoid introducing distortion into the
#'  processed data, use caution when applying these operations.
#'
#' @return a numeric matrix with the same dimension as the fullmatrix
#' @author Shuye Pu
#'
#' @examples
#' fullMatrix <- matrix(rlnorm(100), ncol = 10)
#' for (i in 5:6) {
#'     fullMatrix[i, 4:7] <- NaN
#'     fullMatrix[i + 1, 4:7] <- NA
#'     fullMatrix[i + 2, 4:7] <- -Inf
#'     fullMatrix[i - 1, 4:7] <- 0
#'     fullMatrix[i - 2, 1:3] <- Inf
#' }
#' fullMatrix[9, 4:7] <- runif(4) + 90
#'
#' wo <- process_scoreMatrix(fullMatrix, rmOutlier = 3, verbose = TRUE)
#' tf <- process_scoreMatrix(fullMatrix,
#'     rmOutlier = 0, transform = "log2", verbose = TRUE
#' )
#' scaled <- process_scoreMatrix(fullMatrix, scale = TRUE, verbose = TRUE)
#'
#' @export process_scoreMatrix
#'
#'
process_scoreMatrix <- function(fullmatrix,
                                scale = FALSE,
                                rmOutlier = 0,
                                transform = NA,
                                verbose = FALSE) {
    stopifnot(is.numeric(rmOutlier))
    stopifnot(is.logical(scale))
    stopifnot(transform %in% c("log", "log2", "log10", NA))

    if(verbose) message("[process_scoreMatrix] started ...\n")
    if(is.null(nrow(fullmatrix))) fullmatrix <- matrix(fullmatrix, nrow = 1)

    inspect_matrix(fullmatrix, verbose = verbose)

    fullmatrix[is.infinite(fullmatrix)] <- 0
    fullmatrix[is.na(fullmatrix)] <- 0

    ## remove outliers from reference regions, using Hampel filter with
    ## rmOutlier * mad instead of 3 * mad, which is generally used for normal
    ## distribution.
    ## if outliers are detected, replace the outliers with up bound
    if (rmOutlier > 0) {
        fullmatrix <- rm_outlier(fullmatrix, verbose = verbose,
                                 multiplier = rmOutlier)
        if(is.null(nrow(fullmatrix))) fullmatrix <- matrix(fullmatrix, nrow = 1)
    }

    if (!is.na(transform)) {

        fullmatrix <- impute_hm(fullmatrix, verbose)
        if(is.null(nrow(fullmatrix))) fullmatrix <- matrix(fullmatrix, nrow = 1)

        ## impute to avoid taking log of zero, also to avoid distortion by
        ## adding pseudocount 1 when the vast majority of values of the matrix
        ## are less than 1, like in the case of a ratio matrix
        if (min(fullmatrix) < 0) {
            message("Negative values are found in the matrix, log transformation
                    cannot be applied!\n")
        } else if (transform == "log") {
            fullmatrix <- log(fullmatrix)
        } else if (transform == "log2") {
            fullmatrix <- log2(fullmatrix)
        } else if (transform == "log10") {
            fullmatrix <- log10(fullmatrix)
        }
    }

    if (scale) {
        smc <- t(apply(fullmatrix, 1, scales::rescale)) # rescale to 0:1 range
        # smc <- t(base::scale(t(fullmatrix))) ## rescale to zscore by row
        fullmatrix <- as.matrix(smc)

        fullmatrix[is.na(fullmatrix)] <- 0
        allSame <- apply(fullmatrix, 1, function(x) all(x == mean(x)))
        fullmatrix[allSame, ] <- 0
        ## rescale will set the entire row to 0.5 if all values are 0,
        ## which will distort the downstream analysis
        count_allSame_rows <- sum(allSame)
        if (count_allSame_rows > 0 && verbose) {
            message(count_allSame_rows, " rows have only one distinct value in
                    the entire row after rescale!\n")
        }
    }
    fullmatrix[is.na(fullmatrix)] <- 0
    if(is.null(nrow(fullmatrix))) fullmatrix <- matrix(fullmatrix, nrow = 1)

    if(verbose) message("[process_scoreMatrix] finished!\n")
    invisible(fullmatrix)
}

#' @title Remove outliers from scoreMatrix
#'
#' @description This is a helper function for dealing with excessively high
#' values using Hampel filter. If outliers are detected, replace the outliers
#' with the up bound = median(rowmax) + multiplier*mad(rowmax). This function
#' is experimental. For data with normal distribution, the multiplier is
#' usually set at 3. As the read counts data distribution is highly skewed, it
#' is difficult to define a boundary for outliers, try the multiplier values
#' between 10 to 1000.
#'
#' @param fullmatrix a numeric matrix, with bins in columns and genomic windows
#'  in rows
#' @param verbose logical, whether to output the outlier information to the
#'  console
#' @param multiplier a numeric value to multiple the 'mad', default 1000, maybe
#'  adjusted based on data
#'
#' @return a numeric matrix
#' @author Shuye Pu
#'
#' @examples
#' fullmatrix <- matrix(rnorm(100), ncol = 10)
#' maxm <- max(fullmatrix)
#' fullmatrix[3, 9] <- maxm + 1000
#' fullmatrix[8, 1] <- maxm + 500
#' rm_outlier(fullmatrix, verbose = TRUE, multiplier = 100)
#' rm_outlier(fullmatrix, verbose = TRUE, multiplier = 1000)
#'
#' @export rm_outlier
#'

rm_outlier <- function(fullmatrix, verbose = FALSE, multiplier = 1000) {
    fullmatrix[is.na(fullmatrix)] <- 0
    rowmax <- apply(fullmatrix, 1, max)
    k <- 1.4826
    ## k is the scaling constant for estimating rolling standard deviation
    ## using median absolute deviation, its value is 1.4826 most of the time
    M <- mad(rowmax) * k

    if (M > 0) {
        up_bound <- median(rowmax) + multiplier * M
    } else {
        up_bound <- max(rowmax) * 0.99
        ## for extremely skewed data, use 99th percentile of the maximum
    }

    fullmatrix <- as.matrix(fullmatrix)
    fn <- ecdf(fullmatrix)
    percentile <- fn(up_bound)

    if (length(which(rowmax > up_bound)) > 0) {
        outliers <- fullmatrix[fullmatrix > up_bound]

        if (verbose) {
            message("Outlier detected!!!\n")
            message("Maximum of the matrix: ", max(fullmatrix))
            message("\nMedian of row max: ", median(rowmax))
            message("\nMedian absolute deviation (mad) of row max: ", M)
            message("\nMulitplier of mad: ", multiplier)
            message("\nUp_bound and replace value: ", up_bound)
            message("\nPercentile of up_bound: ", percentile)
            message("\nNumber of outlier rows: ",
                    length(which(rowmax > up_bound)))
            message("\nNumber of outliers: ", length(outliers))
            message("\nFraction of outliers: ",
                    length(outliers) / (nrow(fullmatrix) * ncol(fullmatrix)))
            message("\nValues of outliers:\n")
            message(paste(outliers, collapse = " "), "\n")
        }

        fullmatrix[fullmatrix > up_bound] <- up_bound
    } else {
        if (verbose) {
            message("Outlier not detected:\n")
            message("Maximum of the matrix: ", max(rowmax))
            message("\nMedian of row max: ", median(rowmax))
            message("\nMedian absolute deviation (mad) of row max: ", M)
            message("\nMulitplier of mad: ", multiplier)
            message("\nUp_bound and replace value: ", up_bound)
            message("\nPercentile of up_bound: ", percentile, "\n")
        }
    }

    invisible(fullmatrix)
}

#' @title Perform one-way ANOVA and post hoc TukeyHSD tests
#'
#' @description This is a helper function for performing one-way ANOVA analysis
#' and post hoc Tukey's Honest Significant Differences tests
#'
#' @param df a dataframe
#' @param xc a string denoting column name for grouping
#' @param yc a string denoting column name for numeric data to be plotted
#' @param op output prefix for statistical analysis results
#' @param verbose logical, to indicate whether a file should be produced to
#'  save the test results
#'
#' @return a list of two elements, the first is the p-value of ANOVA test and
#' the second is a matrix of the output of TukeyHSD tests
#' @author Shuye Pu
#'
#' @note used in \code{plot_locus}
#'
#' @examples
#' stat_df <- data.frame(
#'     Feature = rep(c("A", "B"), c(20, 30)),
#'     Intensity = c(rnorm(20, mean = 2, sd = 1), rnorm(30, mean = 3, sd = 1))
#' )
#'
#' out <- aov_TukeyHSD(stat_df, xc = "Feature")
#' out
#'
#' @export aov_TukeyHSD
#'

aov_TukeyHSD <- function(df, xc = "Group", yc = "Intensity", op = NULL,
                         verbose = FALSE) {
    stopifnot(c(xc, yc) %in% colnames(df))
    if (verbose) {
        sink(file = paste0(op, "_TukeyHSD.txt"), append = TRUE, split = TRUE)
        cat(paste("Performing one-way ANOVA analysis for", op, "\n"))
    }
    formu <- as.formula(paste(yc, "~", xc))
    res.aov <- aov(formu, data = df)
    s <- summary(res.aov)
    p <- s[[1]][1, "Pr(>F)"]
    if (verbose) {
        print(s) ## anova summary
        cat("\nPost hoc Tukey Honest Significant Differences test\n")
    }
    v <- TukeyHSD(res.aov, which = xc)
    if (verbose) {
        print(v)
        sink()
    }
    invisible(list("ANOVA" = p, "HSD" = v[[xc]]))
}

#' @title Convert GRanges to dataframe
#' @description Convert a GRanges object with meta data columns to a dataframe,
#' with the first 6 columns corresponding those of BED6 format, and the meta
#' data as additional columns
#'
#' @param gr a GRanges object
#'
#' @return a dataframe
#'
#' @author Shuye Pu
#'
#' @examples
#' gr2 <- GenomicRanges::GRanges(c("chr1", "chr1"),
#'     IRanges::IRanges(c(7, 13), width = 3),
#'     strand = c("+", "-")
#' )
#' GenomicRanges::mcols(gr2) <- data.frame(
#'     score = c(0.3, 0.9),
#'     cat = c(TRUE, FALSE)
#' )
#' df2 <- gr2df(gr2)
#'
#' @export gr2df

gr2df <- function(gr) {
    chr <- as.vector(seqnames(gr))
    start <- start(gr) - 1 # convert to 0-based for bed
    end <- end(gr)
    width <- width(gr)
    strand <- as.vector(strand(gr))
    meta <- mcols(gr, us.names = TRUE)

    if ("name" %in% colnames(meta)) {
        name <- gr$name
        meta <- meta[, !colnames(meta) %in% c("name"), drop = FALSE]
    } else if ("id" %in% colnames(meta)) {
        name <- gr$id
        meta <- meta[, !colnames(meta) %in% c("id"), drop = FALSE]
    } else if (!is.null(names(gr))) {
        name <- names(gr)
    } else {
        name <- seq_along(gr)
    }

    if ("score" %in% colnames(meta)) {
        score <- gr$score
        meta <- meta[, !colnames(meta) %in% c("score"), drop = FALSE]
    } else {
        score <- width
    }

    df <- as.data.frame(cbind(chr, start, end, name, score, strand, meta))

    return(df)
}

#' @title compute ratio over input
#'
#' @description compute enrichment of IP samples over Input samples
#'
#' @param IP a numerical matrix
#' @param Input another numerical matrix with same dimensions as the IP matrix
#' @param verbose logical, whether to output additional information
#'
#' @return a numerical matrix with same dimensions as the IP matrix
#'
#' @author Shuye Pu
#'
#' @examples
#'
#' IP <- matrix(rlnorm(100), ncol = 10)
#' Input <- matrix(runif(100), ncol = 10)
#'
#' ratio <- GenomicPlot:::ratio_over_input(IP, Input, verbose = TRUE)
#'
#' @keywords internal
#'
ratio_over_input <- function(IP, Input, verbose = FALSE) {
    if (!identical(dim(IP), dim(Input))) stop("IP matrix and Input matrix must
                                              have same dimensions")
    if (min(IP) < 0 || min(Input) < 0) stop("IP matrix and Input matrix cannot
                                            have negative values")

    ## regularize the matrices to avoid unreasonably high ratios.
    ## The maximum of IP determines the size of the regularizing term which is
    ## a pseudo number. As a pseudo number is essentially a noise, if it is too
    ## large, it will mask signals.

    reg <- 0
    numerator <- max(IP)
    if (numerator > 10) {
        reg <- 1
    } else if (numerator > 1) {
        reg <- 0.1
    } else if (numerator > 0.1) {
        reg <- 0.01
    } else if (numerator > 0.01) {
        reg <- 0.001
    } else if (numerator > 0.001) {
        reg <- 0.0001
    }

    if (verbose) {
        message("Maximum if IP matrix: ", max(IP))
        message("Minimum of Input matrix: ", min(Input))
        message("Pseudo number added: ", reg)
    }

    ratio <- (IP + reg) / (Input + reg)

    return(ratio)
}
