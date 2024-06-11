#' @title Plot promoter, 5'UTR, CDS, 3'UTR and TTS
#'
#' @description Plot reads or peak Coverage/base/gene of samples given in the
#' query files around genes. The upstream and downstream windows flanking genes
#' can be given separately, metagene plots are generated with 5'UTR, CDS and
#' 3'UTR segments. The length of each segments are prorated according to the
#' median length of each segments. If Input files are provided, ratio over Input
#' is computed and displayed as well.
#'
#' @param queryFiles a vector of sample file names. The file should be in .bam,
#'  .bed, .wig or .bw format, mixture of formats is allowed
#' @param gFeatures_list a list of genomic features as output of the function
#'  \code{\link{prepare_5parts_genomic_features}}
#' @param inputFiles a vector of input sample file names. The file should be in
#'  .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param importParams a list of parameters for \code{handle_input}
#' @param stranded logical, indicating whether the strand of the feature should
#'  be considered
#' @param scale logical, indicating whether the score matrix should be scaled to
#'  the range 0:1, so that samples with different baseline can be compared
#' @param smooth logical, indicating whether the line should smoothed with a
#'  spline smoothing algorithm
#' @param heatmap logical, indicating whether a heatmap of the score matrix
#'  should be generated
#' @param heatRange a numeric vector with three elements, defining custom range for
#'  color ramp, default=NULL, i.e. the range is defined automatically based on
#'  the c(minimun, median, maximum) of a data matrix
#' @param rmOutlier a numeric value serving as a multiplier of the MAD in Hampel
#'  filter for outliers identification, 0 indicating not removing outliers.
#'  For Gaussian distribution, use 3, adjust based on data distribution.
#' @param transform logical, whether to log2 transform the matrix
#' @param outPrefix a string specifying output file prefix for plots
#'  (outPrefix.pdf)
#' @param verbose logical, indicating whether to output additional information
#'  (data used for plotting or statistical test results)
#' @param Ylab a string for y-axis label
#' @param hw a vector of two elements specifying the height and width of the
#'  output figures
#' @param nc integer, number of cores for parallel processing
#'
#' @return a dataframe containing the data used for plotting
#' @author Shuye Pu
#'
#' @examples
#'
#' data(gf5_meta)
#' queryfiles <- system.file("extdata", "treat_chr19.bam",
#'                           package = "GenomicPlot")
#' names(queryfiles) <- "clip_bam"
#' inputfiles <- system.file("extdata", "input_chr19.bam",
#'                           package = "GenomicPlot")
#' names(inputfiles) <- "clip_input"
#'
#' bamimportParams <- setImportParams(
#'     offset = -1, fix_width = 0, fix_point = "start", norm = TRUE,
#'     useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' plot_5parts_metagene(
#'     queryFiles = queryfiles,
#'     gFeatures_list = list("metagene" = gf5_meta),
#'     inputFiles = inputfiles,
#'     scale = FALSE,
#'     verbose = FALSE,
#'     transform = NA,
#'     smooth = TRUE,
#'     stranded = TRUE,
#'     outPrefix = NULL,
#'     importParams = bamimportParams,
#'     heatmap = TRUE,
#'     rmOutlier = 0,
#'     nc = 2
#' )
#'
#' @export plot_5parts_metagene
#'

plot_5parts_metagene <- function(queryFiles,
                                 gFeatures_list,
                                 inputFiles = NULL,
                                 importParams = NULL,
                                 verbose = FALSE,
                                 transform = NA,
                                 smooth = FALSE,
                                 scale = FALSE,
                                 stranded = TRUE,
                                 outPrefix = NULL,
                                 heatmap = FALSE,
                                 heatRange = NULL,
                                 rmOutlier = 0,
                                 Ylab = "Coverage/base/gene",
                                 hw = c(10, 10),
                                 nc = 2) {
    stopifnot(is.numeric(c(nc, hw, rmOutlier)))
    stopifnot(transform %in% c("log", "log2", "log10", NA))
    stopifnot(all(file.exists(queryFiles)))
    if (is.null(names(queryFiles)) || any(names(queryFiles) == ""))
        stop("Each file must have a name attribute!")

    if (verbose) message("[plot_5parts_metagene] started ...\n")
    functionName <- as.character(match.call()[[1]])
    params <- plot_named_list(as.list(environment()))
    force(params)

    if (!is.null(outPrefix)) {
        while (!is.null(dev.list())) {
            dev.off()
        }
        pdf(paste(outPrefix, "pdf", sep = "."), height = hw[1], width = hw[2])
    }

    if (is.null(inputFiles)) {
        inputLabels <- NULL
        queryInputs <- handle_input(inputFiles = queryFiles, importParams,
                                    verbose = verbose, nc = nc)
    } else {
        inputLabels <- names(inputFiles)
        queryLabels <- names(queryFiles)
        if (length(queryFiles) == length(inputFiles)) {
            queryInputs <- handle_input(inputFiles = c(queryFiles, inputFiles),
                                        importParams, verbose = verbose,
                                        nc = nc)
        } else if (length(inputFiles) == 1) {
            queryInputs <- handle_input(inputFiles = c(queryFiles, inputFiles),
                                        importParams, verbose = verbose,
                                        nc = nc)
            queryInputs <- queryInputs[c(queryLabels,
                                         rep(inputLabels,
                                             length(queryLabels)))]
            ## make each inputLabels unique
            inputLabels <- paste0(names(inputFiles), seq_along(queryFiles))
            names(queryInputs) <- c(queryLabels, inputLabels)
        } else {
            stop("the number of inputFiles must be 1 or equal to the number of
                 queryFiles!")
        }
    }
    queryLabels <- names(queryInputs)

    mplot_dfs <- NULL
    mplot_dfs_ratio <- NULL
    heatmap_list <- list()
    heatmap_list_ratio <- list()
    for (aFeature in names(gFeatures_list)) {
        if (verbose) message("Computing coverage for query files in ", aFeature)
        gFeatures <- gFeatures_list[[aFeature]]

        windowRs <- gFeatures$windowRs
        featureNames <- names(windowRs)
        if (verbose) {
            message("Number of features:\n")
            message(paste(vapply(windowRs, length, numeric(1)), collapse = " "))
        }

        nbins <- gFeatures$nbins
        scaled_bins <- gFeatures$scaled_bins
        meta <- gFeatures$meta
        fiveP <- gFeatures$fiveP
        threeP <- gFeatures$threeP

        if (verbose) {
            message("Number of scaled bins:\n")
            message(paste(scaled_bins, collapse = " "), "\n")
        }

        scoreMatrix_list <- list()

        for (queryLabel in queryLabels) {
            if (verbose) message(queryLabel)
            Input <- queryInputs[[queryLabel]]
            libsize <- Input$size
            queryRegions <- Input$query
            fileType <- Input$type
            weight_col <- Input$weight

            for (w in featureNames) {
                if (verbose) message("Feature name: ", w, "\n")
                windowR <- windowRs[[w]]
                bin_num <- scaled_bins[w]

                bin_op <- "mean"
                if (bin_num > 0) {

                    fullMatrix <- parallel_scoreMatrixBin(
                        queryRegions, windowR, bin_num, bin_op, weight_col,
                        stranded, nc = nc)
                    scoreMatrix_list[[queryLabel]][[w]] <- fullMatrix
                } else {
                    scoreMatrix_list[[queryLabel]][[w]] <- NULL
                }
            }
        }

        if (verbose) message("Preparing data for individual plotting...\n")

        ## x axis points for vlines that demarcate the genomic features
        vx <- c(1, cumsum(scaled_bins[seq_len((length(scaled_bins) - 1))]) + 1)
        names(vx) <- featureNames

        processed_matrix <- list()
        mplot_df <- NULL
        Ylab <- ifelse(!is.na(transform) && is.null(inputFiles),
                       paste0(transform, " (", Ylab, ")"), Ylab)

        for (queryLabel in queryLabels) {
            plot_df <- NULL
            dims <- vapply(scoreMatrix_list[[queryLabel]], dim, numeric(2))

            if (any(dims[1, ] != dims[1, 1])) {
                message(paste(dims[1, ], collapse = " "), "\n")
                stop("Number of genes are not equal among features, make sure
                     all feature windows are within chromosome lengths of query
                     regions, as genomation will remvove all feature windows
                     outside chromosome boundaries")
            } else {
                featureMatrix <- as.matrix(bind_cols(
                    scoreMatrix_list[[queryLabel]]))
                if (is.null(inputFiles)) {
                    featureMatrix <- process_scoreMatrix(featureMatrix, scale,
                                                         rmOutlier,
                                                         transform = transform,
                                                         verbose = verbose)
                } else {
                    featureMatrix <- process_scoreMatrix(featureMatrix,
                                                         scale = FALSE,
                                                         rmOutlier = rmOutlier,
                                                         transform = NA,
                                                         verbose = verbose)
                }
                processed_matrix[[queryLabel]] <- featureMatrix

                colm <- apply(featureMatrix, 2, mean)
                if(nrow(featureMatrix) == 1){
                    colsd <- rep(0, ncol(featureMatrix))
                }else{
                    colsd <- apply(featureMatrix, 2, sd)
                }
                colse <- colsd / sqrt(nrow(featureMatrix))
                querybed <- rep(queryLabel, ncol(featureMatrix))
                collabel <- list()
                featuretype <- list()
                for (w in featureNames) {
                    # if(verbose) print(w)
                    if (scaled_bins[w] > 0) {
                        bin_num <- scaled_bins[w]
                        collabel[[w]] <- seq(vx[w], vx[w] + bin_num - 1)
                        featuretype[[w]] <- rep(w, bin_num)
                    }
                }
                collabel <- unlist(collabel)
                featuretype <- unlist(featuretype)
                names(collabel) <- featuretype

                if (heatmap) {
                    dataname <- paste(Ylab, queryLabel, aFeature, sep = ":")
                    rgs <- NULL
                    if(is.null(inputFiles)) rgs <- heatRange
                    heatmap_list[[dataname]] <- draw_matrix_heatmap(
                        featureMatrix, dataName = dataname,
                        labels_col = collabel, levels_col = featureNames,
                        ranges = rgs, verbose = verbose)
                }
                plot_df <- data.frame("Intensity" = colm, "sd" = colsd,
                                      "se" = colse, "Position" = collabel,
                                      "Query" = paste(querybed,
                                                      aFeature, sep = ":"),
                                      "Feature" = featuretype)
            }

            if (smooth) {
                plot_df$Intensity <- as.vector(
                    smooth.spline(plot_df$Intensity,
                                  df = as.integer(nbins / 5))$y)
                plot_df$se <- as.vector(
                    smooth.spline(plot_df$se, df = as.integer(nbins / 5))$y)
            }

            mplot_df <- rbind(mplot_df, plot_df)
        }

        mplot_df <- mutate(mplot_df, lower = Intensity - se,
                           upper = Intensity + se)
        mplot_dfs <- rbind(mplot_dfs, mplot_df)

        ## if inputFiles are provided, plot ratio over input

        if (!is.null(inputFiles)) {
            if (verbose) message("Preparing data for ratio plotting...\n")
            Ylabr <- ifelse(is.na(transform), "Ratio-over-Input",
                            paste0(transform, " (Ratio-over-Input)"))

            ratiolabels <- queryLabels[!queryLabels %in% inputLabels]
            inputMatrix_list <- processed_matrix[inputLabels]
            ratioMatrix_list <- processed_matrix[ratiolabels]

            for (i in seq_along(ratiolabels)) {
                fullMatrix <- ratio_over_input(
                    ratioMatrix_list[[ratiolabels[i]]],
                    inputMatrix_list[[inputLabels[i]]], verbose)
                fullMatrix <- process_scoreMatrix(
                    fullMatrix, scale, rmOutlier, transform = transform,
                    verbose = verbose)

                ratioMatrix_list[[ratiolabels[i]]] <- fullMatrix
            }

            mplot_df <- NULL

            for (ratiolabel in ratiolabels) {
                plot_df <- NULL
                featureMatrix <- as.matrix(ratioMatrix_list[[ratiolabel]])

                colm <- apply(featureMatrix, 2, mean)
                if(nrow(featureMatrix) == 1){
                    colsd <- rep(0, ncol(featureMatrix))
                }else{
                    colsd <- apply(featureMatrix, 2, sd)
                }
                colse <- colsd / sqrt(nrow(featureMatrix))
                querybed <- rep(ratiolabel, ncol(featureMatrix))
                collabel <- list()
                featuretype <- list()
                for (w in featureNames) {
                    if (scaled_bins[w] > 0) {
                        bin_num <- scaled_bins[w]
                        collabel[[w]] <- seq(vx[w], vx[w] + bin_num - 1)
                        featuretype[[w]] <- rep(w, bin_num)
                    }
                }
                collabel <- unlist(collabel)
                featuretype <- unlist(featuretype)
                names(collabel) <- featuretype

                if (heatmap) {
                    dataname <- paste(Ylabr, ratiolabel, aFeature, sep = ":")
                    heatmap_list_ratio[[dataname]] <- draw_matrix_heatmap(
                        featureMatrix, dataName = dataname,
                        labels_col = collabel, levels_col = featureNames,
                        ranges = heatRange, verbose = verbose)
                }
                plot_df <- data.frame("Intensity" = colm, "sd" = colsd,
                                      "se" = colse, "Position" = collabel,
                                      "Query" = paste(querybed, aFeature,
                                                      sep = ":"),
                                      "Feature" = featuretype)

                if (smooth) {
                    plot_df$Intensity <- as.vector(
                        smooth.spline(plot_df$Intensity,
                                      df = as.integer(nbins / 5))$y)
                    plot_df$se <- as.vector(
                        smooth.spline(plot_df$se, df = as.integer(nbins / 5))$y)
                }

                mplot_df <- rbind(mplot_df, plot_df)
            }

            mplot_df <- mutate(mplot_df, lower = Intensity - se,
                               upper = Intensity + se)
            mplot_dfs_ratio <- rbind(mplot_dfs_ratio, mplot_df)
        }
    }

    if (verbose) message("Start plotting...\n")

    xmax <- max(mplot_dfs$Position)
    pp <- draw_region_landmark(featureNames, vx, xmax)
    ppp <- draw_region_name(featureNames, scaled_bins, xmax)

    ## plot individual sample lines with error band
    plot_list <- list()
    for (aFeature in names(gFeatures_list)) {
        for (queryLabel in queryLabels) {
            aplot_df <- mplot_dfs %>%
                filter(Query == paste(queryLabel, aFeature, sep = ":"))
            p <- draw_region_profile(plot_df = aplot_df, cn = "Query",
                                     vx = vx, Ylab = Ylab)
            outp <- plot_grid(p, pp, ppp, ncol = 1, align = "v", axis = "lr",
                              rel_heights = c(20, 1, 2.5))
            plot_list[[paste(queryLabel, aFeature, sep = ":")]] <- outp
        }
    }
    draw_stacked_plot(plot_list, heatmap_list)

    ## plot multi-sample lines with error band
    if (length(queryLabels) * length(gFeatures_list) > 1) {
        p <- draw_region_profile(plot_df = mplot_dfs, cn = "Query", vx = vx,
                                 Ylab = Ylab)
        outp <- plot_grid(p, pp, ppp, ncol = 1, align = "v", axis = "lr",
                          rel_heights = c(20, 1, 2.5))
        print(outp)
    }

    if (!is.null(inputFiles)) {
        xmax <- max(mplot_dfs_ratio$Position)
        pp <- draw_region_landmark(featureNames, vx, xmax)
        ppp <- draw_region_name(featureNames, scaled_bins, xmax)

        ## plot individual sample lines with error band
        plot_list <- list()
        for (aFeature in names(gFeatures_list)) {
            for (ratiolabel in ratiolabels) {
                aplot_df <- mplot_dfs_ratio %>%
                    filter(Query == paste(ratiolabel, aFeature, sep = ":"))
                p <- draw_region_profile(plot_df = aplot_df, cn = "Query",
                                         vx = vx, Ylab = Ylabr)
                outp <- plot_grid(p, pp, ppp, ncol = 1, align = "v",
                                  axis = "lr", rel_heights = c(20, 1, 2.5))
                plot_list[[paste(ratiolabel, aFeature, sep = ":")]] <- outp
            }
        }

        draw_stacked_plot(plot_list, heatmap_list_ratio)

        ## plot multi-sample lines with error band
        if (length(ratiolabels) * length(gFeatures_list) > 1) {
            p <- draw_region_profile(plot_df = mplot_dfs_ratio, cn = "Query",
                                     vx = vx, Ylab = Ylabr)
            outp <- plot_grid(p, pp, ppp, ncol = 1, align = "v", axis = "lr",
                              rel_heights = c(20, 1, 2.5))
            print(outp)
        }
    }

    if (!is.null(outPrefix)) {
        print(params)
        on.exit(dev.off(), add = TRUE)
        if (verbose) {
            message("[plot_5parts_metagene] runs successfully!\n")
        }
    }

    invisible(mplot_dfs)
}
