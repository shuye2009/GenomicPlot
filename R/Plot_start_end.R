#' @title  Plot signals around the start and the end of genomic features
#
#' @description   Plot reads or peak Coverage/base/gene of samples given in the
#' query files around start and end of custom features. The upstream and
#' downstream windows can be given separately, within the window, a smaller
#' window can be defined to highlight region of interest. A line plot will be
#' displayed for both start and end of feature. If Input files are provided,
#' ratio over Input is computed and displayed as well.
#'
#' @param queryFiles a vector of sample file names. The file should be in .bam,
#'  .bed, .wig or .bw format, mixture of formats is allowed
#' @param centerFiles  bed files that define the custom features, or features in
#'  c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene"), multiple
#'  features are allowed.
#' @param txdb a TxDb object defined in the GenomicFeatures package. Default
#'  NULL, needed only when genomic features are used in the place of centerFiles.
#' @param inputFiles a vector of input sample file names. The file should be in
#'  .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param importParams a list of parameters for \code{\link{handle_input}}
#' @param binSize an integer defines bin size for intensity calculation
#' @param ext a vector of four integers defining upstream and downstream
#'  boundaries of the plot window, flanking the start and end of features
#' @param hl a vector of four integers defining upstream and downstream
#'  boundaries of the highlight window, flanking the start and end of features
#' @param insert an integer specifies the length of the center regions to be
#'  included, in addition to the start and end of the feature
#' @param stranded logical, indicating whether the strand of the feature should
#'  be considered
#' @param scale logical, indicating whether the score matrix should be scaled to
#'  the range 0:1, so that samples with different baseline can be compared
#' @param smooth logical, indicating whether the line should smoothed with a
#'  spline smoothing algorithm
#' @param rmOutlier a numeric value serving as a multiplier of the MAD in Hampel
#'  filter for outliers identification, 0 indicating not removing outliers. For
#'  Gaussian distribution, use 3, adjust based on data distribution
#' @param outPrefix a string specifying output file prefix for plots
#' (outPrefix.pdf)
#' @param transform a string in c("log", "log2", "log10"), default = NA,
#'  indicating no transformation of data matrix
#' @param verbose logical, whether to output additional information (including
#'  data used for plotting or statistical test results)
#' @param Ylab a string for y-axis label
#' @param shade logical indicating whether to place a shaded rectangle around
#'  the point of interest
#' @param hw a vector of two elements specifying the height and width of the
#'  output figures
#' @param nc integer, number of cores for parallel processing
#'
#' @return a list of two objects, the first is a GRanges object, the second is
#'  a GRangesList object
#' @author Shuye Pu
#'
#' @examples
#'
#' gtfFile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#'     package = "GenomicPlot"
#' )
#'
#' txdb <- custom_TxDb_from_GTF(gtfFile, genome = "hg19")
#' bamQueryFiles <- system.file("extdata", "treat_chr19.bam",
#'     package = "GenomicPlot")
#' names(bamQueryFiles) <- "clip_bam"
#' bamInputFiles <- system.file("extdata", "input_chr19.bam",
#'                              package = "GenomicPlot")
#' names(bamInputFiles) <- "clip_input"
#'
#' bamimportParams <- setImportParams(
#'     offset = -1, fix_width = 0, fix_point = "start", norm = TRUE,
#'     useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' plot_start_end(
#'     queryFiles = bamQueryFiles,
#'     inputFiles = bamInputFiles,
#'     txdb = txdb,
#'     centerFiles = "intron",
#'     binSize = 10,
#'     importParams = bamimportParams,
#'     ext = c(-500, 200, -200, 500),
#'     hl = c(-100, 100, -100, 100),
#'     insert = 100,
#'     stranded = TRUE,
#'     scale = FALSE,
#'     smooth = TRUE,
#'     transform = "log2",
#'     outPrefix = NULL,
#'     nc = 2
#' )
#' @export plot_start_end
#'

plot_start_end <- function(queryFiles,
                           inputFiles = NULL,
                           centerFiles,
                           txdb = NULL,
                           importParams = NULL,
                           binSize = 10,
                           insert = 0,
                           verbose = FALSE,
                           ext = c(-500, 100, -100, 500),
                           hl = c(-50, 50, -50, 50),
                           stranded = TRUE,
                           scale = FALSE,
                           smooth = FALSE,
                           rmOutlier = 0,
                           outPrefix = NULL,
                           transform = NA,
                           shade = TRUE,
                           Ylab = "Coverage/base/gene",
                           hw = c(8, 8),
                           nc = 2) {
    stopifnot(is.numeric(c(binSize, insert, ext, hl, nc, hw, rmOutlier)))
    stopifnot(transform %in% c("log", "log2", "log10", NA))
    stopifnot(all(file.exists(queryFiles)))
    if (is.null(names(queryFiles)) || any(names(queryFiles) == ""))
        stop("Each file must have a name attribute!")

    if (verbose) message("[plot_start_end] started ...\n")
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
            queryInputs <- queryInputs[c(queryLabels, rep(inputLabels,
                                                          length(queryLabels)))]
            ## expand the list

            inputLabels <- paste0(names(inputFiles), seq_along(queryFiles))
            ## make each inputLabels unique
            names(queryInputs) <- c(queryLabels, inputLabels)
        } else {
            stop("the number of inputFiles must be 1 or equal to the number of
                 queryFiles!")
        }
    }
    queryLabels <- names(queryInputs)

    bedparam <- importParams
    bedparam$offset <- 0
    bedparam$fix_width <- 0
    bedparam$fix_point <- "start"
    bedparam$useScore <- FALSE
    bedparam$outRle <- FALSE
    bedparam$norm <- FALSE
    bedparam$useSizeFactor <- FALSE

    features <- list()
    minimal_width <- ext[2] - ext[3] + insert

    for (featureName in centerFiles) {
        if (featureName %in% c("utr3", "utr5", "cds", "intron", "exon",
                               "transcript", "gene")) {
            featureGR <- get_genomic_feature_coordinates(txdb, featureName,
                                                         longest = TRUE,
                                            protein_coding = TRUE)[["GRanges"]]
            featureGR <- featureGR[width(featureGR) > minimal_width]
            feature <- list("query" = featureGR)
            features[[featureName]] <- feature
        } else if (file.exists(featureName)) {
            names(featureName) <- names(centerFiles)[centerFiles == featureName]
            feature <- handle_input(featureName, bedparam, verbose = verbose,
                                    nc = nc)
            featureGR <- feature[[1]]$query
            featureGR <- featureGR[width(featureGR) > minimal_width]
            feature[[1]]$query <- featureGR
            features[[names(feature)[1]]] <- feature[[1]]
        } else {
            stop("featureName is not supported!")
        }
    }

    featureNames <- names(features)

    ext[2] <- ext[2] - (ext[2] - ext[1]) %% binSize
    ## to avoid binSize inconsistency, as the final binSize depends on bin_num
    bin_num_s <- round((ext[2] - ext[1]) / binSize)
    ext[4] <- ext[4] - (ext[4] - ext[3]) %% binSize
    bin_num_e <- round((ext[4] - ext[3]) / binSize)
    bin_num_c <- round(insert / binSize)

    if (verbose) message("Bin numbers: ", bin_num_s, " ", bin_num_c, " ",
                         bin_num_e, "\n")
    scoreMatrix_lists <- list()
    mat_lists <- list()
    plot_df <- NULL
    for (featureName in featureNames) {
        feature <- features[[featureName]]$query
        nf <- length(feature)
        if (verbose) message("Number of features: ", featureName, " ", nf, "\n")

        fs <- promoters(resize(feature, width = 1, fix = "start"),
                        upstream = -ext[1], downstream = ext[2])
        fe <- promoters(resize(feature, width = 1, fix = "end"),
                        upstream = -ext[3], downstream = ext[4])
        fc <- promoters(resize(feature, width = 1, fix = "center"),
                        upstream = round(insert / 2),
                        downstream = round(insert / 2))

        mat_list <- list()
        mat_list[["Start"]] <- list("window" = fs, s = ext[1], e = ext[2],
                                    "xmin" = hl[1], "xmax" = hl[2],
                                    "bin_num" = bin_num_s)
        mat_list[["Center"]] <- list("window" = fc, s = -round(insert / 2),
                                     e = round(insert / 2), "xmin" = 0,
                                     "xmax" = 0, "bin_num" = bin_num_c)
        mat_list[["End"]] <- list("window" = fe, s = ext[3], e = ext[4],
                                  "xmin" = hl[3], "xmax" = hl[4],
                                  "bin_num" = bin_num_e)

        mat_lists[[featureName]] <- mat_list

        scoreMatrix_list <- list()

        for (locus in names(mat_list)) {
            windowR <- mat_list[[locus]]$window
            bin_num <- mat_list[[locus]]$bin_num
            if (verbose) message("Locus: ", locus, " ", bin_num, " ",
                                 length(windowR), "\n")
            if (bin_num <= 0) next

            for (queryLabel in queryLabels) {
                if (verbose) message("Query label: ", queryLabel, "\n")
                queryRegions <- queryInputs[[queryLabel]]$query
                libsize <- queryInputs[[queryLabel]]$size

                bin_op <- "mean"
                weight_col <- queryInputs[[queryLabel]]$weight

                fullMatrix <- parallel_scoreMatrixBin(queryRegions, windowR,
                                                      bin_num, bin_op,
                                                      weight_col, stranded,
                                                      nc = nc)
                if (is.null(inputFiles)) {
                    fullMatrix <- process_scoreMatrix(fullMatrix, scale,
                                                      rmOutlier,
                                                      transform = transform,
                                                      verbose = verbose)
                } else {
                    fullMatrix <- process_scoreMatrix(fullMatrix, scale = FALSE,
                                                      rmOutlier = rmOutlier,
                                                      transform = NA,
                                                      verbose = verbose)
                }

                scoreMatrix_list[[queryLabel]][[locus]] <- fullMatrix
            }
        }

        scoreMatrix_lists[[featureName]] <- scoreMatrix_list

        for (locus in names(mat_list)) {
            xmin <- mat_list[[locus]]$xmin
            xmax <- mat_list[[locus]]$xmax
            bin_num <- mat_list[[locus]]$bin_num
            start <- mat_list[[locus]]$s
            end <- mat_list[[locus]]$e

            if (bin_num <= 0) next
            for (queryLabel in queryLabels) {
                if (verbose) message("Query label: ", queryLabel, "\n")

                fullMatrix <- scoreMatrix_list[[queryLabel]][[locus]]

                colm <- apply(fullMatrix, 2, mean)
                if(nrow(fullMatrix) == 1){
                    colsd <- rep(0, ncol(fullMatrix))
                }else{
                    colsd <- apply(fullMatrix, 2, sd)
                }
                colse <- colsd / sqrt(nrow(fullMatrix))
                collabel <- seq(start, (end - binSize), binSize)
                querybed <- as.factor(rep(queryLabel, ncol(fullMatrix)))
                location <- as.factor(rep(locus, ncol(fullMatrix)))
                levels(location) <- rev(levels(location))
                featurename <- as.factor(rep(featureName, ncol(fullMatrix)))
                Xmin <- rep(xmin, ncol(fullMatrix))
                Xmax <- rep(xmax, ncol(fullMatrix))
                halfmin <- min(fullMatrix)
                intervals <- apply(fullMatrix, 2, function(x)
                    length(x[x > halfmin]))

                sub_df <- NULL
                sub_df <- data.frame("Intensity" = colm, "sd" = colsd,
                                     "se" = colse, "Interval" = intervals,
                                     "Position" = collabel, "Query" = querybed,
                                     "Location" = location,
                                     "Feature" = featurename)
                if (smooth) {
                    sub_df$Intensity <- as.vector(smooth.spline(
                        sub_df$Intensity, df = as.integer(bin_num / 5))$y)
                    sub_df$se <- as.vector(smooth.spline(
                        sub_df$se, df = as.integer(bin_num / 5))$y)
                    sub_df$Interval <- as.vector(smooth.spline(
                        sub_df$Interval, df = as.integer(bin_num / 5))$y)
                }
                sub_df <- mutate(sub_df, lower = Intensity - se,
                                 upper = Intensity + se)

                plot_df <- rbind(plot_df, sub_df)
            }
        }
    }

    Ylab <- ifelse(!is.na(transform) && is.null(inputFiles),
                   paste0(transform, " (", Ylab, ")"), Ylab)
    ## plot multi feature lines for one query
    for (query in unique(plot_df$Query)) {
        qplot_df <- plot_df %>%
            filter(Query == query)

        plots <- draw_stacked_profile(plot_df = qplot_df, cn = "Feature",
                                      ext = ext, hl = hl, atitle = query,
                                      insert = insert, Ylab = Ylab,
                                      shade = shade)

        print(plots)
    }
    ## plot multi query lines for one feature
    for (feature in unique(plot_df$Feature)) {
        fplot_df <- plot_df %>%
            filter(Feature == feature)

        plots <- draw_stacked_profile(plot_df = fplot_df, cn = "Query",
                                      ext = ext, hl = hl, atitle = feature,
                                      insert = insert, Ylab = Ylab,
                                      shade = shade)

        print(plots)
    }

    ## compute and plot ratio over input
    if (!is.null(inputFiles)) {
        Ylab <- ifelse(is.na(transform), "Ratio-over-Input",
                       paste0(transform, " (Ratio-over-Input)"))

        plot_df <- NULL
        for (featureName in featureNames) {
            feature <- features[[featureName]]$query

            scoreMatrix_list <- scoreMatrix_lists[[featureName]]
            mat_list <- mat_lists[[featureName]]

            ratiolabels <- queryLabels[!queryLabels %in% inputLabels]
            inputMatrix_list <- scoreMatrix_list[inputLabels]
            ratioMatrix_list <- scoreMatrix_list[ratiolabels]


            for (locus in names(mat_list)) {
                bin_num <- mat_list[[locus]]$bin_num
                if (bin_num <= 0) next
                for (i in seq_along(ratiolabels)) {
                    rm <- ratioMatrix_list[[ratiolabels[i]]][[locus]]
                    im <- inputMatrix_list[[inputLabels[i]]][[locus]]
                    minrow <- min(nrow(rm), nrow(im))

                    fullMatrix <- ratio_over_input(rm[seq_len(minrow), ],
                                                   im[seq_len(minrow), ],
                                                   verbose)
                    fullMatrix <- process_scoreMatrix(fullMatrix, scale,
                                                      rmOutlier, transform,
                                                      verbose = verbose)

                    ratioMatrix_list[[ratiolabels[i]]][[locus]] <- fullMatrix
                }
            }

            for (locus in names(mat_list)) {
                xmin <- mat_list[[locus]]$xmin
                xmax <- mat_list[[locus]]$xmax
                bin_num <- mat_list[[locus]]$bin_num
                start <- mat_list[[locus]]$s
                end <- mat_list[[locus]]$e

                if (bin_num <= 0) next
                for (ratiolabel in ratiolabels) {
                    if (verbose) message("Ratio label: ", ratiolabel, "\n")

                    fullMatrix <- ratioMatrix_list[[ratiolabel]][[locus]]

                    colm <- apply(fullMatrix, 2, mean)
                    if(nrow(fullMatrix) == 1){
                        colsd <- rep(0, ncol(fullMatrix))
                    }else{
                        colsd <- apply(fullMatrix, 2, sd)
                    }
                    colse <- colsd / sqrt(nrow(fullMatrix))
                    collabel <- seq(start, (end - binSize), binSize)
                    ratiobed <- as.factor(rep(ratiolabel, ncol(fullMatrix)))
                    location <- as.factor(rep(locus, ncol(fullMatrix)))
                    featurename <- as.factor(rep(featureName, ncol(fullMatrix)))
                    levels(location) <- rev(levels(location))
                    halfmin <- min(fullMatrix)

                    intervals <- apply(fullMatrix, 2, function(x)
                        length(x[x > halfmin]))

                    sub_df <- NULL
                    sub_df <- data.frame("Intensity" = colm, "sd" = colsd,
                                         "se" = colse, "Interval" = intervals,
                                         "Position" = collabel,
                                         "Query" = ratiobed,
                                         "Location" = location,
                                         "Feature" = featurename)
                    if (smooth) {
                        sub_df$Intensity <- as.vector(smooth.spline(
                            sub_df$Intensity, df = as.integer(bin_num / 5))$y)
                        sub_df$se <- as.vector(smooth.spline(sub_df$se,
                                                df = as.integer(bin_num / 5))$y)
                        sub_df$Interval <- as.vector(smooth.spline(
                            sub_df$Interval, df = as.integer(bin_num / 5))$y)
                    }
                    sub_df <- mutate(sub_df, lower = Intensity - se,
                                     upper = Intensity + se)
                    plot_df <- rbind(plot_df, sub_df)
                }
            }
        }

        ## plot multi feature lines for one query
        for (query in unique(plot_df$Query)) {
            qplot_df <- plot_df %>%
                filter(Query == query)

            plots <- draw_stacked_profile(plot_df = qplot_df, cn = "Feature",
                                          ext = ext, hl = hl, atitle = query,
                                          insert = insert, Ylab = Ylab,
                                          shade = shade)

            print(plots)
        }
        ## plot multi query lines for one feature
        for (feature in unique(plot_df$Feature)) {
            fplot_df <- plot_df %>%
                filter(Feature == feature)

            plots <- draw_stacked_profile(plot_df = fplot_df, cn = "Query",
                                          ext = ext, hl = hl, atitle = feature,
                                          insert = insert, Ylab = Ylab,
                                          shade = shade)

            print(plots)
        }
    }

    if (!is.null(outPrefix)) {
        print(params)
        on.exit(dev.off(), add = TRUE)
    }

    if (verbose) message("[plot_start_end] finished!\n")
    invisible(plot_df)
}
