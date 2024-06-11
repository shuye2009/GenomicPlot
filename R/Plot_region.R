#' @title Plot signal inside as well as around custom genomic regions
#'
#' @description Plot reads or peak Coverage/base/gene of samples given in the
#' query files inside regions defined in the centerFiles. The upstream and
#' downstream flanking windows can be given separately. If Input files are
#' provided, ratio over Input is computed and displayed as well.
#'
#' @param queryFiles a named vector of sample file names. The file should be in
#'  .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param centerFiles a named vector of reference file names or genomic features
#'  in  c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene"). The
#'  file should be in .bed format only
#' @param txdb a TxDb object defined in the GenomicFeatures package. Default
#'  NULL, needed only when genomic features are used as centerFiles.
#' @param inputFiles a named vector of input sample file names. The file should
#'  be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param nbins an integer defines the total number of bins
#' @param fiveP an integer, indicating extension out or inside of the 5'
#'  boundary of gene by negative or positive number
#' @param threeP an integer, indicating extension out or inside of the 5'
#'  boundary of gene by positive or negative number
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
#'  filter for outliers identification, 0 indicating not removing outliers. For
#'  Gaussian distribution, use 3, adjust based on data distribution
#' @param importParams a list of parameters for \code{handle_input}
#' @param outPrefix a string specifying output file prefix for plots
#'  (outPrefix.pdf)
#' @param regionName a string specifying the name of the center region in the
#'  plots
#' @param Ylab a string for y-axis label
#' @param transform a string in c("log", "log2", "log10"), default = NA
#'  indicating no transformation of data matrix
#' @param statsMethod a string in c("wilcox.test", "t.test"), for pair-wise
#'  group comparisons
#' @param verbose logical, indicating whether to output additional information
#'  (data used for plotting or statistical test results)
#' @param hw a vector of two elements specifying the height and width of the
#'  output figures
#' @param nc integer, number of cores for parallel processing
#'
#' @return a dataframe containing the data used for plotting
#' @author Shuye Pu
#'
#' @examples
#' centerfiles <- system.file("extdata", "test_chip_peak_chr19.narrowPeak",
#' package = "GenomicPlot")
#' names(centerfiles) <- c("NarrowPeak")
#' queryfiles <- c(
#'   system.file("extdata", "chip_treat_chr19.bam", package = "GenomicPlot"))
#' names(queryfiles) <- c("chip_bam")
#' inputfiles <- c(
#'   system.file("extdata", "chip_input_chr19.bam", package = "GenomicPlot"))
#' names(inputfiles) <- c("chip_input")
#'
#' chipimportParams <- setImportParams(
#'   offset = 0, fix_width = 150, fix_point = "start", norm = TRUE,
#'   useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19",
#'   chr = c("chr19"))
#'
#' plot_region(
#'   queryFiles = queryfiles,
#'   centerFiles = centerfiles,
#'   inputFiles = inputfiles,
#'   nbins = 100,
#'   heatmap = TRUE,
#'   scale = FALSE,
#'   regionName = "narrowPeak",
#'   importParams = chipimportParams,
#'   verbose = FALSE,
#'   fiveP = -500,
#'   threeP = 500,
#'   smooth = TRUE,
#'   transform = "log2",
#'   stranded = TRUE,
#'   outPrefix = NULL,
#'   Ylab = "Coverage/base/peak",
#'   rmOutlier = 0,
#'   nc = 2
#' )
#'
#' @export plot_region

plot_region <- function(queryFiles,
                        centerFiles,
                        txdb = NULL,
                        regionName = "region",
                        inputFiles = NULL,
                        nbins = 100,
                        importParams = NULL,
                        verbose = FALSE,
                        scale = FALSE,
                        heatmap = FALSE,
                        fiveP = -1000,
                        threeP = 1000,
                        smooth = FALSE,
                        stranded = TRUE,
                        transform = NA,
                        outPrefix = NULL,
                        rmOutlier = 0,
                        heatRange = NULL,
                        Ylab = "Coverage/base/gene",
                        statsMethod = "wilcox.test",
                        hw = c(8, 8),
                        nc = 2) {
    stopifnot(is.numeric(c(nbins, fiveP, threeP, nc, hw, rmOutlier)))
    stopifnot(transform %in% c("log", "log2", "log10", NA))
    stopifnot(all(file.exists(queryFiles)))
    if (is.null(names(queryFiles)) || any(names(queryFiles) == ""))
        stop("Each file must have a name attribute!")

    if (verbose) message("[plot_region] started ...\n")
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
    bedparam$useScore <- FALSE
    bedparam$outRle <- FALSE
    bedparam$useSizeFactor <- FALSE

    centerInputs <- list()

    for (featureName in centerFiles) {
        if (featureName %in% c("utr3", "utr5", "cds", "intron", "exon",
                               "transcript", "gene")) {
            featureGR <- get_genomic_feature_coordinates(txdb, featureName,
                           longest = TRUE, protein_coding = TRUE)[["GRanges"]]
            feature <- list("query" = featureGR)
            centerInputs[[featureName]] <- feature
        } else if (file.exists(featureName)) {
            names(featureName) <- names(centerFiles)[centerFiles == featureName]
            feature <- handle_input(featureName, bedparam, verbose = verbose,
                                    nc = nc)
            centerInputs[[names(feature)[1]]] <- feature[[1]]
        } else {
            stop("featureName is not supported!")
        }
    }
    centerLabels <- names(centerInputs)

    five <- fiveP / 1000
    five <- paste0(five, "Kb")
    fiveL <- -fiveP
    if (fiveP >= 0) {
        fiveL <- 0
    }
    three <- threeP / 1000
    three <- paste0(three, "Kb")
    threeL <- threeP
    if (threeP <= 0) {
        threeL <- 0
    }

    featureNames <- c(five, regionName, three)
    # featureNames <- featureNames[featureNames != ""]

    all_regions <- unlist(as(lapply(centerInputs, function(x) x$query),
                             "GRangesList"))
    regionLen <- median(width(all_regions))
    lens <- c(fiveL, regionLen, threeL)

    scaled_bins <- round(lens * nbins / sum(lens))
    # scaled_bins <- scaled_bins[scaled_bins > 0]
    names(scaled_bins) <- featureNames

    if (verbose) message("Computing coverage for sample...\n")
    sml <- list() # scoreMatrix_list

    for (queryLabel in queryLabels) {
        if (verbose) message("Query label: ", queryLabel, "\n")

        myInput <- queryInputs[[queryLabel]]
        libsize <- myInput$size
        queryRegions <- myInput$query
        fileType <- myInput$type
        weight_col <- myInput$weight

        for (centerLabel in centerLabels) {
            if (verbose) message("Center label: ", centerLabel, "\n")
            centerInput <- centerInputs[[centerLabel]]
            centerGr <- centerInput$query
            centerGr <- check_constraints(centerGr, importParams$genome,
                                          queryRegions)

            if (fiveP > 0 && threeP < 0) {
               # to avoid generating negative width in the narrow function below
                centerGr <- centerGr[width(centerGr) > (abs(fiveP) +
                                                           abs(threeP))]
            } else if (fiveP > 0) {
                centerGr <- centerGr[width(centerGr) > abs(fiveP)]
            } else if (threeP < 0) {
                centerGr <- centerGr[width(centerGr) > abs(threeP)]
            }
            centerGrPlus <- centerGr[as.vector(strand(centerGr)) %in%
                                        c("+", "*")]
            centerGrMinus <- centerGr[as.vector(strand(centerGr)) == "-"]

            upstreamGr <- flank(centerGr, width = fiveL, start = TRUE,
                                both = FALSE, use.names = TRUE,
                                ignore.strand = !stranded)
            downstreamGr <- flank(centerGr, width = threeL, start = FALSE,
                                  both = FALSE, use.names = TRUE,
                                  ignore.strand = !stranded)
            upstreamGr <- check_constraints(upstreamGr, importParams$genome,
                                            queryRegions)
            downstreamGr <- check_constraints(downstreamGr, importParams$genome,
                                              queryRegions)

            if (fiveP > 0) {
                centerGrPlus <- narrow(centerGrPlus, start = fiveP, end = NA)
                centerGrMinus <- narrow(centerGrMinus, start = NA, end = fiveP)
            }
            if (threeP < 0) {
                centerGrPlus <- narrow(centerGrPlus, start = NA, end = threeP)
                centerGrMinus <- narrow(centerGrMinus, start = -threeP,
                                        end = NA)
            }
            centerGr <- c(centerGrPlus, centerGrMinus)
            centerGr <- centerGr[width(centerGr) >= scaled_bins[regionName]]


            if (fiveP < 0 && threeP > 0) {
                commonNames <- Reduce(intersect, list(
                     names(upstreamGr), names(downstreamGr),names(centerGr)))

                upstreamGr <- upstreamGr[commonNames]
                downstreamGr <- downstreamGr[commonNames]
                centerGr <- centerGr[commonNames]
            }

            windowRegions <- split(centerGr, f = factor(names(centerGr)))
            windowUp <- split(upstreamGr, f = factor(names(upstreamGr)))
            windowDown <- split(downstreamGr, f = factor(names(downstreamGr)))
            windowRs <- list(windowUp, windowRegions, windowDown)
            names(windowRs) <- featureNames

            for (w in featureNames) {
                # w <- "utr5"
                if (verbose) message("Feature name: ", w, "\n")
                if (scaled_bins[w] > 0) {
                    windowR <- as(windowRs[[w]], "GRangesList")
                    bin_num <- scaled_bins[w]
                    bin_op <- "mean"
                    fullMatrix <- parallel_scoreMatrixBin(queryRegions, windowR,
                                                          bin_num, bin_op,
                                                          weight_col, stranded,
                                                          nc = nc)
                    rownames(fullMatrix) <- names(windowR)
                    sml[[queryLabel]][[centerLabel]][[w]] <- fullMatrix
                }
            }
        }
    }

    mplot_df <- list()
    stat_df <- list()
    vx <- c(1, cumsum(scaled_bins[seq_len(length(scaled_bins) - 1)]) + 1)
    ## x-axis points for vlines that demarcate the genomic features
    names(vx) <- featureNames

    heatmap_list <- list()
    Ylab <- ifelse(!is.na(transform) && is.null(inputFiles),
                   paste0(transform, " (", Ylab, ")"), Ylab)

    processed_matrix <- list()
    processed_region_matrix <- list()
    if (verbose) message("Plotting coverage profiles...\n")
    for (queryLabel in queryLabels) {
        if (verbose) message("Query label: ", queryLabel, "\n")
        for (centerLabel in centerLabels) {
            if (verbose) message("Center label: ", centerLabel, "\n")
            plot_df <- NULL

            dims <- vapply(sml[[queryLabel]][[centerLabel]], dim,
                           numeric(2))

            if (any(dims[1, ] != dims[1, 1])) {
                message(paste(dims[1, ], collapse = " "), "\n")
                stop("Number of genes are not equal among features, make sure
                     all feature windows are within chromosome lengths of query
                     regions, as genomation will remvove all feature windows
                     outside chromosome boundaries")
            } else {
                featureMatrix <- as.matrix(
                   bind_cols(sml[[queryLabel]][[centerLabel]]))
                rownames(featureMatrix) <- rownames(
                   sml[[queryLabel]][[centerLabel]][[1]])
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
                processed_matrix[[queryLabel]][[centerLabel]] <- featureMatrix

                colm <- apply(featureMatrix, 2, mean)
                if(nrow(featureMatrix) == 1){
                    colsd <- rep(0, ncol(featureMatrix))
                }else{
                    colsd <- apply(featureMatrix, 2, sd)
                }
                colse <- colsd / sqrt(nrow(featureMatrix))
                querybed <- rep(queryLabel, ncol(featureMatrix))
                centerbed <- rep(centerLabel, ncol(featureMatrix))
                collabel_list <- list()
                featuretype <- list()
                for (w in featureNames) {
                    if (verbose) message("Feature name: ", w, "\n")
                    if (scaled_bins[w] > 0) {
                        bin_num <- scaled_bins[w]
                        collabel_list[[w]] <- seq(vx[w], vx[w] + bin_num - 1)
                        featuretype[[w]] <- rep(w, bin_num)
                    }
                }
                collabel <- unlist(collabel_list)
                featuretype <- unlist(featuretype)
                names(collabel) <- featuretype

                if (heatmap) {
                    dataname <- paste(Ylab, queryLabel, centerLabel, sep = ":")

                    heatmap_list[[dataname]] <- draw_matrix_heatmap(
                       featureMatrix, dataName = dataname,
                       labels_col = collabel,
                       levels_col = names(scaled_bins[scaled_bins > 0]),
                       ranges = heatRange, verbose = verbose)
                }

                plot_df <- data.frame("Intensity" = colm, "sd" = colsd,
                                      "se" = colse, "Position" = collabel,
                                      "Query" = querybed,
                                      "Reference" = centerbed,
                                      "Feature" = featuretype)
            }

            if (smooth) {
                plot_df$Intensity <- as.vector(smooth.spline(
                      plot_df$Intensity, df = as.integer(nbins / 5))$y)
                plot_df$se <- as.vector(smooth.spline(
                   plot_df$se, df = as.integer(nbins / 5))$y)
            }

            mplot_df[[paste(queryLabel, centerLabel, sep = ":")]] <- plot_df

            regionMatrix <- featureMatrix[, collabel_list[[regionName]], drop = FALSE]

            processed_region_matrix[[queryLabel]][[centerLabel]] <- regionMatrix

            Intensity <- as.numeric(rowMeans(regionMatrix))

            Query <- as.factor(rep(queryLabel, length(Intensity)))
            Reference <- as.factor(rep(centerLabel, length(Intensity)))
            subdf <- data.frame(Intensity, Query, Reference)
            stat_df[[paste(queryLabel, centerLabel, sep = ":")]] <- subdf
        }
    }

    mplot_df <- bind_rows(mplot_df) %>%
        mutate(Group = as.factor(paste0(Query, ":", Reference)),
               .keep = "all") %>%
        mutate(lower = Intensity - se, upper = Intensity + se)
    mstat_df <- bind_rows(stat_df) %>%
        mutate(Group = as.factor(paste0(Query, ":", Reference)), .keep = "all")

    plot_list <- list()
    xmax <- max(mplot_df$Position)
    pp <- draw_region_landmark(featureNames, vx, xmax)
    ppp <- draw_region_name(featureNames, scaled_bins, xmax)
    marker <- plot_grid(pp, ppp, ncol = 1, align = "v", axis = "lr",
                        rel_heights = c(1, 3))

    for (i in seq_along(queryLabels)) {
        for (beds in combn(queryLabels, i, simplify = FALSE)) {
            for (j in seq_along(centerLabels)) {
                for (centers in combn(centerLabels, j, simplify = FALSE)) {
                    if (verbose) {
                        message("Beds: ", paste(beds, collapse = " "), "\n")
                        message("Centers: ",
                                paste(centers, collapse = " "), "\n")
                    }

                    aplot_df <- mplot_df %>%
                        filter(Query %in% beds & Reference %in% centers) %>%
                        mutate(Group = paste(Query, Reference, sep = ":"),
                               .keep = "all")

                    ## plot multi-sample lines with error band
                    p <- draw_region_profile(plot_df = aplot_df, cn = "Group",
                                             vx = vx, Ylab = Ylab)
                    outp <- plot_grid(p, marker, ncol = 1, align = "v",
                                      axis = "lr", rel_heights = c(10, 2))

                    astat_df <- mstat_df %>%
                        filter(Query %in% beds & Reference %in% centers)
                    ## unify order of factors to get consistent color mapping
                    astat_df <- astat_df %>%
                        mutate(Query = factor(
                           Query,levels = sort(unique(Query)))) %>%
                        mutate(Reference = factor(
                           Reference, levels = sort(unique(Reference))))


                    if ((i == 1 && j > 1) || (i > 1 && j == 1)) {
                        if (j > 1) {
                            p <- draw_region_profile(plot_df = aplot_df,
                                                     cn = "Reference",
                                                     sn = "Query", vx = vx,
                                                     Ylab = Ylab)
                            outp <- plot_grid(p, marker, ncol = 1, align = "v",
                                              axis = "lr",
                                              rel_heights = c(10, 2))

                            comp <- combn(seq_along(
                               centers), 2, simplify = FALSE)
                            combo <- draw_combo_plot(stat_df = astat_df,
                                                     xc = "Reference",
                                                     yc = "Intensity",
                                                     comp = comp,
                                                     stats = statsMethod,
                                                     Ylab = Ylab)
                        } else {
                            p <- draw_region_profile(plot_df = aplot_df,
                                                     cn = "Query",
                                                     sn = "Reference",
                                                     vx = vx, Ylab = Ylab)
                            outp <- plot_grid(p, marker, ncol = 1, align = "v",
                                              axis = "lr",
                                              rel_heights = c(10, 2))
                            comp <- combn(seq_along(beds), 2, simplify = FALSE)

                            combo <- draw_combo_plot(stat_df = astat_df,
                                                     xc = "Query",
                                                     yc = "Intensity",
                                                     comp = comp,
                                                     stats = statsMethod,
                                                     Ylab = Ylab)
                        }

                        print(outp)
                        print(combo)
                    } else if (i == 1 && j == 1) {
                        plot_list[[paste(Ylab, beds,
                                         centers, sep = ":")]] <- outp
                    } else if (i == length(queryLabels) &&
                               j == length(centerLabels)) {
                        print(outp)
                    }
                }
            }
        }
    }

    draw_stacked_plot(plot_list, heatmap_list)

    if (!is.null(inputFiles)) {
        Ylab <- ifelse(is.na(transform), "Ratio-over-Input",
                       paste0(transform, " (Ratio-over-Input)"))

        ratiolabels <- queryLabels[!queryLabels %in% inputLabels]
        inputMatrix_list <- processed_matrix[inputLabels]
        ratioMatrix_list <- processed_matrix[ratiolabels]

        for (centerLabel in centerLabels) {
            for (i in seq_along(ratiolabels)) {
                rm <- ratioMatrix_list[[ratiolabels[i]]][[centerLabel]]
                im <- inputMatrix_list[[inputLabels[i]]][[centerLabel]]
                commonrow <- intersect(rownames(rm), rownames(im))

                fullMatrix <- ratio_over_input(rm[commonrow, ], im[commonrow, ],
                                               verbose)
                fullMatrix <- process_scoreMatrix(fullMatrix, scale, rmOutlier,
                                                  transform = transform,
                                                  verbose = verbose)

                ratioMatrix_list[[ratiolabels[i]]][[centerLabel]] <- fullMatrix
            }
        }

        imrl <- processed_region_matrix[inputLabels] # inputMatrix_region_list
        rmrl <- processed_region_matrix[ratiolabels] # ratioMatrix_region_list
        for (centerLabel in centerLabels) {
            for (i in seq_along(ratiolabels)) {
                rm <- rmrl[[ratiolabels[i]]][[centerLabel]]
                im <- imrl[[inputLabels[i]]][[centerLabel]]
                commonrow <- intersect(rownames(rm), rownames(im))

                fullMatrix <- ratio_over_input(rm[commonrow, ], im[commonrow, ],
                                               verbose)
                fullMatrix <- process_scoreMatrix(fullMatrix, scale, rmOutlier,
                                                  transform = transform,
                                                  verbose = verbose)

                rmrl[[ratiolabels[i]]][[centerLabel]] <- fullMatrix
            }
        }

        mplot_df <- list()
        stat_df <- list()
        heatmap_list <- list()
        if (verbose) message("Plotting coverage profiles...\n")
        for (ratiolabel in ratiolabels) {
            if (verbose) message("Ratio label: ", ratiolabel, "\n")
            for (centerLabel in centerLabels) {
                if (verbose) message("Center label: ", centerLabel, "\n")
                plot_df <- NULL

                featureMatrix <- ratioMatrix_list[[ratiolabel]][[centerLabel]]

                colm <- apply(featureMatrix, 2, mean)
                if(nrow(featureMatrix) == 1){
                    colsd <- rep(0, ncol(featureMatrix))
                }else{
                    colsd <- apply(featureMatrix, 2, sd)
                }
                colse <- colsd / sqrt(nrow(featureMatrix))
                ratiobed <- rep(ratiolabel, ncol(featureMatrix))
                centerbed <- rep(centerLabel, ncol(featureMatrix))
                collabel <- list()
                featuretype <- list()
                for (w in featureNames) {
                    if (verbose) message("Feature name: ", w, "\n")
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
                    dataname <- paste(Ylab, ratiolabel, centerLabel, sep = ":")

                    heatmap_list[[dataname]] <- draw_matrix_heatmap(
                       featureMatrix, dataName = dataname,
                       labels_col = collabel,
                       levels_col = names(scaled_bins[scaled_bins > 0]),
                       ranges = heatRange, verbose = verbose)
                }
                plot_df <- data.frame("Intensity" = colm, "sd" = colsd,
                                      "se" = colse, "Position" = collabel,
                                      "Query" = ratiobed,
                                      "Reference" = centerbed,
                                      "Feature" = featuretype)

                if (smooth) {
                    plot_df$Intensity <- as.vector(
                       smooth.spline(plot_df$Intensity,
                                     df = as.integer(nbins / 5))$y)
                    plot_df$se <- as.vector(
                       smooth.spline(plot_df$se, df = as.integer(nbins / 5))$y)
                }

                mplot_df[[paste(ratiolabel, centerLabel, sep = ":")]] <- plot_df

                regionMatrix <- as.matrix(
                   rmrl[[ratiolabel]][[centerLabel]])

                Intensity <- as.numeric(rowMeans(regionMatrix))

                Query <- as.factor(rep(ratiolabel, length(Intensity)))
                Reference <- as.factor(rep(centerLabel, length(Intensity)))
                subdf <- data.frame(Intensity, Query, Reference)
                stat_df[[paste(ratiolabel, centerLabel, sep = ":")]] <- subdf
            }
        }

        mplot_df <- bind_rows(mplot_df) %>%
            mutate(Group = as.factor(paste0(Query, ":", Reference)),
                   .keep = "all") %>%
            mutate(lower = Intensity - se, upper = Intensity + se)
        mstat_df <- bind_rows(stat_df) %>%
            mutate(Group = as.factor(paste0(Query, ":", Reference)),
                   .keep = "all")

        plot_list <- list()
        for (i in seq_along(ratiolabels)) {
            for (beds in combn(ratiolabels, i, simplify = FALSE)) {
                for (j in seq_along(centerLabels)) {
                    for (centers in combn(centerLabels, j, simplify = FALSE)) {
                        if (verbose) {
                            message("Beds: ", paste(beds, collapse = " "))
                            message("Centers: ", paste(centers, collapse = " "))
                        }

                        aplot_df <- mplot_df %>%
                            filter(Query %in% beds & Reference %in% centers) %>%
                            mutate(Group = paste(Query, Reference, sep = ":"),
                                   .keep = "all")

                        ## plot multi-sample lines with error band
                        p <- draw_region_profile(plot_df = aplot_df,
                                                 cn = "Group", vx = vx,
                                                 Ylab = Ylab)
                        outp <- plot_grid(p, marker, ncol = 1, align = "v",
                                          axis = "lr", rel_heights = c(10, 2))

                        astat_df <- mstat_df %>%
                            filter(Query %in% beds & Reference %in% centers)
                        ## unify order of factors to get consistent color
                        astat_df <- astat_df %>%
                            mutate(Query = factor(
                               Query, levels = sort(unique(Query)))) %>%
                            mutate(Reference = factor(
                               Reference, levels = sort(unique(Reference))))


                        if ((i == 1 && j > 1) || (i > 1 && j == 1)) {
                            if (j > 1) {
                                p <- draw_region_profile(plot_df = aplot_df,
                                                         cn = "Reference",
                                                         sn = "Query", vx = vx,
                                                         Ylab = Ylab)
                                outp <- plot_grid(p, marker, ncol = 1,
                                                  align = "v", axis = "lr",
                                                  rel_heights = c(10, 2))

                                comp <- combn(
                                   seq_along(centers), 2, simplify = FALSE)
                                combo <- draw_combo_plot(stat_df = astat_df,
                                                         xc = "Reference",
                                                         yc = "Intensity",
                                                         comp = comp,
                                                         stats = statsMethod,
                                                         Ylab = Ylab)
                            } else {
                                p <- draw_region_profile(plot_df = aplot_df,
                                                         cn = "Query",
                                                         sn = "Reference",
                                                         vx = vx, Ylab = Ylab)
                                outp <- plot_grid(p, marker, ncol = 1,
                                                  align = "v", axis = "lr",
                                                  rel_heights = c(10, 2))

                                comp <- combn(
                                   seq_along(beds), 2, simplify = FALSE)
                                combo <- draw_combo_plot(stat_df = astat_df,
                                                         xc = "Query",
                                                         yc = "Intensity",
                                                         comp = comp,
                                                         stats = statsMethod,
                                                         Ylab = Ylab)
                            }

                            print(outp)
                            print(combo)
                        } else if (i == 1 && j == 1) {
                            plot_list[[paste(Ylab, beds, centers,
                                             sep = ":")]] <- outp
                        } else if (i == length(queryLabels) &&
                                   j == length(centerLabels)) {
                            print(outp)
                        }
                    }
                }
            }
        }


        draw_stacked_plot(plot_list, heatmap_list)
    }

    if (!is.null(outPrefix)) {
        print(params)
        on.exit(dev.off(), add = TRUE)
    }

    if (verbose) message("[plot_region] finished!\n")
    invisible(mplot_df)
}
