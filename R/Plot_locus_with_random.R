#' @title Plot signal around custom genomic loci and random loci for comparison
#
#' @description Plot reads or peak Coverage/base/gene of samples given in the
#' query files around reference locus defined in the centerFiles. The upstream
#' and downstream windows flanking loci can be given separately, a smaller
#' window can be defined to allow statistical comparisons between reference and
#' random loci. The loci are further divided into sub-groups that are
#' overlapping with c("5'UTR", "CDS", "3'UTR"), "unrestricted" means all loci
#' regardless of overlapping.
#'
#' @param queryFiles a vector of sample file names. The file should be in .bam,
#'  .bed, .wig or .bw format, mixture of formats is allowed
#' @param centerFiles a vector of reference file names. The file should be .bed
#'  format only
#' @param inputFiles a vector of input sample file names. The file should be in
#'  .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param txdb a TxDb object defined in the `GenomicFeatures` package
#' @param importParams a list of parameters for \code{\link{handle_input}}
#' @param ext a vector of two integers defining upstream and downstream
#'  boundaries of the plot window, flanking the reference locus
#' @param hl a vector of two integers defining upstream and downstream
#'  boundaries of the highlight window, flanking the reference locus
#' @param stranded logical, indicating whether the strand of the feature should
#'  be considered
#' @param scale logical, indicating whether the score matrix should be scaled to
#'  the range 0:1, so that samples with different baseline can be compared
#' @param smooth logical, indicating whether the line should smoothed with a
#'  spline smoothing algorithm
#' @param rmOutlier a numeric value serving as a multiplier of the MAD in Hampel
#'  filter for outliers identification, 0 indicating not removing outliers.
#'  For Gaussian distribution, use 3, adjust based on data distribution
#' @param outPrefix a string specifying output file prefix for plots
#' (outPrefix.pdf)
#' @param refPoint a string in c("start", "center", "end")
#' @param Xlab a string denotes the label on x-axis
#' @param Ylab a string for y-axis label
#' @param shade logical indicating whether to place a shaded rectangle around
#'  the point of interest
#' @param binSize an integer defines bin size for intensity calculation
#' @param transform a string in c("log", "log2", "log10"), default = NA i
#'  ndicating no transformation of data matrix
#' @param n_random an integer denotes the number of randomization should be
#'  performed
#' @param statsMethod a string in c("wilcox.test", "t.test"), for pair-wise
#'  groups comparisons
#' @param verbose logical, indicating whether to output additional information
#'  (data used for plotting or statistical test results)
#' @param hw a vector of two elements specifying the height and width of the
#'  output figures
#' @param nc integer, number of cores for parallel processing
#' @param detailed logical, indicating whether to plot each parts of gene.
#'
#' @return a dataframe containing the data used for plotting
#' @author Shuye Pu
#'
#' @examples
#'
#' gtfFile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#'     package = "GenomicPlot"
#' )
#'
#' txdb <- custom_TxDb_from_GTF(gtfFile, genome = "hg19")
#' bedQueryFiles <- c(
#' system.file("extdata", "test_chip_peak_chr19.narrowPeak",
#'             package = "GenomicPlot"),
#' system.file("extdata", "test_chip_peak_chr19.bed", package = "GenomicPlot"),
#' system.file("extdata", "test_clip_peak_chr19.bed", package = "GenomicPlot")
#' )
#' names(bedQueryFiles) <- c("NarrowPeak", "SummitPeak", "iCLIPPeak")
#'
#' bamQueryFiles <- system.file("extdata", "treat_chr19.bam",
#'                              package = "GenomicPlot")
#' names(bamQueryFiles) <- "clip_bam"
#' bamInputFiles <- system.file("extdata", "input_chr19.bam",
#'                              package = "GenomicPlot")
#' names(bamInputFiles) <- "clip_input"
#'
#' bamImportParams <- setImportParams(
#'   offset = -1, fix_width = 0, fix_point = "start", norm = TRUE,
#'   useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#' )
#' plot_locus_with_random(
#'     queryFiles = bamQueryFiles,
#'     centerFiles = bedQueryFiles[3],
#'     txdb = txdb,
#'     ext = c(-200, 200),
#'     hl = c(-50, 50),
#'     shade = TRUE,
#'     importParams = bamImportParams,
#'     verbose = FALSE,
#'     smooth = TRUE,
#'     transform = NA,
#'     binSize = 10,
#'     refPoint = "center",
#'     Xlab = "Center",
#'     Ylab = "Coverage/base/peak",
#'     inputFiles = bamInputFiles,
#'     stranded = TRUE,
#'     scale = FALSE,
#'     outPrefix = NULL,
#'     rmOutlier = 0,
#'     n_random = 1,
#'     hw = c(8, 8),
#'     detailed = FALSE,
#'     statsMethod = "wilcox.test",
#'     nc = 2)
#'
#' @export plot_locus_with_random

plot_locus_with_random <- function(queryFiles,
                                   centerFiles,
                                   txdb,
                                   ext = c(-200, 200),
                                   hl = c(-100, 100),
                                   shade = FALSE,
                                   importParams = NULL,
                                   verbose = FALSE,
                                   smooth = FALSE,
                                   transform = NA,
                                   binSize = 10,
                                   refPoint = "center",
                                   Xlab = "Center",
                                   Ylab = "Coverage/base/gene",
                                   inputFiles = NULL,
                                   stranded = TRUE,
                                   scale = FALSE,
                                   outPrefix = NULL,
                                   rmOutlier = 0,
                                   n_random = 1,
                                   hw = c(8, 8),
                                   detailed = FALSE,
                                   statsMethod = "wilcox.test",
                                   nc = 2) {
    stopifnot(is.numeric(c(ext, hl, binSize, nc, hw, rmOutlier)))
    stopifnot(transform %in% c("log", "log2", "log10", NA))
    stopifnot(all(file.exists(queryFiles)))
    if (is.null(names(queryFiles)) || any(names(queryFiles) == ""))
        stop("Each file must have a name attribute!")

    if (verbose) message("[plot_locus_with_random] started ...\n")
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
            queryInputs <- handle_input(
                inputFiles = c(queryFiles, inputFiles), importParams,
                verbose = verbose, nc = nc)
        } else if (length(inputFiles) == 1) {
            queryInputs <- handle_input(
                inputFiles = c(queryFiles, inputFiles), importParams,
                verbose = verbose, nc = nc)
            queryInputs <- queryInputs[c(queryLabels,
                                         rep(inputLabels, length(queryLabels)))]

            inputLabels <- paste0(names(inputFiles), seq_along(queryFiles))
            names(queryInputs) <- c(queryLabels, inputLabels)
        } else {
            stop("the number of inputFiles must be 1 or equal to the number of
                 queryFiles!")
        }
    }
    queryLabels <- names(queryInputs)

    ## to avoid binSize inconsistency, as the final binSize depends bin_num
    ext[2] <- ext[2] - (ext[2] - ext[1]) %% binSize
    colLabel <- seq(ext[1], (ext[2] - binSize), binSize)
    bin_num <- round((ext[2] - ext[1]) / binSize)
    bin_op <- "mean"

    ## ranges for overlap count
    rb <- ext[2] # broad
    rn <- hl[2] # narrow

    if (detailed) {

        ## get protein-coding genes features

        if (verbose) message("Collecting protein_coding gene features...\n")
        exons <- get_genomic_feature_coordinates(txdb, "exon", longest = TRUE)
        utr5 <- get_genomic_feature_coordinates(txdb, "utr5", longest = TRUE)
        utr3 <- get_genomic_feature_coordinates(txdb, "utr3", longest = TRUE)
        cds <- get_genomic_feature_coordinates(txdb, "cds", longest = TRUE)
        gene <- get_genomic_feature_coordinates(txdb, "gene", longest = TRUE)
        region_list <- list(
            "Transcript" = exons$GRanges,
            "5'UTR" = utr5$GRanges,
            "CDS" = cds$GRanges,
            "3'UTR" = utr3$GRanges,
            "Gene" = gene$GRanges,
            "unrestricted" = NULL
        )
    } else {
        region_list <- list("unrestricted" = NULL)
    }

    if (verbose) {
        message("Region length: ",
                paste(lapply(region_list, length), collapse = " "), "\n")
        message("Computing coverage for Sample...\n")
    }
    sml <- list() # scoreMatrix_list
    smlr <- list() # scoreMatrix_list_random

    bedparam <- importParams
    bedparam$offset <- 0
    bedparam$fix_width <- 0
    bedparam$useScore <- FALSE
    bedparam$outRle <- FALSE
    bedparam$useSizeFactor <- FALSE

    centerInputs <- handle_input(
        centerFiles, bedparam, verbose = verbose, nc = nc)
    centerLabels <- names(centerInputs)

    for (queryLabel in queryLabels) {
        myInput <- queryInputs[[queryLabel]]
        libSize <- myInput$size
        queryRegions <- myInput$query
        fileType <- myInput$type
        weight_col <- myInput$weight

        if (verbose) {
            message("Query label: ", queryLabel)
            message("Library size: ", libSize, "\n")
        }

        for (centerLabel in centerLabels) {

            centerInput <- centerInputs[[centerLabel]]
            centerSize <- centerInput$size
            windowRegionsALL <- centerInput$query
            if (verbose) {
                message("Center label: ", centerLabel)
                message("Center size: ", centerSize, "\n")
            }

            for (regionName in names(region_list)) {
                if (verbose) message("Processing genomic region: ", regionName)

                if (regionName == "unrestricted") {
                    windowRegions <- windowRegionsALL
                } else {
                    region <- region_list[[regionName]]
                    if (any(unique(as.vector(
                        strand(windowRegionsALL))) %in% c("*", ".", ""))) {
                        windowRegions <- plyranges::filter_by_overlaps(
                            windowRegionsALL, region)
                        if (verbose) message("The center file is Unstranded\n")
                    } else {
                        windowRegions <- filter_by_overlaps_stranded(
                            windowRegionsALL, region, ignore.order = FALSE)
                        if (verbose) message("The center file is stranded\n")
                    }
                }

                windowRs <- resize(windowRegions, width = 1, fix = refPoint)
                windowRs <- promoters(windowRs, upstream = -ext[1],
                                      downstream = ext[2])
                windowRs <- as(split(windowRs, f = factor(seq_along(windowRs))),
                               "GRangesList")

                fullMatrix <- parallel_scoreMatrixBin(
                    queryRegions, windowRs, bin_num, bin_op, weight_col,
                    stranded, nc = nc)
                if (is.null(inputFiles)) {
                    fullMatrix <- process_scoreMatrix(
                        fullMatrix, scale, rmOutlier, transform = transform,
                        verbose = verbose)
                } else {
                    fullMatrix <- process_scoreMatrix(
                        fullMatrix, scale = FALSE, rmOutlier = rmOutlier,
                        transform = NA, verbose = verbose)
                }

                sml[[queryLabel]][[centerLabel]][[regionName]] <-
                    fullMatrix

                ## create randomized centers, repeat n_random times
                windowRs <- list()
                for (i in seq(1, n_random)) {
                    if (verbose) message("Creating random centers ... ", i)
                    random_points <- sample(
                        ext[1]:ext[2], length(windowRegions), replace = TRUE)
                    rwindowRegions <- shift(
                        windowRegions, shift = random_points, use.names = TRUE)

                    windowR <- resize(rwindowRegions, width = 1, fix = refPoint)
                    windowR <- promoters(
                        windowR, upstream = -ext[1], downstream = ext[2])
                    windowR <- as(split(windowR,
                                        f = factor(seq_along(windowR))),
                                  "GRangesList")
                    windowRs[[i]] <- windowR
                }

                random_matricies <- list()
                for (i in seq_along(windowRs)) {
                    if (verbose) message("Computing for random centers ... ", i)
                    mat <- parallel_scoreMatrixBin(queryRegions, windowRs[[i]],
                                                   bin_num, bin_op, weight_col,
                                                   stranded, nc = nc)
                    random_matricies[[i]] <- mat
                }

                ## unify the dimensions of random_matrices
                if (verbose) message("Averaging random matricies ...\n")
                dims <- vapply(random_matricies, dim, numeric(2))
                minrows <- min(dims[1, ])
                random_matricies <- lapply(random_matricies, function(x)
                    x[seq_len(minrows), ])

                ## turn the list of matrices into a 3-d array
                a3d <- array(unlist(random_matricies),
                             c(dim(random_matricies[[1]]),
                               length(random_matricies)))
                ## take average of all matrices
                fullMatrix <- apply(a3d, seq_len(2), mean)

                fullMatrix <- process_scoreMatrix(
                    fullMatrix, scale = scale, rmOutlier = rmOutlier,
                    transform = transform, verbose = verbose)

                smlr[[queryLabel]][[centerLabel]][[regionName]] <- fullMatrix
            }
        }
    }


    ## plot reference center and random center for each bed

    Ylab <- ifelse(!is.na(transform) && is.null(inputFiles),
                   paste0(transform, " (", Ylab, ")"), Ylab)
    for (queryLabel in queryLabels) {
        if (verbose) message("Processing query: ", queryLabel)

        for (centerLabel in centerLabels) {
            if (verbose) message("Processing reference: ", centerLabel)

            for (regionName in names(region_list)) {
                if (verbose) message("Processing genomic region: ", regionName)
                plot_df <- NULL
                stat_df <- NULL
                countOverlap_df <- NULL

                fullMatrix_list <- list(
                    sml[[queryLabel]][[centerLabel]][[regionName]],
                    smlr[[queryLabel]][[centerLabel]][[regionName]]
                )
                names(fullMatrix_list) <- c(centerLabel, "Random")

                refsize <- nrow(fullMatrix_list[[centerLabel]])

                for (alabel in c(centerLabel, "Random")) {
                    fullMatrix <- fullMatrix_list[[alabel]]
                    colm <- apply(fullMatrix, 2, mean)
                    if(nrow(fullMatrix) == 1){
                        colsd <- rep(0, ncol(fullMatrix))
                    }else{
                        colsd <- apply(fullMatrix, 2, sd)
                    }
                    colse <- colsd / sqrt(apply(fullMatrix, 2, length))
                    collabel <- colLabel
                    querybed <- as.factor(rep(queryLabel, length(colm)))
                    refbed <- as.factor(rep(alabel, length(colm)))

                    sub_df <- data.frame(colm, colsd, colse, collabel,
                                         querybed, refbed)
                    colnames(sub_df) <- c("Intensity", "sd", "se", "Position",
                                          "Query", "Reference")

                    if (smooth) {
                        sub_df$Intensity <- as.vector(
                            smooth.spline(sub_df$Intensity,
                                          df = as.integer(bin_num / 5))$y)
                        sub_df$se <- as.vector(
                            smooth.spline(sub_df$se,
                                          df = as.integer(bin_num / 5))$y)
                    }

                    sub_df <- mutate(sub_df, lower = Intensity - se,
                                     upper = Intensity + se)
                    plot_df <- rbind(plot_df, sub_df)

                    if (hl[2] > hl[1]) {
                        xmin <- which(colLabel == hl[1])
                        xmax <- which(colLabel == hl[2])
                        if (length(xmax) == 0) xmax <- length(colLabel)
                        submatrix <- fullMatrix[, xmin:xmax, drop = FALSE]
                        submatrix[is.na(submatrix)] <- 0
                        Intensity <- as.numeric(rowMeans(submatrix))

                        Query <- as.factor(rep(queryLabel, length(Intensity)))
                        Reference <- as.factor(rep(alabel, length(Intensity)))
                        subdf <- data.frame(Intensity, Query, Reference)

                        stat_df <- rbind(stat_df, subdf)
                    }
                }

                plot_df <- plot_df %>%
                    mutate(Query = factor(
                        Query, levels = sort(unique(Query)))) %>%
                    mutate(Reference = factor(
                        Reference, levels = sort(unique(Reference))))

                p <- draw_locus_profile(
                    plot_df, cn = "Reference", sn = "Query", Xlab = Xlab,
                    Ylab = Ylab, shade = shade, hl = hl) +
                    ggtitle(paste("Feature:", regionName, "\nReference size:",
                                  refsize, "\nSample name:", queryLabel))

                if (hl[2] > hl[1]) {
                    stat_df <- stat_df %>%
                        mutate(Query = factor(
                            Query, levels = sort(unique(Query)))) %>%
                        mutate(Reference = factor(
                            Reference, levels = sort(unique(Reference))))
                    comp <- list(c(1, 2))

                    combo <- draw_combo_plot(
                        stat_df = stat_df, xc = "Reference", yc = "Intensity",
                        comp = comp, stats = statsMethod, Ylab = Ylab)
                    print(p)
                    print(combo)
                } else {
                    print(p)
                }
            }
        }
    }

    if (!is.null(inputFiles)) {
        Ylab <- ifelse(is.na(transform), "Ratio-over-Input",
                       paste0(transform, " (Ratio-over-Input)"))

        rLs <- queryLabels[!queryLabels %in% inputLabels] # ratioLabels
        iml <- sml[inputLabels] # inputMatrix_list
        rml <- sml[rLs] # ratioMatrix_list

        imlr <- smlr[inputLabels] # inputMatrix_list_random
        rmlr <- smlr[rLs] # ratioMatrix_list_random

        for (centerLabel in centerLabels) {
            for (regionName in names(region_list)) {
                for (i in seq_along(rLs)) {
                    rm <- rml[[rLs[i]]][[centerLabel]][[regionName]]
                    im <- iml[[inputLabels[i]]][[centerLabel]][[regionName]]

                    fullMatrix <- ratio_over_input(rm, im, verbose)

                    fullMatrix <- process_scoreMatrix(
                        fullMatrix, scale, rmOutlier, transform, verbose)

                    rml[[rLs[i]]][[centerLabel]][[regionName]] <- fullMatrix

                    ## for random centers
                    rmr <- rmlr[[rLs[i]]][[centerLabel]][[regionName]]
                    imr <- imlr[[inputLabels[i]]][[centerLabel]][[regionName]]
                    minrowr <- min(nrow(rmr), nrow(imr))

                    fullMatrix <- ratio_over_input(
                        rmr[seq_len(minrowr), ], imr[seq_len(minrowr), ])

                    fullMatrix <- process_scoreMatrix(
                        fullMatrix, scale, rmOutlier, transform, verbose)

                    rmlr[[rLs[i]]][[centerLabel]][[regionName]] <- fullMatrix
                }
            }
        }

        for (ratiolabel in rLs) {
            if (verbose) message("Processing ratio for query: ", ratiolabel)
            for (centerLabel in centerLabels) {
                if (verbose) message("Processing ratio for reference: ",
                                     centerLabel)
                for (regionName in names(region_list)) {
                    if (verbose) message(
                        "Processing ratio for genomic region: ", regionName)
                    plot_df <- NULL
                    stat_df <- NULL
                    countOverlap_df <- NULL
                    refsize <- NULL

                    fullMatrix_list <- list(
                        rml[[ratiolabel]][[centerLabel]][[regionName]],
                        rmlr[[ratiolabel]][[centerLabel]][[regionName]]
                    )
                    names(fullMatrix_list) <- c(centerLabel, "Random")

                    refsize <- nrow(fullMatrix_list[[centerLabel]])

                    for (alabel in c(centerLabel, "Random")) {

                        fullMatrix <- fullMatrix_list[[alabel]]

                        colm <- apply(fullMatrix, 2, mean)
                        if(nrow(fullMatrix) == 1){
                            colsd <- rep(0, ncol(fullMatrix))
                        }else{
                            colsd <- apply(fullMatrix, 2, sd)
                        }
                        colse <- colsd / sqrt(apply(fullMatrix, 2, length))
                        collabel <- colLabel
                        querybed <- as.factor(rep(queryLabel, length(colm)))
                        refbed <- as.factor(rep(alabel, length(colm)))

                        sub_df <- data.frame(colm, colsd, colse, collabel,
                                             querybed, refbed)
                        colnames(sub_df) <- c("Intensity", "sd", "se",
                                              "Position", "Query", "Reference")

                        if (smooth) {
                            sub_df$Intensity <- as.vector(
                                smooth.spline(sub_df$Intensity,
                                              df = as.integer(bin_num / 5))$y)
                            sub_df$se <- as.vector(
                                smooth.spline(sub_df$se,
                                              df = as.integer(bin_num / 5))$y)
                        }

                        sub_df <- mutate(sub_df, lower = Intensity - se,
                                         upper = Intensity + se)
                        plot_df <- rbind(plot_df, sub_df)

                        if (hl[2] > hl[1]) {
                            xmin <- which(colLabel == hl[1])
                            xmax <- which(colLabel == hl[2])
                            if (length(xmax) == 0) xmax <- length(colLabel)
                            submatrix <- fullMatrix[, xmin:xmax, drop = FALSE]
                            submatrix[is.na(submatrix)] <- 0
                            Intensity <- as.numeric(rowMeans(submatrix))
                            Query <- as.factor(rep(queryLabel,
                                                   length(Intensity)))
                            Reference <- as.factor(rep(alabel,
                                                       length(Intensity)))
                            subdf <- data.frame(Intensity, Query, Reference)

                            stat_df <- rbind(stat_df, subdf)
                        }
                    }

                    plot_df <- plot_df %>%
                        mutate(Query = factor(
                            Query, levels = sort(unique(Query)))) %>%
                        mutate(Reference = factor(
                            Reference, levels = sort(unique(Reference))))

                    p <- draw_locus_profile(
                        plot_df, cn = "Reference", sn = "Query", Xlab = Xlab,
                        Ylab = Ylab, shade = shade, hl = hl) +
                        ggtitle(paste("Feature:", regionName,
                                      "\nReference size:", refsize,
                                      "\nSample name:", ratiolabel))

                    if (hl[2] > hl[1]) {
                        stat_df <- stat_df %>%
                            mutate(Query = factor(
                                Query, levels = sort(unique(Query)))) %>%
                            mutate(Reference = factor(
                                Reference, levels = sort(unique(Reference))))
                        comp <- list(c(1, 2))

                        combo <- draw_combo_plot(
                            stat_df = stat_df, xc = "Reference",
                            yc = "Intensity", comp = comp,
                            stats = statsMethod, Ylab = Ylab)

                        print(p)
                        print(combo)
                    } else {
                        print(p)
                    }
                }
            }
        }
    }

    if (!is.null(outPrefix)) {
        print(params)
        on.exit(dev.off(), add = TRUE)
    }
    if (verbose) message("[plot_locus_with_random] finished!\n")
}
