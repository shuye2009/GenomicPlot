#' @title Plot correlation of bam files
#'
#' @description Plot correlation in reads coverage distributions along the genome for bam files. Generates a fingerprint plot, a heatmap of correlation coefficients with hierarchical clustering, a pairwise correlation plot and a PCA plot.
#'
#' @param bamFiles a named vector of strings denoting file names
#' @param binSize an integer denoting the tile width for tiling the genome, default 1000000
#' @param outPrefix a string denoting output file name in pdf format
#' @param importParams a list of parameters for \code{handle_input}
#' @param grouping a named vector for bamFiles group assignment
#' @param verbose logical, indicating whether to output additional information
#' @param hw a vector of two elements specifying the height and width of the output figures
#' @param nc integer, number of cores for parallel processing
#'
#' @return a dataframe of read counts per bin per sample
#'
#' @examples
#'
#' bamQueryFiles <- c(
#'     system.file("extdata", "chip_input_chr19.bam", package = "GenomicPlot"),
#'     system.file("extdata", "chip_treat_chr19.bam", package = "GenomicPlot")
#' )
#' grouping <- c(1, 2)
#' names(bamQueryFiles) <- names(grouping) <- c("chip_input", "chip_treat")
#'
#' bamImportParams <- setImportParams(
#'     offset = 0, fix_width = 150, fix_point = "start", norm = FALSE,
#'     useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' plot_bam_correlation(
#'     bamFiles = bamQueryFiles, binSize = 100000, outPrefix = NULL,
#'     importParams = bamImportParams, nc = 2, verbose = FALSE
#' )
#'
#' @export plot_bam_correlation
#'
plot_bam_correlation <- function(bamFiles,
                                 binSize = 1e6,
                                 outPrefix = NULL,
                                 importParams = NULL,
                                 grouping = NULL,
                                 verbose = FALSE,
                                 hw = c(8, 8),
                                 nc = 2) {
    stopifnot(is.numeric(c(binSize, nc, hw)))
    stopifnot(all(file.exists(bamFiles)))
    if (is.null(names(bamFiles)) || any(names(bamFiles) == ""))
        stop("Each file must have a name attribute!")

    if(verbose) message("[plot_bam_correlation] started ...\n")
    functionName <- as.character(match.call()[[1]])
    params <- plot_named_list(as.list(environment()))
    force(params)

    if (!is.null(outPrefix)) {
        while (!is.null(dev.list())) {
            dev.off()
        }
        pdf(paste0(outPrefix, ".pdf"), width = hw[2], height = hw[1])
    }

    bamlabels <- names(bamFiles)

    if (is.null(importParams)) {
        importParams <- setImportParams(outRle = FALSE)
    } else {
        importParams$outRle <- FALSE # force imported data to be GRanges
    }

    if (verbose) message("Computing bam correlation...\n")
    outlist <- handle_input(inputFiles = bamFiles, importParams, nc = nc)

    seqi <- set_seqinfo(importParams$genome)

    tileBins <- tileGenome(seqi, tilewidth = binSize,
                           cut.last.tile.in.chrom = TRUE)

    grange_list <- lapply(outlist, function(x) x$query)

    score_list <- parallel_countOverlaps(grange_list, tileBins, nc = nc)

    bins_df <- data.frame(chr = as.vector(seqnames(tileBins)),
                          start = start(tileBins), end = end(tileBins),
                          strand = strand(tileBins))
    bins <- do.call(paste, c(bins_df, sep = "_"))

    count_mat <- data.matrix(bind_cols(score_list))
    rownames(count_mat) <- bins
    count_mat[is.na(count_mat)] <- 0
    count_mat <- count_mat[apply(count_mat, 1, sum) > 0, ]

    norm_factor <- vapply(outlist, function(x) x$size / 1e6, numeric(1))

    ## convert to counts per million (CPM)
    df <- as.data.frame(t(t(count_mat) / norm_factor))
    colnames(df) <- bamlabels

    long_df <- pivot_longer(
        df, cols = colnames(df), names_to = "Sample", values_to = "Count") %>%
        mutate(Sample = as.factor(Sample)) %>%
        group_by(Sample) %>%
        arrange(Count) %>%
        mutate(cumCount = cumsum(Count)) %>%
        mutate(Rank = order(cumCount)) %>%
        mutate(Fraction = cumCount / max(cumCount), Rank = Rank / max(Rank))

    p1 <- ggplot(data = long_df, aes(x = Rank, y = Fraction, color = Sample)) +
        geom_line() +
        ggtitle(paste("Binned read counts distribution: bin size =", binSize)) +
        labs(x = "Rank(Count)", y = paste("Fraction over highest coverage"))
    print(p1)

    p2 <- ggplot(data = long_df, aes(x = Count, color = Sample)) +
        stat_ecdf() +
        ggtitle(paste("Binned read counts distribution: bin size =", binSize)) +
        labs(x = expression(paste(log[2], " (Count)")), y = paste("Percentage"))
    # print(p2)


    ## code from pairs example of base R
    ## put (absolute) correlations on the upper panels,
    ## with size proportional to the correlations.
    panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
        # usr <- par("usr"); on.exit(par(usr))
        old <- par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y))
        txt <- format(c(r, 0.123456789), digits = digits)[1]
        txt <- paste0(prefix, txt)
        if (missing(cex.cor)) cex.cor <- 0.8 / strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor * r)
        on.exit(par(old))
    }

    ## put histograms on the diagonal
    panel.hist <- function(x, ...) {
        usr <- par("usr")
        old <- par(usr = c(usr[seq_len(2)], 0, 1.5))
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks
        nB <- length(breaks)
        y <- h$counts
        y <- y / max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
        on.exit(par(old))
    }

    ## END code from pairs example

    mat <- cor(log2(df + 1))
    mat_long <- pivot_longer(as.data.frame(mat), cols = seq_len(ncol(mat)),
                             names_to = "X", values_to = "correlation") %>%
        mutate(Y = rep(rownames(mat), each = ncol(mat)))
    g <- ggplot(mat_long, aes(X, Y)) +
        geom_tile(aes(fill = correlation)) +
        geom_text(aes(label = round(correlation, digits = 2),
                      color = "white")) +
        scale_fill_viridis(discrete = FALSE) +
        theme_minimal() +
        theme(
            legend.position = "none",
            axis.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank()
        )
    # print(g)

    ph <- ComplexHeatmap::pheatmap(mat, col = colorRamp2(
        range(mat), viridis(2)), display_numbers = TRUE,
        heatmap_legend_param = list(title = "correlation"))
    draw(ph)

    if (length(bamFiles) <= 6)
        pairs(log2(df + 1), lower.panel = panel.smooth, upper.panel = panel.cor,
              diag.panel = panel.hist,
              main = paste("log2(CPM/bin), bin size =", binSize))

    # PCA analysis plot
    if (is.null(grouping)) {
        grouping <- seq_along(bamlabels)
        names(grouping) <- bamlabels
    }
    res <- prcomp(t(log2(df + 1)), scale. = TRUE)
    var_explained <- res$sdev^2 / sum(res$sdev^2)
    pca <- as.data.frame(res$x) %>%
        mutate(Samples = rownames(.)) %>%
        mutate(Groups = as.factor(grouping[Samples]))

    pcaplot <- ggplot(pca, aes(x = PC1, y = PC2, color = Groups)) +
        geom_text(label = pca$Samples) +
        theme_bw(base_size = 12) +
        scale_color_npg() +
        labs(
            x = paste0("PC1: ", round(var_explained[1] * 100, 1), "%"),
            y = paste0("PC2: ", round(var_explained[2] * 100, 1), "%")
        ) +
        xlim(c(min(pca$PC1) * 1.2), max(pca$PC1) * 1.2)

    print(pcaplot)

    if (!is.null(outPrefix)) {
        print(params)
        on.exit(dev.off(), add = TRUE)
    }

    if(verbose) message("[plot_bam_correlation] finished!\n")
    invisible(df)
}
