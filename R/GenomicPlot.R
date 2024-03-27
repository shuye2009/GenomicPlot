
#' @title GenomicPlot-package
#' @description An R package for efficient and flexible visualization of
#' genome-wide NGS coverage profiles
#'
#' @details The goal of `GenomicPlot` is to provide an efficient visualization
#' tool for next generation sequencing (NGS) data with rich functionality and
#' flexibility. `GenomicPlot` enables plotting of NGS data in various formats
#' (bam, bed, wig and bigwig); both coverage and enrichment over input can be
#' computed and displayed with respect to genomic features (such as UTR, CDS,
#' enhancer), and user defined genomic loci or regions. Statistical tests on
#' signal intensity within user defined regions of interest can be performed
#' and presented as box plots or pie charts. Parallel processing is enabled to
#' speed up computation on multi-core platforms.
#' Main functions are as follows:
#' \itemize{
#'  \item   \code{\link{plot_5parts_metagene}} generates genomic (with introns)
#'      or metagenomic (without introns) plots around gene body and its upstream
#'      and downstream regions, the gene body can be further segmented into
#'      5'UTR, CDS and 3'UTR.
#'  \item   \code{\link{plot_start_end}} plots genomic profiles around the start
#'      and end of genomic features (like exons or introns), or user defined
#'      genomic regions. A center region with user defined width can be plotted
#'      simultaneously.
#'  \item   \code{\link{plot_locus}} plots distance between sample peaks and
#'      genomic features, or distance from one set of peaks to another set of
#'      peaks.
#'  \item   \code{\link{plot_region}} plots signal profiles within and around
#'      genomic features, or user defined genomic regions.
#'  \item   \code{\link{plot_peak_annotation}} plots peak annotation statistics
#'      (distribution in different type of genes, and in different parts of
#'      genes).
#'  \item   \code{\link{plot_overlap_bed}} plots peak overlaps as Venn diagrams.
#'  \item   Random features can be generated and plotted to serve as contrast to
#'      real features in \code{plot_locus_with_random} and
#'      \code{\link{plot_start_end_with_random}}.
#'  \item   All profile line plots have error bands.
#'  \item   Statistical analysis results on user defined regions of interest are
#'      plotted along with the profile plots in \code{\link{plot_region}},
#'      \code{\link{plot_locus}} and \code{\link{plot_locus_with_random}}.
#' }
#'
#' @author Shuye Pu
#'
#' _PACKAGE
#' @name GenomicPlot
NULL


#' @title Plot Venn diagrams depicting overlap of genomic regions
#' @description This function takes a list of up to 4 bed file names, and
#' produce a Venn diagram
#' @param bedList a named list of bed files, with list length = 2, 3 or 4
#' @param outPrefix a string for plot file name
#' @param importParams a list of parameters for \code{handle_input}
#' @param stranded logical, indicating whether the feature is stranded. For
#'  nonstranded feature, only "*" is accepted as strand
#' @param pairOnly logical, indicating whether only pair-wise overlap is
#'  desirable
#' @param verbose logical, indicating whether to output additional information
#' @param hw a vector of two elements specifying the height and width of the
#'  output figures
#'
#' @return a ggplot object
#' @author Shuye Pu
#'
#' @examples
#'
#' queryFiles <- c(
#'     system.file("extdata", "test_chip_peak_chr19.narrowPeak",
#'         package = "GenomicPlot"
#'     ),
#'     system.file("extdata", "test_chip_peak_chr19.bed",
#'         package = "GenomicPlot"
#'     ),
#'     system.file("extdata", "test_clip_peak_chr19.bed",
#'         package = "GenomicPlot"
#'     )
#' )
#' names(queryFiles) <- c("narrowPeak", "summitPeak", "clipPeak")
#'
#' bedimportParams <- setImportParams(
#'     offset = 0, fix_width = 100, fix_point = "center", norm = FALSE,
#'     useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' plot_overlap_bed(
#'     bedList = queryFiles, importParams = bedimportParams, pairOnly = FALSE,
#'     stranded = FALSE, outPrefix = NULL
#' )
#'
#' @export plot_overlap_bed


plot_overlap_bed <- function(bedList,
                             outPrefix = NULL,
                             importParams = NULL,
                             pairOnly = TRUE,
                             stranded = TRUE,
                             hw = c(8, 8),
                             verbose = FALSE) {
    stopifnot(all(file.exists(bedList)))
    stopifnot(!is.null(names(bedList)))
    if (!is.null(outPrefix)) {
        pdf(paste0(outPrefix, ".pdf"), width = hw[2], height = hw[1])
    }

    functionName <- as.character(match.call()[[1]])
    params <- plot_named_list(as.list(environment()))
    force(params)

    if (is.null(importParams)) {
        importParams <- setImportParams(outRle = FALSE)
    } else {
        importParams$outRle <- FALSE # force imported data to be GRanges
    }

    inputList <- handle_input(bedList, importParams)
    names(inputList) <- names(bedList)
    grList <- lapply(inputList, function(x) x$query)
    sizeList <- lapply(inputList, function(x) x$size)
    if (verbose) message("Sizes: ", paste(sizeList, collapse = " "), "\n")

    # get all pair-wise overlap counts into a matrix, display as a heatmap

    fil <- ifelse(stranded,
                  filter_by_overlaps_stranded,
                  filter_by_overlaps_nonstranded)

    counts <- vapply(grList, function(i)
        vapply(grList, function(j)
            length(fil(i, j, ignore.order = FALSE)), integer(1)),
        integer(length(grList))) |> t()

    rownames(counts) <- colnames(counts) <- names(bedList)
    counts_long <- tidyr::pivot_longer(as.data.frame(counts),
                                cols = seq_len(ncol(counts)),
                                names_to = "X", values_to = "count") %>%
        mutate(Y = rep(rownames(counts), each = ncol(counts)))

    pairs <- combn(grList, 2, simplify = FALSE)

    g <- ggplot(counts_long, aes(X, Y)) +
        geom_tile(aes(fill = count)) +
        geom_text(aes(label = count, color = "white", size = 10)) +
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
    print(g)
    grid.newpage()

    lapply(pairs, overlap_pair, fil)
    if (!pairOnly) {
        if (length(grList) > 2) {
            triples <- combn(grList, 3, simplify = FALSE)
            lapply(triples, overlap_triple, fil)
            if (length(grList) > 3) {
                quads <- combn(grList, 4, simplify = FALSE)
                lapply(quads, overlap_quad, fil)
            }
        }
    }

    if (!is.null(outPrefix)) {
        print(params)
        on.exit(dev.off(), add = TRUE)
    }
}

#' @title Plot Venn diagrams depicting overlap of gene lists
#' @description This function takes a list of (at most 4) tab-delimited file
#' names, and produce a Venn diagram
#' @param fileList, a named list of tab-delimited files
#' @param columnList a vector of integers denoting the columns that have gene
#'  names in the list of files
#' @param outPrefix, a string for plot file name
#' @param pairOnly, logical, indicating whether only pair-wise overlap is
#'  desirable
#' @param hw a vector of two elements specifying the height and width of the
#'  output figures
#'
#' @return a list of vectors of gene names
#' @author Shuye Pu
#'
#' @examples
#' testfile1 <- system.file("extdata", "test_file1.txt",
#'     package = "GenomicPlot"
#' )
#' testfile2 <- system.file("extdata", "test_file2.txt",
#'     package = "GenomicPlot"
#' )
#' testfile3 <- system.file("extdata", "test_file3.txt",
#'     package = "GenomicPlot"
#' )
#' testfile4 <- system.file("extdata", "test_file4.txt",
#'     package = "GenomicPlot"
#' )
#' testfiles <- c(testfile1, testfile2, testfile3, testfile4)
#' names(testfiles) <- c("test1", "test2", "test3", "test4")
#'
#' plot_overlap_genes(testfiles, c(3, 2, 1, 1), pairOnly = FALSE)
#'
#' @export plot_overlap_genes

plot_overlap_genes <- function(fileList,
                               columnList,
                               pairOnly = TRUE,
                               hw = c(8, 8),
                               outPrefix = NULL) {
    stopifnot(all(file.exists(fileList)))
    stopifnot(is.numeric(columnList))

    functionName <- as.character(match.call()[[1]])
    params <- plot_named_list(as.list(environment()))
    force(params)

    geneList <- mapply(x = fileList, y = columnList, function(x, y) {
        df <- read.delim(x, header = TRUE, sep = "\t")
        genes <- unique(df[, y])
        genes
    })


    if (!is.null(outPrefix)) {
        pdf(paste0(outPrefix, ".pdf"), width = hw[2], height = hw[1])
    }
    pairs <- combn(geneList, 2, simplify = FALSE)

    lapply(pairs, overlap_pair, intersect)
    if (!pairOnly) {
        if (length(geneList) > 2) {
            triples <- combn(geneList, 3, simplify = FALSE)
            lapply(triples, overlap_triple, intersect)
            if (length(geneList) > 3) {
                quads <- combn(geneList, 4, simplify = FALSE)
                lapply(quads, overlap_quad, intersect)
            }
        }
    }

    if (!is.null(outPrefix)) {
        print(params)
        on.exit(dev.off(), add = TRUE)
    }

    invisible(geneList)
}

#' @title plot a named list as a figure
#' @description This is a helper function for displaying function arguments for
#' a plotting function. If the runtime value of the argument is a small object,
#' its values is displayed, otherwise, only the name of the value of the
#' argument is displayed.
#'
#' @param params a list produced by as.list(environment()), with names being the
#'  arguments and values being the runtime values when the function is called.
#'
#' @return a ggplot object
#'
#' @author Shuye Pu
#'
#' @examples
#'
#' data(gf5_genomic)
#'
#' gtfFile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#'     package = "GenomicPlot"
#' )
#'
#' txdb <- custom_TxDb_from_GTF(gtfFile, genome = "hg19")
#'
#' queryfiles <- system.file("extdata", "treat_chr19.bam",
#'     package = "GenomicPlot"
#' )
#' names(queryfiles) <- "query"
#'
#' inputfiles <- system.file("extdata", "input_chr19.bam",
#'     package = "GenomicPlot"
#' )
#' names(inputfiles) <- "input"
#'
#' bamimportParams <- setImportParams(
#'     offset = -1, fix_width = 0, fix_point = "start", norm = TRUE,
#'     useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' alist <- list(
#'     "txdb" = txdb, "treat" = queryfiles, "control" = inputfiles,
#'     "feature" = gf5_genomic, "param" = bamimportParams
#' )
#'
#' GenomicPlot:::plot_named_list(alist)
#'
#' @keywords  internal

plot_named_list <- function(params) {
    s <- "Plotting parameters:\n"
    for (aname in names(params)) {
        osize <- object.size(params[[aname]])
        if (osize > 2048) { # if (value_length > 10) {
            value <- deparse1(substitute(params[[aname]]))
        } else {
            value <- deparse1(params[[aname]])
        }
        s <- paste(s, paste0(aname, ": ", paste(strwrap(value, width = 80),
            collapse = "\n"
        )), sep = "\n")
    }
    p <- ggplot() +
        annotate("text",
            x = 0.5,
            y = 0.5,
            label = s
        ) +
        theme_void()

    return(p)
}
