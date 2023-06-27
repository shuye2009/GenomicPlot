
#' @title GenomicPlot-package
#' @description An R package for efficient and flexible visualization of NGS 
#' coverage profiles
#' 
#' @details The goal of `GenomicPlot` is to provide an efficient visualization
#' tool for next generation sequencing (NGS) data with rich functionality and
#' flexibility. `GenomicPlot` enables plotting of NGS data in various formats 
#' (bam, bed, wig and bigwig); both coverage and enrichment over input can be
#' computed and displayed with respect to genomic features (such as UTR, CDS, 
#' enhancer), and user defined genomic loci or regions. Statistical tests on 
#' signal intensity within user defined regions of interest can be performed 
#' and presented as box plots or pie charts. Parallel processing is enabled to 
#' speed up computation on multi-core platforms. Main functions are as follows:
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
#' @docType package
#' @name GenomicPlot
NULL

#' 
#' @title Toy data for examples and testing
#' @description The data files in the extdata directory contain data for next 
#' generation sequencing read alignments, MACS2 peaks and gene annotation, which 
#' are used to test the package and generate plots in the package vignettes. 
#' Except for the gtf file, all other files are derived from experimental data 
#' produced in-house. To meet the package file size limit, all data are 
#' restricted to chr19:58000-507000 of the human genome version hg19. Details 
#' for each file are as follows.
#' @details
#' \itemize{
#'  \item "gencode.v19.annotation_chr19.gtf" is an excerpt of the gene 
#'      annotation file downloaded from 
#'      \url{https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz}. 
#'      "gencode.v19.annotation_chr19.gtf.granges.rds" is a GRanges object as a 
#'      result of importing the gtf file.
#'  \item "chip_treat_chr19.bam(.bai)" and "chip_input_chr19.bam(.bai)" are 
#'      paired-end read alignment data from ChIPseq experiments.
#'  \item "treat_chr19.bam(.bai)" and "input_chr19.bam(.bai)" are single-end 
#'      read alignment data from iCLIP experiments.
#'  \item "test_wig_chr19_+(-).wig", "test_wig_chr19_+(-).bw" are iCLIP 
#'      alignment data in WIG and BIGWIG format, respectively; '+' and '-' 
#'      represent forward and reverse strand, respectively. 
#'  \item "test_clip_peak_chr19.bed" contains strand-specific iCLIP peak in BED 
#'      format. 
#'  \item "test_chip_peak_chr19.bed" and "test_chip_peak_chr19.narrowPeak" 
#'      contain ChIPseq peaks generated with MACS2, in summit peak and narrow 
#'      peak format, respectively. 
#'  \item "test_file1.txt", "test_file2.txt", "test_file3.txt" and 
#'      "test_file4.txt" are tab-delimited text files,  each contains various 
#'      human gene names in different columns. 
#' }
#' @returns These files will be read by handle_input function, which will 
#'      produce GRanges object or RleList object, depending on import parameters
#' @author Shuye Pu
#' @docType data
#' @name extdata
#' @keywords datasets
#' 
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

    counts <- matrix(rep(0L, length(grList)^2), nrow = length(grList))
    for (i in seq_along(grList)) {
        for (j in seq_along(grList)) {
            if (stranded) {
                counts[i, j] <- length(filter_by_overlaps_stranded(grList[[i]], 
                                                                   grList[[j]]))
            } else {
                counts[i, j] <- length(filter_by_overlaps(grList[[i]], 
                                                          grList[[j]]))
            }
        }
    }
    rownames(counts) <- colnames(counts) <- names(bedList)
    counts_long <- pivot_longer(as.data.frame(counts), 
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
    if (stranded) {
        lapply(pairs, overlap_pair, filter_by_overlaps_stranded)
        if (!pairOnly) {
            if (length(grList) > 2) {
                triples <- combn(grList, 3, simplify = FALSE)
                lapply(triples, overlap_triple, filter_by_overlaps_stranded)
                if (length(grList) > 3) {
                    quads <- combn(grList, 4, simplify = FALSE)
                    lapply(quads, overlap_quad, filter_by_overlaps_stranded)
                }
            }
        }
    } else {
        lapply(pairs, overlap_pair, filter_by_overlaps)
        if (!pairOnly) {
            if (length(grList) > 2) {
                triples <- combn(grList, 3, simplify = FALSE)
                lapply(triples, overlap_triple, filter_by_overlaps)
                if (length(grList) > 3) {
                    quads <- combn(grList, 4, simplify = FALSE)
                    lapply(quads, overlap_quad, filter_by_overlaps)
                }
            }
        }
    }
    if (!is.null(outPrefix)) {
        on.exit(dev.off(), add = TRUE)
    }

    return(g)
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

    if (!is.null(outPrefix)) on.exit(dev.off(), add = TRUE)

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
#' gtffile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#'     package = "GenomicPlot"
#' )
#' gff <- suppressMessages(RCAS::importGtf(saveObjectAsRds = TRUE, 
#'     filePath = gtffile))
#' txdb <- makeTxDbFromGRanges(gff)
#' txdb$user_genome <- "hg19"
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
#' gfeatures <- prepare_5parts_genomic_features(txdb,
#'     meta = TRUE, nbins = 100, fiveP = -1000,
#'     threeP = 1000, longest = TRUE, verbose = FALSE
#' )
#'
#' bamimportParams <- setImportParams(
#'     offset = -1, fix_width = 0, fix_point = "start", norm = TRUE,
#'     useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' alist <- list(
#'     "txdb" = txdb, "treat" = queryfiles, "control" = inputfiles,
#'     "feature" = gfeatures, "param" = bamimportParams
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
