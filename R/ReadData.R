#' @title set parameters for \code{handle_input} function
#'
#' @description This function save as a template for setting up import
#' parameters for reading NGS data, it provides default values for each
#' parameter.
#' @param offset an integer, -1 indicating the bam reads should be shrunk to the
#'  -1 position at the 5'end of the reads, which corresponds to the cross link
#'  site in iCLIP.
#' @param fix_width an integer, for bam file, defines how long the reads should
#' be extended from the start position, ignored when offset is not 0; for bed
#' files, defines the width of each interval centering on the `fix_point`.
#' @param fix_point a string in c("start", "end", "center") denoting the anchor
#'  point for extension, ignored when offset is not 0.
#' @param useScore logical, indicating whether the 'score' column of the bed
#'  file should be used in calculation of coverage.
#' @param outRle logical, indicating whether the output should be RleList
#'  objects or GRanges objects.
#' @param norm logical, indicating whether the output RleList should be
#'  normalized to RPM using library sizes.
#' @param genome a string denoting the genome name and version.
#' @param useSizeFactor logical, indicating whether the library size should be
#'  adjusted with a size factor, using the 'calcNormFactors' function in the
#'  edgeR package, only applicable to ChIPseq data.
#' @param saveRds logical, indicating whether the results of handle_input should
#'  be saved for fast reloading
#' @param val integer, indicating the column that will be used as score/value.
#' default 4 for bedGraph.
#' @param skip integer, indicating how many rows will be skipped before reading
#' in data, default 0.
#' @param chr a vector of string, denoting chromosomes to be included, like
#' c("chr1", "chr2", "chrX"), default NULL indicating all chromosomes will be
#' included.
#' @return a list of nine elements
#'
#' @author Shuye Pu
#'
#' @examples
#' importParams1 <- setImportParams()
#' importParams2 <- setImportParams(offset = -1, saveRds = TRUE)
#'
#' @export setImportParams
#'
setImportParams <- function(
        offset = 0,
        fix_width = 0,
        fix_point = "start",
        norm = FALSE,
        useScore = FALSE,
        outRle = TRUE,
        useSizeFactor = FALSE,
        saveRds = FALSE,
        genome = "hg19",
        val = 4,
        skip = 0,
        chr = NULL) {
    stopifnot(is.numeric(c(offset, fix_width)))
    stopifnot(fix_point %in% c("start", "center", "end"))
    stopifnot(is.logical(c(norm, useScore, outRle, useSizeFactor, saveRds)))
    stopifnot(is.character(genome))
    if(grepl("GRCh37", genome)) genome <- "hg19"
    if(grepl("GRCh38", genome)) genome <- "hg38"
    return(list(
        offset = offset, fix_width = fix_width, fix_point = fix_point,
        norm = norm, useScore = useScore, outRle = outRle,
        useSizeFactor = useSizeFactor, saveRds = saveRds,
        genome = genome, val = val, skip = skip, chr = chr
    ))
}

#' @title Handle import of NGS data with various formats
#'
#' @description This is a wrapper function for read NGS data in different file
#' formats, store the input data in a list of GRanges objects or RleList
#' objects. File names end in bed|bam|bw|bigwig|bigWig|BigWig|BW|BIGWIG are
#' recognized, and a named list of files with mixed formats are allowed.
#'
#' @param inputFiles a vector of strings denoting file names
#' @param importParams a list with the 9 elements: list(offset, fix_width,
#'  fix_point, useScore, outRle, norm, genome, useSizeFactor). Details are
#'  described in the documentation of \code{\link{setImportParams}} function
#' @param nc integer, number of cores for parallel processing
#' @param verbose logical, whether to output additional information
#'
#' @details when 'useScore' is TRUE, the score column of the bed file will be
#' used in the metadata column 'score' of the GRanges object, or the 'Values'
#' field of the RleList object. Otherwise the value 1 will be used instead. When
#' the intended use of the input bed is a reference feature, both 'useScore' and
#' 'outRle' should be set to FALSE.
#'
#' @return a list object with four elements, 'query' is a list GRanges objects
#' or RleList objects, 'size' is the library size, 'type' is the input file
#' type, 'weight' is the name of the metadata column to be used as weight for
#' coverage calculation
#'
#' @author Shuye Pu
#'
#' @examples
#' queryFiles1 <- system.file("extdata", "treat_chr19.bam",
#'     package = "GenomicPlot"
#' )
#' names(queryFiles1) <- "query"
#'
#' inputFiles1 <- system.file("extdata", "input_chr19.bam",
#'     package = "GenomicPlot"
#' )
#' names(inputFiles1) <- "input"
#'
#' bamimportParams <- setImportParams(
#'     offset = -1, fix_width = 0, fix_point = "start", norm = TRUE,
#'     useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' out_list <- handle_input(
#'     inputFiles = c(queryFiles1, inputFiles1),
#'     importParams = bamimportParams, verbose = TRUE, nc = 2
#' )
#'
#' queryFiles2 <- system.file("extdata", "test_wig_chr19_+.wig",
#'     package = "GenomicPlot"
#' )
#' names(queryFiles2) <- "test_wig"
#'
#' wigimportParams <- setImportParams(
#'     offset = 0, fix_width = 0, fix_point = "start", norm = FALSE,
#'     useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' out <- handle_input(queryFiles2, wigimportParams, verbose = TRUE)
#'
#' queryFiles3 <- system.file("extdata", "test_wig_chr19_+.bw",
#'     package = "GenomicPlot"
#' )
#' names(queryFiles3) <- "test_bw"
#'
#' out <- handle_input(c(queryFiles1, queryFiles2, queryFiles3),
#'     wigimportParams,
#'     verbose = TRUE
#' )
#'
#' @export handle_input

handle_input <- function(inputFiles,
                         importParams = NULL,
                         verbose = FALSE,
                         nc = 2) {
    stopifnot(all(file.exists(inputFiles)))
    if (is.null(names(inputFiles)) || any(names(inputFiles) == ""))
        stop("Each file must have a name attribute!")
    if (is.null(importParams)) importParams <- setImportParams()

    if (verbose) message("[handle_input] started ...\n")
    original_outRle <- importParams$outRle
    if (importParams$useSizeFactor && length(inputFiles) > 1) {
        importParams$outRle <- FALSE
        ## effective_size only take granges_list which produce more accurate
        ## normFactor than RleList
    }

    inputFUN <- function(funName, inputFile, importParams, verbose) {
        fileName <- basename(inputFile)
        dirName <- dirname(inputFile)
        if (file.exists(file.path(dirName, paste0(fileName, ".rds")))) {
            temp <- readRDS(file.path(dirName, paste0(fileName, ".rds")))
            if (identical(temp$param, importParams)) {
                out <- temp$Rle
                if (verbose) message("Cached .rds file is used\n")
            } else {
                if (!file.exists(inputFile)) stop("file does not exist,
                                       please check your file name and path!")
                out <- funName(inputFile = inputFile, importParams, verbose)
                if (importParams$saveRds && file.access(dirName,
                                                        mode = 2)[1] == 0) {
                    saveRDS(list(param = importParams, Rle = out),
                            file.path(dirName, paste0(fileName, ".rds")))
                    if (verbose) message("Cached .rds file is modified using
                                         new input parameters\n")
                }
            }
        } else {
            if (!file.exists(inputFile)) stop("file does not exist,
                                     please check your file name and path!")
            out <- funName(inputFile = inputFile, importParams, verbose)
            if (importParams$saveRds && file.access(dirName,
                                                    mode = 2)[1] == 0) {
                saveRDS(list(param = importParams, Rle = out),
                        file.path(dirName, paste0(fileName, ".rds")))
                if (verbose) message("Input data is cached as .rds file for
                                     fast reloading\n")
            }
        }
        return(out)
    }

    bed_suffix <- paste(
        c(paste0("\\.", c("bed", "BED", "Bed", "narrowPeak", "broadPeak"), "$"),
        paste0("\\.", c("bed", "BED", "Bed", "narrowPeak", "broadPeak"), "\\.gz$")),
        collapse = "|")
    bedGraph_suffix <- paste(
        c(paste0("\\.", c("bedGraph", "bedgraph", "BedGraph", "BEDGRAPH", "bg",
                        "BG"), "$"),
          paste0("\\.", c("bedGraph", "bedgraph", "BedGraph", "BEDGRAPH", "bg",
                                                     "BG"), "\\.gz$")),
        collapse = "|")
    bigwig_suffix <- paste(
        paste0("\\.", c("bw", "bigwig", "Bigwig", "bigWig", "BigWig", "BIGWIG",
                        "BW"), "$"), collapse = "|")
    outlist <- lapply(names(inputFiles), function(aname) {
        inputFile <- inputFiles[aname]
        if (grepl(bed_suffix, inputFile)) {
            fileType <- "bed"
            if (verbose) message("Reading ", fileType, "file: ", inputFile)
            out <- inputFUN(handle_bed, inputFile = inputFile, importParams,
                            verbose)
        }else if (grepl(bedGraph_suffix, inputFile)) {
            fileType <- "bedGraph"
            if (verbose) message("Reading ", fileType, "file: ", inputFile)
            out <- inputFUN(handle_bedGraph, inputFile = inputFile, importParams,
                            verbose)
        } else if (grepl("\\.bam$|\\.BAM$|\\.Bam$", inputFile)) {
            fileType <- "bam"
            if (verbose) message("Reading ", fileType, "file: ", inputFile)
            out <- inputFUN(handle_bam, inputFile = inputFile, importParams,
                            verbose)
        } else if (grepl("\\.wig$|\\.WIG$|\\.Wig$", inputFile)) {
            fileType <- "wig"
            if (verbose) message("Reading ", fileType, "file: ", inputFile)
            out <- inputFUN(handle_wig, inputFile = inputFile, importParams,
                            verbose)
        } else if (grepl(bigwig_suffix, inputFile)) {
            fileType <- "bw"
            if (verbose) message("Reading ", fileType, "file: ", inputFile)
            out <- inputFUN(handle_bw, inputFile = inputFile, importParams,
                            verbose)
        } else {
            stop("The format of file is not supported, please convert it to one
                 of the following format: bed, bam, wig, bigwig")
        }

        out
    })

    names(outlist) <- names(inputFiles)
    importParams$outRle <- original_outRle # if modified, restore

    ## modify library size
    if (importParams$useSizeFactor && (length(inputFiles) > 1)) {
        outlist <- effective_size(outlist = outlist,
                                  outRle = importParams$outRle,
                                  genome = importParams$genome, nc = nc)
    }

    ## compute RPM
    if (importParams$norm) {
        outlist <- lapply(outlist, function(x) {
            if (importParams$outRle) {
                x$query <- x$query * 1e6 / x$size
            } else {
                score(x$query) <- x$query$score * 1e6 / x$size
            }
            x
        })
    }
    if (verbose) message("[handle_input] finished!\n")
    invisible(outlist)
}


#' @title Normalize sample library size to effective size
#'
#' @description This is a helper function for handle_input.
#' edgeR::calcNormFactors function is used to estimate normalizing factors,
#' which is used to multiply library sizes.
#'
#' @param outlist a list of list objects with four elements, 'query' is a
#'  GRanges object, 'size' is the library size, 'type' is the input file type,
#'  'weight' is the name of the metadata column
#' @param outRle logical, indicating whether the 'query' element of the output
#'  should be an RleList object or a GRanges object
#' @param genome a string denoting the genome name and version
#' @param nc integer, number of cores for parallel processing
#' @param verbose logical, whether to output additional information
#'
#' @return a list of list objects with four elements ('query', 'size', 'type',
#'  'weight'), with the 'size' element modified.
#'
#' @author Shuye Pu
#'
#' @examples
#' queryFiles <- system.file("extdata", "chip_treat_chr19.bam",
#'     package = "GenomicPlot"
#' )
#' names(queryFiles) <- "query"
#'
#' inputFiles <- system.file("extdata", "chip_input_chr19.bam",
#'     package = "GenomicPlot"
#' )
#' names(inputFiles) <- "input"
#'
#' chipImportParams <- setImportParams(
#'     offset = 0, fix_width = 150, fix_point = "start", norm = TRUE,
#'     useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' out_list <- handle_input(
#'     inputFiles = c(queryFiles, inputFiles),
#'     importParams = chipImportParams, verbose = TRUE, nc = 2
#' )
#'
#' out <- effective_size(out_list, outRle = TRUE)
#'
#' @export effective_size
#'
effective_size <- function(outlist,
                           outRle,
                           genome = "hg19",
                           nc = 2,
                           verbose = FALSE) {
    if (verbose) message("Estimating size factor\n")

    seqi <- set_seqinfo(genome)

    grange_list <- lapply(outlist, function(x) x$query)

    tilewidth <- 100000
    tileBins <- tileGenome(seqi, tilewidth = tilewidth,
                           cut.last.tile.in.chrom = TRUE)

    score_list <- parallel_countOverlaps(grange_list, tileBins, nc = nc)

    mat <- data.matrix(bind_cols(score_list))
    mat[is.na(mat)] <- 0
    mat <- round(mat * tilewidth)
    mat <- mat[apply(mat, 1, sum) > 0, ]

    # normFactor <- DESeq2::estimateSizeFactorsForMatrix(mat)
    lib.size <- vapply(outlist, function(x) x$size, numeric(1))
    normFactor <- edgeR::calcNormFactors(mat, lib.size = lib.size)

    ## update library size with normFactor
    normlist <- mapply(x = outlist, y = normFactor, function(x, y) {
        x$size <- as.integer(x$size * y)
        x
    }, SIMPLIFY = FALSE, USE.NAMES = TRUE)

    lib.size.adj <- vapply(normlist, function(x) x$size, numeric(1))

    if (outRle) {
        normlist <- lapply(normlist, function(x) {
            x$query <- coverage(x$query, weight = x$weight)
            x
        })
    }

    names(normFactor) <- NULL

    if (verbose) {
        message("Library size before adjusting: ", lib.size)
        message("Library size after adjusting: ", lib.size.adj)
        message("Library size normalizing factors: ", normFactor)
    }

    invisible(normlist)
}

#' @title Handle files in bed|narrowPeak|broadPeak format
#'
#' @description This is a function for read peaks data in bed format, store the
#' input data in a list of GRanges objects or RleList objects.
#'
#' @param inputFile a string denoting path to the input file
#' @param importParams a list of parameters, refer to \code{\link{handle_input}}
#'  for details
#' @param verbose logical, whether to output additional information
#'
#' @return a list object with four elements, 'query' is a list GRanges objects
#' or RleList objects, 'size' is the library size, 'type' is the input file
#' type, 'weight' is the name of the metadata column to be used as weight for
#' coverage calculation
#'
#' @author Shuye Pu
#'
#' @examples
#' queryFiles <- system.file("extdata", "test_chip_peak_chr19.narrowPeak",
#'     package = "GenomicPlot"
#' )
#' names(queryFiles) <- "narrowPeak"
#'
#' bedimportParams <- setImportParams(
#'     offset = 0, fix_width = 100, fix_point = "center", norm = FALSE,
#'     useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' out <- handle_bed(queryFiles, bedimportParams, verbose = TRUE)
#' lapply(out$query, sum)
#'
#' @export handle_bed

handle_bed <- function(inputFile,
                       importParams = NULL,
                       verbose = FALSE) {
    stopifnot(file.exists(inputFile))
    if (is.null(names(inputFile)) || names(inputFile) == "")
        stop("Each file must have a name attribute!")
    if (verbose) message("[handle_bed] started ...\n")
    beddata <- read.delim2(inputFile, header = FALSE, comment.char = "#")

    nco <- ncol(beddata)
    if (nco > 6 && verbose) {
        message("The input file ", inputFile, " have more than 6 columns,
                only the first 6 columns will be used!\n")
    } else if (nco < 6 && verbose) {
        message("The input file ", inputFile, " have less than 6 columns!\n")
    }

    beddata <- type.convert(beddata[, seq_len(min(6, nco))], as.is = TRUE)
    ## ignore extra columns, which cause problem in import.bed()
    standard_col <- c("chr", "start", "end", "name", "score", "strand")
    colnames(beddata) <- standard_col[seq_len(min(6, ncol(beddata)))]

    if(!is.null(importParams$chr)){
        beddata <- beddata %>%
            dplyr::filter(chr %in% importParams$chr)
    }

    queryRegions <- makeGRangesFromDataFrame(beddata, keep.extra.columns = TRUE,
                                             starts.in.df.are.0based = TRUE)

    if("name" %in% colnames(beddata)) {
        if (sum(duplicated(beddata$name)) > 0) {
            ## if the names are not unique, force them to be unique
            message("bed interval names are modified as they are not unique!")
            names(queryRegions) <- paste(beddata$name,
                                         seq_along(beddata$name), sep = "_")
        }else{
            names(queryRegions) <- beddata$name
        }
    }else{
        names(queryRegions) <- paste("peak", seq_along(beddata[,1]), sep = "_")
    }


    if (ncol(beddata) < 6) strand(queryRegions) <- "*"

    if (importParams$fix_width > 0) {
        queryRegions <- resize(queryRegions,
                               width = 1,
                               fix = rep(importParams$fix_point, length(queryRegions)),
                               ignore.strand = FALSE
        )
        queryRegions <- flank(queryRegions,
                              width = as.integer(importParams$fix_width/2),
                              start = TRUE,
                              both = TRUE,
                              ignore.strand = FALSE
        )
    }
    weight_col <- "score"

    if (!importParams$useScore || ncol(beddata) < 5) {
        score(queryRegions) <- 1
    }

    ## make input comply with GenomeInfoDb
    if("NCBI" %in% seqlevelsStyle(queryRegions)){
        seqlevelsStyle(queryRegions) <- "UCSC"
    }
    seqInfo <- set_seqinfo(importParams$genome)

    queryRegions <- queryRegions[as.vector(seqnames(queryRegions))
                                 %in% seqnames(seqInfo)]
    seqlevels(queryRegions) <- seqlevels(seqInfo)
    seqinfo(queryRegions) <- seqInfo

    libsize <- sum(queryRegions$score, na.rm = TRUE)

    if (importParams$outRle) {
        queryRegions <- coverage(queryRegions, weight = weight_col)
        seqinfo(queryRegions) <- seqInfo
    }

    if (verbose) message("[handle_bed] finished!\n")
    invisible(list("query" = queryRegions, "size" = libsize, "type" = "bed",
                   "weight" = weight_col))
}

#' @title Handle files in bedGraph format
#'
#' @description This is a function for read peaks data in bedGraph format,
#' store the input data in a list of GRanges objects or RleList objects.
#'
#' @param inputFile a string denoting path to the input file
#' @param importParams a list of parameters, refer to \code{\link{handle_input}}
#'  for details
#' @param verbose logical, whether to output additional information
#'
#' @return a list object with four elements, 'query' is a list GRanges objects
#' or RleList objects, 'size' is the library size, 'type' is the input file
#' type, 'weight' is the name of the metadata column to be used as weight for
#' coverage calculation
#'
#' @author Shuye Pu
#'
#' @examples
#' queryFiles <- system.file("extdata", "test_chr19.bedGraph",
#'     package = "GenomicPlot"
#' )
#' names(queryFiles) <- "chipPeak"
#'
#' importParams <- setImportParams(
#'     offset = 0, fix_width = 0, fix_point = "start", norm = FALSE,
#'     useScore = TRUE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19",
#'     val = 4, skip = 1
#' )
#'
#' out <- handle_bedGraph(queryFiles, importParams, verbose = TRUE)
#' out$query
#'
#' @export handle_bedGraph

handle_bedGraph <- function(inputFile,
                            importParams = NULL,
                            verbose = FALSE) {
    stopifnot(file.exists(inputFile))
    if (is.null(names(inputFile)) || names(inputFile) == "")
        stop("Each file must have a name attribute!")

    if (verbose) message("[handle_bedGraph] started ...\n")
    beddata <- read.delim2(inputFile, header = FALSE, comment.char = "#",
                           skip = importParams$skip)

    nco <- ncol(beddata)
    if (nco > 4 && verbose) {
        message("The input file ", inputFile, " have more than 4 columns,
                score column has to be specified!\n")
    }
    beddata <- type.convert(beddata[, c(1,2,3, importParams$val)], as.is = TRUE)
    ## ignore extra columns, which cause problem in import.bed()
    standard_col <- c("chr", "start", "end", "score")
    colnames(beddata) <- standard_col

    if(!is.null(importParams$chr)){
        beddata <- beddata %>%
            dplyr::filter(chr %in% importParams$chr)
    }

    queryRegions <- makeGRangesFromDataFrame(beddata, keep.extra.columns = TRUE,
                                             starts.in.df.are.0based = TRUE,
                                             ignore.strand = TRUE)

    names(queryRegions) <- paste("region", seq_along(beddata[,1]), sep = "_")

    if (importParams$fix_width > 0) {
        queryRegions <- resize(queryRegions,
                               width = 1,
                               fix = rep(importParams$fix_point, length(queryRegions)),
                               ignore.strand = FALSE
        )
        queryRegions <- flank(queryRegions,
                              width = as.integer(importParams$fix_width/2),
                              start = TRUE,
                              both = TRUE,
                              ignore.strand = FALSE
        )
    }
    weight_col <- "score"

    if (!importParams$useScore) {
        score(queryRegions) <- 1
    }

    ## make input comply with GenomeInfoDb
    if("NCBI" %in% seqlevelsStyle(queryRegions)){
        seqlevelsStyle(queryRegions) <- "UCSC"
    }
    seqInfo <- set_seqinfo(importParams$genome)

    queryRegions <- queryRegions[as.vector(seqnames(queryRegions))
                                 %in% seqnames(seqInfo)]
    seqlevels(queryRegions) <- seqlevels(seqInfo)
    seqinfo(queryRegions) <- seqInfo

    libsize <- sum(queryRegions$score, na.rm = TRUE)

    if (importParams$outRle) {
        queryRegions <- coverage(queryRegions, weight = weight_col)
        seqinfo(queryRegions) <- seqInfo
    }

    if (verbose) message("[handle_bedGraph] finished!\n")
    invisible(list("query" = queryRegions, "size" = libsize, "type" = "bedGraph",
                   "weight" = weight_col))
}


#' @title Handle files in bam format
#'
#' @description This is a function for read NGS reads data in bam format, store
#' the input data in a list of GRanges objects or RleList objects. For
#' paired-end reads, only take the second read in a pair, assuming which is the
#' sense read for strand-specific RNAseq.
#'
#' @param inputFile a string denoting path to the input file
#' @param importParams a list of parameters, refer to \code{\link{handle_input}}
#'  for details
#' @param verbose logical, whether to output additional information
#'
#' @details The reads are filtered using mapq score >= 10 by default, only
#'  mapped reads are counted towards library size.
#'
#' @return a list object with four elements, 'query' is a list GRanges objects
#' or RleList objects, 'size' is the library size, 'type' is the input file
#' type, weight' is the name of the metadata column to be used as weight for
#' coverage calculation
#'
#' @author Shuye Pu
#'
#' @examples
#'
#' queryFiles <- system.file("extdata", "treat_chr19.bam",
#'     package = "GenomicPlot"
#' )
#' names(queryFiles) <- "query"
#'
#' bamimportParams <- setImportParams(
#'     offset = -1, fix_width = 0, fix_point = "start", norm = TRUE,
#'     useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' out <- handle_bam(
#'     inputFile = queryFiles, importParams = bamimportParams, verbose = TRUE
#' )
#'
#' @export handle_bam
#'
handle_bam <- function(inputFile, importParams = NULL, verbose = FALSE) {

    stopifnot(file.exists(inputFile))
    if (is.null(names(inputFile)) || names(inputFile) == "")
        stop("Each file must have a name attribute!")

    if (verbose) message("[handle_bam] started ...\n")
    paired.end <- testPairedEndBam(inputFile)
    if (paired.end) {
        flag <- scanBamFlag(isPaired = TRUE, isProperPair = TRUE,
                            isUnmappedQuery = FALSE, hasUnmappedMate = FALSE,
                            isSecondMateRead = TRUE)
        ## assume the second read in a pair is the sense read for
        ## strand-specific RNAseq
        param <- ScanBamParam(mapqFilter = 10, flag = flag)
        if (verbose) message("Paired end bam file detected!\n")
    } else {
        param <- ScanBamParam(mapqFilter = 10)
    }
    if (verbose) message("\nBam file ", inputFile, " is loaded\n")

    ga <- readGAlignments(inputFile, use.names = TRUE, param = param)
    libsize <- sum(idxstatsBam(inputFile)$mapped, na.rm = TRUE)
    if (importParams$offset != 0) {
        ##for iCLIP, use offset = -1 to get the 5'-end -1 position of the reads,
        ##which is the crosslink sites for iCLIP reads
        queryRegions <- flank(granges(ga),
                               width = -1 * importParams$offset,
                               start = TRUE, both = FALSE,
                               ignore.strand = FALSE)
        queryRegions <- resize(queryRegions, width = 1, fix = "start",
                               ignore.strand = FALSE)
        score(queryRegions) <- 1
    } else if (importParams$fix_width > 0) {
        queryRegions <- resize(granges(ga), width = importParams$fix_width,
                               fix = "start", ignore.strand = FALSE)
        score(queryRegions) <- 1
    } else {
        queryRegions <- unlist(grglist(ga))
        score(queryRegions) <- 1
    }

    if(!is.null(importParams$chr)){
        queryRegions <- queryRegions %>%
            plyranges::filter(seqnames %in% importParams$chr)
    }

    weight_col <- "score"

    ## make input comply with GenomeInfoDb, use cached chromInfo to avoid
    ## dependency on UCSC web service
    if("NCBI" %in% seqlevelsStyle(queryRegions)){
        seqlevelsStyle(queryRegions) <- "UCSC"
    }
    seqInfo <- set_seqinfo(importParams$genome)
    queryRegions <- queryRegions[as.vector(seqnames(queryRegions))
                                 %in% seqnames(seqInfo)]
    seqlevels(queryRegions) <- seqlevels(seqInfo)
    seqinfo(queryRegions) <- seqInfo


    if (importParams$outRle) {
        queryRegions <- coverage(queryRegions, weight = weight_col)
        seqinfo(queryRegions) <- seqInfo
    }

    if (verbose) message("[handle_bam] finished!\n")
    invisible(list("query" = queryRegions, "size" = libsize, "type" = "bam",
                   "weight" = weight_col))
}


#' @title Handle files in bw|bigwig|bigWig|BigWig|BW|BIGWIG format
#'
#' @description This is a function for read NGS coverage data in bigwig format,
#' store the input data in a list of GRanges objects or RleList objects. The
#' input bw file can be stranded or non-stranded. Library size is calculate as
#' the sum of all coverage.
#'
#' @param inputFile a string denoting path to the input file
#' @param importParams a list of parameters, refer to \code{\link{handle_input}}
#'  for details
#' @param verbose logical, whether to output additional information
#'
#' @details For stranded files, forward and reverse strands are stored in
#' separate files, with '+' or 'p' in the forward strand file name and '-' or
#' 'm' in the reverse strand  file name.
#'
#' @return a list object with four elements, 'query' is a list GRanges objects
#' or RleList objects, 'size' is the estimated library size, 'type' is the
#' input file type, weight' is the name of the metadata column to be used as
#' weight for coverage calculation
#'
#' @author Shuye Pu
#'
#' @examples
#'
#' queryFiles <- system.file("extdata", "test_wig_chr19_+.bw",
#'     package = "GenomicPlot"
#' )
#' names(queryFiles) <- "test_bw"
#'
#' wigimportParams <- setImportParams(
#'     offset = 0, fix_width = 0, fix_point = "start", norm = FALSE,
#'     useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' out <- handle_bw(queryFiles, wigimportParams, verbose = TRUE)
#'
#' @export handle_bw
#'
handle_bw <- function(inputFile, importParams, verbose = FALSE) {
    stopifnot(file.exists(inputFile))
    if (is.null(names(inputFile)) || names(inputFile) == "")
        stop("Each file must have a name attribute!")

    if (verbose) message("[handle_bw] started ...\n")
    weight_col <- "score"
    libsize <- NULL

    neg_file <- find_mate(inputFile, verbose)
    stranded <- ifelse(!is.null(neg_file), TRUE, FALSE)

    if (stranded) {
        posGR <- import.bw(inputFile)
        strand(posGR) <- "+"

        negGR <- import.bw(neg_file)
        strand(negGR) <- "-"

        queryRegions <- append(posGR, negGR)
    } else {
        iGR <- import.bw(inputFile)

        queryRegions <- iGR
        strand(queryRegions) <- "*"
    }

    if(!is.null(importParams$chr)){
        queryRegions <- queryRegions %>%
            plyranges::filter(seqnames %in% importParams$chr)
    }

    ## make input comply with GenomeInfoDb
    if("NCBI" %in% seqlevelsStyle(queryRegions)){
        seqlevelsStyle(queryRegions) <- "UCSC"
    }
    seqInfo <- set_seqinfo(importParams$genome)
    queryRegions <- queryRegions[as.vector(seqnames(queryRegions))
                                 %in% seqnames(seqInfo)]
    GenomeInfoDb::seqlevels(queryRegions) <- GenomeInfoDb::seqlevels(seqInfo)
    GenomeInfoDb::seqinfo(queryRegions) <- seqInfo

    libsize <- sum(score(queryRegions) * width(queryRegions),
                              na.rm = TRUE)/100 # assuming read length is 100

    if (importParams$outRle) {
        queryRegions <- coverage(queryRegions, weight = weight_col)
        GenomeInfoDb::seqinfo(queryRegions) <- seqInfo
    }

    if (verbose) message("[handle_bw] finished!\n")
    invisible(list("query" = queryRegions, "size" = libsize, "type" = "bw",
                   "weight" = weight_col))
}

#' @title Handle files in wig format
#'
#' @description This is a function for read NGS coverage data in wig format,
#' store the input data in a list of GRanges objects or RleList objects.
#' The input wig file can be stranded or non-stranded. Library size is calculate
#' as the sum of all coverage.
#'
#' @param inputFile a string denoting path to the input file
#' @param importParams a list of parameters, refer to \code{\link{handle_input}}
#'  for details
#' @param verbose logical, whether to output additional information
#'
#' @details For stranded files, forward and reverse strands are stored in
#' separate files, with '+' or 'p' in the forward strand file name and '-' or
#' 'm' in the reverse strand file name.
#'
#' @return a list object with four elements, 'query' is a list GRanges objects
#' or RleList objects, 'size' is the library size, 'type' is the input file
#' type, 'weight' is the name of the metadata column to be used as weight for
#' coverage calculation
#'
#' @author Shuye Pu
#'
#' @examples
#'
#' queryFiles <- system.file("extdata", "test_wig_chr19_+.wig",
#'     package = "GenomicPlot"
#' )
#' names(queryFiles) <- "test_wig"
#'
#' wigimportParams <- setImportParams(
#'     offset = 0, fix_width = 0, fix_point = "start", norm = FALSE,
#'     useScore = FALSE, outRle = TRUE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' out <- handle_wig(queryFiles, wigimportParams, verbose = TRUE)
#'
#' @export handle_wig
#'
handle_wig <- function(inputFile,
                       importParams,
                       verbose = FALSE) {
    stopifnot(file.exists(inputFile))
    if (is.null(names(inputFile)) || names(inputFile) == "")
        stop("Each file must have a name attribute!")

    if (verbose) message("[handle_wig] started ...\n")
    neg_file <- find_mate(inputFile, verbose)
    stranded <- ifelse(!is.null(neg_file), TRUE, FALSE)

    seqinfo <- set_seqinfo(importParams$genome)
    wigToBigWig(inputFile, seqinfo, clip = TRUE)

    if (stranded) {
        wigToBigWig(neg_file, seqinfo, clip = TRUE)
    }

    bwfile <- gsub("\\.wig$|\\.WIG$|\\.Wig$", "\\.bw", inputFile)

    out <- handle_bw(bwfile, importParams, verbose)

    if (verbose) message("[handle_wig] finished!\n")
    invisible(out)
}


#' @title Find wig/bw file for the negative strand
#' @description Find the file name of the negative strand, if a .wig/bw file for
#' positive strand if provided, by looking for file names with one character
#' difference. If no negative strand file is found, assume the input .wig/bw
#' file is non-stranded
#'
#' @param inputFile path to a .wig/bw file, presumably for positive strand
#' @param verbose logical, whether to output additional information
#'
#' @return path to the negative .wig/bw file or NULL
#'
#' @author Shuye Pu
#'
#' @examples
#'
#' queryFile <- system.file("extdata", "test_wig_chr19_+.wig",
#'     package = "GenomicPlot"
#' )
#' names(queryFile) <- "test_wig"
#'
#' out <- GenomicPlot:::find_mate(inputFile = queryFile, verbose = TRUE)
#'
#' @keywords internal


find_mate <- function(inputFile,
                      verbose = FALSE) {
    fileName <- basename(inputFile)
    dirName <- dirname(inputFile)

    fch_v <- unlist(strsplit(fileName, split = ""))
    otherFiles <- list.files(dirName)

    mate <- NULL
    for (afile in otherFiles) {
        och_v <- unlist(strsplit(afile, split = ""))
        if (length(fch_v) == length(och_v) && sum(fch_v != och_v) == 1) {
            mate <- file.path(dirName, afile)
            diff <- base::setdiff(och_v, fch_v)
            if (verbose) message("Mate found:\n", inputFile, "\n", mate, "\n",
                                 diff, "\n")
        }
    }
    invisible(mate)
}
