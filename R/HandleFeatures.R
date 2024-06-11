#' @title Extract the longest transcript for each protein-coding genes
#
#' @description Gene level computations require selecting one transcript per
#' gene to avoid bias by genes with multiple isoforms. In ideal case, the most
#' abundant transcript (principal or canonical isoform) should be chosen.
#' However, the most abundant isoform may vary depending on tissue type or
#' physiological condition, the longest transcript is usually the principal
#' isoform, and alternatively spliced isoforms are not. This method get the
#' longest transcript for each gene. The longest transcript is defined as the
#' isoform that has the longest transcript length. In case of tie, the one
#' with longer CDS is selected. If the lengths of CDS tie again, the transcript
#' with smaller id is selected arbitrarily.
#'
#' @param txdb a TxDb object defined in the GenomicFeatures package
#'
#' @return a dataframe of transcript information with the following columns:
#'  "tx_id tx_name gene_id nexon tx_len cds_len utr5_len utr3_len"
#' @author Shuye Pu
#'
#' @examples
#'
#' gtfFile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#'     package = "GenomicPlot"
#' )
#'
#' txdb <- custom_TxDb_from_GTF(gtfFile, genome = "hg19")
#' longestTx <- extract_longest_tx(txdb)
#'
#' @export extract_longest_tx
#'
#'
extract_longest_tx <- function(txdb) {
    message("[extract_longest_tx]")
    tl <- GenomicFeatures::transcriptLengths(txdb, with.utr5_len = TRUE,
                                             with.cds_len = TRUE,
                                             with.utr3_len = TRUE)
    pc <- tl[tl$cds_len > 0, ] # pc stands for protein-coding

    ## choose tx with longest length for each gene
    longest_tx <- aggregate(pc$tx_len, list(pc$gene_id), max)
    colnames(longest_tx) <- c("gene_id", "tx_len")
    tx_longest <- merge(pc, longest_tx)

    ## for two tx of the same gene, if tx_len are same, choose longer cds_len
    longest_cds <- aggregate(tx_longest$cds_len, list(tx_longest$gene_id,
                                                      tx_longest$tx_len), max)
    colnames(longest_cds) <- c("gene_id", "tx_len", "cds_len")
    longest_cdstx <- merge(tx_longest, longest_cds)

    dup_genes <- longest_cdstx$gene_id[duplicated(longest_cdstx$gene_id)]
    dup_genestx <- longest_cdstx[longest_cdstx$gene_id %in% dup_genes, ]

    ## for two tx of the same gene, if both cds_len and tx_len are the same,
    ## choose smaller tx_id (this is arbitrary)
    longest_cdstx_id <- aggregate(longest_cdstx$tx_id,
                                  list(longest_cdstx$gene_id,
                                       longest_cdstx$cds_len,
                                       longest_cdstx$tx_len), min)
    colnames(longest_cdstx_id) <- c("gene_id", "cds_len", "tx_len", "tx_id")
    longest_cdstxid <- merge(longest_cdstx, longest_cdstx_id)
    length(longest_cdstxid$gene_id[duplicated(longest_cdstxid$gene_id)])
    summary(longest_cdstxid)
    head(longest_cdstxid)

    npc <- tl %>%
        filter(!gene_id %in% longest_cdstxid$gene_id)
    ## npc stands for non-protein-coding

    ## choose tx with longest length for each non-protein-coding gene
    longest_npc <- aggregate(npc$tx_len, list(npc$gene_id), max)
    colnames(longest_npc) <- c("gene_id", "tx_len")
    tx_longest <- merge(npc, longest_npc)
    tx_longest_id <- aggregate(tx_longest$tx_id, list(tx_longest$gene_id,
                                                      tx_longest$tx_len), min)
    colnames(tx_longest_id) <- c("gene_id", "tx_len", "tx_id")
    tx_longestid <- merge(tx_longest, tx_longest_id)

    longest_txid <- rbind(longest_cdstxid[, colnames(tl)],
                          tx_longestid[, colnames(tl)])

    invisible(longest_txid)
}

#'
#'
#' @title Set standard chromosome size of model organisms
#'
#' @description This is a helper function for making Seqinfo objects, which is
#' a components of GRanges and TxDb objects. It also serves to unify
#' seqlevels between GRanges and TxDb objects. Mitochondrial chromosome is not
#' included.
#'
#' @param genome a string denoting the genome name and version
#'
#' @return a Seqinfo object defined in the GenomeInfoDb package.
#'
#' @author Shuye Pu
#'
#' @examples
#'
#' out <- set_seqinfo(genome = "hg19")
#'
#' @export set_seqinfo
#'
set_seqinfo <- function(genome = "hg19") {
    message("[set_seqinfo]")
    if(grepl("GRCh37", genome)) genome <- "hg19"
    if(grepl("GRCh38", genome)) genome <- "hg38"

    chromInfo <- circlize::read.chromInfo(species = genome)$df
    seqi <- Seqinfo(seqnames = chromInfo$chr, seqlengths = chromInfo$end,
                    isCircular = rep(FALSE, nrow(chromInfo)),
                    genome = genome)
    return(seqi)
}

#' @title Make custom TxDb object from a GTF/GFF file
#'
#' @description This is a helper function for creating custom TxDb object from
#' a GTF/GFF file. Mitochondrial chromosome is excluded.
#'
#' @param gtfFile path to a gene annotation gtf file
#' @param genome a string denoting the genome name and version
#'
#' @return a TxDb object defined in the GenomicFeatures package.
#'
#' @author Shuye Pu
#'
#' @examples
#'
#' gtfFile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#'     package = "GenomicPlot"
#' )
#'
#' txdb <- custom_TxDb_from_GTF(gtfFile, genome = "hg19")
#'
#' @export custom_TxDb_from_GTF
#'
custom_TxDb_from_GTF <- function(gtfFile, genome = "hg19") {
    if(grepl("GRCh37", genome)) genome <- "hg19"
    if(grepl("GRCh38", genome)) genome <- "hg38"

    message("[custom_TxDb_from_GTF]")
    chromInfo <- set_seqinfo(genome)
    gff <- RCAS::importGtf(saveObjectAsRds = TRUE, filePath = gtfFile)
    gff <- gff[as.vector(seqnames(gff)) %in% seqlevels(chromInfo)]
    seqlevels(gff) <- seqlevels(chromInfo)
    seqinfo(gff) <- chromInfo
    txdb <- txdbmaker::makeTxDbFromGRanges(gff)

    return(txdb)
}

#' @title Extract genomic features from TxDb object
#
#' @description Extract genomic coordinates and make bed or bed 12 files from a
#' TxDb object for a variety of annotated genomic features. The output of this
#' function is a list. The first element of the list is a GRanges object that
#' provide the start and end information of the feature. The second element is a
#' GRangesList providing information for sub-components. The third element is
#' the name of a bed file.
#' @details For "utr3", "utr5", "cds" and "transcript", the GRanges object
#' denotes the start and end of the feature in one transcript, and the range is
#' named by the transcript id and may span introns; the GrangesList object is a
#' list of exons comprising each feature and indexed on transcript id. The bed
#' file is in bed12 format.
#'      For "exon" and "intron", the GRanges object denotes unnamed ranges of
#' individual exon and intron, and the GrangesList object is a list of exons or
#' introns belonging to one transcript and indexed on transcript id. The bed
#' file is in bed6 format.
#'      For "gene", both GRanges object and GRangesList object have the same
#' ranges and names. The bed file is in bed6 format.
#'
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param featureName one of the genomic feature in c("utr3", "utr5", "cds",
#'  "intron", "exon", "transcript", "gene")
#' @param featureSource the name of the gtf/gff3 file or the online database
#'  from which txdb is derived, used as name of output file
#' @param export logical, indicating if the bed file should be produced
#' @param longest logical, indicating whether the output should be limited to
#'  the longest transcript of each gene
#' @param protein_coding logical, indicating whether to limit to protein_coding
#'  genes
#'
#' @return a list of three objects, the first is a GRanges object, the second is
#'  a GRangesList object, the last is the output file name if export is TRUE.
#' @author Shuye Pu
#'
#' @examples
#'
#' gtfFile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#'     package = "GenomicPlot"
#' )
#'
#' txdb <- custom_TxDb_from_GTF(gtfFile, genome = "hg19")
#'
#' output <- get_genomic_feature_coordinates(txdb,
#'     featureName = "cds", featureSource = "gencode",
#'     export = FALSE, longest = TRUE, protein_coding = TRUE
#' )
#'
#' @export get_genomic_feature_coordinates
#'
get_genomic_feature_coordinates <- function(txdb,
                                            featureName,
                                            featureSource = NULL,
                                            export = FALSE,
                                            longest = FALSE,
                                            protein_coding = FALSE) {
    stopifnot(featureName %in% c("utr3", "utr5", "cds", "intron", "exon",
                                 "transcript", "gene"))
    message("[get_genomic_feature_coordinates] for ", featureName)
    if (is.null(featureSource)) {
        featureSource <- "unknown_source"
    }
    outfile <- NULL

    seqInfo <- seqinfo(txdb)

    tl <- transcriptLengths(txdb, with.utr5_len = TRUE, with.cds_len = TRUE,
                            with.utr3_len = TRUE)
    tl_protein_coding <- tl[tl$cds_len > 0, ]

    feature <- NULL

    if (featureName == "utr3") {
        feature <- threeUTRsByTranscript(txdb, use.name = TRUE) # grl
    } else if (featureName == "utr5") {
        feature <- fiveUTRsByTranscript(txdb, use.name = TRUE) # grl
    } else if (featureName == "intron") {
        feature <- intronsByTranscript(txdb, use.name = TRUE) # grl
    } else if (featureName == "exon") {
        feature <- exonsBy(txdb, by = "tx", use.name = TRUE) # grl
    } else if (featureName == "cds") {
        feature <- cdsBy(txdb, by = "tx", use.name = TRUE) # grl
    } else if (featureName == "transcript") {
        feature <- exonsBy(txdb, by = "tx", use.name = TRUE) # grl
    } else if (featureName == "gene") {
        feature <- genes(txdb) # gr
    } else {
        stop("Feature is not defined!")
    }

    if (featureName == "gene") {
       if (protein_coding) feature <- feature[unique(tl_protein_coding$gene_id)]
    } else {
        if (protein_coding) feature <- feature[names(feature) %in%
                                                   tl_protein_coding$tx_name]
    }

    seqinfo(feature) <- seqInfo
    if (longest) {
        longest_tx <- extract_longest_tx(txdb)
        if (featureName == "gene") {
            feature_longest <- feature[names(feature) %in% longest_tx$gene_id]
            seqinfo(feature_longest) <- seqInfo
        } else {
            feature_longest <- feature[names(feature) %in% longest_tx$tx_name]
            seqinfo(feature_longest) <- seqInfo
        }

        if (featureName %in% c("cds", "utr5", "utr3", "transcript")) {
            gr_feature_longest <- rtracklayer::asBED(feature_longest)
            ## convert each element of Grangeslist to one Grange with blocks
            ## info as metadata
            names(gr_feature_longest) <- gr_feature_longest$name
            seqinfo(gr_feature_longest) <- seqInfo
            if (export) {
                outfile <- paste(featureSource, "_", featureName,
                                 "_longest.bed12", sep = "")
                export.bed(gr_feature_longest, outfile)
            }
        }
        if (featureName %in% c("gene", "intron", "exon")) {
            if (featureName == "gene") {
                gr_feature_longest <- feature_longest # gr
                feature_longest <- as(split(feature_longest,
                                            f = feature_longest$gene_id),
                                      "GRangesList") # grl
            } else {
                gr_feature_longest <- unlist(feature_longest, use.names = TRUE)
                ## convert each element of Grangeslist to multiple Granges
                ## if the names are not unique, force them to be unique
                if (sum(duplicated(names(gr_feature_longest))) > 0) {
                   names(gr_feature_longest) <- paste(names(gr_feature_longest),
                        seq_along(names(gr_feature_longest)), sep = "_")
                }
            }

            seqinfo(gr_feature_longest) <- seqInfo
            if (export) {
                outfile <- paste(featureSource, "_", featureName,
                                 "_longest.bed", sep = "")
                export.bed(gr_feature_longest, outfile)
            }
        }
        message("[get_genomic_feature_coordinates] for ", featureName, " finished!\n")
        invisible(list("GRanges" = gr_feature_longest,
                       "GRangesList" = feature_longest, "Output" = outfile))
    } else {
        if (featureName %in% c("cds", "utr5", "utr3", "transcript")) {
            gr_feature <- rtracklayer::asBED(feature)
            ## convert each element of Grangeslist to one Grange with blocks
            ## info as metadata
            names(gr_feature) <- gr_feature$name
            seqinfo(gr_feature) <- seqInfo
            if (export) {
                outfile <- paste(featureSource, "_", featureName, ".bed12",
                                 sep = "")
                export.bed(gr_feature, outfile)
            }
        }
        if (featureName %in% c("gene", "intron", "exon")) {
            if (featureName == "gene") {
                gr_feature <- feature # gr
                feature <- as(split(feature, f = feature$gene_id),
                              "GRangesList") # grl
            } else {
                gr_feature <- unlist(feature, use.names = TRUE)
                ## convert Grangeslist object to GRanges object
                ## if the names are not unique, force them to be unique
                if (sum(duplicated(names(gr_feature))) > 0) {
                    names(gr_feature) <- paste(names(gr_feature),
                                        seq_along(names(gr_feature)), sep = "_")
                }
            }

            seqinfo(gr_feature) <- seqInfo
            if (export) {
                outfile <- paste(featureSource, "_", featureName, ".bed",
                                 sep = "")
                export.bed(gr_feature, outfile)
            }
        }
        message("[get_genomic_feature_coordinates] for ", featureName, " finished!\n")
        invisible(list("GRanges" = gr_feature, "GRangesList" = feature,
                       "Output" = outfile))
    }
}


#' @title Demarcate genes into promoter, gene body  and TTS features
#
#' @description This is a helper function for 'plot_3parts_metagene', used to
#' speed up plotting of multiple data sets with the same configuration. Use
#' featureName='transcript' and meta=FALSE and longest=TRUE for genes.
#'
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param featureName one of the gene feature in c("utr3", "utr5", "cds",
#' "transcript")
#' @param meta logical, indicating whether a metagene (intron excluded) or
#'  genomic (intron included) plot should be produced
#' @param nbins an integer defines the total number of bins
#' @param fiveP extension out of the 5' boundary of gene
#' @param threeP extension out of the 3' boundary of gene
#' @param verbose logical, whether to output additional information
#' @param longest logical, indicating whether the output should be limited to
#'  the longest transcript of each gene
#' @param protein_coding logical, indicating whether to limit to protein_coding
#'  genes
#'
#' @return a named list with the elements c("windowRs", "nbins", "scaled_bins",
#'  "fiveP", "threeP", "meta", "longest")
#' @author Shuye Pu
#'
#' @examples
#' gtfFile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#'     package = "GenomicPlot"
#' )
#'
#' txdb <- custom_TxDb_from_GTF(gtfFile, genome = "hg19")
#'
#' gf <- prepare_3parts_genomic_features(txdb,
#'     meta = FALSE, nbins = 100, fiveP = -1000, threeP = 1000,
#'     longest = FALSE
#' )
#'
#' @export prepare_3parts_genomic_features

prepare_3parts_genomic_features <- function(txdb,
                                            featureName = "transcript",
                                            meta = TRUE,
                                            nbins = 100,
                                            fiveP = -1000,
                                            threeP = 1000,
                                            longest = TRUE,
                                            protein_coding = TRUE,
                                            verbose = FALSE) {
    ## prepare transcripts
    stopifnot(featureName %in% c("utr3", "utr5", "cds", "transcript"))
    stopifnot(is.numeric(c(nbins, fiveP, threeP)))
    stopifnot(fiveP <= 0 && threeP >= 0)
    if (verbose)
        message("[prepare_3parts_genomic_features] started ...\n")

    ## prepare transcripts that are suitable for overlap
    if (featureName %in% c("utr5", "cds", "utr3")) {
        featureName <- ifelse(meta, featureName,
                              paste0(featureName, "(with intron)"))
    } else if (featureName == "transcript") {
        featureName <- ifelse(meta, featureName, "gene")
    }

    five <- fiveP / 1000
    five <- paste0(five, "Kb")
    if (fiveP == 0) five <- "-0Kb"
    three <- threeP / 1000
    three <- paste0(three, "Kb")
    if (threeP == 0) three <- "0Kb"

    featureNames <- c(five, featureName, three)

    gn <- get_genomic_feature_coordinates(txdb, featureName, longest = longest,
                                          protein_coding = protein_coding)

    if (meta) {
        feature <- gn$GRangesList
    } else {
        feature <- as(split(gn$GRanges, as.factor(names(gn$GRanges))),
                      "GRangesList")
    }

    wf <- vapply(as.list(width(feature)), sum, numeric(1))
    means <- c(promoter = -fiveP, median(wf), TTS = threeP)
    scaled_bins <- round(means * nbins / sum(means))
    sel_tx <- names(feature[wf > scaled_bins[2]]) # selected_tx

    promoter <- flank(gn$GRanges, width = -fiveP, both = FALSE, start = TRUE,
                      ignore.strand = FALSE)
    TTS <- flank(gn$GRanges, width = threeP, both = FALSE, start = FALSE,
                 ignore.strand = FALSE)

    promoter <- check_constraints(promoter, genome = txdb$user_genome[1])
    TTS <- check_constraints(TTS, genome = txdb$user_genome[1])

    sel_tx <- sel_tx[sel_tx %in% intersect(names(promoter), names(TTS))]
    windowRs <- list(
        as(split(promoter, as.factor(names(promoter))),
           "GRangesList")[sel_tx],
        feature[sel_tx],
        as(split(TTS, as.factor(names(TTS))), "GRangesList")[sel_tx]
    )

    names(windowRs) <- featureNames
    names(scaled_bins) <- featureNames

    if (verbose) {
        message("Median sizes for features: ", paste(means, collase = " "))
        message("Bin sizes for features: ", paste(scaled_bins, collapse = " "))
        message(
            "Number of transcripts: ",
            paste(vapply(windowRs, length, numeric(1)), collapse = " "), "\n"
        )
        message("[prepare_3parts_genomic_features] finished!\n")
    }

    invisible(list(
        "windowRs" = windowRs, "nbins" = nbins, "scaled_bins" = scaled_bins,
        "fiveP" = fiveP, "threeP" = threeP, "meta" = meta, "longest" = longest
    ))
}

#' @title Demarcate genes into promoter, 5'UTR, CDS, 3'UTR and TTS features
#
#' @description This is a helper function for 'plot_5parts_metagene', used to
#' speed up plotting of multiple data sets with the same configuration. Only
#' protein-coding genes are considered.
#'
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param meta logical, indicating whether a metagene (intron excluded) or gene
#'  (intron included) plot should be produced
#' @param nbins an integer defines the total number of bins
#' @param fiveP extension out of the 5' boundary of gene
#' @param threeP extension out of the 3' boundary of gene
#' @param verbose logical, whether to output additional information
#' @param longest logical, indicating whether the output should be limited to
#'  the longest transcript of each gene
#' @param subsetTx a vector of transcript names (eg. ENST00000587541.1) for
#'  subsetting the genome
#'
#' @return a named list with the elements c("windowRs", "nbins", "scaled_bins",
#'  "fiveP", "threeP", "meta", "longest")
#' @author Shuye Pu
#'
#' @examples
#' gtfFile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#'     package = "GenomicPlot"
#' )
#'
#' txdb <- custom_TxDb_from_GTF(gtfFile, genome = "hg19")
#'
#' gf <- prepare_5parts_genomic_features(txdb,
#'     meta = TRUE, nbins = 100, fiveP = -0, threeP = 0,
#'     longest = TRUE
#' )
#'
#' @export prepare_5parts_genomic_features
#'
prepare_5parts_genomic_features <- function(txdb,
                                            meta = TRUE,
                                            nbins = 100,
                                            fiveP = -1000,
                                            threeP = 1000,
                                            longest = TRUE,
                                            verbose = FALSE,
                                            subsetTx = NULL) {
    ## prepare transcripts
    stopifnot(is.numeric(c(nbins, fiveP, threeP)))
    stopifnot(fiveP <= 0 && threeP >= 0)
    if (verbose)
        message("[prepare_5parts_genomic_features] started ...\n")
    five <- fiveP / 1000
    five <- paste0(five, "Kb")
    if (fiveP == 0) five <- "-0Kb"
    three <- threeP / 1000
    three <- paste0(three, "Kb")
    if (threeP == 0) three <- "0Kb"
    featureNames <- c(five, "5'UTR", "CDS", "3'UTR", three)

    if (!meta) longest <- TRUE
    ## always use the longest transcript to represent the gene

    utr5 <- get_genomic_feature_coordinates(txdb, "utr5", longest = longest,
                                            protein_coding = TRUE)
    utr3 <- get_genomic_feature_coordinates(txdb, "utr3", longest = longest,
                                            protein_coding = TRUE)
    cds <- get_genomic_feature_coordinates(txdb, "cds", longest = longest,
                                           protein_coding = TRUE)
    gn <- get_genomic_feature_coordinates(txdb, "transcript", longest = longest,
                                          protein_coding = TRUE)

    promoter <- flank(gn$GRanges, width = -fiveP, both = FALSE, start = TRUE,
                      ignore.strand = FALSE)
    TTS <- flank(gn$GRanges, width = threeP, both = FALSE, start = FALSE,
                 ignore.strand = FALSE)

    promoter <- check_constraints(promoter, genome = txdb$user_genome[1])
    TTS <- check_constraints(TTS, genome = txdb$user_genome[1])

    if (meta) {
        utr5_grl <- utr5$GRangesList
        cds_grl <- cds$GRangesList
        utr3_grl <- utr3$GRangesList
    } else {
        utr5_grl <- as(split(utr5$GRanges, as.factor(names(utr5$GRanges))),
                       "GRangesList")
        cds_grl <- as(split(cds$GRanges, as.factor(names(cds$GRanges))),
                      "GRangesList")
        utr3_grl <- as(split(utr3$GRanges, as.factor(names(utr3$GRanges))),
                       "GRangesList")
    }

    grls <- list("5'UTR" = utr5_grl, "CDS" = cds_grl, "3'UTR" = utr3_grl)
    l3 <- vapply(grls, function(feature)
        median(vapply(as.list(width(feature)), sum, numeric(1))), numeric(1))

    means <- c(promoter = -fiveP, l3, TTS = threeP)
    names(means) <- featureNames
    scaled_bins <- round(means * nbins / sum(means))
    names(scaled_bins) <- featureNames

    selected_tx <- lapply(names(grls), function(x) {
        len <- vapply(as.list(width(grls[[x]])), sum, numeric(1))
        ## len is named vector, where names are the tx_ids
        y <- names(len)[which(len >= scaled_bins[x])]
        ## since 'check_constraints' may eliminate some promoter or TTS
        ## make sure the transcript is common to all features
        y <- y[y %in% intersect(names(promoter), names(TTS))]
    })

    names(selected_tx) <- names(grls)

    if (!is.null(subsetTx)) {
        selected_tx[["custom"]] <- subsetTx
    }
    selected_tx <- Reduce(intersect, selected_tx)

    windowRs <- list(
        as(split(promoter, as.factor(names(promoter))),
           "GRangesList")[selected_tx],
        utr5_grl[selected_tx],
        cds_grl[selected_tx],
        utr3_grl[selected_tx],
        as(split(TTS, as.factor(names(TTS))), "GRangesList")[selected_tx]
    )

    names(windowRs) <- featureNames

    if (verbose) {
        message("Median sizes for features: ", paste(means, collase = " "))
        message("Bin sizes for features: ", paste(scaled_bins, collapse = " "))
        message("Number of transcripts: ",
            paste(vapply(windowRs, length, numeric(1)), collapse = " "), "\n"
        )
        message("[prepare_5parts_genomic_features] finished!\n")
    }

    invisible(list("windowRs" = windowRs, "nbins" = nbins,
                   "scaled_bins" = scaled_bins, "fiveP" = fiveP,
                   "threeP" = threeP, "meta" = meta, "longest" = longest))
}

#' @title Get genomic coordinates of features of protein-coding genes
#
#' @description Get genomic coordinates of promoter, 5'UTR, CDS, 3'UTR, TTS and
#' intron for the longest transcript of protein-coding genes. The range of
#' promoter is defined by fiveP and dsTSS upstream and downstream TSS,
#' respectively, the TTS ranges from the 3' end of the gene to threeP
#' downstream, or the start of a downstream gene, whichever is closer.
#'
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param fiveP extension upstream of the 5' boundary of genes
#' @param threeP extension downstream of the 3' boundary of genes
#' @param dsTSS range of promoter extending downstream of TSS
#' @param nc number of cores for parallel processing
#'
#' @return a GRangesList object
#' @author Shuye Pu
#'
#' @examples
#' gtfFile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#'     package = "GenomicPlot"
#' )
#'
#' txdb <- custom_TxDb_from_GTF(gtfFile, genome = "hg19")
#'
#' f <- get_txdb_features(txdb, dsTSS = 100, fiveP = -100, threeP = 100)
#'
#' @export get_txdb_features
#'
get_txdb_features <- function(txdb,
                              fiveP = -1000,
                              dsTSS = 300,
                              threeP = 1000,
                              nc = 2) {
    stopifnot(is.numeric(c(dsTSS, fiveP, threeP)))
    stopifnot(fiveP <= 0 && threeP >= 0)
    message("[get_txdb_features]")
    utr5 <- get_genomic_feature_coordinates(txdb, "utr5", longest = TRUE,
                                            protein_coding = TRUE)
    utr3 <- get_genomic_feature_coordinates(txdb, "utr3", longest = TRUE,
                                            protein_coding = TRUE)
    cds <- get_genomic_feature_coordinates(txdb, "cds", longest = TRUE,
                                           protein_coding = TRUE)
    intron <- get_genomic_feature_coordinates(txdb, "intron", longest = TRUE,
                                              protein_coding = TRUE)

    gene_gr <- get_genomic_feature_coordinates(txdb, "transcript",
                                               longest = TRUE,
                                               protein_coding = TRUE)$GRanges

    features <- GRangesList(
        "5'UTR" = unlist(utr5$GRangesList, use.names = TRUE),
        "3'UTR" = unlist(utr3$GRangesList, use.names = TRUE),
        "CDS" = unlist(cds$GRangesList, use.names = TRUE),
        "Intron" = unlist(intron$GRangesList, use.names = TRUE),
        compress = FALSE
    )
    if (fiveP <= 0) {
        promoter <- GenomicRanges::promoters(gene_gr, upstream = -fiveP,
                                             downstream = dsTSS,
                                             use.names = TRUE)
        promoter <- check_constraints(promoter, genome = txdb$user_genome[1])

        features[["Promoter"]] <- promoter
    } else {
        promoter <- NULL
    }

    if (threeP >= 0) {
        TTS <- GenomicRanges::flank(gene_gr, width = threeP, both = FALSE,
                                    start = FALSE, ignore.strand = FALSE)
        TTS <- check_constraints(TTS, genome = txdb$user_genome[1])

        ol <- GenomicRanges::findOverlaps(TTS, gene_gr)
        queries <- unique(ol@from)

        if (length(queries > 0)) {
            cl <- start_parallel(nc = nc)

            parallel::clusterExport(cl, varlist = c("TTS", "gene_gr", "ol"),
                                    envir = environment())
            parallel::clusterExport(cl, varlist = c("nearest", "runValue",
                                                    "strand", "start", "end",
                                                    "start<-", "end<-"),
                                    envir = environment())

            out <- parallel::parLapply(cl, queries, function(tts_idx) {
                gene_idx <- ol@to[ol@from == tts_idx]
                ogenes <- gene_gr[gene_idx]
                ogene <- ogenes[nearest(TTS[tts_idx], ogenes,
                                        select = "arbitrary")]

                if (runValue(strand(TTS[tts_idx])) == "+") {
                    dist <- start(ogene) - start(TTS[tts_idx])
                    if (dist > 0) {
                        end(TTS[tts_idx]) <- start(TTS[tts_idx]) + dist
                    } else {
                        end(TTS[tts_idx]) <- start(TTS[tts_idx])
                    }
                } else if (runValue(strand(TTS[tts_idx])) == "-") {
                    dist <- end(TTS[tts_idx]) - end(ogene)
                    if (dist > 0) {
                        start(TTS[tts_idx]) <- end(TTS[tts_idx]) - dist
                    } else {
                        start(TTS[tts_idx]) <- end(TTS[tts_idx])
                    }
                } else {
                    stop("invalid strand")
                }
                return(c(tts_idx, start(TTS[tts_idx]), end(TTS[tts_idx])))
            })
            stop_parallel(cl)

            out1 <- as.data.frame(Reduce(rbind, out))
            start(TTS[out1[, 1]]) <- out1[, 2]
            end(TTS[out1[, 1]]) <- out1[, 3]
        }

        features[["TTS"]] <- TTS
    } else {
        TTS <- NULL
    }

    for (aname in names(features)) {
        f <- features[[aname]]
        mcols(f) <- NULL
        f$tx_name <- names(f)
        features[[aname]] <- f
    }

    invisible(features)
}

#' @title Get the number of peaks overlapping each feature of all protein-coding
#' genes
#
#' @description Annotate each peak with genomic features based on overlap, and
#' produce summary statistics for distribution of peaks in features of
#' protein-coding genes. If a peak overlap multiple features, a feature is
#' assigned to the peak in the following order of precedence: "5'UTR",
#' "3'UTR", "CDS", "Intron", "Promoter", "TTS".
#'
#' @param peak a GRanges object defining query ranges
#' @param features a GRangesList object representing genomic features
#' @param stranded logical, indicating whether the overlap should be
#'  strand-specific
#'
#' @return a list object
#' @author Shuye Pu
#'
#' @examples
#'
#' gtfFile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#'     package = "GenomicPlot"
#' )
#'
#' txdb <- custom_TxDb_from_GTF(gtfFile, genome = "hg19")
#' f <- get_txdb_features(txdb, dsTSS = 100, fiveP = 0, threeP = 1000)
#'
#' p <- RCAS::importBed(system.file("extdata", "test_chip_peak_chr19.bed",
#'     package = "GenomicPlot"
#' ))
#' ann <- get_targeted_genes(peak = p, features = f, stranded = FALSE)
#'
#' @note used in \code{plot_peak_annotation}
#' @export get_targeted_genes

get_targeted_genes <- function(peak,
                               features,
                               stranded = TRUE) {
    precedence <- function(terms) {
        scope <- c("5'UTR", "3'UTR", "CDS", "Intron", "Promoter", "TTS")
        rank <- seq(6)
        names(rank) <- scope
        pre <- names(rank[rank == min(rank[terms])])
        return(pre)
    }

    message("[get_targeted_genes]")
    num_peaks <- length(peak)
    num_genes <- length(unique(features$CDS$tx_name))

    annot_list <- lapply(names(features), function(feature) {
        featureGr <- features[[feature]]
        featureOverlaps <- findOverlaps(peak, featureGr,
                                        ignore.strand = !stranded)
        peak_df <- gr2df(peak[featureOverlaps@from])
        if (is.null(peak_df$strand)) peak_df$strand <- rep("*", nrow(peak_df))
        peak_df <- peak_df %>%
            mutate(
                chrPeak = as.character(chr),
                startPeak = as.integer(start) - 1,
                endPeak = as.integer(end),
                idPeak = as.character(name),
                scorePeak = as.integer(score),
                strandPeak = as.character(strand),
                .keep = "unused"
            )
        feature_df <- featureGr[featureOverlaps@to]
        names(feature_df) <- seq_along(feature_df)
        feature_df <- gr2df(feature_df) %>%
            mutate(
                chrfeature = as.character(chr),
                startfeature = as.integer(start) - 1,
                endfeature = as.integer(end),
                widthfeature = as.integer(score),
                strandfeature = as.character(strand),
                .keep = "unused"
            )

        ot <- cbind(peak_df, feature_df) %>%
            select(chrPeak, startPeak, endPeak, idPeak, scorePeak, strandPeak,
                   chrfeature, startfeature, endfeature, widthfeature,
                   strandfeature, tx_name) %>%
            mutate(feature_name = feature)
        return(ot)
    })

    annot_table <- bind_rows(annot_list)
    annot_table <- annot_table %>%
        group_by(chrPeak, startPeak, endPeak, strandPeak) %>%
        filter(n() == 1 | feature_name == precedence(unique(feature_name)))
    ## if the peak is assigned to only one feature, associate that feature with
    ## the peak, else if the peak is assigned to multiple features, associate
    ## the feature with the best precedence order with the peak.

    annot_count <- as.data.frame(annot_table %>%
        group_by(tx_name) %>%
        count(feature_name))

    overlap_genes <- length(unique(annot_table$tx_name))
    overlap_peaks <- length(unique(annot_table$idPeak))


    annot_df <- data.frame(tx_name = unique(features$CDS$tx_name),
                           Promoter = 0, `5'UTR` = 0, CDS = 0, `3'UTR` = 0,
                           TTS = 0, Intron = 0)
    colnames(annot_df) <- c("tx_name", "Promoter", "5'UTR", "CDS", "3'UTR",
                            "TTS", "Intron")
    rownames(annot_df) <- annot_df$tx_name


    for (i in seq_len(nrow(annot_count))) {
        txi <- annot_count[i, 1]
        fn <- annot_count[i, 2]
        count <- annot_count[i, 3]
        annot_df[txi, fn] <- count
    }

    annot_df <- annot_df %>%
        mutate(Exon = `5'UTR` + CDS + `3'UTR`) %>%
        mutate(Transcript = Intron + Exon)

    feature_counts <- apply(annot_df[, c("Promoter", "5'UTR", "CDS", "3'UTR",
                                        "TTS", "Intron", "Exon", "Transcript")],
                            2, sum)

    invisible(list(gene_table = annot_df, peak_table = annot_table,
                   num_peak = num_peaks, num_gene = num_genes,
                   feature_count = feature_counts, overlap_peak = overlap_peaks,
                   overlap_gene = overlap_genes))
}

#' @title Make TxDb object from a GTF file for a subset of genes
#
#' @description Make a partial TxDb object given a GTF file and a list of gene
#' names in a file or in a character vector.
#'
#' @param gtfFile path to a GTF file
#' @param genome version of genome, like "hg19"
#' @param geneList path to a tab-delimited text file with one gene name on each
#'  line, or a character vector of gene names
#' @param geneCol the position of the column that containing gene names in the
#'  case that geneList is a file
#'
#' @return a TxDb object
#' @author Shuye Pu
#'
#' @examples
#'
#' gtfFile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#'     package = "GenomicPlot"
#' )
#' genes <- c("RPRD1A", "RPAP2", "RPRD1B", "RPRD2", "ZNF281", "YTHDF2")
#'
#' txdb <- make_subTxDb_from_GTF(gtfFile = gtfFile, geneList = genes)
#'
#' @export make_subTxDb_from_GTF

make_subTxDb_from_GTF <- function(gtfFile,
                                  genome = "hg19",
                                  geneList,
                                  geneCol = 1) {
    message("[make_subTxDb_from_GTF] ", gtfFile)
    gff <- RCAS::importGtf(saveObjectAsRds = TRUE, filePath = gtfFile)
    if (length(geneList) == 1) {
        if (file.exists(geneList)) {
            aList <- read.delim2(geneList, comment.char = "#")
            geneList <- as.character(aList[, geneCol])
        } else {
            stop("Gene list file does not exist!")
        }
    }

    subgff <- gff[gff$gene_name %in% geneList]
    maploss <- length(geneList) - length(subgff)
    message("In make_subTxDb_from_GTF, number of gene symbols failed to map: ",
            maploss, "\n")

    txdb <- txdbmaker::makeTxDbFromGRanges(subgff,
                                   metadata = data.frame(name = "genome",
                                                         value = genome))

    return(txdb)
}

#' @title Translate gene names to transcript ids using a GTF file for a subset
#' of genes
#
#' @description Given a list of gene names in a file or in a character vector,
#' turn them into a vector of transcript ids.
#'
#' @param gtfFile path to a GTF file
#' @param geneList path to a tab-delimited text file with one gene name on each
#'  line, or a character vector of gene names (eg. RPRD1B)
#' @param geneCol the position of the column that containing gene names in the
#'  case that geneList is a file
#'
#' @return a vector of transcript ids (eg. ENST00000577222.1)
#' @author Shuye Pu
#' @examples
#'
#' gtfFile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#'     package = "GenomicPlot"
#' )
#' genes <- c("RPRD1A", "RPAP2", "RPRD1B", "RPRD2", "ZNF281", "YTHDF2")
#'
#' tx <- gene2tx(gtfFile = gtfFile, geneList = genes)
#'
#' @export gene2tx

gene2tx <- function(gtfFile,
                    geneList,
                    geneCol = 1) {
    message("[gene2tx]")
    gff <- RCAS::importGtf(saveObjectAsRds = TRUE,
                           readFromRds = FALSE, filePath = gtfFile)

    if (length(geneList) == 1) {
        if (file.exists(geneList)) {
            aList <- read.delim2(geneList, comment.char = "#")
            geneList <- as.character(aList[, geneCol])
        } else {
            stop("Gene list file does not exist!")
        }
    }
    subgff <- gff[gff$gene_name %in% geneList]

    tx <- as.character(unique(subgff$transcript_id))

    return(tx)
}

#' @title Check constraints of genomic ranges
#' @description Make sure the coordinates of GRanges are within the boundaries
#' of chromosomes, and trim anything that goes beyond. Also, remove entries
#' whose seqname is not in the seqname of a query GRanges.
#'
#' @param gr a GenomicRanges object
#' @param genome genomic version name such as "hg19"
#' @param queryRle a RleList object used as a query against gr
#'
#' @return a GRanges object
#' @author Shuye Pu
#'
#' @examples
#' subject <- GRanges("chr19",
#'     IRanges(rep(c(10, 15), 2), width = c(1, 20, 400, 2e+8)),
#'     strand = c("+", "+", "-", "-")
#' )
#'
#' g <- check_constraints(gr = subject, genome = "hg19")
#' identical(g, subject)
#'
#' subject1 <- GRanges("chr19",
#'     IRanges(rep(c(10, 15), 2), width = c(1, 20, 400, 28)),
#'     strand = c("+", "+", "-", "-")
#' )
#'
#' g1 <- check_constraints(gr = subject1, genome = "hg19")
#' identical(g1, subject1)
#'
#' @export check_constraints
#'
check_constraints <- function(gr, genome, queryRle = NULL) {
    stopifnot(is.character(genome))

    message("[check_constraints]")
    seqInfo <- set_seqinfo(genome)
    len <- seqlengths(seqInfo)

    # limit gr to chromosomes in chromInfo
    gr <- gr[as.vector(seqnames(gr)) %in% seqlevels(seqInfo)]

    if (any(start(gr) < 1)) {
        start(gr)[start(gr) < 1] <- 1
    }
    too_long <- end(gr) > len[as.vector(seqnames(gr))]

    if (any(too_long)) {
        end(gr)[too_long] <- len[as.vector(seqnames(gr))][too_long]
    }

    if (!is.null(queryRle)) {
        gr <- gr[as.vector(seqnames(gr)) %in% names(queryRle)]
    }

    return(gr)
}

#' @title Filter GRanges by overlaps in a stranded way
#' @description This function reports all query GRanges that have overlaps in
#' subject GRanges. Strand information is used to define overlap.
#' @param query a GRanges object
#' @param subject a GRanges object
#' @param maxgap an integer denoting the distance that define overlap
#' @param minoverlap The minimum amount of overlap between intervals as a
#' single integer greater than 0. If you modify this argument, maxgap must be
#' held fixed.
#' @param ignore.order logical, indicating whether the order of query and
#' subject can be switched, default = TRUE. Overlaps in query and subject often
#' have different sizes. This parameter will make the function use whichever is
#' smaller to avoid errors when making Venn diagrams.
#'
#' @return a GRanges object
#' @author Shuye Pu
#'
#' @examples
#'
#' query <- GRanges("chr19",
#'     IRanges(rep(c(10, 15), 2), width = c(1, 20, 40, 50)),
#'     strand = c("+", "+", "-", "-")
#' )
#'
#' subject <- GRanges("chr19",
#'     IRanges(rep(c(13, 150), 2), width = c(10, 14, 20, 28)),
#'     strand = c("+", "-", "-", "+")
#' )
#'
#' res <- filter_by_overlaps_stranded(query, subject)
#' res
#' resf <- filter_by_overlaps_stranded(query, subject, ignore.order = FALSE)
#' resf
#'
#' @export filter_by_overlaps_stranded

filter_by_overlaps_stranded <- function(query, subject, maxgap = -1L,
                                        minoverlap = 0L, ignore.order = TRUE)
    {
    message("[filter_by_overlaps_stranded]")
    plus_query <- query[strand(query) == "+"]
    minus_query <- query[strand(query) == "-"]
    plus_subject <- subject[strand(subject) == "+"]
    minus_subject <- subject[strand(subject) == "-"]

    overlap_plus <- filter_by_overlaps(plus_query, plus_subject,
                                    maxgap = maxgap, minoverlap = minoverlap)
    overlap_minus <- filter_by_overlaps(minus_query, minus_subject,
                                    maxgap = maxgap, minoverlap = minoverlap)

    overlaps <- c(overlap_plus, overlap_minus)

    if (ignore.order) {
        overlap_plus_r <- filter_by_overlaps(plus_subject, plus_query,
                                           maxgap = maxgap,
                                           minoverlap = minoverlap)
        overlap_minus_r <- filter_by_overlaps(minus_subject, minus_query,
                                            maxgap = maxgap,
                                            minoverlap = minoverlap)

        overlaps_r <- c(overlap_plus_r, overlap_minus_r)

        if(length(overlaps_r) < length(overlaps)){
            overlaps <- overlaps_r
            message("The overlaps is from the Subject instead of the Query!\n")
        }

    }

    invisible(overlaps)
}

#' @title Filter GRanges by nonoverlaps in a stranded way
#' @description This function reports all query GRanges that do not overlaps
#' GRanges in subject. Strand information is used to define overlap.
#' @param query a GRanges object
#' @param subject a GRanges object
#' @param maxgap an integer denoting the distance that define overlap
#' @param minoverlap The minimum amount of overlap between intervals as a
#' single integer greater than 0. If you modify this argument, maxgap must be
#' held fixed.
#' @param ignore.order logical, indicating whether the order of query and
#' subject can be switched, default = TRUE. This parameter is used to avoid
#' the situation that the size of overlaps is bigger than the size of subject,
#' which will produce an error when plotting Venn diagrams.
#'
#' @return a GRanges object
#' @author Shuye Pu
#' @examples
#'
#' query <- GRanges("chr19",
#'     IRanges(rep(c(10, 15), 2), width = c(1, 20, 40, 50)),
#'     strand = c("+", "+", "-", "-")
#' )
#'
#' subject <- GRanges("chr19",
#'     IRanges(rep(c(13, 150), 2), width = c(10, 14, 20, 28)),
#'     strand = c("+", "-", "-", "+")
#' )
#'
#' res <- filter_by_nonoverlaps_stranded(query, subject)
#' res
#'
#' @export filter_by_nonoverlaps_stranded
#'
filter_by_nonoverlaps_stranded <- function(query, subject,  maxgap = -1L,
                                           minoverlap = 0L,
                                           ignore.order = TRUE) {
    message("[filter_by_nonoverlaps_stranded]")
    overlaps <- filter_by_overlaps_stranded(query, subject, maxgap = maxgap,
                                            minoverlap = minoverlap,
                                            ignore.order = ignore.order)

    nonoverlaps <- GenomicRanges::setdiff(query, overlaps)
    invisible(nonoverlaps)
}

#' @title Filter GRanges by overlaps in a nonstranded way
#' @description This function reports all query GRanges that have overlaps in
#' subject GRanges. Strand information is not required.
#' @param query a GRanges object
#' @param subject a GRanges object
#' @param maxgap an integer denoting the distance that define overlap
#' @param minoverlap The minimum amount of overlap between intervals as a
#' single integer greater than 0. If you modify this argument, maxgap must be
#' held fixed.
#' @param ignore.order logical, indicating whether the order of query and
#' subject can be switched, default = TRUE. This parameter is used to avoid
#' the situation that the size of overlaps is bigger than the size of subject,
#' which will produce an error when plotting Venn diagrams.
#'
#' @return a GRanges object
#' @author Shuye Pu
#'
#' @examples
#'
#' query <- GRanges("chr19",
#'     IRanges(rep(c(10, 15), 2), width = c(1, 20, 40, 50)),
#'     strand = c("+", "+", "-", "-")
#' )
#'
#' subject <- GRanges("chr19",
#'     IRanges(rep(c(13, 150), 2), width = c(10, 14, 20, 28)),
#'     strand = c("+", "-", "-", "+")
#' )
#'
#' res <- filter_by_overlaps_nonstranded(query, subject, ignore.order = TRUE)
#' res
#'
#' @export filter_by_overlaps_nonstranded

filter_by_overlaps_nonstranded <- function(query, subject, maxgap = -1L,
                                        minoverlap = 0L, ignore.order = TRUE)
{
    message("[filter_by_overlaps_nonstranded]\n")
    overlaps <- filter_by_overlaps(query, subject, maxgap = maxgap,
                                   minoverlap = minoverlap)

    if (ignore.order) {
        overlaps_r <- filter_by_overlaps(subject, query, maxgap = maxgap,
                                           minoverlap = minoverlap)
        if(length(overlaps_r) < length(overlaps)){
            overlaps <- overlaps_r
            message("The overlaps is from the Subject instead of the Query!\n")
        }
    }

    invisible(overlaps)
}

