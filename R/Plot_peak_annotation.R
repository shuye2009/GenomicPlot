#' @title Annotate peaks with genomic features and genes
#'
#' @description  Produce a table of transcripts targeted by peaks, and generate
#' plots for target gene types, and peak distribution in genomic features
#'
#' @param peakFile a string denoting the peak file name, only .bed format is
#'  allowed
#' @param gtfFile path to a gene annotation gtf file with gene_biotype field
#' @param importParams a list of parameters for \code{handle_input}
#' @param fiveP extension out of the 5' boundary of genes for defining promoter:
#'  fiveP TSS + dsTSS
#' @param dsTSS extension downstream of TSS for defining promoter: fiveP TSS +
#'  dsTSS
#' @param threeP extension out of the 3' boundary of genes for defining
#'  termination region: -0 TTS + threeP
#' @param outPrefix a string denoting output file name in pdf format
#' @param simple logical, indicating whether 5'UTR, CDS and 3'UTR are annotated
#'  in the gtfFile
#' @param verbose, logical, to indicate whether to write the annotation results
#'  to a file
#' @param hw a vector of two elements specifying the height and width of the
#'  output figures
#' @param nc number of cores for parallel processing
#'
#' @return a list of three dataframes, 'annotation' is the annotation of peaks
#'  into gene types, 'stat' is the summary stats for pie chart, 'simplified' is
#'  the summary stats excluding intron
#'
#' @author Shuye Pu
#' @examples
#'
#' gtfFile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#'     package = "GenomicPlot"
#' )
#'
#' centerFile <- system.file("extdata", "test_chip_peak_chr19.bed",
#'     package = "GenomicPlot"
#' )
#' names(centerFile) <- c("summitPeak")
#'
#' bedimportParams <- setImportParams(
#'     offset = 0, fix_width = 100, fix_point = "center", norm = FALSE,
#'     useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19"
#' )
#'
#' plot_peak_annotation(
#'     peakFile = centerFile, gtfFile = gtfFile, importParams = bedimportParams,
#'     fiveP = -2000, dsTSS = 200, threeP = 2000, simple = FALSE
#' )
#'
#' @export plot_peak_annotation
#'
plot_peak_annotation <- function(peakFile,
                                 gtfFile,
                                 importParams = NULL,
                                 fiveP = -1000,
                                 dsTSS = 300,
                                 threeP = 1000,
                                 simple = FALSE,
                                 outPrefix = NULL,
                                 verbose = FALSE,
                                 hw = c(8, 8),
                                 nc = 2) {
    stopifnot(is.numeric(c(fiveP, dsTSS, threeP, nc, hw)))
    stopifnot(all(file.exists(peakFile)))
    if (is.null(names(peakFile)) || any(names(peakFile) == ""))
        stop("Each file must have a name attribute!")

    if (verbose) message("[plot_peak_annotation] started ...\n")
    functionName <- as.character(match.call()[[1]])
    params <- plot_named_list(as.list(environment()))
    force(params)

    peakLabel <- names(peakFile)
    importParams$useScore <- FALSE
    importParams$outRle <- FALSE
    bedin <- handle_input(inputFiles = peakFile, importParams = importParams)
    peak <- bedin[[peakLabel]]$query

    stranded <- TRUE
    if (any(runValue(strand(peak)) %in% c("*", ".", " ", ""))) {
        stranded <- FALSE
    }

    if (!is.null(outPrefix)) {
        while (!is.null(dev.list())) {
            dev.off()
        }
        pdf(paste0(outPrefix, ".pdf"), height = hw[1], width = hw[2])
    }

    chromInfo <- set_seqinfo(importParams$genome)
    gff <- RCAS::importGtf(saveObjectAsRds = TRUE, filePath = gtfFile)
    gff <- gff[as.vector(seqnames(gff)) %in% seqlevels(chromInfo)]
    seqlevels(gff) <- seqlevels(chromInfo)
    seqinfo(gff) <- chromInfo
    txdb <- txdbmaker::makeTxDbFromGRanges(gff)

    if (simple) {
        ## the 5'UTR and 3'UTR are not annotted
        exon <- get_genomic_feature_coordinates(txdb, "exon", longest = FALSE)
        intron <- get_genomic_feature_coordinates(txdb, "intron",
                                                  longest = FALSE)
        gene <- get_genomic_feature_coordinates(txdb, "gene",
                                                longest = FALSE)$GRanges

        promoter <- promoters(gene, upstream = -fiveP, downstream = 300,
                              use.names = FALSE)
        names(promoter) <- NULL
        TTS <- flank(gene, width = threeP, both = FALSE, start = FALSE,
                     ignore.strand = FALSE)
        names(TTS) <- NULL

        pg <- mergeByOverlaps(peak, gene)
        targeted_gene <- cbind(gr2df(pg$peak), gr2df(pg$gene))

        pp <- mergeByOverlaps(peak, promoter)
        targeted_promoter <- cbind(gr2df(pp$peak), gr2df(pp$promoter))


        if (verbose) {
            write.table (targeted_gene,
                paste(peakLabel, "_targeted_annotated_gene.tab", sep = ""),
                sep = "\t", row.names = FALSE, quote = FALSE)
            write.table(targeted_promoter,
                        paste(peakLabel, "_targeted_promoter.tab", sep = ""),
                        sep = "\t", row.names = FALSE, quote = FALSE)
        }

        features <- GRangesList(
            "Promoter" = promoter,
            "TTS" = TTS,
            "Exon" = unlist(exon$GRangesList, use.names = FALSE),
            "Intron" = unlist(intron$GRangesList, use.names = FALSE),
            compress = FALSE
        )

        featureNames <- c("Promoter", "Exon", "TTS", "Intron")
        annot <- annotateWithFeatures(peak, features, strand.aware = TRUE,
                                      intersect.chr = FALSE)
        if (verbose) print(annot)
        precedence_count <- annot@num.precedence
        lengths <- vapply(features, FUN = function(x) sum(width(x)),
                          FUN.VALUE = numeric(1))

        df <- data.frame(count = precedence_count, len = lengths) %>%
            mutate(percent = count / sum(count)) %>%
            mutate(feature = rownames(.)) %>%
            mutate(labels = percent(percent, accuracy = 0.1)) %>%
            mutate(percent = round(percent * 100, digits = 1)) %>%
            mutate(norm_count = round(count * 1000 / len, digits = 3)) %>%
            mutate(norm_percent = norm_count / sum(norm_count)) %>%
            mutate(norm_labels = percent(norm_percent, accuracy = 0.1)) %>%
            mutate(norm_percent = round(norm_percent * 100, digits = 1)) %>%
            mutate(feature = factor(feature, levels = featureNames))

        df2 <- df %>%
            mutate(
                csum = rev(cumsum(rev(percent))),
                pos = percent / 2 + lead(csum, 1),
                pos = ifelse(is.na(pos), percent / 2, pos),
                norm_csum = rev(cumsum(rev(norm_percent))),
                norm_pos = norm_percent / 2 + lead(norm_csum, 1),
                norm_pos = ifelse(is.na(norm_pos), norm_percent / 2, norm_pos)
            )
        if (verbose) print(df)

        ap1 <- ggplot(df, aes(x = "", y = percent, fill = feature)) +
            geom_col(color = "white") +
            scale_y_continuous(breaks = df2$pos, labels = df2$labels) +
            guides(fill = guide_legend(nrow = 2)) +
            scale_fill_manual(values = viridis(n = nrow(df))) +
            coord_polar(theta = "y") +
            ggtitle("Absolute count") +
            theme(
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                axis.text = element_text(size = 12),
                legend.position = "top",
                legend.title = element_blank(),
                panel.background = element_rect(fill = "white")
            )

        ap2 <- ggplot(df, aes(x = "", y = norm_percent, fill = feature)) +
            geom_col(color = "white") +
            scale_y_continuous(breaks = df2$norm_pos,
                               labels = df2$norm_labels) +
            guides(fill = guide_legend(nrow = 2)) +
            scale_fill_manual(values = viridis(n = nrow(df))) +
            coord_polar(theta = "y") +
            ggtitle("Length-normalized count") +
            theme(
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                axis.text = element_text(size = 12),
                legend.position = "top",
                legend.title = element_blank(),
                panel.background = element_rect(fill = "white")
            )

        print(plot_grid(ap1, ap2))
        if (!is.null(outPrefix)) {
            print(params)
            on.exit(dev.off(), add = TRUE)
            if (verbose) sink()
        }
        if (verbose) message("[plot_peak_annotation] finished!\n")
        invisible(list(annotation = NULL, stat = df, simplified = NULL))
    } else {
        if (verbose) message("Collecting gene info...\n")

        overlaps <- gr2df(RCAS::queryGff(queryRegions = peak, gffData = gff))

        geneType <- NULL
        if ("gene_type" %in% colnames(overlaps)) {
            geneType <- "gene_type"
        } else if ("gene_biotype" %in% colnames(overlaps)) {
            geneType <- "gene_biotype"
        }
        if (is.null(geneType)) {
            stop("The annotation gtf file must have a field
                 'gene_type' or 'gene_biotype'")
        }
        gene_info_table <- gr2df(gff) %>%
            filter(type == "transcript") %>%
            select(chr, start, end, name, strand, transcript_id, gene_id,
                   gene_name, all_of(geneType))

        ## To find out the distribution of the query regions across gene types:
        if (verbose) message("Computing barchart of gene types...\n")

        df <- overlaps %>%
            filter(type == "gene") %>%
            group_by(.data[[geneType]]) %>%
            summarize(count = dplyr::n_distinct(queryIndex)) %>%
            rename_with(~ c("feature", "count"))

        intergenic_count <- length(peak) - sum(df$count)
        if (intergenic_count < 0) intergenic_count <- 0
        intergenic <- data.frame("feature" = "intergenic",
                                 "count" = intergenic_count)

        selected_feature <- c("protein_coding", "lincRNA", "antisense",
                              "pseudogene", "snRNA", "snoRNA", "rRNA")
        selected_count <- rep(0, length(selected_feature))
        selected <- data.frame("feature" = selected_feature,
                               "count" = selected_count)

        for (ft in selected_feature) {
            if (ft %in% df$feature) {
                selected[selected$feature == ft, "count"] <-
                    df[df$feature == ft, "count"]
            }
        }

        other_df <- filter(df, !feature %in% selected_feature)
        other <- data.frame("feature" = "other", "count" = sum(other_df$count))

        dfs <- as.data.frame(rbind(selected, other, intergenic)) %>%
            mutate(percent = round(count * 100 / sum(count), 1)) %>%
            arrange(desc(count))
        rownames(dfs) <- as.character(dfs$feature)

        pbar <- ggplot(dfs, aes(x = reorder(feature, -percent), y = percent)) +
            geom_bar(stat = "identity", aes(fill = feature)) +
            scale_fill_viridis(discrete = TRUE) +
            geom_text(aes(y = percent + 2), label = dfs$count) +
            labs(
                x = "", y = "Percent overlap (%)",
                title = "Annotation of peaks to all type of genes",
                subtitle = paste0("(Total number of peaks = ",
                                  length(peak), ")")
            ) +
            theme_bw(base_size = 14) +
            theme(
                axis.text.x = element_text(angle = 90),
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5)
            )


        ## get targeted genes table
        ## utr5, utr3 and cds here do not contain introns
        features <- get_txdb_features(
            txdb, fiveP = fiveP, dsTSS = dsTSS, threeP = threeP, nc = nc)
        dt <- get_targeted_genes(peak, features, stranded = stranded)

        if (verbose) write.table(
            dt$peak_table, paste(peakLabel, "_peak_annotations.tab", sep = ""),
            sep = "\t", row.names = FALSE, quote = FALSE)

        dt_gene <- left_join(dt$gene_table, gene_info_table,
                             by = c("tx_name" = "transcript_id"))

        dim(dt_gene)
        head(dt_gene)
        tail(dt_gene)

        ## filter based on peak type, for ChIPseq peak, only output genes
        ## targeted in promoters, for CLIPseq peaks only output genes targeted
        ##  in transcripts(5'UTR, CDS, 3'UTR, intron)
        if (fiveP < 0) {
            dt_gene <- dt_gene %>%
                arrange(desc(Promoter))
        } else {
            dt_gene <- dt_gene %>%
                arrange(desc(Transcript))
        }

        if (verbose) write.table(
            dt_gene, paste(peakLabel, "_targeted_annotated_gene.tab", sep = ""),
            sep = "\t", row.names = FALSE, quote = FALSE)

        # Plotting overlap counts between query regions and transcript features
        # here utr5, utr3 and cds do not contain introns

        if (verbose) message("Computing annotation stats...\n")
        lengths <- vapply(features, FUN = function(x) sum(width(x)),
                          FUN.VALUE = numeric(1))
        feature_counts <- dt$feature_count
        feature_counts <- feature_counts[!names(feature_counts) %in% c("Exon",
                                                                "Transcript")]

        precedence_c <- feature_counts[feature_counts > 0]

        lengths <- lengths[names(precedence_c)]

        if (verbose) message("Plotting piechart...\n")

        featureNames <- c("Promoter", "5'UTR", "CDS", "3'UTR", "TTS", "Intron")
        dfa <- data.frame(count = precedence_c, len = lengths) %>%
            mutate(feature = rownames(.)) %>%
            mutate(percent = count / sum(count)) %>%
            mutate(labels = scales::percent(percent, accuracy = 0.1)) %>%
            mutate(percent = round(percent * 100, digits = 1)) %>%
            mutate(norm_count = round(count * 1000 / len, digits = 3)) %>%
            mutate(norm_percent = norm_count / sum(norm_count)) %>%
            mutate(norm_labels = scales::percent(norm_percent,
                                                 accuracy = 0.1)) %>%
            mutate(norm_percent = round(norm_percent * 100, digits = 1)) %>%
            mutate(feature = factor(feature, levels = featureNames))

        dfa2 <- dfa %>%
            mutate(
                csum = rev(cumsum(rev(percent))),
                pos = percent / 2 + lead(csum, 1),
                pos = ifelse(is.na(pos), percent / 2, pos),
                norm_csum = rev(cumsum(rev(norm_percent))),
                norm_pos = norm_percent / 2 + lead(norm_csum, 1),
                norm_pos = ifelse(is.na(norm_pos), norm_percent / 2, norm_pos)
            )
        if (verbose) print(dfa2)

        fil <- viridis(n = nrow(dfa))
        apa1 <- ggplot(dfa, aes(x = "", y = percent, fill = feature)) +
            geom_col(color = "white" ) +
            coord_polar(theta = "y") +
            scale_y_continuous(breaks = dfa2$pos, labels = dfa2$labels) +
            guides(fill = guide_legend(nrow = 2)) +
            scale_fill_manual(values = fil) +
            theme(
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                axis.text = element_text(size = 12),
                legend.position = "top",
                legend.title = element_blank(),
                panel.background = element_rect(fill = "white")
            )
        apa2 <- ggplot(dfa, aes(x = "", y = norm_percent, fill = feature)) +
            geom_col(color = "white") +
            coord_polar(theta = "y") +
            scale_y_continuous(breaks = dfa2$norm_pos,
                               labels = dfa2$norm_labels) +
            guides(fill = guide_legend(nrow = 2)) +
            scale_fill_manual(values = fil) +
            theme(
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                axis.text = element_text(size = 12),
                legend.position = "top",
                legend.title = element_blank(),
                panel.background = element_rect(fill = "white")
            )

        dfb <- data.frame(count = precedence_c, len = lengths) %>%
            mutate(feature = rownames(.)) %>%
            filter(!feature %in% c("Intron", "NoFeature")) %>%
            mutate(percent = count / sum(count)) %>%
            mutate(labels = scales::percent(percent, accuracy = 0.1)) %>%
            mutate(percent = round(percent * 100, digits = 1)) %>%
            mutate(norm_count = round(count * 1000 / len, digits = 3)) %>%
            mutate(norm_percent = norm_count / sum(norm_count)) %>%
            mutate(norm_labels = scales::percent(norm_percent,
                                                 accuracy = 0.1)) %>%
            mutate(norm_percent = round(norm_percent * 100, digits = 1)) %>%
            mutate(feature = factor(feature, levels = featureNames))

        dfb2 <- dfb %>%
            mutate(
                csum = rev(cumsum(rev(percent))),
                pos = percent / 2 + lead(csum, 1),
                pos = ifelse(is.na(pos), percent / 2, pos),
                norm_csum = rev(cumsum(rev(norm_percent))),
                norm_pos = norm_percent / 2 + lead(norm_csum, 1),
                norm_pos = ifelse(is.na(norm_pos), norm_percent / 2, norm_pos)
            )

        apb1 <- ggplot(dfb, aes(x = "", y = percent, fill = feature)) +
            geom_col(color = "white") +
            scale_y_continuous(breaks = dfb2$pos, labels = dfb2$labels) +
            guides(fill = guide_legend(nrow = 2)) +
            scale_fill_manual(values = fil[1:nrow(dfb)]) +
            coord_polar(theta = "y") +
            theme(
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                axis.text = element_text(size = 12),
                legend.position = "top",
                legend.title = element_blank(),
                panel.background = element_rect(fill = "white")
            )

        apb2 <- ggplot(dfb, aes(x = "", y = norm_percent, fill = feature)) +
            geom_col(color = "white") +
            scale_y_continuous(breaks = dfb2$norm_pos,
                               labels = dfb2$norm_labels) +
            guides(fill = guide_legend(nrow = 2)) +
            scale_fill_manual(values = fil[1:nrow(dfb)]) +
            coord_polar(theta = "y") +
            theme(
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                axis.text = element_text(size = 12),
                legend.position = "top",
                legend.title = element_blank(),
                panel.background = element_rect(fill = "white")
            )

        print(pbar)
        print(plot_grid(apa1, apb1, nrow = 1,
                        labels = "Absolute counts", label_y = 1))
        print(plot_grid(apa2, apb2, nrow = 1,
                        labels = "Length-normalized counts", label_y = 1))

        if (!is.null(outPrefix)) {
            print(params)
            on.exit(dev.off(), add = TRUE)
        }
        if (verbose) message("[plot_peak_annotation] finished!\n")
        invisible(list(annotation = dfs, stat = dfa, simplified = dfb))
    }
}
