
#' @title Extract the longest transcript for each protein-coding genes
#
#' @description Gene level computations require selecting one transcript per gene to avoid bias by genes with multiple isoforms.
#' In ideal case, the most abundant transcript (principal or canonical isoform) should be chosen. However, the
#' most abundant isoform may vary depending on tissue type or physiological condition, the longest transcript is usually
#' the principal isoform, and alternatively spliced isoforms are not. This method get the longest transcript for each
#' gene. The longest transcript is defined as the isoform that has the longest CDS. In case of tie in the lengths of CDS,
#' the one with longer transcript is selected. If the lengths of transcripts tie again, the transcript with smaller id
#' is selected arbitrarily.
#'
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param plot logical, indicating whether feature length plots should be generated
#' @return a vector of transcript ids
#' @author Shuye Pu
#'
#' @examples
#' txdb <- AnnotationDbi::loadDb(system.file("data", "txdb.sql", package="GenomicPlotData"))
#' longestTx <- extract_longest_tx(txdb, plot=FALSE)
#'
#' @export extract_longest_tx
#'
#'
extract_longest_tx <- function(txdb, plot=FALSE){
   tl <- GenomicFeatures::transcriptLengths(txdb, with.utr5_len = TRUE, with.cds_len = TRUE, with.utr3_len = TRUE)
   pc <- tl[tl$cds_len > 0, ]  ## pc stands for protein-coding
   ## choose tx with longest cds for each gene
   longest_cds <- aggregate(pc$cds_len, list(pc$gene_id), max)
   colnames(longest_cds) <- c("gene_id", "cds_len")
   tx_longest_cds <- merge(pc, longest_cds)

   ## for two tx of the same gene, if cds_len are same, choose longer tx_len
   longest_tx <- aggregate(tx_longest_cds$tx_len, list(tx_longest_cds$gene_id, tx_longest_cds$cds_len), max)
   colnames(longest_tx) <- c("gene_id", "cds_len", "tx_len")
   longest_cdstx <- merge(longest_tx, tx_longest_cds)

   dup_genes <- longest_cdstx$gene_id[duplicated(longest_cdstx$gene_id)]
   dup_genestx <- longest_cdstx[longest_cdstx$gene_id %in% dup_genes,]
   ## for two tx of the same gene, if both cds_len and tx_len are the same, choose smaller tx_id (this is arbitrary)
   longest_cdstx_id <- aggregate(longest_cdstx$tx_id, list(longest_cdstx$gene_id, longest_cdstx$cds_len, longest_cdstx$tx_len), min)
   colnames(longest_cdstx_id) <- c("gene_id", "cds_len", "tx_len", "tx_id")
   longest_cdstxid <- merge(longest_cdstx, longest_cdstx_id)
   length(longest_cdstxid$gene_id[duplicated(longest_cdstxid$gene_id)])
   summary(longest_cdstxid)
   head(longest_cdstxid)

   npc <- tl %>% 
      filter(!gene_id %in% longest_cdstxid$gene_id) ## npc stands for non-protein-coding
   
   ## choose tx with longest length for each gene
   longest_npc <- aggregate(npc$tx_len, list(npc$gene_id), max)
   colnames(longest_npc) <- c("gene_id", "tx_len")
   tx_longest <- merge(npc, longest_npc)
   tx_longest_id <- aggregate(tx_longest$tx_id, list(tx_longest$gene_id, tx_longest$tx_len), min)
   colnames(tx_longest_id) <- c("gene_id", "tx_len", "tx_id")
   tx_longestid <- merge(tx_longest, tx_longest_id)

   longest_txid <- rbind(longest_cdstxid[, colnames(tl)], tx_longestid[, colnames(tl)])

   if(plot){

      ## plot intron, exon
      feature_list <- list("Exon"= exonsBy(txdb, by="tx", use.name=TRUE), "Intron"=intronsByTranscript(txdb, use.name=TRUE))
      plot_list <- list()
      pdf("Exon_intron_length_distribution.pdf", width=12, height=8)
      for(featureName in names(feature_list)){
         feature <- feature_list[[featureName]]
         feature_gr <- unlist(feature)
         feature_length <- width(feature_gr)
         length_df <- data.frame(tx=names(feature_gr), feature_length)

         p1 <- ggplot(length_df, aes(x=feature_length)) +
            geom_density() +
            scale_x_log10() +
            stat_ecdf(color="red") +
            scale_y_continuous(name = "Density", sec.axis = sec_axis( trans=~.*1, name="Probability")) +
            geom_vline(xintercept = c(mean(feature_length), median(feature_length)), linetype="dotted", color = c("brown4", "blue"), size=0.5) +
            annotate(
               "text", label = round(c(mean(feature_length), median(feature_length))),
               x = c(mean(feature_length), median(feature_length)), y = c(0, 0.1), size = 6, color = c("brown4", "blue")
            ) +
            xlab(paste(featureName, "length"))

         #print(p1)
         plot_list[[featureName]] <- p1
      }
      outp <- plot_grid(plot_list[[1]], plot_list[[2]], ncol = 2)
      print(outp)
      on.exit(dev.off(), add=TRUE)

      ## plot 5'utr, cds, 3'utr
      pdf("comparision_of_scaled_bins.pdf", width=12, height=8)
      p1 <- ggplot(tl, aes(x=cds_len)) +
         geom_density() +
         scale_x_log10() +
         stat_ecdf(color="red") +
         scale_y_continuous(name = "Density", sec.axis = sec_axis( trans=~.*1, name="Probability")) +
         geom_vline(xintercept = c(mean(tl$cds_len), median(tl$cds_len), 305, 22), linetype="dotted", color = c("green", "blue", "brown4", "magenta"), size=0.5) +
         annotate(
            "text", label = round(c(mean(tl$cds_len), median(tl$cds_len), 305, 22)),
            x = c(mean(tl$cds_len), median(tl$cds_len), 305, 22), y = c(0, 0.05, 0.1, 0.15), size = 4, color = c("green", "blue", "brown4", "magenta")
         )

      #p1
      p2 <- ggplot(longest_cdstxid, aes(x=cds_len)) +
         geom_density() +
         scale_x_log10() +
         stat_ecdf(color="red") +
         scale_y_continuous(name = "Density", sec.axis = sec_axis( trans=~.*1, name="Probability")) +
         geom_vline(xintercept = c(mean(longest_cdstxid$cds_len), median(longest_cdstxid$cds_len), 305, 31), linetype="dotted", color = c("green", "blue", "brown4", "magenta"), size=0.5) +
         annotate(
            "text", label = round(c(mean(longest_cdstxid$cds_len), median(longest_cdstxid$cds_len), 305, 31)),
            x = c(mean(longest_cdstxid$cds_len), median(longest_cdstxid$cds_len), 305, 31), y = c(0, 0.05, 0.1, 0.15), size = 4, color = c("green", "blue", "brown4", "magenta")
         )

      #p2

      p3 <- ggplot(tl, aes(x=utr5_len)) +
         geom_density() +
         scale_x_log10() +
         stat_ecdf(color="red") +
         scale_y_continuous(name = "Density", sec.axis = sec_axis( trans=~.*1, name="Probability")) +
         geom_vline(xintercept = c(mean(tl$utr5_len), median(tl$utr5_len), 45, 4), linetype="dotted", color = c("green", "blue", "brown4", "magenta"), size=0.5) +
         annotate(
            "text", label = round(c(mean(tl$utr5_len), median(tl$utr5_len), 45, 4)),
            x = c(mean(tl$utr5_len), median(tl$utr5_len), 45, 4), y = c(0, 0.05, 0.1, 0.15), size = 4, color = c("green", "blue", "brown4", "magenta")
         )
      #p3
      p4 <- ggplot(longest_cdstxid, aes(x=utr5_len)) +
         geom_density() +
         scale_x_log10() +
         stat_ecdf(color="red") +
         scale_y_continuous(name = "Density", sec.axis = sec_axis( trans=~.*1, name="Probability")) +
         geom_vline(xintercept = c(mean(longest_cdstxid$utr5_len), median(longest_cdstxid$utr5_len), 45, 4), linetype="dotted", color = c("green", "blue", "brown4", "magenta"), size=0.5) +
         annotate(
            "text", label = round(c(mean(longest_cdstxid$utr5_len), median(longest_cdstxid$utr5_len), 45, 4)),
            x = c(mean(longest_cdstxid$utr5_len), median(longest_cdstxid$utr5_len), 45, 4), y = c(0, 0.05, 0.1, 0.15), size = 4, color = c("green", "blue", "brown4", "magenta")
         )
      #p4

      p5 <- ggplot(tl, aes(x=utr3_len)) +
         geom_density() +
         scale_x_log10() +
         stat_ecdf(color="red") +
         scale_y_continuous(name = "Density", sec.axis = sec_axis( trans=~.*1, name="Probability")) +
         geom_vline(xintercept = c(mean(tl$utr3_len), median(tl$utr3_len), 155, 10), linetype="dotted", color = c("green", "blue", "brown4", "magenta"), size=0.5) +
         annotate(
            "text", label = round(c(mean(tl$utr3_len), median(tl$utr3_len), 155, 10)),
            x = c(mean(tl$utr3_len), median(tl$utr3_len), 155, 10), y = c(0, 0.05, 0.1, 0.15), size = 4, color = c("green", "blue", "brown4", "magenta")
         )
      #p5
      p6 <- ggplot(longest_cdstxid, aes(x=utr3_len)) +
         geom_density() +
         scale_x_log10() +
         stat_ecdf(color="red") +
         scale_y_continuous(name = "Density", sec.axis = sec_axis( trans=~.*1, name="Probability")) +
         geom_vline(xintercept = c(mean(longest_cdstxid$utr3_len), median(longest_cdstxid$utr3_len), 155, 16), linetype="dotted", color = c("green", "blue", "brown4", "magenta"), size=0.5) +
         annotate(
            "text", label = round(c(mean(longest_cdstxid$utr3_len), median(longest_cdstxid$utr3_len), 155, 16)),
            x = c(mean(longest_cdstxid$utr3_len), median(longest_cdstxid$utr3_len), 155, 16), y = c(0, 0.05, 0.1, 0.15), size = 4, color = c("green", "blue", "brown4", "magenta")
         )
      #p6

      outp <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)
      print(outp)
      on.exit(dev.off(), add=TRUE)
   }

   invisible(longest_txid)
}


#' @title Extract genomic features from TxDb object
#
#' @description Extract genomic coordinates and make bed and bed 12 files from a TxDb object for a variety of annotated genomic features.
#' The output of this function is a list. The first element of the list is a GRanges object that provide the start and end information.
#' The second element is a GRangesList providing information for subcomponents. For "utr3", "utr5", "cds" and "transcript", the returned GRanges
#' object has only one range for each feature, denoting the start and end of the feature in one transcript, and the range may contain introns;
#' the returned GrangesList object is a list of exons, belonging to one transcript and indexed on transcript id. The third element is a bed12 file.
#' For "exon", "intron" and "gene", the GRanges object denotes ranges of individual exon or intron, and the GrangesList object is a list of exons or
#' introns belonging to one transcript and indexed on transcript id. The third element is a bed6 file.
#'
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param featureName one of the gene feature in c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene")
#' @param featureSource the name of the gtf/gff3 file or the online database from which txdb is derived, used as name of output file
#' @param export logical, indicating if the bed file should be produced
#' @param longest logical, indicating whether the output should be limited to the longest transcript of each gene
#' @param protein_coding logical, indicating whether to limit to protein_coding genes
#'
#' @return a list of three objects, the first is a GRanges object, the second is a GRangesList object, the last is the output file name
#' @author Shuye Pu
#'
#' @examples
#' txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
#' output <- get_genomic_feature_coordinates(txdb, featureName="cds", featureSource="gencode",
#' export=FALSE, longest=TRUE, protein_coding=TRUE)
#' @export get_genomic_feature_coordinates
#'
#'
get_genomic_feature_coordinates <- function(txdb, featureName, featureSource=NULL, export=FALSE, longest=FALSE, protein_coding=FALSE){

   if(is.null(featureSource)){featureSource <- "unknown_source"}
   outfile <- NULL

   seqInfo <- seqinfo(txdb)

   tl <- transcriptLengths(txdb, with.utr5_len = TRUE, with.cds_len = TRUE, with.utr3_len = TRUE)
   tl_protein_coding <- tl[tl$cds_len > 0, ]

   feature <- NULL

   if(featureName == "utr3"){
      feature <- threeUTRsByTranscript(txdb, use.name=FALSE) # grl
   }else if(featureName == "utr5"){
      feature <- fiveUTRsByTranscript(txdb, use.name=FALSE) # grl
   }else if(featureName == "intron"){
      feature <- intronsByTranscript(txdb, use.name=FALSE) #grl
      if(protein_coding) feature <- feature[tl_protein_coding$tx_id]
   }else if(featureName == "exon"){
      feature <- exonsBy(txdb, by="tx", use.name=FALSE) #grl
      if(protein_coding) feature <- feature[tl_protein_coding$tx_id]
   }else if(featureName == "cds"){
      feature <- cdsBy(txdb, by="tx", use.name=FALSE) # grl
   }else if(featureName == "transcript"){
      feature <- exonsBy(txdb, by="tx", use.name=FALSE) #grl
      if(protein_coding) feature <- feature[tl_protein_coding$tx_id]
   }else if(featureName == "gene"){
      feature <- genes(txdb)                       # gr
      feature <- as(split(feature, f=feature$gene_id), "GRangesList") # grl
      if(protein_coding) feature <- feature[unique(tl_protein_coding$gene_id)]
   }else {
      stop("Feature is not defined!")
   }

   seqinfo(feature) <- seqInfo
   if(longest){
      longest_tx <- extract_longest_tx(txdb)
      if(featureName == "gene"){
         feature_longest <- feature[names(feature) %in% longest_tx$gene_id] # protein-coding
         seqinfo(feature_longest) <- seqInfo
      }else{
         feature_longest <- feature[names(feature) %in% longest_tx$tx_id]
         seqinfo(feature_longest) <- seqInfo
      }

      if(featureName %in% c("cds", "utr5", "utr3", "transcript")){
         gr_feature_longest <- rtracklayer::asBED(feature_longest) ## convert each element of Grangeslist to one Grange with blocks info as metadata
         names(gr_feature_longest) <- gr_feature_longest$name
         seqinfo(gr_feature_longest) <- seqInfo
         if(export){
            outfile <- paste(featureSource, "_", featureName, "_longest.bed12", sep="")
            export.bed(gr_feature_longest, outfile)
         }
      }
      if(featureName %in% c("gene", "intron", "exon")){
         gr_feature_longest <- unlist(feature_longest, use.names=FALSE) ## convert each element of Grangeslist to multiple Granges
         seqinfo(gr_feature_longest) <- seqInfo
         if(export){
            outfile <- paste(featureSource, "_", featureName, "_longest.bed", sep="")
            export.bed(gr_feature_longest, outfile)
         }
      }
      invisible(list("GRanges"=gr_feature_longest, "GRangesList"=feature_longest, "Output"=outfile))
   }else{
      if(featureName %in% c("cds", "utr5", "utr3", "transcript")){
         gr_feature <- asBED(feature) ## convert each element of Grangeslist to one Grange with blocks info as metadata
         names(gr_feature) <- gr_feature$name
         seqinfo(gr_feature) <- seqInfo
         if(export){
            outfile <- paste(featureSource, "_", featureName, ".bed12", sep="")
            export.bed(gr_feature, outfile)
         }
      }
      if(featureName %in% c("gene", "intron", "exon")){
         gr_feature <- unlist(feature, use.names=FALSE) ## convert Grangeslist object to GRanges object
         seqinfo(gr_feature) <- seqInfo
         if(export){
            outfile <- paste(featureSource, "_", featureName, ".bed", sep="")
            export.bed(gr_feature, outfile)
         }
      }
      invisible(list("GRanges"=gr_feature, "GRangesList"=feature, "Output"=outfile))
   }
}


#' @title Demarcate genes into promoter, gene body  and TTS features
#
#' @description This is a helper function for 'plot_3parts_metagene', used to speed up plotting of multiple data sets with the same configuration. Use featureName='transcript' and meta=FALSE and longest=TRUE for genes.
#'
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param featureName one of the gene feature in c("utr3", "utr5", "cds", "intron", "exon", "transcript")
#' @param meta logical, indicating whether a metagene (intron excluded) or gene (intron included) plot should be produced
#' @param nbins an integer defines the total number of bins
#' @param fiveP extension out of the 5' boundary of gene
#' @param threeP extension out of the 3' boundary of gene
#' @param longest logical, indicating whether the output should be limited to the longest transcript of each gene
#' @param protein_coding logical, indicating whether to limit to protein_coding genes
#'
#' @return a named list with the elements c("windowRs", "nbins", "scaled_bins", "fiveP", "threeP", "meta", "longest")
#' @author Shuye Pu
#'
#' @examples
#' txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
#'
#' gf <- prepare_3parts_genomic_features(txdb, meta=FALSE, nbins=100, fiveP=1000, threeP=1000,
#' longest=FALSE)
#'
#' @export prepare_3parts_genomic_features

prepare_3parts_genomic_features <- function(txdb, featureName="transcript", meta=TRUE, nbins=100, fiveP=1000, threeP=1000, longest=TRUE, protein_coding=TRUE){
   ## prepare transcripts

   print("Preparing features ...")

   ## prepare transcripts that are suitable for overlap
   if(featureName %in% c("utr5", "cds", "utr3")){
      featureName <- ifelse(meta, featureName, paste0(featureName,"(with intron)"))
   }else if(featureName == "transcript"){
      featureName <- ifelse(meta, featureName, "gene")
   }

   five <- -fiveP/1000
   five <- paste0(five, "K")
   if(fiveP==0) five=""
   three <- threeP/1000
   three <- paste0(three, "K")
   if(threeP==0) three=""

   featureNames <-  c(five, featureName, three)

   gn <- get_genomic_feature_coordinates(txdb, featureName, longest=longest, protein_coding=protein_coding)

   if(meta){
      feature <- gn$GRangesList
   }else{
      feature <- as(split(gn$GRanges, as.factor(names(gn$GRanges))), "GRangesList")
   }

   wf <- vapply(as.list(width(feature)), sum, numeric(1))
   means <- c(promoter=fiveP, median(wf), TTS=threeP)
   scaled_bins <- round(means*nbins/sum(means))
   selected_tx <- names(feature[wf > scaled_bins[2]])

   promoter <- flank(gn$GRanges, width=fiveP, both=FALSE, start=TRUE, ignore.strand=FALSE)
   TTS <- flank(gn$GRanges, width=threeP, both=FALSE, start=FALSE, ignore.strand=FALSE)

   windowRs <- list(as(split(promoter, as.factor(names(promoter))), "GRangesList")[selected_tx],
                    feature[selected_tx],
                    as(split(TTS, as.factor(names(TTS))), "GRangesList")[selected_tx])

   names(windowRs) <- featureNames
   names(scaled_bins) <- featureNames

   print("median sizes for features")
   print(means)
   print("bin sizes for features")
   print(scaled_bins)
   print("number of transcripts")
   print(vapply(windowRs, length, numeric(1)))

   invisible(list("windowRs"=windowRs, "nbins"=nbins, "scaled_bins"=scaled_bins, "fiveP"=fiveP, "threeP"=threeP, "meta"=meta, "longest"=longest))
}

#' @title Demarcate genes into promoter, 5'UTR, CDS, 3'UTR and TTS features
#
#' @description This is a helper function for 'plot_5parts_metagene', used to speed up plotting of multiple data sets with the same configuration.
#'
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param meta logical, indicating whether a metagene (intron excluded) or gene (intron included) plot should be produced
#' @param nbins an integer defines the total number of bins
#' @param fiveP extension out of the 5' boundary of gene
#' @param threeP extension out of the 3' boundary of gene
#' @param longest logical, indicating whether the output should be limited to the longest transcript of each gene
#' @param subsetTx a vector of transcript ids for subseting the genome
#'
#' @return a named list with the elements c("windowRs", "nbins", "scaled_bins", "fiveP", "threeP", "meta", "longest")
#' @author Shuye Pu
#'
#' @examples
#' txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
#'
#' gf <- prepare_5parts_genomic_features(txdb, meta=TRUE, nbins=100, fiveP=1000, threeP=1000,
#' longest=TRUE)
#'
#' @export prepare_5parts_genomic_features
#'
prepare_5parts_genomic_features <- function(txdb, meta=TRUE, nbins=100, fiveP=1000, threeP=1000, longest=TRUE, subsetTx=NULL){
   ## prepare transcripts
   
   print("Preparing genomic features ... ")
   five <- -fiveP/1000
   five <- paste0(five, "K")
   if(fiveP==0) five=""
   three <- threeP/1000
   three <- paste0(three, "K")
   if(threeP==0) three=""
   featureNames <-  c(five, "5'UTR", "CDS", "3'UTR", three)

   if(!meta) longest <- TRUE # always use the longest transcript to represent the gene

   utr5 <- get_genomic_feature_coordinates(txdb, "utr5", longest=longest, protein_coding=TRUE)
   utr3 <- get_genomic_feature_coordinates(txdb, "utr3", longest=longest, protein_coding=TRUE)
   cds <- get_genomic_feature_coordinates(txdb, "cds", longest=longest, protein_coding=TRUE)
   gn <- get_genomic_feature_coordinates(txdb, "transcript", longest=longest, protein_coding=TRUE)

   #feature_coordinates <- parallel_feature_coordinates(txdb, featureNames=featurelist, longest=longest, protein_coding = TRUE)

   promoter <- flank(gn$GRanges, width=fiveP, both=FALSE, start=TRUE, ignore.strand=FALSE)
   TTS <- flank(gn$GRanges, width=threeP, both=FALSE, start=FALSE, ignore.strand=FALSE)

   if(meta){
      utr5_grl <- utr5$GRangesList
      cds_grl <- cds$GRangesList
      utr3_grl <- utr3$GRangesList
   }else{
      utr5_grl <- as(split(utr5$GRanges, as.factor(names(utr5$GRanges))), "GRangesList")
      cds_grl <- as(split(cds$GRanges, as.factor(names(cds$GRanges))), "GRangesList")
      utr3_grl <- as(split(utr3$GRanges, as.factor(names(utr3$GRanges))), "GRangesList")
   }

   grls <- list("5'UTR"=utr5_grl, "CDS"=cds_grl, "3'UTR"=utr3_grl)
   l3 <- vapply(grls, function(feature)median(vapply(as.list(width(feature)), sum, numeric(1))), numeric(1))

   means <- c(promoter=fiveP, l3, TTS=threeP)
   names(means) <- featureNames
   scaled_bins <- round(means*nbins/sum(means))
   names(scaled_bins) <- featureNames

   selected_tx <- lapply(names(grls), function(x){
      len <- vapply(as.list(width(grls[[x]])),sum, numeric(1)) # len is named vector, where names are the tx_ids
      y <- names(len)[which(len >= scaled_bins[x])]
   })

   names(selected_tx) <- names(grls)
   if(!is.null(subsetTx)){
      tl <- transcriptLengths(txdb, with.utr5_len = TRUE, with.cds_len = TRUE, with.utr3_len = TRUE)
      txids <- tl %>%
         filter(tx_name %in% subsetTx) %>%
         select(tx_id)
      selected_tx[["custom"]] <- as.character(txids[,1])
   }
   selected_tx <- Reduce(intersect, selected_tx)

   windowRs <- list(as(split(promoter, as.factor(names(promoter))), "GRangesList")[selected_tx],
                    utr5_grl[selected_tx],
                    cds_grl[selected_tx],
                    utr3_grl[selected_tx],
                    as(split(TTS, as.factor(names(TTS))), "GRangesList")[selected_tx])

   names(windowRs) <- featureNames

   print("median sizes for features")
   print(means)
   print("bin sizes for features")
   print(scaled_bins)
   print("number of transcripts")
   print(vapply(windowRs, length, numeric(1)))

   invisible(list("windowRs"=windowRs, "nbins"=nbins, "scaled_bins"=scaled_bins, "fiveP"=fiveP, "threeP"=threeP, "meta"=meta, "longest"=longest))
}

#' @title Get genomic coordinates of features of protein-coding genes
#
#' @description Get genomic coordinates of promoter, 5'UTR, CDS, 3'UTR, TTS and intron for the longest transcript of protein-coding genes. The range of promoter is defined by fiveP and dsTSS upstream and downstream TSS, respectively, the TTS ranges from the 3' end of the gene to threeP downstream.
#'
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param fiveP extension upstream of the 5' boundary of genes
#' @param threeP extension downstream of the 3' boundary of genes
#' @param dsTSS range of promoter extending downstream of TSS
#'
#' @return a GRangesList object
#' @author Shuye Pu
#'
#' @examples
#' txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
#'
#' f <- get_txdb_features(txdb, dsTSS=100, fiveP=0, threeP=0)
#'
#' @export get_txdb_features
#'
get_txdb_features <- function(txdb, fiveP=1000, dsTSS=300, threeP=1000){
   
   utr5 <- get_genomic_feature_coordinates(txdb, "utr5", longest=TRUE, protein_coding=TRUE)
   utr3 <- get_genomic_feature_coordinates(txdb, "utr3", longest=TRUE, protein_coding=TRUE)
   cds <- get_genomic_feature_coordinates(txdb, "cds", longest=TRUE, protein_coding=TRUE)
   intron <- get_genomic_feature_coordinates(txdb, "intron", longest=TRUE, protein_coding=TRUE)
   
   gene <- get_genomic_feature_coordinates(txdb, "transcript", longest=TRUE, protein_coding=TRUE)$GRanges
   
   features <- GRangesList("5'UTR"=unlist(utr5$GRangesList, use.names=TRUE),
                           "3'UTR"=unlist(utr3$GRangesList, use.names=TRUE),
                           "CDS"=unlist(cds$GRangesList, use.names=TRUE),
                           "Intron"=unlist(intron$GRangesList, use.names=TRUE),
                           compress=FALSE)
   if(fiveP > 0){
      promoter <- GenomicRanges::promoters(gene, upstream=fiveP, downstream=dsTSS, use.names=TRUE)
      features[["Promoter"]] <- promoter
   }else{
      promoter <- NULL
   }
   
   if(threeP > 0){
      TTS <- GenomicRanges::flank(gene, width=threeP, both=FALSE, start=FALSE, ignore.strand=FALSE)
      features[["TTS"]] <- TTS
   }else{
      TTS <- NULL
   }
   
   tx <- transcripts(txdb)
   txmap <- tx$tx_name
   names(txmap) <- tx$tx_id
   
   for(aname in names(features)){
      f <- features[[aname]]
      mcols(f) <- NULL
      f$tx_name <- txmap[names(f)]
      features[[aname]] <- f
   }
   
   invisible(features)
}

#' @title Get the number of peaks overlapping each feature of all protein-coding genes.
#
#' @description Annotate each peak with genomic features based on overlap, and produce summary statistics for distribution of peaks in features of protein-coding genes.
#'
#' @param peak a GRanges object defining query ranges
#' @param features a GRangesList object representing genomic features
#' @param stranded logical, indicating whether the overlap should be strand-specific
#'
#' @return a list object
#' @author Shuye Pu
#'
#' @examples
#' txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
#'
#' f <- get_txdb_features(txdb, dsTSS=100, fiveP=0, threeP=1000)
#' p <- RCAS::importBed(system.file("data", "test_chip_peak_chr19.bed", package="GenomicPlotData"))
#' 
#' ann <- get_targeted_genes(peak=p, features=f, stranded=FALSE)
#' 
#' pp <- RCAS::importBed(system.file("data", "test_clip_peak_chr19.bed", package="GenomicPlotData"))
#' ann <- get_targeted_genes(peak=pp, features=f, stranded=TRUE)
#'
#' @export get_targeted_genes

get_targeted_genes <- function(peak, features, stranded=TRUE){
   
   num_peaks <- length(peak)
   num_genes <- length(unique(features$CDS$tx_name))
   
   annot_list <- lapply(names(features), function(feature){
      featureGr <- features[[feature]]
      featureOverlaps <- findOverlaps(peak, featureGr, ignore.strand=!stranded)
      peak_df <- gr2df(peak[featureOverlaps@from])
      if(is.null(peak_df$strand)) peak_df$strand <- rep("*", nrow(peak_df))
      peak_df <- peak_df %>%
         mutate(chrPeak=as.character(chr),
                startPeak=as.integer(start)-1,
                endPeak=as.integer(end),
                idPeak=as.character(id),
                widthPeak=as.integer(width),
                strandPeak=as.character(strand),
                .keep="unused")
      feature_df <- featureGr[featureOverlaps@to]
      names(feature_df) <- seq_along(feature_df)
      feature_df <- gr2df(feature_df) %>%
         mutate(chrfeature=as.character(chr),
                startfeature=as.integer(start)-1,
                endfeature=as.integer(end),
                widthfeature=as.integer(width),
                strandfeature=as.character(strand),
                .keep="unused")
      
      ot <- cbind(peak_df, feature_df) %>%
         select(chrPeak, startPeak, endPeak, idPeak, widthPeak, strandPeak, chrfeature, startfeature, endfeature, widthfeature, strandfeature, tx_name) %>%
         mutate(feature_name=feature)
      return(ot)
   })
   
   annot_table <- bind_rows(annot_list)
   
   annot_count <- as.data.frame(annot_table %>%
      group_by(tx_name) %>%
      count(feature_name))
   
   overlap_genes <- length(unique(annot_table$tx_name))
   overlap_peaks <- length(unique(annot_table$idPeak))
   
   
   annot_df <- data.frame(tx_name=unique(features$CDS$tx_name), Promoter=0, `5'UTR`=0, CDS=0, `3'UTR`=0, TTS=0, Intron=0) %>% rename("5'UTR" = "X5.UTR", "3'UTR" = "X3.UTR")
   rownames(annot_df) <- annot_df$tx_name
   
   system.time(
   for(i in 1:nrow(annot_count)){
      txi <- annot_count[i,1]
      fn <- annot_count[i,2]
      count <- annot_count[i,3]
      annot_df[txi, fn] <- count
   })
   
   annot_df <- annot_df %>%
      mutate(Exon = `5'UTR` + CDS + `3'UTR`) %>%
      mutate(Transcript = Intron + Exon)
   
   feature_counts <- apply(annot_df[, c("Promoter", "5'UTR", "CDS", "3'UTR", "TTS", "Intron", "Exon", "Transcript")], 2, sum)
   
   invisible(list(gene_table=annot_df, peak_table=annot_table, num_peak=num_peaks, num_gene=num_genes, feature_count=feature_counts, overlap_peak=overlap_peaks, overlap_gene=overlap_genes))
}

#' @title Make TxDb object from a GTF file for a subset of genes.
#
#' @description Make a partial TxDb object given a GTF file and a list of gene names in a file or in a character vector.
#'
#' @param gtfFfile path to a GTF file
#' @param geneList path to a line-delimited text file with one gene name on each line, or a character vector of gene names
#'
#' @return a TxDb object
#' @author Shuye Pu
#' 
#' @export makeSubTxDbFromGtf

makeSubTxDbFromGtf <- function(gtfFile, geneList, geneCol=1){
   gff <- RCAS::importGtf(saveObjectAsRds = TRUE, filePath = gtfFile)
   if(length(geneList == 1)){
      aList <- read.delim2(geneList, comment.char = "#")
      geneList <- as.character(aList[, geneCol])
   }
   subgff <- gff[gff$gene_name %in% geneList]
   TxDb <- makeTxDbFromGRanges(subgff)
   
   return(TxDb)
}

#' @title Translate gene names to transcript ids using a GTF file for a subset of genes.
#
#' @description Given a list of gene names in a file or in a character vector, turn them into a vector of transcript ids
#'
#' @param gtfFfile path to a GTF file
#' @param geneList path to a line-delimited text file with one gene name on each line, or a character vector of gene names
#'
#' @return a vector
#' @author Shuye Pu
#' 
#' @export gene2tx

gene2tx <- function(gtfFile, geneList, geneCol=1){
   gff <- RCAS::importGtf(saveObjectAsRds = TRUE, filePath = gtfFile)
   if(length(geneList) == 1){
      aList <- read.delim2(geneList, comment.char = "#")
      geneList <- as.character(aList[, geneCol])
   }
   subgff <- gff[gff$gene_name %in% geneList]
   
   tx <- as.character(unique(subgff$transcript_id))
   
   return(tx)
}

