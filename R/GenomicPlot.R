
#' @import methods
NULL

#library(Matrix)
library(parallel)
library(data.table)
#library(MatrixGenerics)
library(tidyr)
library(dplyr)
library(DESeq2)
library(magrittr)
library(cowplot)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(genomation)
library(GenomicRanges)
library(plyranges)
library(GenomicFeatures)
library(GenomicAlignments)
library(gridExtra)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(R.utils)
library(hrbrthemes)
library(pheatmap)
library(scales)
library(RMariaDB)
library(RCAS)
library(tictoc)
#library(gg.layers)



#' start_parallel
#'
#' Method for starting a virtual cluster needed for parallel processing
#'
#' @param nc, a positive integer greater than 1
#' @return an object of class c("SOCKcluster", "cluster"), depending on platform
#'
#' @examples
#' cl <- start_parallel(5L)
#'
#' @export start_parallel
#'
start_parallel <- function(nc){

  n.cores <- parallel::detectCores()
  fnc <- min(nc, n.cores-1)

  switch(Sys.info()[['sysname']],
         Windows = {my.cluster <- parallel::makeCluster(fnc, type = "PSOCK")},
         Linux = {my.cluster <- parallel::makeCluster(fnc, type = "FORK")}
  )

  print(paste("Using", fnc, "out of", n.cores, "cores!"))

  my.cluster
}

#' stop_parallel
#
#' Method for stopping a virtual cluster needed for parallel processing
#'
#' @param cl, a cluster or SOCKcluster object depending on platform
#' @return NULL
#'
#' @examples
#' cl <- start_parallel(5L)
#' stop_parallel(cl)
#'
#' @export stop_parallel
#'

stop_parallel <- function(cl){
  parallel::stopCluster(cl = cl)
}


#' extract_longest_tx
#
#' Gene level computations require selecting one transcript per gene to avoid bias by genes with multiple isoforms.
#' In ideal case, the most abundant transcript (principal or canonical isoform) should be chosen. However, the
#' most abundant isoform may vary depending on tissue type or physiological condition, the longest transcript is usually
#' the principal isoform, and alternatively spliced isoforms are not. This method get the longest transcript for each
#' gene. The longest transcript is defined as the isoform that has the longest CDS. In case of tie in the lengths of CDS,
#' the one with longer transcript is selected. If the lengths of transcripts tie again, the transcript with smaller id
#' is selected arbitrarily.
#'
#' @param txdb, a TxDb object defined in GenomicFeatures package
#' @param plot, a boolean object indicating whether feature length plots should be generated
#' @return a vector of transcript ids
#'
#' @examples
#' txdb <- makeTxDbFromEnsembl("Saccharomyces cerevisiae", server="useastdb.ensembl.org")
#' longestTx <- extract_longest_tx(txdb, plot=T)
#' @export extract_longest_tx
#'
#'
extract_longest_tx <- function(txdb, plot=FALSE){
  tl <- transcriptLengths(txdb, with.utr5_len = T, with.cds_len = T, with.utr3_len = T)
  tl <- tl[tl$cds_len > 0, ]
  ## choose tx with longest cds for each gene
  longest_cds <- aggregate(tl$cds_len, list(tl$gene_id), max)
  colnames(longest_cds) <- c("gene_id", "cds_len")
  tx_longest_cds <- merge(tl, longest_cds)

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

  if(plot){
    ## plot intron, exon
    feature_list <- list("Exon"= exonsBy(txdb, by="tx", use.name=T), "Intron"=intronsByTranscript(txdb, use.name=T))
    plot_list <- list()
    pdf("Exon_intron_length_distribution.pdf", width=12, height=8)
    for(featureName in names(feature_list)){
      feature <- feature_list[[featureName]]
      feature_gr <- stack(feature)
      feature_length <- end(feature_gr) - start(feature_gr) + 1

      length_df <- data.frame(tx=feature_gr$name, feature_length)

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
    dev.off()

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
    dev.off()
  }

  return(longest_cdstxid)
}


#' gtf_to_bed_longest_tx
#
#' Make bed and bed 12 files from gtf file for protein coding genes. For "utr3", "utr5", "cds" and "transcript",  the output granges object
#' gives start and end of the entire feature, while the grangslist object gives start and end of each exonic segment. For "exon"
#' and "intron", the output granges object gives start and end of each feature, while the grangslist object also provides
#' transcript ID in addition to these information. For "gene", the granges object and the grangeslist object are equivalent.
#' If 'longest'=TRUE, only the longest transcript of each gene is considered, the exons, introns, 3'utr and 5'utr are limited to
#' the longest transcript of protein-coding genes.
#'
#' @param txdb, a TxDb object defined in GenomicFeatures package
#' @param feaureNmae, one of the gene feature in c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene")
#' @param featureSource, the name of the gtf/gff3 file or the online database from which txdb is derived, used as name of output file
#' @param export, a boolean object to indicate if the bed file should be produced
#' @param longest, a boolean object to indicate whether the output should be limited to the longest transcript of each gene
#' @return a list of two objects, the first is a GRanges object, the second is a GRangesList object
#'
#' @examples
#' txdb <- makeTxDbFromEnsembl("Saccharomyces cerevisiae", server="useastdb.ensembl.org")
#' output <- gtf_to_bed_longest_tx(txdb, featureName="utr3", featureSource="Ensembl", export=F, longest=T)
#' @export gtf_to_bed_longest_tx
#'
#'
gtf_to_bed_longest_tx <- function(txdb, featureName, featureSource=NULL, export=FALSE, longest=TRUE){

 # featureName <- "transcript"
  longest_tx <- extract_longest_tx(txdb)
  if(!is.null(featureSource)){featureSource <- gsub("\\.gtf", "", featureSource)}
  feature <- NULL

  if(featureName == "utr3"){
    feature <- threeUTRsByTranscript(txdb, use.name=F) # grl
  }else if(featureName == "utr5"){
    feature <- fiveUTRsByTranscript(txdb, use.name=F) # grl
  }else if(featureName == "intron"){
    feature <- intronsByTranscript(txdb, use.name=F) #grl
    tl <- transcriptLengths(txdb, with.utr5_len = T, with.cds_len = T, with.utr3_len = T)
    tl_protein_coding <- tl[tl$cds_len > 0, ]
    feature <- feature[tl_protein_coding$tx_id]
  }else if(featureName == "exon"){
    feature <- exonsBy(txdb, by="tx", use.name=F) #grl
    tl <- transcriptLengths(txdb, with.utr5_len = T, with.cds_len = T, with.utr3_len = T)
    tl_protein_coding <- tl[tl$cds_len > 0, ]
    feature <- feature[tl_protein_coding$tx_id]
  }else if(featureName == "cds"){
    feature <- cdsBy(txdb, by="tx", use.name=F) # grl
  }else if(featureName == "transcript"){
    feature <- exonsBy(txdb, by="tx", use.name=F) #grl
    tl <- transcriptLengths(txdb, with.utr5_len = T, with.cds_len = T, with.utr3_len = T)
    tl_protein_coding <- tl[tl$cds_len > 0, ]
    feature <- feature[tl_protein_coding$tx_id]
  }else if(featureName == "gene"){
    feature <- genes(txdb)                       # gr
    feature <- split(feature, f=feature$gene_id)
    tl <- transcriptLengths(txdb, with.utr5_len = T, with.cds_len = T, with.utr3_len = T)
    tl_protein_coding <- tl[tl$cds_len > 0, ]
    feature <- feature[unique(tl_protein_coding$gene_id)]
  }else {
    stop("Feature is not defined!")
  }

  if(longest){
    if(featureName == "gene"){
      feature_longest <- feature[names(feature) %in% longest_tx$gene_id]
    }else{
      feature_longest <- feature[names(feature) %in% longest_tx$tx_id]
    }

    if(featureName %in% c("cds", "utr5", "utr3", "transcript")){
      gr_feature_longest <- asBED(feature_longest) ## convert each element of Grangeslist to one Grange with blocks info as metadata
      if(export){
        export.bed(gr_feature_longest, paste(featureSource, "_", featureName, "_longest.bed12", sep=""))
      }
    }
    if(featureName %in% c("gene", "intron", "exon")){
      gr_feature_longest <- stack(feature_longest) ## convert each element of Grangeslist to multiple Granges
      if(export){
        export.bed(gr_feature_longest, paste(featureSource, "_", featureName, "_longest.bed", sep=""))
      }
    }
    return(list("GRanges"=gr_feature_longest, "GRangesList"=feature_longest))
  }else{
    if(featureName %in% c("cds", "utr5", "utr3", "transcript")){
      gr_feature <- asBED(feature) ## convert each element of Grangeslist to one Grange with blocks info as metadata
      if(export){
        export.bed(gr_feature, paste(featureSource, "_", featureName, ".bed12", sep=""))
      }
    }
    if(featureName %in% c("gene", "intron", "exon")){
      gr_feature <- stack(feature) ## convert each element of Grangeslist to multiple Granges
      if(export){
        export.bed(gr_feature, paste(featureSource, "_", featureName, ".bed", sep=""))
      }
    }
    return(list("GRanges"=gr_feature, "GRangesList"=feature))

  }

}

#' @title  plot_start_end_feature
#
#'@describeIn   Plot reads or peak signal intensity of samples in the query files around stat and end of genomic features. The upstream and downstream windows
#' can be given separately, within the window, a smaller window can be defined to test statistical difference between the signal intensity
#' among samples. A line plot and a boxplot are displayed side by side for both start and end of feature.
#'
#' @param queryfiles, a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param querylabels, a vector of charactor strings serving as short labels of the queryfiles, will be used as sample labels in the plots
#' @param txdb, a TxDb object defined in GenomicFeatures package
#' @param feaureNmae, one of the gene feature in c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene")
#' @param CLIP_reads, a boolean object to indicate if the bam reads should be shifted to the -1 positon at the 5' of the reads
#' @param binsize, an integer defines bin size for intensity calculation
#' @param fix_width, an integer defines how long should the reads should be extended to
#' @param norm, an boolean indicates if the intensities should be normalized by the sample library sizes
#' @param longest, a boolean object to indicate whether the output should be limited to the longest transcript of each gene
#' @param ext, a vector of four integers defining upstream and downstream boundaries of the plot window, flanking the start and end of features
#' @param hl, a vector of four integers defining upstream and downstream boundaries of the highlight window, flanking the start and end of features
#' @param randomize, a boolean indicate if randomized feature should generated and used as a contrast to the real feature
#' @param stranded, a boolean indicate whether the strand of the feature should be considered
#' @param scale, a boolean indicates whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smo, a boolean indicates whether the line should smoothed with a spline smoothing algorithm
#' @param heatmap, a boolean indicates whether a heatmap of the score matrix should be generated
#' @param rm.outlier, a boolean indicates whether a row with abnormally high values in the score matrix should be removed
#' @param genome, genome of the features, the program is mostly concerned with human, if non-human genome is used, certain features may not work
#' @param outPrefix, a character string specifying output file prefix for plots (outPrefix.pdf)
#'
#' @return a list of two objects, the first is a GRanges object, the second is a GRangesList object
#'
#' @examples
#' txdb <- makeTxDbFromEnsembl("Homo sapiens", release=75)
#' queryfiles <- "YTHDF2.merged.thUni.bam"
#' querylabels <- "YTHDF2_merged"
#' inputfiles <- "Input_YTHDF2.thUni.bam"
#' inputlabels <- "YTHDF2_Input"
#' ext <- c(-500, 200, -200, 500)
#' hl <- c(-50, 50, -50, 50)
#' op <- "ratio_YTHDF2_merged_at_intron_boundary"
#'
#' plot_start_end_feature(queryfiles=queryfiles, querylabels=querylabels, inputfiles=inputfiles, inputlabels=inputlabels, txdb=txdb, featureName="intron", CLIP_read=F, binsize=10, fix_width=0,
#' longest=T, ext=ext, hl=hl, randomize=T, stranded=T, norm=F, scale=F, smo=F, heatmap=F, rm.outlier=F, genome="hg19", outPrefix=op)
#'
#' @export plot_start_end_feature
#'
#'

plot_start_end_feature <- function(queryfiles, inputfiles=NULL, txdb, featureName, CLIP_reads=F, binsize=10, fix_width=0,
                               longest=T, ext=c(-500, 200, -200, 500), hl=c(-50, 50, -50, 50), randomize=FALSE, stranded=T, norm=F, scale=F, smo=F, heatmap=F,
                               rm.outlier=F, genome="hg19", outPrefix="plots", useScore=FALSE, useSizeFactor=FALSE){

  querylabels <- names(queryfiles)
  names(querylabels) <- queryfiles
  inputlabels <- names(inputfiles)
  names(inputlabels) <- inputfiles
  if(!is.null(inputfiles)){
    if(length(queryfiles) != length(inputfiles)){
      stop("the number of queryfiles and inputfiles must be the same!")
    }
    names(inputlabels) <- inputfiles
    queryfiles <- c(queryfiles, inputfiles)
    querylabels <- c(querylabels, inputlabels)
  }

  feature <- rfeature <- NULL
  fs <- fe <- rfs <- rfe <- NULL

  feature <- gtf_to_bed_longest_tx(txdb, featureName, longest=longest)[["GRanges"]]
  minimal_width <- ext[2]-ext[3]
  feature <- feature[width(feature)>minimal_width]
  nf <- length(feature)

  print(paste("number of features: ", featureName, nf))

  fs <- promoters(resize(feature,width=1,fix="start"), upstream=-ext[1],downstream=ext[2])
  fe <- promoters(resize(feature,width=1,fix="end"), upstream=-ext[3],downstream=ext[4])

  if(randomize){
    #rfeature <- randomizeFeature(disjoin(feature))
    random_points <- sample(ext[1]:ext[4], length(feature), replace=T)
    rfeature <- shift(feature, shift=random_points)
    rfs <- promoters(resize(rfeature,width=1,fix="start"), upstream=-ext[1],downstream=ext[2])
    rfe <- promoters(resize(rfeature,width=1,fix="end"), upstream=-ext[3],downstream=ext[4])
  }

  ext[2] <- ext[2] - (ext[2]-ext[1])%%binsize ## to avoid binsize inconsistency, as the final binsize is dictated by bin_num
  bin_num_s <- round((ext[2]-ext[1])/binsize)
  ext[4] <- ext[4] - (ext[4]-ext[3])%%binsize ## to avoid binsize inconsistency, as the final binsize is dictated by bin_num
  bin_num_e <- round((ext[4]-ext[3])/binsize)

  mat_list <- NULL
  mat_list[[paste("Start of", featureName)]] <- list("window"=fs, "rwindow"=rfs, s=ext[1], e=ext[2], "xmin"=hl[1], "xmax"=hl[2], "bin_num"=bin_num_s)
  mat_list[[paste("End of", featureName)]] <- list("window"=fe, "rwindow"=rfe,  s=ext[3], e=ext[4], "xmin"=hl[3], "xmax"=hl[4], "bin_num"=bin_num_e)


  bedInputs <- handle_input(queryfiles, CLIP_reads=CLIP_reads, fix_width=fix_width,  useScore=useScore, outRle=TRUE, useSizeFactor=useSizeFactor, genome=genome)

  scoreMatrix_list <- list()
  scoreMatrix_list_random <- list()
  for(locus in names(mat_list)){
    windowR <- mat_list[[locus]]$window
    rwindowR <- mat_list[[locus]]$rwindow
    xmin <- mat_list[[locus]]$xmin
    xmax <- mat_list[[locus]]$xmax
    bin_num <- mat_list[[locus]]$bin_num
    start <- mat_list[[locus]]$s
    end <- mat_list[[locus]]$e

    for(queryfile in queryfiles){

      querylabel <- querylabels[queryfile]
      print(querylabel)
      queryRegions <- bedInputs[[queryfile]]$query
      libsize <- bedInputs[[queryfile]]$size

      bin_op <- "mean"
      weight_col <- bedInputs[[queryfile]]$weight

      fullmatrix <- parallel_scoreMatrixBin(queryRegions, windowR, bin_num, bin_op, weight_col, stranded)
      fullmatrix <- process_scoreMatrix(fullmatrix, libsize, norm, scale, heatmap, rm.outlier=(rm.outlier&&is.null(inputfiles)))

      scoreMatrix_list[[querylabel]][[locus]] <- fullmatrix

      if(randomize){

        rfullmatrix <- parallel_scoreMatrixBin(queryRegions, rwindowR, bin_num, bin_op, weight_col, stranded)
        rfullmatrix <- process_scoreMatrix(rfullmatrix, libsize, norm, scale, heatmap, rm.outlier=(rm.outlier&&is.null(inputfiles)))
        scoreMatrix_list_random[[querylabel]][[locus]] <- rfullmatrix
      }
    }
  }

  plot_df <- NULL
  for(locus in names(mat_list)){

    xmin <- mat_list[[locus]]$xmin
    xmax <- mat_list[[locus]]$xmax
    bin_num <- mat_list[[locus]]$bin_num
    start <- mat_list[[locus]]$s
    end <- mat_list[[locus]]$e

    for(queryfile in queryfiles){

      querylabel <- querylabels[queryfile]
      print(querylabel)

      fullmatrix <- scoreMatrix_list[[querylabel]][[locus]]

      colm <- apply(fullmatrix, 2, mean)
      colsd <- apply(fullmatrix, 2, sd)
      colse <- colsd/sqrt(nrow(fullmatrix))
      collabel <- seq(start, (end-binsize), binsize)
      querybed <- as.factor(rep(querylabel, ncol(fullmatrix)))
      location <- as.factor(rep(locus, ncol(fullmatrix)))
      levels(location) <- rev(levels(location))

      sub_df <- NULL
      sub_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=querybed, "Feature"=location)
      if(smo){
        sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df=as.integer(bin_num/5))$y)
        sub_df$se <- as.vector(smooth.spline(sub_df$se, df=as.integer(bin_num/5))$y)
      }
      sub_df <- mutate(sub_df, lower=Intensity-se, upper=Intensity+se)
      plot_df <- rbind(plot_df, sub_df)

      if(randomize){

        rfullmatrix <- scoreMatrix_list_random[[querylabel]][[locus]]

        rcolm <- apply(rfullmatrix, 2, mean)
        rcolsd <- apply(rfullmatrix, 2, sd)
        rcolse <- rcolsd/sqrt(nrow(rfullmatrix))
        rcollabel <- seq(start, (end-binsize), binsize)
        rquerybed <- as.factor(rep(paste0(querylabel,":Random"), ncol(rfullmatrix)))
        location <- as.factor(rep(locus, ncol(fullmatrix)))
        levels(location) <- rev(levels(location))

        rsub_df <- NULL
        rsub_df <- data.frame("Intensity"=rcolm, "sd"=rcolsd, "se"=rcolse, "Position"=rcollabel, "Query"=rquerybed, "Feature"=location)
        if(smo){
          rsub_df$Intensity <- as.vector(smooth.spline(rsub_df$Intensity, df=as.integer(bin_num/5))$y)
          rsub_df$se <- as.vector(smooth.spline(rsub_df$se, df=as.integer(bin_num/5))$y)
        }
        rsub_df <- mutate(rsub_df, lower=Intensity-se, upper=Intensity+se)

        plot_df <- rbind(plot_df, rsub_df)

      }
    }
  }

  if(!is.null(outPrefix)){
    while(!is.null(dev.list())){
      dev.off()
    }
    pdf(paste0(outPrefix, ".pdf"), width=10, height=8)
  }
  ## plot multi bed lines for one feature
  p <- ggplot(plot_df, aes(x=Position, y=Intensity, color=Query)) +
    geom_line(size=1) + #geom_point(color="grey30", size=2) +
    geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=1) +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=Query), linetype=0, alpha=0.3) +
    annotate("rect", xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill="grey", color="grey", alpha=0.3) +
    theme_classic() + theme(legend.position="top") + xlab("") + ylab("Signal intensity") +
    ggtitle(featureName) +
    facet_wrap(~Feature, scales="free_x")
  print(p)

  if(!is.null(inputfiles)){
    Ylab <- "Ratio-over-Input"
    logYlab <- expression(paste(Log[10], " Ratio-over-Input"))
    pseudo_count <- 0.1

    ratiofiles <- queryfiles[!queryfiles %in% inputfiles]
    ratiolabels <- querylabels[!querylabels %in% inputlabels]
    inputMatrix_list <- scoreMatrix_list[inputlabels]
    ratioMatrix_list <- scoreMatrix_list[ratiolabels]

    inputMatrix_list_random <- scoreMatrix_list_random[inputlabels]
    ratioMatrix_list_random <- scoreMatrix_list_random[ratiolabels]

    for(locus in names(mat_list)){
      if(0)locus <- names(mat_list)[1]
      for(i in seq_along(ratiolabels)){
        if(0) i <- 1
        rm <- ratioMatrix_list[[ratiolabels[i]]][[locus]]
        im <- inputMatrix_list[[inputlabels[i]]][[locus]]
        minrow <- min(nrow(rm), nrow(im))

        fullmatrix <- (rm[1:minrow,] + pseudo_count)/(im[1:minrow,] + pseudo_count)
        fullmatrix[is.na(fullmatrix)] <- 0

        if(rm.outlier){
          fullmatrix  <- rmOutlier(fullmatrix)
        }
        ratioMatrix_list[[ratiolabels[i]]][[locus]] <- fullmatrix

        ## for random feature
        if(randomize){
          rmr <- ratioMatrix_list_random[[ratiolabels[i]]][[locus]]
          imr <- inputMatrix_list_random[[inputlabels[i]]][[locus]]
          minrowr <- min(nrow(rmr), nrow(imr))

          fullmatrix <- (rmr[1:minrowr,] + pseudo_count)/(imr[1:minrowr,] + pseudo_count)
          fullmatrix[is.na(fullmatrix)] <- 0

          if(rm.outlier){
            fullmatrix  <- rmOutlier(fullmatrix)
          }
          ratioMatrix_list_random[[ratiolabels[i]]][[locus]] <- fullmatrix
        }
      }
    }


    plot_df <- NULL
    for(locus in names(mat_list)){

      xmin <- mat_list[[locus]]$xmin
      xmax <- mat_list[[locus]]$xmax
      bin_num <- mat_list[[locus]]$bin_num
      start <- mat_list[[locus]]$s
      end <- mat_list[[locus]]$e

      for(ratiofile in ratiofiles){

        ratiolabel <- ratiolabels[ratiofile]
        print(ratiolabel)

        fullmatrix <- ratioMatrix_list[[ratiolabel]][[locus]]

        colm <- apply(fullmatrix, 2, mean)
        colsd <- apply(fullmatrix, 2, sd)
        colse <- colsd/sqrt(nrow(fullmatrix))
        collabel <- seq(start, (end-binsize), binsize)
        ratiobed <- as.factor(rep(ratiolabel, ncol(fullmatrix)))
        location <- as.factor(rep(locus, ncol(fullmatrix)))
        levels(location) <- rev(levels(location))

        sub_df <- NULL
        sub_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=ratiobed, "Feature"=location)
        if(smo){
          sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df=as.integer(bin_num/5))$y)
          sub_df$se <- as.vector(smooth.spline(sub_df$se, df=as.integer(bin_num/5))$y)
        }
        sub_df <- mutate(sub_df, lower=Intensity-se, upper=Intensity+se)
        plot_df <- rbind(plot_df, sub_df)

        if(randomize){

          rfullmatrix <- ratioMatrix_list_random[[ratiolabel]][[locus]]

          rcolm <- apply(rfullmatrix, 2, mean)
          rcolsd <- apply(rfullmatrix, 2, sd)
          rcolse <- rcolsd/sqrt(nrow(rfullmatrix))
          rcollabel <- seq(start, (end-binsize), binsize)
          rratiobed <- as.factor(rep(paste0(ratiolabel,":Random"), ncol(rfullmatrix)))
          location <- as.factor(rep(locus, ncol(fullmatrix)))
          levels(location) <- rev(levels(location))

          rsub_df <- NULL
          rsub_df <- data.frame("Intensity"=rcolm, "sd"=rcolsd, "se"=rcolse, "Position"=rcollabel, "Query"=rratiobed, "Feature"=location)
          if(smo){
            rsub_df$Intensity <- as.vector(smooth.spline(rsub_df$Intensity, df=as.integer(bin_num/5))$y)
            rsub_df$se <- as.vector(smooth.spline(rsub_df$se, df=as.integer(bin_num/5))$y)
          }
          rsub_df <- mutate(rsub_df, lower=Intensity-se, upper=Intensity+se)

          plot_df <- rbind(plot_df, rsub_df)
        }
      }
    }

    ## plot multi bed lines for one feature
    p <- ggplot(plot_df, aes(x=Position, y=Intensity, color=Query)) +
      geom_line(size=1) + #geom_point(color="grey30", size=2) +
      geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=1) +
      geom_ribbon(aes(ymin=lower, ymax=upper, fill=Query), linetype=0, alpha=0.3) +
      annotate("rect", xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill="grey", color="grey", alpha=0.3) +
      theme_classic() + theme(legend.position="top") + xlab("") + ylab(Ylab) +
      ggtitle(featureName) +
      facet_wrap(~Feature, scales="free_x")
    print(p)
  }
  if(!is.null(outPrefix)){
    dev.off()
  }

  return(plot_df)
}

#' parallel_apply_scoreMatrixBin
#'
#' Function for parallel computation of scoreMatrixBin. The 'windows' parameter of the scoreMatrixBin method is provided as a list of
#' GRangesList objects, allowing simultaneous processing of multiple windows objects.
#'
#' @param windowRs, a list of GRangesList objects.
#' @param queryRegions, a RleList object or Granges object providing input for the 'target' parameter of the scoreMatrixBin method
#' @param bin_num, number of bins the windows should be divided into
#' @param bin_op, operation on the signals in a bin, a string in c("mean", "max", "min", "median", "sum") is accepted.
#' @param weight_col, if the queryRegions is a GRanges object, a numeric column in meta data part can be used as weights.
#' @param stranded, a boolean object to indicate if the strand of the windows should be considered to determine upstream and downstream
#'
#' @return a list of numeric matrices
#'
#' @example
#'
#' @export parallel_apply_scoreMatrixBin
#'
#'
parallel_apply_scoreMatrixBin <- function(windowRs, queryRegions, bin_num, bin_op, weight_col, stranded){
  nc <- length(windowRs)
  cl <- start_parallel(nc)
  clusterExport(cl, varlist=c("ScoreMatrixBin"), envir=.GlobalEnv)
  clusterExport(cl, varlist=c("queryRegions", "bin_op", "bin_num", "weight_col", "stranded"), envir=environment())
  call_scoreMatrixBin <- function(windowR){
    smc <- ScoreMatrixBin(target = queryRegions, windows = windowR, bin.num=bin_num, bin.op=bin_op, weight.col=weight_col, strand.aware = stranded)
    smc@.Data
  }

  smcList <- parLapply(cl, windowRs, call_scoreMatrixBin)
  stop_parallel(cl)

  smcList
}

#' parallel_scoreMatrixBin
#'
#' Function for parallel computation of scoreMatrixBin. The 'windows' parameter of the scoreMatrixBin method is split into 5 chunks,
#' and scoreMatrixBin is called on each chunk simultaneously to speed up the computation.
#'
#' @param windowRs, a GRangesList object.
#' @param queryRegions, a RleList object or Granges object providing input for the 'target' parameter of the scoreMatrixBin method
#' @param bin_num, number of bins the windows should be divided into
#' @param bin_op, operation on the signals in a bin, a string in c("mean", "max", "min", "median", "sum") is accepted.
#' @param weight_col, if the queryRegions is a GRanges object, a numeric column in meta data part can be used as weights.
#' @param stranded, a boolean object to indicate if the strand of the windows should be considered to determine upstream and downstream
#'
#' @return a numeric matrices
#'
#' @example
#'
#' @export parallel_scoreMatrixBin
#'
#'
parallel_scoreMatrixBin <- function(queryRegions, windowRs, bin_num, bin_op, weight_col, stranded){

  nc <- 5
  cl <- start_parallel(nc)
  clusterExport(cl, varlist=c("ScoreMatrixBin"), envir=.GlobalEnv)
  clusterExport(cl, varlist=c("queryRegions", "bin_num", "bin_op", "weight_col", "stranded"), envir=environment())
  call_scoreMatrixBin <- function(windowR){

    ScoreMatrixBin(target = queryRegions, windows = windowR, bin.num=bin_num, bin.op=bin_op, weight.col=weight_col, strand.aware = stranded)
  }

  wRs <- split(windowRs, factor(cut(seq_along(windowRs), breaks=nc)))
  smc <- parLapply(cl, wRs, call_scoreMatrixBin)
  stop_parallel(cl)

  sm <- lapply(smc, function(x)data.table(x@.Data))

  smdt <- rbindlist(sm)
  smdt
}

#' plot_3parts_metagene
#
#' Plot reads or peak signal intensity of samples in the query files around genes. The upstream and downstream windows flanking genes
#' can be given separately, the parameter 'meta' controls if gene or metagene plots are generated.
#'
#' @param queryfiles, a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param querylabels, a vector of character strings serving as short labels of the queryfiles, will be used as sample labels in the plots
#' @param txdb, a TxDb object defined in GenomicFeatures package
#' @param meta, a boolean object indicating whether a metagene (intron excluded) or gene (intron included) plot should be produced
#' @param inputfiles, a vector of input sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param inputlabels, a vector of character strings serving as short labels of the inputfiles, will be used as input labels in the plots
#' @param CLIP_reads, a boolean object to indicate if the bam reads should be shifted to the -1 position at the 5' of the reads
#' @param nbins, an integer defines the total number of bins
#' @param fiveP, extension out of the 5' boundary of gene
#' @param threeP, extension out of the 3' boundary of gene
#' @param fix_width, an integer defines how long should the reads should be extended to
#' @param norm, an boolean indicates if the intensities should be normalized by the sample library sizes
#' @param longest, a boolean object to indicate whether the output should be limited to the longest transcript of each gene
#' @param stranded, a boolean indicate whether the strand of the feature should be considered
#' @param scale, a boolean indicates whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smo, a boolean indicates whether the line should smoothed with a spline smoothing algorithm
#' @param heatmap, a boolean indicates whether a heatmap of the score matrix should be generated
#' @param rm.outlier, a boolean indicates whether a row with abnormally high values in the score matrix should be removed
#' @param genome, genome of the features, the program is mostly concerned with human, if non-human genome is used, certain features may not work
#' @param outPrefix, a character string specifying output file prefix for plots (outPrefix.pdf)
#'
#' @return a dataframe containing the data used for plotting
#'
#' @examples
#'
#'
#' @export plot_3parts_metagene


plot_3parts_metagene <- function(queryfiles, txdb, meta=T, inputfiles=NULL, nbins=100, norm=FALSE, longest=TRUE, scale=FALSE,
                                 fix_width=0, CLIP_reads=FALSE, fiveP=1000, threeP=1000, smo=FALSE, stranded=TRUE, outPrefix="plots", genome="hg19", useScore=FALSE,
                                 heatmap=FALSE, rm.outlier=FALSE, useSizeFactor=FALSE){


  querylabels <- names(queryfiles)
  names(querylabels) <- queryfiles
  inputlabels <- names(inputfiles)
  names(inputlabels) <- inputfiles

  if(!is.null(inputfiles)){
    if(length(queryfiles) != length(inputfiles)){
      stop("the number of queryfiles and inputfiles must be the same!")
    }
    names(inputlabels) <- inputfiles
    queryfiles <- c(queryfiles, inputfiles)
    querylabels <- c(querylabels, inputlabels)
  }

  featureName <- ifelse(meta, "transcript", "gene")
  ## prepare transcripts that are suitable for overlap
  featureNames <-  c("promoter", featureName, "TTS")

  tl <- transcriptLengths(txdb, with.cds_len = T)
  ## selecte the longest transcript for each gene
  if(longest){
    longest_tx <- extract_longest_tx(txdb)
    tl <- tl %>% filter(tx_id %in% longest_tx$tx_id)
  }

  if(meta){
    means <- c(promoter=fiveP, median(tl$tx_len), TTS=threeP)
    scaled_bins <- round(means*nbins/sum(means))

    tl_selected <- tl %>% filter(tx_len > scaled_bins[2])

    feature <- exonsBy(txdb, by="tx")
    feature <- feature[tl_selected$tx_id]

    transcript <- transcripts(txdb)
    transcript <- transcript[transcript$tx_id %in% tl_selected$tx_id, ]

    promoter <- flank(transcript, width=fiveP, both=F, start=T, ignore.strand=FALSE)
    TTS <- flank(transcript, width=threeP, both=F, start=F, ignore.strand=FALSE)

    windowRs <- list(split(promoter, as.factor(promoter$tx_id)),
                     feature,
                     split(TTS, as.factor(TTS$tx_id)))

  }else{
    feature <- genes(txdb)
    feature <- feature[feature$gene_id %in% tl$gene_id]
    means <- c(promoter=fiveP, median(width(feature)), TTS=threeP)
    scaled_bins <- round(means*nbins/sum(means))
    feature <- feature[width(feature) > scaled_bins[2]]

    promoter <- flank(feature, width=fiveP, both=F, start=T, ignore.strand=FALSE)
    TTS <- flank(feature, width=threeP, both=F, start=F, ignore.strand=FALSE)

    windowRs <- list(split(promoter, as.factor(promoter$gene_id)),
                     split(feature, as.factor(feature$gene_id)),
                     split(TTS, as.factor(TTS$gene_id)))
  }



  names(windowRs) <- featureNames
  names(scaled_bins) <- featureNames

  print("median sizes for features")
  print(means)
  print("bin sizes for features")
  print(scaled_bins)
  print("Number of features included in analysis")
  print(unlist(lapply(windowRs, length)))
  ## start overlapping

  Ylab <- "Signal intensity"
  logYlab <- expression(paste(Log[10], " Signal intensity"))


  color_store <- c("#00AFBB", "#E7B800", "#A0BDE0", "#0020C2", "#64E986", "#F5DEB3", "#C19A6B", "#E8A317", "#8E7618", "#A0522D", "#990012", "#CB6D51")

  print("computing coverage for queryfiles")
  scoreMatrix_list <- list()

  bedInputs <- handle_input(queryfiles, CLIP_reads=CLIP_reads, fix_width=fix_width, useScore=useScore, outRle=TRUE, useSizeFactor=useSizeFactor, genome=genome)

  for(queryfile in queryfiles){

    if(0)queryfile <- queryfiles[1]
    querylabel <- querylabels[queryfile]

    myInput <- bedInputs[[queryfile]]
    libsize <- myInput$size
    queryRegions <- myInput$query
    fileType <- myInput$type
    weight_col <- myInput$weight

    for(w in featureNames){
      #w <- "utr5"
      print(w)
      windowR <- windowRs[[w]]
      bin_num <- scaled_bins[w]

      bin_op <- "mean"
      fullmatrix <- parallel_scoreMatrixBin(queryRegions, windowR, bin_num, bin_op, weight_col, stranded)

      fullmatrix <- process_scoreMatrix(fullmatrix, libsize=libsize, norm=norm, scale=scale, heatmap=heatmap, rm.outlier=(rm.outlier&&is.null(inputfiles)))
      scoreMatrix_list[[querylabel]][[w]] <- fullmatrix
    }
  }

  print("plotting coverage for queryfiles")
  mplot_df <- NULL
  vx <- c(1, cumsum(scaled_bins[1:(length(scaled_bins)-1)])+1) ## x axis points for vlines that demarcate the genomic features
  names(vx) <- featureNames

  for(queryfile in queryfiles){

    if(0)queryfile <- queryfiles[1]
    querylabel <- querylabels[queryfile]

    plot_df <- NULL

    for(w in featureNames){
      #w <- "utr5"
      print(w)

      bin_num <- scaled_bins[w]

      fullmatrix <- scoreMatrix_list[[querylabel]][[w]]

      colm <- apply(fullmatrix, 2, mean)
      colsd <- apply(fullmatrix, 2, sd)
      colse <- colsd/sqrt(nrow(fullmatrix))
      collabel <- seq(vx[w], vx[w]+bin_num-1)
      querybed <- rep(querylabel, bin_num)
      featuretype <- rep(w, bin_num)

      sub_df <- NULL
      sub_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=querybed, "Feature"=featuretype)
      plot_df <- rbind(plot_df, sub_df)
    }

    if(smo){
      plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df=as.integer(nbins/5))$y)
      plot_df$se <- as.vector(smooth.spline(plot_df$se, df=as.integer(nbins/5))$y)
    }

    mplot_df <- rbind(mplot_df, plot_df)

  }

  mplot_df <- mutate(mplot_df, lower=Intensity-se, upper=Intensity+se)


  if(!is.null(outPrefix)){
    pdf(paste(outPrefix, "pdf", sep="."), height=8, width=12)
  }


  ## plot multi-sample lines with error band
  p <- ggplot(mplot_df, aes(x=Position, y=Intensity, color=Query)) + scale_fill_manual(values=color_store[1:length(queryfiles)]) +
    geom_line(size=1) + #geom_point(color="grey30", size=2) +
    geom_vline(xintercept = vx[2:length(vx)], linetype="dotted", color = "blue", size=0.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=Query), linetype=0, alpha=0.3) +
    theme_classic() +
    theme(legend.position="top",
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank()) +
    ylab(Ylab)


  values <- data.frame(id=featureNames, value=c(1.25, 1.75, 1.25))
  xmax <- max(mplot_df$Position)
  positions <- data.frame(
    id = rep(featureNames, each = 4),
    x = c(vx[2], vx[1], vx[1], vx[2], vx[3], vx[2], vx[2], vx[3], xmax, vx[3], vx[3], xmax),
    y = c(3, 3, 4, 4, 2.5, 2.5, 4.5, 4.5, 3, 3, 4, 4) - 2
  )


  datapoly <- merge(values, positions, by = c("id"))

  pp <- ggplot(datapoly, aes(x = x, y = y)) +
    geom_polygon(aes(fill = value, group = id)) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())

  annotx <- scaled_bins/2
  for(i in 2:length(scaled_bins)){
    annotx[i] <- annotx[i] + sum(scaled_bins[1:(i-1)])
  }

  five <- fiveP/1000
  five <- paste0(five, "kb")
  three <- threeP/1000
  three <- paste0(three, "kb")
  annot <- data.frame(
    fn = c(five, featureName, three),
    x = annotx,
    y = 0
  )
  ppp <- ggplot(annot, aes(x = x, y = y, label=fn)) +
    geom_text(size=4) +
    coord_cartesian(xlim = c(1, xmax)) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())


  outp <- plot_grid(p, pp, ppp, ncol = 1, align = 'v', axis="b", rel_heights = c(20,1,1))
  print(outp)

  if(!is.null(inputfiles)){
    print("computing coverage for ratio over input")
    Ylab <- "Ratio-over-Input"
    logYlab <- expression(paste(Log[10], " Ratio-over-Input"))
    pseudo_count <- 0.1

    inputMatrix_list <- scoreMatrix_list[inputlabels]
    ratiofiles <- queryfiles[!queryfiles %in% inputfiles]
    ratiolabels <- querylabels[!querylabels %in% inputlabels]
    ratioMatrix_list <- scoreMatrix_list[ratiolabels]
    for(w in featureNames){
      for(i in seq_along(querylabels)){
        fullmatrix <- (ratioMatrix_list[[ratiolabels[i]]][[w]] + pseudo_count)/(inputMatrix_list[[inputlabels[i]]][[w]] + pseudo_count)
        fullmatrix[is.na(fullmatrix)] <- 0

        if(rm.outlier){
          fullmatrix <- rmOutlier(fullmatrix)
        }

        ratioMatrix_list[[ratiolabels[i]]][[w]] <- fullmatrix
      }
    }

    print("plotting coverage for ratio over input")
    mplot_df <- NULL
    for(ratiofile in ratiofiles){

      if(0)ratiofile <- ratiofiles[1]
      ratiolabel <- ratiolabels[ratiofile]

      plot_df <- NULL

      for(w in featureNames){
        #w <- "utr5"
        print(w)

        bin_num <- scaled_bins[w]

        fullmatrix <- scoreMatrix_list[[ratiolabel]][[w]]

        colm <- apply(fullmatrix, 2, mean)
        colsd <- apply(fullmatrix, 2, sd)
        colse <- colsd/sqrt(nrow(fullmatrix))
        collabel <- seq(vx[w], vx[w]+bin_num-1)
        querybed <- rep(ratiolabel, bin_num)
        featuretype <- rep(w, bin_num)

        sub_df <- NULL
        sub_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=querybed, "Feature"=featuretype)
        plot_df <- rbind(plot_df, sub_df)
      }

      if(smo){
        plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df=as.integer(nbins/5))$y)
        plot_df$se <- as.vector(smooth.spline(plot_df$se, df=as.integer(nbins/5))$y)
      }

      mplot_df <- rbind(mplot_df, plot_df)

    }

    mplot_df <- mutate(mplot_df, lower=Intensity-se, upper=Intensity+se)

    ## plot multi-sample lines with error band
    p <- ggplot(mplot_df, aes(x=Position, y=Intensity, color=Query)) + scale_fill_manual(values=color_store[1:length(ratiofiles)]) +
      geom_line(size=1) + #geom_point(color="grey30", size=2) +
      geom_vline(xintercept = vx[2:length(vx)], linetype="dotted", color = "blue", size=0.5) +
      geom_ribbon(aes(ymin=lower, ymax=upper, fill=Query), linetype=0, alpha=0.3) +
      theme_classic() +
      theme(legend.position="top",
            axis.title.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.x=element_blank()) +
      ylab(Ylab)

    outp <- plot_grid(p, pp, ppp, ncol = 1, align = 'v', axis="b", rel_heights = c(20,1,1))
    print(outp)

  }

  if(!is.null(outPrefix)){
    dev.off()
  }

  return(mplot_df)
}

#' plot_reference_region
#
#' Plot reads or peak signal intensity of samples in the query files inside regions defined in the centerfiles. The upstream and downstream windows flanking genes
#' can be given separately, the parameter 'meta' controls if gene or metagene plots are generated.
#'
#' @param queryfiles, a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param querylabels, a vector of character strings serving as short labels of the queryfiles, will be used as sample labels in the plots
#' @param centerfiles, a vector of reference file names. The file should be .bed format only
#' @param centerlabels, a vector of character strings serving as short labels of the centerfiles, will be used as sample labels in the plots
#' @param inputfiles, a vector of input sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param inputlabels, a vector of character strings serving as short labels of the inputfiles, will be used as input labels in the plots
#' @param CLIP_reads, a boolean object to indicate if the bam reads should be shifted to the -1 position at the 5' of the reads
#' @param nbins, an integer defines the total number of bins
#' @param fiveP, extension out of the 5' boundary of gene
#' @param threeP, extension out of the 3' boundary of gene
#' @param fix_width, an integer defines how long should the reads should be extended to
#' @param norm, an boolean indicates if the intensities should be normalized by the sample library sizes
#' @param longest, a boolean object to indicate whether the output should be limited to the longest transcript of each gene
#' @param stranded, a boolean indicate whether the strand of the feature should be considered
#' @param scale, a boolean indicates whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smo, a boolean indicates whether the line should smoothed with a spline smoothing algorithm
#' @param heatmap, a boolean indicates whether a heatmap of the score matrix should be generated
#' @param rm.outlier, a boolean indicates whether a row with abnormally high values in the score matrix should be removed
#' @param genome, genome of the features, the program is mostly concerned with human, if non-human genome is used, certain features may not work
#' @param outPrefix, a character string specifying output file prefix for plots (outPrefix.pdf)
#'
#' @return a dataframe containing the data used for plotting
#'
#' @examples
#'
#'
#' @export plot_reference_region

plot_reference_region <- function(queryfiles, centerfiles, inputfiles=NULL, nbins=100, norm=FALSE, useSizeFactor=FALSE,
                                  fix_width=0, CLIP_reads=FALSE, scale=FALSE, heatmap=FALSE, longest=TRUE, useScore=FALSE,
                                  fiveP=1000, threeP=1000, smo=FALSE, stranded=TRUE, outPrefix="plots", genome="hg19", rm.outlier){

  if(!is.null(outPrefix)){
    pdf(paste(outPrefix, "pdf", sep="."), height=8, width=12)
  }

  querylabels <- names(queryfiles)
  names(querylabels) <- queryfiles
  inputlabels <- names(inputfiles)
  names(inputlabels) <- inputfiles
  centerlabels <- names(centerfiles)
  names(centerlabels) <- centerfiles

  if(!is.null(inputfiles)){
    if(length(queryfiles) != length(inputfiles)){
      stop("the number of queryfiles and inputfiles must be the same!")
    }
    names(inputlabels) <- inputfiles
    queryfiles <- c(queryfiles, inputfiles)
    querylabels <- c(querylabels, inputlabels)
  }

  ## prepare transcripts that are suitable for overlap
  featureNames <-  c("upstream", "region", "downstream")

  means <- c(upstream=fiveP, region=fiveP+threeP, downstream=threeP)
  scaled_bins <- round(means*nbins/sum(means))
  names(scaled_bins) <- featureNames

  bedInputs <- handle_input(queryfiles, CLIP_reads=CLIP_reads, fix_width=fix_width, useScore=useScore, outRle=TRUE, useSizeFactor=useSizeFactor, genome=genome)
  centerInputs <- handle_input(centerfiles, CLIP_reads=FALSE, fix_width=0, useScore=FALSE, outRle=FALSE, useSizeFactor=FALSE, genome=genome)

  print("computing coverage for Sample")
  scoreMatrix_list <- list()

  for(queryfile in queryfiles){
    if(0)queryfile <- queryfiles[1]
    querylabel <- querylabels[queryfile]
    print(querylabel)

    myInput <- bedInputs[[queryfile]]
    libsize <- myInput$size
    queryRegions <- myInput$query
    fileType <- myInput$type
    weight_col <- myInput$weight

    for(centerfile in centerfiles){
      if(0)centerfile <- centerfiles[1]
      centerlabel <- centerlabels[centerfile]
      print(centerlabel)
      centerInput <- centerInputs[[centerfile]]
      centerGr <- centerInput$query

      centerGr <- centerGr[width(centerGr) >= scaled_bins["region"]]

      upstreamGr <- flank(centerGr, width=fiveP, start=T, both=F, use.names=T, ignore.strand=FALSE)
      downstreamGr <- flank(centerGr, width=threeP, start=F, both=F, use.names=T, ignore.strand=FALSE)
      windowRegions <- split(centerGr, f=factor(seq(1:length(centerGr))))
      windowUp <- split(upstreamGr, f=factor(seq(1:length(upstreamGr))))
      windowDown <- split(downstreamGr, f=factor(seq(1:length(downstreamGr))))

      windowRs <- list(windowUp, windowRegions, windowDown)
      names(windowRs) <- featureNames

      for(w in featureNames){
        #w <- "utr5"
        print(w)
        windowR <- windowRs[[w]]
        bin_num <- scaled_bins[w]
        bin_op <- "mean"
        fullmatrix <- parallel_scoreMatrixBin(queryRegions, windowR, bin_num, bin_op, weight_col, stranded)
        fullmatrix <- process_scoreMatrix(fullmatrix, libsize=libsize, norm=norm, scale=scale, heatmap=heatmap, rm.outlier=(rm.outlier&&is.null(inputfiles)))
        scoreMatrix_list[[querylabel]][[centerlabel]][[w]] <- fullmatrix
      }
    }
  }



  mplot_df <- NULL
  vx <- c(1, cumsum(scaled_bins[1:(length(scaled_bins)-1)])+1) ## x axis points for vlines that demarcate the genomic features
  names(vx) <- featureNames
  color_store <- c("#00AFBB", "#E7B800", "#A0BDE0", "#0020C2", "#64E986", "#F5DEB3", "#C19A6B", "#E8A317", "#8E7618", "#A0522D", "#990012", "#CB6D51")

  print("plotting coverage profiles")
  for(queryfile in queryfiles){
    if(0)queryfile <- queryfiles[1]
    querylabel <- querylabels[queryfile]
    print(querylabel)
    for(centerfile in centerfiles){
      if(0)centerfile <- centerfiles[1]
      centerlabel <- centerlabels[centerfile]
      print(centerlabel)
      plot_df <- NULL
      for(w in featureNames){
        print(w)
        bin_num <- scaled_bins[w]
        fullmatrix <- scoreMatrix_list[[querylabel]][[centerlabel]][[w]]

        colm <- apply(fullmatrix, 2, mean)
        colsd <- apply(fullmatrix, 2, sd)
        colse <- colsd/sqrt(nrow(fullmatrix))
        collabel <- seq(vx[w], vx[w]+bin_num-1)
        querybed <- rep(querylabel, bin_num)
        centerbed <- rep(centerlabel, bin_num)
        featuretype <- rep(w, bin_num)

        sub_df <- NULL
        sub_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=querybed, "Reference"=centerbed, "Feature"=featuretype)
        plot_df <- rbind(plot_df, sub_df)
      }

      if(smo){
        plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df=as.integer(nbins/5))$y)
        plot_df$se <- as.vector(smooth.spline(plot_df$se, df=as.integer(nbins/5))$y)
      }

      mplot_df <- rbind(mplot_df, plot_df)
    }
  }

  mplot_df <- mplot_df %>%
    mutate(lower=Intensity-se, upper=Intensity+se)


  values <- data.frame(id=featureNames, value=c(1.25, 1.75, 1.25))
  xmax <- max(mplot_df$Position)
  positions <- data.frame(
    id = rep(featureNames, each = 4),
    x = c(vx[2], vx[1], vx[1], vx[2], vx[3], vx[2], vx[2], vx[3], xmax, vx[3], vx[3], xmax),
    y = c(3, 3, 4, 4, 2.5, 2.5, 4.5, 4.5, 3, 3, 4, 4) - 2
  )


  datapoly <- merge(values, positions, by = c("id"))

  pp <- ggplot(datapoly, aes(x = x, y = y)) +
    geom_polygon(aes(fill = value, group = id)) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())

  annotx <- scaled_bins/2
  for(i in 2:length(scaled_bins)){
    annotx[i] <- annotx[i] + sum(scaled_bins[1:(i-1)])
  }

  five <- fiveP/1000
  five <- paste0(five, "kb")
  three <- threeP/1000
  three <- paste0(three, "kb")
  annot <- data.frame(
    fn = c(five, "region", three),
    x = annotx,
    y = 0
  )
  ppp <- ggplot(annot, aes(x = x, y = y, label=fn)) +
    geom_text(size=4) +
    coord_cartesian(xlim = c(1, xmax)) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())


  for(i in seq_along(querylabels)){
    for(beds in combn(querylabels, i, simplify = F)){
      for(j in seq_along(centerlabels)){
        for(centers in combn(centerlabels, j, simplify = F)){
          print(beds)
          print(centers)

          aplot_df <- mplot_df %>%
            filter(Query %in% beds & Reference %in% centers) %>%
            mutate(Group=paste(Query,Reference,sep=":"), .keep="all")

          ## plot multi-sample lines with error band
          p <- ggplot(aplot_df, aes(x=Position, y=Intensity, color=Group)) + scale_fill_manual(values=color_store[1:length(queryfiles)]) +
            geom_line(size=1) + #geom_point(color="grey30", size=2) +
            geom_vline(xintercept = vx[2:length(vx)], linetype="dotted", color = "blue", size=0.5) +
            geom_ribbon(aes(ymin=lower, ymax=upper, fill=Group), linetype=0, alpha=0.3) +
            theme_classic() +
            theme(legend.position="top",
                  axis.title.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.text.x=element_blank()) +
            ylab("Average signal +/- SE")

          outp <- plot_grid(p, pp, ppp, ncol = 1, align = 'v', axis="b", rel_heights = c(20,1,1))
          print(outp)
        }
      }
    }
  }


  if(!is.null(inputfiles)){
    Ylab <- "Ratio-over-Input"
    logYlab <- expression(paste(Log[10], " Ratio-over-Input"))
    pseudo_count <- 0.1

    ratiofiles <- queryfiles[!queryfiles %in% inputfiles]
    ratiolabels <- querylabels[!querylabels %in% inputlabels]
    inputMatrix_list <- scoreMatrix_list[inputlabels]
    ratioMatrix_list <- scoreMatrix_list[ratiolabels]

    for(centerfile in centerfiles){
      centerlabel <- centerlabels[centerfile]
      for(regionName in featureNames){
        for(i in seq_along(ratiolabels)){
          rm <- ratioMatrix_list[[ratiolabels[i]]][[centerlabel]][[regionName]]
          im <- inputMatrix_list[[inputlabels[i]]][[centerlabel]][[regionName]]
          minrow <- min(nrow(rm), nrow(im))

          fullmatrix <- (rm[1:minrow,] + pseudo_count)/(im[1:minrow,] + pseudo_count)
          fullmatrix[is.na(fullmatrix)] <- 0

          if(rm.outlier){
            fullmatrix <- rmOutlier(fullmatrix)
          }

          ## scaled to the range of 0:1 to allow comparison with random ratio
          smc <- t(apply(fullmatrix, 1, scales::rescale)) ## rescale to 0:1 range
          fullmatrix <- as.data.frame(smc)

          ratioMatrix_list[[ratiolabels[i]]][[centerlabel]][[regionName]] <- fullmatrix
        }
      }
    }

    mplot_df <- NULLprint("plotting coverage profiles")
    for(ratiofile in ratiofiles){
      if(0)ratiofile <- ratiofiles[1]
      ratiolabel <- ratiolabels[ratiofile]
      print(ratiolabel)
      for(centerfile in centerfiles){
        if(0)centerfile <- centerfiles[1]
        centerlabel <- centerlabels[centerfile]
        print(centerlabel)
        plot_df <- NULL
        for(w in featureNames){
          print(w)
          bin_num <- scaled_bins[w]
          fullmatrix <- scoreMatrix_list[[ratiolabel]][[centerlabel]][[w]]

          colm <- apply(fullmatrix, 2, mean)
          colsd <- apply(fullmatrix, 2, sd)
          colse <- colsd/sqrt(nrow(fullmatrix))
          collabel <- seq(vx[w], vx[w]+bin_num-1)
          ratiobed <- rep(ratiolabel, bin_num)
          centerbed <- rep(centerlabel, bin_num)
          featuretype <- rep(w, bin_num)

          sub_df <- NULL
          sub_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "ratio"=ratiobed, "Reference"=centerbed, "Feature"=featuretype)
          plot_df <- rbind(plot_df, sub_df)
        }

        if(smo){
          plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df=as.integer(nbins/5))$y)
          plot_df$se <- as.vector(smooth.spline(plot_df$se, df=as.integer(nbins/5))$y)
        }

        mplot_df <- rbind(mplot_df, plot_df)
      }
    }

    mplot_df <- mplot_df %>%
      mutate(lower=Intensity-se, upper=Intensity+se)

    for(i in seq_along(ratiolabels)){
      for(beds in combn(ratiolabels, i, simplify = F)){
        for(j in seq_along(centerlabels)){
          for(centers in combn(centerlabels, j, simplify = F)){
            print(beds)
            print(centers)

            aplot_df <- mplot_df %>%
              filter(Query %in% beds & Reference %in% centers) %>%
              mutate(Group=paste(Query,Reference,sep=":"), .keep="all")

            ## plot multi-sample lines with error band
            p <- ggplot(aplot_df, aes(x=Position, y=Intensity, color=Group)) + scale_fill_manual(values=color_store[1:length(ratiofiles)]) +
              geom_line(size=1) + #geom_point(color="grey30", size=2) +
              geom_vline(xintercept = vx[2:length(vx)], linetype="dotted", color = "blue", size=0.5) +
              geom_ribbon(aes(ymin=lower, ymax=upper, fill=Group), linetype=0, alpha=0.3) +
              theme_classic() +
              theme(legend.position="top",
                    axis.title.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.text.x=element_blank()) +
              ylab(Ylab)

            outp <- plot_grid(p, pp, ppp, ncol = 1, align = 'v', axis="b", rel_heights = c(20,1,1))
            print(outp)
          }
        }
      }
    }

  }

  if(!is.null(outPrefix)){
    dev.off()
  }

  return(mplot_df)
}

#' plot_5parts_metagene
#
#' Plot reads or peak signal intensity of samples in the query files around genes. The upstream and downstream windows flanking genes
#' can be given separately, metagene plots are generated with 5'UTR, CDS and 3'UTR segments. The length of each segments are prorated
#' according to the median length of each segments.
#'
#' @param queryfiles, a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param querylabels, a vector of character strings serving as short labels of the queryfiles, will be used as sample labels in the plots
#' @param txdb, a TxDb object defined in GenomicFeatures package
#' @param meta, a boolean object indicating whether metagene (intron excluded) or gene (intron included) plots should be produced
#' @param inputfiles, a vector of input sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param inputlabels, a vector of character strings serving as short labels of the inputfiles, will be used as input labels in the plots
#' @param CLIP_reads, a boolean object to indicate if the bam reads should be shifted to the -1 position at the 5' of the reads
#' @param nbins, an integer defines the total number of bins
#' @param fiveP, extension out of the 5' boundary of gene
#' @param threeP, extension out of the 3' boundary of gene
#' @param fix_width, an integer defines how long should the reads should be extended to
#' @param norm, an boolean indicates if the intensities should be normalized by the sample library sizes
#' @param longest, a boolean object to indicate whether the output should be limited to the longest transcript of each gene
#' @param stranded, a boolean indicate whether the strand of the feature should be considered
#' @param scale, a boolean indicates whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smo, a boolean indicates whether the line should smoothed with a spline smoothing algorithm
#' @param heatmap, a boolean indicates whether a heatmap of the score matrix should be generated
#' @param rm.outlier, a boolean indicates whether a row with abnormally high values in the score matrix should be removed
#' @param genome, genome of the features, the program is mostly concerned with human, if non-human genome is used, certain features may not work
#' @param outPrefix, a character string specifying output file prefix for plots (outPrefix.pdf)
#' @param useIntron, a boolean indicates whether intron instead of TTS should be plotted
#' @param useScore, a boolean indicates whether to use score column of the bed file for computation of signal intensity
#'
#' @return a dataframe containing the data used for plotting
#'
#' @examples
#'
#'
#' @export plot_5parts_metagene


plot_5parts_metagene <- function(queryfiles, gFeatures, inputfiles=NULL, norm=FALSE,
                                 fix_width=0, CLIP_reads=FALSE, verbose=FALSE,
                                 smo=FALSE, scale=FALSE, stranded=TRUE, outPrefix=NULL, genome="hg19",
                                 heatmap=F, rm.outlier=F, useScore=FALSE, useSizeFactor=FALSE){

  querylabels <- names(queryfiles)
  names(querylabels) <- queryfiles
  inputlabels <- names(inputfiles)
  names(inputlabels) <- inputfiles

  if(!is.null(inputfiles)){
    if(length(queryfiles) != length(inputfiles)){
      stop("the number of queryfiles and inputfiles must be the same!")
    }
    names(inputlabels) <- inputfiles
    queryfiles <- c(queryfiles, inputfiles)
    querylabels <- c(querylabels, inputlabels)
  }

  windowRs <- gFeatures$windowRs
  featureNames <- names(windowRs)
  print(featureNames)
  print(lapply(windowRs, length))

  nbins <- gFeatures$nbins
  scaled_bins  <- gFeatures$scaled_bins
  meta <- gFeatures$meta
  fiveP <- gFeatures$fiveP
  threeP <- gFeatures$threeP
  useIntron <- gFeatures$useIntron

  ## start overlapping

  scoreMatrix_list <- list()

  bedInputs <- handle_input(inputFiles=queryfiles, CLIP_reads=CLIP_reads, fix_width=fix_width, useScore=useScore, outRle=TRUE, useSizeFactor=useSizeFactor, genome=genome)

  print("computing coverage for query files")

  for(queryfile in queryfiles){
    querylabel <- querylabels[queryfile]
    print(queryfile)
    bedInput <- bedInputs[[queryfile]]
    libsize <- bedInput$size
    queryRegions <- bedInput$query
    fileType <- bedInput$type
    weight_col <- bedInput$weight

    for(w in featureNames){
      #w <- "utr5"
      print(w)
      windowR <- windowRs[[w]]
      bin_num <- scaled_bins[w]

      bin_op <- "mean"
      if(bin_num > 0){
         fullmatrix <- parallel_scoreMatrixBin(queryRegions, windowR, bin_num, bin_op, weight_col, stranded)
         fullmatrix <- process_scoreMatrix(fullmatrix, libsize=libsize, norm=norm, scale=scale, heatmap=heatmap, rm.outlier=(rm.outlier&&is.null(inputfiles)))

         scoreMatrix_list[[querylabel]][[w]] <- fullmatrix
      }else{
         scoreMatrix_list[[querylabel]][[w]] <- NULL
      }

    }
  }



  print("Preparing data for individual plotting")

  vx <- c(1, cumsum(scaled_bins[1:(length(scaled_bins)-1)])+1) ## x axis points for vlines that demarcate the genomic features
  names(vx) <- featureNames

  mplot_df <- NULL
  for(querylabel in querylabels){
    plot_df <- NULL

    for(w in featureNames){
      #w <- "utr5"
      print(w)
      if(scaled_bins[w] > 0){
         bin_num <- scaled_bins[w]
         fullmatrix <- scoreMatrix_list[[querylabel]][[w]]

         fullmatrix[is.na(fullmatrix)] <- 0

         print(dim(fullmatrix))

         colm <- apply(fullmatrix, 2, mean)
         colsd <- apply(fullmatrix, 2, sd)
         colse <- colsd/sqrt(nrow(fullmatrix))
         collabel <- seq(vx[w], vx[w]+bin_num-1)
         querybed <- rep(querylabel, bin_num)
         featuretype <- rep(w, bin_num)

         sub_df <- NULL
         sub_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=querybed, "Feature"=featuretype)
         plot_df <- rbind(plot_df, sub_df)
      }

    }

    if(smo){
      plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df=as.integer(nbins/5))$y)
      plot_df$se <- as.vector(smooth.spline(plot_df$se, df=as.integer(nbins/5))$y)
    }

    mplot_df <- rbind(mplot_df, plot_df)

  }

  mplot_df <- mutate(mplot_df, lower=Intensity-se, upper=Intensity+se)

  if(!is.null(outPrefix)){
     print("Start plotting")
     values <- data.frame(id=featureNames, value=c(1.75, 1.5, 1.25, 1.5, 1.75))
     xmax <- max(mplot_df$Position)
     positions <- data.frame(
       id = rep(featureNames, each = 4),
       x = c(vx[2], vx[1], vx[1], vx[2], vx[3], vx[2], vx[2], vx[3], vx[4], vx[3], vx[3], vx[4], vx[5], vx[4], vx[4], vx[5], xmax, vx[5], vx[5], xmax),
       y = c(3, 3, 4, 4, 2.5, 2.5, 4.5, 4.5, 2, 2, 5, 5, 2.5, 2.5, 4.5, 4.5, 3, 3, 4, 4) - 2
     )


     datapoly <- merge(values, positions, by = c("id"))

     pp <- ggplot(datapoly, aes(x = x, y = y)) +
       geom_polygon(aes(fill = value, group = id)) +
       theme(axis.line=element_blank(),
             axis.text.x=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks=element_blank(),
             axis.title.x=element_blank(),
             axis.title.y=element_blank(),
             legend.position="none",
             panel.background=element_blank(),
             panel.border=element_blank(),
             panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             plot.background=element_blank())

     annotx <- scaled_bins/2
     for(i in 2:length(scaled_bins)){
       annotx[i] <- annotx[i] + sum(scaled_bins[1:(i-1)])
     }

     fiveP <- fiveP/1000
     five <- paste0(fiveP, "kb")
     if(fiveP==0) five=""
     threeP <- threeP/1000
     three <- paste0(threeP, "kb")
     if(threeP==0) three=""
     #if(!meta) five <- three <- ""
     if(useIntron){
       three <- "Intron"
     }
     annot <- data.frame(
       fn = c(five, "5'UTR", "CDS", "3'UTR", three),
       x = annotx,
       y = 0
     )
     ppp <- ggplot(annot, aes(x = x, y = y, label=fn)) +
       geom_text(size=4) +
       coord_cartesian(xlim = c(1, xmax)) +
       theme(axis.line=element_blank(),
             axis.text.x=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks=element_blank(),
             axis.title.x=element_blank(),
             axis.title.y=element_blank(),
             legend.position="none",
             panel.background=element_blank(),
             panel.border=element_blank(),
             panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             plot.background=element_blank())

     color_store <- c("#00AFBB", "#E7B800", "#A0BDE0", "#0020C2", "#64E986", "#F5DEB3", "#C19A6B", "#E8A317", "#8E7618", "#A0522D", "#990012", "#CB6D51")


    pdf(paste(outPrefix, "pdf", sep="."), height=8, width=12)

     for(i in seq_along(querylabels)){
       for(beds in combn(querylabels, i, simplify = F)){
         print(beds)
         if(verbose || i == length(querylabels)){
            aplot_df <- mplot_df %>%
                  filter(Query %in% beds)
            ## plot multi-sample lines with error band
            p <- ggplot(aplot_df, aes(x=Position, y=Intensity, color=Query)) + scale_fill_manual(values=color_store[1:length(beds)]) +
              geom_line(size=1) + #geom_point(color="grey30", size=2) +
              geom_vline(xintercept = vx[2:length(vx)], linetype="dotted", color = "blue", size=0.5) +
              geom_ribbon(aes(ymin=lower, ymax=upper, fill=Query), linetype=0, alpha=0.3) +
              theme_classic() +
              theme(legend.position="top",
                    axis.title.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.text.x=element_blank()) +
              ylab("Signal intensity")

            outp <- plot_grid(p, pp, ppp, ncol = 1, align = 'v', axis="b", rel_heights = c(20,1,1))
            print(outp)
         }

       }
     }
  }

  ## if ratio is required, plot ratio over input

  if(!is.null(inputfiles)){
    print("Preparing data for ratio plotting")

    pseudo_count <- 0.1

    inputMatrix_list <- scoreMatrix_list[inputlabels]
    ratiofiles <- queryfiles[!queryfiles %in% inputfiles]
    ratiolabels <- querylabels[!querylabels %in% inputlabels]
    ratioMatrix_list <- scoreMatrix_list[ratiolabels]
    for(w in featureNames){
      if(scaled_bins[w] > 0){
         for(i in seq_along(ratiolabels)){
           fullmatrix <- (ratioMatrix_list[[ratiolabels[i]]][[w]] + pseudo_count)/(inputMatrix_list[[inputlabels[i]]][[w]] + pseudo_count)

           if(rm.outlier){
             fullmatrix <- rmOutlier(fullmatrix)
           }

           ratioMatrix_list[[ratiolabels[i]]][[w]] <- fullmatrix
         }
      }else{
         ratioMatrix_list[[ratiolabels[i]]][[w]] <- NULL
      }
    }

    mplot_df <- NULL
    for(ratiolabel in ratiolabels){
      plot_df <- NULL
      for(w in featureNames){
        #w <- "utr5"
        print(w)
         if(scaled_bins[w] > 0){
           bin_num <- scaled_bins[w]
           fullmatrix <- ratioMatrix_list[[ratiolabel]][[w]]
           fullmatrix[is.na(fullmatrix)] <- 0

           print(dim(fullmatrix))

           colm <- apply(fullmatrix, 2, mean)
           colsd <- apply(fullmatrix, 2, sd)
           colse <- colsd/sqrt(nrow(fullmatrix))
           collabel <- seq(vx[w], vx[w]+bin_num-1)
           querybed <- rep(ratiolabel, bin_num)
           featuretype <- rep(w, bin_num)

           sub_df <- NULL
           sub_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=querybed, "Feature"=featuretype)
           plot_df <- rbind(plot_df, sub_df)
         }
      }

      if(smo){
        plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df=as.integer(nbins/5))$y)
        plot_df$se <- as.vector(smooth.spline(plot_df$se, df=as.integer(nbins/5))$y)
      }

      mplot_df <- rbind(mplot_df, plot_df)

    }

    mplot_df <- mutate(mplot_df, lower=Intensity-se, upper=Intensity+se)

    if(!is.null(outPrefix)){
       for(i in seq_along(ratiolabels)){
         for(beds in combn(ratiolabels, i, simplify = F)){
           print(beds)
           if(verbose || i == length(ratiolabels)){
              aplot_df <- mplot_df %>%
                filter(Query %in% beds)

              ## plot multi-sample lines with error band
              p <- ggplot(aplot_df, aes(x=Position, y=Intensity, color=Query)) + scale_fill_manual(values=color_store[1:length(beds)]) +
                geom_line(size=1) + #geom_point(color="grey30", size=2) +
                geom_vline(xintercept = vx[2:length(vx)], linetype="dotted", color = "blue", size=0.5) +
                geom_ribbon(aes(ymin=lower, ymax=upper, fill=Query), linetype=0, alpha=0.3) +
                theme_classic() +
                theme(legend.position="top",
                      axis.title.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.text.x=element_blank()) +
                ylab("Ratio over input")

              outp <- plot_grid(p, pp, ppp, ncol = 1, align = 'v', axis="b", rel_heights = c(20,1,1))
              print(outp)
           }
         }
       }
    }
  }

  if(!is.null(outPrefix)){
    dev.off()
  }

  print("plot_5parts_metagene runs successfully!")
  return(mplot_df)
}


#' plot_reference_locus
#
#' Plot reads or peak signal intensity of samples in the query files around reference locus defined in the centerfiles. The upstream and downstream windows
#' flanking genes can be given separately, a smaller window can be defined to allow statistical comparisons between samples for the same reference,
#' or between references for a given sample.
#'
#'
#' @param queryfiles, a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param querylabels, a vector of character strings serving as short labels of the queryfiles, will be used as sample labels in the plots
#' @param centerfiles, a vector of reference file names. The file should be .bed format only
#' @param centerlabels, a vector of character strings serving as short labels of the centerfiles, will be used as sample labels in the plots
#' @param inputfiles, a vector of input sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param inputlabels, a vector of character strings serving as short labels of the inputfiles, will be used as input labels in the plots
#' @param CLIP_reads, a boolean object to indicate if the bam reads should be shifted to the -1 position at the 5' of the reads
#' @param ext, a vector of two integers defining upstream and downstream boundaries of the plot window, flanking the reference locus
#' @param hl, a vector of two integers defining upstream and downstream boundaries of the highlight window, flanking the reference locus
#' @param fix_width, an integer defines how long should the reads should be extended to
#' @param norm, an boolean indicates if the intensities should be normalized by the sample library sizes
#' @param longest, a boolean object to indicate whether the output should be limited to the longest transcript of each gene
#' @param stranded, a boolean indicate whether the strand of the feature should be considered
#' @param scale, a boolean indicates whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smo, a boolean indicates whether the line should smoothed with a spline smoothing algorithm
#' @param heatmap, a boolean indicates whether a heatmap of the score matrix should be generated
#' @param rm.outlier, a boolean indicates whether a row with abnormally high values in the score matrix should be removed
#' @param genome, genome of the features, the program is mostly concerned with human, if non-human genome is used, certain features may not work
#' @param outPrefix, a character string specifying output file prefix for plots (outPrefix.pdf)
#' @param refPoint, a sting in c("start", "center", "end")
#' @param Xlab, a string denotes the label on x-axis
#' @param useScore, a boolean indicates whether to use score column of the bed file for computation of signal intensity
#' @param stats.method, a string in c("wilcox.test", "t.test"), for pair-wise groups comparisons
#'
#' @return a list of two dataframes containing the data used for plotting and tht used for statistical testing
#'
#' @examples
#'
#'
#' @export plot_reference_locus


plot_reference_locus <- function(queryfiles, centerfiles, ext=c(0,0), hl=c(0,0), shade=T, smo=FALSE, CLIP_reads=FALSE, useSizeFactor=FALSE, verbose=FALSE,
                            fix_width=0, norm=FALSE, binsize=10, refPoint="center", Xlab="Center", inputfiles=NULL, stranded=TRUE,
                            heatmap=FALSE, scale=FALSE, outPrefix=NULL, genome="hg19", rm.outlier=F, useScore=F, stats.method="wilcox.test"){


  querylabels <- names(queryfiles)
  names(querylabels) <- queryfiles
  inputlabels <- names(inputfiles)
  names(inputlabels) <- inputfiles
  centerlabels <- names(centerfiles)
  names(centerlabels) <- centerfiles

  if(!is.null(inputfiles)){
    if(length(queryfiles) != length(inputfiles)){
      stop("the number of queryfiles and inputfiles must be the same!")
    }
    names(inputlabels) <- inputfiles
    queryfiles <- c(queryfiles, inputfiles)
    querylabels <- c(querylabels, inputlabels)
  }

  ext[2] <- ext[2] - (ext[2]-ext[1])%%binsize ## to avoid binsize inconsistency, as the final binsize is dictated by bin_num
  bin_num <- round((ext[2]-ext[1])/binsize)
  colLabel <- seq(ext[1], (ext[2]-binsize), binsize)
  Ylab <- "Signal intensity"
  logYlab <- expression(paste(Log[10], " Signal intensity"))

  print("computing coverage for Sample")
  scoreMatrix_list <- list()

  bedInputs <- handle_input(queryfiles, CLIP_reads=CLIP_reads, fix_width=fix_width, useScore=useScore, outRle=TRUE, useSizeFactor=useSizeFactor, genome=genome)
  centerInputs <- handle_input(centerfiles, CLIP_reads=FALSE, fix_width=0, useScore=FALSE, outRle=FALSE, useSizeFactor=FALSE, genome=genome)

  for(queryfile in queryfiles){
    myInput <- bedInputs[[queryfile]]
    libsize <- myInput$size
    queryRegions <- myInput$query
    fileType <- myInput$type
    weight_col <- myInput$weight

    querylabel <- querylabels[queryfile]
    print(querylabel)

    for(centerfile in centerfiles){
      centerlabel <- centerlabels[centerfile]
      print(centerlabel)
      centerInput <- centerInputs[[centerfile]]
      centerGr <- centerInput$query
      windowRegions <- split(centerGr, f=factor(seq(1:length(centerGr))))
      if(refPoint %in% c("center", "start", "end")){
        windowRs <- resize(windowRegions, width = 1, fix = refPoint)
        windowRs <- promoters(windowRs, upstream = -ext[1], downstream = ext[2])
      }else{
        stop("invalid reference point!")
      }

      bin_op <- "mean"

      fullmatrix <- parallel_scoreMatrixBin(queryRegions, windowRs, bin_num, bin_op, weight_col, stranded)

      fullmatrix <- process_scoreMatrix(fullmatrix, libsize=libsize, norm=norm, scale=scale, heatmap=heatmap, rm.outlier=(rm.outlier&&is.null(inputfiles)))
      scoreMatrix_list[[querylabel]][[centerlabel]] <- fullmatrix

      if(verbose){
         print(dim(fullmatrix))
         print(length(unique(get_genomicCoordinates(windowRs))))
         #rownames(fullmatrix) <- get_genomicCoordinates(windowRs)
         write.table(fullmatrix, paste0(querylabel, "_", centerlabel, "_scoreMatrix.tab"), col.names=NA, sep="\t", quote=F)
      }
    }
  }

  color_store <- c("#00AFBB", "#E7B800", "#A0BDE0", "#0020C2", "#64E986", "#F5DEB3", "#C19A6B", "#E8A317", "#8E7618", "#A0522D", "#990012", "#CB6D51")

  plot_df <- list() # per gene, averaged over position
  stat_df <- list() # per gene, averaged over position
  twoFactor_df <- list() # gene by position

  print("collecting coverage data") ## plot multiple bed files on each center

  for(queryfile in queryfiles){

    querylabel <- querylabels[queryfile]
    print(querylabel)
    for(centerfile in centerfiles){

      if(0) centerfile <- centerfiles[1]
      centerlabel <- centerlabels[centerfile]
      print(centerlabel)

      fullmatrix <- scoreMatrix_list[[querylabel]][[centerlabel]]

      colm <- apply(fullmatrix, 2, mean)
      colsd <- apply(fullmatrix, 2, sd)
      colse <- colsd/sqrt(apply(fullmatrix, 2, length))
      collabel <- colLabel
      querybed <- as.factor(rep(querylabel, length(colm)))
      refbed <- as.factor(rep(centerlabel, length(colm)))


      sub_df <- data.frame(colm, colsd, colse, collabel, querybed, refbed)
      colnames(sub_df) <- c("Intensity", "sd", "se", "Position", "Query", "Reference")

      if(smo){
        sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df=as.integer(bin_num/5))$y)
        sub_df$se <- as.vector(smooth.spline(sub_df$se, df=as.integer(bin_num/5))$y)
      }

      plot_df[[paste(querylabel,centerlabel,sep=":")]] <- sub_df

      if(hl[2] > hl[1]){

         xmin <- which(colLabel == hl[1])
         xmax <- which(colLabel == hl[2])
         submatrix <- (fullmatrix[, xmin:xmax])
         submatrix[is.na(submatrix)] <- 0
         Intensity <- as.numeric(rowMeans(submatrix)) + 0.001 ## add 0.001 to allow log transformation of y-axis

         Query <- as.factor(rep(querylabel, length(Intensity)))
         Reference <- as.factor(rep(centerlabel, length(Intensity)))
         subdf <- data.frame(Intensity, Query, Reference)
         stat_df[[paste(querylabel,centerlabel,sep=":")]] <- subdf

         if(verbose){
               tictoc::tic()

            df_list <- lapply(seq(1:nrow(submatrix)), function(j){
              df <- data.frame(t(submatrix[j,]), as.factor(colLabel[xmin:xmax]),
                                                       as.factor(rep(querylabel, ncol(submatrix))), as.factor(rep(centerlabel, ncol(submatrix))))
            })

            factor_df <- dplyr::bind_rows(df_list)

            colnames(factor_df) <- c("Intensity", "Position", "Query", "Reference")
            twoFactor_df[[paste(querylabel,centerlabel,sep=":")]] <-  factor_df
            tictoc::toc()
         }
      }
    }
  }


  print("Transforming lists to dataframes")
  mplot_dt <- dplyr::bind_rows(plot_df) %>%
      mutate(Group=paste0(Query, ":", Reference), .keep="all") %>%
      mutate(lower=Intensity-se, upper=Intensity+se, .keep="all")

  mstat_dt <- NULL
  mtwoFactor_dt <- NULL
  if(hl[2] > hl[1]){
     mstat_dt <- dplyr::bind_rows(stat_df) %>%
       mutate(Group=as.factor(paste0(Query, ":", Reference)), .keep="all")

     if(verbose) mtwoFactor_dt <- dplyr::bind_rows(twoFactor_df) %>%
       mutate(Group=as.factor(paste0(Query, ":", Reference)), .keep="all")
  }

  if(!is.null(outPrefix)){
    while(!is.null(dev.list())){
      dev.off()
    }
    pdf(paste(outPrefix, "pdf", sep="."), height=8, width=12)
    if(verbose){
       if(hl[2] > hl[1]){
          write.table(mstat_dt, paste(outPrefix, "_multiRef.tsv", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
          sink(paste0(outPrefix, "_TukeyHSD.txt"))
       }
       write.table(mplot_dt, paste(outPrefix, "_multiRef_bin.tsv", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
    }
  }

  print("Plotting profile and boxplot")
  for(i in seq_along(querylabels)){
    for(beds in combn(querylabels, i, simplify = F)){
      for(j in seq_along(centerlabels)){
        for(centers in combn(centerlabels, j, simplify = F)){
          print(beds)
          print(centers)

          aplot_df <- mplot_dt %>%
            filter(Query %in% beds & Reference %in% centers)

          p <- ggplot(aplot_df, aes(x=Position, y=Intensity, color=Group)) + scale_fill_manual(values=color_store[1:(i*j)]) +
            geom_line(size=2) + geom_point(color="grey30", size=1) +
            geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
            geom_ribbon(aes(ymin=lower, ymax=upper, fill=Group), linetype=0, alpha=0.3) +
            theme_classic() + theme(legend.position="top") + xlab(Xlab) + ylab(Ylab) +
            theme(axis.title.x = element_text(face="bold", size=10), axis.title.y = element_text(face="bold", size=10))

          if(shade) p <- p + annotate("rect", xmin=hl[1], xmax=hl[2], ymin=-Inf, ymax=Inf, fill="grey", color="grey", alpha=0.3)

          if((i == 1 && j > 1) || (i > 1 && j == 1)){
            if(hl[2] > hl[1]){
               astat_df <- mstat_dt %>%
                 filter(Query %in% beds & Reference %in% centers)

               if(verbose){
                  atwoFactor_df <- mtwoFactor_dt %>%
                     filter(Query %in% beds & Reference %in% centers)
                  aovTukeyHSD(df=atwoFactor_df)
               }

               if(j > 1){
                 comp <- combn(seq_along(centers),2, simplify=F)
                  if(0){
                 ps1old <- ggplot(astat_df, aes(x=Reference, y=Intensity, fill=Reference)) + scale_fill_manual(values=color_store[1:j]) +
                   geom_boxplot2(width = 0.8, width.errorbar = 0.5) +
                   theme_classic() +
                   theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                   theme(legend.position = "none") +
                   labs(y=Ylab) +
                   theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                   stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
                   #geom_signif(comparisons = comp, test=stats.method, map_signif_level=TRUE, step_increase = 0.1) +
                   scale_y_continuous(trans='log10')


                  }

                 ps1 <- ggplot(astat_df, aes(x=Reference, y=Intensity, fill=Reference)) + scale_fill_manual(values=color_store[1:j]) +
                    geom_boxplot(notch=TRUE) +
                    theme_classic() +
                    theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                    theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                    theme(legend.position = "none") +
                    labs(y=Ylab) + #expression(paste(Log[10], "(Signal intensity)"))) +
                    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
                    geom_signif(comparisons = comp, test=stats.method, map_signif_level=TRUE) +
                    scale_y_continuous(trans='log10')

                 ps1_wo_outlier <- ggplot(astat_df, aes(x=Reference, y=Intensity, fill=Reference)) + scale_fill_manual(values=color_store[1:j]) +
                    geom_boxplot2(width = 0.8, width.errorbar = 0.5) +
                    theme_classic() +
                    theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                    theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                    theme(legend.position = "none") +
                    labs(y=Ylab) +
                    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black")
                 #geom_signif(comparisons = comp, test=stats.method, map_signif_level=TRUE, step_increase = 0.1) +
                 #scale_y_continuous(trans='log10')

                 means_se <- astat_df %>%
                    group_by(Reference) %>%
                    summarize(mean_Intensity=mean(Intensity),
                              sd_Intensity=sd(Intensity),
                              N_Intensity=length(Intensity),
                              se=sd_Intensity/sqrt(N_Intensity),
                              upper_limit=mean_Intensity+se,
                              lower_limit=mean_Intensity-se
                    )

                 ps1_mean_se <- ggplot(means_se, aes(x=Reference, y=mean_Intensity, fill=Reference)) + scale_fill_manual(values=color_store[1:j]) +
                    geom_bar(stat="identity") +
                    geom_errorbar(aes(ymin=mean_Intensity, ymax=upper_limit)) +
                    theme_classic() +
                    theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                    theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                    theme(legend.position = "none") +
                    labs(y=Ylab)


               }else{
                 comp <- combn(seq_along(beds),2, simplify=F)

                 ps1 <- ggplot(astat_df, aes(x=Query, y=Intensity, fill=Query)) + scale_fill_manual(values=color_store[1:i]) +
                    geom_boxplot(notch=TRUE) +
                    theme_classic() +
                    theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                    theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                    theme(legend.position = "none") +
                    labs(y=Ylab) + #expression(paste(Log[10], "(Signal intensity)"))) +
                    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
                    geom_signif(comparisons = comp, test=stats.method, map_signif_level=TRUE) +
                    scale_y_continuous(trans='log10')

                 ps1_wo_outlier <- ggplot(astat_df, aes(x=Query, y=Intensity, fill=Query)) + scale_fill_manual(values=color_store[1:i]) +
                    geom_boxplot2(width = 0.8, width.errorbar = 0.5) +
                    theme_classic() +
                    theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                    theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                    theme(legend.position = "none") +
                    labs(y=Ylab) +
                    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black")
                 #geom_signif(comparisons = comp, test=stats.method, map_signif_level=TRUE, step_increase = 0.1) +
                 #scale_y_continuous(trans='log10')

                 means_se <- astat_df %>%
                    group_by(Query) %>%
                    summarize(mean_Intensity=mean(Intensity),
                              sd_Intensity=sd(Intensity),
                              N_Intensity=length(Intensity),
                              se=sd_Intensity/sqrt(N_Intensity),
                              upper_limit=mean_Intensity+se,
                              lower_limit=mean_Intensity-se
                    )

                 ps1_mean_se <- ggplot(means_se, aes(x=Query, y=mean_Intensity, fill=Query)) + scale_fill_manual(values=color_store[1:i]) +
                    geom_bar(stat="identity") +
                    geom_errorbar(aes(ymin=mean_Intensity, ymax=upper_limit)) +
                    theme_classic() +
                    theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                    theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                    theme(legend.position = "none") +
                    labs(y=Ylab)


               }


               print(plot_grid(p, ps1, ps1_wo_outlier, ps1_mean_se, ncol = 2, rel_widths = c(1,1)))
            }else{
               print(p)
            }
          }else if((i == 1 && j == 1) || (i == length(querylabels) && j == length(centerlabels))){
            print(p)
          }

        }
      }
    }
  }

  if(!is.null(inputfiles)){
    print("Computing Ratio over input")
    Ylab <- "Ratio-over-Input"
    logYlab <- expression(paste(Log[10], " Ratio-over-Input"))
    pseudo_count <- 0.1

    inputMatrix_list <- scoreMatrix_list[inputlabels]
    ratiofiles <- queryfiles[!queryfiles %in% inputfiles]
    ratiolabels <- querylabels[!querylabels %in% inputlabels]
    ratioMatrix_list <- scoreMatrix_list[ratiolabels]
    for(centerfile in centerfiles){
      centerlabel <- centerlabels[centerfile]
      for(i in seq_along(querylabels)){
        fullmatrix <- (ratioMatrix_list[[ratiolabels[i]]][[centerlabel]] + pseudo_count)/(inputMatrix_list[[inputlabels[i]]][[centerlabel]] + pseudo_count)
        fullmatrix[is.na(fullmatrix)] <- 0

        if(rm.outlier){
          fullmatrix <- rmOutlier(fullmatrix)
        }

        ratioMatrix_list[[ratiolabels[i]]][[centerlabel]] <- fullmatrix
      }
    }

    plot_df <- list()
    stat_df <- list()
    twoFactor_df <- list()

    print("collecting ratio data") ## plot multiple bed files on each center

    for(ratiolabel in ratiolabels){

      print(ratiolabel)
      for(centerfile in centerfiles){

        if(0) centerfile <- centerfiles[1]
        centerlabel <- centerlabels[centerfile]
        print(centerlabel)

        fullmatrix <- ratioMatrix_list[[ratiolabel]][[centerlabel]]

        colm <- apply(fullmatrix, 2, mean)
        colsd <- apply(fullmatrix, 2, sd)
        colse <- colsd/sqrt(apply(fullmatrix, 2, length))
        collabel <- colLabel
        querybed <- as.factor(rep(ratiolabel, length(colm)))
        refbed <- as.factor(rep(centerlabel, length(colm)))


        sub_df <- data.frame(colm, colsd, colse, collabel, querybed, refbed)
        colnames(sub_df) <- c("Intensity", "sd", "se", "Position", "Query", "Reference")

        if(smo){
          sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df=as.integer(bin_num/5))$y)
          sub_df$se <- as.vector(smooth.spline(sub_df$se, df=as.integer(bin_num/5))$y)
        }

        plot_df[[paste(ratiolabel,centerlabel,sep=":")]] <- sub_df

        if(hl[2] > hl[1]){
           xmin <- which(colLabel == hl[1])
           xmax <- which(colLabel == hl[2])
           submatrix <- (fullmatrix[, xmin:xmax])
           submatrix[is.na(submatrix)] <- 0
           Intensity <- as.numeric(rowMeans(submatrix)) + 0.001 ## add 0.001 to allow log transformation of y-axis

           Query <- as.factor(rep(ratiolabel, length(Intensity)))
           Reference <- as.factor(rep(centerlabel, length(Intensity)))
           subdf <- data.frame(Intensity, Query, Reference)
           stat_df[[paste(ratiolabel,centerlabel,sep=":")]] <- subdf

           if(verbose){
               tictoc::tic()

              df_list <- lapply(seq(1:nrow(submatrix)), function(j){
                 df <- data.frame(t(submatrix[j,]), as.factor(colLabel[xmin:xmax]),
                                  as.factor(rep(ratiolabel, ncol(submatrix))), as.factor(rep(centerlabel, ncol(submatrix))))
              })

              factor_df <- dplyr::bind_rows(df_list)

              colnames(factor_df) <- c("Intensity", "Position", "Query", "Reference")
              twoFactor_df[[paste(ratiolabel,centerlabel,sep=":")]] <-  factor_df
              tictoc::toc()
           }
        }
      }
    }

    print("Transforming ratio data")
    mplot_dt <- dplyr::bind_rows(plot_df) %>%
      mutate(Group=paste0(Query, ":", Reference), .keep="all") %>%
      mutate(lower=Intensity-se, upper=Intensity+se, .keep="all")

    mstat_dt <- NULL
    mtwoFactor_dt <- NULL
    if(hl[2] > hl[1]){
       mstat_dt <- dplyr::bind_rows(stat_df) %>%
         mutate(Group=as.factor(paste0(Query, ":", Reference)), .keep="all")

       if(verbose) mtwoFactor_dt <- dplyr::bind_rows(twoFactor_df) %>%
         mutate(Group=as.factor(paste0(Query, ":", Reference)), .keep="all")
    }

    if(!is.null(outPrefix) && verbose){
       if(hl[2] > hl[1]) write.table(mstat_dt, paste(outPrefix, "_multiRef_ratio.tsv", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
      write.table(mplot_dt, paste(outPrefix, "_multiRef_bin_ratio.tsv", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
    }

    print("Plotting ratio profile and boxplot")
    for(i in seq_along(ratiolabels)){
      for(beds in combn(ratiolabels, i, simplify = F)){
        for(j in seq_along(centerlabels)){
          for(centers in combn(centerlabels, j, simplify = F)){
            print(beds)
            print(centers)

            aplot_df <- mplot_dt %>%
              filter(Query %in% beds & Reference %in% centers)

            p <- ggplot(aplot_df, aes(x=Position, y=Intensity, color=Group)) + scale_fill_manual(values=color_store[1:(i*j)]) +
              geom_line(size=2) + geom_point(color="grey30", size=1) +
              geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
              geom_ribbon(aes(ymin=lower, ymax=upper, fill=Group), linetype=0, alpha=0.3) +
              annotate("rect", xmin=hl[1], xmax=hl[2], ymin=-Inf, ymax=Inf, fill="grey", color="grey", alpha=0.3) +
              theme_classic() + theme(legend.position="top") + xlab(Xlab) + ylab(Ylab) +
              theme(axis.title.x = element_text(face="bold", size=10), axis.title.y = element_text(face="bold", size=10))

            if((i == 1 && j > 1) || (i > 1 && j == 1)){
               if(hl[2] > hl[1]){
                 astat_df <- mstat_dt %>%
                   filter(Query %in% beds & Reference %in% centers)

                 if(verbose){
                    atwoFactor_df <- mtwoFactor_dt %>%
                       filter(Query %in% beds & Reference %in% centers)
                    aovTukeyHSD(df=atwoFactor_df)
                 }


                 if(j > 1){
                    comp <- combn(seq_along(centers),2, simplify=F)

                    ps1 <- ggplot(astat_df, aes(x=Reference, y=Intensity, fill=Reference)) + scale_fill_manual(values=color_store[1:j]) +
                       geom_boxplot(notch=TRUE) +
                       theme_classic() +
                       theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                       theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                       theme(legend.position = "none") +
                       labs(y=Ylab) + #expression(paste(Log[10], "(Signal intensity)"))) +
                       stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
                       geom_signif(comparisons = comp, test=stats.method, map_signif_level=TRUE) +
                       scale_y_continuous(trans='log10')

                    ps1_wo_outlier <- ggplot(astat_df, aes(x=Reference, y=Intensity, fill=Reference)) + scale_fill_manual(values=color_store[1:j]) +
                       geom_boxplot2(width = 0.8, width.errorbar = 0.5) +
                       theme_classic() +
                       theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                       theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                       theme(legend.position = "none") +
                       labs(y=Ylab) +
                       stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black")
                    #geom_signif(comparisons = comp, test=stats.method, map_signif_level=TRUE, step_increase = 0.1) +
                    #scale_y_continuous(trans='log10')

                    means_se <- astat_df %>%
                       group_by(Reference) %>%
                       summarize(mean_Intensity=mean(Intensity),
                                 sd_Intensity=sd(Intensity),
                                 N_Intensity=length(Intensity),
                                 se=sd_Intensity/sqrt(N_Intensity),
                                 upper_limit=mean_Intensity+se,
                                 lower_limit=mean_Intensity-se
                       )

                    ps1_mean_se <- ggplot(means_se, aes(x=Reference, y=mean_Intensity, fill=Reference)) + scale_fill_manual(values=color_store[1:j]) +
                       geom_bar(stat="identity") +
                       geom_errorbar(aes(ymin=mean_Intensity, ymax=upper_limit)) +
                       theme_classic() +
                       theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                       theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                       theme(legend.position = "none") +
                       labs(y=Ylab)



                 }else{
                    comp <- combn(seq_along(beds),2, simplify=F)

                    ps1 <- ggplot(astat_df, aes(x=Query, y=Intensity, fill=Query)) + scale_fill_manual(values=color_store[1:i]) +
                       geom_boxplot(notch=TRUE) +
                       theme_classic() +
                       theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                       theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                       theme(legend.position = "none") +
                       labs(y=Ylab) + #expression(paste(Log[10], "(Signal intensity)"))) +
                       stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
                       geom_signif(comparisons = comp, test=stats.method, map_signif_level=TRUE) +
                       scale_y_continuous(trans='log10')

                    ps1_wo_outlier <- ggplot(astat_df, aes(x=Query, y=Intensity, fill=Query)) + scale_fill_manual(values=color_store[1:i]) +
                       geom_boxplot2(width = 0.8, width.errorbar = 0.5) +
                       theme_classic() +
                       theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                       theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                       theme(legend.position = "none") +
                       labs(y=Ylab) +
                       stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black")
                    #geom_signif(comparisons = comp, test=stats.method, map_signif_level=TRUE, step_increase = 0.1) +
                    #scale_y_continuous(trans='log10')

                    means_se <- astat_df %>%
                       group_by(Query) %>%
                       summarize(mean_Intensity=mean(Intensity),
                                 sd_Intensity=sd(Intensity),
                                 N_Intensity=length(Intensity),
                                 se=sd_Intensity/sqrt(N_Intensity),
                                 upper_limit=mean_Intensity+se,
                                 lower_limit=mean_Intensity-se
                       )

                    ps1_mean_se <- ggplot(means_se, aes(x=Query, y=mean_Intensity, fill=Query)) + scale_fill_manual(values=color_store[1:i]) +
                       geom_bar(stat="identity") +
                       geom_errorbar(aes(ymin=mean_Intensity, ymax=upper_limit)) +
                       theme_classic() +
                       theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                       theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                       theme(legend.position = "none") +
                       labs(y=Ylab)

                 }


                 print(plot_grid(p, ps1, ps1_wo_outlier, ps1_mean_se, ncol = 2, rel_widths = c(1,1)))
               }else{
                  print(p)
               }
            }else if((i == 1 && j == 1) || (i == length(ratiolabels) && j == length(centerlabels))){
              print(p)
            }

          }
        }
      }
    }

  }

  if(!is.null(outPrefix)){
    dev.off()
    if((hl[2] > hl[1]) && verbose) sink()
  }

  return(list("plot"=mplot_dt, "stat"=mstat_dt))
}

#' plot_reference_locus_random
#
#' Plot reads or peak signal intensity of samples in the query files around reference locus defined in the centerfiles. The upstream and downstream windows flanking genes
#' can be given separately, a smaller window can be defined to allow statistical comparisons between reference and random loci.
#'
#' @param queryfiles, a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param querylabels, a vector of character strings serving as short labels of the queryfiles, will be used as sample labels in the plots
#' @param centerfiles, a vector of reference file names. The file should be .bed format only
#' @param centerlabels, a vector of character strings serving as short labels of the centerfiles, will be used as sample labels in the plots
#' @param inputfiles, a vector of input sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param inputlabels, a vector of character strings serving as short labels of the inputfiles, will be used as input labels in the plots
#' @param CLIP_reads, a boolean object to indicate if the bam reads should be shifted to the -1 position at the 5' of the reads
#' @param ext, a vector of two integers defining upstream and downstream boundaries of the plot window, flanking the reference locus
#' @param hl, a vector of two integers defining upstream and downstream boundaries of the highlight window, flanking the reference locus
#' @param fix_width, an integer defines how long should the reads should be extended to
#' @param norm, an boolean indicates if the intensities should be normalized by the sample library sizes
#' @param longest, a boolean object to indicate whether the output should be limited to the longest transcript of each gene
#' @param stranded, a boolean indicate whether the strand of the feature should be considered
#' @param scale, a boolean indicates whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smo, a boolean indicates whether the line should smoothed with a spline smoothing algorithm
#' @param heatmap, a boolean indicates whether a heatmap of the score matrix should be generated
#' @param rm.outlier, a boolean indicates whether a row with abnormally high values in the score matrix should be removed
#' @param genome, genome of the features, the program is mostly concerned with human, if non-human genome is used, certain features may not work
#' @param outPrefix, a character string specifying output file prefix for plots (outPrefix.pdf)
#' @param refPoint, a sting in c("start", "center", "end")
#' @param Xlab, a string denotes the label on x-axis
#' @param n_random, an integer denotes the number of randomization should be formed
#' @param stats.method, a string in c("wilcox.test", "t.test"), for pair-wise groups comparisons
#'
#' @return a dataframe containing the data used for plotting
#'
#' @examples
#'
#'
#' @export plot_reference_locus_random

plot_reference_locus_with_random <- function(queryfiles, centerfiles, txdb, ext=c(0,0), hl=c(0,0), shade=F, useSizeFactor=FALSE,
                                              smo=FALSE,  CLIP_reads=FALSE,  fix_width=0, norm=FALSE, binsize=10, refPoint="center", Xlab="Center",
                                              inputfiles=NULL, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=NULL, genome="hg19",
                                              rm.outlier=F, n_random=1, stats.method="wilcox.test", useScore=FALSE){



  querylabels <- names(queryfiles)
  names(querylabels) <- queryfiles
  inputlabels <- names(inputfiles)
  names(inputlabels) <- inputfiles
  centerlabels <- names(centerfiles)
  names(centerlabels) <- centerfiles

  if(!is.null(inputfiles)){
    if(length(queryfiles) != length(inputfiles)){
      stop("the number of queryfiles and inputfiles must be the same!")
    }
    names(inputlabels) <- inputfiles
    queryfiles <- c(queryfiles, inputfiles)
    querylabels <- c(querylabels, inputlabels)
  }

  ext[2] <- ext[2] - (ext[2]-ext[1])%%binsize ## to avoid binsize inconsistency, as the final binsize is dictated by bin_num
  bin_num <- round((ext[2]-ext[1])/binsize)
  bin_op <- "mean"
  colLabel <- seq(ext[1], (ext[2]-binsize), binsize)
  Ylab <- "Signal intensity"
  logYlab <- expression(paste(Log[10], " Signal intensity"))

  ## ranges for overlap count
  rb <- ext[2] #broad
  rn <- hl[2]  #narrow
  color_store <- c("#00AFBB", "#E7B800", "#A0BDE0", "#0020C2", "#64E986", "#F5DEB3", "#C19A6B", "#E8A317", "#8E7618", "#A0522D", "#990012", "#CB6D51")

  ## get protein-coding genes features

  print("Collecting protein_coding gene features")
  exons <- gtf_to_bed_longest_tx(txdb, "exon", longest=TRUE)
  utr5 <- gtf_to_bed_longest_tx(txdb, "utr5", longest=TRUE)
  utr3 <- gtf_to_bed_longest_tx(txdb, "utr3", longest=TRUE)
  cds <- gtf_to_bed_longest_tx(txdb, "cds", longest=TRUE)
  gene <- gtf_to_bed_longest_tx(txdb, "transcript", longest=TRUE)


  region_list <- list(#"Transcript" = exons$GRanges,
                      "5'UTR" = utr5$GRanges,
                      "CDS" = cds$GRanges,
                      "3'UTR" = utr3$GRanges,
                      #"Gene" = gene$GRanges,
                      "unrestricted" = NULL)

  print(lapply(region_list, length))

  print("computing coverage for Sample")
  scoreMatrix_list <- list()
  quantile_list <- list()
  scoreMatrix_list_random <- list()
  quantile_list_random <- list()

  bedInputs <- handle_input(queryfiles, CLIP_reads=CLIP_reads, fix_width=fix_width,  useScore=useScore, outRle=FALSE, useSizeFactor=useSizeFactor, genome=genome)
  centerInputs <- handle_input(centerfiles, CLIP_reads=FALSE, fix_width=0, useScore=FALSE, outRle=FALSE, useSizeFactor=FALSE, genome=genome)

  for(queryfile in queryfiles){
    if(0){
      queryfile <- queryfiles[1]
    }
    myInput <- bedInputs[[queryfile]]
    libsize <- myInput$size
    queryRegions <- myInput$query
    fileType <- myInput$type
    weight_col <- myInput$weight

    querylabel <- querylabels[queryfile]
    print(querylabel)
    print(libsize)

    for(centerfile in centerfiles){
      if(0){
        centerfile <- centerfiles[1]
      }
      centerlabel <- centerlabels[centerfile]
      print(centerlabel)
      centerInput <- centerInputs[[centerfile]]
      #centerGr <- centerInput$query
      #windowRegions <- split(centerGr, f=factor(seq(1:length(centerGr))))

      windowRegionsALL <- centerInput$query


      for(regionName in names(region_list)){
        if(0){
          regionName <- names(region_list)[3]
        }
        print(paste("Processing genomic region", regionName))

        if(regionName == "unrestricted"){
          windowRegions <- windowRegionsALL
        }else{
          region <- region_list[[regionName]]
          windowRegions <- filter_by_overlaps_stranded(windowRegionsALL, region)
        }
        refsize <- length(windowRegions)

        if(refPoint %in% c("center", "start", "end")){
          windowRs <- resize(windowRegions, width = 1, fix = refPoint)
          windowRs <- promoters(windowRs, upstream = -ext[1], downstream = ext[2])

          windowsbroad <- resize(windowRegions, width = 2*rb, fix = refPoint)
          windowsnarrow <- resize(windowRegions, width = 2*rn, fix = refPoint)

          broad <- countOverlaps(windowsbroad, queryRegions)
          qbroad <- quantile(broad, probs=seq(0,1, 0.01))

          narrow <- countOverlaps(windowsnarrow, queryRegions)
          qnarrow <- quantile(narrow, probs=seq(0,1, 0.01))


          quantile_list[[querylabel]][[centerlabel]][[regionName]] <- list("broad"=broad,
                                                                         "narrow"=narrow,
                                                                         "qbroad"=qbroad,
                                                                         "qnarrow"=qnarrow,
                                                                         "refsize"=refsize)
        }else{
          stop("invalid reference point!")
        }

        queryRegionsRle <- coverage(queryRegions, weight=weight_col)
        windowRs <- split(windowRs, f=factor(seq(1:length(windowRs))))

        fullmatrix <- parallel_scoreMatrixBin(queryRegionsRle, windowRs, bin_num, bin_op, weight_col, stranded)
        fullmatrix <- process_scoreMatrix(fullmatrix, libsize=libsize, norm=norm, scale=scale, heatmap=heatmap, rm.outlier=(rm.outlier&&is.null(inputfiles)))

        scoreMatrix_list[[querylabel]][[centerlabel]][[regionName]] <- fullmatrix


        ## create randomized centers, repeat n_random times

        random_results <- lapply(seq(1, n_random), function(i){
          random_points <- sample(ext[1]:ext[2], length(windowRegions), replace=T)
          rwindowRegions <- shift(windowRegions, shift=random_points)

          windowsbroad <- resize(rwindowRegions, width = 2*rb, fix = refPoint)
          windowsnarrow <- resize(rwindowRegions, width = 2*rn, fix = refPoint)

          rbroad <- countOverlaps(windowsbroad, queryRegions)
          rqbroad <- quantile(rbroad, probs=seq(0,1, 0.01))

          rnarrow <- countOverlaps(windowsnarrow, queryRegions)
          rqnarrow <- quantile(rnarrow, probs=seq(0,1, 0.01))

          windowR <- resize(rwindowRegions, width = 1, fix = refPoint)
          windowR <- promoters(windowR, upstream = -ext[1], downstream = ext[2])
          windowR <- split(windowR, f=factor(seq(1:length(windowR))))
          return(list("windowR"=windowR, "broad"=rbroad, "qbroad"=rqbroad, "narrow"=rnarrow, "qnarrow"=rqnarrow))
        })

        windowRs <- lapply(random_results, function(x) x$windowR)

        print(system.time(
        random_matricies <- lapply(windowRs, function(x){
          parallel_scoreMatrixBin(queryRegions, x, bin_num, bin_op, weight_col, stranded)
        })))

        ## unify the dimensions of random_matrices
        dims <- sapply(random_matricies, dim)
        minrows <- min(dims[1,])
        random_matricies <- lapply(random_matricies, function(x)x[1:minrows,])

        a3d <- array(unlist(random_matricies), c(dim(random_matricies[[1]]), length(random_matricies))) ## turn the list of matrices into a 3-d array
        fullmatrix <- apply(a3d, 1:2, mean)

        fullmatrix <- process_scoreMatrix(fullmatrix, libsize=libsize, norm=norm, scale=scale, heatmap=heatmap, rm.outlier=(rm.outlier&&is.null(inputfiles)))

        scoreMatrix_list_random[[querylabel]][[centerlabel]][[regionName]] <- fullmatrix

        random_count_list <- lapply(c("broad", "qbroad", "narrow", "qnarrow"), function(v){
          alist <- lapply(random_results, function(x) x[[v]])
          amatrix <- do.call(cbind, alist)
          amatrix_m <- apply(amatrix, 1, mean)
          amatrix_se <- apply(amatrix, 1, sd)/sqrt(n_random)
          adf <- data.frame("mean"=amatrix_m, "se"=amatrix_se) %>%
            arrange(mean)

          return(as.list(adf))
        })

        names(random_count_list) <- c("broad", "qbroad", "narrow", "qnarrow")
        quantile_list_random[[querylabel]][[centerlabel]][[regionName]] <- random_count_list
      }
    }
  }


  ## plot reference center and random center for each bed
  if(!is.null(outPrefix)){
    while(!is.null(dev.list())){
      dev.off()
    }
    pdf(paste(outPrefix, "pdf", sep="."), height=8, width=12)
  }

  for(queryfile in queryfiles){
    if(0){
      queryfile <- queryfiles[1]
    }
    print(paste("Processing bed", queryfile))
    querylabel <- querylabels[queryfile]

    for(centerfile in centerfiles){
      if(0){
        centerfile <- centerfiles[1]
      }
      print(paste("Processing reference", centerfile))
      centerlabel <- centerlabels[centerfile]

      for(regionName in names(region_list)){
        if(0){
          regionName <- names(region_list)[3]
        }
        print(paste("Processing genomic region", regionName))
        plot_df <- NULL
        stat_df <- NULL
        countOverlap_df <- NULL

        fullmatrix_list <- list(scoreMatrix_list[[querylabel]][[centerlabel]][[regionName]],
                                scoreMatrix_list_random[[querylabel]][[centerlabel]][[regionName]])
        names(fullmatrix_list) <- c(centerlabel, "Random")

        for(alabel in c(centerlabel, "Random")){
          fullmatrix <- fullmatrix_list[[alabel]]
          colm <- apply(fullmatrix, 2, mean)
          colsd <- apply(fullmatrix, 2, sd)
          colse <- colsd/sqrt(apply(fullmatrix, 2, length))
          collabel <- colLabel
          querybed <- as.factor(rep(querylabel, length(colm)))
          refbed <- as.factor(rep(alabel, length(colm)))

          sub_df <- data.frame(colm, colsd, colse, collabel, querybed, refbed)
          colnames(sub_df) <- c("Intensity", "sd", "se", "Position", "Query", "Reference")

          if(smo){
            sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df=as.integer(bin_num/5))$y)
            sub_df$se <- as.vector(smooth.spline(sub_df$se, df=as.integer(bin_num/5))$y)
          }

          sub_df <- mutate(sub_df, lower=Intensity-se, upper=Intensity+se)
          plot_df <- rbind(plot_df, sub_df)

          if(hl[2] > hl[1]){
             xmin <- which(colLabel == hl[1])
             xmax <- which(colLabel == hl[2])
             submatrix <- (fullmatrix[, xmin:xmax])
             submatrix[is.na(submatrix)] <- 0
             Intensity <- as.numeric(rowMeans(submatrix)) + 0.001 ## add 0.001 to allow log transformation of y-axis

             Query <- as.factor(rep(querylabel, length(Intensity)))
             Reference <- as.factor(rep(alabel, length(Intensity)))
             subdf <- data.frame(Intensity, Query, Reference)

             stat_df <- rbind(stat_df, subdf)
          }
        }

        refsize <- quantile_list[[querylabel]][[centerlabel]][[regionName]]$refsize

        qbroad <- quantile_list[[querylabel]][[centerlabel]][[regionName]]$qbroad
        qnarrow <-quantile_list[[querylabel]][[centerlabel]][[regionName]]$qnarrow
        random_qbroad_m <- quantile_list_random[[querylabel]][[centerlabel]][[regionName]]$qbroad$mean
        random_qnarrow_m <- quantile_list_random[[querylabel]][[centerlabel]][[regionName]]$qnarrow$mean
        random_qbroad_se <- quantile_list_random[[querylabel]][[centerlabel]][[regionName]]$qbroad$se
        random_qnarrow_se <- quantile_list_random[[querylabel]][[centerlabel]][[regionName]]$qnarrow$se

        percent <- seq(100,0,-1)
        quantileOverlap_df <- data.frame(Percentage=rep(percent,4),
                                       Group=as.factor(c(rep(paste("-",rb,":",rb,sep=""), length(qbroad)), rep(paste("-",rn,":",rn,sep=""), length(qnarrow)), rep(paste("random-",rb,":",rb,sep=""), length(qbroad)), rep(paste("random-",rn,":",rn,sep=""), length(qnarrow)))),
                                       OverlapCount=c(qbroad, qnarrow, random_qbroad_m, random_qnarrow_m),
                                       se=c(rep(0, length(qbroad)), rep(0, length(qnarrow)), random_qbroad_se, random_qnarrow_se)) %>%
          mutate(lower=OverlapCount-se, upper=OverlapCount+se)

        broad <- quantile_list[[querylabel]][[centerlabel]][[regionName]]$broad
        narrow <-quantile_list[[querylabel]][[centerlabel]][[regionName]]$narrow
        random_broad_m <- quantile_list_random[[querylabel]][[centerlabel]][[regionName]]$broad$mean
        random_narrow_m <- quantile_list_random[[querylabel]][[centerlabel]][[regionName]]$narrow$mean
        random_broad_se <- quantile_list_random[[querylabel]][[centerlabel]][[regionName]]$broad$se
        random_narrow_se <- quantile_list_random[[querylabel]][[centerlabel]][[regionName]]$narrow$se

        countOverlap_df <- data.frame(Rank=c(seq_along(broad), seq_along(narrow), seq_along(random_broad_m), seq_along(random_narrow_m)),
                                      Group=as.factor(c(rep(paste("-",rb,":",rb,sep=""), length(broad)), rep(paste("-",rn,":",rn,sep=""), length(narrow)), rep(paste("random-",rb,":",rb,sep=""), length(random_broad_m)), rep(paste("random-",rn,":",rn,sep=""), length(random_narrow_m)))),
                                      OverlapCount=c(sort(broad), sort(narrow), random_broad_m, random_narrow_m),
                                      se=c(rep(0, length(broad)*2), random_broad_se, random_narrow_se)) %>%
          mutate(lower=OverlapCount-se, upper=OverlapCount+se)


        p <- ggplot(plot_df, aes(x=Position, y=Intensity, color=Reference)) + scale_fill_manual(values=c("#00AFBB", "#E7B800")) +
          geom_line(size=1) + #geom_point(color="grey30", size=2) +
          geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
          geom_ribbon(aes(ymin=lower, ymax=upper, fill=Reference), linetype=0, alpha=0.3) +
          theme_classic() + theme(legend.position="bottom") + xlab(Xlab) + ylab(Ylab) +
          theme(axis.title.x = element_text(face="bold", size=10), axis.title.y = element_text(face="bold", size=10)) +
          ggtitle(paste("Feature:", regionName, "\nReference size:", refsize, "\nSample name:", querylabel))

        ps2 <- ggplot(quantileOverlap_df, aes(x=Percentage, y=OverlapCount, color=Group)) + scale_fill_manual(values=color_store[1:4]) +
           geom_line(size=1) +
           geom_ribbon(aes(ymin=lower, ymax=upper, fill=Group), linetype=0, alpha=0.3) +
           theme_classic() + theme(legend.position="bottom") + ylab("Overlap count") + xlab("Percentage") +
           geom_hline(yintercept = 1, linetype="dotted", color = "blue", size=0.5) +
           theme(axis.title.x = element_text(face="bold", size=10), axis.title.y = element_text(face="bold", size=10)) +
           scale_y_continuous(trans='log10')
        #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)) +


        ps3 <- ggplot(countOverlap_df, aes(x=Rank, y=OverlapCount, color=Group)) + scale_fill_manual(values=color_store[1:4]) +
           geom_line(size=1) +
           geom_ribbon(aes(ymin=lower, ymax=upper, fill=Group), linetype=0, alpha=0.3) +
           theme_classic() + theme(legend.position="bottom") + ylab("Overlap count") + xlab("Rank") +
           theme(axis.title.x = element_text(face="bold", size=10), axis.title.y = element_text(face="bold", size=10))

        if(shade) p <- p + annotate("rect", xmin=hl[1], xmax=hl[2], ymin=-Inf, ymax=Inf, fill="grey", color="grey", alpha=0.3)

        if(hl[2] > hl[1]){

           ps1 <- ggplot(stat_df, aes(x=Reference, y=Intensity, fill=Reference)) + scale_fill_manual(values=c("#00AFBB", "#E7B800")) +
             geom_boxplot(notch=TRUE) +
             theme_classic() +
             theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
             theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
             theme(legend.position = "none") +
             labs(y="Signal intensity") + #expression(paste(Log[10], "(Signal intensity)"))) +
             scale_x_discrete(limits=c("Random", centerlabel)) + # labels=c("TSN" = expression(paste(m^6,"Am")), featureName = expression(paste("5'-UTR ", m^6, "A"))))
             stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
             geom_signif(comparisons = list(c(1, 2)), test=stats.method, map_signif_level=TRUE) +
             scale_y_continuous(trans='log10')

           ps1_wo_outlier <- ggplot(stat_df, aes(x=Reference, y=Intensity, fill=Reference)) + scale_fill_manual(values=c("#00AFBB", "#E7B800")) +
              geom_boxplot2(width = 0.8, width.errorbar = 0.5) +
              theme_classic() +
              theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
              theme(legend.position = "none") +
              labs(y=Ylab) +
              theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
              stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black")
              #geom_signif(comparisons = comp, test=stats.method, map_signif_level=TRUE, step_increase = 0.1) +
              #scale_y_continuous(trans='log10')

           means_se <- stat_df %>%
              group_by(Reference) %>%
              summarize(mean_Intensity=mean(Intensity),
                        sd_Intensity=sd(Intensity),
                        N_Intensity=length(Intensity),
                        se=sd_Intensity/sqrt(N_Intensity),
                        upper_limit=mean_Intensity+se,
                        lower_limit=mean_Intensity-se
              )

           ps1_mean_se <- ggplot(means_se, aes(x=Reference, y=mean_Intensity, fill=Reference)) + scale_fill_manual(values=c("#00AFBB", "#E7B800")) +
              geom_bar(stat="identity") +
              geom_errorbar(aes(ymin=mean_Intensity, ymax=upper_limit)) +
              theme_classic() +
              theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
              theme(legend.position = "none") +
              labs(y=Ylab) +
              theme(axis.text.x = element_text(face="bold", size=10, color="black"))


           print(plot_grid(p, ps1, ps1_wo_outlier, ps1_mean_se, ncol = 2, rel_widths = c(1,1)))
        }else{
           print(p)
        }
      }
    }
  }

  if(!is.null(inputfiles)){
    Ylab <- "Ratio-over-Input"
    logYlab <- expression(paste(Log[10], " Ratio-over-Input"))
    pseudo_count <- 0.1

    ratiofiles <- queryfiles[!queryfiles %in% inputfiles]
    ratiolabels <- querylabels[!querylabels %in% inputlabels]
    inputMatrix_list <- scoreMatrix_list[inputlabels]
    ratioMatrix_list <- scoreMatrix_list[ratiolabels]

    inputMatrix_list_random <- scoreMatrix_list_random[inputlabels]
    ratioMatrix_list_random <- scoreMatrix_list_random[ratiolabels]

    for(centerfile in centerfiles){
      centerlabel <- centerlabels[centerfile]
      for(regionName in names(region_list)){
        for(i in seq_along(ratiolabels)){
          rm <- ratioMatrix_list[[ratiolabels[i]]][[centerlabel]][[regionName]]
          im <- inputMatrix_list[[inputlabels[i]]][[centerlabel]][[regionName]]
          minrow <- min(nrow(rm), nrow(im))

          fullmatrix <- (rm[1:minrow,] + pseudo_count)/(im[1:minrow,] + pseudo_count)
          fullmatrix[is.na(fullmatrix)] <- 0

          if(rm.outlier){
            fullmatrix <- rmOutlier(fullmatrix)
          }

          ## scaled to the range of 0:1 to allow comparison with random ratio
          smc <- t(apply(fullmatrix, 1, scales::rescale)) ## rescale to 0:1 range
          fullmatrix <- as.data.frame(smc)

          ratioMatrix_list[[ratiolabels[i]]][[centerlabel]][[regionName]] <- fullmatrix


          ## for random centers
          rmr <- ratioMatrix_list_random[[ratiolabels[i]]][[centerlabel]][[regionName]]
          imr <- inputMatrix_list_random[[inputlabels[i]]][[centerlabel]][[regionName]]
          minrowr <- min(nrow(rmr), nrow(imr))

          fullmatrix <- (rmr[1:minrowr,] + pseudo_count)/(imr[1:minrowr,] + pseudo_count)
          fullmatrix[is.na(fullmatrix)] <- 0

          if(rm.outlier){
            fullmatrix <- rmOutlier(fullmatrix)
          }

          ## scaled to the range of 0:1 to allow comparison with random ratio
          smc <- t(apply(fullmatrix, 1, scales::rescale)) ## rescale to 0:1 range
          fullmatrix <- as.data.frame(smc)
          ratioMatrix_list_random[[ratiolabels[i]]][[centerlabel]][[regionName]] <- fullmatrix
        }
      }
    }

    for(ratiofile in ratiofiles){
      if(0){
        ratiofile <- ratiofiles[1]
      }
      print(paste("Processing ratio for query", ratiofile))
      ratiolabel <- ratiolabels[ratiofile]

      for(centerfile in centerfiles){
        if(0){
          centerfile <- centerfiles[1]
        }
        print(paste("Processing ratio for reference", centerfile))
        centerlabel <- centerlabels[centerfile]

        for(regionName in names(region_list)){
          if(0){
            regionName <- names(region_list)[3]
          }
          print(paste("Processing ratio for genomic region", regionName))
          plot_df <- NULL
          stat_df <- NULL
          countOverlap_df <- NULL
          refsize <- NULL

          fullmatrix_list <- list(ratioMatrix_list[[ratiolabel]][[centerlabel]][[regionName]],
                                  ratioMatrix_list_random[[ratiolabel]][[centerlabel]][[regionName]])
          names(fullmatrix_list) <- c(centerlabel, "Random")

          for(alabel in c(centerlabel, "Random")){

            fullmatrix <- fullmatrix_list[[alabel]]
            refsize <- nrow(fullmatrix)
            colm <- apply(fullmatrix, 2, mean)
            colsd <- apply(fullmatrix, 2, sd)
            colse <- colsd/sqrt(apply(fullmatrix, 2, length))
            collabel <- colLabel
            querybed <- as.factor(rep(querylabel, length(colm)))
            refbed <- as.factor(rep(alabel, length(colm)))

            sub_df <- data.frame(colm, colsd, colse, collabel, querybed, refbed)
            colnames(sub_df) <- c("Intensity", "sd", "se", "Position", "Query", "Reference")

            if(smo){
              sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df=as.integer(bin_num/5))$y)
              sub_df$se <- as.vector(smooth.spline(sub_df$se, df=as.integer(bin_num/5))$y)
            }

            sub_df <- mutate(sub_df, lower=Intensity-se, upper=Intensity+se)
            plot_df <- rbind(plot_df, sub_df)

            if(hl[2] > hl[1]){
               xmin <- which(colLabel == hl[1])
               xmax <- which(colLabel == hl[2])
               submatrix <- (fullmatrix[, xmin:xmax])
               submatrix[is.na(submatrix)] <- 0
               Intensity <- as.numeric(rowMeans(submatrix)) + 0.001 ## add 0.001 to allow log transformation of y-axis

               Query <- as.factor(rep(querylabel, length(Intensity)))
               Reference <- as.factor(rep(alabel, length(Intensity)))
               subdf <- data.frame(Intensity, Query, Reference)

               stat_df <- rbind(stat_df, subdf)
            }
          }


          p <- ggplot(plot_df, aes(x=Position, y=Intensity, color=Reference)) + scale_fill_manual(values=c("#00AFBB", "#E7B800")) +
            geom_line(size=1) + #geom_point(color="grey30", size=2) +
            geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
            geom_ribbon(aes(ymin=lower, ymax=upper, fill=Reference), linetype=0, alpha=0.3) +
            theme_classic() + theme(legend.position="bottom") + xlab(Xlab) + ylab(Ylab) +
            theme(axis.title.x = element_text(face="bold", size=10), axis.title.y = element_text(face="bold", size=10)) +
            ggtitle(paste("Feature:", regionName, "\nReference size:", refsize, "\nSample name:", ratiolabel))

          if(shade) p <- p + annotate("rect", xmin=hl[1], xmax=hl[2], ymin=-Inf, ymax=Inf, fill="grey", color="grey", alpha=0.3)

          if(hl[2] > hl[1]){

             ps1 <- ggplot(stat_df, aes(x=Reference, y=Intensity, fill=Reference)) + scale_fill_manual(values=c("#00AFBB", "#E7B800")) +
               geom_boxplot(notch=FALSE) +
               theme_classic() +
               theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
               theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
               theme(legend.position = "none") +
               labs(y=Ylab) + #expression(paste(Log[10], "(Signal intensity)"))) +
               scale_x_discrete(limits=c("Random", centerlabel)) + # labels=c("TSN" = expression(paste(m^6,"Am")), featureName = expression(paste("5'-UTR ", m^6, "A"))))
               stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
               geom_signif(comparisons = list(c(1, 2)), test=stats.method, map_signif_level=TRUE) +
               scale_y_continuous(trans='log10')

             ps1_wo_outlier <- ggplot(stat_df, aes(x=Reference, y=Intensity, fill=Reference)) + scale_fill_manual(values=c("#00AFBB", "#E7B800")) +
                geom_boxplot2(width = 0.8, width.errorbar = 0.5) +
                theme_classic() +
                theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                theme(legend.position = "none") +
                labs(y=Ylab) +
                theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black")
             #geom_signif(comparisons = comp, test=stats.method, map_signif_level=TRUE, step_increase = 0.1) +
             #scale_y_continuous(trans='log10')

             means_se <- stat_df %>%
                group_by(Reference) %>%
                summarize(mean_Intensity=mean(Intensity),
                          sd_Intensity=sd(Intensity),
                          N_Intensity=length(Intensity),
                          se=sd_Intensity/sqrt(N_Intensity),
                          upper_limit=mean_Intensity+se,
                          lower_limit=mean_Intensity-se
                )

             ps1_mean_se <- ggplot(means_se, aes(x=Reference, y=mean_Intensity, fill=Reference)) + scale_fill_manual(values=c("#00AFBB", "#E7B800")) +
                geom_bar(stat="identity") +
                geom_errorbar(aes(ymin=mean_Intensity, ymax=upper_limit)) +
                theme_classic() +
                theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                theme(legend.position = "none") +
                labs(y=Ylab) +
                theme(axis.text.x = element_text(face="bold", size=10, color="black"))


             print(plot_grid(p, ps1, ps1_wo_outlier, ps1_mean_se, ncol = 2, rel_widths = c(1,1)))
          }else{
             print(p)
          }
        }
      }
    }
  }

  if(!is.null(outPrefix)) dev.off()
}


#' @title handle_input
#'
#' @descirption This is a wrapper function for read NGS data in different file formats, store the input data in a list of GRanges objects or RleList objects.
#' File names end in bed|bam|bw|bigwig|bigWig|BigWig|BW|BIGWIG are recognized, and a list of files with mixed formats are allowed.
#'
#' @param inputFiles, a vector of strings denoting file names
#' @param CLIP_reads, a boolean object to indicate if the bam reads should be shifted to the -1 position at the 5' of the reads
#' @param fix_width, an integer defines how long should the reads should be extended to
#' @param useScore, a boolean object indicating whether the 'score' column of the bed file should be used in the output data structure
#' @param outRle, a boolean object indicating whether the output should be a list of RleList objects or GRanges objects
#' @param norm, a boolean object indicating whether the library size should be adjusted with a size factor
#' @param genome, a string denoting the genome version
#'
#' @return a list object with four elements, 'query' is a list GRanges objects or RleList objects, 'size' is the library size, 'type' is the input file type,
#' 'weight' is the name of the metadata column to be used as weight for coverage calculation
#'
#' @author Shuye Pu
#' @examples
#'
#' @export handle_input

handle_input <- function(inputFiles, CLIP_reads=FALSE, fix_width=0, useScore=FALSE, outRle=TRUE, useSizeFactor=FALSE, genome="hg19"){

  if(0){
    inputFiles <- queryfiles
    CLIP_reads=FALSE
    fix_width=0
    useScore=FALSE
    outRle=TRUE
    useSizeFactor=FALSE
    genome="tetrahymena"

  }

  outlist <- lapply(inputFiles, function(inputFile){

    print(paste("Reading file:", inputFile))

    if(grepl("\\.bed$", inputFile)){
      fileType <- "bed"
      out <- handle_bed(inputFile, fix_width, useScore, outRle, genome)
    }else if(grepl("\\.bam$", inputFile)){
      fileType <- "bam"
      out <- handle_bam(inputFile, CLIP_reads, fix_width, outRle, genome)
    }else if(grepl("\\.wig$", inputFile)){
      fileType <- "wig"
      out <- handle_wig(inputFile, outRle, genome)
    }else if(grepl("\\.bw|bigwig|bigWig|BigWig|BW|BIGWIG$", inputFile)){
      fileType <- "bw"
      out <- handle_bw(inputFile, outRle, genome)
    }else{
      stop(paste("The input file format is not supported:", inputFile))
    }

    out
  })

  names(outlist) <- inputFiles

  if(useSizeFactor && length(inputFiles)>1 && genome %in% c("hg19", "hg38")){
    outlist <- effectiveSize(outlist, outRle, genome)
  }

  return(outlist)
}


#' @title effectiveSize
#'
#' @descirption This is a helper function for handle_input. DESeq2::estimateSizeFactorsForMatrix function used to estimate size factor,
#' which is used to multiply library sizes. The function only works for human genome only so far.
#'
#' @param outlist, a list object with four elements, 'query' is a list GRanges objects or RleList objects, 'size' is the library size, 'type' is the input file type,
#' 'weight' is the name of the metadata column
#' @param outRle, a boolean object indicating whether the output is a list of RleList objects or GRanges objects
#' @param genome, a string denoting the genome version
#'
#' @return a list object with four elements, with the 'size' element modified.
#'
#' @author Shuye Pu
#' @examples
#'
#' @export effectiveSize
#'
effectiveSize <- function(outlist, outRle, genome="hg19"){
  print("Estimating size factor")
  ## find common chromosomes
  chr_list <- lapply(outlist, function(x) names(x$query))
  chr_list[["regular_chr"]] <- paste0("chr", c(seq(1,22), c("X", "Y", "M")))
  comchr <- sort(Reduce(intersect, chr_list))

  Rle_list <- outlist
  if(!outRle){
    Rle_list <- lapply(outlist, function(x){
      y <- x
      y$query <- coverage(x$query, weight= x$weight)
      y
    })
    names(Rle_list) <- names(outlist)

    chr_list <- lapply(Rle_list, function(x) names(x$query))
    chr_list[["regular_chr"]] <- paste0("chr", c(seq(1,22), c("X", "Y", "M")))
    comchr <- sort(Reduce(intersect, chr_list))
  }

  seqi <- Seqinfo(genome=genome)
  tilewidth <- 100000
  tileBins <- tileGenome(seqi[comchr], tilewidth=tilewidth, cut.last.tile.in.chrom=TRUE)

  cl <- start_parallel(length(outlist))
  parallel::clusterExport(cl, c("binnedAverage"))
  parallel::clusterExport(cl, c("seqi", "tilewidth", "tileBins", "comchr"), envir=environment())

  score_list <- parLapply(cl, Rle_list, function(x){
    binAverage <- binnedAverage(tileBins, x$query[comchr], varname="binned_score", na.rm=F)
    binAverage$binned_score
  })
  stop_parallel(cl)

  mat <- matrix(unlist(score_list), ncol=length(score_list), byrow=F)
  mat[is.na(mat)] <- 0
  mat <- round(mat * tilewidth) + 1
  sizeFactor <- DESeq2::estimateSizeFactorsForMatrix(mat)

  print(sizeFactor)

  ## update library size with sizeFactor
  outlist <- mapply(x=outlist, y=sizeFactor, function(x,y){
    x$size <- x$size*y
    x
  }, SIMPLIFY=FALSE, USE.NAMES = TRUE)

  #lapply(out, "[[", "size")

  outlist
}

#' @title handle_bed
#'
#' @descirption This is a function for read peaks data in bed format, store the input data in a list of GRanges objects or RleList objects.
#'
#' @param inputFile, a string denoting the bed file name
#' @param fix_width, an integer defines how long should the peaks should be extended to
#' @param useScore, a boolean object indicating whether the 'score' column of the bed file should be used in the output data structure
#' @param outRle, a boolean object indicating whether the output should be a list of RleList objects or GRanges objects
#' @param genome, a string denoting the genome version
#'
#' @details when useScore is TRUE, the score column of the bed file will be used in the metadata column 'score' of the GRanges object,
#' or the 'Values' field of the RleList object. Otherwise the value 1 will be used instead. When the intended use of the input bed is a
#' reference feature, both useScore and outRle should be set to FALSE.
#'
#'
#' @return a list object with four elements, 'query' is a list GRanges objects or RleList objects, 'size' is the library size, 'type' is the input file type,
#' 'weight' is the name of the metadata column to be used as weight for coverage calculation
#'
#' @author Shuye Pu
#' @examples
#'
#' @export handle_bed

handle_bed <- function(inputFile, fix_width=0, useScore=FALSE, outRle=TRUE, fixPoint="center", genome="hg19"){

  beddata <- read.delim2(inputFile, header=F)
  beddata <- beddata[, 1:min(6,ncol(beddata))]  ## ignore extra columns, which cause problem in import.bed()
  colnames(beddata) <- c("chr", "start", "end", "id", "score", "strand")[1:min(6,ncol(beddata))]
  queryRegions <- makeGRangesFromDataFrame(beddata, keep.extra.columns=TRUE, starts.in.df.are.0based=TRUE)
  queryRegions$score <- ifelse(ncol(beddata)==6, as.numeric(queryRegions$score), 1)

  if(genome %in% c("hg19", "hg38")){
    regular_chr <- paste0("chr", c(seq(1,22), c("X", "Y", "M")))
    queryRegions <- queryRegions[seqnames(queryRegions) %in% regular_chr]
    seqlevels(queryRegions) <- regular_chr
  }

  if(!useScore){
    score(queryRegions) <- 1
  }
  libsize <- sum(queryRegions$score, na.rm=T)
  weight_col <- "score"

  if(fix_width > 0) queryRegions <- resize(queryRegions, width=fix_width, fix=fixPoint, ignore.strand=FALSE)

  if(outRle) queryRegions <- coverage(queryRegions, weight=weight_col)  ## when bed is used is windows in ScoreMatrix, do not covert to Rle

  return(list("query"=queryRegions, "size"=libsize, "type"="bed", "weight"=weight_col))
}

#' @title handle_bam
#'
#' @descirption This is a function for read NGS reads data in bam format, store the input data in a list of GRanges objects or RleList objects.
#'
#' @param inputFile, a string denoting the bed file name
#' @param CLIP_reads, a boolean object to indicate if the reads should be shifted to their -1 position at the 5' end and shrink to length 1
#' @param fix_width, an integer defines how long should the reads should be extended to
#' @param outRle, a boolean object indicating whether the output should be a list of RleList objects or GRanges objects
#' @param genome, a string denoting the genome version
#'
#' @details The reads are filtered use mapq score >= 10 by default, only mapped reads are counted towards library size.
#'
#'
#' @return a list object with four elements, 'query' is a list GRanges objects or RleList objects, 'size' is the library size, 'type' is the input file type,
#' 'weight' is the name of the metadata column to be used as weight for coverage calculation
#'
#' @author Shuye Pu
#' @examples
#'
#' @export handle_bam
#'
handle_bam <- function(inputFile, CLIP_reads=FALSE, fix_width=0, outRle=TRUE, genome="hg19"){

  if(0){
    inputFile <- "Greenblatt_001_Strip1_A01_E-H3_3_S1_S34_sorted.bam"
    CLIP_reads=FALSE
    fix_width=0
    outRle=TRUE
    genome="tetrahymena"
  }

  ga <- readGAlignments(inputFile, use.names=T, param=ScanBamParam(mapqFilter=10))
  libsize <- sum(idxstatsBam(inputFile)$mapped)
  if(CLIP_reads){
    ## get the 5'-end -1 position of the reads, which is the crosslink sites for iCLIP reads
    queryRegions <- flank(granges(ga), width=1, both=F, start=T, ignore.strand=FALSE)
    score(queryRegions) <- 1
  }else if(fix_width > 0){
    queryRegions <- resize(granges(ga), width=fix_width, fix="start", ignore.strand=FALSE)
  }else{
    queryRegions <- stack(grglist(ga))
    score(queryRegions) <- 1
  }
  weight_col <- "score"
  if(genome %in% c("hg19", "hg38")){
    regular_chr <- paste0("chr", c(seq(1,22), c("X", "Y", "M")))
    queryRegions <- queryRegions[seqnames(queryRegions) %in% regular_chr]
    seqlevels(queryRegions) <- regular_chr
  }
  if(outRle) queryRegions <- coverage(queryRegions, weight= weight_col)
  return(list("query"=queryRegions, "size"=libsize, "type"="bam", "weight"=weight_col))

}


#' @title handle_bw
#'
#' @descirption This is a function for read NGS coverage data in bigwig format, store the input data in a list of GRanges objects or RleList objects. The input bw
#' file can be stranded or non-stranded. Library size is calculate as the sum of all coverage.
#'
#'
#' @param inputFile, a string denoting the bw file name
#' @param outRle, a boolean object indicating whether the output should be a list of RleList objects or GRanges objects
#' @param genome, a string denoting the genome version
#'
#' @details For stranded files, forward and reverse strands are stored in separate files, with '+' or 'p' in the forward strand file name and '-' or 'm' in the reverse
#' strand  file name.
#'
#'
#' @return a list object with four elements, 'query' is a list GRanges objects or RleList objects, 'size' is the library size, 'type' is the input file type,
#' 'weight' is the name of the metadata column to be used as weight for coverage calculation
#'
#' @author Shuye Pu
#' @examples
#'
#' @export handle_bw
#'
handle_bw <- function(inputFile, outRle=TRUE, genome="hg19"){
  if(0)inputFile <- queryfiles[1]
  weight_col <- "score"
  libsize <- NULL
  regular_chr <- paste0("chr", c(seq(1,22), c("X", "Y", "M")))
  #inputFile <- inputFiles[1]
  stranded <- TRUE
  if(grepl("_\\+_", inputFile)){
    neg_file <- gsub("_\\+_", "_\\-_", inputFile)
  }else if(grepl("\\.p\\.", inputFile)){
    neg_file <- gsub("\\.p\\.", "\\.m\\.", inputFile)
  }else{
    stranded <- FALSE
  }

  if(stranded){
    system.time(
    posGR <- import.bw(inputFile)
    )
    #posGR <- as_ranges(coverage(posGR, weight=weight_col)) ## merge contigous single nucleotide ranges with same score into ranges with score
    posGR <- posGR[seqnames(posGR) %in% regular_chr]
    seqlevels(posGR) <- regular_chr
    strand(posGR) <- "+"
    negGR <- import.bw(neg_file)
    #negGR <- as_ranges(coverage(negGR, weight=weight_col)) ## merge contigous single nucleotide ranges with same score into ranges with score
    negGR <- negGR[seqnames(negGR) %in% regular_chr]
    seqlevels(negGR) <- regular_chr
    strand(negGR) <- "-"

    queryRegions <- append(posGR, negGR)

  }else{
    iGR <- import.bw(inputFile)
    #iGR <- as_ranges(coverage(iGR, weight=weight_col)) ## merge contigous single nucleotide ranges with same score into ranges with score
    queryRegions <- iGR[seqnames(iGR) %in% regular_chr]
    seqlevels(queryRegions) <- regular_chr
    strand(queryRegions) <- "*"
  }


  libsize <- sum(score(queryRegions) * width(queryRegions))
  if(outRle){
    queryRegions <- coverage(queryRegions, weight=weight_col)
    score_list <- mapply(x=runValue(queryRegions), y=runLength(queryRegions), function(x,y) x*y)
    libsize <- sum(sapply(score_list, sum))
  }
  return(list("query"=queryRegions, "size"=libsize, "type"="bw", "weight"=weight_col))

}

#' @title handle_wig
#'
#' @descirption This is a function for read NGS coverage data in wig format, store the input data in a list of GRanges objects or RleList objects. The input wig
#' file can be stranded or non-stranded. Library size is calculate as the sum of all coverage.
#'
#'
#' @param inputFile, a string denoting the wig file name
#' @param outRle, a boolean object indicating whether the output should be a list of RleList objects or GRanges objects
#' @param genome, a string denoting the genome version
#'
#' @details For stranded files, forward and reverse strands are stored in separate files, with '+' or 'p' in the forward strand file name and '-' or 'm' in the reverse
#' strand  file name.
#'
#'
#' @return a list object with four elements, 'query' is a list GRanges objects or RleList objects, 'size' is the library size, 'type' is the input file type,
#' 'weight' is the name of the metadata column to be used as weight for coverage calculation
#'
#' @author Shuye Pu
#' @examples
#'
#' @export handle_wig
#'
handle_wig <- function(inputFile, outRle=TRUE, genome="hg19"){

  stranded <- FALSE
  if(grepl("_\\+_", inputFile)){
    neg_file <- gsub("_\\+_", "_\\-_", inputFile)
    stranded <- TRUE
  }else if(grepl("\\.p\\.wig", inputFile)){
    neg_file <- gsub("\\.p\\.wig", "\\.m\\.wig", inputFile)
    stranded <- TRUE
  }

  wigToBigWig(inputFile, Seqinfo(genome=genome))
  if(stranded){
    wigToBigWig(neg_file, Seqinfo(genome=genome))
  }
  bwfile <- gsub("\\.wig", "\\.bw", inputFile)

  out <- handle_bw(bwfile, outRle, genome)
  out
}




#' @title plot_bam_correlation
#'
#' @descirption plot correlation in reads coverage distributions along the genome for bam files
#'
#' @param bamfiles, a vector of strings denoting file names
#' @param bamlabels, a vector of strings denoting short labels to appeare in the plot
#' @param binsize, an integer denoting the tile width for tiling the genome, default 1000000
#' @param outPrefix, a string denoting output file name in pdf format
#' @param genome, a string denoting the genome version
#'
#' @return NULL
#'
#' @examples
#' bamfiles <- queryfiles
#' genome <- "hg19"
#' binsize <- 10000
#' outPrefix="Bam_correlation"
#' @export plot_bam_correlation
#'
plot_bam_correlation <- function(bamfiles, binsize=1000000, outPrefix="Bam_correlation", genome="hg19"){
  #genome <- "hg19"
  bamlabels <- names(bamfiles)

  print("Computing bam correlation")
  outlist <- handle_input(inputFiles=bamfiles, CLIP_reads=FALSE, fix_width=0, useScore=FALSE, outRle=TRUE, useSizeFactor=FALSE, genome=genome)
  ## find common chromosomes
  chr_list <- lapply(outlist, function(x) names(x$query))
  if(genome %in% c("hg19", "hg38")){
    chr_list[["regular_chr"]] <- paste0("chr", c(seq(1,22), c("X", "Y", "M")))
  }
  comchr <- sort(Reduce(intersect, chr_list))

  seqi <- Seqinfo(genome=genome)

  cl <- start_parallel(length(outlist))
  parallel::clusterExport(cl, c("binnedAverage"))
  tileBins <- tileGenome(seqi[comchr], tilewidth=binsize, cut.last.tile.in.chrom=TRUE)
  parallel::clusterExport(cl, c("seqi", "binsize", "comchr", "tileBins"), envir=environment())
  score_list <- parLapply(cl, outlist, function(x){
    binAverage <- binnedAverage(tileBins, x$query[comchr], varname="binned_score", na.rm=F)
    binAverage$binned_score
  })
  stop_parallel(cl)
  bins_df <- data.frame(chr=seqnames(tileBins), start=start(tileBins),end=end(tileBins),strand=strand(tileBins))
  bins <- do.call(paste, c(bins_df, sep="_"))
  mat <- matrix(unlist(score_list), ncol=length(score_list), byrow=F)
  rownames(mat) <- bins
  mat[is.na(mat)] <- 0
  mat <- mat[apply(mat, 1, sum)>0,]

  norm_factor <- sapply(outlist, function(x) x$size/1000000)

  df <- round(mat * binsize) + 1
  df <- as.data.frame(t(t(df)/norm_factor)) ## convert to counts per million (CPM)

  colnames(df) <- bamlabels

  long_df <- tidyr::pivot_longer(df, cols=colnames(df), names_to="Sample", values_to="Count") %>%
    mutate(Sample=as.factor(Sample)) %>%
     group_by(Sample) %>%
     arrange(Count) %>%
     filter(Count < 10000) %>%
     mutate(Rank=order(Count)) %>%
     mutate(Fraction=Count/max(Count), Rank=Rank/max(Rank))

  pdf(paste0(outPrefix, ".pdf"), width=10, height=8)
  p <- ggplot(data=long_df, aes(x=Rank, y=(Fraction), color=Sample)) +
     geom_line() +
     ggtitle(paste("Binned read counts distribution: bin size =", binsize)) +
     labs(x="Rank(Count)", y=paste("Fraction over highest coverage"))
  print(p)

  p <- ggplot(data=long_df, aes(x=(Count), color=Sample)) +
    stat_ecdf() +
    ggtitle(paste("Binned read counts distribution: bin size =", binsize)) +
    labs(x=expression(paste(log[2], " (Count)")), y=paste("Pencentage"))
  print(p)


  ## code from pairs example
  ## put (absolute) correlations on the upper panels,
  ## with size proportional to the correlations.
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }

  ## put histograms on the diagonal
  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
  }

  ## END code from pairs example

  if(length(bamfiles)>3) pheatmap(cor(log(df+1)), display_numbers = T)
  pairs(log(df+1), lower.panel = panel.smooth, upper.panel=panel.cor, diag.panel = panel.hist, main=paste("log (reads/bin)\nbin size =", binsize))
  dev.off()

  return(df)

}

#' annotate_peaks
#'
#' Produce a table of genes targeted by peaks, and generate plots for target gene types, peak distribution in genomic features
#'
#' @param peakfile, a string denoting the peak file name, only .bed format is allowed
#' @param peaklabel, a string denoting short labels to appeare in the plot
#' @param gtffile, a string denoting the path to the gtf or gff3 file
#' @param fiveP, extension out of the 5' boundary of genes for defining promoter: -fiveP TSS +100
#' @param threeP, extension out of the 3' boundary of genes for defining termination region: -100 TTS +threeP
#' @param simple, a boolean object indicating whether 5'UTR and 3'UTR are annotated in the gtffile
#' @param RNA, a boolean object indicating whether only peaks in transcripts should be considered for pie chart plot
#'
#' @return NULL
#'
#' @examples
#'
#' @export annotate_peaks
#'
annotate_peaks <- function(peakfile, gtffile, genome="hg19", fiveP=1000, threeP=1000, simple=FALSE, RNA=TRUE){
  if(0){
    peakfile <- "m6A_Hek_merged.thUni.crossLink_site.0.05.bed"
    genome <- "hg19"
    gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
    peaklabel <- "m6A_Hek_cits"
    fiveP <- 1000
    threeP <- 1000
    simple <- FALSE
    RNA <- T

    peakfile <- "Tetrahymena_ChIPseq_summits.bed"
    gtffile <- "C:/GREENBLATT/Nabeel/Tetrahymena/2-Genome_GFF3.gff3"
    genome <- "tetrahymena"
    fiveP <- 2000
    threeP <- 1000
    simple <- TRUE
    RNA <- F
    peaklabel <- "Tetrahymena_ChIPseq_summits"
  }

   return_table <- NULL
  peaklabel <- names(peakfile)
  peak <- handle_bed(inputFile=peakfile, fix_width=0, useScore=FALSE, outRle=FALSE, genome=genome)$query

  if(simple){
    txdb <-  makeTxDbFromGFF(gtffile)
    gene <- gtf_to_bed_longest_tx(txdb, "transcript", longest=F)$GRanges
    exon <- gtf_to_bed_longest_tx(txdb, "exon", longest=F)
    intron <- gtf_to_bed_longest_tx(txdb, "intron", longest=F)
    promoter <- promoters(txdb, upstream=fiveP, downstream=100, use.names=FALSE)
    TTS <- promoters(txdb, upstream=100, downstream=threeP, use.names=FALSE)

    targeted_gene <- as.data.frame(mergeByOverlaps(peak, gene)) %>%
      select(colnames(.)[grepl("\\.", colnames(.))])
    targeted_gene_count <- targeted_gene %>%
      filter(!is.na(gene.tx_name)) %>%
      count(gene.tx_name, name="Peak_in_geneBody")

    targeted_promoter <- as.data.frame(mergeByOverlaps(peak, promoter)) %>%
      select(colnames(.)[grepl("\\.", colnames(.))])
    targeted_promoter_count <- targeted_promoter %>%
      filter(!is.na(promoter.tx_name)) %>%
      count(promoter.tx_name, name="Peak_in_promoter")

    summary_table <- full_join(targeted_gene_count, targeted_promoter_count, by=c("gene.tx_name"="promoter.tx_name")) %>%
      arrange(desc(Peak_in_promoter))

    return_table <- targeted_gene
    write.table(targeted_gene, paste(peaklabel, "_targeted_gene.tab", sep=""), sep="\t", row.names=F, quote=F)
    write.table(targeted_promoter, paste(peaklabel, "_targeted_promoter.tab", sep=""), sep="\t", row.names=F, quote=F)
    write.table(summary_table, paste(peaklabel, "_target_summary.tab", sep=""), sep="\t", row.names=F, quote=F)

    features <- GRangesList("Promoter"=promoter,
                            "TTS"=TTS,
                            "Exon"=exon$GRanges,
                            "Intron"=intron$GRanges,
                            compress=F)

    pdf(paste0(peaklabel, "_feature_type_piechart.pdf"))

    annot = annotateWithFeatures(peak, features, strand.aware=TRUE, intersect.chr=FALSE)
    precedence_count <- annot@num.precedence
    precedence_count["NoFeature"] <- length(peak) - sum(precedence_count)

    df <- data.frame(precedence_count) %>%
      mutate(percent = precedence_count / sum(precedence_count)) %>%
      mutate(feature = rownames(.)) %>%
      mutate(labels = scales::percent(percent, accuracy = 0.1))

    p <- ggpie(df, x="percent", label="labels", lab.pos="out", fill="feature", color="white", palette="ucscgb")
    print(p)
    dev.off()

  }else{

    gff <- importGtf(filePath = gtffile)

    overlaps <- as.data.table(queryGff(queryRegions=peak, gffData=gff))

    gene_info_table <- overlaps[, c("transcript_id", "gene_id", "gene_name", "gene_biotype")]

    ## get targeted genes table,
    ## beware that utr5, utr3 and cds in txdbFeatures DO contain introns
    if(0){ #obsolete
       txdbFeatures <- getTxdbFeaturesFromGRanges(gff)

       dt <- getTargetedGenesTable(queryRegions = peak, txdbFeatures = txdbFeatures)
       dt <- dt[order(transcripts, decreasing = TRUE)]

       dt_gene <- unique(merge(dt, gene_info_table, by.x="tx_name", by.y="transcript_id", all.x=T))
       dt_gene <- dt_gene[order(transcripts, decreasing = TRUE)]
       dt_gene <- dt_gene[!is.na(gene_biotype)]
       dim(dt_gene)
       head(dt_gene)
       tail(dt_gene)
    }

    gene <- gtf_to_bed_longest_tx(txdb, "gene", longest=FALSE)
    geneGr <- stack(gene$GRangesList)
    geneOverlaps <- GenomicRanges::findOverlaps(peak, geneGr)
    peak_df <- annoGR2DF(peak[geneOverlaps@from]) %>%
       mutate(chrPeak=as.character(chr),
              startPeak=as.integer(start)-1,
              endPeak=as.integer(end),
              widthPeak=as.integer(width),
              strandPeak=as.character(strand),
              .keep="unused")
    gene_df <- geneGr[geneOverlaps@to]
    names(gene_df) <- seq_along(gene_df)
    gene_df <- annoGR2DF(gene_df) %>%
       mutate(chrGene=as.character(chr),
              startGene=as.integer(start)-1,
              endGene=as.integer(end),
              widthGene=as.integer(width),
              strandGene=as.character(strand),
              .keep="unused")

    ot <- cbind(peak_df, gene_df) %>%
       select(chrPeak, startPeak, endPeak, id, widthPeak, strandPeak, chrGene, startGene, endGene, gene_id, widthGene, strandGene)
    ot_gene <- merge(ot, gene_info_table, by.x="gene_id", by.y="gene_id", all.x=T) %>%
       select(chrPeak, startPeak, endPeak, id, widthPeak, strandPeak, chrGene, startGene, endGene, gene_id, widthGene, strandGene, gene_name, gene_biotype) %>%
       unique()

    return_table <- ot_gene
    write.table(ot_gene, paste(peaklabel, "_targeted_gene.tab", sep=""), sep="\t", row.names=F, quote=F)

    # To find out the distribution of the query regions across gene types:
    biotype_col <- grep('gene_biotype', colnames(overlaps), value = T)
    df <- overlaps[,length(unique(queryIndex)), by = biotype_col] %>%
      rename_with(~ c("gene_type", "count"))
    intergenic_count <- length(peak) - sum(df$count)
    if(intergenic_count < 0) intergenic_count <- 0
    intergenic <- data.frame("gene_type"= "intergenic", "count"= intergenic_count)

    selected <- c("protein_coding", "lincRNA", "antisense", "pseudogene", "snRNA", "snoRNA", "rRNA")
    selected_df <- filter(df, gene_type %in% selected)
    other_df <- filter(df, !gene_type %in% selected)
    other <- data.frame("gene_type"="other", "count"=sum(other_df$count))

    df <- rbind(selected_df, other, intergenic) %>%
      mutate(percent = round(count*100 / sum(count), 1)) %>%
      arrange(desc(count))

    pdf(paste0(peaklabel, "_gene_type_distribution_piechart.pdf"), height=8, width=10)
    p <- ggplot2::ggplot(df, aes(x = reorder(gene_type, -percent), y = percent)) +
      geom_bar(stat = 'identity', aes(fill = gene_type)) +
      geom_label(aes(y = percent + 0.5), label = df$count) +
      labs(x = '', y = 'Percent overlap (%)',
           title="Annotation of peaks to all type of genes",
           subtitle=paste0('(Total number of peaks = ', length(peak), ')')) +
      theme_bw(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust=0.5),
            plot.subtitle=element_text(hjust=0.5))

    print(p)


    # Plotting overlap counts between query regions and transcript features
    # here utr5, utr3 and cds do not contain introns
    protein_coding_gff <- gff[gff$gene_biotype == "protein_coding",]
    txdb <- makeTxDbFromGRanges(protein_coding_gff)


    utr5 <- gtf_to_bed_longest_tx(txdb, "utr5", longest=F)
    utr3 <- gtf_to_bed_longest_tx(txdb, "utr3", longest=F)
    cds <- gtf_to_bed_longest_tx(txdb, "cds", longest=F)
    intron <- gtf_to_bed_longest_tx(txdb, "intron", longest=F)

    promoter <- promoters(txdb, upstream=fiveP, downstream=100, use.names=FALSE)
    TTS <- promoters(txdb, upstream=100, downstream=threeP, use.names=FALSE)

    features <- GRangesList("Promoter"=promoter,
                            "TTS"=TTS,
                            "5'UTR"=stack(utr5$GRangesList),
                            "3'UTR"=stack(utr3$GRangesList),
                            "CDS"=stack(cds$GRangesList),
                            "Intron"=intron$GRanges, compress=F)
    if(RNA){
      features <- GRangesList("5'UTR"=stack(utr5$GRangesList),
                              "3'UTR"=stack(utr3$GRangesList),
                              "CDS"=stack(cds$GRangesList),
                              "Intron"=intron$GRanges, compress=F)
    }

    annot = annotateWithFeatures(peak, features, strand.aware=TRUE, intersect.chr=FALSE)
    precedence_count <- annot@num.precedence
    #precedence_count["NoFeature"] <- length(peak) - sum(precedence_count)

    #plotTargetAnnotation(annot)
    #summary <- summarizeQueryRegions(queryRegions = peak,
    #                                txdbFeatures = features)

    df <- data.frame(precedence_count) %>%
      mutate(feature = rownames(.)) %>%
      mutate(percent = precedence_count / sum(precedence_count)) %>%
      mutate(labels = scales::percent(percent, accuracy = 0.1))

    print(df)

    p <- ggpie(df, x="percent", label="labels", lab.pos="out", fill="feature", color="white", palette="ucscgb")
    print(p)

    df <- data.frame(precedence_count) %>%
      mutate(feature = rownames(.)) %>%
      filter(!feature %in% c("Intron", "NoFeature")) %>%
      mutate(percent = precedence_count / sum(precedence_count)) %>%
      mutate(labels = scales::percent(percent, accuracy = 0.1))

    p <- ggpie(df, x="percent", label="labels", lab.pos="out", fill="feature", color="white", palette="ucscgb")
    print(p)
    dev.off()

    # alternative plots
    if(0){
    ggplot2::ggplot(df, aes(x = reorder(feature, -percent), y = percent)) +
      geom_bar(stat = 'identity', aes(fill = feature)) +
      geom_label(aes(y = percent + 3), label = df$count) +
      labs(x = 'transcript feature', y = paste0('percent overlap (n = ', length(gr), ')')) +
      theme_bw(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90))


    ggplot(df, aes(x = "", y = percent, fill = feature)) +
      geom_col(color = "white") +
      geom_label(aes(label = labels), color = "white",
                 position = position_stack(vjust = 0.5),
                 show.legend = FALSE) +
      guides(fill = guide_legend(title = "Feature")) +
      scale_fill_viridis_d() +
      coord_polar(theta = "y") +
      theme_void()
    }
  }

  return(return_table)
}


#' process_scoreMatrix
#'
#' This is a helper function for manipulate the score matrix produced by ScoreMatrix or ScoreMatrinBin functions defined in the 'genomation' package.
#'
#' @param fullmatrix, a numeric matrix, with bins in columns and genomic windows in rows
#' @param libsize, sample library size in terms of read counts
#' @param norm, a boolean object indicating whether the matrix should be normalized to reads per million (RPM)
#' @param scale, a boolean indicates whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param heatmap, a boolean indicates whether a heatmap of the score matrix should be generated
#' @param rm.outlier, a boolean indicates whether a row with abnormally high values in the score matrix should be removed
#'
#' @return a numeric matrix
#'
#' @examples
#'
#' @export process_scoreMatrix
#'
#'
process_scoreMatrix <- function(fullmatrix, libsize=1000000, norm=F, scale=F, heatmap=F, rm.outlier=F){
  fullmatrix[is.na(fullmatrix)] <- 0

  if(norm && !is.null(libsize)){
    fullmatrix <- fullmatrix*1000000/libsize
  }
  if(scale){
    smc <- t(apply(fullmatrix, 1, scales::rescale)) ## rescale to 0:1 range
    fullmatrix <- as.data.frame(smc)
  }

  if(heatmap){
    fullmatrixm <- log10(fullmatrix+1)
    print(pheatmap(fullmatrixm, cluster_rows = T, cluster_cols = F, fontsize_col=8, fontsize_row=8, labels_col=colLabel, angle_col=0, main=paste(querylabel, centerlabel)))
  }


  ## remove outliers (highest score) from reference regions, using Hampel filter with 10mad instead of 3mad.
  ## if outliers are detected, remove the top outlier only
  if(rm.outlier){
    fullmatrix <- rmOutlier(fullmatrix)
    if(heatmap){
      fullmatrixm <- log10(fullmatrix+1)
      print(pheatmap(fullmatrixm, cluster_rows = T, cluster_cols = F, fontsize_col=8, fontsize_row=8, labels_col=colLabel, angle_col=0, main=paste(querylabel, centerlabel, "after removing outliers")))
    }

  }
  return(fullmatrix)
}

#' rmOutlier
#'
#' This is a helper function for removing rows with excessively high values, using Hampel filter with 10mad instead of 3mad.
#' If outliers are detected, remove the top outlier only, to minimize data loss.
#'
#' @param fullmatrix, a numeric matrix, with bins in columns and genomic windows in rows
#'
#' @return a numeric matrix
#'
#' @examples
#'
#' @export rmOutlier
#'

rmOutlier <- function(fullmatrix){
  fullmatrix[is.na(fullmatrix)] <- 0
  rowm <- apply(fullmatrix, 1, max)
  up_bound <- median(rowm) + 100*mad(rowm)
  if(length(which(rowm > up_bound)) > 0){
    print("Outlier detected:")
    outliers <- fullmatrix[which.max(rowm),]

    print(paste("median of row max", median(rowm), "up_bound", up_bound))
    print(outliers)

    fullmatrix <- fullmatrix[-which.max(rowm),]
  }
  return(fullmatrix)
}


aovTukeyHSD <- function(df){
  #sink(paste0(outfile,"two-way_anova_TukeyHSD.txt"))
  cat("Performing tow-way ANOVA analysis\n")
  res.aov2 <- aov(Intensity ~ Group * Position, data = df)
  print(summary(res.aov2)) ## anova summary
  cat("\nPost hoc Tukey Honest Significant Differences test\n")
  print(TukeyHSD(res.aov2, which = "Group"))

  #sink()

}


#' @param bedList, a named list of bed files, with length = 2 or 3
#' @param outPrefix, a string for plot file name
#'
#' @examples
#' setwd("C:/GREENBLATT/Nabeel/Gio/cits_peaks")
#' bedList <- list("combined_CIMS_0.001_YTHDF2_top10percent.bed", "combined_CITS_0.001_YTHDF2_top10percent.bed", "combined_CITS_0.001_ZNF121_top10percent.bed")
#' names(bedList) <- c("CIMS_YTHDF2", "CITS_YTHDF2", "CITS_ZNF121")
#' outPrefix <- "test"
#'
overlapBed <- function(bedList, outPrefix=NULL, fix_width=0L, fixPoint="center", pairOnly=TRUE, genome="hg19"){

  inputList <- lapply(bedList, handle_bed, fix_width=fix_width, fixPoint=fixPoint, useScore=FALSE, outRle=FALSE, genome=genome)
  grList <- lapply(inputList, function(x)x$query)
  sizeList <- lapply(inputList, function(x)x$size)

  overlap_pair <- function(apair){
    sizes <- sapply(apair, length)
    overlap <- length(Reduce(filter_by_overlaps_stranded, apair))
    venn.plot <- VennDiagram::draw.pairwise.venn(sizes[1], sizes[2], overlap, category=names(apair), lty=rep("blank",2), fill=c("#0020C2", "#64E986"),
                                    cat.just = rep(list(c(0.5, 0)),2), cex = rep(2, 3), cat.pos = c(0, 0))

    grid.draw(venn.plot)
    grid.newpage()
  }

  overlap_triple <- function(atriple){
    sizes <- sort(sapply(atriple, length), decreasing=T)
    atriple <- atriple[names(sizes)] ## sort the gr by decreasing size to avoid n13 < n123

    overlap12 <- length(filter_by_overlaps_stranded(atriple[[1]], atriple[[2]]))
    overlap13 <- length(filter_by_overlaps_stranded(atriple[[1]], atriple[[3]]))
    overlap23 <- length(filter_by_overlaps_stranded(atriple[[2]], atriple[[3]]))
    overlap123 <- length(filter_by_overlaps_stranded(filter_by_overlaps_stranded(atriple[[1]], atriple[[3]]), atriple[[2]]))

    venn.plot <- VennDiagram::draw.triple.venn(sizes[1], sizes[2], sizes[3], overlap12, overlap23, overlap13, overlap123, category=names(atriple),
                                               lty=rep("blank",3), fill=c("#0020C2", "#64E986", "#990012"),
                                                 cat.just = rep(list(c(0.5, 0)),3), cex = rep(2, 7), cat.pos = c(180, 180, 0))

    grid.draw(venn.plot)
    grid.newpage()
  }


  pairs <- combn(grList, 2, simplify = F)

  if(!is.null(outPrefix)){
    pdf(paste0(outPrefix, ".pdf"), width=8, height=8)

    lapply(pairs, overlap_pair)
    if(!pairOnly){
       triples <- combn(grList, 3, simplify = F)
       lapply(triples, overlap_triple)
    }
    dev.off()
  }
}

if(0){

#' @param
#' @exmaples
#'

countFile <- "all_factor.CITS_count.bed"
pcutoff <- 1e-10

system.time(count_table <- data.table::fread(countFile))
system.time(count_table <- as.data.frame(read.delim2(countFile, header=F)))
counts <- count_table[,5]
uniq_counts <- sort(unique(counts))

para_Poisson <- MASS::fitdistr(counts, "Poisson")
#para_NB <- MASS::fitdistr(counts, "negative binomial")

dpois(x=0, lambda=para_Poisson$estimate, log = F)
ppois(0.5, lambda=para_Poisson$estimate, lower.tail=T)
qpois(1e-10, lambda=para_Poisson$estimate, lower.tail=F)

dt <- data.frame(uniq_counts, dPoisson)
cv <- dt[dt[,2] <= pcutoff, ]
cv <- qpois(pcutoff, lambda=para_Poisson$estimate, lower.tail=F)
car::densityPlot(counts)
cv
hotspot <- count_table[count_table[,5] > cv,]

write.table(hotspot, "crosslink_hotspot.bed", col.names=F, row.names=F, sep="\t", quote=F)
write.table(dt,  "Poisson_fitted_pvalue.tab", col.names=T, row.names=F, sep="\t", quote=F)

}

filter_by_overlaps_stranded <- function(query, subject, maxgap=-1L){

  plus_query <- query[strand(query)=="+"]
  minus_query <- query[strand(query)=="-"]
  plus_subject <- subject[strand(subject)=="+"]
  minus_subject <- subject[strand(subject)=="-"]

  overlap_plus <- filter_by_overlaps(plus_query, plus_subject, maxgap=maxgap)
  overlap_minus <- filter_by_overlaps(minus_query, minus_subject, maxgap=maxgap)

  overlaps <- c(overlap_plus, overlap_minus)

  if(length(overlaps) > min(length(query), length(subject))){
     warning("Size of overlap is greater than min(sizeOfQuery, sizeOfSubject!")
  }
  overlaps
}

filter_by_nonoverlaps_stranded <- function(query, subject, maxgap=-1L){

  overlaps <- filter_by_overlaps_stranded(query, subject, maxgap=-1L)
  if(length(overlaps) > min(length(query), length(subject))){
     warning("Size of overlap is greater than min(sizeOfQuery, sizeOfSubject!")
  }
  nonoverlaps <- query[!query %in% overlaps]
  nonoverlaps
}

prepare_5parts_genomic_features <- function(txdb, longest=TRUE, meta=TRUE, nbins=100, useIntron=FALSE, fiveP=1000, threeP=1000){
   ## prepare transcripts that are suitable for overlap
   print("Preparing genomic features ... ")
   featureNames <-  c("promoter", "utr5", "cds", "utr3", "TTS")
   if(useIntron) featureNames <-  c("promoter", "utr5", "cds", "utr3", "Intron")
   tl <- transcriptLengths(txdb, with.utr5_len = T, with.cds_len = T, with.utr3_len = T)
   tl <- tl[tl$cds_len > 0,]  ## choose protein coding genes only

   ## selecte the longest transcript for each gene
   if(longest){
      longest_tx <- extract_longest_tx(txdb)
      tl <- tl[tl$tx_name %in% longest_tx$tx_name, ]
   }

   utr5 <- gtf_to_bed_longest_tx(txdb, "utr5", longest=longest)
   utr3 <- gtf_to_bed_longest_tx(txdb, "utr3", longest=longest)
   cds <- gtf_to_bed_longest_tx(txdb, "cds", longest=longest)
   if(useIntron) intron <- gtf_to_bed_longest_tx(txdb, "intron", longest=longest)

   if(meta){
      means <- c(promoter=fiveP, apply(tl[,c(7,6,8)], 2, median), TTS=threeP)
      scaled_bins <- round(means*nbins/sum(means))
      names(scaled_bins) <- featureNames

      tl_selected <- tl[tl$utr5_len >= scaled_bins["utr5"] & tl$cds_len >= scaled_bins["cds"] & tl$utr3_len >= scaled_bins["utr3"], ]

      print("median sizes for features")
      print(means)
      print("bin sizes for features")
      print(scaled_bins)
      print("Number of genes included in analysis")
      print(nrow(tl_selected))

      gene <- GenomicFeatures::genes(txdb)
      gene <- gene[gene$gene_id %in% tl_selected$gene_id,]

      promoter <- flank(gene, width=fiveP, both=F, start=T, ignore.strand=FALSE)
      TTS <- flank(gene, width=threeP, both=F, start=F, ignore.strand=FALSE)

      windowRs <- list(split(promoter, as.factor(promoter$gene_id)),
                       utr5$GRangesList[as.character(tl_selected$tx_id)],
                       cds$GRangesList[as.character(tl_selected$tx_id)],
                       utr3$GRangesList[as.character(tl_selected$tx_id)],
                       split(TTS, as.factor(TTS$gene_id)))
      if(useIntron){
         intronGr <- intron$GRanges[width(intron$GRanges) >= scaled_bins["Intron"]]
         windowRs <- list(split(promoter, as.factor(promoter$gene_id)),
                          utr5$GRangesList[as.character(tl_selected$tx_id)],
                          cds$GRangesList[as.character(tl_selected$tx_id)],
                          utr3$GRangesList[as.character(tl_selected$tx_id)],
                          split(intronGr, as.factor(seq_along(intronGr))))
      }
   }else{
      if(0){
         cl <- start_parallel(nc=20)
         print("getting utr5 GRanges")
         utr5grl <- parLapply(cl, utr5$GRangesList[names(utr5$GRangesList) %in% as.character(tl$tx_id)], range)
         print("getting cds GRanges")
         cdsgrl <- parLapply(cl, cds$GRangesList[names(cds$GRangesList) %in% as.character(tl$tx_id)], range)
         print("getting utr3 GRanges")
         utr3grl <- parLapply(cl, utr3$GRangesList[names(utr3$GRangesList) %in% as.character(tl$tx_id)], range)
         stop_parallel(cl)
      }

      utr5gr <- utr5$GRanges[utr5$GRanges$name %in% as.character(tl$tx_id)]
      utr5grl <- split(utr5gr, as.factor(utr5gr$name))
      utr3gr <- utr3$GRanges[utr3$GRanges$name %in% as.character(tl$tx_id)]
      utr3grl <- split(utr3gr, as.factor(utr3gr$name))
      cdsgr <- cds$GRanges[cds$GRanges$name %in% as.character(tl$tx_id)]
      cdsgrl <- split(cdsgr, as.factor(cdsgr$name))

      l3 <- sapply(list(utr5grl, cdsgrl, utr3grl), function(x)median(sapply(x, width)))

      means <- c(promoter=fiveP, l3, TTS=threeP)
      scaled_bins <- round(means*nbins/sum(means))
      names(scaled_bins) <- featureNames

      grls <- list("utr5"=utr5grl, "cds"=cdsgrl, "utr3"=utr3grl)
      selected_tx <- lapply(names(grls), function(x){
         len <- sapply(grls[[x]], width) # len is named vector, where names are the tx_ids
         nbin <- scaled_bins[x]
         y <- names(len[which(len >= nbin)])
      })

      selected_tx <- Reduce(intersect, selected_tx)
      selected_gene <- filter(tl, tx_id %in% selected_tx) %>% select(gene_id)

      print("median sizes for features")
      print(means)
      print("bin sizes for features")
      print(scaled_bins)
      print("Number of genes included in analysis")
      print(nrow(selected_gene))


      gene <- genes(txdb)
      gene <- gene[gene$gene_id %in% selected_gene$gene_id,]

      promoter <- flank(gene, width=fiveP, both=F, start=T, ignore.strand=FALSE)
      TTS <- flank(gene, width=threeP, both=F, start=F, ignore.strand=FALSE)

      windowRs <- list(split(promoter, as.factor(promoter$gene_id)),
                       as(utr5grl[selected_tx], "GRangesList"),
                       as(cdsgrl[selected_tx], "GRangesList"),
                       as(utr3grl[selected_tx], "GRangesList"),
                       split(TTS, as.factor(TTS$gene_id)))
   }

   names(windowRs) <- featureNames
   print(lapply(windowRs, length))

   return(list("windowRs"=windowRs, "nbins"=nbins, "scaled_bins"=scaled_bins, "fiveP"=fiveP, "threeP"=threeP, "useIntron"=useIntron, "meta"=meta, "longest"=longest))
}

get_genomicCoordinates <- function(x){
   if(grepl("GRangesList", class(x))) x <- stack(x) ## convert grl to gr

   chr <- seqnames(x) %>% as.character
   start <- start(x) %>% as.numeric
   end <- end(x) %>% as.numeric
   strand <- strand(x) %>% as.character

   out <- paste0(chr, ":", start, "-", end, "(", strand, ")")
   out
}

peak_targeted_gene <- function(peakfile, gtffile, genome="hg19", filterfile=NULL, maxg=-1L){
   if(0){
      peakfile <- "combined_CITS_0.01_m6AHek.merged_DRACH_slopl2r2.bed"
      names(peakfile) <- "Hek_m6A"
      genome <- "hg19"
      gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
      RNA <- T

      filterfile <- "C:/GREENBLATT/Nabeel/ZBTB48/ZBTB48_clip/combined_crosslink_0.01_ZBTB48.recurring.bed"
      names(filterfile) <- "ZBTB48_CLIP"
   }

   return_table <- NULL
   peaklabel <- names(peakfile)
   peak <- handle_bed(inputFile=peakfile, fix_width=0, useScore=FALSE, outRle=FALSE, genome=genome)$query
   if(!is.null(filterfile)){
      filterlabel <- names(filterfile)
      filter <- handle_bed(inputFile=filterfile, fix_width=0, useScore=FALSE, outRle=FALSE, genome=genome)$query
   }

   gff <- importGtf(filePath = gtffile)
   txdb <-  makeTxDbFromGFF(gtffile)

   overlaps <- as.data.table(queryGff(queryRegions=peak, gffData=gff))

   gene_info_table <- overlaps[, c("transcript_id", "gene_id", "gene_name", "gene_biotype")]

   gene <- gtf_to_bed_longest_tx(txdb, "gene", longest=FALSE)
   geneGr <- gene$GRanges

   peak_list <- list()
   if(!is.null(filterfile)){
      overlapped_peak <- filter_by_overlaps_stranded(peak, filter, maxgap=maxg)
      nonoverlapped_peak <- filter_by_nonoverlaps_stranded(peak, filter, maxgap=maxg)
      peak_list[[paste(filterlabel,"overlapped", peaklabel, sep="_")]] <- overlapped_peak
      peak_list[[paste(filterlabel,"nonoverlapped", peaklabel, sep="_")]] <- nonoverlapped_peak
   }else{
      peak_list[[peaklabel]] <- peak
   }

   return_list <- lapply(names(peak_list), function(aname){
      peak <- peak_list[[aname]]
      print(aname)
      print(length(peak))
      geneOverlaps <- GenomicRanges::findOverlaps(peak, geneGr)
      peak_df <- annoGR2DF(peak[geneOverlaps@from]) %>%
         mutate(chrPeak=as.character(chr),
                startPeak=as.integer(start)-1,
                endPeak=as.integer(end),
                widthPeak=as.integer(width),
                strandPeak=as.character(strand),
                .keep="unused")
      gene_df <- geneGr[geneOverlaps@to]
      names(gene_df) <- seq_along(gene_df)
      gene_df <- annoGR2DF(gene_df) %>%
         mutate(chrGene=as.character(chr),
                startGene=as.integer(start)-1,
                endGene=as.integer(end),
                widthGene=as.integer(width),
                strandGene=as.character(strand),
                .keep="unused")

      ot <- cbind(peak_df, gene_df) %>%
         select(chrPeak, startPeak, endPeak, id, widthPeak, strandPeak, chrGene, startGene, endGene, gene_id, widthGene, strandGene)
      ot_gene <- merge(ot, gene_info_table, by.x="gene_id", by.y="gene_id", all.x=T) %>%
         select(chrPeak, startPeak, endPeak, id, widthPeak, strandPeak, chrGene, startGene, endGene, gene_id, widthGene, strandGene, gene_name, gene_biotype) %>%
         unique() %>%
         mutate(overlap_status=aname)

      print(length(unique(ot_gene$gene_name)))
      return(ot_gene)
   })


   return_table <- bind_rows(return_list)
   write.table(return_table, paste(peaklabel, "overlap", filterlabel, "targeted_gene.tab", sep="_"), sep="\t", row.names=F, quote=F)

   return(return_table)
}


peak_targeted_enhancer <- function(peakfile, enhancerfile, enhancerGenefile, maxg=-1L, genome="hg19"){
   if(0){
      peakfile <- queryfiles
      enhancerfile <- "C:/GREENBLATT/resource/EnhancerAtlas2.0/HEK293.bed"
      enhancerGenefile <- "C:/GREENBLATT/resource/EnhancerAtlas2.0/HEK293_EP.txt"
      maxg <- -1L
      genome <- "hg19"
   }

   EG <- read.delim(enhancerGenefile, header=F)
   EG_list <- lapply(EG[,1], function(x){
      str1 <- unlist(strsplit(x, split=":|_|\\$"))
      return(c(str1[1], unlist(strsplit(str1[2], split="-")), str1[3:7]))
   })
   names(EG_list) <- paste0("EG", seq_along(EG_list))
   EG_df <- as.data.frame(do.call(rbind, EG_list))
   colnames(EG_df) <- c("chrEnhancer", "startEnhancer", "endEnhancer", "geneId", "geneName", "chrGene", "coordGene", "strandGene")

   return_table <- NULL
   peaklabel <- names(peakfile)
   peak <- handle_bed(inputFile=peakfile, fix_width=0, useScore=FALSE, outRle=FALSE, genome=genome)$query
   enhancer <- handle_bed(inputFile=enhancerfile, fix_width=0, useScore=FALSE, outRle=FALSE, genome=genome)$query

   print(length(peak))
   print(length(enhancer))
   Overlaps <- GenomicRanges::findOverlaps(peak, enhancer)
   peak_gr <- peak[Overlaps@from]
   peak_df <- annoGR2DF(peak_gr) %>%
      mutate(chrPeak=as.character(chr),
             startPeak=as.integer(start)-1,
             endPeak=as.integer(end),
             widthPeak=as.integer(width),
             strandPeak="*",
             .keep="unused")
   enhancer_gr <- enhancer[Overlaps@to]
   names(enhancer_gr) <- seq_along(enhancer_gr)
   enhancer_df <- annoGR2DF(enhancer_gr) %>%
      mutate(chrEnhancer=as.character(chr),
             startEnhancer=as.integer(start)-1,
             endEnhancer=as.integer(end),
             widthEnhancer=as.integer(width),
             strandEnhancer="*",
             .keep="unused")

   ot <- cbind(peak_df, enhancer_df) %>%
      select(chrPeak, startPeak, endPeak, widthPeak, strandPeak, chrEnhancer, startEnhancer, endEnhancer, widthEnhancer, strandEnhancer)

   return_table <- merge(ot, EG_df)
   return_table <- return_table[,colnames(return_table)[c(4:8, 1:3, 9:15)]]
   uni_peaks <- unique(return_table[, 1:5])
   uni_enhancers <- unique(return_table[, 6:10])
   uni_genes <- unique(return_table[11:15])

   fraction_peaks <- nrow(uni_peaks)/length(peak)
   fraction_enhancers <- nrow(uni_enhancers)/length(enhancer)

   sink(paste(peaklabel, "targeted_enhancer_gene_stats.tab", sep="_"))
   cat(paste0("GenomicPlot.R::peak_targeted_enhancer()\t",date()))
   cat(paste0("\nSummary stats for\t", peaklabel, "_targeted_enhancer_gene.tab"))
   cat(paste0("\nTotal number of summit peaks\t", length(peak)))
   cat(paste0("\nTotal number of enhancers in HEK293\t", length(enhancer)))
   cat(paste0("\nNumber of peaks targeting enhancer\t", nrow(uni_peaks)))
   cat(paste0("\nFraction of peaks targeting enhancer\t", fraction_peaks))
   cat(paste0("\nNumber of enhancers targeted by peaks\t", nrow(uni_enhancers)))
   cat(paste0("\nFraction of enhancer targeted by peaks\t", fraction_enhancers))
   cat(paste0("\nNumber of genes associated with targeted enhancer\t", nrow(uni_genes)))
   sink()

   write.table(return_table, paste(peaklabel, "targeted_enhancer_gene.tab", sep="_"), sep="\t", row.names=F, quote=F)

   return(return_table)
}
