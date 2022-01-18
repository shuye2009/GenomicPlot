
#' @import methods
NULL

#' setup
#'
#' Method for load all packages required by other methods
#'
#' @param NULL
#' @return NULL
setup <- function(){
  list.of.packages <- c(
    "data.table",
    "MatrixGenerics",
    "dplyr",
    "DESeq2",
    "magrittr",
    "cowplot",
    "rtracklayer",
    "TxDb.Hsapiens.UCSC.hg19.knownGene",
    "BSgenome.Hsapiens.UCSC.hg19",
    "genomation",
    "GenomicRanges",
    "plyranges",
    "GenomicFeatures",
    "GenomicAlignments",
    "gridExtra",
    "ggplot2",
    "ggsignif",
    "ggpubr",
    "R.utils",
    "hrbrthemes",
    "pheatmap",
    "scales",
    "RMariaDB"
  )
  for(package.i in list.of.packages){
    suppressPackageStartupMessages(
      library(package.i, character.only = TRUE)
    )
  }
}

setup()

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
#' Make bed and bed 12 files from gtf file for protein coding genes, exon.bed12 will the transcript.bed12. If 'longest'=TRUE, only the
#' longest transcript of each gene is considered, thus the exons, introns, 3'utr and 5'utr are limited to the longest transcript of
#' protein-coding genes.
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
#' output <- gtf_to_bed_longest_tx(txdb, featureName="utr3", feature_source="Ensembl", export=F, longest=T)
#' @export gtf_to_bed_longest_tx
#'
#'
gtf_to_bed_longest_tx <- function(txdb, featureName, featureSource=NULL, export=FALSE, longest=TRUE){

 # featureName <- "transcript"
  longest_tx <- extract_longest_tx(txdb)
  if(!is.null(feature_source)){featureSource <- gsub("\\.gtf", "", feature_source)}
  feature <- NULL

  if(featureName == "utr3"){
    feature <- threeUTRsByTranscript(txdb, use.name=F)
  }else if(featureName == "utr5"){
    feature <- fiveUTRsByTranscript(txdb, use.name=F)
  }else if(featureName == "intron"){
    feature <- intronsByTranscript(txdb, use.name=F)
    tl <- transcriptLengths(txdb, with.utr5_len = T, with.cds_len = T, with.utr3_len = T)
    tl_protein_coding <- tl[tl$cds_len > 0, ]
    feature <- feature[tl_protein_coding$tx_id]
  }else if(featureName == "exon"){
    feature <- exonsBy(txdb, by="tx", use.name=F)
    tl <- transcriptLengths(txdb, with.utr5_len = T, with.cds_len = T, with.utr3_len = T)
    tl_protein_coding <- tl[tl$cds_len > 0, ]
    feature <- feature[tl_protein_coding$tx_id]
  }else if(featureName == "cds"){
    feature <- cdsBy(txdb, by="tx", use.name=F)
  }else if(featureName == "transcript"){
    feature <- transcripts(txdb, use.name=F)
    feature <- split(feature, f=feature$tx_id)
    tl <- transcriptLengths(txdb, with.utr5_len = T, with.cds_len = T, with.utr3_len = T)
    tl_protein_coding <- tl[tl$cds_len > 0, ]
    feature <- feature[tl_protein_coding$tx_id]
  }else if(featureName == "gene"){
    feature <- genes(txdb)
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
    if(export){export.bed(stack(feature_longest), paste(featureSource, "_", featureName, "_longest.bed", sep=""))}

    if(featureName %in% c("cds", "utr5", "utr3", "exon") && export){
      bed12_feature_longest <- asBED(feature_longest) ## convert Grangeslist to Granges with blocks info as metadata
      export.bed(bed12_feature_longest, paste(featureSource, "_", featureName, "_longest.bed12", sep=""))
    }
    return(list("GRanges"=stack(feature_longest), "GRangesList"=feature_longest))
  }else{
    if(export){export.bed(stack(feature), paste(featureSource, "_", featureName, "_all.bed", sep=""))}

    if(featureName %in% c("cds", "utr5", "utr3", "exon") && export){
      bed12_feature_all <- asBED(feature) ## convert Grangeslist to Granges with blocks info as metadata
      export.bed(bed12_feature_all, paste(featureSource, "_", featureName, "_all.bed12", sep=""))
    }
    return(list("GRanges"=stack(feature), "GRangesList"=feature))

  }

}

#' plot_start_end_feature
#
#' Plot reads or peak signal intensity of samples in the query files around stat and end of genomic features. The upstream and downstream windows
#' can be given separately, within the window, a smaller window can be defined to test statistical difference between the signal intensity
#' among samples. A line plot and a boxplot are displayed side by side for both start and end of feature.
#'
#' @param queryfiles, a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param querylabels, a vector of charactor strings serving as short labels of the queryfiles, will be used as sample labels in the plots
#' @param txdb, a TxDb object defined in GenomicFeatures package
#' @param feaureNmae, one of the gene feature in c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene")
#' @param CLIP_reads, a boolean object to indicate if the bam reads should be shifted to the -1 positon at the 5' of the reads
#' @param binsize, an integer defines bin size for intensity calculation
#' @param extend_reads, an integer defines how long should the reads should be extended to
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
#' @param outprefix, a character string specifying output file prefix for plots (outprefix.pdf)
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
#' plot_start_end_feature(queryfiles=queryfiles, querylabels=querylabels, inputfiles=inputfiles, inputlabels=inputlabels, txdb=txdb, featureName="intron", CLIP_read=F, binsize=10, extend_reads=0,
#' longest=T, ext=ext, hl=hl, randomize=T, stranded=T, norm=F, scale=F, smo=F, heatmap=F, rm.outlier=F, genome="hg19", outprefix=op)
#'
#' @export plot_start_end_feature
#'
#'

plot_start_end_feature <- function(queryfiles=NULL, querylabels=NULL, inputfiles=NULL, inputlabels=NULL, txdb=NULL, featureName=NULL, CLIP_reads=F, binsize=10, extend_reads=0,
                               longest=T, ext=c(0,0,0,0), hl=c(0,0,0,0), randomize=FALSE, stranded=T, norm=F, scale=F, smo=F, heatmap=F, rm.outlier=F, genome="hg19", outprefix="plots"){

  names(querylabels) <- queryfiles
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


  bedInputs <- handle_input(queryfiles, CLIP_reads=CLIP_reads, extend_reads=extend_reads,  useScore=TRUE, outRle=TRUE, norm=norm, genome=genome)

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
      fullmatrix <- process_scoreMatrix(fullmatrix, libsize, norm, scale, heatmap, rm.outlier=is.null(inputfiles))

      scoreMatrix_list[[querylabel]][[locus]] <- fullmatrix

      if(randomize){

        rfullmatrix <- parallel_scoreMatrixBin(queryRegions, rwindowR, bin_num, bin_op, weight_col, stranded)
        rfullmatrix <- process_scoreMatrix(rfullmatrix, libsize, norm, scale, heatmap, rm.outlier=is.null(inputfiles))
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
      sub_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=querybed, "Location"=location)
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
        rsub_df <- data.frame("Intensity"=rcolm, "sd"=rcolsd, "se"=rcolse, "Position"=rcollabel, "Query"=rquerybed, "Location"=location)
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
    theme_classic() + theme(legend.position="top") + xlab("") + ylab("Signal intensity") +
    ggtitle(featureName) +
    facet_wrap(~Location, scales="free_x")
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
        sub_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=ratiobed, "Location"=location)
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
          rsub_df <- data.frame("Intensity"=rcolm, "sd"=rcolsd, "se"=rcolse, "Position"=rcollabel, "Query"=rratiobed, "Location"=location)
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
      facet_wrap(~Location, scales="free_x")
    print(p)
  }
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
#' @param extend_reads, an integer defines how long should the reads should be extended to
#' @param norm, an boolean indicates if the intensities should be normalized by the sample library sizes
#' @param longest, a boolean object to indicate whether the output should be limited to the longest transcript of each gene
#' @param stranded, a boolean indicate whether the strand of the feature should be considered
#' @param scale, a boolean indicates whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smo, a boolean indicates whether the line should smoothed with a spline smoothing algorithm
#' @param heatmap, a boolean indicates whether a heatmap of the score matrix should be generated
#' @param rm.outlier, a boolean indicates whether a row with abnormally high values in the score matrix should be removed
#' @param genome, genome of the features, the program is mostly concerned with human, if non-human genome is used, certain features may not work
#' @param outprefix, a character string specifying output file prefix for plots (outprefix.pdf)
#'
#' @return a list of two objects, the first is a GRanges object, the second is a GRangesList object
#'
#' @examples
#'
#'
#' @export plot_3parts_metagene


plot_3parts_metagene <- function(queryfiles=NULL, querylabels=NULL, txdb=NULL, meta=T, inputfiles=NULL, inputlabels=NULL, nbins=100, norm=FALSE, longest=TRUE, scale=FALSE,
                                 extend_reads=0, CLIP_reads=FALSE, fiveP=1000, threeP=1000, smo=FALSE, stranded=TRUE, outPrefix="plots", genome="hg19",
                                 heatmap=FALSE, rm.outlier=FALSE){

  if(0){
    nbins <- 100
    fiveP <- 500
    threeP <- 500
    norm=FALSE
    longest=TRUE
    CLIP_reads=FALSE
    extend_reads=0
    smo=FALSE
    stranded=TRUE
    outPrefix=op
    genome="tetrahymena"
  }


  names(querylabels) <- queryfiles
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

  bedInputs <- handle_input(queryfiles, CLIP_reads=CLIP_reads, extend_reads=extend_reads, useScore=TRUE, outRle=TRUE, norm=norm, genome=genome)

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

      fullmatrix <- process_scoreMatrix(fullmatrix, libsize=libsize, norm=norm, scale=scale, heatmap=heatmap, rm.outlier=is.null(inputfiles))
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


## plot scaled regions for centers given in bed format

plot_reference_region <- function(queryfiles=NULL, querylabels=NULL, centerfiles, centerlabels, nbins=100, norm=FALSE, extend_reads=0, CLIP_reads=FALSE, fiveP=1000, threeP=1000, smo=FALSE, stranded=TRUE, outPrefix=NULL, genome="hg19"){

  if(0){
    nbins <- 100
    fiveP <- 500
    threeP <- 500
    norm=FALSE
    longest=TRUE
    CLIP_reads=FALSE
    extend_reads=0
    smo=FALSE
    stranded=TRUE
    outPrefix=op
    genome="tetrahymena"
  }


  if(!is.null(outPrefix)){
    pdf(paste(outPrefix, "pdf", sep="."), height=8, width=12)
  }

  names(querylabels) <- queryfiles
  names(centerlabels) <- centerfiles

  ## prepare transcripts that are suitable for overlap
  featureNames <-  c("upstream", "region", "downstream")

  means <- c(upstream=fiveP, region=fiveP+threeP, downstream=threeP)
  scaled_bins <- round(means*nbins/sum(means))
  names(scaled_bins) <- featureNames


  mplot_df <- NULL
  vx <- NULL ## x axis point for vlines that demarcate the genomic features
  color_store <- c("#00AFBB", "#E7B800", "#A0BDE0", "#0020C2", "#64E986", "#F5DEB3", "#C19A6B", "#E8A317", "#8E7618", "#A0522D", "#990012", "#CB6D51")

  print("computing coverage for Sample")
  scoreMatrix_list <- list()

  bedInputs <- handle_input(queryfiles, CLIP_reads=CLIP_reads, extend_reads=extend_reads, useScore=TRUE, outRle=TRUE, norm=norm, genome=genome)
  centerInputs <- handle_input(centerfiles, CLIP_reads=FALSE, extend_reads=0, useScore=FALSE, outRle=FALSE, norm=FALSE, genome=genome)

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

      windowRs <- list(upstream=windowUp, region=windowRegions, downstream=windowDown)

      from <- to <- 0
      plot_df <- NULL
      #feature_matrix <- NULL, for scaling across row, still not working, as the smc nrows vary between feature, and cbind fails

      for(w in names(windowRs)){
        #w <- "utr5"
        print(w)
        windowR <- windowRs[[w]]
        bin_num <- scaled_bins[w]
        from <- to + 1
        to <- to + bin_num
        if(length(vx) < length(windowRs)) vx <- c(vx, from)

        #smc = ScoreMatrixBin(target = queryRegions, windows = windowR, bin.num=bin_num, bin.op="mean", weight.col=weight_col, strand.aware = stranded, is.noCovNA=T)

        #fullmatrix <- cov_lst[[w]]

        bin_op <- "mean"
        fullmatrix <- parallel_scoreMatrixBin(queryRegions, windowR, bin_num, bin_op, weight_col, stranded)

        fullmatrix[is.na(fullmatrix)] <- 0
        #feature_matrix <- coverageMatrixBin(queryRegions, windowR, nbins=bin_num, binOp="mean")

        print(dim(fullmatrix))
        ## rpm did not work for the RNAseq bam files, so normalization is performed here
        if(norm && !is.null(libsize)){
          fullmatrix <- fullmatrix/(libsize/1000000)
        }


        colm <- apply(fullmatrix, 2, mean)
        colsd <- apply(fullmatrix, 2, sd)
        colse <- colsd/sqrt(nrow(fullmatrix))
        collabel <- seq(from, to)
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
#' @param extend_reads, an integer defines how long should the reads should be extended to
#' @param norm, an boolean indicates if the intensities should be normalized by the sample library sizes
#' @param longest, a boolean object to indicate whether the output should be limited to the longest transcript of each gene
#' @param stranded, a boolean indicate whether the strand of the feature should be considered
#' @param scale, a boolean indicates whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smo, a boolean indicates whether the line should smoothed with a spline smoothing algorithm
#' @param heatmap, a boolean indicates whether a heatmap of the score matrix should be generated
#' @param rm.outlier, a boolean indicates whether a row with abnormally high values in the score matrix should be removed
#' @param genome, genome of the features, the program is mostly concerned with human, if non-human genome is used, certain features may not work
#' @param outprefix, a character string specifying output file prefix for plots (outprefix.pdf)
#'
#' @return a list of two objects, the first is a GRanges object, the second is a GRangesList object
#'
#' @examples
#'
#'
#' @export plot_5parts_metagene


plot_5parts_metagene <- function(queryfiles=NULL, querylabels=NULL, txdb=NULL, inputfiles=NULL, inputlabels=NULL, meta=TRUE, nbins=100, norm=FALSE,
                                 longest=TRUE, extend_reads=0, CLIP_reads=FALSE, fiveP=1000, threeP=1000, smo=FALSE, scale=FALSE, stranded=TRUE,
                                 outPrefix=NULL, genome="hg19", heatmap=F, rm.outlier=F){

  names(querylabels) <- queryfiles
  if(!is.null(inputfiles)){
    if(length(queryfiles) != length(inputfiles)){
      stop("the number of queryfiles and inputfiles must be the same!")
    }
    names(inputlabels) <- inputfiles
    queryfiles <- c(queryfiles, inputfiles)
    querylabels <- c(querylabels, inputlabels)
  }

  ## prepare transcripts that are suitable for overlap
  featureNames <-  c("promoter", "utr5", "cds", "utr3", "TTS")
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

    gene <- genes(txdb)
    gene <- gene[gene$gene_id %in% tl_selected$gene_id,]

    promoter <- flank(gene, width=fiveP, both=F, start=T, ignore.strand=FALSE)
    TTS <- flank(gene, width=threeP, both=F, start=F, ignore.strand=FALSE)

    windowRs <- list("promoter"=split(promoter, as.factor(promoter$gene_id)),
                     "utr5"=utr5$GRangesList[as.character(tl_selected$tx_id)],
                     "cds"=cds$GRangesList[as.character(tl_selected$tx_id)],
                     "utr3"=utr3$GRangesList[as.character(tl_selected$tx_id)],
                     "TTS"=split(TTS, as.factor(TTS$gene_id)))
  }else{
    utr5grl=lapply(utr5$GRangesList[names(utr5$GRangesList) %in% as.character(tl$tx_id)], range)
    cdsgrl=lapply(cds$GRangesList[as.character(tl$tx_id)], range)
    utr3grl=lapply(utr3$GRangesList[names(utr3$GRangesList) %in% as.character(tl$tx_id)], range)

    l3 <- sapply(list(utr5grl, cdsgrl, utr3grl), function(x)median(sapply(x, width)))

    means <- c(promoter=fiveP, l3, TTS=threeP)
    scaled_bins <- round(means*nbins/sum(means))
    names(scaled_bins) <- featureNames

    grls <- list("utr5"=utr5grl, "cds"=cdsgrl, "utr3"=utr3grl)
    selected_tx <- lapply(names(grs), function(x){
      len <- sapply(grls[[x]], width) # len is named vector, where names are the tx_ids
      nbin <- scaled_bins[x]
      y <- names(len[which(len >= nbin)])
    })

    selected_tx <- Reduce(intersect, selected_tx)
    selected_gene <- filter(tl, tx_id %in% selected_tx) %>% select(gene_id)

    gene <- genes(txdb)
    gene <- gene[gene$gene_id %in% selected_gene$gene_id,]

    promoter <- flank(gene, width=fiveP, both=F, start=T, ignore.strand=FALSE)
    TTS <- flank(gene, width=threeP, both=F, start=F, ignore.strand=FALSE)

    windowRs <- list("promoter"=split(promoter, as.factor(promoter$gene_id)),
                     "utr5"=as(utr5grl[selected_tx], "GRangesList"),
                     "cds"=as(cdsgrl[selected_tx], "GRangesList"),
                     "utr3"=as(utr3grl[selected_tx], "GRangesList"),
                     "TTS"=split(TTS, as.factor(TTS$gene_id)))
  }


  print(lapply(windowRs, length))
  ## start overlapping

  scoreMatrix_list <- list()

  bedInputs <- handle_input(inputFiles=queryfiles, CLIP_reads=CLIP_reads, extend_reads=extend_reads, useScore=TRUE, outRle=TRUE, norm=norm, genome=genome)

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
      fullmatrix <- parallel_scoreMatrixBin(queryRegions, windowR, bin_num, bin_op, weight_col, stranded)
      fullmatrix <- process_scoreMatrix(fullmatrix, libsize=libsize, norm=norm, scale=scale, heatmap=heatmap, rm.outlier=is.null(inputfiles))

      scoreMatrix_list[[querylabel]][[w]] <- fullmatrix
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

    if(smo){
      plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df=as.integer(nbins/5))$y)
      plot_df$se <- as.vector(smooth.spline(plot_df$se, df=as.integer(nbins/5))$y)
    }

    mplot_df <- rbind(mplot_df, plot_df)

  }

  mplot_df <- mutate(mplot_df, lower=Intensity-se, upper=Intensity+se)


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

  five <- fiveP/1000
  five <- paste0(five, "kb")
  three <- threeP/1000
  three <- paste0(three, "kb")
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

  if(!is.null(outPrefix)){
    pdf(paste(outPrefix, "pdf", sep="."), height=8, width=12)
  }

  for(i in seq_along(querylabels)){
    for(beds in combn(querylabels, i, simplify = F)){
      print(beds)

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
        ylab("Signal ntensity")

      outp <- plot_grid(p, pp, ppp, ncol = 1, align = 'v', axis="b", rel_heights = c(20,1,1))
      print(outp)
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

      for(i in seq_along(ratiolabels)){
        fullmatrix <- (ratioMatrix_list[[ratiolabels[i]]][[w]] + pseudo_count)/(inputMatrix_list[[inputlabels[i]]][[w]] + pseudo_count)

        if(rm.outlier){
          fullmatrix <- rmOutlier(fullmatrix)
        }

        ratioMatrix_list[[ratiolabels[i]]][[w]] <- fullmatrix
      }
    }

    mplot_df <- NULL
    for(ratiolabel in ratiolabels){
      plot_df <- NULL
      for(w in featureNames){
        #w <- "utr5"
        print(w)
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

      if(smo){
        plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df=as.integer(nbins/5))$y)
        plot_df$se <- as.vector(smooth.spline(plot_df$se, df=as.integer(nbins/5))$y)
      }

      mplot_df <- rbind(mplot_df, plot_df)

    }

    mplot_df <- mutate(mplot_df, lower=Intensity-se, upper=Intensity+se)

    for(i in seq_along(ratiolabels)){
      for(beds in combn(ratiolabels, i, simplify = F)){
        print(beds)

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
          ylab("Ratio_over_Input")

        outp <- plot_grid(p, pp, ppp, ncol = 1, align = 'v', axis="b", rel_heights = c(20,1,1))
        print(outp)
      }
    }
  }

  if(!is.null(outPrefix)){
    dev.off()
  }

  return(mplot_df)
}


## bed files can be reads (bam), coverage(wig, bw) or peaks (bed6 only)
## center files must be bed6
## if CLIP_reads, use -1 position of 5 prime of the reads

plot_reference_locus <- function(queryfiles=NULL,
                            centerfiles=NULL,
                            ext=c(0,0),
                            hl=c(0,0),
                            querylabels=NULL,
                            centerlabels=NULL,
                            smo=FALSE,
                            CLIP_reads=FALSE,
                            extend_reads=0,
                            norm=FALSE,
                            binsize=10,
                            refPoint="center",
                            Xlab="Center",
                            inputfiles=NULL,
                            inputlabels=NULL,
                            stranded=TRUE,
                            stats=FALSE,
                            heatmap=FALSE,
                            scale=FALSE,
                            outPrefix=NULL,
                            genome="hg19",
                            rm.outlier=F){


  names(querylabels) <- queryfiles
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

  bedInputs <- handle_input(queryfiles, CLIP_reads=CLIP_reads, extend_reads=extend_reads,  useScore=TRUE, outRle=TRUE, norm=norm, genome=genome)
  centerInputs <- handle_input(centerfiles, CLIP_reads=FALSE, extend_reads=0, useScore=FALSE, outRle=FALSE, norm=FALSE, genome=genome)

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
      ## if outliers are detected, remove the top outliers only
      if(is.null(inputfiles) && rm.outlier){
        rowm <- apply(fullmatrix, 1, max)
        up_bound <- median(rowm) + 100*mad(rowm)
        if(length(which(rowm > up_bound)) > 0){
          print("Outlier detected:")
          outliers <- fullmatrix[which.max(rowm),]

          print(paste("median of row max", median(rowm), "up_bound", up_bound))
          print(outliers)

          fullmatrix <- fullmatrix[-which.max(rowm),]

          if(heatmap){
            fullmatrixm <- log10(fullmatrix+1)
            print(pheatmap(fullmatrixm, cluster_rows = T, cluster_cols = F, fontsize_col=8, fontsize_row=8, labels_col=colLabel, angle_col=0, main=paste(querylabel, centerlabel, "after removing outliers")))
          }
        }
      }
      scoreMatrix_list[[querylabel]][[centerlabel]] <- fullmatrix
    }
  }


  color_store <- c("#00AFBB", "#E7B800", "#A0BDE0", "#0020C2", "#64E986", "#F5DEB3", "#C19A6B", "#E8A317", "#8E7618", "#A0522D", "#990012", "#CB6D51")

  plot_df <- list()
  stat_df <- list()
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


      xmin <- which(colLabel == hl[1])
      xmax <- which(colLabel == hl[2])
      submatrix <- (fullmatrix[, xmin:xmax])
      submatrix[is.na(submatrix)] <- 0
      Intensity <- as.numeric(rowMeans(submatrix)) + 0.001 ## add 0.001 to allow log transformation of y-axis

      Query <- as.factor(rep(querylabel, length(Intensity)))
      Reference <- as.factor(rep(centerlabel, length(Intensity)))
      subdf <- data.frame(Intensity, Query, Reference)
      stat_df[[paste(querylabel,centerlabel,sep=":")]] <- subdf

    }
  }

  mplot_dt <- rbindlist(plot_df) %>%
      mutate(Group=paste0(Query, ":", Reference), .keep="all") %>%
      mutate(lower=Intensity-se, upper=Intensity+se, .keep="all")
  mstat_dt <- rbindlist(stat_df) %>%
    mutate(Group=as.factor(paste0(Query, ":", Reference)), .keep="all")


  if(!is.null(outPrefix)){
    while(!is.null(dev.list())){
      dev.off()
    }
    pdf(paste(outPrefix, "pdf", sep="."), height=8, width=12)
    write.table(mstat_dt, paste(outPrefix, "_multiRef.tsv", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
    write.table(mplot_dt, paste(outPrefix, "_multiRef_bin.tsv", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
  }

  for(i in seq_along(querylabels)){
    for(beds in combn(querylabels, i, simplify = F)){
      for(j in seq_along(centerlabels)){
        for(centers in combn(centerlabels, j, simplify = F)){
          print(beds)
          print(centers)

          aplot_df <- mplot_dt %>%
            filter(Query %in% beds & Reference %in% centers)

          p <- ggplot(aplot_df, aes(x=Position, y=Intensity, color=Group)) + scale_fill_manual(values=color_store[1:(i*j)]) +
            geom_line(size=2) +
            geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
            geom_ribbon(aes(ymin=lower, ymax=upper, fill=Group), linetype=0, alpha=0.3) +
            annotate("rect", xmin=hl[1], xmax=hl[2], ymin=-Inf, ymax=Inf, fill="grey", color="grey", alpha=0.3) +
            theme_classic() + theme(legend.position="top") + xlab(Xlab) + ylab(Ylab) +
            theme(axis.title.x = element_text(face="bold", size=10), axis.title.y = element_text(face="bold", size=10))

          if((i == 1 && j > 1) || (i > 1 && j == 1)){

            astat_df <- mstat_dt %>%
              filter(Query %in% beds & Reference %in% centers)

            if(j > 1){
              comp <- combn(seq_along(centers),2, simplify=F)

              ps1 <- ggplot(astat_df, aes(x=Reference, y=Intensity, fill=Reference)) + scale_fill_manual(values=color_store[1:j]) +
                geom_boxplot(notch=FALSE) +
                theme_classic() +
                theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                theme(legend.position = "none") +
                labs(y=Ylab) +
                theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
                geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=TRUE, step_increase = 0.1) +
                scale_y_continuous(trans='log10')
            }else{
              comp <- combn(seq_along(beds),2, simplify=F)

              ps1 <- ggplot(astat_df, aes(x=Query, y=Intensity, fill=Query)) + scale_fill_manual(values=color_store[1:i]) +
                geom_boxplot(notch=FALSE) +
                theme_classic() +
                theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                theme(legend.position = "none") +
                labs(y=Ylab) +
                theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
                geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=TRUE, step_increase = 0.1) +
                scale_y_continuous(trans='log10')
            }


            print(plot_grid(p, ps1, ncol = 2, rel_widths = c(3,2)))
          }else if((i == 1 && j == 1) || (i == length(querylabels) && j == length(centerlabels))){
            print(p)
          }

        }
      }
    }
  }

  if(!is.null(inputfiles)){
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
          rowm <- apply(fullmatrix, 1, max)
          up_bound <- median(rowm) + 100*mad(rowm)  # mad: median absolute deviation

          if(length(which(rowm > up_bound)) > 0){
            print("Outlier detected:")
            outliers <- fullmatrix[which.max(rowm),]

            print(paste("median of row max", median(rowm), "up_bound", up_bound))
            print(outliers)

            fullmatrix <- fullmatrix[-which.max(rowm),]

            if(heatmap){
              fullmatrixm <- log10(fullmatrix+1)
              print(pheatmap(fullmatrixm, cluster_rows = T, cluster_cols = F, fontsize_col=8, fontsize_row=8, labels_col=colLabel, angle_col=0, main=paste(querylabel, centerlabel, "after removing outliers")))
            }
          }
        }

        ratioMatrix_list[[ratiolabels[i]]][[centerlabel]] <- fullmatrix
      }
    }

    plot_df <- list()
    stat_df <- list()
    print("collecting coverage data") ## plot multiple bed files on each center


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


        xmin <- which(colLabel == hl[1])
        xmax <- which(colLabel == hl[2])
        submatrix <- (fullmatrix[, xmin:xmax])
        submatrix[is.na(submatrix)] <- 0
        Intensity <- as.numeric(rowMeans(submatrix)) + 0.001 ## add 0.001 to allow log transformation of y-axis

        Query <- as.factor(rep(ratiolabel, length(Intensity)))
        Reference <- as.factor(rep(centerlabel, length(Intensity)))
        subdf <- data.frame(Intensity, Query, Reference)
        stat_df[[paste(ratiolabel,centerlabel,sep=":")]] <- subdf

      }
    }

    mplot_dt <- rbindlist(plot_df) %>%
      mutate(Group=paste0(Query, ":", Reference), .keep="all") %>%
      mutate(lower=Intensity-se, upper=Intensity+se, .keep="all")
    mstat_dt <- rbindlist(stat_df) %>%
      mutate(Group=as.factor(paste0(Query, ":", Reference)), .keep="all")


    if(!is.null(outPrefix)){
      write.table(mstat_dt, paste(outPrefix, "_multiRef_ratio.tsv", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
      write.table(mplot_dt, paste(outPrefix, "_multiRef_bin_ratio.tsv", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
    }

    for(i in seq_along(ratiolabels)){
      for(beds in combn(ratiolabels, i, simplify = F)){
        for(j in seq_along(centerlabels)){
          for(centers in combn(centerlabels, j, simplify = F)){
            print(beds)
            print(centers)

            aplot_df <- mplot_dt %>%
              filter(Query %in% beds & Reference %in% centers)

            p <- ggplot(aplot_df, aes(x=Position, y=Intensity, color=Group)) + scale_fill_manual(values=color_store[1:(i*j)]) +
              geom_line(size=2) +
              geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
              geom_ribbon(aes(ymin=lower, ymax=upper, fill=Group), linetype=0, alpha=0.3) +
              annotate("rect", xmin=hl[1], xmax=hl[2], ymin=-Inf, ymax=Inf, fill="grey", color="grey", alpha=0.3) +
              theme_classic() + theme(legend.position="top") + xlab(Xlab) + ylab(Ylab) +
              theme(axis.title.x = element_text(face="bold", size=10), axis.title.y = element_text(face="bold", size=10))

            if((i == 1 && j > 1) || (i > 1 && j == 1)){

              astat_df <- mstat_dt %>%
                filter(Query %in% beds & Reference %in% centers)

              if(j > 1){
                comp <- combn(seq_along(centers),2, simplify=F)

                ps1 <- ggplot(astat_df, aes(x=Reference, y=Intensity, fill=Reference)) + scale_fill_manual(values=color_store[1:j]) +
                  geom_boxplot(notch=FALSE) +
                  theme_classic() +
                  theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                  theme(legend.position = "none") +
                  labs(y=Ylab) +
                  theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                  stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
                  geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=TRUE, step_increase = 0.1) +
                  scale_y_continuous(trans='log10')
              }else{
                comp <- combn(seq_along(beds),2, simplify=F)

                ps1 <- ggplot(astat_df, aes(x=Query, y=Intensity, fill=Query)) + scale_fill_manual(values=color_store[1:i]) +
                  geom_boxplot(notch=FALSE) +
                  theme_classic() +
                  theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
                  theme(legend.position = "none") +
                  labs(y=Ylab) +
                  theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
                  stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
                  geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=TRUE, step_increase = 0.1) +
                  scale_y_continuous(trans='log10')
              }


              print(plot_grid(p, ps1, ncol = 2, rel_widths = c(3,2)))
            }else if((i == 1 && j == 1) || (i == length(ratiolabels) && j == length(centerlabels))){
              print(p)
            }

          }
        }
      }
    }

  }

  if(!is.null(outPrefix)) dev.off()

  return(list("plot"=mplot_dt, "stat"=mstat_dt))
}



plot_reference_locus_with_random <- function(txdb=NULL,
                                              queryfiles=NULL,
                                              centerfiles=NULL,
                                              ext=c(0,0),
                                              hl=c(0,0),
                                              querylabels=NULL,
                                              centerlabels=NULL,
                                              smo=FALSE,
                                              CLIP_reads=FALSE,
                                              extend_reads=0,
                                              norm=FALSE,
                                              binsize=10,
                                              refPoint="center",
                                              Xlab="Center",
                                              inputfiles=NULL,
                                              inputlabels=NULL,
                                              stranded=TRUE,
                                              heatmap=FALSE,
                                              scale=FALSE,
                                              outPrefix=NULL,
                                              genome="hg19",
                                              rm.outlier=F,
                                              n_random=1){

  if(0){
    ext=c(-100,100)
    hl=c(-50,50)
    smo=FALSE
    CLIP_reads=FALSE
    extend_reads=0
    norm=FALSE
    binsize=10
    refPoint="center"
    Xlab="Center"
    stranded=TRUE
    heatmap=FALSE
    scale=FALSE
    outPrefix=NULL
    genome="hg19"
    rm.outlier=F
  }


  names(querylabels) <- queryfiles
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

  bedInputs <- handle_input(queryfiles, CLIP_reads=CLIP_reads, extend_reads=extend_reads,  useScore=TRUE, outRle=FALSE, norm=norm, genome=genome)
  centerInputs <- handle_input(centerfiles, CLIP_reads=FALSE, extend_reads=0, useScore=FALSE, outRle=FALSE, norm=FALSE, genome=genome)

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
          windowRegions <- filter_by_overlaps(windowRegionsALL, region)
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
        ## if outliers are detected, remove the top outliers only
        if(is.null(inputfiles) && rm.outlier){
          rowm <- apply(fullmatrix, 1, max)
          up_bound <- median(rowm) + 100*mad(rowm)
          if(length(which(rowm > up_bound)) > 0){
            print("Outlier detected:")
            outliers <- fullmatrix[which.max(rowm),]

            print(paste("median of row max", median(rowm), "up_bound", up_bound))
            print(outliers)

            fullmatrix <- fullmatrix[-which.max(rowm),]

            if(heatmap){
              fullmatrixm <- log10(fullmatrix+1)
              print(pheatmap(fullmatrixm, cluster_rows = T, cluster_cols = F, fontsize_col=8, fontsize_row=8, labels_col=colLabel, angle_col=0, main=paste(querylabel, centerlabel, "after removing outliers")))
            }
          }
        }
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
          }
        )

        windowRs <- lapply(random_results, function(x) x$windowR)

        #system.time(
        #random_matricies <- parallel_apply_scoreMatrixBin(windowRs=windowRs, queryRegions=queryRegions, bin_num=bin_num, bin_op="mean", weight_col=weight_col, stranded=stranded)
        #)

        print(system.time(
        random_matricies <- lapply(windowRs, function(x){
          parallel_scoreMatrixBin(queryRegions, x, bin_num, bin_op, weight_col, stranded)
        })))
        #fullmatrix <- do.call(rbind, random_matricies)
        ## unify the dimensions of random_matrices

        dims <- sapply(random_matricies, dim)
        minrows <- min(dims[1,])
        random_matricies <- lapply(random_matricies, function(x)x[1:minrows,])


        #fullmatrix <- Reduce(`+`, random_matricies) / length(random_matricies) ## since the mean is biased by outliers, median should be used instead
        a3d <- array(unlist(random_matricies), c(dim(random_matricies[[1]]), length(random_matricies)))
        fullmatrix <- apply(a3d, 1:2, mean)


        fullmatrix[is.na(fullmatrix)] <- 0

        if(norm && !is.null(libsize)){
          fullmatrix <- fullmatrix*1000000/libsize
        }
        if(scale){
          smc <- t(apply(fullmatrix, 1, scales::rescale)) ## rescale to 0:1 range
          fullmatrix <- as.data.frame(smc)
        }

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

          xmin <- which(colLabel == hl[1])
          xmax <- which(colLabel == hl[2])
          submatrix <- (fullmatrix[, xmin:xmax])
          submatrix[is.na(submatrix)] <- 0
          Intensity <- as.numeric(rowMeans(submatrix)) + 0.001 ## add 0.001 to allow log transformation of y-axis

          Query <- as.factor(rep(querylabel, length(Intensity)))
          Reference <- as.factor(rep(alabel, length(Intensity)))
          subdf <- data.frame(Intensity, Query, Reference)

          plot_df <- rbind(plot_df, sub_df)
          stat_df <- rbind(stat_df, subdf)
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
          geom_line(size=1) + geom_point(color="grey30", size=2) +
          geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
          geom_ribbon(aes(ymin=lower, ymax=upper, fill=Reference), linetype=0, alpha=0.3) +
          annotate("rect", xmin=hl[1], xmax=hl[2], ymin=-Inf, ymax=Inf, fill="grey", color="grey", alpha=0.3) +
          theme_classic() + theme(legend.position="bottom") + xlab(Xlab) + ylab(Ylab) +
          theme(axis.title.x = element_text(face="bold", size=10), axis.title.y = element_text(face="bold", size=10)) +
          ggtitle(paste("Feature:", regionName, "\nReference size:", refsize, "\nSample name:", querylabel))


        ps1 <- ggplot(stat_df, aes(x=Reference, y=Intensity, fill=Reference)) + scale_fill_manual(values=c("#00AFBB", "#E7B800")) +
          geom_boxplot(notch=TRUE) +
          theme_classic() +
          theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
          theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
          theme(legend.position = "none") +
          labs(y="Signal intensity") + #expression(paste(Log[10], "(Signal intensity)"))) +
          scale_x_discrete(limits=c("Random", centerlabel)) + # labels=c("TSN" = expression(paste(m^6,"Am")), featureName = expression(paste("5'-UTR ", m^6, "A"))))
          stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
          geom_signif(comparisons = list(c(1, 2)), test='wilcox.test', map_signif_level=TRUE) +
          scale_y_continuous(trans='log10')
        #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x))

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


        print(plot_grid(p, ps1, ps2, ps3, ncol = 2, rel_widths = c(1,1)))


        ## the rest of the code is not used
        #ps2 <- add_pval(ps1, pairs = list(c(1, 2)), test='wilcox.test', heights=max(stat_df$Intensity), pval_star=TRUE)

        if(0){
          stat_df <- mutate(stat_df, Ranks=rank(Intensity))
          levels(stat_df$Group)
          RANDOM <- sum(stat_df[stat_df$Group == "Random", "Ranks"])
          CENTER <- sum(stat_df[stat_df$Group != "Random", "Ranks"])
          w <- wilcox.test(Intensity ~ Group, data=stat_df)
          t <- t.test(Intensity ~ Group, data=stat_df)
          k <- kruskal.test(Intensity ~ Group, data=stat_df)
          a <- aov(Intensity ~ Group, data=stat_df)
          h <- TukeyHSD(a)

          px <- ggbarplot(stat_df, x="Group", y="Intensity", add="mean_se", fill="Group", palette = c("#00AFBB", "#E7B800"),
                          position = position_dodge()) +
            stat_compare_means(label.y=mean(stat_df$Intensity)*3)
          #print(psx)
          #ps <- ggboxplot(stat_df, x="Group", y="Intensity", xlab=FALSE, ylab="Signal intensity", notch=TRUE, outlier.shape = 16) +
          #   stat_compare_means(label = "p.signif", method="wilcox.test", comparisons=list(c("random", centerlabel)), hide.ns = FALSE)


          stat_df_SE <- ddply(stat_df, c("Group"), summarise,
                              N    = length(Intensity),
                              mean = mean(Intensity),
                              median = median(Intensity),
                              sd   = sd(Intensity),
                              se   = sd / sqrt(N))

          print(stat_df_SE)
          print(paste("Rank sum Center bed", CENTER, "Random", RANDOM))
          print(k)

          ps4 <- ggplot(stat_df, aes(x=Group, y=Intensity, fill=Group)) +
            geom_bar(position=position_dodge(), stat="summary", fun="mean") +
            theme_classic() +
            theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=14)) +
            theme(legend.position = "none") +
            labs(y="Signal intensity") + #expression(paste(Log[10], "(Signal intensity)"))) +
            scale_x_discrete(limits=c("Random", centerlabel)) + # labels=c("TSN" = expression(paste(m^6,"Am")), featureName = expression(paste("5'-UTR ", m^6, "A")))) +
            theme(axis.text.x = element_text(face="bold", size=14, color="black")) +
            stat_compare_means(comparisons = list(c(1, 2))) +
            scale_y_continuous(trans='log10')


          #print(plot_grid(p, ps3, ncol = 2, rel_widths = c(2,1)))
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
            rowm <- apply(fullmatrix, 1, max)
            up_bound <- median(rowm) + 100*mad(rowm)  # mad: median absolute deviation

            if(length(which(rowm > up_bound)) > 0){
              print("Outlier detected:")
              outliers <- fullmatrix[which.max(rowm),]

              print(paste("median of row max", median(rowm), "up_bound", up_bound))
              print(outliers)

              fullmatrix <- fullmatrix[-which.max(rowm),]

              if(heatmap){
                fullmatrixm <- log10(fullmatrix+1)
                print(pheatmap(fullmatrixm, cluster_rows = T, cluster_cols = F, fontsize_col=8, fontsize_row=8, labels_col=colLabel, angle_col=0, main=paste(querylabel, centerlabel, "after removing outliers")))
              }
            }
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
            rowm <- apply(fullmatrix, 1, max)
            up_bound <- median(rowm) + 100*mad(rowm)  # mad: median absolute deviation

            if(length(which(rowm > up_bound)) > 0){
              print("Outlier detected:")
              outliers <- fullmatrix[which.max(rowm),]

              print(paste("median of row max", median(rowm), "up_bound", up_bound))
              print(outliers)

              fullmatrix <- fullmatrix[-which.max(rowm),]

              if(heatmap){
                fullmatrixm <- log10(fullmatrix+1)
                print(pheatmap(fullmatrixm, cluster_rows = T, cluster_cols = F, fontsize_col=8, fontsize_row=8, labels_col=colLabel, angle_col=0, main=paste(querylabel, centerlabel, "after removing outliers")))
              }
            }
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

            xmin <- which(colLabel == hl[1])
            xmax <- which(colLabel == hl[2])
            submatrix <- (fullmatrix[, xmin:xmax])
            submatrix[is.na(submatrix)] <- 0
            Intensity <- as.numeric(rowMeans(submatrix)) + 0.001 ## add 0.001 to allow log transformation of y-axis

            Query <- as.factor(rep(querylabel, length(Intensity)))
            Reference <- as.factor(rep(alabel, length(Intensity)))
            subdf <- data.frame(Intensity, Query, Reference)

            plot_df <- rbind(plot_df, sub_df)
            stat_df <- rbind(stat_df, subdf)
          }


          p <- ggplot(plot_df, aes(x=Position, y=Intensity, color=Reference)) + scale_fill_manual(values=c("#00AFBB", "#E7B800")) +
            geom_line(size=1) + geom_point(color="grey30", size=2) +
            geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
            geom_ribbon(aes(ymin=lower, ymax=upper, fill=Reference), linetype=0, alpha=0.3) +
            annotate("rect", xmin=hl[1], xmax=hl[2], ymin=-Inf, ymax=Inf, fill="grey", color="grey", alpha=0.3) +
            theme_classic() + theme(legend.position="bottom") + xlab(Xlab) + ylab(Ylab) +
            theme(axis.title.x = element_text(face="bold", size=10), axis.title.y = element_text(face="bold", size=10)) +
            ggtitle(paste("Feature:", regionName, "\nReference size:", refsize, "\nSample name:", ratiolabel))


          ps1 <- ggplot(stat_df, aes(x=Reference, y=Intensity, fill=Reference)) + scale_fill_manual(values=c("#00AFBB", "#E7B800")) +
            geom_boxplot(notch=FALSE) +
            theme_classic() +
            theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
            theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
            theme(legend.position = "none") +
            labs(y=Ylab) + #expression(paste(Log[10], "(Signal intensity)"))) +
            scale_x_discrete(limits=c("Random", centerlabel)) + # labels=c("TSN" = expression(paste(m^6,"Am")), featureName = expression(paste("5'-UTR ", m^6, "A"))))
            stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
            geom_signif(comparisons = list(c(1, 2)), test='wilcox.test', map_signif_level=TRUE) +
            scale_y_continuous(trans='log10')
          #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x))
          print(plot_grid(p, ps1, ncol = 2, rel_widths = c(1,1)))
        }
      }
    }
  }

  if(!is.null(outPrefix)) dev.off()
}



get_grWithA_at_TSS <- function(txdb, longest=TRUE){

  if(longest){
    tl_longest_tx <- extract_longest_tx(txdb)
    tx_longest <- transcripts(txdb, filter=list(tx_name=tl_longest_tx$tx_name), use.name=T)
    tss_tx_longest <- resize(tx_longest, width=1, fix="start", use.names=TRUE, ignore.strand=FALSE)
    seq_tss_tx_longest <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19, tss_tx_longest)
    seq_df_longest <- as.data.frame(seq_tss_tx_longest)
    A_seq_longest <- rownames(seq_df_longest)[seq_df_longest$x == "A"]
    A_tss_tx_longest <- tss_tx_longest[tss_tx_longest$tx_name %in% A_seq_longest]
    return(A_tss_tx_longest)
  }else{
    tl <- transcriptLengths(txdb, with.utr5_len = T, with.cds_len = T, with.utr3_len = T)
    tl_protein_coding <- tl[tl$cds_len > 0, ]
    tx_protein_coding <- transcripts(txdb, filter=list(tx_name=tl_protein_coding$tx_name), use.name=T)
    tss_tx_protein_coding <- resize(tx_protein_coding, width=1, fix="start", use.names=TRUE, ignore.strand=FALSE)
    seq_tss_tx_protein_coding <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19, tss_tx_protein_coding)
    seq_df_protein_coding <- as.data.frame(seq_tss_tx_protein_coding)
    A_seq_protein_coding <- rownames(seq_df_protein_coding)[seq_df_protein_coding$x == "A"]
    A_tss_tx_protein_coding <- tss_tx_protein_coding[tss_tx_protein_coding$tx_name %in% A_seq_protein_coding]
  }
}

count_overlap_at_TSS <- function(queryfiles, querylabels, txdb, ext=0, longest=TRUE){

  names(querylabels) <- queryfiles
  fileType <- "bed"
  if(grepl("\\.bam", queryfiles[1])){fileType <- "bam"}

  print(paste("The input type is", fileType, "file"))

  A_tss <- get_grWithA_at_TSS(txdb, longest=longest)
  gr <- flank(A_tss, width = ext, both=T)

  count_df <- NULL

  for(queryfile in queryfiles){
    #queryfile <- queryfiles[1]
    querylabel <- querylabels[queryfile]
    libsize <- NULL
    if(fileType == "bed"){
      queryRegions <- read_bed(queryfile)
      libsize <- countLines(queryfile)
    }else{

      ## get the 5'-end -1 position of the reads, which is the crosslink sites for iCLIP reads
      bamGR <- readBamFileAsGRanges(queryfile, min.mapq = 10)
      l <- end(bamGR) - start(bamGR)
      summary(l)
      queryRegions =flank(bamGR, width=1, both=F, start=T, ignore.strand=FALSE)
      libsize <- length(queryRegions)

    }

    overlap_count <- countOverlaps(gr, queryRegions, type="any")

    count_df <- cbind(count_df, overlap_count)

  }

  colnames(count_df) <- querylabels

  return(count_df)
}

count_overlap_in_features <- function(queryfiles, querylabels, txdb, longest=T, CLIP_reads=T){

  names(querylabels) <- queryfiles
  fileType <- "bed"
  if(grepl("\\.bam", queryfiles[1])){fileType <- "bam"}

  print(paste("The input type is", fileType, "file"))

  input_list <- list()

  for(queryfile in queryfiles){
    #queryfile <- queryfiles[1]
    querylabel <- querylabels[queryfile]
    libsize <- NULL
    if(fileType == "bed"){
      bedGR <- read_bed(queryfile)
      queryRegions =flank(bedGR, width=1, both=F, start=T, ignore.strand=FALSE)
      queryRegions <- ifelse(CLIP_reads, queryRegions, bedGR)
      libsize <- countLines(queryfile)
    }else{
      ga <- readGAlignments(queryfile, use.names=T, param=ScanBamParam(mapqFilter=10))
      libsize <- length(ga)
      if(CLIP_reads){
        ## get the 5'-end -1 position of the reads, which is the crosslink sites for iCLIP reads
        queryRegions =flank(granges(ga), width=1, both=F, start=T, ignore.strand=FALSE)
      }else{
        queryRegions <- stack(grglist(ga))
      }

    }
    input_list[[querylabel]] <- queryRegions
  }

  feature_count_list <- list()

  for(featureName in c("utr5", "cds", "utr3", "exon")){
    print(featureName)
    #featureName <- "utr5"
    alist <- gtf_to_bed_longest_tx(txdb, featureName, longest=longest)
    original_grl <- alist$GRangesList
    if(featureName %in% c("utr5", "exon")){
      trimmed <- trimTranscripts(original_grl, start=1, end=0)
      feature_grl <- trimmed
    }else{
      feature_grl <- original_grl
    }

    featureCount_df <- NULL

    for(queryfile in queryfiles){
      #queryfile <- queryfiles[1]
      print(queryfile)
      querylabel <- querylabels[queryfile]
      queryRegions <- input_list[[querylabel]]

      overlap_count <- countOverlaps(feature_grl, queryRegions, type="any")

      featureCount_df <- cbind(featureCount_df, overlap_count)
    }

    colnames(featureCount_df) <- querylabels

    feature_count_list[[featureName]] <- featureCount_df
  }
  return(feature_count_list)
}

run_DESeq2 <- function(data_table, s2c, project, contrasts){

  ddsMatrix <- DESeqDataSetFromMatrix(round(data_table),
                                      colData = s2c,
                                      design = ~ groups)
  class(ddsMatrix)

  #sizeFactors(ddsMatrix) <- guide_stat_sum/guide_gomean
  dds <- DESeq(ddsMatrix, fitType = "local") #  , sfType = "iterate"
  #rnames <- resultsNames(dds)
  #rnames

  #normalized_counts <- counts(dds, normalized=TRUE)
  #lfcShrinked <- lfcShrink(dds=dds, coef=3, type="apeglm")
  #lfcShrinked <- lfcShrinked[order(lfcShrinked$padj), ]
  #plotMDS(normalized_counts)
  #pheatmap(cor(normalized_counts))
  print(paste(project, "size factors"))
  print(sizeFactors(dds))

  X <- try(rld <- vst(dds))
  if(class(X) != "try-error"){
    pdf(paste(project, "PCA_plot.pdf", sep="_"))
    print(plotPCA(rld, intgroup=c("groups")))
    dev.off()
  }

  res_list <- NULL

  for(i in 1:nrow(contrasts)){
    #i <- 1
    contrast <- as.vector(as.character(contrasts[i,]))
    print(paste("Processing ", paste(contrast, collapse=" ")))
    res <- results(dds, contrast=contrast, independentFiltering = TRUE, pAdjustMethod = "BH")
    treatSF <- paste(project, contrast[2], "vs", contrast[3], sep="_")

    summary(res)
    res <- na.omit(res)

    mcols(res)

    res <- res[order(-res$stat),]
    res_list[[treatSF]] <- res

    filename <- paste(treatSF, "DESeq2_results_table.tab", sep="_")
    write.table(res, filename, row.names=T, col.names=NA, quote=F, sep="\t")

    pdf(paste(treatSF, "MA_plot.pdf", sep="_"))
    #plotMA(res, ylim=c(-5,5))
    plot(log2(res$baseMean), res$log2FoldChange, ylim=c(-5,5))
    abline(h=0)
    dev.off()

  }
  return(res_list)
}


liftOverBed <- function(chainFile, queryfile){
  outFile <- gsub("\\.bed", "_liftOver\\.bed", queryfile)
  chain <- import.chain(chainFile)
  bedGr <- import.bed(queryfile)
  lifted <- stack(liftOver(bedGr, chain))
  mcols(lifted) <- mcols(lifted)[,2:3] ## remove the name column of numbers
  write_bed(lifted, outFile)
  return(outFile)
}

## v1: large vector to be sub-sampled
## v2: small vector to provide density distribution
## wpv: wilcox.test p-value to indicate similarity between sub-sample and v2
sub_sample <- function(v1, v2, pv){
  #v1 <- out_lens[["_unregulated"]]
  #v2 <- c(out_lens[["_up"]],out_lens[["_down"]])
  #pv <- 0.1


  n = 0
  vsub <- NULL
  ds <- density(v2)
  dv1 <- approx(x=ds$x, y=ds$y, xout=v1, rule=1)
  probs <- dv1$y
  probs[is.na(probs)] <- 0
  summary(probs)

  success <- FALSE
  wtest <- NULL

  ## find subset of vi that has the same distribution as v2
  while(!success){
    n <- n + 1
    vsub <- sample(v1, round(length(probs[probs>0])*0.99), replace=F, prob=probs)
    #vsub <- sample(v1, length(v2), replace=F, prob=probs)
    wtest <- wilcox.test(vsub, v2)
    ktest <- ks.test(vsub, v2)

    print(paste(n, wtest$p.value, ktest$p.value))
    if(wtest$p.value > pv){
    #if(lr < 1){
      success <- TRUE
      print("density distribution of v2")
      print(summary(v2))
      print("subsampling of v1 based on density distribution of v2")
      print(summary(vsub))
    }
    v1 <- vsub
    dv1 <- approx(x=ds$x, y=ds$y, xout=v1, rule=1)
    probs <- dv1$y
    probs[is.na(probs)] <- 0

  }

  ## produce the final subset of v1 that have the same length as v2
  vsub <- sample(v1, length(v2), replace=F, prob=probs)

  return(vsub)
}


handle_input <- function(inputFiles, CLIP_reads=FALSE, extend_reads=0, useScore=FALSE, outRle=TRUE, norm=FALSE, genome="hg19"){

  if(0){
    inputFiles <- queryfiles
    CLIP_reads=FALSE
    extend_reads=0
    useScore=FALSE
    outRle=TRUE
    norm=FALSE
    genome="tetrahymena"

  }

  outlist <- lapply(inputFiles, function(inputFile){

    print(paste("Reading file:", inputFile))

    if(grepl("\\.bed", inputFile)){
      fileType <- "bed"
      out <- handle_bed(inputFile, extend_reads, useScore, outRle, genome)
    }else if(grepl("\\.bam", inputFile)){
      fileType <- "bam"
      out <- handle_bam(inputFile, CLIP_reads, extend_reads, outRle, genome)
    }else if(grepl("\\.wig", inputFile)){
      fileType <- "wig"
      out <- handle_wig(inputFile, outRle, genome)
    }else if(grepl("\\.bw|bigwig|bigWig|BigWig|BW|BIGWIG", inputFile)){
      fileType <- "bw"
      out <- handle_bw(inputFile, outRle, genome)
    }else{
      stop(paste("The input file format is not supported:", inputFile))
    }

    out
  })

  names(outlist) <- inputFiles

  if(norm && length(inputFiles)>1 && genome %in% c("hg19", "hg38")){
    outlist <- effectiveSize(outlist, outRle, genome)
  }

  return(outlist)
}

## outlist is the ouput of handle_input, use DESeq2 estimate size factor, which is used to multiply with libsize
## this fuction works for human only so far
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

handle_bed <- function(inputFile, extend_reads=0, useScore=FALSE, outRle=TRUE, genome="hg19"){

  beddata <- read.delim2(inputFile, header=F)
  beddata <- beddata[, 1:min(6,ncol(beddata))]  ## ignore extra columns, which cause problem in import.bed()
  colnames(beddata) <- c("chr", "start", "end", "id", "score", "strand")[1:min(6,ncol(beddata))]
  queryRegions <- makeGRangesFromDataFrame(beddata, keep.extra.columns=TRUE, starts.in.df.are.0based=TRUE)
  queryRegions$score <- as.numeric(queryRegions$score)

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

  if(extend_reads > 0) queryRegions <- resize(queryRegions, width=extend_reads, fix="start", ignore.strand=FALSE)

  if(outRle) queryRegions <- coverage(queryRegions, weight= weight_col)  ## when bed is used is windows in ScoreMatrix, do not covert to Rle

  return(list("query"=queryRegions, "size"=libsize, "type"="bed", "weight"=weight_col))
}

handle_bam <- function(inputFile, CLIP_reads=FALSE, extend_reads=0, outRle=TRUE, genome="hg19"){

  if(0){
    inputFile <- "Greenblatt_001_Strip1_A01_E-H3_3_S1_S34_sorted.bam"
    CLIP_reads=FALSE
    extend_reads=0
    outRle=TRUE
    genome="tetrahymena"
  }

  ga <- readGAlignments(inputFile, use.names=T, param=ScanBamParam(mapqFilter=10))
  libsize <- sum(idxstatsBam(inputFile)$mapped)
  if(CLIP_reads){
    ## get the 5'-end -1 position of the reads, which is the crosslink sites for iCLIP reads
    queryRegions <- flank(granges(ga), width=1, both=F, start=T, ignore.strand=FALSE)
    score(queryRegions) <- 1
  }else if(extend_reads > 0){
    queryRegions <- resize(granges(ga), width=extend_reads, fix="start", ignore.strand=FALSE)
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

handle_wig <- function(inputFile, outRle=TRUE, norm=TRUE, genome="hg19"){

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

  out <- handle_bw(bwfile, outRle, norm, genome)
  out
}



## find the transcript that a utr3 belongs to
utr3_to_transcript <- function(utr3File, txdb){
  #utr3File <- centerfiles[1]

  transcriptFile <- gsub("\\.bed", "_TSS\\.bed", utr3File)

  longest_tx <- extract_longest_tx(txdb)
  utrInput <- handle_input(utr3File)
  utrGR <- utrInput$query

  w <- width(utrGR)

  feature <- transcripts(txdb, use.name=T)
  feature_longest <- feature[names(feature) %in% longest_tx$tx_name]

  feature_overlap <- subsetByOverlaps(feature_longest, utrGR, type="any", minoverlap=min(w))

  export.bed(feature_overlap, transcriptFile)

  return(transcriptFile)
}

extract_intron_VAST <- function(inputfile, length_filter=NULL){
  inputfile <- "C:/GREENBLATT/Nabeel/VAST-TOOLS/ZBTB7A/Intron_retention_regulated_events.tab"

  in_df <- read.delim(inputfile, header=TRUE)
  outfiles <- list()
  out_beds <- list()
  out_lens <- list()

  for(change in c("_up", "_down", "_unregulated")){
    #change <- "_unregulated"
    print(change)


    outfiles[[change]] <- gsub("\\.tab", paste0(change,"\\.bed"), inputfile)

    sub_df <- filter(in_df, grepl(change, GROUP))

    fullcoor <- lapply(as.list(sub_df$FullCO), function(x){unlist(strsplit(x, split=":", fixed=T))})
    fullcoor_strand <- unlist(lapply(fullcoor, function(x)x[3]))

    introncoor <- lapply(as.list(sub_df$COORD), function(x){unlist(strsplit(x, split=":|-", fixed=F))})
    introncoor_chr <- unlist(lapply(introncoor, function(x)x[1]))
    introncoor_start <- unlist(lapply(introncoor, function(x)x[2]))
    introncoor_end <- unlist(lapply(introncoor, function(x)x[3]))


    out_bed <- data.frame("chr"=introncoor_chr, "start"=as.integer(introncoor_start), "end"=as.integer(introncoor_end), "id"=sub_df$EVENT, "score"=sub_df$MV.dPsi._at_0.95, "strand"=fullcoor_strand)
    #write.table(out_bed, outputfile, row.names=F, col.names=F, sep="\t", quote=F)
    out_beds[[change]] <- out_bed
    out_lens[[change]] <- out_bed$end - out_bed$start + 1
  }

  ## test if up and down have different intron length distribution from the unregulated
  k <- ks.test(out_lens[["_unregulated"]], c(out_lens[["_up"]],out_lens[["_down"]]))

  print(k$p.value)
  print("Performing sub-sampling to match length distribution")
  sub <- sub_sample(out_lens[["_unregulated"]], c(out_lens[["_up"]],out_lens[["_down"]]), 0.9)
  k <- ks.test(sub, c(out_lens[["_up"]],out_lens[["_down"]]))
  print(k$p.value)
  print(summary(sub))

  plot(quantile(c(out_lens[["_unregulated"]]), probs = seq(0, 1, 0.01)), col="cyan")
  points(quantile(c(out_lens[["_up"]], out_lens[["_down"]]), probs = seq(0, 1, 0.01)), col="blue")
  points(quantile(sub, probs = seq(0, 1, 0.01)), col="red")

  qqplot(x=log(sub), y=log(c(out_lens[["_up"]], out_lens[["_down"]])))

  sub_index <- match(sub, out_lens[["_unregulated"]])

  out_beds[["_unregulated"]] <- data.frame(out_beds[["_unregulated"]])[sub_index,]

  ## filter final output using length_filter
  for(change in c("_up", "_down", "_unregulated")){
    if(is.integer(length_filter)){
      out_beds[[change]] <- data.frame(out_beds[[change]]) %>% filter((end-start+1) > length_filter)
    }
    write.table(out_beds[[change]], outfiles[[change]], row.names=F, col.names=F, sep="\t", quote=F)
  }

  unlist(outfiles)
}

plot_bam_correlation <- function(bamfiles, bamlabels, binsize=1000000, outPrefix="Bam_correlation", genome="hg19"){
  #genome <- "hg19"
  print("Computing bam correlation")
  outlist <- handle_input(inputFiles=bamfiles, CLIP_reads=FALSE, extend_reads=0, useScore=FALSE, outRle=TRUE, norm=FALSE, genome=genome)
  ## find common chromosomes
  chr_list <- lapply(outlist, function(x) names(x$query))
  chr_list[["regular_chr"]] <- paste0("chr", c(seq(1,22), c("X", "Y", "M")))
  comchr <- sort(Reduce(intersect, chr_list))

  seqi <- Seqinfo(genome=genome)
  tilewidth <- binsize

  cl <- start_parallel(length(outlist))
  parallel::clusterExport(cl, c("tileGenome", "binnedAverage"))
  parallel::clusterExport(cl, c("seqi", "tilewidth", "comchr"), envir=environment())
  score_list <- parLapply(cl, outlist, function(x){

    tileBins <- tileGenome(seqi[comchr], tilewidth=tilewidth, cut.last.tile.in.chrom=TRUE)
    binAverage <- binnedAverage(tileBins, x$query[comchr], varname="binned_score", na.rm=F)
    binAverage$binned_score
  })
  stop_parallel(cl)

  mat <- matrix(unlist(score_list), ncol=length(score_list), byrow=F)
  mat[is.na(mat)] <- 0
  df <- as.data.frame(round(mat * tilewidth) + 1)

  colnames(df) <- bamlabels

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

  pdf(paste0(op, ".pdf"), width=8, height=8)
  if(length(bamfiles)>3) pheatmap(cor(log(df+1)), display_numbers = T)
  pairs(log(df+1), lower.panel = panel.smooth, upper.panel=panel.cor, diag.panel = panel.hist, main=paste("log (reads/bin)\nbin size =", binsize))
  dev.off()

}

## 'simple' is for annotations withou utrs
## 'RNA' focus on transcripts, not genomic DNA
annotate_peaks <- function(peakfile, peaklabel, gtffile, genome="hg19", fiveP=1000, threeP=1000, simple=FALSE, RNA=TRUE){
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

  peak <- handle_bed(inputFile=peakfile, extend_reads=0, useScore=FALSE, outRle=FALSE, genome=genome)$query

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

    write.table(targeted_gene, paste(peaklabel, "_targeted_gene.tab", sep=""), sep="\t", row.names=F, quote=F)
    write.table(targeted_promoter, paste(peaklabel, "_targeted_promoter.tab", sep=""), sep="\t", row.names=F, quote=F)
    write.table(summary_table, paste(peaklabel, "_target_summary.tab", sep=""), sep="\t", row.names=F, quote=F)

    features <- GRangesList("Promoter"=promoter,
                            "TTS"=TTS,
                            "Exon"=exon$GRanges,
                            "Intron"=intron$GRanges,
                            compress=F)

    pdf(paste(peaklabel, "_feature_type_piechart.pdf"))

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
    txdbFeatures <- getTxdbFeaturesFromGRanges(gff)

    dt <- getTargetedGenesTable(queryRegions = peak, txdbFeatures = txdbFeatures)
    dt <- dt[order(transcripts, decreasing = TRUE)]

    dt_gene <- unique(merge(dt, gene_info_table, by.x="tx_name", by.y="transcript_id", all.x=T))
    dt_gene <- dt_gene[order(transcripts, decreasing = TRUE)]
    dt_gene <- dt_gene[!is.na(gene_biotype)]
    dim(dt_gene)
    head(dt_gene)
    tail(dt_gene)
    write.table(dt_gene, paste(peaklabel, "_targeted_gene.tab", sep=""), sep="\t", row.names=F, quote=F)

    # To find out the distribution of the query regions across gene types:
    biotype_col <- grep('gene_biotype', colnames(overlaps), value = T)
    df <- overlaps[,length(unique(queryIndex)), by = biotype_col] %>%
      rename_with(~ c("gene_type", "count"))
    intergenic <- data.frame("gene_type"= "intergenic", "count"= (length(peak) - sum(df$count)))
    selected <- c("protein_coding", "lincRNA", "antisense", "pseudogene", "snRNA", "snoRNA", "rRNA")
    selected_df <- filter(df, gene_type %in% selected)
    other_df <- filter(df, !gene_type %in% selected)
    other <- data.frame("gene_type"="other", "count"=sum(other_df$count))

    df <- rbind(selected_df, other, intergenic) %>%
      mutate(percent = round(count*100 / sum(count), 1)) %>%
      arrange(desc(count))

    pdf(paste(peaklabel, "_gene_type_distribution_piechart.pdf"), height=8, width=10)
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
                            "5'UTR"=utr5$GRanges,
                            "3'UTR"=utr3$GRanges,
                            "CDS"=cds$GRanges,
                            "Intron"=intron$GRanges, compress=F)
    if(RNA){
      features <- GRangesList("5'UTR"=utr5$GRanges,
                              "3'UTR"=utr3$GRanges,
                              "CDS"=cds$GRanges,
                              "Intron"=intron$GRanges, compress=F)
    }

    annot = annotateWithFeatures(peak, features, strand.aware=TRUE, intersect.chr=FALSE)
    precedence_count <- annot@num.precedence
    precedence_count["NoFeature"] <- length(peak) - sum(precedence_count)
    #plotTargetAnnotation(annot)
    #summary <- summarizeQueryRegions(queryRegions = peak,
    #                                txdbFeatures = features)

    df <- data.frame(precedence_count) %>%
      mutate(feature = rownames(.)) %>%
      mutate(percent = precedence_count / sum(precedence_count)) %>%
      mutate(labels = scales::percent(percent, accuracy = 0.1))

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
}

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
  ## if outliers are detected, remove the top outliers only
  if(rm.outlier){
    fullmatrix <- rmOutlier(fullmatrix)
    if(heatmap){
      fullmatrixm <- log10(fullmatrix+1)
      print(pheatmap(fullmatrixm, cluster_rows = T, cluster_cols = F, fontsize_col=8, fontsize_row=8, labels_col=colLabel, angle_col=0, main=paste(querylabel, centerlabel, "after removing outliers")))
    }

  }
  return(fullmatrix)
}

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
