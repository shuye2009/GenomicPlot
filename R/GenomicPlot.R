


#' @title  Plot signals around the start and the end of genomic features
#'
#' @description   Plot reads or peak signal intensity of samples in the query files around stat, end and center of genomic 
#' features. The upstream and downstream windows can be given separately. If Input files are provided, ratio over Input is
#' computed and displayed as well. A random feature can be generated to serve as a background.
#'
#' @param queryFiles a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param inputFiles a vector of input sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param txdb a TxDb object defined in the GenomicFeatures package. If featureName is "custom", a bed file that define the feature
#' @param featureName one of the gene feature in c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene") or "custom"
#' @param handleInputParams a list of parameters for \code{handle_input}
#' @param binSize an integer defines bin size for intensity calculation
#' @param ext a vector of four integers defining upstream and downstream boundaries of the plot window, flanking the start and end of features
#' @param hl a vector of four integers defining upstream and downstream boundaries of the highlight window, flanking the start and end of features
#' @param insert an integer specifies the length of the center regions to be included, in addition to the start and end of the feature
#' @param randomize logical, indicating if randomized feature should generated and used as a contrast to the real feature
#' @param stranded logical, indicating whether the strand of the feature should be considered
#' @param scale logical, indicating whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param shade logical indicating whether to place a shaded rectangle around the point of interest
#' @param smooth logical, indicating whether the line should smoothed with a spline smoothing algorithm
#' @param rmOutlier logical, indicating whether a row with abnormally high values in the score matrix should be removed
#' @param longest logical, indicating whether the output should be limited to the longest transcript of each gene
#' @param transform logical, whether to log2 transform the data matrix
#' @param outPrefix a string specifying output file prefix for plots (outPrefix.pdf)
#' @param verbose logical, whether to output additional information (data used for plotting or statistical test results)
#' @param nc integer, number of cores for parallel processing
#'
#' @return a list of two objects, the first is a GRanges object, the second is a GRangesList object
#' @author Shuye Pu
#'
#' @examples
#' txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
#' queryFiles <- system.file("data", "test_clip_chr19.bam", package="GenomicPlotData")
#' names(queryFiles) <- "query"
#' inputFiles <- system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData")
#' names(inputFiles) <- "input"
#' ext <- c(-500, 200, -200, 500)
#' hl <- c(-50, 50, -50, 50)
#' op <- NULL
#' handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="start", norm=TRUE, useScore=FALSE,
#'  outRle=TRUE, useSizeFactor=TRUE, genome="hg19")
#' plot_start_end_feature(queryFiles=c(queryFiles,inputFiles), inputFiles=NULL, txdb=txdb, featureName="intron",
#'  binSize=10, handleInputParams=handleInputParams, longest=TRUE, ext=ext, hl=hl, randomize=FALSE,
#'  insert=100, stranded=TRUE, scale=FALSE, smooth=TRUE, outPrefix=op, nc=2)
#'
#' @export plot_start_end_feature
#'
#'

plot_start_end_feature <- function(queryFiles, inputFiles=NULL, txdb, featureName, handleInputParams=NULL, binSize=10,
   insert=0, verbose=FALSE, longest=TRUE, ext=c(-500, 200, -200, 500), hl=c(-50, 50, -50, 50), randomize=FALSE, 
   stranded=TRUE, scale=FALSE, smooth=FALSE, rmOutlier=FALSE, outPrefix="plots", transform=FALSE, shade=TRUE, nc=2,
   Ylab="Signal intensity"){

  queryLabels <- names(queryFiles)
  names(queryLabels) <- queryFiles
  inputLabels <- names(inputFiles)
  names(inputLabels) <- inputFiles

  if(!is.null(inputFiles)){
    if(length(queryFiles) != length(inputFiles)){
      stop("the number of queryFiles and inputFiles must be the same!")
    }
    names(inputLabels) <- inputFiles
    queryFiles <- c(queryFiles, inputFiles)
    queryLabels <- c(queryLabels, inputLabels)
  }

  if(!is.null(outPrefix)){
     while(!is.null(dev.list())){
        dev.off()
     }
     pdf(paste0(outPrefix, ".pdf"), width=10, height=8)
  }

  feature <- rfeature <- NULL
  fs <- fe <- rfs <- rfe <- fc <- rfc <- NULL

  if(featureName %in%  c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene")){
     feature <- get_genomic_feature_coordinates(txdb, featureName, longest=longest, protein_coding=TRUE)[["GRanges"]]
     minimal_width <- ext[2] - ext[3] + insert
     feature <- feature[width(feature)>minimal_width]
  }else if(featureName == "custom"){
     bedparam <- handleInputParams
     # the input bed is going to be used as center, so the handleInputParams has to be modified accordingly
     bedparam$fix_width <- 0
     bedparam$norm <- FALSE
     bedparam$useScore <- FALSE
     bedparam$outRle <- FALSE
     feature <- handle_input(inputFiles=txdb, bedparam, nc=nc)[[1]]$query
     featureName <- names(txdb)
  }else{
     stop(paste(featureName, "is not supported!"))
  }

  if(verbose) print(paste("number of features: ", featureName, length(feature)))

  fs <- trim(promoters(resize(feature,width=1,fix="start"), upstream=-ext[1],downstream=ext[2]))
  fe <- trim(promoters(resize(feature,width=1,fix="end"), upstream=-ext[3],downstream=ext[4]))
  fc <- trim(promoters(resize(feature,width=1,fix="center"), upstream=round(insert/2),downstream=round(insert/2)))

   if(randomize){
    random_points <- sample(ext[1]:ext[4], length(feature), replace=TRUE)
    rfeature <- shift(feature, shift=random_points, use.names=TRUE)
    rfs <- trim(promoters(resize(rfeature,width=1,fix="start"), upstream=-ext[1],downstream=ext[2]))
    rfe <- trim(promoters(resize(rfeature,width=1,fix="end"), upstream=-ext[3],downstream=ext[4]))
    rfc <- trim(promoters(resize(rfeature,width=1,fix="center"), upstream=round(insert/2),downstream=round(insert/2)))
  }

  ext[2] <- ext[2] - (ext[2]-ext[1])%%binSize ## to avoid binSize inconsistency, as the final binSize is dictated by bin_num
  bin_num_s <- round((ext[2]-ext[1])/binSize)
  ext[4] <- ext[4] - (ext[4]-ext[3])%%binSize ## to avoid binSize inconsistency, as the final binSize is dictated by bin_num
  bin_num_e <- round((ext[4]-ext[3])/binSize)
  bin_num_c <- round(insert/binSize)

  mat_list <- NULL
  mat_list[[paste("Start of", featureName)]] <- list("window"=fs, "rwindow"=rfs, s=ext[1], e=ext[2], "xmin"=hl[1], "xmax"=hl[2], "bin_num"=bin_num_s)
  mat_list[[paste("Center of", featureName)]] <- list("window"=fc, "rwindow"=rfc,  s=-round(insert/2), e=round(insert/2), "xmin"=0, "xmax"=0, "bin_num"=bin_num_c)
  mat_list[[paste("End of", featureName)]] <- list("window"=fe, "rwindow"=rfe,  s=ext[3], e=ext[4], "xmin"=hl[3], "xmax"=hl[4], "bin_num"=bin_num_e)

  queryInputs <- handle_input(queryFiles, handleInputParams, nc=nc)

  scoreMatrix_list <- list()
  scoreMatrix_list_random <- list()
  for(locus in names(mat_list)){
    windowR <- mat_list[[locus]]$window
    rwindowR <- mat_list[[locus]]$rwindow

    bin_num <- mat_list[[locus]]$bin_num
    if(bin_num <= 0) next

    for(queryFile in queryFiles){

      queryLabel <- queryLabels[queryFile]
      if(verbose) print(queryLabel)
      queryRegions <- queryInputs[[queryFile]]$query
      libsize <- queryInputs[[queryFile]]$size

      bin_op <- "mean"
      weight_col <- queryInputs[[queryFile]]$weight

     # fullMatrix1 <- parallel_scoreMatrixBin(queryRegions, windowR, bin_num, bin_op, weight_col, stranded)
      fullMatrix <- parallel_scoreMatrixBin(queryRegions, windowR, bin_num, bin_op, weight_col, stranded, nc=nc)
      scoreMatrix_list[[queryLabel]][[locus]] <- fullMatrix

      if(randomize){
       # rfullMatrix1 <- parallel_scoreMatrixBin(queryRegions, rwindowR, bin_num, bin_op, weight_col, stranded)
        rfullMatrix <- parallel_scoreMatrixBin(queryRegions, rwindowR, bin_num, bin_op, weight_col, stranded, nc=nc)
        scoreMatrix_list_random[[queryLabel]][[locus]] <- rfullMatrix
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

    if(bin_num <= 0) next
    for(queryFile in queryFiles){

      queryLabel <- queryLabels[queryFile]
      if(verbose) print(queryLabel)

      fullMatrix <- scoreMatrix_list[[queryLabel]][[locus]]
      fullMatrix <- process_scoreMatrix(fullMatrix, scale, rmOutlier, transform=transform, verbose=verbose)

      colm <- apply(fullMatrix, 2, mean)
      colsd <- apply(fullMatrix, 2, sd)
      colse <- colsd/sqrt(nrow(fullMatrix))
      collabel <- seq(start, (end-binSize), binSize)
      querybed <- as.factor(rep(queryLabel, ncol(fullMatrix)))
      location <- as.factor(rep(locus, ncol(fullMatrix)))
      levels(location) <- rev(levels(location))
      halfmin <- min(fullMatrix)
      intervals <- apply(fullMatrix, 2, function(x) sum(x>halfmin))

      sub_df <- NULL
      sub_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Interval"=intervals, "Position"=collabel, "Query"=querybed, "Location"=location)
      if(smooth){
        sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df=as.integer(bin_num/5))$y)
        sub_df$se <- as.vector(smooth.spline(sub_df$se, df=as.integer(bin_num/5))$y)
      }
      sub_df <- mutate(sub_df, lower=Intensity-se, upper=Intensity+se)
      plot_df <- rbind(plot_df, sub_df)

      if(randomize){

        rfullMatrix <- scoreMatrix_list_random[[queryLabel]][[locus]]
        rfullMatrix <- process_scoreMatrix(rfullMatrix, scale, rmOutlier, transform=transform, verbose=verbose)

        rcolm <- apply(rfullMatrix, 2, mean)
        rcolsd <- apply(rfullMatrix, 2, sd)
        rcolse <- rcolsd/sqrt(nrow(rfullMatrix))
        rcollabel <- seq(start, (end-binSize), binSize)
        rquerybed <- as.factor(rep(paste0(queryLabel,":Random"), ncol(rfullMatrix)))
        location <- as.factor(rep(locus, ncol(fullMatrix)))
        levels(location) <- rev(levels(location))
        halfmin <- min(fullMatrix)
        intervals <- apply(fullMatrix, 2, function(x) sum(x>halfmin))

        rsub_df <- NULL
        rsub_df <- data.frame("Intensity"=rcolm, "sd"=rcolsd, "se"=rcolse, "Interval"=intervals, "Position"=rcollabel, "Query"=rquerybed, "Location"=location)
        if(smooth){
          rsub_df$Intensity <- as.vector(smooth.spline(rsub_df$Intensity, df=as.integer(bin_num/5))$y)
          rsub_df$se <- as.vector(smooth.spline(rsub_df$se, df=as.integer(bin_num/5))$y)
        }
        rsub_df <- mutate(rsub_df, lower=Intensity-se, upper=Intensity+se)

        plot_df <- rbind(plot_df, rsub_df)

      }
    }
  }

  ## plot individual bed line for one feature
 Ylab <- ifelse(transform, paste0("Log2 (", Ylab, ")"), Ylab)
  for (query in unique(plot_df$Query)){
     aplot_df <- plot_df %>%
        filter(Query == query)
     #p <- draw_locus_profile(aplot_df, xc="Position", yc="Intensity", cn="Query", sn="Reference", Xlab="", Ylab=Ylab, shade=shade, hl=c(xmin,xmax)) +
      #  ggtitle(featureName) +
      #  facet_wrap(~Location, scales="free_x")
     #print(p)
     plots <- draw_stacked_profile(plot_df=aplot_df, cn="Query", ext=ext, hl=hl, atitle=featureName, insert=insert, Ylab=Ylab, shade=shade, stack=TRUE)
     
     print(plots)
  }
  ## plot multi bed lines for one feature
  if(length(unique(plot_df$Query))>1){
     #p <- draw_locus_profile(plot_df, xc="Position", yc="Intensity", cn="Query", sn="Reference", Xlab="", Ylab=Ylab, shade=shade, hl=c(xmin,xmax)) +
      #    ggtitle(featureName) +
      # facet_wrap(~Location, scales="free_x")
     #print(p)
     
     plots <- draw_stacked_profile(plot_df=plot_df, cn="Query", ext=ext, hl=hl, atitle=featureName, insert=insert, Ylab=Ylab, shade=shade, stack=TRUE)
     
     print(plots)
  }

  #plots <- draw_stacked_profile(plot_df=plot_df, cn="Query", ext=ext, hl=hl, atitle=featureName, insert=insert, Ylab=Ylab, shade=shade, stack=TRUE)

  #print(plots)

  if(!is.null(inputFiles)){
    Ylab <- ifelse(transform, "Log2 (Ratio-over-Input)", "Ratio-over-Input")

    ratiofiles <- queryFiles[!queryFiles %in% inputFiles]
    ratiolabels <- queryLabels[!queryLabels %in% inputLabels]
    inputMatrix_list <- scoreMatrix_list[inputLabels]
    ratioMatrix_list <- scoreMatrix_list[ratiolabels]

    inputMatrix_list_random <- scoreMatrix_list_random[inputLabels]
    ratioMatrix_list_random <- scoreMatrix_list_random[ratiolabels]

    for(locus in names(mat_list)){

       bin_num <- mat_list[[locus]]$bin_num
       if(bin_num <= 0) next

      for(i in seq_along(ratiolabels)){

        rm <- ratioMatrix_list[[ratiolabels[i]]][[locus]]
        im <- inputMatrix_list[[inputLabels[i]]][[locus]]
        minrow <- min(nrow(rm), nrow(im))

        fullMatrix <- rm[1:minrow,]/im[1:minrow,]
        ratioMatrix_list[[ratiolabels[i]]][[locus]] <- fullMatrix

        ## for random feature
        if(randomize){
          rmr <- ratioMatrix_list_random[[ratiolabels[i]]][[locus]]
          imr <- inputMatrix_list_random[[inputLabels[i]]][[locus]]
          minrowr <- min(nrow(rmr), nrow(imr))

          fullMatrix <- rmr[1:minrowr,]/imr[1:minrowr,]
          ratioMatrix_list_random[[ratiolabels[i]]][[locus]] <- fullMatrix
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

      if(bin_num <= 0) next
      for(ratiofile in ratiofiles){

        ratiolabel <- ratiolabels[ratiofile]
        if(verbose) print(ratiolabel)

        fullMatrix <- ratioMatrix_list[[ratiolabel]][[locus]]
        fullMatrix <- process_scoreMatrix(fullMatrix, scale=FALSE, rmOutlier, transform=transform, verbose=verbose)

        colm <- apply(fullMatrix, 2, mean)
        colsd <- apply(fullMatrix, 2, sd)
        colse <- colsd/sqrt(nrow(fullMatrix))
        collabel <- seq(start, (end-binSize), binSize)
        ratiobed <- as.factor(rep(ratiolabel, ncol(fullMatrix)))
        location <- as.factor(rep(locus, ncol(fullMatrix)))
        levels(location) <- rev(levels(location))
        halfmin <- min(fullMatrix)
        intervals <- apply(fullMatrix, 2, function(x) length(x[x>halfmin]))

        sub_df <- NULL
        sub_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Interval"=intervals, "Position"=collabel, "Query"=ratiobed, "Location"=location)
        if(smooth){
          sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df=as.integer(bin_num/5))$y)
          sub_df$se <- as.vector(smooth.spline(sub_df$se, df=as.integer(bin_num/5))$y)
        }
        sub_df <- mutate(sub_df, lower=Intensity-se, upper=Intensity+se)
        plot_df <- rbind(plot_df, sub_df)

        if(randomize){

          rfullMatrix <- ratioMatrix_list_random[[ratiolabel]][[locus]]
          rfullMatrix <- process_scoreMatrix(rfullMatrix, scale=FALSE, rmOutlier, transform=transform, verbose=verbose)

          rcolm <- apply(rfullMatrix, 2, mean)
          rcolsd <- apply(rfullMatrix, 2, sd)
          rcolse <- rcolsd/sqrt(nrow(rfullMatrix))
          rcollabel <- seq(start, (end-binSize), binSize)
          rratiobed <- as.factor(rep(paste0(ratiolabel,":Random"), ncol(rfullMatrix)))
          location <- as.factor(rep(locus, ncol(fullMatrix)))
          levels(location) <- rev(levels(location))
          halfmin <- min(rfullMatrix)
          intervals <- apply(rfullMatrix, 2, function(x) length(x[x>halfmin]))

          rsub_df <- NULL
          rsub_df <- data.frame("Intensity"=rcolm, "sd"=rcolsd, "se"=rcolse, "Interval"=intervals, "Position"=rcollabel, "Query"=rratiobed, "Location"=location)
          if(smooth){
            rsub_df$Intensity <- as.vector(smooth.spline(rsub_df$Intensity, df=as.integer(bin_num/5))$y)
            rsub_df$se <- as.vector(smooth.spline(rsub_df$se, df=as.integer(bin_num/5))$y)
          }
          rsub_df <- mutate(rsub_df, lower=Intensity-se, upper=Intensity+se)

          plot_df <- rbind(plot_df, rsub_df)
        }
      }
    }

    ## plot individual bed line for one Location
    for (query in unique(plot_df$Query)){
       aplot_df <- plot_df %>%
          filter(Query == query)
       p <- draw_locus_profile(aplot_df, xc="Position", yc="Intensity", cn="Query", sn="Reference", Xlab="", Ylab=Ylab, shade=shade, hl=c(xmin,xmax)) +
          ggtitle(featureName) +
          facet_wrap(~Location, scales="free_x")
       print(p)
    }

    ## plot multi bed lines for one Location
    if(length(unique(plot_df$Query))>1){
       p <- draw_locus_profile(plot_df, xc="Position", yc="Intensity", cn="Query", sn="Reference", Xlab="", Ylab=Ylab, shade=shade, hl=c(xmin,xmax)) +
         ggtitle(featureName) +
         facet_wrap(~Location, scales="free_x")
       print(p)
    }

    plots <- draw_stacked_profile(plot_df=plot_df, cn="Query", ext=ext, hl=hl, atitle=featureName, insert=insert, Ylab=Ylab, shade=shade)

    print(plots)
  }
  if(!is.null(outPrefix)){
    on.exit(dev.off(), add=TRUE)
  }

  invisible(plot_df)
}


#' @title  Plot signals around the start and the end of custom features
#
#' @description   Plot reads or peak signal intensity of samples in the query files around start and end of custom features. The upstream and downstream windows
#' can be given separately, within the window, a smaller window can be defined to highlight region of interest.
#' A line plot will be displayed for both start and end of feature. If Input files are provided, ratio over Input is
#' computed and displayed as well.
#'
#' @param queryFiles a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param centerFiles  bed files that define the custom features
#' @param inputFiles a vector of input sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param handleInputParams a list of parameters for \code{handle_input}
#' @param binSize an integer defines bin size for intensity calculation
#' @param ext a vector of four integers defining upstream and downstream boundaries of the plot window, flanking the start and end of features
#' @param hl a vector of four integers defining upstream and downstream boundaries of the highlight window, flanking the start and end of features
#' @param insert an integer specifies the length of the center regions to be included, in addition to the start and end of the feature
#' @param stranded logical, indicating whether the strand of the feature should be considered
#' @param scale logical, indicating whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smooth logical, indicating whether the line should smoothed with a spline smoothing algorithm
#' @param rmOutlier logical, indicating whether a row with abnormally high values in the score matrix should be removed
#' @param outPrefix a string specifying output file prefix for plots (outPrefix.pdf)
#' @param transform logical, whether to log2 transform the data matrix
#' @param verbose logical, whether to output additional information (data used for plotting or statistical test results)
#' @param shade logical indicating whether to place a shaded rectangle around the point of interest
#' @param nc integer, number of cores for parallel processing
#'
#' @return a list of two objects, the first is a GRanges object, the second is a GRangesList object
#' @author Shuye Pu
#'
#' @examples
#' centerFiles <- c(system.file("data", "test_B.bed", package="GenomicPlotData"),
#' system.file("data", "test_C.bed", package="GenomicPlotData"))
#' names(centerFiles) <- c("TestB", "TestC")
#'
#' @export plot_start_end_reference_region
#'

plot_start_end_reference_region <- function(queryFiles, inputFiles=NULL, centerFiles, handleInputParams=NULL, binSize=10, insert=0, verbose=FALSE, ext=c(-500,100, -100, 500), hl=c(-50, 50, -50, 50), stranded=TRUE, scale=FALSE, smooth=FALSE, rmOutlier=FALSE, outPrefix="plots", transform=TRUE, shade=TRUE, Ylab="Signal intensity", nc=2){

   queryLabels <- names(queryFiles)
   names(queryLabels) <- queryFiles
   inputLabels <- names(inputFiles)
   names(inputLabels) <- inputFiles


   if(!is.null(inputFiles)){
      if(length(queryFiles) != length(inputFiles)){
         stop("the number of queryFiles and inputFiles must be the same!")
      }
      names(inputLabels) <- inputFiles
      queryFiles <- c(queryFiles, inputFiles)
      queryLabels <- c(queryLabels, inputLabels)
   }

   if(!is.null(outPrefix)){
      pdf(paste0(outPrefix, ".pdf"), height=10, width=8)
   }

   queryInputs <- handle_input(inputFiles=queryFiles, handleInputParams, nc=nc)

   bedparam <- handleInputParams
   bedparam$CLIP_reads <- FALSE
   bedparam$fix_width <- 0
   bedparam$useScore <- FALSE
   bedparam$outRle <- FALSE
   bedparam$useSizeFactor <- FALSE
   features <- handle_input(inputFiles=centerFiles, bedparam)
   featureNames <- names(centerFiles)

   ext[2] <- ext[2] - (ext[2]-ext[1])%%binSize ## to avoid binSize inconsistency, as the final binSize is dictated by bin_num
   bin_num_s <- round((ext[2]-ext[1])/binSize)
   ext[4] <- ext[4] - (ext[4]-ext[3])%%binSize ## to avoid binSize inconsistency, as the final binSize is dictated by bin_num
   bin_num_e <- round((ext[4]-ext[3])/binSize)
   bin_num_c <- round(insert/binSize)

   scoreMatrix_lists <- list()
   mat_lists <- list()
   plot_df <- NULL
   for(i in seq_along(featureNames)){
      feature <- features[[i]]$query
      featureName <- featureNames[i]
      nf <- length(feature)
      if(verbose) print(paste("number of features: ", featureName, nf))

      fs <- promoters(resize(feature,width=1,fix="start"), upstream=-ext[1],downstream=ext[2])
      fe <- promoters(resize(feature,width=1,fix="end"), upstream=-ext[3],downstream=ext[4])
      fc <- promoters(resize(feature,width=1,fix="center"), upstream=round(insert/2),downstream=round(insert/2))

      mat_list <- list()
      mat_list[["Start"]] <- list("window"=fs, s=ext[1], e=ext[2], "xmin"=hl[1], "xmax"=hl[2], "bin_num"=bin_num_s)
      mat_list[["Center"]] <- list("window"=fc, s=-round(insert/2), e=round(insert/2), "xmin"=0, "xmax"=0, "bin_num"=bin_num_c)
      mat_list[["End"]] <- list("window"=fe,  s=ext[3], e=ext[4], "xmin"=hl[3], "xmax"=hl[4], "bin_num"=bin_num_e)

      mat_lists[[featureName]] <- mat_list

      scoreMatrix_list <- list()

      for(locus in names(mat_list)){
         windowR <- mat_list[[locus]]$window
         bin_num <- mat_list[[locus]]$bin_num
         if(bin_num <= 0) next

         for(queryFile in queryFiles){

            queryLabel <- queryLabels[queryFile]
            if(verbose) print(queryLabel)
            queryRegions <- queryInputs[[queryFile]]$query
            libsize <- queryInputs[[queryFile]]$size

            bin_op <- "mean"
            weight_col <- queryInputs[[queryFile]]$weight

            fullMatrix <- parallel_scoreMatrixBin(queryRegions, windowR, bin_num, bin_op, weight_col, stranded, nc=nc)

            scoreMatrix_list[[queryLabel]][[locus]] <- fullMatrix

         }
      }

      scoreMatrix_lists[[featureName]] <- scoreMatrix_list

      for(locus in names(mat_list)){

         xmin <- mat_list[[locus]]$xmin
         xmax <- mat_list[[locus]]$xmax
         bin_num <- mat_list[[locus]]$bin_num
         start <- mat_list[[locus]]$s
         end <- mat_list[[locus]]$e

         if(bin_num <= 0) next
         for(queryFile in queryFiles){

            queryLabel <- queryLabels[queryFile]
            if(verbose) print(queryLabel)

            fullMatrix <- scoreMatrix_list[[queryLabel]][[locus]]
            fullMatrix <- process_scoreMatrix(fullMatrix, scale, rmOutlier, transform=transform, verbose=verbose)

            colm <- apply(fullMatrix, 2, mean)
            colsd <- apply(fullMatrix, 2, sd)
            colse <- colsd/sqrt(nrow(fullMatrix))
            collabel <- seq(start, (end-binSize), binSize)
            querybed <- as.factor(rep(queryLabel, ncol(fullMatrix)))
            location <- as.factor(rep(locus, ncol(fullMatrix)))
            levels(location) <- rev(levels(location))
            featurename <- as.factor(rep(featureName, ncol(fullMatrix)))
            Xmin <- rep(xmin, ncol(fullMatrix))
            Xmax <- rep(xmax, ncol(fullMatrix))
            halfmin <- min(fullMatrix)
            intervals <- apply(fullMatrix, 2, function(x) length(x[x>halfmin]))

            sub_df <- NULL
            sub_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Interval"=intervals, "Position"=collabel, "Query"=querybed, "Location"=location, "Feature"=featurename)
            if(smooth){
               sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df=as.integer(bin_num/5))$y)
               sub_df$se <- as.vector(smooth.spline(sub_df$se, df=as.integer(bin_num/5))$y)
            }
            sub_df <- mutate(sub_df, lower=Intensity-se, upper=Intensity+se)

            plot_df <- rbind(plot_df, sub_df)

         }
      }
   }

   Ylab <- ifelse(transform, paste0("Log2 (", Ylab, ")"), Ylab)
   ## plot multi feature lines for one query
   for (query in unique(plot_df$Query)){

      qplot_df <- plot_df %>%
         filter(Query == query)

      plots <- draw_stacked_profile(plot_df=qplot_df, cn="Feature", ext=ext, hl=hl, atitle=query, insert=insert, Ylab=Ylab, shade=shade)

      print(plots)
   }
   ## plot multi query lines for one feature
   for (feature in unique(plot_df$Feature)){

      fplot_df <- plot_df %>%
         filter(Feature == feature)

      plots <- draw_stacked_profile(plot_df=fplot_df, cn="Query", ext=ext, hl=hl, atitle=feature, insert=insert, Ylab=Ylab, shade=shade)

      print(plots)
   }

   ## compute and plot ratio over input
   if(!is.null(inputFiles)){
      Ylab <- ifelse(transform, "Log2 (Ratio-over-Input)", "Ratio-over-Input")

      plot_df <- NULL
      for(i in seq_along(featureNames)){
         feature <- features[[i]]$query
         featureName <- featureNames[i]

         scoreMatrix_list <- scoreMatrix_lists[[featureName]]
         mat_list <- mat_lists[[featureName]]

         ratiofiles <- queryFiles[!queryFiles %in% inputFiles]
         ratiolabels <- queryLabels[!queryLabels %in% inputLabels]
         inputMatrix_list <- scoreMatrix_list[inputLabels]
         ratioMatrix_list <- scoreMatrix_list[ratiolabels]


         for(locus in names(mat_list)){

            bin_num <- mat_list[[locus]]$bin_num
            if(bin_num <= 0) next
            for(i in seq_along(ratiolabels)){
               rm <- ratioMatrix_list[[ratiolabels[i]]][[locus]]
               im <- inputMatrix_list[[inputLabels[i]]][[locus]]
               minrow <- min(nrow(rm), nrow(im))

               fullMatrix <- rm[1:minrow,]/im[1:minrow,]
               fullMatrix <- process_scoreMatrix(fullMatrix, scale=FALSE, rmOutlier, transform=transform, verbose=verbose)

               ratioMatrix_list[[ratiolabels[i]]][[locus]] <- fullMatrix
            }
         }

         for(locus in names(mat_list)){

            xmin <- mat_list[[locus]]$xmin
            xmax <- mat_list[[locus]]$xmax
            bin_num <- mat_list[[locus]]$bin_num
            start <- mat_list[[locus]]$s
            end <- mat_list[[locus]]$e

            if(bin_num <= 0) next
            for(ratiofile in ratiofiles){

               ratiolabel <- ratiolabels[ratiofile]
               if(verbose) print(ratiolabel)

               fullMatrix <- ratioMatrix_list[[ratiolabel]][[locus]]

               colm <- apply(fullMatrix, 2, mean)
               colsd <- apply(fullMatrix, 2, sd)
               colse <- colsd/sqrt(nrow(fullMatrix))
               collabel <- seq(start, (end-binSize), binSize)
               ratiobed <- as.factor(rep(ratiolabel, ncol(fullMatrix)))
               location <- as.factor(rep(locus, ncol(fullMatrix)))
               featurename <- as.factor(rep(featureName, ncol(fullMatrix)))
               levels(location) <- rev(levels(location))
               halfmin <- min(fullMatrix)

               intervals <- apply(fullMatrix, 2, function(x) length(x[x>halfmin]))

               sub_df <- NULL
               sub_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Interval"=intervals,  "Position"=collabel, "Query"=ratiobed, "Location"=location, "Feature"=featurename)
               if(smooth){
                  sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df=as.integer(bin_num/5))$y)
                  sub_df$se <- as.vector(smooth.spline(sub_df$se, df=as.integer(bin_num/5))$y)
               }
               sub_df <- mutate(sub_df, lower=Intensity-se, upper=Intensity+se)
               plot_df <- rbind(plot_df, sub_df)
            }
         }
      }

      ## plot multi feature lines for one query
      for (query in unique(plot_df$Query)){

         qplot_df <- plot_df %>%
            filter(Query == query)

         plots <- draw_stacked_profile(plot_df=qplot_df, cn="Feature", ext=ext, hl=hl, atitle=query, insert=insert, Ylab=Ylab, shade=shade)

         print(plots)
      }
      ## plot multi query lines for one feature
      for (feature in unique(plot_df$Feature)){

         fplot_df <- plot_df %>%
            filter(Feature == feature)

         plots <- draw_stacked_profile(plot_df=fplot_df, cn="Query", ext=ext, hl=hl, atitle=feature, insert=insert, Ylab=Ylab, shade=shade)

         print(plots)
      }
   }
   if(!is.null(outPrefix)){
      on.exit(dev.off(), add=TRUE)
   }

   invisible(plot_df)
}


#' @title Plot promoter, gene body and TTS
#
#' @description Plot reads or peak signal intensity of samples in the query files around genes. The upstream and downstream windows flanking genes
#' can be given separately, the parameter 'meta' controls if gene or metagene plots are generated. If Input files are provided, ratio over Input
#' is computed and displaued as well.
#'
#' @param queryFiles a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param gFeatures genomic features as output of the function 'prepare_3parts_genomic_features'
#' @param inputFiles a vector of input sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param handleInputParams a list of parameters for \code{handle_input}
#' @param stranded logical, indicating whether the strand of the feature should be considered
#' @param scale logical, indicating whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smooth logical, indicating whether the line should smoothed with a spline smoothing algorithm
#' @param heatmap logical, indicating whether a heatmap of the score matrix should be generated
#' @param heatRange a numerical vector of two elements, defining range for heatmap color ramp generation
#' @param rmOutlier logical, indicating whether a row with abnormally high values in the score matrix should be removed
#' @param outPrefix a string specifying output file prefix for plots (outPrefix.pdf)
#' @param transform logical, whether to log2 transform the data matrix
#' @param verbose logical, whether to output additional information (data used for plotting or statistical test results)
#' @param nc integer, number of cores for parallel processing
#'
#' @return a dataframe containing the data used for plotting
#' @author Shuye Pu
#'
#' @examples
#' centerFiles <- c(system.file("data", "test_B.bed", package="GenomicPlotData"),
#' system.file("data", "test_C.bed", package="GenomicPlotData"))
#' names(centerFiles) <- c("TestB", "TestC")
#'
#' @export plot_3parts_metagene


plot_3parts_metagene <- function(queryFiles, gFeatures, inputFiles=NULL, scale=FALSE,  verbose=FALSE, Ylab="Signal intensity", handleInputParams=NULL, smooth=FALSE, stranded=TRUE, outPrefix="plots", heatmap=FALSE, rmOutlier=FALSE, heatRange=NULL, transform=FALSE, nc=2){


  queryLabels <- names(queryFiles)
  names(queryLabels) <- queryFiles
  inputLabels <- names(inputFiles)
  names(inputLabels) <- inputFiles

  if(!is.null(inputFiles)){
    if(length(queryFiles) != length(inputFiles)){
      stop("the number of queryFiles and inputFiles must be the same!")
    }
    names(inputLabels) <- inputFiles
    queryFiles <- c(queryFiles, inputFiles)
    queryLabels <- c(queryLabels, inputLabels)
  }


  if(!is.null(outPrefix)){
     pdf(paste(outPrefix, "pdf", sep="."), height=8, width=12)
  }

  windowRs <- gFeatures$windowRs
  featureNames <- names(windowRs)
  if(verbose) print(featureNames)
  if(verbose) print(vapply(windowRs, length, numeric(1)))

  nbins <- gFeatures$nbins
  scaled_bins  <- gFeatures$scaled_bins
  meta <- gFeatures$meta
  fiveP <- gFeatures$fiveP
  threeP <- gFeatures$threeP

  ## start overlapping

  if(verbose) print("computing coverage for queryFiles")
  scoreMatrix_list <- list()

  queryInputs <- handle_input(queryFiles, handleInputParams, nc=nc)

  for(queryFile in queryFiles){

    queryLabel <- queryLabels[queryFile]

    myInput <- queryInputs[[queryFile]]
    libsize <- myInput$size
    queryRegions <- myInput$query
    fileType <- myInput$type
    weight_col <- myInput$weight

    for(w in featureNames){
      if(verbose) print(w)
      windowR <- as(windowRs[[w]], "GRangesList")
      bin_num <- scaled_bins[w]

      bin_op <- "mean"
      fullMatrix <- parallel_scoreMatrixBin(queryRegions, windowR, bin_num, bin_op, weight_col, stranded, nc=nc)

      scoreMatrix_list[[queryLabel]][[w]] <- fullMatrix
    }
  }

  if(verbose) print("plotting coverage for queryFiles")
  mplot_df <- NULL
  vx <- c(1, cumsum(scaled_bins[1:(length(scaled_bins)-1)])+1) ## x axis points for vlines that demarcate the genomic features
  names(vx) <- featureNames

  Ylab <- ifelse(transform, paste0("Log2 (", Ylab, ")"), Ylab)

  if(heatmap) heatmap_list <- list()

  for(queryFile in queryFiles){
    queryLabel <- queryLabels[queryFile]

    plot_df <- NULL

    dims <- vapply(scoreMatrix_list[[queryLabel]], dim, numeric(2))

    if(any(dims[1,] != dims[1,1])){
       message(paste(dims[1,], collapse = " "))
       stop("Number of genes are not equal among features, make sure all feature windows are within chromosome lengths of query regions,
            as genomation will remvove all feature windows outside chromosome boundaries")
    }else{
       featureMatrix <- as.matrix(bind_cols(scoreMatrix_list[[queryLabel]]))
       featureMatrix <- process_scoreMatrix(featureMatrix, scale=scale, rmOutlier=rmOutlier, transform=transform, verbose=verbose)
       colm <- apply(featureMatrix, 2, mean)
       colsd <- apply(featureMatrix, 2, sd)
       colse <- colsd/sqrt(nrow(featureMatrix))
       querybed <- rep(queryLabel, ncol(featureMatrix))
       collabel <- list()
       featuretype <- list()
       for(w in featureNames){
          if(verbose) print(w)
          if(scaled_bins[w] > 0){
             bin_num <- scaled_bins[w]
             collabel[[w]] <- seq(vx[w], vx[w]+bin_num-1)
             featuretype[[w]] <- rep(w, bin_num)
          }
       }
       collabel <- unlist(collabel)
       featuretype <- unlist(featuretype)
       names(collabel) <- featuretype

       if(heatmap){
          dataname <- paste(Ylab, queryLabel, "gene", sep=":")
          heatmap_list[dataname] <- draw_matrix_heatmap(featureMatrix, dataName=dataname, labels_col=collabel, levels_col=featureNames, ranges=heatRange, verbose=verbose)
       }

       plot_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=querybed, "Feature"=featuretype)
    }

    if(smooth){
      plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df=as.integer(nbins/5))$y)
      plot_df$se <- as.vector(smooth.spline(plot_df$se, df=as.integer(nbins/5))$y)
    }

    mplot_df <- rbind(mplot_df, plot_df)

  }

  mplot_df <- mutate(mplot_df, lower=Intensity-se, upper=Intensity+se)

  xmax <- max(mplot_df$Position)
  pp <- draw_region_landmark(featureNames, vx, xmax)
  ppp <- draw_region_name(featureNames, scaled_bins, xmax)
  marker <- plot_grid(pp, ppp, ncol = 1, align = 'v', axis= "lr", rel_heights = c(1,2))
  ## plot individual sample lines with error band
  plot_list <- list()
  for(queryLabel in queryLabels){
     aplot_df <- mplot_df %>%
        filter(Query == queryLabel)
     p <- draw_region_profile(plot_df=aplot_df, cn="Query", vx=vx, Ylab=Ylab)
     outp <- plot_grid(p, marker, ncol = 1, align = 'v', axis= "lr", rel_heights = c(10,1))
     plot_list[[queryLabel]] <- outp
  }
  rowp <- plot_grid(plotlist=plot_list, nrow = 1, align = 'h')
  #print(rowp)

  if(heatmap){
     groblist <- lapply(heatmap_list, function(x)grid.grabExpr(draw(x, heatmap_legend_side = "top")))
     heatp <- plot_grid(plotlist=groblist, nrow = 1, align = 'h')

     composite <- plot_grid(rowp, heatp, ncol = 1, align = 'v')
     print(composite)
  }else{
     print(rowp)
  }

  ## plot multi-sample lines with error band
  p <- draw_region_profile(plot_df=mplot_df, cn="Query", vx=vx, Ylab=Ylab)
  outp <- plot_grid(p, marker, ncol = 1, align = 'v', axis="lr", rel_heights = c(10,1))
  print(outp)

  if(!is.null(inputFiles)){
    if(verbose) print("computing coverage for ratio over input")
    Ylab <- ifelse(transform, "Log2 (Ratio-over-Input)", "Ratio-over-Input")
    if(heatmap) heatmap_list <- list()

    inputMatrix_list <- scoreMatrix_list[inputLabels]
    ratiofiles <- queryFiles[!queryFiles %in% inputFiles]
    ratiolabels <- queryLabels[!queryLabels %in% inputLabels]
    ratioMatrix_list <- scoreMatrix_list[ratiolabels]
    for(w in featureNames){
      for(i in seq_along(ratiolabels)){
        fullMatrix <- ratioMatrix_list[[ratiolabels[i]]][[w]]/inputMatrix_list[[inputLabels[i]]][[w]]
        ratioMatrix_list[[ratiolabels[i]]][[w]] <- fullMatrix
      }
    }

    if(verbose) print("plotting coverage for ratio over input")
    mplot_df <- NULL
    for(ratiofile in ratiofiles){
      ratiolabel <- ratiolabels[ratiofile]

      plot_df <- NULL
      dims <- vapply(ratioMatrix_list[[ratiolabel]], dim, numeric(2))

      if(any(dims[1,] != dims[1,1])){
         message(paste(dims[1,], collapse = " "))
         stop("Number of genes are not equal among features, make sure all feature windows are within chromosome lengths of query regions,
            as genomation will remvove all feature windows outside chromosome boundaries")
      }else{
         featureMatrix <- as.matrix(bind_cols(ratioMatrix_list[[ratiolabel]]))
         featureMatrix <- process_scoreMatrix(featureMatrix, scale=FALSE, rmOutlier=rmOutlier, transform=transform, verbose=verbose) # do not scale ration over input
         colm <- apply(featureMatrix, 2, mean)
         colsd <- apply(featureMatrix, 2, sd)
         colse <- colsd/sqrt(nrow(featureMatrix))
         querybed <- rep(ratiolabel, ncol(featureMatrix))
         collabel <- list()
         featuretype <- list()
         for(w in featureNames){
            if(verbose) print(w)
            if(scaled_bins[w] > 0){
               bin_num <- scaled_bins[w]
               collabel[[w]] <- seq(vx[w], vx[w]+bin_num-1)
               featuretype[[w]] <- rep(w, bin_num)
            }
         }
         collabel <- unlist(collabel)
         featuretype <- unlist(featuretype)
         names(collabel) <- featuretype

         if(heatmap){
            dataname <- paste(Ylab, ratiolabel, "gene", sep=":")
            heatmap_list[dataname] <- draw_matrix_heatmap(featureMatrix, dataName=dataname, labels_col=collabel, levels_col=featureNames, ranges=heatRange, verbose=verbose)
         }
         plot_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=querybed, "Feature"=featuretype)
      }

      if(smooth){
        plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df=as.integer(nbins/5))$y)
        plot_df$se <- as.vector(smooth.spline(plot_df$se, df=as.integer(nbins/5))$y)
      }

      mplot_df <- rbind(mplot_df, plot_df)

    }

    mplot_df <- mutate(mplot_df, lower=Intensity-se, upper=Intensity+se)

    ## plot individual sample lines with error band
    plot_list <- list()
    for(ratiolabel in ratiolabels){
       aplot_df <- mplot_df %>%
          filter(Query == ratiolabel)
       p <- draw_region_profile(plot_df=aplot_df, cn="Query", vx=vx, Ylab=Ylab)
       outp <- plot_grid(p, marker, ncol = 1, align = 'v', axis= "lr", rel_heights = c(10,1))
       plot_list[[ratiolabel]] <- outp
    }
    rowp <- plot_grid(plotlist=plot_list, nrow = 1, align = 'h')
    #print(rowp)

    if(heatmap){
       groblist <- lapply(heatmap_list, function(x)grid.grabExpr(draw(x, heatmap_legend_side = "top")))
       heatp <- plot_grid(plotlist=groblist, nrow = 1, align = 'h')
       composite <- plot_grid(rowp, heatp, ncol = 1, align = 'v')
       print(composite)
    }else{
       print(rowp)
    }

    ## plot multi-sample lines with error band
    p <- draw_region_profile(plot_df=mplot_df, cn="Query", vx=vx, Ylab=Ylab)
    outp <- plot_grid(p, marker, ncol = 1, align = 'v', axis="lr", rel_heights = c(10,1))
    print(outp)

  }

  if(!is.null(outPrefix)){
    on.exit(dev.off(), add=TRUE)
  }

  invisible(mplot_df)
}

#' @title Plot signal inside as well as around custom genomic regions
#'
#' @description Plot reads or peak signal intensity of samples in the query files inside regions defined in the centerFiles. The upstream and downstream flanking windows
#' can be given separately. If Input files are provided, ratio over Input is computed and displaued as well.
#'
#' @param queryFiles a named vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param centerFiles a named vector of reference file names or genomic features in  c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene"). The file should be in .bed format only
#' @param txdb a TxDb object defined in the GenomicFeatures package.
#' @param inputFiles a named vector of input sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param nbins an integer defines the total number of bins
#' @param fiveP an integer, indicating extension out or inside of the 5' boundary of gene by negative or positive number 
#' @param threeP an integer, indicating extension out or inside of the 5' boundary of gene by positive or negative number
#' @param stranded logical, indicating whether the strand of the feature should be considered
#' @param scale logical, indicating whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smooth logical, indicating whether the line should smoothed with a spline smoothing algorithm
#' @param heatmap logical, indicating whether a heatmap of the score matrix should be generated
#' @param heatRange a numerical vector of two elements, defining range for heatmap color ramp generation
#' @param rmOutlier logical, indicating whether a row with abnormally high values in the score matrix should be removed
#' @param handleInputParams a list of parameters for \code{handle_input}
#' @param outPrefix a string specifying output file prefix for plots (outPrefix.pdf)
#' @param regionName a string specifying the name of the center region in the plots
#' @param transform logical, whether to log2 transform the data matrix
#' @param statsMethod a string in c("wilcox.test", "t.test"), for pair-wise group comparisons
#' @param verbose logical, indicating whether to output additional information (data used for plotting or statistical test results)
#' @param nc integer, number of cores for parallel processing
#'
#' @return a dataframe containing the data used for plotting
#' @author Shuye Pu
#'
#' @examples
#' centerFiles <- c(system.file("data", "test_B.bed", package="GenomicPlotData"),
#' system.file("data", "test_C.bed", package="GenomicPlotData"))
#' names(centerFiles) <- c("TestB", "TestC")
#'
#' @export plot_reference_region

plot_reference_region <- function(queryFiles, centerFiles, txdb=NULL, regionName="region", inputFiles=NULL, nbins=100, handleInputParams=NULL, verbose=FALSE, scale=FALSE, heatmap=FALSE, fiveP=-1000, threeP=1000, smooth=FALSE, stranded=TRUE, transform=FALSE, outPrefix="plots", rmOutlier=FALSE, heatRange=NULL, Ylab="Signal intensity", statsMethod="wilcox.test", nc=2){

  if(!is.null(outPrefix)){
    pdf(paste(outPrefix, "pdf", sep="."), height=8, width=12)
  }

  queryLabels <- names(queryFiles)
  names(queryLabels) <- queryFiles
  inputLabels <- names(inputFiles)
  names(inputLabels) <- inputFiles
  

  if(!is.null(inputFiles)){
    if(length(queryFiles) != length(inputFiles)){
      stop("the number of queryFiles and inputFiles must be the same!")
    }
    names(inputLabels) <- inputFiles
    queryFiles <- c(queryFiles, inputFiles)
    queryLabels <- c(queryLabels, inputLabels)
  }

  queryInputs <- handle_input(queryFiles, handleInputParams, nc=nc)

  bedparam <- handleInputParams
  bedparam$CLIP_reads <- FALSE
  bedparam$fix_width <- 0
  bedparam$useScore <- FALSE
  bedparam$outRle <- FALSE
  bedparam$useSizeFactor <- FALSE
  
  centerInputs <- list()
  centerLabels <- vector()
  
  for(featureName in centerFiles){
     if(featureName %in% c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene")){
        feature <- get_genomic_feature_coordinates(txdb, featureName, longest=TRUE, protein_coding=TRUE)[["GRanges"]]
        featureGR <- list("query"=feature)
        centerInputs[[featureName]] <- featureGR
        
        centerLabels[featureName] <- paste0(featureName, "Region")
        
     }else if(file.exists(featureName)){
        feature <- handle_input(featureName, bedparam, nc=nc)[[1]]
        centerInputs[[featureName]] <- feature
        centerLabels[featureName] <- names(centerFiles[centerFiles==featureName])
     }else{
        stop(paste(featureName, "is not supported!"))
     }
   }
  

  five <- fiveP/1000
  five <- paste0(five, "K")
  fiveL <- -fiveP
  if(fiveP>=0){fiveL <- 0}
  three <- threeP/1000
  three <- paste0(three, "K")
  threeL <- threeP
  if(threeP<=0){threeL <- 0}

  featureNames <-  c(five, regionName, three)
  #featureNames <- featureNames[featureNames != ""]

  all_regions <- unlist(as(lapply(centerInputs, function(x)x$query), "GRangesList"))
  regionLen <- median(width(all_regions))
  lens <- c(fiveL, regionLen, threeL)
  
  scaled_bins <- round(lens*nbins/sum(lens))
  #scaled_bins <- scaled_bins[scaled_bins > 0]
  names(scaled_bins) <- featureNames

  if(verbose) print("computing coverage for Sample")
  scoreMatrix_list <- list()

  for(queryFile in queryFiles){
    queryLabel <- queryLabels[queryFile]
    if(verbose) print(paste("queryLabel", queryLabel))

    myInput <- queryInputs[[queryFile]]
    libsize <- myInput$size
    queryRegions <- myInput$query
    fileType <- myInput$type
    weight_col <- myInput$weight

    for(centerFile in centerFiles){
      centerLabel <- centerLabels[centerFile]
      if(verbose) print(paste("centerLabel", centerLabel))
      centerInput <- centerInputs[[centerFile]]
      centerGr <- centerInput$query
      centerGr <- check_constraints(centerGr, handleInputParams$genome, queryRegions)
      
      if(fiveP > 0 && threeP < 0){ # to avoid generating negative width in the narrow function below
         centerGr <- centerGr[width(centerGr) > (abs(fiveP) + abs(threeP))]
      }else if(fiveP > 0){
         centerGr <- centerGr[width(centerGr) > abs(fiveP)]
      }else if(threeP < 0){
         centerGr <- centerGr[width(centerGr) > abs(threeP)]
      }
      centerGrPlus <- centerGr[as.vector(strand(centerGr)) %in% c("+", "*")]
      centerGrMinus <- centerGr[as.vector(strand(centerGr)) == "-"]

      upstreamGr <- flank(centerGr, width=fiveL, start=TRUE, both=FALSE, use.names=TRUE, ignore.strand=FALSE)
      downstreamGr <- flank(centerGr, width=threeL, start=FALSE, both=FALSE, use.names=TRUE, ignore.strand=FALSE)
      upstreamGr <- check_constraints(upstreamGr, handleInputParams$genome, queryRegions)
      downstreamGr <- check_constraints(downstreamGr, handleInputParams$genome, queryRegions)
      
      if(fiveP > 0){
         centerGrPlus <- narrow(centerGrPlus, start=fiveP, end=NA)
         centerGrMinus <- narrow(centerGrMinus, start=NA, end=fiveP)
      }
      if(threeP < 0){
         centerGrPlus <- narrow(centerGrPlus, start=NA, end=threeP)
         centerGrMinus <- narrow(centerGrMinus, start=-threeP, end=NA)
      }
      centerGr <- c(centerGrPlus, centerGrMinus)
      centerGr <- centerGr[width(centerGr) >= scaled_bins[regionName]]
      
      
      if(fiveP < 0 && threeP > 0){
         commonNames <- Reduce(intersect, list(names(upstreamGr), names(downstreamGr), names(centerGr)))
         
         upstreamGr <- upstreamGr[commonNames]
         downstreamGr <- downstreamGr[commonNames]
         centerGr <- centerGr[commonNames]  
      }
      
      windowRegions <- split(centerGr, f=factor(names(centerGr)))
      windowUp <- split(upstreamGr, f=factor(names(upstreamGr)))
      windowDown <- split(downstreamGr, f=factor(names(downstreamGr)))
      windowRs <- list(windowUp, windowRegions, windowDown)
      names(windowRs) <- featureNames

      for(w in featureNames){
        #w <- "utr5"
        if(verbose) print(paste("featureName", w))
        if(scaled_bins[w] > 0){
           windowR <- as(windowRs[[w]], "GRangesList")
           bin_num <- scaled_bins[w]
           bin_op <- "mean"
           fullMatrix <- parallel_scoreMatrixBin(queryRegions, windowR, bin_num, bin_op, weight_col, stranded, nc=nc)
           rownames(fullMatrix) <- names(windowR)
           scoreMatrix_list[[queryLabel]][[centerLabel]][[w]] <- fullMatrix
        }
      }
    }
  }

  mplot_df <- list()
  stat_df <- list()
  vx <- c(1, cumsum(scaled_bins[1:(length(scaled_bins)-1)])+1) ## x-axis points for vlines that demarcate the genomic features
  names(vx) <- featureNames
  
  if(heatmap)heatmap_list <- list()
  Ylab <- ifelse(transform, paste0("Log2 (", Ylab, ")"), Ylab)

  if(verbose) print("plotting coverage profiles")
  for(queryFile in queryFiles){
    queryLabel <- queryLabels[queryFile]
    if(verbose) print(paste("queryLabel",queryLabel))
    for(centerFile in centerFiles){
      centerLabel <- centerLabels[centerFile]
      if(verbose) print(paste("centerLabel", centerLabel))
      plot_df <- NULL

      dims <- vapply(scoreMatrix_list[[queryLabel]][[centerLabel]], dim, numeric(2))

      if(any(dims[1,] != dims[1,1])){
         message(paste(dims[1,], collapse = " "))
         stop("Number of genes are not equal among features, make sure all feature windows are within chromosome lengths of query regions,
            as genomation will remvove all feature windows outside chromosome boundaries")
      }else{
         featureMatrix <- as.matrix(bind_cols(scoreMatrix_list[[queryLabel]][[centerLabel]]))
         rownames(featureMatrix) <- rownames(scoreMatrix_list[[queryLabel]][[centerLabel]][[1]])
         featureMatrix <- process_scoreMatrix(featureMatrix, scale=scale, rmOutlier=rmOutlier, transform=transform, verbose=verbose)
         colm <- apply(featureMatrix, 2, mean)
         colsd <- apply(featureMatrix, 2, sd)
         colse <- colsd/sqrt(nrow(featureMatrix))
         querybed <- rep(queryLabel, ncol(featureMatrix))
         centerbed <- rep(centerLabel, ncol(featureMatrix))
         collabel <- list()
         featuretype <- list()
         for(w in featureNames){
            if(verbose) print(paste("featureName", w))
            if(scaled_bins[w] > 0){
               bin_num <- scaled_bins[w]
               collabel[[w]] <- seq(vx[w], vx[w]+bin_num-1)
               featuretype[[w]] <- rep(w, bin_num)
            }
         }
         collabel <- unlist(collabel)
         featuretype <- unlist(featuretype)
         names(collabel) <- featuretype

         if(heatmap){
           
            dataname <- paste(Ylab, queryLabel, centerLabel, sep=":")
            heatmap_list[dataname] <- draw_matrix_heatmap(featureMatrix, dataName=dataname, labels_col=collabel, levels_col=names(scaled_bins[scaled_bins>0]), ranges=heatRange, verbose=verbose)
         }

         plot_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=querybed, "Reference"=centerbed, "Feature"=featuretype)
      }

      if(smooth){
        plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df=as.integer(nbins/5))$y)
        plot_df$se <- as.vector(smooth.spline(plot_df$se, df=as.integer(nbins/5))$y)
      }

      mplot_df[[paste(queryLabel,centerLabel,sep=":")]] <- plot_df
      
      regionMatrix <- as.matrix(scoreMatrix_list[[queryLabel]][[centerLabel]][[regionName]])
      regionMatrix <- process_scoreMatrix(regionMatrix, scale=scale, rmOutlier=rmOutlier, transform=transform, verbose=verbose)
      Intensity <- as.numeric(rowMeans(regionMatrix))
      
      Query <- as.factor(rep(queryLabel, length(Intensity)))
      Reference <- as.factor(rep(centerLabel, length(Intensity)))
      subdf <- data.frame(Intensity, Query, Reference)
      stat_df[[paste(queryLabel,centerLabel,sep=":")]] <- subdf
    }
  }

  mplot_df <- bind_rows(mplot_df) %>%
     mutate(Group=as.factor(paste0(Query, ":", Reference)), .keep="all") %>%
     mutate(lower=Intensity-se, upper=Intensity+se)
  mstat_df <- bind_rows(stat_df) %>%
     mutate(Group=as.factor(paste0(Query, ":", Reference)), .keep="all")

  plot_list <- list()
  xmax <- max(mplot_df$Position)
  pp <- draw_region_landmark(featureNames, vx, xmax)
  ppp <- draw_region_name(featureNames, scaled_bins, xmax)
  marker <- plot_grid(pp, ppp, ncol = 1, align = 'v', axis= "lr", rel_heights = c(1,2))

  for(i in seq_along(queryLabels)){
    for(beds in combn(queryLabels, i, simplify = FALSE)){
      for(j in seq_along(centerLabels)){
        for(centers in combn(centerLabels, j, simplify = FALSE)){
          if(verbose) print(paste("beds", beds))
          if(verbose) print(paste("centers", centers))

          aplot_df <- mplot_df %>%
            filter(Query %in% beds & Reference %in% centers) %>%
            mutate(Group=paste(Query,Reference,sep=":"), .keep="all")

          ## plot multi-sample lines with error band
          p <- draw_region_profile(plot_df=aplot_df, cn="Group", vx=vx, Ylab=Ylab)
          outp <- plot_grid(p, marker, ncol = 1, align = 'v', axis="lr", rel_heights = c(10,1))
          
          astat_df <- mstat_df %>%
             filter(Query %in% beds & Reference %in% centers)
          astat_df <- astat_df %>% # unify order of factors to get consistent color mapping
             mutate(Query=factor(Query, levels=sort(unique(Query)))) %>%
             mutate(Reference=factor(Reference, levels=sort(unique(Reference))))
          
          
          if((i == 1 && j > 1) || (i > 1 && j == 1)){
                
             if(j > 1){
                p <- draw_region_profile(plot_df=aplot_df, cn="Reference", sn="Query", vx=vx, Ylab=Ylab)
                outp <- plot_grid(p, marker, ncol = 1, align = 'v', axis="lr", rel_heights = c(10,1))
                
                comp <- combn(seq_along(centers),2, simplify=FALSE)
                
                ps1 <- draw_boxplot_logy(stat_df=astat_df, xc="Reference", yc="Intensity",  comp=comp, stats=statsMethod, Ylab=Ylab, logy="linear")
                ps2 <- draw_boxplot_logy(stat_df=astat_df, xc="Reference", yc="Intensity",  comp=comp, stats=statsMethod, Ylab=Ylab, logy="logy")
                ps1_wo_outlier <- draw_boxplot_wo_outlier(stat_df=astat_df, xc="Reference", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab)
                ps1_mean_se <- draw_mean_se_barplot(stat_df=astat_df, xc="Reference", yc="Intensity", comp=comp, Ylab=Ylab)
                prank <- draw_rank_plot(stat_df=astat_df, xc="Reference", yc="Intensity", Ylab=Ylab, ecdf=TRUE, rank=FALSE)
             }else{
                p <- draw_region_profile(plot_df=aplot_df, cn="Query", sn="Reference", vx=vx, Ylab=Ylab)
                outp <- plot_grid(p, marker, ncol = 1, align = 'v', axis="lr", rel_heights = c(10,1))
                comp <- combn(seq_along(beds),2, simplify=FALSE)
                
                ps1 <- draw_boxplot_logy(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab, logy="linear")
                ps2 <- draw_boxplot_logy(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab, logy="logy")
                ps1_wo_outlier <- draw_boxplot_wo_outlier(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab)
                ps1_mean_se <- draw_mean_se_barplot(stat_df=astat_df, xc="Query", yc="Intensity", Ylab=Ylab, comp=comp)
                prank <- draw_rank_plot(stat_df=astat_df, xc="Query", yc="Intensity", Ylab=Ylab, ecdf=TRUE, rank=FALSE)
             }
             
             if(max(i,j)>3){
                lapply(list(outp, ps1_mean_se, prank, ps1, ps2, ps1_wo_outlier), print)
             }else{
                comp1 <- plot_grid(outp, ps1_mean_se, prank, ncol = 3, rel_widths = c(1,1,1))
                comp2 <- plot_grid(ps1, ps2, ps1_wo_outlier, ncol = 3, align='h', axis='b', rel_widths = c(1,1,1))
                print(plot_grid(comp1, comp2, ncol=1, rel_heights=c(1,1)))
             }
          }else if(i == 1 && j == 1){
             plot_list[[paste(Ylab, beds,centers,sep=":")]] <- outp
          }else if(i == length(queryLabels) && j == length(centerLabels)){
             print(outp)
          }
        }
      }
    }
  }


  ## plot individual sample lines with error band

  rowp <- plot_grid(plotlist=plot_list, nrow = 1, align = 'h')
  #print(rowp)

  if(heatmap){
     groblist <- lapply(heatmap_list, function(x)grid.grabExpr(draw(x, heatmap_legend_side = "top")))
     heatp <- plot_grid(plotlist=groblist, nrow = 1, align = 'h')
     composite <- plot_grid(rowp, heatp, ncol=1, rel_heights=c(1,1))
     print(composite)
  }else{
     print(rowp)
  }


  if(!is.null(inputFiles)){
    Ylab <- ifelse(transform, "Log2 (Ratio-over-Input)", "Ratio-over-Input")

    ratiofiles <- queryFiles[!queryFiles %in% inputFiles]
    ratiolabels <- queryLabels[!queryLabels %in% inputLabels]
    inputMatrix_list <- scoreMatrix_list[inputLabels]
    ratioMatrix_list <- scoreMatrix_list[ratiolabels]

    for(centerFile in centerFiles){
      centerLabel <- centerLabels[centerFile]
      for(featureName in featureNames){
         if(scaled_bins[featureName]> 0){
            for(i in seq_along(ratiolabels)){
               rm <- ratioMatrix_list[[ratiolabels[i]]][[centerLabel]][[featureName]]
               im <- inputMatrix_list[[inputLabels[i]]][[centerLabel]][[featureName]]
               commonrow <- intersect(rownames(rm), rownames(im))
               
               fullMatrix <- rm[commonrow,]/im[commonrow,]
               
               ratioMatrix_list[[ratiolabels[i]]][[centerLabel]][[featureName]] <- fullMatrix
            }
         }
      }
    }

    mplot_df <- list()
    stat_df <- list()
    if(heatmap)heatmap_list <- list()
    if(verbose) print("plotting coverage profiles")
    for(ratiofile in ratiofiles){
      ratiolabel <- ratiolabels[ratiofile]
      if(verbose) print(paste("ratiolabel", ratiolabel))
      for(centerFile in centerFiles){
        centerLabel <- centerLabels[centerFile]
        if(verbose) print(paste("centerlabel", centerLabel))
        plot_df <- NULL
        dims <- vapply(scoreMatrix_list[[ratiolabel]][[centerLabel]], dim, numeric(2))

        if(any(dims[1,] != dims[1,1])){
           message(paste(dims[1,], collapse = " "))
           stop("Number of genes are not equal among features, make sure all feature windows are within chromosome lengths of query regions,
            as genomation will remvove all feature windows outside chromosome boundaries")
        }else{
           featureMatrix <- as.matrix(bind_cols(scoreMatrix_list[[ratiolabel]][[centerLabel]]))
           rownames(featureMatrix) <- rownames(scoreMatrix_list[[ratiolabel]][[centerLabel]][[1]])
           featureMatrix <- process_scoreMatrix(featureMatrix, scale=scale, rmOutlier=rmOutlier, transform=transform, verbose=verbose)
           colm <- apply(featureMatrix, 2, mean)
           colsd <- apply(featureMatrix, 2, sd)
           colse <- colsd/sqrt(nrow(featureMatrix))
           ratiobed <- rep(ratiolabel, ncol(featureMatrix))
           centerbed <- rep(centerLabel, ncol(featureMatrix))
           collabel <- list()
           featuretype <- list()
           for(w in featureNames){
              if(verbose) print(w)
              if(scaled_bins[w] > 0){
                 bin_num <- scaled_bins[w]
                 collabel[[w]] <- seq(vx[w], vx[w]+bin_num-1)
                 featuretype[[w]] <- rep(w, bin_num)
              }
           }
           collabel <- unlist(collabel)
           featuretype <- unlist(featuretype)
           names(collabel) <- featuretype

           if(heatmap){
              dataname <- paste(Ylab, ratiolabel, centerLabel, sep=":")
              heatmap_list[dataname] <- draw_matrix_heatmap(featureMatrix, dataName=dataname, labels_col=collabel, levels_col=names(scaled_bins[scaled_bins>0]), ranges=heatRange, verbose=verbose)
           }
           plot_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=ratiobed, "Reference"=centerbed, "Feature"=featuretype)
        }


        if(smooth){
          plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df=as.integer(nbins/5))$y)
          plot_df$se <- as.vector(smooth.spline(plot_df$se, df=as.integer(nbins/5))$y)
        }

        mplot_df[[paste(ratiolabel,centerLabel,sep=":")]] <- plot_df
        
        regionMatrix <- as.matrix(scoreMatrix_list[[queryLabel]][[centerLabel]][[regionName]])
        regionMatrix <- process_scoreMatrix(regionMatrix, scale=scale, rmOutlier=rmOutlier, transform=transform, verbose=verbose)
        Intensity <- as.numeric(rowMeans(regionMatrix))
        
        Query <- as.factor(rep(ratiolabel, length(Intensity)))
        Reference <- as.factor(rep(centerLabel, length(Intensity)))
        subdf <- data.frame(Intensity, Query, Reference)
        stat_df[[paste(ratiolabel,centerLabel,sep=":")]] <- subdf
      }
    }
    
    mplot_df <- bind_rows(mplot_df) %>%
       mutate(Group=as.factor(paste0(Query, ":", Reference)), .keep="all") %>%
       mutate(lower=Intensity-se, upper=Intensity+se)
    mstat_df <- bind_rows(stat_df) %>%
       mutate(Group=as.factor(paste0(Query, ":", Reference)), .keep="all")
    
    plot_list <- list()
    for(i in seq_along(ratiolabels)){
      for(beds in combn(ratiolabels, i, simplify = FALSE)){
        for(j in seq_along(centerLabels)){
          for(centers in combn(centerLabels, j, simplify = FALSE)){
            if(verbose) print(beds)
            if(verbose) print(centers)

            aplot_df <- mplot_df %>%
               filter(Query %in% beds & Reference %in% centers) %>%
               mutate(Group=paste(Query,Reference,sep=":"), .keep="all")
            
            ## plot multi-sample lines with error band
            p <- draw_region_profile(plot_df=aplot_df, cn="Group", vx=vx, Ylab=Ylab)
            outp <- plot_grid(p, marker, ncol = 1, align = 'v', axis="lr", rel_heights = c(10,1))
            
            astat_df <- mstat_df %>%
               filter(Query %in% beds & Reference %in% centers)
            astat_df <- astat_df %>% # unify order of factors to get consistent color mapping
               mutate(Query=factor(Query, levels=sort(unique(Query)))) %>%
               mutate(Reference=factor(Reference, levels=sort(unique(Reference))))
            
            
            if((i == 1 && j > 1) || (i > 1 && j == 1)){
               
               if(j > 1){
                  p <- draw_region_profile(plot_df=aplot_df, cn="Reference", sn="Query", vx=vx, Ylab=Ylab)
                  outp <- plot_grid(p, marker, ncol = 1, align = 'v', axis="lr", rel_heights = c(10,1))
                  
                  comp <- combn(seq_along(centers),2, simplify=FALSE)
                  
                  ps1 <- draw_boxplot_logy(stat_df=astat_df, xc="Reference", yc="Intensity",  comp=comp, stats=statsMethod, Ylab=Ylab, logy="linear")
                  ps2 <- draw_boxplot_logy(stat_df=astat_df, xc="Reference", yc="Intensity",  comp=comp, stats=statsMethod, Ylab=Ylab, logy="logy")
                  ps1_wo_outlier <- draw_boxplot_wo_outlier(stat_df=astat_df, xc="Reference", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab)
                  ps1_mean_se <- draw_mean_se_barplot(stat_df=astat_df, xc="Reference", yc="Intensity", comp=comp, Ylab=Ylab)
                  prank <- draw_rank_plot(stat_df=astat_df, xc="Reference", yc="Intensity", Ylab=Ylab, ecdf=TRUE, rank=FALSE)
               }else{
                  p <- draw_region_profile(plot_df=aplot_df, cn="Query", sn="Reference", vx=vx, Ylab=Ylab)
                  outp <- plot_grid(p, marker, ncol = 1, align = 'v', axis="lr", rel_heights = c(10,1))
                  comp <- combn(seq_along(beds),2, simplify=FALSE)
                  
                  ps1 <- draw_boxplot_logy(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab, logy="linear")
                  ps2 <- draw_boxplot_logy(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab, logy="logy")
                  ps1_wo_outlier <- draw_boxplot_wo_outlier(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab)
                  ps1_mean_se <- draw_mean_se_barplot(stat_df=astat_df, xc="Query", yc="Intensity", Ylab=Ylab, comp=comp)
                  prank <- draw_rank_plot(stat_df=astat_df, xc="Query", yc="Intensity", Ylab=Ylab, ecdf=TRUE, rank=FALSE)
               }
               
               if(max(i,j)>3){
                  lapply(list(outp, ps1_mean_se, prank, ps1, ps2, ps1_wo_outlier), print)
               }else{
                  comp1 <- plot_grid(outp, ps1_mean_se, prank, ncol = 3, rel_widths = c(1,1,1))
                  comp2 <- plot_grid(ps1, ps2, ps1_wo_outlier, ncol = 3, align='h', axis='b', rel_widths = c(1,1,1))
                  print(plot_grid(comp1, comp2, ncol=1, rel_heights=c(1,1)))
               }
            }else if(i == 1 && j == 1){
               plot_list[[paste(Ylab, beds,centers,sep=":")]] <- outp
            }else if(i == length(queryLabels) && j == length(centerLabels)){
               print(outp)
            }
          }
        }
      }
    }


    ## plot individual sample lines with error band

    rowp <- plot_grid(plotlist=plot_list, nrow = 1, align = 'h')
    #print(rowp)

    if(heatmap){
       groblist <- lapply(heatmap_list, function(x)grid.grabExpr(draw(x, heatmap_legend_side = "top")))
       heatp <- plot_grid(plotlist=groblist, nrow = 1, align = 'h')
       composite <- plot_grid(rowp, heatp, ncol=1, rel_widths=c(2,3))
       print(composite)
    }else{
       print(rowp)
    }
  }

  if(!is.null(outPrefix)){
    on.exit(dev.off(), add=TRUE)
  }

  invisible(mplot_df)
}

#' @title Plot promoter, 5'UTR, CDS, 3'UTR and TTS
#'
#' @description Plot reads or peak signal intensity of samples in the query files around genes. The upstream and downstream windows flanking genes can be given separately, metagene plots are generated with 5'UTR, CDS and 3'UTR segments. The length of each segments are prorated according to the median length of each segments. If Input files are provided, ratio over Input is computed and displayed as well.
#'
#' @param queryFiles a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param gFeatures_list a list of genomic features as output of the function 'prepare_5parts_genomic_features'
#' @param inputFiles a vector of input sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param handleInputParams a list of parameters for \code{handle_input}
#' @param stranded logical, indicating whether the strand of the feature should be considered
#' @param scale logical, indicating whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smooth logical, indicating whether the line should smoothed with a spline smoothing algorithm
#' @param heatmap logical, indicating whether a heatmap of the score matrix should be generated
#' @param heatRange a numerical vector of two elements, defining range for heatmap color ramp generation
#' @param rmOutlier logical, indicating whether a row with abnormally high values in the score matrix should be removed
#' @param transform logical, whether to log2 transform the matrix
#' @param outPrefix a string specifying output file prefix for plots (outPrefix.pdf)
#' @param verbose logical, indicating whether to output additional information (data used for plotting or statistical test results)
#' @param nc integer, number of cores for parallel processing
#'
#' @return a dataframe containing the data used for plotting
#' @author Shuye Pu
#'
#' @examples
#' txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlotData"))
#' gf <- prepare_5parts_genomic_features(txdb, meta=FALSE, nbins=100, fiveP=500, threeP=500, longest=TRUE)
#' queryFiles <- c(system.file("data", "test_chip_peak_chr19.narrowPeak", package="GenomicPlotData"),
#' system.file("data", "test_chip_peak_chr19.bed", package="GenomicPlotData"),
#' system.file("data", "test_clip_peak_chr19.bed", package="GenomicPlotData"))
#' names(queryFiles) <- c("narrowPeak", "summitPeak", "iCLIPPeak")
#' op <- NULL
#'
#' handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="center", norm=FALSE, useScore=FALSE,
#'                           outRle=TRUE, useSizeFactor=FALSE, genome="hg19")
#'
#' plot_5parts_metagene(queryFiles, gFeatures=gf, inputFiles=NULL, handleInputParams=handleInputParams,
#'                      verbose=FALSE, smooth=TRUE, scale=FALSE, stranded=TRUE, outPrefix=op, transform=FALSE,
#'                      heatmap =TRUE,  rmOutlier=TRUE, nc=2)
#'
#' @export plot_5parts_metagene
#'

plot_5parts_metagene <- function(queryFiles, gFeatures_list, inputFiles=NULL, handleInputParams=NULL,
                                 verbose=FALSE, transform=FALSE, smooth=FALSE, scale=FALSE, stranded=TRUE,
                                 outPrefix=NULL, heatmap=FALSE, heatRange=NULL, rmOutlier=FALSE,
                                 Ylab="Signal intensity", nc=2){

  queryLabels <- names(queryFiles)
  names(queryLabels) <- queryFiles
  inputLabels <- names(inputFiles)
  names(inputLabels) <- inputFiles

  if(!is.null(inputFiles)){
    if(length(queryFiles) != length(inputFiles)){
      stop("the number of queryFiles and inputFiles must be the same!")
    }
    names(inputLabels) <- inputFiles
    queryFiles <- c(queryFiles, inputFiles)
    queryLabels <- c(queryLabels, inputLabels)
  }

  queryInputs <- handle_input(inputFiles=queryFiles, handleInputParams, verbose=verbose, nc=nc)
  
  if(!is.null(outPrefix)) pdf(paste(outPrefix, "pdf", sep="."), height=8, width=12)

  mplot_dfs <- NULL
  mplot_dfs_ratio <- NULL
  heatmap_list <- list()
  heatmap_list_ratio <- list()
  for(aFeature in names(gFeatures_list)){
     if(verbose) print(paste0("computing coverage for query files in ", aFeature))
     gFeatures <- gFeatures_list[[aFeature]]
     
     windowRs <- gFeatures$windowRs
     featureNames <- names(windowRs)
     if(verbose) print("Number of features:")
     if(verbose) print(vapply(windowRs, length, numeric(1)))
   
     nbins <- gFeatures$nbins
     scaled_bins  <- gFeatures$scaled_bins
     meta <- gFeatures$meta
     fiveP <- gFeatures$fiveP
     threeP <- gFeatures$threeP

     if(verbose) print("Number of scaled bins")
     if(verbose) print(scaled_bins)
   
     scoreMatrix_list <- list()

     for(queryFile in queryFiles){
       queryLabel <- queryLabels[queryFile]
       if(verbose) print(queryFile)
       Input <- queryInputs[[queryFile]]
       libsize <- Input$size
       queryRegions <- Input$query
       fileType <- Input$type
       weight_col <- Input$weight
   
       for(w in featureNames){
         if(verbose) print(w)
         windowR <- windowRs[[w]]
         bin_num <- scaled_bins[w]
   
         bin_op <- "mean"
         if(bin_num > 0){
            #fullMatrix1 <- parallel_scoreMatrixBin(queryRegions, windowR, bin_num, bin_op, weight_col, stranded)
            fullMatrix <- parallel_scoreMatrixBin(queryRegions, windowR, bin_num, bin_op, weight_col, stranded, nc=nc)
            scoreMatrix_list[[queryLabel]][[w]] <- fullMatrix
         }else{
            scoreMatrix_list[[queryLabel]][[w]] <- NULL
         }
   
       }
     }
   
     if(verbose) print("Preparing data for individual plotting")

     vx <- c(1, cumsum(scaled_bins[1:(length(scaled_bins)-1)])+1) ## x axis points for vlines that demarcate the genomic features
     names(vx) <- featureNames
   
     mplot_df <- NULL
     Ylab <- ifelse(transform, paste0("Log2 (", Ylab, ")"), Ylab)
   
     for(queryLabel in queryLabels){
       plot_df <- NULL
       dims <- vapply(scoreMatrix_list[[queryLabel]], dim, numeric(2))
   
       if(any(dims[1,] != dims[1,1])){
          message(paste(dims[1,], collapse = " "))
          stop("Number of genes are not equal among features, make sure all feature windows are within chromosome lengths of query regions,
               as genomation will remvove all feature windows outside chromosome boundaries")
       }else{
          featureMatrix <- as.matrix(bind_cols(scoreMatrix_list[[queryLabel]]))
          featureMatrix <- process_scoreMatrix(featureMatrix, scale=scale, rmOutlier=rmOutlier, transform=transform, verbose=verbose)
          colm <- apply(featureMatrix, 2, mean)
          colsd <- apply(featureMatrix, 2, sd)
          colse <- colsd/sqrt(nrow(featureMatrix))
          querybed <- rep(queryLabel, ncol(featureMatrix))
          collabel <- list()
          featuretype <- list()
          for(w in featureNames){
             #if(verbose) print(w)
             if(scaled_bins[w] > 0){
                bin_num <- scaled_bins[w]
                collabel[[w]] <- seq(vx[w], vx[w]+bin_num-1)
                featuretype[[w]] <- rep(w, bin_num)
             }
          }
          collabel <- unlist(collabel)
          featuretype <- unlist(featuretype)
          names(collabel) <- featuretype

          if(heatmap){
             dataname <- paste(Ylab, queryLabel, aFeature, sep=":")
             heatmap_list[[dataname]] <- draw_matrix_heatmap(featureMatrix, dataName=dataname, labels_col=collabel, levels_col=featureNames, ranges=heatRange, verbose=verbose)
          }
          plot_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=paste(querybed,aFeature, sep=":"), "Feature"=featuretype)
       }
   
       if(smooth){
         plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df=as.integer(nbins/5))$y)
         plot_df$se <- as.vector(smooth.spline(plot_df$se, df=as.integer(nbins/5))$y)
       }
   
       mplot_df <- rbind(mplot_df, plot_df)
     }

     mplot_df <- mutate(mplot_df, lower=Intensity-se, upper=Intensity+se)
     mplot_dfs <- rbind(mplot_dfs, mplot_df)
  
     ## if inputFiles are provided, plot ratio over input
   
     if(!is.null(inputFiles)){
       if(verbose) print("Preparing data for ratio plotting")
       Ylab <- ifelse(transform, "Log2 (Ratio-over-Input)", "Ratio-over-Input")
   
       inputMatrix_list <- scoreMatrix_list[inputLabels]
       ratiofiles <- queryFiles[!queryFiles %in% inputFiles]
       ratiolabels <- queryLabels[!queryLabels %in% inputLabels]
       ratioMatrix_list <- scoreMatrix_list[ratiolabels]
       for(w in featureNames){
         if(scaled_bins[w] > 0){
            for(i in seq_along(ratiolabels)){
              fullMatrix <- ratioMatrix_list[[ratiolabels[i]]][[w]]/inputMatrix_list[[inputLabels[i]]][[w]]
              ratioMatrix_list[[ratiolabels[i]]][[w]] <- fullMatrix
            }
         }else{
            ratioMatrix_list[[ratiolabels[i]]][[w]] <- NULL
         }
       }

       mplot_df <- NULL
       for(ratiolabel in ratiolabels){
         plot_df <- NULL
   
         dims <- vapply(ratioMatrix_list[[ratiolabel]], dim, numeric(2))
   
         if(any(dims[1,] != dims[1,1])){
            message(paste(dims[1,], collapse = " "))
            stop("Number of genes are not equal among features, make sure all feature windows are within chromosome lengths of query regions,
               as genomation will remvove all feature windows outside chromosome boundaries")
         }else{
            featureMatrix <- as.matrix(bind_cols(ratioMatrix_list[[ratiolabel]]))
            featureMatrix <- process_scoreMatrix(featureMatrix, scale=FALSE, rmOutlier=rmOutlier, transform=transform, verbose=verbose)  # do not scale ratio over input
            colm <- apply(featureMatrix, 2, mean)
            colsd <- apply(featureMatrix, 2, sd)
            colse <- colsd/sqrt(nrow(featureMatrix))
            querybed <- rep(ratiolabel, ncol(featureMatrix))
            collabel <- list()
            featuretype <- list()
            for(w in featureNames){
               #if(verbose) print(w)
               if(scaled_bins[w] > 0){
                  bin_num <- scaled_bins[w]
                  collabel[[w]] <- seq(vx[w], vx[w]+bin_num-1)
                  featuretype[[w]] <- rep(w, bin_num)
               }
            }
            collabel <- unlist(collabel)
            featuretype <- unlist(featuretype)
            names(collabel) <- featuretype
   
            if(heatmap){
               dataname <- paste(Ylab, ratiolabel, aFeature, sep=":")
               heatmap_list_ratio[[dataname]] <- draw_matrix_heatmap(featureMatrix, dataName=dataname, labels_col=collabel, levels_col=featureNames, ranges=heatRange, verbose=verbose)
            }
            plot_df <- data.frame("Intensity"=colm, "sd"=colsd, "se"=colse, "Position"=collabel, "Query"=paste(querybed,aFeature,sep=":"), "Feature"=featuretype)
         }
         if(smooth){
           plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df=as.integer(nbins/5))$y)
           plot_df$se <- as.vector(smooth.spline(plot_df$se, df=as.integer(nbins/5))$y)
         }
   
         mplot_df <- rbind(mplot_df, plot_df)
       }

       mplot_df <- mutate(mplot_df, lower=Intensity-se, upper=Intensity+se)
       mplot_dfs_ratio <- rbind(mplot_dfs_ratio, mplot_df)
       
     }
  }

  if(verbose) print("Start plotting")
  
  xmax <- max(mplot_dfs$Position)
  pp <- draw_region_landmark(featureNames, vx, xmax)
  ppp <- draw_region_name(featureNames, scaled_bins, xmax)
  #marker <- plot_grid(pp, ppp, ncol = 1, align = 'v', axis= "lr", rel_heights = c(1,2))
  
  ## plot individual sample lines with error band
  plot_list <- list()
  for(queryLabel in queryLabels){
     for(aFeature in names(gFeatures_list)){
        aplot_df <- mplot_dfs %>%
           filter(Query == paste(queryLabel,aFeature,sep=":"))
        p <- draw_region_profile(plot_df=aplot_df, cn="Query", vx=vx, Ylab=Ylab)
        outp <- plot_grid(p, pp, ppp, ncol = 1, align = 'v', axis= "lr", rel_heights = c(25,1,2))
        plot_list[[paste(queryLabel,aFeature,sep=":")]] <- outp
     }
  }
  rowp <- plot_grid(plotlist=plot_list, nrow = 1, align = 'h', axis= "tb")
  #print(rowp)
  if(heatmap){
     groblist <- lapply(heatmap_list, function(x)grid.grabExpr(draw(x, heatmap_legend_side = "top")))
     heatp <- plot_grid(plotlist=groblist, nrow = 1, align = 'h')
     composite <- plot_grid(rowp, heatp, ncol=1, align = "v")
     print(composite)
  }else{
     print(rowp)
  }
  
  ## plot multi-sample lines with error band
  if(length(queryLabels)*length(gFeatures_list)>1){
     p <- draw_region_profile(plot_df=mplot_dfs, cn="Query", vx=vx, Ylab=Ylab)
     outp <- plot_grid(p, pp, ppp, ncol = 1, align = 'v', axis="lr", rel_heights = c(25,1,2))
     print(outp)
  }
  
  if(!is.null(inputFiles)){
    xmax <- max(mplot_dfs_ratio$Position)
    pp <- draw_region_landmark(featureNames, vx, xmax)
    ppp <- draw_region_name(featureNames, scaled_bins, xmax)
    #marker <- plot_grid(pp, ppp, ncol = 1, align = 'v', axis= "lr", rel_heights = c(1,2))

    ## plot individual sample lines with error band
    plot_list <- list()
    for(ratiolabel in ratiolabels){
       for(aFeature in names(gFeatures_list)){
          aplot_df <- mplot_dfs_ratio %>%
             filter(Query == paste(ratiolabel, aFeature, sep=":"))
          p <- draw_region_profile(plot_df=aplot_df, cn="Query", vx=vx, Ylab=Ylab)
          outp <- plot_grid(p, pp, ppp, ncol = 1, align = 'v', axis="lr", rel_heights = c(25,1,2))
          plot_list[[paste(ratiolabel, aFeature, sep=":")]] <- outp
       }
    }
    rowp <- plot_grid(plotlist=plot_list, nrow = 1, align = "h", axis= "tb")
    #print(rowp)

    if(heatmap){
       groblist <- lapply(heatmap_list_ratio, function(x)grid.grabExpr(draw(x, heatmap_legend_side = "top")))
       heatp <- plot_grid(plotlist=groblist, nrow = 1, align = "v")
       composite <- plot_grid(rowp, heatp, ncol=1)
       print(composite)
    }else{
       print(rowp)
    }

    ## plot multi-sample lines with error band
    if(length(ratiolabels)*length(gFeatures_list)>1){
       p <- draw_region_profile(plot_df=mplot_dfs_ratio, cn="Query", vx=vx, Ylab=Ylab)
       outp <- plot_grid(p, pp, ppp, ncol = 1, align = 'v', axis="lr", rel_heights = c(25,1,2))
       print(outp)
    }
  }

  if(!is.null(outPrefix)){
    on.exit(dev.off(), add=TRUE)
  }

  if(verbose) print("plot_5parts_metagene runs successfully!")
  invisible(mplot_dfs)
  
}


#' @title Plot signal around custom genomic loci
#' @description  Plot reads or peak signal intensity of samples in the query files around reference locus (start, end or center of a genomic region)
#' defined in the centerFiles. The upstream and downstream windows flanking loci can be given separately, a smaller window can be defined to allow
#' statistical comparisons between samples for the same reference, or between references for a given sample. If Input files are provided, ratio over
#' Input is computed and displaued as well.
#'
#' @param queryFiles a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param centerFiles a named vector of reference file names or genomic features in  c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene"). The file should be in .bed format only
#' @param txdb a TxDb object defined in the GenomicFeatures package.
#' @param inputFiles a vector of input sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param handleInputParams a list of parameters for \code{handle_input}
#' @param ext a vector of two integers defining upstream and downstream boundaries of the plot window, flanking the reference locus
#' @param hl a vector of two integers defining upstream and downstream boundaries of the highlight window, flanking the reference locus
#' @param stranded logical, indicating whether the strand of the feature should be considered
#' @param scale logical, indicating whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smooth logical, indicating whether the line should smoothed with a spline smoothing algorithm
#' @param heatmap logical, indicating whether a heatmap of the score matrix should be generated
#' @param heatRange a numerical vector of two elements, defining range for heatmap color ramp generation
#' @param rmOutlier logical, indicating whether a row with abnormally high values in the score matrix should be removed
#' @param outPrefix a string specifying output file prefix for plots (outPrefix.pdf)
#' @param refPoint a string in c("start", "center", "end")
#' @param Xlab a string denotes the label on x-axis
#' @param shade logical indicating whether to place a shaded rectangle around the point of interest
#' @param binSize an integer defines bin size for intensity calculation
#' @param transform logical, whether to log2 transform the data matrix
#' @param statsMethod a string in c("wilcox.test", "t.test"), for pair-wise group comparisons
#' @param verbose logical, indicating whether to output additional information (data used for plotting or statistical test results)
#' @param nc integer, number of cores for parallel processing
#'
#' @return a list of two dataframes containing the data used for plotting and for statistical testing
#' @author Shuye Pu
#'
#' @examples
#' centerFiles <- c(system.file("data", "test_B.bed", package="GenomicPlotData"),
#' system.file("data", "test_C.bed", package="GenomicPlotData"))
#' names(centerFiles) <- c("TestB", "TestC")
#'
#' @export plot_reference_locus


plot_reference_locus <- function(queryFiles, centerFiles, txdb=NULL, ext=c(-100,100), hl=c(0,0), shade=TRUE, smooth=FALSE,
                                 handleInputParams=NULL, verbose=FALSE, binSize=10, refPoint="center", Xlab="Center", 
                                 Ylab="Signal intensity", inputFiles=NULL, stranded=TRUE, heatmap=TRUE, scale=FALSE,
                                 outPrefix=NULL, rmOutlier=FALSE, transform=FALSE, statsMethod="wilcox.test", heatRange=NULL,
                                 nc=2){

  queryLabels <- names(queryFiles)
  names(queryLabels) <- queryFiles
  inputLabels <- names(inputFiles)
  names(inputLabels) <- inputFiles
  

  if(!is.null(inputFiles)){
    if(length(queryFiles) != length(inputFiles)){
      stop("the number of queryFiles and inputFiles must be the same!")
    }
    names(inputLabels) <- inputFiles
    queryFiles <- c(queryFiles, inputFiles)
    queryLabels <- c(queryLabels, inputLabels)
  }

  five <- ext[1]/1000
  five <- paste0(five, "K")
  if(ext[1]==0) five=""
  three <- ext[2]/1000
  three <- paste0(three, "K")
  if(ext[2]==0) three=""
  featureNames <- c(five, Xlab, three)

  ext[2] <- ext[2] - (ext[2]-ext[1])%%binSize ## to avoid binSize inconsistency, as the final binSize is dictated by bin_num
  bin_num <- round((ext[2]-ext[1])/binSize)
  colLabel <- seq(ext[1], (ext[2]-binSize), binSize)
  names(colLabel) <- rep(featureNames, c(sum(colLabel<0), sum(colLabel==0), sum(colLabel>0)))

  if(!is.null(outPrefix)){
     while(!is.null(dev.list())){
        dev.off()
     }
     pdf(paste(outPrefix, "pdf", sep="."), height=12, width=12)
  }

  scoreMatrix_list <- list()

  queryInputs <- handle_input(queryFiles,  handleInputParams, nc=nc)

  bedparam <- handleInputParams
  bedparam$CLIP_reads <- FALSE
  bedparam$fix_width <- 0
  bedparam$useScore <- FALSE
  bedparam$outRle <- FALSE
  bedparam$useSizeFactor <- FALSE
  
  centerInputs <- list()
  centerLabels <- vector()
  
  for(featureName in centerFiles){
     if(featureName %in% c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene")){
        feature <- get_genomic_feature_coordinates(txdb, featureName, longest=TRUE, protein_coding=TRUE)[["GRanges"]]
        featureGR <- list("query"=feature)
        centerInputs[[featureName]] <- featureGR
        
        centerLabels[featureName] <- paste0(featureName, refPoint)
        
     }else if(file.exists(featureName)){
        feature <- handle_input(featureName, bedparam, nc=nc)[[1]]
        centerInputs[[featureName]] <- feature
        centerLabels[featureName] <- names(centerFiles[centerFiles==featureName])
     }else{
        stop(paste(featureName, "is not supported!"))
     }
  }

  if(verbose) print("computing coverage for Sample")
  for(queryFile in queryFiles){
    myInput <- queryInputs[[queryFile]]
    libsize <- myInput$size
    queryRegions <- myInput$query
    fileType <- myInput$type
    weight_col <- myInput$weight
    if(verbose) print(paste("size of query regions", libsize))

    queryLabel <- queryLabels[queryFile]
    if(verbose) print(queryLabel)

    for(centerFile in centerFiles){
      centerLabel <- centerLabels[centerFile]
      if(verbose) print(centerLabel)
      centerInput <- centerInputs[[centerFile]]
      centerGr <- centerInput$query
      if(verbose) print("preparing window regions")

      if(refPoint %in% c("center", "start", "end")){
         windowRegions <- resize(centerGr, width = 1, fix = refPoint)
         windowRegions <- promoters(windowRegions, upstream = -ext[1], downstream = ext[2])
         windowRegions <- check_constraints(windowRegions, handleInputParams$genome, queryRegions)
      }else{
        stop("invalid reference point! Must be one of c('center', 'start', 'end')")
      }
      windowRs <- as(split(windowRegions, f=factor(names(windowRegions))), "GRangesList")

      if(verbose) print(paste("number of window regions", length(windowRs)))

      bin_op <- "mean"

      fullMatrix <- parallel_scoreMatrixBin(queryRegions, windowRs, bin_num, bin_op, weight_col, stranded, nc=nc)
      colnames(fullMatrix) <- as.character(colLabel)
      rownames(fullMatrix) <- names(windowRegions)

      scoreMatrix_list[[queryLabel]][[centerLabel]] <- fullMatrix

      if(verbose){
         print(dim(fullMatrix))
         print(length(unique(format_genomic_coordinates(windowRs))))
         #rownames(fullMatrix) <- format_genomic_coordinates(windowRs)
         #write.table(fullMatrix, paste0(queryLabel, "_", centerLabel, "_scoreMatrix.tab"), col.names=NA, sep="\t", quote=FALSE)
      }
    }
  }

  plot_df <- list() # per position, averaged over gene
  stat_df <- list() # per gene, averaged over position
  if(heatmap) heatmap_list <- list()
  Ylab <- ifelse(transform, paste0("Log2 (", Ylab, ")"), Ylab)

  if(verbose) print("collecting coverage data") ## plot multiple bed files on each center

  for(queryFile in queryFiles){

    queryLabel <- queryLabels[queryFile]
    if(verbose) print(queryLabel)
    for(centerFile in centerFiles){
      centerLabel <- centerLabels[centerFile]
      if(verbose) print(centerLabel)

      fullMatrix <- scoreMatrix_list[[queryLabel]][[centerLabel]]
      fullMatrix <- process_scoreMatrix(fullMatrix, scale=scale, rmOutlier=rmOutlier, transform=transform, verbose=verbose)

      colm <- apply(fullMatrix, 2, mean)
      colsd <- apply(fullMatrix, 2, sd)
      colse <- colsd/sqrt(apply(fullMatrix, 2, length))
      collabel <- colLabel
      querybed <- as.factor(rep(queryLabel, length(colm)))
      refbed <- as.factor(rep(centerLabel, length(colm)))


      sub_df <- data.frame(colm, colsd, colse, collabel, querybed, refbed)
      colnames(sub_df) <- c("Intensity", "sd", "se", "Position", "Query", "Reference")

      if(heatmap){
         dataname <- paste(Ylab, queryLabel, centerLabel, sep=":")
         heatmap_list[[dataname]] <- draw_matrix_heatmap(fullMatrix, dataName=dataname, labels_col=collabel,
                                                         levels_col=featureNames, ranges=heatRange, verbose=verbose)
      }

      if(smooth){
        sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df=as.integer(bin_num/5))$y)
        sub_df$se <- as.vector(smooth.spline(sub_df$se, df=as.integer(bin_num/5))$y)
      }

      plot_df[[paste(queryLabel,centerLabel,sep=":")]] <- sub_df

      if(hl[2] > hl[1]){

         xmin <- which(colLabel == hl[1])
         xmax <- which(colLabel == hl[2])
         if(length(xmax)==0) xmax=length(colLabel)
         submatrix <- (fullMatrix[, xmin:xmax])
         submatrix[is.na(submatrix)] <- 0
         Intensity <- as.numeric(rowMeans(submatrix))

         Query <- as.factor(rep(queryLabel, length(Intensity)))
         Reference <- as.factor(rep(centerLabel, length(Intensity)))
         subdf <- data.frame(Intensity, Query, Reference)
         stat_df[[paste(queryLabel,centerLabel,sep=":")]] <- subdf

      }
    }
  }

  mplot_dt <- bind_rows(plot_df) %>%
      mutate(Group=paste0(Query, ":", Reference), .keep="all") %>%
      mutate(lower=Intensity-se, upper=Intensity+se, .keep="all")

  mstat_dt <- NULL
  if(hl[2] > hl[1]){
     mstat_dt <- bind_rows(stat_df) %>%
       mutate(Group=as.factor(paste0(Query, ":", Reference)), .keep="all")
  }

  if(verbose) print("Plotting profile and boxplot")

  plot_list <- list()
  queryLabels <- queryLabels[!queryLabels %in% inputLabels]
  for(i in seq_along(queryLabels)){
    for(beds in combn(queryLabels, i, simplify = FALSE)){
      for(j in seq_along(centerLabels)){
        for(centers in combn(centerLabels, j, simplify = FALSE)){
          if(verbose) print(beds)
          if(verbose) print(centers)

          aplot_df <- mplot_dt %>%
            filter(Query %in% beds & Reference %in% centers)
          aplot_df <- aplot_df %>%  # unify order of factors to get consistent color mapping
             mutate(Query=factor(Query, levels=sort(unique(Query)))) %>%
             mutate(Reference=factor(Reference, levels=sort(unique(Reference))))

          p <- draw_locus_profile(plot_df=aplot_df, cn="Group", sn="Group", Xlab=Xlab, Ylab=Ylab, shade=shade, hl=hl)

          if((i == 1 && j > 1) || (i > 1 && j == 1)){
            if(hl[2] > hl[1]){
               astat_df <- mstat_dt %>%
                 filter(Query %in% beds & Reference %in% centers)
               astat_df <- astat_df %>% # unify order of factors to get consistent color mapping
                  mutate(Query=factor(Query, levels=sort(unique(Query)))) %>%
                  mutate(Reference=factor(Reference, levels=sort(unique(Reference))))

               if(j > 1){
                 p <- draw_locus_profile(plot_df=aplot_df, cn="Reference", sn="Query", Xlab=Xlab, Ylab=Ylab, shade=shade, hl=hl)
                 comp <- combn(seq_along(centers),2, simplify=FALSE)

                 ps1 <- draw_boxplot_logy(stat_df=astat_df, xc="Reference", yc="Intensity",  comp=comp, stats=statsMethod, Ylab=Ylab, logy="linear")
                 ps2 <- draw_boxplot_logy(stat_df=astat_df, xc="Reference", yc="Intensity",  comp=comp, stats=statsMethod, Ylab=Ylab, logy="logy")
                 ps1_wo_outlier <- draw_boxplot_wo_outlier(stat_df=astat_df, xc="Reference", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab)
                 ps1_mean_se <- draw_mean_se_barplot(stat_df=astat_df, xc="Reference", yc="Intensity", comp=comp, Ylab=Ylab)
                 prank <- draw_rank_plot(stat_df=astat_df, xc="Reference", yc="Intensity", Ylab=Ylab, ecdf=TRUE, rank=FALSE)
               }else{
                 p <- draw_locus_profile(plot_df=aplot_df, cn="Query", sn="Reference", Xlab=Xlab, Ylab=Ylab, shade=shade, hl=hl)
                 comp <- combn(seq_along(beds),2, simplify=FALSE)

                 ps1 <- draw_boxplot_logy(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab, logy="linear")
                 ps2 <- draw_boxplot_logy(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab, logy="logy")
                 ps1_wo_outlier <- draw_boxplot_wo_outlier(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab)
                 ps1_mean_se <- draw_mean_se_barplot(stat_df=astat_df, xc="Query", yc="Intensity", Ylab=Ylab, comp=comp)
                 prank <- draw_rank_plot(stat_df=astat_df, xc="Query", yc="Intensity", Ylab=Ylab, ecdf=TRUE, rank=FALSE)
               }

               if(max(i,j)>3){
                  lapply(list(p, ps1_mean_se, prank, ps1, ps2, ps1_wo_outlier), print)
               }else{
                  comp1 <- plot_grid(p, ps1_mean_se, prank, ncol = 3, rel_widths = c(1,1,1))
                  comp2 <- plot_grid(ps1, ps2, ps1_wo_outlier, ncol = 3, align='h', axis='b', rel_widths = c(1,1,1))
                  print(plot_grid(comp1, comp2, ncol=1, rel_heights=c(1,1)))
               }

            }else{
               print(p)
            }
          }else if(i == 1 && j == 1){
             plot_list[[paste(Ylab, beds,centers,sep=":")]] <- p
          }else if(i == length(queryLabels) && j == length(centerLabels)){
            print(p)
          }
        }
      }
    }
  }

  ## plot query profile and heatmap side by side

   rowp <- plot_grid(plotlist=plot_list, nrow = 1, align = 'h', axis= "b")

   if(heatmap){
      qheatmap_list <- heatmap_list[names(plot_list)]
      groblist <- lapply(qheatmap_list, function(x)grid.grabExpr(draw(x, heatmap_legend_side = "top")))
      heatp <- plot_grid(plotlist=groblist, nrow = 1, align = 'h')
      composite <- plot_grid(rowp, heatp, ncol=1, align='v', axis='l', rel_heights=c(1,1))
      print(composite)
   }else{
      print(rowp)
   }


  if(!is.null(inputFiles)){
     plot_list <- list()
     for(i in seq_along(inputLabels)){
        for(beds in combn(inputLabels, i, simplify = FALSE)){
           for(j in seq_along(centerLabels)){
              for(centers in combn(centerLabels, j, simplify = FALSE)){
                 if(verbose) print(beds)
                 if(verbose) print(centers)

                 aplot_df <- mplot_dt %>%
                    filter(Query %in% beds & Reference %in% centers)
                 aplot_df <- aplot_df %>%  # unify order of factors to get consistent color mapping
                    mutate(Query=factor(Query, levels=sort(unique(Query)))) %>%
                    mutate(Reference=factor(Reference, levels=sort(unique(Reference))))

                 p <- draw_locus_profile(plot_df=aplot_df, cn="Group", sn="Group", Xlab=Xlab, Ylab=Ylab, shade=shade, hl=hl)

                 if((i == 1 && j > 1) || (i > 1 && j == 1)){
                    if(hl[2] > hl[1]){
                       astat_df <- mstat_dt %>%
                          filter(Query %in% beds & Reference %in% centers)
                       astat_df <- astat_df %>%
                          mutate(Query=factor(Query, levels=sort(unique(Query)))) %>%
                          mutate(Reference=factor(Reference, levels=sort(unique(Reference))))

                       if(j > 1){
                          p <- draw_locus_profile(plot_df=aplot_df, cn="Reference", sn="Query", Xlab=Xlab, Ylab=Ylab, shade=shade, hl=hl)
                          comp <- combn(seq_along(centers),2, simplify=FALSE)

                          ps1 <- draw_boxplot_logy(stat_df=astat_df, xc="Reference", yc="Intensity",  comp=comp, stats=statsMethod, Ylab=Ylab, logy="linear")
                          ps2 <- draw_boxplot_logy(stat_df=astat_df, xc="Reference", yc="Intensity",  comp=comp, stats=statsMethod, Ylab=Ylab, logy="logy")
                          ps1_wo_outlier <- draw_boxplot_wo_outlier(stat_df=astat_df, xc="Reference", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab)
                          ps1_mean_se <- draw_mean_se_barplot(stat_df=astat_df, xc="Reference", yc="Intensity", comp=comp, Ylab=Ylab)
                          prank <- draw_rank_plot(stat_df=astat_df, xc="Reference", yc="Intensity", Ylab=Ylab, ecdf=TRUE, rank=FALSE)
                       }else{
                          p <- draw_locus_profile(plot_df=aplot_df, cn="Query", sn="Reference", Xlab=Xlab, Ylab=Ylab, shade=shade, hl=hl)
                          comp <- combn(seq_along(beds),2, simplify=FALSE)

                          ps1 <- draw_boxplot_logy(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab, logy="linear")
                          ps2 <- draw_boxplot_logy(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab, logy="logy")
                          ps1_wo_outlier <- draw_boxplot_wo_outlier(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab)
                          ps1_mean_se <- draw_mean_se_barplot(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, Ylab=Ylab)
                          prank <- draw_rank_plot(stat_df=astat_df, xc="Query", yc="Intensity", Ylab=Ylab, ecdf=TRUE, rank=FALSE)
                       }

                       if(max(i,j)>3){
                          lapply(list(p, ps1_mean_se, prank, ps1, ps2, ps1_wo_outlier), print)
                       }else{
                          comp1 <- plot_grid(p, ps1_mean_se, prank, ncol = 3, rel_widths = c(1,1,1))
                          comp2 <- plot_grid(ps1, ps2, ps1_wo_outlier, ncol = 3, align='h', axis='b', rel_widths = c(1,1,1))
                          print(plot_grid(comp1, comp2, ncol=1, rel_heights=c(1,1)))
                       }
                       
                    }else{
                       print(p)
                    }
                 }else if(i == 1 && j == 1){
                    plot_list[[paste(Ylab, beds,centers,sep=":")]] <- p
                 }else if(i == length(inputLabels) && j == length(centerLabels)){
                    print(p)
                 }

              }
           }
        }
     }


     ## plot input profile and heatmap side by side
     rowp <- plot_grid(plotlist=plot_list, nrow = 1, align = 'h', axis= "b")
     #print(rowp)
     if(heatmap){
        iheatmap_list <- heatmap_list[names(plot_list)]
        groblist <- lapply(iheatmap_list, function(x)grid.grabExpr(draw(x, heatmap_legend_side = "top")))
        heatp <- plot_grid(plotlist=groblist, nrow = 1, align = 'h')
        composite <- plot_grid(rowp, heatp, ncol=1, align='v', axis='l', rel_widths=c(1,1))
        print(composite)
     }else{
        print(rowp)
     }


    if(verbose) print("Computing Ratio over input")
    Ylab <- ifelse(transform, "Log2 (Ratio-over-Input)", "Ratio-over-Input")

    inputMatrix_list <- scoreMatrix_list[inputLabels]
    ratiofiles <- queryFiles[!queryFiles %in% inputFiles]
    ratiolabels <- queryLabels[!queryLabels %in% inputLabels]
    ratioMatrix_list <- scoreMatrix_list[ratiolabels]
    for(centerFile in centerFiles){
      centerLabel <- centerLabels[centerFile]
      for(i in seq_along(ratiolabels)){
        fullMatrix <- ratioMatrix_list[[ratiolabels[i]]][[centerLabel]]/inputMatrix_list[[inputLabels[i]]][[centerLabel]]
        colnames(fullMatrix) <- as.character(colLabel)

        ratioMatrix_list[[ratiolabels[i]]][[centerLabel]] <- fullMatrix
      }
    }

    plot_df <- list()
    stat_df <- list()
    if(heatmap) heatmap_list <- list()

    if(verbose) print("collecting ratio data") ## plot multiple bed files on each center

    for(ratiolabel in ratiolabels){

      if(verbose) print(ratiolabel)
      for(centerFile in centerFiles){
        centerLabel <- centerLabels[centerFile]
        if(verbose) print(centerLabel)

        fullMatrix <- ratioMatrix_list[[ratiolabel]][[centerLabel]]
        fullMatrix <- process_scoreMatrix(fullMatrix, scale=FALSE, rmOutlier=rmOutlier, transform=transform, verbose=verbose)

        colm <- apply(fullMatrix, 2, mean)
        colsd <- apply(fullMatrix, 2, sd)
        colse <- colsd/sqrt(apply(fullMatrix, 2, length))
        collabel <- colLabel
        querybed <- as.factor(rep(ratiolabel, length(colm)))
        refbed <- as.factor(rep(centerLabel, length(colm)))


        sub_df <- data.frame(colm, colsd, colse, collabel, querybed, refbed)
        colnames(sub_df) <- c("Intensity", "sd", "se", "Position", "Query", "Reference")

        if(smooth){
          sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df=as.integer(bin_num/5))$y)
          sub_df$se <- as.vector(smooth.spline(sub_df$se, df=as.integer(bin_num/5))$y)
        }

        plot_df[[paste(ratiolabel,centerLabel,sep=":")]] <- sub_df

        if(heatmap){
           dataname <- paste(Ylab, ratiolabel, centerLabel, sep=":")
           heatmap_list[[dataname]] <- draw_matrix_heatmap(fullMatrix, dataName=dataname, labels_col=collabel, levels_col=featureNames, ranges=heatRange, verbose=verbose)
        }

        if(hl[2] > hl[1]){
           xmin <- which(colLabel == hl[1])
           xmax <- which(colLabel == hl[2])
           if(length(xmax)==0) xmax=length(colLabel)
           submatrix <- (fullMatrix[, xmin:xmax])
           submatrix[is.na(submatrix)] <- 0
           Intensity <- as.numeric(rowMeans(submatrix))

           Query <- as.factor(rep(ratiolabel, length(Intensity)))
           Reference <- as.factor(rep(centerLabel, length(Intensity)))
           subdf <- data.frame(Intensity, Query, Reference)
           stat_df[[paste(ratiolabel,centerLabel,sep=":")]] <- subdf
        }
      }
    }

    mplot_dt <- bind_rows(plot_df) %>%
      mutate(Group=paste0(Query, ":", Reference), .keep="all") %>%
      mutate(lower=Intensity-se, upper=Intensity+se, .keep="all")

    mstat_dt <- NULL

    if(hl[2] > hl[1]){
       mstat_dt <- bind_rows(stat_df) %>%
         mutate(Group=as.factor(paste0(Query, ":", Reference)), .keep="all")
    }

    if(verbose) print("Plotting ratio profile and boxplot")
    plot_list <- list()
    for(i in seq_along(ratiolabels)){
      for(beds in combn(ratiolabels, i, simplify = FALSE)){
        for(j in seq_along(centerLabels)){
          for(centers in combn(centerLabels, j, simplify = FALSE)){
            if(verbose) print(beds)
            if(verbose) print(centers)

            aplot_df <- mplot_dt %>%
              filter(Query %in% beds & Reference %in% centers)
            aplot_df <- aplot_df %>%  # unify order of factors to get consistent color mapping
               mutate(Query=factor(Query, levels=sort(unique(Query)))) %>%
               mutate(Reference=factor(Reference, levels=sort(unique(Reference))))

            p <- draw_locus_profile(plot_df=aplot_df, cn="Group", sn="Group", Xlab=Xlab, Ylab=Ylab, shade=shade, hl=hl)

            if((i == 1 && j > 1) || (i > 1 && j == 1)){
               if(hl[2] > hl[1]){
                 astat_df <- mstat_dt %>%
                   filter(Query %in% beds & Reference %in% centers)
                 astat_df <- astat_df %>%
                    mutate(Query=factor(Query, levels=sort(unique(Query)))) %>%
                    mutate(Reference=factor(Reference, levels=sort(unique(Reference))))

                 if(j > 1){
                    p <- draw_locus_profile(plot_df=aplot_df, cn="Reference", sn="Query", Xlab=Xlab, Ylab=Ylab, shade=shade, hl=hl)
                    comp <- combn(seq_along(centers),2, simplify=FALSE)

                    ps1 <- draw_boxplot_logy(stat_df=astat_df, xc="Reference", yc="Intensity",  comp=comp, stats=statsMethod, Ylab=Ylab, logy="linear")
                    ps2 <- draw_boxplot_logy(stat_df=astat_df, xc="Reference", yc="Intensity",  comp=comp, stats=statsMethod, Ylab=Ylab, logy="logy")
                    ps1_wo_outlier <- draw_boxplot_wo_outlier(stat_df=astat_df, xc="Reference", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab)
                    ps1_mean_se <- draw_mean_se_barplot(stat_df=astat_df, xc="Reference", yc="Intensity", comp=comp, Ylab=Ylab)
                    prank <- draw_rank_plot(stat_df=astat_df, xc="Reference", yc="Intensity", Ylab=Ylab, ecdf=TRUE, rank=FALSE)
                 }else{
                    p <- draw_locus_profile(plot_df=aplot_df, cn="Query", sn="Reference", Xlab=Xlab, Ylab=Ylab, shade=shade, hl=hl)
                    comp <- combn(seq_along(beds),2, simplify=FALSE)

                    ps1 <- draw_boxplot_logy(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab, logy="linear")
                    ps2 <- draw_boxplot_logy(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab, logy="logy")
                    ps1_wo_outlier <- draw_boxplot_wo_outlier(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab)
                    ps1_mean_se <- draw_mean_se_barplot(stat_df=astat_df, xc="Query", yc="Intensity", comp=comp, Ylab=Ylab)
                    prank <- draw_rank_plot(stat_df=astat_df, xc="Query", yc="Intensity", Ylab=Ylab, ecdf=TRUE, rank=FALSE)
                 }

                 if(max(i,j)>3){
                    lapply(list(p, ps1_mean_se, prank, ps1, ps2, ps1_wo_outlier), print)
                 }else{
                    comp1 <- plot_grid(p, ps1_mean_se, prank, ncol = 3, rel_widths = c(1,1,1))
                    comp2 <- plot_grid(ps1, ps2, ps1_wo_outlier, ncol = 3, align='h', axis='b', rel_widths = c(1,1,1))
                    print(plot_grid(comp1, comp2, ncol=1, rel_heights=c(1,1)))
                 }
                 
               }else{
                  print(p)
               }
            }else if(i == 1 && j == 1){
               plot_list[[paste(Ylab, beds,centers,sep=":")]] <- p
            }else if(i == length(ratiolabels) && j == length(centerLabels)){
               print(p)
            }
          }
        }
      }
    }

    ## plot ratio profile and heatmap side by side
    if(heatmap){
       rheatmap_list <- heatmap_list[names(plot_list)]
       groblist <- lapply(rheatmap_list, function(x)grid.grabExpr(draw(x, heatmap_legend_side = "top")))
       heatp <- plot_grid(plotlist=groblist, nrow = 1, align = 'h')
       #print(heatp)
    }
    rowp <- plot_grid(plotlist=plot_list, nrow = 1, align = 'h', axis='b')
    #print(rowp)
    if(heatmap){
       composite <- plot_grid(rowp, heatp, ncol=1, align='v', axis='l', rel_widths=c(1,1))
       print(composite)
    }else{
       print(rowp)
    }
  }

  if(!is.null(outPrefix)){
    on.exit(dev.off(), add=TRUE)
  }

  invisible(list("plot"=mplot_dt, "stat"=mstat_dt))
}

#' @title Plot signal around custom genomic loci and random loci for comparison
#
#' @description Plot reads or peak signal intensity of samples in the query files around reference locus defined in the centerFiles. The upstream and downstream windows flanking genes
#' can be given separately, a smaller window can be defined to allow statistical comparisons between reference and random loci.
#'
#' @param queryFiles a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param centerFiles a vector of reference file names. The file should be .bed format only
#' @param inputFiles a vector of input sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param handleInputParams a list of parameters for \code{handle_input}
#' @param ext a vector of two integers defining upstream and downstream boundaries of the plot window, flanking the reference locus
#' @param hl a vector of two integers defining upstream and downstream boundaries of the highlight window, flanking the reference locus
#' @param stranded logical, indicating whether the strand of the feature should be considered
#' @param scale logical, indicating whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smooth logical, indicating whether the line should smoothed with a spline smoothing algorithm
#' @param rmOutlier logical, indicating whether a row with abnormally high values in the score matrix should be removed
#' @param outPrefix a string specifying output file prefix for plots (outPrefix.pdf)
#' @param refPoint a string in c("start", "center", "end")
#' @param Xlab a string denotes the label on x-axis
#' @param shade logical indicating whether to place a shaded rectangle around the point of interest
#' @param binSize an integer defines bin size for intensity calculation
#' @param transform logical, whether to log2 transform the data matrix
#' @param n_random an integer denotes the number of randomization should be formed
#' @param statsMethod a string in c("wilcox.test", "t.test"), for pair-wise groups comparisons
#' @param verbose logical, indicating whether to output additional information (data used for plotting or statistical test results)
#' @param nc integer, number of cores for parallel processing
#'
#' @return a dataframe containing the data used for plotting
#' @author Shuye Pu
#'
#' @examples
#' queryFiles <- system.file("data", "test_clip.bam", package="GenomicPlotData")
#' names(queryFiles) <- "query"
#' inputFiles <- system.file("data", "test_clip_input.bam", package="GenomicPlotData")
#' names(inputFiles) <- "input"
#'
#' @export plot_reference_locus_with_random

plot_reference_locus_with_random <- function(queryFiles, centerFiles, txdb, ext=c(0,0), hl=c(0,0), shade=FALSE,
                                             handleInputParams=NULL, verbose=FALSE, smooth=FALSE, transform=FALSE, 
                                             binSize=10, refPoint="center", Xlab="Center", Ylab="Signal intensity",
                                             inputFiles=NULL, stranded=TRUE, scale=FALSE, outPrefix=NULL, 
                                             rmOutlier=FALSE, n_random=1, statsMethod="wilcox.test", nc=2){

  queryLabels <- names(queryFiles)
  names(queryLabels) <- queryFiles
  inputLabels <- names(inputFiles)
  names(inputLabels) <- inputFiles
  centerLabels <- names(centerFiles)
  names(centerLabels) <- centerFiles

  if(!is.null(inputFiles)){
    if(length(queryFiles) != length(inputFiles)){
      stop("the number of queryFiles and inputFiles must be the same!")
    }
    names(inputLabels) <- inputFiles
    queryFiles <- c(queryFiles, inputFiles)
    queryLabels <- c(queryLabels, inputLabels)
  }

  if(!is.null(outPrefix)){
     while(!is.null(dev.list())){
        dev.off()
     }
     pdf(paste(outPrefix, "pdf", sep="."), height=8, width=12)
  }

  ext[2] <- ext[2] - (ext[2]-ext[1])%%binSize ## to avoid binSize inconsistency, as the final binSize is dictated by bin_num
  colLabel <- seq(ext[1], (ext[2]-binSize), binSize)
  bin_num <- round((ext[2]-ext[1])/binSize)
  bin_op <- "mean"

  ## ranges for overlap count
  rb <- ext[2] #broad
  rn <- hl[2]  #narrow

  ## get protein-coding genes features

  if(verbose) print("Collecting protein_coding gene features")
  #exons <- get_genomic_feature_coordinates(txdb, "exon", longest=TRUE)
  utr5 <- get_genomic_feature_coordinates(txdb, "utr5", longest=TRUE)
  utr3 <- get_genomic_feature_coordinates(txdb, "utr3", longest=TRUE)
  cds <- get_genomic_feature_coordinates(txdb, "cds", longest=TRUE)
  #gene <- get_genomic_feature_coordinates(txdb, "transcript", longest=TRUE)


  region_list <- list(#"Transcript" = exons$GRanges,
                      "5'UTR" = utr5$GRanges,
                      "CDS" = cds$GRanges,
                      "3'UTR" = utr3$GRanges,
                      #"Gene" = gene$GRanges,
                      "unrestricted" = NULL)

  if(verbose) print(lapply(region_list, length))

  if(verbose) print("computing coverage for Sample")
  scoreMatrix_list <- list()
  quantile_list <- list()
  scoreMatrix_list_random <- list()
  quantile_list_random <- list()

  queryInputs <- handle_input(queryFiles, handleInputParams, nc=nc)

  bedparam <- handleInputParams
  bedparam$CLIP_reads <- FALSE
  bedparam$fix_width <- 0
  bedparam$useScore <- FALSE
  bedparam$outRle <- FALSE
  bedparam$useSizeFactor <- FALSE
  centerInputs <- handle_input(centerFiles, bedparam)

  for(queryFile in queryFiles){
    myInput <- queryInputs[[queryFile]]
    libsize <- myInput$size
    queryRegions <- myInput$query
    fileType <- myInput$type
    weight_col <- myInput$weight

    queryLabel <- queryLabels[queryFile]
    if(verbose) print(queryLabel)
    if(verbose) print(libsize)

    for(centerFile in centerFiles){
      centerLabel <- centerLabels[centerFile]
      if(verbose) print(centerLabel)
      centerInput <- centerInputs[[centerFile]]

      windowRegionsALL <- centerInput$query

      for(regionName in names(region_list)){
        if(verbose) print(paste("Processing genomic region", regionName))

        if(regionName == "unrestricted"){
          windowRegions <- windowRegionsALL
        }else{
          region <- region_list[[regionName]]
          if(any(unique(as.vector(strand(windowRegionsALL))) %in% c("*", ".", ""))){
             windowRegions <- plyranges::filter_by_overlaps(windowRegionsALL, region)
             if(verbose) print("the center file is Unstranded")
          }else{
             windowRegions <- filter_by_overlaps_stranded(windowRegionsALL, region)
             if(verbose) print("the center file is stranded")
          }
        }

        windowRs <- resize(windowRegions, width = 1, fix = refPoint)
        windowRs <- promoters(windowRs, upstream = -ext[1], downstream = ext[2])
        windowRs <- as(split(windowRs, f=factor(seq(1:length(windowRs)))), "GRangesList")

        fullMatrix <- parallel_scoreMatrixBin(queryRegions, windowRs, bin_num, bin_op, weight_col, stranded, nc=nc)
        fullMatrix <- process_scoreMatrix(fullMatrix, scale=scale, rmOutlier, transform=transform, verbose=verbose)

        scoreMatrix_list[[queryLabel]][[centerLabel]][[regionName]] <- fullMatrix


        ## create randomized centers, repeat n_random times

        windowRs <- random_results <- lapply(seq(1, n_random), function(i){
          random_points <- sample(ext[1]:ext[2], length(windowRegions), replace=TRUE)
          rwindowRegions <- shift(windowRegions, shift=random_points, use.names=TRUE)

          windowR <- resize(rwindowRegions, width = 1, fix = refPoint)
          windowR <- promoters(windowR, upstream = -ext[1], downstream = ext[2])
          windowR <- as(split(windowR, f=factor(seq(1:length(windowR)))), "GRangesList")
          invisible(windowR)
        })

        random_matricies <- lapply(windowRs, function(x){
          parallel_scoreMatrixBin(queryRegions, x, bin_num, bin_op, weight_col, stranded, nc=nc)
        })

        ## unify the dimensions of random_matrices
        dims <- vapply(random_matricies, dim, numeric(2))
        minrows <- min(dims[1,])
        random_matricies <- lapply(random_matricies, function(x)x[1:minrows,])

        a3d <- array(unlist(random_matricies), c(dim(random_matricies[[1]]), length(random_matricies))) ## turn the list of matrices into a 3-d array
        fullMatrix <- apply(a3d, 1:2, mean)  # take average of all matrices

        fullMatrix <- process_scoreMatrix(fullMatrix, scale=scale, rmOutlier=rmOutlier, transform=transform, verbose=verbose)

        scoreMatrix_list_random[[queryLabel]][[centerLabel]][[regionName]] <- fullMatrix

      }
    }
  }


  ## plot reference center and random center for each bed

  Ylab <- ifelse(transform, paste0("Log2 (", Ylab, ")"), Ylab)
  for(queryFile in queryFiles){

    if(verbose) print(paste("Processing query", queryFile))
    queryLabel <- queryLabels[queryFile]

    for(centerFile in centerFiles){
      if(verbose) print(paste("Processing reference", centerFile))
      centerLabel <- centerLabels[centerFile]

      for(regionName in names(region_list)){
        if(verbose) print(paste("Processing genomic region", regionName))
        plot_df <- NULL
        stat_df <- NULL
        countOverlap_df <- NULL

        fullMatrix_list <- list(scoreMatrix_list[[queryLabel]][[centerLabel]][[regionName]],
                                scoreMatrix_list_random[[queryLabel]][[centerLabel]][[regionName]])
        names(fullMatrix_list) <- c(centerLabel, "Random")

        refsize <- nrow(fullMatrix_list[[centerLabel]])

        for(alabel in c(centerLabel, "Random")){
          fullMatrix <- fullMatrix_list[[alabel]]
          colm <- apply(fullMatrix, 2, mean)
          colsd <- apply(fullMatrix, 2, sd)
          colse <- colsd/sqrt(apply(fullMatrix, 2, length))
          collabel <- colLabel
          querybed <- as.factor(rep(queryLabel, length(colm)))
          refbed <- as.factor(rep(alabel, length(colm)))

          sub_df <- data.frame(colm, colsd, colse, collabel, querybed, refbed)
          colnames(sub_df) <- c("Intensity", "sd", "se", "Position", "Query", "Reference")

          if(smooth){
            sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df=as.integer(bin_num/5))$y)
            sub_df$se <- as.vector(smooth.spline(sub_df$se, df=as.integer(bin_num/5))$y)
          }

          sub_df <- mutate(sub_df, lower=Intensity-se, upper=Intensity+se)
          plot_df <- rbind(plot_df, sub_df)

          if(hl[2] > hl[1]){
             xmin <- which(colLabel == hl[1])
             xmax <- which(colLabel == hl[2])
             if(length(xmax)==0) xmax=length(colLabel)
             submatrix <- (fullMatrix[, xmin:xmax])
             submatrix[is.na(submatrix)] <- 0
             Intensity <- as.numeric(rowMeans(submatrix))

             Query <- as.factor(rep(queryLabel, length(Intensity)))
             Reference <- as.factor(rep(alabel, length(Intensity)))
             subdf <- data.frame(Intensity, Query, Reference)

             stat_df <- rbind(stat_df, subdf)
          }
        }

        plot_df <- plot_df %>%  # unify order of factors to get consistent color mapping
           mutate(Query=factor(Query, levels=sort(unique(Query)))) %>%
           mutate(Reference=factor(Reference, levels=sort(unique(Reference))))

        p <- draw_locus_profile(plot_df, cn="Reference", sn="Query", Xlab=Xlab, Ylab=Ylab, shade=shade, hl=hl) +
          ggtitle(paste("Feature:", regionName, "\nReference size:", refsize, "\nSample name:", queryLabel))

        if(hl[2] > hl[1]){
           stat_df <- stat_df %>%  # unify order of factors to get consistent color mapping
              mutate(Query=factor(Query, levels=sort(unique(Query)))) %>%
              mutate(Reference=factor(Reference, levels=sort(unique(Reference))))
           comp <- list(c(1, 2))

           ps1 <- draw_boxplot_logy(stat_df=stat_df, xc="Reference", yc="Intensity",  comp=comp, stats=statsMethod, Ylab=Ylab, logy="linear")
           ps2 <- draw_boxplot_logy(stat_df=stat_df, xc="Reference", yc="Intensity",  comp=comp, stats=statsMethod, Ylab=Ylab, logy="logy")
           ps1_wo_outlier <- draw_boxplot_wo_outlier(stat_df=stat_df, xc="Reference", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab)
           ps1_mean_se <- draw_mean_se_barplot(stat_df=stat_df, xc="Reference", yc="Intensity", comp=comp, Ylab=Ylab)
           prank <- draw_rank_plot(stat_df=stat_df, xc="Reference", yc="Intensity", Ylab=Ylab)

           comp1 <- plot_grid(p, ps1_mean_se, prank, ncol = 3, rel_widths = c(1,1,1))
           comp2 <- plot_grid(ps1, ps2, ps1_wo_outlier, ncol = 3, rel_widths = c(1,1,1))
           print(plot_grid(comp1, comp2, ncol=1, rel_heights=c(1,1)))
        }else{
           print(p)
        }
      }
    }
  }

  if(!is.null(inputFiles)){
    Ylab <- ifelse(transform, "Log2 (Ratio-over-Input)", "Ratio-over-Input")

    ratiofiles <- queryFiles[!queryFiles %in% inputFiles]
    ratiolabels <- queryLabels[!queryLabels %in% inputLabels]
    inputMatrix_list <- scoreMatrix_list[inputLabels]
    ratioMatrix_list <- scoreMatrix_list[ratiolabels]

    inputMatrix_list_random <- scoreMatrix_list_random[inputLabels]
    ratioMatrix_list_random <- scoreMatrix_list_random[ratiolabels]

    for(centerFile in centerFiles){
      centerLabel <- centerLabels[centerFile]
      for(regionName in names(region_list)){
        for(i in seq_along(ratiolabels)){
          rm <- ratioMatrix_list[[ratiolabels[i]]][[centerLabel]][[regionName]]
          im <- inputMatrix_list[[inputLabels[i]]][[centerLabel]][[regionName]]

          fullMatrix <- rm/im
          fullMatrix <- process_scoreMatrix(fullMatrix, scale=FALSE, rmOutlier=rmOutlier, transform=transform, verbose=verbose)

          ratioMatrix_list[[ratiolabels[i]]][[centerLabel]][[regionName]] <- fullMatrix

          ## for random centers
          rmr <- ratioMatrix_list_random[[ratiolabels[i]]][[centerLabel]][[regionName]]
          imr <- inputMatrix_list_random[[inputLabels[i]]][[centerLabel]][[regionName]]
          minrowr <- min(nrow(rmr), nrow(imr))

          fullMatrix <- rmr[1:minrowr,]/imr[1:minrowr,]
          fullMatrix <- process_scoreMatrix(fullMatrix, scale=FALSE, rmOutlier=rmOutlier, transform=transform, verbose=verbose)

          ratioMatrix_list_random[[ratiolabels[i]]][[centerLabel]][[regionName]] <- fullMatrix
        }
      }
    }

    for(ratiofile in ratiofiles){
      if(verbose) print(paste("Processing ratio for query", ratiofile))
      ratiolabel <- ratiolabels[ratiofile]

      for(centerFile in centerFiles){
        if(verbose) print(paste("Processing ratio for reference", centerFile))
        centerLabel <- centerLabels[centerFile]

        for(regionName in names(region_list)){
          if(verbose) print(paste("Processing ratio for genomic region", regionName))
          plot_df <- NULL
          stat_df <- NULL
          countOverlap_df <- NULL
          refsize <- NULL

          fullMatrix_list <- list(ratioMatrix_list[[ratiolabel]][[centerLabel]][[regionName]],
                                  ratioMatrix_list_random[[ratiolabel]][[centerLabel]][[regionName]])
          names(fullMatrix_list) <- c(centerLabel, "Random")

          refsize <- nrow(fullMatrix_list[[centerLabel]])

          for(alabel in c(centerLabel, "Random")){

            fullMatrix <- fullMatrix_list[[alabel]]

            colm <- apply(fullMatrix, 2, mean)
            colsd <- apply(fullMatrix, 2, sd)
            colse <- colsd/sqrt(apply(fullMatrix, 2, length))
            collabel <- colLabel
            querybed <- as.factor(rep(queryLabel, length(colm)))
            refbed <- as.factor(rep(alabel, length(colm)))

            sub_df <- data.frame(colm, colsd, colse, collabel, querybed, refbed)
            colnames(sub_df) <- c("Intensity", "sd", "se", "Position", "Query", "Reference")

            if(smooth){
              sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df=as.integer(bin_num/5))$y)
              sub_df$se <- as.vector(smooth.spline(sub_df$se, df=as.integer(bin_num/5))$y)
            }

            sub_df <- mutate(sub_df, lower=Intensity-se, upper=Intensity+se)
            plot_df <- rbind(plot_df, sub_df)

            if(hl[2] > hl[1]){
               xmin <- which(colLabel == hl[1])
               xmax <- which(colLabel == hl[2])
               if(length(xmax)==0) xmax=length(colLabel)
               submatrix <- (fullMatrix[, xmin:xmax])
               submatrix[is.na(submatrix)] <- 0
               Intensity <- as.numeric(rowMeans(submatrix))
               Query <- as.factor(rep(queryLabel, length(Intensity)))
               Reference <- as.factor(rep(alabel, length(Intensity)))
               subdf <- data.frame(Intensity, Query, Reference)

               stat_df <- rbind(stat_df, subdf)
            }
          }

          plot_df <- plot_df %>%  # unify order of factors to get consistent color mapping
             mutate(Query=factor(Query, levels=sort(unique(Query)))) %>%
             mutate(Reference=factor(Reference, levels=sort(unique(Reference))))

          p <- draw_locus_profile(plot_df, cn="Reference", sn="Query", Xlab=Xlab, Ylab=Ylab, shade=shade, hl=hl) +
            ggtitle(paste("Feature:", regionName, "\nReference size:", refsize, "\nSample name:", ratiolabel))

          if(hl[2] > hl[1]){
             stat_df <- stat_df %>%  # unify order of factors to get consistent color mapping
                mutate(Query=factor(Query, levels=sort(unique(Query)))) %>%
                mutate(Reference=factor(Reference, levels=sort(unique(Reference))))
             comp <- list(c(1, 2))

             ps1 <- draw_boxplot_logy(stat_df=stat_df, xc="Reference", yc="Intensity",  comp=comp, stats=statsMethod, Ylab=Ylab, logy="linear") 
             ps2 <- draw_boxplot_logy(stat_df=stat_df, xc="Reference", yc="Intensity",  comp=comp, stats=statsMethod, Ylab=Ylab, logy="logy") 
             ps1_wo_outlier <- draw_boxplot_wo_outlier(stat_df=stat_df, xc="Reference", yc="Intensity", comp=comp, stats=statsMethod, Ylab=Ylab)
             ps1_mean_se <- draw_mean_se_barplot(stat_df=stat_df, xc="Reference", yc="Intensity", comp=comp, Ylab=Ylab)
             prank <- draw_rank_plot(stat_df=stat_df, xc="Reference", yc="Intensity", Ylab=Ylab)

             comp1 <- plot_grid(p, ps1_mean_se, prank, ncol = 3, rel_widths = c(1,1,1))
             comp2 <- plot_grid(ps1, ps2, ps1_wo_outlier, ncol = 3, rel_widths = c(1,1,1))
             print(plot_grid(comp1, comp2, ncol=1, rel_heights=c(1,1)))

          }else{
             print(p)
          }
        }
      }
    }
  }

  if(!is.null(outPrefix)) on.exit(dev.off(), add=TRUE)
}



#' @title Plot correlation of bam files
#'
#' @description plot correlation in reads coverage distributions along the genome for bam files
#'
#' @param bamfiles a named vector of strings denoting file names
#' @param binSize an integer denoting the tile width for tiling the genome, default 1000000
#' @param outPrefix a string denoting output file name in pdf format
#' @param handleInputParams a list of parameters for \code{handle_input}
#' @param nc integer, number of cores for parallel processing
#'
#' @return NULL
#'
#' @examples
#' queryFiles <- c(system.file("data", "test_chip_chr19.bam", package="GenomicPlotData"),
#' system.file("data", "test_clip_chr19.bam", package="GenomicPlotData"))
#' names(queryFiles) <- c("chip_query", "clip_query")
#' inputFiles <- c(system.file("data", "test_chip_input_chr19.bam", package="GenomicPlotData"),
#'                 system.file("data", "test_clip_input_chr19.bam", package="GenomicPlotData"))
#' names(inputFiles) <- c("chip_input", "clip_input")
#'
#' op=NULL
#' handleInputParams <- list(CLIP_reads=FALSE, fix_width=0, fix_point="start", norm=FALSE, useScore=FALSE,
#'                           outRle=FALSE, useSizeFactor=FALSE, genome="hg19")
#' plot_bam_correlation(bamfiles=c(queryFiles[1], inputFiles[1]), binSize=10000, outPrefix=op,
#'                      handleInputParams=handleInputParams, nc=2)
#'
#'
#' @export plot_bam_correlation
#'
plot_bam_correlation <- function(bamfiles, binSize=1e6, outPrefix=NULL, handleInputParams=NULL, verbose=FALSE, nc=2){

  bamlabels <- names(bamfiles)
  handleInputParams$outRle <- FALSE # force query to be GRanges
  if(verbose) print("Computing bam correlation")
  outlist <- handle_input(inputFiles=bamfiles, handleInputParams, nc=nc)

  seqi <- Seqinfo(genome=handleInputParams$genome)

  tileBins <-  tileGenome(seqi, tilewidth=binSize, cut.last.tile.in.chrom=TRUE)

  grange_list <- lapply(outlist, function(x)x$query)

  score_list <- parallel_countOverlaps(grange_list, tileBins, nc=nc)

  bins_df <- data.frame(chr=as.vector(seqnames(tileBins)), start=start(tileBins),end=end(tileBins),strand=strand(tileBins))
  bins <- do.call(paste, c(bins_df, sep="_"))

  count_mat <- data.matrix(bind_cols(score_list))
  rownames(count_mat) <- bins
  count_mat[is.na(count_mat)] <- 0
  count_mat <- count_mat[apply(count_mat, 1, sum)>0,]

  norm_factor <- vapply(outlist, function(x) x$size/1e6, numeric(1))

  df <- as.data.frame(t(t(count_mat)/norm_factor)) ## convert to counts per million (CPM)
  colnames(df) <- bamlabels

  long_df <- pivot_longer(df, cols=colnames(df), names_to="Sample", values_to="Count") %>%
    mutate(Sample=as.factor(Sample)) %>%
     group_by(Sample) %>%
     arrange(Count) %>%
     mutate(cumCount=cumsum(Count)) %>%
     mutate(Rank=order(cumCount)) %>%
     mutate(Fraction=cumCount/max(cumCount), Rank=Rank/max(Rank))

  if(!is.null(outPrefix)) pdf(paste0(outPrefix, ".pdf"), width=10, height=8)

  p1 <- ggplot(data=long_df, aes(x=Rank, y=Fraction, color=Sample)) +
     geom_line() +
     ggtitle(paste("Binned read counts distribution: bin size =", binSize)) +
     labs(x="Rank(Count)", y=paste("Fraction over highest coverage"))
  print(p1)

  p2 <- ggplot(data=long_df, aes(x=Count, color=Sample)) +
    stat_ecdf() +
    ggtitle(paste("Binned read counts distribution: bin size =", binSize)) +
    labs(x=expression(paste(log[2], " (Count)")), y=paste("Pencentage"))
  #print(p2)


  ## code from pairs example of base R
  ## put (absolute) correlations on the upper panels,
  ## with size proportional to the correlations.
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
  {
    #usr <- par("usr"); on.exit(par(usr))
    old=par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
    on.exit(par(old))
  }

  ## put histograms on the diagonal
  panel.hist <- function(x, ...)
  {
    usr <- par("usr")
    old=par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
    on.exit(par(old))
  }

  ## END code from pairs example

  mat <- cor(log2(df+1))
  mat_long <- pivot_longer(as.data.frame(mat), cols=1:ncol(mat), names_to="X", values_to="correlation") %>%
     mutate(Y=rep(rownames(mat), each=ncol(mat)))
  g <- ggplot(mat_long, aes(X, Y)) +
     geom_tile(aes(fill = correlation)) +
     geom_text(aes(label = round(correlation, digits=2), color="white")) +
     scale_fill_viridis(discrete=FALSE) +
     theme_minimal() +
     theme(legend.position="none",
           axis.title=element_blank(),
           panel.grid.major = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.ticks = element_blank())
  print(g)
  grid.newpage()
  #pheatmap(cor(log2(df+1)), display_numbers = TRUE)
  if(length(bamfiles)<6) pairs(log2(df+1), lower.panel = panel.smooth, upper.panel=panel.cor, diag.panel = panel.hist, main=paste("log2(CPM/bin), bin size =", binSize))

  if(!is.null(outPrefix)) on.exit(dev.off(), add=TRUE)

  invisible(df)
}

#' @title Annotate peaks with genomic features and genes
#'
#' @description  Produce a table of transcripts targeted by peaks, and generate plots for target gene types, and peak distribution in genomic features
#'
#' @param peakFile a string denoting the peak file name, only .bed format is allowed
#' @param gtfFile path to a gene annotation gtf file with gene_biotype field
#' @param handleInputParams a list of parameters for \code{handle_input}
#' @param fiveP extension out of the 5' boundary of genes for defining promoter: -fiveP TSS + dsTSS
#' @param dsTSS extension downstream of TSS for defining promoter: -fiveP TSS + dsTSS
#' @param threeP extension out of the 3' boundary of genes for defining termination region: -0 TTS + threeP
#' @param outPrefix a string denoting output file name in pdf format
#' @param simple logical, indicating whether 5'UTR and 3'UTR are annotated in the gtffile
#' @param verbose, logical, to indicate whether to write the annotation results to a file
#'
#' @return a list of two dataframes, 'annotation' is the annotation per peak, 'stat' is the summary stats for pie chart
#' @author Shuye Pu
#' @examples
#' library(BiocFileCache)
#' bfc <- BiocFileCache()
#' url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
#' gz <- bfcrpath(bfc, url)
#' gtffile <- gsub(".gz", "", gz)
#' if(!file.exists(gtffile)) gtffile <- R.utils::gunzip(gz, remove=FALSE)
#'
#' centerFile <- system.file("data", "test_chip_peak.narrowPeak", package="GenomicPlotData")
#' names(centerFile) <- c("narrowPeak")
#' op <- "test_plot_peak_annotation"
#' handleBedparams <- list(fix_width=0, fix_point="center", useScore=FALSE, outRle=FALSE, CLIP_reads=FALSE,
#'                         norm=FALSE, useSizeFactor=FALSE, genome="hg19")
#'
#' plot_peak_annotation(peakFile=centerFile, gtfFile=gtffile, handleInputParams=handleBedparams, fiveP=0,
#'                      threeP=0, simple=FALSE)
#'
#' @export plot_peak_annotation
#'
plot_peak_annotation <- function(peakFile, gtfFile, handleInputParams=NULL, fiveP=1000, dsTSS=300, threeP=1000, simple=FALSE, outPrefix=NULL, verbose=FALSE){

  peakLabel <- names(peakFile)
  handleInputParams$useScore=FALSE
  handleInputParams$outRle=FALSE
  bedin <- handle_input(inputFiles=peakFile, handleInputParams)
  peak <- bedin[[peakFile]]$query
  
  stranded <- TRUE
  if(any(runValue(strand(peak)) %in% c("*", ".", " ", ""))){
     stranded <- FALSE
  }
  
  suppressWarnings(suppressMessages(txdb <- makeTxDbFromGFF(gtfFile)))

  if(simple){
      ## the 5'UTR and 3'UTR are not annotted
    exon <- get_genomic_feature_coordinates(txdb, "exon", longest=FALSE)
    intron <- get_genomic_feature_coordinates(txdb, "intron", longest=FALSE)
    gene <- get_genomic_feature_coordinates(txdb, "gene", longest=FALSE)$GRanges

    promoter <- promoters(gene, upstream=fiveP, downstream=300, use.names=FALSE)
    names(promoter) <- NULL
    TTS <- flank(gene, width=threeP, both=FALSE, start=FALSE, ignore.strand=FALSE)
    names(TTS) <- NULL

    pg <- mergeByOverlaps(peak, gene)
    targeted_gene <- cbind(gr2df(pg$peak), gr2df(pg$gene))

    pp <- mergeByOverlaps(peak, promoter)
    targeted_promoter <- cbind(gr2df(pp$peak), gr2df(pp$promoter))
    

    if(verbose){
       write.table(targeted_gene, paste(peakLabel, "_targeted_annotated_gene.tab", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
       write.table(targeted_promoter, paste(peakLabel, "_targeted_promoter.tab", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
    }

    features <- GRangesList("Promoter"=promoter,
                            "TTS"=TTS,
                            "Exon"=unlist(exon$GRangesList, use.names=FALSE),
                            "Intron"=unlist(intron$GRangesList, use.names=FALSE),
                            compress=FALSE)

    if(!is.null(outPrefix)) pdf(paste0(outPrefix, ".pdf"), height=8, width=12)

    suppressWarnings(annot <- annotateWithFeatures(peak, features, strand.aware=TRUE, intersect.chr=FALSE))
    if(verbose) print(annot)
    precedence_count <- annot@num.precedence
    lengths <- vapply(features, FUN=function(x)sum(width(x)), FUN.VALUE=numeric(1))

    df <- data.frame(count=precedence_count, len=lengths) %>%
      mutate(percent = count / sum(count)) %>%
      mutate(feature = rownames(.)) %>%
      mutate(labels = percent(percent, accuracy = 0.1)) %>%
      mutate(percent = round(percent*100, digits=1))%>%
      mutate(norm_count=round(count*1000/len, digits=3)) %>%
      mutate(norm_percent = norm_count / sum(norm_count)) %>%
      mutate(norm_labels = percent(norm_percent, accuracy = 0.1)) %>%
      mutate(norm_percent = round(norm_percent*100, digits=1)) %>%
      mutate(feature = forcats::fct_inorder(feature))

    df2 <- df %>%
       mutate(csum = rev(cumsum(rev(percent))),
              pos = percent/2 + lead(csum, 1),
              pos = ifelse(is.na(pos), percent/2, pos),
              norm_csum = rev(cumsum(rev(norm_percent))),
              norm_pos = norm_percent/2 + lead(norm_csum, 1),
              norm_pos = ifelse(is.na(norm_pos), norm_percent/2, norm_pos))
    if(verbose) print(df)

    ap1 <- ggplot(df, aes(x = "", y = percent, fill = feature)) +
       geom_col(color = "white") +
       scale_y_continuous(breaks = df2$pos, labels = df2$labels) +
       guides(fill = guide_legend(title = "Feature", nrow=2)) +
       scale_fill_viridis(discrete=TRUE) +
       coord_polar(theta = "y") +
       ggtitle("Absolute count") +
       theme(axis.ticks = element_blank(),
             axis.title = element_blank(),
             axis.text = element_text(size = 15),
             legend.position = "top",
             panel.background = element_rect(fill = "white"))
    
    ap2 <- ggplot(df, aes(x = "", y = norm_percent, fill = feature)) +
       geom_col(color = "white") +
       coord_polar(theta = "y") +
       scale_y_continuous(breaks = df2$norm_pos, labels = df2$norm_labels) +
       guides(fill = guide_legend(title = "Feature", nrow=2)) +
       scale_fill_viridis(discrete=TRUE) +
       ggtitle("Length-normalized count") +
       theme(axis.ticks = element_blank(),
             axis.title = element_blank(),
             axis.text = element_text(size = 10),
             legend.position = "top",
             panel.background = element_rect(fill = "white"))

    print(plot_grid(ap1, ap2))
    if(!is.null(outPrefix)) on.exit(dev.off(), add=TRUE)
    invisible(list(annotation=NULL, stat=df, simplified=NULL))
  }else{

    if(verbose) print("Collecting gene info")

    suppressWarnings(suppressMessages({gff <- RCAS::importGtf(saveObjectAsRds = TRUE, filePath = gtfFile);
                     overlaps <- gr2df(RCAS::queryGff(queryRegions=peak, gffData=gff))}))
    geneType <- NULL
    if("gene_type" %in% colnames(overlaps)){
       geneType <- "gene_type"
    }else if("gene_biotype" %in% colnames(overlaps)){
       geneType <- "gene_biotype"
    }
    if(is.null(geneType)){
       stop("The annotation gtf file must have a field 'gene_type' or 'gene_biotype'")
    }
    gene_info_table <- gr2df(gff) %>% 
       filter(type == "transcript") %>%
       select(chr, start, end, name, strand, transcript_id, gene_id, gene_name, all_of(geneType))

    # To find out the distribution of the query regions across gene types:
    if(verbose) print("Computing barchart of gene types")

    df <- overlaps %>%
       filter(type == "gene") %>%
       group_by(.data[[geneType]]) %>%
       summarize(count=dplyr::n_distinct(queryIndex)) %>%
       rename_with(~ c("feature", "count"))

    intergenic_count <- length(peak) - sum(df$count)
    if(intergenic_count < 0) intergenic_count <- 0
    intergenic <- data.frame("feature"= "intergenic", "count"= intergenic_count)

    selected_feature <- c("protein_coding", "lincRNA", "antisense", "pseudogene", "snRNA", "snoRNA", "rRNA")
    selected_count <- rep(0, length(selected_feature))
    selected <- data.frame("feature"=selected_feature, "count"=selected_count)

    for(ft in selected_feature){
       if(ft %in% df$feature){
          selected[selected$feature == ft, "count"] <- df[df$feature == ft, "count"]
       }
    }

    other_df <- filter(df, !feature %in% selected_feature)
    other <- data.frame("feature"="other", "count"=sum(other_df$count))

    dfs <- as.data.frame(rbind(selected, other, intergenic)) %>%
      mutate(percent = round(count*100 / sum(count), 1)) %>%
      arrange(desc(count))
    rownames(dfs) <- as.character(dfs$feature)

    if(!is.null(outPrefix)) pdf(paste0(outPrefix, ".pdf"), height=8, width=12)
    pbar <- ggplot(dfs, aes(x = reorder(feature, -percent), y = percent)) +
      geom_bar(stat = 'identity', aes(fill = feature)) +
      scale_fill_viridis(discrete=TRUE) +
      geom_text(aes(y = percent + 2), label = dfs$count) +
      labs(x = '', y = 'Percent overlap (%)',
           title="Annotation of peaks to all type of genes",
           subtitle=paste0('(Total number of peaks = ', length(peak), ')')) +
      theme_bw(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust=0.5),
            plot.subtitle=element_text(hjust=0.5))

    
    ## get targeted genes table
    
    features <- get_txdb_features(txdb, fiveP=fiveP, dsTSS=dsTSS, threeP=threeP) # utr5, utr3 and cds here do not contain introns
    dt <- get_targeted_genes(peak, features, stranded=stranded)

    if(verbose) write.table(dt$peak_table, paste(peakLabel, "_peak_annotations.tab", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
    
    dt_gene <- left_join(dt$gene_table, gene_info_table, by=c("tx_name" = "transcript_id"))
    
    dim(dt_gene)
    head(dt_gene)
    tail(dt_gene)
    
    ## filter based on peak type, for ChIPseq peak, only output genes targeted in promoters,
    ## for CLIPseq peaks only output genes targeted in transcripts(5'UTR, CDS, 3'UTR, intron)
    if(fiveP > 0){
       dt_gene <- dt_gene %>%
          #filter(Promoter > 0) %>%
          arrange(desc(Promoter))
    }else{
       dt_gene <- dt_gene %>%
          #filter(Transcript > 0) %>%
          arrange(desc(Transcript))
    }
    
    write.table(dt_gene, paste(peakLabel, "_targeted_annotated_gene.tab", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

    # Plotting overlap counts between query regions and transcript features
    # here utr5, utr3 and cds do not contain introns

    if(verbose) print("Annotation stats:")
    lengths <- vapply(features, FUN=function(x)sum(width(x)), FUN.VALUE=numeric(1))
    feature_counts <- dt$feature_count
    feature_counts <- feature_counts[!names(feature_counts) %in% c("Exon", "Transcript")]
    
    precedence_c <- feature_counts[feature_counts>0]
   
    lengths <- lengths[names(precedence_c)]

    if(verbose) print("plot piechart")

    dfa <- data.frame(count=precedence_c,len=lengths) %>%
      mutate(feature = rownames(.)) %>%
      mutate(percent = count / sum(count)) %>%
      mutate(labels = scales::percent(percent, accuracy = 0.1)) %>%
      mutate(percent = round(percent*100, digits=1)) %>%
      mutate(norm_count=round(count*1000/len, digits=3)) %>%
      mutate(norm_percent = norm_count / sum(norm_count)) %>%
      mutate(norm_labels = scales::percent(norm_percent, accuracy = 0.1)) %>%
      mutate(norm_percent = round(norm_percent*100, digits=1)) %>%
      mutate(feature = forcats::fct_inorder(feature))


    dfa2 <- dfa %>%
       mutate(csum = rev(cumsum(rev(percent))),
              pos = percent/2 + lead(csum, 1),
              pos = ifelse(is.na(pos), percent/2, pos),
              norm_csum = rev(cumsum(rev(norm_percent))),
              norm_pos = norm_percent/2 + lead(norm_csum, 1),
              norm_pos = ifelse(is.na(norm_pos), norm_percent/2, norm_pos))
    if(verbose) print(dfa2)

    apa1 <- ggplot(dfa, aes(x = "", y = percent, fill = feature)) +
       geom_col(color = "white") +
       coord_polar(theta = "y") +
       scale_y_continuous(breaks = dfa2$pos, labels = dfa2$labels) +
       guides(fill = guide_legend(title = "Feature")) +
       scale_fill_viridis(discrete=TRUE) +
       theme(axis.ticks = element_blank(),
             axis.title = element_blank(),
             axis.text = element_text(size = 10),
             legend.position = "top",
             panel.background = element_rect(fill = "white"))
    apa2 <- ggplot(dfa, aes(x = "", y = norm_percent, fill = feature)) +
       geom_col(color = "white") +
       coord_polar(theta = "y") +
       scale_y_continuous(breaks = dfa2$norm_pos, labels = dfa2$norm_labels) +
       guides(fill = guide_legend(title = "Feature")) +
       scale_fill_viridis(discrete=TRUE) +
       theme(axis.ticks = element_blank(),
             axis.title = element_blank(),
             axis.text = element_text(size = 10),
             legend.position = "top",
             panel.background = element_rect(fill = "white"))

    #p1 <- ggpie(df, x="percent", label="labels", lab.pos="out", fill="feature", color="white", palette="ucscgb")

    dfb <- data.frame(count=precedence_c,len=lengths) %>%
      mutate(feature = rownames(.)) %>%
      filter(!feature %in% c("Intron", "NoFeature")) %>%
      mutate(percent = count / sum(count)) %>%
      mutate(labels = scales::percent(percent, accuracy = 0.1)) %>%
      mutate(percent = round(percent*100, digits=1))%>%
      mutate(norm_count=round(count*1000/len, digits=3)) %>%
      mutate(norm_percent = norm_count / sum(norm_count)) %>%
      mutate(norm_labels = scales::percent(norm_percent, accuracy = 0.1)) %>%
      mutate(norm_percent = round(norm_percent*100, digits=1)) %>%
      mutate(feature = forcats::fct_inorder(feature))

    dfb2 <- dfb %>%
       mutate(csum = rev(cumsum(rev(percent))),
              pos = percent/2 + lead(csum, 1),
              pos = ifelse(is.na(pos), percent/2, pos),
              norm_csum = rev(cumsum(rev(norm_percent))),
              norm_pos = norm_percent/2 + lead(norm_csum, 1),
              norm_pos = ifelse(is.na(norm_pos), norm_percent/2, norm_pos))

    apb1 <- ggplot(dfb, aes(x = "", y = percent, fill = feature)) +
      geom_col(color = "white") +
      scale_y_continuous(breaks = dfb2$pos, labels = dfb2$labels) +
      guides(fill = guide_legend(title = "Feature")) +
      scale_fill_viridis(discrete=TRUE) +
      coord_polar(theta = "y") +
      theme(axis.ticks = element_blank(),
             axis.title = element_blank(),
             axis.text = element_text(size = 10),
            legend.position = "top",
             panel.background = element_rect(fill = "white")) 
    
    apb2 <- ggplot(dfb, aes(x = "", y = norm_percent, fill = feature)) +
       geom_col(color = "white") +
       scale_y_continuous(breaks = dfb2$norm_pos, labels = dfb2$norm_labels) +
       guides(fill = guide_legend(title = "Feature")) +
       scale_fill_viridis(discrete=TRUE) +
       coord_polar(theta = "y") +
       theme(axis.ticks = element_blank(),
             axis.title = element_blank(),
             axis.text = element_text(size = 10),
             legend.position = "top",
             panel.background = element_rect(fill = "white")) 

    print(pbar)
    print(plot_grid(apa1, apb1, nrow=1, labels="Absolute counts", label_y=1))
    print(plot_grid(apa2, apb2, nrow=1, labels="Length-normalized counts", label_y=1))
    
    if(!is.null(outPrefix)) on.exit(dev.off(), add=TRUE)
    invisible(list(annotation=dfs, stat=dfa, simplified=dfb))
  }
}

#' @title Plot two-sets Venn diagram
#'
#' @description This is a helper function for Venn diagram plot. A Venn diagram is plotted as output.
#' @param apair a list of two vectors
#' @param overlap_fun the name of the function that defines overlap, depending on the type of object in the vectors.
#'
#' @return NULL
#' @author Shuye Pu
#'
#'
#' @export overlap_pair

overlap_pair <- function(apair, overlap_fun){
   sizes <- vapply(apair, length, numeric(1))
   overlap <- length(Reduce(overlap_fun, apair))
   jaccard <- round(overlap/(sum(sizes)-overlap), digits=5)
   venn.plot <- draw.pairwise.venn(sizes[1], sizes[2], overlap, category=names(apair),
                                   lty=rep("blank",2), fill=c("#0020C2", "#64E986"),
                                   cat.just = rep(list(c(0.5, 0)),2), cex = rep(2, 3), cat.pos = c(0, 0))

   grid.draw(venn.plot)
   grid.text(paste("Jaccard:", jaccard), unit(0.2, "npc"), unit(0.9, "npc"), draw = TRUE)
   grid.newpage()

   return(NULL)
}

#' @title Plot three-sets Venn diagram
#'
#' @description This is a helper function for Venn diagram plot. A Venn diagram is plotted as output.
#' @param atriple a list of three vectors
#' @param overlap_fun the name of the function that defines overlap
#'
#' @return NULL
#' @author Shuye Pu
#'
#'
#' @export overlap_triple

overlap_triple <- function(atriple, overlap_fun){
   sizes <- sort(vapply(atriple, length, numeric(1)), decreasing=TRUE)
   atriple <- atriple[names(sizes)] ## sort the gr by decreasing size to avoid n13 < n123

   overlap12 <- length(Reduce(overlap_fun, atriple[c(1,2)]))
   overlap13 <- length(Reduce(overlap_fun, atriple[c(1,3)]))
   overlap23 <- length(Reduce(overlap_fun, atriple[c(2,3)]))
   overlap123 <- length(Reduce(overlap_fun, atriple))

   venn.plot <- draw.triple.venn(sizes[1], sizes[2], sizes[3], overlap12, overlap23, overlap13,
                                 overlap123, category=names(atriple), lty=rep("blank",3),
                                 fill=c("#0020C2", "#64E986", "#990012"), cat.just = rep(list(c(0.5, 0)),3),
                                 cex = rep(2, 7), cat.pos = c(0, 0, 180))

   grid.draw(venn.plot)
   grid.newpage()

   return(NULL)
}

#' @title Plot four-sets Venn diagram
#'
#' @description This is a helper function for Venn diagram plot. A Venn diagram is plotted as output.
#' @param aquad a list of four vectors
#' @param overlap_fun the name of the function that defines overlap
#'
#' @return NULL
#' @author Shuye Pu
#'
#'
#' @export overlap_quad
#'
overlap_quad <- function(aquad, overlap_fun){
   sizes <- sort(vapply(aquad, length, numeric(1)), decreasing=TRUE)
   aquad <- aquad[names(sizes)] ## sort the gr by decreasing size to avoid n13 < n123

   overlap12 <- length(Reduce(overlap_fun, aquad[c(1,2)]))
   overlap13 <- length(Reduce(overlap_fun, aquad[c(1,3)]))
   overlap14 <- length(Reduce(overlap_fun, aquad[c(1,4)]))
   overlap23 <- length(Reduce(overlap_fun, aquad[c(2,3)]))
   overlap24 <- length(Reduce(overlap_fun, aquad[c(2,4)]))
   overlap34 <- length(Reduce(overlap_fun, aquad[c(3,4)]))
   overlap123 <- length(Reduce(overlap_fun, aquad[c(1,2,3)]))
   overlap124 <- length(Reduce(overlap_fun, aquad[c(1,2,4)]))
   overlap134 <- length(Reduce(overlap_fun, aquad[c(1,3,4)]))
   overlap234 <- length(Reduce(overlap_fun, aquad[c(2,3,4)]))
   overlap1234 <- length(Reduce(overlap_fun, aquad))

   venn.plot <- draw.quad.venn(sizes[1], sizes[2], sizes[3], sizes[4], overlap12, overlap13,
                               overlap14, overlap23, overlap24, overlap34, overlap123, overlap124,
                               overlap134, overlap234, overlap1234, category=names(aquad),
                               lty=rep("blank",4), fill=c("#0020C2", "#64E986", "#990012", "#c6dcff"),
                               cat.just = rep(list(c(0.5, 0)),4), cex = rep(2, 15), cat.pos = c(0, 0, 0, 0))

   grid.draw(venn.plot)
   grid.newpage()
   return(NULL)
}

#' @title Plot Venn diagrams depicting overlap of genomic regions
#' @description This function takes a list of bed file names, and produce a Venn diagram
#' @param bedList a named list of bed files, with length = 2, 3 or 4
#' @param outPrefix a string for plot file name
#' @param handleInputParams a list of parameters for \code{handle_input}
#' @param stranded logical, indicating whether the feature is stranded. For nonstranded feature, only "*" is accepted as strand
#' @param pairOnly logical, indicating whether only pair-wise overlap is desired
#'
#' @return NULL
#' @author Shuye Pu
#'
#' @examples
#' queryFiles <- c(system.file("data", "test_chip_peak_chr19.narrowPeak", package="GenomicPlotData"),
#' system.file("data", "test_chip_peak_chr19.bed", package="GenomicPlotData"),
#' system.file("data", "test_clip_peak_chr19.bed", package="GenomicPlotData"))
#' names(queryFiles) <- c("narrowPeak", "summitPeak", "iCLIPPeak")
#' op <- NULL
#' handleBedParams <- list(fix_width=100, fix_point="center", useScore=FALSE, outRle=FALSE,
#'                         CLIP_reads=FALSE, norm=FALSE, useSizeFactor=FALSE, genome="hg19")
#'
#' plot_overlap_bed(bedList=queryFiles, outPrefix=op, handleInputParams=handleBedParams, pairOnly=FALSE, stranded=FALSE)
#'
#' @export plot_overlap_bed


plot_overlap_bed <- function(bedList, outPrefix=NULL, handleInputParams=NULL, pairOnly=TRUE, stranded=TRUE){

  inputList <- handle_input(bedList, handleInputParams)
  names(inputList) <- names(bedList)
  grList <- lapply(inputList, function(x)x$query)
  sizeList <- lapply(inputList, function(x)x$size)
  if(verbose) print(sizeList)

  # get all pair-wise overlap counts into a matrix, display as a heatmap

  counts <- matrix(rep(0L, length(grList)^2), nrow=length(grList))
  for(i in seq_along(grList)){
     for(j in seq_along(grList)){
        if(stranded){
           counts[i,j] <- length(filter_by_overlaps_stranded(grList[[i]], grList[[j]]))
        }else{
           counts[i,j] <- length(filter_by_overlaps(grList[[i]], grList[[j]]))
        }
        
     }
  }
  rownames(counts) <- colnames(counts) <- names(bedList)
  counts_long <- pivot_longer(as.data.frame(counts), cols=1:ncol(counts), names_to="X", values_to="count") %>%
     mutate(Y=rep(rownames(counts), each=ncol(counts)))

  pairs <- combn(grList, 2, simplify = FALSE)

  if(!is.null(outPrefix)) pdf(paste0(outPrefix, ".pdf"), width=8, height=8)

  g <- ggplot(counts_long, aes(X, Y)) +
     geom_tile(aes(fill = count)) +
     geom_text(aes(label = count, color="white"))  +
     scale_fill_viridis(discrete=FALSE) +
     theme_minimal() +
     theme(legend.position="none",
           axis.title=element_blank(),
           panel.grid.major = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.ticks = element_blank())
  print(g)
  grid.newpage()
  if(stranded){
     lapply(pairs, overlap_pair, filter_by_overlaps_stranded)
     if(!pairOnly){
        if(length(grList)>2){
           triples <- combn(grList, 3, simplify = FALSE)
           lapply(triples, overlap_triple, filter_by_overlaps_stranded)
           if(length(grList)>3){
              quads <- combn(grList, 4, simplify = FALSE)
              lapply(quads, overlap_quad, filter_by_overlaps_stranded)
           }
        }
     }
  }else{
     lapply(pairs, overlap_pair, filter_by_overlaps)
     if(!pairOnly){
        if(length(grList)>2){
           triples <- combn(grList, 3, simplify = FALSE)
           lapply(triples, overlap_triple, filter_by_overlaps)
           if(length(grList)>3){
              quads <- combn(grList, 4, simplify = FALSE)
              lapply(quads, overlap_quad, filter_by_overlaps)
           }
        }
     }
  }
  if(!is.null(outPrefix)) on.exit(dev.off(), add=TRUE)

  return(NULL)
}

#' @title Plot Venn diagrams depicting overlap of gene lists
#' @description This function takes a list of (at most 3) tab-delimited file names, and produce a Venn diagram
#' @param fileList, a named list of tab-dlimited files
#' @param columnList a vector of integers denoting the columns that have gene names in the list of files
#' @param outPrefix, a string for plot file name
#' @param pairOnly, logical, indicating whether only pair-wise overlap is desired
#'
#' @return a list of vectors of gene names
#' @author Shuye Pu
#'
#'
#' @export plot_overlap_genes

plot_overlap_genes <- function(fileList, columnList, pairOnly=TRUE, outPrefix=NULL){

   geneList <- mapply(x=fileList, y=columnList, function(x, y){
      df <- read.delim(x, header=TRUE, sep="\t")
      genes <- unique(df[, y])
      genes
   })


   if(!is.null(outPrefix)){
      pdf(paste0(outPrefix, ".pdf"), width=8, height=8)
      pairs <- combn(geneList, 2, simplify = FALSE)

      lapply(pairs, overlap_pair, intersect)
      if(!pairOnly){
         if(length(geneList)>2){
            triples <- combn(geneList, 3, simplify = FALSE)
            lapply(triples, overlap_triple, intersect)
            if(length(geneList)>3){
               quads <- combn(geneList, 4, simplify = FALSE)
               lapply(quads, overlap_quad, intersect)
            }
         }
      }

      on.exit(dev.off(), add=TRUE)
   }
   invisible(geneList)
}

#' @title Filter GRanges by overlaps in stranded way
#' @description This function reports all query GRanges that have overlaps in subject GRanges. Strand information is used to define overlap.
#' @param query a GRanges object
#' @param subject a GRanges object
#' @param maxgap an integer denoting the distance that define overlap
#' @param minoverlap The minimum amount of overlap between intervals as a single integer greater than 0.
#' If you modify this argument, maxgap must be held fixed.
#'
#' @return a GRanges object
#' @author Shuye Pu
#'
#'
#' @export filter_by_overlaps_stranded

filter_by_overlaps_stranded <- function(query, subject, maxgap=-1L, minoverlap=0L){

  plus_query <- query[strand(query)=="+"]
  minus_query <- query[strand(query)=="-"]
  plus_subject <- subject[strand(subject)=="+"]
  minus_subject <- subject[strand(subject)=="-"]

  overlap_plus <- filter_by_overlaps(plus_query, plus_subject, maxgap=maxgap, minoverlap=minoverlap)
  overlap_minus <- filter_by_overlaps(minus_query, minus_subject, maxgap=maxgap, minoverlap=minoverlap)

  overlaps <- c(overlap_plus, overlap_minus)

  if(length(overlaps) > min(length(query), length(subject))){
     message("Size of overlap is greater than min(sizeOfQuery, sizeOfSubject)!")
  }
  invisible(overlaps)
}

#' @title Filter GRanges by nonoverlaps in stranded way
#' @description This function reports all query GRanges that do not overlaps GRanges in subject. Strand information is used to define overlap.
#' @param query a GRanges object
#' @param subject a GRanges object
#'
#' @return a GRanges object
#' @author Shuye Pu
#'
#'
#' @export filter_by_nonoverlaps_stranded
#'
filter_by_nonoverlaps_stranded <- function(query, subject){

  overlaps <- filter_by_overlaps_stranded(query, subject, maxgap=-1L)
  if(length(overlaps) > min(length(query), length(subject))){
     message("Size of overlap is greater than min(sizeOfQuery, sizeOfSubject!")
  }
  nonoverlaps <- GenomicRanges::setdiff(query, overlaps)
  invisible(nonoverlaps)
}


#' @title Format genomic coordinates in GRanges or GRrangesList as strings used in igv
#' @description This function takes a GRanges or GRangesList object, and transform each range into a string
#' @param x a GRanges or GRangesList object
#' @return a vector of strings in the format of 'chr:start-end(strand)'
#' @author Shuye Pu
#'
#' @examples
#' gr1 <- GenomicRanges::GRanges("chr2", IRanges::IRanges(3, 6))
#' gr2 <- GenomicRanges::GRanges(c("chr1", "chr1"), IRanges::IRanges(c(7,13), width=3),
#'  strand=c("+", "-"))
#' gr3 <- GenomicRanges::GRanges(c("chr1", "chr2"), IRanges::IRanges(c(1, 4), c(3, 9)),
#'  strand="-")
#'
#' grl <- GenomicRanges::GRangesList(gr1= gr1, gr2=gr2, gr3=gr3)
#' grl
#'
#' out <- format_genomic_coordinates(grl)
#' cat(out)
#'
#' @export format_genomic_coordinates
#'
format_genomic_coordinates <- function(x){
   if(grepl("GRangesList", class(x))) x <- unlist(x) ## convert grl to gr

   chr <- as.vector(seqnames(x)) %>% as.character
   start <- start(x) %>% as.numeric
   end <- end(x) %>% as.numeric
   strand <- strand(x) %>% as.character

   out <- paste0(chr, ":", start, "-", end, "(", strand, ")")
   invisible(out)
}

