#' @title Plot signal around custom genomic loci
#' @description  Plot reads or peak Coverage/base/gene of samples given in the query files around reference locus (start, end or center of a genomic region)
#' defined in the centerFiles. The upstream and downstream windows flanking loci can be given separately, a smaller window can be defined to allow
#' statistical comparisons between samples for the same reference, or between references for a given sample. If Input files are provided, ratio over
#' Input is computed and displayed as well.
#'
#' @param queryFiles a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param centerFiles a named vector of reference file names or genomic features in  c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene"). The file should be in .bed format only
#' @param txdb a TxDb object defined in the GenomicFeatures package. Default NULL, needed only when genomic features are used in the place of centerFiles.
#' @param inputFiles a vector of input sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param importParams a list of parameters for \code{handle_input}
#' @param ext a vector of two integers defining upstream and downstream boundaries of the plot window, flanking the reference locus
#' @param hl a vector of two integers defining upstream and downstream boundaries of the highlight window, flanking the reference locus
#' @param stranded logical, indicating whether the strand of the feature should be considered
#' @param scale logical, indicating whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smooth logical, indicating whether the line should smoothed with a spline smoothing algorithm
#' @param heatmap logical, indicating whether a heatmap of the score matrix should be generated
#' @param heatRange a numerical vector of two elements, defining range for heatmap color ramp generation
#' @param rmOutlier a numeric value serving as a multiplier of the MAD in Hampel filter for outliers identification, 0 indicating not removing outliers. For Gaussian distribution, use 3, adjust based on data distribution.
#' @param outPrefix a string specifying output file prefix for plots (outPrefix.pdf)
#' @param refPoint a string in c("start", "center", "end")
#' @param Xlab a string denotes the label on x-axis
#' @param Ylab a string for y-axis label
#' @param shade logical indicating whether to place a shaded rectangle around the point of interest
#' @param binSize an integer defines bin size for intensity calculation
#' @param transform a string in c("log", "log2", "log10"), default = NA indicating no transformation of data matrix
#' @param statsMethod a string in c("wilcox.test", "t.test"), for pair-wise group comparisons
#' @param verbose logical, indicating whether to output additional information (data used for plotting or statistical test results)
#' @param hw a vector of two elements specifying the height and width of the output figures
#' @param nc integer, number of cores for parallel processing
#'
#' @return a list of two dataframes containing the data used for plotting and for statistical testing
#' @author Shuye Pu
#'
#' @export plot_locus


plot_locus <- function(queryFiles,
                       centerFiles,
                       txdb = NULL,
                       ext = c(-100, 100),
                       hl = c(0, 0),
                       shade = TRUE,
                       smooth = FALSE,
                       importParams = NULL,
                       verbose = FALSE,
                       binSize = 10,
                       refPoint = "center",
                       Xlab = "Center",
                       Ylab = "Coverage/base/gene",
                       inputFiles = NULL,
                       stranded = TRUE,
                       heatmap = TRUE,
                       scale = FALSE,
                       outPrefix = NULL,
                       rmOutlier = 0,
                       transform = NA,
                       statsMethod = "wilcox.test",
                       heatRange = NULL,
                       hw = c(8, 8),
                       nc = 2) {
   functionName <- as.character(match.call()[[1]])
   params <- plot_named_list(as.list(environment()))
   force(params)
   
   if (!is.null(outPrefix)) {
      while (!is.null(dev.list())) {
         dev.off()
      }
      pdf(paste(outPrefix, "pdf", sep = "."), height = hw[1], width = hw[2])
   }
   
   if (is.null(inputFiles)) {
      inputLabels <- NULL
      queryInputs <- handle_input(inputFiles = queryFiles, importParams, verbose = verbose, nc = nc)
   } else {
      inputLabels <- names(inputFiles)
      queryLabels <- names(queryFiles)
      if (length(queryFiles) == length(inputFiles)) {
         queryInputs <- handle_input(inputFiles = c(queryFiles, inputFiles), importParams, verbose = verbose, nc = nc)
      } else if (length(inputFiles) == 1) {
         queryInputs <- handle_input(inputFiles = c(queryFiles, inputFiles), importParams, verbose = verbose, nc = nc)
         queryInputs <- queryInputs[c(queryLabels, rep(inputLabels, length(queryLabels)))] # expand the list
         
         inputLabels <- paste0(names(inputFiles), seq_along(queryFiles)) # make each inputLabels unique
         names(queryInputs) <- c(queryLabels, inputLabels)
      } else {
         stop("the number of inputFiles must be 1 or equal to the number of queryFiles!")
      }
   }
   queryLabels <- names(queryInputs)
   
   five <- ext[1] / 1000
   five <- paste0(five, "K")
   if (ext[1] == 0) five <- "-0K"
   three <- ext[2] / 1000
   three <- paste0(three, "K")
   if (ext[2] == 0) three <- "0K"
   featureNames <- c(five, Xlab, three)
   
   ext[2] <- ext[2] - (ext[2] - ext[1]) %% binSize ## to avoid binSize inconsistency, as the final binSize is dictated by bin_num
   bin_num <- round((ext[2] - ext[1]) / binSize)
   colLabel <- seq(ext[1], (ext[2] - binSize), binSize)
   names(colLabel) <- rep(featureNames, c(sum(colLabel < 0), sum(colLabel == 0), sum(colLabel > 0)))
   
   scoreMatrix_list <- list()
   
   bedparam <- importParams
   bedparam$CLIP_reads <- FALSE
   bedparam$fix_width <- 0
   bedparam$useScore <- FALSE
   bedparam$outRle <- FALSE
   bedparam$useSizeFactor <- FALSE
   
   centerInputs <- list()
   
   for (featureName in centerFiles) {
      if (featureName %in% c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene")) {
         featureGR <- get_genomic_feature_coordinates(txdb, featureName, longest = TRUE, protein_coding = TRUE)[["GRanges"]]
         feature <- list("query" = featureGR)
         centerInputs[[featureName]] <- feature
      } else if (file.exists(featureName)) {
         names(featureName) <- names(centerFiles)[centerFiles == featureName]
         feature <- handle_input(featureName, bedparam, verbose = verbose, nc = nc)
         centerInputs[[names(feature)[1]]] <- feature[[1]]
      } else {
         stop("featureName is not supported or the file does not exist, please check your file name and path!")
      }
   }
   
   if (verbose) message("Preparing centers...\n")
   centerLabels <- names(centerInputs)
   centerList <- list()
   
   for (centerLabel in centerLabels) {
      centerInput <- centerInputs[[centerLabel]]
      centerGr <- centerInput$query
      
      if (verbose) message("Center label: ", centerLabel, "\n")
      
      if (refPoint %in% c("center", "start", "end")) {
         windowRegions <- resize(centerGr, width = 1, fix = refPoint)
         windowRegions <- promoters(windowRegions, upstream = -ext[1], downstream = ext[2])
         windowRegions <- check_constraints(windowRegions, importParams$genome)
      } else {
         stop("invalid reference point! Must be one of c('center', 'start', 'end')")
      }
      windowRs <- as(split(windowRegions, f = factor(names(windowRegions))), "GRangesList")
      centerList[[centerLabel]] <- windowRs
      if (verbose) message("Number of window regions ", length(windowRs), "\n")
   }
   
   
   if (verbose) message("Computing coverage for Sample...\n")
   
   for (queryLabel in queryLabels) {
      myInput <- queryInputs[[queryLabel]]
      queryRegions <- myInput$query
      weight_col <- myInput$weight
      libsize <- myInput$size
      
      if (verbose) {
         message("size of query regions: ", libsize, "\n")
         message("Query label: ", queryLabel, "\n")
      }
      
      for (centerLabel in centerLabels) {
         bin_op <- "mean"
         windowRs <- centerList[[centerLabel]]
         fullMatrix <- parallel_scoreMatrixBin(queryRegions, windowRs, bin_num, bin_op, weight_col, stranded, nc = nc)
         if (is.null(inputFiles)) {
            fullMatrix <- process_scoreMatrix(fullMatrix, scale, rmOutlier, transform = transform, verbose = verbose)
         } else {
            fullMatrix <- process_scoreMatrix(fullMatrix, scale = FALSE, rmOutlier = rmOutlier, transform = NA, verbose = verbose)
         }
         colnames(fullMatrix) <- as.character(colLabel)
         rownames(fullMatrix) <- names(windowRs)
         
         scoreMatrix_list[[queryLabel]][[centerLabel]] <- fullMatrix
         
         if (verbose) {
            message("Dimension of fullMatrix:\n")
            message(paste(dim(fullMatrix), collapse = " "), "\n")
            #print("Number of unique window regions")
            # print(length(unique(format_genomic_coordinates(windowRs))))
            # rownames(fullMatrix) <- format_genomic_coordinates(windowRs)
            # write.table(fullMatrix, paste0(queryLabel, "_", centerLabel, "_scoreMatrix.tab"), col.names=NA, sep="\t", quote=FALSE)
         }
      }
   }
   
   plot_df <- list() # per position, averaged over gene
   stat_df <- list() # per gene, averaged over position
   if (heatmap) heatmap_list <- list()
   Ylab <- ifelse(!is.na(transform) && is.null(inputFiles), paste0(transform, " (", Ylab, ")"), Ylab)
   
   if (verbose) message("Collecting coverage data...\n") ## plot multiple bed files on each center
   
   for (queryLabel in queryLabels) {
      if (verbose) message("Query label: ", queryLabel, "\n")
      for (centerLabel in centerLabels) {
         if (verbose) message("Center label: ", centerLabel, "\n")
         
         fullMatrix <- scoreMatrix_list[[queryLabel]][[centerLabel]]
         
         colm <- apply(fullMatrix, 2, mean)
         colsd <- apply(fullMatrix, 2, sd)
         colse <- colsd / sqrt(apply(fullMatrix, 2, length))
         collabel <- colLabel
         querybed <- as.factor(rep(queryLabel, length(colm)))
         refbed <- as.factor(rep(centerLabel, length(colm)))
         
         
         sub_df <- data.frame(colm, colsd, colse, collabel, querybed, refbed)
         colnames(sub_df) <- c("Intensity", "sd", "se", "Position", "Query", "Reference")
         
         if (heatmap) {
            dataname <- paste(Ylab, queryLabel, centerLabel, sep = ":")
            heatmap_list[[dataname]] <- draw_matrix_heatmap(fullMatrix,
                                                            dataName = dataname, labels_col = collabel, ranking = "Sum",
                                                            levels_col = featureNames, ranges = heatRange, verbose = verbose
            )
         }
         
         if (smooth) {
            sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df = as.integer(bin_num / 5))$y)
            sub_df$se <- as.vector(smooth.spline(sub_df$se, df = as.integer(bin_num / 5))$y)
         }
         
         plot_df[[paste(queryLabel, centerLabel, sep = ":")]] <- sub_df
         
         if (hl[2] > hl[1]) {
            xmin <- which(colLabel == hl[1])
            xmax <- which(colLabel == hl[2])
            if (length(xmax) == 0) xmax <- length(colLabel)
            submatrix <- (fullMatrix[, xmin:xmax])
            submatrix[is.na(submatrix)] <- 0
            Intensity <- as.numeric(rowMeans(submatrix))
            
            Query <- as.factor(rep(queryLabel, length(Intensity)))
            Reference <- as.factor(rep(centerLabel, length(Intensity)))
            subdf <- data.frame(Intensity, Query, Reference)
            stat_df[[paste(queryLabel, centerLabel, sep = ":")]] <- subdf
         }
      }
   }
   
   mplot_dt <- bind_rows(plot_df) %>%
      mutate(Group = paste0(Query, ":", Reference), .keep = "all") %>%
      mutate(lower = Intensity - se, upper = Intensity + se, .keep = "all")
   
   mstat_dt <- NULL
   if (hl[2] > hl[1]) {
      mstat_dt <- bind_rows(stat_df) %>%
         mutate(Group = as.factor(paste0(Query, ":", Reference)), .keep = "all")
   }
   
   if (verbose) message("Plotting profile and boxplot...\n")
   
   plot_list <- list()
   queryLabels <- queryLabels[!queryLabels %in% inputLabels]
   for (i in seq_along(queryLabels)) {
      for (beds in combn(queryLabels, i, simplify = FALSE)) {
         for (j in seq_along(centerLabels)) {
            for (centers in combn(centerLabels, j, simplify = FALSE)) {
               if (verbose) {
                  message("Beds: ", paste(beds, collapse = " "), "\n")
                  message("Centers: ", paste(centers, collapse = " "), "\n")
               }
               
               aplot_df <- mplot_dt %>%
                  filter(Query %in% beds & Reference %in% centers)
               aplot_df <- aplot_df %>% # unify order of factors to get consistent color mapping
                  mutate(Query = factor(Query, levels = sort(unique(Query)))) %>%
                  mutate(Reference = factor(Reference, levels = sort(unique(Reference))))
               
               p <- draw_locus_profile(plot_df = aplot_df, cn = "Group", sn = "Group", Xlab = Xlab, Ylab = Ylab, shade = shade, hl = hl)
               
               if ((i == 1 && j > 1) || (i > 1 && j == 1)) {
                  if (hl[2] > hl[1]) {
                     astat_df <- mstat_dt %>%
                        filter(Query %in% beds & Reference %in% centers)
                     astat_df <- astat_df %>% # unify order of factors to get consistent color mapping
                        mutate(Query = factor(Query, levels = sort(unique(Query)))) %>%
                        mutate(Reference = factor(Reference, levels = sort(unique(Reference))))
                     
                     if (j > 1) {
                        p <- draw_locus_profile(plot_df = aplot_df, cn = "Reference", sn = "Query", Xlab = Xlab, Ylab = Ylab, shade = shade, hl = hl)
                        comp <- combn(seq_along(centers), 2, simplify = FALSE)
                        
                        combo <- draw_combo_plot(stat_df = astat_df, xc = "Reference", yc = "Intensity", comp = comp, stats = statsMethod, Ylab = Ylab)
                     } else {
                        p <- draw_locus_profile(plot_df = aplot_df, cn = "Query", sn = "Reference", Xlab = Xlab, Ylab = Ylab, shade = shade, hl = hl)
                        comp <- combn(seq_along(beds), 2, simplify = FALSE)
                        
                        combo <- draw_combo_plot(stat_df = astat_df, xc = "Query", yc = "Intensity", comp = comp, stats = statsMethod, Ylab = Ylab)
                     }
                     
                     lapply(list(p, combo), print)
                  } else {
                     print(p)
                  }
               } else if (i == 1 && j == 1) {
                  plot_list[[paste(Ylab, beds, centers, sep = ":")]] <- p
               } else if (i == length(queryLabels) && j == length(centerLabels)) {
                  print(p)
               }
            }
         }
      }
   }
   
   ## plot query profile and heatmap side by side
   
   rowp <- plot_grid(plotlist = plot_list, nrow = 1, align = "h", axis = "b")
   
   if (heatmap) {
      qheatmap_list <- heatmap_list[names(plot_list)]
      groblist <- lapply(qheatmap_list, function(x) grid.grabExpr(draw(x, heatmap_legend_side = "left")))
      heatp <- plot_grid(plotlist = groblist, nrow = 1, align = "h")
      composite <- plot_grid(rowp, heatp, ncol = 1, align = "v", axis = "l", rel_heights = c(1, 1))
      print(composite)
   } else {
      print(rowp)
   }
   
   
   if (!is.null(inputFiles)) {
      plot_list <- list()
      for (i in seq_along(inputLabels)) {
         for (beds in combn(inputLabels, i, simplify = FALSE)) {
            for (j in seq_along(centerLabels)) {
               for (centers in combn(centerLabels, j, simplify = FALSE)) {
                  if (verbose) {
                     message("Beds: ", paste(beds, collapse = " "), "\n")
                     message("Centers: ", paste(centers, collapse = " "), "\n")
                  }
                  
                  aplot_df <- mplot_dt %>%
                     filter(Query %in% beds & Reference %in% centers)
                  aplot_df <- aplot_df %>% # unify order of factors to get consistent color mapping
                     mutate(Query = factor(Query, levels = sort(unique(Query)))) %>%
                     mutate(Reference = factor(Reference, levels = sort(unique(Reference))))
                  
                  p <- draw_locus_profile(plot_df = aplot_df, cn = "Group", sn = "Group", Xlab = Xlab, Ylab = Ylab, shade = shade, hl = hl)
                  
                  if ((i == 1 && j > 1) || (i > 1 && j == 1)) {
                     if (hl[2] > hl[1]) {
                        astat_df <- mstat_dt %>%
                           filter(Query %in% beds & Reference %in% centers)
                        astat_df <- astat_df %>%
                           mutate(Query = factor(Query, levels = sort(unique(Query)))) %>%
                           mutate(Reference = factor(Reference, levels = sort(unique(Reference))))
                        
                        if (j > 1) {
                           p <- draw_locus_profile(plot_df = aplot_df, cn = "Reference", sn = "Query", Xlab = Xlab, Ylab = Ylab, shade = shade, hl = hl)
                           comp <- combn(seq_along(centers), 2, simplify = FALSE)
                           
                           combo <- draw_combo_plot(stat_df = astat_df, xc = "Reference", yc = "Intensity", comp = comp, stats = statsMethod, Ylab = Ylab)
                        } else {
                           p <- draw_locus_profile(plot_df = aplot_df, cn = "Query", sn = "Reference", Xlab = Xlab, Ylab = Ylab, shade = shade, hl = hl)
                           comp <- combn(seq_along(beds), 2, simplify = FALSE)
                           
                           combo <- draw_combo_plot(stat_df = astat_df, xc = "Query", yc = "Intensity", comp = comp, stats = statsMethod, Ylab = Ylab)
                        }
                        
                        lapply(list(p, combo), print)
                     } else {
                        print(p)
                     }
                  } else if (i == 1 && j == 1) {
                     plot_list[[paste(Ylab, beds, centers, sep = ":")]] <- p
                  } else if (i == length(inputLabels) && j == length(centerLabels)) {
                     print(p)
                  }
               }
            }
         }
      }
      
      
      ## plot input profile and heatmap side by side
      rowp <- plot_grid(plotlist = plot_list, nrow = 1, align = "h", axis = "b")
      # print(rowp)
      if (heatmap) {
         iheatmap_list <- heatmap_list[names(plot_list)]
         groblist <- lapply(iheatmap_list, function(x) grid.grabExpr(draw(x, heatmap_legend_side = "left")))
         heatp <- plot_grid(plotlist = groblist, nrow = 1, align = "h")
         composite <- plot_grid(rowp, heatp, ncol = 1, align = "v", axis = "l", rel_widths = c(1, 1))
         print(composite)
      } else {
         print(rowp)
      }
      
      
      if (verbose) message("Computing Ratio over input...\n")
      Ylab <- ifelse(is.na(transform), "Ratio-over-Input", paste0(transform, " (Ratio-over-Input)"))
      
      inputMatrix_list <- scoreMatrix_list[inputLabels]
      
      ratiolabels <- queryLabels[!queryLabels %in% inputLabels]
      ratioMatrix_list <- scoreMatrix_list[ratiolabels]
      for (centerLabel in centerLabels) {
         for (i in seq_along(ratiolabels)) {
            fullMatrix <- ratio_over_input(ratioMatrix_list[[ratiolabels[i]]][[centerLabel]], inputMatrix_list[[inputLabels[i]]][[centerLabel]], verbose)
            
            fullMatrix <- process_scoreMatrix(fullMatrix, scale, rmOutlier, transform = transform, verbose = verbose)
            colnames(fullMatrix) <- as.character(colLabel)
            
            ratioMatrix_list[[ratiolabels[i]]][[centerLabel]] <- fullMatrix
         }
      }
      
      plot_df <- list()
      stat_df <- list()
      if (heatmap) heatmap_list <- list()
      
      if (verbose) message("Collecting ratio data...\n") ## plot multiple bed files on each center
      
      for (ratiolabel in ratiolabels) {
         if (verbose) message("Ratio label: ", ratiolabel, "\n")
         for (centerLabel in centerLabels) {
            if (verbose) message("Center label: ", centerLabel, "\n")
            
            fullMatrix <- ratioMatrix_list[[ratiolabel]][[centerLabel]]
            # fullMatrix <- process_scoreMatrix(fullMatrix, scale=FALSE, rmOutlier=rmOutlier, transform=transform, verbose=verbose)
            
            colm <- apply(fullMatrix, 2, mean)
            colsd <- apply(fullMatrix, 2, sd)
            colse <- colsd / sqrt(apply(fullMatrix, 2, length))
            collabel <- colLabel
            querybed <- as.factor(rep(ratiolabel, length(colm)))
            refbed <- as.factor(rep(centerLabel, length(colm)))
            
            
            sub_df <- data.frame(colm, colsd, colse, collabel, querybed, refbed)
            colnames(sub_df) <- c("Intensity", "sd", "se", "Position", "Query", "Reference")
            
            if (smooth) {
               sub_df$Intensity <- as.vector(smooth.spline(sub_df$Intensity, df = as.integer(bin_num / 5))$y)
               sub_df$se <- as.vector(smooth.spline(sub_df$se, df = as.integer(bin_num / 5))$y)
            }
            
            plot_df[[paste(ratiolabel, centerLabel, sep = ":")]] <- sub_df
            
            if (heatmap) {
               dataname <- paste(Ylab, ratiolabel, centerLabel, sep = ":")
               heatmap_list[[dataname]] <- draw_matrix_heatmap(fullMatrix, dataName = dataname, labels_col = collabel, levels_col = featureNames, ranges = heatRange, verbose = verbose)
            }
            
            if (hl[2] > hl[1]) {
               xmin <- which(colLabel == hl[1])
               xmax <- which(colLabel == hl[2])
               if (length(xmax) == 0) xmax <- length(colLabel)
               submatrix <- (fullMatrix[, xmin:xmax])
               submatrix[is.na(submatrix)] <- 0
               Intensity <- as.numeric(rowMeans(submatrix))
               
               Query <- as.factor(rep(ratiolabel, length(Intensity)))
               Reference <- as.factor(rep(centerLabel, length(Intensity)))
               subdf <- data.frame(Intensity, Query, Reference)
               stat_df[[paste(ratiolabel, centerLabel, sep = ":")]] <- subdf
            }
         }
      }
      
      mplot_dt <- bind_rows(plot_df) %>%
         mutate(Group = paste0(Query, ":", Reference), .keep = "all") %>%
         mutate(lower = Intensity - se, upper = Intensity + se, .keep = "all")
      
      mstat_dt <- NULL
      
      if (hl[2] > hl[1]) {
         mstat_dt <- bind_rows(stat_df) %>%
            mutate(Group = as.factor(paste0(Query, ":", Reference)), .keep = "all")
      }
      
      if (verbose) message("Plotting ratio profile and boxplot...\n")
      plot_list <- list()
      for (i in seq_along(ratiolabels)) {
         for (beds in combn(ratiolabels, i, simplify = FALSE)) {
            for (j in seq_along(centerLabels)) {
               for (centers in combn(centerLabels, j, simplify = FALSE)) {
                  if (verbose) {
                     message("Beds: ", paste(beds, collapse = " "), "\n")
                     message("Centers: ", paste(centers, collapse = " "), "\n")
                  }
                  
                  aplot_df <- mplot_dt %>%
                     filter(Query %in% beds & Reference %in% centers)
                  aplot_df <- aplot_df %>% # unify order of factors to get consistent color mapping
                     mutate(Query = factor(Query, levels = sort(unique(Query)))) %>%
                     mutate(Reference = factor(Reference, levels = sort(unique(Reference))))
                  
                  p <- draw_locus_profile(plot_df = aplot_df, cn = "Group", sn = "Group", Xlab = Xlab, Ylab = Ylab, shade = shade, hl = hl)
                  
                  if ((i == 1 && j > 1) || (i > 1 && j == 1)) {
                     if (hl[2] > hl[1]) {
                        astat_df <- mstat_dt %>%
                           filter(Query %in% beds & Reference %in% centers)
                        astat_df <- astat_df %>%
                           mutate(Query = factor(Query, levels = sort(unique(Query)))) %>%
                           mutate(Reference = factor(Reference, levels = sort(unique(Reference))))
                        
                        if (j > 1) {
                           p <- draw_locus_profile(plot_df = aplot_df, cn = "Reference", sn = "Query", Xlab = Xlab, Ylab = Ylab, shade = shade, hl = hl)
                           comp <- combn(seq_along(centers), 2, simplify = FALSE)
                           
                           combo <- draw_combo_plot(stat_df = astat_df, xc = "Reference", yc = "Intensity", comp = comp, stats = statsMethod, Ylab = Ylab)
                        } else {
                           p <- draw_locus_profile(plot_df = aplot_df, cn = "Query", sn = "Reference", Xlab = Xlab, Ylab = Ylab, shade = shade, hl = hl)
                           comp <- combn(seq_along(beds), 2, simplify = FALSE)
                           
                           combo <- draw_combo_plot(stat_df = astat_df, xc = "Query", yc = "Intensity", comp = comp, stats = statsMethod, Ylab = Ylab)
                        }
                        
                        lapply(list(p, combo), print)
                     } else {
                        print(p)
                     }
                  } else if (i == 1 && j == 1) {
                     plot_list[[paste(Ylab, beds, centers, sep = ":")]] <- p
                  } else if (i == length(ratiolabels) && j == length(centerLabels)) {
                     print(p)
                  }
               }
            }
         }
      }
      
      ## plot ratio profile and heatmap side by side
      if (heatmap) {
         rheatmap_list <- heatmap_list[names(plot_list)]
         groblist <- lapply(rheatmap_list, function(x) grid.grabExpr(draw(x, heatmap_legend_side = "left")))
         heatp <- plot_grid(plotlist = groblist, nrow = 1, align = "h")
         # print(heatp)
      }
      rowp <- plot_grid(plotlist = plot_list, nrow = 1, align = "h", axis = "b")
      # print(rowp)
      if (heatmap) {
         composite <- plot_grid(rowp, heatp, ncol = 1, align = "v", axis = "l", rel_widths = c(1, 1))
         print(composite)
      } else {
         print(rowp)
      }
   }
   
   if (!is.null(outPrefix)) {
      print(params)
      on.exit(dev.off(), add = TRUE)
   }
   
   invisible(list("plot" = mplot_dt, "stat" = mstat_dt))
}