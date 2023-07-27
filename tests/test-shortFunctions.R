library(GenomicPlot)
library(testthat)

Sys.setenv("R_TESTS" = "")

test_that("testing start and stop of cluster", {
   cl <- start_parallel(2L)
   stop_parallel(cl)
})

test_that("testing inspect_matrix", {
   fullMatrix <- matrix(rnorm(100), ncol = 10)
   for (i in 5:6) {
     fullMatrix[i, 4:7] <- NaN
     fullMatrix[i+1, 4:7] <- NA
     fullMatrix[i+2, 4:7] <- -Inf
     fullMatrix[i-1, 4:7] <- 0
     fullMatrix[i-2, 1:3] <- Inf
   }
   
   GenomicPlot:::inspect_matrix(fullMatrix, verbose = TRUE)  
})

test_that("testing impute_hm", {
   fullMatrix <- matrix(rlnorm(100), ncol = 10)
   for (i in 5:6) {
     fullMatrix[i-1, 4:7] <- 0
   }
   
   imp <- GenomicPlot:::impute_hm(fullMatrix, verbose = TRUE)
})

test_that("testing ratio_over_input", {
   IP <- matrix(rlnorm(100), ncol = 10)
   Input <- matrix(runif(100), ncol = 10)
   
   ratio <- GenomicPlot:::ratio_over_input(IP, Input, verbose = TRUE)  
})

test_that("testing gr2df", {
   gr2 <- GenomicRanges::GRanges(c("chr1", "chr1"), IRanges::IRanges(c(7, 13),
                                                                     width = 3),
     strand = c("+", "-")
   )
   GenomicRanges::mcols(gr2) <- data.frame(score = c(0.3, 0.9),
                                           cat = c(TRUE, FALSE))
   df2 <- gr2df(gr2)
})

test_that("testing aov_TukeyHSD", {
   stat_df <- data.frame(
     Feature = rep(c("A", "B"), c(20, 30)),
     Intensity = c(rnorm(20, 2), rnorm(30, 3))
   )
   
   out <- aov_TukeyHSD(stat_df, xc="Feature")
})

test_that("testing rm_outlier", {
   fullmatrix <- matrix(rnorm(100), ncol = 10)
   maxm <- max(fullmatrix)
   fullmatrix[3, 9] <- maxm + 1000
   fullmatrix[8, 1] <- maxm + 500
   rm_outlier(fullmatrix, verbose = TRUE, multiplier = 100)
   rm_outlier(fullmatrix, verbose = TRUE, multiplier = 1000)
})

test_that("testing process_scoreMatrix", {
   fullMatrix <- matrix(rlnorm(100), ncol = 10)
   for (i in 5:6) {
     fullMatrix[i, 4:7] <- NaN
     fullMatrix[i+1, 4:7] <- NA
     fullMatrix[i+2, 4:7] <- -Inf
     fullMatrix[i-1, 4:7] <- 0
     fullMatrix[i-2, 1:3] <- Inf
   }
   fullMatrix[9, 4:7] <- runif(4) + 90
   
   wo <- process_scoreMatrix(fullMatrix, rmOutlier = 3, verbose = TRUE)
   tf <- process_scoreMatrix(fullMatrix, 
     rmOutlier = 0, transform = "log2", verbose = TRUE
   )
   scaled <- process_scoreMatrix(fullMatrix, scale = TRUE, verbose = TRUE) 
})

test_that("testing rank_rows", {
   fullMatrix <- matrix(rnorm(100), ncol = 10)
   for (i in 5:8) {
     fullMatrix[i, 4:7] <- runif(4) + i
   }
   apply(fullMatrix, 1, sum)
   ranked <- rank_rows(fullMatrix, ranking = "Sum")   
})

test_that("testing overlap_quad", {
   test_list <- list(A = c(1, 2, 3, 4, 5), B = c(4, 5, 7), C = c(1, 3), D = 6)
   overlap_quad(test_list, intersect)
   
   ## GRanges overlap
   query1 <- GRanges("chr19", 
     IRanges(rep(c(10, 15), 2), width=c(1, 20, 40, 50)), 
     strand=c("+", "+", "-", "-")
   )
   
   query2 <- GRanges("chr19", 
     IRanges(rep(c(1, 15), 2), width=c(1, 20, 40, 50)), 
     strand=c("+", "+", "-", "-")
   )
   
   subject1 <- GRanges("chr19", 
     IRanges(rep(c(13, 150), 2), width=c(10, 14, 20, 28)), 
     strand=c("+", "-", "-", "+")
   )
   
   subject2 <- GRanges("chr19", 
     IRanges(rep(c(13, 50), 2), width=c(10, 14, 20, 21)), 
     strand=c("+", "-", "-", "+")
   )
   
   p <- overlap_quad(list(subject1 = subject1, subject2 = subject2,
                          query1 = query1,
     query2 = query2), filter_by_overlaps_stranded
   )  
})

test_that("testing overlap_triple", {
   test_list <- list(A = c(1, 2, 3, 4, 5), B = c(4, 5, 7), C = c(1, 3))
   overlap_triple(test_list, intersect)
   
   ## GRanges overlap
   query <- GRanges("chr19", 
     IRanges(rep(c(10, 15), 2), width=c(1, 20, 40, 50)), 
     strand=c("+", "+", "-", "-")
   )
   
   subject1 <- GRanges("chr19", 
     IRanges(rep(c(13, 150), 2), width=c(10, 14, 20, 28)), 
     strand=c("+", "-", "-", "+")
   )
   
   subject2 <- GRanges("chr19", 
     IRanges(rep(c(13, 50), 2), width=c(10, 14, 20, 21)), 
     strand=c("+", "-", "-", "+")
   )
   
   p <- overlap_triple(list(subject1 = subject1, subject2 = subject2,
                            query = query), filter_by_overlaps_stranded)
})

test_that("testing overlap_pair", {
   test_list <- list(A = c(1, 2, 3, 4, 5), B = c(4, 5, 7))
   overlap_pair(test_list, intersect)
   
   ## GRanges overlap
   query <- GRanges("chr19", 
     IRanges(rep(c(10, 15), 2), width=c(1, 20, 40, 50)), 
     strand=c("+", "+", "-", "-")
   )
   
   subject <- GRanges("chr19", 
     IRanges(rep(c(13, 150), 2), width=c(10, 14, 20, 28)), 
     strand=c("+", "-", "-", "+")
   )
   
   p <- overlap_pair(list(query = query, subject = subject), 
                     filter_by_overlaps_stranded)   
})

test_that("testing draw_combo_plot", {
   stat_df <- data.frame(
     Feature = rep(c("A", "B"), c(200, 300)),
     Intensity = c(rnorm(200, 2, 5), rnorm(300, 3, 5)),
     Height = c(rnorm(200, 5, 5), rnorm(300, 1, 5))
   )
   stat_df_long <- tidyr::pivot_longer(stat_df,
     cols = c(Intensity, Height),
     names_to = "type", values_to = "value"
   )
   
   p <- draw_combo_plot(stat_df_long,
     xc = "Feature", yc = "value", fc = "type",
     Ylab = "value", comp = list(c(1, 2), c(3, 4), c(1, 3), c(2, 4)), nf = 2
   )
})

test_that("testing draw_rank_plot", {
   stat_df <- data.frame(
     Feature = rep(c("A", "B"), c(20, 30)),
     Intensity = c(rlnorm(20, 5, 5), rlnorm(30, 1, 5))
   )
   stat_df1 <- data.frame(
     Feature = rep(c("A", "B"), c(20, 30)),
     Height = c(rnorm(20, 5, 5), rnorm(30, 1, 5))
   )
   
   p1 <- draw_rank_plot(stat_df, xc = "Feature", yc = "Intensity",
                        Ylab = "Intensity")
   p2 <- draw_rank_plot(stat_df1, xc = "Feature", yc = "Height", 
                        Ylab = "Height")   
})

test_that("testing draw_quantile_plot", {
   stat_df <- data.frame(
     Feature = rep(c("A", "B"), c(20, 30)),
     Intensity = c(rnorm(20, 2, 5), rnorm(30, 3, 5)),
     Height = c(rnorm(20, 5, 5), rnorm(30, 1, 5))
   )
   stat_df_long <- tidyr::pivot_longer(stat_df,
     cols = c(Intensity, Height), names_to = "type",
     values_to = "value"
   )
   
   p1 <- draw_quantile_plot(stat_df, xc = "Feature", yc = "Intensity")
   p2 <- draw_quantile_plot(stat_df, xc = "Feature", yc = "Height")
   p3 <- draw_quantile_plot(stat_df_long, xc = "Feature", yc = "value", 
                            fc = "type", Ylab = "value") 
})

test_that("testing draw_mean_se_barplot", {
   stat_df <- data.frame(
        Feature = rep(c("A", "B"), c(20, 30)),
        Intensity = c(rnorm(20, 2), rnorm(30, 3))
      )
      p <- draw_mean_se_barplot(stat_df, xc = "Feature", yc = "Intensity", 
                                Ylab = "Intensity")  
})

test_that("testing draw_boxplot_wo_outlier", {
   stat_df <- data.frame(Feature = rep(c("A", "B"), c(20, 30)), 
                         Intensity = c(rnorm(20, 2), rnorm(30, 3)))
      p <- draw_boxplot_wo_outlier(stat_df,
        xc = "Feature", yc = "Intensity",
        Ylab = "Signal Intensity"
      ) 
})

test_that("testing draw_boxplot_by_factor", {
   stat_df <- data.frame(
     Feature = rep(c("A", "B"), c(20, 30)),
     Intensity = c(rnorm(20, 2, 0.5), rnorm(30, 3, 0.6))
   )
   p <- draw_boxplot_by_factor(stat_df,
     xc = "Feature", yc = "Intensity",
     Ylab = "Signal Intensity"
   )
      
})

test_that("testing draw_locus_profile", {
   library(dplyr)
   Reference <- rep(rep(c("Ref1", "Ref2"), each = 100), 2)
   Query <- rep(c("Query1", "Query2"), each = 200)
   Position <- rep(seq(-50, 49), 4)
   Intensity <- rlnorm(400)
   se <- runif(400)
   df <- data.frame(Intensity, se, Position, Query, Reference) %>%
     mutate(lower = Intensity - se, upper = Intensity + se) %>%
     mutate(Group = paste(Query, Reference, sep = ":"))
   
   p <- draw_locus_profile(df, cn = "Group", shade = TRUE, hl = c(-10, 20)) 
})

test_that("testing draw_region_profile", {
   library(dplyr)
   Reference <- rep(rep(c("Ref1", "Ref2"), each = 100), 2)
   Query <- rep(c("Query1", "Query2"), each = 200)
   Position <- rep(seq_len(100), 4)
   Intensity <- rlnorm(400)
   se <- runif(400)
   df <- data.frame(Intensity, se, Position, Query, Reference) %>%
     mutate(lower = Intensity - se, upper = Intensity + se) %>%
     mutate(Group = paste(Query, Reference, sep = ":"))
   vx <- c(1, 23, 70)
   
   p <- draw_region_profile(df, cn = "Group", vx = vx)
})

test_that("testing draw_region_name", {
   fn <- c("5'UTR", "CDS", "3'UTR")
   bins <- c(5, 15, 5)
   xmax <- 25
   
   p <- draw_region_name(featureNames = fn, scaled_bins = bins, xmax = xmax)
})

test_that("testing draw_region_landmark", {
   fn <- c("5'UTR", "CDS", "3'UTR")
   mark <- c(1, 5, 20)
   xmax <- 25
   
   p <- draw_region_landmark(featureNames = fn, vx = mark, xmax = xmax)
})

test_that("testing draw_matrix_heatmap", {
   fullMatrix <- matrix(rnorm(10000), ncol = 100)
   for (i in seq_len(80)) {
     fullMatrix[i, 16:75] <- runif(60) + i
   }
   labels_col <- as.character(seq_len(100))
   levels_col <- c("start", "center", "end")
   names(labels_col) <- rep(levels_col, c(15, 60, 25))
   
   p <- draw_matrix_heatmap(fullMatrix, dataName = "test", labels_col,
                            levels_col)
})
