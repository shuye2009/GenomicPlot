#' @title Display matrix as a heatmap
#' @description Make a complex heatmap with column annotations
#'
#' @param fullMatrix a numeric matrix
#' @param dataName the nature of the numeric data
#' @param labels_col a named vector for column annotation
#' @param levels_col factor levels for names of labels_col, specifying the
#'  order of labels_col
#' @param ranking method for ranking the rows of the input matrix, options are
#'  c("Sum", "Max", "Hierarchical", "None")
#' @param ranges a numeric vector with three elements, defining custom range for
#'  color ramp, default=NULL, i.e. the range is defined automatically based on
#'  the c(minimun, median, maximum) of fullMatrix
#' @param verbose logical, whether to output the input matrix for inspection
#'
#' @return a grob object
#'
#' @author Shuye Pu
#'
#' @examples
#' fullMatrix <- matrix(rnorm(10000), ncol = 100)
#' for (i in seq_len(80)) {
#'     fullMatrix[i, 16:75] <- runif(60) + i
#' }
#' labels_col <- as.character(seq_len(100))
#' levels_col <- c("start", "center", "end")
#' names(labels_col) <- rep(levels_col, c(15, 60, 25))
#'
#' draw_matrix_heatmap(fullMatrix, dataName = "test", labels_col, levels_col,
#'   ranges = c(-2, 0, 20))
#'
#' @export draw_matrix_heatmap
#'
#'


draw_matrix_heatmap <- function(fullMatrix,
                                dataName = "geneData",
                                labels_col = NULL,
                                levels_col = NULL,
                                ranking = "Sum",
                                ranges = NULL,
                                verbose = FALSE) {
    if (is.null(labels_col)) {
        labels_col <- seq(1, ncol(fullMatrix))
        names(labels_col) <- rep("column", ncol(fullMatrix))
        levels_col <- "column"
    }
    colnames(fullMatrix) <- labels_col

    showRN <- ifelse(nrow(fullMatrix) <= 30, TRUE, FALSE) # show row name?

    fullMatrix <- rank_rows(fullMatrix, ranking)

    if (verbose) {
        message("[draw_matrix_heatmap] Drawing heatmap ...\n")
        vdataName <- gsub(":|/|,|\\.|\\s", "_", dataName, fixed = FALSE)
        ## replace characters not allowed in file names
        write.table(fullMatrix, paste(vdataName, "_matrix.tab", sep = ""),
                    row.names = TRUE, col.names = TRUE, sep = "\t",
                    quote = FALSE)
    }

    features <- factor(names(labels_col), levels = levels_col)
    cols <- ggsci::pal_npg()(length(levels(features)))
    names(cols) <- levels(features)
    mycols <- cols[features]
    names(mycols) <- features

    ha <- HeatmapAnnotation(
        df = data.frame(feature = features), col = list(feature = mycols),
        which = "column", show_legend = FALSE, annotation_label = ""
    )

    if (verbose) {
        message("quantile(fullMatrix, c(seq(0.9, 1, 0.05)), na.rm=TRUE)\n")
        message(paste(quantile(fullMatrix, c(seq(0.9, 1, 0.05)), na.rm = TRUE),
                      collapse = " "), "\n")
    }
    if (is.null(ranges)) {
        ranges <- quantile(fullMatrix, c(0.25, 0.5, 0.75), na.rm = TRUE)
        if (ranges[1] == ranges[3]) {
            message("75% of values are not unique, heatmap may not show
                    signals effectively\n")
          # use quantile of non-zero values
            ranges <- c(min(fullMatrix), quantile(fullMatrix[fullMatrix != 0],
                                                  c(0.5, 0.95), na.rm = TRUE))
        }
    }

    h <- Heatmap(fullMatrix,
        name = "Value",
        col = colorRamp2(ranges, c("#F5F5F5", "#1F45FC", "#151B54")),
        bottom_annotation = ha,
        heatmap_legend_param = list(legend_direction = "horizontal",
                                    title_position = "leftcenter"),
        show_row_names = showRN,
        show_column_names = FALSE,
        show_row_dend = FALSE,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_split = features,
        column_gap = unit(0, "mm"),
        column_title = "%s",
        column_title_gp = gpar(fontsize = 16, fontface = "plain"),
        column_title_side = "bottom"
    )

    message("[draw_matrix_heatmap] finished!\n")
    return(h)
}


#' @title Plot genomic region landmark indicator
#' @description Plot a gene centered polygon for demarcating gene and its
#' upstream and downstream regions
#'
#' @param featureNames a string vector giving names of sub-regions
#' @param vx a vector on integers denoting the x coordinates of start of each
#'  sub-region
#' @param xmax an integer denoting the left most boundary
#' @return a ggplot object
#' @note used by \code{\link{plot_5parts_metagene}}, \code{\link{plot_region}}
#'
#' @author Shuye Pu
#'
#' @examples
#' fn <- c("5'UTR", "CDS", "3'UTR")
#' mark <- c(1, 5, 20)
#' xmax <- 25
#'
#' p <- draw_region_landmark(featureNames = fn, vx = mark, xmax = xmax)
#'
#' @export draw_region_landmark


draw_region_landmark <- function(featureNames, vx, xmax) {
    nfeatures <- length(featureNames)
    cols <- ggsci::pal_npg()(nfeatures)
    names(cols) <- featureNames

    if (nfeatures == 5) {
        values <- data.frame(fid = featureNames,
                             value = cols)
        positions <- data.frame(
            fid = rep(featureNames, each = 4),
            x = c(vx[2], vx[1], vx[1], vx[2], vx[3], vx[2], vx[2], vx[3], vx[4],
                  vx[3], vx[3], vx[4], vx[5], vx[4], vx[4], vx[5], xmax, vx[5],
                  vx[5], xmax),
            y = c(3, 3, 4, 4, 2.5, 2.5, 4.5, 4.5, 2, 2, 5, 5, 2.5, 2.5, 4.5,
                  4.5, 3, 3, 4, 4) - 2
        )
    } else if (nfeatures == 3) {
        values <- data.frame(fid = factor(featureNames), value = cols)
        positions <- data.frame(
            fid = rep(featureNames, each = 4),
            x = c(vx[2], vx[1], vx[1], vx[2], vx[3], vx[2], vx[2], vx[3], xmax,
                  vx[3], vx[3], xmax),
            y = c(3, 3, 4, 4, 2.5, 2.5, 4.5, 4.5, 3, 3, 4, 4) - 2
        )
    } else {
        stop("Number of feautre names must be 3 or 5! other numbers are not
             supported at this point.")
    }

    datapoly <- merge(values, positions, by = c("fid"))

    pp <- ggplot(datapoly, aes(x = x, y = y)) +
        geom_polygon(aes(group = fid, fill = fid)) +
        scale_fill_manual(values = cols) +
        theme(
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none",
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "mm")
        )

    return(pp)
}

#' @title Plot genomic region names
#' @description Plot sub-region labels under the landmark
#'
#' @param featureNames a string vector giving names of sub-regions
#' @param scaled_bins a vector of integers denoting the lengths of each
#'  sub-region
#' @param xmax an integer denoting the right most boundary
#' @return a ggplot object
#' @note used by \code{\link{plot_5parts_metagene}}, \code{\link{plot_region}}
#'
#' @author Shuye Pu
#' @examples
#' fn <- c("5'UTR", "CDS", "3'UTR")
#' bins <- c(5, 15, 5)
#' xmax <- 25
#'
#' p <- draw_region_name(featureNames = fn, scaled_bins = bins, xmax = xmax)
#'
#' @export draw_region_name

draw_region_name <- function(featureNames,
                             scaled_bins,
                             xmax) {
    stopifnot(is.character(featureNames))
    stopifnot(is.numeric(scaled_bins))
    stopifnot(length(featureNames) == length(scaled_bins))

    annotx <- scaled_bins / 2
    for (i in 2:length(scaled_bins)) {
        annotx[i] <- annotx[i] + sum(scaled_bins[seq_len(i - 1)])
    }

    annot <- data.frame(
        fn <- featureNames,
        x <- annotx,
        y <- 0
    )
    ppp <- ggplot(annot, aes(x = x, y = y, label = fn)) +
        geom_text(size = 16/.pt, check_overlap = TRUE) +
        coord_cartesian(xlim = c(1, xmax)) +
        theme(
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none",
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank(),
            plot.margin = unit(c(1, 0, 1, 0), "mm")
        )

    return(ppp)
}

#' @title Plot signal profile in genomic regions
#' @description Plot lines with standard error as the error band
#'
#' @param plot_df a dataframe with column names c(xc, yc, cn, "lower", "upper")
#' @param xc a string denoting column name for values on x-axis
#' @param yc a string denoting column name for numeric data to be plotted
#' @param vx a vector on integers denoting the x coordinates of start of each
#'  sub-region
#' @param cn column name in plot_df for query samples grouping
#' @param sn column name in plot_df for subject name to be shown in the plot
#'  title
#' @param Ylab a string for Y-axis label
#' @return a ggplot object
#' @note used by \code{\link{plot_5parts_metagene}}, \code{\link{plot_region}}
#'
#' @author Shuye Pu
#'
#' @examples
#' library(dplyr)
#' Reference <- rep(rep(c("Ref1", "Ref2"), each = 100), 2)
#' Query <- rep(c("Query1", "Query2"), each = 200)
#' Position <- rep(seq_len(100), 4)
#' Intensity <- rlnorm(400)
#' se <- runif(400)
#' df <- data.frame(Intensity, se, Position, Query, Reference) %>%
#'     mutate(lower = Intensity - se, upper = Intensity + se) %>%
#'     mutate(Group = paste(Query, Reference, sep = ":"))
#' vx <- c(1, 23, 70)
#'
#' p <- draw_region_profile(df, cn = "Group", vx = vx)
#' p
#'
#' @export draw_region_profile
#'
#'

draw_region_profile <- function(plot_df,
                                xc = "Position",
                                yc = "Intensity",
                                cn = "Query",
                                sn = NULL,
                                Ylab = "Signal Intensity",
                                vx) {
    stopifnot(c(xc, yc, cn) %in% colnames(plot_df))
    message("[draw_region_profile] started ...\n")
    p <- ggplot(plot_df, aes(x = .data[[xc]], y = .data[[yc]],
                             color = .data[[cn]])) +
        scale_fill_npg() +
        scale_color_npg() +
        geom_line(size = 1.25) +
        geom_vline(xintercept = vx[2:length(vx)], linetype = "dotted",
                   color = "blue", size = 0.5) +
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = .data[[cn]]),
                    linetype = 0, alpha = 0.3) +
        theme_classic() +
        ylab(Ylab) +
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(face="bold", size=18),
            axis.text.y = element_text(size=16),
            plot.margin = unit(c(1, 1, 0, 1), "lines")
        )
    if(!is.null(sn))
        p <- p + ggtitle(paste(unique(plot_df[[sn]]), collapse = ", "))

    message("[draw_region_profile] finished!\n")
    return(p)
}

#' @title Plot signal profile around genomic loci
#' @description Plot lines with standard error as the error band
#'
#' @param plot_df a dataframe with column names c(xc, yc, cn, "lower", "upper")
#' @param xc a string denoting column name for values on x-axis
#' @param yc a string denoting column name for numeric data to be plotted
#' @param cn a string denoting column name for sample grouping, like 'Query' or
#'  'Reference'
#' @param sn a string denoting column name for the subject of sample grouping,
#'  if 'cn' is 'Query', then 'sn' will be 'Reference'
#' @param Xlab a string for x-axis label
#' @param Ylab a string for y-axis label
#' @param hl a vector of two integers defining upstream and downstream
#'  boundaries of the rectangle
#' @param shade logical indicating whether to place a shaded rectangle
#'  around the loci bounded by hl
#'
#' @return a ggplot object
#' @note used by \code{\link{plot_locus}}, \code{\link{plot_locus_with_random}}
#' @author Shuye Pu
#'
#' @examples
#' library(dplyr)
#' Reference <- rep(rep(c("Ref1", "Ref2"), each = 100), 2)
#' Query <- rep(c("Query1", "Query2"), each = 200)
#' Position <- rep(seq(-50, 49), 4)
#' Intensity <- rlnorm(400)
#' se <- runif(400)
#' df <- data.frame(Intensity, se, Position, Query, Reference) %>%
#'     mutate(lower = Intensity - se, upper = Intensity + se) %>%
#'     mutate(Group = paste(Query, Reference, sep = ":"))
#'
#' p <- draw_locus_profile(df, cn = "Group", shade = TRUE, hl = c(-10, 20))
#' p
#'
#' @export draw_locus_profile
#'

draw_locus_profile <- function(plot_df,
                               xc = "Position",
                               yc = "Intensity",
                               cn = "Query",
                               sn = NULL,
                               Xlab = "Center",
                               Ylab = "Signal Intensity",
                               shade = FALSE,
                               hl = c(0, 0)) {
    stopifnot(c(xc, yc, cn) %in% colnames(plot_df))

  message("[draw_locus_profile] started ...\n")
    p <- ggplot(plot_df, aes(x = .data[[xc]], y = .data[[yc]],
                             color = .data[[cn]])) +
        scale_fill_npg() +
        scale_color_npg() +
        geom_line(size = 1.25) +
        geom_vline(xintercept = 0, linetype = "dotted", color = "blue",
                   size = 0.5) +
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = .data[[cn]]),
                    linetype = 0, alpha = 0.3) +
        theme_classic() +
        xlab(Xlab) +
        ylab(Ylab) +
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            axis.text = element_text(face = "plain", size = 16),
            axis.title = element_text(face = "bold", size = 18)
        )

    if(!is.null(sn))
        p <- p + ggtitle(paste(unique(plot_df[[sn]]), collapse = ", "))

    if (shade) p <- p + annotate("rect", xmin = hl[1], xmax = hl[2],
                                 ymin = -Inf, ymax = Inf, fill = "grey",
                                 color = "grey", alpha = 0.3)

    message("[draw_locus_profile] finished\n")
    return(p)
}


#' @title draw stacked plot
#' @description Plot profile on top of heatmap, and align feature labels.
#'
#' @param plot_list a list of profile plots
#' @param heatmap_list a list of heatmaps
#'
#' @return a null value
#' @note used by \code{\link{plot_locus}}, \code{\link{plot_5parts_metagene}},
#'    \code{\link{plot_region}}
#' @author Shuye Pu
#'
#' @export draw_stacked_plot
#'

draw_stacked_plot <- function(plot_list, heatmap_list){
  if (length(heatmap_list) > 0) {
    groblist <- lapply(heatmap_list, function(x)
      grid.grabExpr(draw(x, heatmap_legend_side = "bottom")))
    names(groblist) <- names(heatmap_list)
  }else{
    groblist <- NULL
  }

  for(i in seq_along(plot_list)){
    if(!is.null(groblist)){
      composite <- ggdraw() +
        draw_plot(plot_list[[i]], 0, 0.5, 1, 0.5) +
        draw_plot(groblist[[i]], 0.16, 0, 0.8, 0.5, halign = 0.7)
      print(composite)
    }else{
      print(plot_list[[i]])
    }
  }
  return(NULL)
}

#' @title Plot boxplot with two factors
#' @description Plot violin plot with boxplot components for data with one or
#' two factors, p-value significance levels are displayed, "***" = 0.001,
#' "**" = 0.01, "*" = 0.05.
#'
#' @param stat_df a dataframe with column names c(xc, yc)
#' @param xc a string denoting column name for grouping
#' @param fc a string denoting column name for sub-grouping based on an
#'  additional factor
#' @param yc a string denoting column name for numeric data to be plotted
#' @param Xlab a string for x-axis label
#' @param Ylab a string for y-axis label
#' @param comp a list of vectors denoting pair-wise comparisons to be performed
#'  between groups
#' @param stats the name of pair-wise statistical tests, like t.test or
#'  wilcox.test
#' @param nf a integer normalizing factor for correct count of observations when
#'  the data table has two factors, such as those produced by `pivot_longer`,
#'  equals to the number of factors
#'
#' @return a ggplot object
#' @note used by \code{\link{plot_locus}}, \code{\link{plot_locus_with_random}},
#'    \code{\link{plot_region}}
#' @author Shuye Pu
#'
#' @examples
#' stat_df <- data.frame(
#'     Feature = rep(c("A", "B"), c(20, 30)),
#'     Intensity = c(rnorm(20, 2, 0.5), rnorm(30, 3, 0.6))
#' )
#' p <- draw_boxplot_by_factor(stat_df,
#'     xc = "Feature", yc = "Intensity",
#'     Ylab = "Signal Intensity"
#' )
#' p
#' @export draw_boxplot_by_factor
#'

draw_boxplot_by_factor <- function(stat_df,
                                   xc = "Feature",
                                   yc = "Intensity",
                                   fc = xc,
                                   comp = list(c(1, 2)),
                                   stats = "wilcox.test",
                                   Xlab = xc,
                                   Ylab = yc,
                                   nf = 1) {
    stopifnot(c(xc, yc, fc) %in% colnames(stat_df))
    xlabs <- paste(levels(as.factor(stat_df[[xc]])), "\n(",
                   table(stat_df[[xc]]) / nf, ")", sep = "")
    ypos <- rep(max(stat_df[[yc]]), length(comp)) *
      seq(1, 1 + (length(comp) - 1) * 0.1, 0.1)
    outlier.shape <- 19

    if (fc == xc) {
        p <- ggplot(stat_df, aes(x = .data[[xc]], y = .data[[yc]],
                                 fill = .data[[fc]])) +
            geom_violin(width = 0.5) +
            geom_boxplot(width = 0.2, outlier.shape = outlier.shape) +
            scale_fill_npg() +
            scale_color_npg() +
            theme_classic() +
            theme(
                axis.text = element_text(face = "plain", size = 12,
                                         color = "black"),
                # axis.title.x = element_blank(),
                axis.title = element_text(face = "bold", size = 14,
                                          color = "black"),
                legend.position = "none"
            ) +
            labs(y = Ylab, x = Xlab) +
            geom_signif(comparisons = comp, test = stats,
                        map_signif_level = TRUE, y_position = ypos) +
            scale_x_discrete(labels = xlabs)
    } else {
        mid <- function(v) {
            m <- rep(0, (length(v) / 2))
            for (i in seq_along(m)) {
                m[i] <- (v[i * 2 - 1] + v[i * 2]) / 2
            }
            return(m)
        }
        stat_df <- stat_df %>%
            mutate(x2 = as.integer(interaction(.data[[fc]], .data[[xc]])))
        p <- ggplot(stat_df, aes(x = x2, y = .data[[yc]], group = x2,
                                 fill = .data[[fc]])) +
            geom_violin(width = 0.5) +
            geom_boxplot(width = 0.2, outlier.shape = outlier.shape) +
            scale_fill_npg() +
            scale_color_npg() +
            theme_classic() +
            theme(
                axis.text = element_text(face = "plain", size = 12,
                                         color = "black"),
                axis.title.x = element_blank(),
                axis.title = element_text(face = "bold", size = 14,
                                          color = "black"),
                legend.position = "bottom"
            ) +
            labs(y = Ylab) +
            geom_signif(comparisons = comp, test = stats,
                        map_signif_level = TRUE, y_position = ypos) +
            scale_x_continuous(breaks = mid(sort(unique(stat_df$x2))),
                               labels = xlabs)
    }

    return(p)
}

#' @title Plot boxplot without outliers
#' @description Plot boxplot without outliers, useful when outliers have a wide
#' range and the median is squeezed at the bottom of the plot. The p-value
#' significance level is the same as those in
#' \code{\link{draw_boxplot_by_factor}}, but not displayed.
#'
#' @param stat_df a dataframe with column names c(xc, yc)
#' @param xc a string denoting column name for grouping
#' @param fc a string denoting column name for sub-grouping
#' @param yc a string denoting column name for numeric data to be plotted
#' @param Xlab a string for x-axis label
#' @param Ylab a string for y-axis label
#' @param comp a list of vectors denoting pair-wise comparisons to be performed
#'  between groups
#' @param stats the name of pair-wise statistical tests, like t.test or
#'  wilcox.test
#' @param nf a integer normalizing factor for correct count of observations when
#'  the data table has two factors, such as those produced by `pivot_longer`,
#'  equals to the number of factors
#'
#' @return a ggplot object
#' @export draw_boxplot_wo_outlier
#' @examples
#' stat_df <- data.frame(
#'     Feature = rep(c("A", "B"), c(20, 30)),
#'     Intensity = c(rnorm(20, 2), rnorm(30, 3))
#' )
#'
#' p <- draw_boxplot_wo_outlier(stat_df,
#'     xc = "Feature", yc = "Intensity",
#'     Ylab = "Signal Intensity"
#' )
#' p
#'
draw_boxplot_wo_outlier <- function(stat_df,
                                    xc = "Feature",
                                    yc = "Intensity",
                                    fc = xc,
                                    comp = list(c(1, 2)),
                                    stats = "wilcox.test",
                                    Xlab = xc,
                                    Ylab = yc,
                                    nf = 1) {
    stopifnot(c(xc, yc, fc) %in% colnames(stat_df))

    xlabs <- paste(levels(as.factor(stat_df[[xc]])), "\n(",
                   table(stat_df[[xc]]) / nf, ")", sep = "")
    fomu <- as.formula(paste(yc, "~", xc))
    bp <- boxplot(fomu, stat_df, plot = FALSE)

    lim <- c(min(bp$stats), max(bp$stats))

    if (fc == xc) {
        p <- ggplot(stat_df, aes(x = .data[[xc]], y = .data[[yc]],
                                 fill = .data[[fc]])) +
            geom_boxplot(width = 0.2, outlier.shape = NA) +
            scale_fill_npg() +
            scale_color_npg() +
            theme_classic() +
            theme(
                axis.text = element_text(face = "plain", size = 12,
                                         color = "black"),
                axis.title = element_text(face = "bold", size = 14,
                                          color = "black"),
                legend.position = "none"
            ) +
            labs(y = Ylab, x = Xlab) +
            coord_cartesian(ylim = lim) +
            scale_x_discrete(labels = xlabs)
    } else {
        mid <- function(v) {
            m <- rep(0, (length(v) / 2))
            for (i in seq_along(m)) {
                m[i] <- (v[i * 2 - 1] + v[i * 2]) / 2
            }
            return(m)
        }
        stat_df <- stat_df %>%
            mutate(x2 = as.integer(interaction(.data[[fc]], .data[[xc]])))
        p <- ggplot(stat_df, aes(x = x2, y = .data[[yc]], group = x2,
                                 fill = .data[[fc]])) +
            geom_boxplot(width = 0.2, outlier.shape = NA) +
            scale_fill_npg() +
            scale_color_npg() +
            theme_classic() +
            theme(
                axis.text = element_text(face = "plain", size = 12,
                                         color = "black"),
                axis.title.x = element_blank(),
                axis.title = element_text(face = "bold", size = 14,
                                          color = "black"),
                legend.position = "bottom"
            ) +
            labs(y = Ylab) +
            scale_x_continuous(breaks = mid(sort(unique(stat_df$x2))),
                               labels = xlabs) +
            coord_cartesian(ylim = lim)
    }

    return(p)
}
#' @title Plot barplot for mean with standard error bars
#' @description Plot barplot for mean with standard error bars, no p-value
#' significance levels are displayed, but ANOVA p-value is provided as tag and
#' TukeyHSD test are displayed as caption.
#'
#' @param stat_df a dataframe with column names c(xc, yc)
#' @param xc a string denoting column name for grouping
#' @param yc a string denoting column name for numeric data to be plotted
#' @param fc a string denoting column name for sub-grouping based on an
#'  additional factor
#' @param comp a list of vectors denoting pair-wise comparisons to be performed
#'  between groups
#' @param Xlab a string for x-axis label
#' @param Ylab a string for y-axis label
#' @param Ylim a numeric vector of two elements, defining custom limits of
#'  y-axis
#' @param nf a integer normalizing factor for correct count of observations when
#'  the data table has two factors, such as those produced by pivot_longer,
#'  equals to the number of factors
#'
#' @return a ggplot object
#' @note used by \code{\link{plot_locus}}, \code{\link{plot_locus_with_random}}
#' @author Shuye Pu
#'
#' @examples
#' stat_df <- data.frame(
#'     Feature = rep(c("A", "B"), c(20, 30)),
#'     Intensity = c(rnorm(20, 2), rnorm(30, 3))
#' )
#' p <- draw_mean_se_barplot(stat_df,
#'     xc = "Feature", yc = "Intensity",
#'     Ylab = "Intensity"
#' )
#' p
#' @export draw_mean_se_barplot
#'

draw_mean_se_barplot <- function(stat_df,
                                 xc = "Feature",
                                 yc = "Intensity",
                                 fc = xc,
                                 comp = list(c(1, 2)),
                                 Xlab = xc,
                                 Ylab = yc,
                                 Ylim = NULL,
                                 nf = 1) {
    stopifnot(c(xc, yc, fc) %in% colnames(stat_df))
    if (fc == xc) {
        stat_df[[xc]] <- as.factor(stat_df[[xc]])

        stats <- aov_TukeyHSD(stat_df, xc, yc)

        means_se <- stat_df %>%
            group_by(.data[[xc]]) %>%
            summarize(
                mean_Intensity = mean(.data[[yc]]),
                sd_Intensity = sd(.data[[yc]]),
                N_Intensity = length(.data[[yc]]),
                se = sd_Intensity / sqrt(N_Intensity),
                upper_limit = mean_Intensity + se,
                lower_limit = mean_Intensity - se
            )
        means_se <- means_se %>%
            mutate(labelx = paste0(.data[[xc]], "\n(", N_Intensity / nf, ")"))
        ## now .data is means_se, not stat_df anymore
        levels(means_se[[xc]]) <- levels(stat_df[[xc]])

        p <- ggplot(means_se, aes(x = .data[[xc]], y = mean_Intensity,
                                  fill = .data[[xc]])) +
            scale_fill_npg() +
            scale_color_npg() +
            geom_col(position = "identity") +
            geom_errorbar(aes(ymin = lower_limit, ymax = upper_limit),
                          position = position_dodge(width = 0.2), width = 0.2) +
            theme_classic() +
            theme(
                axis.title = element_text(face = "bold", size = 14,
                                          color = "black", vjust = 0.25),
                axis.title.x = element_blank(),
                axis.text = element_text(face = "plain", size = 12,
                                         color = "black"),
                legend.position = "none"
            ) +
            labs(y = paste0(Ylab, "\n(mean \u00b1 SE)"), x = Xlab) +
            scale_x_discrete(labels = means_se$labelx) +
            ggtitle(label = paste("ANOVA p-value =",
                                  format(stats$ANOVA, digits = 3)))

        if (!is.null(Ylim)) {
            p <- p + coord_cartesian(ylim = Ylim)
        }

        stats$HSD[, seq_len(3)] <- round(stats$HSD[, seq_len(3)], digits = 3)
        stats$HSD[, 4] <- format(stats$HSD[, 4], digits = 3)

        comp_row <- vapply(comp, function(x) {
            arow <- paste0(levels(stat_df[[xc]])[x[2]], "-",
                           levels(stat_df[[xc]])[x[1]])
        }, character(1))
        stats_selected <- as.data.frame(stats$HSD) %>%
            filter(row.names(stats$HSD) %in% comp_row)
        ptable <- tab_add_title(ggtexttable(stats_selected,
                                            theme = ttheme(base_size = 9)),
                                text = "post hoc TukeyHSD test",
                                padding = unit(1.0, "line"), just = "left")

        outp <- cowplot::plot_grid(p, ptable, ncol = 1, align = "l",
                                   rel_heights = c(3, 2))
    } else {
        stat_df <- stat_df %>%
            mutate(x2 = as.factor(as.integer(interaction(.data[[fc]],
                                                         .data[[xc]]))))
        stats <- aov_TukeyHSD(stat_df, "x2", yc)

        means_se <- stat_df %>%
            group_by(.data[["x2"]]) %>%
            summarize(
                mean_Intensity = mean(.data[[yc]]),
                sd_Intensity = sd(.data[[yc]]),
                N_Intensity = length(.data[[yc]]),
                se = sd_Intensity / sqrt(N_Intensity),
                upper_limit = mean_Intensity + se,
                lower_limit = mean_Intensity - se
            )

        levels(means_se[["x2"]]) <- levels(stat_df[["x2"]])
        means_se <- means_se %>%
            mutate(labelx = paste0(.data[["x2"]], "\n(", N_Intensity, ")"))
        ## now .data is means_se, not stat_df anymore

        p <- ggplot(means_se, aes(x = x2, y = mean_Intensity, group = x2,
                                  fill = x2)) +
            scale_fill_npg() +
            scale_color_npg() +
            geom_col(position = "identity") +
            geom_errorbar(aes(ymin = lower_limit, ymax = upper_limit),
                          position = position_dodge(width = 0.2), width = 0.2) +
            theme_classic() +
            theme(
                axis.title = element_text(face = "bold", size = 14,
                                          color = "black", vjust = 0.25),
                axis.title.x = element_blank(),
                axis.text = element_text(face = "plain", size = 12,
                                         color = "black"),
                legend.position = "none"
            ) +
            labs(y = paste0(Ylab, "(mean \u00b1 SE)"), x = Xlab) +
            scale_x_discrete(labels = means_se$labelx) +
            ggtitle(label = paste("ANOVA p-value =",
                                  format(stats$ANOVA, digits = 3)))

        if (!is.null(Ylim)) {
            p <- p + coord_cartesian(ylim = Ylim)
        }

        stats$HSD[, seq_len(3)] <- round(stats$HSD[, seq_len(3)], digits = 3)
        stats$HSD[, 4] <- format(stats$HSD[, 4], digits = 3)

        comp_row <- vapply(comp, function(x) {
            arow <- paste0(levels(stat_df[["x2"]])[x[2]], "-",
                           levels(stat_df[["x2"]])[x[1]])
        }, character(1))
        stats_selected <- as.data.frame(stats$HSD) %>%
            filter(row.names(stats$HSD) %in% comp_row)
        ptable <- tab_add_title(ggtexttable(stats_selected,
                                            theme = ttheme(base_size = 9)),
                                text = "post hoc TukeyHSD test",
                                padding = unit(1.0, "line"), just = "left")

        outp <- cowplot::plot_grid(p, ptable, ncol = 1, align = "l",
                                   rel_heights = c(3, 2))
    }

    return(outp)
}

#' @title Plot quantile over value
#' @description Plot quantiles as y-axis, and values as x-axis. Same as
#' `geom_ecdf`, but allows sub-grouping by a second factor.
#'
#' @param stat_df a dataframe with column names c(xc, yc)
#' @param xc a string denoting column name for grouping
#' @param yc a string denoting column name for numeric data to be plotted
#' @param fc a string denoting column name for sub-grouping based on an
#'  additional factor
#' @param Ylab a string for y-axis label
#'
#' @return a ggplot object
#' @note used by \code{\link{plot_locus}}, \code{\link{plot_locus_with_random}}
#' @author Shuye Pu
#'
#' @export draw_quantile_plot
#'
#' @examples
#' stat_df <- data.frame(
#'     Feature = rep(c("A", "B"), c(20, 30)),
#'     Intensity = c(rnorm(20, 2, 5), rnorm(30, 3, 5)),
#'     Height = c(rnorm(20, 5, 5), rnorm(30, 1, 5))
#' )
#' stat_df_long <- tidyr::pivot_longer(stat_df,
#'     cols = c(Intensity, Height), names_to = "type",
#'     values_to = "value"
#' )
#'
#' print(draw_quantile_plot(stat_df, xc = "Feature", yc = "Intensity"))
#' print(draw_quantile_plot(stat_df, xc = "Feature", yc = "Height"))
#' print(draw_quantile_plot(stat_df_long,
#'     xc = "Feature", yc = "value",
#'     fc = "type", Ylab = "value"
#' ))
#'
draw_quantile_plot <- function(stat_df,
                               xc = "Feature",
                               yc = "Intensity",
                               Ylab = yc,
                               fc = xc) {
    stopifnot(c(xc, yc, fc) %in% colnames(stat_df))
    if (fc == xc) {
        long_df <- stat_df %>%
            group_by(.data[[xc]]) %>%
            mutate(Fraction = ecdf(.data[[yc]])(.data[[yc]]))

        p <- ggplot(data = long_df, aes(x = .data[[yc]], y = Fraction,
                                        color = .data[[xc]])) +
            scale_color_npg() +
            geom_line(size = 1.25) +
            labs(x = Ylab, y = "Cumulative fraction") +
            theme_classic() +
            theme(
                legend.position = "top",
                legend.title = element_blank(),
                axis.text = element_text(angle = 0, size = 12, vjust = 0),
                axis.title = element_text(face = "bold", color = "black",
                                          size = 14, vjust = 0.25)
            )
    } else {
        stat_df <- stat_df %>%
            mutate(x2 = as.integer(interaction(.data[[fc]], .data[[xc]])))
        long_df <- stat_df %>%
            group_by(x2) %>%
            mutate(Fraction = ecdf(.data[[yc]])(.data[[yc]]))

        p <- ggplot(data = long_df, aes(x = .data[[yc]], y = Fraction,
                                        color = as.factor(x2))) +
            scale_color_npg() +
            geom_line(size = 1.25) +
            labs(x = Ylab, y = "Cumulative fraction") +
            theme_classic() +
            theme(
                legend.position = "top",
                legend.title = element_blank(),
                axis.text = element_text(angle = 0, size = 12, vjust = 0),
                axis.title = element_text(face = "bold", color = "black",
                                          size = 14, vjust = 0.25)
            )
    }

    return(p)
}

#' @title Plot fraction of cumulative sum over rank
#' @description Plot cumulative sum over rank as line plot, both cumulative sum
#' and rank are scaled between 0 and 1. This is the same as the fingerprint plot
#' of the deepTools.
#'
#' @param stat_df a dataframe with column names c(xc, yc)
#' @param xc a string denoting column name for grouping
#' @param yc a string denoting column name for numeric data to be plotted
#' @param Ylab a string for y-axis label
#'
#' @return a ggplot object
#'
#' @author Shuye Pu
#'
#' @export draw_rank_plot
#'
#' @examples
#' stat_df <- data.frame(
#'     Feature = rep(c("A", "B"), c(20, 30)),
#'     Intensity = c(rlnorm(20, 5, 5), rlnorm(30, 1, 5))
#' )
#' stat_df1 <- data.frame(
#'     Feature = rep(c("A", "B"), c(20, 30)),
#'     Height = c(rnorm(20, 5, 5), rnorm(30, 1, 5))
#' )
#'
#' print(draw_rank_plot(stat_df,
#'     xc = "Feature", yc = "Intensity",
#'     Ylab = "Intensity"
#' ))
#' print(draw_rank_plot(stat_df1,
#'     xc = "Feature", yc = "Height",
#'     Ylab = "Height"
#' ))
#'
draw_rank_plot <- function(stat_df,
                           xc = "Feature",
                           yc = "Intensity",
                           Ylab = yc) {
    stopifnot(c(xc, yc) %in% colnames(stat_df))
    long_df <- stat_df %>%
        group_by(.data[[xc]]) %>%
        arrange(.data[[yc]]) %>%
        mutate(cumCount = cumsum(.data[[yc]])) %>%
        mutate(Rank = rank(cumCount)) %>%
        mutate(Fraction = cumCount / max(cumCount), Rank = Rank / max(Rank))

    p <- ggplot(data = long_df, aes(x = Rank, y = Fraction,
                                    color = .data[[xc]])) +
        scale_color_npg() +
        geom_line(size = 1.25) +
        labs(x = paste0("Rank (", Ylab, ")"), y = "Cumulative sum fraction") +
        theme_classic() +
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            axis.text = element_text(angle = 0, size = 12, vjust = 0),
            axis.title = element_text(face = "bold", color = "black", size = 14,
                                      vjust = 0.25)
        ) +
        ggtitle(paste("Cumulative sum fraction of", Ylab))

    return(p)
}

#' @title Make combo plot for statistics plots
#' @description Place violin plot, boxplot without outliers, mean+se barplot and
#' quantile plot on the same page
#'
#' @param stat_df a dataframe with column names c(xc, yc)
#' @param xc a string denoting column name for grouping
#' @param yc a string denoting column name for numeric data to be plotted
#' @param fc a string denoting column name for sub-grouping based on an
#'  additional factor
#' @param comp a list of vectors denoting pair-wise comparisons to be performed
#'  between groups
#' @param stats the name of pair-wise statistical tests, like t.test or
#'  wilcox.test
#' @param Xlab a string for x-axis label
#' @param Ylab a string for y-axis label
#' @param Ylim a numeric vector of two elements, defining custom limits of y-axis
#' @param title a string for plot title
#' @param nf a integer normalizing factor for correct count of observations when
#'  the data table has two factors, such as those produced by pivot_longer,
#'  equals to the number of factors
#'
#' @return a ggplot object
#'
#' @author Shuye Pu
#'
#' @examples
#'
#' stat_df <- data.frame(
#'     Feature = rep(c("A", "B"), c(200, 300)),
#'     Intensity = c(rnorm(200, 2, 5), rnorm(300, 3, 5)),
#'     Height = c(rnorm(200, 5, 5), rnorm(300, 1, 5))
#' )
#' stat_df_long <- tidyr::pivot_longer(stat_df,
#'     cols = c(Intensity, Height),
#'     names_to = "type", values_to = "value"
#' )
#'
#' print(draw_combo_plot(stat_df_long,
#'     xc = "Feature", yc = "value", fc = "type",
#'     Ylab = "value", comp = list(c(1, 2), c(3, 4), c(1, 3), c(2, 4)), nf = 2
#' ))
#'
#' @export draw_combo_plot

draw_combo_plot <- function(stat_df,
                            xc = "Feature",
                            yc = "Intensity",
                            comp = list(c(1, 2)),
                            Xlab = xc,
                            Ylab = yc,
                            stats = "wilcox.test",
                            fc = xc,
                            Ylim = NULL,
                            title = "",
                            nf = 1) {
    stopifnot(c(xc, yc, fc) %in% colnames(stat_df))
    bbf <- draw_boxplot_by_factor(stat_df, xc = xc, yc = yc, comp = comp,
                                  stats = stats, Xlab = Xlab, Ylab = Ylab,
                                  fc = fc, nf = nf)
    bwo <- draw_boxplot_wo_outlier(stat_df, xc = xc, yc = yc, comp = comp,
                                   stats = stats, Xlab = Xlab, Ylab = Ylab,
                                   fc = fc, nf = nf)
    msb <- draw_mean_se_barplot(stat_df, xc = xc, yc = yc, comp = comp,
                                Xlab = Xlab, Ylab = Ylab, fc = fc, Ylim = Ylim,
                                nf = nf)
    q <- draw_quantile_plot(stat_df, xc = xc, yc = yc, Ylab = Ylab, fc = fc)

    combo <- cowplot::plot_grid(bbf, bwo, q, msb, nrow = 2, axis = "b",
                                rel_widths = c(1, 1), rel_heights = c(1, 1))
    title <- ggdraw() + draw_label(title)
    combo <- plot_grid(title, combo, ncol=1, rel_heights=c(0.1, 1))


    return(combo)
}

#' @title Plot signal profile around start, center, and end of genomic regions
#' @description Plot lines with standard error as the error band, also plots
#' number of regions having non-zero signals
#'
#' @param plot_df a dataframe with column names c(xc, yc, cn, "Interval",
#'  "lower", "upper")
#' @param xc a string denoting column name for values on x-axis
#' @param yc a string denoting column name for numeric data to be plotted
#' @param cn a string denoting column name for grouping
#' @param ext a vector of 4 integers denoting upstream and downstream extension
#'  around start and end, the range of extensions must be within the range of
#'  `xc` of the `plot_df`
#' @param hl a vector of 4 integers defining upstream and downstream boundaries
#'  of the rectangle for start and end
#' @param atitle a string for the title of the plot
#' @param insert a integer denoting the width of the center region
#' @param Ylab a string for y-axis label
#' @param shade logical, indicating whether to place a shaded rectangle around
#'  the point of interest
#' @param stack logical, indicating whether to plot the number of valid
#'  (non-zero) data points in each bin
#'
#' @return a ggplot object
#' @note used by \code{\link{plot_start_end}},
#'  \code{\link{plot_start_end_with_random}}
#' @author Shuye Pu
#'
#' @examples
#' library(dplyr)
#' Reference <- rep(rep(c("Ref1", "Ref2"), each = 100), 2)
#' Query <- rep(c("Query1", "Query2"), each = 200)
#' Position <- rep(seq(-50, 49), 4)
#' Intensity <- rlnorm(400)
#' se <- runif(400)
#' start_df <- data.frame(Intensity, se, Position, Query, Reference) %>%
#'     mutate(lower = Intensity - se, upper = Intensity + se) %>%
#'     mutate(Group = paste(Query, Reference, sep = ":")) %>%
#'     mutate(Location = rep("Start", 400)) %>%
#'     mutate(Interval = sample.int(1000, 400))
#' Intensity <- rlnorm(400, meanlog = 1.5)
#' se <- runif(400)
#' center_df <- data.frame(Intensity, se, Position, Query, Reference) %>%
#'     mutate(lower = Intensity - se, upper = Intensity + se) %>%
#'     mutate(Group = paste(Query, Reference, sep = ":")) %>%
#'     mutate(Location = rep("Center", 400)) %>%
#'     mutate(Interval = sample.int(600, 400))
#' Intensity <- rlnorm(400, meanlog = 2)
#' se <- runif(400)
#' end_df <- data.frame(Intensity, se, Position, Query, Reference) %>%
#'     mutate(lower = Intensity - se, upper = Intensity + se) %>%
#'     mutate(Group = paste(Query, Reference, sep = ":")) %>%
#'     mutate(Location = rep("End", 400)) %>%
#'     mutate(Interval = sample.int(2000, 400))
#'
#' df <- rbind(start_df, center_df, end_df)
#' p <- draw_stacked_profile(df, cn = "Group", shade = TRUE,
#'     ext = c(-50, 50, -50, 50),
#'     hl = c(-20, 20, -25, 25), insert = 100)
#' p
#'
#'
#' @export draw_stacked_profile
#'

draw_stacked_profile <- function(plot_df,
                                 xc = "Position",
                                 yc = "Intensity",
                                 cn = "Query",
                                 ext = c(0, 0, 0, 0),
                                 hl = c(0, 0, 0, 0),
                                 atitle = "title",
                                 insert = 0,
                                 Ylab = "Signal Intensity",
                                 shade = FALSE,
                                 stack = TRUE) {
    stopifnot(c(xc, yc, cn) %in% colnames(plot_df))
    stopifnot(is.numeric(c(ext, hl, insert)))

    ylimits <- c(min(plot_df$lower), max(plot_df$upper))
    ylimits_intervals <- c(0, max(plot_df$Interval))

    aplot_df_start <- plot_df %>%
        filter(grepl("Start", Location))

    ps <- ggplot(aplot_df_start, aes(x = .data[[xc]], y = .data[[yc]],
                                     color = .data[[cn]])) +
        scale_fill_npg() +
        scale_color_npg() +
        geom_line(size = 1.25) +
        ylim(ylimits) +
        geom_vline(xintercept = 0, linetype = "dotted", color = "blue",
                   size = 0.5) +
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = .data[[cn]]),
                    linetype = 0, alpha = 0.3) +
        theme_classic() +
        theme(legend.position = "none") +
        scale_x_continuous(limits = c(ext[1], ext[2])) +
        ylab(Ylab) +
        ggtitle(atitle) +
        theme(
            plot.margin = unit(c(1.2, 0.5, 0, 1.2), "lines"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.line.x = element_blank()
        )
    if (shade) ps <- ps + annotate("rect", xmin = hl[1], xmax = hl[2],
                                   ymin = -Inf, ymax = Inf, fill = "grey",
                                   alpha = 0.3)

    psi <- ggplot(aplot_df_start, aes(x = .data[[xc]], y = Interval,
                                      color = .data[[cn]])) +
        scale_fill_npg() +
        scale_color_npg() +
        geom_line(size = 1.25) +
        ylim(ylimits_intervals) +
        geom_vline(xintercept = 0, linetype = "dotted", color = "blue",
                   size = 0.5) +
        theme_classic() +
        theme(
            plot.margin = unit(c(0, 0.5, 1.2, 1.2), "lines"),
            legend.position = "none"
        ) +
        xlab("5'") +
        scale_x_continuous(limits = c(ext[1], ext[2])) +
        ylab("n")
    if (shade) psi <- psi + annotate("rect", xmin = hl[1], xmax = hl[2],
                                     ymin = -Inf, ymax = Inf, fill = "grey",
                                     alpha = 0.3)
    ps_stack <- plot_grid(ps, psi, ncol = 1, align = "v",
                          rel_heights = c(1, 0.25))

    if (insert > 0) {
        aplot_df_center <- plot_df %>%
            filter(grepl("Center", Location))
        pc <- ggplot(aplot_df_center, aes(x = .data[[xc]], y = .data[[yc]],
                                          color = .data[[cn]])) +
            scale_fill_npg() +
            scale_color_npg() +
            geom_line(size = 1.25) +
            ylim(ylimits) + # geom_point(color="grey30", size=2) +
            geom_vline(xintercept = 0, linetype = "dotted", color = "blue",
                       size = 0.5) +
            geom_ribbon(aes(ymin = lower, ymax = upper, fill = .data[[cn]]),
                        linetype = 0, alpha = 0.3) +
            theme_classic() +
            theme(legend.position = "none") +
            ggtitle("") +
            theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                axis.line = element_blank(),
                plot.margin = unit(c(1.2, 0.5, 0, 0.5), "lines")
            )

        pci <- ggplot(aplot_df_center, aes(x = .data[[xc]], y = Interval,
                                           color = .data[[cn]])) +
            scale_fill_npg() +
            scale_color_npg() +
            geom_line(size = 1.25) +
            ylim(ylimits_intervals) +
            geom_vline(xintercept = 0, linetype = "dotted", color = "blue",
                       size = 0.5) +
            theme_classic() +
            theme(legend.position = "none") +
            xlab("Center") +
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                axis.line.y = element_blank(),
                plot.margin = unit(c(0, 0.5, 1.2, 0.5), "lines")
            )
        pc_stack <- plot_grid(pc, pci, ncol = 1, align = "v",
                              rel_heights = c(1, 0.25))
    } else {
        pc_stack <- NULL
    }

    aplot_df_end <- plot_df %>%
        filter(grepl("End", Location))
    pe <- ggplot(aplot_df_end, aes(x = .data[[xc]], y = .data[[yc]],
                                   color = .data[[cn]])) +
        scale_fill_npg() +
        scale_color_npg() +
        geom_line(size = 1.25) +
        ylim(ylimits) +
        geom_vline(xintercept = 0, linetype = "dotted", color = "blue",
                   size = 0.5) +
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = .data[[cn]]),
                    linetype = 0, alpha = 0.3) +
        theme_classic() +
        theme(
            legend.position = c(0.7, 0.9),
            legend.title = element_blank(),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.spacing.y = unit(1, "mm")
        ) +
        scale_x_continuous(limits = c(ext[3], ext[4])) +
        ggtitle("") +
        guides(fill = "none", color = guide_legend(keyheight = 0.75)) +
        theme(
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            plot.margin = unit(c(1.2, 1.2, 0, 0.5), "lines")
        )

    if (shade) pe <- pe + annotate("rect", xmin = hl[3], xmax = hl[4],
                                   ymin = -Inf, ymax = Inf, fill = "grey",
                                   alpha = 0.3)

    pei <- ggplot(aplot_df_end, aes(x = .data[[xc]], y = Interval,
                                    color = .data[[cn]])) +
        scale_fill_npg() +
        scale_color_npg() +
        geom_line(size = 1.25) +
        ylim(ylimits_intervals) +
        theme_classic() +
        theme(legend.position = "none") +
        xlab("3'") +
        geom_vline(xintercept = 0, linetype = "dotted", color = "blue",
                   size = 0.5) +
        scale_x_continuous(limits = c(ext[3], ext[4])) +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = unit(c(0, 1.2, 1.2, 0.5), "lines")
        )
    if (shade) pei <- pei + annotate("rect", xmin = hl[3], xmax = hl[4],
                                     ymin = -Inf, ymax = Inf, fill = "grey",
                                     alpha = 0.3)

    pe_stack <- plot_grid(pe, pei, ncol = 1, align = "v",
                          rel_heights = c(1, 0.25))

    if (insert > 0) {
        if (stack) {
            p <- plot_grid(ps_stack, pc_stack, pe_stack, ncol = 3, align = "h",
                           rel_widths = c(1, 0.5, 0.9))
        } else {
            p <- plot_grid(ps, pc, pe, ncol = 3, align = "h",
                           rel_widths = c(1, 0.5, 0.9))
        }
    } else {
        if (stack) {
            p <- plot_grid(ps_stack, pe_stack, ncol = 2, align = "h",
                           rel_widths = c(1, 0.9))
        } else {
            p <- plot_grid(ps, pe, ncol = 2, align = "h",
                           rel_widths = c(1, 0.9))
        }
    }
    return(p)
}

#' @title Plot two-sets Venn diagram
#'
#' @description This is a helper function for Venn diagram plot. A Venn diagram
#' is plotted as output. For GRanges, as A overlap B may not be the same as B
#' overlap A, the order of GRanges in a list matters, certain order may produce
#' an error.
#' @param apair a list of two vectors
#' @param overlap_fun the name of the function that defines overlap, depending
#'  on the type of object in the vectors. For GRanges, use
#'  \code{\link{filter_by_overlaps_stranded}} or
#'  \code{\link{filter_by_nonoverlaps_stranded}}, for gene names, use intersect.
#' @param title main title of the figure
#'
#' @return a VennDiagram object
#' @author Shuye Pu
#'
#' @examples
#' test_list <- list(A = c(1, 2, 3, 4, 5), B = c(4, 5, 7))
#' overlap_pair(test_list, intersect, title = "test")
#'
#' ## GRanges overlap
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
#' overlap_pair(
#'     list(query = query, subject = subject),
#'     filter_by_overlaps_stranded
#' )
#'
#' @export overlap_pair

overlap_pair <- function(apair, overlap_fun, title = NULL) {
    sizes <- vapply(apair, length, numeric(1))
    overlap <- length(Reduce(overlap_fun, apair))
    jaccard <- round(overlap / (sum(sizes) - overlap), digits = 5)
    venn.plot <- draw.pairwise.venn(sizes[1], sizes[2], overlap,
        category = names(apair),
        lty = rep("blank", 2), fill = c("#0020C2", "#64E986"),
        cat.just = rep(list(c(0.5, 0)), 2), cex = rep(2, 3), cat.pos = c(0, 0)
    )

    grid.draw(venn.plot)
    grid.text(paste("Jaccard:", jaccard), unit(0.8, "npc"), unit(0.025, "npc"),
              draw = TRUE)
    grid.text(title, unit(0.25, "npc"), unit(0.975, "npc"), draw = TRUE)
    grid.newpage()

    invisible(venn.plot)
}

#' @title Plot three-sets Venn diagram
#'
#' @description This is a helper function for Venn diagram plot. A Venn diagram
#' is plotted as output. For GRanges, as A overlap B may not be the same as B
#' overlap A, the order of GRanges in a list matters, certain order may produce
#' an error.
#' @param atriple a list of three vectors
#' @param overlap_fun the name of the function that defines overlap, depending
#'  on the type of object in the vectors. For GRanges, use
#'  \code{\link{filter_by_overlaps_stranded}} or
#'  \code{\link{filter_by_nonoverlaps_stranded}}, for gene names, use intersect.
#' @param title main title of the figure
#' @return a VennDiagram object
#' @author Shuye Pu
#'
#' @examples
#'
#' test_list <- list(A = c(1, 2, 3, 4, 5), B = c(4, 5, 7), C = c(1, 3))
#' overlap_triple(test_list, intersect, title = "test")
#'
#' ## GRanges overlap
#' query <- GRanges("chr19",
#'     IRanges(rep(c(10, 15), 2), width = c(1, 20, 40, 50)),
#'     strand = c("+", "+", "-", "-")
#' )
#'
#' subject1 <- GRanges("chr19",
#'     IRanges(rep(c(13, 150), 2), width = c(10, 14, 20, 28)),
#'     strand = c("+", "-", "-", "+")
#' )
#'
#' subject2 <- GRanges("chr19",
#'     IRanges(rep(c(13, 50), 2), width = c(10, 14, 20, 21)),
#'     strand = c("+", "-", "-", "+")
#' )
#'
#' overlap_triple(
#'     list(subject1 = subject1, subject2 = subject2, query = query),
#'     filter_by_overlaps_stranded
#' )
#'
#' @export overlap_triple

overlap_triple <- function(atriple, overlap_fun, title = NULL) {
    ## sort the gr by decreasing size to avoid n13 < n123
    sizes <- sort(vapply(atriple, length, numeric(1)), decreasing = TRUE)
    atriple <- atriple[names(sizes)]

    overlap12 <- length(Reduce(overlap_fun, atriple[c(1, 2)]))
    overlap13 <- length(Reduce(overlap_fun, atriple[c(1, 3)]))
    overlap23 <- length(Reduce(overlap_fun, atriple[c(2, 3)]))
    overlap123 <- length(Reduce(overlap_fun, atriple))

    # check consistency if bed overlaps produce inconsistent values
    if(overlap123 > min(overlap12, overlap13, overlap23)){
      overlap123 <- min(overlap12, overlap13, overlap23)
    }

    venn.plot <- draw.triple.venn(sizes[1], sizes[2], sizes[3], overlap12,
                                  overlap23, overlap13, overlap123,
                                  category = names(atriple),
                                  lty = rep("blank", 3),
                                  fill = c("#0020C2", "#64E986", "#990012"),
                                  cat.just = rep(list(c(0.5, 0)), 3),
                                  cex = rep(2, 7), cat.pos = c(0, 0, 180)
    )

    grid.draw(venn.plot)
    grid.text(title, unit(0.25, "npc"), unit(0.975, "npc"), draw = TRUE)
    grid.newpage()

    invisible(venn.plot)
}

#' @title Plot four-sets Venn diagram
#'
#' @description This is a helper function for Venn diagram plot. A Venn diagram
#' is plotted as output. For GRanges, as A overlap B may not be the same as B
#' overlap A, the order of GRanges in a list matters, certain order may produce
#' an error.
#' @param aquad a list of four vectors
#' @param overlap_fun the name of the function that defines overlap, depending
#'  on the type of object in the vectors. For GRanges, use
#'  \code{\link{filter_by_overlaps_stranded}} or
#'  \code{\link{filter_by_nonoverlaps_stranded}}, for gene names, use intersect.
#' @param title main title of the figure
#'
#' @return a VennDiagram object
#' @author Shuye Pu
#'
#' @examples
#'
#' test_list <- list(A = c(1, 2, 3, 4, 5), B = c(4, 5, 7), C = c(1, 3), D = 6)
#' overlap_quad(test_list, intersect)
#'
#' ## GRanges overlap
#' query1 <- GRanges("chr19",
#'     IRanges(rep(c(10, 15), 2), width = c(1, 20, 40, 50)),
#'     strand = c("+", "+", "-", "-")
#' )
#'
#' query2 <- GRanges("chr19",
#'     IRanges(rep(c(1, 15), 2), width = c(1, 20, 40, 50)),
#'     strand = c("+", "+", "-", "-")
#' )
#'
#' subject1 <- GRanges("chr19",
#'     IRanges(rep(c(13, 150), 2), width = c(10, 14, 20, 28)),
#'     strand = c("+", "-", "-", "+")
#' )
#'
#' subject2 <- GRanges("chr19",
#'     IRanges(rep(c(13, 50), 2), width = c(10, 14, 20, 21)),
#'     strand = c("+", "-", "-", "+")
#' )
#'
#' overlap_quad(list(
#'     subject1 = subject1, subject2 = subject2, query1 = query1,
#'     query2 = query2
#' ), filter_by_overlaps_stranded)
#'
#' @export overlap_quad
#'
overlap_quad <- function(aquad, overlap_fun, title = NULL) {
    ## sort the gr by decreasing size to avoid n13 < n123
    sizes <- sort(vapply(aquad, length, numeric(1)), decreasing = TRUE)
    aquad <- aquad[names(sizes)]

    overlap12 <- length(Reduce(overlap_fun, aquad[c(1, 2)]))
    overlap13 <- length(Reduce(overlap_fun, aquad[c(1, 3)]))
    overlap14 <- length(Reduce(overlap_fun, aquad[c(1, 4)]))
    overlap23 <- length(Reduce(overlap_fun, aquad[c(2, 3)]))
    overlap24 <- length(Reduce(overlap_fun, aquad[c(2, 4)]))
    overlap34 <- length(Reduce(overlap_fun, aquad[c(3, 4)]))
    overlap123 <- length(Reduce(overlap_fun, aquad[c(1, 2, 3)]))
    overlap124 <- length(Reduce(overlap_fun, aquad[c(1, 2, 4)]))
    overlap134 <- length(Reduce(overlap_fun, aquad[c(1, 3, 4)]))
    overlap234 <- length(Reduce(overlap_fun, aquad[c(2, 3, 4)]))
    overlap1234 <- length(Reduce(overlap_fun, aquad))

    # check consistency if bed overlaps produce inconsistent values
    if(overlap123 > min(overlap12, overlap13, overlap23)){
      overlap123 <- min(overlap12, overlap13, overlap23)
    }
    if(overlap124 > min(overlap12, overlap14, overlap24)){
      overlap124 <- min(overlap12, overlap14, overlap24)
    }
    if(overlap134 > min(overlap13, overlap14, overlap34)){
      overlap134 <- min(overlap13, overlap14, overlap34)
    }
    if(overlap234 > min(overlap23, overlap24, overlap34)){
      overlap234 <- min(overlap23, overlap24, overlap34)
    }

    if(overlap1234 > min(overlap123, overlap124, overlap134, overlap234)){
      overlap1234 <- min(overlap123, overlap124, overlap134, overlap234)
    }

    venn.plot <- draw.quad.venn(sizes[1], sizes[2], sizes[3], sizes[4],
                                overlap12, overlap13, overlap14, overlap23,
                                overlap24, overlap34, overlap123, overlap124,
                                overlap134, overlap234, overlap1234,
                                category = names(aquad),
                                lty = rep("blank", 4),
                                fill = c("#0020C2", "#64E986", "#990012",
                                         "#c6dcff"),
                                cat.just = rep(list(c(0.5, 0)), 4),
                                cex = rep(2, 15), cat.pos = c(0, 0, 0, 0)
    )

    grid.draw(venn.plot)
    grid.text(title, unit(0.25, "npc"), unit(0.975, "npc"), draw = TRUE)
    grid.newpage()
    invisible(venn.plot)
}
