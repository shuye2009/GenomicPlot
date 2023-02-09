
#' @title Display matrix as a heatmap
#' @description Make a complex heatmap with column annotations
#'
#' @param fullmatrix a numeric matrix
#' @param dataName the nature of the numeric data
#' @param labels_col a vector make column annotation
#' @param levels_col factor levels for labels_col, specifying the order of labels_col
#' @param ranking method for ranking the rows of the input matrix
#' @param verbose logical, whether to output the input matrix for inspection
#' @return NULL
#'
#' @author Shuye Pu
#'
#' @examples
#' fullmatrix <- matrix(rnorm(10000), ncol=100)
#' labels_col <- as.character(seq(1:100))
#' levels_col <- c("start", "center", "end")
#' names(labels_col) <- rep(levels_col, c(15, 60, 25))
#'
#' draw_matrix_heatmap(fullmatrix, dataName="test", labels_col, levels_col, ranking="Sum")
#'
#' @export draw_matrix_heatmap
#'
#'

draw_matrix_heatmap <- function(fullmatrix, dataName="geneData", labels_col=NULL, levels_col=NULL, ranking="Sum", verbose=FALSE){
   
   message("drawing heatmap")
   #inspect_matrix(fullmatrix, verbose)
   ## reduce the size by removing rows with all 0s
   #fullmatrix[is.na(fullmatrix)] <- 0
   #all0 <- apply(fullmatrix, 1, function(x)all(x == sum(x)/length(x)))
   #fullmatrix <- data.matrix(fullmatrix[!all0,])
   #count_all0_rows <- sum(all0)
   #if(count_all0_rows > 0){
   #   message(count_all0_rows, " all 0 rows are removed!")
   #}

   if(verbose){
      vdataName <- gsub(":", "_", dataName, fixed=TRUE)
      write.table(fullmatrix, paste(vdataName, "_matrix.tab", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
   }

   if(is.null(labels_col)){
      labels_col <- seq(1, ncol(fullmatrix))
      names(labels_col) <- rep("column", ncol(fullmatrix))
      levels_col <- "column"
   }
   colnames(fullmatrix) <- labels_col

   fullmatrix <- rank_rows(fullmatrix, ranking)

   features <- factor(names(labels_col), levels=levels_col)
   cols <- ggsci::pal_npg()(length(levels(features)))
   names(cols) <- levels(features)
   mycols <- cols[features]
   names(mycols) <- features

   ha <- HeatmapAnnotation(df = data.frame(feature = features), col=list(feature=mycols), which="column", show_legend=FALSE, annotation_label = "")
   #y <- matrix(as.vector(fullmatrix), ncol=1)
   if(verbose) print(quantile(fullmatrix, c(seq(0.9, 1, 0.005)), na.rm=TRUE))
   ranges <- quantile(fullmatrix, c(0.025, 0.975), na.rm=TRUE)
   if(ranges[1] == ranges[2]){
      message("97.5% of values are not unique, heatmap may not show signals effectively")
      
      ranges <- quantile(fullmatrix, c(0, 0.995), na.rm=TRUE) ## Need to have a better way for determine the upper bound
      ranges[2] <- ranges[2]*2
   }

   h <- Heatmap(fullmatrix,
                name = unlist(strsplit(dataName, split=":", fixed=TRUE))[1],
                col = colorRamp2(ranges, viridis(2)),
                bottom_annotation = ha,
                heatmap_legend_param=list(legend_direction = "horizontal"),
                show_row_names = FALSE,
                show_column_names = FALSE,
                show_row_dend = FALSE,
                cluster_columns=FALSE,
                cluster_rows=FALSE,
                column_split=features,
                column_gap = unit(0, "mm"),
                column_title="%s",
                column_title_gp = gpar(fontsize = 10, fontface = "plain"),
                column_title_side = "bottom")

   return(h)
}


#' @title Plot genomic region landmark indicator
#' @description Plot a gene centered polygon for demarcating gene and its upstream and downstream regions
#'
#' @param featureNames a string vector giving names of sub-regions
#' @param vx a vector on integers denoting the x coordinates of start of each sub-region
#' @param xmax an integer denoting the left most boundary
#' @return a ggplot object
#' @note used by \code{plot_3parts_metagene}, \code{plot_5parts_metagene}, \code{plot_reference_region}
#'
#' @author Shuye Pu
#'
#' @export draw_region_landmark


draw_region_landmark <- function(featureNames, vx, xmax){
   nfeatures <- length(featureNames)
   if(nfeatures == 5){
      values <- data.frame(id=featureNames, value=c(1.75, 1.5, 1.25, 1.5, 1.75))
      positions <- data.frame(
         id = rep(featureNames, each = 4),
         x = c(vx[2], vx[1], vx[1], vx[2], vx[3], vx[2], vx[2], vx[3], vx[4], vx[3], vx[3], vx[4], vx[5], vx[4], vx[4], vx[5], xmax, vx[5], vx[5], xmax),
         y = c(3, 3, 4, 4, 2.5, 2.5, 4.5, 4.5, 2, 2, 5, 5, 2.5, 2.5, 4.5, 4.5, 3, 3, 4, 4) - 2
      )
   }else if(nfeatures == 3){
      values <- data.frame(id=featureNames, value=c(1.25, 1.75, 1.25))
      positions <- data.frame(
         id = rep(featureNames, each = 4),
         x = c(vx[2], vx[1], vx[1], vx[2], vx[3], vx[2], vx[2], vx[3], xmax, vx[3], vx[3], xmax),
         y = c(3, 3, 4, 4, 2.5, 2.5, 4.5, 4.5, 3, 3, 4, 4) - 2
      )
   }else{
      stop("Number of feautre names must be 3 or 5! other numbers are not supported at this point.")
   }

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
            plot.background=element_blank(),
            plot.margin=unit(c(0, 0, 0, 0), "mm"))

   return(pp)
}

#' @title Plot genomic region names
#' @description Plot sub-region labels under the landmark
#'
#' @param featureNames a string vector giving names of sub-regions
#' @param scaled_bins a vector on integers denoting the length of each sub-region
#' @param xmax an integer denoting the left most boundary
#' @return a ggplot object
#' @note used by \code{plot_3parts_metagene}, \code{plot_5parts_metagene}, \code{plot_reference_region}
#'
#' @author Shuye Pu
#'
#' @export draw_region_name

draw_region_name <- function(featureNames, scaled_bins, xmax){

   annotx <- scaled_bins/2
   for(i in 2:length(scaled_bins)){
      annotx[i] <- annotx[i] + sum(scaled_bins[1:(i-1)])
   }

   annot <- data.frame(
      fn <- featureNames,
      x <- annotx,
      y <- 0
   )
   ppp <- ggplot(annot, aes(x = x, y = y, label=fn)) +
      geom_text(size=3, check_overlap=TRUE) +
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
            plot.background=element_blank(),
            plot.margin=unit(c(1, 0, 1, 0), "mm"))

   return(ppp)
}

#' @title Plot signal profile in genomic regions
#' @description Plot lines with standard error as the error band
#'
#' @param plot_df a dataframe with column names c(xc, yc, cn, "lower", "upper")
#' @param xc a string denoting column name for values on x-axis
#' @param yc a string denoting column name for numeric data to be plotted
#' @param vx a vector on integers denoting the x coordinates of start of each sub-region
#' @param cn column name for grouping
#' @param Ylab a string for Y-axis label
#' @return a ggplot object
#' @note used by \code{plot_3parts_metagene}, \code{plot_5parts_metagene}, \code{plot_reference_region}
#'
#' @author Shuye Pu
#'
#' @export draw_region_profile
#'
#'

draw_region_profile <- function(plot_df, xc="Position", yc="Intensity", cn="Query", Ylab="Signal Intensity", vx){

   p <- ggplot(plot_df, aes(x=.data[[xc]], y=.data[[yc]], color=.data[[cn]])) + scale_fill_npg() + scale_color_npg() +
      geom_line(size=2) + #geom_point(color="grey30", size=2) +
      geom_vline(xintercept = vx[2:length(vx)], linetype="dotted", color = "blue", size=0.5) +
      geom_ribbon(aes(ymin=lower, ymax=upper, fill=.data[[cn]]), linetype=0, alpha=0.3) +
      theme_classic() + ylab(Ylab) +
      theme(legend.position="top",
            axis.title.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.x=element_blank(),
            plot.margin=unit(c(1, 1, 0, 1), "lines"))

   return(p)
}

#' @title Plot signal profile around genomic loci
#' @description Plot lines with standard error as the error band
#'
#' @param plot_df a dataframe with column names c(xc, yc, cn, "lower", "upper")
#' @param xc a string denoting column name for values on x-axis
#' @param yc a string denoting column name for numeric data to be plotted
#' @param cn a string denoting column name for sample grouping, like 'Query' or 'Reference'
#' @param sn a string denoting column name for the subject of sample grouping, if 'cn' is 'Query', then 'sn' will be 'Reference'
#' @param Xlab a string for x-axis label
#' @param Ylab a string for y-axis label
#' @param shade logical indicating whether to place a shaded rectangle around the loci
#' @param hl a vector of two integers defining upstream and downstream boundaries of the rectangle
#'
#' @return a ggplot object
#' @note used by \code{plot_reference_locus}, \code{plot_reference_locus_with_random}
#' @author Shuye Pu
#' @export draw_locus_profile
#'

draw_locus_profile <- function(plot_df, xc="Position", yc="Intensity", cn="Query", sn="Reference", Xlab="Center", Ylab="Signal Intensity", shade=FALSE, hl=c(0,0)){

   p <- ggplot(plot_df, aes(x=.data[[xc]], y=.data[[yc]], color=.data[[cn]])) +
      scale_fill_npg() + scale_color_npg() +
      geom_line(size=2) + #geom_point(color="grey30", size=1) +
      geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
      geom_ribbon(aes(ymin=lower, ymax=upper, fill=.data[[cn]]), linetype=0, alpha=0.3) +
      theme_classic() + xlab(Xlab) + ylab(Ylab) +
      theme(legend.position="top", legend.title=element_blank(),
            axis.text = element_text(face="plain", size=14),
            axis.title = element_text(face="bold", size=16)) +
      ggtitle(unique(plot_df[[sn]]))

   if(shade) p <- p + annotate("rect", xmin=hl[1], xmax=hl[2], ymin=-Inf, ymax=Inf, fill="grey", color="grey", alpha=0.3)

   return(p)
}

#' @title Plot boxplot with log scale y-axis
#' @description Plot boxplot with linear or log scale y-axis, with p-value significance levels displayed
#'
#' @param stat_df a dataframe with column names c(xc, yc)
#' @param xc a string denoting column name for grouping
#' @param fc a string denoting column name for sub-grouping
#' @param yc a string denoting column name for numeric data to be plotted
#' @param Ylab a string for y-axis label
#' @param comp a list of vectors denoting pair-wise comparisons to be performed between groups
#' @param stats the name of pair-wise statistical tests, like t.test or wilcox.test
#' @param logy character string in c(linear, logy, wo), "wo" = without outliers
#' @param nf a integr normalizing factor for correct count of observations when the data table is produced by pivot_longer, equals to the number of columns in cols argument of pivot_longer
#'
#' @return a ggplot object
#' @note used by \code{plot_reference_locus}, \code{plot_reference_locus_with_random}
#' @author Shuye Pu
#'
#' @export draw_boxplot_logy
#'

draw_boxplot_logy <- function(stat_df, xc="Feature", yc="Intensity", fc=xc, comp=list(c(1,2)), stats="wilcox.test", Xlab=xc, Ylab=yc, logy="linear", nf=1){
   xlabs <- paste(levels(as.factor(stat_df[[xc]])),"\n(",table(stat_df[[xc]])/nf,")",sep="")
   ypos <- rep(max(stat_df[[yc]]), length(comp))*seq(1, 1+(length(comp)-1)*0.1, 0.1)
   outlier.shape = 19
   if(logy == "wo"){
      fomu <- as.formula(paste(yc, "~", xc))
      bp <- boxplot(fomu, stat_df, plot=FALSE)
      lim <- c(min(bp$stats), max(bp$stats) + abs(max(bp$stats))*0.25)
      ypos <- rep(lim[2], length(comp))*seq(1, 1+(length(comp)-1)*0.1, 0.1)
      lim[2] <- max(ypos)*1.25
      outlier.shape = NA
   }
   
   if(fc == xc){
      p <- ggplot(stat_df, aes(x=.data[[xc]], y=.data[[yc]], fill=.data[[fc]])) +
         geom_violin(width=0.5) +
         geom_boxplot(width=0.2, outlier.shape = outlier.shape) +
         scale_fill_npg() +
         scale_color_npg() +
         theme_classic() +
         theme(axis.text = element_text(face="plain", size=14, color="black"),
               #axis.title.x = element_blank(),
               axis.title = element_text(face="bold", size=16, color="black"),
               legend.position = "bottom") +
         labs(y=Ylab, x=Xlab) +
         geom_signif(comparisons = comp, test=stats, map_signif_level=TRUE, y_position=ypos)+
         scale_x_discrete(labels = xlabs) 
   }else{
      mid <- function(v){
         m <- rep(0, (length(v)/2))
         for(i in 1:length(m)){
            m[i] <- (v[i*2-1]+v[i*2])/2
         }
         return(m)
      }
      stat_df <- stat_df %>%
         mutate(x2 = as.integer(interaction(.data[[fc]], .data[[xc]])))
      p <- ggplot(stat_df, aes(x=x2, y=.data[[yc]], group=x2, fill=.data[[fc]])) +
         geom_violin(width=0.5) +
         geom_boxplot(width=0.2, outlier.shape = outlier.shape) +
         scale_fill_npg() +
         scale_color_npg() +
         theme_classic() +
         theme(axis.text = element_text(face="plain", size=14, color="black"),
               #axis.title.x = element_blank(),
               axis.title = element_text(face="bold", size=16, color="black"),
               legend.position = "bottom") +
         labs(y=Ylab, x=Xlab) +
         geom_signif(comparisons = comp, test=stats, map_signif_level=TRUE, y_position=ypos) +
         scale_x_continuous(breaks = mid(sort(unique(stat_df$x2))), labels = xlabs) 
   }
   
   if(logy == "logy"){
      p <- p + scale_y_continuous(trans='log10', labels = scales::scientific)
      message("log scale")
   }else if(logy == "linear"){
      message("linear scale")
   }else if(logy == "wo"){
      p <- p + coord_cartesian(ylim=lim)
      message("without outliers")
   }else{
      stop("plot type is not supported, use one of c(linear, logy, wo) only")
   }
   
   return(p)
}

#' @title Plot boxplot without outliers
#' @description Plot boxplot without outliers, with p-value significance levels displayed
#'
#' @param stat_df a dataframe with column names c(xc, yc)
#' @param xc a string denoting column name for grouping
#' @param fc a string denoting column name for sub-grouping
#' @param yc a string denoting column name for numeric data to be plotted
#' @param Ylab a string for y-axis label
#' @param comp a list of vectors denoting pair-wise comparisons to be performed between groups
#' @param stats the name of pair-wise statistical tests, like t.test or wilcox.test
#' @param nf a integr normalizing factor for correct count of observations when the data table is produced by pivot_longer, equals to the number of columns in cols argument of pivot_longer
#'
#' @return a ggplot object
#' @export draw_boxplot_wo_outlier
#' @examples 
#' stat_df <- data.frame(Feature=rep(c("A", "B"), c(20, 30)), Intensity=rlnorm(50))
#' p <- draw_boxplot_wo_outlier(stat_df, xc="Feature", yc="Intensity", Ylab="Signal Intensity")
#' p
#' 
draw_boxplot_wo_outlier <- function(stat_df, xc="Feature", yc="Intensity", fc=xc, comp=list(c(1,2)), stats="wilcox.test", Xlab=xc, Ylab=yc, nf=1){
   xlabs <- paste(levels(as.factor(stat_df[[xc]])),"\n(",table(stat_df[[xc]])/nf,")",sep="")
   fomu <- as.formula(paste(yc, "~", xc))
   bp <- boxplot(fomu, stat_df, plot=FALSE)
   lim <- c(min(bp$stats) - abs(min(bp$stats))*0.25, max(bp$stats) + abs(max(bp$stats))*0.25)
   ypos <- rep(lim[2], length(comp))*seq(1, 1+(length(comp)-1)*0.1, 0.1)
   lim[2] <- max(ypos)+0.05
   violin <- FALSE
   if(max(stat_df[[yc]])/lim[2] < 10) violin <- TRUE
   
   if(fc == xc){
      p <- ggplot(stat_df, aes(x=.data[[xc]], y=.data[[yc]], fill=.data[[fc]])) +
         geom_boxplot(width=0.2, outlier.shape = NA) +
         scale_fill_npg() +
         scale_color_npg() +
         theme_classic() +
         theme(axis.text = element_text(face="plain", size=14, color="black"),
               #axis.title.x = element_blank(),
               axis.title = element_text(face="bold", size=16, color="black"),
               legend.position = "bottom") +
         labs(y=Ylab, x=Xlab) +
         geom_signif(comparisons=comp, test=stats, map_signif_level=TRUE, y_position=ypos, tip_length=(lim[2]-lim[1])*0.03) +
         coord_cartesian(ylim=lim) +
         scale_x_discrete(labels = xlabs) 
      
   }else{
      mid <- function(v){
         m <- rep(0, (length(v)/2))
         for(i in 1:length(m)){
            m[i] <- (v[i*2-1]+v[i*2])/2
         }
         return(m)
      }
      stat_df <- stat_df %>%
         mutate(x2 = as.integer(interaction(.data[[fc]], .data[[xc]])))
      p <- ggplot(stat_df, aes(x=x2, y=.data[[yc]], group=x2, fill=.data[[fc]])) +
         geom_boxplot(width=0.2, outlier.shape = NA) +
         scale_fill_npg() +
         scale_color_npg() +
         theme_classic() +
         theme(axis.text = element_text(face="plain", size=14, color="black"),
               #axis.title.x = element_blank(),
               axis.title = element_text(face="bold", size=16, color="black"),
               legend.position = "bottom") +
         labs(y=Ylab, x=Xlab) +
         geom_signif(comparisons=comp, test=stats, map_signif_level=TRUE, y_position=ypos, tip_length=(lim[2]-lim[1])*0.03) +
         scale_x_continuous(breaks = mid(sort(unique(stat_df$x2))), labels = xlabs) +
         coord_cartesian(ylim=lim) 
         
   }
   
   if(violin) p <- p + geom_violin(width=0.5)
      
   return(p)
}
#' @title Plot barplot for mean with standard error bars
#' @description Plot barplot for mean with standard error bars, no p-value significance levels are displayed, but ANOVA p-value is provided as tag and TukeyHSD test are displayed as caption.
#'
#' @param stat_df a dataframe with column names c(xc, yc)
#' @param xc a string denoting column name for grouping
#' @param yc a string denoting column name for numeric data to be plotted
#' @param comp a list of vectors denoting pair-wise comparisons to be performed between groups
#' @param Ylab a string for y-axis label
#'
#' @return a ggplot object
#' @note used by \code{plot_reference_locus}, \code{plot_reference_locus_with_random}
#' @author Shuye Pu
#'
#' @export draw_mean_se_barplot
#'

draw_mean_se_barplot <- function(stat_df, xc="Feature", yc="Intensity", comp=list(c(1,2)), Xlab=xc, Ylab=yc){
   stat_df[[xc]] <- as.factor(stat_df[[xc]])
   
   stats <- aov_TukeyHSD(stat_df, xc, yc)

   means_se <- stat_df %>%
      group_by(.data[[xc]]) %>%
      summarize(mean_Intensity=mean(.data[[yc]]),
                sd_Intensity=sd(.data[[yc]]),
                N_Intensity=length(.data[[yc]]),
                se=sd_Intensity/sqrt(N_Intensity),
                upper_limit=mean_Intensity+se,
                lower_limit=mean_Intensity-se
      )
   means_se <- means_se %>%
      mutate(labelx=paste0(.data[[xc]], "\n(", N_Intensity, ")")) ## now .data is means_se, not stat_df anymore
   levels(means_se[[xc]]) <- levels(stat_df[[xc]])

   p <- ggplot(means_se, aes(x=.data[[xc]], y=mean_Intensity, fill=.data[[xc]])) +
      scale_fill_npg() + scale_color_npg() +
      geom_col(stat="identity") +
      geom_errorbar(aes(ymin=lower_limit, ymax=upper_limit), position=position_dodge(width = 0.2), width=0.2) +
      theme_classic() +
      theme(axis.title = element_text(face="bold", size=14, color="black", vjust=0.25),
            axis.text = element_text(face="plain", size=16, color="black"),
            legend.position = "none") +
      labs(y=Ylab, x=Xlab, caption="post hoc TukeyHSD test") +
      scale_x_discrete(labels = means_se$labelx) +
      ggtitle(label="Mean + SE" , subtitle=paste("ANOVA p-value =",format(stats$ANOVA, digits=3)))
   
   stats$HSD[,1:3] <- round(stats$HSD[,1:3], digits=3)
   stats$HSD[,4] <- format(stats$HSD[,4], digits=3)
   
   comp_row <- sapply(comp, function(x){arow <- paste0(levels(stat_df[[xc]])[x[2]], "-", levels(stat_df[[xc]])[x[1]])})
   stats_selected <- as.data.frame(stats$HSD) %>%
      filter(row.names(stats$HSD) %in% comp_row)
   ptable <- ggtexttable(stats_selected)
   
   outp <- cowplot::plot_grid(p, ptable, ncol=1, rel_heights = c(3,1))

   return(outp)
}

#' @title Plot cumulative sum or quantile over rank
#' @description Plot cumulative sum over rank as line plot, both cumulative sum and rank are scaled between 0 and 1. This is the same as the fingerprint plot of the deepTools. Quantiles can also be used as y-axis, and values can also be used as x-axis. If the curve is skewed toward ends, the x-axis is truncated for better visualization.
#'
#' @param stat_df a dataframe with column names c(xc, yc)
#' @param xc a string denoting column name for grouping
#' @param yc a string denoting column name for numeric data to be plotted
#' @param Ylab a string for y-axis label
#' @param ecdf logical, indicating using quantile instead of cumulative sum as y-axis
#' @param rank logical, indicating using rank of values instead of value itself as x-axis
#' 
#' @return a ggplot object
#' @note used by \code{plot_reference_locus}, \code{plot_reference_locus_with_random}
#' @author Shuye Pu
#'
#' @export draw_rank_plot
#'
#' @examples 
#' stat_df <- data.frame(Feature=rep(c("A", "B"), c(20, 30)), Intensity=c(rlnorm(20, 5, 5), rlnorm(30, 1, 5)))
#' stat_df1 <- data.frame(Feature=rep(c("A", "B"), c(20, 30)), Intensity=c(rnorm(20, 5, 5), rnorm(30, 1, 5)))
#' 
#' for(e in c(TRUE, FALSE)){
#'    for (r in c(TRUE, FALSE)){
#'       draw_rank_plot(stat_df, xc="Feature", yc="Intensity", Ylab="Signal Intensity", ecdf=e, rank=r)
#'       draw_rank_plot(stat_df1, xc="Feature", yc="Intensity", Ylab="Signal Intensity", ecdf=e, rank=r)
#'    }
#' }
#'
#' 
draw_rank_plot <- function(stat_df, xc="Feature", yc="Intensity", Ylab="Signal Intensity", ecdf=TRUE, rank=FALSE){
   
   if(ecdf){
      long_df <- stat_df %>%
         group_by(.data[[xc]]) %>%
         arrange(.data[[yc]]) %>%
         mutate(Rank=rank(.data[[yc]])) %>%
         mutate(Fraction=ecdf(.data[[yc]])(.data[[yc]]), Rank=Rank/max(Rank))
      
      if(rank){
         p <- ggplot(data=long_df, aes(x=Rank, y=Fraction, color=.data[[xc]])) +
            scale_color_npg() +
            geom_line(size=2) +
            labs(x=paste0("Rank (",Ylab,")"), y="Cumulative fraction") +
            theme_classic() +
            theme(legend.position="top", 
                  legend.title=element_blank(),
                  axis.text=element_text(angle=0, size=14, vjust=0),
                  axis.title=element_text(face="bold", color="black", size=16, vjust=0.25) 
            ) +
            ggtitle(paste("Cumulative fraction of ", Ylab))
      }else{
         
         p <- ggplot(data=long_df, aes(x=.data[[yc]], y=Fraction, color=.data[[xc]])) +
            scale_color_npg() +
            geom_line(size=2) +
            labs(x=Ylab, y="Cumulative fraction") +
            theme_classic() +
            theme(legend.position="top", 
                  legend.title=element_blank(),
                  axis.text=element_text(angle=0, size=14, vjust=0),
                  axis.title=element_text(face="bold", color="black", size=16, vjust=0.25) 
            ) +
            ggtitle(paste("Cumulative fraction of ", Ylab))
         
         if(max(stat_df[[yc]])/quantile(stat_df[[yc]], 0.9) > 100 ){
            p <- p + coord_cartesian(xlim=quantile(stat_df[[yc]], c(0, 0.9)))
         }
      }
      
   }else{
      long_df <- stat_df %>%
         group_by(.data[[xc]]) %>%
         arrange(.data[[yc]]) %>%
         mutate(cumCount=cumsum(.data[[yc]])) %>%
         mutate(Rank=rank(cumCount)) %>%
         mutate(Fraction=cumCount/max(cumCount), Rank=Rank/max(Rank))
      if(rank){
         p <- ggplot(data=long_df, aes(x=Rank, y=Fraction, color=.data[[xc]])) +
            scale_color_npg() +
            geom_line(size=2) +
            labs(x=paste0("Rank (",Ylab,")"), y="Cumulative sum fraction") +
            theme_classic() +
            theme(legend.position="top", 
                  legend.title=element_blank(),
                  axis.text=element_text(angle=0, size=14, vjust=0),
                  axis.title=element_text(face="bold", color="black", size=16, vjust=0.25)
            ) +
            ggtitle(paste("Cumulative sum fraction of", Ylab)) 
         
         if(max(stat_df[[yc]])/quantile(stat_df[[yc]], 0.9) > 100 ){
            p <- p + coord_cartesian(xlim=c(0.9, 1))
         }
      }else{
         p <- ggplot(data=long_df, aes(x=.data[[yc]], y=Fraction, color=.data[[xc]])) +
            scale_color_npg() +
            geom_line(size=2) +
            labs(x=Ylab, y="Cumulative sum fraction") +
            theme_classic() +
            theme(legend.position="top", 
                  legend.title=element_blank(),
                  axis.text=element_text(angle=0, size=14, vjust=0),
                  axis.title=element_text(face="bold", color="black", size=16, vjust=0.25)
            ) +
            ggtitle(paste("Cumulative sum fraction of", Ylab))
         
         if(max(stat_df[[yc]])/quantile(stat_df[[yc]], 0.9) > 100 ){
            p <- p + coord_cartesian(xlim=quantile(stat_df[[yc]], c(0.9, 1)))
         }
      }
       
   }

   return(p)
}


#' @title Plot signal profile around start, center, and end of genomic regions
#' @description Plot lines with standard error as the error band, also plots number of regions having non-zero signals
#'
#' @param plot_df a dataframe with column names c(xc, yc, cn, "Interval", "lower", "upper")
#' @param xc a string denoting column name for values on x-axis
#' @param yc a string denoting column name for numeric data to be plotted
#' @param cn a string denoting column name for grouping
#' @param ext a vector of 4 integers denoting upstream and downstream extension around start and end
#' @param hl a vector of 4 integers defining upstream and downstream boundaries of the rectangle for start and end
#' @param atitle a string for the title of the plot
#' @param insert a integer denoting the width of the center region
#' @param Ylab a string for y-axis label
#' @param shade logical indicating whether to place a shaded rectangle around the point of interest
#'
#' @return a ggplot object
#' @note used by \code{plot_start_end_feature}, \code{plot_start_end_reference_region}
#' @author Shuye Pu
#' @export draw_stacked_profile
#'

draw_stacked_profile <- function(plot_df, xc="Position", yc="Intensity", cn="Query", ext=c(0,0,0,0), hl=c(0,0,0,0), atitle="title", insert=0, Ylab="Signal Intensity", shade=FALSE){

   ylimits <- c(min(plot_df$lower), max(plot_df$upper))
   ylimits_intervals <- c(0, max(plot_df$Interval))

   aplot_df_start <- plot_df %>%
      filter(grepl("Start", Location))

   ps <- ggplot(aplot_df_start, aes(x=.data[[xc]], y=.data[[yc]], color=.data[[cn]])) + scale_fill_npg() + scale_color_npg() +
      geom_line(size=2) + ylim(ylimits) + #geom_point(color="grey30", size=2) +
      geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
      geom_ribbon(aes(ymin=lower, ymax=upper, fill=.data[[cn]]), linetype=0, alpha=0.3) +
      theme_classic() + theme(legend.position="none") +
      scale_x_continuous(limits = c(ext[1], ext[2])) +
      ylab(Ylab) +
      ggtitle(atitle) +
      theme(plot.margin = unit(c(1.2, 0.5, 0, 1.2), "lines"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.line.x = element_blank())
   if(shade) ps <- ps  + annotate("rect", xmin=hl[1], xmax=hl[2], ymin=-Inf, ymax=Inf, fill="grey", alpha=0.3)

   psi <- ggplot(aplot_df_start, aes(x=.data[[xc]], y=Interval, color=.data[[cn]])) + scale_fill_npg() + scale_color_npg() +
      geom_line(size=1) + ylim(ylimits_intervals) +
      geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
      theme_classic() +
      theme(plot.margin = unit(c(0, 0.5, 1.2, 1.2), "lines"),
            legend.position="none") + xlab("5'") +
      scale_x_continuous(limits = c(ext[1], ext[2])) +
      ylab("n")
   if(shade) psi <- psi  + annotate("rect", xmin=hl[1], xmax=hl[2], ymin=-Inf, ymax=Inf, fill="grey", alpha=0.3)
   ps_stack <-  plot_grid(ps, psi, ncol = 1, align = "v", rel_heights = c(1, 0.25))

   if(insert > 0){
      aplot_df_center <- plot_df %>%
         filter(grepl("Center", Location))
      pc <- ggplot(aplot_df_center, aes(x=.data[[xc]], y=.data[[yc]], color=.data[[cn]])) + scale_fill_npg() + scale_color_npg() +
         geom_line(size=2) + ylim(ylimits) + #geom_point(color="grey30", size=2) +
         geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
         geom_ribbon(aes(ymin=lower, ymax=upper, fill=.data[[cn]]), linetype=0, alpha=0.3) +
         #annotate("rect", xmin=unique(Xmin), xmax=unique(Xmax), ymin=-Inf, ymax=Inf, fill="grey", alpha=0.3) +
         theme_classic() + theme(legend.position="none") +
         ggtitle("") +
         theme(axis.text = element_blank(),
               axis.ticks = element_blank(),
               axis.title = element_blank(),
               axis.line = element_blank(),
               plot.margin = unit(c(1.2, 0.5, 0, 0.5), "lines"))

      pci <- ggplot(aplot_df_center, aes(x=.data[[xc]], y=Interval, color=.data[[cn]])) + scale_fill_npg() + scale_color_npg() +
         geom_line(size=1) + ylim(ylimits_intervals) +
         geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
         theme_classic() + theme(legend.position="none") + xlab("Center") +
         theme(axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               axis.title.y = element_blank(),
               axis.line.y = element_blank(),
               plot.margin = unit(c(0, 0.5, 1.2, 0.5), "lines"))
      pc_stack <- plot_grid(pc, pci, ncol = 1, align = "v", rel_heights = c(1, 0.25))
   }else{
      pc_stack <- NULL
   }

   aplot_df_end <- plot_df %>%
      filter(grepl("End", Location))
   pe <- ggplot(aplot_df_end, aes(x=.data[[xc]], y=.data[[yc]])) + scale_fill_npg() + scale_color_npg() +
      geom_line(aes(color=.data[[cn]]), size=2) + ylim(ylimits) + #geom_point(color="grey30", size=2) +
      geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
      geom_ribbon(aes(ymin=lower, ymax=upper, fill=.data[[cn]]), linetype=0, alpha=0.3) +
      theme_classic() +
      theme(legend.position=c(0.7, 0.9),
            legend.title=element_blank(),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.spacing.y = unit(1, "mm")) +
      scale_x_continuous(limits = c(ext[3], ext[4])) +
      ggtitle("") + guides(fill="none", color=guide_legend(keyheight = 0.75)) + # guides(fill="none") removes legend for fill
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            plot.margin = unit(c(1.2, 1.2, 0, 0.5), "lines"))

   if(shade) pe <- pe  + annotate("rect", xmin=hl[3], xmax=hl[4], ymin=-Inf, ymax=Inf, fill="grey", alpha=0.3)

   pei <- ggplot(aplot_df_end, aes(x=.data[[xc]], y=Interval, color=.data[[cn]])) + scale_fill_npg() + scale_color_npg() +
      geom_line(size=1) + ylim(ylimits_intervals) + #scale_y_continuous(position = "right") +
      theme_classic() + theme(legend.position="none") + xlab("3'") +
      geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=0.5) +
      scale_x_continuous(limits = c(ext[3], ext[4])) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = unit(c(0, 1.2, 1.2, 0.5), "lines"))
   if(shade) pei <- pei  + annotate("rect", xmin=hl[3], xmax=hl[4], ymin=-Inf, ymax=Inf, fill="grey", alpha=0.3)

   pe_stack <- plot_grid(pe, pei, ncol = 1, align = "v", rel_heights = c(1, 0.25))

   if(insert >0){
      p <- plot_grid(ps_stack, pc_stack, pe_stack, ncol = 3, align = "h", rel_widths = c(1, 0.5, 0.9))
   }else{
      p <- plot_grid(ps_stack, pe_stack, ncol = , align = "h", rel_widths = c(1, 0.9))
   }
   return(p)

}

