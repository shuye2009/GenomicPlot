# managing namespace

#' @import methods
#' @import ggplot2
#' @import Rsamtools
#' @import GenomicFeatures
#' @import GenomicRanges
#' @importFrom parallel detectCores clusterExport parLapply makeCluster stopCluster
#' @importFrom dplyr bind_cols bind_rows mutate filter arrange %>% group_by select count full_join left_join desc inner_join
#' @importFrom dplyr rename_with summarize n_distinct lead n
#' @importFrom plyranges filter_by_overlaps
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomeInfoDb Seqinfo seqinfo seqlevels seqnames seqinfo<- seqlevels<- seqnames<-
#' @importFrom rtracklayer export.bed import.bw wigToBigWig asBED
#' @importFrom RCAS importGtf queryGff
#' @importFrom ggpubr ggtexttable
#' @importFrom ggplotify as.grob
#' @importFrom ggsci pal_npg scale_color_npg scale_fill_npg
#' @importFrom ggsignif geom_signif
#' @importFrom cowplot plot_grid
#' @importFrom forcats fct_inorder
#' @importFrom edgeR calcNormFactors
#' @importFrom genomation ScoreMatrixBin annotateWithFeatures
#' @importFrom BiocGenerics strand strand<- score score<- append start start<- end end<- Reduce
#' @importFrom grid grid.grabExpr grid.draw grid.text grid.newpage gpar
#' @importFrom tidyr pivot_longer
#' @importFrom scales percent
#' @importFrom IRanges mergeByOverlaps width
#' @importFrom GenomicRanges trim shift resize promoters flank countOverlaps coverage granges grglist tileGenome setdiff narrow
#' @importFrom GenomicRanges binnedAverage findOverlaps GRanges GRangesList mcols mcols<- makeGRangesFromDataFrame
#' @importFrom GenomicFeatures transcriptLengths exonsBy intronsByTranscript threeUTRsByTranscript fiveUTRsByTranscript
#' @importFrom GenomicFeatures cdsBy genes transcripts makeTxDbFromGRanges makeTxDbFromGFF makeTxDbFromEnsembl
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap draw
#' @importFrom circlize colorRamp2
#' @importFrom viridis viridis scale_fill_viridis
#' @importFrom VennDiagram draw.pairwise.venn draw.triple.venn draw.quad.venn
#' @importFrom grDevices dev.list dev.off pdf
#' @importFrom graphics abline hist pairs panel.smooth par points rect strwidth text
#' @importFrom utils combn head read.delim read.delim2 read.table tail type.convert write.table
#' @importFrom stats TukeyHSD aggregate aov approx
#' @importFrom stats as.formula cor density dist ecdf hclust
#' @importFrom stats ks.test mad median na.omit qqplot quantile
#' @importFrom stats reorder sd smooth.spline wilcox.test
NULL
