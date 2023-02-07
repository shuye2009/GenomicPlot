
# GenomicPlot

<!-- badges: start -->
<!-- badges: end -->

The goal of GenomicPlot is to create a efficient visualization tool for next generation sequencing data with rich functionality and flexibility. 'GenomicPlot' facilitates plotting of NGS data in various formats (bam, bed, wig and bigwig); both coverage and enrichment over input can be computed and displayed with respect to genomic features (such as UTR, CDS, enhancer), and user defined genomic loci or regions. Statistical tests on signal intensity within user defined regions of interest can be performed and presented as boxplots or pie charts. Parallel processing is used to speed up computation on multi-core platforms. In addition to genomic plots which is suitable for displaying of coverage of genomic DNA (such as ChIP-seq data), metagenomic (without introns) plots can also be made for RNA-seq or CLIP-seq data as well. Peak annotation in intron-free regions is allowed.

## Installation

The following packages are prerequisites: 

GenomicRanges (>= 1.46.1), GenomicFeatures, Rsamtools, ggplot2 (>= 3.3.5), tidyr, rtracklayer (>= 1.54.0), plyranges (>= 1.14.0), dplyr (>= 1.0.8), cowplot (>= 1.1.1), VennDiagram, ggplotify, GenomeInfoDb, IRanges, ComplexHeatmap, RCAS (>= 1.20.0), scales (>= 1.2.0), GenomicAlignments (>= 1.30.0), edgeR, forcats, circlize, viridis, ggsignif (>= 0.6.3), ggsci (>= 2.9), genomation (>= 1.26.0)

You can install the development version of GenomicPlot from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("shuye2009/GenomicPlot")
```

## Example

This is a basic example which shows you how to visualize your ChIP-seq peaks in different parts of genes:

``` r
library(GenomicPlot)

## basic example code

txdb <- AnnotationDbi::loadDb(system.file("data", "txdb_chr19.sql", package="GenomicPlot"))
gf <- prepare_5parts_genomic_features(txdb, meta=TRUE, nbins=100, fiveP=1000, threeP=1000, longest=TRUE)

queryfiles <- c(system.file("data", "test_chip_peak_chr19.narrowPeak", package="GenomicPlot"),
                system.file("data", "test_chip_peak_chr19.bed", package="GenomicPlot"),
                system.file("data", "test_clip_peak_chr19.bed", package="GenomicPlot"))
names(queryfiles) <- c("narrowPeak", "summitPeak", "iCLIPPeak")

op <- "test_plot_5parts_metagene3"

handleInputParams <- list(CLIP_reads=FALSE, fix_width=150, fix_point="center", norm=FALSE, useScore=FALSE,
                          outRle=TRUE, useSizeFactor=FALSE, genome="hg19")

plot_5parts_metagene(queryfiles, gFeatures=list(metaF=gf), inputfiles=NULL, handleInputParams=handleInputParams,
                     verbose=FALSE, smooth=TRUE, scale=FALSE, stranded=TRUE, outPrefix=op, transform=FALSE, heatmap=TRUE,
                     rmOutlier=TRUE, nc=4)
```

The resulting figure is 'test_plot_5parts_metagene3.pdf'.

