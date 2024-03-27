#'
#' @title Toy data for examples and testing of the `GenomicPlot` package
#' @description The data files in the extdata directory contain data for next
#' generation sequencing read alignments, MACS2 peaks and gene annotation, which
#' are used to test the package and generate plots in the package vignettes.
#'     To meet the package file size limit, all data are restricted to
#' chr19:58000-507000 of the human genome version hg19. Details for each file
#' are as follows.
#' @details
#' \itemize{
#'  \item "gencode.v19.annotation_chr19.gtf" is an excerpt of a gene
#'      annotation file by limiting to chr19:58000-507000 of the human genome.
#'  \item "gencode.v19.annotation_chr19.gtf.granges.rds" is a GRanges object
#'      produced by importing the above gtf file using RCAS::importGtf.
#'  \item "chip_treat_chr19.bam(.bai)" and "chip_input_chr19.bam(.bai)" are
#'      paired-end read alignment data from ChIPseq experiments.
#'  \item "treat_chr19.bam(.bai)" and "input_chr19.bam(.bai)" are single-end
#'      read alignment data from iCLIP experiments.
#'  \item "test_wig_chr19_+(-).wig", "test_wig_chr19_+(-).bw" are iCLIP
#'      alignment data in WIG and BIGWIG format, respectively; '+' and '-'
#'      represent forward and reverse strand, respectively.
#'  \item "test_clip_peak_chr19.bed" contains strand-specific iCLIP peak in BED
#'      format.
#'  \item "test_chip_peak_chr19.bed" and "test_chip_peak_chr19.narrowPeak"
#'      contain ChIPseq peaks generated with MACS2, in summit peak and narrow
#'      peak format, respectively. "test_chr19.bedGraph" contains the same data
#'      in bedGraph format.
#'  \item "test_file1.txt", "test_file2.txt", "test_file3.txt" and
#'      "test_file4.txt" are tab-delimited text files,  each contains various
#'      human gene names in different columns.
#' }
#' @source The original gene annotation (gtf) file is downloaded from
#'      \href{https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz}{https://www.gencodegenes.org/human/}. \cr
#'      Except for the gtf file, all other files are derived from experimental
#'      data produced in-house at the \href{https://thedonnellycentre.utoronto.ca/faculty/jack-greenblatt}{Greenblatt Lab, University of Toronto, Canada}.
#' @return Various files used as inputs to run examples and tests
#' @aliases test_file1.txt test_file2.txt test_file3.txt test_file4.txt
#'     gencode.v19.annotation_chr19.gtf chip_treat_chr19.bam treat_chr19.bam
#'     chip_input_chr19.bam input_chr19.bam test_wig_chr19_+.wig
#'     test_wig_chr19_+.bw test_clip_peak_chr19.bed test_chip_peak_chr19.bed
#'     test_chip_peak_chr19.narrowPeak test_chr19.bedGraph
#' @author Shuye Pu
#' @docType data
#' @name extdata
#' @keywords datasets
#'
NULL


#'
#' @title Toy data for examples and testing of the `GenomicPlot` package
#' @description A tiny TxDb object holding genomic feature coordinates of 72
#' transcripts in hg19.
#' @source The data is produced by running the following code: \cr
#' gtffile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
#' package = "GenomicPlot") \cr
#' txdb <- custom_TxDb_from_GTF(gtffile, genome = "hg19") \cr
#' AnnotationDbi::saveDb(txdb, "./inst/extdata/txdb.sql")
#'
#' @return
#' A SQLlite database
#'
#' @author Shuye Pu
#' @docType data
#' @name txdb.sql
#' @keywords datasets
#'
NULL

#'
#' @title Toy data for examples and testing of the `GenomicPlot` package
#' @description  Genomic  coordinates of 72 transcripts in hg19 for genomic
#' features promoter, 5'UTR, CDS, 3'UTR, TTS, as well as user inputs for
#' processing these features.
#' See \code{\link{prepare_5parts_genomic_features}} for details.
#'
#' @source The data is produced by running the following code: \cr
#' txdb <- AnnotationDbi::loadDb(system.file("extdata", "txdb.sql",
#'     package = "GenomicPlot")) \cr
#' gf5_genomic <- GenomicPlot::prepare_5parts_genomic_features(txdb,
#'     meta = FALSE, nbins = 100, fiveP = -2000, threeP = 1000, longest = TRUE)
#'
#' @return
#' A named list with the following elements:
#' \describe{
#'  \item{windowRs}{a list of 5 GrangesList objects for the 5 genomic features}
#'  \item{nbins}{a positive integer}
#'  \item{scaled_bins}{a vector of 5 integers}
#'  \item{fiveP}{a negative integer}
#'  \item{threeP}{a positive integer}
#'  \item{meta}{logical}
#'  \item{longest}{logical}
#' }
#'
#' @author Shuye Pu
#' @docType data
#' @name gf5_genomic
#' @keywords datasets
#'
NULL

#'
#' @title Toy data for examples and testing of the `GenomicPlot` package
#' @description  Metagenomic coordinates of 72 transcripts in hg19 for genomic
#' features promoter, 5'UTR, CDS, 3'UTR, TTS, as well as user inputs for
#' processing these features.
#' See \code{\link{prepare_5parts_genomic_features}} for details.
#'
#' @source The data is produced by running the following code: \cr
#' txdb <- AnnotationDbi::loadDb(system.file("extdata", "txdb.sql",
#'     package = "GenomicPlot")) \cr
#' gf5_meta <- GenomicPlot::prepare_5parts_genomic_features(txdb, meta = TRUE,
#'     nbins = 100, fiveP = -2000, threeP = 1000, longest = TRUE)
#'
#' @return
#' A named list with the following elements:
#' \describe{
#'  \item{windowRs}{a list of 5 GrangesList objects for the 5 genomic features}
#'  \item{nbins}{a positive integer}
#'  \item{scaled_bins}{a vector of 5 integers}
#'  \item{fiveP}{a negative integer}
#'  \item{threeP}{a positive integer}
#'  \item{meta}{logical}
#'  \item{longest}{logical}
#' }
#'
#' @author Shuye Pu
#' @docType data
#' @name gf5_meta
#' @keywords datasets
#'
NULL
