
#' 
#' @title Toy data for examples and testing of the `GenomicPlot` package
#' @description The data files in the extdata directory contain data for next 
#' generation sequencing read alignments, MACS2 peaks and gene annotation, which 
#' are used to test the package and generate plots in the package vignettes. 
#' Except for the gtf file, all other files are derived from experimental data 
#' produced in-house at the \href{https://thedonnellycentre.utoronto.ca/faculty/jack-greenblatt}{Greenblatt Lab, University of Toronto, Canada}. 
#'     To meet the package file size limit, all data are restricted to 
#' chr19:58000-507000 of the human genome version hg19. Details for each file 
#' are as follows.
#' @details
#' \itemize{
#'  \item "gencode.v19.annotation_chr19.gtf" is an excerpt of the gene 
#'      annotation file downloaded from 
#'      \url{https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz}. 
#'      "gencode.v19.annotation_chr19.gtf.granges.rds" is a GRanges object as a 
#'      result of importing the gtf file.
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
#'      peak format, respectively. 
#'  \item "test_file1.txt", "test_file2.txt", "test_file3.txt" and 
#'      "test_file4.txt" are tab-delimited text files,  each contains various 
#'      human gene names in different columns. 
#' }
#' @aliases test_file1.txt test_file2.txt test_file3.txt test_file4.txt 
#'     gencode.v19.annotation_chr19.gtf chip_treat_chr19.bam treat_chr19.bam
#'     chip_input_chr19.bam input_chr19.bam test_wig_chr19_+.wig 
#'     test_wig_chr19_+.bw test_clip_peak_chr19.bed test_chip_peak_chr19.bed
#'     test_chip_peak_chr19.narrowPeak 
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
#' @format 
#' SQLlite database  
#' @source
#' The data is produced by running the following code:
#' \describe{
#' gtffile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf", 
#' package = "GenomicPlot")
#' }
#' \describe{
#' gff <- RCAS::importGtf(saveObjectAsRds = TRUE, readFromRds = FALSE, 
#'     filePath = gtffile)
#' }
#' \describe{
#' md <- data.frame(name = "Genome", value = "hg19")
#' }
#' \describe{
#' txdb <- GenomicFeatures::makeTxDbFromGRanges(gff, metadata = md)
#' }
#' \describe{
#' AnnotationDbi::saveDb(txdb, "./inst/extdata/txdb.sql")
#' }
#' @author Shuye Pu
#' @docType data
#' @name txdb.sql
#' @keywords datasets
#' 
NULL

#' 
#' @title Toy data for examples and testing of the `GenomicPlot` package
#' @description  Genomic  coordinates of 72 transcripts in hg19 for genomic 
#' features promoter, 5'UTR, CDS, 3'UTR, TTS. 
#' See \code{\link{prepare_5parts_genomic_features}} for details.
#' @format 
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
#' @source
#' The data is produced by running the following code:
#' \describe{
#' txdb <- AnnotationDbi::loadDb(system.file("extdata", "txdb.sql", 
#'     package = "GenomicPlot"))
#' }
#' \describe{
#' gf5_genomic <- GenomicPlot::prepare_5parts_genomic_features(txdb, 
#'     meta = FALSE, nbins = 100, fiveP = -2000, threeP = 1000, longest = TRUE)
#' }
#' @author Shuye Pu
#' @docType data
#' @name gf5_genomic
#' @keywords datasets
#' 
NULL

#' 
#' @title Toy data for examples and testing of the `GenomicPlot` package
#' @description  Metagenomic coordinates of 72 
#' transcripts in hg19 for genomic features promoter, 5'UTR, CDS, 3'UTR, TTS. 
#' See \code{\link{prepare_5parts_genomic_features}} for details.
#' @format 
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
#' @source
#' The data is produced by running the following code:
#' \describe{
#' txdb <- AnnotationDbi::loadDb(system.file("extdata", "txdb.sql", 
#'     package = "GenomicPlot"))
#' }
#' \describe{
#' gf5_meta <- GenomicPlot::prepare_5parts_genomic_features(txdb, meta = TRUE, 
#'     nbins = 100, fiveP = -2000, threeP = 1000, longest = TRUE)
#' }
#' 
#' @author Shuye Pu
#' @docType data
#' @name gf5_meta
#' @keywords datasets
#' 
NULL