#' @title Hand input of NGS data with various formats
#'
#' @description This is a wrapper function for read NGS data in different file formats, store the input data in a list of GRanges objects or RleList objects. File names end in bed|bam|bw|bigwig|bigWig|BigWig|BW|BIGWIG are recognized, and a list of files with mixed formats are allowed.
#'
#' @param inputFiles a vector of strings denoting file names
#' @param handleInputParams a list with the following elements:
#' 'CLIP_reads' logical, indicating if the bam reads should be shifted to the -1 position at the 5' of the reads.
#' 'fix_width' an integer defines how long should the reads should be extended to.
#' 'fix_point' a string in c("start", "end", "center") denoting the anchor point for extension.
#' 'useScore' logical, indicating whether the 'score' column of the bed file should be used in calculation of coverage.
#' 'outRle' logical, indicating whether the output should be RleList objects or GRanges objects.
#' 'norm' logical, indicating whether the output RleList should be normalized to RPM using library sizes.
#' 'genome' a string denoting the genome name and version.
#' 'useSizeFactor' logical, indicating whether the library size should be adjusted with a size factor, using the 'calcNormFactors' function in the edgeR package
#' @param nc integer, number of cores for parallel processing
#' @param verbose logical, whether to output additional information
#'
#' @details when 'useScore' is TRUE, the score column of the bed file will be used in the metadata column 'score' of the GRanges object, or the 'Values' field of the RleList object. Otherwise the value 1 will be used instead. When the intended use of the input bed is a reference feature, both 'useScore' and 'outRle' should be set to FALSE.
#'
#' @return a list object with four elements, 'query' is a list GRanges objects or RleList objects, 'size' is the library size, 'type' is the input file type,
#' 'weight' is the name of the metadata column to be used as weight for coverage calculation
#'
#' @author Shuye Pu
#' @examples
#' centerfiles <- c(system.file("data", "test_B.bed", package="GenomicPlotData"),
#' system.file("data", "test_C.bed", package="GenomicPlotData"))
#' names(centerfiles) <- c("TestB", "TestC")
#'
#' @export handle_input

handle_input <- function(inputFiles, handleInputParams=NULL, verbose=FALSE, nc=2){

   if(any(is.null(names(inputFiles)))) stop("Each file must have a name attribute!")
   if(is.null(handleInputParams)) handleInputParams=list(CLIP_reads=FALSE, fix_width=0, fix_point="center", useScore=FALSE, outRle=TRUE, norm=TRUE, useSizeFactor=FALSE, genome="hg19")

   original_outRle <- handleInputParams$outRle
   if(handleInputParams$useSizeFactor && length(inputFiles)>1){
      handleInputParams$outRle <- FALSE ## effective_size only take granges_list which produce more accurate normFactor than RleList
   }

   inputFUN <- function(funName, inputFile, handleInputParams, verbose){
      fileName <- basename(inputFile)
      dirName <- dirname(inputFile)
      if(file.exists(file.path(dirName, paste0(fileName, ".rds")))){
         temp <- readRDS(file.path(dirName, paste0(fileName, ".rds")))
         if(identical(temp$param, handleInputParams)){
            out <- temp$Rle
            if(verbose)message("cached .rds file is used")
         }else{
            out <- funName(inputFile=inputFile, handleInputParams, verbose)
            saveRDS(list(param=handleInputParams, Rle=out), file.path(dirName, paste0(fileName, ".rds")))
            if(verbose)message("cached .rds file is modified using new input parameters")
         }
      }else{
         out <- funName(inputFile=inputFile, handleInputParams, verbose)
         saveRDS(list(param=handleInputParams, Rle=out), file.path(dirName, paste0(fileName, ".rds")))
         if(verbose)message("input data is cached as .rds file for fast reloading")
      }
      return(out)
   }
   
   outlist <- lapply(inputFiles, function(inputFile){
      if(!file.exists(inputFile)) stop("file does not exist, please check your file name and path!")
      if(grepl("\\.bed|BED|Bed|narrowPeak|broadPeak$", inputFile)){
         fileType <- "bed"
         if(verbose) print(paste("Reading", fileType, "file:", inputFile))
         out <- inputFUN(handle_bed, inputFile=inputFile, handleInputParams, verbose)
      }else if(grepl("\\.bam|BAM|Bam$", inputFile)){
         fileType <- "bam"
         if(verbose) print(paste("Reading", fileType, "file:", inputFile))
         out <- inputFUN(handle_bam, inputFile=inputFile, handleInputParams, verbose)
      }else if(grepl("\\.wig|WIG|Wig$", inputFile)){
         fileType <- "wig"
         if(verbose) print(paste("Reading", fileType, "file:", inputFile))
         out <- inputFUN(handle_wig,inputFile=inputFile, handleInputParams, verbose)
      }else if(grepl("\\.bw|bigwig|bigWig|BigWig|BW|BIGWIG$", inputFile)){
         fileType <- "bw"
         if(verbose) print(paste("Reading", fileType, "file:", inputFile))
         out <- inputFUN(handle_bw, inputFile=inputFile, handleInputParams, verbose)
      }else{
         stop("The format of file is not supported, please convert it to one of the following format: bed, bam, wig, bigwig")
      }

      out
   })

   names(outlist) <- names(inputFiles)
   handleInputParams$outRle <- original_outRle #if modified, restore

   ## modify library size
   if(handleInputParams$useSizeFactor && length(inputFiles)>1){
      outlist <- effective_size(outlist=outlist, outRle=handleInputParams$outRle, genome=handleInputParams$genome, nc=nc)
   }

   ## compute RPM
   if(handleInputParams$norm){
      outlist <- lapply(outlist, function(x){
         if(handleInputParams$outRle){
            x$query <- x$query*1e6/x$size
         }else{
            score(x$query) <- x$query$score*1e6/x$size
         }
         x
      })
   }
   invisible(outlist)
}


#' @title Normalize sample library size to effective size
#'
#' @description This is a helper function for handle_input. edgeR::calcNormFactors function is used to estimate normalizing factors,
#' which is used to multiply library sizes. The function only works for human genome only at present.
#'
#' @param outlist a list object with four elements, 'query' is a list GRanges objects or RleList objects, 'size' is the library size, 'type' is the input file type,
#' 'weight' is the name of the metadata column
#' @param outRle logical, indicating whether the output is a list of RleList objects or GRanges objects
#' @param genome a string denoting the genome name and version
#' @param nc integer, number of cores for parallel processing
#' @param verbose logical, whether to output additional information
#'
#' @return a list object with four elements ('query', 'size', 'type', 'weight'), with the 'size' element modified.
#'
#' @author Shuye Pu
#'
#' @export effective_size
#'
effective_size <- function(outlist, outRle, genome="hg19", nc=2, verbose=FALSE){
   if(verbose) print("Estimating size factor")

   seqi <- Seqinfo(genome=genome)

   grange_list <- lapply(outlist, function(x)x$query)

   tilewidth <- 100000
   tileBins <- tileGenome(seqi, tilewidth=tilewidth, cut.last.tile.in.chrom=TRUE)

   #score_list1 <-  biocParallel_binAverage(Rle_list, tileBins)
   #score_list <- parallel_binnedAverage(Rle_list, tileBins, nc=nc)
   score_list <- parallel_countOverlaps(grange_list, tileBins, nc=nc)

   mat <- data.matrix(bind_cols(score_list))
   mat[is.na(mat)] <- 0
   mat <- round(mat * tilewidth)
   mat <- mat[apply(mat, 1, sum)>0,]

   #normFactor <- DESeq2::estimateSizeFactorsForMatrix(mat)
   lib.size <- vapply(outlist, function(x)x$size, numeric(1))
   normFactor <- edgeR::calcNormFactors(mat, lib.size=lib.size)

   ## update library size with normFactor
   normlist <- mapply(x=outlist, y=normFactor, function(x,y){
      x$size <- as.integer(x$size*y)
      x
   }, SIMPLIFY=FALSE, USE.NAMES = TRUE)

   if(outRle){
      normlist <- lapply(normlist, function(x){
         x$query <- coverage(x$query, weight=x$weight)
         x
      })
   }

   names(normFactor) <- NULL
   if(verbose) print("Library size normalizing factors:")
   if(verbose) print(normFactor)
   invisible(normlist)
}

#' @title Handle files in bed format
#'
#' @description This is a function for read peaks data in bed format, store the input data in a list of GRanges objects or RleList objects.
#'
#' @param inputFile a string denoting path to the input file
#' @param handleInputParams a list of parameters for \code{handle_input}
#' @param verbose logical, whether to output additional information
#'
#' @return a list object with four elements, 'query' is a list GRanges objects or RleList objects, 'size' is the library size, 'type' is the input file type, 'weight' is the name of the metadata column to be used as weight for coverage calculation
#'
#' @author Shuye Pu
#'
#' @export handle_bed

handle_bed <- function(inputFile, handleInputParams=NULL, verbose=FALSE){

   beddata <- read.delim2(inputFile, header=FALSE, comment.char = "#")

   nco <- ncol(beddata)
   if(nco > 6 && verbose){
      message(paste("The input file", inputFile, "have more than 6 columns, only the first 6 columns will be used!"))
   }else if(nco < 6 && verbose){
      message(paste("The input file", inputFile, "have only", nco, "columns!"))
   }

   beddata <- type.convert(beddata[, 1:min(6,nco)], as.is=TRUE)  ## ignore extra columns, which cause problem in import.bed()
   colnames(beddata) <- c("chr", "start", "end", "name", "score", "strand")[1:min(6,ncol(beddata))]
   queryRegions <- makeGRangesFromDataFrame(beddata, keep.extra.columns=TRUE, starts.in.df.are.0based=TRUE)
   names(queryRegions) <- beddata$name
   if(sum(duplicated(names(queryRegions))) > 0){ ## if the names are not unique, force them to be unique
      names(queryRegions) <- paste(names(queryRegions), seq_along(names(queryRegions)), sep="_")  
   }
   if(ncol(beddata)<6) strand(queryRegions) <- "*"

   if(handleInputParams$fix_width > 0) queryRegions <- resize(queryRegions, width=handleInputParams$fix_width,
                                                              fix=rep(handleInputParams$fix_point, length(queryRegions)), ignore.strand=FALSE)
   weight_col <- "score"

   if(!handleInputParams$useScore || ncol(beddata)<5){
      score(queryRegions) <- 1
   }

   ## make input comply with GenomeInfoDb
   seqInfo <- Seqinfo(genome=handleInputParams$genome)
   seqInfo <- keepStandardChromosomes(seqInfo)
   
   if(c("chr1") %in% as.vector(seqnames(queryRegions)) && seqnames(seqInfo)[1] == "1"){
      seqnames(seqInfo) <- paste0("chr", seqnames(seqInfo))
      seqnames(seqInfo)[25] = "chrM"
   }else if(c("1") %in% as.vector(seqnames(queryRegions)) && seqnames(seqInfo)[1] == "chr1"){
      seqnames(seqInfo) <- gsub("chr", "", seqnames(seqInfo))
      seqnames(seqInfo)[25] = "MT"
   }  
   queryRegions <- queryRegions[as.vector(seqnames(queryRegions)) %in% seqnames(seqInfo)]
   seqlevels(queryRegions) <- seqlevels(seqInfo)
   seqinfo(queryRegions) <- seqInfo

   libsize <- sum(queryRegions$score, na.rm=TRUE)
  
   if("name" %in% colnames(mcols(queryRegions))) names(queryRegions) <- queryRegions$name
  
   if(handleInputParams$outRle){
      queryRegions <- coverage(queryRegions, weight=weight_col)
      seqinfo(queryRegions) <- seqInfo
   }

   invisible(list("query"=queryRegions, "size"=libsize, "type"="bed", "weight"=weight_col))

}

#' @title Handle files in bam format
#'
#' @description This is a function for read NGS reads data in bam format, store the input data in a list of GRanges objects or RleList objects. For paired-end reads, only take the second read in a pair, assuming which is the sense read for strand-specific RNAseq.
#'
#' @param inputFile a string denoting path to the input file
#' @param handleInputParams a list of parameters for \code{handle_input}
#' @param verbose logical, whether to output additional information
#'
#' @details The reads are filtered using mapq score >= 10 by default, only mapped reads are counted towards library size.
#'
#' @return a list object with four elements, 'query' is a list GRanges objects or RleList objects, 'size' is the library size, 'type' is the input file type, weight' is the name of the metadata column to be used as weight for coverage calculation
#'
#' @author Shuye Pu
#'
#' @export handle_bam
#'
handle_bam <- function(inputFile, handleInputParams=NULL, verbose=FALSE){

   paired.end <- testPairedEndBam(inputFile)
   if(paired.end){
      flag <- scanBamFlag(isPaired=TRUE, isProperPair=TRUE, isUnmappedQuery=FALSE, hasUnmappedMate=FALSE, isSecondMateRead=TRUE)  ## assume the second read in a pair is the sense read for strand-specific RNAseq
      param <- ScanBamParam(mapqFilter=10, flag=flag)
      if(verbose) message("Paired end bam file detected!")
   }else{
      param <- ScanBamParam(mapqFilter=10)
   }
   if(verbose) message(paste("\nbam file", inputFile, "is loaded\n"))
   
   ga <- readGAlignments(inputFile, use.names=TRUE, param=param)
   libsize <- sum(idxstatsBam(inputFile)$mapped)
   if(handleInputParams$CLIP_reads){
      ## get the 5'-end -1 position of the reads, which is the crosslink sites for iCLIP reads
      queryRegions <- flank(granges(ga), width=1, both=FALSE, start=TRUE, ignore.strand=FALSE)
      score(queryRegions) <- 1
   }else if(handleInputParams$fix_width > 0){
      queryRegions <- resize(granges(ga), width=handleInputParams$fix_width, fix=rep("start", length(granges(ga))), ignore.strand=FALSE)
      score(queryRegions) <- 1
   }else{
      queryRegions <- unlist(grglist(ga))
      score(queryRegions) <- 1
   }
   weight_col <- "score"

   ## make input comply with GenomeInfoDb
   seqInfo <- Seqinfo(genome=handleInputParams$genome)
   seqInfo <- keepStandardChromosomes(seqInfo)
   queryRegions <- queryRegions[as.vector(seqnames(queryRegions)) %in% seqnames(seqInfo)]
   seqlevels(queryRegions) <- seqlevels(seqInfo)
   seqinfo(queryRegions) <- seqInfo
   

   if(handleInputParams$outRle){
      queryRegions <- coverage(queryRegions, weight= weight_col)
      seqinfo(queryRegions) <- seqInfo
   }

   invisible(list("query"=queryRegions, "size"=libsize, "type"="bam", "weight"=weight_col))

}


#' @title Handle files in bw|bigwig|bigWig|BigWig|BW|BIGWIG format
#'
#' @description This is a function for read NGS coverage data in bigwig format, store the input data in a list of GRanges objects or RleList objects. The input bw file can be stranded or non-stranded. Library size is calculate as the sum of all coverage.
#'
#' @param inputFile a string denoting path to the input file
#' @param handleInputParams a list of parameters for \code{handle_input}
#' @param verbose logical, whether to output additional information
#'
#' @details For stranded files, forward and reverse strands are stored in separate files, with '+' or 'p' in the forward strand file name and '-' or 'm' in the reverse strand  file name.
#'
#' @return a list object with four elements, 'query' is a list GRanges objects or RleList objects, 'size' is the estimated library size, 'type' is the input file type, weight' is the name of the metadata column to be used as weight for coverage calculation
#'
#' @author Shuye Pu
#'
#' @export handle_bw
#'
handle_bw <- function(inputFile, handleInputParams, verbose=FALSE){

   weight_col <- "score"
   libsize <- NULL

   neg_file <- find_mate(inputFile, verbose)
   stranded <- ifelse(!is.null(neg_file), TRUE, FALSE)

   if(stranded){
      posGR <- import.bw(inputFile)
      strand(posGR) <- "+"

      negGR <- import.bw(neg_file)
      strand(negGR) <- "-"

      queryRegions <- append(posGR, negGR)
   }else{
      iGR <- import.bw(inputFile)

      queryRegions <- iGR
      strand(queryRegions) <- "*"
   }

   ## make input comply with GenomeInfoDb
   seqInfo <- Seqinfo(genome=handleInputParams$genome)
   seqInfo <- keepStandardChromosomes(seqInfo)
   queryRegions <- queryRegions[as.vector(seqnames(queryRegions)) %in% seqnames(seqInfo)]
   seqlevels(queryRegions) <- seqlevels(seqInfo)
   seqinfo(queryRegions) <- seqInfo

   fragmentLength <- 100 ## this is an assumption
   libsize <- as.integer(sum(score(queryRegions) * width(queryRegions))/fragmentLength)

   if(handleInputParams$outRle){
      queryRegions <- coverage(queryRegions, weight=weight_col)
      seqinfo(queryRegions) <- seqInfo
   }

   invisible(list("query"=queryRegions, "size"=libsize, "type"="bw", "weight"=weight_col))
}

#' @title Handle files in wig format
#'
#' @description This is a function for read NGS coverage data in wig format, store the input data in a list of GRanges objects or RleList objects. The input wig file can be stranded or non-stranded. Library size is calculate as the sum of all coverage.
#'
#' @param inputFile a string denoting path to the input file
#' @param handleInputParams a list of parameters for \code{handle_input}
#' @param verbose logical, whether to output additional information
#'
#' @details For stranded files, forward and reverse strands are stored in separate files, with '+' or 'p' in the forward strand file name and '-' or 'm' in the reverse strand file name.
#'
#' @return a list object with four elements, 'query' is a list GRanges objects or RleList objects, 'size' is the library size, 'type' is the input file type, 'weight' is the name of the metadata column to be used as weight for coverage calculation
#'
#' @author Shuye Pu
#'
#' @export handle_wig
#'
handle_wig <- function(inputFile, handleInputParams, verbose=FALSE){

   neg_file <- find_mate(inputFile, verbose)
   stranded <- ifelse(!is.null(neg_file), TRUE, FALSE)

   seqinfo <- Seqinfo(genome=handleInputParam$genome)
   wigToBigWig(inputFile, seqinfo)

   if(stranded){
      wigToBigWig(neg_file, seqinfo)
   }

   bwfile <- gsub("\\.wig", "\\.bw", inputFile)

   out <- handle_bw(bwfile, handleInputParams, verbose)

   invisible(out)
}


#' @title Find wig/bw file for the negative strand
#' @description Find the file name of the negative strand, if a .wig/bw file for positive strand if provided, by looking for file names with one character difference.
#' If no negative strand file is found, assume the input .wig/bw file is non-stranded
#'
#' @param inputFile path to a .wig/bw file, presumably for positive strand
#' @param verbose logical, whether to output additional information
#'
#' @return path to the negative .wig/bw file or NULL
#'
#' @keywords internal


find_mate <- function(inputFile, verbose=FALSE){
   fileName <- basename(inputFile)
   dirName <- dirname(inputFile)

   fch_v <- unlist(strsplit(fileName,split=""))
   otherFiles <- list.files(dirName)

   mate <- NULL
   for(afile in otherFiles){
      och_v <- unlist(strsplit(afile,split=""))
      if(length(fch_v) == length(och_v) && sum(fch_v != och_v) == 1){
         mate <- file.path(dirName, afile)
         diff <- base::setdiff(och_v, fch_v)
         if(verbose) print(diff)
      }
   }
   invisible(mate)
}


