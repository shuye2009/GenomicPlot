% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenomicPlot.R
\name{handle_bw}
\alias{handle_bw}
\title{handle_bw}
\usage{
handle_bw(inputFile, outRle = TRUE, genome = "hg19")
}
\arguments{
\item{inputFile, }{a string denoting the bw file name}

\item{outRle, }{a boolean object indicating whether the output should be a list of RleList objects or GRanges objects}

\item{genome, }{a string denoting the genome version}
}
\value{
a list object with four elements, 'query' is a list GRanges objects or RleList objects, 'size' is the library size, 'type' is the input file type,
'weight' is the name of the metadata column to be used as weight for coverage calculation
}
\description{
handle_bw
}
\details{
For stranded files, forward and reverse strands are stored in separate files, with '+' or 'p' in the forward strand file name and '-' or 'm' in the reverse
strand  file name.
}
\author{
Shuye Pu
}
