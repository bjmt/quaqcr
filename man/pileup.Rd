% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pileup.R
\name{pileup}
\alias{pileup}
\title{Generate read pileups from BAMs with quaqc.}
\usage{
pileup(target.regions, bam.files, bkg.regions = NULL, normalize = c("rpm",
  "bkg", "no"), region.size = 5001, qlen = 0, verbose = 0, ...)
}
\arguments{
\item{target.regions}{Either (1) a filename of a BED file containing
target region positions, or (2) a \code{GRanges} object from the
\pkg{GenomicRanges} package.}

\item{bam.files}{Character vector of BAM file names. Must be coordinate
sorted. If no index file can be found, \code{quaqc} will generate them.}

\item{bkg.regions}{Either (1) a filename of a BED file containing
target region positions, or (2) a \code{GRanges} object from the
\pkg{GenomicRanges} package.}

\item{normalize}{"bkg": Converts read density into values
relative to the background (the first 25\% of the window).
"rpm": Conver to reads per million. "no": Return as the average
number of reads per window.}

\item{region.size}{The input regions will be uniformly resized
to a single size.}

\item{qlen}{The size of the reads when they are included in the
pileup. A \code{qlen} of 0 means preserving the original read sizes;
otherwise the reads are resized from their 5-prime ends.}

\item{verbose}{Integer, a value from 0 to 2 for the level of program verbosity.}

\item{...}{See \code{\link[=quaqc]{quaqc()}}.}
}
\value{
A \code{data.frame} containing read pileup data.
}
\description{
quaqc maintains an internal TSS pileup in order to calcualte a TSS
enrichment score. This function takes advantage of this feature to
instead generate read pileups for arbitrary sets of regions.
}
\examples{
\dontrun{
peaks <- system.file("extdata", "peaks.bed.gz", package = "quaqcr")
pileup(peaks, "Sample.bam")
}

}
\seealso{
\code{\link[=quaqc]{quaqc()}}, \code{\link[=melt_reports]{melt_reports()}}
}
\author{
Benjamin Jean-Marie Tremblay, \email{benjmtremblay@gmail.com}
}
