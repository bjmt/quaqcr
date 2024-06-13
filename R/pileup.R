#' Generate read pileups from BAMs with quaqc.
#'
#' quaqc maintains an internal TSS pileup in order to calcualte a TSS
#' enrichment score. This function takes advantage of this feature to
#' instead generate read pileups for arbitrary sets of regions.
#'
#' @param target.regions Either (1) a filename of a BED file containing
#' target region positions, or (2) a `GRanges` object from the
#' \pkg{GenomicRanges} package.
#' @param bkg.regions Either (1) a filename of a BED file containing
#' target region positions, or (2) a `GRanges` object from the
#' \pkg{GenomicRanges} package.
#' @param normalize "bkg": Converts read density into values
#' relative to the background (the first 25% of the window).
#' "rpm": Conver to reads per million. "no": Return as the average
#' number of reads per window.
#' @param region.size The input regions will be uniformly resized
#' to a single size.
#' @param qlen The size of the reads when they are included in the
#' pileup. A `qlen` of 0 means preserving the original read sizes;
#' otherwise the reads are resized from their 5-prime ends.
#' @param ... See [quaqc()].
#'
#' @return A `data.frame` containing read pileup data.
#'
#' @examples
#' \dontrun{
#' peaks <- system.file("extdata", "peaks.bed.gz", package = "quaqcr")
#' pileup(peaks, "Sample.bam")
#' }
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjmtremblay@@gmail.com}
#' @seealso [quaqcr::quaqc()], [quaqcr::melt_reports()]
#' @inheritParams quaqc
#' @export
pileup <- function(target.regions, bam.files, bkg.regions = NULL,
  normalize = c("rpm", "bkg", "no"), region.size = 5001, qlen = 0,
  verbose = 0, ...) {

  normalize <- match.arg(normalize)
  args <- list(...)
  if ("json" %in% names(args) && is.null(args$json)) stop("'json' cannot be NULL")
  if ("tss" %in% names(args)) warning("'tss' cannot be set")

  target <- quaqc(bam.files = bam.files, tss.size = region.size, tss.qlen = qlen,
    fast = TRUE, tss = target.regions, verbose = verbose, ...)
  target.pileup <- melt_reports(target, "tss_pileup", normalize.tss = normalize)
  target.pileup$Region <- "Target"

  # TODO: Loop through samples and divide by # of passing alignments

  if (!is.null(bkg.regions)) {
    bkg <- quaqc(bam.files = bam.files, tss.size = region.size, tss.qlen = qlen,
      fast = TRUE, tss = bkg.regions, verbose = verbose, ...)
    bkg.pileup <- melt_reports(bkg, "tss_pileup", normalize.tss = normalize)
    bkg.pileup$Region <- "Background"
    target.pileup <- rbind(target.pileup, bkg.pileup)
  }

  colnames(target.pileup) <- c("Sample", "Position", "Signal", "Target")

  target.pileup[, c(1, 4, 2, 3)]

}

