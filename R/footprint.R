#' Get transcription factor footprints with quaqc.
#'
#' The TSS pileup feature of quaqc can instead be used to get single
#' base resolution transcription factor footprints from ATAC-seq data.
#' This function provides a convenient wrapper around such
#' functionality. Target regions can be compared to unbound or background
#' regions. (Note that no Tn5 bias correction is applied.)
#'
#' @param target.motifs Either (1) a filename of a BED file containing
#' target motif positions, or (2) a `GRanges` object from the
#' \pkg{GenomicRanges} package.
#' @param bkg.motifs Either (1) a filename of a BED file containing
#' target motif positions, or (2) a `GRanges` object from the
#' \pkg{GenomicRanges} package. (Optional.)
#' @param normalize Converts read density into values
#' relative to the background (the first 25% of the window).
#' @param ... See [quaqc()].
#'
#' @return A `data.frame` containing read pileup data.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjmtremblay@@gmail.com}
#' @seealso [base::system2()], [quaqcr::quaqc()]
#' @inheritParams quaqc
#' @export
footprint <- function(target.motifs, bam.files, bkg.motifs = NULL,
  normalize = c("bkg", "rpm", "no"), tss.size = 501, tss.qlen = 1, tss.tn5 = TRUE,
  nfr = TRUE, verbose = 0, ...) {

  normalize <- match.arg(normalize)
  args <- list(...)
  if ("json" %in% names(args) && is.null(json)) stop("'json' cannot be NULL")
  if ("tss" %in% names(args)) warning("'tss' cannot be set")

  target <- quaqc(bam.files = bam.files, tss.size = tss.size, tss.qlen = tss.qlen,
    tss.tn5 = tss.tn5, fast = TRUE, nfr = nfr, tss = target.motifs,
    verbose = verbose, ...)
  target.pileup <- melt_reports(target, "tss_pileup",
    normalize.tss = normalize)
  target.pileup$Region <- "Target"

  if (!is.null(bkg.motifs)) {
    bkg <- quaqc(bam.files = bam.files, tss.size = tss.size, tss.qlen = tss.qlen,
      tss.tn5 = tss.tn5, fast = TRUE, nfr = nfr, tss = bkg.motifs,
      verbose = verbose, ...)
    bkg.pileup <- melt_reports(bkg, "tss_pileup",
      normalize.tss = normalize)
    bkg.pileup$Region <- "Background"
    target.pileup <- rbind(target.pileup, bkg.pileup)
  }

  # TODO: divide the pileup values by the number of total reads?

  colnames(target.pileup) <- c("Sample", "Distance", "Frequency", "Target")

  target.pileup[, c(1, 4, 2, 3)]

}

