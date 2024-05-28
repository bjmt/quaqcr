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
#' @param ... See [quaqc()].
#'
#' @return A `data.frame` containing read pileup data.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjmtremblay@@gmail.com}
#' @seealso [base::system2()], [quaqcr::quaqc()]
#' @inheritParams quaqc
#' @export
footprint <- function(target.motifs, bam.files,
  bkg.motifs = NULL, tss.size = 501, tss.qlen = 1, tss.tn5 = TRUE,
  fast = TRUE, nfr = TRUE, verbose = 0, ...) {

  args <- list(...)
  if ("json" %in% names(args) && is.null(json)) stop("'json' cannot be NULL")
  if ("tss" %in% names(args)) warning("'tss' cannot be set")

  get_res <- function(x) {
    data.frame(row.names = NULL, Target = TRUE, Sample = x$sample,
      Position = x$filtered$nuclear$tss$pileup$x,
      Signal = x$filtered$nuclear$tss$pileup$y)
  }

  target <- quaqc(bam.files = bam.files, tss.size = tss.size, tss.qlen = tss.qlen,
    tss.tn5 = tss.tn5, fast = fast, nfr = nfr, tss = target.motifs,
    verbose = verbose, ...)
  target.tss <- lapply(target$reports, get_res)
  target.tss <- do.call(rbind, target.tss)

  bkg <- NULL
  if (!is.null(bkg.motifs)) {
    bkg <- quaqc(bam.files = bam.files, tss.size = tss.size, tss.qlen = tss.qlen,
      tss.tn5 = tss.tn5, fast = fast, nfr = nfr, tss = bkg.motifs,
      verbose = verbose, ...)
    bkg.tss <- lapply(bkg$reports, get_res)
    bkg.tss <- do.call(rbind, bkg.tss)
    bkg.tss$Target <- FALSE
  }

  if (!is.null(bkg.tss)) target.tss <- rbind(target.tss, bkg.tss)

  target.tss

}

