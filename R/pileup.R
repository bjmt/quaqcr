
pileup <- function(target.regions, bam.files, bkg.regions = NULL,
  normalize = c("rpm", "bkg", "no"), region.size = 5001, qlen = 0,
  verbose = 0, ...) {

  normalize <- match.arg(normalize)
  args <- list(...)
  if ("json" %in% names(args) && is.null(json)) stop("'json' cannot be NULL")
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

