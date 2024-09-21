#' Melt sections of a quaqc report into a data.frame.
#' 
#' The quaqc report class type in R is divided into lists of lists,
#' which can require additional manipulation. This function will
#' "melt" these individual sections into \code{data.frame} objects.
#'
#' @param report A \code{quaqc} object.
#' @param section The \code{quaqc} object subsection to melt.
#' @param use.basename Whether to use the [base::basename()] function
#' on the sample names.
#' @param normalize.tss How to normalize the TSS pileup. "no": Keep
#' the signal as the average number of reads per window. "bkg":
#' Calculate the signal relative to the background (the first 25%
#' of the window). "rpm": Convert to reads per million.
#' @param normalize.hist How to normalize the alignment size,
#' fragment size, GC percent, and read depth histograms.
#' "no": Keep as the total number of reads per bin. "proportion":
#' Divide by the sum of reads across all windows. "max": Divide
#' by the max bin count.
#'
#' @return A \code{data.frame} with varying columns based on the
#' section being melted.
#'
#' @examples
#' report.file <- system.file("extdata", "report.json.gz", package = "quaqcr")
#' report <- parse_quaqc_file(report.file)
#' melt_reports(report, "overview_filt")
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjmtremblay@@gmail.com}
#' @seealso [quaqcr::quaqc()], [quaqcr::parse_quaqc_file()]
#' @export
melt_reports <- function(report, section = c("bam_stats",
    "overview_unfilt", "overview_filt", "nucl_stats", "nucl_addn",
    "peak_stats", "tss_stats", "tss_pileup", "aln_hist", "frag_hist", "gc_hist",
    "depth_hist", "genome"),
  use.basename = FALSE, normalize.tss = c("no", "bkg", "rpm"),
  normalize.hist = c("no", "proportion", "max")) {

  if (!is(report, "quaqc")) stop("'report' must be a 'quaqc' class object")

  section <- match.arg(section)
  normalize.tss <- match.arg(normalize.tss)
  normalize.hist <- match.arg(normalize.hist)

  reports <- report$reports
  sample.names <- sapply(reports, function(x) x$sample)

  out <- switch(section,
    "bam_stats" = get_bam_stats(reports, sample.names),
    "overview_unfilt" = get_overview_unfilt(reports, sample.names),
    "overview_filt" = get_overview_filt(reports, sample.names),
    "nucl_stats" = get_nucl_stats(reports, sample.names),
    "nucl_addn" = get_nucl_addn(reports, sample.names),
    "peak_stats" = get_peak_stats(reports, sample.names),
    "tss_stats" = get_tss_stats(reports, sample.names),
    "tss_pileup" = get_tss_pileup(reports, sample.names, normalize.tss),
    "aln_hist" = get_aln_hist(reports, sample.names, normalize.hist),
    "frag_hist" = get_frag_hist(reports, sample.names, normalize.hist),
    "gc_hist" = get_gc_hist(reports, sample.names, normalize.hist),
    "depth_hist" = get_depth_hist(reports, sample.names, normalize.hist),
    "genome" = get_genome(reports, sample.names)
  )

  if (use.basename) out$Sample <- basename(out$Sample)

  out

}

get_depth_hist <- function(r, sn, normhist) {
  nRow <- sapply(r, function(x) length(x$filtered$nuclear$histograms$depth$x))
  o <- data.frame(row.names = NULL,
    Sample = rep(sn, times = nRow),
    ReadDepth = NA,
    Count = NA
  )
  for (i in seq_along(r)) {
    i_i <- (1:(nRow[i])) + sum(nRow[0:(i - 1)])
    o$ReadDepth[i_i] <- r[[i]]$filtered$nuclear$histograms$depth$x
    o$Count[i_i] <- r[[i]]$filtered$nuclear$histograms$depth$y
    if (normhist == "proportion") {
      o$Count[i_i] <- o$Count[i_i] / sum(o$Count[i_i])
    } else if (normhist == "max") {
      o$Count[i_i] <- o$Count[i_i] / max(o$Count[i_i])
    }
  }
  o
}

get_gc_hist <- function(r, sn, normhist) {
  nRow <- sapply(r, function(x) length(x$filtered$nuclear$histograms$gc$x))
  o <- data.frame(row.names = NULL,
    Sample = rep(sn, times = nRow),
    GCPercent = NA,
    Count = NA
  )
  for (i in seq_along(r)) {
    i_i <- (1:(nRow[i])) + sum(nRow[0:(i - 1)])
    o$GCPercent[i_i] <- r[[i]]$filtered$nuclear$histograms$gc$x
    o$Count[i_i] <- r[[i]]$filtered$nuclear$histograms$gc$y
    if (normhist == "proportion") {
      o$Count[i_i] <- o$Count[i_i] / sum(o$Count[i_i])
    } else if (normhist == "max") {
      o$Count[i_i] <- o$Count[i_i] / max(o$Count[i_i])
    }
  }
  o
}

get_frag_hist <- function(r, sn, normhist) {
  nRow <- sapply(r, function(x) length(x$filtered$nuclear$histograms$fragment$x))
  o <- data.frame(row.names = NULL,
    Sample = rep(sn, times = nRow),
    FragSize = NA,
    Count = NA
  )
  for (i in seq_along(r)) {
    i_i <- (1:(nRow[i])) + sum(nRow[0:(i - 1)])
    o$FragSize[i_i] <- r[[i]]$filtered$nuclear$histograms$fragment$x
    o$Count[i_i] <- r[[i]]$filtered$nuclear$histograms$fragment$y
    if (normhist == "proportion") {
      o$Count[i_i] <- o$Count[i_i] / sum(o$Count[i_i])
    } else if (normhist == "max") {
      o$Count[i_i] <- o$Count[i_i] / max(o$Count[i_i])
    }
  }
  o
}

get_aln_hist <- function(r, sn, normhist) {
  nRow <- sapply(r, function(x) length(x$filtered$nuclear$histograms$alignment$x))
  o <- data.frame(row.names = NULL,
    Sample = rep(sn, times = nRow),
    AlnSize = NA,
    Count = NA
  )
  for (i in seq_along(r)) {
    i_i <- (1:(nRow[i])) + sum(nRow[0:(i - 1)])
    o$AlnSize[i_i] <- r[[i]]$filtered$nuclear$histograms$alignment$x
    o$Count[i_i] <- r[[i]]$filtered$nuclear$histograms$alignment$y
    if (normhist == "proportion") {
      o$Count[i_i] <- o$Count[i_i] / sum(o$Count[i_i])
    } else if (normhist == "max") {
      o$Count[i_i] <- o$Count[i_i] / max(o$Count[i_i])
    }
  }
  o
}

get_tss_pileup <- function(r, sn, normtss) {
  nRow <- sapply(r, function(x) length(x$filtered$nuclear$tss$pileup$x))
  o <- data.frame(row.names = NULL,
    Sample = rep(sn, times = nRow),
    Coordinate = NA,
    Depth = NA
  )
  for (i in seq_along(r)) {
    i_i <- (1:(nRow[i])) + sum(nRow[0:(i - 1)])
    o$Coordinate[i_i] <- r[[i]]$filtered$nuclear$tss$pileup$x
    o$Depth[i_i] <- r[[i]]$filtered$nuclear$tss$pileup$y
    if (normtss == "bkg") {
      bRange <- r[[i]]$filtered$nuclear$tss$pileup$range
      bRange <- seq(bRange[1], bRange[2])
      bRange <- bRange[seq_len(length(bRange) * 0.25)]
      bDepths <- o$Depth[i_i][o$Coordinate[i_i] %in% bRange]
      bDepths <- c(bDepths, rep(0, length(bRange) - length(bDepths)))
      o$Depth[i_i] <- o$Depth[i_i] / mean(bDepths)
    } else if (normtss == "rpm") {
      # TODO: What if some ranges are outside the effective area?
      nranges <- r[[i]]$params$integer["tss_bed_n"]
      aln.total <- r[[i]]$filtered$overview[1, 1]
      o$Depth[i_i] <- (o$Depth[i_i] * (1e6 / aln.total)) / nranges
    } else {
      nranges <- r[[i]]$params$integer["tss_bed_n"]
      o$Depth[i_i] <- o$Depth[i_i] / nranges
    }
  }
  o
}

get_tss_stats <- function(r, sn) {
  nRow <- 1
  o <- data.frame(row.names = NULL,
    Sample = rep(sn, each = nRow),
    TSSCount = NA,
    TES = NA
  )
  for (i in seq_along(r)) {
    i_i <- (1:nRow) + nRow * (i - 1)
    o$TSSCount[i_i] <- r[[i]]$filtered$nuclear$tss$stats["n"]
    o$TES[i_i] <- r[[i]]$filtered$nuclear$tss$stats["tss_enrichment_score"]
  }
  o
}

get_peak_stats <- function(r, sn) {
  nRow <- 1
  o <- data.frame(row.names = NULL,
    Sample = rep(sn, each = nRow),
    PeakCount = NA,
    PeakGenomeCov = NA,
    FRIP = NA
  )
  for (i in seq_along(r)) {
    i_i <- (1:nRow) + nRow * (i - 1)
    o$PeakCount[i_i] <- r[[i]]$filtered$nuclear$peaks["n"]
    o$PeakGenomeCov[i_i] <- r[[i]]$filtered$nuclear$peaks["coverage"]
    o$FRIP[i_i] <- r[[i]]$filtered$nuclear$peaks["fraction_of_reads_in_peaks"]
  }
  o
}

get_nucl_addn <- function(r, sn) {
  nRow <- 1
  o <- data.frame(row.names = NULL,
    Sample = rep(sn, each = nRow),
    GenomeCoverage = NA,
    AlnNoMAPQ = NA
  )
  for (i in seq_along(r)) {
    i_i <- (1:nRow) + nRow * (i - 1)
    o$GenomeCoverage[i_i] <- r[[i]]$filtered$nuclear$addn.stats["genome_cov"]
    o$AlnNoMAPQ[i_i] <- r[[i]]$filtered$nuclear$addn.stats["alignments_without_mapq"]
  }
  o
}

get_nucl_stats <- function(r, sn) {
  nRow <- 6
  o <- data.frame(row.names = NULL,
    Sample = rep(sn, each = nRow),
    Reads = rep(c("Min", "1stPctile", "Average", "SD", "99thPctile", "Max"),
      length(sn)),
    AlignmentSize = NA,
    FragmentSize = NA,
    MAPQ = NA,
    ReadDepth = NA,
    GCPercent = NA
  )
  for (i in seq_along(r)) {
    i_i <- (1:nRow) + nRow * (i - 1)
    o$AlignmentSize[i_i] <- r[[i]]$filtered$nuclear$stats[, "alignments"]
    o$FragmentSize[i_i] <- r[[i]]$filtered$nuclear$stats[, "fragments"]
    o$MAPQ[i_i] <- r[[i]]$filtered$nuclear$stats[, "mapq"]
    o$ReadDepth[i_i] <- r[[i]]$filtered$nuclear$stats[, "read_depth"]
    o$GCPercent[i_i] <- r[[i]]$filtered$nuclear$stats[, "gc_percent"]
  }
  o
}

get_overview_filt <- function(r, sn) {
  nRow <- 7
  o <- data.frame(row.names = NULL,
    Sample = rep(sn, each = nRow),
    Reads = rep(c("AlnPassFilters", "AlnSizeAvg", "ReadDepthAvg",
        "FragPassFilters", "FragSizeAvg", "MAPQAvg", "GCPctAvg"),
      length(sn)),
    Nuclear = NA,
    Mitochondrial = NA,
    Plastidic = NA
  )
  for (i in seq_along(r)) {
    i_i <- (1:nRow) + nRow * (i - 1)
    o$Nuclear[i_i] <- r[[i]]$filtered$overview[, "nuclear"]
    o$Mitochondrial[i_i] <- r[[i]]$filtered$overview[, "mitochondrial"]
    o$Plastidic[i_i] <- r[[i]]$filtered$overview[, "plastidic"]
  }
  o
}

get_overview_unfilt <- function(r, sn) {
  nRow <- 9
  o <- data.frame(row.names = NULL,
    Sample = rep(sn, each = nRow),
    Reads = rep(c("Total", "Duplicated", "SE", "PE", "ProperMate",
        "Primary", "PrimaryDuplicated", "Secondary", "Supplementary"),
      length(sn)),
    Nuclear = NA_integer_,
    Mitochondrial = NA_integer_,
    Plastidic = NA_integer_
  )
  for (i in seq_along(r)) {
    i_i <- (1:nRow) + nRow * (i - 1)
    o$Nuclear[i_i] <- r[[i]]$unfiltered$effective[, "nuclear"]
    o$Mitochondrial[i_i] <- r[[i]]$unfiltered$effective[, "mitochondrial"]
    o$Plastidic[i_i] <- r[[i]]$unfiltered$effective[, "plastidic"]
  }
  o
}

get_bam_stats <- function(r, sn) {
  nRow <- 4
  o <- data.frame(row.names = NULL,
    Sample = rep(sn, each = nRow),
    Reads = rep(c("Total", "Unmapped", "Mapped", "Blacklisted"),
      length(sn)),
    Count = NA_integer_
  )
  for (i in seq_along(r)) {
    i_i <- (1:nRow) + nRow * (i - 1)
    o$Count[i_i] <- r[[i]]$unfiltered$total
  }
  o
}

get_genome <- function(r, sn) {
  nRow <- 5
  o <- data.frame(row.names = NULL,
    Sample = rep(sn, each = nRow),
    Sequence = rep(c("total", "nuclear", "mitochondrial", "plastidic", "effective"),
      length(sn)),
    Count = NA_integer_,
    Size = NA_integer_
  )
  for (i in seq_along(r)) {
    i_i <- (1:nRow) + nRow * (i - 1)
    o$Count[i_i] <- r[[i]]$genome[, "sequences"]
    o$Size[i_i] <- r[[i]]$genome[, "size"]
  }
  o
}

