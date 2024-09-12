#' Parse a JSON quaqc report file.
#'
#' Parse the output of \code{quaqc --json} into a easier to use
#' \code{quaqc}-class object within R.
#'
#' @param json.file The JSON report filename.
#' @param json.text The JSON report as a character vector.
#'
#' @return A \code{quaqc}-class object.
#'
#' @details
#'
#' A \code{quaqc} object is a higher level format encompassing the
#' quaqc run parameters (accessible via `$metadata`) and the actual
#' individual reports for each sample (accessible via `$reports`). The
#' reports are themselves \code{quaqc_report}-class objects with
#' multiple list slots, including:
#'
#' * `$sample`: The filename of the sample.
#' * `$success`: Whether the sample was successfully analyzed.
#' * `$params`: Values for all quaqc parameters used to analyze this sample.
#' * `$genome`: Data about the genome taken from the BAM header.
#' * `$unfiltered`: Basic stats about the total number of reads before filtering.
#' * `$filtered`: Contains the majority of the data output by quaqc.
#'
#' This final `$filtered` slot itself is broken down into several sub-lists:
#'
#' * `$overview`: Average values for several stats such as fragment size.
#' * `$nuclear$stats`: Further breakdown of the previous stats for nuclear reads.
#' * `$nuclear$stats.warn`: Whether any quaqc parameters prevented it from accurately collecting some data.
#' * `$nuclear$addn.stats`: Genome coverage and the number of alignments without a MAPQ score.
#' * `$nuclear$histograms`: Raw histogram data for alignment size, fragment size, GC percent, and read depth.
#' * `$nuclear$peaks`: The number of peaks, the fraction of the effective genome covered by them, and the FRIP score.
#' * `$nuclear$tss`: The read pileup around TSSs as well as the TSS enrichment score.
#' 
#' Note that the word 'effective' refers to reads which are visible to
#' quaqc within target regions or outside blacklisted regions, as well
#' as reads associated with any specified target read groups. 
#'
#' @examples
#' report.file <- system.file("extdata", "report.json.gz", package = "quaqcr")
#'
#' ## Option 1: parse a report already read into R
#' f <- gzfile(report.file, "rt")
#' json <- jsonlite::fromJSON(readLines(f), simplifyDataFrame = FALSE)
#' close(f)
#' report <- parse_quaqc(json)
#'
#' ## Option 2: parse a report directly from a file
#' report <- parse_quaqc_file(report.file)
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjmtremblay@@gmail.com}
#' @seealso [quaqcr::quaqc()]
#' @rdname parse_quaqc
#' @export
parse_quaqc <- function(json.text) {
  json <- json.text
  stopifnot(is.list(json))
  stopifnot(names(json) %in% c("quaqc_version", "quaqc_run_title", "quaqc_args",
    "quaqc_time_start", "quaqc_params", "quaqc_reports", "quaqc_time_end",
    "quaqc_max_bytes"))
  meta <- list(
    version = json$quaqc_version,
    title = json$quaqc_run_title,
    date = json$quaqc_time_end,
    args = json$quaqc_args,
    samples = json$quaqc_params$sample_n
    # runtime = as.integer(difftime(json$quaqc_time_end, json$quaqc_time_start, units = "secs")),
    # memory = json$quaqc_max_bytes,
  )
  params.nume <- unlist(json$quaqc_params[c(
      "target_seqs_n", "peak_bed_n", "tss_bed_n", "target_list_bed_n",
      "blacklist_bed_n", "read_groups_n", "mapq_min", "alignment_size_min",
      "alignment_size_max", "fragment_size_min", "fragment_size_max",
      "alignment_histogram_max", "fagment_histogram_max",
      "depth_histogram_max", "tss_histogram_max", "tss_histogram_size",
      "tss_read_size"
  )])
  params.logi <- unlist(json$quaqc_params[c(
    "tss_tn5_shift", "use_secondary_alignments", "use_supplementary_alignments",
    "use_improper_mates", "use_duplicates", "no_se", "use_dovetails", "use_all"
  )])
  reports <- lapply(json$quaqc_reports, function(x) parse_single_report(x, params.nume, params.logi))
  structure(list(metadata = meta, reports = reports), class = "quaqc")
}

#' @rdname parse_quaqc
#' @export
parse_quaqc_file <- function(json.file) {
  json <- json.file
  if (!is.character(json) && length(json) != 1) {
    stop("'json.file' must be a length 1 character vector")
  }
  f <- gzfile(normalizePath(path.expand(json))[1], "rt")
  json <- fromJSON(readLines(f), simplifyDataFrame = FALSE)
  close(f)
  parse_quaqc(json)
}

parse_single_report <- function(report, params.nume, params.logi) {
  if (!report$status_success) {
    structure(list(
      sample = report$sample,
      success = FALSE,
      params = list(integer = params.nume, boolean = params.logi),
      genome = NULL,
      unfiltered = NULL,
      filtered = NULL
    ), class = c("quaqc_report"))
  } else {
    sample <- report$sample
    report <- report$report
    genome <- list(
      total = c(
        sequences = report$genome_stats$total$genome_wide$seq_n,
        size = report$genome_stats$total$genome_wide$seq_size
      ),
      nuclear = c(
        sequences = report$genome_stats$total$nuclear$seq_n,
        size = report$genome_stats$total$nuclear$seq_size
      ),
      mitochondrial = c(
        sequences = report$genome_stats$total$mitochondrial$seq_n,
        size = report$genome_stats$total$mitochondrial$seq_size
      ),
      plastidic = c(
        sequences = report$genome_stats$total$plastidic$seq_n,
        size = report$genome_stats$total$plastidic$seq_size
      ),
      effective = c(
        sequences = report$genome_stats$effective$nuclear$seq_n,
        size = report$genome_stats$effective$nuclear$seq_size
      )
    )
    genome <- do.call(rbind, genome)
    unfiltered <- list(total = unlist(report$unfiltered_read_stats$total))
    unfiltered.stats <- list(
      nuclear = unlist(report$unfiltered_read_stats$effective$nuclear),
      mitochondrial = unlist(report$unfiltered_read_stats$effective$mitochondrial),
      plastidic = unlist(report$unfiltered_read_stats$effective$plastidic)
    )
    for (i in seq_along(unfiltered.stats)) {
      unfiltered.stats[[i]] <- unfiltered.stats[[i]][
        names(unfiltered.stats[[i]]) != "non_redundant_fraction"
      ]
    }
    unfiltered$effective <- do.call(cbind, unfiltered.stats)
    overview.nuc <- c(
      alignments_passing_filters = report$filtered_read_stats$nuclear$alignment$passing_filters,
      alignment_size_average = report$filtered_read_stats$nuclear$alignment$size_average,
      read_depth_average = nullToNA(report$filtered_read_stats$nuclear$depth$read_depth_average),
      fragments_passing_filters = nullToNA(report$filtered_read_stats$nuclear$fragment$passing_filters),
      fragment_size_average = nullToNA(report$filtered_read_stats$nuclear$fragment$size_average),
      mapq_score_average = nullToNA(report$filtered_read_stats$nuclear$mapq$score_average),
      gc_pct_average = nullToNA(report$filtered_read_stats$nuclear$gc$pct_average)
    )
    filtered <- list(
      nuclear = overview.nuc,
      mitochondrial = unlist(report$filtered_read_stats$mitochondrial),
      plastidic = unlist(report$filtered_read_stats$plastidic)
    )
    nuclear.stats.names1 <- c("size_min", "size_1st_pctile", "size_average",
      "size_sd", "size_99th_pctile", "size_max")
    nuclear.stats.names2 <- c("score_min", "score_1st_pctile", "score_average",
      "score_sd", "score_99th_pctile", "score_max")
    nuclear.stats.names3 <- c("depths_min", "depths_1st_pctile", "read_depth_average",
      "depths_sd", "depths_99th_pctile", "depths_max")
    nuclear.stats.names4 <- c("pct_min", "pct_1st_pctile", "pct_average",
      "pct_sd", "pct_99th_pctile", "pct_max")
    nuclear.stats <- data.frame(check.names = FALSE,
      row.names = c("min", "1st_pctile", "average", "sd", "99th_pctile", "max"),
      alignments = unlist(report$filtered_read_stats$nuclear$alignment[nuclear.stats.names1]),
      fragments = nullToNA(unlist(report$filtered_read_stats$nuclear$fragment[nuclear.stats.names1])),
      mapq = unlist(report$filtered_read_stats$nuclear$mapq[nuclear.stats.names2]),
      read_depth = nullToNA(unlist(report$filtered_read_stats$nuclear$depth[nuclear.stats.names3])),
      gc_percent = nullToNA(unlist(report$filtered_read_stats$nuclear$gc[nuclear.stats.names4]))
    )
    nuclear.stats.warn = c(
      alignment = report$filtered_read_stats$nuclear$alignment$addn_stats_are_suspect,
      fragment = report$filtered_read_stats$nuclear$fragment$addn_stats_are_suspect,
      depth = report$filtered_read_stats$nuclear$depth$addn_stats_are_suspect
    )
    peaks <- if (is.null(report$filtered_read_stats$nuclear$peaks)) NULL else {
      unlist(report$filtered_read_stats$nuclear$peaks)
    }
    tss <- if (is.null(report$filtered_read_stats$nuclear$tss)) NULL else {
      list(
        stats = unlist(report$filtered_read_stats$nuclear$tss[c("n", "tss_enrichment_score")]),
        pileup = report$filtered_read_stats$nuclear$tss$tss_pileup
      )
    }
    nuclear <- list(
      stats = nuclear.stats,
      stats.warn = nuclear.stats.warn,
      addn.stats = c("genome_cov" = report$filtered_read_stats$nuclear$depth$genome_cov,
        "alignments_without_mapq" = report$filtered_read_stats$nuclear$mapq$no_mapq_n),
      histograms = list(
        alignment = report$filtered_read_stats$nuclear$alignment$size_histogram,
        fragment = report$filtered_read_stats$nuclear$fragment$size_histogram,
        mapq = report$filtered_read_stats$nuclear$mapq$score_histogram,
        gc = report$filtered_read_stats$nuclear$gc$pct_histogram,
        depth = report$filtered_read_stats$nuclear$depth$depths_histogram
      ),
      peaks = peaks, tss = tss
    )
    filtered <- list(
      overview = do.call(cbind, filtered),
      nuclear = nuclear
    )
    structure(list(sample = sample, success = TRUE,
      params = list(integer = params.nume, boolean = params.logi),
      genome = genome, unfiltered = unfiltered,
      filtered = filtered), class = c("quaqc_report"))
  }
}

identical_params <- function(x, y) {
  if (!is(x, "quaqc_report") || !is(y, "quaqc_report"))
    stop("'x' and 'y' must be 'quaqc_report' objects")
  if (!validate_quaqc_report(x) || !validate_quaqc_report(y))
    warning("'quaqc_report' object(s) may be invalid", call. = FALSE, immediate. = TRUE)
  if (any(x$params$integer != y$params$integer) ||
      any(x$params$boolean != y$params$boolean)) {
    int_i <- which(x$params$integer != y$params$integer)
    int_b <- which(x$params$boolean != y$params$boolean)
    if (length(int_i) + length(int_b)) message("Found the following differences:")
    for (i in int_i) {
      message("  ", names(x$params$integer)[i], ": x=", x$params$integer[i], " y=", y$params$integer[i])
    }
    for (i in int_b) {
      message("  ", names(x$params$boolean)[i], ": x=", x$params$boolean[i], " y=", y$params$boolean[i])
    }
    FALSE
  } else {
    TRUE
  }
}

println <- function(...) cat(..., "\n", sep = "")

#' @export
print.quaqc <- function(x, ...) {
  if (!validate_quaqc(x))
    warning("'quaqc' object may be invalid", call. = FALSE, immediate. = TRUE)
  println("quaqc v", x$metadata$version)
  println("Run title: ", if (!nchar(x$metadata$title)) "---" else x$metadata$title)
  println("Run date: ", x$metadata$date)
  println("Number of samples: ", x$metadata$samples)
  println("Number of fails:   ", sum(vapply(x$reports, function(y) !y$success, logical(1))))
  println("Reports:")
  nmax <- length(x$reports)
  for (i in seq_len(min(nmax, 3))) {
    status <- if (x$reports[[i]]$success) "[SUCCESS]" else "[FAILURE]"
    println("    ", status, " ", x$reports[[i]]$sample)
  }
  if (nmax > 3) println("    ...")
  println("To examine run data: $metadata")
  println("To examine reports:  $reports")
  println("See ?melt_reports or ?parse_quaqc for ways to access data.")
  invisible(x)
}

#' @export
print.quaqc_report <- function(x, ...) {
  if (!validate_quaqc_report(x))
    warning("'quaqc_report' object may be invalid", call. = FALSE, immediate. = TRUE)
  println("Sample: ", x$sample)
  println("Status: ", if (x$success) "success" else "fail")
  if (x$success) {
    println("Processed ", x$genome["effective", "sequences"], " sequences (",
      x$genome["effective", "size"], " bp).")
    println("Reads: ", x$unfiltered$total["reads"], " total; ",
      x$filtered$overview["alignments_passing_filters", "nuclear"], " passing.")
    println(
      "GC histo: ", if (!is.null(x$filtered$nuclear$histograms$gc)) "yes" else "no", ". ",
      "Depths histo: ", if (!is.null(x$filtered$nuclear$histograms$depth)) "yes" else "no", ". ",
      "\nPeak stats: ", if (!is.null(x$filtered$nuclear$peaks)) "yes" else "no", ". ",
      "TSS pileup: ", if (!is.null(x$filtered$nuclear$tss)) "yes" else "no", "."
    )
  }
  println("Available slots:")
  println("$sample      -- Sample filename")
  println("$success     -- Run status")
  println("$params      -- Run parameters")
  println("$genome      -- Genome stats")
  println("$unfiltered  -- Pre-filter read stats")
  println("$filtered    -- Post-filter read stats")
  invisible(x)
}

validate_quaqc_report <- function(x) {
  if (!is(x, "quaqc_report")) stop("Object is not of class 'quaqc_report'")
  # TODO:
  TRUE
}

validate_quaqc <- function(x) {
  if (!is(x, "quaqc")) stop("Object is not of class 'quaqc'")
  # TODO:
  TRUE
}

remove_fails <- function(x) {
  if (!validate_quaqc(x))
    warning("'quaqc_report' object may be invalid", call. = FALSE, immediate. = TRUE)
  x$reports <- x$reports[sapply(x$reports, function(y) y$success)]
  x$metadata$samples <- length(x$reports)
  x
}

nullToNA <- function(x, n = 1) if (is.null(x)) rep(NA, n) else x

