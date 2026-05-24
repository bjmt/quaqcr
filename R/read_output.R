#' Read quaqc output files associated with a report.
#'
#' These helpers locate the per-sample output files that \code{quaqc}
#' produced for a given run and read them into \code{data.frame}s.
#' Output paths are reconstructed from each sample's BAM filename
#' (stored in the report) combined with the directory and extension
#' that \code{quaqc} would have used.
#'
#' @param report A \code{quaqc} object (as returned by [quaqc()] or
#' [parse_quaqc_file()]) or a single \code{quaqc_report} object.
#' @param dir Optional directory containing the output files. If
#' \code{NULL} (the default), the directory of each input BAM is used,
#' matching the default behaviour of \code{quaqc} when no \code{*-dir}
#' option is provided.
#' @param ext Filename extension that \code{quaqc} appended to each
#' sample's basename (after stripping the \code{.bam}/\code{.cram}
#' suffix). The defaults match those of \code{quaqc}.
#' @param file Path to the quantification TSV produced by
#' \code{quaqc --quant}. This file path is not stored in the JSON
#' report, so it must be supplied explicitly.
#' @param reader Function used to read each file. The default
#' [utils::read.table()] handles gzip and BGZF transparently. Replace
#' with e.g. \code{data.table::fread} for large files.
#'
#' @return When \code{report} is a \code{quaqc} object, a named list of
#' \code{data.frame}s (one per successful sample, named by the BAM
#' basename without the \code{.bam} suffix). When \code{report} is a
#' single \code{quaqc_report}, a single \code{data.frame}.
#' [read_quant()] always returns a single \code{data.frame}.
#'
#' @details
#' Column names are assigned based on each format:
#'
#' * `read_bedgraph` / `read_qscore`: \code{chrom}, \code{start},
#'   \code{end}, \code{value} (or \code{qscore}).
#' * `read_bed`: BED6 by default (\code{chrom}, \code{start},
#'   \code{end}, \code{name}, \code{score}, \code{strand}), or BED3
#'   (\code{chrom}, \code{start}, \code{end}) when \code{--bed-ins}
#'   was used. The format is detected from the boolean parameters in
#'   the report.
#' * `read_narrowpeak`: standard narrowPeak columns
#'   (\code{chrom}, \code{start}, \code{end}, \code{name},
#'   \code{score}, \code{strand}, \code{signalValue}, \code{pValue},
#'   \code{qValue}, \code{peak}).
#' * `read_quant`: read directly from the TSV header.
#'
#' Reading the filtered BAM output (\code{--keep}) is not supported
#' here; use a dedicated BAM-reading package such as \pkg{Rsamtools}.
#'
#' @examples
#' \donttest{
#' bam <- "Sample.bam"
#' if (nzchar(Sys.which("quaqc")) && file.exists(bam)) {
#'   res <- quaqc(bam, bedgraph = TRUE, bed = TRUE, call.peaks = TRUE)
#'   bgs <- read_bedgraph(res)
#'   beds <- read_bed(res)
#'   peaks <- read_narrowpeak(res)
#' }
#' }
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjmtremblay@@gmail.com}
#' @seealso [quaqc()], [parse_quaqc_file()]
#' @name read_quaqc_output
NULL

# Build the output path for a single sample.
quaqc_output_path <- function(sample, dir, ext) {
  base <- sub("\\.(bam|cram)$", "", basename(sample), ignore.case = TRUE)
  d <- if (is.null(dir)) dirname(sample) else dir
  if (!nzchar(d)) d <- "."
  file.path(d, paste0(base, ext))
}

# Apply a per-sample reader across a quaqc or quaqc_report object.
read_quaqc_outputs <- function(report, dir, ext, reader_fn, label) {
  if (is(report, "quaqc_report")) {
    if (!report$success) {
      stop("Sample '", report$sample, "' did not complete successfully")
    }
    path <- quaqc_output_path(report$sample, dir, ext)
    if (!file.exists(path)) {
      stop("Cannot find ", label, " file for sample '", report$sample,
        "' at '", path, "'")
    }
    return(reader_fn(path, report))
  }
  if (!is(report, "quaqc")) {
    stop("'report' must be a 'quaqc' or 'quaqc_report' object")
  }
  reports <- report$reports
  ok <- vapply(reports, function(x) isTRUE(x$success), logical(1))
  if (!any(ok)) stop("No successful samples found in report")
  if (!all(ok))
    message("Skipping ", sum(!ok), " failed report(s).")
  reports <- reports[ok]
  paths <- vapply(reports, function(r) quaqc_output_path(r$sample, dir, ext),
    character(1))
  exists <- file.exists(paths)
  if (!all(exists)) {
    missing <- vapply(reports[!exists], function(r) r$sample, character(1))
    warning("Missing ", label, " file(s) for: ",
      paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- vector("list", sum(exists))
  reports <- reports[exists]
  paths <- paths[exists]
  for (i in seq_along(paths)) {
    out[[i]] <- reader_fn(paths[i], reports[[i]])
  }
  names(out) <- vapply(reports,
    function(r) sub("\\.(bam|cram)$", "", basename(r$sample), ignore.case = TRUE),
    character(1))
  out
}

#' @rdname read_quaqc_output
#' @export
read_bedgraph <- function(report, dir = NULL, ext = ".bedGraph.gz",
  reader = utils::read.table) {
  reader_fn <- function(path, r) {
    df <- reader(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
      col.names = c("chrom", "start", "end", "value"),
      colClasses = c("character", "integer", "integer", "numeric"))
    as.data.frame(df)
  }
  read_quaqc_outputs(report, dir, ext, reader_fn, "bedGraph")
}

#' @rdname read_quaqc_output
#' @export
read_bed <- function(report, dir = NULL, ext = ".bed.gz",
  reader = utils::read.table) {
  reader_fn <- function(path, r) {
    ins <- isTRUE(unname(r$params$boolean["bed_ins"]) == TRUE) ||
           isTRUE(unname(r$params$boolean["bed_ins"]) == 1)
    if (ins) {
      df <- reader(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
        col.names = c("chrom", "start", "end"),
        colClasses = c("character", "integer", "integer"))
    } else {
      df <- reader(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
        col.names = c("chrom", "start", "end", "name", "score", "strand"),
        colClasses = c("character", "integer", "integer",
          "character", "numeric", "character"))
    }
    as.data.frame(df)
  }
  read_quaqc_outputs(report, dir, ext, reader_fn, "bed")
}

#' @rdname read_quaqc_output
#' @export
read_narrowpeak <- function(report, dir = NULL, ext = ".narrowPeak.gz",
  reader = utils::read.table) {
  reader_fn <- function(path, r) {
    df <- reader(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
      col.names = c("chrom", "start", "end", "name", "score", "strand",
        "signalValue", "pValue", "qValue", "peak"),
      colClasses = c("character", "integer", "integer", "character",
        "numeric", "character", "numeric", "numeric", "numeric", "integer"))
    as.data.frame(df)
  }
  read_quaqc_outputs(report, dir, ext, reader_fn, "narrowPeak")
}

#' @rdname read_quaqc_output
#' @export
read_qscore <- function(report, dir = NULL, ext = ".qscore.bedGraph.gz",
  reader = utils::read.table) {
  reader_fn <- function(path, r) {
    df <- reader(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
      col.names = c("chrom", "start", "end", "qscore"),
      colClasses = c("character", "integer", "integer", "numeric"))
    as.data.frame(df)
  }
  read_quaqc_outputs(report, dir, ext, reader_fn, "qscore bedGraph")
}

#' @rdname read_quaqc_output
#' @export
read_quant <- function(report, file, reader = utils::read.table) {
  if (!is(report, "quaqc") && !is(report, "quaqc_report"))
    stop("'report' must be a 'quaqc' or 'quaqc_report' object")
  if (missing(file) || is.null(file))
    stop("'file' must be supplied; the quant TSV path is not stored in the JSON report")
  file <- as.character(file)[1]
  if (!file.exists(file))
    stop("Cannot find quant file at '", file, "'")
  df <- reader(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
    check.names = FALSE)
  as.data.frame(df)
}
