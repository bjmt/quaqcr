#' Run quaqc from within R.
#'
#' Interactive wrapper for running quaqc from within R. For a detailed
#' description of the program, see the manual: execute \code{man quaqc} from the
#' command line or open `doc/quaqc.1.md` in the program folder. For a brief
#' description of the command parameters, as well as to see default values,
#' call `quaqc()` without any arguments.
#'
#' @param bam.files Character vector of BAM file names. Must be coordinate
#' sorted. If no index file can be found, `quaqc` will generate them.
#' @param mitochondria Character vector of mitochondria names. Provide `""`
#' to clear the defaults.
#' @param plastids Character vector of plastid names. Provide `""` to clear
#' the defaults.
#' @param peaks Either (1) a filename of a BED file containing peaks, or (2) a
#' `GRanges` object from the \pkg{GenomicRanges} package.
#' @param tss Either (1) a filename of a BED file containing TSSs, or (2) a
#' `GRanges` object from the \pkg{GenomicRanges} package.
#' @param target.names Character vector of sequence names to restrict `quaqc`.
#' @param target.list Either (1) a filename of a BED file containing ranges
#' to restrict `quaqc`, or (2) a `GRanges` object from the \pkg{GenomicRanges}
#' package.
#' @param blacklist Either (1) a filename of a BED file containing blacklist
#' ranges, or (2) a `GRanges` object from the \pkg{GenomicRanges}
#' package.
#' @param rg.names A character vector of read group (RG) names to restrict
#' `quaqc`.
#' @param rg.list Filename of a text file containing read group (RG) names to
#' restrict `quaqc`, one name per line.
#' @param use.secondary Logical, allow secondary alignments.
#' @param use.nomate Logical, allow PE reads when the mate does not align
#' properly.
#' @param use.dups Logical, allow duplicate reads.
#' @param use.chimeric Logical, allow supplemental or chimeric alignments.
#' @param use.dovetails Logical, allow dovetailing PE reads.
#' @param no.se Logical, discard SE reads.
#' @param mapq Integer, min MAPQ score.
#' @param min.qlen Integer, min alignment length.
#' @param min.flen Integer, min fragment length.
#' @param max.qlen Integer, max alignment length.
#' @param max.flen Integer, max fragment length.
#' @param use.all Logical, discard all filters and keep all reads.
#' @param max.depth Integer, max base depth for read depth histogram.
#' @param max.qhist Integer, max alignment length for histogram.
#' @param max.fhist Integer, max fragment length for histogram.
#' @param tss.size Integer, size of the TSS region for pileup.
#' @param tss.qlen Integer, resize reads (centered on the 5-prime end for pileup.
#' @param tss.tn5 Logical, shift 5-prime end coordinates +4/-5 bases for pileup.
#' @param omit.gc Logical, omit calculation of read GC content.
#' @param omit.depth Logical, omit calculation of read depths.
#' @param fast Logical, turn on fast mode.
#' @param lenient Logical, turn on lenient mode.
#' @param nfr Logical, turn on NFR mode.
#' @param nbr Logical, turn on NBR mode.
#' @param footprint Logical, turn on footprinting mode.
#' @param chip Logical, turn on ChIP-seq mode.
#' @param output.dir Name of directory to save QC report if not that of input.
#' @param output.ext Filename extension for output files.
#' @param no.output Logical, suppress creation of output QC reports. Note that
#' option is turned on by default when run from \pkg{quaqcr}.
#' @param json Filename of JSON file to save combined QC results to. Set to
#' `NULL` to suppress this. The default is to pipe the JSON output directly to
#' R and not save to a file.
#' @param keep Logical, save passing nuclear reads to a new BAM file.
#' @param keep.dir Directory name to save filtered BAMs.
#' @param keep.ext Extension of filtered BAMs.
#' @param threads Integer, number of worker threads. Max one per sample.
#' @param title Assign a title to run.
#' @param continue Logical, do not return an error and instead continue running
#' if samples trigger program errors.
#' @param verbose Integer, a value from 0 to 2 for the level of program verbosity.
#' @param timeout Integer, number of seconds before stopping quaqc. By default
#' it is allowed to run indefinitely. See [base::system2()].
#' @param env A character vector of environment variables to set when
#' running quaqc. See [base::system2()].
#' @param stderr.file Filename to save quaqc messages. By default they are
#' printed in the console.
#' @param bin Path to quaqc binary. Alternatively, set `options(quaqc.bin)`.
#' If the binary is present in the current working directory, provide `"./quaqc"`.
#' The default, set when \pkg{quaqcr} is loaded, is to assume the binary is
#' present in your PATH.
#'
#' @return If nothing is provided to `bam.files`, then the help message is
#' printed to the console and returned as a character vector, invisibly.
#' Alternatively: if `json = NULL` then `NULL`, otherwise the JSON
#' output as parsed by \pkg{jsonlite}.
#'
#' @examples
#' \dontrun{
#' ## To check that you are properly linking to the binary and view help:
#' quaqc()
#' }
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjmtremblay@@gmail.com}
#' @seealso [base::system2()]
#' @export
quaqc <- function(bam.files, mitochondria = NULL, plastids = NULL, peaks = NULL,
  tss = NULL, target.names = NULL, target.list = NULL, blacklist = NULL,
  rg.names = NULL, rg.list = NULL, use.secondary = FALSE, use.nomate = FALSE,
  use.dups = FALSE, use.chimeric = FALSE, use.dovetails = FALSE, no.se = FALSE,
  mapq = NULL, min.qlen = NULL, min.flen = NULL, max.qlen = NULL,
  max.flen = NULL, use.all = FALSE, max.depth = NULL, max.qhist = NULL,
  max.fhist = NULL, tss.size = NULL, tss.qlen = NULL, tss.tn5 = FALSE,
  omit.gc = FALSE, omit.depth = FALSE, fast = FALSE, lenient = FALSE,
  nfr = FALSE, nbr = FALSE, footprint = FALSE, chip = FALSE, output.dir = NULL,
  output.ext = NULL, no.output = TRUE, json = "-", keep = FALSE,
  keep.dir = NULL, keep.ext = NULL, threads = NULL, title = NULL, continue = FALSE,
  verbose = 1, timeout = 0, env = character(), stderr.file = "",
  bin = getOption("quaqc.bin")) {

  if (is.null(bin)) stop("Please specify the location of the quaqc binary")

  run_quaqc <- function(x) {
    res <- try(suppressWarnings(
      system2(path.expand(bin), args = x, stdout = TRUE, stderr = stderr.file,
        env = env, timeout = timeout)), silent = TRUE)
    if (!is.null(attr(res, "status")) && attr(res, "status") > 0) {
      stop("quaqc encountered an error while running, see stderr")
    } else if (is(res, "try-error")) {
      stop("Failed to execute quaqc program, make sure `bin` is correct")
    }
    res
  }

  if (missing(bam.files)) {
    help <- run_quaqc("-h")
    cat(help, sep = "\n")
    return(invisible(c(version, help)))
  }

  args <- as.character(bam.files)

  # BED/GRanges args

  if (!is.null(peaks)) {
    if (is.character(peaks)) {
      args <- c("--peaks", normalizePath(path.expand(peaks), mustWork = TRUE)[1], args)
    } else if (is(peaks, "GRanges")) {
      peaks <- as.data.frame(peaks)
      peaks <- cbind(peaks$seqnames, peaks$start, peaks$end, ".", 1, peaks$strand)
      peaks.tmp <- tempfile("peaks", fileext = ".bed")
      on.exit(unlink(peaks.tmp))
      write.table(peaks, peaks.tmp, sep = "\t", row.names = FALSE, col.names = FALSE,
        quote = FALSE)
      args <- c("--peaks", peaks.tmp, args)
    } else {
      stop("Expected a filename or GRanges object for peaks")
    }
  }
  if (!is.null(tss)) {
    if (is.character(tss)) {
      args <- c("--tss", normalizePath(path.expand(tss), mustWork = TRUE)[1], args)
    } else if (is(tss, "GRanges")) {
      tss <- as.data.frame(tss)
      tss <- cbind(tss$seqnames, tss$start, tss$end, ".", 1, tss$strand)
      tss.tmp <- tempfile("tss", fileext = ".bed")
      on.exit(unlink(tss.tmp))
      write.table(tss, tss.tmp, sep = "\t", row.names = FALSE, col.names = FALSE,
        quote = FALSE)
      args <- c("--tss", tss.tmp, args)
    } else {
      stop("Expected a filename or GRanges object for tss")
    }
  }
  if (!is.null(target.list)) {
    if (is.character(target.list)) {
      args <- c("--target-list", normalizePath(path.expand(target.list), mustWork = TRUE)[1], args)
    } else if (is(target.list, "GRanges")) {
      target.list <- as.data.frame(target.list)
      target.list <- cbind(target.list$seqnames, target.list$start, target.list$end, ".", 1, target.list$strand)
      target.list.tmp <- tempfile("target.list", fileext = ".bed")
      on.exit(unlink(target.list.tmp))
      write.table(target.list, target.list.tmp, sep = "\t", row.names = FALSE, col.names = FALSE,
        quote = FALSE)
      args <- c("--target-list", target.list.tmp, args)
    } else {
      stop("Expected a filename or GRanges object for target.list")
    }
  }
  if (!is.null(blacklist)) {
    if (is.character(blacklist)) {
      args <- c("--blacklist", normalizePath(path.expand(blacklist), mustWork = TRUE)[1], args)
    } else if (is(blacklist, "GRanges")) {
      blacklist <- as.data.frame(blacklist)
      blacklist <- cbind(blacklist$seqnames, blacklist$start, blacklist$end, ".", 1, blacklist$strand)
      blacklist.tmp <- tempfile("blacklist", fileext = ".bed")
      on.exit(unlink(blacklist.tmp))
      write.table(blacklist, blacklist.tmp, sep = "\t", row.names = FALSE, col.names = FALSE,
        quote = FALSE)
      args <- c("--blacklist", blacklist.tmp, args)
    } else {
      stop("Expected a filename or GRanges object for blacklist")
    }
  }

  # Logical args

  if (use.secondary) args <- c("--use-secondary", args)
  if (use.nomate) args <- c("--use-nomate", args)
  if (use.dups) args <- c("--use-dups", args)
  if (use.chimeric) args <- c("--use-chimeric", args)
  if (use.dovetails) args <- c("--use-dovetails", args)
  if (no.se) args <- c("--no-se", args)
  if (use.all) args <- c("--use-all", args)
  if (tss.tn5) args <- c("--tss-tn5", args)
  if (omit.gc) args <- c("--omit-gc", args)
  if (omit.depth) args <- c("--omit-depth", args)
  if (fast) args <- c("--fast", args)
  if (lenient) args <- c("--lenient", args)
  if (nfr) args <- c("--nfr", args)
  if (nbr) args <- c("--nbr", args)
  if (footprint) args <- c("--footprint", args)
  if (chip) args <- c("--chip", args)
  if (no.output) args <- c("--no-output", args)
  if (keep) args <- c("--keep", args)
  if (continue) args <- c("--continue", args)

  # Character args

  if (!is.null(mitochondria)) args <- c("--mitochondria", paste0(mitochondria, collapse = ","), args)
  if (!is.null(plastids)) args <- c("--plastids", paste0(plastids, collapse = ","), args)
  if (!is.null(target.names)) args <- c("--target-names", paste0(target.names, collapse = ","), args)
  if (!is.null(rg.names)) args <- c("--rg-names", paste0(rg.names, collapse = ","), args)
  if (!is.null(rg.list)) args <- c("--rg-list", normalizePath(path.expand(rg.list), mustWork = TRUE)[1], args)
  if (!is.null(output.dir)) args <- c("--output-dir", normalizePath(path.expand(output.dir), mustWork = TRUE)[1], args)
  if (!is.null(output.ext)) args <- c("--output-ext", output.ext[1], args)
  if (!is.null(keep.dir)) args <- c("--keep-dir", normalizePath(path.expand(keep.dir), mustWork = TRUE)[1], args)
  if (!is.null(keep.ext)) args <- c("--keep-ext", keep.ext[1], args)
  if (!is.null(title)) args <- c("--title", paste0(title, collapse = " "), args)

  isPiped <- FALSE
  if (!is.null(json)) {
    json <- as.character(json)[1]
    if (json == "-" || json == "/dev/stdin") isPiped <- TRUE
    args <- c("--json", json, args)
  }

  # Integer args

  if (!is.null(mapq)) args <- c("--mapq", as.integer(mapq)[1], args)
  if (!is.null(min.qlen)) args <- c("--min-qlen", as.integer(min.qlen)[1], args)
  if (!is.null(min.flen)) args <- c("--min-flen", as.integer(min.flen)[1], args)
  if (!is.null(max.qlen)) args <- c("--max-qlen", as.integer(max.qlen)[1], args)
  if (!is.null(max.flen)) args <- c("--max-flen", as.integer(max.flen)[1], args)
  if (!is.null(max.depth)) args <- c("--max-depth", as.integer(max.depth)[1], args)
  if (!is.null(max.qhist)) args <- c("--max-qhist", as.integer(max.qhist)[1], args)
  if (!is.null(max.fhist)) args <- c("--max-fhist", as.integer(max.fhist)[1], args)
  if (!is.null(tss.size)) args <- c("--tss-size", as.integer(tss.size)[1], args)
  if (!is.null(tss.qlen)) args <- c("--tss-qlen", as.integer(tss.qlen)[1], args)
  if (!is.null(threads)) args <- c("--threads", as.integer(threads)[1], args)

  # Program args

  if (verbose) {
    if (verbose == 1) args <- c("-v", args)
    else args <- c("-vv", args)
  }

  message("quaqc ", paste0(args, collapse = " "))
  res <- run_quaqc(args)

  if (isPiped) {
    parse_quaqc(fromJSON(res, simplifyDataFrame = FALSE))
  } else if (!is.null(json)) {
    f <- gzfile(normalizePath(path.expand(json))[1], "rt")
    parse_quaqc(fromJSON(readLines(f), simplifyDataFrame = FALSE))
    close(f)
  } else {
    invisible(NULL)
  }

}

# extfun <- function(PKG, FUN, env = parent.frame()) {
#   if (requireNamespace(PKG, quietly = TRUE)) {
#     eval(substitute(FUN), envir = env)
#   } else {
#     stop("The ", PKG, " package must be installed", call. = FALSE)
#   }
# }

