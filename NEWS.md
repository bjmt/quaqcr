# quaqcr 1.0.4

## CRAN submission preparation

* `DESCRIPTION`: removed redundant "in R" from the title.
* `DESCRIPTION`: expanded acronyms (ATAC-seq, NGS, ChIP-seq) and added a
  reference to the methods paper in the form
  `<doi:10.1093/bioinformatics/btae649>`.
* `quaqc()`: replaced the `cat()` call that printed the program help with a
  `message()` call so that it can be suppressed by the user.
* `footprint()` and `pileup()` examples: switched from `\dontrun{}` to a
  guarded `\donttest{}` that no-ops when the `quaqc` binary or a sample BAM
  is unavailable.

# quaqcr 1.0.3

## CRAN submission preparation

* Reverted the `footprint()` and `pileup()` examples to `\dontrun{}` because
  they additionally require a user-supplied BAM file that the user must obtain
  themselves; the previous `\donttest{}` wrapper caused
  `R CMD check --run-donttest` to fail. The `quaqc()` example continues to use
  `\donttest{ if (nzchar(Sys.which("quaqc"))) }`.
* Added `cran-comments.md` to `.Rbuildignore`.

# quaqcr 1.0.2

## CRAN submission preparation

* Added `SystemRequirements: quaqc` to DESCRIPTION.
* Moved `options(quaqc.bin = ...)` initialisation into `.onLoad()` in `R/zzz.R`.
* Migrated `inst/CITATION` from deprecated `citEntry()`/`personList()` to
  `bibentry()`/`c(person(), person())`.
* Updated `.Rbuildignore` to exclude development files (`CLAUDE.md`, `.claude/`).
* Removed empty `merge_reports()` stub (`R/merge.R`).
* Converted `\dontrun{}` to `\donttest{ if (nzchar(Sys.which("quaqc"))) }` in
  examples for `quaqc()`, `footprint()`, and `pileup()`.
* Fixed typo `genereal` → `general` in DESCRIPTION.

# quaqcr 1.0.1

## Bug fixes

* `parse_quaqc()`: input validation now checks that all required JSON keys are
  present and reports which are missing, rather than accepting any subset of
  valid names.
* `parse_quaqc()`: non-list input now produces a clear error message
  (`'json.text' must be a list`).
* `parse_quaqc_file()`: the length/type guard used `&&` instead of `||`,
  allowing numeric length-1 inputs or character length->1 inputs to bypass
  validation and crash inside `gzfile()`.
* `parse_quaqc()`: fixed typo `fagment_histogram_max` → `fragment_histogram_max`
  in the params key vector; the misspelled key was silently absent from
  `$params$integer` on every parsed report.
* `parse_quaqc()`: removed `tss_histogram_max` from the params key vector; this
  key is not emitted by quaqc (only `tss_histogram_size` is), so it was
  silently missing from every parsed report.
* `quaqc()`: calling `quaqc()` with no arguments (to print help) would error
  with `object 'version' not found`; fixed by returning just the help text.
* `quaqc()`: when `json` is a file path, the parsed report was discarded and
  `NULL` returned; the file connection was also leaked on error. Both issues
  are fixed.
* `quaqc()`: when `peaks`, `tss`, `target.list`, or `blacklist` was supplied as
  a `GRanges` object, the temporary BED file's `on.exit(unlink(...))` handler
  lacked `add = TRUE`; passing more than one GRanges argument would leak all
  but the last temp file.
* `quaqc()`: GRanges objects are 1-based inclusive; the BED files written from
  them now correctly use 0-based half-open coordinates (`start - 1L`).
* `quaqc()`: `verbose = NULL` caused an error (`argument is of length zero`);
  the verbosity block now guards against `NULL`.
* `melt_reports()`: calling with a `quaqc` object containing failed reports
  caused a crash on NULL subscript-assignment inside the `get_*` helper
  functions; failed reports are now filtered out automatically (with a message).
* `melt_reports()`: histogram and TSS-pileup helpers used `1:nRow[i]` to build
  row indices, which produces `c(1, 0)` when a sample has zero rows; replaced
  with `seq_len(nRow[i])` throughout.
* `melt_reports()`: TSS background normalisation could call `rep(0, n)` with a
  negative `n`; clamped to `max(0L, n)`.
* `melt_reports()`: dividing TSS depth by `nranges` when `nranges` is 0 or NA
  produced `Inf`/`NaN`; guarded with a fallback to 1.
* `melt_reports()`: `use.basename = TRUE` failed on R < 4.0 because
  `data.frame()` converts character columns to factors by default in those
  versions; the `Sample` column is now explicitly coerced to character before
  `basename()` is called.
* `pileup()`: passing `tss` in `...` previously issued a warning then crashed
  with a confusing "formal argument matched by multiple actual arguments" error
  from R; it now stops with a clear message.
* `footprint()`: the unreachable `tss`-in-`...` warning (R errors before the
  body runs due to ambiguous partial matching against `tss.size`/`tss.qlen`/
  `tss.tn5`) has been removed.

## Other

* Added a `testthat` test suite covering parsing, all 13 `melt_reports()`
  sections, validators, print snapshots, and input-validation paths.
