## Submission

This is a resubmission of quaqcr addressing the reviewer feedback on version
1.0.3:

* Removed the redundant "in R" from the package title.
* Expanded the acronyms in the Description field (ATAC-seq, NGS, ChIP-seq)
  and added a reference to the methods paper in the form
  `<doi:10.1093/bioinformatics/btae649>`.
* Replaced the `cat()` call in `R/quaqc.R` that printed the program help
  with a `message()` call, so the output can be suppressed via
  `suppressMessages()`. All other console output in the package already
  uses `message()` or is confined to S3 `print()` methods.
* Replaced `\dontrun{}` with `\donttest{}` in the `footprint()` and
  `pileup()` examples. Because both functions also require a user-supplied
  BAM file that the package cannot bundle, the example body is guarded by
  `if (nzchar(Sys.which("quaqc")) && file.exists(bam))` so that
  `R CMD check --run-donttest` runs the example as a no-op and does not
  error on machines without the binary or sample data.

## Test environments

* local macOS, R release
* GitHub Actions:
  * macOS-latest, R release
  * windows-latest, R release
  * windows-latest, R oldrel-4
  * ubuntu-latest, R devel
  * ubuntu-latest, R release
  * ubuntu-latest, R oldrel-1, oldrel-2, oldrel-3, oldrel-4
  * ubuntu-latest, R 3.6.3 (the declared minimum)
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## SystemRequirements

The package wraps the external command-line program `quaqc`
(<https://github.com/bjmt/quaqc>), declared in `SystemRequirements`. The
program is not bundled with the package and must be installed separately by
the user.

* The `quaqc()` example is wrapped in
  `\donttest{ if (nzchar(Sys.which("quaqc"))) { ... } }` so it runs under
  `R CMD check --run-donttest` only on machines where the program is on
  `PATH`.
* The `footprint()` and `pileup()` examples are wrapped in
  `\donttest{ if (nzchar(Sys.which("quaqc")) && file.exists(bam)) { ... } }`
  because, in addition to needing the `quaqc` binary, they require a
  user-supplied BAM file the package cannot bundle. The guard makes the
  example a no-op under `R CMD check --run-donttest` on machines without
  the binary or sample data.
* The test suite uses the same `nzchar(Sys.which("quaqc"))` guard so it
  silently skips when the program is unavailable.

## Downstream dependencies

There are currently no downstream dependencies.
