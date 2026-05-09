## Submission

This is the first submission of quaqcr to CRAN.

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
* The `footprint()` and `pileup()` examples are wrapped in `\dontrun{}`
  because they additionally require a user-supplied BAM file (sequencing
  data the user must obtain themselves) and therefore cannot be executed
  in any check environment.
* The test suite uses the same `nzchar(Sys.which("quaqc"))` guard so it
  silently skips when the program is unavailable.

## Downstream dependencies

There are currently no downstream dependencies.
