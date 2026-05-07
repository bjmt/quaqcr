.onLoad <- function(libname, pkgname) {
  if (is.null(getOption("quaqc.bin"))) {
    options(quaqc.bin = "quaqc")
  }
  invisible()
}
