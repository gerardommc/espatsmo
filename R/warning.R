.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    packageStartupMessage(
      "Package 'INLA' is necessary but is not installed.",
      "\n Please install it by running:",
      "\n install.packages(\"INLA\",repos=c(getOption(\"repos\"),INLA=\"https://inla.r-inla-download.org/R/stable\"), dep=TRUE)"
    )
  }
}