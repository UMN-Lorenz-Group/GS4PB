#' Install all GS4PB dependencies
#'
#' Convenience wrapper that sources the bundled dependency installer.
#' Run this once after installing the GS4PB package.  It installs CRAN,
#' Bioconductor, and GitHub packages, then sets up the Python virtualenv
#' for AlphaPlantImpute.
#'
#' @param cran_mirror Character. CRAN mirror URL.
#'   Default: \code{"https://cloud.r-project.org"}.
#' @param python_venv Character. Name of the reticulate virtualenv to create
#'   for AlphaPlantImpute.  Default: \code{"gs4pb-python"}.
#' @param install_python Logical. Whether to set up the Python virtualenv.
#'   Set \code{FALSE} if you do not need the AlphaPlantImpute imputation
#'   method.  Default: \code{TRUE}.
#'
#' @return Invisible \code{NULL}. Messages are printed to the console.
#'
#' @examples
#' \dontrun{
#' # Install everything (R + Python)
#' GS4PB::install_gs4pb_deps()
#'
#' # Skip Python setup
#' GS4PB::install_gs4pb_deps(install_python = FALSE)
#' }
#'
#' @export
install_gs4pb_deps <- function(
    cran_mirror    = "https://cloud.r-project.org",
    python_venv    = "gs4pb-python",
    install_python = TRUE
) {
  script <- system.file("setup", "install_dependencies.R", package = "GS4PB")
  if (!nzchar(script)) stop("Installer script not found in package.")

  # Source the script's internal function without auto-running it,
  # then call with user-supplied arguments.
  env <- new.env(parent = baseenv())
  sys.source(script, envir = env, keep.source = FALSE)
  env$.gs4pb_install(
    cran_mirror    = cran_mirror,
    python_venv    = python_venv,
    install_python = install_python
  )
}
