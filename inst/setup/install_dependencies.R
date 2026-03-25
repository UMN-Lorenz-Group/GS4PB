###############################################################################
##  GS4PB — Dependency Installer
##
##  Run this script ONCE after installing the GS4PB package to pull in
##  all required R packages, Bioconductor packages, GitHub packages, and
##  the optional Python environment (AlphaPlantImpute).
##
##  Usage (from an R console):
##    source(system.file("setup", "install_dependencies.R", package = "GS4PB"))
##  OR equivalently:
##    GS4PB::install_gs4pb_deps()
##
##  Prerequisites:
##    - R >= 4.1
##    - Java JDK >= 11 (required by rTASSEL / rJava)
##        Windows : https://adoptium.net/
##        Linux   : sudo apt install default-jdk && R CMD javareconf
##    - Python >= 3.9 (only needed for AlphaPlantImpute imputation)
##        Windows : https://www.python.org/downloads/
##        Linux   : sudo apt install python3 python3-venv python3-pip
###############################################################################

.gs4pb_install <- function(
    cran_mirror   = "https://cloud.r-project.org",
    python_venv   = "gs4pb-python",       # reticulate virtualenv name
    install_python = TRUE,                # set FALSE to skip Python setup
    ask            = interactive()        # set FALSE for non-interactive use
) {

  ## ── helpers ──────────────────────────────────────────────────────────────
  .installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)

  .try_install <- function(pkgs, installer, ...) {
    need <- pkgs[!vapply(pkgs, .installed, logical(1))]
    if (length(need) == 0) {
      message("  [ok] already installed: ", paste(pkgs, collapse = ", "))
      return(invisible(NULL))
    }
    message("  --> installing: ", paste(need, collapse = ", "))
    installer(need, ...)
  }

  .header <- function(txt) {
    message("\n", paste(rep("─", 70), collapse = ""), "\n  ", txt,
            "\n", paste(rep("─", 70), collapse = ""))
  }

  ## ── 0. Bootstrap installers ───────────────────────────────────────────────
  .header("Step 0 · Bootstrap remotes + BiocManager")

  if (!.installed("remotes")) {
    message("  --> installing remotes")
    utils::install.packages("remotes", repos = cran_mirror, quiet = TRUE)
  } else {
    message("  [ok] remotes")
  }

  if (!.installed("BiocManager")) {
    message("  --> installing BiocManager")
    utils::install.packages("BiocManager", repos = cran_mirror, quiet = TRUE)
  } else {
    message("  [ok] BiocManager")
  }

  ## ── 1. CRAN packages ──────────────────────────────────────────────────────
  .header("Step 1 · CRAN packages")

  cran_pkgs <- c(
    # Core pipeline
    "BGGE", "BGLR", "bWGR", "NAM", "rrBLUP", "STPGA", "ASRgenomics",
    # Data handling
    "dplyr", "plyr", "reshape2", "vcfR",
    # Parallel
    "callr", "doParallel", "foreach", "parallel", "future", "promises",
    # Geo / environment
    "nasapower", "geodata", "terra",
    # Shiny / UI
    "shiny", "shinyBS", "shinyFiles", "shinyjs", "shinyWidgets",
    "shinydashboard", "DT", "htmlwidgets", "plotly", "ggplot2", "gplots",
    # Python bridge
    "reticulate",
    # Misc
    "config", "fs", "golem", "grDevices", "graphics", "renv",
    "stats", "tools", "utils"
  )

  .try_install(cran_pkgs,
               utils::install.packages,
               repos = cran_mirror,
               quiet = TRUE)

  ## ── 2. Bioconductor packages ──────────────────────────────────────────────
  .header("Step 2 · Bioconductor packages")

  bioc_pkgs <- c("SummarizedExperiment")

  need_bioc <- bioc_pkgs[!vapply(bioc_pkgs, .installed, logical(1))]
  if (length(need_bioc) > 0) {
    message("  --> installing via BiocManager: ", paste(need_bioc, collapse = ", "))
    BiocManager::install(need_bioc, ask = FALSE, update = FALSE)
  } else {
    message("  [ok] ", paste(bioc_pkgs, collapse = ", "))
  }

  ## ── 3. GitHub packages ────────────────────────────────────────────────────
  .header("Step 3 · GitHub packages")

  # rTASSEL requires Java >= 11 and rJava to be working.
  # If rJava fails to load, run  R CMD javareconf  (Linux/Mac) or reinstall
  # the JDK and restart R (Windows) before retrying.

  github_pkgs <- list(
    rJava    = NULL,                                          # CRAN fallback
    rTASSEL  = "maize-genetics-genomics-lab/rTASSEL",
    EnvRtype = "allogamous/EnvRtype"
  )

  # rJava — try CRAN first; GitHub version rarely needed
  if (!.installed("rJava")) {
    message("  --> installing rJava from CRAN")
    utils::install.packages("rJava", repos = cran_mirror, quiet = TRUE)
    if (!.installed("rJava")) {
      message("  [!] rJava failed. Ensure a JDK >= 11 is installed and",
              "\n      run  R CMD javareconf  (Linux/Mac) or restart R (Windows).")
    }
  } else {
    message("  [ok] rJava")
  }

  for (pkg in names(github_pkgs)) {
    repo <- github_pkgs[[pkg]]
    if (is.null(repo)) next          # handled above
    if (!.installed(pkg)) {
      message("  --> remotes::install_github('", repo, "')")
      remotes::install_github(repo, upgrade = "never", quiet = TRUE)
    } else {
      message("  [ok] ", pkg, " (", repo, ")")
    }
  }

  ## ── 4. FactoMineR (Genomic Kinship PCA) ──────────────────────────────────
  .header("Step 4 · FactoMineR (Genomic Kinship PCA)")

  .try_install("FactoMineR",
               utils::install.packages,
               repos = cran_mirror,
               quiet = TRUE)

  ## ── 5. Python virtualenv (AlphaPlantImpute) ───────────────────────────────
  if (install_python) {
    .header("Step 5 · Python virtualenv for AlphaPlantImpute")

    if (!.installed("reticulate")) {
      message("  [!] reticulate not available — skipping Python setup.")
    } else {

      # Locate the bundled .whl file shipped with the package
      whl_path <- system.file("python",
                               "alphaplantimpute2-1.5.3-py3-none-any.whl",
                               package = "GS4PB")

      if (!nzchar(whl_path)) {
        message("  [!] alphaplantimpute2 wheel not found in package.",
                "\n      Place alphaplantimpute2-1.5.3-py3-none-any.whl in",
                "\n      inst/python/ and rebuild the package.")
      } else {

        venv_exists <- tryCatch(
          { reticulate::virtualenv_exists(python_venv); TRUE },
          error = function(e) FALSE
        )

        if (!venv_exists) {
          message("  --> creating virtualenv '", python_venv, "'")
          reticulate::virtualenv_create(python_venv)
        } else {
          message("  [ok] virtualenv '", python_venv, "' already exists")
        }

        # Core numeric stack (order matters: llvmlite before numba)
        pip_pkgs <- c("setuptools", "wheel", "numpy", "llvmlite", "numba")
        message("  --> pip install: ", paste(pip_pkgs, collapse = ", "))
        reticulate::virtualenv_install(python_venv, pip_pkgs,
                                       ignore_installed = FALSE,
                                       pip_options      = "--quiet")

        # AlphaPlantImpute from local wheel
        message("  --> pip install alphaplantimpute2 from bundled wheel")
        reticulate::virtualenv_install(python_venv, whl_path,
                                       ignore_installed = FALSE,
                                       pip_options      = "--quiet")

        message("  [ok] Python virtualenv '", python_venv, "' ready.")
        message("  To activate in future sessions:\n",
                "      reticulate::use_virtualenv('", python_venv, "', required=TRUE)")
      }
    }
  } else {
    message("\n  [skipped] Python setup (install_python = FALSE)")
  }

  ## ── 6. Java check ─────────────────────────────────────────────────────────
  .header("Step 6 · Java check (required for rTASSEL)")

  java_ok <- tryCatch({
    requireNamespace("rJava", quietly = TRUE) &&
      { rJava::.jinit(); TRUE }
  }, error = function(e) FALSE)

  if (java_ok) {
    jv <- tryCatch(rJava::.jcall("java/lang/System", "S", "getProperty",
                                  "java.version"), error = function(e) "unknown")
    message("  [ok] Java detected — version: ", jv)
  } else {
    message("  [!] Java / rJava not working.",
            "\n      Install JDK >= 11 from https://adoptium.net/",
            "\n      Linux: sudo apt install default-jdk && R CMD javareconf",
            "\n      Windows: reinstall JDK then restart R.")
  }

  ## ── Done ──────────────────────────────────────────────────────────────────
  .header("GS4PB dependency setup complete")
  message("  Start TASSEL in every session with:  rTASSEL::initializeTASSEL()")
  message("  Launch the app with:                 GS4PB::run_app()\n")
  invisible(NULL)
}

# Run immediately when sourced directly
.gs4pb_install()
