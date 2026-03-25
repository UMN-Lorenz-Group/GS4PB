### .onLoad function 

.onLoad <- function(libname, pkgname) {

  envname <- "GS4PB_CondaEnv"
  # Use Conda — non-fatal if Conda or the env is absent (e.g. during R CMD check)
  tryCatch({
    if (envname %in% reticulate::conda_list()$name) {
      reticulate::use_condaenv(envname, required = FALSE)
    } else {
      packageStartupMessage(
        "Conda environment '", envname, "' not found. ",
        "Run GS4PB:::setup_python_env() to create it."
      )
    }
  }, error = function(e) {
    packageStartupMessage(
      "Conda not available or errored on load. ",
      "Run GS4PB:::setup_python_env() after installing Miniconda."
    )
  })

}



utils::globalVariables(c("env.id", "lat", "lon", "Loc")) 

utils::globalVariables(c("POS", "REF", "ALT", "ID","CHROM")) 

utils::globalVariables(c("IDColsMEList","tstIndices2","DT_tstIndices2"))

utils::globalVariables("FN")

utils::globalVariables(c("ID","LocationME","X","YearME","fixME","nFact","strainGeno","varEnv"))
