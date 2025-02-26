### .onLoad function 

.onLoad <- function(libname, pkgname) {
  library(reticulate)

  envname <- "GS4PB_CondaEnv"
  ## Activate Python virtual environment on package load
  #if (virtualenv_exists(envname)) {
   # use_virtualenv(envname, required = TRUE)
  #} else {
   # warning("Python virtual environment not found. Run setup_python_env() to create it.")
  #}

  # Use Conda instead of Virtualenv
  if (envname %in% reticulate::conda_list()$name) {
    reticulate::use_condaenv(envname, required = TRUE)
  } else {
    warning("Conda environment not found. Run GS4PB:::setup_python_env() to create it.")
  }


}



utils::globalVariables(c("env.id", "lat", "lon", "Loc")) 

utils::globalVariables(c("POS", "REF", "ALT", "ID","CHROM")) 

utils::globalVariables(c("IDColsMEList","tstIndices2","DT_tstIndices2"))

utils::globalVariables("FN")

utils::globalVariables(c("ID","LocationME","X","YearME","fixME","nFact","strainGeno","varEnv"))
