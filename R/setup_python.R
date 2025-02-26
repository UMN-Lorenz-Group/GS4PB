
setup_python_env <- function(envname = "GS4PB_CondaEnv",python_version = "3.12", packages = c("numba","libstdcxx-ng","numpy", "pandas")){
  library(reticulate)
   

  # Check if Conda is installed
  if (is.null(reticulate::conda_binary())) {
    stop("Conda not found. Install Anaconda or Miniconda first.")
  }

  # Create a Conda environment if it does not exist
  if (!envname %in% reticulate::conda_list()$name) {
    message("Creating Conda environment: ", envname)
    reticulate::conda_create(envname,python_version= "3.12")
  }

  # Install missing Python packages in Conda environment
  installed_pkgs <- reticulate::py_list_packages()
  missing_pkgs <- setdiff(packages, installed_pkgs$package)

  if (length(missing_pkgs) > 0) {
    message("Installing missing Python packages: ", paste(missing_pkgs, collapse = ", "))
    reticulate::conda_install(envname, missing_pkgs)
  }
  
   # Install custom package from .whl or .tar.gz
  python_path <- system.file("python", package = "GS4PB")
  custom_package <- list.files(python_path, pattern = "\\.whl|\\.tar\\.gz$", full.names = TRUE)
  
  if (length(custom_package) > 0) {
    message("Installing custom Python package: ", custom_package)
    system(paste(shQuote(reticulate::py_exe()), "-m pip install", shQuote(custom_package)))
  }
  
  
  new_path <- paste(Sys.getenv("PATH"), "~/.local/bin", sep = ":")
  Sys.setenv(PATH = new_path)

  message("Conda Python environment setup complete.")
}



