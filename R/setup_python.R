
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

  # Install missing Python packages in Conda environment.
  # envname must be passed: with no argument reticulate queries whatever python
  # is currently active, and at this point none is, so it fails with
  # "invalid 'path' argument" from path.expand() inside python_info(). That
  # aborted setup_python_env() after the environment was created but before any
  # package -- including the alphaplantimpute2 wheel -- was installed into it.
  installed_pkgs <- reticulate::py_list_packages(envname = envname)
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
    # conda_python(envname) rather than py_exe(): py_exe() reports whichever
    # interpreter reticulate has initialised, which is not necessarily -- and at
    # this point in a fresh session is not at all -- the environment just
    # created. The wheel has to land in `envname`, since that is what
    # use_condaenv("GS4PB_CondaEnv") loads at import time.
    py <- reticulate::conda_python(envname)

    # alphaplantimpute2.py imports pkg_resources, which setuptools stopped
    # shipping at version 81. A fresh conda env installs current setuptools, so
    # the wheel installs cleanly and then fails at import with
    # "ModuleNotFoundError: No module named 'pkg_resources'" -- reaching R as
    # the far less obvious "TypeError: 'ModuleNotFoundError' object does not
    # support item assignment". Pin setuptools so the module stays importable.
    # Remove this pin once alphaplantimpute2 drops pkg_resources.
    system(paste(shQuote(py), "-m pip install", shQuote("setuptools<81")))

    system(paste(shQuote(py), "-m pip install", shQuote(custom_package)))
  }
  
  
  new_path <- paste(Sys.getenv("PATH"), "~/.local/bin", sep = ":")
  Sys.setenv(PATH = new_path)

  message("Conda Python environment setup complete.")
}



