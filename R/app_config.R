#' Access files in the current app
#'
#' NOTE: If you manually change your package name in the DESCRIPTION,
#' don't forget to change it here too, and in the config file.
#' For a safer name change mechanism, use the `golem::set_golem_name()` function.
#'
#' @param ... character vectors, specifying subdirectory and file(s)
#' within your package. The default, none, returns the root of the app.
#'
#' @noRd

app_sys <- function(...) {
  # Ensure the resource path is set
  golem::add_resource_path("www", system.file("app/www", package = "GS4PB"))
  
  # Retrieve the file path within the package
  # If system.file() doesn't find the file, construct the path relative to inst/
  file_path <- file.path(system.file(package = "GS4PB"),...)
  
  return(file_path)
}



#' Set Global Options
#'
#' Configures necessary global options for the application to function properly.
#' @importFrom config get
#' @export
#'
set_global_options <- function() {
  config_file <- app_sys("golem-config.yml")

  if (!file.exists(config_file)) {
    stop("Config file not found at: ", config_file)
  }

  config <- get(file = config_file)

  options(java.parameters = config$java.parameters)
  options(error = eval(parse(text = config$error)))
  options(shiny.maxRequestSize = eval(parse(text=config$maxRequestSize)))

}

#' Initialize App Environment
#'
#' Configures the Python virtual environment and Shiny options.
#' @importFrom config get
#' @importFrom reticulate use_virtualenv
#' @export
initialize_environment <- function() {
  config_file <- app_sys("golem-config.yml")

  if (!file.exists(config_file)) {
    stop("Config file not found at: ", config_file)
  }

  config <- get(file = config_file)

   #  use_virtualenv(
   #   config$virtualenv,
   #  required = TRUE
   #)
}


#' Read App Config
#'
#' @param value Value to retrieve from the config file.
#' @param config GOLEM_CONFIG_ACTIVE value. If unset, R_CONFIG_ACTIVE.
#' If unset, "default".
#' @param use_parent Logical, scan the parent directory for config file.
#' @param file Location of the config file
#' @importFrom config get
#'
#' @noRd
get_golem_config <- function(
    value,
    config = Sys.getenv(
      "GOLEM_CONFIG_ACTIVE",
      Sys.getenv(
        "R_CONFIG_ACTIVE",
        "default"
      )
    ),
    use_parent = TRUE,
    # Modify this if your config file is somewhere else
    file = app_sys("golem-config.yml")
) {
  get(
    value = value,
    config = config,
    file = file,
    use_parent = use_parent
  )
}
