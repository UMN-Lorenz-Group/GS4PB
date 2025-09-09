#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options

run_app <- function(test_mode = FALSE,  
           onStart = NULL,
           options = list(),
           enableBookmarking = NULL,
           uiPattern = "/",
  ...
) {
	tryCatch({
	  set_global_options()
	  #initialize_environment()
          #setup_python()
	}, error = function(e) {
	  stop("Failed to initialize app environment: ", e$message)
	})

  with_golem_options(
    # app = shinyApp(
      # ui = app_ui,
      # server = app_server,
      # onStart = onStart,
      # options = options,
      # enableBookmarking = enableBookmarking,
      # uiPattern = uiPattern
     # ),

    app = shinyApp(
      ui = app_ui(),
      server = function(input, output, session) {
        if (test_mode) {
          options(shiny.testmode = TRUE)  # Enable test mode
        }
        app_server(input, output, session)
      },
	  onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
	),
   golem_opts = list(...)
  )
}
