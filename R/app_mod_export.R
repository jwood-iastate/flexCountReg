# Export page: download the session script (everything the app ran) with
# the settings of each fit, a manifest and, optionally, the uploaded data.

mod_export_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      shiny::checkboxInput(ns("include_data"),
                           "Include uploaded data files", FALSE),
      help_text("Uploaded data are left out unless you include them; the ",
                "script then names the files it expects. Example datasets ",
                "are loaded from flexCountReg."),
      shiny::downloadButton(ns("download"), "Download the analysis",
                            class = "btn-primary")),
    bslib::card(bslib::card_header("Contents"), fill = FALSE,
                shiny::tableOutput(ns("manifest"))),
    bslib::card(bslib::card_header("analysis.R (the session script)"),
                fill = FALSE,
                shiny::verbatimTextOutput(ns("script"))))
}

# Write the analysis bundle of a session into `dir`; returns the file names
# written. `state` is read without reactivity (call inside isolate()).
write_bundle <- function(state, dir, include_data = FALSE) {
  files <- character(0)
  writeLines(app_script(state), file.path(dir, "analysis.R"),
             useBytes = TRUE)
  files <- c(files, "analysis.R")
  fits <- state$fits
  if (length(fits)) {
    specs <- lapply(fits, function(b) unclass(b$spec))
    names(specs) <- vapply(fits, `[[`, "", "fit_name")
    writeLines(jsonlite::toJSON(specs, auto_unbox = TRUE, null = "null",
                                pretty = TRUE, digits = NA),
               file.path(dir, "settings.json"), useBytes = TRUE)
    files <- c(files, "settings.json")
  }
  ups <- state$uploads
  if (include_data) {
    for (u in ups) {
      file.copy(u$path, file.path(dir, u$name), overwrite = TRUE)
      files <- c(files, u$name)
    }
  }
  info <- c(
    paste("Created:", format(Sys.time())),
    paste("Versions:", paste(names(fcr_versions()), fcr_versions(),
                             collapse = ", ")),
    if (length(ups)) paste("Uploaded data:",
                           toString(vapply(ups, `[[`, "", "name")),
                           if (include_data) "(included)" else
                             "(not included; analysis.R names the files)"),
    "", "Fits:",
    unlist(lapply(fits, function(b) {
      f <- get(b$fit_name, envir = state$env)
      m <- f$model
      c(paste0("  ", b$fit_name, ": ", b$label),
        paste0("    rows used: ", m$nobs %||% NA,
               " | log-likelihood: ", format(m$maximum, digits = 10),
               " | converged: ", isTRUE(m$converged)),
        if (length(b$warnings)) paste0("    warning: ", b$warnings))
    })),
    "", "Session:", utils::capture.output(utils::sessionInfo()))
  writeLines(info, file.path(dir, "manifest.txt"), useBytes = TRUE)
  c(files, "manifest.txt")
}

mod_export_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    output$script <- shiny::renderText({
      if (!length(state$log)) return("Nothing has been run yet.")
      app_script(state)
    })
    output$manifest <- shiny::renderTable({
      ups <- state$uploads
      data.frame(
        file = c("analysis.R", if (length(state$fits)) "settings.json",
                 "manifest.txt",
                 if (isTRUE(input$include_data))
                   vapply(ups, `[[`, "", "name")),
        contents = c("R code of everything the app ran, in order",
                     if (length(state$fits))
                       "settings of each fit (JSON)",
                     "fits, sample sizes, warnings, versions, session",
                     if (isTRUE(input$include_data))
                       rep("uploaded data", length(ups))))
    })
    output$download <- shiny::downloadHandler(
      filename = function() paste0("flexCountReg-analysis-",
                                   format(Sys.time(), "%Y%m%d-%H%M%S"),
                                   ".zip"),
      content = function(file) {
        dir <- tempfile("bundle-")
        dir.create(dir)
        on.exit(unlink(dir, recursive = TRUE))
        files <- shiny::isolate(write_bundle(state, dir,
                                             isTRUE(input$include_data)))
        zip::zip(file, files, root = dir)
      })
  })
}
