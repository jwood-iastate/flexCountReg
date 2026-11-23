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
# Uploaded files are written to their own paths in the bundle
# (data/data_1.csv, data/newdata_1.csv, ...; see .app_upload_file()), the
# paths the session script reads, so files uploaded under the same name
# stay distinct. Each copy is checked against the MD5 checksum recorded
# when the file was loaded.
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
      dest <- file.path(dir, u$file)
      dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
      ok <- file.exists(u$path) && file.copy(u$path, dest)
      if (!isTRUE(ok)) {
        stop("Could not copy the uploaded file \"", u$name, "\" to ",
             u$file, ".", call. = FALSE)
      }
      if (!identical(unname(tools::md5sum(dest)), u$md5)) {
        stop("The copy of the uploaded file \"", u$name, "\" (", u$file,
             ") does not match the file that was loaded.", call. = FALSE)
      }
      files <- c(files, u$file)
    }
  }
  info <- c(
    paste("Created:", format(Sys.time())),
    paste("Versions:", paste(names(fcr_versions()), fcr_versions(),
                             collapse = ", ")),
    if (length(ups)) c(
      paste("Uploaded data", if (include_data) "(included):" else
        "(not included; put each file at the path analysis.R reads):"),
      unlist(lapply(names(ups), function(obj) {
        u <- ups[[obj]]
        paste0("  ", u$file, " <- uploaded as \"", u$name, "\" (object ",
               obj, "; MD5 ", u$md5, ")")
      }))),
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
    output_manifest(state$log),
    "", "Session:", utils::capture.output(utils::sessionInfo()))
  writeLines(info, file.path(dir, "manifest.txt"), useBytes = TRUE)
  c(files, "manifest.txt")
}

# Manifest lines of the outputs in a session log (predictions, intervals,
# diagnostics, comparisons): each object with its title, the details of how
# it was obtained (uncertainty method, bootstrap refits, seed) and the
# warnings raised while it ran
output_manifest <- function(log) {
  outs <- Filter(function(e) !e$kind %in% c("data", "transform", "fit",
                                            "plot_file"), log)
  if (!length(outs)) return(character(0))
  c("", "Outputs:",
    unlist(lapply(outs, function(e) {
      c(paste0("  ", e$ref, ": ", gsub("[[:space:]]+", " ", e$title)),
        if (length(e$details)) paste0("    ", e$details),
        if (length(e$warnings)) paste0("    warning: ", e$warnings))
    })))
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
                   vapply(ups, `[[`, "", "file")),
        contents = c("R code of everything the app ran, in order",
                     if (length(state$fits))
                       "settings of each fit (JSON)",
                     paste("fits, sample sizes, outputs, warnings, seeds,",
                           "checksums, versions, session"),
                     if (isTRUE(input$include_data))
                       paste0("uploaded data (", vapply(ups, `[[`, "",
                                                         "name"), ")")))
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
