# Export page: a downloadable analysis bundle for the selected fit.

mod_export_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      shiny::selectInput(ns("fit"), "Fit to export", character()),
      shiny::checkboxGroupInput(
        ns("parts"), "Include in the script",
        c("Rootogram data" = "rootogram",
          "Quantile residuals" = "quantile_residuals",
          "Predictions for the scenarios on the Predict page" = "scenarios"),
        selected = "rootogram"),
      shiny::checkboxInput(ns("include_data"),
                           "Include the uploaded data file", FALSE),
      help_text("Uploaded data are left out unless you include them; the ",
                "script then names the file it expects. Bundled examples ",
                "are loaded from flexCountReg."),
      shiny::downloadButton(ns("download"), "Download analysis bundle",
                            class = "btn-primary")),
    bslib::card(bslib::card_header("Bundle contents"),
                shiny::tableOutput(ns("manifest"))),
    bslib::card(bslib::card_header("analysis.R"),
                shiny::verbatimTextOutput(ns("script"))))
}

# Write an analysis bundle for a fit bundle into directory `dir`; returns
# the file names written.
write_bundle <- function(b, dir, options = list(), data = NULL,
                         include_data = FALSE) {
  spec <- b$spec
  files <- character(0)
  script <- analysis_script(spec, options, fit_label = b$label,
                            fit_id = b$id, versions = b$versions)
  writeLines(script, file.path(dir, "analysis.R"), useBytes = TRUE)
  spec_to_json(spec, file.path(dir, "settings.json"))
  files <- c(files, "analysis.R", "settings.json")
  if (include_data && spec$data$source != "bundled" && !is.null(data)) {
    utils::write.csv(data, file.path(dir, spec$data$name), row.names = FALSE)
    files <- c(files, spec$data$name)
  }
  f <- b$fit
  co <- as.data.frame(generics_tidy(f))
  utils::write.csv(co, file.path(dir, "coefficients.csv"), row.names = FALSE)
  files <- c(files, "coefficients.csv")
  m <- f$model
  info <- c(
    paste("Fit:", b$label), paste("Request:", b$id),
    paste("Created:", format(b$created)),
    paste("Data:", b$data_name, if (spec$data$source != "bundled")
      if (include_data) "(included)" else "(not included; see analysis.R)"),
    paste("Rows used:", m$nobs, "| rows removed for missing values:",
          m$n_dropped %||% 0),
    paste("Log-likelihood:", format(m$maximum, digits = 10)),
    paste("Converged:", isTRUE(m$converged)),
    if (!is.null(m$seed)) paste("Seed:", m$seed),
    paste("Versions:", paste(names(b$versions), b$versions, collapse = ", ")),
    if (length(b$warnings)) c("Warnings:", paste(" -", b$warnings)),
    "", "Session:", utils::capture.output(utils::sessionInfo()))
  writeLines(info, file.path(dir, "manifest.txt"), useBytes = TRUE)
  c(files, "manifest.txt")
}

generics_tidy <- function(f) flexCountReg::tidy(f, conf.int = TRUE)

mod_export_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observe({
      fits <- state$fits
      shiny::updateSelectInput(session, "fit",
                               choices = label_choices(fits),
                               selected = state$active)
    })
    bundle <- shiny::reactive({
      shiny::req(nzchar(input$fit %||% ""))
      state$fits[[input$fit]]
    })
    opts <- shiny::reactive({
      p <- input$parts %||% character()
      list(rootogram = "rootogram" %in% p,
           quantile_residuals = "quantile_residuals" %in% p,
           scenarios = if ("scenarios" %in% p) state$scenarios)
    })
    output$script <- shiny::renderText({
      b <- bundle()
      analysis_script(b$spec, opts(), fit_label = b$label, fit_id = b$id,
                      versions = b$versions)
    })
    output$manifest <- shiny::renderTable({
      b <- bundle()
      up <- b$spec$data$source != "bundled"
      data.frame(
        file = c("analysis.R", "settings.json", "coefficients.csv",
                 "manifest.txt", if (up && isTRUE(input$include_data))
                   b$spec$data$name),
        contents = c("R code that reproduces the fit and selected outputs",
                     "versioned app settings (can be reloaded)",
                     "estimates, standard errors and confidence limits",
                     "fit identity, sample, warnings, versions, session",
                     if (up && isTRUE(input$include_data)) "the uploaded data"))
    })
    output$download <- shiny::downloadHandler(
      filename = function() paste0("flexCountReg-analysis-", bundle()$id,
                                   ".zip"),
      content = function(file) {
        dir <- tempfile("bundle-")
        dir.create(dir)
        on.exit(unlink(dir, recursive = TRUE))
        files <- write_bundle(bundle(), dir, opts(), state$data,
                              isTRUE(input$include_data))
        zip::zip(file, files, root = dir)
      })
  })
}
