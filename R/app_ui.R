# Packages the app needs; they are suggested, not imported, so the
# statistical functions of flexCountReg never require Shiny.
.app_packages <- c("shiny", "bslib", "DT", "callr", "jsonlite", "zip")

.app_check_packages <- function() {
  rlang::check_installed(.app_packages,
                         reason = "to use the flexCountReg app.")
}

#' Interactive analysis app
#'
#' `launch_app()` opens an interactive interface to flexCountReg in a web
#' browser. In it you can load data (a bundled example or a CSV file),
#' specify and estimate count regression, finite-mixture and latent-class
#' models, check their numerical status, diagnose them (rootograms,
#' randomized quantile residuals, CURE plots), compare fits, compute
#' predictions for scenarios, and download the analysis.
#'
#' Every fit is defined by one specification that is compiled into R code:
#' the code that runs and the script you download are the same, so the
#' downloaded `analysis.R` reproduces the fit without the app (the
#' "R code" panel at the bottom shows the draft and the fitted code).
#' Models are estimated in background R processes, so the interface stays
#' responsive and long fits can be cancelled. Data, fits and jobs belong to
#' one browser session and are removed when it ends; uploaded data are read
#' as CSV only and stay on your computer.
#'
#' The app needs the packages shiny, bslib, DT, callr, jsonlite and zip
#' (suggested, not required, by flexCountReg); `launch_app()` offers to
#' install any that are missing.
#'
#' @param ... arguments passed to [shiny::runApp()], such as `port` or
#'   `launch.browser`.
#' @param max_upload_mb largest CSV file that can be uploaded, in
#'   megabytes.
#' @returns `launch_app()` runs the app (it returns when the app is
#'   closed). `app_object()` returns the app as a `shiny.appobj`, for
#'   deployment or testing.
#' @examples
#' if (interactive()) {
#'   launch_app()
#' }
#' @export
launch_app <- function(..., max_upload_mb = 50) {
  .app_check_packages()
  shiny::runApp(app_object(max_upload_mb = max_upload_mb), ...)
}

#' @rdname launch_app
#' @export
app_object <- function(max_upload_mb = 50) {
  .app_check_packages()
  shiny::shinyApp(
    ui = app_ui(),
    server = app_server,
    onStart = function() {
      old <- options(shiny.maxRequestSize = max_upload_mb * 1024^2)
      shiny::onStop(function() options(old))
    })
}

app_ui <- function() {
  theme <- bslib::bs_theme(version = 5, preset = "shiny")
  bslib::page_navbar(
    title = "flexCountReg",
    id = "page",
    theme = theme,
    fillable = FALSE,
    header = shiny::div(
      class = "container-fluid py-1 small border-bottom bg-light",
      role = "status", `aria-live` = "polite",
      shiny::textOutput("status_bar", inline = TRUE)),
    bslib::nav_panel("Data", mod_data_ui("data")),
    bslib::nav_panel("Build", mod_build_ui("build")),
    bslib::nav_panel("Fit", mod_fit_ui("fit")),
    bslib::nav_panel("Diagnose", mod_diagnose_ui("diagnose")),
    bslib::nav_panel("Compare", mod_compare_ui("compare")),
    bslib::nav_panel("Predict", mod_predict_ui("predict")),
    bslib::nav_panel("Export", mod_export_ui("export")),
    footer = bslib::accordion(
      open = FALSE, class = "mt-3",
      bslib::accordion_panel(
        "R code", icon = NULL,
        bslib::navset_tab(
          bslib::nav_panel("Draft code",
                           shiny::p(class = "small text-muted mt-2",
                                    "The analysis defined by the current ",
                                    "controls (not yet estimated)."),
                           shiny::verbatimTextOutput("draft_code")),
          bslib::nav_panel("Fitted code",
                           shiny::p(class = "small text-muted mt-2",
                                    "The exact code of the selected ",
                                    "completed fit."),
                           shiny::verbatimTextOutput("fitted_code"))))))
}

# Choices for a list of fits or jobs: ids named by their labels
label_choices <- function(x) {
  if (!length(x)) return(character())
  stats::setNames(names(x), vapply(x, function(e) e$label, character(1)))
}

# A small helper for help text under a control
help_text <- function(...) shiny::p(class = "form-text small", ...)
