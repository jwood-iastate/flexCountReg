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
#' browser, with a tab for each step: load data (a bundled example or a CSV
#' file) and create new variables (logs, powers, standardized values,
#' indicators, products, ratios, categories, factors); specify and estimate
#' count regression ([countreg()]), random-parameter ([countreg_rp()]),
#' hierarchical ([countreg_hier()]: random intercepts and slopes, crossed
#' or nested levels, normal or non-normal effects), fixed-effects
#' ([countreg_fe()]), random-effects panel ([renb()], [poisLind.re()],
#' [repois()], [negmultinom()]), finite-mixture and latent-class
#' ([countreg_mix()], [countreg_lc()]) models, with smooth (nonlinear)
#' effects `s(x)`, correlated random effects (`cre`), cluster-robust
#' standard errors and sampling weights where the model allows them; check
#' their numerical status; diagnose them (rootograms, randomized quantile
#' residuals, CURE plots, smooth effects); compare fits;
#' compute predictions, confidence or prediction intervals and tolerance
#' intervals for a table of scenarios or for uploaded new data (with the
#' variables created on the Data tab, using the means, standard deviations
#' and categories of the training data); save plots as JPEG files at a chosen
#' size and resolution; and download the analysis. The Diagnose and Predict
#' tabs work on any completed model.
#'
#' Everything the app runs is compiled into R code that calls the package's
#' functions, evaluated from that code, and appended to a session script:
#' data loading, new variables, fits, diagnostics, predictions and
#' comparisons, in order. The "R code" panel at the bottom shows the
#' session script (and the code of the model being built), and the Export
#' tab downloads it as `analysis.R`, which reproduces the whole session
#' without the app. Column names and typed values enter the code only as
#' names, numbers or strings, never as code.
#' Models are estimated in background R processes, so the interface stays
#' responsive and long fits can be canceled. Data, fits and jobs belong to
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
  bslib::page_fluid(
    title = "flexCountReg",
    theme = theme,
    shiny::div(
      class = "d-flex flex-wrap align-items-baseline gap-3 pt-2",
      shiny::h4("flexCountReg", class = "m-0"),
      shiny::div(class = "small text-muted", role = "status",
                 `aria-live` = "polite",
                 shiny::textOutput("status_bar", inline = TRUE))),
    # The pages are tabs across the top (they scroll sideways on a narrow
    # screen instead of collapsing into a menu)
    shiny::tags$style(shiny::HTML(
      "#page.nav-tabs { flex-wrap: nowrap; overflow-x: auto;
                        overflow-y: hidden; }
       #page.nav-tabs .nav-link { white-space: nowrap; }")),
    bslib::navset_tab(
      id = "page",
      bslib::nav_panel("Data", mod_data_ui("data")),
      bslib::nav_panel("Build", mod_build_ui("build")),
      bslib::nav_panel("Fit", mod_fit_ui("fit")),
      bslib::nav_panel("Diagnose", mod_diagnose_ui("diagnose")),
      bslib::nav_panel("Compare", mod_compare_ui("compare")),
      bslib::nav_panel("Predict", mod_predict_ui("predict")),
      bslib::nav_panel("Export", mod_export_ui("export"))),
    bslib::accordion(
      open = FALSE, class = "mt-3 mb-3",
      bslib::accordion_panel(
        "R code", icon = NULL,
        bslib::navset_tab(
          bslib::nav_panel(
            "Session script",
            shiny::p(class = "small text-muted mt-2",
                     "Everything the app has run, in order: data, fits, ",
                     "diagnostics, predictions and comparisons. Export ",
                     "downloads this script."),
            shiny::verbatimTextOutput("session_code")),
          bslib::nav_panel(
            "Draft model",
            shiny::p(class = "small text-muted mt-2",
                     "The code the Build page would add for the next fit ",
                     "(not yet run)."),
            shiny::verbatimTextOutput("draft_code"))))))
}

# Choices for a list of fits or jobs: ids named by their labels
label_choices <- function(x) {
  if (!length(x)) return(character())
  stats::setNames(names(x), vapply(x, function(e) e$label, character(1)))
}

# A small helper for help text under a control
help_text <- function(...) shiny::p(class = "form-text small", ...)

# What an output result `r` reports beside it: how it was obtained
# (r$details: uncertainty method, bootstrap refits accepted, seed), the
# warnings raised while it ran (r$warnings) and its messages (r$messages).
# Warnings are shown as on the Fit page's Messages tab, not only in the
# console.
output_notes <- function(r) {
  if (is.null(r)) return(NULL)
  li <- function(x) shiny::tags$ul(lapply(x, shiny::tags$li))
  shiny::tagList(
    if (length(r$details)) {
      shiny::p(class = "small text-muted",
               paste(r$details, collapse = "; "), ".")
    },
    if (length(r$warnings)) {
      shiny::div(class = "alert alert-warning small", role = "alert",
                 shiny::strong("Warnings"), li(r$warnings))
    },
    if (length(r$messages)) {
      shiny::div(class = "small text-muted", shiny::strong("Messages"),
                 li(r$messages))
    })
}
