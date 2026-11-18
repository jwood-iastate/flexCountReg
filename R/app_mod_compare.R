# Compare page: information criteria of saved fits, with the conditions for
# comparing them stated next to the numbers.

mod_compare_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 300,
      shiny::checkboxGroupInput(ns("fits"), "Fits to compare", character()),
      help_text("Information criteria compare fits of the same response ",
                "on the same rows with the same weights.")),
    bslib::card(bslib::card_header("Comparison"),
                shiny::tableOutput(ns("table")),
                shiny::uiOutput(ns("notes"))))
}

# One row per fit, and the reasons the fits are not directly comparable
compare_fits <- function(bundles) {
  rows <- lapply(bundles, function(b) {
    f <- b$fit
    m <- f$model
    data.frame(fit = b$label, model = class(f)[1],
               family = paste(unique(m$family), collapse = "+"),
               classes = m$K %||% 1L, logLik = as.numeric(stats::logLik(f)),
               df = length(stats::coef(f)), AIC = stats::AIC(f),
               BIC = stats::BIC(f),
               ICL = m$ICL %||% NA_real_, rows = m$nobs,
               units = attr(stats::logLik(f), "nobs"),
               converged = isTRUE(m$converged), stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  why <- character(0)
  key <- function(b) {
    s <- b$spec
    list(data = b$data_name, response = s$recipe$response,
         weights = s$recipe$weights %||% "", rows = b$fit$model$rows_used,
         units = attr(stats::logLik(b$fit), "nobs"))
  }
  ks <- lapply(bundles, key)
  if (length(ks) > 1) {
    same <- function(nm) all(vapply(ks[-1], function(k)
      identical(k[[nm]], ks[[1]][[nm]]), logical(1)))
    if (!same("data")) why <- c(why, "The fits use different datasets.")
    if (!same("response")) why <- c(why, "The responses differ.")
    if (!same("weights")) why <- c(why, "The frequency weights differ.")
    if (!same("rows")) why <- c(why, paste(
      "The fits use different rows (missing values in different",
      "columns); refit them on a common sample before ranking."))
    if (!same("units")) why <- c(why, paste(
      "The number of independent units differs (rows versus panels), so",
      "BIC penalties are not on the same scale."))
  }
  list(table = tab, not_comparable = why)
}

mod_compare_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observe({
      fits <- state$fits
      shiny::updateCheckboxGroupInput(
        session, "fits",
        choices = label_choices(fits),
        selected = shiny::isolate(input$fits))
    })
    cmp <- shiny::reactive({
      shiny::req(length(input$fits) >= 1)
      compare_fits(state$fits[input$fits])
    })
    output$table <- shiny::renderTable(cmp()$table, digits = 3)
    output$notes <- shiny::renderUI({
      cp <- cmp()
      mix <- any(cp$table$model %in% c("countreg_mix", "countreg_lc"))
      shiny::tagList(
        if (length(cp$not_comparable)) shiny::div(
          class = "alert alert-warning", role = "alert",
          shiny::strong("Not directly comparable:"),
          shiny::tags$ul(lapply(cp$not_comparable, shiny::tags$li)))
        else if (nrow(cp$table) > 1) shiny::p(
          "Same data, response, rows and weights: AIC and BIC (smaller is ",
          "better) can be compared."),
        if (mix) shiny::p(
          "No likelihood-ratio p-value is shown for numbers of classes: ",
          "its usual chi-square distribution does not hold. Use AIC/BIC, ",
          "start stability and held-out predictive scores. ICL (smaller ",
          "is better) also rewards well-separated classes."))
    })
  })
}
