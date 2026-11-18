# Compare page: information criteria of saved fits (computed by logged code)
# with the conditions for comparing them stated next to the numbers.

mod_compare_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      shiny::checkboxGroupInput(ns("fits"), "Fits to compare", character()),
      help_text("Information criteria compare fits of the same response ",
                "on the same rows with the same weights."),
      shiny::actionButton(ns("go"), "Compare fits", class = "btn-primary")),
    bslib::card(bslib::card_header(shiny::textOutput(ns("title"))),
                fill = FALSE,
                shiny::tableOutput(ns("table")),
                shiny::uiOutput(ns("notes"))),
    bslib::card(bslib::card_header("Code that produced it"), fill = FALSE,
                shiny::verbatimTextOutput(ns("code"))))
}

# Reasons the fits (bundles with their fitted objects) are not directly
# comparable
comparability <- function(bundles) {
  why <- character(0)
  if (length(bundles) < 2) return(why)
  key <- function(b) {
    s <- b$spec
    list(data = s$data$object, response = s$recipe$response,
         weights = s$recipe$weights %||% "",
         rows = b$fit$model$rows_used,
         units = attr(stats::logLik(b$fit), "nobs"))
  }
  ks <- lapply(bundles, key)
  same <- function(nm) all(vapply(ks[-1], function(k)
    identical(k[[nm]], ks[[1]][[nm]]), logical(1)))
  if (!same("data")) why <- c(why, "The fits use different datasets.")
  if (!same("response")) why <- c(why, "The responses differ.")
  if (!same("weights")) why <- c(why, "The frequency weights differ.")
  if (!same("rows")) why <- c(why, paste(
    "The fits use different rows (missing values in different columns);",
    "refit them on a common sample before ranking."))
  if (!same("units")) why <- c(why, paste(
    "The number of independent units differs (rows versus panels), so the",
    "BIC penalties are not on the same scale."))
  why
}

mod_compare_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    current <- shiny::reactiveVal(NULL)
    shiny::observe({
      fits <- state$fits
      shiny::updateCheckboxGroupInput(
        session, "fits", choices = label_choices(fits),
        selected = shiny::isolate(input$fits))
    })
    shiny::observeEvent(input$go, app_try(what = "The comparison", {
      ids <- input$fits
      if (length(ids) < 1) stop("Choose at least one fit.", call. = FALSE)
      bl <- lapply(ids, function(id) {
        b <- state$fits[[id]]
        b$fit <- app_fit_object(state, id)
        b
      })
      req <- list(type = "compare",
                  fits = vapply(bl, `[[`, "", "fit_name"),
                  labels = vapply(bl, `[[`, "", "label"))
      res <- app_run_output(state, req)
      res$why <- comparability(bl)
      res$mixture <- any(vapply(bl, function(b)
        inherits(b$fit, "countreg_mix"), logical(1)))
      current(res)
    }))
    output$title <- shiny::renderText({
      r <- current()
      if (is.null(r)) "No comparison yet" else paste("Comparison", r$out)
    })
    output$table <- shiny::renderTable({
      r <- current()
      shiny::req(r)
      r$value
    }, digits = 3)
    output$code <- shiny::renderText({
      r <- current()
      shiny::req(r)
      paste(vapply(r$exprs, .deparse, character(1)), collapse = "\n")
    })
    output$notes <- shiny::renderUI({
      r <- current()
      shiny::req(r)
      shiny::tagList(
        if (length(r$why)) shiny::div(
          class = "alert alert-warning", role = "alert",
          shiny::strong("Not directly comparable:"),
          shiny::tags$ul(lapply(r$why, shiny::tags$li)))
        else if (nrow(r$value) > 1) shiny::p(
          "Same data, response, rows and weights: AIC and BIC (smaller is ",
          "better) can be compared."),
        if (isTRUE(r$mixture)) shiny::p(
          "No likelihood-ratio p-value is shown for numbers of classes: ",
          "its usual chi-square distribution does not hold. Use AIC/BIC, ",
          "start stability and held-out predictive scores."))
    })
  })
}
