# Predict page: an editable table of scenarios; predictions, confidence or
# prediction intervals, and tolerance intervals of the selected fit. Each
# calculation runs its compiled code, which is added to the session script.

mod_predict_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340,
      shiny::uiOutput(ns("fit_note")),
      bslib::accordion(
        multiple = FALSE,
        bslib::accordion_panel(
          "Predictions",
          shiny::uiOutput(ns("method_ui")),
          shiny::selectInput(ns("interval"), "Interval",
                             c("None" = "none",
                               "Confidence interval for the expected count" =
                                 "confidence",
                               "Prediction interval for a new count" =
                                 "prediction")),
          shiny::numericInput(ns("level"), "Level", 0.95, min = 0.5,
                              max = 0.999, step = 0.01),
          help_text("A confidence interval concerns the expected count; a ",
                    "prediction interval concerns one future count."),
          shiny::actionButton(ns("go"), "Calculate predictions",
                              class = "btn-primary")),
        bslib::accordion_panel(
          "Tolerance intervals",
          help_text("Limits that contain at least a proportion (content) ",
                    "of the count distribution of each scenario, with a ",
                    "stated confidence."),
          shiny::numericInput(ns("content"), "Content", 0.90, min = 0.5,
                              max = 0.999, step = 0.01),
          shiny::numericInput(ns("confidence"), "Confidence", 0.95,
                              min = 0.5, max = 0.999, step = 0.01),
          shiny::selectInput(ns("side"), "Side",
                             c("Two-sided" = "two-sided",
                               "Upper limit" = "upper",
                               "Lower limit" = "lower")),
          shiny::selectInput(ns("scope"), "Scope",
                             c("Pointwise (each scenario)" = "pointwise",
                               "Simultaneous (all scenarios)" =
                                 "simultaneous")),
          shiny::selectInput(ns("tol_method"), "Method",
                             c("auto", "joint", "envelope", "exact")),
          help_text("See ?flexCountReg::tolerance_interval: the joint and ",
                    "envelope methods evaluate a finite set of parameter ",
                    "points; exact applies to intercept-only Poisson ",
                    "models."),
          shiny::actionButton(ns("go_tol"), "Calculate tolerance intervals",
                              class = "btn-primary")))),
    bslib::card(
      bslib::card_header("Scenarios"), fill = FALSE,
      shiny::div(
        class = "d-flex gap-2 mb-2",
        shiny::actionButton(ns("reset"), "Start from typical values",
                            class = "btn-outline-secondary btn-sm"),
        shiny::actionButton(ns("add"), "Add a row",
                            class = "btn-outline-secondary btn-sm")),
      help_text("Double-click a cell to edit it. Numeric columns start at ",
                "their mean, categorical ones at their first level."),
      DT::DTOutput(ns("scen"))),
    bslib::card(bslib::card_header(shiny::textOutput(ns("res_title"))),
                fill = FALSE,
                shiny::uiOutput(ns("res_note")),
                shiny::tableOutput(ns("res"))),
    bslib::card(bslib::card_header("Code that produced it"), fill = FALSE,
                shiny::verbatimTextOutput(ns("code"))))
}

# Columns a fit needs in new data
needed_columns <- function(spec) {
  r <- spec$recipe
  m <- spec$model
  unique(c(r$predictors, r$exposure, r$log_offset, m$dis1, m$mixing,
           m$rpar))
}

# One row of typical values (mean of numeric columns, first level or most
# common value otherwise), after the recipe's factor conversions
typical_row <- function(spec, data) {
  for (f in spec$recipe$factors) data[[f]] <- factor(data[[f]])
  cols <- needed_columns(spec)
  if (!length(cols)) return(data.frame(row.names = 1L))
  out <- lapply(cols, function(cn) {
    v <- data[[cn]]
    if (is.factor(v)) {
      factor(levels(v)[1], levels = levels(v))
    } else if (is.numeric(v)) {
      signif(mean(v, na.rm = TRUE), 4)
    } else {
      tab <- sort(table(v), decreasing = TRUE)
      if (length(tab)) names(tab)[1] else NA
    }
  })
  stats::setNames(as.data.frame(out, check.names = FALSE,
                                stringsAsFactors = FALSE), cols)
}

mod_predict_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    scen <- shiny::reactiveVal(NULL)
    current <- shiny::reactiveVal(NULL)
    bundle <- shiny::reactive({
      b <- state$fits[[state$active %||% ""]]
      shiny::req(b)
      b$fit <- app_fit_object(state, b$id)
      b
    })
    output$fit_note <- shiny::renderUI({
      b <- state$fits[[state$active %||% ""]]
      help_text(if (is.null(b)) "Estimate a model first (Build page)." else
        paste0("Predictions of the selected fit: ", b$label, " (",
               b$fit_name, ")."))
    })
    output$method_ui <- shiny::renderUI({
      m <- app_predict_methods(bundle()$fit)
      if (is.null(m)) return(NULL)
      shiny::selectInput(session$ns("method"), "Prediction method", m)
    })
    shiny::observeEvent(list(bundle(), input$reset),
                        app_try(what = "Building the scenario table", {
      b <- bundle()
      data <- get(b$spec$data$object, envir = state$env)
      scen(typical_row(b$spec, data))
    }))
    shiny::observeEvent(input$add, {
      s <- scen()
      if (!is.null(s) && nrow(s)) scen(rbind(s, s[nrow(s), , drop = FALSE]))
    })
    output$scen <- DT::renderDT({
      s <- scen()
      shiny::req(s)
      DT::datatable(s, rownames = FALSE, editable = "cell",
                    fillContainer = FALSE,
                    options = list(dom = "t", pageLength = 50,
                                   scrollX = TRUE))
    })
    shiny::observeEvent(input$scen_cell_edit, app_try(what = "Editing", {
      s <- scen()
      e <- input$scen_cell_edit
      j <- e$col + 1
      v <- s[[j]]
      if (is.factor(v) && !e$value %in% levels(v)) {
        stop("Unknown level \"", e$value, "\"; levels: ",
             toString(levels(v)), call. = FALSE)
      }
      s[e$row, j] <- if (is.numeric(v)) as.numeric(e$value) else e$value
      scen(s)
    }))
    shiny::observe(state$scenarios <- scen())

    run <- function(req) {
      b <- bundle()
      s <- scen()
      if (is.null(s) || !nrow(s)) {
        stop("The scenario table is empty.", call. = FALSE)
      }
      req$scenarios <- s
      res <- app_run_output(state, req, b$id)
      res$type <- req$type
      res$fit$fit <- b$fit
      current(res)
    }
    shiny::observeEvent(input$go, app_try(what = "The predictions", {
      f <- bundle()$fit
      interval <- input$interval
      if (interval != "none") {
        sp <- app_supports(f, "interval")
        if (!isTRUE(sp$ok)) stop(sp$reason, call. = FALSE)
      }
      run(list(type = "predict", interval = interval, level = input$level,
               method = if (!is.null(app_predict_methods(f))) input$method,
               mixture = inherits(f, "countreg_mix")))
    }))
    shiny::observeEvent(input$go_tol, app_try(what = "The tolerance intervals", {
      f <- bundle()$fit
      sp <- app_supports(f, "tolerance")
      if (!isTRUE(sp$ok)) stop(sp$reason, call. = FALSE)
      run(list(type = "tolerance", content = input$content,
               confidence = input$confidence, side = input$side,
               scope = input$scope, tol_method = input$tol_method))
    }))
    output$res_title <- shiny::renderText({
      r <- current()
      if (is.null(r)) "No predictions yet" else
        paste0(if (identical(r$type, "tolerance")) "Tolerance intervals"
               else "Predictions", " - ", r$fit$label,
          " (", r$out, ")")
    })
    output$res_note <- shiny::renderUI({
      r <- current()
      shiny::req(r)
      v <- r$value
      ti <- attr(v, "tolerance")
      if (!is.null(ti)) {
        return(help_text(
          "Content ", ti$content %||% "", ", confidence ",
          ti$confidence %||% "", "; method used: ", ti$method %||% "",
          ". min_content is the smallest content over the evaluated ",
          "parameter points, not a certified global minimum."))
      }
      if (inherits(r$fit$fit, "countreg_mix")) help_text(
        "expected_count is the mixture mean (class probabilities from the ",
        "membership variables; no outcomes used); mean.*: class expected ",
        "counts; prob.*: class probabilities.")
    })
    output$res <- shiny::renderTable({
      r <- current()
      shiny::req(r)
      as.data.frame(r$value)
    }, digits = 4)
    output$code <- shiny::renderText({
      r <- current()
      shiny::req(r)
      paste(vapply(r$exprs, .deparse, character(1)), collapse = "\n")
    })
  })
}
