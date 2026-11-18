# Predict page: an editable table of scenarios and the predictions of the
# selected fit for them.

mod_predict_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      shiny::actionButton(ns("reset"), "Start from typical values",
                          class = "btn-outline-secondary btn-sm"),
      shiny::actionButton(ns("add"), "Add a row",
                          class = "btn-outline-secondary btn-sm"),
      help_text("Double-click a cell to edit it. Numeric columns start at ",
                "their mean, categorical ones at their first level."),
      shiny::selectInput(ns("interval"), "Interval",
                         c("None" = "none",
                           "Confidence interval for the expected count" =
                             "confidence",
                           "Prediction interval for a new count" =
                             "prediction")),
      shiny::numericInput(ns("level"), "Level", 0.95, min = 0.5, max = 0.999,
                          step = 0.01),
      help_text("A confidence interval concerns the expected count; a ",
                "prediction interval concerns one future count. For ",
                "population content use flexCountReg::tolerance_interval()."),
      shiny::actionButton(ns("go"), "Calculate predictions",
                          class = "btn-primary")),
    bslib::card(bslib::card_header("Scenarios"), DT::DTOutput(ns("scen"))),
    bslib::card(bslib::card_header(shiny::textOutput(ns("res_title"))),
                shiny::uiOutput(ns("res_note")),
                DT::DTOutput(ns("res"))))
}

# Columns a fit needs in new data
needed_columns <- function(spec) {
  r <- spec$recipe
  m <- spec$model
  unique(c(r$predictors, r$exposure, r$log_offset, m$dis1, m$mixing))
}

# One row of typical values (mean of numeric columns, first level or most
# common value otherwise), after the recipe's factor conversions
typical_row <- function(spec, data) {
  for (f in spec$recipe$factors) data[[f]] <- factor(data[[f]])
  cols <- needed_columns(spec)
  out <- lapply(cols, function(cn) {
    v <- data[[cn]]
    if (is.factor(v)) factor(levels(v)[1], levels = levels(v))
    else if (is.numeric(v)) mean(v, na.rm = TRUE)
    else names(sort(table(v), decreasing = TRUE))[1]
  })
  if (!length(out)) return(data.frame(row = 1)[, 0, drop = FALSE][1, , drop = FALSE])
  stats::setNames(as.data.frame(out, check.names = FALSE), cols)
}

mod_predict_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    scen <- shiny::reactiveVal(NULL)
    bundle <- shiny::reactive({
      shiny::req(state$active)
      state$fits[[state$active]]
    })
    shiny::observeEvent(list(bundle(), input$reset), {
      b <- bundle()
      scen(typical_row(b$spec, state$data))
    })
    shiny::observeEvent(input$add, {
      s <- scen()
      if (!is.null(s) && nrow(s)) scen(rbind(s, s[nrow(s), , drop = FALSE]))
    })
    output$scen <- DT::renderDT({
      DT::datatable(scen(), rownames = FALSE, editable = "cell",
                    options = list(dom = "t", pageLength = 50))
    })
    shiny::observeEvent(input$scen_cell_edit, {
      s <- scen()
      e <- input$scen_cell_edit
      j <- e$col + 1
      v <- s[[j]]
      s[e$row, j] <- if (is.factor(v)) {
        if (!e$value %in% levels(v)) {
          shiny::showNotification(paste0("Unknown level \"", e$value,
                                         "\"; levels: ",
                                         toString(levels(v))),
                                  type = "error")
          return()
        }
        e$value
      } else if (is.numeric(v)) as.numeric(e$value) else e$value
      scen(s)
    })
    shiny::observe(state$scenarios <- scen())

    res <- shiny::eventReactive(input$go, {
      b <- bundle()
      f <- b$fit
      s <- scen()
      mixture <- inherits(f, "countreg_mix")
      if (input$interval != "none" && mixture) {
        shiny::showNotification(fcr_capabilities("countreg_mix",
                                                 "POISSON")$intervals_reason,
                                type = "warning")
      }
      out <- s
      if (input$interval != "none" && !mixture) {
        p <- stats::predict(f, newdata = s, interval = input$interval,
                            level = input$level)
        out <- cbind(s, p)
        attr(out, "note") <- paste0(
          if (input$interval == "confidence") "Confidence" else "Prediction",
          " intervals at level ", input$level, " (", attr(p, "interval")$uncertainty %||%
            "normal", " parameter uncertainty).")
      } else {
        out$expected_count <- stats::predict(f, newdata = s)
        if (mixture) {
          comp <- stats::predict(f, newdata = s, type = "component")
          colnames(comp) <- paste0("mean_", colnames(comp))
          pr <- flexCountReg::component_probs(f, newdata = s)
          colnames(pr) <- paste0("prob_", colnames(pr))
          out <- cbind(out, comp, pr)
          attr(out, "note") <- paste(
            "expected_count is the mixture mean (prior class probabilities",
            "from the membership variables; no outcomes used). mean_*: class",
            "expected counts; prob_*: class probabilities.")
        }
      }
      list(data = out, fit = b)
    })
    output$res_title <- shiny::renderText(paste("Predictions - fit:",
                                                res()$fit$label))
    output$res_note <- shiny::renderUI({
      n <- attr(res()$data, "note")
      if (!is.null(n)) shiny::p(class = "small", n)
    })
    output$res <- DT::renderDT({
      d <- res()$data
      num <- vapply(d, is.numeric, logical(1))
      d[num] <- lapply(d[num], signif, 5)
      DT::datatable(d, rownames = FALSE, options = list(dom = "tip"))
    })
  })
}
