# Predict page: an editable table of scenarios; predictions, confidence or
# prediction intervals, and tolerance intervals of the selected fit. Each
# calculation runs its compiled code, which is added to the session script.

mod_predict_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340,
      shiny::selectInput(ns("fit"), "Model", character()),
      shiny::uiOutput(ns("fit_note")),
      shiny::radioButtons(ns("rows"), "Rows to predict",
                          c("Scenario table (edit below)" = "table",
                            "New data (CSV file)" = "newdata")),
      shiny::conditionalPanel(
        "input.rows == 'newdata'", ns = ns,
        shiny::fileInput(ns("nd_file"), "Upload new data (CSV)",
                         accept = c(".csv", "text/csv")),
        shiny::actionButton(ns("nd_load"), "Load new data",
                            class = "btn-outline-primary btn-sm"),
        shiny::selectInput(ns("newdata"), "New data", character()),
        help_text("Variables you created on the Data page are created in ",
                  "the new data in the same way, with the means, standard ",
                  "deviations and categories of the training data.")),
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
                shiny::uiOutput(ns("running")),
                shiny::uiOutput(ns("res_note")),
                shiny::uiOutput(ns("res_warnings")),
                shiny::tableOutput(ns("res")),
                shiny::downloadButton(ns("download"), "Download results",
                                      class = "btn-sm")),
    bslib::card(bslib::card_header("Code that produced it"), fill = FALSE,
                shiny::verbatimTextOutput(ns("code"))))
}

# Columns a fit needs in new data
needed_columns <- function(spec) {
  r <- spec$recipe
  m <- spec$model
  # Hierarchical and fixed-effects models also need the grouping columns
  grp <- if (m$estimator %in% c("countreg_hier", "countreg_fe"))
    c(r$group, r$group2)
  unique(c(r$predictors, r$exposure, r$log_offset, m$dis1, m$mixing,
           m$rpar, grp))
}

# One row of typical values (mean of numeric columns, first level or most
# common value otherwise), after the recipe's factor conversions
typical_row <- function(spec, data) {
  for (f in spec$recipe$factors) data[[f]] <- factor(data[[f]])
  cols <- needed_columns(spec)
  # Grouping columns identify groups by their labels: a character value
  # (the most common group), not a factor with every level
  grp <- intersect(c(spec$recipe$group, spec$recipe$group2), cols)
  if (spec$model$estimator %in% c("countreg_hier", "countreg_fe")) {
    for (g in grp) data[[g]] <- as.character(data[[g]])
  }
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
    app_fit_selector(input, session, state, "fit")
    bundle <- shiny::reactive({
      b <- state$fits[[input$fit %||% ""]]
      shiny::req(b)
      b$fit <- app_fit_object(state, b$id)
      b
    })
    output$fit_note <- shiny::renderUI({
      b <- state$fits[[input$fit %||% ""]]
      help_text(if (is.null(b)) "Estimate a model first (Build page)." else
        paste0("Object ", b$fit_name, "; data ", b$spec$data$object, "."))
    })
    # New data
    shiny::observeEvent(input$nd_load, app_try(what = "Loading new data", {
      f <- input$nd_file
      if (is.null(f)) stop("Choose a CSV file first.", call. = FALSE)
      keep <- file.path(state$session_dir,
                        paste0("newdata-", length(state$newdata) + 1, ".csv"))
      if (!isTRUE(file.copy(f$datapath, keep, overwrite = TRUE))) {
        stop("Could not keep a copy of the uploaded file.", call. = FALSE)
      }
      obj <- app_load_newdata(state, f$name, keep)
      shiny::updateSelectInput(session, "newdata", choices = state$newdata,
                               selected = obj)
      shiny::showNotification(paste("Loaded", f$name, "as", obj),
                              type = "message")
    }))
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
      if (identical(input$rows, "newdata")) {
        nd <- input$newdata
        if (is.null(nd) || !nzchar(nd)) {
          stop("Load new data first.", call. = FALSE)
        }
        req$newdata <- nd
      } else {
        s <- scen()
        if (is.null(s) || !nrow(s)) {
          stop("The scenario table is empty.", call. = FALSE)
        }
        req$scenarios <- s
      }
      # Predictions and intervals run in the background (a parametric
      # bootstrap can take minutes); the result is picked up when the job
      # finishes
      if (!identical(req$interval %||% "none", "none") ||
          identical(req$type, "tolerance")) {
        if (app_needs_bootstrap(b$fit)) {
          shiny::showNotification(paste(
            "The covariance matrix of this model is not positive definite",
            "(e.g. a parameter at a boundary), so the interval uses a",
            "parametric bootstrap with 200 refits. This can take several",
            "minutes; it runs in the background and can be canceled."),
            type = "warning", duration = 15)
        }
      }
      job <- app_start_output_job(state, req, b$id)
      pending(job$id)
    }
    pending <- shiny::reactiveVal(NULL)
    # Pick up the result (or the failure) of the pending job
    shiny::observe({
      id <- pending()
      if (is.null(id)) return()
      res <- state$outputs[[id]]
      if (!is.null(res)) {
        res$fit$fit <- app_fit_object(state, res$fit$id)
        current(res)
        pending(NULL)
        for (m in res$messages) {
          shiny::showNotification(m, type = "message", duration = 10)
        }
        if (length(res$warnings)) {
          shiny::showNotification(paste(
            length(res$warnings), "warning(s); see below the results."),
            type = "warning", duration = 10)
        }
        return()
      }
      j <- state$jobs[[id]]
      if (!is.null(j) && j$status %in% c("failed", "canceled")) pending(NULL)
    })
    output$running <- shiny::renderUI({
      id <- pending()
      if (is.null(id)) return(NULL)
      j <- state$jobs[[id]]
      shiny::div(
        class = "alert alert-info d-flex align-items-center gap-3",
        role = "status",
        shiny::span(paste0("Running ", j$label, " in the background...")),
        shiny::actionButton(session$ns("cancel"), "Cancel",
                            class = "btn-sm btn-outline-danger"))
    })
    shiny::observeEvent(input$cancel, app_try(what = "Canceling", {
      id <- pending()
      if (is.null(id)) return()
      jl <- state$jobs
      jl[[id]] <- cancel_fit_job(jl[[id]])
      state$jobs <- jl
      pending(NULL)
      shiny::showNotification("Canceled.", type = "warning")
    }))
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
      n <- nrow(as.data.frame(v))
      if (n > 200) return(help_text("Showing the first 200 of ", n,
                                    " rows; download for all of them."))
      if (inherits(r$fit$fit, "countreg_mix")) help_text(
        "expected_count is the mixture mean (class probabilities from the ",
        "membership variables; no outcomes used); mean.*: class expected ",
        "counts; prob.*: class probabilities.")
    })
    # Warnings, messages and the uncertainty settings stay beside the
    # result (and are written to the manifest of the export)
    output$res_warnings <- shiny::renderUI(output_notes(current()))
    output$res <- shiny::renderTable({
      r <- current()
      shiny::req(r)
      utils::head(as.data.frame(r$value), 200)
    }, digits = 4)
    output$download <- shiny::downloadHandler(
      filename = function() paste0(current()$out %||% "predictions", ".csv"),
      content = function(file) {
        utils::write.csv(as.data.frame(current()$value), file,
                         row.names = FALSE)
      })
    output$code <- shiny::renderText({
      r <- current()
      shiny::req(r)
      paste(vapply(r$exprs, .deparse, character(1)), collapse = "\n")
    })
  })
}
