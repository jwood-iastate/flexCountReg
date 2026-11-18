# Build page: the model specification. Editing controls changes only the
# draft; estimation starts only when "Estimate model" is pressed.

mod_build_ui <- function(id) {
  ns <- shiny::NS(id)
  ft <- fcr_family_table()
  fam_choices <- stats::setNames(ft$id, paste0(ft$name, " (", ft$id, ")"))
  none <- c("(none)" = "")
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      shiny::selectInput(ns("estimator"), "Model", .app_estimators),
      shiny::selectizeInput(ns("family"), "Count distribution (family)",
                            fam_choices, selected = "NB2"),
      shiny::conditionalPanel(
        "input.estimator != 'countreg'", ns = ns,
        shiny::numericInput(ns("components"), "Number of classes", 2,
                            min = 1, max = 6, step = 1),
        shiny::selectizeInput(ns("family_each"),
                              "Different family in each class (optional)",
                              fam_choices, multiple = TRUE,
                              options = list(placeholder =
                                               "same family in every class")),
        help_text("Give one family per class, in class order, to mix ",
                  "families; leave empty to use the family above.")),
      shiny::selectInput(ns("response"), "Response (count)", character()),
      shiny::selectizeInput(ns("predictors"), "Predictors", character(),
                            multiple = TRUE),
      shiny::selectizeInput(ns("factors"), "Treat as categorical (factor)",
                            character(), multiple = TRUE),
      shiny::selectInput(ns("exposure"), "Exposure (logged: offset(log(x)))",
                         none),
      shiny::selectInput(ns("log_offset"),
                         "Offset already on the log scale", none),
      shiny::selectInput(ns("weights"), "Frequency weights", none),
      help_text("Frequency weights count each row that many times; they ",
                "are not survey weights."),
      shiny::uiOutput(ns("dis1_ui")),
      shiny::conditionalPanel(
        "input.estimator != 'countreg'", ns = ns,
        shiny::selectizeInput(ns("mixing"), "Membership predictors",
                              character(), multiple = TRUE),
        help_text("Variables that change the class probabilities ",
                  "(a multinomial logit), separate from the count ",
                  "predictors."),
        shiny::selectizeInput(ns("shared_mean"),
                              "Same coefficient in every class",
                              character(), multiple = TRUE),
        shiny::checkboxInput(ns("shared_dis1"),
                             "Same distribution parameter in every class"),
        shiny::selectInput(ns("group"), "Panel identifier (one class per panel)",
                           none)),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Advanced settings",
          shiny::numericInput(ns("n_starts"), "Starts (mixture models)", 10,
                              min = 1, step = 1),
          shiny::numericInput(ns("seed"), "Random seed", 1, step = 1),
          shiny::selectInput(ns("method"), "Optimizer",
                             c("BFGS", "NR", "BHHH", "NM", "CG")),
          shiny::numericInput(ns("max_iters"), "Iteration limit", 1000,
                              min = 1, step = 1),
          shiny::selectInput(ns("stderr"), "Standard errors",
                             c("Observed information (Hessian)" = "normal",
                               "Robust (sandwich)" = "robust")))),
      shiny::textInput(ns("label"), "Fit label", placeholder = "e.g. NB2 base"),
      shiny::actionButton(ns("estimate"), "Estimate model",
                          class = "btn-primary")),
    bslib::card(bslib::card_header("Model"),
                shiny::uiOutput(ns("summary"))),
    bslib::card(bslib::card_header("Checks before fitting"),
                shiny::uiOutput(ns("checks"))))
}

mod_build_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    none <- c("(none)" = "")
    # New data: reset the column choices
    shiny::observeEvent(state$data_meta, {
      d <- state$data
      cols <- names(d)
      cl <- column_summary(d)
      counts <- cl$column[cl$count_like]
      m <- state$data_meta
      resp <- if (!is.null(m$response) && m$response %in% cols) m$response else
        counts[1]
      shiny::updateSelectInput(session, "response", choices = cols,
                               selected = resp)
      for (nm in c("predictors", "factors", "mixing", "shared_mean")) {
        shiny::updateSelectizeInput(session, nm, choices = cols,
                                    selected = character())
      }
      num <- cl$column[cl$type %in% c("numeric", "integer")]
      for (nm in c("exposure", "log_offset", "weights")) {
        shiny::updateSelectInput(session, nm, choices = c(none, num),
                                 selected = "")
      }
      shiny::updateSelectInput(session, "group", choices = c(none, cols),
                               selected = "")
    })
    # Shared terms can only be chosen among the predictors
    shiny::observeEvent(input$predictors, ignoreNULL = FALSE, {
      shiny::updateSelectizeInput(
        session, "shared_mean",
        choices = c("(Intercept)", input$predictors),
        selected = intersect(input$shared_mean,
                             c("(Intercept)", input$predictors)))
    })

    families <- shiny::reactive({
      if (input$estimator != "countreg" && length(input$family_each)) {
        input$family_each
      } else input$family
    })
    caps <- shiny::reactive(fcr_capabilities(input$estimator, families()))

    output$dis1_ui <- shiny::renderUI({
      cp <- caps()
      if (!isTRUE(cp$dis1)) return(NULL)
      shiny::tagList(
        shiny::selectizeInput(session$ns("dis1"),
                              paste0("Predictors of ", cp$dis1_label),
                              names(state$data), multiple = TRUE,
                              selected = shiny::isolate(input$dis1)),
        help_text("Optional: model the distribution parameter (",
                  cp$dis1_label, ") as a function of variables."))
    })
    # A family without a distribution parameter clears its predictors,
    # visibly (no hidden setting stays in the model)
    shiny::observeEvent(caps(), {
      if (!isTRUE(caps()$dis1) && length(input$dis1)) {
        shiny::showNotification(
          paste("The", toString(unique(families())), "family has no",
                "distribution parameter; its predictors were cleared."),
          type = "warning")
      }
    })

    spec <- shiny::reactive({
      shiny::req(state$data_meta)
      opt <- function(v) if (is.null(v) || !nzchar(v)) NULL else v
      m <- state$data_meta
      app_spec(
        data = list(source = m$source, name = m$name),
        recipe = list(response = opt(input$response),
                      predictors = input$predictors %||% character(),
                      factors = input$factors %||% character(),
                      exposure = opt(input$exposure),
                      log_offset = opt(input$log_offset),
                      weights = opt(input$weights),
                      group = if (input$estimator != "countreg")
                        opt(input$group)),
        model = list(
          estimator = input$estimator, family = families(),
          components = if (input$estimator != "countreg" &&
                           length(input$family_each))
            length(input$family_each) else input$components,
          dis1 = if (isTRUE(caps()$dis1)) input$dis1 %||% character() else
            character(),
          mixing = if (input$estimator != "countreg")
            input$mixing %||% character() else character(),
          shared_mean = if (input$estimator != "countreg")
            input$shared_mean %||% character() else character(),
          shared_dis1 = input$estimator != "countreg" &&
            isTRUE(input$shared_dis1) && isTRUE(caps()$dis1),
          n_starts = input$n_starts, seed = input$seed,
          method = input$method, max_iters = input$max_iters,
          stderr = input$stderr))
    })
    shiny::observe(state$spec <- spec())

    output$summary <- shiny::renderUI({
      s <- spec()
      ft <- fcr_family_table()
      fam <- ft[match(s$model$family, ft$id), ]
      f <- tryCatch(paste(deparse(compile_formula(s)), collapse = " "),
                    error = function(e) "(incomplete)")
      shiny::tagList(
        shiny::p(shiny::strong("Formula: "), shiny::code(f)),
        shiny::p(shiny::strong("Family: "),
                 paste(unique(paste0(fam$name, " (", fam$id, ")")),
                       collapse = "; ")),
        if (any(!is.na(fam$param1))) shiny::p(
          "Distribution parameter(s): ",
          toString(unique(stats::na.omit(c(fam$param1, fam$param2)))),
          ". Labels follow ", shiny::code("?flexCountReg::count_families"),
          ": not every parameter is a log dispersion."),
        if (isTRUE(caps()$mixture)) shiny::p(
          "Classes are a modeling device: a good fit with several classes ",
          "does not show that real subpopulations exist."),
        if (isTRUE(caps()$experimental)) shiny::div(
          class = "alert alert-warning", role = "alert",
          "Mixtures of this family are experimental in flexCountReg ",
          "(recovery not yet validated)."))
    })
    output$checks <- shiny::renderUI({
      s <- spec()
      d <- state$data
      p <- validate_spec(s, d)
      used <- intersect(c(s$recipe$response, s$recipe$predictors,
                          s$recipe$exposure, s$recipe$log_offset,
                          s$recipe$weights, s$recipe$group, s$model$dis1,
                          s$model$mixing), names(d))
      n_cc <- sum(stats::complete.cases(d[, used, drop = FALSE]))
      shiny::tagList(
        shiny::p("Analysis sample: ", n_cc, " of ", nrow(d), " rows have ",
                 "no missing values in the columns this model uses."),
        if (length(p)) shiny::div(
          class = "alert alert-danger", role = "alert",
          shiny::strong("Fix before estimating:"),
          shiny::tags$ul(lapply(p, shiny::tags$li))) else
            shiny::p(class = "text-success", "Ready to estimate."))
    })

    shiny::observeEvent(input$estimate, {
      s <- spec()
      p <- validate_spec(s, state$data)
      if (length(p)) {
        shiny::showNotification("Fix the problems listed under 'Checks' first.",
                                type = "error")
        return()
      }
      n <- length(state$jobs) + 1
      label <- if (nzchar(input$label)) input$label else
        paste0("Fit ", n, ": ", paste(unique(s$model$family), collapse = "+"),
               if (s$model$estimator != "countreg")
                 paste0(", K = ", s$model$components))
      job <- start_fit_job(s, data = if (s$data$source != "bundled")
        state$data, dir = state$session_dir, label = label)
      jl <- state$jobs
      jl[[job$id]] <- job
      state$jobs <- jl
      state$last_submitted <- job$id
      shiny::showNotification(paste("Estimating", label, "..."),
                              type = "message")
    })
  })
}
