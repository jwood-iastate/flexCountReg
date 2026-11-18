# Build page: the model specification. Editing controls changes only the
# draft; estimation starts only when "Estimate model" is pressed.

mod_build_ui <- function(id) {
  ns <- shiny::NS(id)
  ft <- fcr_family_table()
  fam_choices <- stats::setNames(ft$id, paste0(ft$name, " (", ft$id, ")"))
  none <- c("(none)" = "")
  js_in <- function(ids) {
    paste0("['", paste(ids, collapse = "','"), "'].indexOf(input.estimator) >= 0")
  }
  family_models <- .app_estimators$id[.app_estimators$uses_family]
  mixture_models <- c("countreg_mix", "countreg_lc")
  panel_models <- .app_estimators$id[.app_estimators$group != "none"]
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 380,
      shiny::selectInput(ns("estimator"), "Model", .app_estimator_choices()),
      shiny::uiOutput(ns("estimator_help")),
      shiny::conditionalPanel(
        js_in(family_models), ns = ns,
        shiny::selectizeInput(ns("family"), "Count distribution (family)",
                              fam_choices, selected = "NB2")),
      shiny::conditionalPanel(
        "input.estimator == 'negmultinom'", ns = ns,
        shiny::selectInput(ns("nm_family"), "Negative multinomial form",
                           c("NM2" = "NM2", "NM1" = "NM1", "NM-P" = "NMP"))),
      shiny::conditionalPanel(
        js_in(mixture_models), ns = ns,
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
      help_text("Frequency weights count each row (or panel) that many ",
                "times; they are not survey weights."),
      shiny::conditionalPanel(
        js_in(panel_models), ns = ns,
        shiny::selectInput(ns("group"), "Panel identifier", none),
        shiny::uiOutput(ns("group_help"))),
      shiny::uiOutput(ns("dis1_ui")),
      shiny::conditionalPanel(
        "input.estimator == 'countreg_rp'", ns = ns,
        shiny::selectizeInput(ns("rpar"), "Random parameters", character(),
                              multiple = TRUE),
        help_text("Predictors whose coefficients vary across observations ",
                  "(or panels). They are taken out of the fixed part."),
        shiny::checkboxInput(ns("rpar_intercept"), "Random intercept"),
        shiny::selectInput(ns("rpar_dist"), "Distribution of the random parameters",
                           c("Normal" = "n", "Lognormal" = "ln",
                             "Triangular" = "t", "Uniform" = "u",
                             "Gamma" = "g")),
        shiny::checkboxInput(ns("correlated"),
                             "Correlated random parameters (normal only)"),
        shiny::numericInput(ns("ndraws"), "Halton draws", 500, min = 10,
                            step = 50)),
      shiny::conditionalPanel(
        js_in(mixture_models), ns = ns,
        shiny::selectizeInput(ns("mixing"), "Membership predictors",
                              character(), multiple = TRUE),
        help_text("Variables that change the class probabilities ",
                  "(a multinomial logit), separate from the count ",
                  "predictors."),
        shiny::selectizeInput(ns("shared_mean"),
                              "Same coefficient in every class",
                              character(), multiple = TRUE),
        shiny::checkboxInput(ns("shared_dis1"),
                             "Same distribution parameter in every class")),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Advanced settings",
          shiny::numericInput(ns("n_starts"), "Starts (mixture models)", 10,
                              min = 1, step = 1),
          shiny::numericInput(ns("seed"), "Random seed (mixture models)", 1,
                              step = 1),
          shiny::selectInput(ns("method"), "Optimizer",
                             c("BFGS", "NR", "BHHH", "NM", "CG")),
          shiny::numericInput(ns("max_iters"), "Iteration limit", 1000,
                              min = 1, step = 1),
          shiny::selectInput(ns("stderr"), "Standard errors",
                             c("Observed information (Hessian)" = "normal",
                               "Robust (sandwich)" = "robust")))),
      shiny::textInput(ns("label"), "Fit label (optional)",
                       placeholder = "e.g. NB2 base"),
      shiny::actionButton(ns("estimate"), "Estimate model",
                          class = "btn-primary")),
    bslib::card(bslib::card_header("Model"), fill = FALSE,
                shiny::uiOutput(ns("summary"))),
    bslib::card(bslib::card_header("Checks before fitting"), fill = FALSE,
                shiny::uiOutput(ns("checks"))))
}

mod_build_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    none <- c("(none)" = "")
    last_object <- NULL
    # New data: reset the column choices. A new variable in the same dataset
    # only extends the choices and keeps the current selections.
    shiny::observeEvent(state$data_meta, app_try(what = "Updating choices", {
      d <- state$data
      cols <- names(d)
      cl <- column_summary(d)
      counts <- cl$column[cl$count_like]
      m <- state$data_meta
      same <- identical(m$object, last_object)
      last_object <<- m$object
      keep <- function(v, default) {
        if (same) intersect(v %||% character(), cols) else default
      }
      resp <- if (same && isTRUE(input$response %in% cols)) input$response else
        if (!is.null(m$response) && m$response %in% cols) m$response else
          counts[1]
      shiny::updateSelectInput(session, "response", choices = cols,
                               selected = resp)
      for (nm in c("predictors", "factors", "mixing", "rpar")) {
        shiny::updateSelectizeInput(session, nm, choices = cols,
                                    selected = keep(input[[nm]], character()))
      }
      num <- cl$column[cl$type %in% c("numeric", "integer")]
      for (nm in c("exposure", "log_offset", "weights")) {
        shiny::updateSelectInput(session, nm, choices = c(none, num),
                                 selected = keep(input[[nm]], ""))
      }
      g <- if (!is.null(m$group) && !is.na(m$group) && m$group %in% cols)
        m$group else ""
      shiny::updateSelectInput(session, "group", choices = c(none, cols),
                               selected = if (same) input$group %||% "" else g)
    }))
    # Shared terms and random parameters are chosen among the predictors
    shiny::observeEvent(input$predictors, ignoreNULL = FALSE, {
      shiny::updateSelectizeInput(
        session, "shared_mean",
        choices = c("(Intercept)", input$predictors),
        selected = intersect(input$shared_mean,
                             c("(Intercept)", input$predictors)))
      shiny::updateSelectizeInput(
        session, "rpar", choices = input$predictors,
        selected = intersect(input$rpar, input$predictors))
    })

    families <- shiny::reactive({
      if (.app_is_mixture(input$estimator) && length(input$family_each)) {
        input$family_each
      } else input$family
    })
    caps <- shiny::reactive(fcr_capabilities(input$estimator, families()))

    output$estimator_help <- shiny::renderUI({
      txt <- switch(input$estimator,
        countreg = "Cross-sectional regression with any count family.",
        countreg_rp = paste(
          "Coefficients that vary randomly across observations (or panels),",
          "estimated by maximum simulated likelihood."),
        renb = , poisLind.re = , repois = paste(
          "A random effect shared by the rows of each panel."),
        negmultinom = paste(
          "Panel model with a gamma random effect shared within each panel",
          "(negative multinomial)."),
        countreg_mix = paste(
          "A mixture of count distributions or regressions; components",
          "are a flexible description of the outcome."),
        countreg_lc = paste(
          "The same mixture model, read as latent classes of units."))
      help_text(txt)
    })
    output$group_help <- shiny::renderUI({
      g <- .app_estimator(input$estimator)$group
      help_text(if (identical(g, "required"))
        "Required: the column that identifies the panels." else
          if (.app_is_mixture(input$estimator))
            "Optional: one class per panel." else
              "Optional: each panel shares one draw of the random parameters.")
    })
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
          "This model has no distribution parameter to model; its predictors were cleared.",
          type = "warning")
      }
    })

    spec <- shiny::reactive({
      shiny::req(state$data_meta)
      opt <- function(v) if (is.null(v) || !nzchar(v)) NULL else v
      num <- function(v, default) if (is.null(v) || is.na(v)) default else v
      m <- state$data_meta
      est <- input$estimator
      mixture <- .app_is_mixture(est)
      grp <- if (.app_estimator(est)$group != "none") opt(input$group)
      app_spec(
        data = list(source = m$source, name = m$name, object = m$object),
        recipe = list(response = opt(input$response),
                      predictors = input$predictors %||% character(),
                      factors = input$factors %||% character(),
                      exposure = opt(input$exposure),
                      log_offset = opt(input$log_offset),
                      weights = opt(input$weights), group = grp),
        model = list(
          estimator = est,
          family = if (caps()$uses_family) families() else "NB2",
          components = if (mixture && length(input$family_each))
            length(input$family_each) else num(input$components, 2),
          dis1 = if (isTRUE(caps()$dis1)) input$dis1 %||% character() else
            character(),
          mixing = if (mixture) input$mixing %||% character() else
            character(),
          shared_mean = if (mixture) input$shared_mean %||% character() else
            character(),
          shared_dis1 = mixture && isTRUE(input$shared_dis1) &&
            isTRUE(caps()$dis1),
          rpar = if (est == "countreg_rp") input$rpar %||% character() else
            character(),
          rpar_intercept = est == "countreg_rp" &&
            isTRUE(input$rpar_intercept),
          rpar_dist = input$rpar_dist %||% "n",
          correlated = est == "countreg_rp" && isTRUE(input$correlated),
          ndraws = num(input$ndraws, 500), nm_family = input$nm_family %||%
            "NM2",
          n_starts = num(input$n_starts, 10), seed = num(input$seed, 1),
          method = input$method, max_iters = num(input$max_iters, 1000),
          stderr = input$stderr))
    })
    shiny::observe(state$spec <- tryCatch(spec(), error = function(e) NULL))

    output$summary <- shiny::renderUI({
      s <- spec()
      ft <- fcr_family_table()
      fam <- ft[match(s$model$family, ft$id), ]
      f <- tryCatch(paste(deparse(compile_formula(s)), collapse = " "),
                    error = function(e) "(choose a response)")
      shiny::tagList(
        shiny::p(shiny::strong("Model: "),
                 .app_estimator(s$model$estimator)$label),
        shiny::p(shiny::strong("Formula: "), shiny::code(f)),
        if (s$model$estimator == "countreg_rp") shiny::p(
          shiny::strong("Random parameters: "),
          toString(c(if (isTRUE(s$model$rpar_intercept)) "(Intercept)",
                     s$model$rpar))),
        if (caps()$uses_family) shiny::p(
          shiny::strong("Family: "),
          paste(unique(paste0(fam$name, " (", fam$id, ")")),
                collapse = "; ")),
        if (caps()$uses_family && any(!is.na(fam$param1))) shiny::p(
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
        shiny::showModal(shiny::modalDialog(
          title = "The model cannot be estimated yet",
          shiny::tags$ul(lapply(p, shiny::tags$li)), easyClose = TRUE))
        return()
      }
      job <- app_try(app_submit_fit(state, s, input$label),
                     "Starting the fit")
      if (is.null(job)) return()
      shiny::updateTextInput(session, "label", value = "")
      shiny::showNotification(paste("Estimating", job$label, "..."),
                              type = "message")
      app_go_to(session, "Fit")
    })
  })
}
