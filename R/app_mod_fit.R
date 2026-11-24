# Fit page: jobs, saved fits and the results of the selected fit, including
# its marginal effects (computed on request, as a background job).

mod_fit_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      shiny::selectInput(ns("active"), "Selected fit", character()),
      help_text("Results on every page refer to this fit and its data."),
      shiny::h6("Jobs"),
      shiny::tableOutput(ns("jobs")),
      shiny::selectInput(ns("cancel_id"), "Running job", character()),
      shiny::actionButton(ns("cancel"), "Cancel job",
                          class = "btn-outline-danger btn-sm")),
    bslib::navset_card_tab(
      bslib::nav_panel("Summary", shiny::verbatimTextOutput(ns("summary"))),
      bslib::nav_panel("Coefficients", DT::DTOutput(ns("coefs"))),
      bslib::nav_panel("Numerical checks", shiny::uiOutput(ns("checks"))),
      bslib::nav_panel("Classes", shiny::uiOutput(ns("classes"))),
      bslib::nav_panel("Marginal effects", margeff_ui(ns)),
      bslib::nav_panel("Messages", shiny::uiOutput(ns("messages")))))
}

# Controls and results of the marginal effects of the selected fit (see
# ?flexCountReg::margEffTable)
margeff_ui <- function(ns) {
  bslib::layout_columns(
    col_widths = c(4, 8),
    shiny::div(
      shiny::selectInput(ns("me_measure"), "Measure",
                         c("Automatic (AME; indicators as below)" = "auto",
                           "Average marginal effect (AME)" = "ame",
                           "Elasticity" = "elasticity")),
      shiny::selectInput(ns("me_indicator"), "Indicator (0/1) variables",
                         c("Pseudo-elasticity (%)" = "pseudo_elasticity",
                           "Discrete change" = "discrete_change")),
      shiny::selectInput(ns("me_se"), "Standard errors",
                         c("Delta method" = "delta",
                           "Bootstrap" = "bootstrap")),
      shiny::conditionalPanel(
        "input.me_se == 'bootstrap'", ns = ns,
        shiny::numericInput(ns("me_bootstraps"), "Bootstrap resamples", 200,
                            min = 10, step = 50),
        help_text("Each request gets its own seed, written into the code, ",
                  "so the result can be reproduced.")),
      shiny::selectizeInput(ns("me_vars"), "Variables", character(),
                            multiple = TRUE,
                            options = list(placeholder =
                                             "all variables of the model")),
      help_text("Effects on the expected count, averaged over the rows ",
                "used in estimation. See ",
                shiny::code("?flexCountReg::margEffTable"), "."),
      shiny::actionButton(ns("me_go"), "Compute marginal effects",
                          class = "btn-primary")),
    shiny::div(
      shiny::uiOutput(ns("me_status")),
      shiny::uiOutput(ns("me_title")),
      DT::DTOutput(ns("me_table")),
      shiny::uiOutput(ns("me_notes")),
      shiny::downloadButton(ns("me_download"), "Download table (CSV)",
                            class = "btn-sm"),
      shiny::verbatimTextOutput(ns("me_code"))))
}

# Variables a marginal effect can be asked for: the columns the fit uses,
# except offsets (which margEffTable() also leaves out by default)
margeff_var_choices <- function(spec) {
  setdiff(needed_columns(spec), c(spec$recipe$exposure,
                                  spec$recipe$log_offset, spec$recipe$group,
                                  spec$recipe$group2))
}

# A marginal-effects request from the page's settings
margeff_request <- function(measure, indicator, se, bootstraps = 200,
                            vars = NULL) {
  if (identical(se, "bootstrap") &&
      (!is.numeric(bootstraps) || length(bootstraps) != 1 ||
       is.na(bootstraps) || bootstraps < 1)) {
    stop("The number of bootstrap resamples must be at least 1.",
         call. = FALSE)
  }
  list(type = "margeff", measure = measure, indicator = indicator, se = se,
       bootstraps = if (identical(se, "bootstrap")) as.integer(bootstraps),
       vars = if (length(vars)) vars)
}

# Numerical status of a fit, with unavailable checks shown as such
fit_checks <- function(fit) {
  m <- fit$model
  num <- m$numerical
  if (is.null(num)) {
    H <- m$hessian
    ev <- if (is.matrix(H) && all(is.finite(H)))
      eigen(-(H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
    num <- list(
      max_abs_gradient = if (!is.null(m$gradient)) max(abs(m$gradient)) else
        NA_real_,
      min_eigenvalue = if (length(ev)) min(ev) else NA_real_,
      information = if (!length(ev)) "unavailable" else
        if (min(ev) > sqrt(.Machine$double.eps) * max(abs(ev)))
          "positive definite" else "not positive definite",
      condition = if (length(ev) && min(ev) > 0) max(ev) / min(ev) else
        NA_real_)
  }
  # Standard errors from the covariance matrix when the object has none;
  # negative or non-finite variances give NA (unavailable)
  se <- m$se %||% .fcr_se(m$vcov, length(m$estimate))
  data.frame(
    check = c("Optimizer reported convergence", "Optimizer message",
              "Iterations", "Log-likelihood", "Largest absolute gradient",
              "Observed information", "Smallest information eigenvalue",
              "Condition number", "Unavailable standard errors"),
    value = c(format(isTRUE(m$converged)), trimws(m$message %||% "unavailable"),
              format(m$iterations %||% NA), format(round(m$maximum, 4)),
              format(signif(num$max_abs_gradient, 3)), num$information,
              format(signif(num$min_eigenvalue, 3)),
              if (is.finite(num$condition)) format(signif(num$condition, 3))
              else "not meaningful (information not positive definite)",
              if (any(is.na(se))) toString(names(m$estimate)[is.na(se)]) else
                "none"),
    stringsAsFactors = FALSE)
}

mod_fit_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observe({
      fits <- state$fits
      ch <- label_choices(fits)
      shiny::updateSelectInput(session, "active", choices = ch,
                               selected = state$active)
    })
    shiny::observeEvent(input$active, {
      if (nzchar(input$active)) state$active <- input$active
    })
    shiny::observe({
      run <- Filter(function(j) identical(j$status, "running"), state$jobs)
      shiny::updateSelectInput(session, "cancel_id",
                               choices = label_choices(run))
    })
    shiny::observeEvent(input$cancel, app_try(what = "Canceling", {
      id <- input$cancel_id
      if (!nzchar(id %||% "")) return()
      jl <- state$jobs
      jl[[id]] <- cancel_fit_job(jl[[id]])
      state$jobs <- jl
      shiny::showNotification(paste("Canceled", jl[[id]]$label),
                              type = "warning")
    }))
    output$jobs <- shiny::renderTable({
      jobs <- state$jobs
      if (!length(jobs)) return(data.frame(job = "none yet"))
      data.frame(
        job = vapply(jobs, `[[`, "", "label"),
        status = vapply(jobs, function(j) {
          if (identical(j$status, "failed")) {
            paste("failed:", j$result$error %||% "")
          } else j$status
        }, ""),
        seconds = vapply(jobs, function(j) {
          end <- j$finished %||% Sys.time()
          round(as.numeric(difftime(end, j$started, units = "secs")), 1)
        }, numeric(1)))
    })

    bundle <- shiny::reactive({
      shiny::req(state$active)
      b <- state$fits[[state$active]]
      shiny::req(b)
      b$fit <- app_fit_object(state, b$id)
      b
    })
    output$summary <- shiny::renderPrint({
      b <- bundle()
      cat(b$label, "| object", b$fit_name, "| data:", b$data_name, "\n\n")
      print(get(paste0(b$fit_name, "_summary"), envir = state$env))
    })
    output$coefs <- DT::renderDT({
      td <- flexCountReg::tidy(bundle()$fit, conf.int = TRUE)
      num <- vapply(td, is.numeric, logical(1))
      td[num] <- lapply(td[num], signif, 4)
      DT::datatable(as.data.frame(td), rownames = FALSE,
                    options = list(pageLength = 25, scrollX = TRUE))
    })
    output$checks <- shiny::renderUI({
      b <- bundle()
      tab <- fit_checks(b$fit)
      shiny::tagList(
        shiny::p("Convergence alone is not enough: check the gradient, the ",
                 "observed information and standard errors that are not ",
                 "available. See ", shiny::code(
                   "vignette(\"numerical-reliability\", package = \"flexCountReg\")"),
                 "."),
        shiny::renderTable(tab))
    })
    output$classes <- shiny::renderUI({
      f <- bundle()$fit
      if (!inherits(f, "countreg_mix")) {
        return(shiny::p("Class results are shown for finite-mixture and ",
                        "latent-class fits."))
      }
      cl <- f$model$classes
      st <- f$model$starts
      shiny::tagList(
        shiny::h6("Class sizes"),
        shiny::p(class = "small", "prior_share: mean membership probability ",
                 "(from covariates); posterior_share: mean posterior ",
                 "probability (given the outcomes); effective_size: sum of ",
                 "posterior probabilities; modal_count: units whose most ",
                 "probable class it is (an uncertain classification)."),
        shiny::renderTable(cl, digits = 4),
        shiny::p("Entropy: ", round(f$model$entropy, 3),
                 " | relative entropy (1 = perfect classification): ",
                 round(f$model$relative_entropy, 3)),
        shiny::h6("Starts"),
        shiny::p(class = "small", sum(st$reached_best), " of ", nrow(st),
                 " starts reached the selected log-likelihood. A solution ",
                 "found by one start only deserves caution."),
        shiny::renderTable(st, digits = 4))
    })
    # Marginal effects: a background job of the selected fit. Its result,
    # or the error it stopped with, is shown in the panel.
    me_pending <- shiny::reactiveVal(NULL)
    me_result <- shiny::reactiveVal(NULL)
    me_error <- shiny::reactiveVal(NULL)
    shiny::observeEvent(bundle(), {
      ch <- margeff_var_choices(bundle()$spec)
      shiny::updateSelectizeInput(session, "me_vars", choices = ch,
                                  selected = intersect(input$me_vars, ch))
    })
    shiny::observeEvent(input$me_go, {
      me_error(NULL)
      tryCatch({
        b <- bundle()
        req <- margeff_request(input$me_measure, input$me_indicator,
                               input$me_se, input$me_bootstraps,
                               input$me_vars)
        job <- app_start_output_job(state, req, b$id)
        me_result(NULL)
        me_pending(job$id)
      }, error = function(e) me_error(conditionMessage(e)))
    })
    shiny::observe({
      id <- me_pending()
      if (is.null(id)) return()
      res <- state$outputs[[id]]
      if (!is.null(res)) {
        me_result(res)
        me_pending(NULL)
        return()
      }
      j <- state$jobs[[id]]
      if (!is.null(j) && j$status %in% c("failed", "canceled")) {
        me_error(if (identical(j$status, "failed"))
          j$result$error %||% "The computation failed." else "Canceled.")
        me_pending(NULL)
      }
    })
    shiny::observeEvent(input$me_cancel, {
      id <- me_pending()
      if (is.null(id)) return()
      jl <- state$jobs
      jl[[id]] <- cancel_fit_job(jl[[id]])
      state$jobs <- jl
      me_error("Canceled.")
      me_pending(NULL)
    })
    output$me_status <- shiny::renderUI({
      id <- me_pending()
      err <- me_error()
      if (!is.null(id)) {
        return(shiny::div(
          class = "alert alert-info d-flex align-items-center gap-3",
          role = "status",
          shiny::span(paste0("Running ", state$jobs[[id]]$label,
                             " in the background...")),
          shiny::actionButton(session$ns("me_cancel"), "Cancel",
                              class = "btn-sm btn-outline-danger")))
      }
      if (!is.null(err)) {
        shiny::div(class = "alert alert-danger", role = "alert",
                   shiny::strong("Marginal effects not computed: "), err)
      }
    })
    output$me_title <- shiny::renderUI({
      r <- me_result()
      if (is.null(r)) return(shiny::p("No marginal effects yet."))
      shiny::h6(paste0("Marginal effects - ", r$fit$label, " (", r$out, ")"))
    })
    output$me_table <- DT::renderDT({
      r <- me_result()
      shiny::req(r)
      # Rounded for display only; the object and the download keep every
      # digit
      td <- as.data.frame(r$value)
      num <- vapply(td, is.numeric, logical(1))
      td[num] <- lapply(td[num], signif, 4)
      DT::datatable(td, rownames = FALSE,
                    options = list(dom = "t", pageLength = 100,
                                   scrollX = TRUE))
    })
    output$me_notes <- shiny::renderUI(output_notes(me_result()))
    output$me_code <- shiny::renderText({
      r <- me_result()
      shiny::req(r)
      paste(vapply(r$exprs, .deparse, character(1)), collapse = "\n")
    })
    output$me_download <- shiny::downloadHandler(
      filename = function() paste0(me_result()$out %||% "margeff", ".csv"),
      content = function(file) {
        utils::write.csv(as.data.frame(me_result()$value), file,
                         row.names = FALSE)
      })
    output$messages <- shiny::renderUI({
      b <- bundle()
      shiny::tagList(
        shiny::h6("Warnings"),
        if (length(b$warnings)) shiny::tags$ul(lapply(b$warnings, shiny::tags$li))
        else shiny::p("None."),
        shiny::h6("Messages"),
        if (length(b$messages)) shiny::tags$ul(lapply(b$messages, shiny::tags$li))
        else shiny::p("None."),
        shiny::p(class = "small text-muted", "Elapsed: ",
                 round(b$elapsed, 2), " s; versions: ",
                 paste(names(b$versions), b$versions, collapse = ", ")))
    })
  })
}
