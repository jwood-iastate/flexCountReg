# Fit page: jobs, saved fits and the results of the selected fit.

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
      bslib::nav_panel("Messages", shiny::uiOutput(ns("messages")))))
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
  se <- m$se %||% sqrt(pmax(diag(m$vcov), NA))
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
    shiny::observeEvent(input$cancel, {
      id <- input$cancel_id
      if (!nzchar(id %||% "")) return()
      jl <- state$jobs
      jl[[id]] <- cancel_fit_job(jl[[id]])
      state$jobs <- jl
      shiny::showNotification(paste("Cancelled", jl[[id]]$label),
                              type = "warning")
    })
    output$jobs <- shiny::renderTable({
      jobs <- state$jobs
      if (!length(jobs)) return(data.frame(job = "none yet"))
      data.frame(
        job = vapply(jobs, `[[`, "", "label"),
        status = vapply(jobs, `[[`, "", "status"),
        seconds = vapply(jobs, function(j) {
          end <- j$finished %||% Sys.time()
          round(as.numeric(difftime(end, j$started, units = "secs")), 1)
        }, numeric(1)))
    })

    bundle <- shiny::reactive({
      shiny::req(state$active)
      state$fits[[state$active]]
    })
    output$summary <- shiny::renderPrint({
      b <- bundle()
      cat("Fit:", b$label, "| request", b$id, "| data:", b$data_name, "\n\n")
      print(summary(b$fit))
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
