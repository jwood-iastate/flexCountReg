# Module composition and the session's background-job polling. All state
# lives in one reactiveValues object per session (app_session.R), never in
# global variables, so concurrent sessions cannot see each other's data,
# jobs or fits.

app_server <- function(input, output, session) {
  state <- new_app_state()
  dir.create(shiny::isolate(state$session_dir), showWarnings = FALSE)

  mod_data_server("data", state)
  mod_build_server("build", state)
  mod_fit_server("fit", state)
  mod_diagnose_server("diagnose", state)
  mod_compare_server("compare", state)
  mod_predict_server("predict", state)
  mod_export_server("export", state)

  # Poll running jobs
  shiny::observe({
    jobs <- state$jobs
    running <- vapply(jobs, function(j) identical(j$status, "running"),
                      logical(1))
    if (!any(running)) return()
    shiny::invalidateLater(500)
    for (id in names(jobs)[running]) {
      j <- poll_fit_job(jobs[[id]])
      if (identical(j$status, "running")) next
      app_try(shiny::isolate(record_finished_job(state, j)),
              paste("Recording", j$label))
      if (identical(j$status, "completed")) {
        shiny::showNotification(paste(j$label, "completed."),
                                type = "message")
      } else {
        shiny::showNotification(
          paste0(j$label, " ", j$status, ": ", j$result$error %||% ""),
          type = "error", duration = NULL)
      }
    }
  })

  output$status_bar <- shiny::renderText({
    dm <- state$data_meta
    act <- state$fits[[state$active %||% ""]]
    n_run <- sum(vapply(state$jobs, function(j)
      identical(j$status, "running"), logical(1)))
    paste0("Data: ", if (is.null(dm)) "none" else
      paste0(dm$name, " (", dm$object, ")"),
      " | Selected fit: ", if (is.null(act)) "none" else act$label,
      " | Running jobs: ", n_run)
  })

  output$session_code <- shiny::renderText({
    if (!length(state$log)) {
      return("Nothing has been run yet. Load data on the Data page.")
    }
    app_script(state)
  })
  output$draft_code <- shiny::renderText({
    spec <- state$spec
    if (is.null(spec)) return("Choose data and a model on the Build page.")
    probs <- validate_spec(spec, state$data)
    code <- tryCatch(draft_script(spec, "fit_next"),
                     error = function(e) "(incomplete specification)")
    if (length(probs)) {
      paste0("# Unresolved problems:\n", paste0("#  - ", probs,
                                                collapse = "\n"),
             "\n\n", code)
    } else code
  })
  # The code drawer is collapsed at first; keep its (cheap) outputs current
  shiny::outputOptions(output, "session_code", suspendWhenHidden = FALSE)
  shiny::outputOptions(output, "draft_code", suspendWhenHidden = FALSE)

  session$onSessionEnded(function() {
    jobs <- shiny::isolate(state$jobs)
    for (j in jobs) cancel_fit_job(j)
    unlink(shiny::isolate(state$session_dir), recursive = TRUE)
  })
}

# Switch the app to a page from inside a module
app_go_to <- function(session, page) {
  bslib::nav_select("page", page, session = session$rootScope())
}
