# Session state and module composition. All state lives in one
# reactiveValues object per session (never in global variables), so
# concurrent sessions cannot see each other's data, jobs or fits.

new_app_state <- function(session_dir = tempfile("fcr-app-")) {
  shiny::reactiveValues(
    data = NULL,            # source data frame
    data_meta = NULL,       # list(source, name, title)
    spec = NULL,            # draft app_spec (from the Build page)
    jobs = list(),          # running and finished jobs, by request id
    fits = list(),          # completed fit bundles, by request id
    active = NULL,          # id of the selected fit
    last_submitted = NULL,  # id of the most recent submission
    scenarios = NULL,       # scenario table of the Predict page
    session_dir = session_dir)
}

# Record a finished job: store it, keep a completed fit as an immutable
# bundle, and make it the active fit only if it is the most recent
# submission (so a slow earlier job never replaces a newer result).
record_finished_job <- function(state, j) {
  jl <- state$jobs
  jl[[j$id]] <- j
  state$jobs <- jl
  if (identical(j$status, "completed")) {
    fl <- state$fits
    fl[[j$id]] <- make_fit_bundle(j, state$data_meta$name)
    state$fits <- fl
    if (identical(j$id, state$last_submitted)) state$active <- j$id
  }
  invisible(j$status)
}

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

  # Poll running jobs; a finished job becomes the active fit only if it is
  # still the most recent submission (a slow earlier job does not replace a
  # newer result)
  shiny::observe({
    jobs <- state$jobs
    running <- vapply(jobs, function(j) identical(j$status, "running"),
                      logical(1))
    if (!any(running)) return()
    shiny::invalidateLater(500)
    for (id in names(jobs)[running]) {
      j <- poll_fit_job(jobs[[id]])
      if (!identical(j$status, "running")) {
        shiny::isolate(record_finished_job(state, j))
        if (identical(j$status, "completed")) {
          shiny::showNotification(paste("Fit", j$label, "completed."),
                                  type = "message")
        } else {
          shiny::showNotification(
            paste0("Fit ", j$label, " ", j$status, ": ",
                   j$result$error %||% ""), type = "error", duration = NULL)
        }
      }
    }
  })

  output$status_bar <- shiny::renderText({
    dm <- state$data_meta
    act <- state$fits[[state$active %||% ""]]
    n_run <- sum(vapply(state$jobs, function(j)
      identical(j$status, "running"), logical(1)))
    paste0("Data: ", if (is.null(dm)) "none" else dm$name,
           " | Active fit: ", if (is.null(act)) "none" else act$label,
           " | Running jobs: ", n_run)
  })

  output$draft_code <- shiny::renderText({
    spec <- state$spec
    if (is.null(spec)) return("Choose data and a model on the Build page.")
    probs <- validate_spec(spec, state$data)
    code <- analysis_script(spec)
    if (length(probs)) {
      paste0("# Unresolved problems:\n", paste0("#  - ", probs,
                                                collapse = "\n"),
             "\n\n", code)
    } else code
  })
  output$fitted_code <- shiny::renderText({
    b <- state$fits[[state$active %||% ""]]
    if (is.null(b)) "No completed fit is selected." else b$script
  })

  # The code drawer is collapsed at first; keep its (cheap) outputs current
  shiny::outputOptions(output, "draft_code", suspendWhenHidden = FALSE)
  shiny::outputOptions(output, "fitted_code", suspendWhenHidden = FALSE)

  session$onSessionEnded(function() {
    jobs <- shiny::isolate(state$jobs)
    for (j in jobs) cancel_fit_job(j)
    unlink(shiny::isolate(state$session_dir), recursive = TRUE)
  })
}
