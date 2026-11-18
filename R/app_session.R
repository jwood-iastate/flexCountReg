# Session state and the operations the pages perform. Everything the app
# runs is compiled to R code (app_compile.R), evaluated in the session
# environment `state$env` (where the objects of the session script live:
# data_1, fit_1, rootogram_1, ...) and appended to `state$log`, from which
# the session script is written. All state belongs to one browser session.

new_app_state <- function(session_dir = tempfile("fcr-app-")) {
  shiny::reactiveValues(
    env = new.env(parent = globalenv()),  # objects of the session script
    log = list(),           # what ran, in order (the session script)
    counters = list(),      # numbering of data, fits and outputs
    data = NULL,            # the current data frame
    data_meta = NULL,       # list(source, name, title, object, ...)
    uploads = list(),       # uploaded files, by data object
    spec = NULL,            # draft specification (Build page)
    jobs = list(),          # running and finished jobs, by request id
    fits = list(),          # completed fits, by request id
    active = NULL,          # id of the selected fit
    last_submitted = NULL,  # id of the most recent submission
    scenarios = NULL,       # scenario table of the Predict page
    session_dir = session_dir)
}

# Next number for a kind of object (data, fit, rootogram, ...)
.app_next <- function(state, kind) {
  cn <- state$counters
  n <- (cn[[kind]] %||% 0L) + 1L
  cn[[kind]] <- n
  state$counters <- cn
  n
}

# Append an entry to the session log
app_log <- function(state, title, display, kind, ref = NULL) {
  state$log <- c(state$log, list(list(title = title, display = display,
                                      kind = kind, ref = ref,
                                      time = Sys.time())))
  invisible(NULL)
}

# The session script so far
app_script <- function(state) {
  session_script(state$log)
}

# Load a dataset (bundled name, or an uploaded file with its local path)
app_load_data <- function(state, source, name, upload_path = NULL,
                          title = name, response = NULL, group = NA) {
  k <- .app_next(state, "data")
  obj <- paste0("data_", k)
  code <- compile_data_load(source, name, obj)
  env <- state$env
  if (!identical(source, "bundled")) {
    assign(paste0(obj, "_file"), upload_path, envir = env)
    up <- state$uploads
    up[[obj]] <- list(path = upload_path, name = name)
    state$uploads <- up
  }
  run_compiled(code$run, env)
  state$data <- env[[obj]]
  state$data_meta <- list(source = source, name = name, title = title,
                          object = obj, response = response, group = group)
  app_log(state, paste0("Data ", k, ": ", name), code$display, "data", obj)
  invisible(obj)
}

# Submit a fit of `spec` (its data object must be in the session); returns
# the job
app_submit_fit <- function(state, spec, label = NULL) {
  obj <- spec$data$object
  p <- validate_spec(spec, state$env[[obj]])
  if (length(p)) stop(paste(p, collapse = " "), call. = FALSE)
  j <- .app_next(state, "fit")
  fit_name <- paste0("fit_", j)
  label <- if (is.null(label) || !nzchar(label)) {
    paste0("Fit ", j, ": ", .app_estimator(spec$model$estimator)$label,
           if (.app_estimator(spec$model$estimator)$uses_family)
             paste0(", ", paste(unique(spec$model$family), collapse = "+")),
           if (.app_is_mixture(spec$model$estimator))
             paste0(", K = ", spec$model$components))
  } else label
  job <- start_fit_job(spec, fit_name,
                       objects = stats::setNames(list(state$env[[obj]]), obj),
                       dir = state$session_dir, label = label)
  jl <- state$jobs
  jl[[job$id]] <- job
  state$jobs <- jl
  state$last_submitted <- job$id
  job
}

# Record a finished job: keep it; for a completed fit, put its objects in
# the session environment, log its code, and make it the active fit only if
# it is the most recent submission (a slow earlier job never replaces a
# newer result).
record_finished_job <- function(state, j) {
  jl <- state$jobs
  jl[[j$id]] <- j
  state$jobs <- jl
  if (identical(j$status, "completed")) {
    env <- state$env
    for (nm in names(j$result$objects)) {
      assign(nm, j$result$objects[[nm]], envir = env)
    }
    res <- j$result
    fl <- state$fits
    fl[[j$id]] <- list(id = j$id, label = j$label, fit_name = j$fit_name,
                       spec = j$spec, warnings = res$warnings,
                       messages = res$messages, elapsed = res$elapsed,
                       versions = res$versions,
                       data_name = j$spec$data$name,
                       created = j$finished %||% Sys.time())
    state$fits <- fl
    app_log(state, j$label, j$exprs, "fit", j$id)
    if (identical(j$id, state$last_submitted)) state$active <- j$id
  }
  invisible(j$status)
}

# The fitted model object of a fit id
app_fit_object <- function(state, id) {
  b <- state$fits[[id]]
  if (is.null(b)) return(NULL)
  get(b$fit_name, envir = state$env)
}

# Short object-name stems of the outputs
.app_out_stem <- c(rootogram = "rootogram", quantile = "qresid",
                   cure = "cure", predict = "pred", tolerance = "tol",
                   compare = "comparison")

# Run an output request on a fit (or a comparison of fits): compile it,
# evaluate it in the session environment (plots go to a null device; the
# page redraws them), log it, and return its object name and code.
app_run_output <- function(state, req, fit_id = NULL) {
  b <- if (!is.null(fit_id)) state$fits[[fit_id]]
  fit_name <- if (!is.null(b)) b$fit_name else "fit"
  n <- .app_next(state, req$type)
  out <- paste0(.app_out_stem[[req$type]], "_", n)
  req$title <- req$title %||% (if (!is.null(b)) b$label)
  exprs <- compile_output(req, fit_name, out)
  env <- state$env
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  run_compiled(exprs, env)
  title <- switch(req$type,
                  rootogram = "Rootogram", quantile = "Quantile residuals",
                  cure = "CURE plot", predict = "Predictions",
                  tolerance = "Tolerance intervals",
                  compare = "Model comparison")
  app_log(state, paste0(title, " ", n, if (!is.null(b))
    paste0(" (", b$label, ")")), exprs, req$type, out)
  list(out = out, exprs = exprs, value = get(out, envir = env), fit = b)
}

# Run code from a page, turning errors into a visible notification instead
# of ending the session
app_try <- function(expr, what = "This action") {
  tryCatch(expr, error = function(e) {
    shiny::showNotification(paste0(what, " failed: ", conditionMessage(e)),
                            type = "error", duration = NULL)
    NULL
  })
}
