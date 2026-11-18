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
    transforms = list(),    # variables created, by data object
    newdata = character(),  # new-data objects for predictions
    outputs = list(),       # results of background output jobs, by id
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

# Load an uploaded CSV file of new data (for predictions) as object
# newdata_<k>; the current analysis data do not change
app_load_newdata <- function(state, name, upload_path) {
  k <- .app_next(state, "newdata")
  obj <- paste0("newdata_", k)
  code <- compile_data_load("upload", name, obj)
  env <- state$env
  assign(paste0(obj, "_file"), upload_path, envir = env)
  up <- state$uploads
  up[[obj]] <- list(path = upload_path, name = name)
  state$uploads <- up
  run_compiled(code$run, env)
  state$newdata <- c(state$newdata, stats::setNames(obj, paste0(
    obj, ": ", name, " (", nrow(get(obj, envir = env)), " rows)")))
  app_log(state, paste0("New data ", k, ": ", name), code$display, "data",
          obj)
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
    if (identical(j$kind, "output")) {
      app_log(state, j$label, j$exprs, j$type, j$out)
      ol <- state$outputs
      ol[[j$id]] <- list(out = j$out, exprs = j$exprs, type = j$type,
                         value = get(j$out, envir = env),
                         fit = state$fits[[j$fit_id]],
                         messages = j$result$messages,
                         warnings = j$result$warnings)
      state$outputs <- ol
      return(invisible(j$status))
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
  o <- app_prepare_output(state, req, fit_id)
  env <- state$env
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  run_compiled(o$exprs, env)
  app_log(state, o$title, o$exprs, req$type, o$out)
  list(out = o$out, exprs = o$exprs, value = get(o$out, envir = env),
       fit = o$fit, type = req$type)
}

# Compile an output request: its object name, code and log title. Problems
# that can be seen before running (missing columns of new data) are errors
# here.
app_prepare_output <- function(state, req, fit_id = NULL) {
  b <- if (!is.null(fit_id)) state$fits[[fit_id]]
  fit_name <- if (!is.null(b)) b$fit_name else "fit"
  if (!is.null(req$newdata)) {
    # New data: the columns the fit needs must be there (or be created by
    # the fit's variable creations)
    tr <- state$transforms[[b$spec$data$object]] %||% list()
    nd <- get(req$newdata, envir = state$env)
    miss <- setdiff(raw_columns_needed(b$spec, tr), names(nd))
    if (length(miss)) {
      stop("The new data lack column(s) the model needs: ", toString(miss),
           ".", call. = FALSE)
    }
  }
  n <- .app_next(state, req$type)
  out <- paste0(.app_out_stem[[req$type]], "_", n)
  req$title <- req$title %||% (if (!is.null(b)) b$label)
  if (!is.null(req$newdata)) {
    req$prep <- compile_prediction_data(b$spec, fit_name, tr, req$newdata,
                                        paste0(out, "_data"))
  }
  exprs <- compile_output(req, fit_name, out)
  title <- switch(req$type,
                  rootogram = "Rootogram", quantile = "Quantile residuals",
                  cure = "CURE plot", predict = "Predictions",
                  tolerance = "Tolerance intervals",
                  compare = "Model comparison")
  list(out = out, exprs = exprs, fit = b, type = req$type,
       title = paste0(title, " ", n, if (!is.null(b))
         paste0(" (", b$label, ")")))
}

# Will an interval of this fit need the parametric bootstrap (refits)
# because its covariance matrix is not positive definite?
app_needs_bootstrap <- function(fit) {
  V <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  !is.matrix(V) || !.fcr_is_pd(V)
}

# Run an output request in a background R process (predictions and
# intervals can take minutes, e.g. a parametric bootstrap); returns the
# job. The objects its code uses are sent to the worker; the objects it
# creates come back when it finishes (record_finished_job()).
app_start_output_job <- function(state, req, fit_id) {
  o <- app_prepare_output(state, req, fit_id)
  env <- state$env
  used <- unique(unlist(lapply(o$exprs, all.vars)))
  objects <- mget(intersect(used, ls(env)), envir = env)
  job <- start_job(o$exprs, objects, dir = state$session_dir,
                   label = o$title)
  job$kind <- "output"
  job$out <- o$out
  job$type <- o$type
  job$fit_id <- fit_id
  jl <- state$jobs
  jl[[job$id]] <- job
  state$jobs <- jl
  job
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

# Keep a page's model selector (input `id`) in step with the completed fits:
# it lists every fit and moves to each newly completed fit; in between, it
# stays where the user put it.
app_fit_selector <- function(input, session, state, id) {
  shiny::observe({
    ch <- label_choices(state$fits)
    cur <- shiny::isolate(input[[id]])
    sel <- if (!is.null(cur) && cur %in% ch) cur else utils::tail(ch, 1)
    shiny::updateSelectInput(session, id, choices = ch,
                             selected = if (length(ch)) sel)
  })
  shiny::observeEvent(state$active, {
    shiny::updateSelectInput(session, id, selected = state$active)
  })
}

# Code that writes the plot of an output (its plotting expressions) to a
# JPEG file
compile_plot_file <- function(exprs, filename, width, height, res) {
  c(list(call("<-", .sym("plot_file"), filename),
         bquote(grDevices::jpeg(plot_file, width = .(width),
                                height = .(height), units = "in",
                                res = .(res), quality = 95))),
    exprs, list(quote(grDevices::dev.off())))
}

# Save the plot of a diagnostic result `r` (from app_run_output()) as a
# JPEG file at `path`, and log the code (with the file name <out>.jpg)
app_save_plot <- function(state, r, path, width = 7, height = 5, res = 300) {
  ok <- function(v, lo, hi) is.numeric(v) && length(v) == 1 && !is.na(v) &&
    v >= lo && v <= hi
  if (!ok(width, 1, 30) || !ok(height, 1, 30)) {
    stop("Width and height must be between 1 and 30 inches.", call. = FALSE)
  }
  if (!ok(res, 50, 2400)) {
    stop("The resolution must be between 50 and 2400 dpi.", call. = FALSE)
  }
  fname <- paste0(r$out, ".jpg")
  # Draw in a child of the session environment (objects stay unchanged)
  env <- new.env(parent = state$env)
  run_compiled(compile_plot_file(r$exprs, path, width, height, res), env)
  app_log(state, paste0("Plot file ", fname, " (", width, " x ", height,
                        " in, ", res, " dpi)"),
          compile_plot_file(r$exprs, fname, width, height, res),
          "plot_file", fname)
  invisible(path)
}
