# Background estimation. Each fit runs in its own R process (callr::r_bg),
# from an immutable request saved in a session-specific directory. The
# worker loads the installed flexCountReg, evaluates the compiled
# expressions and saves the fit bundle; the app polls the process, never
# blocking the session. Results carry their request id, so a slow earlier
# job cannot overwrite a newer one.

.new_request_id <- function() {
  paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "-",
         paste(sample(c(letters, 0:9), 6, replace = TRUE), collapse = ""))
}

# The worker (runs in the child process; uses only base R and
# flexCountReg, so it does not depend on the app's global state)
.fit_worker <- function(request_path, result_path) {
  req <- readRDS(request_path)
  warnings <- character(0)
  messages <- character(0)
  t0 <- proc.time()[["elapsed"]]
  suppressPackageStartupMessages(requireNamespace("flexCountReg"))
  env <- new.env(parent = globalenv())
  if (!is.null(req$data_file)) assign("data_file", req$data_file, envir = env)
  status <- "completed"
  err <- NULL
  set.seed(if (is.null(req$seed)) 1L else req$seed)
  withCallingHandlers(
    tryCatch(for (e in req$exprs) eval(e, env),
             error = function(e) {
               status <<- "failed"
               err <<- conditionMessage(e)
             }),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      messages <<- c(messages, trimws(conditionMessage(m)))
      invokeRestart("muffleMessage")
    })
  fit <- if (exists("fit", envir = env, inherits = FALSE)) env$fit
  # The data are already kept in the app; drop them from the stored fit
  # only when asked (they are needed for fitted values and diagnostics)
  saveRDS(list(id = req$id, status = status, error = err, fit = fit,
               warnings = unique(warnings), messages = unique(messages),
               elapsed = proc.time()[["elapsed"]] - t0,
               versions = c(R = paste(R.version$major, R.version$minor,
                                      sep = "."),
                            flexCountReg = as.character(
                              utils::packageVersion("flexCountReg")))),
          result_path)
  invisible(status)
}

# Start a background fit
#
start_fit_job <- function(spec, data = NULL, dir = tempfile("fcr-app-"),
                          label = NULL) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  id <- .new_request_id()
  jdir <- file.path(dir, id)
  dir.create(jdir)
  data_file <- NULL
  if (!identical(spec$data$source, "bundled")) {
    data_file <- file.path(jdir, "data.csv")
    utils::write.csv(data, data_file, row.names = FALSE)
  }
  a <- compile_analysis(spec)
  req <- list(id = id, spec = spec, exprs = c(a$prep, a$fit),
              data_file = data_file, seed = spec$model$seed,
              label = label %||% id, created = Sys.time())
  req_path <- file.path(jdir, "request.rds")
  res_path <- file.path(jdir, "result.rds")
  saveRDS(req, req_path)
  proc <- callr::r_bg(.fit_worker, args = list(req_path, res_path),
                      libpath = .libPaths(), supervise = TRUE,
                      stdout = file.path(jdir, "stdout.txt"),
                      stderr = file.path(jdir, "stderr.txt"))
  list(id = id, label = req$label, spec = spec, proc = proc, dir = jdir,
       result_path = res_path, started = Sys.time(), status = "running")
}

# Status of a background fit
poll_fit_job <- function(job) {
  if (!identical(job$status, "running")) return(job)
  if (job$proc$is_alive()) return(job)
  if (file.exists(job$result_path)) {
    job$result <- readRDS(job$result_path)
    job$status <- job$result$status
  } else {
    err <- tryCatch(paste(readLines(file.path(job$dir, "stderr.txt"),
                                    warn = FALSE), collapse = "\n"),
                    error = function(e) "")
    job$result <- list(id = job$id, status = "failed",
                       error = paste("The worker stopped without a result.",
                                     err))
    job$status <- "failed"
  }
  job$finished <- Sys.time()
  job
}

cancel_fit_job <- function(job) {
  if (identical(job$status, "running")) {
    try(job$proc$kill_tree(), silent = TRUE)
    try(job$proc$kill(), silent = TRUE)
    job$status <- "canceled"
    job$finished <- Sys.time()
  }
  job
}

# A fit bundle: the immutable record of one completed fit
make_fit_bundle <- function(job, data_name) {
  res <- job$result
  list(id = job$id, label = job$label, spec = job$spec,
       fit = res$fit, warnings = res$warnings, messages = res$messages,
       elapsed = res$elapsed, versions = res$versions,
       data_name = data_name, created = job$finished %||% Sys.time(),
       script = analysis_script(job$spec, fit_label = job$label,
                                fit_id = job$id))
}
