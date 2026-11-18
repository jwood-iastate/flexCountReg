# Background estimation. Each fit runs in its own R process (callr::r_bg)
# from an immutable request saved in the session directory: the data object
# of the session and the compiled fit expressions. The worker evaluates the
# expressions and saves the fitted objects; the app polls the process and
# never blocks. Results carry their request id, so a slow earlier job never
# replaces a newer one.

.new_request_id <- function() {
  paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "-",
         paste(sample(c(letters, 0:9), 6, replace = TRUE), collapse = ""))
}

# The worker (runs in the child process; base R and flexCountReg only, so it
# does not depend on the app's state)
.fit_worker <- function(request_path, result_path) {
  req <- readRDS(request_path)
  warnings <- character(0)
  messages <- character(0)
  t0 <- proc.time()[["elapsed"]]
  suppressPackageStartupMessages(requireNamespace("flexCountReg"))
  env <- new.env(parent = globalenv())
  for (nm in names(req$objects)) assign(nm, req$objects[[nm]], envir = env)
  status <- "completed"
  err <- NULL
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
  created <- setdiff(ls(env), names(req$objects))
  saveRDS(list(id = req$id, status = status, error = err,
               objects = mget(created, envir = env),
               warnings = unique(warnings), messages = unique(messages),
               elapsed = proc.time()[["elapsed"]] - t0,
               versions = c(R = paste(R.version$major, R.version$minor,
                                      sep = "."),
                            flexCountReg = as.character(
                              utils::packageVersion("flexCountReg")))),
          result_path)
  invisible(status)
}

# Start a background fit. `objects` is a named list of the objects the
# expressions need (the session's data object).
start_fit_job <- function(spec, fit_name, objects, dir = tempfile("fcr-app-"),
                          label = NULL) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  id <- .new_request_id()
  jdir <- file.path(dir, id)
  dir.create(jdir)
  exprs <- compile_fit(spec, fit_name)
  req <- list(id = id, spec = spec, fit_name = fit_name, exprs = exprs,
              objects = objects, label = label %||% fit_name,
              created = Sys.time())
  req_path <- file.path(jdir, "request.rds")
  res_path <- file.path(jdir, "result.rds")
  saveRDS(req, req_path)
  proc <- callr::r_bg(.fit_worker, args = list(req_path, res_path),
                      libpath = .libPaths(), supervise = TRUE,
                      stdout = file.path(jdir, "stdout.txt"),
                      stderr = file.path(jdir, "stderr.txt"))
  list(id = id, label = req$label, spec = spec, fit_name = fit_name,
       exprs = exprs, proc = proc, dir = jdir, result_path = res_path,
       started = Sys.time(), status = "running")
}

# Update a job's status ("running", "completed", "failed" or "canceled")
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
