# App: fixes from the external review of the app. Uploads with the same
# file name stay distinct in the export; warnings of predictions and
# diagnostics are kept, shown and exported; intervals that may use the
# parametric bootstrap carry a frozen seed; standard errors fall back to the
# covariance matrix correctly.

skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("DT")
skip_if_not_installed("callr")
skip_if_not_installed("jsonlite")
skip_if_not_installed("zip")

wait_job <- function(job, timeout = 180) {
  t0 <- Sys.time()
  while (identical((job <- poll_fit_job(job))$status, "running")) {
    if (difftime(Sys.time(), t0, units = "secs") > timeout) {
      stop("job timed out")
    }
    Sys.sleep(0.2)
  }
  job
}

new_state <- function() {
  st <- new_app_state()
  dir.create(shiny::isolate(st$session_dir))
  st
}

test_that("uploads with the same file name export to distinct files", {
  skip_on_cran()
  state <- new_state()
  sdir <- shiny::isolate(state$session_dir)
  set.seed(5)
  train <- data.frame(y = rpois(120, 2), x = rnorm(120))
  nd <- data.frame(x = c(-1, 0, 1.5))
  p1 <- file.path(sdir, "upload-1.csv")
  p2 <- file.path(sdir, "newdata-1.csv")
  utils::write.csv(train, p1, row.names = FALSE)
  utils::write.csv(nd, p2, row.names = FALSE)
  shiny::isolate({
    # Training data and new data both uploaded as "data.csv"
    app_load_data(state, "upload", "data.csv", upload_path = p1)
    s <- app_spec(data = list(source = "upload", name = "data.csv",
                              object = "data_1"),
                  recipe = list(response = "y", predictors = "x"),
                  model = list(family = "POISSON"))
    job <- app_submit_fit(state, s)
    record_finished_job(state, wait_job(job))
    obj <- app_load_newdata(state, "data.csv", p2)
    r <- app_run_output(state, list(type = "predict", newdata = obj,
                                    interval = "none"), job$id)
    ll <- as.numeric(logLik(state$env$fit_1))
    dir <- tempfile("bundle-")
    dir.create(dir)
    files <- write_bundle(state, dir, include_data = TRUE)
    script <- app_script(state)
  })
  expect_true(all(c("data/data_1.csv", "data/newdata_1.csv") %in% files))
  expect_match(script, "data_1_file <- \"data/data_1.csv\"", fixed = TRUE)
  expect_match(script, "newdata_1_file <- \"data/newdata_1.csv\"",
               fixed = TRUE)
  # Each copy is the file that was loaded, and the manifest says so
  md5 <- unname(tools::md5sum(c(p1, p2)))
  expect_identical(unname(tools::md5sum(file.path(
    dir, c("data/data_1.csv", "data/newdata_1.csv")))), md5)
  manifest <- readLines(file.path(dir, "manifest.txt"))
  expect_match(manifest, paste0("data/data_1.csv <- uploaded as \"data.csv\"",
                                " (object data_1; MD5 ", md5[1], ")"),
               fixed = TRUE, all = FALSE)
  expect_match(manifest, paste0("data/newdata_1.csv <- uploaded as ",
                                "\"data.csv\" (object newdata_1; MD5 ",
                                md5[2], ")"),
               fixed = TRUE, all = FALSE)
  # The zipped bundle replays both datasets in a clean R process
  zf <- tempfile(fileext = ".zip")
  zip::zip(zf, files, root = dir)
  dir2 <- tempfile("unzipped-")
  utils::unzip(zf, exdir = dir2)
  out <- callr::r(function(dir) {
    setwd(dir)
    env <- new.env()
    sys.source("analysis.R", envir = env)
    list(ll = as.numeric(stats::logLik(env$fit_1)), n = nrow(env$data_1),
         nd = env$newdata_1, pred = env$pred_1$expected_count)
  }, args = list(dir2), libpath = .libPaths())
  expect_equal(out$ll, ll, tolerance = 1e-10)
  expect_identical(out$n, 120L)
  expect_equal(out$nd, nd)
  expect_equal(out$pred, r$value$expected_count, tolerance = 1e-10)
  # A copy that cannot be made, or differs from what was loaded, is an error
  writeLines("x\n99", p2)
  dir3 <- tempfile("bundle-")
  dir.create(dir3)
  expect_error(shiny::isolate(write_bundle(state, dir3, include_data = TRUE)),
               "does not match")
  unlink(p1)
  dir4 <- tempfile("bundle-")
  dir.create(dir4)
  expect_error(shiny::isolate(write_bundle(state, dir4, include_data = TRUE)),
               "Could not copy")
  # Without the data the bundle can still be written
  dir5 <- tempfile("bundle-")
  dir.create(dir5)
  expect_false(any(grepl("^data/", shiny::isolate(
    write_bundle(state, dir5, include_data = FALSE)))))
})

test_that("warnings of a successful output are kept, shown and exported", {
  skip_on_cran()
  state <- new_state()
  set.seed(2)
  n <- 300
  d <- data.frame(x = runif(n), z = rnorm(n))
  d$y <- rpois(n, exp(0.5 + 0.4 * d$x + 0.3 * d$z))
  path <- file.path(shiny::isolate(state$session_dir), "upload-1.csv")
  utils::write.csv(d, path, row.names = FALSE)
  # The exact prediction at x = 1e5 overflows: a warning on a result that
  # is otherwise returned
  sc <- data.frame(z = c(0, 0), x = c(1, 1e5))
  req <- list(type = "predict", scenarios = sc, interval = "none",
              method = "Exact")
  shiny::isolate({
    app_load_data(state, "upload", "rp.csv", upload_path = path)
    s <- app_spec(data = list(source = "upload", name = "rp.csv",
                              object = "data_1"),
                  recipe = list(response = "y", predictors = c("z", "x")),
                  model = list(estimator = "countreg_rp", family = "POISSON",
                               rpar = "x", ndraws = 50L))
    fj <- app_submit_fit(state, s)
    record_finished_job(state, wait_job(fj))
    # In the session (as the Diagnose and Compare pages run their outputs)
    r <- app_run_output(state, req, fj$id)
    # In a background job (as the Predict page runs them)
    j <- wait_job(app_start_output_job(state, req, fj$id))
    expect_identical(j$status, "completed")
    record_finished_job(state, j)
    res <- state$outputs[[j$id]]
    log <- state$log
  })
  expect_match(r$warnings, "overflow", all = FALSE)
  expect_true(is.infinite(r$value$expected_count[2]))
  expect_match(res$warnings, "overflow", all = FALSE)
  expect_match(log[[length(log)]]$warnings, "overflow", all = FALSE)
  expect_match(log[[length(log) - 1]]$warnings, "overflow", all = FALSE)
  # Shown beside the result
  html <- as.character(output_notes(res))
  expect_match(html, "alert-warning", fixed = TRUE)
  expect_match(html, "overflow", fixed = TRUE)
  # In the manifest of the export, under each output
  dir <- tempfile("bundle-")
  dir.create(dir)
  shiny::isolate(write_bundle(state, dir))
  manifest <- readLines(file.path(dir, "manifest.txt"))
  expect_true("Outputs:" %in% manifest)
  expect_match(manifest, "^  pred_2: Predictions 2", all = FALSE)
  expect_equal(sum(grepl("^    warning: .*overflow", manifest)), 2)

  # The Predict page shows the warnings of the result it picks up
  shiny::testServer(mod_predict_server, args = list(state = state), {
    session$setInputs(fit = fj$id, rows = "table", interval = "none",
                      level = 0.95, method = "Exact")
    # x is the second column of the scenario table
    session$setInputs(scen_cell_edit = list(row = 1, col = 1,
                                            value = "1e5"))
    expect_equal(state$scenarios$x, 1e5)
    session$setInputs(go = 1)
    ids <- names(state$jobs)
    job <- wait_job(state$jobs[[ids[length(ids)]]])
    record_finished_job(state, job)
    session$flushReact()
    expect_match(output$res_title, "pred_3", fixed = TRUE)
    expect_match(output$res_warnings$html, "overflow", fixed = TRUE)
  })
})

test_that("a bootstrap interval replays exactly from its frozen seed", {
  skip_on_cran()
  state <- new_state()
  set.seed(11)
  d <- data.frame(x = rnorm(150))
  d$y <- rpois(150, exp(0.3 + 0.5 * d$x))
  path <- file.path(shiny::isolate(state$session_dir), "upload-1.csv")
  utils::write.csv(d, path, row.names = FALSE)
  sc <- data.frame(x = c(-1, 1))
  req <- list(type = "predict", scenarios = sc, interval = "confidence",
              level = 0.9)
  shiny::isolate({
    app_load_data(state, "upload", "boot.csv", upload_path = path)
    s <- app_spec(data = list(source = "upload", name = "boot.csv",
                              object = "data_1"),
                  recipe = list(response = "y", predictors = "x"),
                  model = list(family = "POISSON"))
    fj <- app_submit_fit(state, s)
    record_finished_job(state, wait_job(fj))
    # A covariance matrix that is not positive definite: the interval falls
    # back to the parametric bootstrap
    state$env$fit_1$model$vcov[] <- 0
    expect_true(app_needs_bootstrap(state$env$fit_1))
    fit <- state$env$fit_1
    set.seed(1)
    r <- app_run_output(state, req, fj$id)
    # The same request in a background job gets its own frozen seed
    j <- app_start_output_job(state, req, fj$id)
    j <- wait_job(j)
    record_finished_job(state, j)
    res <- state$outputs[[j$id]]
    dir <- tempfile("bundle-")
    dir.create(dir)
    write_bundle(state, dir)
  })
  expect_true(is.integer(r$seed) && length(r$seed) == 1)
  expect_true(is.integer(res$seed))
  # The seed is in the code that ran, which is the code in the script
  code <- paste(vapply(r$exprs, .deparse, character(1)), collapse = "\n")
  expect_match(code, paste0("seed = ", r$seed, "L"), fixed = TRUE)
  expect_match(paste(vapply(res$exprs, .deparse, character(1)),
                     collapse = "\n"),
               paste0("seed = ", res$seed, "L"), fixed = TRUE)
  # Uncertainty method, refits and seed beside the result and in the
  # manifest
  expect_match(r$details, "parametric bootstrap", all = FALSE)
  expect_match(r$details, "^Bootstrap refits accepted: [0-9]+ of 200$",
               all = FALSE)
  expect_true(paste0("Seed: ", r$seed) %in% r$details)
  expect_match(r$messages, "parametric bootstrap", all = FALSE)
  manifest <- readLines(file.path(dir, "manifest.txt"))
  expect_true(paste0("    Seed: ", r$seed) %in% manifest)
  expect_true(paste0("    Seed: ", res$seed) %in% manifest)
  expect_true("    Uncertainty method: parametric bootstrap" %in% manifest)
  # The worker used the seed of its request: the same limits as running
  # that code here
  env <- new.env()
  env$fit_1 <- fit
  set.seed(123)
  run_compiled(res$exprs, env)
  expect_equal(env$pred_2, res$value)
  # Exact replay in a clean R process that starts from another RNG state
  out <- callr::r(function(fit, exprs) {
    # As in an exported script, the package's methods must be registered
    loadNamespace("flexCountReg")
    set.seed(999)
    env <- new.env()
    env$fit_1 <- fit
    for (e in exprs) eval(e, env)
    env$pred_1
  }, args = list(fit, r$exprs), libpath = .libPaths())
  expect_identical(out, r$value)
  # Tolerance intervals carry their seed too
  tol <- compile_output(list(type = "tolerance", scenarios = sc,
                             content = 0.9, confidence = 0.9,
                             side = "upper", scope = "pointwise",
                             tol_method = "joint", seed = 7L),
                        "fit_1", "tol_1")
  expect_match(.deparse(tol[[length(tol)]]), "seed = 7L", fixed = TRUE)
  # Requests that draw no random numbers get no seed
  shiny::isolate({
    o <- app_prepare_output(state, list(type = "predict", scenarios = sc,
                                        interval = "none"), fj$id)
  })
  expect_null(o$seed)
})

test_that("unavailable standard errors come from the covariance matrix", {
  skip_on_cran()
  d <- load_bundled_data("washington_roads")
  f <- countreg(Total_crashes ~ lnaadt, data = d, family = "NB2")
  unavailable <- function(f) {
    tab <- fit_checks(f)
    tab$value[tab$check == "Unavailable standard errors"]
  }
  # An otherwise valid object without $se: every standard error available
  f$model$se <- NULL
  expect_identical(unavailable(f), "none")
  # A negative variance is unavailable, the others are not
  f$model$vcov[2, 2] <- -1
  expect_identical(unavailable(f), names(f$model$estimate)[2])
  # No covariance matrix: all are unavailable, without an error
  f$model$vcov <- NULL
  expect_identical(unavailable(f), toString(names(f$model$estimate)))
})
