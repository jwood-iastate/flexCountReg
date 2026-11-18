# App: background jobs, the session log and script, and its reproduction
# in a clean R session.

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

test_that("a whole session is logged and reproduced by its script", {
  skip_on_cran()
  state <- new_state()
  shiny::isolate({
    app_load_data(state, "bundled", "washington_roads")
    app_add_variable(state, list(type = "log", x = "AADT",
                                 name = "log aadt"))
    base <- list(source = "bundled", name = "washington_roads",
                 object = "data_1")
    rec <- list(response = "Total_crashes",
                predictors = c("log aadt", "speed50"), exposure = "Length")
    s1 <- app_spec(data = base, recipe = rec, model = list(family = "NB2"))
    s2 <- app_spec(data = base, recipe = c(rec, list(group = "ID")),
                   model = list(estimator = "renb"))
    j1 <- app_submit_fit(state, s1, "nb2")
    j2 <- app_submit_fit(state, s2)
    record_finished_job(state, wait_job(j1))
    record_finished_job(state, wait_job(j2))
    expect_identical(unname(vapply(state$fits, `[[`, "", "fit_name")),
                     c("fit_1", "fit_2"))
    id1 <- j1$id
    sc <- typical_row(s1, state$data)
    sc <- rbind(sc, sc)
    sc[["log aadt"]][2] <- 9
    r_root <- app_run_output(state, list(type = "rootogram"), id1)
    r_q <- app_run_output(state, list(type = "quantile"), id1)
    r_cure <- app_run_output(state, list(type = "cure", var = "log aadt"),
                             id1)
    r_pred <- app_run_output(state, list(type = "predict", scenarios = sc,
                                         interval = "prediction",
                                         level = 0.9), id1)
    r_tol <- app_run_output(state, list(type = "tolerance", scenarios = sc,
                                        content = 0.9, confidence = 0.9,
                                        side = "upper", scope = "pointwise",
                                        tol_method = "auto"), id1)
    r_cmp <- app_run_output(state, list(type = "compare",
                                        fits = c("fit_1", "fit_2"),
                                        labels = c("A", "B")))
    expect_identical(r_root$out, "rootogram_1")
    expect_identical(names(r_pred$value), c("log aadt", "speed50",
                                            "Length", "fit", "lwr", "upr"))
    expect_identical(r_cmp$value$fit, c("A", "B"))
    kinds <- vapply(state$log, `[[`, "", "kind")
    expect_identical(kinds, c("data", "transform", "fit", "fit",
                              "rootogram", "quantile", "cure", "predict",
                              "tolerance", "compare"))
    script <- app_script(state)
    env <- state$env
    local_vals <- list(
      ll = vapply(c("fit_1", "fit_2"), function(n)
        as.numeric(logLik(get(n, envir = env))), 0),
      pred = env$pred_1, tol = as.data.frame(env$tol_1),
      cmp = env$comparison_1, root = as.data.frame(unclass(env$rootogram_1)),
      q = env$qresid_1)
  })
  f <- tempfile(fileext = ".R")
  writeLines(script, f)
  out <- callr::r(function(f) {
    env <- new.env()
    grDevices::pdf(NULL)
    sys.source(f, envir = env)
    list(ll = vapply(c("fit_1", "fit_2"), function(n)
      as.numeric(stats::logLik(env[[n]])), 0),
      pred = env$pred_1, tol = as.data.frame(env$tol_1),
      cmp = env$comparison_1, root = as.data.frame(unclass(env$rootogram_1)),
      q = env$qresid_1)
  }, args = list(f), libpath = .libPaths())
  expect_equal(out$ll, local_vals$ll, tolerance = 1e-10)
  expect_equal(out$pred, local_vals$pred, tolerance = 1e-10)
  expect_equal(out$tol$upper, local_vals$tol$upper)
  expect_equal(out$cmp, local_vals$cmp, tolerance = 1e-10)
  expect_equal(out$root$expected, local_vals$root$expected,
               tolerance = 1e-10)
  expect_equal(out$q, local_vals$q, tolerance = 1e-10)
})

test_that("uploaded data: the script names the file, the bundle can hold it", {
  skip_on_cran()
  state <- new_state()
  set.seed(3)
  d <- data.frame(`count y` = c(rpois(150, 1), rpois(100, 6)),
                  `x 1` = rnorm(250), check.names = FALSE)
  up <- file.path(shiny::isolate(state$session_dir), "u.csv")
  utils::write.csv(d, up, row.names = FALSE)
  shiny::isolate({
    app_load_data(state, "upload", "my data.csv", upload_path = up)
    s <- app_spec(data = list(source = "upload", name = "my data.csv",
                              object = "data_1"),
                  recipe = list(response = "count y", predictors = "x 1"),
                  model = list(estimator = "countreg_mix",
                               family = "POISSON", components = 2L,
                               n_starts = 2L, seed = 4L))
    record_finished_job(state, wait_job(app_submit_fit(state, s)))
    ll <- as.numeric(logLik(state$env$fit_1))
    dir <- tempfile("bundle-")
    dir.create(dir)
    files <- write_bundle(state, dir, include_data = TRUE)
    dir2 <- tempfile("bundle-")
    dir.create(dir2)
    files2 <- write_bundle(state, dir2, include_data = FALSE)
  })
  expect_setequal(files, c("analysis.R", "settings.json", "my data.csv",
                           "manifest.txt"))
  expect_false("my data.csv" %in% files2)
  expect_match(readLines(file.path(dir2, "analysis.R")),
               "data_1_file <- \"my data.csv\"", all = FALSE, fixed = TRUE)
  out <- callr::r(function(dir) {
    setwd(dir)
    env <- new.env()
    sys.source("analysis.R", envir = env)
    as.numeric(stats::logLik(env$fit_1))
  }, args = list(dir), libpath = .libPaths())
  expect_equal(out, ll, tolerance = 1e-10)
  # The settings file rebuilds the same fit
  sj <- jsonlite::fromJSON(file.path(dir, "settings.json"),
                           simplifyVector = TRUE)
  s2 <- .spec_from_list(sj$fit_1)
  expect_identical(compile_fit(s2, "fit_1"), compile_fit(s, "fit_1"))
})

test_that("failures and cancellations are reported, not hidden", {
  skip_on_cran()
  state <- new_state()
  shiny::isolate({
    app_load_data(state, "bundled", "abalone_rings")
    bad <- app_spec(recipe = list(response = "rings"),
                    model = list(family = "NOT_A_FAMILY"))
    expect_error(app_submit_fit(state, bad), "Unknown family")
    expect_length(state$jobs, 0)
    # A fit that fails in the worker
    s_err <- app_spec(recipe = list(response = "rings",
                                    predictors = "not_a_column"),
                      model = list(family = "NB2"))
    jerr <- start_fit_job(s_err, "fit_y", list(data_1 = state$data),
                          dir = state$session_dir)
    jerr <- wait_job(jerr)
    expect_identical(jerr$status, "failed")
    expect_match(jerr$result$error, "not_a_column")
    record_finished_job(state, jerr)
    expect_length(state$fits, 0)
    # Cancellation of a long job
    set.seed(1)
    big <- data.frame(y = rpois(20000, 3), x = rnorm(20000))
    s2 <- app_spec(recipe = list(response = "y", predictors = "x"),
                   model = list(estimator = "countreg_mix", family = "COM",
                                components = 3L, n_starts = 20L))
    j2 <- start_fit_job(s2, "fit_z", list(data_1 = big),
                        dir = state$session_dir)
    Sys.sleep(1)
    j2 <- cancel_fit_job(j2)
    expect_identical(j2$status, "canceled")
    Sys.sleep(0.5)
    expect_false(j2$proc$is_alive())
  })
})

test_that("a finished older job does not replace a newer active fit", {
  skip_on_cran()
  state <- new_state()
  shiny::isolate({
    app_load_data(state, "bundled", "washington_roads")
    s <- app_spec(data = list(source = "bundled", name = "washington_roads",
                              object = "data_1"),
                  recipe = list(response = "Total_crashes",
                                predictors = "lnaadt"),
                  model = list(family = "POISSON"))
    j_old <- app_submit_fit(state, s, "old")
    j_new <- app_submit_fit(state, s, "new")
    # The newer job finishes first, the older one afterwards
    record_finished_job(state, wait_job(j_new))
    record_finished_job(state, wait_job(j_old))
    expect_identical(state$active, j_new$id)
    expect_setequal(names(state$fits), c(j_old$id, j_new$id))
    expect_identical(state$fits[[j_old$id]]$label, "old")
    expect_identical(state$fits[[j_old$id]]$fit_name, "fit_1")
  })
})

test_that("operations report what a model does not support", {
  skip_on_cran()
  d <- load_bundled_data("washington_roads")
  env <- new.env()
  env$data_1 <- d
  s <- app_spec(recipe = list(response = "Total_crashes",
                              predictors = "lnaadt", group = "ID"),
                model = list(estimator = "renb"))
  run_compiled(compile_fit(s, "fit_1"), env)
  f <- env$fit_1
  expect_false(app_supports(f, "rootogram")$ok)
  expect_match(app_supports(f, "tolerance")$reason, "countreg")
  expect_true(app_supports(f, "cure")$ok)
})
