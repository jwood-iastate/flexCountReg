# Background jobs, stale results, cancellation and clean-session export.

skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("DT")
skip_if_not_installed("callr")
skip_if_not_installed("jsonlite")
skip_if_not_installed("zip")

wait_job <- function(job, timeout = 120) {
  t0 <- Sys.time()
  while (identical((job <- poll_fit_job(job))$status, "running")) {
    if (difftime(Sys.time(), t0, units = "secs") > timeout) {
      stop("job timed out")
    }
    Sys.sleep(0.2)
  }
  job
}

test_that("a background fit equals the same fit in this session", {
  skip_on_cran()
  s <- app_spec(data = list(source = "bundled", name = "washington_roads"),
                recipe = list(response = "Total_crashes",
                              predictors = c("lnaadt", "speed50"),
                              exposure = "Length"),
                model = list(family = "NB2"))
  job <- wait_job(start_fit_job(s, label = "nb2"))
  expect_identical(job$status, "completed")
  env <- run_compiled(c(compile_analysis(s)$prep, compile_analysis(s)$fit))
  expect_equal(as.numeric(logLik(job$result$fit)),
               as.numeric(logLik(env$fit)))
  b <- make_fit_bundle(job, "washington_roads")
  expect_identical(b$id, job$id)
  expect_match(b$script, job$id, fixed = TRUE)
})

test_that("failures and cancellations are reported, not hidden", {
  skip_on_cran()
  d <- data.frame(y = c(rpois(50, 2), NA), x = c(rnorm(50), NA))
  s <- app_spec(data = list(source = "upload", name = "d.csv"),
                recipe = list(response = "y", predictors = "x"),
                model = list(family = "NOT_A_FAMILY"))
  job <- wait_job(start_fit_job(s, data = d))
  expect_identical(job$status, "failed")
  expect_match(job$result$error, "Unknown")
  # Cancellation of a long job
  set.seed(1)
  big <- data.frame(y = rpois(20000, 3), x = rnorm(20000))
  s2 <- app_spec(data = list(source = "upload", name = "big.csv"),
                 recipe = list(response = "y", predictors = "x"),
                 model = list(estimator = "countreg_mix", family = "COM",
                              components = 3L, n_starts = 20L))
  job2 <- start_fit_job(s2, data = big)
  Sys.sleep(1)
  job2 <- cancel_fit_job(job2)
  expect_identical(job2$status, "canceled")
  Sys.sleep(0.5)
  expect_false(job2$proc$is_alive())
  expect_identical(poll_fit_job(job2)$status, "canceled")
})

test_that("the exported bundle reproduces the fit in a clean R session", {
  skip_on_cran()
  set.seed(3)
  d <- data.frame(`count y` = c(rpois(150, 1), rpois(100, 6)),
                  `x 1` = rnorm(250), grp = sample(c("a", "b"), 250, TRUE),
                  check.names = FALSE)
  s <- app_spec(data = list(source = "upload", name = "my data.csv"),
                recipe = list(response = "count y",
                              predictors = c("x 1", "grp"),
                              factors = "grp"),
                model = list(estimator = "countreg_mix", family = "POISSON",
                             components = 2L, n_starts = 3L, seed = 4L))
  job <- wait_job(start_fit_job(s, data = d, label = "mix"))
  expect_identical(job$status, "completed")
  b <- make_fit_bundle(job, "my data.csv")
  dir <- tempfile("bundle-")
  dir.create(dir)
  scen <- data.frame(`x 1` = c(0, 1), grp = factor(c("a", "b")),
                     check.names = FALSE)
  files <- write_bundle(b, dir, list(rootogram = TRUE, scenarios = scen),
                        data = d, include_data = TRUE)
  expect_setequal(files, c("analysis.R", "settings.json", "my data.csv",
                           "coefficients.csv", "manifest.txt"))
  out <- callr::r(function(dir) {
    setwd(dir)
    env <- new.env()
    sys.source("analysis.R", envir = env)
    list(ll = as.numeric(stats::logLik(env$fit)), coef = stats::coef(env$fit),
         pred = env$scenario_predictions$expected_count,
         rows = env$fit$model$rows_used)
  }, args = list(dir), libpath = .libPaths())
  f <- b$fit
  expect_equal(out$ll, as.numeric(logLik(f)), tolerance = 1e-10)
  expect_equal(out$coef, coef(f), tolerance = 1e-8)
  expect_equal(out$rows, f$model$rows_used)
  expect_equal(out$pred, predict(f, newdata = scen), tolerance = 1e-8)
  # The settings file rebuilds the same analysis
  s2 <- spec_from_json(path = file.path(dir, "settings.json"))
  expect_identical(compile_analysis(s2), compile_analysis(s))
  # Without the data, the script names the file it expects
  dir2 <- tempfile("bundle-")
  dir.create(dir2)
  files2 <- write_bundle(b, dir2, data = d, include_data = FALSE)
  expect_false("my data.csv" %in% files2)
  expect_match(readLines(file.path(dir2, "analysis.R")),
               "data_file <- \"my data.csv\"", all = FALSE, fixed = TRUE)
  expect_match(readLines(file.path(dir2, "manifest.txt")), "not included",
               all = FALSE)
})

test_that("comparisons state when fits are not comparable", {
  skip_on_cran()
  d <- load_bundled_data("washington_roads")
  mk <- function(fam, preds) {
    s <- app_spec(data = list(source = "bundled", name = "washington_roads"),
                  recipe = list(response = "Total_crashes",
                                predictors = preds),
                  model = list(family = fam))
    env <- run_compiled(c(compile_analysis(s)$prep, compile_analysis(s)$fit))
    list(id = fam, label = fam, spec = s, fit = env$fit,
         data_name = "washington_roads")
  }
  a <- mk("POISSON", "lnaadt")
  b <- mk("NB2", "lnaadt")
  cp <- compare_fits(list(a, b))
  expect_length(cp$not_comparable, 0)
  expect_equal(cp$table$AIC, c(AIC(a$fit), AIC(b$fit)))
  b2 <- b
  b2$spec$recipe$response <- "Animal"
  expect_match(compare_fits(list(a, b2))$not_comparable, "responses",
               all = FALSE)
})
