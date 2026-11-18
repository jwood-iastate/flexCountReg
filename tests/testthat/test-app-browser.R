# App: browser workflow in a real headless browser (shinytest2): data ->
# new variable -> two fits -> diagnostic -> predictions -> comparison ->
# download, and the downloaded script reproduces the session.

skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("DT")
skip_if_not_installed("callr")
skip_if_not_installed("jsonlite")
skip_if_not_installed("zip")

test_that("a user can run an analysis and download a working script", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if(is.null(tryCatch(chromote::find_chrome(), error = function(e) NULL)),
          "no Chrome")
  app <- shinytest2::AppDriver$new(app_object(), name = "workflow",
                                   load_timeout = 60000, timeout = 60000)
  on.exit(app$stop(), add = TRUE)
  status <- function() app$get_value(output = "status_bar")
  wait_for_fit <- function(label) {
    app$wait_for_js(paste0(
      "document.getElementById('status_bar').innerText",
      ".includes('Selected fit: ", label, "')"), timeout = 120000)
  }

  app$set_inputs(`data-source` = "bundled",
                 `data-bundled` = "washington_roads")
  app$click("data-use")
  app$wait_for_idle()
  expect_match(status(), "washington_roads \\(data_1\\)")

  # A new variable
  app$set_inputs(`data-tf_type` = "indicator", `data-tf_x` = "AADT")
  app$set_inputs(`data-tf_op` = ">", `data-tf_value` = "10000",
                 `data-tf_name` = "aadt10k")
  app$click("data-tf_add")
  app$wait_for_idle()

  # Two fits, one after the other
  app$set_inputs(page = "Build")
  app$set_inputs(`build-response` = "Total_crashes")
  app$set_inputs(`build-predictors` = c("lnaadt", "aadt10k"))
  app$set_inputs(`build-exposure` = "Length")
  app$set_inputs(`build-label` = "NB2 A")
  app$click("build-estimate")
  wait_for_fit("NB2 A")
  app$set_inputs(page = "Build")
  app$set_inputs(`build-family` = "POISSON", `build-label` = "Poisson B")
  app$click("build-estimate")
  wait_for_fit("Poisson B")
  app$wait_for_idle()
  expect_match(app$get_value(output = "fit-summary"), "Poisson B")

  # A diagnostic and predictions of the selected fit
  app$set_inputs(page = "Diagnose")
  app$set_inputs(`diagnose-diag` = "cure", wait_ = FALSE)
  app$click("diagnose-go")
  app$wait_for_idle()
  expect_match(app$get_value(output = "diagnose-title"), "cure_1")
  app$set_inputs(page = "Predict")
  app$wait_for_idle()
  app$set_inputs(`predict-interval` = "confidence", wait_ = FALSE)
  app$click("predict-go")
  app$wait_for_idle()
  expect_match(app$get_value(output = "predict-res_title"), "pred_1")
  app$click("predict-go_tol")
  app$wait_for_idle()
  expect_match(app$get_value(output = "predict-res_title"), "tol_1")

  # A comparison of the two fits
  app$set_inputs(page = "Compare")
  app$wait_for_idle()
  app$run_js(paste0(
    "document.querySelectorAll('#compare-fits input').forEach(",
    "e => { e.checked = true; $(e).trigger('change'); });"))
  app$wait_for_idle()
  app$click("compare-go")
  app$wait_for_idle()
  expect_match(app$get_value(output = "compare-title"), "comparison_1")

  # Everything that ran is in the session script
  script <- app$get_value(output = "session_code")
  for (pat in c("data_1 <- local", "data_1[[\"aadt10k\"]] <-",
                "fit_1 <- flexCountReg::countreg", "fit_2 <-",
                "cure_1 <-", "pred_1 <-", "tol_1 <-", "comparison_1 <-")) {
    expect_match(script, pat, fixed = TRUE)
  }

  # The download reproduces it without the app
  app$set_inputs(page = "Export")
  app$wait_for_idle()
  zipfile <- app$get_download("export-download")
  dir <- tempfile("dl-")
  utils::unzip(zipfile, exdir = dir)
  expect_true(all(c("analysis.R", "settings.json", "manifest.txt") %in%
                    list.files(dir)))
  out <- callr::r(function(dir) {
    setwd(dir)
    env <- new.env()
    grDevices::pdf(NULL)
    sys.source("analysis.R", envir = env)
    c(as.numeric(stats::logLik(env$fit_1)),
      as.numeric(stats::logLik(env$fit_2)), nrow(env$comparison_1))
  }, args = list(dir), libpath = .libPaths())
  d <- load_bundled_data("washington_roads")
  d$aadt10k <- as.integer(d$AADT > 10000)
  ref <- countreg(Total_crashes ~ lnaadt + aadt10k + offset(log(Length)),
                  data = d, family = "NB2", method = "BFGS")
  expect_equal(out[1], as.numeric(logLik(ref)), tolerance = 1e-8)
  expect_equal(out[3], 2)
})
