# Browser workflow: data -> model -> fit -> diagnostics -> download, in a
# real headless browser (shinytest2). Skipped where no browser is available.

skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("DT")
skip_if_not_installed("callr")
skip_if_not_installed("jsonlite")
skip_if_not_installed("zip")

test_that("a user can fit a model and download a working analysis", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if(is.null(tryCatch(chromote::find_chrome(), error = function(e) NULL)),
          "no Chrome")
  app <- shinytest2::AppDriver$new(app_object(), name = "workflow",
                                   load_timeout = 60000, timeout = 60000)
  on.exit(app$stop(), add = TRUE)

  app$set_inputs(`data-source` = "bundled",
                 `data-bundled` = "washington_roads")
  app$click("data-use")
  app$wait_for_idle()
  expect_match(app$get_value(output = "status_bar"), "washington_roads")

  app$set_inputs(page = "Build")
  app$set_inputs(`build-response` = "Total_crashes")
  app$set_inputs(`build-predictors` = c("lnaadt", "speed50"))
  app$set_inputs(`build-exposure` = "Length")
  app$set_inputs(`build-label` = "NB2 browser")
  app$wait_for_idle()
  draft <- app$get_value(output = "draft_code")
  expect_match(draft, "offset(log(Length))", fixed = TRUE)
  app$click("build-estimate")
  # The fit runs in a background process; wait until it is the active fit
  app$wait_for_js(paste0("document.getElementById('status_bar').innerText",
                         ".includes('Active fit: NB2 browser')"),
                  timeout = 60000)
  app$set_inputs(page = "Fit")
  app$wait_for_idle()
  expect_match(app$get_value(output = "fit-summary"), "NB2 browser")
  expect_match(app$get_value(output = "status_bar"), "NB2 browser")
  fitted <- app$get_value(output = "fitted_code")
  expect_match(fitted, "flexCountReg::countreg(", fixed = TRUE)

  app$set_inputs(page = "Diagnose")
  app$set_inputs(`diagnose-diag` = "rootogram", wait_ = FALSE)
  app$click("diagnose-go")
  app$wait_for_idle()
  expect_match(app$get_value(output = "diagnose-title"), "Rootogram")

  app$set_inputs(page = "Export")
  app$wait_for_idle()
  zipfile <- app$get_download("export-download")
  dir <- tempfile("dl-")
  utils::unzip(zipfile, exdir = dir)
  expect_true(all(c("analysis.R", "settings.json", "manifest.txt",
                    "coefficients.csv") %in% list.files(dir)))
  ll <- callr::r(function(dir) {
    setwd(dir)
    env <- new.env()
    sys.source("analysis.R", envir = env)
    as.numeric(stats::logLik(env$fit))
  }, args = list(dir), libpath = .libPaths())
  ref <- flexCountReg::countreg(
    Total_crashes ~ lnaadt + speed50 + offset(log(Length)),
    data = load_bundled_data("washington_roads"), family = "NB2",
    method = "BFGS")
  expect_equal(ll, as.numeric(logLik(ref)), tolerance = 1e-8)
})
