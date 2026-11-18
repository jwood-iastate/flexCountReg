# Module servers, with shiny::testServer().

skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("DT")
skip_if_not_installed("callr")
skip_if_not_installed("jsonlite")
skip_if_not_installed("zip")

test_that("loading a dataset updates the session state only", {
  state <- shiny::isolate(new_app_state())
  shiny::testServer(mod_data_server, args = list(state = state), {
    session$setInputs(source = "bundled", bundled = "stream_salamanders")
    session$setInputs(use = 1)
    expect_identical(state$data_meta$name, "stream_salamanders")
    expect_equal(nrow(state$data), 644)
  })
  # A second session starts empty
  other <- shiny::isolate(new_app_state())
  expect_null(shiny::isolate(other$data))
})

test_that("editing the build controls changes the draft, not the fits", {
  state <- new_app_state()
  shiny::isolate({
    state$data <- load_bundled_data("abalone_rings")
    state$data_meta <- list(source = "bundled", name = "abalone_rings",
                            response = "rings")
  })
  shiny::testServer(mod_build_server, args = list(state = state), {
    session$setInputs(estimator = "countreg", family = "NB2",
                      family_each = NULL, components = 2,
                      response = "rings", predictors = "shell_weight",
                      factors = NULL, exposure = "", log_offset = "",
                      weights = "", group = "", mixing = NULL,
                      shared_mean = NULL, shared_dis1 = FALSE, n_starts = 3,
                      seed = 1, method = "BFGS", max_iters = 500,
                      stderr = "normal", label = "")
    s <- state$spec
    expect_identical(s$recipe$predictors, "shell_weight")
    expect_identical(s$model$family, "NB2")
    expect_length(state$jobs, 0)
    # A family without a distribution parameter drops its predictors
    session$setInputs(dis1 = "length")
    expect_identical(state$spec$model$dis1, "length")
    session$setInputs(family = "POISSON")
    expect_length(state$spec$model$dis1, 0)
    # An invalid draft cannot be submitted
    session$setInputs(response = "length")
    session$setInputs(estimate = 1)
    expect_length(state$jobs, 0)
  })
})

test_that("a finished older job does not replace a newer active fit", {
  skip_on_cran()
  s <- app_spec(data = list(source = "bundled", name = "washington_roads"),
                recipe = list(response = "Total_crashes",
                              predictors = "lnaadt"),
                model = list(family = "POISSON"))
  wait <- function(j) {
    while (identical((j <- poll_fit_job(j))$status, "running")) Sys.sleep(0.2)
    j
  }
  j_old <- start_fit_job(s, label = "old")
  j_new <- start_fit_job(s, label = "new")
  state <- new_app_state()
  shiny::isolate({
    state$data_meta <- list(source = "bundled", name = "washington_roads")
    state$last_submitted <- j_new$id
    # The newer job finishes first, the older one afterwards
    record_finished_job(state, wait(j_new))
    record_finished_job(state, wait(j_old))
    expect_identical(state$active, j_new$id)
    expect_setequal(names(state$fits), c(j_old$id, j_new$id))
    # Each bundle keeps its own request id and label
    expect_identical(state$fits[[j_old$id]]$label, "old")
  })
})

test_that("typical scenario rows contain every column the model needs", {
  d <- load_bundled_data("abalone_rings")
  s <- app_spec(recipe = list(response = "rings",
                              predictors = c("shell_weight", "sex")),
                model = list(estimator = "countreg_lc", mixing = "length"))
  r <- typical_row(s, d)
  expect_setequal(names(r), c("shell_weight", "sex", "length"))
  expect_true(is.factor(r$sex))
  expect_equal(r$shell_weight, mean(d$shell_weight))
})
