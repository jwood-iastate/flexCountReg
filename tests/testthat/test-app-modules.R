# App: module servers, with shiny::testServer().

skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("DT")
skip_if_not_installed("callr")
skip_if_not_installed("jsonlite")
skip_if_not_installed("zip")

test_that("loading data and creating variables are logged", {
  state <- new_app_state()
  shiny::testServer(mod_data_server, args = list(state = state), {
    session$setInputs(source = "bundled", bundled = "stream_salamanders")
    session$setInputs(use = 1)
    expect_identical(state$data_meta$object, "data_1")
    expect_equal(nrow(state$data), 644)
    session$setInputs(tf_type = "indicator_level", tf_x = "mined",
                      tf_level = "yes", tf_name = "mined_yes",
                      tf_op = ">", tf_value = "0", tf_constant = "1",
                      tf_y = "", tf_breaks = "")
    session$setInputs(tf_add = 1)
    expect_true("mined_yes" %in% names(state$data))
    expect_equal(sum(state$data$mined_yes),
                 sum(state$data$mined == "yes"))
    # An invalid request is refused with a notification, not an error
    session$setInputs(tf_type = "log", tf_x = "cover", tf_name = "lc")
    session$setInputs(tf_add = 2)
    expect_false("lc" %in% names(state$data))
    expect_identical(vapply(state$log, `[[`, "", "kind"),
                     c("data", "transform"))
    # A second dataset gets its own object
    session$setInputs(bundled = "abalone_rings")
    session$setInputs(use = 2)
    expect_identical(state$data_meta$object, "data_2")
  })
  other <- new_app_state()
  expect_length(shiny::isolate(other$log), 0)
})

test_that("editing the build controls changes the draft, not the fits", {
  state <- new_app_state()
  shiny::isolate(app_load_data(state, "bundled", "abalone_rings",
                               response = "rings"))
  shiny::testServer(mod_build_server, args = list(state = state), {
    session$setInputs(estimator = "countreg", family = "NB2",
                      family_each = NULL, components = 2, nm_family = "NM2",
                      response = "rings", predictors = "shell_weight",
                      factors = NULL, exposure = "", log_offset = "",
                      weights = "", group = "", mixing = NULL,
                      shared_mean = NULL, shared_dis1 = FALSE, rpar = NULL,
                      rpar_intercept = FALSE, rpar_dist = "n",
                      correlated = FALSE, ndraws = 500, n_starts = 3,
                      seed = 1, method = "BFGS", max_iters = 500,
                      stderr = "normal", label = "")
    s <- state$spec
    expect_identical(s$recipe$predictors, "shell_weight")
    expect_identical(s$data$object, "data_1")
    expect_length(state$jobs, 0)
    session$setInputs(dis1 = "length")
    expect_identical(state$spec$model$dis1, "length")
    session$setInputs(family = "POISSON")
    expect_length(state$spec$model$dis1, 0)
    # Panel estimators take the group; count regression ignores it
    session$setInputs(group = "sex")
    expect_null(state$spec$recipe$group)
    session$setInputs(estimator = "renb")
    expect_identical(state$spec$recipe$group, "sex")
    # An invalid draft cannot be submitted
    session$setInputs(estimator = "countreg", response = "length")
    session$setInputs(estimate = 1)
    expect_length(state$jobs, 0)
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
  expect_equal(r$shell_weight, signif(mean(d$shell_weight), 4))
  d$txt <- NA_character_
  s2 <- app_spec(recipe = list(response = "rings", predictors = "txt"))
  expect_silent(typical_row(s2, d))
})
