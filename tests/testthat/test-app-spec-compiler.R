# Specification, validation, settings files and the code compiler.

skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("DT")
skip_if_not_installed("callr")
skip_if_not_installed("jsonlite")
skip_if_not_installed("zip")

ab_spec <- function(...) {
  app_spec(data = list(source = "bundled", name = "abalone_rings"),
           recipe = list(response = "rings",
                         predictors = c("shell_weight", "sex")),
           model = list(family = "NB2", ...))
}

test_that("the family table matches the installed flexCountReg", {
  ft <- fcr_family_table()
  ids <- utils::getFromNamespace(".fcr_families", "flexCountReg")
  expect_identical(ft$id, ids)
  expect_true(all(nzchar(ft$name)))
  expect_identical(ft$param1[ft$id == "NB2"], "ln(alpha)")
  expect_true(is.na(ft$param1[ft$id == "POISSON"]))
  expect_identical(ft$param1[ft$id == "DELAPORTE"], "logit(rho)")
  expect_setequal(ft$id[ft$mixture_status == "validated"],
                  c("POISSON", "NB1", "NB2"))
})

test_that("every exposed family maps to a call that fits", {
  skip_on_cran()
  d <- load_bundled_data("washington_roads")
  d <- d[seq_len(300), ]
  for (fam in c("POISSON", "NB2", "PIG2", "DPO")) {
    s <- app_spec(data = list(source = "bundled", name = "washington_roads"),
                  recipe = list(response = "Total_crashes",
                                predictors = "lnaadt"),
                  model = list(family = fam))
    expect_length(validate_spec(s, d), 0)
    env <- new.env()
    env$analysis_data <- d
    eval(compile_fit(s), env)
    expect_s3_class(env$fit, "countreg")
    expect_identical(env$fit$model$family, fam)
  }
})

test_that("unusual column names are compiled as symbols", {
  d <- data.frame(`my y` = rpois(30, 2), `x "q"` = rnorm(30),
                  `a\`b` = runif(30, 1, 2), check.names = FALSE)
  s <- app_spec(data = list(source = "upload", name = "d.csv"),
                recipe = list(response = "my y", predictors = "x \"q\"",
                              exposure = "a`b"),
                model = list(family = "POISSON"))
  expect_length(validate_spec(s, d), 0)
  f <- compile_formula(s)
  expect_identical(all.vars(f), c("my y", "x \"q\"", "a`b"))
  env <- new.env()
  env$analysis_data <- d
  eval(compile_fit(s), env)
  ref <- flexCountReg::countreg(
    stats::as.formula(paste0("`my y` ~ `x \"q\"` + offset(log(`a\\`b`))")),
    data = d, family = "POISSON", method = "BFGS")
  expect_equal(as.numeric(logLik(env$fit)), as.numeric(logLik(ref)))
  # The script deparses the same expressions (and parses back)
  txt <- analysis_script(s)
  expect_silent(parse(text = txt))
})

test_that("the compiler never builds code by pasting user text", {
  # A column name that is R code stays a name
  evil <- "system('echo hacked')"
  d <- data.frame(y = rpois(20, 1), x = rnorm(20))
  names(d)[2] <- evil
  s <- app_spec(data = list(source = "upload", name = "d.csv"),
                recipe = list(response = "y", predictors = evil),
                model = list(family = "POISSON"))
  f <- compile_formula(s)
  expect_identical(all.vars(f), c("y", evil))
  expect_false(any(grepl("^system$", all.names(compile_fit(s)))))
})

test_that("validation reports problems before fitting", {
  d <- load_bundled_data("abalone_rings")
  expect_length(validate_spec(ab_spec(), d), 0)
  p <- validate_spec(app_spec(recipe = list(response = "length")), d)
  expect_match(p, "whole numbers", all = FALSE)
  p <- validate_spec(app_spec(recipe = list(response = "rings",
                                            predictors = "nope")), d)
  expect_match(p, "not in the data", all = FALSE)
  d$e <- 1
  d$le <- 0
  p <- validate_spec(app_spec(recipe = list(response = "rings",
                                            exposure = "e",
                                            log_offset = "le")), d)
  expect_match(p, "log it twice", all = FALSE)
  d$e0 <- 0
  p <- validate_spec(app_spec(recipe = list(response = "rings",
                                            exposure = "e0")), d)
  expect_match(p, "positive", all = FALSE)
  p <- validate_spec(app_spec(recipe = list(response = "rings"),
                              model = list(family = "POISSON",
                                           dis1 = "length")), d)
  expect_match(p, "no distribution parameter", all = FALSE)
  p <- validate_spec(app_spec(recipe = list(response = "rings",
                                            group = "sex")), d)
  expect_match(p, "mixture estimators", all = FALSE)
  p <- validate_spec(app_spec(recipe = list(response = "rings"),
                              model = list(estimator = "countreg_lc",
                                           family = c("POISSON", "NB2"),
                                           components = 3)), d)
  expect_match(p, "one per class", all = FALSE)
  expect_error(app_spec(model = list(familly = "NB2")), "Unknown")
})

test_that("settings round-trip through JSON and reject unknown fields", {
  s <- app_spec(data = list(source = "bundled", name = "abalone_rings"),
                recipe = list(response = "rings",
                              predictors = c("shell_weight", "sex"),
                              factors = "sex"),
                model = list(estimator = "countreg_lc", family = "NB2",
                             components = 3L, mixing = "length",
                             shared_mean = "shell_weight",
                             shared_dis1 = TRUE, seed = 7L))
  s2 <- spec_from_json(spec_to_json(s))
  expect_identical(compile_analysis(s2), compile_analysis(s))
  expect_identical(analysis_script(s2), analysis_script(s))
  j <- jsonlite::fromJSON(spec_to_json(s), simplifyVector = FALSE)
  j$run_this <- "system('x')"
  expect_error(spec_from_json(jsonlite::toJSON(j, auto_unbox = TRUE)),
               "Unknown top-level")
  j$run_this <- NULL
  j$schema <- 99
  expect_error(spec_from_json(jsonlite::toJSON(j, auto_unbox = TRUE)),
               "schema")
  j$schema <- 1
  j$model$extra <- 1
  expect_error(spec_from_json(jsonlite::toJSON(j, auto_unbox = TRUE)),
               "Unknown specification")
})

test_that("mixture specifications compile to the latent-class call", {
  s <- app_spec(data = list(source = "bundled", name = "stream_salamanders"),
                recipe = list(response = "count", predictors = "cover",
                              group = "site"),
                model = list(estimator = "countreg_lc", family = "POISSON",
                             components = 2L, mixing = "mined",
                             shared_mean = "cover", n_starts = 3L,
                             seed = 5L))
  cl <- compile_fit(s)[[3]]
  expect_identical(cl[[1]], quote(flexCountReg::countreg_lc))
  expect_identical(cl$classes, 2L)
  expect_identical(cl$group_var, "site")
  expect_identical(cl$membership_formula, quote(~mined))
  expect_identical(cl$shared, quote(list(mean = "cover")))
  # Different family per class: no `classes` argument (taken from family)
  s$model$family <- c("POISSON", "NB2")
  cl <- compile_fit(s)[[3]]
  expect_null(cl$classes)
  expect_identical(cl$family, quote(c("POISSON", "NB2")))
})
