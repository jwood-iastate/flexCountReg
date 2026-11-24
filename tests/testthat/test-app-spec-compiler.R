# App: model specification, validation, settings files and code compiler.

skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("DT")
skip_if_not_installed("callr")
skip_if_not_installed("jsonlite")
skip_if_not_installed("zip")

eval_fit <- function(spec, data) {
  env <- new.env()
  assign(spec$data$object, data, envir = env)
  run_compiled(compile_fit(spec, "fit_1"), env)
  env
}

test_that("the family table matches the package registry", {
  ft <- fcr_family_table()
  expect_identical(ft$id, .fcr_families)
  expect_true(all(nzchar(ft$name)))
  expect_identical(ft$param1[ft$id == "NB2"], "ln(alpha)")
  expect_true(is.na(ft$param1[ft$id == "POISSON"]))
  expect_identical(ft$param1[ft$id == "DELAPORTE"], "logit(rho)")
  expect_setequal(ft$id[ft$mixture_status == "validated"],
                  c("POISSON", "NB1", "NB2"))
  expect_true(all(.app_estimators$id %in% getNamespaceExports("flexCountReg")))
})

test_that("unusual column names are compiled as symbols", {
  d <- data.frame(`my y` = rpois(30, 2), `x "q"` = rnorm(30),
                  `a\`b` = runif(30, 1, 2), check.names = FALSE)
  s <- app_spec(data = list(source = "upload", name = "d.csv",
                            object = "data_1"),
                recipe = list(response = "my y", predictors = "x \"q\"",
                              exposure = "a`b"),
                model = list(family = "POISSON"))
  expect_length(validate_spec(s, d), 0)
  expect_identical(all.vars(compile_formula(s)), c("my y", "x \"q\"", "a`b"))
  env <- eval_fit(s, d)
  ref <- countreg(
    stats::as.formula("`my y` ~ `x \"q\"` + offset(log(`a\\`b`))"),
    data = d, family = "POISSON", method = "BFGS")
  expect_equal(as.numeric(logLik(env$fit_1)), as.numeric(logLik(ref)))
  expect_silent(parse(text = draft_script(s)))
})

test_that("column names that look like code stay names", {
  evil <- "system('echo hacked')"
  d <- data.frame(y = rpois(20, 1), x = rnorm(20))
  names(d)[2] <- evil
  s <- app_spec(recipe = list(response = "y", predictors = evil),
                model = list(family = "POISSON"))
  expect_identical(all.vars(compile_formula(s)), c("y", evil))
  expect_false("system" %in% all.names(compile_fit(s, "fit_1")$fit))
})

test_that("every model type compiles to a call that fits", {
  skip_on_cran()
  d <- load_bundled_data("washington_roads")
  base <- list(source = "bundled", name = "washington_roads",
               object = "data_1")
  rec <- list(response = "Total_crashes", predictors = c("lnaadt",
                                                         "speed50"))
  cases <- list(
    list(model = list(family = "PIG2"), cls = "countreg"),
    list(model = list(estimator = "countreg_rp", family = "POISSON",
                      rpar = "speed50", ndraws = 50), cls = "countreg_rp"),
    list(model = list(estimator = "countreg_rp", family = "POISSON",
                      rpar_intercept = TRUE, ndraws = 50,
                      method = "BFGS"), cls = "countreg_rp"),
    list(model = list(estimator = "renb"), group = "ID", cls = "renb"),
    list(model = list(estimator = "repois"), group = "ID", cls = "repois"),
    list(model = list(estimator = "poisLind.re"), group = "ID",
         cls = "poislind_re"),
    list(model = list(estimator = "negmultinom", nm_family = "NM1"),
         group = "ID", cls = "negmultinom"),
    list(model = list(estimator = "countreg_mix", family = "POISSON",
                      components = 2, n_starts = 2), cls = "countreg_mix"))
  for (cs in cases) {
    s <- app_spec(data = base, recipe = c(rec, list(group = cs$group)),
                  model = cs$model)
    expect_length(validate_spec(s, d), 0)
    env <- suppressWarnings(eval_fit(s, d))
    expect_s3_class(env$fit_1, cs$cls)
    expect_s3_class(env$fit_1_summary, "summary.flexCountReg")
  }
})

test_that("the hierarchical, fixed-effects and smooth options compile and fit", {
  skip_on_cran()
  skip_if_not_installed("mgcv")
  skip_if_not_installed("reformulas")
  d <- load_bundled_data("washington_roads")
  d$Year_f <- d$Year
  base <- list(source = "bundled", name = "washington_roads",
               object = "data_1")
  rec <- list(response = "Total_crashes", predictors = c("lnaadt",
                                                         "speed50"),
              log_offset = "lnlength")
  cases <- list(
    # Hierarchical: random intercept and slope, crossed years
    list(recipe = list(group = "ID", group2 = "Year"),
         model = list(estimator = "countreg_hier", family = "POISSON",
                      re_slopes = "lnaadt", re_correlated = FALSE),
         cls = "countreg_hier",
         formula = quote(Total_crashes ~ lnaadt + speed50 +
                           offset(lnlength) + (1 + lnaadt || ID) +
                           (1 | Year))),
    # Hierarchical with gamma effects
    list(recipe = list(group = "ID"),
         model = list(estimator = "countreg_hier", family = "POISSON",
                      re_dist = "gamma"),
         cls = "countreg_hier", args = list(re_dist = "gamma")),
    # Fixed effects in two dimensions, cluster-robust by default
    list(recipe = list(group = "ID", group2 = "Year"),
         model = list(estimator = "countreg_fe", family = "POISSON",
                      stderr = "cluster"),
         cls = "countreg_fe", args = list(fixef = quote(~ID + Year),
                                          stderr = "cluster")),
    # Smooth effect, cluster-robust SEs and sampling weights
    list(recipe = list(smooth = "lnaadt", cluster = "ID",
                       weights = "Length"),
         model = list(family = "NB2", stderr = "cluster",
                      weight_type = "sampling"),
         cls = "countreg",
         formula = quote(Total_crashes ~ s(lnaadt) + speed50 +
                           offset(lnlength)),
         args = list(weight_type = "sampling", cluster = quote(~ID))),
    # Correlated random effects in a random-effects panel model
    list(recipe = list(group = "ID"),
         model = list(estimator = "repois", cre = "lnaadt"),
         cls = "repois", args = list(cre = quote(~lnaadt))))
  for (cs in cases) {
    s <- app_spec(data = base, recipe = utils::modifyList(rec, cs$recipe),
                  model = cs$model)
    expect_length(validate_spec(s, d), 0)
    cl <- compile_fit(s, "fit_1")$fit[[3]]
    if (!is.null(cs$formula)) expect_identical(cl$formula, cs$formula)
    for (a in names(cs$args)) expect_identical(cl[[a]], cs$args[[a]])
    env <- suppressMessages(suppressWarnings(eval_fit(s, d)))
    expect_s3_class(env$fit_1, cs$cls)
    expect_silent(parse(text = draft_script(s)))
  }
  # Nested second level
  s <- app_spec(data = base, recipe = utils::modifyList(
    rec, list(group = "ID", group2 = "Year")),
    model = list(estimator = "countreg_hier", group2_type = "nested"))
  expect_identical(compile_formula(s)[[3]][[3]],
                   quote((1 | ID:Year)))
})

test_that("new options are validated", {
  d <- load_bundled_data("washington_roads")
  chk <- function(pattern, recipe = list(), model = list()) {
    s <- app_spec(recipe = utils::modifyList(
      list(response = "Total_crashes", predictors = c("lnaadt", "speed50")),
      recipe), model = model)
    expect_match(validate_spec(s, d), pattern, all = FALSE)
  }
  chk("numeric predictors", list(smooth = "speed50", factors = "speed50"))
  chk("at least 10 distinct", list(smooth = "speed50"))
  chk("must be predictors", list(smooth = "Length"))
  chk("hierarchical and fixed-effects", list(group2 = "Year"))
  chk("must differ", list(group = "ID", group2 = "ID"),
      list(estimator = "countreg_hier"))
  chk("Random slopes apply", model = list(re_slopes = "lnaadt"))
  chk("single grouping column", list(group = "ID", group2 = "Year"),
      list(estimator = "countreg_hier", re_dist = "uniform"))
  chk("random intercept only", list(group = "ID"),
      list(estimator = "countreg_hier", re_dist = "gamma",
           re_slopes = "lnaadt"))
  chk("not available for this model", list(group = "ID"),
      list(estimator = "renb", stderr = "cluster"))
  chk("cluster column", model = list(stderr = "cluster"))
  chk("only with cluster-robust", list(cluster = "ID"))
  chk("need a weights column", model = list(weight_type = "sampling"))
  chk("Sampling weights are available", list(group = "ID",
                                             weights = "Length"),
      list(estimator = "renb", weight_type = "sampling"))
  chk("Correlated random effects .* not available", list(group = "ID"),
      list(estimator = "poisLind.re", cre = "lnaadt"))
  chk("need the panel identifier", model = list(estimator = "countreg_rp",
                                                rpar = "speed50",
                                                cre = "lnaadt"))
  # Settings files keep the new fields, and schema 3 files still load
  s <- app_spec(recipe = list(response = "y", predictors = c("a", "b"),
                              smooth = "a", group = "g", group2 = "h"),
                model = list(estimator = "countreg_hier",
                             re_slopes = c("a", "b"), cre = "b",
                             re_dist = "uniform", group2_type = "nested"))
  r <- spec_from_json(spec_to_json(s))
  expect_identical(r$recipe$smooth, "a")
  expect_identical(r$recipe$group2, "h")
  expect_identical(r$model$re_slopes, c("a", "b"))
  expect_identical(r$model$re_dist, "uniform")
  old <- jsonlite::fromJSON(spec_to_json(app_spec(
    recipe = list(response = "y"))), simplifyVector = TRUE)
  old$schema <- 3L
  old$recipe[c("smooth", "group2", "cluster")] <- NULL
  old$model[c("weight_type", "cre", "re_slopes", "re_correlated",
              "group2_type", "re_dist", "fe_jackknife")] <- NULL
  r3 <- .spec_from_list(old)
  expect_identical(r3$schema, .spec_schema)
  expect_identical(r3$model$re_dist, "normal")
  expect_identical(r3$recipe$smooth, character())
})

test_that("random parameters leave the fixed part of the formula", {
  s <- app_spec(recipe = list(response = "y", predictors = c("a", "b")),
                model = list(estimator = "countreg_rp", rpar = "b"))
  cl <- compile_fit(s, "fit_1")$fit[[3]]
  expect_identical(cl$formula, quote(y ~ a))
  expect_identical(cl$rpar_formula, quote(~-1 + b))
  s$model$rpar_intercept <- TRUE
  cl <- compile_fit(s, "fit_1")$fit[[3]]
  expect_identical(cl$formula, quote(y ~ a - 1))
  expect_identical(cl$rpar_formula, quote(~b))
})

test_that("validation reports problems before fitting", {
  d <- load_bundled_data("abalone_rings")
  ok <- app_spec(recipe = list(response = "rings",
                               predictors = c("shell_weight", "sex")))
  expect_length(validate_spec(ok, d), 0)
  chk <- function(pattern, recipe = list(), model = list()) {
    s <- app_spec(recipe = utils::modifyList(list(response = "rings"),
                                             recipe), model = model)
    expect_match(validate_spec(s, d), pattern, all = FALSE)
  }
  chk("whole numbers", list(response = "length"))
  chk("not in the data", list(predictors = "nope"))
  d$e <- 1
  d$le <- 0
  d$e0 <- 0
  chk("log it twice", list(exposure = "e", log_offset = "le"))
  chk("positive", list(exposure = "e0"))
  chk("no distribution parameter", model = list(family = "POISSON",
                                                dis1 = "length"))
  chk("does not use a panel", list(group = "sex"))
  chk("needs a panel identifier", model = list(estimator = "renb"))
  chk("one per class", model = list(estimator = "countreg_lc",
                                    family = c("POISSON", "NB2"),
                                    components = 3))
  chk("at least one random parameter",
      model = list(estimator = "countreg_rp"))
  chk("must be normal", list(predictors = "length"),
      model = list(estimator = "countreg_rp", rpar = "length",
                   rpar_dist = "ln", correlated = TRUE))
  expect_error(app_spec(model = list(familly = "NB2")), "Unknown")
})

test_that("settings round-trip through JSON and reject unknown fields", {
  s <- app_spec(recipe = list(response = "rings",
                              predictors = c("shell_weight", "sex"),
                              factors = "sex"),
                model = list(estimator = "countreg_lc", family = "NB2",
                             components = 3L, mixing = "length",
                             shared_mean = "shell_weight",
                             shared_dis1 = TRUE, seed = 7L))
  s2 <- spec_from_json(spec_to_json(s))
  expect_identical(compile_fit(s2, "f"), compile_fit(s, "f"))
  j <- jsonlite::fromJSON(spec_to_json(s), simplifyVector = FALSE)
  j$run_this <- "system('x')"
  expect_error(spec_from_json(jsonlite::toJSON(j, auto_unbox = TRUE)),
               "Unknown top-level")
  j$run_this <- NULL
  j$schema <- 99
  expect_error(spec_from_json(jsonlite::toJSON(j, auto_unbox = TRUE)),
               "schema")
  j$schema <- .spec_schema
  j$model$extra <- 1
  expect_error(spec_from_json(jsonlite::toJSON(j, auto_unbox = TRUE)),
               "Unknown specification")
})

test_that("mixture specifications compile to the latent-class call", {
  s <- app_spec(recipe = list(response = "count", predictors = "cover",
                              group = "site"),
                model = list(estimator = "countreg_lc", family = "POISSON",
                             components = 2L, mixing = "mined",
                             shared_mean = "cover", n_starts = 3L,
                             seed = 5L))
  cl <- compile_fit(s, "fit_1")$fit[[3]]
  expect_identical(cl[[1]], quote(flexCountReg::countreg_lc))
  expect_identical(cl$classes, 2L)
  expect_identical(cl$group_var, "site")
  expect_identical(cl$membership_formula, quote(~mined))
  expect_identical(cl$shared, quote(list(mean = "cover")))
  s$model$family <- c("POISSON", "NB2")
  cl <- compile_fit(s, "fit_1")$fit[[3]]
  expect_null(cl$classes)
  expect_identical(cl$family, quote(c("POISSON", "NB2")))
})

test_that("new variables are compiled from structured requests", {
  d <- data.frame(x = c(1, 2, 4), z = c(2, 2, 1), g = c("a", "b", "a"))
  apply_t <- function(t) {
    expect_length(check_transform(t, d), 0)
    env <- new.env()
    env$data_1 <- d
    run_compiled(list(compile_transform(t, "data_1")), env)
    env$data_1[[t$name]]
  }
  expect_equal(apply_t(list(type = "log", x = "x", name = "lx")), log(d$x))
  expect_equal(apply_t(list(type = "log_c", x = "x", constant = "1",
                            name = "l1")), log(d$x + 1))
  expect_equal(apply_t(list(type = "power", x = "x", constant = "2",
                            name = "x2")), d$x^2)
  expect_equal(apply_t(list(type = "standardize", x = "x", name = "zx")),
               as.vector(scale(d$x)))
  expect_equal(apply_t(list(type = "indicator", x = "x", op = ">=",
                            value = "2", name = "big")), c(0L, 1L, 1L))
  expect_equal(apply_t(list(type = "indicator_level", x = "g", value = "a",
                            name = "is_a")), c(1L, 0L, 1L))
  expect_equal(apply_t(list(type = "product", x = "x", y = "z",
                            name = "xz")), d$x * d$z)
  expect_equal(apply_t(list(type = "ratio", x = "x", y = "z", name = "r")),
               d$x / d$z)
  expect_equal(as.character(apply_t(list(type = "cut", x = "x",
                                         value = "2, 3", name = "cx"))),
               c("[-Inf,2)", "[2,3)", "[3, Inf)"))
  expect_s3_class(apply_t(list(type = "factor", x = "z", name = "fz")),
                  "factor")
  # Problems are reported instead of producing bad columns
  expect_match(check_transform(list(type = "log", x = "z", name = "x"), d),
               "already exists", all = FALSE)
  d$w <- c(0, 1, 2)
  expect_match(check_transform(list(type = "log", x = "w", name = "lw"), d),
               "positive", all = FALSE)
  expect_match(check_transform(list(type = "indicator", x = "x", op = ">",
                                    value = "system('x')", name = "i"), d),
               "number", all = FALSE)
  expect_match(check_transform(list(type = "log", x = "g", name = "lg"), d),
               "not numeric", all = FALSE)
  # A value typed by the user is data, not code
  e <- compile_transform(list(type = "indicator_level", x = "g",
                              value = "a\") ; system(\"x", name = "i"),
                         "data_1")
  expect_false("system" %in% all.names(e))
})
