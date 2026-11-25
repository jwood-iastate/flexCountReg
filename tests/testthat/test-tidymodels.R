# tidymodels integration: engines, workflows, tuning, metrics.

roads_small <- function() {
  ids <- unique(washington_roads$ID)[1:150]
  d <- tibble::as_tibble(washington_roads[washington_roads$ID %in% ids, ])
  set.seed(4)
  d$c1 <- stats::rbinom(nrow(d), d$Total_crashes, 0.3)
  d$c2 <- d$Total_crashes - d$c1
  d
}

engine_cases <- function() {
  f <- Total_crashes ~ lnaadt
  re <- list(group_var = "ID", method = "BFGS")
  list(
    flexCountReg = list(args = list(method = "BFGS", offset = "lnlength"),
                        formula = f, fun = "countreg"),
    countreg_rp = list(
      args = list(rpar_formula = ~ -1 + lnlength, ndraws = 25,
                  max.iters = 30), formula = f),
    countreg_mix = list(args = list(n_starts = 2, seed = 1), formula = f),
    countreg_lc = list(args = list(n_starts = 2, seed = 1, group_var = "ID"),
                       formula = f),
    countreg_fe = list(args = list(fixef = ~ ID), formula = f,
                       predict_args = list(newgroup = "mean_fe")),
    countreg_hier = list(args = list(),
                         formula = Total_crashes ~ lnaadt + (1 | ID),
                         predict_args = list(marginal = TRUE)),
    countreg_alloc = list(args = list(), formula = cbind(c1, c2) ~ lnaadt,
                          outcomes = c("c1", "c2")),
    countreg_mv = list(args = list(), formula = cbind(c1, c2) ~ lnaadt,
                       outcomes = c("c1", "c2")),
    renb = list(args = re, formula = f),
    repois = list(args = re, formula = f),
    poislind_re = list(args = re, formula = f),
    negmultinom = list(args = list(group_var = "ID"), formula = f))
}

quiet <- function(x) suppressMessages(suppressWarnings(x))

test_that("every estimator is a registered engine", {
  skip_if_not_installed("parsnip")
  count_reg()
  eng <- parsnip::show_engines("count_reg")$engine
  expect_setequal(eng, names(engine_cases()))
  # Each estimator named in the registry exists and is exported
  funs <- vapply(.fcr_parsnip_engines(), `[[`, "", "fun")
  expect_true(all(funs %in% getNamespaceExports("flexCountReg")))
})

for (eng in names(engine_cases())) {
  test_that(paste("engine", eng, "follows the tidymodels contract"), {
    skip_if_not_installed("parsnip")
    skip_if_not_installed("workflows")
    k <- engine_cases()[[eng]]
    d <- roads_small()
    y <- k$outcomes %||% "Total_crashes"
    spec <- do.call(parsnip::set_engine,
                    c(list(count_reg(), eng), k$args,
                      list(predict_args = k$predict_args)))
    wf <- workflows::workflow() |>
      workflows::add_variables(tidyselect::all_of(!!y),
                               c(lnaadt, lnlength, ID)) |>
      workflows::add_model(spec, formula = k$formula)
    w <- quiet(parsnip::fit(wf, d))
    m <- hardhat::extract_fit_engine(w)
    direct <- quiet(do.call(k$fun %||% eng,
                            c(list(k$formula, data = d), k$args)))
    expect_equal(coef(m), coef(direct), tolerance = 1e-8)

    new <- d[1:5, setdiff(names(d), y)]
    new$lnaadt[2] <- NA
    p <- predict(w, new)
    expect_named(p, if (is.null(k$outcomes)) ".pred" else
      paste0(".pred_", k$outcomes))
    expect_equal(nrow(p), 5)
    expect_true(is.na(p[[1]][2]))
    expect_false(anyNA(p[[1]][-2]))
    expect_equal(as.matrix(p), unname(as.matrix(do.call(predict, c(
      list(direct, newdata = new), k$predict_args)))),
      ignore_attr = TRUE)
    expect_equal(predict(w, new[1, ]), p[1, ])
    p0 <- expect_silent(predict(w, new[0, ]))
    expect_equal(dim(p0), c(0L, ncol(p)))
    expect_named(p0, names(p))
    expect_no_error(predict(w, new, type = "raw"))

    expect_s3_class(tidy(w), "tbl_df")
    expect_equal(nrow(glance(w)), 1)
    expect_equal(nrow(augment(w, d[1:4, ])), 4)
    expect_equal(nrow(augment(hardhat::extract_fit_parsnip(w), d[1:4, ])), 4)
    expect_equal(nrow(augment(m)), nobs(m))

    # Without a model formula, columns with their own role are not
    # predictors of the mean
    if (eng != "countreg_hier" && is.null(k$outcomes)) {
      vars <- unique(c("lnaadt", "ID", if (eng %in% c("flexCountReg",
                                                      "countreg_rp"))
        "lnlength"))
      w2 <- quiet(parsnip::fit(
        workflows::workflow() |>
          workflows::add_variables(Total_crashes,
                                   tidyselect::all_of(!!vars)) |>
          workflows::add_model(spec),
        if (eng %in% c("flexCountReg", "countreg_rp", "countreg_mix")) {
          transform(d, ID = 1 + (seq_along(ID) %% 2))
        } else d))
      trm <- attr(stats::terms(formula(hardhat::extract_fit_engine(w2))),
                  "term.labels")
      expect_true("lnaadt" %in% trm)
      expect_false("lnlength" %in% trm)
      if (!eng %in% c("flexCountReg", "countreg_rp", "countreg_mix")) {
        expect_false("ID" %in% trm)
      }
    }
  })
}

test_that("family can be tuned and updated", {
  skip_if_not_installed("tune")
  skip_if_not_installed("rsample")
  skip_if_not_installed("workflows")
  d <- roads_small()
  spec <- update(count_reg(), family = "NB1")
  expect_equal(rlang::eval_tidy(spec$args$family), "NB1")
  spec <- count_reg(family = hardhat::tune()) |>
    parsnip::set_engine("flexCountReg", method = "BFGS")
  wf <- workflows::workflow() |>
    workflows::add_formula(Total_crashes ~ lnaadt) |>
    workflows::add_model(spec)
  set.seed(1)
  cv <- rsample::group_vfold_cv(d, group = ID, v = 2)
  g <- tune::tune_grid(wf, cv, grid = data.frame(family = c("POISSON", "NB2")))
  expect_equal(nrow(tune::collect_notes(g)), 0)
  expect_setequal(tune::collect_metrics(g)$family, c("POISSON", "NB2"))
  f <- parsnip::fit(tune::finalize_workflow(
    wf, tune::select_best(g, metric = "rmse")), d)
  expect_s3_class(hardhat::extract_fit_engine(f), "countreg")
})

test_that("engine arguments have tuning parameters", {
  skip_if_not_installed("tune")
  skip_if_not_installed("dials")
  par_of <- function(spec, i = 1) {
    hardhat::extract_parameter_set_dials(spec)$object[[i]]
  }
  s <- count_reg() |> parsnip::set_engine("countreg_lc",
                                          classes = hardhat::tune())
  expect_equal(par_of(s)$range, list(lower = 2L, upper = 4L))
  s <- count_reg(family = hardhat::tune()) |>
    parsnip::set_engine("negmultinom", group_var = "ID")
  expect_equal(par_of(s)$values, c("NM1", "NM2", "NMP"))
  s <- count_reg(family = hardhat::tune()) |>
    parsnip::set_engine("flexCountReg", zero_type = hardhat::tune())
  expect_equal(par_of(s, 2)$values, c("inflated", "hurdle"))
  for (p in list(count_zero_type(), count_correlated(), count_re_dist())) {
    expect_s3_class(p, "qual_param")
  }
})

test_that("intervals and raw predictions match predict()", {
  skip_if_not_installed("parsnip")
  d <- roads_small()
  f <- parsnip::fit(count_reg() |>
                      parsnip::set_engine("flexCountReg", method = "BFGS"),
                    Total_crashes ~ lnaadt, d)
  nd <- d[1:4, ]
  for (ty in c("confidence", "prediction")) {
    p <- predict(f, nd, type = c(confidence = "conf_int",
                                 prediction = "pred_int")[[ty]], level = 0.8)
    q <- predict(f$fit, nd, interval = ty, level = 0.8)
    expect_named(p, c(".pred_lower", ".pred_upper"))
    expect_equal(p$.pred_lower, q$lwr)
    expect_equal(p$.pred_upper, q$upr)
  }
  expect_equal(predict(f, nd, type = "raw",
                       opts = list(type = "prob", max_count = 3)),
               predict(f$fit, nd, type = "prob", max_count = 3))
})

test_that("engine arguments are checked", {
  d <- roads_small()
  expect_error(count_reg_fit(Total_crashes ~ lnaadt, d, rpar_formula = ~ ID),
               "not used by countreg")
  expect_error(count_reg_fit(Total_crashes ~ lnaadt, d, family = "NB2",
                             model = "renb", group_var = "ID"),
               "fixed distribution")
  expect_error(count_reg_fit(Total_crashes ~ lnaadt, d, predict_args = 1),
               "must be a list")
  expect_warning(
    count_reg_fit(Total_crashes ~ lnaadt, d, method = "BFGS",
                  weights = rep(c(0.5, 1), length.out = nrow(d))),
    "sampling")
  expect_no_warning(
    count_reg_fit(Total_crashes ~ lnaadt, d, method = "BFGS",
                  weights = rep(c(0.5, 1), length.out = nrow(d)),
                  weight_type = "sampling"))
  # `engine` is an argument of set_engine() itself
  m <- count_reg_fit(Total_crashes ~ lnaadt + (1 | ID), d,
                     model = "countreg_hier", fit_engine = "laplace")
  expect_equal(m$model$engine, "laplace")
})

test_that("`.` leaves out columns with a role, and keeps the rest", {
  d <- data.frame(..y = 1, x = 1, g = 1, `a b` = 1, off = 1,
                  check.names = FALSE)
  f <- .fcr_parsnip_formula(..y ~ ., d, list(group_var = "g",
                                             offset = "off"), "renb")
  expect_equal(f[[3]], quote(x + `a b`))
  f <- .fcr_parsnip_formula(..y ~ . + (1 | g), d, list(offset = "off"),
                            "countreg_hier")
  expect_equal(f[[3]], quote((x + `a b`) + (1 | g)))
  f <- .fcr_parsnip_formula(..y ~ ., d, list(
    rpar_formula = ~ x, fixef = ~ g, offset = c("off", "a b")),
    "countreg_rp")
  expect_equal(f[[3]], 1)
  f <- ..y ~ x
  expect_identical(.fcr_parsnip_formula(f, d, list(), "countreg"), f)
})

test_that("extracted models work with update() and bootstrap refits", {
  skip_if_not_installed("parsnip")
  d <- roads_small()
  d$lnaadt[3] <- NA
  w <- rep(1:2, length.out = nrow(d))
  fit <- function(eng, f, ...) quiet(count_reg_fit(f, d, model = eng, ...))
  for (m in list(
    fit("countreg", Total_crashes ~ lnaadt, method = "BFGS", weights = w),
    fit("countreg_alloc", cbind(c1, c2) ~ lnaadt, weights = w),
    fit("countreg_mv", cbind(c1, c2) ~ lnaadt, weights = w),
    fit("repois", Total_crashes ~ lnaadt, group_var = "ID", method = "BFGS",
        cre = ~ lnaadt))) {
    expect_equal(coef(quiet(update(m))), coef(m), tolerance = 1e-4)
    expect_equal(coef(quiet(update(update(m)))), coef(m), tolerance = 1e-4)
  }
  m <- fit("countreg", Total_crashes ~ lnaadt, method = "BFGS")
  expect_equal(update(m, family = "NB1")$model$family, "NB1")
  expect_true("lnlength" %in% names(coef(update(m, . ~ . + lnlength))))
  expect_equal(nobs(update(m, data = d[1:300, ])), 299)
  # Models fitted directly still use stats::update() semantics
  m2 <- countreg(Total_crashes ~ lnaadt, d, method = "BFGS")
  expect_equal(update(m2, family = "NB1")$model$family, "NB1")
  expect_equal(names(coef(update(m2, . ~ . + lnlength)))[3], "lnlength")
  expect_true(is.call(update(m2, family = "NB1", evaluate = FALSE)))
  # Where the package is not attached
  eng <- m
  environment(eng$formula) <- baseenv()
  expect_s3_class(quiet(margEffTable(eng, se = "boot", bootstraps = 3)),
                  "flexEffectTable")
})

test_that("rmse() and mae() work as yardstick metrics", {
  y <- c(1, 2, 3)
  mu <- c(1.1, 1.9, 3.5)
  expect_equal(rmse(y, mu), sqrt(mean((y - mu)^2)))
  expect_equal(mae(y = y, mu = mu), mean(abs(y - mu)))
  skip_if_not_installed("yardstick")
  d <- data.frame(y = y, mu = mu)
  expect_equal(rmse(d, y, mu), yardstick::rmse(d, y, mu))
  expect_equal(mae(d, truth = y, estimate = mu), yardstick::mae(d, y, mu))
  expect_equal(rmse(data = d, truth = y, estimate = mu)$.estimate,
               rmse(y, mu))
  ms <- yardstick::metric_set(rmse, mae)
  expect_equal(ms(d, truth = y, estimate = mu)$.estimate,
               c(rmse(y, mu), mae(y, mu)))
  skip_if_not_installed("tune")
  skip_if_not_installed("rsample")
  r <- roads_small()
  set.seed(1)
  z <- tune::fit_resamples(
    count_reg() |> parsnip::set_engine("flexCountReg", method = "BFGS"),
    Total_crashes ~ lnaadt, rsample::vfold_cv(r, v = 2), metrics = ms)
  expect_equal(nrow(tune::collect_metrics(z)), 2)
})

test_that("predict() of every class honors new_data", {
  d <- roads_small()
  ms <- list(
    quiet(countreg_alloc(cbind(c1, c2) ~ lnaadt, d)),
    quiet(countreg_mv(cbind(c1, c2) ~ lnaadt, d)),
    quiet(countreg_hier(Total_crashes ~ lnaadt + (1 | ID), d)))
  for (m in ms) {
    expect_equal(predict(m, new_data = d[3:5, ]),
                 predict(m, newdata = d[3:5, ]))
  }
})
