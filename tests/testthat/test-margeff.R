# margEffTable(): marginal effects, elasticities and their standard errors,
# checked against closed forms for a log-linear mean.

nb2 <- function() {
  countreg(Total_crashes ~ lnaadt + speed50, washington_roads,
           offset = "lnlength", method = "BFGS")
}

test_that("AMEs and pseudo-elasticities match their closed forms", {
  m <- nb2()
  b <- coef(m)
  mu <- predict(m, newdata = washington_roads)
  me <- margEffTable(m, digits = 8)
  expect_s3_class(me, "flexEffectTable")
  expect_equal(me$term, c("lnaadt", "speed50"))
  expect_equal(me$variable_type, c("continuous", "indicator"))
  # d mu / d x = beta * mu for a log-linear mean
  expect_equal(me$estimate[1], unname(b["lnaadt"] * mean(mu)),
               tolerance = 1e-5)
  # Pseudo-elasticity of an indicator: 100 (exp(beta) - 1)
  expect_equal(me$estimate[2], unname(100 * (exp(b["speed50"]) - 1)),
               tolerance = 1e-6)
  # Its delta-method standard error: 100 exp(beta) se(beta)
  se_b <- sqrt(vcov(m)["speed50", "speed50"])
  expect_equal(me$std_error[2], unname(100 * exp(b["speed50"]) * se_b),
               tolerance = 1e-4)
  expect_equal(me$lower_ci, me$estimate - qnorm(0.975) * me$std_error,
               tolerance = 1e-6)
  expect_equal(me$n_obs, rep(nrow(washington_roads), 2))
})

test_that("elasticities and discrete changes", {
  m <- nb2()
  b <- coef(m)
  d <- washington_roads
  el <- margEffTable(m, measure = "elasticity", digits = 8)
  # Elasticity of mu with respect to x is beta * x
  expect_equal(el$estimate[1], unname(b["lnaadt"] * mean(d$lnaadt)),
               tolerance = 1e-5)
  dc <- margEffTable(m, vars = "speed50", indicator = "discrete_change",
                     digits = 8)
  d1 <- d0 <- d
  d1$speed50 <- 1
  d0$speed50 <- 0
  expect_equal(dc$estimate,
               mean(predict(m, newdata = d1) - predict(m, newdata = d0)),
               tolerance = 1e-6)
  expect_equal(dc$effect_metric, "Discrete change")
})

test_that("factor and logical indicators", {
  d <- washington_roads
  d$fast <- factor(ifelse(d$speed50 == 1, "yes", "no"), levels = c("no", "yes"))
  d$wide <- d$ShouldWidth04 == 0
  m <- countreg(Total_crashes ~ lnaadt + fast + wide, d, method = "BFGS")
  me <- margEffTable(m, vars = c("fast", "wide"), digits = 8)
  expect_equal(me$variable_type, c("indicator", "indicator"))
  expect_equal(me$estimate,
               unname(100 * (exp(coef(m)[c("fastyes", "wideTRUE")]) - 1)),
               tolerance = 1e-6)
  d$three <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  m3 <- countreg(Total_crashes ~ lnaadt + three, d, method = "BFGS")
  expect_error(margEffTable(m3, vars = "three"), "more than two levels")
})

test_that("variables from auxiliary formulas are included and validated", {
  m <- countreg(Total_crashes ~ lnaadt, washington_roads,
                underreport_formula = ~ ShouldWidth04, method = "BFGS")
  me <- margEffTable(m)
  expect_equal(me$term, c("lnaadt", "ShouldWidth04"))
  # The reporting probability changes the expected observed count
  expect_gt(abs(me$estimate[2]), 0)
  expect_error(margEffTable(m, vars = "not_a_variable"), "not in `data`")
  expect_error(margEffTable(m, confint_level = 2), "confint_level")
  expect_error(margEffTable(m, bootstraps = 0), "bootstraps")
  expect_error(margEffTable(m, digits = -1), "digits")
})

test_that("bootstrap standard errors resample rows or clusters", {
  m <- nb2()
  set.seed(1)
  delta <- margEffTable(m, vars = "lnaadt", digits = 8)
  boot <- margEffTable(m, vars = "lnaadt", se = "bootstrap", bootstraps = 15,
                       digits = 8)
  expect_equal(boot$estimate, delta$estimate)
  expect_equal(boot$successful_bootstraps, 15)
  # Bootstrap and delta-method SEs are of the same size
  expect_lt(abs(log(boot$std_error / delta$std_error)), log(2.5))
  set.seed(2)
  cl <- margEffTable(m, vars = "lnaadt", se = "bootstrap", bootstraps = 10,
                     cluster_var = "ID", digits = 8)
  expect_equal(cl$successful_bootstraps, 10)
  expect_true(is.finite(cl$std_error))
  # Rows of a resampled cluster stay together: each sampled ID appears a
  # multiple of its original number of rows
  set.seed(3)
  d30 <- washington_roads[1:30, ]
  rs <- .margEff_resample_data(d30, cluster_var = "ID")
  n_orig <- table(as.character(d30$ID))
  n_boot <- table(as.character(rs$ID))
  expect_true(all(n_boot %% n_orig[names(n_boot)] == 0))
  expect_equal(length(unique(as.character(rs$ID))) <= length(n_orig), TRUE)
})

test_that("random-parameter models use Exact predictions by default", {
  rp <- countreg_rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                    data = washington_roads, ndraws = 50, method = "BFGS",
                    max.iters = 100)
  me <- margEffTable(rp, vars = "lnaadt")
  expect_equal(attr(me, "effect_info")$pred_method, "Exact")
  expect_output(print(me), "Prediction method: Exact")
})

test_that("gt and LaTeX output", {
  m <- nb2()
  l <- margEffTable(m, tableType = "latex")
  expect_match(as.character(l), "tabular")
  skip_if_not_installed("gt")
  skip_if_not(inherits(gt::gt(data.frame(a = 1)), "gt_tbl"),
              "a working gt installation is required")
  g <- margEffTable(m, tableType = "gt")
  expect_s3_class(g, "gt_tbl")
})

# Evaluate `expr`, expecting one or more warnings that all match `pattern`
# (margEffTable() warns once per variable).
expect_warnings <- function(expr, pattern) {
  ws <- character()
  out <- withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_gt(length(ws), 0)
  expect_match(ws, pattern, all = TRUE)
  out
}

test_that("delta-method SEs fall back to the Hessian or the stored SEs", {
  m <- nb2()
  ref <- margEffTable(m, digits = 8)
  # Without the covariance matrix: the inverse negative Hessian
  m_h <- m
  m_h$model$vcov <- NULL
  m_h$model$hessian <- -solve(vcov(m))
  expect_equal(margEffTable(m_h, digits = 8)$std_error, ref$std_error,
               tolerance = 1e-6)
  # Without the Hessian: a diagonal approximation, with a warning
  m_d <- m_h
  m_d$model$hessian <- NULL
  me <- expect_warnings(margEffTable(m_d, digits = 8), "diagonal covariance")
  # The pseudo-elasticity depends on one coefficient only, so it is unchanged
  expect_equal(me$std_error[2], ref$std_error[2], tolerance = 1e-6)
  # Nothing at all: NA standard errors with a warning
  m_n <- m_d
  m_n$model$se <- NULL
  me <- expect_warnings(margEffTable(m_n), "No covariance matrix")
  expect_true(all(is.na(me$std_error)))
  # Mismatched dimensions
  m_x <- m
  m_x$model$vcov <- diag(2)
  me <- expect_warnings(margEffTable(m_x), "dimensions do not match")
  expect_true(all(is.na(me$std_error)))
})

test_that("margEffTable needs data and at least one covariate", {
  m <- countreg(Total_crashes ~ lnaadt, washington_roads, method = "BFGS",
                keep_data = FALSE)
  expect_error(margEffTable(m), "No data supplied")
  expect_equal(margEffTable(m, data = washington_roads)$term, "lnaadt")
  m0 <- countreg(Total_crashes ~ 1, washington_roads, offset = "lnlength",
                 method = "BFGS")
  expect_error(margEffTable(m0), "No eligible variables")
})
