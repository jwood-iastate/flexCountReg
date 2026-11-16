# S3 methods, summary object, class structure and argument harmonization.

nb2_fit <- function() {
  countreg(Total_crashes ~ lnaadt + lnlength + speed50, washington_roads,
           method = "BFGS")
}

test_that("standard methods match MASS::glm.nb", {
  m <- nb2_fit()
  mm <- MASS::glm.nb(Total_crashes ~ lnaadt + lnlength + speed50,
                     washington_roads)
  expect_s3_class(m, c("countreg", "flexCountReg"))
  expect_equal(unname(coef(m)[1:4]), unname(coef(mm)), tolerance = 1e-4)
  expect_equal(AIC(m), AIC(mm), tolerance = 1e-6)
  expect_equal(BIC(m), BIC(mm), tolerance = 1e-6)
  expect_equal(nobs(m), 1501L)
  expect_equal(as.numeric(fitted(m)), as.numeric(fitted(mm)),
               tolerance = 1e-4)
  expect_equal(as.numeric(residuals(m, "pearson")),
               as.numeric(residuals(mm, "pearson")), tolerance = 1e-3)
  expect_equal(residuals(m),
               as.numeric(washington_roads$Total_crashes) - fitted(m))
  expect_identical(rownames(vcov(m)), names(coef(m)))
  expect_equal(AIC(m), -2 * as.numeric(logLik(m)) + 2 * length(coef(m)))
  ci <- confint(m, level = 0.9)
  se <- sqrt(diag(vcov(m)))
  expect_equal(unname(ci[, 2]), unname(coef(m) + qnorm(0.95) * se))
  expect_equal(deparse(formula(m)),
               "Total_crashes ~ lnaadt + lnlength + speed50")
  expect_output(print(m), "family: NB2")
  u <- update(m, . ~ . - speed50)
  expect_false("speed50" %in% names(coef(u)))
})

test_that("summary returns an object with a print method and tibble path", {
  m <- nb2_fit()
  s <- summary(m)
  expect_s3_class(s, "summary.flexCountReg")
  expect_equal(colnames(s$coefficients)[1:4],
               c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_equal(unname(s$coefficients[, "Estimate"]), unname(coef(m)))
  expect_equal(s$AIC, AIC(m))
  expect_equal(s$vcov_type, "hessian")
  expect_false(is.null(s$lr_test))
  expect_true(s$lr_test$boundary)
  expect_output(print(s), "LR test vs Poisson")
  tb <- tibble::as_tibble(s)
  expect_s3_class(tb, "tbl_df")
  expect_equal(tb$parameter, names(coef(m)))
  mo <- countreg(Total_crashes ~ lnaadt, washington_roads,
                 offset = "lnlength", method = "BFGS")
  expect_output(print(summary(mo)), "Offset: lnlength")
  expect_false(any(grepl("Offset", rownames(summary(mo)$coefficients))))
})

test_that("Pearson residuals with underreporting use the model's variance", {
  m <- countreg(Total_crashes ~ lnaadt + lnlength, washington_roads,
                underreport_formula = ~ ShouldWidth04, method = "BFGS")
  p <- predict(m, newdata = washington_roads, type = "prob_reported")
  mu <- predict(m, newdata = washington_roads, type = "true")
  alpha <- exp(coef(m)["ln(alpha)"])
  # The observed count is NB2 with mean mu * p (as in the likelihood)
  v <- mu * p + alpha * (mu * p)^2
  expect_equal(as.numeric(residuals(m, "pearson")),
               (as.numeric(washington_roads$Total_crashes) - mu * p) /
                 sqrt(v))
})

test_that("simulate draws from the fitted distribution", {
  m <- nb2_fit()
  s1 <- simulate(m, nsim = 3, seed = 4)
  s2 <- simulate(m, nsim = 3, seed = 4)
  expect_identical(s1, s2)
  expect_equal(dim(s1), c(1501L, 3L))
  sims <- simulate(m, nsim = 300, seed = 1)
  expect_equal(mean(as.matrix(sims)), mean(fitted(m)), tolerance = 0.02)
})

test_that("every estimator returns a subclass with working methods", {
  d <- washington_roads
  fits <- list(
    countreg_rp = countreg_rp(Total_crashes ~ lnaadt,
                              rpar_formula = ~ -1 + speed50, data = d,
                              ndraws = 30, method = "BFGS", max.iters = 50),
    renb = renb(Total_crashes ~ lnaadt, group_var = "ID", data = d,
                method = "BFGS"),
    poislind_re = poislind_re(Total_crashes ~ lnaadt, group_var = "ID",
                              data = d, method = "BFGS"),
    repois = repois(Total_crashes ~ lnaadt, group_var = "ID", data = d,
                    method = "BFGS"),
    negmultinom = negmultinom(Total_crashes ~ lnaadt, group_var = "ID",
                              data = d))
  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_s3_class(f, c(nm, "flexCountReg"))
    expect_true(is.numeric(coef(f)) && !is.null(names(coef(f))))
    expect_equal(dim(vcov(f)), rep(length(coef(f)), 2))
    expect_true(is.finite(AIC(f)) && is.finite(BIC(f)))
    expect_length(fitted(f), nobs(f))
    expect_s3_class(summary(f), "summary.flexCountReg")
    expect_null(f$model$data)
  }
  expect_equal(attr(logLik(fits$renb), "nobs"), length(unique(d$ID)))
  expect_equal(predict(fits$repois, newdata = d[1:5, ]),
               as.numeric(fitted(fits$repois))[1:5], tolerance = 1e-6)
})

test_that("old argument names still work with a deprecation warning", {
  d <- washington_roads
  rlang::reset_warning_verbosity("fcr_renb_print.level")
  expect_warning(renb(Total_crashes ~ lnaadt, group_var = "ID", data = d,
                      method = "BFGS", print.level = 0), "deprecated")
  rlang::reset_warning_verbosity("fcr_countreg.rp_panel_id")
  expect_warning(m <- countreg.rp(Total_crashes ~ lnaadt,
                                  rpar_formula = ~ -1 + speed50, data = d,
                                  ndraws = 20, method = "BFGS",
                                  max.iters = 20, panel_id = "ID"),
                 "deprecated")
  expect_equal(m$model$n_groups, length(unique(d$ID)))
  # bootstraps alone still requests bootstrapped standard errors
  set.seed(1)
  b <- poislind_re(Total_crashes ~ lnaadt, group_var = "ID", data = d,
                   method = "BFGS", bootstraps = 5)
  expect_equal(b$model$vcov_type, "bootstrap")
})

test_that("robust standard errors are available for panel estimators", {
  d <- washington_roads
  for (fn in list(renb, poislind_re, repois)) {
    f <- fn(Total_crashes ~ lnaadt, group_var = "ID", data = d,
            method = "BFGS", stderr = "robust")
    expect_equal(f$model$vcov_type, "robust")
    expect_true(all(is.finite(sqrt(diag(vcov(f))))))
  }
})

test_that("keep_data = FALSE drops the data but keeps prediction working", {
  d <- washington_roads
  m <- countreg(Total_crashes ~ lnaadt, d, method = "BFGS", keep_data = FALSE)
  expect_null(m$data)
  expect_length(predict(m, newdata = d[1:3, ]), 3)
  expect_error(fitted(m), "keep_data")
})

test_that("regCompTest uses the same rows, offset and nesting rules", {
  d <- washington_roads
  m1 <- countreg(Total_crashes ~ lnaadt, d, offset = "lnlength",
                 method = "BFGS")
  g <- glm(Total_crashes ~ lnaadt + offset(lnlength), poisson, d)
  r <- regCompTest(m1, basemodel = "Poisson", variables = TRUE)
  expect_equal(r$LLbase, as.numeric(logLik(g)), tolerance = 1e-8)
  m0 <- countreg(Total_crashes ~ lnaadt, d, offset = "lnlength",
                 family = "Poisson", method = "BFGS")
  a <- regCompTest(m0, model2 = m1)
  b <- regCompTest(m1, model2 = m0)
  expect_equal(a$LR, b$LR)
  expect_gt(a$LR, 0)
  expect_equal(a$base_index, 1L)
  expect_equal(b$base_index, 2L)
  expect_error(regCompTest(m1, model2 = countreg(Total_crashes ~ lnaadt,
                                                 d[1:500, ], method = "BFGS")),
               "different numbers of observations")
})

test_that("cor2cov validates its input", {
  C <- matrix(c(1, 0.2, 0.2, 1), 2)
  expect_equal(cor2cov(C, c(2, 3)), matrix(c(4, 1.2, 1.2, 9), 2))
  expect_error(cor2cov(matrix(c(1, 0.2, 0.3, 1), 2), c(1, 1)), "symmetric")
  expect_error(cor2cov(matrix(c(2, 0.2, 0.2, 1), 2), c(1, 1)), "diagonal")
})
