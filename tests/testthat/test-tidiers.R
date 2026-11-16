# tidy(), glance() and augment() methods (generics/broom conventions).

test_that("tidy follows broom conventions with a component column", {
  m <- countreg(Total_crashes ~ lnaadt + lnlength, washington_roads,
                underreport_formula = ~ ShouldWidth04, method = "BFGS")
  td <- tidy(m)
  expect_named(td, c("term", "component", "estimate", "std.error",
                     "statistic", "p.value"))
  expect_equal(nrow(td), length(coef(m)))
  expect_equal(td$component,
               c("mean", "mean", "mean", "distribution", "underreporting",
                 "underreporting"))
  expect_equal(td$std.error, unname(.fcr_se(vcov(m))))
  expect_true(all(is.na(td$std.error) | td$std.error > 0))
  ci <- tidy(m, conf.int = TRUE, conf.level = 0.9)
  expect_equal(ci$conf.low, unname(confint(m, level = 0.9)[, 1]))
  ex <- tidy(m, exponentiate = TRUE)
  expect_equal(ex$estimate[1:3], unname(exp(coef(m)[1:3])))
  expect_equal(ex$estimate[4], unname(coef(m)[4]))
})

test_that("tidy components for random-parameter models", {
  m <- countreg_rp(Total_crashes ~ lnlength, rpar_formula = ~ -1 + speed50,
                   het_mean_formula = ~ ShouldWidth04,
                   data = washington_roads, ndraws = 30, method = "BFGS",
                   max.iters = 50)
  td <- tidy(m)
  expect_equal(td$component[td$term == "speed50:Mean"], "random_mean")
  expect_equal(td$component[td$term == "speed50:St.Dev"], "random_scale")
  expect_equal(td$component[grepl("^HetMean", td$term)], "heterogeneity")
  expect_equal(td$component[td$term == "ln(alpha)"], "distribution")
})

test_that("tidy works for every estimator", {
  d <- washington_roads
  for (f in list(
    renb(Total_crashes ~ lnaadt, group_var = "ID", data = d, method = "BFGS"),
    poislind_re(Total_crashes ~ lnaadt, group_var = "ID", data = d,
                method = "BFGS"),
    negmultinom(Total_crashes ~ lnaadt, group_var = "ID", data = d))) {
    td <- tidy(f)
    expect_equal(nrow(td), length(coef(f)))
    expect_false(anyNA(td$estimate))
  }
})

test_that("glance returns one row of model statistics", {
  m <- countreg(Total_crashes ~ lnaadt, washington_roads, method = "BFGS")
  g <- glance(m)
  expect_equal(nrow(g), 1)
  expect_equal(g$AIC, AIC(m))
  expect_equal(g$nobs, 1501L)
  expect_true(g$converged)
  expect_equal(g$vcov_type, "hessian")
})

test_that("augment adds fitted values and residuals", {
  d <- washington_roads
  d$lnaadt[3] <- NA
  m <- suppressMessages(countreg(Total_crashes ~ lnaadt, d, method = "BFGS"))
  a <- augment(m)
  expect_equal(nrow(a), nrow(d) - 1)
  expect_equal(a$.fitted, fitted(m))
  expect_equal(a$.resid, residuals(m))
  nd <- d[1:5, ]
  a2 <- augment(m, new_data = nd)
  expect_equal(nrow(a2), 5)
  expect_true(is.na(a2$.fitted[3]))
  expect_equal(a2$.fitted, predict(m, new_data = nd))
  a3 <- augment(m, new_data = nd[, setdiff(names(nd), "Total_crashes")])
  expect_false(".resid" %in% names(a3))
})
