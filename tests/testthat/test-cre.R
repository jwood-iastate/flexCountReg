# Correlated random effects (cre =) and cre_test().

cre_panel <- function(corr, G = 200, Tn = 5, seed = 1) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(G), each = Tn))
  xbar <- stats::rnorm(G)
  d$x <- xbar[d$id] + stats::rnorm(nrow(d), 0, 0.7)
  u <- stats::rgamma(G, 2, 2) * exp(corr * xbar - corr^2 / 2)
  d$y <- stats::rpois(nrow(d), u[d$id] * exp(0.3 + 0.5 * d$x))
  d
}

test_that("cre adds training group means and corrects the RE slope", {
  d <- cre_panel(0.8)
  plain <- repois(y ~ x, group_var = "id", data = d, method = "BFGS")
  m <- repois(y ~ x, group_var = "id", data = d, method = "BFGS",
              cre = ~ x, stderr = "robust")
  expect_identical(unname(m$model$cre$cols), "mean_x")
  expect_true("mean_x" %in% names(coef(m)))
  expect_equal(m$data$mean_x, stats::ave(d$x, d$id))
  fe <- suppressMessages(countreg_fe(y ~ x, fixef = ~ id, data = d))
  # Correlated effects bias the plain RE slope; the CRE slope is close to FE
  expect_gt(abs(coef(plain)[["x"]] - coef(fe)[["x"]]),
            3 * abs(coef(m)[["x"]] - coef(fe)[["x"]]))
  t <- cre_test(m)
  expect_s3_class(t, "htest")
  expect_lt(t$p.value, 0.01)
  b <- coef(m)[["mean_x"]]
  expect_equal(unname(t$statistic), b^2 / vcov(m)["mean_x", "mean_x"])
  # Uncorrelated effects: no evidence against RE
  u <- repois(y ~ x, group_var = "id", data = cre_panel(0, seed = 3),
              method = "BFGS", cre = ~ x, stderr = "robust")
  expect_gt(cre_test(u)$p.value, 0.01)
  expect_error(cre_test(plain), "not fitted with `cre`")
})

test_that("predictions use the training means of known groups", {
  d <- cre_panel(0.5)
  m <- repois(y ~ x, group_var = "id", data = d, method = "BFGS",
              cre = ~ x)
  nd <- d[d$id == 1, ]
  p0 <- predict(m, newdata = nd)
  nd2 <- nd
  nd2$x[1] <- nd2$x[1] + 1
  p1 <- predict(m, newdata = nd2)
  # Only the within effect: the group's mean stays at its training value
  expect_equal(p1[1] / p0[1], exp(coef(m)[["x"]]))
  expect_equal(p1[-1], p0[-1])
  # A new group gets the mean of its own rows
  nd3 <- nd
  nd3$id <- 999
  xb <- coef(m)[["(Intercept)"]] + coef(m)[["x"]] * nd3$x +
    coef(m)[["mean_x"]] * mean(nd3$x)
  expect_equal(predict(m, newdata = nd3), exp(xb), tolerance = 1e-8)
  nd4 <- nd
  nd4$x <- NULL
  expect_error(predict(m, newdata = nd4), "needs the column")
})

test_that("cre works with renb() and countreg_rp() and checks inputs", {
  skip_on_cran()
  d <- cre_panel(0.5, G = 150)
  nb <- suppressWarnings(renb(y ~ x, group_var = "id", data = d,
                              cre = ~ x))
  expect_true("mean_x" %in% names(coef(nb)))
  rp <- suppressWarnings(countreg_rp(
    y ~ x - 1, rpar_formula = ~ 1, data = d, group_var = "id",
    family = "POISSON", ndraws = 30, cre = ~ x))
  expect_true("mean_x" %in% names(coef(rp)))
  expect_true(is.finite(cre_test(rp)$statistic))
  expect_error(repois(y ~ x, group_var = "id", data = d, cre = "x"),
               "one-sided formula")
  expect_error(repois(y ~ x, group_var = "id", data = d, cre = ~ nope),
               "not found")
  d$f <- factor(d$id %% 2)
  expect_error(repois(y ~ x, group_var = "id", data = d, cre = ~ f),
               "numeric")
})
