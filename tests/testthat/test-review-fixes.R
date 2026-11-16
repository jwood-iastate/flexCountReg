# Regression tests for issues found in code review (PR #3).

test_that("repois() bootstrap covariance is the coefficient covariance", {
  set.seed(1)
  b <- repois(Total_crashes ~ lnaadt, group_var = "ID",
              data = washington_roads, method = "BFGS", bootstraps = 8)
  expect_equal(dim(vcov(b)), c(3L, 3L))
  expect_identical(rownames(vcov(b)), names(coef(b)))
  expect_true(all(is.finite(tidy(b)$std.error)))
  expect_equal(b$model$n_boot_ok, 8)
  expect_equal(unname(tidy(b)$std.error), unname(sqrt(diag(vcov(b)))))
})

test_that("repois() records rows used and keeps a prediction blueprint", {
  d <- washington_roads
  d$lnaadt[1:3] <- NA
  r <- suppressMessages(repois(Total_crashes ~ lnaadt, group_var = "ID",
                               data = d, method = "BFGS"))
  expect_equal(r$model$rows_used, 4:nrow(d))
  a <- augment(r)
  expect_equal(nrow(a), nrow(d) - 3)
  expect_equal(a$.fitted, as.numeric(fitted(r)))
  k <- repois(Total_crashes ~ lnaadt, group_var = "ID",
              data = washington_roads, method = "BFGS", keep_data = FALSE)
  expect_null(k$data)
  expect_length(predict(k, newdata = washington_roads[1:4, ]), 4)
  # Several offsets and formula offsets now work as in the other estimators
  a1 <- repois(Total_crashes ~ lnaadt + offset(lnlength), group_var = "ID",
               data = washington_roads, method = "BFGS")
  a2 <- repois(Total_crashes ~ lnaadt, group_var = "ID", offset = "lnlength",
               data = washington_roads, method = "BFGS")
  expect_equal(coef(a1), coef(a2), tolerance = 1e-6)
})

test_that("infinite values are not counts", {
  expect_equal(.fcr_is_count(c(1, Inf, -Inf, NaN, NA, 2.5)),
               c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
  expect_warning(p <- dplind(c(Inf, 1), mean = 1, theta = 1), "whole number")
  expect_equal(p[1], 0)
  expect_equal(suppressWarnings(dpinvgaus(Inf, 1, 1)), 0)
})

test_that("dot formulas are expanded against the data", {
  d <- washington_roads[, c("Total_crashes", "lnaadt", "lnlength", "speed50")]
  a <- countreg(Total_crashes ~ ., d, method = "BFGS")
  b <- countreg(Total_crashes ~ lnaadt + lnlength + speed50, d,
                method = "BFGS")
  expect_equal(coef(a), coef(b))
  expect_equal(predict(a, newdata = d[1:5, ]), predict(b, newdata = d[1:5, ]))
  # Also in an auxiliary formula
  d2 <- d[, c("Total_crashes", "lnaadt", "speed50")]
  m <- countreg(Total_crashes ~ lnaadt, d2,
                dis_param_formula_1 = ~ . - Total_crashes - lnaadt,
                method = "BFGS")
  expect_true("ln(alpha):speed50" %in% names(coef(m)))
})

test_that("auxiliary-formula contrasts are kept for prediction", {
  d <- washington_roads[, c("Total_crashes", "lnaadt", "speed50")]
  d$f <- factor(ifelse(d$speed50 == 1, "hi", "lo"))
  contrasts(d$f) <- contr.sum(2)
  m <- countreg(Total_crashes ~ lnaadt, d, dis_param_formula_1 = ~ f,
                method = "BFGS")
  expect_true("ln(alpha):f1" %in% names(coef(m)))
  nd <- d
  attr(nd$f, "contrasts") <- NULL
  expect_no_warning(p1 <- predict(m, newdata = nd))
  expect_equal(p1, predict(m, newdata = d))
  expect_equal(as.numeric(fitted(m)), p1)
})

test_that("regCompTable() without R^2 works when data are not kept", {
  a <- countreg(Total_crashes ~ lnaadt, washington_roads, method = "BFGS",
                keep_data = FALSE)
  b <- countreg(Total_crashes ~ lnaadt + speed50, washington_roads,
                method = "BFGS", keep_data = FALSE)
  tab <- regCompTable(list(a = a, b = b), RSquare = FALSE)
  expect_false("Pseudo-R-Sq." %in% tab$Parameter)
  expect_equal(as.numeric(tab$a[tab$Parameter == "AIC"]), round(AIC(a), 3))
})

test_that("tidy() components come from the model, not covariate names", {
  d <- data.frame(y = washington_roads$Total_crashes,
                  p = washington_roads$lnaadt,
                  gamma = washington_roads$lnlength)
  m <- countreg(y ~ p + gamma, d, method = "BFGS")
  td <- tidy(m, exponentiate = TRUE)
  expect_equal(td$component, c("mean", "mean", "mean", "distribution"))
  expect_equal(td$estimate[2:3], unname(exp(coef(m)[2:3])))
  expect_equal(td$estimate[4], unname(coef(m)[4]))
})

test_that("augment() residuals use the transformed response", {
  m <- countreg(I(Total_crashes + 1) ~ lnaadt, washington_roads,
                method = "BFGS")
  nd <- washington_roads[1:4, ]
  a <- augment(m, new_data = nd)
  expect_equal(a$.resid, nd$Total_crashes + 1 - a$.fitted)
  expect_equal(augment(m)$.resid, as.numeric(residuals(m)))
  expect_error(augment(m, new_data = nd, type.residuals = "pearson"),
               "response residuals")
})
