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

# Second Copilot review (PR #5) -------------------------------------------------

test_that("fitted objects keep the expanded `.` formulas", {
  d <- washington_roads[, c("Total_crashes", "lnaadt", "lnlength")]
  m <- countreg(Total_crashes ~ ., d, method = "BFGS")
  expect_equal(attr(terms(m$formula), "term.labels"), c("lnaadt", "lnlength"))
  # Methods that re-read the formula without data work
  expect_equal(margEffTable(m)$term, c("lnaadt", "lnlength"))
  # Bootstrap refits of countreg.rp() must not expand `.` against the
  # bootstrap panel column (it would add a coefficient)
  set.seed(1)
  r <- countreg.rp(Total_crashes ~ ., rpar_formula = ~ -1 + lnlength,
                   data = d, ndraws = 20, method = "BFGS", stderr = "boot",
                   bootstraps = 3, max.iters = 50)
  expect_equal(dim(vcov(r)), rep(length(coef(r)), 2))
  expect_equal(r$model$n_boot_ok, 3)
})

test_that("factor responses are rejected instead of fit as level codes", {
  d <- washington_roads
  d$Total_crashes <- factor(d$Total_crashes)
  msg <- "response must be a numeric count vector, not factor"
  expect_error(countreg(Total_crashes ~ lnaadt, d), msg)
  expect_error(repois(Total_crashes ~ lnaadt, d, group_var = "ID"), msg)
  expect_error(renb(Total_crashes ~ lnaadt, d, group_var = "ID"), msg)
  expect_error(negmultinom(Total_crashes ~ lnaadt, "ID", d), msg)
})

test_that("COM upper tails are summed directly (no negative probabilities)", {
  # nu = 1 is the Poisson distribution
  up <- pcom(c(5, 25, 40), lambda = 1, nu = 1, lower.tail = FALSE)
  expect_true(all(up > 0))
  expect_equal(up, ppois(c(5, 25, 40), 1, lower.tail = FALSE),
               tolerance = 1e-10)
  expect_equal(pcom(0:40, lambda = 4, nu = 1, lower.tail = FALSE,
                    log.p = TRUE),
               ppois(0:40, 4, lower.tail = FALSE, log.p = TRUE),
               tolerance = 1e-10)
  lo <- pcom(0:30, lambda = 3, nu = 0.7)
  expect_lte(max(lo), 1)
  expect_equal(lo + pcom(0:30, lambda = 3, nu = 0.7, lower.tail = FALSE),
               rep(1, 31), tolerance = 1e-12)
})

test_that("missing parameters give NA quantiles, not an error", {
  expect_identical(is.na(qplind(c(0.5, 0.5), mean = c(NA, 2))),
                   c(TRUE, FALSE))
  expect_true(is.na(qplindLnorm(0.5, mean = NA_real_, theta = 1.5,
                                sigma = 0.5, ndraws = 100)))
})

test_that("sigma = 0 is the Poisson-Lindley limit in qplindLnorm()", {
  p <- c(0.1, 0.5, 0.9)
  expect_equal(qplindLnorm(p, mean = 2, theta = 1.5, sigma = 0),
               qplind(p, mean = 2, theta = 1.5))
  expect_no_warning(rplindLnorm(5, mean = 2, theta = 1.5, sigma = 0))
  expect_warning(qplindLnorm(0.5, mean = 2, theta = 1.5, sigma = -1),
                 "non-negative")
})

test_that("method = \"Nelder-Mead\" is accepted as an alias of \"NM\"", {
  expect_identical(.fcr_check_method("Nelder-Mead"), "NM")
  m <- repois(Total_crashes ~ lnaadt, washington_roads, group_var = "ID",
              method = "Nelder-Mead")
  expect_identical(m$model$method, "NM")
})

test_that("augment() evaluates response transformations where the model was fit", {
  fit_local <- function() {
    tf <- function(y) 2 * y
    countreg(tf(Total_crashes) ~ lnaadt, washington_roads, method = "BFGS")
  }
  m <- fit_local()
  nd <- washington_roads[1:4, ]
  a <- augment(m, new_data = nd)
  expect_equal(a$.resid, 2 * nd$Total_crashes - a$.fitted)
})
