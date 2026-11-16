# Input checks of the shared data preparation, the bootstrap helper, and the
# less common branches of the standard methods.

test_that("prepare_model_data() checks offsets, weights and variables", {
  d <- washington_roads
  f <- Total_crashes ~ lnaadt
  expect_error(countreg(f, d, offset = 1), "character vector of column names")
  expect_error(countreg(f, d, offset = "nope"),
               "Offset column\\(s\\) not found in `data`: nope")
  expect_error(countreg(f, d, weights = c("a", "b")),
               "single column name or a numeric vector")
  expect_error(countreg(f, d, weights = 1:3), "one value per row")
  expect_error(countreg(f, d, weights = list(1)),
               "column name or a numeric vector")
  expect_error(countreg(Total_crashes ~ lnaadt + nope, d),
               "Variable\\(s\\) not found in `data`: nope")
  d_na <- d
  d_na$lnaadt <- NA
  expect_error(countreg(f, d_na), "No complete cases")
  d_inf <- d
  d_inf$len0 <- log(d$Length - d$Length)   # log(0) = -Inf
  expect_error(countreg(f, d_inf, offset = "len0"), "Offsets must be finite")
  expect_error(countreg(f, d, weights = rep(-1, nrow(d))), "non-negative")
  # Missing numeric weights drop their rows
  w <- rep(1, nrow(d))
  w[1:3] <- NA
  m <- countreg(f, d, weights = w, method = "BFGS")
  expect_equal(nobs(m), nrow(d) - 3)
  expect_error(countreg(f, d, start.vals = c("a", "b", "c")),
               "`start.vals` must be numeric")
})

test_that("the bootstrap helper drops failed replicates", {
  set.seed(1)
  calls <- 0
  refit <- function(idx) {
    calls <<- calls + 1
    if (calls %% 2 == 0) stop("no convergence")
    c(a = mean(idx), b = 1)
  }
  expect_warning(dr <- .fcr_bootstrap(10, 10, c("a", "b"), refit),
                 "5 of 10 bootstrap replicates failed")
  expect_equal(sum(stats::complete.cases(dr)), 5)
  # No usable replicate: a warning, and NA standard errors
  expect_warning(dr <- .fcr_bootstrap(10, 5, c("a", "b"),
                                      function(idx) stop("x")),
                 "Fewer than 2 bootstrap replicates succeeded \\(error: 5\\)")
  expect_true(all(is.na(.fcr_boot_cov(dr))))
  # Replicates missing an estimate count as failures
  expect_warning(.fcr_bootstrap(10, 5, c("a", "b"), function(idx) c(a = 1)),
                 "non-finite estimates: 5")
})

test_that("verbose = TRUE prints the optimizer trace", {
  expect_output(countreg(Total_crashes ~ lnaadt, washington_roads,
                         method = "BFGS", verbose = TRUE))
})

test_that("standard methods handle objects with fewer stored parts", {
  d <- washington_roads
  m <- countreg(Total_crashes ~ lnaadt + lnlength, d, method = "BFGS")
  # terms() without a blueprint falls back to the formula
  m2 <- m
  m2$model$blueprint <- NULL
  expect_equal(attr(terms(m2), "term.labels"), c("lnaadt", "lnlength"))
  # confint() for chosen parameters
  ci <- confint(m, parm = "lnaadt", level = 0.9)
  expect_equal(dim(ci), c(1L, 2L))
  expect_equal(colnames(ci), c("5 %", "95 %"))
  # residuals() recover the response from the data when y is not stored
  m3 <- m
  m3$model$y <- NULL
  expect_equal(residuals(m3), residuals(m))
  # augment() with the full training data uses only the rows in the fit
  d_na <- d
  d_na$lnaadt[1:4] <- NA
  mna <- suppressMessages(countreg(Total_crashes ~ lnaadt, d_na,
                                   method = "BFGS"))
  a <- augment(mna, data = d_na)
  expect_equal(nrow(a), nrow(d) - 4)
  expect_equal(a$.fitted, fitted(mna))
  # tidy() and summary() when no covariance matrix is available
  m4 <- m
  m4$model$vcov <- NULL
  expect_true(all(is.na(tidy(m4)$std.error)))
  expect_true(all(is.na(summary(m4)$coefficients[, "Std. Error"])))
})

test_that("Pearson residuals and simulate() are for countreg() models", {
  re <- renb(Total_crashes ~ lnaadt, washington_roads, group_var = "ID",
             method = "BFGS")
  expect_error(residuals(re, type = "pearson"), "available for countreg")
  expect_error(simulate(re), "available for countreg")
  out <- capture.output(print(re))
  expect_true(any(grepl("panels =", out)))
  s <- capture.output(print(summary(re)))
  expect_true(any(grepl("panels:", s)))
})

test_that("printing reports convergence problems, dropped rows and the call", {
  d <- washington_roads
  d$lnaadt[1:2] <- NA
  m <- suppressMessages(countreg(Total_crashes ~ lnaadt, d, method = "BFGS"))
  s <- capture.output(print(summary(m)))
  expect_true(any(grepl("rows removed \\(missing values\\): 2", s)))
  m$model$converged <- FALSE
  expect_output(print(m), "did not converge")
  expect_output(print(summary(m)), "did not report convergence")
  m$call <- NULL
  expect_output(print(summary(m)), "Total_crashes ~ lnaadt")
})

test_that("the LR test against Poisson is skipped when it does not apply", {
  d <- washington_roads
  pois <- countreg(Total_crashes ~ lnaadt, d, family = "POISSON",
                   method = "BFGS")
  expect_null(summary(pois)$lr_test)
  nb <- countreg(Total_crashes ~ lnaadt, d, method = "BFGS", keep_data = FALSE)
  expect_null(summary(nb)$lr_test)
  re <- renb(Total_crashes ~ lnaadt, d, group_var = "ID", method = "BFGS")
  expect_null(summary(re)$lr_test)
  nb <- countreg(Total_crashes ~ lnaadt, d, method = "BFGS")
  lr <- summary(nb)$lr_test
  expect_true(lr$boundary)
  # Closed-form check of the statistic
  p <- glm(Total_crashes ~ lnaadt, family = poisson, data = d)
  expect_equal(lr$statistic, 2 * (nb$model$maximum - as.numeric(logLik(p))),
               tolerance = 1e-6)
})
