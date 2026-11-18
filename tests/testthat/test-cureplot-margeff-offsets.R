# cureplot() without cureplots' internal-name messages; margEffTable()
# without offset-only variables.

offset_data <- function(n = 200, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), e = stats::runif(n, 0.5, 2))
  d$le <- log(d$e)
  d$y <- stats::rpois(n, d$e * exp(0.3 + 0.4 * d$x))
  d
}

test_that("margEffTable() leaves out variables used only in offset()", {
  d <- offset_data()
  m <- countreg(y ~ x + offset(log(e)), data = d, family = "POISSON",
                method = "BFGS")
  expect_identical(margEffTable(m)$term, "x")
  m2 <- countreg(y ~ x + offset(le), data = d, family = "POISSON",
                 method = "BFGS")
  expect_identical(margEffTable(m2)$term, "x")
  # A variable that is also a predictor stays a candidate
  m3 <- countreg(y ~ x + log(e) + offset(log(e)), data = d,
                 family = "POISSON", method = "BFGS")
  expect_setequal(margEffTable(m3)$term, c("x", "e"))
  # Offsets given with the `offset` argument are excluded as before
  m4 <- countreg(y ~ x, data = d, family = "POISSON", method = "BFGS",
                 offset = "le")
  expect_identical(margEffTable(m4)$term, "x")
  # They can still be requested explicitly
  expect_setequal(margEffTable(m, vars = c("x", "e"))$term, c("x", "e"))
  # The same rule for mixture models (a degenerate two-component fit of
  # Poisson data; only the list of variables matters here)
  mx <- suppressWarnings(countreg_mix(y ~ x + offset(log(e)), data = d,
                                      components = 2, n_starts = 1,
                                      seed = 1))
  expect_identical(suppressWarnings(margEffTable(mx))$term, "x")
})

test_that("cureplot() does not report cureplots' internal names", {
  d <- offset_data()
  m <- countreg(y ~ x + offset(log(e)), data = d, family = "POISSON",
                method = "BFGS")
  expect_no_message(p <- cureplot(m, indvar = "x"))
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$x, "x")
  expect_no_message(p2 <- cureplot(m))
  expect_identical(p2$labels$x, "Predicted Values")
  # The plotted data are the cumulative residuals in covariate order
  r <- d$y - predict(m, newdata = d)
  o <- order(d$x)
  expect_equal(p$data$cumres, cumsum(r[o]))
  # Other messages are not muffled
  expect_message(flexCountReg:::.fcr_muffle_cure_messages(message("other")),
                 "other")
})
