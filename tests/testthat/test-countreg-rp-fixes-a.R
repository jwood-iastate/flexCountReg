# Split by run time so the parts run in parallel (see Config/testthat/start-first).
# Regression tests for countreg_rp(): draw conventions, rp_mgf(), Exact vs
# Simulated predictions, heterogeneity in means, input validation and
# compatibility with fits from earlier versions.

rp_data <- function() {
  d <- washington_roads
  d$x <- d$lnaadt - mean(d$lnaadt)
  d
}

test_that("Exact and Simulated agree with heterogeneity, for each distribution", {
  skip_on_cran()  # slow (stochastic recovery); run on CI and locally
  d <- rp_data()
  for (dist in c("n", "u", "t", "g")) {
    m <- suppressWarnings(countreg_rp(
      Total_crashes ~ lnlength, rpar_formula = ~ -1 + x, data = d,
      rpardists = c(x = dist), ndraws = 100, method = "BFGS",
      het_mean_formula = ~ speed50, het_var_formula = ~ ShouldWidth04))
    m$model$ndraws <- 5000
    ps <- predict(m, newdata = d, method = "Simulated")
    pe <- predict(m, newdata = d, method = "Exact")
    expect_lt(max(abs(ps / pe - 1)), 0.01)
  }
})
