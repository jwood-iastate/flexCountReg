# Regression tests for the corrections from the development review of
# flexCountReg 0.1.3 (findings F01-F28).

test_that("F01: PGE and PW estimate only the shape parameter", {
  h <- randtoolbox::halton(200)
  # The scale of the mixing distribution cancels out of the mean-parameterized
  # PMF, so it cannot be estimated
  expect_equal(sapply(c(0.1, 1, 10, 100), function(s)
    dpge(3, mean = 2, shape = 1.5, scale = s, haltons = h)),
    rep(dpge(3, mean = 2, shape = 1.5, scale = 1, haltons = h), 4))
  expect_equal(sapply(c(0.1, 1, 10, 100), function(s)
    dpWeib_cpp(3, mean = 2, alpha = 1.5, sigma = s, h = h)),
    rep(dpWeib_cpp(3, mean = 2, alpha = 1.5, sigma = 1, h = h), 4))
  d <- washington_roads
  for (fam in c("PGE", "PW")) {
    m <- countreg(Total_crashes ~ lnaadt, d, family = fam, offset = "lnlength",
                  ndraws = 100, method = "BFGS")
    expect_length(coef(m), 3)
    expect_equal(attr(logLik(m), "df"), 3)
    expect_true(all(is.finite(sqrt(diag(vcov(m))))))
    expect_error(countreg(Total_crashes ~ lnaadt, d, family = fam,
                          dis_param_formula_2 = ~ speed50),
                 "cancels out of the likelihood")
    expect_error(countreg.rp(Total_crashes ~ lnaadt,
                             rpar_formula = ~ -1 + speed50, data = d,
                             family = fam, dis_param_formula_2 = ~ speed50),
                 "cancels out of the likelihood")
  }
})
