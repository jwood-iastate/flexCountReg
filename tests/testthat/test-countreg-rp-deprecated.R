# countreg.rp() is the deprecated former name of countreg_rp().

test_that("countreg.rp() warns once and gives the countreg_rp() fit", {
  data("washington_roads", package = "flexCountReg")
  d <- washington_roads[1:300, ]
  new <- countreg_rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                     data = d, family = "POISSON", ndraws = 50,
                     method = "BFGS")
  rlang::reset_warning_verbosity("fcr_countreg.rp")
  expect_warning(
    old <- countreg.rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                       data = d, family = "POISSON", ndraws = 50,
                       method = "BFGS"),
    "countreg.rp() is deprecated; use countreg_rp()", fixed = TRUE,
    class = "lifecycle_warning_deprecated")
  expect_equal(coef(old), coef(new))
  # the stored call has the caller's own expressions, so refits work
  expect_identical(old$call$data, quote(d))
  expect_identical(old$call$family, "POISSON")
})
