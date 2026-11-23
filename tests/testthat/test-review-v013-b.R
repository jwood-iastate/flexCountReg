# Split by run time so the parts run in parallel (see Config/testthat/start-first).
# Regression tests for the corrections from the development review of
# flexCountReg 0.1.2 (findings F01-F28).

test_that("F05: two-model LR tests check nesting and the sample", {
  d <- washington_roads
  nb1 <- countreg(Total_crashes ~ lnaadt, d, family = "NB1", method = "BFGS")
  nbp <- countreg(Total_crashes ~ lnaadt, d, family = "NBP", method = "BFGS")
  pl <- countreg(Total_crashes ~ lnaadt, d, family = "PL", method = "BFGS")
  plg <- countreg(Total_crashes ~ lnaadt, d, family = "PLG", method = "BFGS")
  # NB1 is NBP with P = 1 (interior)
  r <- regCompTest(nb1, model2 = nbp)
  expect_equal(r$LR_type, "interior")
  expect_equal(r$LR_pvalue, pchisq(r$LR, 1, lower.tail = FALSE))
  # PL is PLG with alpha = 0 (boundary)
  expect_equal(.fcr_nesting(.fcr_model_spec(pl), .fcr_model_spec(plg))$type,
               "boundary")
  # Poisson within NB2, two-model form: boundary test
  po <- countreg(Total_crashes ~ lnaadt, d, family = "POISSON",
                 method = "BFGS")
  nb2 <- countreg(Total_crashes ~ lnaadt, d, method = "BFGS")
  r <- regCompTest(nb2, model2 = po)
  expect_equal(r$LR_type, "boundary")
  expect_equal(r$LR_pvalue, 0.5 * pchisq(r$LR, 1, lower.tail = FALSE))
  # PL and NB1 have the same number of parameters; PL vs NBP is not nested
  expect_warning(r <- regCompTest(pl, model2 = nbp),
                 "not nested|not a subset")
  expect_true(is.na(r$LR_pvalue))
  # Same number of rows but different observations
  d2 <- d
  d2$Total_crashes <- rev(d2$Total_crashes)
  other <- countreg(Total_crashes ~ lnaadt, d2, family = "NBP",
                    method = "BFGS")
  expect_error(regCompTest(nb1, model2 = other), "same observations")
})
