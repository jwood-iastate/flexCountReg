test_that("rpnb Tests", {
  
  data("washington_roads")
  
  # Optimization may produce "NaNs produced" warnings during intermediate steps.
  # This is expected behavior for maxLik with complex likelihoods.
  # We suppress warnings to ensure the check passes if the model fits successfully.
  expect_error(
    suppressWarnings(
      nbp.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
                     rpar_formula = ~ speed50,
                     data = washington_roads,
                     ndraws = 100,
                     correlated = FALSE,
                     rpardists = c(intercept="n", speed50="n"),
                     form = 'nbp',
                     method = "bfgs")
    ), 
    NA
  )
  
  # Ensure the object is actually returned and is of the correct class
  expect_true(exists("nbp.rp"))
  expect_s3_class(nbp.rp, "flexCountReg")
  
})