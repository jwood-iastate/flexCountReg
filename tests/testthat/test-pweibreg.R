test_that("Poisson-Weibull model runs and returns correct output", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- pwiebreg(Total_crashes ~ lnlength + lnaadt,
                    data = washington_roads,
                    method="NM",
                    ndraws = 100)
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
  expect_named(model$model$estimate, c("lnlength", "lnaadt", "ln(alpha)", "ln(sigma)"))  # Check names of estimates
  expect_true(model$model$LL < 0)  # Log-likelihood should be negative
})

# Test the Poisson-Weibull model with alpha_formula
test_that("Poisson-Weibull model runs correctly with alpha_formula", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- pwiebreg(Total_crashes ~ lnlength + lnaadt,
                    alpha_formula = ~ lnaadt,
                    data = washington_roads,
                    ndraws = 100)
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
  expect_named(model$model$estimate, c("lnlength", "lnaadt", "ln(alpha):(Intercept)", "ln(alpha):lnaadt", "ln(sigma)"))
})

# Test the Poisson-Weibull model with sigma_formula
test_that("Poisson-Weibull model runs correctly with sigma_formula", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- pwiebreg(Total_crashes ~ lnlength + lnaadt,
                    sigma_formula = ~ speed50,
                    data = washington_roads,
                    ndraws = 100)
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
  expect_named(model$model$estimate, c("(Intercept)","lnlength", "lnaadt", "ln(alpha)", "ln(sigma):(Intercept)", "ln(sigma):speed50"))
})

# Test the Poisson-Weibull model with both alpha_formula and sigma_formula
test_that("Poisson-Weibull model runs correctly with both alpha_formula and sigma_formula", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- pwiebreg(Total_crashes ~ lnlength + lnaadt,
                    alpha_formula = ~ lnaadt,
                    sigma_formula = ~ speed50,
                    data = washington_roads,
                    ndraws = 100)
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
  expect_named(model$model$estimate, c("(Intercept)","lnlength", "lnaadt", "ln(alpha):(Intercept)", "ln(alpha):lnaadt", "ln(sigma):(Intercept)", "ln(sigma):speed50"))
})

# Test the Poisson-Weibull model with bootstrapping
test_that("Poisson-Weibull model runs correctly with bootstrapping", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- pwiebreg(Total_crashes ~ lnlength + lnaadt,
                    data = washington_roads,
                    ndraws = 100,
                    bootstraps = 10)
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
  expect_named(model$model$bootstrapped_se, names(model$estimate))  # Check bootstrapped SEs are returned
  expect_true(all(model$bootstrapped_se > 0))  # SEs should be positive
})
