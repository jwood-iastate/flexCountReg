test_that("NB-2 model runs and returns correct output", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k,
               data = washington_roads, form = 'nb2', method = 'BHHH', max.iters = 3000)
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
  expect_named(model$model$estimate, c("(Intercept)", "lnaadt", "lnlength", "speed50", "ShouldWidth04", "AADTover10k", "ln(alpha)"))  # Check names of estimates
  expect_true(model$model$LL < 0)  # Log-likelihood should be negative
})

# Test the generalized NB-2 model with ln.alpha.formula
test_that("Generalized NB-2 model runs correctly with ln.alpha.formula", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k,
               data = washington_roads, form = 'nb2', method = 'BHHH', max.iters = 3000,
               ln.alpha.formula = ~ 1 + lnlength)
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
  expect_named(model$model$estimate, c("(Intercept)", "lnaadt", "lnlength", "speed50", "ShouldWidth04", "AADTover10k", "ln(alpha): (Intercept)", "ln(alpha): lnlength"))
})

# Test the NB-P model
test_that("NB-P model runs and returns correct output", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k,
               data = washington_roads, form = 'nbp', method = 'BHHH', max.iters = 3000)
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
  expect_named(model$model$estimate, c("(Intercept)", "lnaadt", "lnlength", "speed50", "ShouldWidth04", "AADTover10k", "ln(alpha)", "P"))
})

# Test the NB-1 model with bootstrapping
test_that("NB-1 model runs correctly with bootstrapping", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k,
               data = washington_roads, form = 'nb1', method = 'NM', max.iters = 3000, bootstraps = 10)
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
  expect_named(model$model$bootstrapped_se, names(model$model$estimate))  # Check bootstrapped SEs are returned
  expect_true(all(model$model$bootstrapped_se > 0))  # SEs should be positive
})

# Test error handling for invalid form
test_that("Error is thrown for invalid form input", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  expect_error(nbg(Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k,
                   data = washington_roads, form = 'invalid', method = 'BHHH', max.iters = 3000),
               "Invalid form specified") 
})

# Test that weights parameter works correctly
test_that("Weighted regression runs correctly", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  washington_roads$weights <- runif(nrow(washington_roads), min=0.5, max=1.5)  # Add random weights
  model <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k,
               data = washington_roads, form = 'nb2', method = 'BHHH', max.iters = 3000, weights = "weights")
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
  expect_named(model$model$estimate, c("(Intercept)", "lnaadt", "lnlength", "speed50", "ShouldWidth04", "AADTover10k", "ln(alpha)"))
})

# Test that the model runs with a different optimization method
test_that("Model runs correctly with different optimization method", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k,
               data = washington_roads, form = 'nb2', method = 'NM', max.iters = 3000)  # Use Newton-Raphson method
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})

# Test the NB-1 model with generalized options
test_that("Generalized NB-1 model runs correctly with ln.alpha.formula", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k,
               data = washington_roads, form = 'nb1', method = 'BHHH', max.iters = 3000,
               ln.alpha.formula = ~ lnlength)
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
  expect_named(model$model$estimate, c("(Intercept)", "lnaadt", "lnlength", "speed50", "ShouldWidth04", "AADTover10k", "ln(alpha): (Intercept)", "ln(alpha): lnlength"))
})
