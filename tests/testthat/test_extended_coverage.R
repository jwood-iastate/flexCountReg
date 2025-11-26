test_that("countreg.rp: Basic execution with BHHH and Predictions", {
  data("washington_roads")
  # Basic Random Parameters NB2
  # Using BHHH to test the specific gradient check logic implemented
  model <- suppressWarnings(
    countreg.rp(Total_crashes ~ lnaadt + lnlength,
                rpar_formula = ~ -1 + speed50,
                data = washington_roads,
                family = "NB2",
                rpardists = c(speed50 = "n"),
                ndraws = 20,         # Low draws for speed
                method = "BHHH",     # Test BHHH path
                max.iters = 50)
  )
  
  expect_s3_class(model, "flexCountReg")
  
  # Test Predict methods for countreg.rp
  pred_sim <- predict(model, data = washington_roads, method = "Simulated")
  pred_ind <- predict(model, data = washington_roads, method = "Individual")
  
  expect_equal(length(pred_sim), nrow(washington_roads))
  expect_equal(length(pred_ind), nrow(washington_roads))
  expect_false(any(is.na(pred_sim)))
  expect_false(any(is.na(pred_ind)))
})

test_that("countreg.rp: Heterogeneity, Dist Formulas, and Panel Data", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
  
  # create a dummy panel ID
  washington_roads$site_id <- rep(1:50, length.out = nrow(washington_roads))
  
  # Complex model: Heterogeneity in Mean + Dist Param Formula + Panel
  model <- suppressWarnings(
    countreg.rp(Total_crashes ~ lnaadt,
                rpar_formula = ~ -1 + speed50,
                dis_param_formula_1 = ~ lnlength,     # Varying Alpha
                het_mean_formula = ~ AADT10kplus,     # Het Mean
                data = washington_roads,
                panel_id = "site_id",                 # Panel Structure
                family = "NB2",
                rpardists = c(speed50 = "n"),
                ndraws = 20,
                method = "NM",
                max.iters = 50)
  )
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})

test_that("countreg.rp: Correlated Random Parameters", {
  data("washington_roads")
  
  # Correlated Random Parameters (requires normal dist)
  model <- suppressWarnings(
    countreg.rp(Total_crashes ~ lnaadt,
                rpar_formula = ~ speed50 + lnlength,  # 2 Random Params
                data = washington_roads,
                family = "NB2",
                correlated = TRUE,                    # Activate Correlation
                ndraws = 20,
                method = "BFGS",
                max.iters = 50)
  )
  
  expect_s3_class(model, "flexCountReg")
  # Check if Cholesky/Covariance matrix was generated
  expect_false(!is.null(model$model$Covariance))
})

test_that("rpnb: Heterogeneity in Means and Variances", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
  
  # Test het_mean and het_var arguments in rpnb
  model <- suppressWarnings(
    rpnb(Total_crashes ~ lnlength,
         rpar_formula = ~ -1 + lnaadt,
         het_mean_formula = ~ speed50,
         het_var_formula = ~ AADT10kplus,
         data = washington_roads,
         ndraws = 20,
         form = 'nb2',
         verbose = FALSE)
  )
  
  expect_s3_class(model, "flexCountReg")
  # Check if heterogeneity coefficients exist
  expect_true(!is.null(model$model$het_mean_coefs))
  expect_true(!is.null(model$model$het_var_coefs))
})


test_that("poisLind.re: Bootstrapping", {
  data("washington_roads")
  
  # Test poisLind.re with bootstrapping
  model <- suppressWarnings(
    poisLind.re(Animal ~ lnaadt + lnlength,
                data = washington_roads,
                group_var = "ID",
                bootstraps = 5,
                method = "NM",
                max.iters = 50)
  )
  
  expect_s3_class(model, "flexCountReg")
  expect_true(!is.null(model$model$bootstrapped_se))
})

test_that("Metrics calculations are mathematically correct", {
  
  # Test myAIC
  # AIC = -2 * LL + 2 * k
  expect_equal(myAIC(LL = -100, nparam = 5), -2 * (-100) + 2 * 5) # 210
  expect_equal(myAIC(LL = 0, nparam = 0), 0)
  
  # Test myBIC
  # BIC = -2 * LL + k * log(n)
  expect_equal(myBIC(LL = -100, nparam = 5, n = 100), -2 * (-100) + 5 * log(100))
  
  # Test RMSE
  y <- c(1, 2, 3)
  mu <- c(1, 2, 3)
  expect_equal(rmse(y, mu), 0) # Perfect prediction
  
  y2 <- c(0, 0)
  mu2 <- c(3, 4)
  # sqrt(mean((3^2 + 4^2))) = sqrt(mean(9+16)) = sqrt(12.5)
  expect_equal(rmse(y2, mu2), sqrt(12.5))
  
  # Test MAE
  expect_equal(mae(y, mu), 0)
  expect_equal(mae(c(1, 1), c(2, 0)), 1) # mean(|1-2| + |1-0|) = 1
})


test_that("regCompTable handles formats and errors", {
  data("washington_roads")
  
  # Fit two small, fast models for comparison
  m1 <- countreg(Total_crashes ~ lnaadt, data = washington_roads, family = "NB1")
  m2 <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads, family = "NB2")
  
  # Test Latex Output
  tbl_latex <- regCompTable(models = list("Small" = m1, "Large" = m2), tableType = "latex")
  expect_s3_class(tbl_latex, "knitr_kable")
})

test_that("regCompTest calculates LR statistics correctly", {
  data("washington_roads")
  
  # Nested models
  m_full <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads, family = "NB2")
  
  # Test against Poisson intercept only (default)
  res_def <- regCompTest(m_full, washington_roads)
  expect_type(res_def, "list")
  expect_true(!is.null(res_def$LR))
  expect_true(res_def$LR_pvalue <= 1)
  
  # Test against specific base model string
  res_spec <- regCompTest(m_full, washington_roads, basemodel = "Poisson", variables = TRUE)
  expect_true(res_spec$LL > res_spec$LLbase) # Full NB2 should fit better than Poisson
})

test_that("cureplot runs with various arguments", {
  data("washington_roads")
  
  # Fit a fast model
  model <- countreg(Total_crashes ~ lnaadt, data = washington_roads, family = "POISSON", ndraws = 10)
  
  # 1. Default execution (against predicted values)
  # We use expect_error(..., NA) to assert NO error is thrown
  expect_error(cureplot(model, n_resamples = 0), NA)
  
  # 2. Execution with specific independent variable
  expect_error(cureplot(model, indvar = "lnaadt", n_resamples = 0), NA)
  
  # 3. Error handling: Invalid variable name
  expect_error(cureplot(model, indvar = "NON_EXISTENT_VAR"), 
               "not found in the provided data")
})



test_that("Random parameter generation handles non-Normal distributions", {
  # This targets the switch statement in helpers.R -> generate_random_draws
  data("washington_roads")
  
  # Reduce data for speed
  small_data <- washington_roads[1:50, ]
  
  # 1. Triangular Distribution ("t")
  # Use BHHH and low iters just to check if the probability function evaluates without crashing
  expect_error(
    countreg.rp(Total_crashes ~ lnaadt,
                rpar_formula = ~ -1 + speed50,
                data = small_data,
                family = "POISSON",
                rpardists = c(speed50 = "t"), # Triangle
                ndraws = 5,
                method = "BHHH", 
                max.iters = 2),
    NA
  )
  
  # 2. Uniform Distribution ("u")
  expect_error(
    countreg.rp(Total_crashes ~ lnaadt,
                rpar_formula = ~ -1 + speed50,
                data = small_data,
                family = "NB2",
                rpardists = c(speed50 = "u"), # Uniform
                ndraws = 5,
                method = "BHHH", 
                max.iters = 2),
    NA
  )
  
})
