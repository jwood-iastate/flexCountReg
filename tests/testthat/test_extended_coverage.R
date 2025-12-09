test_that("countreg.rp: Basic execution with BHHH and Predictions", {
  data("washington_roads")
  # Use smaller subset for speed
  small_data <- washington_roads[1:300, ]
  
  model <- suppressWarnings(
    countreg.rp(Total_crashes ~ lnaadt + lnlength,
                rpar_formula = ~ -1 + speed50,
                data = small_data,  # Use smaller data
                family = "NB2",
                rpardists = c(speed50 = "n"),
                ndraws = 15,         # Reduced from 20
                method = "BHHH",
                max.iters = 30)      # Reduced from 50
  )
  
  expect_s3_class(model, "flexCountReg")
  
  # Test Predict methods for countreg.rp
  pred_sim <- predict(model, data = small_data, method = "Simulated")
  pred_ind <- predict(model, data = small_data, method = "Individual")
  
  expect_equal(length(pred_sim), nrow(small_data))
  expect_equal(length(pred_ind), nrow(small_data))
  expect_false(any(is.na(pred_sim)))
  expect_false(any(is.na(pred_ind)))
})

test_that("countreg.rp: Heterogeneity, Dist Formulas, and Panel Data", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
  
  # Use smaller subset
  small_data <- washington_roads[1:200, ]
  small_data$site_id <- rep(1:20, length.out = nrow(small_data))
  
  model <- suppressWarnings(
    countreg.rp(Total_crashes ~ lnaadt,
                rpar_formula = ~ -1 + speed50,
                dis_param_formula_1 = ~ lnlength,
                het_mean_formula = ~ AADT10kplus,
                data = small_data,
                panel_id = "site_id",
                family = "NB2",
                rpardists = c(speed50 = "n"),
                ndraws = 15,       # Reduced
                method = "NM",
                max.iters = 30)    # Reduced
  )
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})

test_that("countreg.rp: Correlated Random Parameters", {
  data("washington_roads")
  small_data <- washington_roads[1:200, ]
  
  model <- suppressWarnings(
    countreg.rp(Total_crashes ~ lnaadt,
                rpar_formula = ~ speed50 + lnlength,
                data = small_data,
                family = "NB2",
                correlated = TRUE,
                ndraws = 15,       # Reduced
                method = "BFGS",
                max.iters = 30)    # Reduced
  )
  
  expect_s3_class(model, "flexCountReg")
  expect_false(!is.null(model$model$Covariance))
})

test_that("poisLind.re: Bootstrapping", {
  skip_on_cran()  # Skip on CRAN - bootstrapping is slow
  
  data("washington_roads")
  small_data <- washington_roads[1:300, ]
  
  model <- suppressWarnings(
    poisLind.re(Animal ~ lnaadt + lnlength,
                data = small_data,
                group_var = "ID",
                bootstraps = 3,     # Reduced from 5
                method = "NM",
                max.iters = 30)     # Reduced
  )
  
  expect_s3_class(model, "flexCountReg")
  expect_true(!is.null(model$model$bootstrapped_se))
})

test_that("Metrics calculations are mathematically correct", {
  # Test myAIC
  expect_equal(myAIC(LL = -100, nparam = 5), -2 * (-100) + 2 * 5)
  expect_equal(myAIC(LL = 0, nparam = 0), 0)
  
  # Test myBIC
  expect_equal(myBIC(
    LL = -100, nparam = 5, n = 100), -2 * (-100) + 5 * log(100))
  
  # Test RMSE
  y <- c(1, 2, 3)
  mu <- c(1, 2, 3)
  expect_equal(rmse(y, mu), 0)
  
  y2 <- c(0, 0)
  mu2 <- c(3, 4)
  expect_equal(rmse(y2, mu2), sqrt(12.5))
  
  # Test MAE
  expect_equal(mae(y, mu), 0)
  expect_equal(mae(c(1, 1), c(2, 0)), 1)
})


test_that("regCompTable handles formats and errors", {
  data("washington_roads")
  small_data <- washington_roads[1:300, ]
  
  m1 <- countreg(
    Total_crashes ~ lnaadt, data = small_data, family = "NB1", max.iters = 100)
  m2 <- countreg(
    Total_crashes ~ lnaadt + lnlength, 
    data = small_data, 
    family = "NB2", 
    max.iters = 100)
  
  tbl_latex <- regCompTable(
    models = list("Small" = m1, "Large" = m2), tableType = "latex")
  expect_s3_class(tbl_latex, "knitr_kable")
})

test_that("regCompTest calculates LR statistics correctly", {
  data("washington_roads")
  small_data <- washington_roads[1:300, ]
  
  m_full <- countreg(
    Total_crashes ~ lnaadt + lnlength, 
    data = small_data, 
    family = "NB2", 
    max.iters = 100)
  
  res_def <- regCompTest(m_full, small_data)
  expect_type(res_def, "list")
  expect_true(!is.null(res_def$LR))
  expect_true(res_def$LR_pvalue <= 1)
  
  res_spec <- regCompTest(
    m_full, small_data, basemodel = "Poisson", variables = TRUE)
  expect_true(res_spec$LL > res_spec$LLbase)
})

test_that("cureplot runs with various arguments", {
  data("washington_roads")
  small_data <- washington_roads[1:200, ]
  
  model <- countreg(
    Total_crashes ~ lnaadt, 
    data = small_data, 
    family = "POISSON", 
    max.iters = 50)
  
  expect_error(cureplot(model, n_resamples = 0), NA)
  expect_error(cureplot(model, indvar = "lnaadt", n_resamples = 0), NA)
})


test_that("Random parameter generation handles non-Normal distributions", {
  data("washington_roads")
  small_data <- washington_roads[1:30, ]  # Very small for distribution tests
  
  # 1. Triangular Distribution ("t")
  expect_error(
    countreg.rp(Total_crashes ~ lnaadt,
                rpar_formula = ~ -1 + speed50,
                data = small_data,
                family = "POISSON",
                rpardists = c(speed50 = "t"),
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
                rpardists = c(speed50 = "u"),
                ndraws = 5,
                method = "BHHH", 
                max.iters = 2),
    NA
  )
})

###################################
## dsichel and related - STREAMLINED
###################################

test_that("dsichel handles normal cases", {
  p <- dsichel(0:10, mu = 5, sigma = 1, gamma = -0.5)
  expect_true(all(is.finite(p)))
  expect_true(all(p >= 0))
  expect_true(all(p <= 1))
  # Reduced range for sum check
  expect_equal(
    sum(dsichel(0:100, mu = 5, sigma = 1, gamma = -0.5)), 1, tolerance = 0.01)
})

test_that("dsichel handles extreme parameters", {
  # Combined extreme tests for speed
  p1 <- dsichel(0:5, mu = 5, sigma = 0.01, gamma = -0.5)
  p2 <- dsichel(0:5, mu = 5, sigma = 100, gamma = -0.5)
  p3 <- dsichel(0:5, mu = 5, sigma = 1, gamma = 50)
  p4 <- dsichel(0:5, mu = 5, sigma = 1, gamma = -50)
  
  expect_true(all(is.finite(p1)))
  expect_true(all(is.finite(p2)))
  expect_true(all(is.finite(p3)))
  expect_true(all(is.finite(p4)))
})

test_that("dsichel log option works correctly", {
  p <- dsichel(5, mu = 5, sigma = 1, gamma = -0.5)
  log_p <- dsichel(5, mu = 5, sigma = 1, gamma = -0.5, log = TRUE)
  expect_equal(log(p), log_p, tolerance = 1e-10)
})

test_that("psichel CDF is monotonic", {
  cdf <- psichel(0:15, mu = 5, sigma = 1, gamma = -0.5)  # Reduced from 0:20
  expect_true(all(diff(cdf) >= 0))
  expect_true(all(cdf >= 0 & cdf <= 1))
})

test_that("qsichel is inverse of psichel", {
  p_vals <- c(0.25, 0.5, 0.75)  # Reduced test points
  q_vals <- qsichel(p_vals, mu = 5, sigma = 1, gamma = -0.5)
  p_back <- psichel(q_vals, mu = 5, sigma = 1, gamma = -0.5)
  
  expect_true(all(p_back >= p_vals - 0.01))
})
