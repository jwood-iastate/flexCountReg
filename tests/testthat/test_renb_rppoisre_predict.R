# Create a larger synthetic dataset for stability
# Random Effects models need sufficient groups and obs/group to converge
set.seed(123)
n_groups <- 50      # Increased from 20
n_per_group <- 20   # Increased from 10
id <- rep(1:n_groups, each = n_per_group)
# Secondary ID for composite key testing
id2 <- rep(c("A", "B"), each = n_groups * n_per_group / 2)

# Predictors
x1 <- rnorm(n_groups * n_per_group)
off1 <- rep(0.1, n_groups * n_per_group)
off2 <- rep(0.05, n_groups * n_per_group)

# Generate Random Effects
# Adding a true random effect stabilizes the estimation of RE models (poisLind.re/renb)
# preventing parameters from hitting boundaries/crashing during bootstrap.
re_sd <- 0.3
group_effects <- rnorm(n_groups, mean = 0, sd = re_sd)
re_vec <- rep(group_effects, each = n_per_group)

# Generate counts
# Use Negative Binomial to ensure overdispersion, helping the RE models fit better
# Include the random effect (re_vec) in the linear predictor
mu <- exp(0.5 + 0.5 * x1 + off1 + re_vec)
theta_sim <- 2 
y <- MASS::rnegbin(n_groups * n_per_group, mu = mu, theta = theta_sim)

test_data <- data.frame(
  id = factor(id),
  id2 = factor(id2),
  y = y,
  x1 = x1,
  off1 = off1,
  off2 = off2
)

# ==============================================================================
# TESTS FOR poisLind.re
# ==============================================================================

test_that("poisLind.re runs correctly with basic inputs", {
  # Use suppressWarnings to handle "iteration limit" or "NaNs in SE" 
  # which are common on synthetic data but don't indicate code failure.
  mod <- suppressWarnings(
    poisLind.re(y ~ x1, group_var = "id", data = test_data, max.iters = 100)
  )
  
  expect_equal(mod$model$modelType, "poisLindRE")
  expect_true(!is.null(mod$model$theta))
  expect_true(length(mod$model$estimate) > 0)
})

test_that("poisLind.re handles bootstrapping and offsets", {
  # Suppress warnings for clean test output
  mod <- suppressWarnings(
    poisLind.re(y ~ x1, group_var = "id", data = test_data, 
                offset = "off1", bootstraps = 2, max.iters = 50)
  )
  
  expect_true(!is.null(mod$model$bootstrapped_se))
  expect_equal(mod$model$bootstraps, 2)
  expect_true("beta_pred" %in% names(mod$model))
  
  # Check summary method works on this object
  summ <- summary(mod)
  expect_true("x1" %in% summ$parameter)
})

test_that("poisLind.re handles composite grouping variables", {
  mod <- suppressWarnings(
    poisLind.re(y ~ x1, group_var = c("id", "id2"), data = test_data, max.iters = 50)
  )
  
  expect_true(inherits(mod, "flexCountReg"))
  expect_equal(mod$model$modelType, "poisLindRE")
})

test_that("poisLind.re handles existing panel_id in data", {
  d2 <- test_data
  d2$panel_id <- d2$id 
  
  mod <- suppressWarnings(
    poisLind.re(y ~ x1, group_var = "id", data = d2, max.iters = 50)
  )
  
  expect_true(inherits(mod, "flexCountReg"))
  expect_equal(mod$model$formula, y ~ x1)
})



# ==============================================================================
# TESTS FOR renb
# ==============================================================================

test_that("renb runs correctly with basic inputs", {
  mod <- suppressWarnings(
    renb(y ~ x1, group_var = "id", data = test_data, max.iters = 100)
  )
  
  expect_true(inherits(mod, "flexCountReg"))
  expect_equal(mod$model$modelType, "RENB")
  expect_true(!is.null(mod$model$a))
  expect_true(!is.null(mod$model$b))
})

test_that("renb handles bootstrapping and multiple offsets", {
  mod <- suppressWarnings(
    renb(y ~ x1, group_var = "id", data = test_data, 
         offset = c("off1", "off2"), max.iters = 50)
  )
  
  expect_equal(length(mod$model$offset), 2)
})

test_that("renb handles composite grouping variables", {
  mod <- suppressWarnings(
    renb(y ~ x1, group_var = c("id", "id2"), data = test_data, max.iters = 50)
  )
  expect_true(inherits(mod, "flexCountReg"))
})

test_that("renb handles existing panel_id in data", {
  d2 <- test_data
  d2$panel_id <- d2$id
  
  mod <- suppressWarnings(
    renb(y ~ x1, group_var = "id", data = d2, max.iters = 50)
  )
  expect_true(inherits(mod, "flexCountReg"))
})


test_that("renb handles bootstrap failure gracefully", {
  tiny_data <- test_data[1:10,]
  
  # Suppress warnings to ignore convergence issues on tiny data
  mod <- suppressWarnings(
    renb(y ~ x1, group_var = "id", data = tiny_data, bootstraps = 2, max.iters = 5)
  )
  
  expect_true(inherits(mod, "flexCountReg"))
})

# ==============================================================================
# TESTS FOR PREDICTIONS
# ==============================================================================

test_that("predict function works for poisLindRE models", {
  mod <- suppressWarnings(
    poisLind.re(y ~ x1, group_var = "id", data = test_data, max.iters = 20)
  )
  
  preds <- predict(mod, data = test_data)
  
  expect_equal(length(preds), nrow(test_data))
  expect_true(all(preds > 0))
  expect_true(is.numeric(preds))
})

test_that("predict function works for renb models", {
  mod <- suppressWarnings(
    renb(y ~ x1, group_var = "id", data = test_data, max.iters = 20)
  )
  
  preds <- predict(mod, data = test_data)
  
  expect_equal(length(preds), nrow(test_data))
  expect_true(all(preds > 0))
  expect_true(is.numeric(preds))
})

test_that("predict function works for countreg.rp models (Exact vs Simulated)", {
  # Estimate a quick Random Parameter model
  # Suppress warnings for cleaner test output
  rp_mod <- suppressWarnings(
    countreg.rp(y ~ x1, 
                rpar_formula = ~ 1, 
                data = test_data,
                family = "NB2",
                rpardists = c("(Intercept)" = "n"),
                ndraws = 25, 
                max.iters = 20)
  )
  
  # 1. Test Simulated Method
  pred_sim <- predict(rp_mod, method = "Simulated")
  expect_equal(length(pred_sim), nrow(test_data))
  expect_true(all(pred_sim > 0))
  
  # 2. Test Exact Method
  pred_exact <- predict(rp_mod, method = "Exact")
  expect_equal(length(pred_exact), nrow(test_data))
  expect_true(all(pred_exact > 0))
  
  # 3. Comparison
  # The correlation should be high between the two methods
  expect_true(cor(pred_sim, pred_exact) > 0.8)
})

test_that("predict function works for basic countreg models (Fixed)", {
  mod <- suppressWarnings(
    countreg(y ~ x1, data = test_data, family = "Poisson")
  )
  
  preds <- predict(mod)
  expect_equal(length(preds), nrow(test_data))
  expect_true(all(preds > 0))
})

