# Split by run time so the parts run in parallel (see Config/testthat/start-first).
# Regression tests for the corrections from the development review of
# flexCountReg 0.1.2 (findings F01-F28).

test_that("F22: inner mixing draws are separate, scrambled draws reproducible", {
  d <- washington_roads
  t_small <- system.time(
    a <- countreg_rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                     data = d, family = "PLN", ndraws = 30, mix_draws = 10,
                     method = "BFGS", max.iters = 10))["elapsed"]
  expect_equal(a$model$mix_draws, 10)
  expect_equal(countreg_rp(Total_crashes ~ lnaadt,
                           rpar_formula = ~ -1 + speed50, data = d,
                           family = "PLN", ndraws = 300, max.iters = 0,
                           method = "BFGS")$model$mix_draws, 200)
  # Scrambled draws: predictions reuse the fit's draws and leave the
  # caller's random-number stream unchanged
  set.seed(4)
  s1 <- countreg_rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                    data = d, ndraws = 50, scrambled = TRUE, method = "BFGS",
                    max.iters = 30)
  set.seed(4)
  s2 <- countreg_rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                    data = d, ndraws = 50, scrambled = TRUE, method = "BFGS",
                    max.iters = 30)
  expect_equal(coef(s1), coef(s2))
  set.seed(10)
  u1 <- runif(1)
  set.seed(10)
  p1 <- predict(s1, newdata = d[1:20, ], method = "Simulated")
  u2 <- runif(1)
  expect_equal(u1, u2)
  expect_equal(predict(s1, newdata = d[1:20, ], method = "Simulated"), p1)
})
