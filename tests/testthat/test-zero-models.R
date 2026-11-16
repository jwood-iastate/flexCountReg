# Zero-inflated and hurdle models (zero_formula / zero_type in countreg()).

zi_data <- function(n = 3000, seed = 5) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rbinom(n, 1, 0.5)
  mu <- exp(0.5 + 0.4 * x)
  y <- ifelse(runif(n) < plogis(-1 + 1.2 * z), 0,
              rnbinom(n, size = 2, mu = mu))
  data.frame(y, x, z)
}

test_that("zero-inflated and hurdle models match pscl", {
  skip_if_not_installed("pscl")
  d <- zi_data()
  zp <- countreg(y ~ x, d, family = "Poisson", zero_formula = ~ z,
                 method = "BFGS")
  ref <- pscl::zeroinfl(y ~ x | z, d, dist = "poisson")
  expect_equal(as.numeric(logLik(zp)), as.numeric(logLik(ref)),
               tolerance = 1e-6)
  expect_equal(unname(coef(zp)), unname(coef(ref)), tolerance = 1e-3)
  expect_equal(as.numeric(fitted(zp)), as.numeric(fitted(ref)),
               tolerance = 1e-4)
  
  zn <- countreg(y ~ x, d, family = "NB2", zero_formula = ~ z, method = "NR")
  refn <- pscl::zeroinfl(y ~ x | z, d, dist = "negbin")
  expect_equal(as.numeric(logLik(zn)), as.numeric(logLik(refn)),
               tolerance = 1e-6)
  
  hn <- countreg(y ~ x, d, family = "NB2", zero_formula = ~ z,
                 zero_type = "hurdle", method = "BFGS")
  refh <- pscl::hurdle(y ~ x | z, d, dist = "negbin")
  expect_equal(as.numeric(logLik(hn)), as.numeric(logLik(refh)),
               tolerance = 1e-6)
  expect_equal(as.numeric(fitted(hn)), as.numeric(fitted(refh)),
               tolerance = 1e-2)
})

test_that("zero-inflated NB2 recovers simulated parameters", {
  d <- zi_data(n = 4000)
  m <- countreg(y ~ x, d, family = "NB2", zero_formula = ~ z, method = "NR")
  cf <- coef(m)
  se <- sqrt(diag(vcov(m)))
  truth <- c(0.5, 0.4, log(0.5), -1, 1.2)
  expect_true(all(abs(cf - truth) < 4 * se))
  expect_equal(tidy(m)$component,
               c("mean", "mean", "distribution", "zero", "zero"))
})

test_that("predictions, probabilities and residuals include the zero model", {
  d <- zi_data()
  m <- countreg(y ~ x, d, family = "NB2", zero_formula = ~ z, method = "BFGS")
  nd <- d[1:5, ]
  pi0 <- predict(m, newdata = nd, type = "zero")
  mu <- predict(m, newdata = nd, type = "true")
  expect_equal(predict(m, newdata = nd), (1 - pi0) * mu)
  P <- predict(m, newdata = nd, type = "prob", max_count = 300)
  expect_equal(rowSums(P), rep(1, 5), tolerance = 1e-8)
  size <- 1 / exp(coef(m)["ln(alpha)"])
  expect_equal(P[, "0"], pi0 + (1 - pi0) * dnbinom(0, size = size, mu = mu))
  expect_equal(as.numeric(P %*% (0:300)), predict(m, newdata = nd),
               tolerance = 1e-6)
  set.seed(1)
  r <- residuals(m, type = "quantile")
  expect_equal(sd(r), 1, tolerance = 0.06)
  # Pearson variance from the mixture moments
  v <- rowSums(P * outer(rep(1, 5), (0:300)^2)) - (P %*% (0:300))^2
  pr <- residuals(m, type = "pearson")[1:5]
  expect_equal(as.numeric(pr),
               as.numeric((d$y[1:5] - predict(m, newdata = nd)) / sqrt(v)),
               tolerance = 1e-6)
  expect_null(summary(m)$lr_test)
})

test_that("hurdle predictions and simulation", {
  d <- zi_data()
  m <- countreg(y ~ x, d, family = "Poisson", zero_formula = ~ z,
                zero_type = "hurdle", method = "BFGS")
  nd <- d[1:5, ]
  P <- predict(m, newdata = nd, type = "prob", max_count = 100)
  expect_equal(P[, "0"], predict(m, newdata = nd, type = "zero"))
  expect_equal(rowSums(P), rep(1, 5), tolerance = 1e-10)
  expect_equal(as.numeric(P %*% (0:100)), predict(m, newdata = nd),
               tolerance = 1e-8)
  sims <- as.matrix(simulate(m, nsim = 100, seed = 2))
  # The zero share is matched exactly by the zero model; the simulated mean
  # matches the model's fitted mean
  expect_equal(mean(sims == 0), mean(d$y == 0), tolerance = 0.03)
  expect_equal(mean(sims), mean(fitted(m)), tolerance = 0.02)
})

test_that("zero models validate their arguments", {
  d <- zi_data(500)
  expect_error(countreg(y ~ x, d, zero_formula = ~ z,
                        underreport_formula = ~ z), "cannot be combined")
  expect_error(countreg(y ~ x, d, zero_formula = ~ z, zero_type = "hurdel"))
  m <- countreg(y ~ x, d, method = "BFGS")
  expect_error(predict(m, newdata = d, type = "zero"), "zero_formula")
})
