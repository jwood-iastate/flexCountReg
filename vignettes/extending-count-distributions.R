## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")


## ----setup--------------------------------------------------------------------
library(flexCountReg)


## ----nb2-entry----------------------------------------------------------------
nb2 <- flexCountReg:::.fcr_family_registry$NB2
names(nb2)
nb2$params
names(formals(nb2$dens))
body(nb2$dens)


## ----nb2-check----------------------------------------------------------------
y <- 0:6
all.equal(nb2$dens(y, predicted = 2.5, alpha = 0.4, sigma = NULL),
          dnbinom(y, size = 1 / 0.4, mu = 2.5))


## ----capabilities-------------------------------------------------------------
reg <- flexCountReg:::.fcr_family_registry
nest <- flexCountReg:::.fcr_family_nesting
label <- function(x) if (is.null(x)) "--" else x
capabilities <- data.frame(
  family = names(reg),
  param_1 = vapply(reg, function(f) label(f$params[[1]]), ""),
  param_2 = vapply(reg, function(f) label(f$params[[2]]), ""),
  variance = vapply(reg, function(f) !is.null(f$variance), TRUE),
  cdf = vapply(reg, function(f) !is.null(f$cdf), TRUE),
  own_mean = vapply(reg, function(f) !is.null(f$mean), TRUE),
  thinning = names(reg) %in% flexCountReg:::.fcr_thinning_families,
  lr_vs_poisson = names(reg) %in% nest$full[nest$restricted == "POISSON"],
  row.names = NULL)
capabilities


## ----dpo-source---------------------------------------------------------------
body(ddpo)
body(qdpo)


## ----contract-----------------------------------------------------------------
mu <- 3; rho <- 0.4; k <- 2
x <- 0:200
p <- ddelaporte(x, mu = mu, rho = rho, k = k)
c(total = sum(p),
  mean = sum(x * p),
  cdf_matches_pmf = isTRUE(all.equal(cumsum(p)[1:20],
                                     pdelaporte(0:19, mu, rho, k))),
  log_tail_direct = pdelaporte(150, mu, rho, k, lower.tail = FALSE,
                               log.p = TRUE))
# a quantile is the smallest count whose CDF reaches the probability
q <- qdelaporte(0.9, mu, rho, k)
c(q = q, cdf_q = pdelaporte(q, mu, rho, k),
  cdf_q_minus_1 = pdelaporte(q - 1, mu, rho, k))
# invalid parameters and non-integer counts
suppressWarnings(c(ddelaporte(2, mu = -1, rho = rho, k = k),
                   ddelaporte(2.5, mu = mu, rho = rho, k = k)))


## ----trace--------------------------------------------------------------------
set.seed(4)
d <- data.frame(x = rnorm(300))
d$y <- rnbinom(300, size = 2, mu = exp(0.5 + 0.4 * d$x))
fit <- countreg(y ~ x, data = d, family = "NB2", method = "BFGS")
mu_hat <- predict(fit)
alpha_hat <- exp(coef(fit)[["ln(alpha)"]])
# predict(type = "prob") uses the registry density
all.equal(unname(predict(fit, type = "prob", max_count = 3)[1, ]),
          nb2$dens(0:3, predicted = mu_hat[1], alpha = alpha_hat,
                   sigma = NULL))
# simulate() uses the registry's random draws
sims <- simulate(fit, nsim = 4000, seed = 1)
c(mean_ratio = mean(rowMeans(sims) / mu_hat),
  var_ratio = mean(apply(sims, 1, var) /
                     nb2$variance(mu_hat, alpha_hat, NULL)))
# quantile residuals use the cdf hook
summary(residuals(fit, type = "quantile"))


## ----test-example, eval = requireNamespace("testthat", quietly = TRUE)--------
testthat::test_that("countreg() recovers NB2 parameters", {
  set.seed(11)
  n <- 4000
  d <- data.frame(x = rnorm(n), t = runif(n, 0.5, 2))
  d$y <- rnbinom(n, size = 1 / 0.5, mu = d$t * exp(0.5 + 0.4 * d$x))
  m <- countreg(y ~ x + offset(log(t)), data = d, family = "NB2",
                method = "BFGS")
  testthat::expect_true(m$model$converged)
  est <- coef(m)
  se <- summary(m)$coefficients[, 2]
  truth <- c(0.5, 0.4, log(0.5))
  # within four standard errors: fails rarely by chance, often if wrong
  testthat::expect_true(all(abs(est - truth) < 4 * se))
})

