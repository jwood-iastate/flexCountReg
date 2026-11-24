# Multivariate count models: countreg_mv()

sim_mv <- function(n = 800, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), z = stats::rbinom(n, 1, 0.5),
                  len = stats::runif(n, 0.5, 2), site = rep(seq_len(n / 4),
                                                           each = 4))
  d$a <- stats::rnbinom(n, size = 2, mu = exp(0.8 + 0.5 * d$x) * d$len)
  d$b <- stats::rpois(n, exp(-0.2 + 0.3 * d$x - 0.4 * d$z) * d$len)
  d$c <- stats::rnbinom(n, size = 0.7, mu = exp(-1 + 0.6 * d$z) * d$len)
  d
}

test_that("independent equations equal separate countreg() fits", {
  d <- sim_mv()
  m <- countreg_mv(cbind(a, b, c) ~ x + z + offset(log(len)), data = d,
                   family = c(a = "NB2", b = "POISSON", c = "NB1"),
                   dis_param_formula_1 = list(a = ~ z))
  expect_s3_class(m, "countreg_mv")
  expect_true(m$converged)
  r <- list(
    countreg(a ~ x + z + offset(log(len)), data = d, family = "NB2",
             dis_param_formula_1 = ~ z, method = "BFGS"),
    countreg(b ~ x + z + offset(log(len)), data = d, family = "POISSON",
             method = "BFGS"),
    countreg(c ~ x + z + offset(log(len)), data = d, family = "NB1",
             method = "BFGS"))
  expect_equal(as.numeric(logLik(m)),
               sum(vapply(r, function(f) as.numeric(logLik(f)), 0)),
               tolerance = 1e-7)
  expect_equal(unname(coef(m)), unname(unlist(lapply(r, coef))),
               tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(vcov(m)))),
               unname(unlist(lapply(r, function(f) sqrt(diag(vcov(f)))))),
               tolerance = 1e-3)
  expect_equal(names(coef(m))[1:4], c("a:(Intercept)", "a:x", "a:z",
                                      "a:ln(alpha):(Intercept)"))
  expect_equal(attr(logLik(m), "df"), length(coef(m)))
})

test_that("the analytic gradient matches numeric differences", {
  skip_if_not_installed("numDeriv")
  d <- sim_mv(n = 200)
  E <- countreg_mv(cbind(a, b, c) ~ x + z, data = d,
                   family = c("NB1", "PLN", "NB2"),
                   dis_param_formula_1 = ~ x, .setup_only = TRUE)
  set.seed(2)
  p <- flexCountReg:::.fcr_mv_start(E) + stats::runif(length(E$names),
                                                      -0.1, 0.1)
  f <- function(q) sum(flexCountReg:::.fcr_mv_ll(q, E)$ll)
  expect_equal(colSums(flexCountReg:::.fcr_mv_ll(p, E, grad = TRUE)$G),
               numDeriv::grad(f, p), tolerance = 1e-5)
})

test_that("formulas per outcome, weights and missing rows", {
  d <- sim_mv()
  d$b[3] <- NA
  expect_message(
    m <- countreg_mv(list(first = a ~ x, second = b ~ z), data = d,
                     family = "POISSON"), "1 row")
  expect_equal(m$outcomes, c("first", "second"))
  expect_equal(names(coef(m)), c("first:(Intercept)", "first:x",
                                 "second:(Intercept)", "second:z"))
  expect_equal(m$nobs, nrow(d) - 1)
  # Frequency weights equal the expanded data
  d2 <- sim_mv(n = 200)
  d2$w <- rep(1:2, 100)
  mw <- countreg_mv(cbind(a, b) ~ x, data = d2, family = "NB2",
                    weights = "w")
  me <- countreg_mv(cbind(a, b) ~ x, data = d2[rep(seq_len(200), d2$w), ],
                    family = "NB2")
  expect_equal(coef(mw), coef(me), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(mw)), as.numeric(logLik(me)),
               tolerance = 1e-8)
  mc <- countreg_mv(cbind(a, b) ~ x, data = d2, family = "NB2",
                    group_var = "site", stderr = "cluster")
  expect_false(isTRUE(all.equal(vcov(mc), vcov(mw))))
})

test_that("predictions, simulation and methods", {
  d <- sim_mv()
  m <- countreg_mv(cbind(a, b) ~ x + z + offset(log(len)), data = d,
                   family = c("NB2", "POISSON"))
  mu <- predict(m, newdata = d[1:5, ])
  expect_equal(dim(mu), c(5, 2))
  expect_equal(unname(mu[, 2]),
               exp(coef(m)["b:(Intercept)"] + coef(m)["b:x"] * d$x[1:5] +
                     coef(m)["b:z"] * d$z[1:5]) * d$len[1:5],
               ignore_attr = TRUE)
  expect_equal(predict(m, newdata = d[1:5, ], type = "variance")[, 2],
               mu[, 2])
  pr <- predict(m, newdata = d[1:3, ], type = "prob", max_count = 200)
  expect_named(pr, c("a", "b"))
  expect_equal(unname(rowSums(pr$a)), rep(1, 3), tolerance = 1e-8)
  expect_equal(as.vector(pr$b %*% (0:200)), unname(mu[1:3, 2]),
               tolerance = 1e-8)
  expect_equal(dim(fitted(m)), c(nrow(d), 2))
  expect_equal(dim(residuals(m, type = "pearson")), c(nrow(d), 2))
  s <- simulate(m, nsim = 2, seed = 1)
  expect_equal(dim(s$sim_1), c(nrow(d), 2))
  expect_equal(colMeans(s$sim_1), colMeans(m$y), tolerance = 0.15)
  td <- tidy(m, conf.int = TRUE)
  expect_equal(unique(td$outcome), c("a", "b"))
  expect_equal(sum(td$component == "distribution"), 1)
  expect_equal(glance(m)$outcomes, 2)
  expect_output(print(summary(m)), "a \\(NB2\\)")
  expect_output(print(m), "2 outcomes")
})

test_that("input errors", {
  d <- sim_mv(n = 40)
  expect_error(countreg_mv(a ~ x, data = d), "two or more outcomes")
  expect_error(countreg_mv(cbind(a, b) ~ x, data = d,
                           family = c(a = "NB2", q = "PIG")), "Unknown: q")
  expect_error(countreg_mv(cbind(a, b) ~ x, data = d,
                           family = c("NB2", "PIG", "NB1")), "one per outcome")
  expect_error(countreg_mv(cbind(a, b) ~ x, data = d,
                           dis_param_formula_1 = list(q = ~ x)),
               "Unknown outcome")
  expect_error(countreg_mv(cbind(a, b) ~ s(x), data = d), "Smooth")
  expect_error(countreg_mv(cbind(a, b) ~ x, data = d, stderr = "cluster"),
               "group_var")
})
