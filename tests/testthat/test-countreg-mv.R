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

# ---- Random parameters across equations --------------------------------------

sim_mvrp <- function(G = 300, Tn = 4, seed = 7, sdv = c(0.5, 0.4, 0.6),
                     cr = c(0.6, -0.4, 0.3), nb = TRUE) {
  set.seed(seed)
  id <- rep(seq_len(G), each = Tn)
  n <- G * Tn
  Rm <- diag(3)
  Rm[lower.tri(Rm)] <- cr
  Rm[upper.tri(Rm)] <- t(Rm)[upper.tri(Rm)]
  Sig <- diag(sdv) %*% Rm %*% diag(sdv)
  b <- (matrix(stats::rnorm(3 * G), G) %*% chol(Sig))[id, ]
  d <- data.frame(id = id, x = stats::rnorm(n), z = stats::rbinom(n, 1, 0.5))
  d$a <- stats::rpois(n, exp(0.5 + 0.4 * d$x + b[, 1]))
  mu_b <- exp(-0.3 + 0.3 * d$z + b[, 2])
  d$b <- if (nb) stats::rnbinom(n, size = 2, mu = mu_b) else
    stats::rpois(n, mu_b)
  d$c <- stats::rpois(n, exp(-0.8 + 0.5 * d$x + b[, 3]))
  d
}

test_that("the simulated likelihood has the right analytic gradient", {
  skip_if_not_installed("numDeriv")
  d <- sim_mvrp(G = 60)
  cases <- list(
    list(family = "POISSON", rpar_formula = ~ 1, group_var = NULL),
    list(family = c("POISSON", "NB2", "PIG"), rpar_formula = ~ 1),
    list(family = "NB2", rpar_formula = list(a = ~ 1 + x, b = ~ z),
         correlation = "within"),
    list(family = "NB1", rpar_formula = ~ 1, correlation = "none",
         dis_param_formula_1 = ~ z),
    list(family = "POISSON", rpar_formula = ~ 1,
         weights = rep(1:2, nrow(d) / 2)))
  for (cs in cases) {
    args <- utils::modifyList(list(
      formula = cbind(a, b, c) ~ x + z, group_var = "id", rpar_draws = 15,
      seed = 3, .setup_only = TRUE), cs)
    args$data <- d
    E <- do.call(countreg_mv, args)
    set.seed(2)
    p <- flexCountReg:::.fcr_mv_start(E) +
      stats::runif(length(E$names), -0.1, 0.1)
    f <- function(q) sum(flexCountReg:::.fcr_mv_sim_ll(q, E)$ll)
    expect_equal(colSums(flexCountReg:::.fcr_mv_sim_ll(p, E, grad = TRUE)$G),
                 numDeriv::grad(f, p), tolerance = 1e-5,
                 info = paste(names(cs), collapse = ","))
  }
})

test_that("the simulated likelihood approaches the exact integral", {
  # Bivariate Poisson-lognormal with correlated random intercepts per row:
  # exact by two-dimensional Gauss-Hermite quadrature
  set.seed(4)
  n <- 200
  d <- data.frame(x = stats::rnorm(n))
  b <- matrix(stats::rnorm(2 * n), n) %*% chol(matrix(c(0.25, 0.15, 0.15,
                                                        0.36), 2))
  d$a <- stats::rpois(n, exp(0.3 + 0.4 * d$x + b[, 1]))
  d$b <- stats::rpois(n, exp(-0.2 + b[, 2]))
  E <- countreg_mv(cbind(a, b) ~ x, data = d, family = "POISSON",
                   rpar_formula = list(a = ~ 1, b = ~ 1), rpar_draws = 2000,
                   seed = 1, .setup_only = TRUE)
  p <- stats::setNames(c(0.3, 0.4, -0.2, 0, 0.5, 0.3, 0.52), E$names)
  sim <- sum(flexCountReg:::.fcr_mv_sim_ll(p, E)$ll)
  L <- flexCountReg:::.fcr_alloc_chol(p, E$RS)
  gh <- flexCountReg:::.fcr_gh_grid(30, 2)
  B <- gh$Z %*% t(L)
  exact <- sum(vapply(seq_len(n), function(i) {
    log(sum(gh$w * stats::dpois(d$a[i], exp(0.3 + 0.4 * d$x[i] + B[, 1])) *
              stats::dpois(d$b[i], exp(-0.2 + B[, 2]))))
  }, numeric(1)))
  expect_lt(abs(sim - exact), 0.2)
})

test_that("correlated random intercepts are recovered", {
  skip_on_cran()
  d <- sim_mvrp(G = 400, Tn = 6, nb = FALSE)
  m <- countreg_mv(cbind(a, b, c) ~ x + z, data = d, family = "POISSON",
                   rpar_formula = ~ 1, group_var = "id", rpar_draws = 100,
                   seed = 2)
  expect_true(m$converged)
  cf <- coef(m)
  expect_true(all(cf[grep("^sd:", names(cf))] > 0))
  vc <- flexCountReg:::.fcr_alloc_varcorr(m)
  expect_equal(unname(sqrt(diag(vc))), c(0.5, 0.4, 0.6), tolerance = 0.2)
  cr <- stats::cov2cor(vc)
  expect_lt(max(abs(c(cr[2, 1], cr[3, 1], cr[3, 2]) - c(0.6, -0.4, 0.3))),
            0.2)
  m0 <- countreg_mv(cbind(a, b, c) ~ x + z, data = d, family = "POISSON")
  expect_lt(AIC(m), AIC(m0) - 100)
  expect_output(print(summary(m)), "Standard deviations and correlations")
  expect_true("random" %in% tidy(m)$component)
})

test_that("joint and conditional probabilities", {
  d <- sim_mvrp(G = 80)
  m <- countreg_mv(cbind(a, b) ~ x, data = d, family = "POISSON",
                   rpar_formula = ~ 1, group_var = "id", rpar_draws = 50,
                   seed = 1)
  nd <- d[1:2, ]
  grid <- expand.grid(a = 0:60, b = 0:40)
  tot <- vapply(1:2, function(i) sum(predict(
    m, newdata = nd[rep(i, nrow(grid)), ], type = "joint", y = grid)), 0)
  expect_equal(tot, c(1, 1), tolerance = 1e-6)
  cp <- predict(m, newdata = nd, type = "conditional", outcome = "b",
                given = c(a = 3), max_count = 60)
  expect_equal(unname(rowSums(cp)), c(1, 1), tolerance = 1e-6)
  # P(b = 1 | a = 3) = P(a = 3, b = 1) / P(a = 3)
  j31 <- predict(m, newdata = nd, type = "joint", y = c(a = 3, b = 1))
  pa3 <- predict(m, newdata = nd, type = "prob", max_count = 3)$a[, "3"]
  expect_equal(unname(cp[, "1"]), j31 / pa3, tolerance = 1e-8)
  # Without random parameters the conditional equals the marginal
  mi <- countreg_mv(cbind(a, b) ~ x, data = d, family = "POISSON")
  expect_equal(unname(predict(mi, newdata = nd, type = "conditional",
                              outcome = "b", given = c(a = 3),
                              max_count = 5)),
               unname(predict(mi, newdata = nd, type = "prob",
                              max_count = 5)$b))
  # Posterior (Empirical Bayes) predictions for estimated units
  pm <- predict(m, newdata = d[1:4, ])
  pp <- predict(m, newdata = d[1:4, ], random = "posterior")
  expect_false(isTRUE(all.equal(pm, pp)))
  expect_error(predict(m, newdata = d[1:4, "x", drop = FALSE],
                       random = "posterior"), "grouping column")
  expect_error(predict(m, newdata = nd, type = "joint", y = c(a = 1)),
               "every outcome")
  expect_error(predict(m, newdata = nd, type = "conditional", outcome = "b",
                       given = c(b = 1)), "outcome itself")
  s <- simulate(m, nsim = 1, seed = 1)$sim_1
  expect_equal(dim(s), c(nrow(d), 2))
})

test_that("random-parameter input errors", {
  d <- sim_mvrp(G = 20)
  expect_error(countreg_mv(cbind(a, b) ~ x, data = d,
                           rpar_formula = ~ offset(x)), "offset")
  expect_error(countreg_mv(cbind(a, b) ~ x, data = d,
                           rpar_formula = list(q = ~ 1)), "Unknown outcome")
})

test_that("countreg_alloc() makes the intercept random only when asked", {
  set.seed(3)
  n <- 120
  d <- data.frame(x = stats::rnorm(n), site = rep(1:40, each = 3))
  Y <- t(vapply(stats::rpois(n, 3), function(k)
    as.numeric(stats::rmultinom(1, k, c(0.6, 0.3, 0.1))), numeric(3)))
  d$a <- Y[, 1]
  d$b <- Y[, 2]
  d$c <- Y[, 3]
  E <- countreg_alloc(cbind(a, b, c) ~ x, data = d, family = "POISSON",
                      alloc_formula = ~ x, group_var = "site",
                      rpar_total = ~ x, rpar_draws = 5, .setup_only = TRUE)
  expect_equal(E$RS$labels, "Total:x")
  E1 <- countreg_alloc(cbind(a, b, c) ~ x, data = d, family = "POISSON",
                       alloc_formula = ~ x, group_var = "site",
                       rpar_total = ~ 1 + x, rpar_draws = 5,
                       .setup_only = TRUE)
  expect_equal(E1$RS$labels, c("Total:(Intercept)", "Total:x"))
})

# ---- Copulas -----------------------------------------------------------------

test_that("the bivariate normal distribution function is exact", {
  ref <- function(a, b, r) stats::integrate(function(x) stats::dnorm(x) *
      stats::pnorm((b - r * x) / sqrt(1 - r^2)), -Inf, a,
    rel.tol = 1e-13, abs.tol = 1e-15)$value
  for (r in c(-0.99, -0.93, -0.5, 0, 0.2, 0.8, 0.95)) {
    for (a in c(-2, 0.3, 1.7)) for (b in c(-1, 0, 2.2)) {
      expect_equal(flexCountReg:::.fcr_pbvnorm(a, b, r), ref(a, b, r),
                   tolerance = 1e-10)
    }
  }
  expect_equal(flexCountReg:::.fcr_pbvnorm(Inf, 0.3, 0.5), stats::pnorm(0.3))
  expect_equal(flexCountReg:::.fcr_pbvnorm(-Inf, 1, 0.5), 0)
})

test_that("rectangle probabilities are a distribution with the right margins", {
  rect <- flexCountReg:::.fcr_cop_rect
  F1 <- function(y) stats::ppois(y, 1.3)
  F2 <- function(y) stats::pnbinom(y, size = 2, mu = 0.8)
  F3 <- function(y) stats::ppois(y, 0.5)
  g <- expand.grid(a = 0:25, b = 0:25)
  hi <- cbind(F1(g$a), F2(g$b))
  lo <- cbind(F1(g$a - 1), F2(g$b - 1))
  indep <- stats::dpois(g$a, 1.3) * stats::dnbinom(g$b, size = 2, mu = 0.8)
  for (cf in list(c("gaussian", 0.6), c("frank", -4), c("clayton", 2),
                  c("gumbel", 2))) {
    P <- rect(hi, lo, cf[1], as.numeric(cf[2]))
    expect_equal(sum(P), 1, tolerance = 1e-10, info = cf[1])
    expect_gt(min(P), -1e-12)
    # The margins are the outcomes' own distributions
    expect_equal(as.vector(tapply(P, g$a, sum)), stats::dpois(0:25, 1.3),
                 tolerance = 1e-10, info = cf[1])
  }
  expect_equal(rect(hi, lo, "gaussian", 0), indep, tolerance = 1e-12)
  expect_equal(rect(hi, lo, "clayton", 1e-9), indep, tolerance = 1e-8)
  expect_equal(rect(hi, lo, "frank", 0), indep, tolerance = 1e-12)
  g3 <- expand.grid(a = 0:15, b = 0:15, c = 0:12)
  for (fam in c("clayton", "gumbel", "frank")) {
    P3 <- rect(cbind(F1(g3$a), F2(g3$b), F3(g3$c)),
               cbind(F1(g3$a - 1), F2(g3$b - 1), F3(g3$c - 1)), fam, 1.5)
    # (the grid leaves out about 3e-8 of the probability)
    expect_equal(sum(P3), 1, tolerance = 1e-6, info = fam)
  }
  # An outcome left out: the other's own probabilities
  expect_equal(rect(cbind(F1(0:20), NA), cbind(F1(-1:19), NA), "gumbel", 2),
               stats::dpois(0:20, 1.3), tolerance = 1e-12)
})

test_that("copula draws have the copula's Kendall's tau", {
  set.seed(1)
  for (cf in list(list("gaussian", 2, 0.6), list("frank", 2, -4),
                  list("frank", 3, 5), list("clayton", 3, 2),
                  list("gumbel", 3, 2))) {
    U <- flexCountReg:::.fcr_rcopula(3000, cf[[2]], cf[[1]], cf[[3]])
    expect_lt(abs(stats::cor(U[, 1], U[, 2], method = "kendall") -
                    flexCountReg:::.fcr_cop_tau(cf[[3]], cf[[1]])), 0.04)
    expect_lt(abs(mean(U) - 0.5), 0.02)
  }
})

sim_cop <- function(n, family, theta, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), z = stats::rbinom(n, 1, 0.5))
  U <- flexCountReg:::.fcr_rcopula(n, 2, family, theta)
  d$a <- stats::qnbinom(U[, 1], size = 1.5, mu = exp(0.6 + 0.5 * d$x))
  d$b <- stats::qpois(U[, 2], exp(-0.2 + 0.4 * d$z))
  d
}

test_that("copula models recover the dependence and keep the margins", {
  d <- sim_cop(1500, "clayton", 2, seed = 2)
  m <- countreg_mv(cbind(a, b) ~ x + z, data = d,
                   family = c("NB2", "POISSON"), copula = "clayton")
  expect_true(m$converged)
  expect_equal(unname(m$copula$theta), 2, tolerance = 0.15)
  expect_equal(unname(coef(m)[c("a:x", "b:z")]), c(0.5, 0.4),
               tolerance = 0.2)
  expect_output(print(summary(m)), "Kendall's tau")
  expect_true("copula" %in% tidy(m)$component)
  mi <- countreg_mv(cbind(a, b) ~ x + z, data = d,
                    family = c("NB2", "POISSON"))
  expect_gt(as.numeric(logLik(m) - logLik(mi)), 50)
  nd <- d[1:3, ]
  # The margins are the equations': expected counts are those of the
  # margins, and the joint probabilities of all pairs sum to 1
  expect_equal(predict(m, newdata = nd)[, "b"],
               exp(coef(m)["b:(Intercept)"] + coef(m)["b:x"] * nd$x +
                     coef(m)["b:z"] * nd$z),
               ignore_attr = TRUE)
  g <- expand.grid(a = 0:80, b = 0:20)
  expect_equal(sum(predict(m, newdata = nd[rep(1, nrow(g)), ], type = "joint",
                           y = g)), 1, tolerance = 1e-6)
  cp <- predict(m, newdata = nd, type = "conditional", outcome = "b",
                given = c(a = 4), max_count = 20)
  expect_equal(unname(rowSums(cp)), rep(1, 3), tolerance = 1e-8)
  j <- predict(m, newdata = nd, type = "joint", y = c(a = 4, b = 2))
  pa <- predict(m, newdata = nd, type = "prob", max_count = 4)$a[, "4"]
  expect_equal(unname(cp[, "2"]), j / pa, tolerance = 1e-8)
  s <- simulate(m, nsim = 1, seed = 2)$sim_1
  expect_equal(colMeans(s), colMeans(m$y), tolerance = 0.1)
  expect_gt(stats::cor(s[, 1], s[, 2], method = "kendall"), 0.2)
  # The Gaussian copula with a negative dependence
  dg <- sim_cop(1500, "gaussian", -0.5, seed = 4)
  mg <- countreg_mv(cbind(a, b) ~ x + z, data = dg,
                    family = c("NB2", "POISSON"), copula = "gaussian")
  expect_equal(unname(mg$copula$theta), -0.5, tolerance = 0.15)
})

test_that("copula input errors", {
  d <- sim_mv(n = 40)
  expect_error(countreg_mv(cbind(a, b, c) ~ x, data = d,
                           copula = "gaussian"), "two outcomes")
  expect_error(countreg_mv(cbind(a, b) ~ x, data = d, copula = "clayton",
                           rpar_formula = ~ 1), "cannot yet be combined")
})
