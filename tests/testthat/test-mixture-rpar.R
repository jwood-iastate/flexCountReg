# Random parameters within the classes of countreg_mix() / countreg_lc()

sim_lcrp <- function(G = 300, Tn = 4, seed = 1, sd = c(0.4, 0.4), pi1 = 0.6) {
  set.seed(seed)
  id <- rep(seq_len(G), each = Tn)
  cls <- ifelse(stats::runif(G) < pi1, 1, 2)[id]
  x <- stats::rnorm(G * Tn)
  b0 <- c(0.2, 1.5)[cls]
  b1 <- c(0.5, -0.3)[cls] + stats::rnorm(G)[id] * sd[cls]
  data.frame(id = id, x = x, y = stats::rpois(G * Tn, exp(b0 + b1 * x)),
             z = id %% 2)
}

# The training problem of a call, without fitting
lcrp_spec <- function(data, ..., seed = 7) {
  a <- utils::modifyList(list(
    entry = "countreg_lc", formula = y ~ x, family = "POISSON",
    components = 2L, components_missing = FALSE, mixing_formula = ~ 1,
    fixed_mixing = NULL, shared = NULL, dis_param_formula_1 = NULL,
    dis_param_formula_2 = NULL, offset = NULL, weights = NULL,
    group_var = "id", method = "BFGS", max.iters = 100, n_starts = 1,
    seed = 1, ndraws = 1500, engine = "quadrature",
    control = mixture_control(), rpar_formula = ~ x, rpardists = NULL,
    correlated = FALSE, rpar_draws = 20), list(...))
  st <- flexCountReg:::.fcr_fm_settings(a)
  st$rp$seed <- seed
  spec <- flexCountReg:::.fcr_fm_setup(st, data)
  spec$st <- st
  spec
}

test_that("the likelihood and the EM objective have the right gradients", {
  skip_if_not_installed("numDeriv")
  d <- sim_lcrp(G = 40)
  cases <- list(
    list(),
    list(shared = list(rpar = TRUE)),
    list(family = "NB2", rpar_formula = ~ 1 + x, correlated = TRUE),
    list(rpardists = "ln", mixing_formula = ~ z),
    list(family = c("POISSON", "NB2")),
    list(group_var = NULL),
    list(formula = y ~ 1, rpar_formula = ~ 1))
  for (cs in cases) {
    spec <- do.call(lcrp_spec, c(list(data = d), cs))
    set.seed(3)
    th <- flexCountReg:::.fcr_fm_pooled(spec)$theta
    th <- th + stats::rnorm(length(th), 0, 0.1)
    f <- function(p) sum(flexCountReg:::.fcr_fm_loglik_units(p, spec))
    expect_equal(colSums(flexCountReg:::.fcr_fm_grad_units(th, spec)),
                 numDeriv::grad(f, th), tolerance = 1e-5,
                 info = paste(names(cs), collapse = ","))
    tau <- flexCountReg:::.fcr_fm_eval(th, spec)$tau
    idx <- spec$map$comp_par
    q <- function(p) flexCountReg:::.fcr_fm_rp_Q(p, idx, th, spec, tau, 1:2)
    expect_equal(flexCountReg:::.fcr_fm_rp_Q(th[idx], idx, th, spec, tau, 1:2,
                                             grad = TRUE),
                 numDeriv::grad(q, th[idx]), tolerance = 1e-5,
                 info = paste(names(cs), collapse = ","))
  }
})

test_that("the simulated likelihood approaches the exact integral", {
  # Two classes with a random intercept each: exact by Gauss-Hermite
  set.seed(5)
  G <- 150
  id <- rep(seq_len(G), each = 4)
  cls <- ifelse(stats::runif(G) < 0.5, 1, 2)[id]
  d <- data.frame(id = id, x = stats::rnorm(4 * G))
  # Class-specific scales 0.6 and 0.3, as evaluated below (at scales far
  # from the data's, panels deep in a class's tail need far more draws)
  d$y <- stats::rpois(4 * G, exp(c(0, 1.2)[cls] + stats::rnorm(G)[id] *
                                   c(0.6, 0.3)[cls] + 0.4 * d$x))
  spec <- lcrp_spec(d, rpar_formula = ~ 1, rpar_draws = 1000)
  nm <- spec$map$names
  th <- stats::setNames(numeric(length(nm)), nm)
  th[c("Component1:Mean:(Intercept)", "Component1:Mean:x",
       "Component1:SD:(Intercept)", "Component2:Mean:(Intercept)",
       "Component2:Mean:x", "Component2:SD:(Intercept)",
       "Mixing:Component1:(Intercept)")] <- c(0, 0.4, 0.6, 1.2, 0.4, 0.3, 0)
  sim <- sum(flexCountReg:::.fcr_fm_loglik_units(th, spec))
  gh <- flexCountReg:::.fcr_gh_rule(60)
  A <- vapply(1:2, function(k) {
    b0 <- c(0, 1.2)[k]
    s <- c(0.6, 0.3)[k]
    L <- vapply(gh$z, function(z) rowsum(stats::dpois(
      d$y, exp(b0 + s * z + 0.4 * d$x), log = TRUE), d$id)[, 1], numeric(G))
    mx <- apply(L, 1, max)
    mx + log(colSums(t(exp(L - mx)) * gh$w))
  }, numeric(G))
  exact <- sum(log(exp(A) %*% c(0.5, 0.5)))
  # Simulation error of 150 panels' log-likelihoods with 1000 draws each
  expect_lt(abs(sim - exact), 0.5)
})

test_that("classes and random parameters are recovered", {
  skip_on_cran()
  d <- sim_lcrp(G = 400, Tn = 5, seed = 11, sd = c(0.5, 0.3))
  m <- countreg_lc(y ~ x, data = d, family = "POISSON", classes = 2,
                   group_var = "id", rpar_formula = ~ x, rpar_draws = 60,
                   n_starts = 2, seed = 1)
  cf <- coef(m)
  expect_equal(unname(cf[c("Component1:Mean:(Intercept)", "Component1:Mean:x",
                           "Component2:Mean:(Intercept)",
                           "Component2:Mean:x")]),
               c(0.2, 0.5, 1.5, -0.3), tolerance = 0.15)
  expect_equal(unname(cf[c("Component1:SD:x", "Component2:SD:x")]),
               c(0.5, 0.3), tolerance = 0.25)
  expect_true(all(cf[grep("SD:", names(cf))] > 0))
  m0 <- countreg_lc(y ~ x, data = d, family = "POISSON", classes = 2,
                    group_var = "id", n_starts = 2, seed = 1)
  expect_lt(AIC(m), AIC(m0) - 100)
})

test_that("methods integrate over the random parameters", {
  d <- sim_lcrp(G = 80, Tn = 4, seed = 3)
  m <- countreg_lc(y ~ x, data = d, family = "POISSON", classes = 2,
                   group_var = "id", rpar_formula = ~ x, rpar_draws = 30,
                   n_starts = 2, seed = 1)
  expect_output(print(summary(m)), "SD:x")
  expect_true("random" %in% tidy(m)$component)
  expect_true(all(coef(m)[grep("SD:", names(coef(m)))] > 0))
  # The component mean is the mean over the random slope: for a normal
  # slope, exactly exp(b0 + b1 x + s^2 x^2 / 2)
  cf <- coef(m)
  M <- predict(m, newdata = d[1:5, ], type = "component")
  ex <- exp(cf["Component1:Mean:(Intercept)"] +
              cf["Component1:Mean:x"] * d$x[1:5] +
              cf["Component1:SD:x"]^2 * d$x[1:5]^2 / 2)
  expect_equal(unname(M[, 1]), unname(ex), tolerance = 1e-10)
  pr <- predict(m, newdata = d[1:3, ], type = "prob", max_count = 150)
  expect_equal(rowSums(pr), rep(1, 3), tolerance = 1e-6)
  # Probabilities (simulated with 500 draws) and the exact expected counts
  # agree up to simulation error
  expect_equal(as.vector(pr %*% (0:150)),
               predict(m, newdata = d[1:3, ]), tolerance = 0.04)
  post <- component_probs(m, newdata = d[d$id %in% 1:3, ],
                          conditioning = "posterior")
  expect_equal(unname(post[1, ]), unname(post[4, ]))
  h <- d[d$id == 1, ][1, ]
  nd <- d[d$id == 1, ][-1, ]
  expect_length(predict(m, newdata = nd, conditioning = "history",
                        history = h), 3)
  q <- residuals(m, type = "quantile")
  expect_true(all(is.finite(q)))
  s <- simulate(m, nsim = 3, seed = 1)
  expect_equal(mean(colMeans(s)), mean(d$y), tolerance = 0.15)
  mm <- countreg_mix(y ~ x, data = d, family = "POISSON", components = 2,
                     group_var = "id", rpar_formula = ~ x, rpar_draws = 30,
                     n_starts = 2, seed = 1)
  expect_equal(as.numeric(logLik(mm)), as.numeric(logLik(m)))
})

test_that("random-parameter input errors", {
  d <- sim_lcrp(G = 20)
  expect_error(countreg_lc(y ~ 1, data = d, rpar_formula = ~ x,
                           group_var = "id"), "also be in `formula`")
  expect_error(countreg_lc(y ~ x, data = d, rpardists = "ln"),
               "need `rpar_formula`")
  expect_error(countreg_lc(y ~ x, data = d, rpar_formula = ~ x,
                           correlated = TRUE, rpardists = "ln"),
               "Correlated random parameters are normal")
  expect_error(countreg_lc(y ~ x, data = d, rpar_formula = ~ offset(x)),
               "offset")
  expect_error(countreg_lc(y ~ x, data = d, rpar_formula = ~ s(x)),
               "Smooth")
  expect_error(countreg_lc(y ~ x, data = d, shared = list(rpar = 1)),
               "TRUE or FALSE")
  expect_error(countreg_lc(y ~ x, data = d, shared = "mean"),
               "`mean`, `dis1`, `dis2` and `rpar`")
})

test_that("re-signed Cholesky columns keep the fit at the simulated maximum", {
  # A start with a negative diagonal in class 1's Cholesky factor: the fit
  # reports a positive diagonal, with that class's draws re-signed, so the
  # likelihood at the reported coefficients is the reported maximum
  d <- sim_lcrp(G = 60, sd = c(0.5, 0.5))
  args <- list(y ~ x, data = d, classes = 2, group_var = "id",
               rpar_formula = ~ 1 + x, correlated = TRUE, rpar_draws = 20,
               n_starts = 1, seed = 2)
  m0 <- do.call(countreg_lc, args)
  sv <- coef(m0)
  j <- grep("^Component1:Chol:.*_\\(Intercept\\)$", names(sv))
  expect_length(j, 2)
  sv[j] <- -sv[j]
  m <- do.call(countreg_lc, c(args, list(start.vals = sv)))
  cf <- coef(m)
  expect_true(all(cf[grep("Chol:(\\(Intercept\\)_\\(Intercept\\)|x_x)",
                          names(cf))] > 0))
  expect_true(any(unlist(m$model$rpar$signs) == -1))
  spec <- flexCountReg:::.fcr_fm_setup(m$model$settings, d)
  expect_equal(sum(flexCountReg:::.fcr_fm_loglik_units(cf, spec)),
               m$model$maximum, tolerance = 1e-10)
  g <- colSums(flexCountReg:::.fcr_fm_grad_units(cf, spec))
  expect_lt(max(abs(g)), 1e-3)
  expect_equal(max(abs(g)), m$model$numerical$max_abs_gradient,
               tolerance = 1e-6)
  expect_equal(unname(flexCountReg:::.fcr_fm_eval(cf, spec)$tau),
               unname(m$model$posterior), tolerance = 1e-10)
})

test_that("bootstrap refits report positive standard deviations", {
  skip_on_cran()
  d <- sim_lcrp(G = 60, sd = c(0.5, 0))
  m <- countreg_lc(y ~ x, data = d, classes = 2, group_var = "id",
                   rpar_formula = ~ x, rpar_draws = 15, n_starts = 1,
                   seed = 2, max.iters = 200)
  spec <- flexCountReg:::.fcr_fm_setup(m$model$settings, d)
  spec$st <- m$model$settings
  b <- suppressWarnings(flexCountReg:::.fcr_fm_bootstrap(
    spec, coef(m), m$model$hessian, d, 4, 5))
  sds <- b$draws[, grep("SD:", colnames(b$draws)), drop = FALSE]
  sds <- sds[stats::complete.cases(sds), , drop = FALSE]
  expect_gt(nrow(sds), 0)
  expect_true(all(sds >= 0))
})

test_that("component means over random parameters are exact or Inf", {
  d <- sim_lcrp(G = 40, seed = 4)
  d$x <- stats::runif(nrow(d))
  nd <- data.frame(x = c(-1, -0.5, 0.5, 1, 2), id = 1:5, y = 1)
  b0 <- 0.2
  mu <- 0.5
  s <- 0.3
  for (dist in c("n", "ln", "u", "t", "g")) {
    m <- suppressWarnings(countreg_lc(
      y ~ x, data = d, classes = 2, group_var = "id", rpar_formula = ~ x,
      rpardists = dist, rpar_draws = 10, n_starts = 1, seed = 2,
      max.iters = 3))
    th <- coef(m)
    th[c("Component1:Mean:(Intercept)", "Component1:Mean:x",
         "Component1:SD:x")] <- c(b0, mu, s)
    m$model$estimate <- th
    exact <- vapply(nd$x, function(x) switch(
      dist,
      n = exp(b0 + mu * x + 0.5 * s^2 * x^2),
      ln = if (x > 0) Inf else exp(b0) * stats::integrate(function(z)
        exp(x * exp(mu + s * z)) * stats::dnorm(z), -Inf, Inf,
        rel.tol = 1e-12)$value,
      g = exp(b0) * (1 - x * s^2 / mu)^(-mu^2 / s^2),
      u = exp(b0 + mu * x) * sinh(x * s) / (x * s),
      t = exp(b0 + x * (mu - s)) * expm1(x * s)^2 / (x * s)^2),
      numeric(1))
    if (dist == "ln") {
      expect_warning(M <- predict(m, newdata = nd, type = "component"),
                     "lognormal random parameter with a positive x")
    } else {
      M <- predict(m, newdata = nd, type = "component")
    }
    expect_equal(unname(M[, 1]), exact, tolerance = 1e-8, info = dist)
  }
  # Gamma: the mean exists for x sd^2 / mean < 1, the variance for
  # 2 x sd^2 / mean < 1
  th[c("Component1:Mean:x", "Component1:SD:x")] <- c(0.5, 0.8)
  m$model$estimate <- th
  expect_warning(M <- predict(m, newdata = nd, type = "component"),
                 "infinite")
  expect_equal(is.infinite(M[, 1]), nd$x * 0.64 / 0.5 >= 1)
  V <- suppressWarnings(predict(m, newdata = nd, type = "variance"))
  expect_equal(is.infinite(V), nd$x * 2 * 0.64 / 0.5 >= 1)
  expect_true(all(V[is.finite(V)] > 0))
  # A lognormal slope on a positive variable: no expected counts
  expect_warning(
    ml <- countreg_lc(y ~ x, data = d, classes = 2, group_var = "id",
                      rpar_formula = ~ x, rpardists = "ln", rpar_draws = 10,
                      n_starts = 1, seed = 2, max.iters = 20),
    "mean_count")
  expect_true(all(is.infinite(ml$model$classes$mean_count)))
  expect_true(all(is.infinite(suppressWarnings(fitted(ml)))))
  expect_error(suppressWarnings(residuals(ml, type = "pearson")),
               "second moment")
})

test_that("the lognormal moment generating function is accurate", {
  f <- flexCountReg:::.fcr_fm_ln_log_mgf
  for (p in list(c(-0.5, 0.3, 0.8), c(-3, 1, 1.5), c(-0.01, -2, 0.1),
                 c(-50, 2, 2), c(-1e-6, 0, 3))) {
    ex <- log(stats::integrate(function(z) exp(p[1] * exp(p[2] + p[3] * z)) *
                                 stats::dnorm(z), -Inf, Inf,
                               rel.tol = 1e-12)$value)
    expect_equal(exp(f(p[1], p[2], p[3])), exp(ex), tolerance = 1e-8)
  }
  expect_equal(f(c(0, 1, -1), 0.5, c(0.3, 0.3, 0)), c(0, Inf, -exp(0.5)))
})

test_that("posterior probabilities of training units match the fit", {
  d <- sim_lcrp(G = 50, seed = 6)
  m <- countreg_lc(y ~ x, data = d, classes = 2, group_var = "id",
                   rpar_formula = ~ x, rpar_draws = 15, n_starts = 1,
                   seed = 1)
  p_fit <- component_probs(m, conditioning = "posterior")
  p_new <- component_probs(m, newdata = d, conditioning = "posterior")
  expect_equal(as.vector(p_new[!duplicated(d$id), ]), as.vector(p_fit),
               tolerance = 1e-10)
  expect_equal(predict(m, conditioning = "posterior"),
               predict(m, newdata = d, conditioning = "posterior"))
  # The same panels in a different row order
  o <- rev(seq_len(nrow(d)))
  p_rev <- component_probs(m, newdata = d[o, ], conditioning = "posterior")
  expect_equal(as.vector(p_rev[order(o), ]), as.vector(p_new),
               tolerance = 1e-10)
  # Cross-sectional: the stored training rows (predict() without newdata)
  dc <- d[1:120, ]
  mc <- countreg_lc(y ~ x, data = dc, classes = 2, rpar_formula = ~ x,
                    rpar_draws = 15, n_starts = 1, seed = 1)
  P <- exp(flexCountReg:::.fcr_fm_log_class(
    mc, flexCountReg:::.fcr_fm_rows(mc, dc), dc, "posterior",
    train = TRUE))
  expect_equal(unname(P), unname(mc$model$posterior), tolerance = 1e-10)
})
