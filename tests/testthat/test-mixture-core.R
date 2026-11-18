# Finite mixtures and latent classes: likelihood, estimation engine and
# parameter bookkeeping. Reference values are computed independently with
# base R densities.

mix_pois_data <- function(n1 = 300, n2 = 200, l1 = 1, l2 = 6, seed = 1) {
  set.seed(seed)
  data.frame(y = c(stats::rpois(n1, l1), stats::rpois(n2, l2)))
}

mix_reg_data <- function(n = 600, seed = 2) {
  set.seed(seed)
  x <- stats::rnorm(n)
  z <- stats::rbinom(n, 1, 0.5)
  e <- stats::runif(n, 0.5, 2)
  cl <- ifelse(stats::runif(n) < stats::plogis(-0.5 + 1.2 * z), 1, 2)
  mu <- ifelse(cl == 1, exp(0.2 + 0.3 * x), exp(1.5 - 0.4 * x)) * e
  data.frame(y = stats::rnbinom(n, size = 1 / 0.3, mu = mu), x = x, z = z,
             e = e)
}

# Numerical gradient by central differences
num_grad <- function(f, p, h = 1e-5) {
  vapply(seq_along(p), function(j) {
    a <- p
    b <- p
    a[j] <- a[j] + h
    b[j] <- b[j] - h
    (f(a) - f(b)) / (2 * h)
  }, numeric(1))
}

fit_spec <- function(fit, data) {
  st <- fit$model$settings
  sp <- flexCountReg:::.fcr_fm_setup(st, data)
  sp$st <- st
  sp
}

test_that("two-component Poisson mixture matches a hand-computed likelihood", {
  d <- mix_pois_data()
  m <- countreg_mix(y ~ 1, data = d, family = "POISSON", components = 2,
                    n_starts = 3, seed = 42)
  cf <- coef(m)
  p1 <- stats::plogis(cf[["Mixing:Component1:(Intercept)"]])
  l1 <- exp(cf[["Component1:Mean:(Intercept)"]])
  l2 <- exp(cf[["Component2:Mean:(Intercept)"]])
  f1 <- p1 * stats::dpois(d$y, l1)
  f2 <- (1 - p1) * stats::dpois(d$y, l2)
  expect_equal(as.numeric(logLik(m)), sum(log(f1 + f2)), tolerance = 1e-10)
  # Posterior probabilities
  post <- component_probs(m, conditioning = "posterior")
  expect_equal(unname(post[, 1]), f1 / (f1 + f2), tolerance = 1e-8)
  # Recovery and display order (increasing component means)
  expect_lt(l1, l2)
  expect_equal(c(l1, l2), c(1, 6), tolerance = 0.15)
  expect_equal(p1, 0.6, tolerance = 0.1)
  expect_true(m$model$converged)
  expect_lt(m$model$numerical$max_abs_gradient, 1e-3)
  expect_identical(m$model$numerical$information, "positive definite")
})

test_that("countreg_mix() and countreg_lc() give the same fit", {
  d <- mix_reg_data()
  a <- countreg_mix(y ~ x, data = d, family = "NB2", components = 2,
                    mixing_formula = ~ z, n_starts = 3, seed = 9)
  b <- countreg_lc(y ~ x, data = d, family = "NB2", classes = 2,
                   membership_formula = ~ z, n_starts = 3, seed = 9)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)))
  expect_equal(coef(a), coef(b))
  expect_equal(vcov(a), vcov(b))
  expect_equal(predict(a, newdata = d[1:5, ]), predict(b, newdata = d[1:5, ]))
  expect_equal(unname(component_probs(a)), unname(class_probs(b)))
  expect_s3_class(b, c("countreg_lc", "countreg_mix", "flexCountReg"),
                  exact = TRUE)
  expect_s3_class(a, c("countreg_mix", "flexCountReg"), exact = TRUE)
  expect_identical(b$model$labels, c("Class1", "Class2"))
})

test_that("one component reproduces countreg()", {
  d <- mix_reg_data(300)
  for (fam in c("POISSON", "NB2", "NB1", "PIG2")) {
    c1 <- countreg(y ~ x + offset(log(e)), data = d, family = fam,
                   method = "BFGS")
    m1 <- suppressMessages(countreg_mix(y ~ x + offset(log(e)), data = d,
                                        family = fam, components = 1))
    expect_equal(as.numeric(logLik(m1)), as.numeric(logLik(c1)),
                 tolerance = 1e-6, label = fam)
    expect_equal(unname(coef(m1)), unname(coef(c1)), tolerance = 1e-3,
                 label = fam)
    expect_equal(predict(m1, newdata = d[1:4, ]),
                 predict(c1, newdata = d[1:4, ]), tolerance = 1e-3,
                 label = fam)
    expect_false(any(grepl("^Mixing", names(coef(m1)))))
  }
})

test_that("analytic scores agree with numerical derivatives", {
  d <- mix_reg_data(300)
  for (sh in list(NULL, list(mean = "x", dis1 = TRUE))) {
    m <- countreg_mix(y ~ x + offset(log(e)), data = d, family = "NB2",
                      components = 2, mixing_formula = ~ z, shared = sh,
                      n_starts = 2, seed = 1)
    sp <- fit_spec(m, d)
    th <- coef(m) + 0.03
    ga <- colSums(flexCountReg:::.fcr_fm_grad_units(th, sp))
    gn <- num_grad(function(p)
      sum(flexCountReg:::.fcr_fm_loglik_units(p, sp)), th)
    expect_equal(unname(ga), gn, tolerance = 1e-5)
  }
  # Poisson + NB1 components, panel structure
  set.seed(3)
  dp <- data.frame(id = rep(1:40, each = 3), x = stats::rnorm(120),
                   y = stats::rpois(120, 2))
  m <- countreg_mix(y ~ x, data = dp, family = c("POISSON", "NB1"),
                    group_var = "id", n_starts = 2, seed = 1)
  sp <- fit_spec(m, dp)
  th <- coef(m) + 0.02
  ga <- colSums(flexCountReg:::.fcr_fm_grad_units(th, sp))
  gn <- num_grad(function(p)
    sum(flexCountReg:::.fcr_fm_loglik_units(p, sp)), th)
  expect_equal(unname(ga), gn, tolerance = 1e-5)
})

test_that("one EM step matches the closed-form Poisson update with exposure", {
  set.seed(4)
  d <- data.frame(t = stats::runif(200, 0.5, 3))
  d$y <- stats::rpois(200, d$t * ifelse(stats::runif(200) < 0.5, 0.5, 3))
  d$w <- sample(1:3, 200, replace = TRUE)
  m <- countreg_mix(y ~ 1 + offset(log(t)), data = d, components = 2,
                    weights = "w", n_starts = 1, seed = 1)
  sp <- fit_spec(m, d)
  th0 <- c(log(0.8), log(2), 0.3)
  names(th0) <- names(coef(m))
  e <- flexCountReg:::.fcr_fm_eval(th0, sp)
  tau <- e$tau
  w <- d$w
  lam <- colSums(w * tau * d$y) / colSums(w * tau * d$t)
  pk <- colSums(w * tau) / sum(w)
  th1 <- flexCountReg:::.fcr_fm_mstep_mix(th0, sp, tau, 50L)
  th1 <- flexCountReg:::.fcr_fm_mstep_comp(th1, sp, tau, 200L)
  expect_equal(unname(exp(th1[1:2])), unname(lam), tolerance = 1e-6)
  expect_equal(unname(th1[3]), log(pk[1] / pk[2]), tolerance = 1e-10)
})

test_that("frequency tables and expanded rows give the same fit", {
  d <- mix_pois_data(seed = 5)
  tb <- as.data.frame(table(y = d$y))
  tb$y <- as.numeric(as.character(tb$y))
  a <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 3, seed = 1)
  b <- countreg_mix(y ~ 1, data = tb, weights = "Freq", components = 2,
                    n_starts = 3, seed = 1)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-8)
  expect_equal(coef(a), coef(b), tolerance = 1e-5)
  expect_equal(BIC(a), BIC(b), tolerance = 1e-8)
  expect_equal(attr(logLik(b), "nobs"), nrow(d))
  expect_equal(a$model$entropy, b$model$entropy, tolerance = 1e-5)
})

test_that("panel likelihood multiplies within a panel before mixing", {
  set.seed(6)
  G <- 60
  id <- rep(seq_len(G), each = 4)
  cl <- ifelse(stats::runif(G) < 0.5, 1, 2)[id]
  d <- data.frame(id = id, y = stats::rpois(4 * G, ifelse(cl == 1, 1, 5)))
  p <- countreg_mix(y ~ 1, data = d, components = 2, group_var = "id",
                    n_starts = 3, seed = 2)
  cf <- coef(p)
  p1 <- stats::plogis(cf[["Mixing:Component1:(Intercept)"]])
  l <- exp(cf[1:2])
  ll_panel <- sum(vapply(split(d$y, d$id), function(y)
    log(p1 * prod(stats::dpois(y, l[1])) +
          (1 - p1) * prod(stats::dpois(y, l[2]))), numeric(1)))
  ll_rows <- sum(log(p1 * stats::dpois(d$y, l[1]) +
                       (1 - p1) * stats::dpois(d$y, l[2])))
  expect_equal(as.numeric(logLik(p)), ll_panel, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(ll_panel, ll_rows)))
  expect_equal(attr(logLik(p), "nobs"), G)
  expect_equal(p$model$n_groups, G)
  expect_equal(nrow(component_probs(p)), G)
  # Duplicating every panel equals a frequency weight of 2
  d2 <- rbind(d, transform(d, id = id + G))
  dw <- transform(d, w = 2)
  a <- countreg_mix(y ~ 1, data = d2, components = 2, group_var = "id",
                    n_starts = 2, seed = 2)
  b <- countreg_mix(y ~ 1, data = dw, components = 2, group_var = "id",
                    weights = "w", n_starts = 2, seed = 2)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-7)
  expect_equal(coef(a), coef(b), tolerance = 1e-4)
})

test_that("panel inputs must be constant within panels", {
  d <- data.frame(id = rep(1:10, each = 2), y = rpois(20, 2),
                  z = rep(0:1, 10), w = rep(1:2, 10))
  expect_error(countreg_lc(y ~ 1, data = d, membership_formula = ~ z,
                           group_var = "id", n_starts = 1),
               "constant within each panel")
  expect_error(countreg_lc(y ~ 1, data = d, weights = "w", group_var = "id",
                           n_starts = 1),
               "constant within each panel")
})

test_that("free-parameter counts follow the specification", {
  d <- mix_pois_data(seed = 7)
  d$z <- rep(0:1, length.out = nrow(d))
  # (Some of these over-specified models are at a boundary; only the
  # parameter count matters here.)
  np <- function(...) length(coef(suppressWarnings(suppressMessages(
    countreg_mix(y ~ 1, data = d, n_starts = 1, seed = 1, ...)))))
  for (K in 1:3) {
    expect_equal(np(components = K), 2 * K - 1)
    expect_equal(np(components = K, family = "NB2"), 3 * K - 1)
  }
  expect_equal(np(components = 3, family = "NB2",
                  shared = list(dis1 = TRUE)), 2 * 3)
  expect_equal(np(components = 2, fixed_mixing = c(0.3, 0.7)), 2)
  expect_equal(np(components = 3, mixing_formula = ~ z), 3 + 2 * 2)
  m <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 1, seed = 1)
  expect_equal(attr(logLik(m), "df"), 3)
  expect_equal(AIC(m), -2 * as.numeric(logLik(m)) + 2 * 3)
})

test_that("relabeling components leaves the likelihood and predictions unchanged", {
  d <- mix_reg_data(300)
  m <- countreg_mix(y ~ x, data = d, family = "NB2", components = 2,
                    mixing_formula = ~ z, n_starts = 2, seed = 1)
  sp <- fit_spec(m, d)
  th <- coef(m)
  thp <- flexCountReg:::.fcr_fm_permute(th, m$model$map, c(2L, 1L))
  expect_equal(sum(flexCountReg:::.fcr_fm_loglik_units(thp, sp)),
               as.numeric(logLik(m)), tolerance = 1e-10)
  mp <- m
  mp$model$estimate <- thp
  expect_equal(predict(mp, newdata = d[1:6, ]), predict(m, newdata = d[1:6, ]))
  # Permuting twice is the identity; the matrix form agrees
  expect_equal(flexCountReg:::.fcr_fm_permute(thp, m$model$map, c(2L, 1L)),
               th)
  A <- flexCountReg:::.fcr_fm_perm_matrix(m$model$map, c(2L, 1L))
  expect_equal(as.vector(A %*% th), unname(thp))
  # Components of different families cannot be swapped
  mf <- countreg_mix(y ~ 1, data = d, family = c("POISSON", "NB2"),
                     n_starts = 1, seed = 1)
  expect_error(flexCountReg:::.fcr_fm_permute(coef(mf), mf$model$map,
                                              c(2L, 1L)), "same family")
})

test_that("mixed-family mixtures match a hand-computed likelihood", {
  d <- mix_reg_data(400)
  m <- countreg_mix(y ~ x, data = d,
                    family = c(base = "POISSON", tail = "NB2"),
                    n_starts = 3, seed = 3)
  expect_identical(m$model$labels, c("base", "tail"))
  expect_identical(m$model$family, c("POISSON", "NB2"))
  cf <- coef(m)
  expect_setequal(names(cf), c(
    "Component1:Mean:(Intercept)", "Component1:Mean:x",
    "Component2:Mean:(Intercept)", "Component2:Mean:x",
    "Component2:ln(alpha)", "Mixing:Component1:(Intercept)"))
  p1 <- stats::plogis(cf[["Mixing:Component1:(Intercept)"]])
  m1 <- exp(cf[["Component1:Mean:(Intercept)"]] +
              cf[["Component1:Mean:x"]] * d$x)
  m2 <- exp(cf[["Component2:Mean:(Intercept)"]] +
              cf[["Component2:Mean:x"]] * d$x)
  a <- exp(cf[["Component2:ln(alpha)"]])
  ll <- sum(log(p1 * stats::dpois(d$y, m1) +
                  (1 - p1) * stats::dnbinom(d$y, size = 1 / a, mu = m2)))
  expect_equal(as.numeric(logLik(m)), ll, tolerance = 1e-10)
  expect_error(countreg_mix(y ~ x, data = d, family = c("POISSON", "NB2"),
                            shared = list(dis1 = TRUE)),
               "not every component")
  expect_error(countreg_mix(y ~ x, data = d, family = c("NB1", "NB2"),
                            shared = list(dis1 = TRUE)),
               "same family")
})

test_that("shared terms are estimated once", {
  d <- mix_reg_data(300)
  d$f <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
  m <- countreg_mix(y ~ x + f, data = d, family = "NB2", components = 2,
                    shared = list(mean = c("x", "f"), dis1 = TRUE),
                    n_starts = 2, seed = 1)
  cf <- names(coef(m))
  expect_true(all(c("Shared:Mean:x", "Shared:Mean:fb", "Shared:Mean:fc",
                    "Shared:ln(alpha)") %in% cf))
  expect_false(any(grepl("^Component[12]:Mean:(x|f)", cf)))
  expect_error(countreg_mix(y ~ x + f, data = d, components = 2,
                            shared = list(mean = "fb")),
               "not a term of the formula")
  expect_error(countreg_mix(y ~ x, data = d, components = 2,
                            shared = list(mean = c("(Intercept)", "x"))),
               "identical")
})

test_that("known mixing probabilities stay fixed", {
  d <- mix_pois_data(seed = 8)
  m <- countreg_mix(y ~ 1, data = d, components = 2,
                    fixed_mixing = c(0.7, 0.3), n_starts = 3, seed = 1)
  expect_false(any(grepl("^Mixing", names(coef(m)))))
  expect_equal(unname(m$model$prior[1, ]), c(0.7, 0.3))
  expect_equal(unname(component_probs(m, newdata = d[1:2, , drop = FALSE])[1, ]),
               c(0.7, 0.3))
  # The order of fixed components with different probabilities is kept
  cf <- coef(m)
  expect_equal(as.numeric(logLik(m)),
               sum(log(0.7 * stats::dpois(d$y, exp(cf[[1]])) +
                         0.3 * stats::dpois(d$y, exp(cf[[2]])))),
               tolerance = 1e-10)
  expect_error(countreg_mix(y ~ 1, data = d, fixed_mixing = c(0.5, 0.6)),
               "sum to 1")
  expect_error(countreg_mix(y ~ 1, data = d, fixed_mixing = c(1, 0)),
               "zero probability")
  expect_error(countreg_mix(y ~ 1, data = d, fixed_mixing = 1),
               "2 finite")
})

test_that("invalid specifications are rejected before fitting", {
  d <- mix_pois_data(seed = 9)
  expect_error(countreg_mix(y ~ 1, data = d, family = "POISSON",
                            dis_param_formula_1 = ~ 1),
               "no component family")
  expect_error(countreg_mix(y ~ 1, data = d, components = 0), "at least 1")
  expect_error(countreg_mix(y ~ 1, data = d, family = c("POISSON", "NB2"),
                            components = 3), "3 components")
  expect_error(countreg_mix(y ~ 1, data = d, family = "FOO"), "Unknown")
  expect_error(countreg_mix(y ~ 1, data = d, mixing_formula = y ~ 1),
               "one-sided")
  expect_error(countreg_mix(y ~ 1, data = d, shared = list(slopes = "x")),
               "Unknown element")
  expect_error(countreg_mix(y ~ 1, data = d, n_starts = 0), "n_starts")
  expect_error(countreg_mix(y ~ 1, data = d, control = list()),
               "mixture_control")
  expect_error(mixture_control(em_tol = -1), "positive")
})

test_that("the start table records every start", {
  d <- mix_reg_data(300)
  m <- countreg_mix(y ~ x, data = d, components = 2, n_starts = 4,
                    seed = 11)
  st <- m$model$starts
  expect_equal(nrow(st), 4)
  expect_equal(sum(st$selected), 1)
  expect_true(st$reached_best[st$selected])
  expect_true(all(c("strategy", "seed", "em_iterations", "em_status",
                    "loglik_start", "loglik", "converged") %in% names(st)))
  # Reproducible from the seed, without changing the caller's stream
  set.seed(123)
  r1 <- stats::runif(1)
  set.seed(123)
  m2 <- countreg_mix(y ~ x, data = d, components = 2, n_starts = 4,
                     seed = 11)
  expect_equal(stats::runif(1), r1)
  expect_equal(coef(m2), coef(m))
  expect_identical(m2$model$starts$seed, st$seed)
  # User starting values are the first start
  m3 <- countreg_mix(y ~ x, data = d, components = 2, n_starts = 2,
                     seed = 1, start.vals = coef(m))
  expect_identical(m3$model$starts$strategy[1], "user")
  expect_equal(as.numeric(logLik(m3)), as.numeric(logLik(m)),
               tolerance = 1e-6)
  expect_error(countreg_mix(y ~ x, data = d, components = 2,
                            start.vals = c(1, 2)), "start.vals")
})

test_that("the direct algorithm reaches the same maximum", {
  d <- mix_pois_data(seed = 12)
  a <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 3, seed = 1)
  b <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 3, seed = 1,
                    control = mixture_control(algorithm = "direct"))
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-6)
})

test_that("robust and bootstrap standard errors are available", {
  skip_on_cran()
  d <- mix_pois_data(n1 = 150, n2 = 100, seed = 13)
  m <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 2, seed = 1)
  r <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 2, seed = 1,
                    stderr = "robust")
  expect_identical(r$model$vcov_type, "robust")
  expect_equal(coef(r), coef(m))
  expect_true(all(is.finite(r$model$se)))
  b <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 2, seed = 1,
                    stderr = "boot", bootstraps = 10)
  expect_identical(b$model$vcov_type, "bootstrap")
  expect_true(b$model$n_boot_ok >= 8)
  expect_true(all(is.finite(b$model$se)))
  # Bootstrap SEs are of the same order as the Hessian SEs
  expect_true(all(b$model$se / m$model$se > 0.3 & b$model$se /
                    m$model$se < 3))
})

test_that("nearly identical components are detected", {
  d <- mix_pois_data(seed = 14)
  m <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 1, seed = 1)
  sp <- fit_spec(m, d)
  th <- coef(m)
  expect_length(flexCountReg:::.fcr_fm_identical(th, sp), 0)
  th[2] <- th[1]
  expect_identical(flexCountReg:::.fcr_fm_identical(th, sp), "1-2")
  # A single-component Poisson sample fitted with two components
  set.seed(15)
  d1 <- data.frame(y = stats::rpois(300, 3))
  f <- suppressWarnings(countreg_mix(y ~ 1, data = d1, components = 2,
                                     n_starts = 2, seed = 1))
  if (length(f$model$identical_pairs)) {
    expect_output(print(summary(f)), "nearly identical")
  }
  expect_true(is.character(f$model$identical_pairs))
})

test_that("one-component likelihood equals countreg() for every family", {
  skip_on_cran()
  set.seed(11)
  n <- 120
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$y <- stats::rnbinom(n, size = 2, mu = exp(1 + 0.4 * d$x))
  # PLG is left out only because its density is slow at these settings
  for (f in setdiff(flexCountReg:::.fcr_families, "PLG")) {
    f1 <- if (!is.null(flexCountReg:::get_params(f)[[1]])) ~ z
    c0 <- suppressWarnings(countreg(y ~ x, data = d, family = f,
                                    method = "BFGS", max.iters = 0,
                                    dis_param_formula_1 = f1))
    th <- coef(c0) + 0.05
    c1 <- suppressWarnings(countreg(y ~ x, data = d, family = f,
                                    method = "BFGS", max.iters = 0,
                                    dis_param_formula_1 = f1,
                                    start.vals = th))
    st <- flexCountReg:::.fcr_fm_settings(list(
      entry = "countreg_mix", formula = y ~ x, family = f, components = 1,
      components_missing = FALSE, mixing_formula = ~ 1,
      fixed_mixing = NULL, shared = NULL, dis_param_formula_1 = f1,
      dis_param_formula_2 = NULL, offset = NULL, weights = NULL,
      group_var = NULL, method = "BFGS", max.iters = 0, n_starts = 1,
      seed = 1, ndraws = 1500, engine = "quadrature",
      control = mixture_control()))
    sp <- flexCountReg:::.fcr_fm_setup(st, d)
    sp$st <- st
    ll <- sum(flexCountReg:::.fcr_fm_loglik_units(unname(th), sp))
    expect_equal(ll, c1$model$maximum, tolerance = 1e-10, label = f)
  }
})
