# Finite mixtures: fixes from the external review (undefined GP2
# components, panel membership designs, tail probabilities, quantile
# residuals in log space, relabeling by optimal assignment).

# Training problem of a two-component mixture with fixed 50:50 mixing, so
# that theta holds only the component parameters.
rf_spec <- function(d, formula, families, f1 = NULL, weights = NULL) {
  st <- flexCountReg:::.fcr_fm_settings(list(
    entry = "countreg_mix", formula = formula, family = families,
    components = length(families), components_missing = FALSE,
    mixing_formula = ~ 1, fixed_mixing = c(0.5, 0.5), shared = NULL,
    dis_param_formula_1 = f1, dis_param_formula_2 = NULL, offset = NULL,
    weights = weights, group_var = NULL, method = "BFGS", max.iters = 0,
    n_starts = 1, seed = 1, ndraws = 1500, engine = "quadrature",
    control = mixture_control()))
  sp <- flexCountReg:::.fcr_fm_setup(st, d)
  sp$st <- st
  sp
}

rf_theta <- function(sp, values) {
  th <- setNames(rep(NA_real_, sp$map$n_par), sp$map$names)
  th[names(values)] <- values
  stopifnot(!anyNA(th))
  th
}

rf_ll <- function(th, sp) {
  sum(flexCountReg:::.fcr_fm_loglik_units(unname(th), sp))
}

test_that("an undefined GP2 component cannot hide behind a valid one", {
  # GP2 at internal mu = 5, alpha = 0.5 has 1 + (alpha - 1) mu = -1.5 < 0:
  # no distribution. Its registry density is quietly 0, so a nominal 50:50
  # mixture with a Poisson sums to 0.5 over all counts.
  dens <- flexCountReg:::get_probFunc("GP2")
  y <- 0:200
  f <- 0.5 * dens(y, predicted = rep(5, 201), alpha = rep(0.5, 201),
                  sigma = NULL) + 0.5 * dpois(y, 2)
  expect_equal(sum(f), 0.5)
  expect_false(flexCountReg:::.fcr_family_valid("GP2", 5, 0.5))

  d <- data.frame(y = c(0, 1, 2, 3, 4, 5, 2, 1))
  sp <- rf_spec(d, y ~ 1, c("GP2", "POISSON"))
  bad <- rf_theta(sp, c("Component1:Mean:(Intercept)" = log(5),
                        "Component1:ln(alpha+1)" = log(0.5),
                        "Component2:Mean:(Intercept)" = log(2)))
  expect_null(flexCountReg:::.fcr_fm_eval(unname(bad), sp))
  expect_true(all(is.na(
    flexCountReg:::.fcr_fm_loglik_units(unname(bad), sp))))
  # The EM objective rejects it too
  tau <- matrix(0.5, nrow(d), 2)
  expect_true(is.na(flexCountReg:::.fcr_fm_Q(
    unname(bad), seq_along(bad), unname(bad), sp, tau, 1:2)))
})

test_that("a valid underdispersed GP2 keeps a zero for impossible counts", {
  # alpha = 0.9 (a = -0.1) at mu = 5: 1 + a mu = 0.5 > 0, support y < 10.
  # The count 15 is impossible under the GP2 component but the component is
  # defined, so the mixture is valid and the Poisson explains that row.
  d <- data.frame(y = c(0, 3, 5, 6, 4, 15))
  sp <- rf_spec(d, y ~ 1, c("GP2", "POISSON"))
  th <- rf_theta(sp, c("Component1:Mean:(Intercept)" = log(5),
                       "Component1:ln(alpha+1)" = log(0.9),
                       "Component2:Mean:(Intercept)" = log(6)))
  e <- flexCountReg:::.fcr_fm_eval(unname(th), sp)
  expect_false(is.null(e))
  expect_identical(e$LP[6, 1], -Inf)
  expect_equal(e$tau[6, ], c(0, 1))
  dens <- flexCountReg:::get_probFunc("GP2")
  g <- dens(d$y, predicted = rep(5, 6), alpha = rep(0.9, 6), sigma = NULL)
  expect_equal(rf_ll(th, sp), sum(log(0.5 * g + 0.5 * dpois(d$y, 6))))
  expect_true(is.finite(rf_ll(th, sp)))
})

test_that("GP2 validity is checked only on rows with positive weight", {
  # Mean regression: x = 1 gives mu = 5 (undefined with alpha = 0.5),
  # x = 0 gives mu = 1 (1 - 0.5 = 0.5 > 0)
  d <- data.frame(y = c(0, 1, 2, 1, 4, 3), x = c(0, 0, 0, 0, 1, 1),
                  w = c(1, 1, 2, 1, 0, 0))
  vals <- c("Component1:Mean:(Intercept)" = 0,
            "Component1:Mean:x" = log(5),
            "Component1:ln(alpha+1)" = log(0.5),
            "Component2:Mean:(Intercept)" = log(2),
            "Component2:Mean:x" = 0)
  sp0 <- rf_spec(d, y ~ x, c("GP2", "POISSON"), weights = "w")
  th <- rf_theta(sp0, vals)
  ll0 <- rf_ll(th, sp0)
  expect_true(is.finite(ll0))
  # Zero-weight rows leave the log-likelihood unchanged
  dens <- flexCountReg:::get_probFunc("GP2")
  i <- 1:4
  f <- 0.5 * dens(d$y[i], predicted = rep(1, 4), alpha = rep(0.5, 4),
                  sigma = NULL) + 0.5 * dpois(d$y[i], 2)
  expect_equal(ll0, sum(d$w[i] * log(f)))
  d1 <- d
  d1$w[5] <- 1
  sp1 <- rf_spec(d1, y ~ x, c("GP2", "POISSON"), weights = "w")
  expect_true(is.na(rf_ll(th, sp1)))
})

test_that("GP2 validity follows a dispersion regression", {
  # alpha = exp(d0 + d1 z): 0.9 when z = 0 (defined at mu = 5), 0.5 when
  # z = 1 (undefined)
  d <- data.frame(y = c(3, 5, 6, 4, 2, 7), z = c(0, 0, 0, 0, 1, 1),
                  w = c(1, 1, 1, 1, 0, 0))
  vals <- c("Component1:Mean:(Intercept)" = log(5),
            "Component1:ln(alpha+1):(Intercept)" = log(0.9),
            "Component1:ln(alpha+1):z" = log(0.5 / 0.9),
            "Component2:Mean:(Intercept)" = log(4))
  sp0 <- rf_spec(d, y ~ 1, c("GP2", "POISSON"), f1 = ~ z, weights = "w")
  th <- rf_theta(sp0, vals)
  expect_true(is.finite(rf_ll(th, sp0)))
  d$w[] <- 1
  sp1 <- rf_spec(d, y ~ 1, c("GP2", "POISSON"), f1 = ~ z, weights = "w")
  expect_true(is.na(rf_ll(th, sp1)))
})

test_that("predictions report NA where a GP2 component is undefined", {
  set.seed(3)
  d <- data.frame(x = runif(80))
  d$y <- rpois(80, exp(0.2 + 0.3 * d$x))
  m <- suppressWarnings(suppressMessages(countreg_mix(
    y ~ x, data = d, family = c("GP2", "POISSON"),
    fixed_mixing = c(0.5, 0.5), n_starts = 1, seed = 1)))
  est <- m$model$estimate
  est[] <- 0
  est["Component1:Mean:x"] <- log(1.5)
  est["Component1:ln(alpha+1)"] <- log(0.5)
  m$model$estimate <- est
  # x = 0: mu = 1, defined; x = 3: mu = 3.375, 1 - 0.5 * 3.375 < 0
  nd <- data.frame(x = c(0, 3), y = c(1, 1))
  P <- predict(m, newdata = nd, type = "prob", max_count = 40,
               overflow = TRUE)
  expect_equal(sum(P[1, ]), 1)
  expect_true(all(is.na(P[2, ])))
  post <- component_probs(m, newdata = nd, conditioning = "posterior")
  expect_true(all(is.finite(post[1, ])))
  expect_true(all(is.na(post[2, ])))
  expect_warning(q <- flexCountReg:::.fcr_fm_quantile_residuals(
    m, nd, nd$y, u = c(0.5, 0.5)), "invalid fitted probabilities")
  expect_true(is.finite(q[1]))
  expect_true(is.na(q[2]))
})

# ---- Panel membership designs ------------------------------------------------

rf_panel <- function() {
  set.seed(31)
  G <- 60
  id <- rep(seq_len(G), each = 4)
  z <- rbinom(G, 1, 0.5)
  cl <- ifelse(runif(G) < plogis(-0.5 + 1.5 * z), 1, 2)
  d <- data.frame(id = id, t = rep(1:4, G), z = z[id])
  d$y <- rpois(nrow(d), ifelse(cl[id] == 1, 0.7, 4))
  m <- countreg_lc(y ~ 1, data = d, classes = 2, group_var = "id",
                   membership_formula = ~ z, n_starts = 2, seed = 1)
  list(d = d, m = m)
}

test_that("contradictory membership variables in a panel are an error", {
  f <- rf_panel()
  d <- f$d
  m <- f$m
  targ <- d[d$t > 2, ]
  hist <- d[d$t <= 2 & d$id <= 30, ]
  bad <- targ
  i <- which(bad$id == 1)[1]
  bad$z[i] <- 1 - bad$z[i]
  # Before the fix the first row of each panel set its prior, so the
  # result depended on the row order; now both orders are refused
  for (nd in list(bad, bad[rev(seq_len(nrow(bad))), ])) {
    expect_error(predict(m, newdata = nd, conditioning = "posterior"),
                 "constant within each panel \\(in the new data\\)")
    expect_error(predict(m, newdata = nd),
                 "constant within each panel \\(in the new data\\)")
    expect_error(predict(m, newdata = nd, conditioning = "history",
                         history = hist), "in the new data")
  }
  # Inconsistent only within the history rows of a panel
  hbad <- hist
  j <- which(hbad$id == 2)[1]
  hbad$z[j] <- 1 - hbad$z[j]
  expect_error(predict(m, newdata = targ, conditioning = "history",
                       history = hbad), "\\(in `history`\\)")
  expect_error(component_probs(m, newdata = targ, conditioning = "history",
                               history = hbad), "in `history`")
  # Consistent within each, but history and target disagree
  hb2 <- hist
  hb2$z[hb2$id == 3] <- 1 - hb2$z[hb2$id == 3]
  expect_error(predict(m, newdata = targ, conditioning = "history",
                       history = hb2), "between the new data and `history`")
  # Fitting uses the same check
  expect_error(countreg_lc(y ~ 1, data = bad, classes = 2, group_var = "id",
                           membership_formula = ~ z, n_starts = 1),
               "constant within each panel")
  # Prior predictions without the panel identifiers treat rows separately
  expect_length(predict(m, newdata = bad[, c("y", "z")]), nrow(bad))
})

test_that("panel class probabilities do not depend on the row order", {
  f <- rf_panel()
  d <- f$d
  m <- f$m
  targ <- d[d$t > 2, ]
  hist <- d[d$t <= 2 & d$id <= 30, ]
  set.seed(5)
  o1 <- sample(nrow(targ))
  o2 <- sample(nrow(hist))
  for (cond in c("prior", "posterior", "history")) {
    a <- predict(m, newdata = targ, conditioning = cond, history = hist)
    b <- predict(m, newdata = targ[o1, ], conditioning = cond,
                 history = hist[o2, ])
    expect_equal(b[order(o1)], a, label = cond)
    pa <- component_probs(m, newdata = targ, conditioning = cond,
                          history = hist)
    pb <- component_probs(m, newdata = targ[o1, ], conditioning = cond,
                          history = hist[o2, ])
    expect_equal(pb[order(o1), , drop = FALSE], pa[, , drop = FALSE],
                 ignore_attr = TRUE, label = cond)
  }
})

test_that("missing outcomes and missing history are handled explicitly", {
  f <- rf_panel()
  d <- f$d
  m <- f$m
  targ <- d[d$t > 2, ]
  hist <- d[d$t <= 2 & d$id <= 30, ]
  pp <- predict(m, newdata = targ)
  ph <- predict(m, newdata = targ, conditioning = "history", history = hist)
  # A panel whose history rows all lack the outcome keeps its prior; a
  # panel with one missing history row is conditioned on the other row
  h2 <- hist
  h2$y[h2$id == 1] <- NA
  h2$y[which(h2$id == 2)[1]] <- NA
  expect_warning(ph2 <- predict(m, newdata = targ, conditioning = "history",
                                history = h2), "3 history row")
  expect_equal(ph2[targ$id == 1], pp[targ$id == 1])
  ph3 <- predict(m, newdata = targ, conditioning = "history",
                 history = h2[!is.na(h2$y), ])
  expect_equal(ph2, ph3)
  expect_equal(ph2[targ$id > 2], ph[targ$id > 2])
  # In the posterior target a row with a missing outcome adds nothing
  t2 <- targ
  t2$y[which(t2$id == 4)[1]] <- NA
  a <- predict(m, newdata = t2, conditioning = "posterior")
  b <- predict(m, newdata = t2[!is.na(t2$y), ], conditioning = "posterior")
  expect_equal(a[!is.na(t2$y)], b)
  expect_equal(a[t2$id == 4][1], a[t2$id == 4][2])
})

# ---- Tail probabilities and quantile residuals --------------------------------

rf_two_poisson <- function(y, means, p1) {
  d <- data.frame(y = y)
  m <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 1, seed = 1)
  m$model$estimate[] <- c(log(means), log(p1 / (1 - p1)))
  m
}

test_that("the overflow probability keeps its relative accuracy", {
  set.seed(7)
  m <- rf_two_poisson(c(rpois(60, 1), rpois(40, 5)), c(1, 1), 0.6)
  P <- predict(m, newdata = data.frame(y = 0), type = "prob",
               max_count = 20, overflow = TRUE)
  # Relative accuracy (expect_equal's tolerance is absolute for values
  # this small)
  rel <- function(a, b) abs(a / b - 1)
  expect_lt(rel(P[1, ">20"], ppois(20, 1, lower.tail = FALSE)), 1e-12)
  expect_lt(rel(P[1, ">20"], 7.542625e-21), 1e-6)
  m2 <- rf_two_poisson(c(rpois(60, 1), rpois(40, 5)), c(1, 1.5), 0.6)
  P2 <- predict(m2, newdata = data.frame(y = 0), type = "prob",
                max_count = 25, overflow = TRUE)
  expect_lt(rel(P2[1, ">25"], 0.6 * ppois(25, 1, lower.tail = FALSE) +
                  0.4 * ppois(25, 1.5, lower.tail = FALSE)), 1e-12)
})

test_that("quantile residuals are computed in log space", {
  set.seed(8)
  y <- c(rpois(60, 1), rpois(39, 5), 80)
  m <- rf_two_poisson(y, c(1, 5), 0.6)
  # Independent computation at fixed u
  lse <- function(a, b) {
    m <- pmax(a, b)
    ifelse(m == -Inf, -Inf, m + log1p(exp(-abs(a - b))))
  }
  ref <- function(y, u) {
    l1 <- log(0.6)
    l2 <- log(0.4)
    lp <- lse(l1 + dpois(y, 1, log = TRUE), l2 + dpois(y, 5, log = TRUE))
    lo <- lse(l1 + ppois(y - 1, 1, log.p = TRUE),
              l2 + ppois(y - 1, 5, log.p = TRUE))
    s <- lse(l1 + ppois(y, 1, lower.tail = FALSE, log.p = TRUE),
             l2 + ppois(y, 5, lower.tail = FALSE, log.p = TRUE))
    lu <- lse(lo, log(u) + lp)
    l1mu <- lse(s, log1p(-u) + lp)
    low <- lu <= log(0.5)
    out <- qnorm(pmin(l1mu, 0), lower.tail = FALSE, log.p = TRUE)
    out[low] <- qnorm(lu[low], log.p = TRUE)
    out
  }
  set.seed(9)
  q <- residuals(m, type = "quantile")
  set.seed(9)
  u <- runif(length(y))
  expect_true(all(is.finite(q)))
  expect_equal(q, ref(y, u), tolerance = 1e-10)
  # The last observation is far in the upper tail: exponentiating the log
  # distribution function would have given qnorm(1) = Inf
  expect_gt(q[100], 15)
  expect_identical(exp(log(0.6 + 0.4 * ppois(79, 5))), 1)
  qf <- flexCountReg:::.fcr_fm_quantile_residuals(
    m, data.frame(y = c(0, 3, 80)), c(0, 3, 80), u = c(0, 0.25, 1))
  expect_equal(qf, ref(c(0, 3, 80), c(0, 0.25, 1)), tolerance = 1e-10)
})

# ---- Relabeling by optimal assignment ----------------------------------------

# Exhaustive reference: all permutations of 1..K (rows)
rf_all_perms <- function(K) {
  if (K == 1L) return(matrix(1L, 1, 1))
  sub <- rf_all_perms(K - 1L)
  do.call(rbind, lapply(seq_len(K), function(i) {
    cbind(i, matrix(setdiff(seq_len(K), i)[sub], ncol = K - 1L))
  }))
}

test_that("the assignment solver matches exhaustive enumeration", {
  set.seed(11)
  for (K in 1:6) {
    P <- rf_all_perms(K)
    for (rep in 1:15) {
      C <- matrix(rexp(K * K), K)
      if (rep > 10 && K > 2) {
        # Forbidden pairs (different families), keeping the identity allowed
        C[sample(which(row(C) != col(C)), K)] <- Inf
      }
      cost <- apply(P, 1, function(p) sum(C[cbind(seq_len(K), p)]))
      a <- flexCountReg:::.fcr_fm_lap(C, second = TRUE)
      expect_equal(a$cost, min(cost))
      expect_equal(sum(C[cbind(seq_len(K), a$perm)]), min(cost))
      expect_equal(a$second, if (K > 1) sort(cost)[2] else Inf)
    }
  }
})

rf_map <- function(families, n_z = 2) {
  X <- stats::model.matrix(~ x, data.frame(x = c(0.1, 0.5, 0.9)))
  tt <- stats::terms(~ x)
  Z <- if (n_z > 0) cbind("(Intercept)" = 1, z = c(0, 1, 0))[, seq_len(n_z),
                                                          drop = FALSE]
  D1 <- matrix(1, 3, 1, dimnames = list(NULL, "(Intercept)"))
  flexCountReg:::.fcr_fm_map(families, X, tt, D1, NULL, FALSE, FALSE, Z,
                             flexCountReg:::.fcr_fm_check_shared(
                               list(mean = "x")))
}

test_that("label alignment equals the exhaustive search over relabelings", {
  set.seed(12)
  for (fams in list(rep("POISSON", 4), c("POISSON", "NB2", "POISSON", "NB2"),
                    rep("NB2", 3))) {
    for (n_z in 0:2) {
      map <- rf_map(fams, n_z)
      K <- map$K
      P <- rf_all_perms(K)
      P <- P[apply(P, 1, function(p) identical(fams[p], fams)), ,
             drop = FALSE]
      for (rep in 1:5) {
        theta <- rnorm(map$n_par)
        est <- rnorm(map$n_par)
        scale <- runif(map$n_par, 0.5, 2)
        dist <- apply(P, 1, function(p) sum(((flexCountReg:::.fcr_fm_permute(
          est, map, p) - theta) / scale)^2))
        al <- flexCountReg:::.fcr_fm_align(est, theta, map, scale)
        expect_equal(al$dist, min(dist))
        expect_equal(al$perm, unname(P[which.min(dist), ]))
        expect_equal(al$second, if (nrow(P) > 1) sort(dist)[2] else Inf)
      }
    }
  }
  # Fixed mixing: only components with equal probabilities are exchangeable
  map <- rf_map(rep("POISSON", 3), 0)
  theta <- c(0, 1, 2, 0.3)
  est <- c(2, 1, 0, 0.3)
  al <- flexCountReg:::.fcr_fm_align(est, theta, map, rep(1, 4),
                                     fixed = c(0.25, 0.5, 0.25))
  expect_equal(al$perm, c(3L, 2L, 1L))
  expect_equal(al$dist, 0)
  expect_equal(map$n_par, 4L)
  al <- flexCountReg:::.fcr_fm_align(est, theta, map, rep(1, 4),
                                     fixed = c(0.2, 0.5, 0.3))
  expect_equal(al$perm, 1:3)
})

test_that("eight components are relabeled correctly", {
  set.seed(13)
  map <- rf_map(rep("POISSON", 8), 2)
  theta <- rnorm(map$n_par, sd = 3)
  true <- sample(8L)
  est <- flexCountReg:::.fcr_fm_permute(theta, map, true) +
    rnorm(map$n_par, sd = 0.01)
  t0 <- proc.time()[["elapsed"]]
  al <- flexCountReg:::.fcr_fm_align(est, theta, map, rep(1, map$n_par))
  expect_lt(proc.time()[["elapsed"]] - t0, 5)
  expect_equal(flexCountReg:::.fcr_fm_permute(est, map, al$perm), theta,
               tolerance = 0.05)
  expect_equal(al$perm, order(true))
  expect_gt(al$second, 1.5 * al$dist)
})

test_that("a bootstrap of an eight-component mixture aligns its labels", {
  skip_on_cran()
  # Relabeling used to enumerate all K! permutations and stopped above
  # seven components, after the fit had been computed
  set.seed(14)
  mu <- c(0.3, 1.5, 4, 9, 16, 28, 45, 70)
  d <- data.frame(y = rpois(480, rep(mu, each = 60)))
  m <- suppressWarnings(countreg_mix(
    y ~ 1, data = d, components = 8, n_starts = 1, seed = 1,
    stderr = "boot", bootstraps = 2,
    control = mixture_control(boot_starts = 1, em_maxit = 20)))
  expect_identical(m$model$K, 8L)
  expect_identical(m$model$n_boot_ok, 2L)
  expect_true(all(is.finite(sqrt(diag(vcov(m)))[1:8])))
})
