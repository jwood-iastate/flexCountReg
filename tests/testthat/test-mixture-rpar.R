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
  # slope, exp(b0 + b1 x + s^2 x^2 / 2), up to simulation error
  cf <- coef(m)
  M <- predict(m, newdata = d[1:5, ], type = "component")
  ex <- exp(cf["Component1:Mean:(Intercept)"] +
              cf["Component1:Mean:x"] * d$x[1:5] +
              cf["Component1:SD:x"]^2 * d$x[1:5]^2 / 2)
  expect_equal(unname(M[, 1]), unname(ex), tolerance = 0.02)
  pr <- predict(m, newdata = d[1:3, ], type = "prob", max_count = 150)
  expect_equal(rowSums(pr), rep(1, 3), tolerance = 1e-6)
  # Probabilities and expected counts agree
  expect_equal(as.vector(pr %*% (0:150)),
               predict(m, newdata = d[1:3, ]), tolerance = 1e-4)
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
})
