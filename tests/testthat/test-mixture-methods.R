# Finite mixtures and latent classes: predictions, probabilities,
# residuals, simulation and integration with the package's utilities.

fm_data <- function(n = 500, seed = 21) {
  set.seed(seed)
  x <- stats::rnorm(n)
  z <- stats::rbinom(n, 1, 0.5)
  cl <- ifelse(stats::runif(n) < stats::plogis(-0.3 + 1.5 * z), 1, 2)
  mu <- ifelse(cl == 1, exp(0.1 + 0.4 * x), exp(1.6 - 0.2 * x))
  data.frame(y = stats::rpois(n, mu), x = x, z = z)
}

test_that("mixture mean and variance match the analytic oracle", {
  # Poisson means 1 and 5 with probabilities 0.6 and 0.4: mean 2.6,
  # variance 0.6 * (1 + 1) + 0.4 * (5 + 25) - 2.6^2 = 6.44
  d <- data.frame(y = c(rpois(60, 1), rpois(40, 5)))
  m <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 1, seed = 1)
  m$model$estimate[] <- c(log(1), log(5), log(0.6 / 0.4))
  nd <- data.frame(y = 0)
  expect_equal(predict(m, newdata = nd), 2.6)
  expect_equal(predict(m, newdata = nd, type = "variance"), 6.44)
  expect_equal(unname(predict(m, newdata = nd, type = "component")[1, ]),
               c(1, 5))
  P <- predict(m, newdata = nd, type = "prob", max_count = 30)
  expect_equal(unname(P[1, ]), 0.6 * dpois(0:30, 1) + 0.4 * dpois(0:30, 5))
  P2 <- predict(m, newdata = nd, type = "prob", max_count = 3,
                overflow = TRUE)
  expect_equal(sum(P2), 1)
  expect_identical(colnames(P2), c("0", "1", "2", "3", ">3"))
  expect_equal(unname(predict(m, newdata = nd, type = "link")[1, ]),
               log(c(1, 5)))
})

test_that("prior predictions never use the outcomes in newdata", {
  d <- fm_data()
  m <- countreg_lc(y ~ x, data = d, classes = 2, membership_formula = ~ z,
                   n_starts = 2, seed = 1)
  nd <- d[1:20, ]
  nd2 <- nd
  nd2$y <- nd2$y + 50
  expect_equal(predict(m, newdata = nd), predict(m, newdata = nd2))
  expect_equal(component_probs(m, newdata = nd),
               component_probs(m, newdata = nd2))
  # The posterior target does use them, explicitly
  expect_false(isTRUE(all.equal(
    predict(m, newdata = nd, conditioning = "posterior"),
    predict(m, newdata = nd2, conditioning = "posterior"))))
  # Batch invariance: predicting rows separately or together
  one <- vapply(1:5, function(i) predict(m, newdata = nd[i, ]), numeric(1))
  expect_equal(one, predict(m, newdata = nd[1:5, ]))
  # Outcome not required for prior predictions
  expect_length(predict(m, newdata = nd[, c("x", "z")]), 20)
  # Mixture mean is the probability-weighted component mean
  P <- component_probs(m, newdata = nd)
  M <- predict(m, newdata = nd, type = "component")
  expect_equal(predict(m, newdata = nd), unname(rowSums(P * M)))
})

test_that("history conditioning updates only panels with history", {
  set.seed(22)
  G <- 80
  id <- rep(seq_len(G), each = 5)
  cl <- ifelse(stats::runif(G) < 0.5, 1, 2)[id]
  d <- data.frame(id = id, t = rep(1:5, G),
                  y = stats::rpois(5 * G, ifelse(cl == 1, 0.8, 4)))
  m <- countreg_lc(y ~ 1, data = d, classes = 2, group_var = "id",
                   n_starts = 2, seed = 1)
  hist <- d[d$t <= 3 & d$id <= 40, ]
  targ <- d[d$t > 3, ]
  ph <- predict(m, newdata = targ, conditioning = "history", history = hist)
  pp <- predict(m, newdata = targ)
  with_hist <- targ$id <= 40
  expect_equal(ph[!with_hist], pp[!with_hist])
  expect_false(isTRUE(all.equal(ph[with_hist], pp[with_hist])))
  # The target outcomes are not used
  targ2 <- targ
  targ2$y <- 0
  expect_equal(ph, predict(m, newdata = targ2, conditioning = "history",
                           history = hist))
  # Hand calculation for one panel
  cf <- coef(m)
  p1 <- stats::plogis(cf[["Mixing:Component1:(Intercept)"]])
  l <- exp(cf[1:2])
  yh <- hist$y[hist$id == 1]
  a <- c(p1 * prod(dpois(yh, l[1])), (1 - p1) * prod(dpois(yh, l[2])))
  expect_equal(ph[targ$id == 1][1], sum(a / sum(a) * l))
  expect_error(predict(m, newdata = targ, conditioning = "history"),
               "needs `history`")
  # No history conditioning without panels
  m0 <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 1, seed = 1)
  expect_error(predict(m0, newdata = targ, conditioning = "history",
                       history = hist), "panel model")
})

test_that("fitted values, residuals and probabilities are consistent", {
  d <- fm_data(300)
  m <- countreg_mix(y ~ x, data = d, components = 2, mixing_formula = ~ z,
                    n_starts = 2, seed = 1)
  mu <- fitted(m)
  expect_equal(residuals(m), d$y - mu)
  v <- predict(m, newdata = d, type = "variance")
  expect_equal(residuals(m, type = "pearson"), (d$y - mu) / sqrt(v))
  set.seed(1)
  q <- residuals(m, type = "quantile")
  expect_length(q, nrow(d))
  expect_true(all(is.finite(q)))
  # Quantile residuals lie between the mixture CDF limits
  P <- predict(m, newdata = d, type = "prob", max_count = max(d$y))
  cdf <- t(apply(P, 1, cumsum))
  up <- cdf[cbind(seq_len(nrow(d)), d$y + 1)]
  lo <- ifelse(d$y == 0, 0, cdf[cbind(seq_len(nrow(d)), pmax(d$y, 1))])
  u <- stats::pnorm(q)
  expect_true(all(u >= lo - 1e-8 & u <= up + 1e-8))
  # Rootogram uses the mixture probabilities of every row
  rt <- rootogram(m, plot = FALSE, overflow = TRUE)
  expect_equal(sum(rt$expected), nrow(d), tolerance = 1e-8)
  expect_equal(sum(rt$observed), nrow(d))
  # augment() and glance()
  au <- augment(m)
  expect_equal(au$.fitted, unname(mu))
  g <- glance(m)
  expect_equal(nrow(g), 1)
  expect_equal(g$components, 2L)
  expect_equal(g$df, length(coef(m)))
  td <- tidy(m, conf.int = TRUE, exponentiate = TRUE)
  expect_identical(td$component[td$term == "Mixing:Component1:z"], "mixing")
  expect_equal(td$estimate[td$term == "Mixing:Component1:z"],
               unname(coef(m)["Mixing:Component1:z"]))
  expect_equal(td$estimate[td$term == "Component1:Mean:x"],
               exp(unname(coef(m)["Component1:Mean:x"])))
})

test_that("Pearson residuals are refused when the variance is infinite", {
  d <- fm_data(300)
  m <- suppressWarnings(suppressMessages(countreg_mix(
    y ~ 1, data = d, family = "GW", components = 1)))
  # Force rho <= 2 (infinite variance)
  m$model$estimate[["Component1:ln(rho-1)"]] <- log(0.5)
  expect_true(all(is.infinite(predict(m, newdata = d[1:3, ],
                                      type = "variance"))))
  expect_error(residuals(m, type = "pearson"), "not available")
})

test_that("simulate() draws one class per unit", {
  set.seed(23)
  G <- 50
  d <- data.frame(id = rep(seq_len(G), each = 4), y = rpois(4 * G, 2))
  m <- countreg_mix(y ~ 1, data = d, components = 2, group_var = "id",
                    n_starts = 1, seed = 1)
  s <- simulate(m, nsim = 3, seed = 1)
  expect_equal(dim(s), c(4 * G, 3))
  expect_equal(dim(attr(s, "classes")), c(G, 3))
  expect_equal(s, simulate(m, nsim = 3, seed = 1), ignore_attr = TRUE)
  # Simulated mixture moments agree with the analytic moments
  m$model$estimate[] <- c(log(1), log(5), log(0.6 / 0.4))
  m$model$prior[] <- rep(c(0.6, 0.4), each = G)
  sims <- unlist(simulate(m, nsim = 400, seed = 2))
  expect_equal(mean(sims), 2.6, tolerance = 0.05)
  # Within-panel dependence from the shared class: Cov = 0.24 * 16 = 3.84
  ss <- simulate(m, nsim = 2000, seed = 3)
  r1 <- unlist(ss[seq(1, 4 * G, by = 4), ])
  r2 <- unlist(ss[seq(2, 4 * G, by = 4), ])
  expect_equal(stats::cov(r1, r2), 3.84, tolerance = 0.1)
})

test_that("marginal effects include the membership pathway", {
  d <- fm_data(400)
  m <- countreg_lc(y ~ x, data = d, classes = 2, membership_formula = ~ z,
                   n_starts = 2, seed = 1)
  me <- margEffTable(m)
  expect_setequal(me$term, c("x", "z"))
  # z enters only the membership model; compare the discrete change
  me_z <- margEffTable(m, vars = "z", indicator = "discrete_change")
  d0 <- d
  d0$z <- 0
  d1 <- d
  d1$z <- 1
  dc <- mean(predict(m, newdata = d1) - predict(m, newdata = d0))
  expect_equal(me_z$estimate, round(dc, 3))
  expect_true(is.finite(me_z$std_error))
  expect_error(margEffTable(m, se = "bootstrap"), "delta")
})

test_that("unsupported operations explain themselves", {
  d <- fm_data(200)
  m <- countreg_mix(y ~ x, data = d, components = 2, n_starts = 1, seed = 1)
  expect_error(predict(m, newdata = d, interval = "confidence",
                       type = "variance"), "type = \"response\"")
  expect_error(tolerance_interval(m, d), "countreg")
  expect_error(regCompTest(m), "does not test finite-mixture")
  expect_error(predict(m, method = "Exact"), "conditioning")
  # keep_data = FALSE
  k <- countreg_mix(y ~ x, data = d, components = 2, n_starts = 1, seed = 1,
                    keep_data = FALSE)
  expect_null(k$data)
  expect_equal(predict(k, newdata = d[1:3, ]), predict(m, newdata = d[1:3, ]))
  expect_error(fitted(k), "keep_data")
  expect_error(predict(k), "keep_data")
})

test_that("summary() reports classes, starts and numerical checks", {
  d <- fm_data(300)
  m <- countreg_lc(y ~ x, data = d, classes = 2, membership_formula = ~ z,
                   n_starts = 3, seed = 1)
  s <- summary(m)
  expect_s3_class(s, "summary.countreg_mix")
  expect_equal(nrow(s$classes), 2)
  expect_equal(sum(s$classes$posterior_share), 1)
  expect_equal(sum(s$classes$prior_share), 1)
  expect_equal(s$classes$effective_size,
               unname(colSums(component_probs(m, conditioning =
                                                "posterior"))))
  out <- capture.output(print(s))
  expect_true(any(grepl("Latent-class count model with 2 classes", out)))
  expect_true(any(grepl("No likelihood-ratio test", out)))
  expect_true(any(grepl("Mixing \\(multinomial logit\\)", out)))
  # Entropy definition
  tau <- component_probs(m, conditioning = "posterior")
  expect_equal(m$model$entropy, -sum(tau * log(tau)))
  expect_equal(m$model$ICL, BIC(m) + 2 * m$model$entropy)
  expect_output(print(m), "Latent-class count model")
})
