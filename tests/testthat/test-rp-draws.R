# Per-unit simulation draws in countreg_rp() (review item N4).

test_that("each unit gets its own consecutive block of Halton draws", {
  u <- .fcr_unit_halton(c(2, 1, 2, 3), R = 5, K = 2)
  expect_s3_class(u, "fcr_unit_draws")
  expect_length(u, 2)
  expect_identical(dim(u[[1]]), c(4L, 5L))
  # Rows of the same unit share its draws; different units differ
  expect_identical(u[[1]][1, ], u[[1]][3, ])
  expect_false(isTRUE(all.equal(u[[1]][1, ], u[[1]][2, ])))
  # Unit k uses points 10 + (k - 1) R + 1, ..., 10 + k R of the sequence
  h <- randtoolbox::halton(10 + 3 * 5, 2)
  expect_equal(u[[2]][4, ], h[10 + 11:15, 2])
  expect_equal(u[[1]][2, ], h[10 + 1:5, 1])
  # Rows can be taken without changing the draws
  r <- .fcr_draw_rows(u, c(4, 1))
  expect_identical(r[[2]], u[[2]][c(4, 1), ])
})

test_that("prepared draws match the distribution quantiles", {
  h <- matrix(c(0.1, 0.5, 0.8, 0.3, 0.6, 0.95), 3)
  d <- .fcr_prepare_draws(h, c("t", "g"), n_obs = 2)
  expect_equal(d$Z[[1]][1, ], qtri(c(0.1, 0.5, 0.8), 0, 1))
  expect_equal(d$Z[[1]][2, ], d$Z[[1]][1, ])
  expect_null(d$Z[[2]])
  expect_equal(d$U[[2]][2, ], c(0.3, 0.6, 0.95))
  expect_identical(.fcr_prepare_draws(d, c("t", "g"), 2), d)
})

rp_fit <- function(...) {
  data("washington_roads", package = "flexCountReg", envir = environment())
  suppressWarnings(countreg_rp(
    Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
    data = washington_roads, family = "NB2", group_var = "ID",
    ndraws = 40, ...))
}

test_that("new fits use per-unit draws; shared draws remain available", {
  m <- rp_fit()
  data("washington_roads", package = "flexCountReg")
  expect_identical(m$model$draws, "per_unit")
  expect_identical(m$model$draw_units,
                   levels(factor(washington_roads$ID)))
  s <- rp_fit(draws = "shared")
  expect_identical(s$model$draws, "shared")
  expect_null(s$model$draw_units)
  expect_false(isTRUE(all.equal(coef(m), coef(s))))
  expect_error(rp_fit(draws = "bogus"), "should be one of")
})

test_that("Empirical Bayes regenerates each panel's own draws", {
  m <- rp_fit()
  data("washington_roads", package = "flexCountReg")
  w <- washington_roads
  eb <- predict(m, newdata = w, method = "Individual")
  # Matched by panel id: reordering or subsetting rows changes nothing
  set.seed(2)
  o <- sample(nrow(w))
  expect_equal(predict(m, newdata = w[o, ], method = "Individual"), eb[o],
               ignore_attr = TRUE)
  keep <- w$ID %in% unique(w$ID)[c(5, 9, 40)]
  expect_equal(predict(m, newdata = w[keep, ], method = "Individual"),
               eb[keep], ignore_attr = TRUE)
  # The draws are the estimation draws: a panel's EB coefficient is the
  # posterior mean over the draws used in the likelihood
  st <- local({
    prep <- prepare_model_data(m$formula, w,
                               aux_formulas = list(rpar = m$rpar_formula),
                               group_var = "ID")
    list(group = prep$group)
  })
  u <- .fcr_unit_halton(as.integer(st$group), 40, 1)
  b <- coef(m)
  beta <- b[["speed50:Mean"]] + abs(b[["speed50:St.Dev"]]) *
    stats::qnorm(u[[1]])
  mu <- exp(b[["(Intercept)"]] + b[["lnaadt"]] * w$lnaadt +
              w$speed50 * beta)
  lp <- matrix(stats::dnbinom(w$Total_crashes, size = exp(-b[["ln(alpha)"]]),
                              mu = mu, log = TRUE), nrow(w))
  L <- rowsum(lp, st$group)
  W <- exp(L - apply(L, 1, max))
  W <- W / rowSums(W)
  eb_coef <- attr(eb, "eb_coefficients")[, "speed50"]
  first <- !duplicated(w$ID)
  # Draws of each panel, in the order of the panels' levels (rows of W)
  beta_unit <- beta[match(levels(st$group), st$group), , drop = FALSE]
  by_panel <- rowSums(W * beta_unit)
  expect_equal(eb_coef[first],
               by_panel[match(as.character(w$ID[first]), levels(st$group))],
               tolerance = 1e-10, ignore_attr = TRUE)
  # A panel not in the estimation data gets new draws
  w2 <- w[w$ID == w$ID[1], ]
  w2$ID <- "new-site"
  p_new <- predict(m, newdata = w2, method = "Individual")
  expect_true(all(is.finite(p_new)))
})

test_that("fits without panels use one block of draws per row", {
  data("washington_roads", package = "flexCountReg")
  w <- washington_roads[1:200, ]
  m <- suppressWarnings(countreg_rp(
    Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50, data = w,
    family = "NB2", ndraws = 30))
  expect_length(m$model$draw_units, 200)
  eb <- predict(m, newdata = w, method = "Individual")
  expect_true(all(is.finite(eb)))
})

test_that("scrambled draws follow their seed and keep R's random numbers", {
  a <- .fcr_halton(5, 2, TRUE, seed = 1)
  expect_identical(.fcr_halton(5, 2, TRUE, seed = 1), a)
  expect_false(isTRUE(all.equal(.fcr_halton(5, 2, TRUE, seed = 2), a)))
  # Randomly shifted plain Halton points, one shift per dimension
  set.seed(1)
  sh <- stats::runif(2)
  expect_equal(a, sweep(as.matrix(randtoolbox::halton(5, 2)), 2, sh, `+`) %% 1)
  # "mixed" with no seed: randtoolbox's default scrambling (seed 0), as used
  # by fits with draw_version < 3
  randtoolbox::setSeed(0)
  def <- as.matrix(randtoolbox::halton(5, 2, mixed = TRUE))
  expect_identical(.fcr_halton(5, 2, TRUE, seed = 9), .fcr_halton(5, 2, TRUE, 9))
  expect_identical(.fcr_halton(5, 2, TRUE, scramble = "mixed"), def)
  set.seed(5)
  r1 <- stats::runif(2)
  set.seed(5)
  invisible(.fcr_halton(5, 2, TRUE, seed = 3))
  expect_identical(stats::runif(2), r1)

  set.seed(1)
  m1 <- rp_fit(scrambled = TRUE)
  set.seed(2)
  m2 <- rp_fit(scrambled = TRUE)
  expect_identical(m1$model$draw_version, 4L)
  expect_false(isTRUE(all.equal(coef(m1), coef(m2))))
  data("washington_roads", package = "flexCountReg")
  expect_identical(predict(m1, newdata = washington_roads,
                           method = "Individual"),
                   predict(m1, newdata = washington_roads,
                           method = "Individual"))
})

test_that("scrambled per-unit draws shift each unit's block by its own amount", {
  u <- .fcr_unit_halton(c(2, 1, 2, 3), R = 5, K = 2, scrambled = TRUE,
                        seed = 4)
  p <- .fcr_unit_halton(c(2, 1, 2, 3), R = 5, K = 2)
  sh <- matrix(.fcr_seeded_runif(3 * 2, 4), 3, 2)
  for (k in 1:2) expect_equal(u[[k]], (p[[k]] + sh[c(2, 1, 2, 3), k]) %% 1)
  expect_identical(u[[1]][1, ], u[[1]][3, ])
  # The old "mixed" scrambling is kept for fits with draw_version < 4
  m <- .fcr_unit_halton(c(2, 1, 2, 3), R = 5, K = 2, scrambled = TRUE,
                        seed = 4, scramble = "mixed")
  randtoolbox::setSeed(4)
  h <- randtoolbox::halton(10 + 3 * 5, 2, mixed = TRUE)
  expect_equal(m[[2]][4, ], h[10 + 11:15, 2])
})

test_that("fits with draw_version 3 predict with their mixed draws", {
  set.seed(1)
  m <- rp_fit(scrambled = TRUE)
  data("washington_roads", package = "flexCountReg")
  w <- washington_roads[1:50, ]
  p4 <- predict(m, newdata = w, method = "Individual")
  old <- m
  old$model$draw_version <- 3L
  p3 <- predict(old, newdata = w, method = "Individual")
  expect_false(isTRUE(all.equal(p3, p4)))
  # The two scramblings give different draws from the same seed
  hd <- .fcr_unit_halton(.fcr_rp_draw_units(m$model, w), m$model$ndraws, 1,
                         TRUE, m$model$draw_seed, scramble = "mixed")
  hn <- .fcr_unit_halton(.fcr_rp_draw_units(m$model, w), m$model$ndraws, 1,
                         TRUE, m$model$draw_seed)
  expect_false(isTRUE(all.equal(hd[[1]], hn[[1]])))
})

test_that("scrambled draws converge to the exact random-intercept likelihood", {
  skip_on_cran()
  # Poisson random intercept: the exact maximum log-likelihood, -1069.365,
  # is from adaptive Gauss-Hermite quadrature (countreg_hier(), equal to
  # lme4::glmer(nAGQ = 9)). Scrambled per-unit draws from randtoolbox's
  # "mixed" sequence gave -1069.054 at 2000 draws, a bias that more draws
  # did not remove.
  data("washington_roads", package = "flexCountReg")
  exact <- -1069.365
  fit <- function(R, seed) {
    set.seed(seed)
    suppressWarnings(countreg_rp(
      Total_crashes ~ -1 + lnaadt + speed50 + offset(lnlength),
      rpar_formula = ~ 1, group_var = "ID", family = "POISSON",
      rpardists = c("(Intercept)" = "n"), data = washington_roads,
      ndraws = R, scrambled = TRUE))
  }
  err <- sapply(1:2, function(s) as.numeric(logLik(fit(2000, s))) - exact)
  expect_true(all(abs(err) < 0.06), info = paste(err, collapse = ", "))
  # Fewer draws: larger simulation error
  err100 <- as.numeric(logLik(fit(100, 1))) - exact
  expect_true(abs(err100) < 1.5)

})
