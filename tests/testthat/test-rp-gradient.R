# Analytic gradient of the countreg_rp() log-likelihood (rp_gradient.R).

# Per-panel scores at the starting values, with the analytic gradient and
# with maxLik's numeric gradient.
rp_scores <- function(...) {
  one <- function(opt) {
    old <- options(flexCountReg.rp_gradient = opt)
    on.exit(options(old))
    m <- suppressWarnings(suppressMessages(countreg_rp(
      ..., ndraws = 25, max.iters = 0, method = "BHHH")))
    m$model$gradientObs
  }
  list(analytic = one("analytic"), numeric = one("numeric"))
}

expect_scores_match <- function(s, tol = 1e-6) {
  a <- s$analytic
  n <- s$numeric
  expect_identical(dim(a), dim(n))
  expect_lt(max(abs(a - n)) / max(1, abs(n)), tol)
  expect_lt(max(abs(colSums(a) - colSums(n)) / pmax(1, abs(colSums(n)))),
            tol)
}

rp_data <- function() {
  data("washington_roads", package = "flexCountReg", envir = environment())
  w <- washington_roads[1:300, ]
  set.seed(4)
  w$z <- stats::rnorm(nrow(w))
  w$s <- stats::runif(nrow(w))
  w
}

test_that("the analytic RP score matches the numeric one for every draw type", {
  w <- rp_data()
  f <- Total_crashes ~ lnaadt
  expect_scores_match(rp_scores(
    f, rpar_formula = ~ -1 + speed50 + lnlength, data = w, family = "NB2",
    group_var = "ID", start.vals = c(-6, 0.7, 0.1, 0.8, 0.3, 0.2, -1)))
  expect_scores_match(rp_scores(
    f, rpar_formula = ~ -1 + speed50 + lnlength + z, data = w,
    family = "POISSON", group_var = "ID",
    rpardists = c(speed50 = "ln", lnlength = "u", z = "t"),
    start.vals = c(-6, 0.7, -1, 0.8, 0.1, 0.2, 0.3, 0.2)))
  expect_scores_match(rp_scores(
    f, rpar_formula = ~ -1 + lnlength, data = w, family = "NB1",
    rpardists = c(lnlength = "g"), start.vals = c(-6, 0.7, 0.8, 0.3, -1)))
  # Draws shared by every panel (the default is one set per panel)
  expect_scores_match(rp_scores(
    f, rpar_formula = ~ -1 + speed50 + lnlength + z, data = w,
    family = "NB2", group_var = "ID", draws = "shared",
    rpardists = c(speed50 = "ln", lnlength = "g", z = "t"),
    start.vals = c(-6, 0.7, -1, 0.8, 0.1, 0.2, 0.3, 0.2, -1)))
  expect_scores_match(rp_scores(
    f, rpar_formula = ~ -1 + speed50 + lnlength, data = w, family = "NB2",
    correlated = TRUE, draws = "shared",
    start.vals = c(-6, 0.7, 0.1, 0.8, 0.3, 0.1, 0.2, -1)))
})

test_that("the analytic RP score covers correlation and heterogeneity", {
  w <- rp_data()
  f <- Total_crashes ~ lnaadt
  expect_scores_match(rp_scores(
    f, rpar_formula = ~ -1 + speed50 + lnlength, data = w, family = "NBP",
    correlated = TRUE, group_var = "ID",
    start.vals = c(-6, 0.7, 0.1, 0.8, 0.3, 0.1, 0.2, -1, 0.5)))
  expect_scores_match(rp_scores(
    f, rpar_formula = ~ -1 + speed50 + lnlength, data = w, family = "NB2",
    rpardists = c(speed50 = "ln", lnlength = "n"),
    het_mean_formula = ~ -1 + z, het_var_formula = ~ -1 + s,
    dis_param_formula_1 = ~ z,
    start.vals = c(-6, 0.7, -1, 0.8, 0.3, 0.2, 0.1, -0.1, 0.2, -1, 0.1)))
  expect_scores_match(rp_scores(
    f, rpar_formula = ~ -1 + lnlength, data = w, family = "NB2",
    rpardists = c(lnlength = "g"), het_mean_formula = ~ -1 + z,
    het_var_formula = ~ -1 + s,
    start.vals = c(-6, 0.7, 0.8, 0.3, 0.1, 0.2, -1)))
  expect_scores_match(rp_scores(
    f, rpar_formula = ~ -1 + speed50 + lnlength, data = w, family = "NB2",
    correlated = TRUE, het_mean_formula = ~ -1 + z,
    het_var_formula = ~ -1 + s, group_var = "ID",
    start.vals = c(-6, 0.7, 0.1, 0.8, 0.3, 0.1, 0.2, 0.1, -0.1, 0.2, -1)))
})

test_that("families without a closed-form score use differences in log mu", {
  w <- rp_data()
  f <- Total_crashes ~ lnaadt
  expect_scores_match(rp_scores(
    f, rpar_formula = ~ -1 + lnlength, data = w, family = "PLN",
    start.vals = c(-6, 0.7, 0.8, 0.3, -1)))
  expect_scores_match(rp_scores(
    f, rpar_formula = ~ -1 + lnlength, data = w, family = "GP1",
    start.vals = c(-6, 0.7, 0.8, 0.3, 0.1)))
  expect_scores_match(rp_scores(
    f, rpar_formula = ~ -1 + lnlength, data = w, family = "SI",
    group_var = "ID", start.vals = c(-6, 0.7, 0.8, 0.3, 0, 0.2)))
})

test_that("weighted panels have weighted scores, and zero weights none", {
  w <- rp_data()
  ids <- unique(w$ID)
  set.seed(1)
  wt <- sample(0:3, length(ids), replace = TRUE)
  w$wt <- wt[match(w$ID, ids)]
  s <- rp_scores(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + lnlength,
                 data = w, family = "NB2", group_var = "ID", weights = "wt",
                 start.vals = c(-6, 0.7, 0.8, 0.3, -1))
  expect_scores_match(s)
  zero <- wt[order(ids)] == 0
  expect_true(all(s$analytic[zero, ] == 0))
})

test_that("closed-form family scores match differences of the log-density", {
  y <- c(0, 1, 3, 10, 40)
  mu <- c(0.5, 2, 3, 8, 30)
  h <- 1e-5
  for (fam in c("POISSON", "NB1", "NB2", "NBP")) {
    np <- lengths(get_params(fam)) > 0
    a <- if (np[1]) c(0.3, 0.8, 1.5, 0.05, 2)
    s <- if (np[2]) c(1.2, 1.5, 1.8, 1.1, 1.6)
    dens <- get_probFunc(fam)
    lf <- function(m, a, s) dens(y, m, a, s, log = TRUE)
    d <- function(up, dn) (up - dn) / (2 * h)
    sc <- .fcr_family_score(fam, y, mu, a, s)
    e <- exp(h)
    expect_equal(sc$eta, d(lf(mu * e, a, s), lf(mu / e, a, s)),
                 tolerance = 1e-7)
    if (np[1]) expect_equal(sc$alpha, d(lf(mu, a * e, s), lf(mu, a / e, s)),
                            tolerance = 1e-7)
    if (np[2]) expect_equal(sc$sigma, d(lf(mu, a, s * e), lf(mu, a, s / e)),
                            tolerance = 1e-7)
  }
  # Near the Poisson limit the NBP score keeps its accuracy
  sc <- .fcr_nb_score("NBP", y, mu, 1e-14, 1.5)
  expect_true(all(is.finite(unlist(sc))))
  expect_equal(sc$eta, y - mu, tolerance = 1e-10)
})

test_that("the fallback score is one-sided where a step leaves the support", {
  # Underdispersed GP1 (phi = -0.25, mu = 1) has support y < 4; a larger
  # mean extends it, so at y = 4 only the upward step has probability
  s <- .fcr_family_score("GP1", c(1, 3), c(1, 1), c(0.75, 0.75))
  expect_true(all(is.finite(s$eta)))
  s_pl <- .fcr_family_score("PLN", c(0, 2), c(1, 1), c(0.5, 0.5),
                            dens_args = list(haltons = NULL,
                                             normed_haltons = NULL,
                                             engine = "quadrature"))
  lf <- function(m) dpLnorm(c(0, 2), m, 0.5, log = TRUE)
  expect_equal(s_pl$eta, (lf(exp(1e-5)) - lf(exp(-1e-5))) / 2e-5,
               tolerance = 1e-6)
})

test_that("fits agree with the analytic and the numeric gradient", {
  skip_on_cran()
  data("washington_roads", package = "flexCountReg")
  fit <- function(opt, method) {
    old <- options(flexCountReg.rp_gradient = opt)
    on.exit(options(old))
    suppressWarnings(countreg_rp(
      Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
      data = washington_roads, family = "NB2", group_var = "ID",
      ndraws = 50, method = method))
  }
  for (method in c("BHHH", "NR", "BFGS")) {
    a <- fit("analytic", method)
    n <- fit("numeric", method)
    expect_equal(as.numeric(logLik(a)), as.numeric(logLik(n)),
                 tolerance = 1e-6)
    # The likelihood depends on |sigma|, so a standard deviation's sign is
    # not identified
    sd_abs <- function(b) {
      sd <- grepl("St.Dev", names(b))
      b[sd] <- abs(b[sd])
      b
    }
    expect_equal(sd_abs(coef(a)), sd_abs(coef(n)), tolerance = 1e-3)
  }
})

test_that("the Newton decrement separates a maximum from a stalled optimizer", {
  H <- -diag(c(4, 1))
  expect_equal(.fcr_newton_decrement(list(gradient = c(2, 1), hessian = H)),
               2)
  expect_equal(.fcr_newton_decrement(list(gradient = c(0, 0), hessian = H)),
               0)
  # A flat direction (singular Hessian) uses the generalized inverse
  expect_equal(.fcr_newton_decrement(
    list(gradient = c(2, 0), hessian = -diag(c(4, 0)))), 1)
  expect_true(is.na(.fcr_newton_decrement(
    list(gradient = c(NA, 1), hessian = H))))
  data("washington_roads", package = "flexCountReg")
  m <- suppressWarnings(countreg_rp(
    Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
    data = washington_roads, family = "NB2", group_var = "ID", ndraws = 50,
    method = "BHHH"))
  expect_true(m$model$converged)
  expect_lt(m$model$newton_decrement, .fcr_decrement_tol)
})
