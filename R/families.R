# Family registry for countreg() and countreg.rp().
#
# Each count family is defined once, here. An entry has:
#   params    list(<parameter 1 label>, <parameter 2 label>); NULL when the
#             family does not have that parameter. Labels are the
#             coefficient names (e.g. "ln(alpha)"); each parameter is
#             exp(W delta) of its linear predictor.
#   start     starting value(s) on the linear-predictor scale for the
#             distribution parameters (recycled).
#   dens      function(y, predicted, alpha, sigma, haltons, normed_haltons,
#             engine) returning P(Y = y), with predicted = exp(X beta +
#             offset) (times the reporting probability) and alpha, sigma the
#             exponentiated distribution parameters.
#   variance  function(mu, alpha, sigma): Var(Y) at the mean parameter mu,
#             or NULL when no closed form is available.
#   rand      function(n, mu, alpha, sigma): random draws.
#   mean_adj  function(alpha, sigma): E[Y] / mu (1 except PLN and PLL).
#   mean      optional function(mu, alpha, sigma): E[Y], for families whose
#             mean is not proportional to mu (default mu * mean_adj).
#
# Adding a family means adding one entry below.

.fcr_family <- function(params, dens, variance = NULL, rand = NULL,
                        mean_adj = function(alpha, sigma) 1, start = 0,
                        mean = NULL) {
  list(params = params, dens = dens, variance = variance, rand = rand,
       mean_adj = mean_adj, start = start, mean = mean)
}

# Density values for the likelihood: parameter values outside the valid
# range (tried by the optimizer) have probability 0, without warnings.
.fcr_quiet_dens <- function(expr) {
  p <- suppressWarnings(expr)
  p[is.na(p)] <- 0
  p
}

# Exact mean and variance of GP1(mu, phi) as implemented. For phi >= 0 these
# are the textbook moments mu and mu (1 + phi)^2. For phi < 0 the support is
# truncated and the PMF renormalized, so mu is a location parameter and the
# moments are summed from the PMF (each distinct (mu, phi) pair once).
.gp1_moments <- function(mu, phi) {
  r <- .fcr_recycle(mu, phi)
  mu <- r[[1]]; phi <- r[[2]]
  m <- mu
  v <- mu * (1 + phi)^2
  neg <- which(!is.na(mu) & !is.na(phi) & phi < 0 & phi > -1 & mu > 0)
  for (idx in split(neg, ceiling(seq_along(neg) / 5000))) {
    lp <- .gp1_trunc_lraw(mu[idx], phi[idx])
    nk <- nrow(lp)
    w <- exp(lp - rep(.fcr_col_max(lp), each = nk))
    w <- w / rep(colSums(w), each = nk)
    y <- 0:(nk - 1)
    m1 <- colSums(y * w)
    m[idx] <- m1
    v[idx] <- colSums(y^2 * w) - m1^2
  }
  list(mean = m, var = v)
}

.fcr_family_registry <- list(
  POISSON = .fcr_family(
    params = list(NULL, NULL),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) stats::dpois(y, predicted),
    variance = function(mu, alpha, sigma) mu,
    rand = function(n, mu, alpha, sigma) stats::rpois(n, mu)),

  GP1 = .fcr_family(
    params = list("ln(phi+1)", NULL),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) .fcr_quiet_dens(dgp1(y, mu = predicted,
                                                  phi = alpha - 1)),
    # phi = alpha - 1; for phi < 0, mu is a location parameter (see
    # .gp1_moments())
    variance = function(mu, alpha, sigma) .gp1_moments(mu, alpha - 1)$var,
    mean = function(mu, alpha, sigma) .gp1_moments(mu, alpha - 1)$mean,
    rand = function(n, mu, alpha, sigma) rgp1(n, mu = mu, phi = alpha - 1)),

  GP2 = .fcr_family(
    params = list("ln(alpha+1)", NULL),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) .fcr_quiet_dens(dgp2(y, mu = predicted,
                                                  alpha = alpha - 1)),
    # GP2(mu, a) is GP1(mu, phi = a * mu), with a = alpha - 1
    variance = function(mu, alpha, sigma)
      .gp1_moments(mu, (alpha - 1) * mu)$var,
    mean = function(mu, alpha, sigma) .gp1_moments(mu, (alpha - 1) * mu)$mean,
    rand = function(n, mu, alpha, sigma) rgp2(n, mu = mu, alpha = alpha - 1)),

  NB1 = .fcr_family(
    params = list("ln(alpha)", NULL),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) .fcr_quiet_dens(
                      stats::dnbinom(y, size = predicted / alpha,
                                     mu = predicted)),
    variance = function(mu, alpha, sigma) mu * (1 + alpha),
    rand = function(n, mu, alpha, sigma)
      stats::rnbinom(n, size = mu / alpha, mu = mu)),

  NB2 = .fcr_family(
    params = list("ln(alpha)", NULL),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) .fcr_quiet_dens(
                      stats::dnbinom(y, size = 1 / alpha, mu = predicted)),
    variance = function(mu, alpha, sigma) mu + alpha * mu^2,
    rand = function(n, mu, alpha, sigma)
      stats::rnbinom(n, size = 1 / alpha, mu = mu)),

  NBP = .fcr_family(
    params = list("ln(alpha)", "ln(p)"),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) .fcr_quiet_dens(stats::dnbinom(
                      y, size = (predicted^(2 - sigma)) / alpha,
                      mu = predicted)),
    variance = function(mu, alpha, sigma) mu + alpha * mu^sigma,
    rand = function(n, mu, alpha, sigma)
      stats::rnbinom(n, size = mu^(2 - sigma) / alpha, mu = mu)),

  PLN = .fcr_family(
    params = list("ln(sigma)", NULL),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) dpLnorm(x = y, mean = predicted, sigma = alpha,
                                    hdraws = haltons, engine = engine),
    variance = function(mu, alpha, sigma) {
      m <- mu * exp(alpha^2 / 2)
      m + (exp(alpha^2) - 1) * m^2
    },
    rand = function(n, mu, alpha, sigma) rpLnorm(n, mean = mu, sigma = alpha),
    mean_adj = function(alpha, sigma) exp(alpha^2 / 2)),

  # PGE and PW: with the mean parameterization the scale of the mixing
  # distribution cancels out of the PMF, so only the shape is estimated and
  # the scale is fixed at 1.
  PGE = .fcr_family(
    params = list("ln(shape)", NULL),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) dpge(y, mean = predicted, shape = alpha,
                                 scale = 1, haltons = haltons),
    variance = function(mu, alpha, sigma) {
      mu + mu^2 * (trigamma(1) - trigamma(alpha + 1)) /
        (digamma(alpha + 1) - digamma(1))^2
    },
    rand = function(n, mu, alpha, sigma)
      rpge(n, mean = mu, shape = alpha, scale = 1)),

  PIG1 = .fcr_family(
    params = list("ln(eta)", NULL),
    dens = function(y, predicted, alpha, sigma, ...)
      dpinvgaus(y, mu = predicted, eta = alpha),
    variance = function(mu, alpha, sigma) mu * (1 + alpha),
    rand = function(n, mu, alpha, sigma) rpinvgaus(n, mu = mu, eta = alpha)),

  PIG2 = .fcr_family(
    params = list("ln(eta)", NULL),
    dens = function(y, predicted, alpha, sigma, ...)
      dpinvgaus(y, mu = predicted, eta = alpha, form = "Type 2"),
    variance = function(mu, alpha, sigma) mu + alpha * mu^2,
    rand = function(n, mu, alpha, sigma)
      rpinvgaus(n, mu = mu, eta = alpha, form = "Type 2")),

  PIG = .fcr_family(
    params = list("ln(eta)", NULL),
    dens = function(y, predicted, alpha, sigma, ...)
      dpinvgamma(y, mu = predicted, eta = alpha),
    variance = function(mu, alpha, sigma) mu + alpha * mu^2,
    rand = function(n, mu, alpha, sigma) rpinvgamma(n, mu = mu, eta = alpha)),

  PL = .fcr_family(
    params = list("ln(theta)", NULL),
    dens = function(y, predicted, alpha, sigma, ...)
      dplind(y, mean = predicted, theta = alpha),
    variance = function(mu, alpha, sigma)
      mu + mu^2 * (1 - 2 / (alpha + 2)^2),
    rand = function(n, mu, alpha, sigma) rplind(n, mean = mu, theta = alpha)),

  PLG = .fcr_family(
    params = list("ln(theta)", "ln(alpha)"),
    dens = function(y, predicted, alpha, sigma, ...)
      dplindGamma(x = y, mean = predicted, theta = alpha, alpha = sigma),
    variance = function(mu, alpha, sigma)
      mu + ((1 + sigma) * (2 - 2 / (alpha + 2)^2) - 1) * mu^2,
    rand = function(n, mu, alpha, sigma)
      rplindGamma(n, mean = mu, theta = alpha, alpha = sigma)),

  PLL = .fcr_family(
    params = list("ln(theta)", "ln(sigma)"),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) dplindLnorm(x = y, mean = predicted, theta = alpha,
                                        sigma = sigma, hdraws = haltons),
    variance = function(mu, alpha, sigma) {
      m <- mu * exp(sigma^2 / 2)
      m + ((2 - 2 / (alpha + 2)^2) * exp(sigma^2) - 1) * m^2
    },
    rand = function(n, mu, alpha, sigma)
      rplindLnorm(n, mean = mu, theta = alpha, sigma = sigma),
    mean_adj = function(alpha, sigma) exp(sigma^2 / 2)),

  PW = .fcr_family(
    params = list("ln(alpha)", NULL),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) dpWeib_cpp(y, mean = predicted, alpha = alpha,
                                       sigma = 1, h = haltons),
    variance = function(mu, alpha, sigma)
      mu + (gamma(1 + 2 / alpha) / gamma(1 + 1 / alpha)^2 - 1) * mu^2,
    rand = function(n, mu, alpha, sigma)
      rpoisweibull(n, lambda = mu / gamma(1 + 1 / alpha),
                   alpha = alpha, sigma = 1)),

  SI = .fcr_family(
    params = list("gamma", "ln(sigma)"),
    # gamma can be negative, so it is passed as log(alpha) = W delta
    dens = function(y, predicted, alpha, sigma, ...)
      dsichel(x = y, mu = predicted, sigma = sigma, gamma = log(alpha)),
    variance = function(mu, alpha, sigma) {
      gam <- log(alpha)
      cc <- besselK(1 / sigma, gam + 1) / besselK(1 / sigma, gam)
      mu + mu^2 * (2 * sigma * (gam + 1) / cc + 1 / cc^2 - 1)
    },
    rand = function(n, mu, alpha, sigma)
      rsichel(n, mu = mu, sigma = sigma, gamma = log(alpha)),
    start = 0.1),

  GW = .fcr_family(
    params = list("ln(k)", "ln(rho-1)"),
    dens = function(y, predicted, alpha, sigma, ..., log = FALSE)
      dgwar(y, mu = predicted, k = alpha, rho = 1 + sigma, log = log),
    variance = function(mu, alpha, sigma) {
      k <- alpha; rho <- 1 + sigma; a <- mu * (rho - 1) / k
      ifelse(rho > 2, a * k * (a + rho - 1) * (k + rho - 1) /
               ((rho - 1)^2 * (rho - 2)), Inf)
    },
    rand = function(n, mu, alpha, sigma)
      rgwar(n, mu = mu, k = alpha, rho = 1 + sigma),
    start = 0.1),

  COM = .fcr_family(
    params = list("ln(nu)", NULL),
    dens = function(y, predicted, alpha, sigma, ...)
      dcom(x = y, mu = predicted, nu = alpha),
    variance = function(mu, alpha, sigma) {
      r <- .fcr_recycle(mu, alpha)
      vapply(seq_along(r[[1]]), function(i) {
        yy <- 0:max(100, ceiling(20 * r[[1]][i] + 50))
        p <- dcom(yy, mu = r[[1]][i], nu = r[[2]][i])
        p <- p / sum(p)
        sum(yy^2 * p) - sum(yy * p)^2
      }, numeric(1))
    },
    rand = function(n, mu, alpha, sigma) rcom(n, mu = mu, nu = alpha))
)

# ---- Accessors used throughout the package ----------------------------------

.fcr_families <- names(.fcr_family_registry)

.fcr_get_family <- function(family) {
  f <- .fcr_family_registry[[family]]
  if (is.null(f)) {
    stop("Unknown `family` \"", family, "\". Must be one of: ",
         toString(.fcr_families), ".", call. = FALSE)
  }
  f
}

# Parameter labels: list(<parameter 1>, <parameter 2>)
# Errors for distribution-parameter formulas a family cannot use.
.fcr_check_dis_formulas <- function(family, params, f1, f2) {
  if (is.null(params[[1]]) && !is.null(f1)) {
    stop("`dis_param_formula_1` was given, but family \"", family,
         "\" has no distribution parameter to model.", call. = FALSE)
  }
  if (is.null(params[[2]]) && !is.null(f2)) {
    why <- if (family %in% c("PGE", "PW")) {
      paste0(" With the mean parameterization the scale of its mixing ",
             "distribution cancels out of the likelihood, so it is not ",
             "estimated (it was estimated, without effect, in ",
             "flexCountReg <= 0.1.2).")
    } else ""
    stop("`dis_param_formula_2` was given, but family \"", family,
         "\" has no second distribution parameter to model.", why,
         call. = FALSE)
  }
  invisible(NULL)
}

get_params <- function(family) {
  .fcr_family_registry[[family]]$params
}

# Probability function P(Y = y | ...)
get_probFunc <- function(family) {
  .fcr_family_registry[[family]]$dens
}

# Variance at mean parameter mu (NULL when not available)
.fcr_family_variance <- function(family, mu, alpha = NULL, sigma = NULL) {
  v <- .fcr_family_registry[[family]]$variance
  if (is.null(v)) NULL else v(mu, alpha, sigma)
}

# Random draws given the mean parameter and distribution parameters
.fcr_family_r <- function(family, n, mu, alpha = NULL, sigma = NULL) {
  r <- .fcr_get_family(family)$rand
  if (is.null(r)) {
    stop("simulate() is not available for family \"", family, "\".",
         call. = FALSE)
  }
  r(n, mu, alpha, sigma)
}

# ---- Underreporting (binomial thinning) -------------------------------------
#
# Observed counts are Y_obs | Y_true ~ Binomial(Y_true, p). For Poisson, NB2,
# PIG2, PIG and the mixed-Poisson families whose mixing distribution does not
# depend on the mean (PLN, PL, PLG, PLL, PGE, PW, SI), Y_obs has the same
# family with mean p * mu and the same distribution parameters. For NB1 and
# PIG1 (dispersion proportional to the mean) the observed dispersion is
# p * alpha, and for NBP it is alpha * p^(2 - P). GP1, GP2, COM and GW are not
# closed under thinning, so underreporting models are not available for them.
.fcr_thinning_families <- c("POISSON", "NB1", "NB2", "NBP", "PLN", "PGE",
                            "PIG1", "PIG2", "PIG", "PL", "PLG", "PLL", "PW",
                            "SI")

.fcr_check_thinning <- function(family) {
  if (!family %in% .fcr_thinning_families) {
    stop("Underreporting models are not available for family \"", family,
         "\": the distribution of the reported counts (binomial thinning of ",
         "the true counts) is not of the same family. Available families: ",
         toString(.fcr_thinning_families), ".", call. = FALSE)
  }
  invisible(NULL)
}

# Distribution parameters of the reported counts, given those of the true
# counts and the reporting probability p.
.fcr_thin_params <- function(family, alpha, sigma, p) {
  if (is.null(p) || all(p == 1)) return(list(alpha = alpha, sigma = sigma))
  if (family %in% c("NB1", "PIG1")) alpha <- alpha * p
  if (family == "NBP") alpha <- alpha * p^(2 - sigma)
  list(alpha = alpha, sigma = sigma)
}

# Mean of the family at mean parameter mu: E[Y] (mu times the mean
# adjustment, or the family's own mean function)
.fcr_family_mean <- function(family, mu, alpha = NULL, sigma = NULL) {
  f <- .fcr_family_registry[[family]]
  if (!is.null(f$mean)) return(f$mean(mu, alpha, sigma))
  mu * .fcr_mean_adjust(family, alpha, sigma)
}

# Is E[Y] proportional to mu (so it can be applied after averaging over
# random parameters)? Not for GP1/GP2 with a negative dispersion.
.fcr_mean_is_linear <- function(family, alpha = NULL) {
  if (!family %in% c("GP1", "GP2") || is.null(alpha)) return(TRUE)
  all(alpha >= 1, na.rm = TRUE)   # alpha = dispersion + 1
}

# Ratio of the mean to exp(X beta)
.fcr_mean_adjust <- function(family, vec_1, vec_2) {
  f <- .fcr_family_registry[[family]]
  if (is.null(f)) return(1)
  a <- f$mean_adj(vec_1, vec_2)
  if (length(a) == 0) 1 else a
}
