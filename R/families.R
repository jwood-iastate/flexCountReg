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
#
# Adding a family means adding one entry below.

.fcr_family <- function(params, dens, variance = NULL, rand = NULL,
                        mean_adj = function(alpha, sigma) 1, start = 0) {
  list(params = params, dens = dens, variance = variance, rand = rand,
       mean_adj = mean_adj, start = start)
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
                    engine) dgp1(y, mu = predicted, phi = (alpha - 1)),
    variance = function(mu, alpha, sigma) alpha^2 * mu,   # phi = alpha - 1
    rand = function(n, mu, alpha, sigma) rgp1(n, mu = mu, phi = alpha - 1)),

  GP2 = .fcr_family(
    params = list("ln(alpha+1)", NULL),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) dgp2(y, mu = predicted, alpha = (alpha - 1)),
    variance = function(mu, alpha, sigma) mu * (1 + (alpha - 1) * mu)^2,
    rand = function(n, mu, alpha, sigma) rgp2(n, mu = mu, alpha = alpha - 1)),

  NB1 = .fcr_family(
    params = list("ln(alpha)", NULL),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) stats::dnbinom(y, size = predicted / alpha,
                                           mu = predicted),
    variance = function(mu, alpha, sigma) mu * (1 + alpha),
    rand = function(n, mu, alpha, sigma)
      stats::rnbinom(n, size = mu / alpha, mu = mu)),

  NB2 = .fcr_family(
    params = list("ln(alpha)", NULL),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) stats::dnbinom(y, size = 1 / alpha,
                                           mu = predicted),
    variance = function(mu, alpha, sigma) mu + alpha * mu^2,
    rand = function(n, mu, alpha, sigma)
      stats::rnbinom(n, size = 1 / alpha, mu = mu)),

  NBP = .fcr_family(
    params = list("ln(alpha)", "ln(p)"),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) stats::dnbinom(
                      y, size = (predicted^(2 - sigma)) / alpha,
                      mu = predicted),
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

  PGE = .fcr_family(
    params = list("ln(shape)", "ln(scale)"),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) dpge(y, mean = predicted, shape = alpha,
                                 scale = sigma, haltons = haltons),
    variance = function(mu, alpha, sigma) {
      mu + mu^2 * (trigamma(1) - trigamma(alpha + 1)) /
        (digamma(alpha + 1) - digamma(1))^2
    },
    rand = function(n, mu, alpha, sigma)
      rpge(n, mean = mu, shape = alpha, scale = sigma)),

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
    params = list("ln(alpha)", "ln(sigma)"),
    dens = function(y, predicted, alpha, sigma, haltons, normed_haltons,
                    engine) dpWeib_cpp(y, mean = predicted, alpha = alpha,
                                       sigma = sigma, h = haltons),
    variance = function(mu, alpha, sigma)
      mu + (gamma(1 + 2 / alpha) / gamma(1 + 1 / alpha)^2 - 1) * mu^2,
    rand = function(n, mu, alpha, sigma)
      rpoisweibull(n, lambda = mu / (sigma * gamma(1 + 1 / alpha)),
                   alpha = alpha, sigma = sigma)),

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

# Ratio of the mean to exp(X beta)
.fcr_mean_adjust <- function(family, vec_1, vec_2) {
  f <- .fcr_family_registry[[family]]
  if (is.null(f)) return(1)
  a <- f$mean_adj(vec_1, vec_2)
  if (length(a) == 0) 1 else a
}
