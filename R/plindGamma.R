# Parameter convention:
#   alpha = 1/r; G ~ Gamma(shape = 1/alpha, rate = 1/alpha).
#   E(Y) = mean; Var(Y) = mean + [(1+alpha)*(2-2/(theta+2)^2)-1]*mean^2.
#
# Computational representation:
#   J = 1 with probability theta/(1+theta), otherwise J = 2;
#   X | J ~ Gamma(J, rate=1); b = (theta+2)/(theta+1); W = X/b;
#   Y | W ~ NB(size=1/alpha, mu=mean*W).
# We integrate each gamma component separately after t=log(X), scaling
# the integrand at its mode. The CDF uses pnbinom(), not a sum of PMFs.
# No gsl dependency is required by these replacements.
#
# Extensions: mean=0 is degenerate at zero; alpha=0 is the exact
# Poisson--Lindley limit. Positive alpha is never silently replaced by zero.

.plg_logadd <- function(a, b) {
  if (is.na(a) || is.na(b)) return(NA_real_)
  if (a == -Inf) return(b)
  if (b == -Inf) return(a)
  m <- max(a, b)
  m + log1p(exp(min(a, b) - m))
}

.plg_log1mexp <- function(a) {
  if (is.na(a) || a > 0) stop("Invalid log probability in complement")
  if (a == -Inf) return(0)
  if (a == 0) return(-Inf)
  if (a < -log(2)) log1p(-exp(a)) else log(-expm1(a))
}

.plg_control <- function(rel.tol, subdivisions) {
  if (length(rel.tol) != 1L || !is.finite(rel.tol) ||
      rel.tol < 1e-12 || rel.tol > 1e-3)
    stop("'rel.tol' must be a finite scalar between 1e-12 and 1e-3")
  if (length(subdivisions) != 1L || !is.finite(subdivisions) ||
      subdivisions < 50 || subdivisions > .Machine$integer.max ||
      subdivisions != floor(subdivisions))
    stop("'subdivisions' must be an integer of at least 50")
}

.plg_flag <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.logical(x))
    stop("'", name, "' must be TRUE or FALSE")
}

.plg_inputs <- function(x, mean, theta, alpha) {
  a <- list(x = x, mean = mean, theta = theta, alpha = alpha)
  if (!all(vapply(a, is.numeric, logical(1))))
    stop("Count and distribution arguments must be numeric")
  if (any(lengths(a) == 0L)) return(NULL)
  n <- max(lengths(a))
  lapply(a, rep_len, length.out = n)
}

# Conditional log probability at log(mean * W). For the CDF the requested
# tail is passed directly to the NB/Poisson routine.
.plg_conditional <- function(v, logmu, r, kind, lower.tail) {
  ans <- numeric(length(logmu))
  normal <- is.finite(logmu) & logmu >= log(.Machine$double.xmin) &
    logmu <= log(.Machine$double.xmax)
  if (any(normal)) {
    mu <- exp(logmu[normal])
    if (kind == "pmf") {
      ans[normal] <- if (is.infinite(r))
        stats::dpois(v, lambda = mu, log = TRUE) else
          stats::dnbinom(v, size = r, mu = mu, log = TRUE)
    } else {
      ans[normal] <- if (is.infinite(r))
        stats::ppois(v, lambda = mu, lower.tail = lower.tail, log.p = TRUE) else
          stats::pnbinom(v, size = r, mu = mu, lower.tail = lower.tail, log.p = TRUE)
    }
  }
  # Values beyond ordinary floating-point means are not automatically zero.
  # PMFs admit a log-mean formula; CDFs use the NB probability parameter when
  # representable. An unresolved calculation raises a numerical failure.
  other <- which(!normal)
  for (i in other) {
    lm <- logmu[i]
    if (!is.finite(lm)) stop("Nonfinite conditional log mean")
    if (kind == "pmf") {
      if (is.infinite(r)) {
        ans[i] <- if (lm > log(.Machine$double.xmax)) -Inf else
          v * lm - exp(lm) - lgamma(v + 1)
      } else {
        lr <- log(r)
        softplus <- function(z) if (z > 0) z + log1p(exp(-z)) else log1p(exp(z))
        lc <- if (v == 0) 0 else -lbeta(r, v + 1) - log(r + v)
        ans[i] <- lc - r * softplus(lm - lr) - v * softplus(lr - lm)
      }
    } else if (lm < log(.Machine$double.xmin)) {
      # The survival series starts at k=v+1. Bound all successive PMF
      # ratios by rho. If rho <= machine epsilon, its first term gives
      # the tail to floating-point relative precision, even when mu
      # itself underflows. Never round a tiny log survival to -Inf.
      k <- v + 1
      if (is.infinite(r)) {
        logrho <- lm - log(k + 1)
      } else {
        lr <- log(r)
        z <- lr - lm
        sp <- if (z > 0) z + log1p(exp(-z)) else log1p(exp(z))
        logrho <- max(0, .plg_logadd(log(k), lr) - log(k + 1)) - sp
      }
      if (logrho > log(.Machine$double.eps))
        stop("Conditional small-mean tail requires higher precision")
      lt <- .plg_conditional(k, lm, r, "pmf", TRUE)
      ans[i] <- if (lower.tail) .plg_log1mexp(lt) else lt
    } else if (!is.infinite(r)) {
      prob <- stats::plogis(log(r) - lm)
      if (prob <= 0 || prob >= 1)
        stop("Conditional CDF exceeds representable NB probability range")
      ans[i] <- stats::pnbinom(v, size = r, prob = prob,
                               lower.tail = lower.tail, log.p = TRUE)
    } else {
      stop("Conditional Poisson CDF exceeds representable mean range")
    }
  }
  if (anyNA(ans) || any(ans == Inf)) stop("Conditional probability evaluation failed")
  ans
}

# Integrate one log-scaled gamma component on the entire real t axis.
.plg_component <- function(v, mean, theta, alpha, shape, kind,
                           lower.tail, rel.tol, subdivisions) {
  r <- if (alpha == 0) Inf else 1 / alpha
  if (alpha > 0 && !is.finite(r))
    stop("1/alpha overflowed; use alpha=0 only for the exact limit")
  b <- 1 + 1 / (theta + 1) # avoids (theta+2)/(theta+1) overflow
  logscale <- log(mean) - log(b)
  integration_peak <- NULL
  lf <- function(t) {
    out <- rep(-Inf, length(t))
    # At +Inf the exp(-exp(t)) factor vanishes. At -Inf both gamma
    # components vanish on the log-X integration scale.
    active <- is.finite(t) & t <= log(.Machine$double.xmax)
    if (any(active)) {
      u <- t[active]
      prior <- shape * u - exp(u) - lgamma(shape)
      # Far beyond the density's right tail, prior=-Inf needs no NB call.
      keep <- is.finite(prior)
      # Since conditional probabilities are <=1, this bounds the scaled
      # integrand below the smallest normal double. It also avoids asking
      # the conditional routine to resolve irrelevant extreme means.
      if (!is.null(integration_peak))
        keep <- keep & prior >= integration_peak + log(.Machine$double.xmin)
      val <- rep(-Inf, length(u))
      if (any(keep)) val[keep] <- prior[keep] +
        .plg_conditional(v, logscale + u[keep], r, kind, lower.tail)
      out[active] <- val
    }
    if (anyNA(out) || any(out == Inf)) stop("Invalid log integrand")
    out
  }
  
  # Locate an interior mode without assuming it is near X=1. Nonunit means
  # and upper-tail counts can move it substantially.
  found <- FALSE
  for (span in c(8, 16, 32, 64, 128, 256, 512, 1024)) {
    grid <- seq(-span, min(span, 700), by = 2)
    values <- lf(grid)
    j <- which.max(values)
    if (is.finite(values[j]) && j > 1 && j < length(grid) &&
        values[1] < values[j] - 35 && tail(values, 1) < values[j] - 35) {
      found <- TRUE
      break
    }
  }
  if (!found) stop("Could not bracket the integration mode")
  opt <- stats::optimize(lf, c(grid[j-1], grid[j+1]), maximum = TRUE,
                         tol = 1e-9)
  mode <- opt$maximum
  peak <- opt$objective
  if (!is.finite(peak)) stop("Nonfinite integration scale")
  integration_peak <- peak
  
  # Locate a one-log-unit drop on each side. Separate scales prevent a
  # narrow or asymmetric integrand from being missed by infinite-interval
  # quadrature after transformation to u >= 0.
  side_scale <- function(direction) {
    h <- 0.25
    for (j in seq_len(20L)) {
      if (lf(mode + direction*h) <= peak - 1) break
      h <- h * 2
    }
    if (lf(mode + direction*h) > peak - 1) stop("Integration scale not bracketed")
    f <- function(d) {
      z <- lf(mode + direction*d) - peak + 1
      if (z == -Inf) -.Machine$double.xmax else z
    }
    width <- stats::uniroot(f, c(0, h), tol = 1e-10)$root
    if (!is.finite(width) || width <= 0) stop("Integration width unresolved")
    width
  }
  integrate_side <- function(direction) {
    width <- side_scale(direction)
    g <- function(u) {
      z <- lf(mode + direction*width*u) - peak
      ans <- exp(z)
      if (any(!is.finite(ans))) stop("Scaled integrand overflow")
      ans
    }
    ans <- stats::integrate(g, 0, Inf, subdivisions = subdivisions,
                            rel.tol = rel.tol, abs.tol = 0,
                            stop.on.error = TRUE)
    if (!identical(ans$message, "OK") || !is.finite(ans$value) || ans$value <= 0)
      stop("Quadrature failed or returned unresolved mass")
    if (ans$abs.error > 5*rel.tol*ans$value)
      stop("Quadrature error estimate exceeds tolerance")
    log(width) + log(ans$value)
  }
  result <- peak + .plg_logadd(integrate_side(-1), integrate_side(1))
  if (!is.finite(result)) stop("Integrated probability is not resolved")
  result
}

.plg_logprob <- function(v, mean, theta, alpha, kind, lower.tail,
                         rel.tol, subdivisions) {
  a <- .plg_component(v, mean, theta, alpha, 1, kind, lower.tail,
                      rel.tol, subdivisions)
  b <- .plg_component(v, mean, theta, alpha, 2, kind, lower.tail,
                      rel.tol, subdivisions)
  ans <- .plg_logadd(log(theta) - log1p(theta) + a,
                     -log1p(theta) + b)
  if (!is.finite(ans) || ans > 5*rel.tol)
    stop("Integrated probability outside [0,1]")
  # Only remove small positive roundoff within the integration tolerance.
  min(ans, 0)
}

.plg_evaluate <- function(v, mean, theta, alpha, kind, lower.tail,
                          rel.tol, subdivisions) {
  ans <- .plg_logprob(v, mean, theta, alpha, kind, lower.tail,
                      rel.tol, subdivisions)
  # A probability near one is represented more accurately through its small
  # complementary tail. Each small tail is itself integrated directly.
  if (kind == "cdf" && ans > -log(2)) {
    opposite <- .plg_logprob(v, mean, theta, alpha, "cdf", !lower.tail,
                             rel.tol, subdivisions)
    ans <- .plg_log1mexp(opposite)
  }
  ans
}

.plg_vector <- function(x, mean, theta, alpha, kind, lower.tail,
                        give.log, rel.tol, subdivisions) {
  a <- .plg_inputs(x, mean, theta, alpha)
  if (is.null(a)) return(numeric())
  n <- length(a$x)
  out <- rep(NA_real_, n)
  invalid <- noninteger <- FALSE
  failures <- character()
  for (i in seq_len(n)) {
    vals <- c(a$x[i], a$mean[i], a$theta[i], a$alpha[i])
    if (any(is.na(vals) & !is.nan(vals))) next
    if (any(is.nan(vals))) { out[i] <- NaN; next }
    if (!is.finite(a$mean[i]) || a$mean[i] < 0 ||
        !is.finite(a$theta[i]) || a$theta[i] <= 0 ||
        !is.finite(a$alpha[i]) || a$alpha[i] < 0) {
      out[i] <- NaN; invalid <- TRUE; next
    }
    v <- a$x[i]
    if (kind == "pmf") {
      if (!is.finite(v) || v < 0) { out[i] <- -Inf; next }
      if (v != floor(v)) { out[i] <- -Inf; noninteger <- TRUE; next }
      if (a$mean[i] == 0) { out[i] <- if (v == 0) 0 else -Inf; next }
    } else {
      if (v < 0) { out[i] <- if (lower.tail) -Inf else 0; next }
      if (v == Inf || a$mean[i] == 0) {
        out[i] <- if (lower.tail) 0 else -Inf; next
      }
      v <- floor(v) # deliberately no integer coercion
    }
    out[i] <- tryCatch(
      .plg_evaluate(v, a$mean[i], a$theta[i], a$alpha[i], kind, lower.tail,
                    rel.tol, subdivisions),
      error = function(e) {
        failures <<- c(failures, paste0("element ", i, ": ", conditionMessage(e)))
        NaN
      }
    )
  }
  if (invalid) warning("Invalid distribution parameter(s): NaNs produced", call. = FALSE)
  if (noninteger) warning("Noninteger 'x': probability zero", call. = FALSE)
  if (length(failures)) warning(
    "PLG numerical evaluation failed; NaNs produced. ",
    paste(utils::head(failures, 3), collapse = "; "),
    if (length(failures) > 3) paste0("; and ", length(failures)-3, " more"),
    call. = FALSE)
  if (!give.log) out <- exp(out)
  if (length(x) == n && !is.null(names(x))) names(out) <- names(x)
  out
}


#' Poisson-Lindley-Gamma (Negative Binomial-Lindley) Distribution
#'
#' These functions provide density, distribution function, quantile
#' function, and random number generation for the Poisson-Lindley-Gamma
#' (PLG) Distribution
#'
#' The Poisson-Lindley-Gamma is a count distribution that captures high
#' densities for small integer values and provides flexibility for heavier
#' tails.
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mean numeric value or vector of mean values for the distribution
#'   (the values have to be greater than 0).
#' @param theta single value or vector of values for the theta parameter of
#'   the distribution (the values have to be greater than 0).
#' @param alpha single value or vector of values for the `alpha` parameter
#'   of the gamma distribution in the special case that the mean = 1 and
#'   the variance = `alpha` (the values for `alpha` have to be greater
#'   than 0).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]}
#'   otherwise, \eqn{P[X>x]}.
#' @param rel.tol Relative numerical integration tolerance; default 1e-8.
#' @param subdivisions Maximum subintervals per component-side integral.
#'
#' @details
#' \code{dplindGamma} computes the density (PDF) of the
#' Poisson-Lindley-Gamma Distribution.
#'
#' \code{pplindGamma} computes the CDF of the Poisson-Lindley-Gamma
#' Distribution.
#'
#' \code{qplindGamma} computes the quantile function of the
#' Poisson-Lindley-Gamma Distribution.
#'
#' \code{rplindGamma} generates random numbers from the
#' Poisson-Lindley-Gamma Distribution.
#'
#' The compound Probability Mass Function (PMF) for the
#' Poisson-Lindley-Gamma (PLG) distribution is:
#' \deqn{
#' f(x|\mu,\theta,\alpha)=
#' \frac{
#'   (\theta+2)^2\Gamma(x+1/\alpha)
#' }{
#'   \alpha\mu^2(\theta+1)^3\Gamma(1/\alpha)
#' }
#' \left(
#'   \frac{\mu\theta(\theta+1)}{\theta+2}
#'   U\left(
#'     x+1,2-1/\alpha,\frac{1/\alpha(\theta+2)}{\mu(\theta+1)}
#'   \right)
#'   + 1/\alpha(x+1)
#'   U\left(
#'     x+2,3-1/\alpha,\frac{1/\alpha(\theta+2)}{\mu(\theta+1)}
#'   \right)
#' \right)
#' }
#'
#' Where \eqn{\theta} is a distribution parameter from the Poisson-Lindley
#' distribution with the restrictions that \eqn{\theta>0}, \eqn{\alpha} is
#' a parameter for the gamma distribution with the restriction
#' \eqn{\alpha>0}, \eqn{\mu} is the mean value, and \eqn{x} is a
#' non-negative integer, and \deqn{U(a,b,z)} is the Tricomi’s confluent 
#' hypergeometric function to - also known as the confluent
#' hypergeometric function of the second kind
#'
#' The expected value of the distribution is:
#' \deqn{E[x]=\mu}
#'
#' The variance is:
#' \deqn{\sigma^2=\mu+\left(\left(1+\alpha\right)\left(2-\frac{2}
#' {(\theta+2)^2}\right)-1\right)\mu^2}
#'
#' 
#' @returns dplindGamma gives the density, pplindGamma gives the distribution 
#'  function, qplindGamma gives the quantile function, and rplindGamma generates
#'  random  deviates.
#' 
#'  The length of the result is determined by n for rplindGamma, and is the 
#'  maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#' dplindGamma(0, mean=0.75, theta=7, alpha=2)
#' pplindGamma(c(0,1,2,3,5,7,9,10), mean=0.75, theta=3, alpha=0.5)
#' qplindGamma(c(0.1,0.3,0.5,0.9,0.95), mean=1.67, theta=0.5, alpha=0.5)
#' rplindGamma(30, mean=0.5, theta=0.5, alpha=2)
#'
#' @importFrom stats runif
#' @importFrom gsl hyperg_U
#' @useDynLib flexCountReg
#' @name NegativeBinomialLindley
#'
#' @rdname NegativeBinomialLindley
#' @export
dplindGamma <- function(x, mean = 1, theta = 1, alpha = 1, log = FALSE,
                        rel.tol = 1e-8, subdivisions = 200L) {
  .plg_flag(log, "log")
  .plg_control(rel.tol, subdivisions)
  .plg_vector(x, mean, theta, alpha, "pmf", TRUE, log,
              rel.tol, as.integer(subdivisions))
}


##' @rdname NegativeBinomialLindley
#' @export
pplindGamma <- function(q, mean = 1, theta = 1, alpha = 1,
                        lower.tail = TRUE, log.p = FALSE,
                        rel.tol = 1e-8, subdivisions = 200L) {
  .plg_flag(lower.tail, "lower.tail")
  .plg_flag(log.p, "log.p")
  .plg_control(rel.tol, subdivisions)
  .plg_vector(q, mean, theta, alpha, "cdf", lower.tail, log.p,
              rel.tol, as.integer(subdivisions))
}

#' @rdname NegativeBinomialLindley
#' @export
qplindGamma <- Vectorize(function(p, mean=1, theta=1, alpha=1) {
  if(p < 0)
    warning("The value of `p` must be a value greater than 0 and less than 1.")
  if(is.na(p)) warning("The value of `p` cannot be an `NA` value")
  
  if(mean<=0 || theta<=0 || alpha<=0)
    warning(paste(
      "The values of `mean`, `theta`, and `alpha` all have to have",
      "values greater than 0."
    ))
  
  y <- 0
  p_value <- max(
    pplindGamma(y, mean, theta, alpha=alpha),
    .Machine$double.xmin
  )
  while(p_value < p){
    y <- y + 1
    p_value_new <- max(
      pplindGamma(y, mean, theta, alpha=alpha),
      .Machine$double.xmin
    )
    if (!is.na(p_value_new)) p_value <- p_value_new else break
  }
  return(y)
})

#' @rdname NegativeBinomialLindley
#' @export
rplindGamma <- function(n, mean=1, theta=1, alpha=1) {
  
  if(mean<=0 || theta<=0  || alpha<=0)
    warning(paste('The values of `mean`, `theta`, and `alpha` all", 
                  "have to have values greater than 0.'))
  
  u <- runif(n)
  y <- lapply(u, function(p) qplindGamma(p, mean, theta, alpha=alpha))
  return(unlist(y))
}

# 
# Optional local checks. This helper is not exported and does not execute
# automatically when the file is sourced. Use the full precision tests before
# adopting this implementation; they exercise numerical behavior, not just syntax.
# plg_check_stable <- function() {
#   close <- function(a, b, tol = 2e-6) {
#     stopifnot(length(a) == length(b), all(is.finite(a)), all(is.finite(b)),
#               max(abs(a-b)) < tol)
#   }
#   close(dplindGamma(c(0, 3), mean = 1, theta = 1, alpha = .5),
#         c(.5381304934467996, .0547430609355798))
#   stopifnot(is.finite(dplindGamma(0, mean = 1, theta = 1, alpha = .005, log = TRUE)))
#   stopifnot(is.finite(dplindGamma(1000, mean = .1, theta = 1, alpha = .5, log = TRUE)))
#   stopifnot(is.finite(dplindGamma(10000, mean = .1, theta = 1, alpha = .5, log = TRUE)))
#   stopifnot(is.finite(pplindGamma(10000, mean = .1, theta = 1, alpha = .5,
#                                   lower.tail = FALSE, log.p = TRUE)))
#   stopifnot(identical(pplindGamma(c(-Inf, Inf)), c(0, 1)))
#   stopifnot(is.na(pplindGamma(NA_real_)), is.nan(pplindGamma(NaN)))
#   stopifnot(identical(pplindGamma(c(-Inf, Inf), lower.tail = FALSE), c(1, 0)))
#   stopifnot(identical(dplindGamma(c(-1, Inf)), c(0, 0)))
#   stopifnot(suppressWarnings(dplindGamma(.5)) == 0)
#   stopifnot(length(dplindGamma(numeric())) == 0L)
#   stopifnot(identical(dplindGamma(c(0, 1), mean = 0), c(1, 0)))
#   stopifnot(identical(pplindGamma(c(-1, 0, 1), mean = 0), c(0, 1, 1)))
#   q <- c(0, 1, 4, 10)
#   F <- pplindGamma(q, mean = 2, theta = .5, alpha = .2)
#   S <- pplindGamma(q, mean = 2, theta = .5, alpha = .2, lower.tail = FALSE)
#   close(F+S, rep(1, length(q)))
#   stopifnot(all(diff(F) >= 0))
#   close(exp(pplindGamma(q, mean = 2, theta = .5, alpha = .2, log.p = TRUE)), F)
#   close(pplindGamma(q+.9, mean = 2, theta = .5, alpha = .2), F)
#   close(pplindGamma(q, mean = 2, theta = .5, alpha = .2) -
#           pplindGamma(q-1, mean = 2, theta = .5, alpha = .2),
#         dplindGamma(q, mean = 2, theta = .5, alpha = .2))
#   close(dplindGamma(c(0, 1), mean = c(.2, 2), theta = c(.5, 3), alpha = c(.1, 2)),
#         c(dplindGamma(0, .2, .5, .1), dplindGamma(1, 2, 3, 2)))
#   close(sum(dplindGamma(0:60, mean = 1, theta = 1, alpha = .5)) +
#           pplindGamma(60, mean = 1, theta = 1, alpha = .5, lower.tail = FALSE), 1)
#   # Exact alpha=0 Poisson--Lindley boundary has an elementary mixed-Poisson PMF.
#   mu <- .7; theta <- 2; b <- 1+1/(theta+1); rate <- b/mu
#   y <- 0:8
#   reference <- theta/(theta+1)*stats::dnbinom(y, size=1, prob=rate/(rate+1)) +
#     1/(theta+1)*stats::dnbinom(y, size=2, prob=rate/(rate+1))
#   close(dplindGamma(y, mu, theta, alpha=0), reference)
#   message("PLG numerical and boundary checks passed.")
#   invisible(TRUE)
# }
# 
