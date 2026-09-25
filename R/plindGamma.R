# Parameter convention:
#   alpha = 1/r; G ~ Gamma(shape = 1/alpha, rate = 1/alpha).
#   E(Y) = mean; Var(Y) = mean + [(1+alpha)*(2-2/(theta+2)^2)-1]*mean^2.
#
# Computational representation:
#   J = 1 with probability theta/(1+theta), otherwise J = 2;
#   X | J ~ Gamma(J, rate=1); b = (theta+2)/(theta+1); W = X/b;
#   Y | W ~ NB(size=1/alpha, mu=mean*W).
# The PMF first tries a guarded GSL Tricomi-U evaluation. On rejection,
# integrate each gamma component after t=log(X), scaling at its mode.
# The CDF always integrates pnbinom()/ppois(), not a sum of PMFs.
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


# Keep the external numerical call separate so fallback tests can substitute
# a failing backend in an isolated environment without changing gsl's namespace.
.plg_gsl_U <- function(a, b, z) {
  gsl::hyperg_U(a, b, z, give = TRUE, strict = TRUE)
}

# Return a log PMF, or NULL to request quadrature. This helper is called
# only for a finite nonnegative integer v and valid positive parameters.
.plg_gsl_logpmf <- function(v, mean, theta, alpha, rel.tol) {
  if (alpha == 0) return(NULL)
  r <- 1 / alpha
  a <- c(v + 1, v + 2)
  b <- c(2 - r, 3 - r)
  logz <- log(r) + log1p(1 / (theta + 1)) - log(mean)
  
  # The historical GSL report identifies a > 8 and z close to zero.
  # 1 is our conservative cutoff, NOT a boundary established by GSL.
  # Large-parameter and representability guards also avoid expensive or
  # poorly resolved recurrences before entering the C library.
  if (!is.finite(r) || r > 100 || any(a > 1e4) ||
      !is.finite(logz) || logz < log(.Machine$double.xmin) ||
      logz > log(.Machine$double.xmax)) return(NULL)
  z <- exp(logz)
  if (!is.finite(z) || z <= 0 || (any(a > 8) && z < 1)) return(NULL)
  
  # give=TRUE exposes val, err and status. Reject warnings as well as
  # errors; a successful-looking numeric value alone is not sufficient.
  u <- tryCatch(
    .plg_gsl_U(a, b, z),
    warning = function(w) NULL,
    error = function(e) NULL
  )
  if (!is.list(u)) return(NULL)
  fields <- c("val", "err", "status")
  valid <- vapply(fields, function(nm) {
    x <- u[[nm]]
    is.numeric(x) && length(x) == 2L && all(is.finite(x))
  }, logical(1))
  if (!all(valid) || any(u$status != 0) ||
      any(u$val < .Machine$double.xmin) || any(u$err < 0)) return(NULL)
  relative.error <- u$err / u$val
  if (any(!is.finite(relative.error)) ||
      any(relative.error > rel.tol / 4)) return(NULL)
  
  # Gamma(v+r)/Gamma(r) = Gamma(v)/B(r,v), for v > 0.
  # This avoids subtracting two nearly equal lgamma values when r >> v.
  lg <- if (v == 0) 0 else lgamma(v) - lbeta(r, v)
  lu <- log(u$val)
  lz <- c(logz, 2 * logz + log1p(v))
  components <- lg + lz + lu
  weights <- c(log(theta) - log1p(theta), -log1p(theta))
  ans <- .plg_logadd(weights[1] + components[1],
                     weights[2] + components[2])
  
  # Conservative roundoff screening for cancellation in the log formula.
  # Neither this heuristic nor GSL's error estimate is a certified bound.
  rounding <- 32 * .Machine$double.eps *
    (1 + abs(lg) + max(abs(lz)) + max(abs(lu)) + max(abs(weights)))
  if (!is.finite(rounding) || rounding > rel.tol / 4 ||
      any(!is.finite(components)) || any(components > 0) ||
      !is.finite(ans) || ans > 0) return(NULL)
  ans
}

.plg_evaluate <- function(v, mean, theta, alpha, kind, lower.tail,
                          rel.tol, subdivisions, method = "auto") {
  if (kind == "pmf" && method == "auto") {
    fast <- .plg_gsl_logpmf(v, mean, theta, alpha, rel.tol)
    if (!is.null(fast)) return(fast)
  }
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
                        give.log, rel.tol, subdivisions, method = "auto") {
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
                    rel.tol, subdivisions, method),
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
#' Probability mass, cumulative distribution, quantiles, and random generation
#' for the mean-parameterized Poisson-Lindley-Gamma (PLG), also called the
#' negative binomial-Lindley (NB-L), distribution. The PMF uses a guarded GSL
#' Tricomi-U evaluation with adaptive quadrature as a fallback.
#'
#' @param x Numeric vector of counts. Negative, infinite, or noninteger counts
#'   have probability zero; noninteger counts produce a warning.
#' @param q Numeric vector of quantiles. Finite values are rounded down.
#' @param p Numeric vector of probabilities in \eqn{[0,1]}.
#' @param n Number of observations to generate: a nonnegative integer scalar.
#' @param mean Marginal mean \eqn{\mu}, finite and nonnegative. May be a vector
#'   for the density, CDF, and quantile functions; must be scalar for random
#'   generation. A zero mean gives a point mass at zero.
#' @param theta Lindley parameter, finite and strictly positive. May be a
#'   vector except for random generation.
#' @param alpha Gamma dispersion, finite and nonnegative. For positive
#'   \code{alpha}, the mean-one gamma factor has shape and rate
#'   \eqn{1/\alpha}, and variance \eqn{\alpha}. Zero selects the exact
#'   Poisson-Lindley limit. May be a vector except for random generation.
#' @param log Logical; return log probabilities from \code{dplindGamma}.
#' @param log.p Logical; return log probabilities from \code{pplindGamma}.
#' @param lower.tail Logical; if \code{TRUE}, return \eqn{P(Y\le q)};
#'   otherwise return \eqn{P(Y>q)}.
#' @param rel.tol Relative tolerance, a scalar between \code{1e-12} and
#'   \code{1e-3}; default \code{1e-8}. Controls adaptive quadrature and, for
#'   the PMF, screening of GSL error estimates and log-formula roundoff.
#'   It is a numerical target, not a certified bound on the final error.
#' @param subdivisions Maximum number of subintervals for each side of each
#'   component integral; an integer of at least 50, default 200.
#' @param method PMF evaluation method. \code{"auto"} (default) attempts GSL
#'   outside the guarded regions and uses quadrature if any check fails.
#'   \code{"quadrature"} always uses adaptive quadrature and is useful for
#'   validation. There is deliberately no option to force an unchecked GSL
#'   value. This argument applies only to \code{dplindGamma}.
#'
#' @details
#' \strong{Parameterization and probability mass.}
#' Let \eqn{L} have Lindley density
#' \deqn{f_L(l)=\frac{\theta^2}{1+\theta}(1+l)e^{-\theta l},\quad l>0,}
#' and let \eqn{W=L/E(L)}, where
#' \eqn{E(L)=(\theta+2)/[\theta(\theta+1)]}. Conditional on \eqn{W},
#' \eqn{Y} is negative binomial with size \eqn{1/\alpha} and mean
#' \eqn{\mu W}. Equivalently, \eqn{Y\mid W,G} is Poisson with mean
#' \eqn{\mu WG}, with independent
#' \eqn{G\sim\mathrm{Gamma}(1/\alpha,1/\alpha)}, using the shape-rate
#' parameterization.
#' Thus \code{mean} is the marginal mean, not the unnormalized Lindley scale.
#'
#' The compound probability mass function (PMF) for the
#' Poisson-Lindley-Gamma distribution is
#' \deqn{
#' f(x\mid\mu,\theta,\alpha)=
#' \frac{
#'   (\theta+2)^2\Gamma(x+1/\alpha)
#' }{
#'   \alpha\mu^2(\theta+1)^3\Gamma(1/\alpha)
#' }
#' \left(
#'   \frac{\mu\theta(\theta+1)}{\theta+2}
#'   U\left(
#'     x+1,2-1/\alpha,\frac{\theta+2}{\alpha\mu(\theta+1)}
#'   \right)
#'   +\frac{x+1}{\alpha}
#'   U\left(
#'     x+2,3-1/\alpha,\frac{\theta+2}{\alpha\mu(\theta+1)}
#'   \right)
#' \right).
#' }
#' This expression applies for \eqn{\mu>0}, \eqn{\theta>0},
#' \eqn{\alpha>0}, and nonnegative integer \eqn{x}. Here \eqn{U} is
#' Tricomi's confluent hypergeometric function of the second kind, with
#' positive-integrand representation
#' \deqn{U(a,b,z)=\frac{1}{\Gamma(a)}\int_0^\infty
#'   e^{-zt}t^{a-1}(1+t)^{b-a-1}\,dt,\quad a>0,\ z>0.}
#' This integral exists for every real \eqn{b}; integer or nonpositive
#' \eqn{b} is not, by itself, a mathematical singularity of \eqn{U}.
#' The expected value of the distribution is
#' \deqn{E[x]=\mu.}
#' The variance is
#' \deqn{\sigma^2=\mu+\left((1+\alpha)
#' \left(2-\frac{2}{(\theta+2)^2}\right)-1\right)\mu^2.}
#'
#' \strong{GSL instability and automatic fallback.}
#' In the PMF above, both Tricomi terms have third argument
#' \deqn{z=\frac{\theta+2}{\alpha\mu(\theta+1)}.}
#' The ROOT discussion cited below reports GSL series-convergence failures
#' for \eqn{U(a,b,z)} when \eqn{a>8} and \eqn{z} is close to zero; an example
#' is \eqn{U(9.50606,0.5,0.000160903)}. The final correction in that thread
#' says \eqn{a>8}, superseding the earlier statement \eqn{a<10}. These are
#' historical implementation failures, not singularities of the function.
#' The discussion also concerns Kummer's first-kind function at negative
#' arguments; those examples are not the basis of this PMF's guard, since
#' its Tricomi argument \eqn{z} is positive.
#'
#' Additional checks against positive-integrand quadrature with GSL 2.7.1
#' found underestimated error for large negative second arguments, e.g.
#' \code{x = 0, mean = 1000, theta = 1, alpha = 0.001}
#' (\eqn{1/\alpha=1000}, \eqn{z=1.5}), where the relative PMF discrepancy was
#' about \eqn{4\times10^{-8}} despite successful GSL status. The guard
#' \eqn{0<\alpha<0.01} conservatively excludes this large-negative-b regime;
#' it too is an implementation policy, not a mathematical boundary.
#'
#' In \code{method = "auto"}, both U calls are bypassed whenever either
#' first argument exceeds 8 and \eqn{z<1}. Here the first arguments are
#' \eqn{x+1} and \eqn{x+2}, so this guard applies to \eqn{x\ge7}.
#' The cutoff \code{1} is a conservative implementation choice around
#' the reported small-positive-argument failures, not a published universal
#' boundary. Tests with GSL 2.7.1 also found inaccurate, successful-status
#' results outside \eqn{z<0.01} (including \eqn{z} near 0.12), motivating
#' the wider \eqn{z<1} guard. Behavior depends on the linked GSL version
#' and on all three arguments. For numerical and computational safeguards,
#' quadrature is also used directly when \eqn{0<\alpha<0.01}, either first
#' argument exceeds \eqn{10^4}, or \eqn{z} cannot be represented as a positive
#' normal double.
#'
#' Otherwise, the two U terms are evaluated together with
#' \code{gsl::hyperg_U(..., give = TRUE, strict = TRUE)}. The complete PMF
#' is recomputed by quadrature if GSL throws an error or warning; returns
#' \code{NULL}, missing/malformed fields, \code{NA}, \code{NaN}, or infinite
#' values; reports nonzero status; returns a nonpositive or subnormal U
#' value; or supplies an invalid error estimate. Each estimated relative
#' U error must be at most \code{rel.tol / 4}. A conservative estimate of
#' roundoff in the log formula must also be at most \code{rel.tol / 4}.
#' Invalid component probabilities or a final probability outside
#' \eqn{(0,1]} trigger the same fallback. GSL errors are not exposed as
#' warnings when quadrature subsequently succeeds. An unavailable GSL
#' namespace also causes fallback when this file is sourced independently.
#'
#' Gamma ratios, powers, and the sum of the two positive terms are evaluated
#' on the log scale. A zero U result is treated as possible underflow, never
#' as proof of a zero probability. Neither success status nor an error
#' estimate proves accuracy for every argument. Use \code{method =
#' "quadrature"} for independent checks in the parameter range of an
#' application, especially after changing the linked GSL version.
#'
#' \strong{Adaptive quadrature and cumulative probabilities.}
#' Write \eqn{X=\theta L}. With probability \eqn{\theta/(1+\theta)},
#' \eqn{X\sim\mathrm{Gamma}(1,1)}; otherwise,
#' \eqn{X\sim\mathrm{Gamma}(2,1)}. With
#' \eqn{c=(\theta+2)/(\theta+1)}, \eqn{W=X/c}. Quadrature integrates
#' the conditional NB PMF separately against these two positive gamma
#' densities. Each integral uses \eqn{t=\log X}, is scaled at its mode,
#' and is integrated on both sides of the mode with separate width scales.
#' Component probabilities are combined using log-sum-exp. The fallback
#' does not use a difference of hypergeometric functions, truncate the
#' count support, or substitute an arbitrary positive probability floor.
#'
#' \code{pplindGamma} always uses quadrature of the conditional NB CDF or
#' survival function. Upper tails are integrated directly; probabilities
#' near one are obtained from the opposite small tail. At \code{alpha = 0},
#' the conditional routines are Poisson and GSL is bypassed. Small positive
#' \code{alpha} is never silently replaced by zero. If quadrature itself
#' fails, the affected result is \code{NaN} with a warning; extreme inputs
#' need not be resolvable in double precision. On the ordinary probability
#' scale, sufficiently small probabilities can still underflow to zero;
#' use \code{log = TRUE} or \code{log.p = TRUE} for such tails.
#'
#' \code{qplindGamma} uses sequential inversion of \code{pplindGamma} with
#' its default integration controls, returning the smallest nonnegative
#' integer with CDF at least \code{p}. The endpoints are zero for
#' \code{p = 0} and infinity for \code{p = 1} when \code{mean > 0}; for
#' \code{mean = 0}, every valid quantile is zero. \code{rplindGamma} retains
#' inverse-CDF generation using uniform draws and this quantile function.
#' Quantiles and simulation therefore do not benefit from the GSL PMF path
#' and can be slow for large means or probabilities near one.
#'
#' Numeric arguments of \code{dplindGamma} and \code{pplindGamma} are
#' recycled to their maximum length; any zero-length numeric argument
#' produces \code{numeric(0)}. Their logical flags and numerical controls
#' must be scalar. Missing values propagate. Invalid distribution parameters
#' produce \code{NaN} with a warning. The quantile function uses
#' \code{Vectorize} recycling and simplification; random generation requires
#' scalar parameters. Names on \code{x} or \code{q} are retained when their
#' length equals the output length.
#'
#' @returns \code{dplindGamma} returns probability masses or log masses;
#'   \code{pplindGamma} returns cumulative or survival probabilities, possibly
#'   on the log scale; \code{qplindGamma} returns quantiles; and
#'   \code{rplindGamma} returns \code{n} random counts.
#'
#' @references
#' ROOT Forum, \emph{Problem in Hypergeometric functions with GSL}, including
#' the April 19, 2010 correction and October 3, 2014 follow-up:
#' \url{https://root-forum.cern.ch/t/problem-in-hypergeometric-functions-with-gsl/9507}.
#'
#' GNU Scientific Library Reference Manual, Hypergeometric Functions:
#' \url{https://www.gnu.org/software/gsl/doc/html/specfunc.html}.
#'
#' @seealso \code{\link[gsl:Hyperg]{hyperg_U}}, \code{\link[stats]{integrate}},
#'   \code{\link[stats]{dnbinom}}
#'
#' @examples
#' dplindGamma(0:5, mean = 0.75, theta = 7, alpha = 2)
#'
#' # Compare the hybrid evaluator with the quadrature reference.
#' x <- 0:10
#' fast <- dplindGamma(x, mean = 2, theta = 1, alpha = 0.5)
#' reference <- dplindGamma(x, mean = 2, theta = 1, alpha = 0.5,
#'                         method = "quadrature")
#' max(abs(fast - reference))
#'
#' # x >= 7 and z = (theta+2)/(alpha*mean*(theta+1)) < 1:
#' # the small-z guard sends this PMF directly to quadrature.
#' dplindGamma(8, mean = 1e4, theta = 1, alpha = 2/3, log = TRUE)
#'
#' # Log probabilities remain useful when ordinary probabilities underflow.
#' dplindGamma(10000, mean = 0.1, theta = 1, alpha = 0.5, log = TRUE)
#' pplindGamma(10, mean = 0.75, theta = 3, alpha = 0.5,
#'             lower.tail = FALSE, log.p = TRUE)
#' dplindGamma(0:3, mean = 0.75, theta = 2, alpha = 0)
#' qplindGamma(c(0, 0.5, 0.9, 1), mean = 1.67, theta = 0.5, alpha = 0.5)
#' rplindGamma(5, mean = 0.5, theta = 0.5, alpha = 2)
#'
#' @importFrom stats runif
#' @importFrom gsl hyperg_U
#' @useDynLib flexCountReg
#' @name NegativeBinomialLindley
#' @rdname NegativeBinomialLindley
#' @export
dplindGamma <- function(x, mean = 1, theta = 1, alpha = 1, log = FALSE,
                        rel.tol = 1e-8, subdivisions = 200L,
                        method = c("auto", "quadrature")) {
  method <- match.arg(method)
  .plg_flag(log, "log")
  .plg_control(rel.tol, subdivisions)
  .plg_vector(x, mean, theta, alpha, "pmf", TRUE, log,
              rel.tol, as.integer(subdivisions), method)
}


#' @rdname NegativeBinomialLindley
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
  vals <- list(p, mean, theta, alpha)
  if (!all(vapply(vals, is.numeric, logical(1))))
    stop("Probability and distribution arguments must be numeric")
  vals <- unlist(vals, use.names = FALSE)
  if (any(is.na(vals) & !is.nan(vals))) return(NA_real_)
  if (any(is.nan(vals))) return(NaN)
  if (any(!is.finite(vals)) || p < 0 || p > 1 ||
      mean < 0 || theta <= 0 || alpha < 0) {
    warning("Invalid probability or distribution parameter(s): NaNs produced",
            call. = FALSE)
    return(NaN)
  }
  if (mean == 0 || p == 0) return(0)
  if (p == 1) return(Inf)
  
  y <- 0
  p_value <- pplindGamma(y, mean, theta, alpha = alpha)
  while (is.finite(p_value) && p_value < p) {
    if (y + 1 == y) {
      warning("Quantile exceeds consecutive integer precision: NaN produced",
              call. = FALSE)
      return(NaN)
    }
    y <- y + 1
    p_value <- pplindGamma(y, mean, theta, alpha = alpha)
  }
  if (is.finite(p_value)) y else NaN
})

#' @rdname NegativeBinomialLindley
#' @export
rplindGamma <- function(n, mean=1, theta=1, alpha=1) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 0 || n != floor(n))
    stop("'n' must be a nonnegative integer scalar")
  parameters <- list(mean, theta, alpha)
  if (!all(vapply(parameters, function(x)
    is.numeric(x) && length(x) == 1L && is.finite(x), logical(1))) ||
    mean < 0 || theta <= 0 || alpha < 0)
    stop("'mean' and 'alpha' must be finite nonnegative scalars; ",
         "'theta' must be a finite positive scalar")
  if (n == 0) return(numeric())
  u <- stats::runif(n)
  y <- lapply(u, function(p) qplindGamma(p, mean, theta, alpha=alpha))
  return(unlist(y))
}
