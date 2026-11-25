#' Poisson-Weibull Distribution Functions
#'
#' These functions provide density, distribution function, quantile function,
#' and random number generation for the Poisson-Weibull Distribution, which is
#' specified either by its shape and scale parameters or by its mean and
#' standard deviation.
#'
#' The Poisson-Weibull distribution uses the Weibull distribution as a mixing
#' distribution for a Poisson process. It is useful for modeling overdispersed
#' count data. The density function (probability mass function) for the
#' Poisson-Weibull distribution is given by:
#' \deqn{P(y|\lambda,\alpha,\sigma) = 
#' \int_0^\infty \frac{e^{-\lambda x} \lambda^y x^y }{y!} 
#'   \left(\frac{\alpha}{\sigma}\right) 
#'   \left(\frac{x}{\sigma}\right)^{\alpha-1}
#'   e^{-\left(\frac{x}{\sigma}\right)^\alpha} dx}
#' where \eqn{f(x| \alpha, \sigma)} is the PDF of the Weibull distribution and
#' \eqn{\lambda} is the mean of the Poisson distribution.
#'
#' @param x A numeric value or vector of values for which the PDF or CDF is
#'   calculated.
#' @param q Quantile or a vector of quantiles.
#' @param p A numeric value or vector of probabilities for the quantile
#'   function.
#' @param n The number of random samples to generate.
#' @param alpha Shape parameter of the Weibull distribution (optional if mean
#'   and sd are provided).
#' @param sigma Scale parameter of the Weibull distribution (optional if mean
#'   and sd are provided).
#' @param mean_value Mean of the Weibull distribution (optional if alpha and
#'   sigma are provided).
#' @param sd_value Standard deviation of the Weibull distribution (optional if
#'   alpha and sigma are provided).
#' @param lambda Mean value of the Poisson distribution.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param log.p Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE, probabilities are \eqn{P[X <= x]}, otherwise
#'   \eqn{P[X > x]}.
#' @param ndraws the number of Halton draws (only with \code{engine = "halton"}).
#' @param engine the numerical method: \code{"quadrature"} (default) computes
#'   the log-probability by adaptive quadrature over the log of the mixing
#'   variable, accurate for every count including far tails;
#'   \code{"halton"} averages the Poisson probabilities over Halton draws
#'   (the method of flexCountReg 0.1.2 and earlier, kept for
#'   reproducibility; badly biased far in the tail).
#'
#' @details
#' \code{dpoisweibull} computes the density of the Poisson-Weibull distribution.
#'
#' \code{ppoisweibull} computes the distribution function of the Poisson-Weibull
#' distribution.
#'
#' \code{qpoisweibull} computes the quantile function of the Poisson-Weibull
#' distribution.
#'
#' \code{rpoisweibull} generates random numbers following the Poisson-Weibull
#' distribution.
#'
#' The shape and scale parameters directly define the Weibull distribution,
#' whereas the mean and standard deviation are used to compute these parameters
#' indirectly.
#' 
#' @returns dpoisweibull gives the density, ppoisweibull gives the distribution 
#'  function, qpoisweibull gives the quantile function, and rpoisweibull 
#'  generates random  deviates.
#' 
#'  The length of the result is determined by n for rpoisweibull, and is the 
#'  maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#' dpoisweibull(4, lambda=1.5, mean_value=1.5, sd_value=0.5, ndraws=10)
#' ppoisweibull(4, lambda=1.5, mean_value=1.5, sd_value=2, ndraws=10)
#' qpoisweibull(0.95, lambda=1.5, mean_value=1.5, sd_value=2, ndraws=10)
#' rpoisweibull(10, lambda=1.5, mean_value=1.5, sd_value=2, ndraws=10)
#'
#' @importFrom stats runif
#' @importFrom randtoolbox halton
#' @export
#' @name PoissonWeibull
#' 
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @rdname PoissonWeibull
#' @export
dpoisweibull <- function(x, lambda = NULL, alpha = NULL, sigma = NULL, 
                         mean_value = NULL, sd_value = NULL, 
                         ndraws = 1500, log = FALSE,
                         engine = c("quadrature", "halton")) {
  engine <- match.arg(engine)
  
  # --- Input Parameter Logic ---
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    sigma <- params[2]
  }
  
  if (is.null(alpha) || is.null(sigma)) 
    warning("Parameters alpha and sigma are required")
  
  from_mean <- is.null(lambda) && !is.null(mean_value)
  if (is.null(lambda)) {
    if (!is.null(mean_value)) {
      lambda <- mean_value / (sigma * gamma(1 + 1 / alpha))
    } else {
      warning("Parameter lambda or mean_value is required")
    }
  }
  
  if (engine == "quadrature") {
    lp <- .pw_quad_logpmf(x, if (from_mean) mean_value else lambda, alpha,
                          sigma, from_mean)
    .fcr_warn_nan(is.nan(lp))
    return(if (log) lp else exp(lp))
  }

  # --- Vectorization Setup ---
  # Ensure parameters are vectors of appropriate length matches x
  n <- length(x)
  if(length(lambda) == 1) lambda <- rep(lambda, n)
  if(length(alpha) == 1) alpha <- rep(alpha, n)
  if(length(sigma) == 1) sigma <- rep(sigma, n)
  
  # --- Optimization ---
  # Generate Halton draws ONCE for the whole batch
  h <- randtoolbox::halton(ndraws)
  
  # Assuming dpWeib_cpp is vectorized or we map over it. 
  # Ideally, dpWeib_cpp should accept vectors. If strictly scalar C++, 
  # we must use mapply here, but passing 'h' prevents regenerating it.
  
  # The C++ kernel is parameterized by the mean of Y,
  # E[Y] = lambda * sigma * Gamma(1 + 1/alpha)
  mean_y <- lambda * sigma * gamma(1 + 1 / alpha)
  p <- dpWeib_cpp(x, mean_y, alpha, sigma, h)
  
  if (log) return(log(p))
  else return(p)
}

#' @rdname PoissonWeibull
#' @export
ppoisweibull <- function(q, lambda = NULL, alpha = NULL, sigma = NULL,
                         mean_value = NULL, sd_value = NULL, ndraws = 1500,
                         lower.tail = TRUE, log.p = FALSE,
                         engine = c("quadrature", "halton")) {
  engine <- match.arg(engine)
  if (engine == "quadrature") {
    if (!is.null(mean_value) && !is.null(sd_value)) {
      params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
      alpha <- params[1]
      sigma <- params[2]
    }
    from_mean <- is.null(lambda) && !is.null(mean_value)
    if (from_mean) lambda <- mean_value
    if (is.null(lambda) || is.null(alpha) || is.null(sigma)) {
      stop("`lambda` (or `mean_value`), `alpha` and `sigma` are required.",
           call. = FALSE)
    }
    # `lambda` holds the mean of Y when from_mean is TRUE
    core <- function(q, lambda, alpha, sigma, lower.tail, log.p) {
      .pw_quad_cdf(q, lambda, alpha, sigma, from_mean, lower.tail, log.p)
    }
    return(.fcr_p_contract(core, q,
                           list(lambda = lambda, alpha = alpha,
                                sigma = sigma), lower.tail, log.p,
                           logdens = function(x, lambda, alpha, sigma)
                             .pw_quad_logpmf(x, lambda, alpha, sigma,
                                             from_mean)))
  }
  .fcr_p_contract(.ppoisweibull_core, q,
                  list(lambda = lambda, alpha = alpha, sigma = sigma,
                       mean_value = mean_value, sd_value = sd_value),
                  lower.tail, log.p, ndraws = ndraws)
}

.ppoisweibull_core <- Vectorize(function(
    q, lambda=NULL, alpha = NULL, sigma = NULL, 
    mean_value = NULL, sd_value = NULL, 
    ndraws=1500, lower.tail = TRUE, log.p = FALSE) {
  
  # Parameter logic repeated ensures safety inside Vectorize
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    sigma <- params[2]
  }
  if (is.null(lambda) && !is.null(mean_value)){
    lambda <- mean_value/(sigma*gamma(1+1/alpha))
  }
  
  if (q < 0) return(if(log.p) -Inf else 0)
  
  # Efficiency: dpoisweibull is now faster
  y <- 0:floor(q)
  
  # Note: dpoisweibull now handles the heavy lifting
  probs <- 
    dpoisweibull(y, lambda=lambda, alpha=alpha, sigma=sigma, ndraws=ndraws,
                 engine = "halton")
  p <- sum(probs)
  
  if(!lower.tail) p <- 1-p
  if (log.p) return(log(p))
  else return(p)
})

#' @rdname PoissonWeibull
#' @export
qpoisweibull <- function(p, lambda = NULL, alpha = NULL, sigma = NULL,
                         mean_value = NULL, sd_value = NULL, ndraws = 1500,
                         engine = c("quadrature", "halton")) {
  engine <- match.arg(engine)
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    sigma <- params[2]
  }
  if (is.null(alpha) || is.null(sigma)) 
    stop("Parameters alpha and sigma are required")
  if (is.null(lambda) && !is.null(mean_value)){
    lambda <- mean_value/(sigma*gamma(1+1/alpha))
  } else if (is.null(lambda)) {
    stop("Parameter lambda or mean_value is required")
  }
  .fcr_qcount(p, function(y, lambda, alpha, sigma) {
    cumsum(suppressWarnings(dpoisweibull(y, lambda = lambda, alpha = alpha,
                                         sigma = sigma, ndraws = ndraws,
                                         engine = engine)))
  }, lambda, alpha, sigma)
}

#' @rdname PoissonWeibull
#' @export
rpoisweibull <- function(n, lambda=NULL, alpha = NULL, sigma = NULL, 
                         mean_value = NULL, sd_value = NULL, ndraws=1500) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    sigma <- params[2]
  }
  if (is.null(alpha) || is.null(sigma)) 
    warning("Parameters alpha and sigma are required")
  if (is.null(lambda) && !is.null(mean_value)){
    lambda <- mean_value/(sigma*gamma(1+1/alpha))
  } else if (is.null(lambda) && is.null(mean_value)) {
    warning("Parameter lambda or mean_value is required")
  }
  
  if (is.null(alpha) || is.null(sigma)) 
    warning("Parameters alpha and sigma are required")
  
  # Use lambda, alpha, and sigma to generate random numbers using weibull random
  # numbers and Poisson distribution
  eta <- stats::rweibull(n, shape = alpha, scale = sigma)
  y <- stats::rpois(n, lambda = lambda * eta)
  return(y)
  # u <- stats::runif(n)
  # y <- 
  #   qpoisweibull(u, lambda=lambda, alpha = alpha, sigma = sigma, ndraws=ndraws)
  # return(y)
}


# log P(Y = x) by quadrature (src/mixtureQuad.cpp), recycled. With
# from_mean = TRUE, `m` is the mean of Y (the Weibull scale sigma cancels but
# must be valid); otherwise `m` is lambda and the mixing variable is
# lambda * sigma * Weibull(alpha, 1), passed to the C++ code as a scale
# (kind 5) so that shapes for which Gamma(1 + 1 / alpha) overflows (alpha
# below about 0.006) still work. Invalid parameters give NaN, missing values
# NA.
.pw_quad_logpmf <- function(x, m, alpha, sigma, from_mean) {
  r <- .fcr_recycle(x, m, alpha, sigma)
  if (length(r[[1]]) == 0L) return(numeric(0))
  m <- r[[2]]; al <- r[[3]]; sg <- r[[4]]
  scale <- if (from_mean) m else m * sg
  lp <- mix_logpmf_cpp(as.numeric(r[[1]]), as.numeric(scale),
                       as.numeric(al), numeric(length(al)),
                       if (from_mean) 2L else 5L)
  bad <- !is.na(m) & !is.na(al) & !is.na(sg) &
    !(is.finite(m) & m > 0 & is.finite(al) & al > 0 & is.finite(sg) &
        sg > 0 & is.finite(scale))
  lp[bad & !is.na(r[[1]])] <- NaN
  lp
}

# P(Y <= q) (or P(Y > q)) for finite q >= 0 by summing the quadrature PMF
# over 0, ..., floor(q), once per distinct parameter set (as .fcr_mix_cdf())
.pw_quad_cdf <- function(q, m, alpha, sigma, from_mean, lower.tail = TRUE,
                         log.p = FALSE) {
  r <- .fcr_recycle(q, m, alpha, sigma)
  q <- floor(r[[1]])
  out <- rep(NA_real_, length(q))
  key <- paste(format(r[[2]], digits = 17), format(r[[3]], digits = 17),
               format(r[[4]], digits = 17), sep = "\r")
  for (k in unique(key)) {
    idx <- which(key == k)
    i1 <- idx[1]
    lp <- .pw_quad_logpmf(0:max(q[idx]), r[[2]][i1], r[[3]][i1], r[[4]][i1],
                          from_mean)
    if (anyNA(lp)) {
      out[idx] <- if (any(is.nan(lp))) NaN else NA_real_
      next
    }
    top <- max(lp)
    lcdf <- if (top == -Inf) rep(-Inf, length(lp)) else
      top + log(cumsum(exp(lp - top)))
    lcdf <- pmin(lcdf[q[idx] + 1], 0)
    out[idx] <- if (lower.tail) {
      if (log.p) lcdf else exp(lcdf)
    } else {
      if (log.p) .fcr_log1mexp(lcdf) else -expm1(lcdf)
    }
  }
  .fcr_warn_nan(is.nan(out))
  out
}

# Weibull shape (alpha) and scale (sigma) with a given mean and standard
# deviation. The coefficient of variation cv = sd / mean depends only on the
# shape: log1p(cv^2) = lgamma(1 + 2 / alpha) - 2 lgamma(1 + 1 / alpha), which
# decreases in alpha. It is solved in log(alpha) with an expanding bracket,
# so any positive cv is reached (a search restricted to alpha in [0.1, 5]
# returned 5 for small cv, for example sd 0.229 instead of 0.01 at mean 1).
calculate_params <- function(mean_value = NULL, sd_value = NULL) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    if (!is.numeric(mean_value) || !is.numeric(sd_value) ||
        length(mean_value) != 1L || length(sd_value) != 1L ||
        !is.finite(mean_value) || !is.finite(sd_value) ||
        mean_value <= 0 || sd_value <= 0) {
      stop("`mean_value` and `sd_value` must be single positive numbers.",
           call. = FALSE)
    }
    target <- log1p((sd_value / mean_value)^2)
    f <- function(la) {
      a <- exp(la)
      lgamma(1 + 2 / a) - 2 * lgamma(1 + 1 / a) - target
    }
    root <- stats::uniroot(f, c(-1, 1), extendInt = "downX",
                           tol = 1e-12, maxiter = 1000)
    alpha <- exp(root$root)
    if (!is.finite(alpha) || abs(root$f.root) > 1e-8 * max(1, target)) {
      stop("No Weibull shape matches the mean and standard deviation.",
           call. = FALSE)
    }
    sigma <- mean_value / gamma(1 + 1/alpha)
    c(alpha, sigma)
  } else {
    msg <- paste(
      "Insufficient parameters: Provide either mean and standard deviation,",
      "or alpha and sigma."
    )
    warning(msg)
  }
}