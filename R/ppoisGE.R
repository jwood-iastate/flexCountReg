#' Poisson-Generalized-Exponential Distribution
#'
#' These functions provide density, distribution function, quantile function,
#' and random number generation for the Poisson-Generalized-Exponential (PGE)
#' Distribution
#'
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mean numeric value or vector of mean values for the distribution (the
#'   values have to be greater than 0). This is NOT the value of \eqn{\lambda}.
#' @param shape numeric value or vector of shape values for the shape parameter
#'   of the generalized exponential distribution (the values have to be greater
#'   than 0).
#' @param scale single value or vector of values for the rate \eqn{b} of the
#'   generalized exponential mixing distribution (greater than 0; despite
#'   its name, the mixing variable is generated as \eqn{X/b} with
#'   \eqn{X \sim GE(a, 1)}). Because the distribution is parameterized by
#'   its mean, \eqn{b} cancels and does not change any result; it is only
#'   checked for validity.
#' @param ndraws the number of Halton draws (only with \code{engine = "halton"}).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]}
#'   otherwise, \eqn{P[X>x]}.
#' @param haltons an optional vector of Halton draws to use instead of ndraws
#'   (only with \code{engine = "halton"}).
#' @param engine the numerical method: \code{"quadrature"} (default) computes
#'   the log-probability by adaptive quadrature over the log of the mixing
#'   variable, accurate for every count including far tails;
#'   \code{"halton"} averages the Poisson probabilities over Halton draws
#'   (the method of flexCountReg 0.1.2 and earlier, kept for
#'   reproducibility; badly biased far in the tail; the default when
#'   \code{haltons} are given).
#'
#' @details
#' \code{dpge} computes the density (PDF) of the PGE Distribution.
#'
#' \code{ppge} computes the CDF of the PGE Distribution.
#'
#' \code{qpge} computes the quantile function of the PGE Distribution.
#'
#' \code{rpge} generates random numbers from the PGE Distribution.
#' 
#' The mixing variable \eqn{X} has a generalized exponential distribution
#' (Gupta and Kundu, 2007) with shape \eqn{a > 0} and rate \eqn{b > 0}:
#' \deqn{f(x \mid a, b) = a b \left(1 - e^{-bx}\right)^{a-1} e^{-bx},
#'   \quad x > 0,}
#' with
#' \deqn{E[X] = \frac{A}{b}, \quad Var[X] = \frac{B}{b^2}, \quad
#'   A = \psi(a+1) - \psi(1), \quad B = \psi'(1) - \psi'(a+1),}
#' where \eqn{\psi} and \eqn{\psi'} are the digamma and trigamma
#' functions. Given \eqn{X}, \eqn{Y \sim Poisson(\lambda X)}, so
#' \deqn{P(Y = y) = \int_0^\infty \frac{(\lambda x)^y e^{-\lambda x}}{y!}
#'   f(x \mid a, b)\, dx.}
#' To make \eqn{E[Y] = \mu} (\code{mean}, \eqn{e^{X\beta}} in a regression),
#' \eqn{\lambda = \mu b / A}. Then
#' \deqn{E[Y] = \mu, \qquad Var[Y] = \mu + \frac{B}{A^2}\mu^2.}
#' Since \eqn{\lambda X = \mu (bX) / A} and \eqn{bX \sim GE(a, 1)}, the
#' distribution of \eqn{Y} depends on \eqn{\mu} and \eqn{a} only: the rate
#' (or scale) cancels. This is why PGE regression models estimate only the
#' shape (\code{ln(shape)}); a scale parameter would not be identified.
#' With \eqn{a = 1} the mixing distribution is exponential and \eqn{Y} is
#' geometric with mean \eqn{\mu}.
#'
#' The integral has no closed form. By default (\code{engine =
#' "quadrature"}) it is computed on the log scale by adaptive quadrature over
#' \eqn{\log X}, centered at the mode of the integrand, which is accurate for
#' every count. \code{engine = "halton"} averages the Poisson probabilities
#' over Halton draws of \eqn{X}, as flexCountReg 0.1.2 and earlier did.
#'
#' @references
#' Gupta, R. D., & Kundu, D. (2007). Generalized exponential distribution:
#' Existing results and some recent developments. Journal of Statistical
#' planning and inference, 137(11), 3537-3547.
#' 
#' @returns dpge gives the density, ppge gives the distribution 
#'  function, qpge gives the quantile function, and rpge generates
#'  random  deviates.
#' 
#'  The length of the result is determined by n for rpge, and is the 
#'  maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#' dpge(0, mean=0.75, shape=2, scale=1, ndraws=2000)
#' ppge(c(0,1,2,3,4,5,6), mean=0.75, shape=2, scale=1, ndraws=500)
#' qpge(c(0.1,0.3,0.5,0.9,0.95), mean=0.75, shape=2, scale=1, ndraws=500)
#' rpge(30, mean=0.75,  shape=2, scale=1, ndraws=500)
#'
#' @importFrom stats runif
#' @importFrom randtoolbox halton
#' @name PoissonGeneralizedExponential
NULL

# Internal: PGE probability for one observation (not exported)
.dpge_one <- function(x, mean, shape, scale, haltons) {
  
  lambda <- mean * scale /
    (digamma(shape + 1) - digamma(1))
  
  gedist <- -log1p(-haltons^(1 / shape)) / scale
  lambda_i <- lambda * gedist
  
  base::mean(
    stats::dpois(
      x = x,
      lambda = lambda_i
    )
  )
}

#' @rdname PoissonGeneralizedExponential
#' @export
dpge <- function(
    x,
    mean = 1,
    shape = 1,
    scale = 1,
    ndraws = 1500,
    log = FALSE,
    haltons = NULL,
    engine = c("quadrature", "halton")
) {
  # Draws given without an engine: the draw-based method
  if (missing(engine) && !is.null(haltons)) engine <- "halton"
  engine <- match.arg(engine)
  if (engine == "quadrature") {
    # The scale cancels from the mean-parameterized PMF; it must still be
    # valid
    r <- .fcr_recycle(x, mean, shape, scale)
    if (length(r[[1]]) == 0L) return(numeric(0))
    sc <- r[[4]]
    bad_sc <- !is.na(sc) & !(is.finite(sc) & sc > 0)
    lp <- suppressWarnings(.fcr_mix_logpmf(r[[1]], r[[2]], r[[3]], 0, "PGE"))
    lp[is.na(sc) & !is.na(r[[1]])] <- NA
    lp[bad_sc & !is.na(lp)] <- NaN
    .fcr_warn_nan(is.nan(lp) & !is.na(r[[1]]))
    return(if (isTRUE(log)) lp else exp(lp))
  }
  
  # Validate observation-specific parameters
  if (any(!is.finite(mean) | mean <= 0)) {
    stop("mean must be finite and greater than 0")
  }
  
  if (any(!is.finite(shape) | shape <= 0)) {
    stop("shape must be finite and greater than 0")
  }
  
  if (any(!is.finite(scale) | scale <= 0)) {
    stop("scale must be finite and greater than 0")
  }
  
  if (length(log) != 1L || is.na(log)) {
    stop("log must be a single TRUE or FALSE value")
  }
  
  # Validate or generate integration draws
  if (is.null(haltons)) {
    
    if (
      length(ndraws) != 1L ||
      !is.finite(ndraws) ||
      ndraws < 1 ||
      ndraws != floor(ndraws)
    ) {
      stop("ndraws must be a positive integer")
    }
    
    haltons <- randtoolbox::halton(n = as.integer(ndraws))
    
  } else {
    
    if (
      length(haltons) == 0L ||
      any(!is.finite(haltons) | haltons <= 0 | haltons >= 1)
    ) {
      stop("haltons must be nonempty, finite, and strictly between 0 and 1")
    }
  }
  
  # Preserve ordinary R-style recycling across observations.
  parameter_lengths <- c(
    length(x),
    length(mean),
    length(shape),
    length(scale)
  )
  
  if (any(parameter_lengths == 0L)) {
    return(numeric(0))
  }
  
  n <- max(parameter_lengths)
  
  x_vec     <- rep_len(x, n)
  mean_vec  <- rep_len(mean, n)
  shape_vec <- rep_len(shape, n)
  scale_vec <- rep_len(scale, n)
  
  probabilities <- vapply(
    seq_len(n),
    FUN = function(i) {
      .dpge_one(
        x = x_vec[i],
        mean = mean_vec[i],
        shape = shape_vec[i],
        scale = scale_vec[i],
        haltons = haltons
      )
    },
    FUN.VALUE = numeric(1)
  )
  
  if (isTRUE(log)) {
    base::log(probabilities)
  } else {
    probabilities
  }
}

#' @rdname PoissonGeneralizedExponential
#' @export
ppge <- function(q, mean = 1, shape = 1, scale = 1, ndraws = 1500,
                 lower.tail = TRUE, log.p = FALSE, haltons = NULL,
                 engine = c("quadrature", "halton")) {
  if (missing(engine) && !is.null(haltons)) engine <- "halton"
  engine <- match.arg(engine)
  .fcr_p_contract(.ppge_core, q, list(mean = mean, shape = shape, scale = scale),
                  lower.tail, log.p, ndraws = ndraws, haltons = haltons,
                  engine = engine,
                  logdens = if (engine == "quadrature") function(x, ...)
                    dpge(x, ..., log = TRUE))
}

.ppge_core <- function(q, mean = 1, shape = 1, scale = 1, ndraws = 1500,
                 lower.tail = TRUE, log.p = FALSE, haltons = NULL,
                 engine = "quadrature") {
  if (engine == "quadrature") {
    r <- .fcr_recycle(q, mean, shape, scale)
    out <- .fcr_mix_cdf(r[[1]], r[[2]], r[[3]], 0, "PGE", lower.tail, log.p)
    out[!is.na(r[[4]]) & !(is.finite(r[[4]]) & r[[4]] > 0)] <- NaN
    return(out)
  }
  
  
  # --- Input Validation ---
  
  if (any(mean <= 0, na.rm = TRUE)) warning("'mean' must be positive")
  if (any(shape <= 0, na.rm = TRUE)) warning("'shape' must be positive")
  if (any(scale <= 0, na.rm = TRUE))  warning("'scale' must be positive")
  
  
  # --- Vectorization Setup ---
  
  n <- max(length(q), length(mean), length(shape), length(scale))
  q <- rep_len(as.integer(floor(q)), n)
  mean <- rep_len(mean, n)
  shape <- rep_len(shape, n)
  scale <- rep_len(scale, n)
  
  
  # --- Initialize Result ---
  
  cdf <- rep(NA_real_, n)
  
  
  # --- Handle Special Cases ---
  
  invalid_q <- is.na(q) | q < 0
  cdf[invalid_q] <- 0
  
  
  # --- Pre-compute Halton draws (shared across all evaluations) ---
  
  if (is.null(haltons)) {
    haltons <- randtoolbox::halton(ndraws)
  }
  
  
  # --- Compute CDF for Valid Cases ---
  
  bad_par <- !(is.finite(mean) & mean > 0 & is.finite(shape) & shape > 0 &
                is.finite(scale) & scale > 0)
  cdf[bad_par & !invalid_q] <- NaN
  valid <- !invalid_q & !bad_par
  
  if (any(valid)) {
    # Group by unique parameter combinations
    param_df <- data.frame(
      idx = which(valid),
      q = q[valid],
      mean = mean[valid],
      shape = shape[valid],
      scale = scale[valid]
    )
    
    unique_params <- unique(param_df[, c("mean", "shape", "scale")])
    
    for (j in seq_len(nrow(unique_params))) {
      m_j <- unique_params$mean[j]
      sh_j <- unique_params$shape[j]
      sc_j <- unique_params$scale[j]
      
      # Find all indices with these parameters
      mask <- 
        param_df$mean == m_j & param_df$shape == sh_j & param_df$scale == sc_j
      indices <- param_df$idx[mask]
      q_vals <- param_df$q[mask]
      max_q <- max(q_vals)
      
      # Compute lambda for this parameter set
      lambda_j <- m_j * sc_j / (digamma(sh_j + 1) - digamma(1))
      
      # Generate GE quantiles from Halton draws
      # Quantile function of GE: Q(p) = -log(1 - p^(1/shape)) / scale
      ge_quantiles <- -log(1 - haltons^(1 / sh_j)) / sc_j
      
      # Compute mu values for each Halton draw
      mu_draws <- lambda_j * ge_quantiles
      
      # Compute PMF for 0:max_q by averaging over Halton draws
      pmf <- numeric(max_q + 1)
      
      for (y in 0:max_q) {
        # Average Poisson PMF over all draws
        pmf[y + 1] <- mean(stats::dpois(y, mu_draws))
      }
      
      # Compute cumulative sums
      cum_pmf <- cumsum(pmf)
      
      # Assign CDF values
      for (k in seq_along(indices)) {
        idx_k <- indices[k]
        q_k <- q_vals[k]
        cdf[idx_k] <- cum_pmf[q_k + 1]
      }
    }
  }
  
  
  # --- Apply Tail and Log Options ---
  
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  
  return(cdf)
}
#' @rdname PoissonGeneralizedExponential
#' @export
qpge <- function(p, mean = 1, shape = 1, scale = 1, ndraws = 1500,
                 engine = c("quadrature", "halton")) {
  engine <- match.arg(engine)
  if (any(mean <= 0 | scale <= 0 | shape <= 0, na.rm = TRUE)) {
    warning("The values of `mean`, `shape`, and `scale` have to have ",
            "values greater than 0.")
  }
  h <- randtoolbox::halton(ndraws)
  .fcr_qcount(p, function(y, mean, shape, scale) {
    cumsum(suppressWarnings(dpge(y, mean = mean, shape = shape,
                                 scale = scale, haltons = h,
                                 engine = engine)))
  }, mean, shape, scale)
}


#' @rdname PoissonGeneralizedExponential
#' @export
rpge <- function(n, mean = 1, shape = 1, scale = 1, ndraws = 1500) {
  n <- .fcr_rn(n)
  mean <- rep_len(mean, n); shape <- rep_len(shape, n)
  scale <- rep_len(scale, n)
  bad <- !(is.finite(mean) & mean > 0 & is.finite(shape) & shape > 0 &
             is.finite(scale) & scale > 0)
  .fcr_warn_na(bad, "`mean`, `shape` and `scale` must be greater than 0.")
  out <- rep(NA_real_, n); ok <- !bad; k <- sum(ok)
  # Generalized-exponential mixing variable by inverse CDF, scaled so that
  # E[Y] = mean (as in dpge()).
  g <- -log1p(-stats::runif(k)^(1 / shape[ok])) / scale[ok]
  lambda <- mean[ok] * scale[ok] / (digamma(shape[ok] + 1) - digamma(1))
  out[ok] <- stats::rpois(k, lambda * g)
  out
}
