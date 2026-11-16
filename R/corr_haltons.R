#' Generate Correlated Random Variables Using Halton or Scrambled Halton Draws
#'
#' This function generates \code{N} correlated random variables using Halton or
#' scrambled Halton draws. The function supports normal and truncated normal
#' distributions.
#'
#' @param means A numeric vector of means for each variable.
#' @param cholesky A Cholesky decomposition matrix to introduce correlation.
#' @param stdev A numeric vector of standard deviations for each variable. If
#'   provided, the function will use these values instead of the Cholesky
#'   decomposition matrix (must also provide a correlation matrix if providing
#'   standard deviations). Default is NULL.
#' @param correlations A correlation matrix to introduce correlation. If
#'   provided, the function will use these values instead of the Cholesky
#'   decomposition matrix (must also provide standard deviations). Default is
#'   NULL.
#' @param hdraws A matrix of Halton or scrambled Halton draws. If provided, the
#'   function will use these draws instead of generating new ones. Default is
#'   NULL.
#' @param ndraws An integer specifying the number of values to simulate for each
#'   variable. Default is 500.
#' @param scrambled A logical value indicating whether to use scrambled Halton
#'   draws. Default is FALSE.
#' @param dist A character string specifying the distribution type. Options are
#'   "normal" and "truncated_normal". Default is "normal".
#' @param lower A numeric value specifying the lower bound for truncated normal
#'   distribution. Default is -Inf.
#' @param upper A numeric value specifying the upper bound for truncated normal
#'   distribution. Default is Inf.
#'
#' @returns A matrix with \code{N} columns and \code{ndraws} rows containing the
#'   simulated values for the correlated random variables.
#'   
#' @importFrom randtoolbox halton
#' @importFrom truncnorm qtruncnorm
#' @importFrom stats qlnorm pnorm
#' @include cor2cov.R
#' @examples
#' # Define mean, correlation, and standard deviations
#' means <- c(3, 2, 0.9)
#' sdevs <- c(0.25,1.5,0.8)
#' CORR <- matrix(c(1, -0.3, 0.5, -0.3, 1, -0.2, 0.5, -0.2, 1), 3, 3)
#'
#' # Create the Cholesky decomposition matrix and set values for ndraws, etc.
#' ndraws <- 5000
#' scrambled <- TRUE
#' dist <- "normal"
#'
#' # simulated the data
#' simulated_data <- corr_haltons(means, stdev=sdevs, correlations=CORR,
#'                                 ndraws=ndraws, scrambled=scrambled,
#'                                 dist=dist)
#'
#' # look at the mean, standard deviation, and correlation of the simulated data
#' apply(simulated_data, 2, mean)
#' apply(simulated_data, 2, sd)
#' cor(simulated_data)
#'
#' # providing a cholesky decomposition matrix
#' dist <- "normal"
#' cholesky <- chol(cor2cov(CORR, sdevs))
#' simulated_data <- corr_haltons(means, cholesky=cholesky, ndraws=ndraws,
#'                                 scrambled=scrambled, dist=dist)
#' apply(simulated_data, 2, mean)
#' apply(simulated_data, 2, sd)
#' cor(simulated_data)
#'
#' # Truncated normal
#' dist <- "truncated_normal"
#' lower <- 0
#' upper <- 30
#' simulated_data <- corr_haltons(means, cholesky=cholesky, ndraws=ndraws,
#'                                 scrambled=scrambled, dist=dist,
#'                                 lower=lower, upper=upper)
#' apply(simulated_data, 2, mean)
#' apply(simulated_data, 2, sd)
#' cor(simulated_data)
#'
#' @export
corr_haltons <- function(
    means, cholesky=NULL, stdev=NULL, correlations=NULL, hdraws=NULL, 
    ndraws=500, scrambled=FALSE, dist="normal", lower=-Inf, upper=Inf) {
  N <- length(means)
  
  if (!(dist %in% c("normal", "truncated_normal"))) {
    stop("Unsupported distribution type (", dist, "). Valid options include ",
         "'normal' and 'truncated_normal'.", call. = FALSE)
  }

  if (!is.null(hdraws)) {
    hdraws <- as.matrix(hdraws)
    if (ncol(hdraws) != length(means)) {
      stop("The number of columns in `hdraws` (", ncol(hdraws),
           ") must be the same as the length of `means` (", length(means),
           ").", call. = FALSE)
    }
    halton_seq <- hdraws
  } else {
    # Generate Halton or scrambled Halton sequence (standard normal)
    halton_seq <- as.matrix(halton(ndraws, dim = N, normal = TRUE,
                                   mixed = scrambled))
  }

  if (!is.null(correlations) && !is.null(stdev)) {
    cholesky <- chol(cor2cov(correlations, stdev))
  } else if (is.null(cholesky)) {
    stop("You must provide either a Cholesky decomposition matrix ",
         "(`cholesky`) or standard deviations and a correlation matrix ",
         "(`stdev` and `correlations`).", call. = FALSE)
  }
  cholesky <- as.matrix(cholesky)
  if (!identical(dim(cholesky), c(N, N))) {
    stop("`cholesky` must be a ", N, " x ", N, " matrix (one row and column ",
         "per element of `means`).", call. = FALSE)
  }
  sdevs <- sqrt(diag(t(cholesky) %*% cholesky))

  # Apply the Cholesky factor to introduce the correlation: each column is
  # normal with mean 0 and standard deviation sdevs[j]
  correlated_seq <- halton_seq %*% cholesky

  if (dist == "normal") {
    result <- sweep(correlated_seq, 2, means, "+")
  } else {
    # Truncated normal: map each column through its own normal CDF and the
    # truncated-normal quantile function with that column's mean and SD
    lower <- rep_len(lower, N)
    upper <- rep_len(upper, N)
    result <- correlated_seq
    for (j in seq_len(N)) {
      u <- stats::pnorm(correlated_seq[, j], sd = sdevs[j])
      result[, j] <- qtruncnorm(u, a = lower[j], b = upper[j],
                                mean = means[j], sd = sdevs[j])
    }
  }

  return(result)
}
