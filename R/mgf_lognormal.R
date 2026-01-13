#' Moment Generating Function for a Lognormal Distribution
#'
#' Computes the value of the moment generating function (MGF) for a lognormal
#' distribution at a given point through numerical integration. This function is
#' particularly useful for distributions where the MGF does not have a
#' closed-form solution. The lognormal distribution is specified by its log-mean
#' (\eqn{\mu}) and log-standard deviation (\eqn{\sigma}).
#'
#' @param mu The mean of the log-transformed variable, corresponding to
#'   \eqn{\mu} in the lognormal distribution's parameters.
#' @param sigma The standard deviation of the log-transformed variable,
#'   corresponding to \eqn{\sigma} in the lognormal distribution's parameters.
#' @param n The point at which to evaluate the MGF, often denoted as \eqn{t} in
#'   the definition of the MGF. This parameter essentially specifies the order
#'   of the moment generating function.
#' @returns The estimated value of the moment generating function (MGF) for the
#'   specified lognormal distribution at the given point.
#' @importFrom stats integrate
#' @details
#' The moment generating function (MGF) for the lognormal distribution does not
#' have a closed form solution. The MGF is defined as:
#' \deqn{
#' M_x(n) = \int_0^\infty e^{nx}\frac{1}{x\sigma\sqrt{2\pi}}
#'          e^{-\frac{(\ln(x)-\mu)^2}{2\sigma^2}}\,dx
#' }
#' 
#' The MGF for the lognormal distribution is useful for adjusting the
#' predictions of generalized linear mixed models (GLMMs) that have parameters
#' that follow a lognormal distribution and use a log link function. The
#' adjustment for the mean value is the MGF with \eqn{n=1} or
#' \eqn{E[e^x]=M_x(n=1)}. The variance for the lognormal random parameter is:
#' \deqn{Var(e^x)=E[e^{2x}]-E[e^x]^2=M_x(n=2)-M_x(n=1)^2}
#' 
#' @examples
#' mu <- 0
#' sigma <- 1
#' n <- 1
#' mgf_value <- mgf_lognormal(mu, sigma, n)
#' print(mgf_value)
#' @export
#' @name mgf_lognormal
mgf_lognormal <- function(mu, sigma, n) {
  # Check if necessary package is available
  if (!requireNamespace("stats", quietly = TRUE)) 
    warning("Package 'stats' is required but is not installed.")
  
  # Define the integrand function for the MGF
  integrand <- function(y, n, mu, sigma) {
    density <-  
      1 / (y * sigma * sqrt(2 * pi)) * exp(-((log(y) - mu)^2 / (2 * sigma^2)))
    return(density * exp(n * y))
  }
  upper.limit <- max(mu+sigma*10,100)
  
  # Perform numerical integration using the integrate() function
  result <- integrate(
    integrand, lower = 0, upper = upper.limit, n = n, mu = mu, sigma = sigma)
  
  if (result$message != "OK") {
    msg <- paste("Numerical integration may not be accurate.",
                 "Consider adjusting the limits.")
    warning(msg)
  }
  
  return(result$value)
}