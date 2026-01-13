#' Generate pseudo-random draws from specified distributions using Halton draws
#'
#' @param dist The distribution type to use. The distribution options include
#'   normal ("n"), lognormal ("ln"), triangular ("t"), uniform ("u"), and gamma
#'   ("g").
#' @param mean The mean value for the random draws.
#' @param sdev The standard deviation value for the random draws.
#' @param hdraw An optional vector of Halton draws to convert to the specified
#'   distribution. If not provided, the function will generate Halton draws.
#' @param ndraws The number of random draws to generate. This is only used if
#'   `hdraw` is not provided.
#'
#' @importFrom randtoolbox halton
#' @importFrom stats qlnorm qgamma qnorm
#' @include tri.R
#' @details
#' This function is used to convert Halton draws to the specified distribution.
#' The function can be used to generate random draws for use in random parameter
#' models, generating Halton-based pseudo-random draws for specified
#' distributions, etc.
#'
#' The distributions generated all use the `mean` (\deqn{\mu}) and `sdev`
#' (\deqn{\sigma}) parameters to generate the random draws. The density
#' functions for the distributions are as follows:
#' The Normal distribution is:
#' \eqn{f(x) = \frac{1}{\sqrt{2\pi\sigma^2}} 
#'   \exp\left(-\frac{(x - \mu)^2}{2\sigma^2}\right)}
#'
#' The Lognormal distribution is:
#' \eqn{f(x) = 
#'   \frac{1}{x\sigma\sqrt{2\pi}} 
#'   \exp\left(-\frac{(\log(x) - \mu)^2}{2\sigma^2}\right)}
#'
#' The Triangular distribution is (note that this is a symmetrical triangular
#' distribution where \deqn{\mu} is the median and \deqn{\sigma} is the
#' half-width):
#' \eqn{f(x) = \begin{cases} 
#'   \frac{(x - \mu + \sigma)}{\sigma^2}, & \text{for } \mu - 
#'     \sigma \leq x \leq \mu \\
#'   \frac{(\mu + \sigma - x)}{\sigma^2}, & \text{for } 
#'     \mu < x \leq \mu + \sigma \\0, & \text{otherwise}\end{cases}}
#'
#' The Uniform distribution is (note that \eqn{\mu} is the midpoint and
#' \eqn{\sigma} is the half-width):
#' \eqn{f(x) = \frac{1}{(\beta_{\mu}+\beta_{\sigma}) - 
#'   (\beta_{\mu}-\beta_{\sigma})}=\frac{1}{2\beta_{\sigma}}}
#'
#' The Gamma distribution is based on \deqn{\mu = \frac{\alpha}{\beta}} and
#' \deqn{\sigma^2 = \frac{\alpha}{\beta^2}}:
#' \eqn{f(x) = 
#'   \frac{\left(\frac{\mu}{\sigma^2}\right)^
#'     {\frac{\mu^2}{\sigma^2}}}{\Gamma\left(\frac{\mu^2}{\sigma^2}\right)} 
#'     x^{\frac{\mu^2}{\sigma^2} - 1} e^{-\frac{\mu}{\sigma^2} x}}
#'
#'
#' @returns A vector of psudo-random draws from the specified distribution, based 
#'  on Halton draws.
#' 
#' @examples
#' # Generate 500 random draws from a normal distribution 
#' halton_dists(dist="n", mean=3, sdev=2, ndraws=500)
#'
#' # Generate 500 random draws from a lognormal distribution 
#' halton_dists(dist="ln", mean=2, sdev=1.5, ndraws=500)
#'
#' # Generate 500 random draws from a triangular distribution
#' halton_dists(dist="t", mean=1, sdev=0.5, ndraws=500)
#'
#' # Generate 500 random draws from a uniform distribution
#' halton_dists(dist="u", mean=8, sdev=3, ndraws=500)
#'
#' # Generate 500 random draws from a gamma distribution
#' halton_dists(dist="g", mean=0.5, sdev=1.5, ndraws=500)
#'
#' @export
halton_dists <- function(dist, mean, sdev, hdraw=NULL, ndraws=500) {
  if (is.null(hdraw)) { # If Halton draws are not provided, generate them
    hdraw <- randtoolbox::halton(ndraws, 2)
  }
  switch(dist,
         "ln" = stats::qlnorm(hdraw, mean, abs(sdev)),
         "t" = qtri(hdraw, mean, abs(sdev)),
         "u" = mean + (hdraw - 0.5) * abs(sdev),
         "g" = stats::qgamma(hdraw, 
                             shape = mean^2 / sdev^2, 
                             rate = mean / sdev^2),
         # default case for normal distribution
         "n" = stats::qnorm(hdraw, mean, abs(sdev)), 
         warning("Invalid distribution type")
  )
}

