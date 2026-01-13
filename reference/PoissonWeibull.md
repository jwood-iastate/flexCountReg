# Poisson-Weibull Distribution Functions

These functions provide density, distribution function, quantile
function, and random number generation for the Poisson-Weibull
Distribution, which is specified either by its shape and scale
parameters or by its mean and standard deviation.

## Usage

``` r
dpoisweibull(
  x,
  lambda = NULL,
  alpha = NULL,
  sigma = NULL,
  mean_value = NULL,
  sd_value = NULL,
  ndraws = 1500,
  log = FALSE
)

ppoisweibull(
  q,
  lambda = NULL,
  alpha = NULL,
  sigma = NULL,
  mean_value = NULL,
  sd_value = NULL,
  ndraws = 1500,
  lower.tail = TRUE,
  log.p = FALSE
)

qpoisweibull(
  p,
  lambda = NULL,
  alpha = NULL,
  sigma = NULL,
  mean_value = NULL,
  sd_value = NULL,
  ndraws = 1500
)

rpoisweibull(
  n,
  lambda = NULL,
  alpha = NULL,
  sigma = NULL,
  mean_value = NULL,
  sd_value = NULL,
  ndraws = 1500
)
```

## Arguments

- x:

  A numeric value or vector of values for which the PDF or CDF is
  calculated.

- lambda:

  Mean value of the Poisson distribution.

- alpha:

  Shape parameter of the Weibull distribution (optional if mean and sd
  are provided).

- sigma:

  Scale parameter of the Weibull distribution (optional if mean and sd
  are provided).

- mean_value:

  Mean of the Weibull distribution (optional if alpha and sigma are
  provided).

- sd_value:

  Standard deviation of the Weibull distribution (optional if alpha and
  sigma are provided).

- ndraws:

  the number of Halton draws to use for the integration.

- log:

  Logical; if TRUE, probabilities p are given as log(p).

- q:

  Quantile or a vector of quantiles.

- lower.tail:

  Logical; if TRUE, probabilities are P\[X \<= x\], otherwise P\[X \>
  x\].

- log.p:

  Logical; if TRUE, probabilities p are given as log(p).

- p:

  A numeric value or vector of probabilities for the quantile function.

- n:

  The number of random samples to generate.

## Value

dpoisweibull gives the density, ppoisweibull gives the distribution
function, qpoisweibull gives the quantile function, and rpoisweibull
generates random deviates.

The length of the result is determined by n for rpoisweibull, and is the
maximum of the lengths of the numerical arguments for the other
functions.

## Details

The Poisson-Weibull distribution uses the Weibull distribution as a
mixing distribution for a Poisson process. It is useful for modeling
overdispersed count data. The density function (probability mass
function) for the Poisson-Weibull distribution is given by:
\$\$P(y\|\lambda,\alpha,\sigma) = \int_0^\infty \frac{e^{-\lambda x}
\lambda^y x^y }{y!} \left(\frac{\alpha}{\sigma}\right)
\left(\frac{x}{\sigma}\right)^{\alpha-1}
e^{-\left(\frac{x}{\sigma}\right)^\alpha} dx\$\$ where \\f(x\| \alpha,
\sigma)\\ is the PDF of the Weibull distribution and \\\lambda\\ is the
mean of the Poisson distribution.

`dpoisweibull` computes the density of the Poisson-Weibull distribution.

`ppoisweibull` computes the distribution function of the Poisson-Weibull
distribution.

`qpoisweibull` computes the quantile function of the Poisson-Weibull
distribution.

`rpoisweibull` generates random numbers following the Poisson-Weibull
distribution.

The shape and scale parameters directly define the Weibull distribution,
whereas the mean and standard deviation are used to compute these
parameters indirectly.

## Examples

``` r
dpoisweibull(4, lambda=1.5, mean_value=1.5, sd_value=0.5, ndraws=10)
#> [1] 0.04381581
ppoisweibull(4, lambda=1.5, mean_value=1.5, sd_value=2, ndraws=10)
#> [1] 0.9693094
qpoisweibull(0.95, lambda=1.5, mean_value=1.5, sd_value=2, ndraws=10)
#> [1] 4
rpoisweibull(10, lambda=1.5, mean_value=1.5, sd_value=2, ndraws=10)
#>  [1] 0 1 0 2 2 3 6 0 0 1
```
