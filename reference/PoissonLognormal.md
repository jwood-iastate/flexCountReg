# Poisson-Lognormal Distribution

These functions provide density, distribution function, quantile
function, and random number generation for the Poisson-Lognormal (PLogN)
Distribution

## Usage

``` r
dpLnorm(x, mean = 1, sigma = 1, ndraws = 1500, log = FALSE, hdraws = NULL)

ppLnorm(
  q,
  mean = 1,
  sigma = 1,
  ndraws = 1500,
  lower.tail = TRUE,
  log.p = FALSE
)

qpLnorm(p, mean = 1, sigma = 1, ndraws = 1500)

rpLnorm(n, mean = 1, sigma = 1, ndraws = 1500)
```

## Arguments

- x:

  numeric value or a vector of values.

- mean:

  numeric value or vector of mean values for the distribution (the
  values have to be greater than 0).

- sigma:

  single value or vector of values for the sigma parameter of the
  lognormal distribution (the values have to be greater than 0).

- ndraws:

  the number of Halton draws to use for the integration.

- log:

  logical; if TRUE, probabilities p are given as log(p).

- hdraws:

  and optional vector of Halton draws to use for the integration.

- q:

  quantile or a vector of quantiles.

- lower.tail:

  logical; if TRUE, probabilities p are \\P\[X\leq x\]\\ otherwise,
  \\P\[X\>x\]\\.

- log.p:

  logical; if TRUE, probabilities p are given as log(p).

- p:

  probability or a vector of probabilities.

- n:

  the number of random numbers to generate.

## Value

dpLnorm gives the density, ppLnorm gives the distribution function,
qpLnorm gives the quantile function, and rpLnorm generates random
deviates.

The length of the result is determined by n for rpLnorm, and is the
maximum of the lengths of the numerical arguments for the other
functions.

## Details

`dpLnorm` computes the density (PDF) of the Poisson-Lognormal
Distribution.

`ppLnorm` computes the CDF of the Poisson-Lognormal Distribution.

`qpLnorm` computes the quantile function of the Poisson-Lognormal
Distribution.

`rpLnorm` generates random numbers from the Poisson-Lognormal
Distribution.

The compound Probability Mass Function (PMF) for the Poisson-Lognormal
distribution is: \$\$f(y\|\mu,\theta,\alpha)= \int_0^\infty \frac{\mu^y
x^y e^{-\mu x}}{y!} \frac{exp\left(-\frac{ln^2(x)}{2\sigma^2}
\right)}{x\sigma\sqrt{2\pi}}dx\$\$

Where \\\sigma\\ is a parameter for the lognormal distribution with the
restriction \\\sigma\>0\\, and \\y\\ is a non-negative integer.

The expected value of the distribution is:
\$\$E\[y\]=e^{X\beta+\sigma^2/2} = \mu e^{\sigma^2/2}\$\$ Halton draws
are used to perform simulation over the lognormal distribution to solve
the integral.

## Examples

``` r
dpLnorm(0, mean=0.75, sigma=2, ndraws=10)
#> [1] 0.5271804
ppLnorm(c(0,1,2,3,5,7,9,10), mean=0.75, sigma=2, ndraws=10)
#> [1] 0.5271804 0.7245134 0.8159930 0.8646682 0.9165510 0.9516880 0.9777184
#> [8] 0.9863234
qpLnorm(c(0.1,0.3,0.5,0.9,0.95), mean=0.75, sigma=2, ndraws=10)
#> [1] 0 0 0 5 7
rpLnorm(30, mean=0.75,  sigma=2, ndraws=10)
#>  [1] 1 0 0 3 0 2 0 6 2 2 1 0 1 5 6 0 0 1 0 4 0 1 0 0 0 0 0 7 2 1
```
