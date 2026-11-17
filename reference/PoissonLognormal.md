# Poisson-Lognormal Distribution

These functions provide density, distribution function, quantile
function, and random number generation for the Poisson-Lognormal (PLogN)
Distribution

## Usage

``` r
dpLnorm(
  x,
  mean = 1,
  sigma = 1,
  ndraws = 1500,
  log = FALSE,
  hdraws = NULL,
  engine = c("quadrature", "halton", "poilog")
)

ppLnorm(
  q,
  mean = 1,
  sigma = 1,
  ndraws = 1500,
  lower.tail = TRUE,
  log.p = FALSE,
  engine = c("quadrature", "halton", "poilog")
)

qpLnorm(
  p,
  mean = 1,
  sigma = 1,
  ndraws = 1500,
  engine = c("quadrature", "halton", "poilog")
)

rpLnorm(n, mean = 1, sigma = 1)
```

## Arguments

- x:

  numeric value or a vector of values.

- mean:

  numeric value or vector of mean values (not the expected value). This
  is the mean for the Poisson portion of the distribution. Note that the
  mean values have to be greater than 0.

- sigma:

  single value or vector of values for the sigma parameter of the
  lognormal distribution (the values have to be greater than 0).

- ndraws:

  the number of Halton draws (only with `engine = "halton"`).

- log:

  logical; if TRUE, probabilities p are given as log(p).

- hdraws:

  and optional vector of Halton draws to use for the integration.

- engine:

  the numerical method: `"quadrature"` (default) computes the
  log-probability by adaptive quadrature over the log of the mixing
  variable, accurate for every count including far tails; `"halton"`
  averages the Poisson probabilities over Halton draws (the method of
  flexCountReg 0.1.2 and earlier, kept for reproducibility; badly biased
  far in the tail; the default when `hdraws` are given), or `"poilog"`
  (the poilog package).

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
distribution is: \$\$f(y\|\mu,\sigma)= \int_0^\infty \frac{\mu^y x^y
e^{-\mu x}}{y!} \frac{exp\left(-\frac{ln^2(x)}{2\sigma^2}
\right)}{x\sigma\sqrt{2\pi}}dx\$\$

Where \\\sigma\\ is a parameter for the lognormal distribution with the
restriction \\\sigma\>0\\, and \\y\\ is a non-negative integer.

The expected value of the distribution is:
\$\$E\[y\]=e^{X\beta+\sigma^2/2} = \mu e^{\sigma^2/2}\$\$

The variance for the distribution is:
\$\$V\[Y\]=E\[Y\]+\left(e^{\sigma^2}-1\right)E\[Y\]^2\$\$

Halton draws are used to perform simulation over the lognormal
distribution to solve the integral if the engine is set to "halton". The
poilog package is used to solve the integral if the engine is set to
"poilog".

## Examples

``` r
dpLnorm(0, mean=0.75, sigma=2, ndraws=10)
#> [1] 0.4609574
if (requireNamespace("poilog", quietly = TRUE)) {
  dpLnorm(0, mean=0.75, sigma=2, engine="poilog")
}
#> [1] 0.4609574
ppLnorm(c(0,1,2,3,5,7,9,10), mean=0.75, sigma=2, ndraws=10)
#> [1] 0.4609574 0.6316035 0.7189017 0.7723869 0.8352530 0.8713687 0.8949394
#> [8] 0.9039127
qpLnorm(c(0.1,0.3,0.5,0.9,0.95), mean=0.75, sigma=2, ndraws=10)
#> [1]  0  0  1 10 21
rpLnorm(30, mean=0.75,  sigma=2)
#>  [1]  7  0  7  0  0  2  4  1  2  0  2  8  0  1 31  2 10  0 76  1 13  2  7  1 13
#> [26]  5  7  3  0  8
```
