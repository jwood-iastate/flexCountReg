# Poisson-Lindley-Lognormal Distribution

These functions provide density, distribution, quantile, and random
generation for the Poisson-Lindley-Lognormal (PLL) Distribution.

## Usage

``` r
dplindLnorm(
  x,
  mean = 1,
  theta = 1,
  sigma = 1,
  ndraws = 1500,
  log = FALSE,
  hdraws = NULL,
  engine = c("quadrature", "halton")
)

pplindLnorm(
  q,
  mean = 1,
  theta = 1,
  lambda = NULL,
  sigma = 1,
  ndraws = 1500,
  lower.tail = TRUE,
  log.p = FALSE,
  engine = c("quadrature", "halton")
)

qplindLnorm(
  p,
  mean = 1,
  theta = 1,
  sigma = 1,
  ndraws = 1500,
  lambda = NULL,
  engine = c("quadrature", "halton")
)

rplindLnorm(n, mean = 1, theta = 1, sigma = 1, ndraws = 1500, lambda = NULL)
```

## Arguments

- x:

  numeric value or vector of values.

- mean:

  the Poisson-Lindley mean before the lognormal factor (\>0); the mean
  of the distribution is `mean * exp(sigma^2 / 2)`.

- theta:

  Poisson-Lindley theta parameter (\>0).

- sigma:

  lognormal sigma parameter (\>= 0; 0 gives the Poisson-Lindley
  distribution).

- ndraws:

  the number of Halton draws (only with `engine = "halton"`).

- log:

  return log-density.

- hdraws:

  optional Halton draws (only with `engine = "halton"`).

- engine:

  the numerical method: `"quadrature"` (default) computes the
  log-probability by adaptive quadrature over the log of the mixing
  variable, accurate for every count including far tails; `"halton"`
  averages the Poisson probabilities over Halton draws (the method of
  flexCountReg 0.1.2 and earlier, kept for reproducibility; badly biased
  far in the tail; the default when `hdraws` are given).

- q:

  quantile or vector of quantiles.

- lambda:

  optional Poisson scale parameter (\>0), an alternative to `mean`:
  `mean = lambda * (theta + 2) / (theta * (theta + 1))`.

- lower.tail:

  TRUE returns \\P\[X \leq x\]\\.

- log.p:

  return log-CDF.

- p:

  probability or vector of probabilities.

- n:

  number of random draws.

## Value

dplindLnorm gives the density, pplindLnorm gives the distribution
function, qplindLnorm gives the quantile function, and rplindLnorm
generates random deviates.

The length of the result is determined by n for rplindLnorm, and is the
maximum of the lengths of the numerical arguments for the other
functions.

## Details

The PLL is a 3-parameter count distribution that captures high mass at
small y and allows flexible heavy tails.

`dplindLnorm` computes the PLL density. `pplindLnorm` computes the PLL
CDF. `qplindLnorm` computes quantiles. `rplindLnorm` generates random
draws.

The PMF is: \$\$ f(y\|\mu,\theta,\sigma)=\int_0^\infty
\frac{\theta^2(\lambda x)^y(\theta+\lambda x+y+1)}
{(\theta+1)(\theta+\lambda x)^{y+2}}
\frac{\exp\left(-\frac{\ln^2(x)}{2\sigma^2}\right)}
{x\sigma\sqrt{2\pi}}dx \$\$

Where \\\mu\\ is the `mean` argument (the Poisson-Lindley mean) and
\\\lambda=\mu\theta(\theta+1)/(\theta+2)\\ is the corresponding
Poisson-Lindley scale. The mean of the distribution is: \$\$ E\[y\]=\mu
e^{\sigma^2/2}=\frac{\lambda(\theta+2)e^{\sigma^2/2}} {\theta(\theta+1)}
\$\$

Halton draws are used to evaluate the integral.

## Examples

``` r
dplindLnorm(0, mean=0.75, theta=7, sigma=2, ndraws=10)
#> [1] 0.5424503
pplindLnorm(0:10, mean=0.75, theta=7, sigma=2, ndraws=10)
#>  [1] 0.5424503 0.6931577 0.7653293 0.8086522 0.8379407 0.8592379 0.8755060
#>  [8] 0.8883823 0.8988519 0.9075462 0.9148902
qplindLnorm(c(0.1,0.5,0.9), lambda=4.67, theta=7, sigma=2)
#> [1] 0 0 9
rplindLnorm(5, mean=0.75, theta=7, sigma=2)
#> [1] 0 0 1 4 0
```
