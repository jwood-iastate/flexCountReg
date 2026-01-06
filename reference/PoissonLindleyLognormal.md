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
  hdraws = NULL
)

pplindLnorm(
  q,
  mean = 1,
  theta = 1,
  lambda = NULL,
  sigma = 1,
  ndraws = 1500,
  lower.tail = TRUE,
  log.p = FALSE
)

qplindLnorm(p, mean = 1, theta = 1, sigma = 1, ndraws = 1500, lambda = NULL)

rplindLnorm(n, mean = 1, theta = 1, sigma = 1, ndraws = 1500, lambda = NULL)
```

## Arguments

- x:

  numeric value or vector of values.

- mean:

  mean (\>0).

- theta:

  Poisson-Lindley theta parameter (\>0).

- sigma:

  lognormal sigma parameter (\>0).

- ndraws:

  number of Halton draws.

- log:

  return log-density.

- hdraws:

  optional Halton draws.

- q:

  quantile or vector of quantiles.

- lambda:

  optional lambda parameter (\>0).

- lower.tail:

  TRUE returns P\[X \$\$\leq\$\$ x\].

- log.p:

  return log-CDF.

- p:

  probability or vector of probabilities.

- n:

  number of random draws.

## Details

The PLL is a 3-parameter count distribution that captures high mass at
small y and allows flexible heavy tails.

`dplind` computes the PLL density. `pplind` computes the PLL CDF.
`qplind` computes quantiles. `rplind` generates random draws.

The PMF is: \$\$ f(y\|\mu,\theta,\sigma)=\int_0^\infty
\frac{\theta^2\mu^y x^y(\theta+\mu x+y+1)} {(\theta+1)(\theta+\mu
x)^{y+2}} \frac{\exp\left(-\frac{\ln^2(x)}{2\sigma^2}\right)}
{x\sigma\sqrt{2\pi}}dx \$\$

Mean: \$\$ E\[y\]=\mu=\frac{\lambda(\theta+2)e^{\sigma^2/2}}
{\theta(\theta+1)} \$\$

Halton draws are used to evaluate the integral.

dplindLnorm gives the density, pplindLnorm gives the distribution
function, qplindLnorm gives the quantile function, and rplindLnorm
generates random deviates.

The length of the result is determined by n for rplindLnorm, and is the
maximum of the lengths of the numerical arguments for the other
functions.

## Examples

``` r
dplindLnorm(0, mean=0.75, theta=7, sigma=2, ndraws=10)
#> [1] 0.6072274
pplindLnorm(0:10, mean=0.75, theta=7, sigma=2, ndraws=10)
#>  [1] 0.6072274 0.7717043 0.8450268 0.8859278 0.9117652 0.9294614 0.9422986
#>  [8] 0.9520167 0.9596144 0.9657016 0.9706706
qplindLnorm(c(0.1,0.5,0.9), lambda=4.67, theta=7, sigma=2)
#> [1]  0  3 61
rplindLnorm(5, mean=0.75, theta=7, sigma=2)
#> [1] 6 0 0 0 0
```
