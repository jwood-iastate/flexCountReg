# Delaporte Distribution

Density, distribution function, quantile function, and random generation
for the Delaporte distribution in a mean parameterization.

## Usage

``` r
ddelaporte(x, mu, rho, k, log = FALSE)

pdelaporte(q, mu, rho, k, lower.tail = TRUE, log.p = FALSE)

qdelaporte(p, mu, rho, k, lower.tail = TRUE, log.p = FALSE)

rdelaporte(n, mu, rho, k)
```

## Arguments

- x:

  vector of (non-negative integer) quantiles.

- mu:

  vector of means \\\mu \geq 0\\.

- rho:

  vector of shares \\0 \leq \rho \leq 1\\ of the mean that come from the
  negative binomial (gamma-mixed) component.

- k:

  vector of (positive) negative binomial size (shape) parameters.

- log, log.p:

  logical; if TRUE, probabilities p are given as log(p).

- q:

  vector of quantiles.

- lower.tail:

  logical; if TRUE (default), probabilities are \\P(Y \leq y)\\,
  otherwise \\P(Y \> y)\\.

- p:

  vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

## Value

`ddelaporte` gives the probability mass, `pdelaporte` the distribution
function, `qdelaporte` the quantile function, and `rdelaporte` random
deviates. The length of the result is determined by `n` for
`rdelaporte`, and is the maximum of the lengths of the numerical
arguments for the other functions.

## Details

A Delaporte random variable is the sum of two independent counts, \$\$Y
= A + B, \quad A \sim \mathrm{Poisson}((1 - \rho)\mu), \quad B \sim
\mathrm{NB}(\mathrm{size} = k, \mathrm{mean} = \rho\mu),\$\$ so that
\$\$E\[Y\] = \mu, \qquad \mathrm{Var}\[Y\] = \mu +
\frac{(\rho\mu)^2}{k}.\$\$ Equivalently it is a Poisson distribution
whose rate is a shifted gamma variable. \\\rho = 0\\ gives the Poisson
distribution with mean \\\mu\\ and \\\rho = 1\\ the negative binomial
distribution with mean \\\mu\\ and size \\k\\ (NB2 with \\\alpha =
1/k\\). Both endpoints are accepted here.

In the \\(\alpha, \beta, \lambda)\\ parameterization of the Delaporte
package, \\\alpha = k\\, \\\beta = \rho\mu / k\\ and \\\lambda = (1 -
\rho)\mu\\.

Probabilities are computed by summing the convolution of the two
components on the log scale (all terms are positive). Both tails are
computed directly from the same decomposition, so `lower.tail = FALSE`
is accurate far into the upper tail. Random deviates are generated as
the sum of independent Poisson and negative binomial draws. The work per
probability grows linearly with the count.

## References

Ruohonen, M. (1988). On a model for the claim number process. *ASTIN
Bulletin*, 18(1), 57-68.

Vose, D. (2008). *Risk Analysis: A Quantitative Guide* (3rd ed.). Wiley.

## Examples

``` r
ddelaporte(0:6, mu = 2, rho = 0.5, k = 1.5)
#> [1] 0.17097491 0.27355986 0.23936488 0.15501726 0.08457559 0.04172700 0.01939783
# rho = 0 is Poisson, rho = 1 is negative binomial
all.equal(ddelaporte(0:10, mu = 2, rho = 0, k = 1.5), dpois(0:10, 2))
#> [1] TRUE
all.equal(ddelaporte(0:10, mu = 2, rho = 1, k = 1.5),
          dnbinom(0:10, size = 1.5, mu = 2))
#> [1] TRUE
pdelaporte(5, mu = 2, rho = 0.5, k = 1.5)
#> [1] 0.9652195
qdelaporte(c(0.5, 0.95), mu = 2, rho = 0.5, k = 1.5)
#> [1] 2 5
set.seed(3)
y <- rdelaporte(5000, mu = 4, rho = 0.6, k = 2)
c(mean = mean(y), var = var(y))   # theory: 4 and 4 + 2.4^2 / 2 = 6.88
#>     mean      var 
#> 4.011400 6.978266 
```
