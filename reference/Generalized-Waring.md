# Generalized Waring Distribution

These functions provide the probability mass function, distribution
function, quantile function, and random number generation for the
Generalized Waring distribution.

## Usage

``` r
dgwar(y, mu, k, rho, log = FALSE)

pgwar(q, mu, k, rho, lower.tail = TRUE, log.p = FALSE)

qgwar(p, mu, k, rho)

rgwar(n, mu, k, rho)
```

## Arguments

- y:

  non-negative integer vector of count outcomes.

- mu:

  positive numeric vector of distribution means.

- k:

  positive numeric vector of shape parameters.

- rho:

  numeric vector of shape parameters greater than 1. A finite variance
  requires \\\rho \> 2\\.

- log:

  logical; if TRUE, log probabilities are returned.

- q:

  non-negative integer vector of quantiles.

- lower.tail:

  logical; if TRUE, probabilities are \\P\[X \leq q\]\\; otherwise,
  \\P\[X \> q\]\\.

- log.p:

  logical; if TRUE, probabilities are returned as log probabilities.

- p:

  numeric vector of probabilities.

- n:

  integer number of random numbers to generate.

## Value

`dgwar` gives the probability mass function, `pgwar` gives the
distribution function, `qgwar` gives the quantile function, and `rgwar`
generates random deviates.

The length of the result is determined by `n` for `rgwar` and is the
maximum of the lengths of the numerical arguments for the other
functions.

## Details

The Generalized Waring distribution is a three-parameter count
distribution used to model overdispersed count data.

`dgwar` computes the probability mass function (PMF) of the Generalized
Waring distribution.

`pgwar` computes the CDF of the Generalized Waring distribution.

`qgwar` computes the quantile function of the Generalized Waring
distribution.

`rgwar` generates random numbers from the Generalized Waring
distribution.

The probability mass function for the Generalized Waring distribution
is: \$\$ f(y \mid a_x, k, \rho) = \frac{ \Gamma(a_x + \rho)\Gamma(k +
\rho) (a_x)\_y(k)\_y }{ y!\Gamma(\rho)\Gamma(a_x + k + \rho) (a_x + k +
\rho)\_y } \$\$ where \\(\alpha)\_r =
\frac{\Gamma(\alpha+r)}{\Gamma(\alpha)}\\, and \\a_x \> 0\\, \\k \> 0\\,
and \\\rho \> 0\\.

When \\\rho \> 1\\, the mean is: \$\$ E\[Y\] = \frac{a_x k}{\rho - 1}.
\$\$

Therefore, the distribution can be parameterized in terms of its mean
using: \$\$ a_x = \frac{\mu(\rho - 1)}{k}. \$\$

For a regression model: \$\$ \mu = \exp(X\beta). \$\$

When \\\rho \> 2\\, the variance under this mean parameterization is:
\$\$ \mathrm{Var}(Y) = \frac{\mu(k+\mu)(k+\rho-1)}{k(\rho-2)}. \$\$

## Examples

``` r
dgwar(0, mu = 1, k = 2, rho = 3)
#> [1] 0.6
pgwar(c(0, 1, 2, 3), mu = 1, k = 2, rho = 3)
#> [1] 0.6000000 0.8000000 0.8857143 0.9285714
qgwar(0.8, mu = 1, k = 2, rho = 3)
#> [1] 1
rgwar(10, mu = 1, k = 2, rho = 3)
#>  [1] 0 0 0 1 0 0 1 4 0 0
```
