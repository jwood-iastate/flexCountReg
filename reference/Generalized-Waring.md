# Generalized Waring Distribution

These functions provide density, distribution function, quantile
function, and random number generation for the Generalized Waring
Distribution.

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

  numeric vector of means of the distribution.

- k:

  non-negative numeric parameter of the distribution.

- rho:

  non-negative numeric parameter of the distribution.

- log:

  logical; if TRUE, probabilities p are given as log(p).

- q:

  non-negative integer vector of quantiles.

- lower.tail:

  logical; if TRUE, probabilities p are \\P\[X\leq x\]\\ otherwise,
  \\P\[X\>x\]\\.

- log.p:

  logical; if TRUE, probabilities p are given as log(p).

- p:

  numeric vector of probabilities.

- n:

  integer number of random numbers to generate.

## Details

The Generalized Waring distribution is a 3-parameter count distribution
that is used to model overdispersed count data.

`dgwar` computes the density (PMF) of the Generalized Waring
Distribution.

`pgwar` computes the CDF of the Generalized Waring Distribution.

`qwaring` computes the quantile function of the Generalized Waring
Distribution.

`rwaring` generates random numbers from the Generalized Waring
Distribution.

The Probability Mass Function (PMF) for the Generalized Waring (GW)
distribution is: \$\$f(y\|a_x,k,\rho) =
\frac{\Gamma(a_x+\rho)\Gamma(k+\rho)\left(a_x\right)\_y(k)\_y}
{y!\Gamma(\rho)\Gamma(a_x+k+\rho)(a_x+k+\rho)\_y}\$\$ Where
\\(\alpha)\_r=\frac{\Gamma(\alpha+r)}{\Gamma(\alpha)}\\, and \\a_x, \\
k, \\ \rho)\>0\\.

The mean value is: \$\$E\[Y\]=\frac{a_x K}{\rho-1}\$\$

Thus, we can use: \$\$a_x=\frac{\mu(\rho-1)}{k}\$\$

This results in a regression model where: \$\$\mu=e^{X\beta}\$\$
\$\$\sigma^2 = \mu \left(1-\frac{1}{\alpha+\rho+1} \right) +
\mu^2\frac{(\alpha+\rho)^2}{\alpha\rho(\alpha+\rho+1)}\$\$

dgwar gives the density, pgwar gives the distribution function, qgwar
gives the quantile function, and rgwar generates random deviates.

The length of the result is determined by n for rgwar, and is the
maximum of the lengths of the numerical arguments for the other
functions.

## Examples

``` r
dgwar(0, mu=1, k=2, rho=3)
#> [1] 0.6
pgwar(c(0,1,2,3), mu=1, k=2, rho=3)
#> [1] 0.6000000 0.8000000 0.8857143 0.9285714
qgwar(0.8, mu=1, k=2, rho=3)
#> [1] 1
rgwar(10, mu=1, k=2, rho=3)
#>  [1] 0 0 1 4 0 0 1 0 1 1
```
