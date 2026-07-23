# Generalized Poisson Version 2 Distribution

These functions provide the density function, distribution function,
quantile function, and random number generation for the Generalized
Poisson Version 2 (GP-2) Distribution.

## Usage

``` r
dgp2(x, mu = 1, alpha = 0, log = FALSE)

pgp2(q, mu = 1, alpha = 0, lower.tail = TRUE, log.p = FALSE)

qgp2(p, mu = 1, alpha = 0, lower.tail = TRUE, log.p = FALSE)

rgp2(n, mu = 1, alpha = 0)
```

## Arguments

- x:

  numeric value or a vector of values.

- mu:

  numeric value or vector of mean values for the distribution (the
  values have to be greater than 0).

- alpha:

  numeric value or vector of values for the dispersion parameter of the
  distribution. Values may be negative, zero, or positive, provided that
  \\1 + \alpha \mu \> 0\\.

- log:

  logical; if TRUE, probabilities are given as log(p).

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

dgp2 gives the density, pgp2 gives the distribution function, qgp2 gives
the quantile function, and rgp2 generates random deviates.

The length of the result is determined by n for rgp2, and is the maximum
of the lengths of the numerical arguments for the other functions.

## Details

`dgp2` computes the density (PMF) of the Generalized Poisson Version 2
Distribution.

`pgp2` computes the CDF of the Generalized Poisson Version 2
Distribution.

`qgp2` computes the quantile function of the Generalized Poisson Version
2 Distribution.

`rgp2` generates random numbers from the Generalized Poisson Version 2
Distribution.

The probability mass function (PMF) for the Generalized Poisson Version
2 distribution (GP-2) is: \$\$ f(y\|\mu,\alpha)= \frac{\mu(\mu+\alpha\mu
y)^{y-1} \exp\left(-\frac{\mu+\alpha\mu y}{1+\alpha\mu}\right)}
{(1+\alpha\mu)^y y!} \$\$

where \\\mu\>0\\ is the mean and \\y\\ is a non-negative integer. This
formulation uses the mean directly.

The variance of the GP-2 distribution is:
\$\$\sigma^2=\mu(1+\alpha\mu)^2\$\$

If \\\alpha \> 0\\, the distribution is overdispersed. If \\\alpha =
0\\, the distribution is equidispersed and reduces to the ordinary
Poisson distribution. If \\\alpha \< 0\\ and \\1 + \alpha\mu \> 0\\, the
distribution is underdispersed and has finite support.

Under the underdispersed case, the largest possible value of \\y\\ is
bounded above by \\\left\lceil -1/\alpha \right\rceil - 1\\, i.e., the
largest integer satisfying \\1 + \alpha y \> 0\\.

## References

Consul, PoC, and Felix Famoye. "Generalized Poisson regression model."
Communications in Statistics-Theory and Methods 21.1 (1992): 89-109.

Wang, Weizhen, and Felix Famoye. "Modeling household fertility decisions
with generalized Poisson regression." Journal of Population Economics
10.3 (1997): 273-283.

Yang, Zhao, James W. Hardin, and Cheryl L. Addy. "A score test for
overdispersion in Poisson regression based on the generalized Poisson-2
model." Journal of Statistical Planning and Inference 139.4 (2009):
1514-1521.

Harris, Tammy, Zhao Yang, and James W. Hardin. "Modeling underdispersed
count data with generalized Poisson regression." Stata Journal 12.4
(2012): 736-747.

## Examples

``` r
dgp2(1, mu = 0.75, alpha = 0.3)
#> [1] 0.2762245
pgp2(c(0, 1, 2, 3, 5, 7, 9, 10), mu = 0.75, alpha = 0.25)
#> [1] 0.5317515 0.8185413 0.9345470 0.9771263 0.9973058 0.9996865 0.9999634
#> [8] 0.9999875
qgp2(c(0.1, 0.3, 0.5, 0.9, 0.95), mu = 0.75, alpha = -0.2)
#> [1] 0 0 1 2 2
rgp2(30, mu = 0.75, alpha = 0.2)
#>  [1] 0 0 0 0 0 0 0 1 2 0 1 1 0 0 0 2 0 0 1 0 5 0 0 2 1 0 1 2 0 0
```
