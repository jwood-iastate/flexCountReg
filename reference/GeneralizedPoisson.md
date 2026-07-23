# Generalized Poisson Version 1 Distribution

These functions provide the density function, distribution function,
quantile function, and random number generation for the Generalized
Poisson Version 1 (GP-1) Distribution

## Usage

``` r
.gp1_validate(mu, phi)

dgp1(x, mu = 1, phi = 1, log = FALSE)

pgp1(q, mu = 1, phi = 1, lower.tail = TRUE, log.p = FALSE)

qgp1(p, mu = 1, phi = 1, lower.tail = TRUE, log.p = FALSE)

rgp1(n, mu = 1, phi = 1)
```

## Arguments

- mu:

  numeric value or vector of mean values for the distribution (the
  values have to be greater than 0).

- phi:

  single value or vector of values for the scale parameter of the
  distribution (the values have to be greater than -1).

- x:

  numeric value or a vector of values.

- log:

  logical; if TRUE, probabilities p are given as log(p).

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

dgp1 gives the density, pgp1 gives the distribution function, qgp1 gives
the quantile function, and rgp1 generates random deviates.

The length of the result is determined by n for rgp1, and is the maximum
of the lengths of the numerical arguments for the other functions.

## Details

`dgp1` computes the density (PDF) of the Generalized Poisson
Distribution.

`pgp1` computes the CDF of the Generalized Poisson Distribution.

`qgp1` computes the quantile function of the Generalized Poisson
Distribution.

`rgp1` generates random numbers from the Generalized Poisson
Distribution.

The compound Probability Mass Function (PMF) for the Generalized Poisson
distribution, version 1, (GP-1) is: \$\$
f(y\|\phi,\mu)=\frac{\mu(\mu+\phi y)^{y-1} exp\left(-\frac{\mu+\phi y}
{1+\phi}\right)}{(1+\phi)^y y!} \$\$

Where \\\phi\\ is a scale parameter with the restriction that
\\\eta\>0\\, \\\mu\>0\\ is the mean value, and \\y\\ is a non-negative
integer. This formulation uses the mean directly.

The variance of the GP-1 distribution is: \$\$\sigma^2=(1+\phi)^2
\mu\$\$

If \\\phi\>0\\, the distribution is overdispersed. If \\\phi=0\\, the
distribution is equidispersed. If \\\phi\<0\\, the distribution is
underdispersed.

Furthermore, \\phi\>-1\\ is required for this distribution. When
\\\phi\<0\\, there is also a maximum value of support for the integar
value \\y\\. This is \\y_max=\left\lfloor -\frac{\mu}{\phi}
\right\rfloor\\.

## References

Consul, PoC, and Felix Famoye. "Generalized Poisson regression model."
Communications in Statistics-Theory and Methods 21.1 (1992): 89-109.

Zamani, Hossein, and Noriszura Ismail. "Functional form for the
generalized Poisson regression model." Communications in
Statistics-Theory and Methods 41.20 (2012): 3666-3675.

## Examples

``` r
dgp1(1, mu=0.75, phi=-0.1)
#> [1] 0.4047265
pgp1(c(0,1,2,3,5,7,9,10), mu=0.75, phi=3)
#> [1] 0.8290291 0.9024552 0.9317198 0.9479434 0.9657091 0.9753198 0.9813291
#> [8] 0.9835507
qgp1(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, phi=0.5)
#> [1] 0 0 0 2 3
rgp1(30, mu=0.75, phi=1.5)
#>  [1]  0  0  0  0  0  1  1 11  6  0  0  0  0  0  0  1  0  0  0  0  1  0  0  0  0
#> [26]  0  4  0  0  0
```
