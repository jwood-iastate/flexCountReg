# Discrete Weibull Distribution

Density, distribution function, quantile function, and random generation
for the (type I) discrete Weibull distribution of Nakagawa and Osaki
(1975).

## Usage

``` r
ddweib(x, theta, shape, log = FALSE)

pdweib(q, theta, shape, lower.tail = TRUE, log.p = FALSE)

qdweib(p, theta, shape, lower.tail = TRUE, log.p = FALSE)

rdweib(n, theta, shape)
```

## Arguments

- x:

  vector of (non-negative integer) quantiles.

- theta:

  vector of values of the parameter \\\theta\\, with \\0 \< \theta \<
  1\\. It equals \\P(Y \geq 1)\\.

- shape:

  vector of values of the shape parameter \\\beta \> 0\\.

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

`ddweib` gives the probability mass, `pdweib` the distribution function,
`qdweib` the quantile function, and `rdweib` random deviates. The length
of the result is determined by `n` for `rdweib`, and is the maximum of
the lengths of the numerical arguments for the other functions (zero if
any is zero).

## Details

The discrete Weibull distribution has survival function \$\$P(Y \geq y)
= \theta^{y^\beta}, \quad y = 0, 1, 2, \ldots\$\$ so that \$\$P(Y = y) =
\theta^{y^\beta} - \theta^{(y+1)^\beta}, \qquad P(Y \leq y) = 1 -
\theta^{(y+1)^\beta}.\$\$ It is the distribution of \\\lfloor X
\rfloor\\ where \\X\\ has a continuous Weibull distribution with shape
\\\beta\\ and survival function \\\exp(-c x^\beta)\\, \\c =
-\log\theta\\; `rdweib()` uses that representation. With \\\beta = 1\\
it is the geometric distribution with success probability \\1 - \theta\\
(`dgeom(x, 1 - theta)`).

The moments have no closed form: \$\$E\[Y\] = \sum\_{k \geq 1}
\theta^{k^\beta}, \qquad E\[Y^2\] = \sum\_{k \geq 1} (2k -
1)\theta^{k^\beta}.\$\$ The distribution can be under- or overdispersed
relative to the Poisson, which makes it a flexible model for counts
(Klakattawi et al., 2018).

This is a discrete response distribution; it is unrelated to the
Poisson-Weibull mixture of
[`dpoisweibull`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonWeibull.md).

All computations use the logarithms of the tail probabilities, so the
upper tail `lower.tail = FALSE` keeps its relative accuracy far into the
tail. `qdweib()` inverts the survival function exactly and returns the
smallest \\y\\ with \\P(Y \leq y) \geq p\\.

Invalid parameter values give `NaN` with a warning, and non-integer `x`
has probability 0 (with a warning), as in base R.

## References

Nakagawa, T. and Osaki, S. (1975). The discrete Weibull distribution.
*IEEE Transactions on Reliability*, R-24(5), 300-301.

Klakattawi, H., Vinciotti, V. and Yu, K. (2018). A simple and adaptive
dispersion regression model for count data. *Entropy*, 20(2), 142.

## Examples

``` r
ddweib(0:5, theta = 0.7, shape = 1.5)
#> [1] 0.30000000 0.33535430 0.20793242 0.09906527 0.03910657 0.01325193
# shape = 1 gives the geometric distribution
all.equal(ddweib(0:10, theta = 0.6, shape = 1), dgeom(0:10, prob = 0.4))
#> [1] TRUE
pdweib(3, theta = 0.7, shape = 1.5)
#> [1] 0.942352
pdweib(60, theta = 0.7, shape = 1.5, lower.tail = FALSE, log.p = TRUE)
#> [1] -169.9289
qdweib(c(0.1, 0.5, 0.9), theta = 0.7, shape = 1.5)
#> [1] 0 1 3
set.seed(1)
y <- rdweib(1000, theta = 0.9, shape = 0.8)
c(mean = mean(y), var = var(y))   # overdispersed
#>     mean      var 
#>  18.9820 549.8495 
```
