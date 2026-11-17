# Negative Binomial-Generalized Exponential Distribution

Density, distribution function, quantile function, and random generation
for the negative binomial-generalized exponential (NBGE) mixture: a
negative binomial (NB2) count whose mean is multiplied by a generalized
exponential random effect with mean one.

## Usage

``` r
dnbge(x, mu, shape, alpha, log = FALSE)

pnbge(q, mu, shape, alpha, lower.tail = TRUE, log.p = FALSE)

qnbge(p, mu, shape, alpha, lower.tail = TRUE, log.p = FALSE)

rnbge(n, mu, shape, alpha)
```

## Arguments

- x:

  vector of (non-negative integer) quantiles.

- mu:

  vector of (positive) means.

- shape:

  vector of (positive) shape parameters \\k\\ of the generalized
  exponential mixing distribution.

- alpha:

  vector of (non-negative) NB2 dispersion parameters; `alpha = 0` gives
  the Poisson-generalized exponential distribution
  ([`dpge`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonGeneralizedExponential.md)).

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

`dnbge` gives the probability mass, `pnbge` the distribution function,
`qnbge` the quantile function, and `rnbge` random deviates. The length
of the result is determined by `n` for `rnbge`, and is the maximum of
the lengths of the numerical arguments for the other functions.

## Details

The distribution is defined hierarchically: \$\$X \sim \mathrm{GE}(k,
1), \quad F_X(x) = (1 - e^{-x})^k, \quad G = X / A, \quad A =
\psi(k + 1) - \psi(1),\$\$ \$\$Y \mid G \sim \mathrm{NB2}(\mathrm{mean}
= \mu G, \mathrm{dispersion} = \alpha),\$\$ where \\\psi\\ is the
digamma function, so that \\E\[G\] = 1\\ and the conditional variance is
\\\mu G + \alpha (\mu G)^2\\. With \\B = \psi'(1) - \psi'(k + 1)\\ and
\\M_2 = 1 + B / A^2\\ (\\= E\[G^2\]\\), \$\$E\[Y\] = \mu, \qquad
\mathrm{Var}\[Y\] = \mu + \left((1 + \alpha) M_2 - 1\right)\mu^2.\$\$
The scale of the generalized exponential distribution cancels in this
mean parameterization and is not a parameter. \\\alpha = 0\\ gives the
Poisson-generalized exponential distribution and, as the
generalized-exponential heterogeneity vanishes, the distribution
approaches NB2. The two dispersion mechanisms have confounded effects on
the variance, so \\k\\ and \\\alpha\\ are only separately identified
through the shape of the distribution; expect weak identification in
small samples.

This construction is specified here; a distribution with a similar name
elsewhere may be defined differently.

Probabilities are computed on the log scale by adaptive quadrature over
the logarithm of the mixing variable (the same deterministic,
mode-centred rule as
[`dpge`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonGeneralizedExponential.md)),
so tail probabilities keep their relative accuracy. Upper tails are
summed directly from the probabilities. Random deviates are generated
from the definition.

## References

Gupta, R. D. and Kundu, D. (1999). Generalized exponential
distributions. *Australian & New Zealand Journal of Statistics*, 41(2),
173-188.

Greene, W. (2008). Functional forms for the negative binomial model for
count data. *Economics Letters*, 99(3), 585-590.

## Examples

``` r
dnbge(0:6, mu = 2, shape = 1.5, alpha = 0.5)
#> [1] 0.36375462 0.22612339 0.13793589 0.08665409 0.05615242 0.03739404 0.02549160
# alpha = 0 gives the Poisson-generalized exponential distribution
all.equal(dnbge(0:10, mu = 2, shape = 1.5, alpha = 0),
          dpge(0:10, mean = 2, shape = 1.5, scale = 1))
#> [1] TRUE
pnbge(5, mu = 2, shape = 1.5, alpha = 0.5)
#> [1] 0.9080144
qnbge(c(0.5, 0.95), mu = 2, shape = 1.5, alpha = 0.5)
#> [1] 1 7
set.seed(8)
y <- rnbge(5000, mu = 2, shape = 1.5, alpha = 0.5)
c(mean = mean(y), var = var(y))
#>    mean     var 
#> 2.08040 9.06975 
```
