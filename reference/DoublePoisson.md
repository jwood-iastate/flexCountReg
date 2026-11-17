# Double Poisson Distribution

Density, distribution function, quantile function, and random generation
for the double Poisson distribution of Efron (1986), normalized exactly.

## Usage

``` r
ddpo(x, lambda, phi, log = FALSE)

pdpo(q, lambda, phi, lower.tail = TRUE, log.p = FALSE)

qdpo(p, lambda, phi, lower.tail = TRUE, log.p = FALSE)

rdpo(n, lambda, phi)
```

## Arguments

- x:

  vector of (non-negative integer) quantiles.

- lambda:

  vector of (positive) location parameters \\\lambda\\.

- phi:

  vector of (positive) precision parameters \\\phi\\. \\\phi = 1\\ is
  the Poisson distribution, \\\phi \< 1\\ gives overdispersion and
  \\\phi \> 1\\ underdispersion.

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

`ddpo` gives the probability mass, `pdpo` the distribution function,
`qdpo` the quantile function, and `rdpo` random deviates. The length of
the result is determined by `n` for `rdpo`, and is the maximum of the
lengths of the numerical arguments for the other functions.

## Details

The double Poisson distribution has probabilities proportional to the
weights (with \\0 \log 0 = 0\\) \$\$w(y) = \exp\left\\\phi y
\log\lambda + (1 - \phi) y \log y + (\phi - 1) y - \log y!\right\\,
\quad y = 0, 1, 2, \ldots\$\$ which is Efron's (1986) expression with
the factors that do not depend on \\y\\ removed. Here the normalizing
constant \\\sum\_{j \geq 0} w(j)\\ is computed numerically (by summing
the weights around their mode until the remaining tail, bounded by a
geometric series, is negligible), so the probabilities sum to one.
Efron's closed-form approximation of the constant is not used: it is
least accurate for small \\\lambda\\ and strong dispersion.

\\\lambda\\ is a location parameter: the mean is close to, but not
exactly, \\\lambda\\, and the variance is approximately \\\lambda /
\phi\\. The exact moments are computed from the normalized probabilities
when needed. Some references (e.g. the `DPO` family of gamlss.dist) use
the dispersion \\\sigma = 1 / \phi\\ instead of the precision.

Upper tail probabilities (`lower.tail = FALSE`) are summed directly,
never computed as one minus the distribution function. `rdpo()` uses
inversion of the distribution function.

## References

Efron, B. (1986). Double exponential families and their use in
generalized linear regression. *Journal of the American Statistical
Association*, 81(395), 709-721.

## Examples

``` r
ddpo(0:6, lambda = 3, phi = 0.5)
#> [1] 0.15271589 0.16043444 0.16854311 0.15333984 0.12400709 0.09103218 0.06158654
# phi = 1 gives the Poisson distribution
all.equal(ddpo(0:10, lambda = 3, phi = 1), dpois(0:10, 3))
#> [1] TRUE
sum(ddpo(0:200, lambda = 3, phi = 2.5))   # exactly normalized
#> [1] 1
pdpo(c(2, 5), lambda = 3, phi = 2.5)
#> [1] 0.3353030 0.9835871
qdpo(c(0.1, 0.5, 0.9), lambda = 3, phi = 0.5)
#> [1] 0 3 6
set.seed(2)
y <- rdpo(2000, lambda = 5, phi = 3)
c(mean = mean(y), var = var(y))   # underdispersed
#>     mean      var 
#> 5.016000 1.775632 
```
