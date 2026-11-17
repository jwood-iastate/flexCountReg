# Hyper-Poisson Distribution

Density, distribution function, quantile function, and random generation
for the hyper-Poisson distribution of Bardwell and Crow (1964).

## Usage

``` r
dhpois(x, lambda, gamma, log = FALSE)

phpois(q, lambda, gamma, lower.tail = TRUE, log.p = FALSE)

qhpois(p, lambda, gamma, lower.tail = TRUE, log.p = FALSE)

rhpois(n, lambda, gamma)
```

## Arguments

- x:

  vector of (non-negative integer) quantiles.

- lambda:

  vector of (positive) rate parameters \\\lambda\\.

- gamma:

  vector of (positive) dispersion parameters \\\gamma\\. \\\gamma = 1\\
  is the Poisson distribution, \\\gamma \> 1\\ gives overdispersion and
  \\\gamma \< 1\\ underdispersion.

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

`dhpois` gives the probability mass, `phpois` the distribution function,
`qhpois` the quantile function, and `rhpois` random deviates. The length
of the result is determined by `n` for `rhpois`, and is the maximum of
the lengths of the numerical arguments for the other functions.

## Details

The hyper-Poisson distribution has probabilities \$\$P(Y = y) =
\frac{\lambda^y}{(\gamma)\_y \\ {}\_1F_1(1; \gamma; \lambda)}, \quad y =
0, 1, 2, \ldots\$\$ where \\(\gamma)\_y = \Gamma(\gamma + y) /
\Gamma(\gamma)\\ is the rising factorial and \\{}\_1F_1(1; \gamma;
\lambda) = \sum\_{j \geq 0} \lambda^j / (\gamma)\_j\\ the confluent
hypergeometric function, the normalizing constant. Successive
probabilities satisfy \\P(Y = y + 1) / P(Y = y) = \lambda / (\gamma +
y)\\.

\\\lambda\\ is not the mean. The mean is \$\$E\[Y\] = \lambda -
(\gamma - 1)(1 - P(Y = 0)),\$\$ and moments are computed from the
normalized probabilities when needed. Regression models usually
reparameterize the distribution by its mean (Saez-Castillo and
Conde-Sanchez, 2013); symbols in other packages can differ from those
used here.

The normalizing constant is computed by summing the series around its
mode until the remaining tail (bounded by a geometric series, because
the ratio of successive terms decreases) is negligible, so the
computations do not overflow for large \\\lambda\\. Upper tail
probabilities are summed directly. `rhpois()` uses inversion of the
distribution function.

## References

Bardwell, G. E. and Crow, E. L. (1964). A two-parameter family of
hyper-Poisson distributions. *Journal of the American Statistical
Association*, 59(305), 133-141.

Saez-Castillo, A. J. and Conde-Sanchez, A. (2013). A hyper-Poisson
regression model for overdispersed and underdispersed count data.
*Computational Statistics & Data Analysis*, 61, 148-157.

## Examples

``` r
dhpois(0:6, lambda = 3, gamma = 2)
#> [1] 0.15718709 0.23578063 0.23578063 0.17683548 0.10610129 0.05305064 0.02273599
# gamma = 1 gives the Poisson distribution
all.equal(dhpois(0:10, lambda = 3, gamma = 1), dpois(0:10, 3))
#> [1] TRUE
phpois(4, lambda = 3, gamma = 0.5)
#> [1] 0.7404145
qhpois(c(0.1, 0.5, 0.9), lambda = 3, gamma = 2)
#> [1] 0 2 4
set.seed(5)
y <- rhpois(2000, lambda = 4, gamma = 0.3)
c(mean = mean(y), var = var(y))   # underdispersed
#>     mean      var 
#> 4.714500 4.059019 
```
