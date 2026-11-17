# Polya-Aeppli Distribution

Density, distribution function, quantile function, and random generation
for the Polya-Aeppli (Poisson-geometric) distribution.

## Usage

``` r
dpaeppli(x, lambda, prob, log = FALSE)

ppaeppli(q, lambda, prob, lower.tail = TRUE, log.p = FALSE)

qpaeppli(p, lambda, prob, lower.tail = TRUE, log.p = FALSE)

rpaeppli(n, lambda, prob)
```

## Arguments

- x:

  vector of (non-negative integer) quantiles.

- lambda:

  vector of (non-negative) Poisson rates \\\lambda\\ of clusters.

- prob:

  vector of geometric parameters \\0 \leq q \< 1\\: each cluster has \\c
  \geq 1\\ events with probability \\(1 - q) q^{c - 1}\\.

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

`dpaeppli` gives the probability mass, `ppaeppli` the distribution
function, `qpaeppli` the quantile function, and `rpaeppli` random
deviates. The length of the result is determined by `n` for `rpaeppli`,
and is the maximum of the lengths of the numerical arguments for the
other functions.

## Details

A Polya-Aeppli random variable is a Poisson number of clusters with
geometric cluster sizes: \$\$Y = \sum\_{i = 1}^{M} C_i, \quad M \sim
\mathrm{Poisson}(\lambda), \quad P(C_i = c) = (1 - q) q^{c - 1}, \\ c =
1, 2, \ldots\$\$ so that \\P(Y = 0) = e^{-\lambda}\\ and, for \\y \geq
1\\, \$\$P(Y = y) = e^{-\lambda} \sum\_{j = 1}^{y} \frac{\lambda^j}{j!}
\binom{y - 1}{j - 1} q^{y - j} (1 - q)^j.\$\$ The mean and variance are
\$\$E\[Y\] = \frac{\lambda}{1 - q}, \qquad \mathrm{Var}\[Y\] =
\frac{\lambda (1 + q)}{(1 - q)^2}.\$\$ \\q = 0\\ gives the Poisson
distribution with mean \\\lambda\\. The parameterization is that of the
polyaAeppli package.

The Polya-Aeppli distribution is the Poisson-Tweedie distribution with
power parameter 1.5: a count with mean \\\mu\\ and variance \\\mu +
\phi\mu^{1.5}\\ is Polya-Aeppli with \\q = d / (2 + d)\\ and \\\lambda =
\mu (1 - q)\\, where \\d = \phi\sqrt{\mu}\\.

Probabilities are computed by summing the decomposition over the number
of clusters on the log scale (all terms are positive); both tails are
computed directly. Random deviates are generated from the definition.

## References

Johnson, N. L., Kemp, A. W. and Kotz, S. (2005). *Univariate Discrete
Distributions* (3rd ed.), Section 9.7. Wiley.

Bonat, W. H., Jorgensen, B., Kokonendji, C. C., Hinde, J. and Demetrio,
C. G. B. (2018). Extended Poisson-Tweedie: properties and regression
models for count data. *Statistical Modelling*, 18(1), 24-49.

## Examples

``` r
dpaeppli(0:6, lambda = 2, prob = 0.3)
#> [1] 0.13533528 0.18946940 0.18946940 0.15852273 0.11829206 0.08134174 0.05255308
# prob = 0 gives the Poisson distribution
all.equal(dpaeppli(0:10, lambda = 2, prob = 0), dpois(0:10, 2))
#> [1] TRUE
ppaeppli(4, lambda = 2, prob = 0.3)
#> [1] 0.7910889
qpaeppli(c(0.5, 0.95), lambda = 2, prob = 0.3)
#> [1] 2 7
set.seed(7)
y <- rpaeppli(5000, lambda = 2, prob = 0.3)
c(mean = mean(y), var = var(y))   # theory: 2.857 and 5.306
#>     mean      var 
#> 2.860200 5.362128 
```
