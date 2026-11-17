# Generalized Hermite Distribution

Density, distribution function, quantile function, and random generation
for the (generalized) Hermite distribution of order \\m\\, in a mean
parameterization.

## Usage

``` r
dghermite(x, mu, rho, m = 2, log = FALSE)

pghermite(q, mu, rho, m = 2, lower.tail = TRUE, log.p = FALSE)

qghermite(p, mu, rho, m = 2, lower.tail = TRUE, log.p = FALSE)

rghermite(n, mu, rho, m = 2)
```

## Arguments

- x:

  vector of (non-negative integer) quantiles.

- mu:

  vector of means \\\mu \geq 0\\.

- rho:

  vector of shares \\0 \leq \rho \leq 1\\ of the mean that come from
  events occurring in batches of size \\m\\.

- m:

  order: the (fixed, whole) batch size, \\m \geq 1\\. The default
  `m = 2` is the classical Hermite distribution.

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

`dghermite` gives the probability mass, `pghermite` the distribution
function, `qghermite` the quantile function, and `rghermite` random
deviates. The length of the result is determined by `n` for `rghermite`,
and is the maximum of the lengths of the numerical arguments for the
other functions.

## Details

A generalized Hermite random variable counts single events and events
that occur in batches of \\m\\: \$\$Y = N_1 + m N_m, \quad N_1 \sim
\mathrm{Poisson}((1 - \rho)\mu), \quad N_m \sim \mathrm{Poisson}(\rho\mu
/ m),\$\$ independently, so that \$\$E\[Y\] = \mu, \qquad
\mathrm{Var}\[Y\] = \mu(1 + (m - 1)\rho).\$\$ The variance-to-mean ratio
(dispersion index) \\d = 1 + (m - 1)\rho\\ lies between 1 and \\m\\: the
distribution is overdispersed but with bounded dispersion, and it can be
multimodal. \\\rho = 0\\ (or \\m = 1\\) gives the Poisson distribution;
with \\\rho = 1\\ all of the mass is on multiples of \\m\\, so other
counts have probability zero. In the hermite package parameterization
(order \\m\\, mean, dispersion index) the dispersion index is \\d\\.

\\m\\ is a fixed model choice, not a continuous parameter. Probabilities
are computed by summing the convolution of the two Poisson components on
the log scale; both tails are computed directly. Random deviates are
generated from the definition.

## References

Kemp, C. D. and Kemp, A. W. (1965). Some properties of the "Hermite"
distribution. *Biometrika*, 52(3/4), 381-394.

Puig, P. (2003). Characterizing additively closed discrete models by a
property of their maximum likelihood estimators, with an application to
generalized Hermite distributions. *Journal of the American Statistical
Association*, 98(463), 687-692.

Morina, D., Higueras, M., Puig, P. and Oliveira, M. (2015). Generalized
Hermite distribution modelling with the R package hermite. *The R
Journal*, 7(2), 263-274.

## Examples

``` r
dghermite(0:6, mu = 2, rho = 0.4)
#> [1] 0.20189652 0.24227582 0.22612410 0.15505653 0.09174178 0.04682707 0.02159765
# rho = 0 is Poisson
all.equal(dghermite(0:10, mu = 2, rho = 0), dpois(0:10, 2))
#> [1] TRUE
# rho = 1: only multiples of m
dghermite(0:6, mu = 2, rho = 1, m = 3)
#> [1] 0.5134171 0.0000000 0.0000000 0.3422781 0.0000000 0.0000000 0.1140927
pghermite(4, mu = 2, rho = 0.4, m = 3)
#> [1] 0.8953356
qghermite(c(0.5, 0.95), mu = 2, rho = 0.4)
#> [1] 2 5
set.seed(4)
y <- rghermite(5000, mu = 3, rho = 0.5, m = 3)
c(mean = mean(y), var = var(y))   # theory: 3 and 3 * (1 + 2 * 0.5) = 6
#>     mean      var 
#> 2.987000 5.854002 
```
