# Gamma Count Distribution

Density, distribution function, quantile function, and random generation
for the gamma renewal count distribution (Winkelmann, 1995): the number
of events in a time window when the times between events are independent
gamma variables.

## Usage

``` r
dgcount(x, shape, rate, time = 1, log = FALSE)

pgcount(q, shape, rate, time = 1, lower.tail = TRUE, log.p = FALSE)

qgcount(p, shape, rate, time = 1, lower.tail = TRUE, log.p = FALSE)

rgcount(n, shape, rate, time = 1)
```

## Arguments

- x:

  vector of (non-negative integer) quantiles.

- shape:

  vector of (positive) shape parameters \\a\\ of the gamma inter-arrival
  times. \\a = 1\\ is the Poisson distribution, \\a \> 1\\ gives
  underdispersion and \\a \< 1\\ overdispersion.

- rate:

  vector of (positive) rate parameters \\b\\ of the gamma inter-arrival
  times.

- time:

  vector of (positive) lengths \\t\\ of the observation window. Default
  1.

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

`dgcount` gives the probability mass, `pgcount` the distribution
function, `qgcount` the quantile function, and `rgcount` random
deviates. The length of the result is determined by `n` for `rgcount`,
and is the maximum of the lengths of the numerical arguments for the
other functions.

## Details

Let the times between events be independent Gamma(shape = \\a\\, rate =
\\b\\) variables and let \\N(t)\\ be the number of events completed by
time \\t\\ (the process starts with an event at time 0). With \\G(s;
c)\\ the distribution function of a Gamma(\\c\\, 1) variable at \\s =
bt\\ (`pgamma(s, c)`), \$\$P(N(t) \leq y) = 1 - G(bt; (y + 1)a), \qquad
P(N(t) = y) = G(bt; ya) - G(bt; (y + 1)a),\$\$ with \\G(bt; 0) = 1\\.
Only the product \\bt\\ enters the distribution. With \\a = 1\\ the
inter-arrival times are exponential and \\N(t)\\ is Poisson with mean
\\bt\\.

The mean is \\E\[N(t)\] = \sum\_{k \geq 1} G(bt; ka)\\, which is
approximately (but not exactly) \\bt / a\\ for large \\bt\\; the
expected count is generally not proportional to the window length \\t\\,
so lengthening the window is not the same as a log-exposure offset.

Probabilities are differences of gamma distribution functions, taken
between the pair of tail probabilities that is not close to one. When
the shape \\a\\ is very small the two terms are nearly equal and some
relative accuracy (about \\\log\_{10}(1/a)\\ digits) is lost. Both tails
of the distribution function are computed directly with
[`pgamma()`](https://rdrr.io/r/stats/GammaDist.html). `rgcount()` uses
inversion of the distribution function.

## References

Winkelmann, R. (1995). Duration dependence and dispersion in count-data
models. *Journal of Business & Economic Statistics*, 13(4), 467-474.

Kharrat, T., Boshnakov, G. N., McHale, I. and Baker, R. (2019). Flexible
regression models for count data based on renewal processes: the Countr
package. *Journal of Statistical Software*, 90(13), 1-35.

## Examples

``` r
dgcount(0:6, shape = 2, rate = 6)
#> [1] 0.01735127 0.13385262 0.29447576 0.29830012 0.17209622 0.06383205 0.01646347
# shape = 1 gives the Poisson distribution with mean rate * time
all.equal(dgcount(0:10, shape = 1, rate = 1.5, time = 2), dpois(0:10, 3))
#> [1] TRUE
pgcount(3, shape = 0.5, rate = 1.5)
#> [1] 0.5578254
qgcount(c(0.1, 0.5, 0.9), shape = 2, rate = 6)
#> [1] 1 3 4
set.seed(6)
y <- rgcount(2000, shape = 3, rate = 12)
c(mean = mean(y), var = var(y))   # underdispersed
#>     mean      var 
#> 3.676500 1.428562 
```
