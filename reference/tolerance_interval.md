# Tolerance Intervals for Count Regression Models

Content tolerance intervals for the observed count of a
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
model: integer intervals intended to contain at least a fraction
`content` of the count distribution at given covariate values, with
confidence `confidence` over repeated samples.

## Usage

``` r
tolerance_interval(
  object,
  newdata = NULL,
  content = 0.95,
  confidence = 0.95,
  side = c("two-sided", "upper", "lower"),
  scope = c("pointwise", "simultaneous"),
  method = c("auto", "joint", "envelope", "exact"),
  uncertainty = c("normal", "bootstrap"),
  npoints = 400,
  bootstraps = 200,
  seed = NULL,
  max_count = 1e+05
)
```

## Arguments

- object:

  a model fitted with
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md).

- newdata:

  a data frame of covariate profiles (scenarios). Defaults to the
  training rows.

- content:

  the required fraction \\p\\ of the count distribution (e.g. 0.95).

- confidence:

  the confidence level \\\gamma\\ (e.g. 0.95).

- side:

  `"two-sided"` (default), `"upper"` (interval `[0, U]`) or `"lower"`
  (interval `[L, Inf)`).

- scope:

  `"pointwise"` (default): each row's interval has the stated confidence
  on its own; `"simultaneous"`: the confidence holds jointly for all
  covariate profiles (based on a confidence set for all parameters, so
  wider).

- method:

  `"auto"` (default): `"exact"` when it applies, otherwise `"joint"`.
  `"joint"`: the narrowest integer interval whose content is at least
  `content` for every parameter value in the confidence set.
  `"envelope"`: equal-tail quantiles taken at the worst-case parameter
  values of each tail (simpler, more conservative). `"exact"`: the
  common-rate Poisson construction, available for intercept-only Poisson
  models without auxiliary components.

- uncertainty:

  how the confidence set is formed: `"normal"` (default), the asymptotic
  normal (Wald) confidence ellipsoid of the estimates, or `"bootstrap"`,
  a parametric-bootstrap confidence region. `"normal"` switches to the
  bootstrap, with a message, when the covariance matrix is not positive
  definite.

- npoints:

  number of points of the confidence set at which the distribution is
  evaluated (normal method). They are deterministic (based on Halton
  sequences), so results are reproducible.

- bootstraps:

  number of parametric-bootstrap refits.

- seed:

  random seed for the parametric bootstrap (the caller's random-number
  stream is not changed). `NULL` uses the current stream.

- max_count:

  largest count used when resolving the upper tail.

## Value

A data frame with one row per row of `newdata`: `fit` (the expected
count), `lower` and `upper` (integer limits; `upper` is `Inf` for
`side = "lower"`), `min_content` (the smallest content of the interval
over the evaluated points of the confidence set) and `status` (`"ok"`,
or the reason a limit is unavailable). The settings are stored in the
attribute `"tolerance"`.

## Details

**What the interval means.** For covariates \\x\\, an integer interval
\\\[L, U\]\\ computed from the estimation sample \\D\\ is a content
tolerance interval if \$\$\Pr_D\\F\_\theta(U \mid x) - F\_\theta(L - 1
\mid x) \ge p\\ \ge \gamma,\$\$ where \\F\_\theta\\ is the distribution
function of the observed count under the true parameters, \\p\\ the
content and \\\gamma\\ the confidence. This differs from a prediction
interval (which covers one future count on average; see
[`predict.flexCountReg()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md))
and from a confidence interval for the mean. The statement is
conditional on the model being correct.

**Construction.** The distribution of the count at \\x\\ depends on the
parameters only through the linear predictors of the model components at
\\x\\ (mean, distribution-parameter formulas, reporting and zero
models). A confidence set for those linear predictors is formed (Wald:
an ellipsoid with \\\chi^2\\ cutoff, with as many degrees of freedom as
linear predictors for pointwise intervals, or as parameters for
simultaneous intervals; bootstrap: the bootstrap draws closest to the
estimate in Mahalanobis distance, a fraction \\\gamma\\ of them). The
distribution is evaluated at `npoints` points of the set (its boundary,
interior and center; deterministic Halton-based points) and the interval
is chosen so that its content is at least \\p\\ at every evaluated
point. When the true parameter lies in the confidence set the interval
has the required content, so the confidence is approximately at least
\\\gamma\\; it is approximate because the confidence set is asymptotic
(or bootstrap-based) and the extrema are found on a finite set of
points. Discrete counts generally give more than the requested content;
endpoints are never interpolated.

**Exact Poisson method.** For an intercept-only Poisson model (a common
rate \\\lambda\\, with exposures from an offset), an exact confidence
interval for \\\lambda\\ from
[`stats::poisson.test()`](https://rdrr.io/r/stats/poisson.test.html) is
mapped to Poisson quantiles at the scenario's exposure, giving an
interval with at least the stated content and confidence (conservative
because of discreteness).

In a Monte Carlo study of Poisson and NB2 regressions (n = 50 and 200,
1000 samples each; content 0.90, confidence 0.95), the joint method
achieved confidence 0.99 to 1.00 and plug-in quantile intervals (which
ignore estimation uncertainty) as little as 0.54. The interval targets
the observed (reported) count, including zero-inflated, hurdle and
underreporting components. Random-parameter and panel models are not yet
supported.

## References

Krishnamoorthy, K. and Mathew, T. (2009). *Statistical Tolerance
Regions: Theory, Applications, and Computation*. Wiley.

Young, D. S. (2010). tolerance: An R package for estimating tolerance
intervals. *Journal of Statistical Software*, 36(5), 1-39.

## See also

[`predict.flexCountReg()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md)
for confidence and prediction intervals.

## Examples

``` r
data("washington_roads")
fit <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads,
                family = "NB2", method = "BFGS")
scen <- data.frame(lnaadt = log(c(2000, 10000)), lnlength = log(1))
tolerance_interval(fit, scen, content = 0.90, confidence = 0.95)
#> Tolerance intervals: content 0.9, confidence 0.95, two-sided, pointwise
#> Method: joint (normal)
#>         fit lower upper min_content status
#> 1 0.4817677     0     2   0.9605276     ok
#> 2 2.9030273     0     8   0.9208954     ok

# Exact common-rate Poisson method
set.seed(1)
d <- data.frame(y = rpois(30, 2))
p0 <- countreg(y ~ 1, data = d, family = "POISSON", method = "BFGS")
tolerance_interval(p0, data.frame(x = 1), content = 0.95)
#> Tolerance intervals: content 0.95, confidence 0.95, two-sided, pointwise
#> Method: exact (exact)
#>   fit lower upper min_content status
#> 1   2     0     6   0.9836322     ok
```
