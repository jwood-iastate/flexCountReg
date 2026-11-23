# Predictions for flexCountReg models

Generates predictions for the expected count (lambda) for observations.

For **countreg_rp** (Random Parameters) models, three methods are
available:

- **Simulated**: Uses Halton draws to simulate the random parameters and
  averages the outcomes. This is a simulation-based approximation.

- **Individual** (alias `EB`): Empirical Bayes observation-specific
  predictions (Wood and Gayah, 2025, Eqs. 5-7). Each simulation draw is
  weighted by the likelihood of the unit's observed outcomes (all of a
  panel's rows when `panel_id` was used), giving conditional
  coefficients \\E\[\beta_i \| y_i\]\\. By default these are used in the
  mean function, as in the paper; pass `eb_type = "mean"` to instead
  return the conditional mean \\E\[\mu\_{it} \| y_i\]\\. The conditional
  coefficients and their conditional standard deviations are returned as
  the attributes `"eb_coefficients"` and `"eb_sd"`. Requires the outcome
  (and the `panel_id` column, if used) in `newdata`, so it is meant for
  data with observed outcomes.

- **Exact**: Uses the moment generating functions (MGFs) of the random
  parameter distributions (Wood and Gayah, 2025, Appendix A), including
  heterogeneity in means and variances. This is exact for normal,
  uniform, triangular and gamma parameters and has no simulation error.
  For lognormal parameters with \\x \< 0\\ the lognormal MGF has no
  closed form, and the Laplace approximation of Asmussen, Jensen and
  Rojas-Nandayapa (2016) is used: accurate, but approximate despite the
  method's name (it needs the lamW package). For \\x \> 0\\ with a
  lognormal parameter (and for gamma parameters with \\x\sigma^2/\mu \ge
  1\\) the expected value is infinite, and `Inf` is returned with a
  warning, for both `Exact` and `Simulated` predictions. Whether the
  expected value exists is decided from the parameters alone, so
  `Simulated` predictions do not need lamW. A finite expected value too
  large to represent is reported as an overflow.

For **countreg**, **poisLindRE**, and **RENB** models, the function
calculates the expected value \\\mu = \exp(X\beta)\\ (with appropriate
adjustments for specific families like PLN, PLL or underreporting).

## Usage

``` r
# S3 method for class 'flexCountReg'
predict(
  object,
  newdata = NULL,
  ...,
  type = c("response", "true", "prob_reported", "link", "prob", "zero"),
  new_data = NULL,
  interval = c("none", "confidence", "prediction"),
  level = 0.95
)
```

## Arguments

- object:

  a model object estimated using this R package.

- newdata:

  optional data frame for which to generate predictions. It must contain
  every variable used by the model's formulas (including `rpar_formula`,
  `dis_param_formula_1/2`, `underreport_formula` and heterogeneity
  formulas) and any offset columns. The outcome is not required. If
  omitted, predictions are made for the data used to fit the model.

- ...:

  optional arguments passed to the function. This includes `method`.
  Passing the data as `data =` is deprecated; use `newdata`.

- type:

  the type of prediction. `"response"` (default) is the expected
  observed count. For models with underreporting this is \\\mu \times
  p\\, the true mean times the probability of being reported. `"true"`
  is the true expected count \\\mu\\ (excluding the underreporting
  probability), `"prob_reported"` is the probability of being reported
  \\p\\, and `"link"` is the linear predictor \\X\beta\\ plus the
  offset. For models without underreporting, `"true"` equals
  `"response"` and `"prob_reported"` is 1. For zero-inflated and hurdle
  models (`zero_formula`), `"response"` includes the zero model,
  `"true"` is the mean of the count distribution alone, and `"zero"` is
  the probability from the zero model: the probability of a structural
  zero (zero-inflated) or of a zero count (hurdle). `"prob"` returns a
  matrix of probabilities \\P(Y = k)\\ with one row per row of `newdata`
  and one column per count \\k = 0, \dots,\\ `max_count` (passed in
  `...`; default the largest count in the training data), for rootograms
  and for scoring models out of sample. With `overflow = TRUE` (in
  `...`), a last column holds the tail mass \\P(Y \> \\
  `max_count`\\)\\, so each row sums to 1. Types other than `"response"`
  are available for `countreg` models only.

- new_data:

  the same as `newdata`, under the name used by tidymodels. Give one or
  the other.

- interval:

  `"none"` (default), `"confidence"` or `"prediction"` (`countreg`
  models). See the section *Intervals*.

- level:

  the confidence or prediction level of the interval.

## Value

A numeric vector of predicted expected counts for each observation in
the provided data. If no data is provided, the predictions for the data
used in estimating the model are provided. With
`interval = "confidence"` or `"prediction"`, a data frame with columns
`fit`, `lwr` and `upr` (see the section Intervals).

## Details

The returned vector always has one element per row of `newdata`, in the
same order. Rows with missing values in any variable the estimate needs
get `NA`: a row with a missing distribution-parameter, zero-model or
reporting covariate has no count distribution, so its probabilities
(`type = "prob"`) and intervals are `NA`, while its expected count is
returned when it does not depend on that parameter (for example the NB2
mean) and is `NA` otherwise. Factor levels and contrasts from the
training data are used, so `newdata` may contain only a subset of the
levels. An error is raised if an offset column is missing from
`newdata`.

## Note

optional parameter `method`: Only valid for random parameters models
(`countreg_rp`). Options include `Simulated` (default), `Individual` (or
`EB`), or `Exact`.

optional parameter `eb_type`: for `method = "Individual"`, either
`"coefficients"` (default) or `"mean"`.

## Intervals

For `countreg` models, `interval = "confidence"` gives a confidence
interval for the expected count (`type = "response"`, the expected
observed count, or `"true"`), and `interval = "prediction"` a prediction
interval for one future observed count (`type = "response"`). The result
is then a data frame with columns `fit`, `lwr` and `upr`, with the
settings in the attribute `"interval"`. For content tolerance intervals
see
[`tolerance_interval`](https://jwood-iastate.github.io/flexCountReg/reference/tolerance_interval.md).

Estimation uncertainty is represented by parameter vectors chosen by the
optional argument `uncertainty` (in `...`):

- `"normal"` (default): the asymptotic normal distribution of the
  estimates, \\\hat\theta + R^\top z\\ with \\R^\top R\\ the covariance
  matrix and \\z\\ the standard normal quantiles of `ndraws`
  (default 500) Halton points. Halton sequences are deterministic, so
  the intervals are reproducible and have less simulation error than
  pseudo-random draws; for models with many parameters,
  `scrambled = TRUE` (with a `seed`) uses scrambled Halton points, which
  avoid the correlation between high dimensions of the standard
  sequence. When the covariance matrix is not positive definite (e.g. a
  dispersion parameter at its boundary) the parametric bootstrap is used
  instead, with a message.

- `"bootstrap"`: a parametric bootstrap. `bootstraps` (default 200)
  response vectors are simulated from the fitted model on the training
  rows and the model is refit to each with its frozen specification;
  refits that fail or do not converge are rejected and counted. Use
  `seed` for reproducibility. Not available for models with weights.

The confidence interval takes the \\(1 - level)/2\\ and \\(1 +
level)/2\\ quantiles of the expected count over these parameter vectors.
The prediction interval uses the predictive distribution of the count,
the probability mass function of the observed count averaged over the
parameter vectors, and returns its integer quantiles: the smallest
counts \\L\\ and \\U\\ with \\P(Y \le L) \ge (1 - level)/2\\ and \\P(Y
\le U) \ge (1 + level)/2\\. Because counts are discrete, the predictive
coverage is usually above `level`; the limits are not interpolated. Both
include zero-inflated, hurdle and underreporting components. `max_count`
(in `...`, default 1e5) bounds the counts used to resolve the upper
tail.

In a Monte Carlo study of Poisson and NB2 regressions (n = 30 and 200),
prediction intervals at level 0.90 covered 0.93 to 0.97 of future
counts, while plug-in intervals that ignore estimation uncertainty fell
to 0.88 for the smaller sample.

## References

Wood, J.S., Gayah, V. (2025). Out-of-sample prediction and
interpretation for random parameter generalized linear models. *Accident
Analysis and Prevention*, 220, 108147.

## Examples

``` r
# \donttest{
# Load data and create a dummy variable
data("washington_roads")
washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)

# =========================================================================
# 1. Fixed Parameter Model (countreg)
# =========================================================================
nb2_fixed <- countreg(Total_crashes ~ lnaadt + lnlength + speed50,
                      data = washington_roads, 
                      family = "NB2")
pred_fixed <- predict(nb2_fixed, newdata = washington_roads)

# =========================================================================
# 2. Random Parameters Model (countreg_rp)
# =========================================================================
rp_nb2 <- countreg_rp(Total_crashes ~ lnaadt + lnlength,
                      rpar_formula = ~ -1 + speed50,
                      data = washington_roads,
                      family = "NB2",
                      rpardists = c(speed50 = "n"),
                      ndraws = 100)

# Method A: Simulated (Default)
pred_sim <- predict(rp_nb2, newdata = washington_roads, method = "Simulated")

# Method B: Exact (Analytical MGF)
pred_exact <- predict(rp_nb2, newdata = washington_roads, method = "Exact")

# Method C: Individual (Empirical Bayes, requires the observed outcome)
pred_eb <- predict(rp_nb2, newdata = washington_roads, method = "Individual")

# =========================================================================
# 3. Random Effects Models (poisLindRE / RENB)
# =========================================================================
pl_re <- poisLind.re(Total_crashes ~ lnaadt + lnlength,
                     data = washington_roads,
                     group_var = "ID")
pred_pl_re <- predict(pl_re, newdata = washington_roads)
# }
```
