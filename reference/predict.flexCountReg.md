# Predictions for flexCountReg models

Generates predictions for the expected count (lambda) for observations.

For **countreg.rp** (Random Parameters) models, three methods are
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
  new_data = NULL
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

## Value

A numeric vector of predicted expected counts for each observation in
the provided data. If no data is provided, the predictions for the data
used in estimating the model are provided.

## Details

The returned vector always has one element per row of `newdata`, in the
same order. Rows with missing values in any variable the model uses get
`NA`. Factor levels and contrasts from the training data are used, so
`newdata` may contain only a subset of the levels. An error is raised if
an offset column is missing from `newdata`.

## Note

optional parameter `method`: Only valid for random parameters models
(`countreg.rp`). Options include `Simulated` (default), `Individual` (or
`EB`), or `Exact`.

optional parameter `eb_type`: for `method = "Individual"`, either
`"coefficients"` (default) or `"mean"`.

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
# 2. Random Parameters Model (countreg.rp)
# =========================================================================
rp_nb2 <- countreg.rp(Total_crashes ~ lnaadt + lnlength,
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
