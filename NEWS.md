# flexCountReg 0.1.3 (development)

## Breaking changes and changed results

* Probit underreporting now uses P(reported) = pnorm(z'g), the same sign
  convention as logit (a positive coefficient means more likely to be
  reported). Previously probit used 1 - pnorm(z'g), so probit coefficients had
  the opposite sign. Refit probit models to update them; `predict()` on probit
  fits made with earlier versions still uses the old convention.
* `countreg()`: giving an offset both with `offset()` in the formula and with
  the `offset` argument is now an error. Before, a formula `offset()` was
  silently ignored.
* `countreg()`: an invalid `method`, `family`, `stderr` or
  `underreport_family` is now an error instead of a warning or a silent
  fallback. `dis_param_formula_1`/`_2` given for a parameter the family does
  not have is an error.
* `countreg()` starting values now come from a Poisson GLM that uses the same
  rows, offset and weights as the model, so results can differ slightly
  (usually closer to the optimum).
* `predict()` errors when an offset column used in fitting is missing from
  `newdata`, instead of silently setting the offset to 0.

## Bug fixes

* `countreg()`: offsets given with `offset()` in the formula are now used in
  estimation and prediction.
* `countreg()`: bootstrapped standard errors are matched to parameters by name.
  They were previously sorted alphabetically and printed against the wrong
  coefficients. Failed bootstrap replicates are counted and reported.
* `countreg()`: all data (main formula, `dis_param_formula_1/2`,
  `underreport_formula`, offset and weights) are taken from one set of complete
  cases, so missing values no longer misalign rows. The number of removed rows
  is reported.
* `countreg()`: `stderr = "robust"` and `stderr = "boot"` now set the stored
  standard errors and covariance matrix (`$model$vcov`, `$model$vcov_type`),
  which `summary()` and `margEffTable()` use. `stderr = "boot"` without
  `bootstraps` uses 200 replicates. The SE type is stored in `$model$stderr`.
* `countreg()`: parameter values where a distribution's density fails no
  longer abort the fit; the likelihood returns `NA` so the optimizer steps
  back.
* `countreg()`: `start.vals` must have one value per parameter, and named
  values are matched by name. `$model$converged` is now set (with a warning
  when the optimizer does not report convergence).
* `predict()` returns one value per row of `newdata`, in order, with `NA` for
  rows with missing values, instead of recycling shorter vectors. Training
  factor levels and contrasts are used, so `newdata` can contain a subset of
  levels. Coefficients are matched by name.
* `predict()` sums multiple offset columns and uses formula offsets.
* `predict()` for PLL models now includes the exp(sigma^2 / 2) lognormal mean
  factor, as it already did for PLN.
* `countreg.rp()` now stores its offset, so `predict()` uses it. Simulated
  predictions now include heterogeneity in means and variances.

## New features

* `predict(method = "Individual")` (alias `"EB"`) for `countreg.rp()` models
  now computes Empirical Bayes observation-specific predictions following
  Wood and Gayah (2025), weighting each draw by the likelihood of the unit's
  observed outcomes. Previously it returned the simulated predictions.
  `eb_type = "mean"` returns the conditional mean instead of plugging in the
  conditional coefficients; the conditional coefficients and their standard
  deviations are returned as attributes.
* `predict()` gains `type = c("response", "true", "prob_reported", "link")`
  for `countreg()` models.
* `countreg()` accepts `weights` as a numeric vector as well as a column name.

## Deprecations

* `predict(object, data = ...)` is deprecated; use `newdata`.

# flexCountReg 0.1.1

* Initial CRAN submission.
