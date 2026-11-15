# flexCountReg 0.1.3 (development)

## Breaking changes and changed results

* `summary()` now returns a `summary.flexCountReg` object with a print method
  (coefficients with z statistics, the standard-error type, log-likelihood,
  AIC, BIC, N and panels, convergence, and a likelihood-ratio test against
  Poisson for `countreg()` models). Use `summary(m)$coefficients` for the
  coefficient matrix, or `tibble::as_tibble(summary(m))` for the tibble that
  `summary()` returned before. The offset is shown as a note, not as a
  coefficient row.
* Fitted objects have class `c(<estimator>, "flexCountReg")`, with
  `<estimator>` one of `countreg`, `countreg_rp`, `renb`, `poislind_re`,
  `repois` and `negmultinom`. The unused S4 class definition is removed. The
  training data are stored once (`obj$data`); `obj$model$data` is gone.
* `regCompTest()` with two models now treats the model with fewer parameters
  as the restricted model (before: the one with the higher log-likelihood),
  errors when the models were fit to different numbers of observations, and
  warns when the restricted model has the higher log-likelihood. Its base and
  null models use the same rows, offset and weights as the fitted model, and
  McFadden's R^2 uses an intercept-only model of the same family for
  `countreg()` models. `regCompTable()` uses the same definitions.
* `cor2cov()` stops (instead of warning and returning `NULL`) on invalid
  input, and checks that `C` is a symmetric correlation matrix.
* `margEffTable()` includes by default the variables of the underreporting,
  random-parameter and heterogeneity formulas (and the dispersion formulas of
  PLN and PLL models), which also change the expected value.

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
* `countreg.rp()`: heterogeneity in means now follows Wood and Gayah (2025,
  Eq. 2). Each random parameter has its own coefficient for each variable in
  `het_mean_formula` (named `HetMean:<parameter>:<variable>`), and the effect
  is additive (on the log scale for lognormal parameters; multiplicative for
  gamma parameters, whose mean must stay positive). Previously a single
  coefficient vector multiplied every random parameter's mean. Estimates from
  models with `het_mean_formula` change; saved fits still predict with their
  original specification.
* `countreg.rp()`: uniform random parameters now span mu +/- St.Dev (draws
  `2u - 1`), matching the paper and the Exact prediction. Previously the draws
  spanned mu +/- St.Dev/2, so estimated uniform St.Dev values are now half as
  large. `predict()` converts saved uniform fits automatically (with a
  one-time message).
* `countreg.rp()`: `rpardists` is matched to the random parameters by name
  (not position) and validated. An intercept in both `formula` and
  `rpar_formula` (an unidentified model) is an error; the intercept implied by
  `rpar_formula` is not a random parameter unless the fixed formula has no
  intercept. An invalid `panel_id`, `family` or `method` is an error.
* `countreg.rp()`: `weights` now weight the panel log-likelihood (frequency
  weights for whole panels) instead of each row's probability.
* `countreg.rp()`: the mixing distribution of PLN, PGE, PLL and PW families
  uses its own Halton dimension; it previously reused the first random
  parameter's draws, which made the two perfectly correlated.
* `predict(method = "Exact")` returns `Inf` (with a warning) where the
  expected value does not exist: lognormal random parameters with a positive
  x, and gamma random parameters beyond the domain of their moment generating
  function. The lognormal formula is documented as a Laplace approximation.
* `renb()` and `poisLind.re()` require `group_var`, and the panels always come
  from it (a `panel_id` column in `data` is no longer used instead).
  Bootstrapped standard errors resample whole panels instead of rows.
* `pplindLnorm()`, `qplindLnorm()` and `rplindLnorm()` convert `lambda` to
  `mean` consistently with `dplindLnorm()` (the extra `exp(sigma^2 / 2)`
  factor is removed).
* `dpLnorm()`, `dpinvgaus()`, `dpinvgamma()`, `dplind()` and `dgwar()` are
  vectorized directly instead of through `Vectorize()`. Off-support `x` gives
  0 with a warning, and invalid parameters give `NaN` with a warning.

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

## Bug fixes (distributions and panel models)

* `dpinvgaus()` and `dpinvgamma()` are computed on the log scale, so they no
  longer return `Inf` or `NaN` at large counts (e.g. y > 85 with mu = 2,
  eta = 0.5).
* r-functions (`rpLnorm()`, `rplind()`, `rplindLnorm()`, `rplindGamma()`,
  `rpge()`, `rcom()`, `rgp1()`, `rgp2()`) accept vector parameters, as needed
  for simulation from fitted models. Mixture families draw the mixing
  variable per observation.
* q-functions no longer loop forever at p = 1 (they return `Inf`), return
  `NaN` for p outside [0, 1], and no longer recompute the CDF at every step.
* `dplind(lambda = )` uses `lambda` (it previously ignored it).
* The RENB likelihood computes log-gamma differences without cancellation for
  large means. Previously the likelihood was overstated far from the optimum,
  and optimizers such as BFGS could converge to meaningless estimates.
* `poisLind.re()` starts from a moment estimate of theta instead of an
  undocumented formula.
* `renb()` and `poisLind.re()` use the same data preparation as `countreg()`:
  missing values are removed consistently, formula offsets work (giving both a
  formula offset and `offset` is an error), and starting values use the
  offset.
* `countreg.rp()` uses the shared data preparation (formula offsets, missing
  values, numeric weights), starts gamma random parameters from a positive
  mean, averages panel likelihoods over draws on the log scale, and reports
  convergence.
* `predict(method = "Exact")` for `countreg.rp()` includes heterogeneity in
  means and variances.

## Documentation

* Corrected the Generalized Poisson parameter restrictions, the parameter
  table in `?countreg` (PGE, Sichel), the Poisson-Lindley-Lognormal mean and
  parameter descriptions, the RENB density, the Sichel reference (Zou, Lord
  and Zhang, 2012) and several typos. `?countreg` no longer claims that an LR
  test is printed; see `regCompTest()`.
* `?countreg.rp` documents each random-parameter distribution's scale, the
  heterogeneity specification and the lognormal caveat.
* `?washington_roads` describes the panel structure; the `Rollover` label in
  the data is corrected. The README function list is corrected and completed.

## Other changes

* `lamW` moves from Suggests to Imports (used by `qlindley()` and lognormal
  predictions). The `poilog` engine is checked with `rlang::check_installed()`.
* `msg1` and `.dpge_one` are no longer exported. The hand-maintained
  `Collate` field is removed.

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
* `print.level` in `renb()`, `poisLind.re()` and `repois()` is deprecated; use
  `verbose`. `panel_id` in `countreg.rp()` is deprecated; use `group_var`.
  The old names still work, with a once-per-session warning.

## Methods and API consistency

* New methods for all models: `coef()`, `vcov()`, `logLik()`, `nobs()`,
  `fitted()`, `residuals()` (response and, for `countreg()`, Pearson),
  `confint()`, `formula()`, `terms()` and `print()`. `AIC()`, `BIC()` and
  `update()` now work. For panel models, `logLik()` reports the number of
  panels as its `nobs`, so `BIC()` uses the number of independent units.
* `simulate()` draws responses from fitted `countreg()` models (all
  families, including underreporting).
* Every estimator stores one covariance matrix (`vcov()`), and `renb()`,
  `poisLind.re()`, `repois()` and `countreg.rp()` gain `stderr = c("normal",
  "robust", "boot")`. Robust standard errors of panel models are clustered by
  panel; bootstraps resample whole panels. `countreg.rp()` gains
  `bootstraps`.
* All estimators gain `keep_data`; `keep_data = FALSE` does not store the
  training data.
* `countreg_rp()` and `poislind_re()` are snake_case names for
  `countreg.rp()` and `poisLind.re()`.
* `repois()` predictions work (its linear-predictor coefficients were
  truncated to one value), and its likelihood returns one value per cluster.
* For fitted models, `AIC()` and `BIC()` replace `myAIC()` and `myBIC()`.

# flexCountReg 0.1.1

* Initial CRAN submission.
