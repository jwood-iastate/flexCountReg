# flexCountReg 0.1.2 (development)

## tidymodels and broom

* `tidy()`, `glance()` and `augment()` methods, using the generics package
  (re-exported, so `library(flexCountReg)` provides them). `tidy()` returns
  the standard broom columns plus a `component` column (mean, distribution,
  underreporting, random-parameter and heterogeneity parameters), with
  `conf.int` and `exponentiate` (incidence rate ratios for the mean model).
  `tidy()` now works for every estimator. broom is no longer a dependency.
* New `count_reg()` parsnip model with a `"flexCountReg"` engine, registered
  when parsnip is loaded (parsnip is suggested, not required). The count
  `family` is a main argument and can be tuned with `count_family()` (a
  dials parameter). Case weights are passed as frequency weights.
* `predict()` accepts `new_data` (the tidymodels name) as well as `newdata`.
* Standard errors are `NA` (instead of 0) when the estimated variance is
  negative, which happens for weakly identified parameters. `regCompTable()`
  shows such standard errors as `NA` without significance stars.

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
  null models use the same rows, offset and weights as the fitted model.
  McFadden's R^2 is computed against an intercept-only Poisson model for
  every model by default (so it is comparable across distributions); a
  different baseline family, with or without the mean covariates, can be
  chosen with `basemodel` and `variables` in both `regCompTest()` and
  `regCompTable()`.
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
* `dgp1()` and `dgp2()` are vectorized, and for underdispersed values the
  normalizing constant is summed only where the probability is not
  negligible. Near phi = 0 (or alpha = 0) from below the old code summed up
  to millions of terms per observation, so GP1/GP2 fits could run for hours
  (for example the GP1 example in `?countreg`, which now takes under a
  second).
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

## Corrections from the development review

* PGE and PW: with the mean parameterization, the scale of the mixing
  distribution cancels out of the likelihood, so it was not identified (the
  PMF is identical for every scale). These families now estimate only the
  shape parameter (one fewer parameter, so AIC and BIC change), and
  `dis_param_formula_2` is an error for them. Predictions from fits made with
  earlier versions are unchanged. The `dpge()`/`dpoisweibull()` families of
  functions keep their scale arguments.
* GP1 and GP2 with a negative dispersion: the PMF is truncated and
  renormalized, so `mu` is a location parameter, not the mean (GP1 with
  mu = 0.75 and phi = -0.5 has mean 0.803). Predictions, Pearson residuals
  (variances), Empirical Bayes and random-parameter predictions now use the
  exact moments of the implemented PMF. Exact random-parameter predictions
  fall back to Simulated for these models, because the mean is not
  proportional to exp(X beta). The GP normalizing constants are computed
  vectorized, which makes GP random-parameter models much faster.
* Underreporting is binomial thinning of the true counts, and the model of
  the reported counts now matches it: for NB1 and PIG1 the reported counts'
  dispersion is p * alpha, and for NBP alpha * p^(2 - P) (before, the
  true-count dispersion was used unchanged, which is not the distribution of
  thinned counts). Poisson, NB2, PIG2, PIG and the mixed-Poisson families
  (PLN, PL, PLG, PLL, PGE, PW, SI) were already exact. GP1, GP2, COM and GW
  are not closed under thinning, so underreporting is now an error for them.
  This changes NB1, NBP and PIG1 underreporting estimates.
* Underreporting identification: an intercept-only `underreport_formula` is
  an error (only exp(b0) * p is identified), and `countreg()` warns when the
  reporting and count coefficients are nearly collinear (correlation above
  0.99), which happens when reporting is separated from the count mean only
  by the curvature of the link function.
* Likelihood-ratio tests are reported only for valid nesting. `summary()`
  no longer reports a test against Poisson for PL, PLG and PLL (which cannot
  reduce to Poisson: their variance is at least mu + mu^2 / 2), for NBP, SI,
  PW, PGE and GW (nonregular restrictions), or for boundary families with a
  dispersion formula. `regCompTest()` uses the same rules (same estimator
  and model components, nested parameters and families), uses a
  50:50 chi-square mixture when one parameter is on the boundary (also with
  extra mean covariates), requires two fitted models to use the same
  observations, outcome, offsets and weights, returns `LR_type`, and gives
  no p-value (instead of 1) when the larger model has the lower
  log-likelihood. The README zero-inflation example now compares models by
  AIC and BIC.
* Likelihoods use log-probabilities computed on the log scale (every family
  density has a `log` argument), so extreme counts keep finite
  log-likelihoods (log P(Y = 1000 | mu = 1) = -5913 instead of log(0)).
  Zero-inflated and hurdle adjustments are done on the log scale, and rows
  or panels with zero weight contribute exactly 0.
* The Conway-Maxwell-Poisson functions use the log normalizing constant (it
  overflowed for lambda = 1000, nu = 1, so `dcom(1000, lambda = 1000,
  nu = 1)` was 0 instead of 0.0126), one support rule for the PMF, CDF,
  quantiles and moments, moments computed from the normalized PMF, and a
  bracketed solver for the rate given the mean (accurate to about 1e-8,
  NaN with a warning when it fails, instead of silently returning an
  unconverged value).
* Correlated random parameters with three or more parameters: the `Chol:`
  labels and the default starting values were built row by row while the
  Cholesky factor was filled column by column, so labels were attached to
  the wrong elements, named starting values went to the wrong elements, and
  the default start was a singular matrix instead of 0.1 * I. Labels,
  starting values, the likelihood and predictions now use one packing.
  Estimates and predictions of earlier fits were internally consistent;
  only their labels were wrong.
* `countreg.rp()` rejects a variable in both `formula` and `rpar_formula`
  when its random parameter is normal, uniform or triangular (only the sum
  of the fixed coefficient and the random mean is identified). Lognormal
  and gamma random parameters with a fixed shift are still allowed.
* `countreg.rp()` starts GP1/GP2 at zero dispersion (the Poisson model)
  instead of a strongly underdispersed value that gave larger counts zero
  probability, so default GP fits no longer fail with "NA in the initial
  gradient".

## Bug fixes found by the extended test suite

* `pcom()` ignored `lower.tail` and `log.p`. Both tails are now summed
  directly on the log scale with a precise normalizing constant, so upper
  tails are never negative and keep their relative precision (e.g.
  `pcom(25, lambda = 1, nu = 1, lower.tail = FALSE)` equals the Poisson
  value 9.5e-28). Recycling of `mu`, `lambda` and
  `nu` in `dcom()`, `pcom()` and `qcom()` follows the base-R rule (length 1
  or the length of the first argument; other lengths are an error).
* `dpoisweibull()` and friends used `lambda` inconsistently between the PMF
  and the mean; parameters given only partly now warn instead of failing.
* The d/p/q functions of the Poisson-Lindley, Poisson-Lindley-Lognormal,
  Poisson-inverse-Gaussian, Poisson-inverse-gamma, Poisson-generalized-
  exponential, Sichel, GP1 and inverse-gamma distributions return `NaN`
  (with one warning) for invalid parameters
  instead of failing, returning 0, or emitting extra warnings. Missing inputs
  give `NA`.
* `qsichel()` returns `Inf` at p = 1 (it returned `NA` with a coercion
  warning) and `NaN` for p outside [0, 1] (it searched 10,000 counts).
  `qgp1()` no longer fails for invalid parameters. `rsichel()` and
  `rinvgamma()` accept a vector `n` (its length is used, as in base R).
* `halton_dists()` draws one dimension by default, uses the stated
  half-width for uniform draws (`mean +/- sdev`), and errors for an unknown
  distribution. `corr_haltons()` gives each column its own mean, standard
  deviation and bounds for truncated normal draws, and errors for an
  unsupported distribution, `hdraws` with the wrong number of columns, or a
  Cholesky factor of the wrong size.
* `margEffTable()` errors for misspelled variables (they were dropped
  silently) and refits bootstrap replicates in the environment of the
  model formula.
* `confint()`, `summary()` and `tidy()` work (with `NA` standard errors) when
  a model has no covariance matrix.
* `bootstraps` is checked up front (a single number of at least 2) instead of
  failing after the model is fit.
* `repois()` uses the shared data preparation and panel bootstrap: the
  bootstrap covariance matrix is computed from the coefficient draws, the
  rows used are recorded, and a prediction blueprint is stored, so
  `predict()` and `augment()` work with `keep_data = FALSE`.
* Data preparation rejects infinite counts and non-numeric (e.g. factor)
  responses, expands `.` in formulas, and keeps the contrasts of factors in
  auxiliary formulas for prediction. Fitted objects store the expanded
  formulas, so `margEffTable()` and bootstrap refits of models written with
  `.` work.
* `method = "Nelder-Mead"` is accepted by every estimator as an alias of
  `"NM"`, and `repois()` accepts all optimizers the other estimators accept.
* `qplind()` and `qplindLnorm()` return `NA` for missing parameters instead
  of failing, and `qplindLnorm()`/`rplindLnorm()` accept `sigma = 0` (the
  Poisson-Lindley limit).
* `regCompTable()` fits the baseline model only when the pseudo-R^2 is
  requested. `tidy()` classifies parameters from the model's bookkeeping
  (covariates named like distribution parameters are no longer
  misclassified), and `augment()` computes residuals from the transformed
  response (e.g. `log(y + 1) ~ x`), evaluated in the environment of the
  model formula.

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

* Dependencies: broom, dplyr, modelr, purrr, stringr and tidyr are no longer
  used; generics is added. Whole-namespace imports are replaced by
  `importFrom()`. lamW stays in Suggests (its uses in `qlindley()` and
  lognormal random-parameter predictions are checked with
  `rlang::check_installed()`); dials and parsnip are added to Suggests.
* `rlindley()` draws from the gamma mixture form directly (no longer needs
  lamW) and no longer warns for every whole-number `n`.
* `src/Makevars` and `src/Makevars.win` are removed (no OpenMP or BLAS is
  used), so `SystemRequirements: GNU make` is dropped.
* Added an R CMD check GitHub Actions workflow (including a job with only
  the hard dependencies installed) and a grouped pkgdown reference index.
  Removed the commented-out `tests/testthat/test-all.R` and a duplicated
  test.

* The `poilog` engine is checked with `rlang::check_installed()`.
* `msg1` and `.dpge_one` are no longer exported. The hand-maintained
  `Collate` field is removed.

## New features

* Zero-inflated and hurdle models for every `countreg()` family, through the
  new `zero_formula` and `zero_type` arguments (checked against
  `pscl::zeroinfl()` and `pscl::hurdle()`).
* `predict(type = "prob")` (probabilities of each count), randomized quantile
  residuals (`residuals(type = "quantile")`) and `rootogram()`.
* Analytic gradients for the Poisson, NB1 and NB2 families (with or without a
  dispersion formula), used by the gradient-based optimizers and by robust
  standard errors.
* Internally, each count family is defined once in a family registry
  (`R/families.R`), which makes adding a family a one-file change.

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
