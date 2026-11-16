# Changelog

## flexCountReg 0.1.3 (development)

### tidymodels and broom

- [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`glance()`](https://generics.r-lib.org/reference/glance.html) and
  [`augment()`](https://generics.r-lib.org/reference/augment.html)
  methods, using the generics package (re-exported, so
  [`library(flexCountReg)`](https://jwood-iastate.github.io/flexCountReg/)
  provides them).
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns the
  standard broom columns plus a `component` column (mean, distribution,
  underreporting, random-parameter and heterogeneity parameters), with
  `conf.int` and `exponentiate` (incidence rate ratios for the mean
  model). [`tidy()`](https://generics.r-lib.org/reference/tidy.html) now
  works for every estimator. broom is no longer a dependency.
- New
  [`count_reg()`](https://jwood-iastate.github.io/flexCountReg/reference/count_reg.md)
  parsnip model with a `"flexCountReg"` engine, registered when parsnip
  is loaded (parsnip is suggested, not required). The count `family` is
  a main argument and can be tuned with
  [`count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/count_family.md)
  (a dials parameter). Case weights are passed as frequency weights.
- [`predict()`](https://rdrr.io/r/stats/predict.html) accepts `new_data`
  (the tidymodels name) as well as `newdata`.
- Standard errors are `NA` (instead of 0) when the estimated variance is
  negative, which happens for weakly identified parameters.
  [`regCompTable()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTable.md)
  shows such standard errors as `NA` without significance stars.

### Breaking changes and changed results

- [`summary()`](https://rdrr.io/r/base/summary.html) now returns a
  `summary.flexCountReg` object with a print method (coefficients with z
  statistics, the standard-error type, log-likelihood, AIC, BIC, N and
  panels, convergence, and a likelihood-ratio test against Poisson for
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  models). Use `summary(m)$coefficients` for the coefficient matrix, or
  `tibble::as_tibble(summary(m))` for the tibble that
  [`summary()`](https://rdrr.io/r/base/summary.html) returned before.
  The offset is shown as a note, not as a coefficient row.

- Fitted objects have class `c(<estimator>, "flexCountReg")`, with
  `<estimator>` one of `countreg`, `countreg_rp`, `renb`, `poislind_re`,
  `repois` and `negmultinom`. The unused S4 class definition is removed.
  The training data are stored once (`obj$data`); `obj$model$data` is
  gone.

- [`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md)
  with two models now treats the model with fewer parameters as the
  restricted model (before: the one with the higher log-likelihood),
  errors when the models were fit to different numbers of observations,
  and warns when the restricted model has the higher log-likelihood. Its
  base and null models use the same rows, offset and weights as the
  fitted model. McFadden’s R^2 is computed against an intercept-only
  Poisson model for every model by default (so it is comparable across
  distributions); a different baseline family, with or without the mean
  covariates, can be chosen with `basemodel` and `variables` in both
  [`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md)
  and
  [`regCompTable()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTable.md).

- [`cor2cov()`](https://jwood-iastate.github.io/flexCountReg/reference/cor2cov.md)
  stops (instead of warning and returning `NULL`) on invalid input, and
  checks that `C` is a symmetric correlation matrix.

- [`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
  includes by default the variables of the underreporting,
  random-parameter and heterogeneity formulas (and the dispersion
  formulas of PLN and PLL models), which also change the expected value.

- Probit underreporting now uses P(reported) = pnorm(z’g), the same sign
  convention as logit (a positive coefficient means more likely to be
  reported). Previously probit used 1 - pnorm(z’g), so probit
  coefficients had the opposite sign. Refit probit models to update
  them; [`predict()`](https://rdrr.io/r/stats/predict.html) on probit
  fits made with earlier versions still uses the old convention.

- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md):
  giving an offset both with
  [`offset()`](https://rdrr.io/r/stats/offset.html) in the formula and
  with the `offset` argument is now an error. Before, a formula
  [`offset()`](https://rdrr.io/r/stats/offset.html) was silently
  ignored.

- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md):
  an invalid `method`, `family`, `stderr` or `underreport_family` is now
  an error instead of a warning or a silent fallback.
  `dis_param_formula_1`/`_2` given for a parameter the family does not
  have is an error.

- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  starting values now come from a Poisson GLM that uses the same rows,
  offset and weights as the model, so results can differ slightly
  (usually closer to the optimum).

- [`predict()`](https://rdrr.io/r/stats/predict.html) errors when an
  offset column used in fitting is missing from `newdata`, instead of
  silently setting the offset to 0.

- [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md):
  heterogeneity in means now follows Wood and Gayah (2025, Eq. 2). Each
  random parameter has its own coefficient for each variable in
  `het_mean_formula` (named `HetMean:<parameter>:<variable>`), and the
  effect is additive (on the log scale for lognormal parameters;
  multiplicative for gamma parameters, whose mean must stay positive).
  Previously a single coefficient vector multiplied every random
  parameter’s mean. Estimates from models with `het_mean_formula`
  change; saved fits still predict with their original specification.

- [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md):
  uniform random parameters now span mu +/- St.Dev (draws `2u - 1`),
  matching the paper and the Exact prediction. Previously the draws
  spanned mu +/- St.Dev/2, so estimated uniform St.Dev values are now
  half as large. [`predict()`](https://rdrr.io/r/stats/predict.html)
  converts saved uniform fits automatically (with a one-time message).

- [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md):
  `rpardists` is matched to the random parameters by name (not position)
  and validated. An intercept in both `formula` and `rpar_formula` (an
  unidentified model) is an error; the intercept implied by
  `rpar_formula` is not a random parameter unless the fixed formula has
  no intercept. An invalid `panel_id`, `family` or `method` is an error.

- [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md):
  `weights` now weight the panel log-likelihood (frequency weights for
  whole panels) instead of each row’s probability.

- [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md):
  the mixing distribution of PLN, PGE, PLL and PW families uses its own
  Halton dimension; it previously reused the first random parameter’s
  draws, which made the two perfectly correlated.

- `predict(method = "Exact")` returns `Inf` (with a warning) where the
  expected value does not exist: lognormal random parameters with a
  positive x, and gamma random parameters beyond the domain of their
  moment generating function. The lognormal formula is documented as a
  Laplace approximation.

- [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)
  and
  [`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)
  require `group_var`, and the panels always come from it (a `panel_id`
  column in `data` is no longer used instead). Bootstrapped standard
  errors resample whole panels instead of rows.

- [`pplindLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md),
  [`qplindLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md)
  and
  [`rplindLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md)
  convert `lambda` to `mean` consistently with
  [`dplindLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md)
  (the extra `exp(sigma^2 / 2)` factor is removed).

- [`dpLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLognormal.md),
  [`dpinvgaus()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGaussian.md),
  [`dpinvgamma()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGamma.md),
  [`dplind()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindley.md)
  and
  [`dgwar()`](https://jwood-iastate.github.io/flexCountReg/reference/Generalized-Waring.md)
  are vectorized directly instead of through
  [`Vectorize()`](https://rdrr.io/r/base/Vectorize.html). Off-support
  `x` gives 0 with a warning, and invalid parameters give `NaN` with a
  warning.

### Bug fixes

- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md):
  offsets given with [`offset()`](https://rdrr.io/r/stats/offset.html)
  in the formula are now used in estimation and prediction.
- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md):
  bootstrapped standard errors are matched to parameters by name. They
  were previously sorted alphabetically and printed against the wrong
  coefficients. Failed bootstrap replicates are counted and reported.
- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md):
  all data (main formula, `dis_param_formula_1/2`,
  `underreport_formula`, offset and weights) are taken from one set of
  complete cases, so missing values no longer misalign rows. The number
  of removed rows is reported.
- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md):
  `stderr = "robust"` and `stderr = "boot"` now set the stored standard
  errors and covariance matrix (`$model$vcov`, `$model$vcov_type`),
  which [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
  use. `stderr = "boot"` without `bootstraps` uses 200 replicates. The
  SE type is stored in `$model$stderr`.
- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md):
  parameter values where a distribution’s density fails no longer abort
  the fit; the likelihood returns `NA` so the optimizer steps back.
- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md):
  `start.vals` must have one value per parameter, and named values are
  matched by name. `$model$converged` is now set (with a warning when
  the optimizer does not report convergence).
- [`predict()`](https://rdrr.io/r/stats/predict.html) returns one value
  per row of `newdata`, in order, with `NA` for rows with missing
  values, instead of recycling shorter vectors. Training factor levels
  and contrasts are used, so `newdata` can contain a subset of levels.
  Coefficients are matched by name.
- [`predict()`](https://rdrr.io/r/stats/predict.html) sums multiple
  offset columns and uses formula offsets.
- [`predict()`](https://rdrr.io/r/stats/predict.html) for PLL models now
  includes the exp(sigma^2 / 2) lognormal mean factor, as it already did
  for PLN.
- [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md)
  now stores its offset, so
  [`predict()`](https://rdrr.io/r/stats/predict.html) uses it. Simulated
  predictions now include heterogeneity in means and variances.

### Bug fixes (distributions and panel models)

- [`dpinvgaus()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGaussian.md)
  and
  [`dpinvgamma()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGamma.md)
  are computed on the log scale, so they no longer return `Inf` or `NaN`
  at large counts (e.g. y \> 85 with mu = 2, eta = 0.5).
- r-functions
  ([`rpLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLognormal.md),
  [`rplind()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindley.md),
  [`rplindLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md),
  [`rplindGamma()`](https://jwood-iastate.github.io/flexCountReg/reference/NegativeBinomialLindley.md),
  [`rpge()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonGeneralizedExponential.md),
  [`rcom()`](https://jwood-iastate.github.io/flexCountReg/reference/COMDistribution.md),
  [`rgp1()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedPoisson.md),
  [`rgp2()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedPoisson2.md))
  accept vector parameters, as needed for simulation from fitted models.
  Mixture families draw the mixing variable per observation.
- q-functions no longer loop forever at p = 1 (they return `Inf`),
  return `NaN` for p outside \[0, 1\], and no longer recompute the CDF
  at every step.
- `dplind(lambda = )` uses `lambda` (it previously ignored it).
- [`dgp1()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedPoisson.md)
  and
  [`dgp2()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedPoisson2.md)
  are vectorized, and for underdispersed values the normalizing constant
  is summed only where the probability is not negligible. Near phi = 0
  (or alpha = 0) from below the old code summed up to millions of terms
  per observation, so GP1/GP2 fits could run for hours (for example the
  GP1 example in
  [`?countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
  which now takes under a second).
- The RENB likelihood computes log-gamma differences without
  cancellation for large means. Previously the likelihood was overstated
  far from the optimum, and optimizers such as BFGS could converge to
  meaningless estimates.
- [`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)
  starts from a moment estimate of theta instead of an undocumented
  formula.
- [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)
  and
  [`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)
  use the same data preparation as
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md):
  missing values are removed consistently, formula offsets work (giving
  both a formula offset and `offset` is an error), and starting values
  use the offset.
- [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md)
  uses the shared data preparation (formula offsets, missing values,
  numeric weights), starts gamma random parameters from a positive mean,
  averages panel likelihoods over draws on the log scale, and reports
  convergence.
- `predict(method = "Exact")` for
  [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md)
  includes heterogeneity in means and variances.

### Documentation

- Corrected the Generalized Poisson parameter restrictions, the
  parameter table in
  [`?countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  (PGE, Sichel), the Poisson-Lindley-Lognormal mean and parameter
  descriptions, the RENB density, the Sichel reference (Zou, Lord and
  Zhang, 2012) and several typos.
  [`?countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  no longer claims that an LR test is printed; see
  [`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md).
- [`?countreg.rp`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md)
  documents each random-parameter distribution’s scale, the
  heterogeneity specification and the lognormal caveat.
- [`?washington_roads`](https://jwood-iastate.github.io/flexCountReg/reference/washington_roads.md)
  describes the panel structure; the `Rollover` label in the data is
  corrected. The README function list is corrected and completed.

### Other changes

- Dependencies: broom, dplyr, modelr, purrr, stringr and tidyr are no
  longer used; generics is added. Whole-namespace imports are replaced
  by `importFrom()`. lamW stays in Suggests (its uses in
  [`qlindley()`](https://jwood-iastate.github.io/flexCountReg/reference/Lindley.md)
  and lognormal random-parameter predictions are checked with
  [`rlang::check_installed()`](https://rlang.r-lib.org/reference/is_installed.html));
  dials and parsnip are added to Suggests.

- [`rlindley()`](https://jwood-iastate.github.io/flexCountReg/reference/Lindley.md)
  draws from the gamma mixture form directly (no longer needs lamW) and
  no longer warns for every whole-number `n`.

- `src/Makevars` and `src/Makevars.win` are removed (no OpenMP or BLAS
  is used), so `SystemRequirements: GNU make` is dropped.

- Added an R CMD check GitHub Actions workflow (including a job with
  only the hard dependencies installed) and a grouped pkgdown reference
  index. Removed the commented-out `tests/testthat/test-all.R` and a
  duplicated test.

- The `poilog` engine is checked with
  [`rlang::check_installed()`](https://rlang.r-lib.org/reference/is_installed.html).

- `msg1` and `.dpge_one` are no longer exported. The hand-maintained
  `Collate` field is removed.

### New features

- Zero-inflated and hurdle models for every
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  family, through the new `zero_formula` and `zero_type` arguments
  (checked against
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html) and
  [`pscl::hurdle()`](https://rdrr.io/pkg/pscl/man/hurdle.html)).

- `predict(type = "prob")` (probabilities of each count), randomized
  quantile residuals (`residuals(type = "quantile")`) and
  [`rootogram()`](https://jwood-iastate.github.io/flexCountReg/reference/rootogram.md).

- Analytic gradients for the Poisson, NB1 and NB2 families (with or
  without a dispersion formula), used by the gradient-based optimizers
  and by robust standard errors.

- Internally, each count family is defined once in a family registry
  (`R/families.R`), which makes adding a family a one-file change.

- `predict(method = "Individual")` (alias `"EB"`) for
  [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md)
  models now computes Empirical Bayes observation-specific predictions
  following Wood and Gayah (2025), weighting each draw by the likelihood
  of the unit’s observed outcomes. Previously it returned the simulated
  predictions. `eb_type = "mean"` returns the conditional mean instead
  of plugging in the conditional coefficients; the conditional
  coefficients and their standard deviations are returned as attributes.

- [`predict()`](https://rdrr.io/r/stats/predict.html) gains
  `type = c("response", "true", "prob_reported", "link")` for
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  models.

- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  accepts `weights` as a numeric vector as well as a column name.

### Deprecations

- `predict(object, data = ...)` is deprecated; use `newdata`.
- `print.level` in
  [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
  [`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)
  and
  [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md)
  is deprecated; use `verbose`. `panel_id` in
  [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md)
  is deprecated; use `group_var`. The old names still work, with a
  once-per-session warning.

### Methods and API consistency

- New methods for all models:
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) (response and,
  for
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
  Pearson), [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`formula()`](https://rdrr.io/r/stats/formula.html),
  [`terms()`](https://rdrr.io/r/stats/terms.html) and
  [`print()`](https://rdrr.io/r/base/print.html).
  [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`update()`](https://rdrr.io/r/stats/update.html) now work. For panel
  models, [`logLik()`](https://rdrr.io/r/stats/logLik.html) reports the
  number of panels as its `nobs`, so
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) uses the number of
  independent units.
- [`simulate()`](https://rdrr.io/r/stats/simulate.html) draws responses
  from fitted
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  models (all families, including underreporting).
- Every estimator stores one covariance matrix
  ([`vcov()`](https://rdrr.io/r/stats/vcov.html)), and
  [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
  [`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md),
  [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md)
  and
  [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md)
  gain `stderr = c("normal", "robust", "boot")`. Robust standard errors
  of panel models are clustered by panel; bootstraps resample whole
  panels.
  [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md)
  gains `bootstraps`.
- All estimators gain `keep_data`; `keep_data = FALSE` does not store
  the training data.
- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md)
  and
  [`poislind_re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)
  are snake_case names for
  [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md)
  and
  [`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md).
- [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md)
  predictions work (its linear-predictor coefficients were truncated to
  one value), and its likelihood returns one value per cluster.
- For fitted models, [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) replace
  [`myAIC()`](https://jwood-iastate.github.io/flexCountReg/reference/myAIC.md)
  and
  [`myBIC()`](https://jwood-iastate.github.io/flexCountReg/reference/myBIC.md).

## flexCountReg 0.1.1

CRAN release: 2026-01-19

- Initial CRAN submission.
