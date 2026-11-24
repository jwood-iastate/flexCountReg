# Changelog

## flexCountReg 0.1.2 (development)

### Total-plus-allocation models

- New
  [`countreg_alloc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_alloc.md)
  fits a total count jointly with its allocation across categories, such
  as crashes by severity: the total follows any count family (with
  distribution-parameter formulas), and the split follows an ordered
  probit, ordered logit or multinomial logit model, with generalized
  thresholds (`threshold_formula`) for the ordered links. The default
  multinomial likelihood is a proper joint likelihood of the category
  counts (likelihood-ratio tests, AIC and BIC are valid); with a Poisson
  total the categories are independent Poissons, and with an NB2 total
  and no covariates the fit equals
  [`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md).
  `likelihood = "fractional"` gives the share-weighted quasi-likelihood
  of Shabab et al. (2026, AMAR 51:100440), with sandwich standard errors
  and no AIC. Zero totals are part of the model (they contribute P(c
  = 0) and nothing about the split); `truncated = TRUE` fits a
  zero-truncated total for data that hold only units with at least one
  event, and a message points to it when the data have no zero totals.
  [`predict()`](https://rdrr.io/r/stats/predict.html) gives the expected
  counts by category, the probabilities, the total, and the expected
  total or counts given at least one event;
  [`simulate()`](https://rdrr.io/r/stats/simulate.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`glance()`](https://generics.r-lib.org/reference/glance.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html) work, and
  [`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
  can be used in every formula. Analytic gradients for all parts.
- [`countreg_alloc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_alloc.md)
  has random parameters in each equation: `rpar_total`, `rpar_alloc`
  (the propensity) and `rpar_threshold` (a common random term in the
  gaps between thresholds) add normal random deviations to the
  coefficients, drawn once per unit of `group_var` and shared by its
  rows (the years of a site), with `correlation = "none"`, `"within"`
  (each equation) or `"across"` (between the total and the allocation).
  The likelihood is simulated with randomly shifted Halton draws
  (`rpar_draws`, a stored `seed`) and has analytic gradients through the
  draws’ posterior weights. Checked against exact Gauss-Hermite
  quadrature (within 0.05 of the log-likelihood with 500 draws) and by
  recovering correlated random intercepts.
  [`predict()`](https://rdrr.io/r/stats/predict.html) gives marginal
  predictions and, with `random = "posterior"`, Empirical Bayes
  predictions for units in the estimation data;
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) draws new
  coefficients for each unit;
  [`summary()`](https://rdrr.io/r/base/summary.html) shows the standard
  deviations and correlations.
- New vignette “Total-plus-allocation models: counts by severity”.

### Time-varying effects (temporal instability)

- New formula term `tv(x, time = year)` gives `x` one effect per period,
  in the formulas of every estimator: the mean, the distribution
  parameters, zero and underreporting models,
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)’s
  `rpar_formula`, and the fixed parts of
  [`countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md),
  [`countreg_fe()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_fe.md),
  [`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md)
  and the random-effects models. `tv(1, time = year)` gives a
  time-varying intercept. Four bases: `"step"` (default; the first
  period’s effect and the change at each later period, so a zero
  coefficient means equal effects in adjacent periods), `"level"` (one
  effect per period), `"hinge"` (the cumulative “spline indicators” of
  Shabab et al. 2026, AMAR 51:100440, where a zero coefficient continues
  the earlier trend instead of keeping the effect equal) and `"smooth"`
  (the varying-coefficient smooth `s(time, by = x)`). The step, level
  and hinge bases give the same fit. The periods are stored for
  prediction; a period not in the estimation data is an error, or uses
  the nearest period with `new_time = "nearest"`.
- New
  [`temporal_effects()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_effects.md)
  gives each period’s effect with its standard error (A V A’) and
  confidence interval, and
  [`temporal_tests()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_tests.md)
  gives Wald tests that the effects are equal in all periods and of each
  change between adjacent periods. Both work for every basis, give the
  same results for the step, level and hinge bases, and accept a robust
  `vcov`.
- [`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
  does not treat the time variable of a
  [`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
  term as a predictor.
- New vignette “Time-varying effects: temporal instability”.

### Sampling weights

- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  and
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  gain `weight_type = c("frequency", "sampling")`. Frequency weights
  (the default) count each row as that many independent observations, as
  before. Sampling (survey or inverse-probability) weights make the
  weighted log-likelihood a pseudo-likelihood: standard errors default
  to the robust sandwich with meat sum((w s)(w s)’) and a G / (G - 1)
  factor, so they do not change when all weights are rescaled (with
  frequency weights, multiplying every weight by 100 shrinks them
  10-fold); the number of units is the number of rows or panels; the
  bootstrap is refused;
  [`summary()`](https://rdrr.io/r/base/summary.html) notes that LR
  tests, AIC and BIC are not valid; and
  [`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md)
  runs no LR test.

### Faster random-parameter estimation

- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  now gives the optimizer the analytic gradient of the simulated
  log-likelihood instead of a numeric one. Each panel’s score is the
  posterior-weighted average over draws of the family’s score in the log
  mean times the derivative of the log mean in each parameter. This
  covers fixed coefficients; normal, lognormal, uniform, triangular and
  gamma random parameters; correlated (Cholesky) parameters;
  heterogeneity in means and in variances; and distribution-parameter
  formulas. Poisson, NB1, NB2 and NBP have closed-form scores. The other
  families use central differences in the log mean and in each
  distribution parameter, which take two density evaluations per
  distribution parameter whatever the number of regression parameters.
  On `washington_roads` with three random parameters (9 parameters), a
  500-draw NB2 fit by BFGS takes 8.6 s instead of 41 s, and a 100-draw
  fit by BHHH takes 6 s instead of 38 s, with the same log-likelihood.
  Newton-Raphson’s Hessian is now the numeric derivative of the analytic
  gradient, which is faster and more accurate than second differences of
  the log-likelihood. `options(flexCountReg.rp_gradient = "numeric")`
  restores the numeric gradient.
- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  no longer reports convergence where the optimizer only stopped because
  the log-likelihood changed little between iterations. BHHH, NR and
  BFGSR stop on that test, and a crawling BHHH also meets it far from
  the maximum. On `washington_roads`, a 100-draw NB2 fit with three
  random parameters started from the matching
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  fit “converged” at -1086.9, 21.6 below the maximum. The fit now checks
  the Newton decrement g’(-H)^-1 g at the estimates, which is twice the
  gain a Newton step predicts. Above 0.1 it sets `converged = FALSE` and
  warns with the predicted gain. Fits with `max.iters = 0` by BFGS,
  which [`optim()`](https://rdrr.io/r/stats/optim.html) reported as
  converged, are flagged the same way. The decrement is stored as
  `model$newton_decrement`.
- Families can declare a `score` function in the internal registry; the
  NB1, NB2 and NBP scores are shared by
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  and
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md).

### Generalized additive count models (smooth terms)

- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  and
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  accept the penalized smooth terms of the mgcv package (`s()`, `te()`,
  `ti()`, `t2()`) in their formulas:
  - In
    [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md):
    the mean, the distribution-parameter formulas, the zero model and
    the underreporting model.
  - In
    [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md):
    the fixed part of the mean, `het_mean_formula`, `het_var_formula`
    and the distribution-parameter formulas.
  - Smooth bases are built once on the estimation rows and rebuilt for
    new data with the same knots.
- Fitting is by penalized maximum likelihood, finished with Newton steps
  on the penalized Hessian so that weakly penalized smooths reach the
  optimum. The smoothing parameters are selected by generalized
  Fellner-Schall updates (Wood and Fasiolo 2017), or fixed with the new
  argument `sp`.
- Standard errors are Bayesian, (H + S_lambda)^-1, by default; robust
  and cluster-robust versions are available.
- [`summary()`](https://rdrr.io/r/base/summary.html) reports each smooth
  by its effective degrees of freedom, with an approximate test (Wood
  2013), instead of its basis coefficients.
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) counts the total
  effective degrees of freedom, so
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) compares smooth and
  parametric models.
- Agreement with mgcv:
  - With fixed smoothing parameters, the estimates, log-likelihood,
    effective degrees of freedom and standard errors equal those of
    [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) for Poisson
    and NB2 (with the dispersion fixed).
  - The selected smoothing parameters were within about 5% of mgcv’s
    REML choice in the tests.
- New
  [`smooth_estimates()`](https://jwood-iastate.github.io/flexCountReg/reference/smooth_estimates.md)
  and
  [`plot_smooths()`](https://jwood-iastate.github.io/flexCountReg/reference/plot_smooths.md)
  evaluate and plot each smooth with pointwise and simultaneous
  confidence bands. Their values equal mgcv’s term predictions.
- The random-effects panel models
  [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
  [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md),
  [`poislind_re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)
  and
  [`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md)
  also accept smooth terms, in the mean and, for
  [`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md),
  in the panel-level dispersion formulas. They have the same penalized
  fitting, `sp` argument and Bayesian standard errors.
  [`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md)
  with family NM2 and
  [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md)
  are the same model, and give the same smooth.
- [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md)
  gains an analytic score for its penalized fits, which made a smooth
  fit about 10 times faster (46 s to 4 s on `washington_roads`).
- [`countreg_fe()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_fe.md)
  accepts smooth terms of variables that change within groups, with all
  three of its estimators. A smooth of a variable that is constant
  within groups is absorbed by the fixed effects and gives an error.
  - The multi-dimension Poisson estimator uses penalized IRLS.
  - The bases are built on the rows used in estimation, after the
    all-zero and single-observation groups are dropped.
  - With the same smoothing parameters, the estimates equal
    [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) with group
    dummies: to 1e-8 for one dimension, 1e-14 for two, and 1e-5 for NB2
    with the dispersion fixed.
- [`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md)
  and
  [`countreg_lc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md)
  accept smooth terms in the mean, the distribution-parameter formulas
  and the mixing (membership) formula.
  - Each component’s copy of a smooth, and each membership logit’s, has
    its own smoothing parameter.
  - The multi-start search is unpenalized; the selected, label-aligned
    solution is then refitted by penalized likelihood.
  - One component reproduces
    [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
    exactly.
  - In a simulated two-component Poisson mixture (a sine-shaped and a
    linear component), the smooths recovered both shapes: a correlation
    of 1.000 with the sine, and a slope of 0.296 for a true 0.3. The
    class shares were 0.386 and 0.614, for true 0.4 and 0.6.
  - Every estimator of the package now accepts smooth terms.
- All mgcv smooth classes can be used, among them `by` smooths (factor
  `by` gives a curve per level; numeric `by` gives a varying
  coefficient), random-effect (`bs = "re"`), factor-smooth (`bs = "fs"`)
  and Markov random field (`bs = "mrf"`, with `xt = list(nb = nb)`)
  smooths. With fixed smoothing parameters they reproduce
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html).
  - Objects in a smooth’s arguments, such as the neighborhood list, are
    taken from the formula’s environment, not from the data.
  - [`smooth_estimates()`](https://jwood-iastate.github.io/flexCountReg/reference/smooth_estimates.md)
    and
    [`plot_smooths()`](https://jwood-iastate.github.io/flexCountReg/reference/plot_smooths.md)
    show every smooth type: curves per level, factor effects with
    intervals, and the coefficient of a numeric `by`.
- New diagnostics for smooth terms:
  - [`smooth_check()`](https://jwood-iastate.github.io/flexCountReg/reference/smooth_check.md):
    a basis-dimension check (a k-index with a permutation p-value)
    computed from randomized quantile residuals.
  - [`concurvity_fcr()`](https://jwood-iastate.github.io/flexCountReg/reference/concurvity_fcr.md):
    the share of each smooth that the other terms can reproduce,
    observed and worst case.
  - [`smooth_derivatives()`](https://jwood-iastate.github.io/flexCountReg/reference/smooth_derivatives.md)
    and
    [`plot_derivatives()`](https://jwood-iastate.github.io/flexCountReg/reference/plot_derivatives.md):
    derivatives with pointwise and simultaneous bands.
  - [`plot_surface()`](https://jwood-iastate.github.io/flexCountReg/reference/plot_surface.md):
    two-variable smooths, masked away from the data.
  - [`calibration_by()`](https://jwood-iastate.github.io/flexCountReg/reference/calibration_by.md):
    observed and predicted means, probabilities of zero and variances by
    quantile groups of a covariate.
- New vignette, “Smooth effects: generalized additive count models”
  ([`vignette("smooth-effects")`](https://jwood-iastate.github.io/flexCountReg/articles/smooth-effects.md)).
- Bootstrap standard errors are not yet available with smooth terms.

### Hierarchical count models

- New
  [`countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md)
  fits count models with normally distributed random effects at any
  number of grouping levels, nested or crossed, written as in lme4:
  `(1 | site)`, `(1 + lnaadt | county)` (correlated intercepts and
  slopes), `(lnaadt || county)` (uncorrelated), `(1 | county/segment)`
  (nested) and `(1 | site) + (1 | year)` (crossed). Every count family
  can be used, including user-defined ones, with distribution-parameter
  formulas, offsets and frequency weights.
- Two integration engines: adaptive Gauss-Hermite quadrature for one
  grouping term (9 points per effect by default), and the Laplace
  approximation with sparse matrices for any number of terms.
  `engine = "auto"` chooses between them. Each covariance matrix is
  estimated through its Cholesky factor with an unconstrained diagonal,
  so a variance can be estimated as zero; such singular fits are
  reported, and the Wald standard errors of that term’s variance
  components are not. The derivatives of the log-density are analytic
  for the Poisson, NB1 and NB2 families and numerical for the others.
- Estimates, log-likelihoods and standard errors equal those of
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) (Poisson,
  Laplace and quadrature) and `glmmTMB` (Poisson, NB1 and NB2; random
  intercepts and slopes, crossed and nested effects), with
  log-likelihood differences of 1e-6 or less.
- Methods: [`summary()`](https://rdrr.io/r/base/summary.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) with `re.form`
  (condition on all, some or none of the levels), `allow.new.levels` and
  `marginal` (integrate the remaining effects out of the mean, exp(eta +
  z’ Sigma z / 2)),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) (with the
  groups’ effects or new ones),
  [`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html),
  [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html)
  (conditional modes with posterior variances),
  [`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html),
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
  [`glance()`](https://generics.r-lib.org/reference/glance.html).
  [`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html),
  [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) and
  [`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) are the nlme
  generics, so they work alongside lme4.
- [`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
  works with
  [`countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md)
  (the grouping variables are not candidates, and delta-method standard
  errors account for the fixed coefficients) and leaves out the
  regressors a fixed-effects model removed as not identified.
- New
  [`re_lrtest()`](https://jwood-iastate.github.io/flexCountReg/reference/re_lrtest.md):
  likelihood-ratio test of a random-effect term, with the
  boundary-corrected (chi-bar-squared) p-value.
- New vignette, “Hierarchical count models”
  ([`vignette("hierarchical-count-models")`](https://jwood-iastate.github.io/flexCountReg/articles/hierarchical-count-models.md)).
- Smooth terms `s()` and `t2()` in the fixed part, fitted in their
  mixed-model form (as in gamm4): the penalized part of each smooth is a
  block of normal coefficients with variance 1 / lambda, so the
  smoothing parameters are estimated by the marginal likelihood with the
  variance components. Smoothing parameters, EDF and curves equal those
  of `mgcv::gam(method = "ML")` with `bs = "re"` group effects. A
  formula with smooths and no grouping term is a GAM fitted by marginal
  likelihood.
- `re_dist = "uniform"`, `"triangular"` or `"gamma"` (one grouping term)
  gives non-normal random effects, fitted by maximum simulated
  likelihood (`engine = "msl"`, `re_draws` randomly shifted Halton draws
  per group). With `"gamma"` the multiplicative group effect is gamma
  with mean 1; a Poisson model with gamma effects reproduces
  [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md).
- `predict(history = )` estimates the effects of the groups from a
  separate set of past counts, so a forecast never uses its own
  outcomes.
- `confint(method = "profile")` gives profile-likelihood intervals for
  the random-effect standard deviations (and fixed coefficients); they
  equal glmmTMB’s.
- `cre` and `cre_group` add the group means of time-varying covariates
  (correlated random effects), tested with
  [`cre_test()`](https://jwood-iastate.github.io/flexCountReg/reference/cre_test.md).
- Not yet available in
  [`countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md):
  non-normal random effects at more than one level, smooth terms in the
  distribution-parameter formulas, and `te()`/`ti()` smooths (use
  `t2()`).

### Fixed-effects count models

- New
  [`countreg_fe()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_fe.md)
  estimates count models with fixed effects for one or more grouping
  variables (`fixef = ~ site` or `~ site + year`). It chooses among
  three estimators:
  - Poisson with one dimension: the conditional likelihood.
  - Poisson with several dimensions: IRLS with alternating projections,
    as in
    [`fixest::fepois()`](https://lrberge.github.io/fixest/reference/feglm.html).
  - Any other family, including families from
    [`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md),
    with one dimension: the dummy-variable likelihood with the group
    intercepts concentrated out.
- Groups with all-zero counts or a single observation are dropped and
  reported, and regressors without variation within groups are removed
  with a warning.
- Standard errors are cluster-robust by the first fixed-effect variable
  by default. The log-likelihood is that of the dummy-variable model,
  with degrees of freedom that count the fixed effects.
- `jackknife = TRUE` applies the split-panel jackknife of Dhaene and
  Jochmans (2015). In simulations with 200 units it removed most of the
  NB2 dispersion bias with 5 or 10 periods, but did not improve the
  slopes, and it is unreliable below 4 periods, where a warning is
  given.
- [`predict()`](https://rdrr.io/r/stats/predict.html) uses each known
  group’s effect. Groups not in the estimation data give `NA` unless
  `newgroup = "mean_fe"`.
- Estimates, log-likelihoods and dropped groups equal those of
  [`fixest::fepois()`](https://lrberge.github.io/fixest/reference/feglm.html)
  (one and two ways) and
  [`fixest::fenegbin()`](https://lrberge.github.io/fixest/reference/femlm.html).
  The cluster-robust standard errors equal fixest’s with
  `ssc(adj = FALSE)`.
- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md),
  [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)
  and
  [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md)
  gain `cre = ~ x1 + x2` for correlated random effects (Mundlak 1978).
  The units’ means of the listed time-varying regressors are added to
  the mean model, computed on the estimation rows and stored with the
  fit. Predictions for known units use their training means; new units
  use the means of their own rows.
- New
  [`cre_test()`](https://jwood-iastate.github.io/flexCountReg/reference/cre_test.md),
  a Wald test that the coefficients of those means are zero: a robust
  test of random against fixed effects. In simulations with unit effects
  correlated with the regressor, the plain random-effects slope was
  biased (0.63 for a true 0.5), the slope with `cre` was 0.505, equal to
  fixed effects, and the test rejected in every replication. With
  uncorrelated effects, it rejected in 7% of replications at the 5%
  level.
- New vignette, “Fixed-effects count models”
  ([`vignette("fixed-effects-count-models")`](https://jwood-iastate.github.io/flexCountReg/articles/fixed-effects-count-models.md)),
  with a section on correlated random effects.

### Cluster-robust standard errors and the sandwich package

- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  and
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  gain `stderr = "cluster"` with `cluster = ~ corridor` (or column
  names, or a vector of ids). Two or more variables
  (`~ corridor + year`) give multiway clustering (Cameron, Gelbach and
  Miller 2011). The scores of the independent units (rows, or panels)
  are summed within clusters, with the factor G / (G - 1). Clusters may
  be coarser than the panels (sites within corridors) but may not split
  a panel. Giving `cluster` alone selects cluster-robust standard
  errors. With frequency weights, a row’s copies are in its cluster.
  [`summary()`](https://rdrr.io/r/base/summary.html) reports the number
  of clusters.
- [`vcov()`](https://rdrr.io/r/stats/vcov.html) gains
  `type = c("fit", "hessian", "robust", "cluster")` and `cluster`, which
  recompute the covariance matrix from the stored scores without
  refitting. This works for
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md),
  [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
  [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md)
  and
  [`poislind_re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)
  fits. With `keep_data = FALSE`, the clusters are given as a vector.
- New `estfun()` and `bread()` methods for the sandwich package, so
  [`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html),
  `vcovHC()`, `sandwich()`,
  [`lmtest::coeftest()`](https://rdrr.io/pkg/lmtest/man/coeftest.html)
  and tools built on them work with flexCountReg fits. One-way and
  two-way cluster covariances equal `sandwich::vcovCL(type = "HC0")` to
  1e-10.

### User-defined count distributions

- New
  [`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md)
  and
  [`count_param()`](https://jwood-iastate.github.io/flexCountReg/reference/count_param.md)
  define a count distribution from its log PMF, with up to two
  parameters, each with a log, logit or identity link. The family object
  can be passed as `family` to
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md),
  [`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md)
  and
  [`countreg_lc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md).
  It works with distribution-parameter formulas, zero-inflated and
  hurdle models, predictions, quantile residuals, rootograms,
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) and robust or
  cluster-robust standard errors.
  - Optional `mean`, `variance`, `cdf`, `rand` and `valid` functions can
    be given. Without them, the moments and the CDF come from summing
    the PMF, and random draws from inverting the CDF.
  - A user NB2 reproduces the built-in NB2 in fixed, random-parameter,
    hurdle and mixture fits. A user zero-inflated NB2 reproduces
    `zero_formula = ~ 1`.
  - Each family is registered under its name and a hash of its
    definition, so two families with the same name do not collide.
  - The fit stores the family, so a model saved with
    [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) works in a new R
    session.
  - Underreporting is not available for user families.
- New
  [`check_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/check_count_family.md)
  checks a family over a grid of means and parameter values. It tests
  that the PMF sums to one, that `mean`, `variance`, `cdf` and `rand`
  agree with the PMF, that the function is vectorized, and that missing
  or invalid values cause no errors.
- New vignette, “Fitting your own count distribution”
  ([`vignette("user-defined-count-distributions")`](https://jwood-iastate.github.io/flexCountReg/articles/user-defined-count-distributions.md)),
  with a one-inflated negative binomial as the worked example.

### Interactive app

- New
  [`launch_app()`](https://jwood-iastate.github.io/flexCountReg/reference/launch_app.md)
  opens an interactive interface to flexCountReg in a web browser, with
  a tab for each step: load a bundled example or a CSV file and create
  new variables (logs, powers, standardized values, indicators,
  products, ratios, categories, factors); specify and estimate count
  regression, random-parameter, random-effects panel (RENB,
  Poisson-Lindley, Poisson, negative multinomial), finite-mixture and
  latent-class models (any family, exposure or log offset, frequency
  weights, distribution-parameter, membership and shared terms, panels);
  inspect the summary, coefficients, numerical checks, class sizes and
  starts; diagnose fits (rootogram, randomized quantile residuals, CURE
  plot); compare them (with the reasons when they are not comparable);
  predict for an editable table of scenarios or for uploaded new data,
  with confidence, prediction and tolerance intervals for count
  regressions; and download the analysis. The Diagnose and Predict tabs
  work on any completed model (a selector, which moves to each new fit),
  plots can be saved as JPEG files at a chosen size and resolution, and
  new data get the variables created on the Data tab that the model
  needs, with the means, standard deviations and categories of the
  training data.
- Everything the app runs (data loading, new variables, fits,
  diagnostics, predictions, intervals and comparisons) is compiled into
  R code that calls the package’s functions, evaluated from that code
  and appended to a session script, shown in the app and downloaded as
  `analysis.R`; the script reproduces the whole session without the app
  (the download also holds the settings of each fit as JSON and a
  manifest with versions and warnings). Column names and typed values
  enter the code only as names, numbers or strings, never as code. An
  action that fails reports its error and leaves the session running.
- Models are estimated in background R processes (callr), so the
  interface stays responsive and fits can be canceled; a slow earlier
  fit never replaces a newer one. Predictions and intervals also run in
  the background: an interval for a model whose covariance matrix is not
  positive definite uses a parametric bootstrap (200 refits), which can
  take many minutes, and the app says so when it starts and offers to
  cancel it. Data, fits and jobs belong to one browser session.
- Models built in the app can include interactions between pairs of
  predictors (entered as `a:b`; main effects stay in the model), and the
  Fit page has a Marginal effects tab running
  [`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
  (average marginal effects or elasticities; pseudo-elasticities or
  discrete changes for indicators; delta-method or bootstrap standard
  errors) in the background, with the code added to the session script
  and the table downloadable as CSV.
- Exports are reproducible: each uploaded file is stored at its own path
  (`data/data_1.csv`, `data/newdata_1.csv`, …) with its MD5 checksum in
  the manifest; warnings and messages of every result are shown beside
  it and listed in the manifest with the uncertainty method, accepted
  bootstrap refits and seed; and stochastic requests (bootstrap
  intervals, tolerance intervals, bootstrap marginal effects) fix their
  seed when they are made, so the exported script replays them exactly.
- The app runs the models added in this release: hierarchical models
  ([`countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md):
  a grouping column, random slopes, correlated or uncorrelated, a second
  level crossed with or nested in the first, and normal, gamma, uniform
  or triangular effects) and fixed-effects models
  ([`countreg_fe()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_fe.md):
  one or two dimensions, the split-panel jackknife, cluster-robust
  standard errors by default). Any model can have smooth (nonlinear)
  effects `s(x)` of numeric predictors, shown on the Diagnose page
  (“Smooth effects”,
  [`plot_smooths()`](https://jwood-iastate.github.io/flexCountReg/reference/plot_smooths.md));
  random-parameter, hierarchical and random-effects models can add group
  means of predictors (correlated random effects, `cre`); count
  regression and random-parameter models offer cluster-robust standard
  errors (with a cluster column) and sampling weights. Predictions for
  hierarchical and fixed-effects models take the group of each scenario,
  and settings files from earlier versions of the app still load.
- The app’s packages (shiny, bslib, DT, callr, jsonlite, zip) are
  suggested, not imported: the statistical functions do not need them,
  and
  [`launch_app()`](https://jwood-iastate.github.io/flexCountReg/reference/launch_app.md)
  offers to install any that are missing.
  [`app_object()`](https://jwood-iastate.github.io/flexCountReg/reference/launch_app.md)
  returns the app for deployment or testing.

### Finite mixtures and latent classes

- New
  [`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md)
  estimates finite mixtures of count distributions and finite-mixture
  count regressions;
  [`countreg_lc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md)
  is the same model under the latent-class vocabulary (`classes`,
  `membership_formula`). Both normalize their arguments into one
  specification and call one engine, so equivalent calls give identical
  fits. Supported specifications:
  - distribution-only mixtures (`y ~ 1`), including fits to exact-count
    frequency tables through `weights`;
  - mixture regressions with component-specific coefficients, or with
    whole terms shared by all components (`shared = list(mean = "x")`)
    and shared distribution-parameter blocks
    (`shared = list(dis1 = TRUE)`);
  - constant mixing proportions, covariate-dependent membership
    probabilities (a multinomial logit, `mixing_formula = ~ z`, last
    component as reference), or known fixed proportions
    (`fixed_mixing`), which add no parameters;
  - any count family of
    [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
    in the components, the same in all of them or different ones
    (e.g. `family = c(base = "POISSON", tail = "NB2")`), with
    distribution-parameter formulas;
  - panel mixtures (`group_var`): one class per panel, with the class
    probability entering once per panel (product of the rows’
    probabilities before the sum over classes);
  - offsets and frequency weights (outside the logarithm of the mixture;
    constant within panels).
- Estimation uses several reproducible starts (`n_starts`, `seed`), EM
  iterations with weighted M-steps (analytic scores for Poisson, NB1 and
  NB2 components), and direct maximization of the observed-data
  log-likelihood to polish the best solutions
  ([`mixture_control()`](https://jwood-iastate.github.io/flexCountReg/reference/mixture_control.md),
  [`lc_control()`](https://jwood-iastate.github.io/flexCountReg/reference/mixture_control.md)).
  Every start is recorded (`fit$model$starts`), including whether
  several starts reached the selected solution. Components of the same
  family are numbered by increasing average expected count.
- Standard errors come from the observed information of the full
  parameter vector (`stderr = "normal"`), unit-level sandwich scores
  (`"robust"`), or a bootstrap of independent units that refits the full
  mixture with several starts and aligns component labels within
  families (`"boot"`). The fit records the gradient, the definiteness
  and conditioning of the observed information, and ambiguous bootstrap
  label matches.
- Methods: [`summary()`](https://rdrr.io/r/base/summary.html)
  (coefficients grouped by component, shared parameters and mixing
  model; prior and posterior class sizes; entropy, relative entropy and
  ICL; start stability; numerical checks),
  [`print()`](https://rdrr.io/r/base/print.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) (free parameters
  only, independent units for BIC),
  [`predict()`](https://rdrr.io/r/stats/predict.html) with explicit
  targets (`type = "response"`, `"component"`, `"variance"`, `"prob"`,
  `"link"`; `conditioning = "prior"`, `"history"` for panel forecasts
  from earlier outcomes only, or the retrospective `"posterior"`),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) (response,
  Pearson with the full mixture variance, and randomized quantile
  residuals from the mixture distribution function),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) (one class per
  row or panel),
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`glance()`](https://generics.r-lib.org/reference/glance.html) and
  [`augment()`](https://generics.r-lib.org/reference/augment.html). New
  [`component_probs()`](https://jwood-iastate.github.io/flexCountReg/reference/component_probs.md)
  /
  [`class_probs()`](https://jwood-iastate.github.io/flexCountReg/reference/component_probs.md)
  return prior, history-conditioned or posterior class probabilities.
- [`rootogram()`](https://jwood-iastate.github.io/flexCountReg/reference/rootogram.md)
  and
  [`cureplot()`](https://jwood-iastate.github.io/flexCountReg/reference/cureplot.md)
  accept mixture models (prior mixture probabilities and means), and
  [`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
  includes variables that enter only the membership model (delta-method
  standard errors with the full covariance matrix).
- Not available yet, with explanatory errors: confidence, prediction and
  tolerance intervals for mixture predictions; bootstrap standard errors
  of marginal effects; likelihood-ratio tests involving mixtures
  ([`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md)
  refuses them, because the chi-square reference distribution does not
  hold when the number of components changes). Zero-inflated, hurdle and
  underreporting components inside mixtures and component-specific
  formulas are not implemented. Mixtures of families other than Poisson,
  NB1 and NB2 are labeled experimental: their component distributions
  are those of
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  (one-component fits agree with it), but recovery of such mixtures has
  not been studied.
- New vignettes: “Finite mixtures of count distributions” and
  “Latent-class and finite-mixture count regression”, and a recovery and
  coverage study, `inst/validation/finite-mixtures/` (its `CORES`
  setting runs the replicates in parallel with the same results). In 200
  replicates of each of seven designs (well and weakly separated Poisson
  mixtures, a latent-class NB2 regression with membership covariates and
  a shared dispersion, a panel latent-class model, and a Poisson + NB2
  mixture), every fit converged; Wald intervals covered 0.89 to 0.985
  (Monte Carlo standard errors about 0.015), and the label-invariant
  mixture mean was nearly unbiased in every design. With weakly
  separated components and n = 300 the low component’s intercept was
  biased (-0.18) and its standard errors unstable, and 40% of the fits
  were reached by only one of five starts: use more starts and treat
  component-specific inference with caution in such data.
- A warning (and a note in
  [`summary()`](https://rdrr.io/r/base/summary.html)) flags components
  of the same family with nearly identical parameters: the fit is then
  effectively one component fewer and its shares are arbitrary.

### New count distributions

- Eight new count distributions, each with density, distribution,
  quantile and random-generation functions following base R conventions
  (recycling, `NA`/`NaN`, `NaN` with a warning for invalid parameters,
  `lower.tail`, `log`/`log.p`):
  - discrete Weibull
    ([`ddweib()`](https://jwood-iastate.github.io/flexCountReg/reference/DiscreteWeibull.md),
    [`pdweib()`](https://jwood-iastate.github.io/flexCountReg/reference/DiscreteWeibull.md),
    [`qdweib()`](https://jwood-iastate.github.io/flexCountReg/reference/DiscreteWeibull.md),
    [`rdweib()`](https://jwood-iastate.github.io/flexCountReg/reference/DiscreteWeibull.md));
  - double Poisson, exactly normalized
    ([`ddpo()`](https://jwood-iastate.github.io/flexCountReg/reference/DoublePoisson.md),
    …);
  - Delaporte in a mean parameterization
    ([`ddelaporte()`](https://jwood-iastate.github.io/flexCountReg/reference/Delaporte.md),
    …);
  - generalized Hermite of order `m`
    ([`dghermite()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedHermite.md),
    …);
  - hyper-Poisson
    ([`dhpois()`](https://jwood-iastate.github.io/flexCountReg/reference/HyperPoisson.md),
    …);
  - gamma renewal count
    ([`dgcount()`](https://jwood-iastate.github.io/flexCountReg/reference/GammaCount.md),
    …);
  - Polya-Aeppli, the Poisson-Tweedie distribution with power 1.5
    ([`dpaeppli()`](https://jwood-iastate.github.io/flexCountReg/reference/PolyaAeppli.md),
    …).
- Probabilities are computed on the log scale, and upper tails
  (`lower.tail = FALSE`) are computed directly rather than as one minus
  the distribution function, so tiny tail probabilities keep their
  relative accuracy. Quantiles are the smallest count whose distribution
  function reaches the requested probability.
- Negative binomial-generalized exponential mixture
  ([`dnbge()`](https://jwood-iastate.github.io/flexCountReg/reference/NBGE.md),
  [`pnbge()`](https://jwood-iastate.github.io/flexCountReg/reference/NBGE.md),
  [`qnbge()`](https://jwood-iastate.github.io/flexCountReg/reference/NBGE.md),
  [`rnbge()`](https://jwood-iastate.github.io/flexCountReg/reference/NBGE.md)):
  NB2 with a mean-one generalized exponential random effect, computed by
  the adaptive quadrature used for PGE.

### New regression families

- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  and
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  gain the families `"DW"`, `"DPO"`, `"DELAPORTE"`, `"HERMITE"` (order
  2), `"HERMITE3"`, `"HERMITE4"`, `"HYPERPOIS"`, `"GAMMACOUNT"`, `"PA"`
  (Poisson-Tweedie with power 1.5) and `"NBGE"`. All are parameterized
  by the mean, so coefficients are effects on the expected count. For
  DW, DPO, HYPERPOIS and GAMMACOUNT, whose means have no closed form,
  the natural parameter of each observation is found by solving the mean
  equation in compiled code; these families are slower to fit and have a
  numerical domain (shape in \[1e-6, 1e6\], mean up to 1e9).
- Shares in \[0, 1\] (Delaporte and Hermite `rho`) use a logit link,
  labeled `logit(rho)`.
- Prediction, Pearson and quantile residuals, simulation, rootograms and
  summaries work for all new families. Likelihood-ratio tests against
  Poisson are reported for DPO, HYPERPOIS and GAMMACOUNT (interior
  restriction). Underreporting (binomial thinning) is available for
  DELAPORTE, NBGE, HERMITE and PA.
- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  uses family-specific starting values for the distribution parameters
  when the registry provides them (the new families start at their
  Poisson or geometric limit rather than at log-scale -2).

### Confidence, prediction and tolerance intervals

- [`predict()`](https://rdrr.io/r/stats/predict.html) gains
  `interval = "confidence"` (for the expected count) and
  `interval = "prediction"` (for one new count), with `level`, for
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  models. Estimation uncertainty is represented by parameter vectors
  from the asymptotic normal distribution of the estimates, evaluated at
  deterministic Halton points (reproducible), or by a parametric
  bootstrap (`uncertainty = "bootstrap"`), which is also used when the
  covariance matrix is not positive definite. Prediction intervals are
  integer quantiles of the predictive distribution (the count PMF
  averaged over the parameter vectors).
- New
  [`tolerance_interval()`](https://jwood-iastate.github.io/flexCountReg/reference/tolerance_interval.md):
  content tolerance intervals (content and confidence levels; two-sided,
  upper or lower; pointwise or simultaneous) for
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  models. The default joint method returns the narrowest integer
  interval whose content is at least the requested fraction at every
  evaluated point of the confidence set; the envelope method uses
  per-tail worst cases; intercept-only Poisson models use the exact
  common-rate construction. Zero-inflated, hurdle and underreporting
  components are included. A reproducible coverage study is in
  `inst/validation/intervals/`.

### Starting values for the dispersion

- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  now starts the dispersion parameter of NB1, NB2, PIG1, PIG2 and PIG at
  a moment estimate from the Poisson starting fit (the auxiliary
  regression of Cameron and Trivedi) instead of `alpha = 1`. With the
  old start, quasi-Newton optimizers could step onto the flat region
  `alpha -> 0` and stop there with a falsely “converged” Poisson fit:
  NB2 with `method = "BFGS"` gave log-likelihood -228.48 instead of
  -195.56 on `warpbreaks`, and -182.29 instead of -180.11 on
  `InsectSprays`. Results now match
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html). Fits
  that already reached the maximum are unchanged.

### Standard errors at flat likelihood directions

- When the Hessian is singular, standard errors still use a generalized
  inverse, but parameters involved in a direction where the likelihood
  is flat now get `NA` variances and standard errors instead of zero (a
  zero standard error read as perfect precision, e.g. the random-effect
  parameter of a Poisson-Lindley random-effects model at its boundary).
  The warning names those parameters.
- When the Hessian is invertible but not negative definite (for example
  an NB2 dispersion drifting toward zero on Poisson-like data), the
  negative variances already gave `NA` standard errors; a warning now
  names those parameters instead of leaving them silently `NA`.

### Vignettes

- Seventeen vignettes, grouped on the package website as Start here,
  Diagnose and interpret, Model structure, Finite mixtures and latent
  classes, Predict and validate, and Advanced. New (besides the two
  finite-mixture vignettes above): getting started; data, exposure and
  weights; diagnosing count models; predictions, marginal effects and
  elasticities; modeling distribution parameters; zeros, hurdles and
  sampling; underreporting and identification; predictive validation;
  prediction for panel data; numerical reliability; simulation studies;
  and adding a count distribution (for contributors). The three earlier
  vignettes are revised.

- `inst/validation/simulation-studies/count-regression-study.R` is a
  parameterized simulation study (Poisson, NB2 and a negative
  multinomial panel) with one random-number stream per replicate; its
  summary and settings are saved alongside it. The per-replicate results
  of this study and of the interval coverage study are kept in the
  source repository but left out of the package build. \## New data

- `abalone_rings`: shell-ring counts and measurements of 4,177 abalone
  from the UCI Machine Learning Repository (CC BY 4.0), a
  non-transportation example with conditionally underdispersed counts.
  The preparation script and source hashes are in `data-raw/` of the
  source repository; data terms are in the installed `DATA-LICENSES.md`.

- `stream_salamanders`: repeated counts of seven stream salamander
  groups at 23 sites (Price et al., 2016; Dryad, CC0), prepared from the
  original deposit. Zero-heavy, overdispersed and grouped by site.

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
- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  and
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  compute the final Hessian, and so every Hessian-based, robust and
  cluster-robust standard error, adaptively. maxLik’s default step is
  1e-6, where rounding error is about 0.2 per Hessian entry for a
  log-likelihood near 1000. That made the standard errors of weakly
  determined parameters wrong. In the README’s generalized NBP model on
  `washington_roads` they were up to 25% too large (ln(p): 0.229, now
  0.183).
  - The Hessian is computed at steps h and h/2 (h = 1e-3 times the size
    of each estimate) and combined by Richardson extrapolation. Larger
    steps (for log-likelihoods computed with quadrature, series or
    Bessel functions) and smaller ones are tried when the two disagree.
  - Families with an analytic gradient (Poisson, NB1, NB2, and
    [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
    in general) difference the gradient; the others use second
    differences of the log-likelihood.
  - The fit stores `model$hessian_accuracy`, the largest relative change
    of a standard error between the two steps, and warns when it exceeds
    1%.
  - Checked against independent Hessians, the standard errors agree to
    about 6e-5 or better for the Poisson, NB1, NB2, NBP, PIG, PIG1, PLN,
    GP1 and DW families. The Sichel family agrees to 0.15%, for its
    weakly identified gamma.
  - BHHH fits keep the outer-product Hessian of the scores.
- Hessian-based, robust and cluster-robust covariances use the symmetric
  part of the numeric Hessian. In poorly conditioned fits the inverse
  amplified the Hessian’s small asymmetry. For example, a
  [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)
  fit by Nelder-Mead on `washington_roads` had a negative robust
  variance for ln(a-1); it is now 0.036. In well-conditioned fits the
  change is negligible.
- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  gives each unit (panel, or row without `group_var`) its own simulation
  draws by default: new argument `draws = c("per_unit", "shared")`. Unit
  u uses draws (u - 1) R + 1 to u R of one Halton sequence, after its
  first 10 points (Train 2009, Sec. 9.3.3). With shared draws, as
  before, draw r is the same for every panel, so simulation errors are
  correlated across panels and do not average out as panels are added.
  On `washington_roads` (NB2, random `speed50` and `lnlength`, 12
  scrambling seeds), per-unit draws cut the error of the
  random-parameter means and standard deviations by a factor of 2 to 5
  at 50 and 200 draws. The intercept and the log-likelihood varied more
  across seeds than with shared draws. Empirical Bayes predictions
  regenerate each estimation panel’s own draws (by panel id), and new
  panels get new draws. `draws = "shared"` reproduces earlier fits.
- Scrambled Halton draws (`scrambled = TRUE` in
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md),
  and in interval simulations) now depend on their seed. `randtoolbox`
  scrambles with its own generator, which
  [`set.seed()`](https://rdrr.io/r/base/Random.html) does not reach. So
  every “scrambled” fit used the same scrambling (randtoolbox’s default,
  seed 0), whatever the stored `draw_seed`, and different seeds gave
  identical results. Fits made before this change still predict with the
  draws they were estimated with.
- `countreg_rp(scrambled = TRUE)` now randomizes the Halton draws with a
  random shift (Cranley-Patterson rotation): each unit’s block of plain
  Halton points is shifted by its own uniform amount per dimension,
  modulo 1, with the shifts generated from the fit’s `draw_seed`. It
  used randtoolbox’s scrambled (“mixed”) Halton sequence, whose
  consecutive per-unit blocks biased the simulated log-likelihood, and
  more draws did not remove the bias. For a Poisson random intercept on
  `washington_roads` (exact maximum -1069.365 by adaptive quadrature),
  scrambled fits gave -1069.406 at 500 draws and -1069.054 at 2000; with
  random shifts the error at fixed parameters averages -0.02 at 500
  draws and -0.008 at 2000 (standard deviation over seeds 0.14 and
  0.027). Scrambled fits therefore give different estimates than before.
  New fits have `draw_version` 4; fits with an earlier `draw_version`
  still predict with the “mixed” draws they were estimated with.
  Prediction intervals with `scrambled = TRUE` keep the “mixed”
  sequence, which uses one point set and no per-unit blocks.
- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  now uses `method = "BFGS"` by default (it was `"BHHH"`). With the
  analytic gradient, BFGS reached the maximum in each of six test
  models. BHHH failed on two of the six. On a 500-draw NB2 model with
  three random parameters it was still at -1093 after 200 iterations,
  while BFGS reached -1065.5 in about 10 s. On the same model with 100
  draws it reported convergence at the Poisson boundary (ln(alpha) =
  -167), 0.08 below the maximum, where no gradient method can recover.
  BHHH remains available and is often faster on simple models. Code that
  relies on the default can give different estimates; pass
  `method = "BHHH"` to reproduce earlier fits.
- Fitted objects no longer keep the optimizer’s log-likelihood function
  (`maxLik::objectiveFn(fit$model)` is unavailable): its closure held
  the data, design matrices and draws, so fits were several times larger
  than the data (30 MB for 4.3 MB of data).
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
- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md):
  heterogeneity in means now follows Wood and Gayah (2025, Eq. 2). Each
  random parameter has its own coefficient for each variable in
  `het_mean_formula` (named `HetMean:<parameter>:<variable>`), and the
  effect is additive (on the log scale for lognormal parameters;
  multiplicative for gamma parameters, whose mean must stay positive).
  Previously a single coefficient vector multiplied every random
  parameter’s mean. Estimates from models with `het_mean_formula`
  change; saved fits still predict with their original specification.
- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md):
  uniform random parameters now span mu +/- St.Dev (draws `2u - 1`),
  matching the paper and the Exact prediction. Previously the draws
  spanned mu +/- St.Dev/2, so estimated uniform St.Dev values are now
  half as large. [`predict()`](https://rdrr.io/r/stats/predict.html)
  converts saved uniform fits automatically (with a one-time message).
- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md):
  `rpardists` is matched to the random parameters by name (not position)
  and validated. An intercept in both `formula` and `rpar_formula` (an
  unidentified model) is an error; the intercept implied by
  `rpar_formula` is not a random parameter unless the fixed formula has
  no intercept. An invalid `panel_id`, `family` or `method` is an error.
- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md):
  `weights` now weight the panel log-likelihood (frequency weights for
  whole panels) instead of each row’s probability.
- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md):
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

- The Poisson-inverse-gamma density
  ([`dpinvgamma()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGamma.md)
  and the PIG family) could take minutes, or effectively hang, when the
  optimizer tried a small `eta`. Its log Bessel function took one
  recurrence step per unit of order, and the order is about 1/eta.
  Orders above 100 now use the uniform asymptotic expansion (accurate to
  about 1e-10 there), so the cost no longer grows with the order. A
  bootstrap refit that took 12 seconds now takes 0.2 seconds.

- [`dpinvgaus()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGaussian.md),
  [`ppinvgaus()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGaussian.md)
  and the PIG1/PIG2 families: P(0) was computed as
  `exp(mu / eta * (1 - sqrt(1 + 2 * eta)))`, which cancels as eta -\> 0
  (at eta = 1e-16 the Type 1 PMF summed to 20, and
  `ppinvgaus(5, 3, 1e-16)` was 18.4), and fits could exploit it to
  report a positive log-likelihood. It is now
  `exp(-2 * mu / (1 + sqrt(1 + 2 * eta)))`, which tends to the Poisson
  value.

- [`dsichel()`](https://jwood-iastate.github.io/flexCountReg/reference/SichelDistribution.md)
  returned -Inf where the log PMF is finite (for example
  `dsichel(500, 500, 1, 1, log = TRUE)`), and
  [`qsichel()`](https://jwood-iastate.github.io/flexCountReg/reference/SichelDistribution.md)
  (and so
  [`rsichel()`](https://jwood-iastate.github.io/flexCountReg/reference/SichelDistribution.md))
  stopped at 10,000 for heavy tails. Bessel functions now come from a
  log-scale order recurrence, and
  [`psichel()`](https://jwood-iastate.github.io/flexCountReg/reference/SichelDistribution.md)
  and
  [`qsichel()`](https://jwood-iastate.github.io/flexCountReg/reference/SichelDistribution.md)
  evaluate the PMF in one recurrence sweep and search the support in
  doubling blocks.

- The inverse-Gaussian draws behind
  [`rpinvgaus()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGaussian.md)
  cancelled to 0 for large dispersions; the root is now computed as m /
  (1 + t + sqrt(t (t + 2))).

- [`dpoisweibull()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonWeibull.md)
  and friends with `mean_value` and `sd_value`: the shape was searched
  only in \[0.1, 5\], so small coefficients of variation were silently
  given shape 5 (SD 0.229 instead of 0.01 at mean 1). The shape is now
  solved from the coefficient of variation with an expanding bracket.

- `regCompTable(coefs = FALSE)` printed the coefficients anyway, and an
  unnamed list of models gave a table without model columns (they are
  now named “Model 1”, “Model 2”, …).

- [`cureplot()`](https://jwood-iastate.github.io/flexCountReg/reference/cureplot.md)
  no longer prints the cureplots messages “Covariate: indvar_values” and
  “CURE data frame was provided…”, which named the wrapper’s internal
  objects; other messages are kept. The documentation of `n_resamples`
  now says that the red reference limits are drawn whatever its value
  (they were described as “no bands” at 0).

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

- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
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
  [`rcom()`](https://jwood-iastate.github.io/flexCountReg/reference/COMDistribution.md))
  accept vector parameters, as needed for simulation from fitted models.
  Mixture families draw the mixing variable per observation.
- q-functions no longer loop forever at p = 1 (they return `Inf`),
  return `NaN` for p outside \[0, 1\], and no longer recompute the CDF
  at every step.
- `dplind(lambda = )` uses `lambda` (it previously ignored it).
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
- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  uses the shared data preparation (formula offsets, missing values,
  numeric weights), starts gamma random parameters from a positive mean,
  averages panel likelihoods over draws on the log scale, and reports
  convergence.
- `predict(method = "Exact")` for
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  includes heterogeneity in means and variances.

### Further bug fixes and changed results

Changed results are marked **(changed)**.

- **Weights for the panel estimators.**
  [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)
  and
  [`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)
  accept `weights`, as
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  and
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  do. In panel models these are frequency weights of the panels
  (constant within a panel): a panel with weight w counts as w
  independent panels in the log-likelihood, the robust covariance, BIC
  and the bootstrap. Fits with weights equal fits to the expanded
  panels.

- **Frequency-weight bootstrap (changed).** With `weights`, bootstrap
  standard errors now draw new frequencies from Multinomial(sum(w), w /
  sum(w)) and refit with them in place of the weights. This equals the
  bootstrap of the expanded observations (or panels). Before, rows were
  resampled with their original weights, which is a different and much
  more variable sampling scheme. Bootstrap standard errors need
  whole-number weights; use `stderr = "robust"` for other weights.

- **Bootstrap acceptance.** A replicate counts only if its refit reports
  convergence (unknown convergence is a failure), has exactly the
  expected coefficients and finite estimates. Failures are stored by
  reason in `model$boot_failures`.

- **Likelihood-ratio tests.** Fits store evidence of their analysis
  sample (rows, response, weights, offsets, panels and the variables
  they use). An LR test between two models requires the same sample
  (shared variables must agree row by row), reported convergence for
  both, nested designs (checked by column space, so equivalent bases
  such as `x` and `2 * x` are accepted) and matching random-parameter
  distributions and correlation. A larger model whose log-likelihood is
  materially below the restricted one’s (LR below -max(0.01, 2e-6
  \|logLik\|)) gives no test, with a reason in `LR_reason`, instead of
  LR = 0. The intercept-only base model of
  [`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md)
  is not nested in a model without an intercept. Objects from older
  versions get no LR test between fitted models (refit them).

- **Mixture families (changed).** The PLN, PGE, PW and PLL probabilities
  are computed on the log scale by deterministic adaptive quadrature
  over the log of the mixing variable (`engine = "quadrature"`, the new
  default of the d/p/q functions,
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  and
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)).
  Averaging over 1500 Halton draws was badly biased in the tails: for
  PGE with shape 1 (geometric with mean 1) and y = 1000 it gave log P =
  -3847 instead of -693.8. The new values agree with exact special cases
  (PGE with shape 1 and PW with alpha 1 are geometric; P(Y = 0) of PGE
  is a beta function) and with independent quadrature to about 1e-9, for
  counts up to 1e5. Fits are faster too. `engine = "halton"` keeps the
  earlier method, and objects without a stored engine predict as before.

- **Moment existence.** Random-parameter predictions decide whether
  E\[exp(x beta)\] exists from the parameters alone, so `Simulated`
  predictions with lognormal parameters no longer need lamW. `Exact`
  uses log-scale MGFs and reports an overflow of a finite expectation as
  such. For lognormal parameters with x \< 0, `Exact` uses a Laplace
  approximation (documented as approximate).

- **Distribution contracts.** Invalid parameters give `NaN` for that
  element, with one warning (base R’s convention), and missing
  parameters `NA`, at every `q` including `q < 0` and `q = Inf` (before,
  those endpoints returned 0 or 1).
  [`dgwar()`](https://jwood-iastate.github.io/flexCountReg/reference/Generalized-Waring.md)/[`pgwar()`](https://jwood-iastate.github.io/flexCountReg/reference/Generalized-Waring.md)
  no longer stop for an invalid element.
  [`plindley()`](https://jwood-iastate.github.io/flexCountReg/reference/Lindley.md)
  computes both tails directly.

- **COM controls.** `find_lambda_cpp()` honors its tolerance (relative
  error of the mean) and iteration limit and reports non-convergence
  instead of accepting a collapsed bracket; `find_lambda_diag_cpp()`
  gives the residual, iterations and status. Unused controls of the
  internal moment functions were removed.

- **[`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)**:
  the log-likelihood is computed as three stable differences of
  log-gamma functions. For large `a` and `b`, terms of size ~1e200
  canceled and the optimizer could report a positive log-likelihood of
  ~1e203 **(changed**, fits that hit this region**)**.

- **[`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)
  standard errors (changed).**
  [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)
  now has analytic scores, so its Hessian is the derivative of an exact
  gradient. Before, it was a second finite difference of the
  log-likelihood, which was inaccurate (it could be indefinite and give
  `NaN` standard errors or delta-method marginal-effect standard errors,
  notably for `ln(a-1)` and `ln(b)`). Estimates are unchanged; standard
  errors can change.

- `rpardists` values are checked before a single value is expanded.

- Documentation: the PGE equations use one rate convention and explain
  why the rate cancels; fitted objects carry a schema version
  (`model$fcr_schema`).

- PGE and PW: with the mean parameterization, the scale of the mixing
  distribution cancels out of the likelihood, so it was not identified
  (the PMF is identical for every scale). These families now estimate
  only the shape parameter (one fewer parameter, so AIC and BIC change),
  and `dis_param_formula_2` is an error for them. Predictions from fits
  made with earlier versions are unchanged. The
  [`dpge()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonGeneralizedExponential.md)/[`dpoisweibull()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonWeibull.md)
  families of functions keep their scale arguments.

- Underreporting is binomial thinning of the true counts, and the model
  of the reported counts now matches it: for NB1 and PIG1 the reported
  counts’ dispersion is p \* alpha, and for NBP alpha \* p^(2 - P)
  (before, the true-count dispersion was used unchanged, which is not
  the distribution of thinned counts). Poisson, NB2, PIG2, PIG and the
  mixed-Poisson families (PLN, PL, PLG, PLL, PGE, PW, SI) were already
  exact. GP1, GP2, COM and GW are not closed under thinning, so
  underreporting is now an error for them. This changes NB1, NBP and
  PIG1 underreporting estimates **(changed)**.

- Underreporting identification: an intercept-only `underreport_formula`
  is an error (only exp(b0) \* p is identified), and
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  warns when the reporting and count coefficients are nearly collinear
  (correlation above 0.99), which happens when reporting is separated
  from the count mean only by the curvature of the link function.

- Likelihood-ratio tests are reported only for valid nesting.
  [`summary()`](https://rdrr.io/r/base/summary.html) no longer reports a
  test against Poisson for PL, PLG and PLL (which cannot reduce to
  Poisson: their variance is at least mu + mu^2 / 2), for NBP, SI, PW,
  PGE and GW (nonregular restrictions), or for boundary families with a
  dispersion formula.
  [`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md)
  uses the same rules (same estimator and model components, nested
  parameters and families), uses a 50:50 chi-square mixture when one
  parameter is on the boundary (also with extra mean covariates),
  requires two fitted models to use the same observations, outcome,
  offsets and weights, returns `LR_type`, and gives no p-value (instead
  of 1) when the larger model has the lower log-likelihood. The README
  zero-inflation example now compares models by AIC and BIC.

- Likelihoods use log-probabilities computed on the log scale (every
  family density has a `log` argument), so extreme counts keep finite
  log-likelihoods (log P(Y = 1000 \| mu = 1) = -5913 instead of log(0)).
  Rows or panels with zero weight contribute exactly 0.

- The Conway-Maxwell-Poisson functions use the log normalizing constant
  (it overflowed for lambda = 1000, nu = 1, so
  `dcom(1000, lambda = 1000, nu = 1)` was 0 instead of 0.0126), one
  support rule for the PMF, CDF, quantiles and moments, moments computed
  from the normalized PMF, and a bracketed solver for the rate given the
  mean (accurate to about 1e-8, NaN with a warning when it fails,
  instead of silently returning an unconverged value).

- Correlated random parameters with three or more parameters: the
  `Chol:` labels and the default starting values were built row by row
  while the Cholesky factor was filled column by column, so labels were
  attached to the wrong elements, named starting values went to the
  wrong elements, and the default start was a singular matrix instead of
  0.1 \* I. Labels, starting values, the likelihood and predictions now
  use one packing. Estimates and predictions of earlier fits were
  internally consistent; only their labels were wrong.

- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  rejects a variable in both `formula` and `rpar_formula` when its
  random parameter is normal, uniform or triangular (only the sum of the
  fixed coefficient and the random mean is identified). Lognormal and
  gamma random parameters with a fixed shift are still allowed.

- Auxiliary formulas (distribution parameters, underreporting, random
  parameters, heterogeneity) keep their trained terms, so data-dependent
  bases such as [`poly()`](https://rdrr.io/r/stats/poly.html),
  [`splines::ns()`](https://rdrr.io/r/splines/ns.html) and
  [`scale()`](https://rdrr.io/r/base/scale.html) are evaluated with the
  training values at prediction time (a row is predicted the same alone
  or in any batch).

- Bootstrapped standard errors use only replicates whose refit converged
  and gave finite estimates for every parameter (failed optimizations
  were counted as successes before). The warning about failed replicates
  gives the reasons, the counts are stored (`boot_failures`), and the
  covariance uses the accepted replicates only (not pairwise deletion).

- Frequency weights now give the same results as expanding the data: the
  robust (sandwich) covariance uses sum(w s s’) instead of sum(w^2 s
  s’), and BIC uses the sum of the weights (of rows, or of panels in
  panel models) as the sample size.
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  weights that vary within a panel are an error (they were silently
  replaced by the first weight of the panel).

- All count p-functions
  ([`pplind()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindley.md),
  [`ppinvgaus()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGaussian.md),
  [`ppinvgamma()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGamma.md),
  [`ppLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLognormal.md),
  [`pplindLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md),
  [`pplindGamma()`](https://jwood-iastate.github.io/flexCountReg/reference/NegativeBinomialLindley.md),
  [`ppge()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonGeneralizedExponential.md),
  [`ppoisweibull()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonWeibull.md),
  [`psichel()`](https://jwood-iastate.github.io/flexCountReg/reference/SichelDistribution.md),
  [`pgwar()`](https://jwood-iastate.github.io/flexCountReg/reference/Generalized-Waring.md),
  [`pcom()`](https://jwood-iastate.github.io/flexCountReg/reference/COMDistribution.md))
  follow one contract: NA gives NA (it gave 0), q \< 0 and -Inf give 0,
  Inf gives 1 (several gave 0 or failed with “result would be too long a
  vector”), the upper tail and log scale follow, and zero-length input
  gives `numeric(0)`. Small upper tails are summed directly on the log
  scale instead of taken as 1 - CDF.

- [`cureplot()`](https://jwood-iastate.github.io/flexCountReg/reference/cureplot.md)
  uses the rows used in estimation by default, builds the outcome,
  predictions and plotting variable with one mask of complete rows (they
  could be misaligned after missing values), errors for an unknown
  `indvar`, and passes `method` only to random-parameter and negative
  multinomial models (the default `NULL` uses the estimator’s default).

- Input checks before optimization, shared by all estimators: the
  response must be finite, non-negative whole numbers (non-integer
  counts are now an error, not a warning); design matrices must be
  finite after formula transformations
  (e.g. [`log()`](https://rdrr.io/r/base/Log.html) of a non-positive
  value) and of full column rank (aliased columns are named); weights
  must be finite with a positive total; `start.vals` must be finite;
  `max.iters` and `ndraws` must be whole numbers.

- Random-parameter predictions: Simulated predictions now report a
  nonexistent expected value (lognormal coefficient with positive x,
  gamma coefficient beyond its MGF domain) like Exact predictions,
  instead of a finite average of draws. `rp_mgf()` handles zero-scale
  lognormal and gamma coefficients exactly and computes uniform and
  triangular MGFs on the log scale.

- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  has a `mix_draws` argument for the inner integration of the PLN, PGE,
  PLL and PW families (default `min(ndraws, 200)`), so the cost is
  proportional to `ndraws * mix_draws` instead of `ndraws^2`. Scrambled
  Halton draws are generated from a seed stored in the fit, so
  predictions reuse the estimation draws and leave the caller’s random
  numbers unchanged.

- When fewer than 2 bootstrap replicates succeed, the fit is returned
  with NA standard errors and a warning, instead of an error.

- Documentation now matches the implementation: the NB2 PMF uses size
  1/alpha (as the variance and the code do); the PLL regression mean is
  exp(X beta + sigma^2 / 2) with the corresponding variance; the PLL
  integral uses the Poisson-Lindley scale lambda, not the mean; the PL
  regression section gives the base PMF evaluated at lambda (the
  expanded formula had the wrong power of theta + 2) and no longer
  describes an LR test against Poisson; the PGE integral is over the
  generalized-exponential distribution. Fitted objects record the
  package version that created them (`model$fcr_version`).

- Packaging: `DESCRIPTION` has `BugReports` and the GitHub URL; a
  `CITATION` file and a contributor guide with a release checklist were
  added; CI also checks the minimum R version; the slowest stochastic
  tests are skipped on CRAN.

- The minimum R version is now 4.5.0, because the current CRAN release
  of the gsl package (an import) requires it. \## Bug fixes
  (distribution functions and inputs)

- [`pcom()`](https://jwood-iastate.github.io/flexCountReg/reference/COMDistribution.md)
  ignored `lower.tail` and `log.p`. Both tails are now summed directly
  on the log scale with a precise normalizing constant, so upper tails
  are never negative and keep their relative precision (e.g.
  `pcom(25, lambda = 1, nu = 1, lower.tail = FALSE)` equals the Poisson
  value 9.5e-28). Recycling of `mu`, `lambda` and `nu` in
  [`dcom()`](https://jwood-iastate.github.io/flexCountReg/reference/COMDistribution.md),
  [`pcom()`](https://jwood-iastate.github.io/flexCountReg/reference/COMDistribution.md)
  and
  [`qcom()`](https://jwood-iastate.github.io/flexCountReg/reference/COMDistribution.md)
  follows the base-R rule (length 1 or the length of the first argument;
  other lengths are an error).

- [`dpoisweibull()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonWeibull.md)
  and friends used `lambda` inconsistently between the PMF and the mean;
  parameters given only partly now warn instead of failing.

- The d/p/q functions of the Poisson-Lindley, Poisson-Lindley-Lognormal,
  Poisson-inverse-Gaussian, Poisson-inverse-gamma, Poisson-generalized-
  exponential, Sichel and inverse-gamma distributions return `NaN` (with
  one warning) for invalid parameters instead of failing, returning 0,
  or emitting extra warnings. Missing inputs give `NA`.

- [`qsichel()`](https://jwood-iastate.github.io/flexCountReg/reference/SichelDistribution.md)
  returns `Inf` at p = 1 (it returned `NA` with a coercion warning) and
  `NaN` for p outside \[0, 1\] (it searched 10,000 counts).
  [`rsichel()`](https://jwood-iastate.github.io/flexCountReg/reference/SichelDistribution.md)
  and
  [`rinvgamma()`](https://jwood-iastate.github.io/flexCountReg/reference/invgamma.md)
  accept a vector `n` (its length is used, as in base R).

- [`halton_dists()`](https://jwood-iastate.github.io/flexCountReg/reference/halton_dists.md)
  draws one dimension by default, uses the stated half-width for uniform
  draws (`mean +/- sdev`), and errors for an unknown distribution.
  [`corr_haltons()`](https://jwood-iastate.github.io/flexCountReg/reference/corr_haltons.md)
  gives each column its own mean, standard deviation and bounds for
  truncated normal draws, and errors for an unsupported distribution,
  `hdraws` with the wrong number of columns, or a Cholesky factor of the
  wrong size.

- [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) work (with
  `NA` standard errors) when a model has no covariance matrix.

- `bootstraps` is checked up front (a single number of at least 2)
  instead of failing after the model is fit.

- Data preparation rejects infinite counts and non-numeric (e.g. factor)
  responses, expands `.` in formulas, and keeps the contrasts of factors
  in auxiliary formulas for prediction. Fitted objects store the
  expanded formulas, so bootstrap refits of models written with `.`
  work.

- `method = "Nelder-Mead"` is accepted by every estimator as an alias of
  `"NM"`.

- [`qplind()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindley.md)
  and
  [`qplindLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md)
  return `NA` for missing parameters instead of failing, and
  [`qplindLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md)/[`rplindLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md)
  accept `sigma = 0` (the Poisson-Lindley limit).

- [`regCompTable()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTable.md)
  fits the baseline model only when the pseudo-R^2 is requested.
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) classifies
  parameters from the model’s bookkeeping (covariates named like
  distribution parameters are no longer misclassified).

### Documentation

- The pkgdown site renders mathematics with MathJax
  (`template: math-rendering: mathjax`). The equations of help pages
  (`\eqn{}`, `\deqn{}`) and articles were shown as raw LaTeX, for
  example in
  [`?countreg_mix`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md).

- The package description (DESCRIPTION) now summarizes the current
  scope: the 28 count families and their distribution functions,
  distribution-parameter, zero, underreporting, random-parameter,
  random-effects and finite-mixture models, the prediction, interval,
  effect and diagnostic tools, and references for latent-class
  regression, randomized quantile residuals and rootograms.

- Corrected the parameter table in
  [`?countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  (PGE, Sichel), the Poisson-Lindley-Lognormal mean and parameter
  descriptions, the RENB density, the Sichel reference (Zou, Lord and
  Zhang, 2012) and several typos.
  [`?countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  no longer claims that an LR test is printed; see
  [`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md).

- [`?countreg_rp`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  documents each random-parameter distribution’s scale, the
  heterogeneity specification and the lognormal caveat.

- [`?renb`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)
  states that the marginal variance exists only for a \> 2 (it is
  infinite for the legal values 1 \< a \<= 2).

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

### Other new features

- New
  [`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md):
  average marginal effects, elasticities and, for indicators,
  pseudo-elasticities or discrete changes, for every estimator, with
  delta-method or bootstrap (independent units, optional `cluster_var`)
  standard errors. Variables that change the expected count only through
  auxiliary formulas (underreporting, zero model, random parameters and
  their heterogeneity, and distribution parameters when the mean depends
  on them, as in hurdle models) are included by default.

- New panel estimators
  [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md)
  (random-effects Poisson) and
  [`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md)
  (negative multinomial NM1, NM2 and NM-P), with weights, robust and
  bootstrap standard errors, predictions and the common methods.

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
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  models now computes Empirical Bayes observation-specific predictions
  following Wood and Gayah (2025), weighting each draw by the likelihood
  of the unit’s observed outcomes. Previously it returned the simulated
  predictions. Panels are matched to the rows of `newdata` by id (a
  `group_var` given as a vector is stored as the column
  `.fcr_panel_id`), and the response expression is evaluated as in
  estimation. `eb_type = "mean"` returns the conditional mean instead of
  plugging in the conditional coefficients; the conditional coefficients
  and their standard deviations are returned as attributes.

- [`predict()`](https://rdrr.io/r/stats/predict.html) gains
  `type = c("response", "true", "prob_reported", "link")` for
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  models.

- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  accepts `weights` as a numeric vector as well as a column name.

### Deprecations

- [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  is deprecated; use
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md),
  now the documented name
  ([`?countreg_rp`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)).
  [`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  still works, with a once-per-session warning, and returns the same
  fit.
- `predict(object, data = ...)` is deprecated; use `newdata`.
- `print.level` in
  [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)
  and
  [`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)
  is deprecated; use `verbose`. `panel_id` in
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
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
  [`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)
  and
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  gain `stderr = c("normal", "robust", "boot")`. Robust standard errors
  of panel models are clustered by panel; bootstraps resample whole
  panels.
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  gains `bootstraps`.
- All estimators gain `keep_data`; `keep_data = FALSE` does not store
  the training data.
- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  and
  [`poislind_re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)
  are snake_case names for
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  and
  [`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md).
- For fitted models, [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) replace
  [`myAIC()`](https://jwood-iastate.github.io/flexCountReg/reference/myAIC.md)
  and
  [`myBIC()`](https://jwood-iastate.github.io/flexCountReg/reference/myBIC.md).

## flexCountReg 0.1.1

- Initial CRAN submission.
