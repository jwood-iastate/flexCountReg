# flexCountReg 0.1.2 (development)

## Interactive app

* New `launch_app()` opens an interactive interface to flexCountReg in a
  web browser, with a tab for each step: load a bundled example or a CSV
  file and create new variables (logs, powers, standardized values,
  indicators, products, ratios, categories, factors); specify and estimate
  count regression, random-parameter, random-effects panel (RENB,
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
  new data get the variables created on the Data tab that the model needs,
  with the means, standard deviations and categories of the training data.
* Everything the app runs (data loading, new variables, fits, diagnostics,
  predictions, intervals and comparisons) is compiled into R code that
  calls the package's functions, evaluated from that code and appended to
  a session script, shown in the app and downloaded as `analysis.R`; the
  script reproduces the whole session without the app (the download also
  holds the settings of each fit as JSON and a manifest with versions and
  warnings). Column names and typed values enter the code only as names,
  numbers or strings, never as code. An action that fails reports its
  error and leaves the session running.
* Models are estimated in background R processes (callr), so the
  interface stays responsive and fits can be canceled; a slow earlier fit
  never replaces a newer one. Predictions and intervals also run in the
  background: an interval for a model whose covariance matrix is not
  positive definite uses a parametric bootstrap (200 refits), which can
  take many minutes, and the app says so when it starts and offers to
  cancel it. Data, fits and jobs belong to one browser
  session.
* The app's packages (shiny, bslib, DT, callr, jsonlite, zip) are
  suggested, not imported: the statistical functions do not need them, and
  `launch_app()` offers to install any that are missing. `app_object()`
  returns the app for deployment or testing.

## Finite mixtures and latent classes

* New `countreg_mix()` estimates finite mixtures of count distributions and
  finite-mixture count regressions; `countreg_lc()` is the same model under
  the latent-class vocabulary (`classes`, `membership_formula`). Both
  normalize their arguments into one specification and call one engine, so
  equivalent calls give identical fits. Supported specifications:
  * distribution-only mixtures (`y ~ 1`), including fits to exact-count
    frequency tables through `weights`;
  * mixture regressions with component-specific coefficients, or with whole
    terms shared by all components (`shared = list(mean = "x")`) and shared
    distribution-parameter blocks (`shared = list(dis1 = TRUE)`);
  * constant mixing proportions, covariate-dependent membership
    probabilities (a multinomial logit, `mixing_formula = ~ z`, last
    component as reference), or known fixed proportions (`fixed_mixing`),
    which add no parameters;
  * any count family of `countreg()` in the components, the same in all of
    them or different ones (e.g. `family = c(base = "POISSON",
    tail = "NB2")`), with distribution-parameter formulas;
  * panel mixtures (`group_var`): one class per panel, with the class
    probability entering once per panel (product of the rows' probabilities
    before the sum over classes);
  * offsets and frequency weights (outside the logarithm of the mixture;
    constant within panels).
* Estimation uses several reproducible starts (`n_starts`, `seed`), EM
  iterations with weighted M-steps (analytic scores for Poisson, NB1 and
  NB2 components), and direct maximization of the observed-data
  log-likelihood to polish the best solutions (`mixture_control()`,
  `lc_control()`). Every start is recorded (`fit$model$starts`), including
  whether several starts reached the selected solution. Components of the
  same family are numbered by increasing average expected count.
* Standard errors come from the observed information of the full parameter
  vector (`stderr = "normal"`), unit-level sandwich scores (`"robust"`), or
  a bootstrap of independent units that refits the full mixture with
  several starts and aligns component labels within families (`"boot"`).
  The fit records the gradient, the definiteness and conditioning of the
  observed information, and ambiguous bootstrap label matches.
* Methods: `summary()` (coefficients grouped by component, shared
  parameters and mixing model; prior and posterior class sizes; entropy,
  relative entropy and ICL; start stability; numerical checks), `print()`,
  `coef()`, `vcov()`, `logLik()` (free parameters only, independent units
  for BIC), `predict()` with explicit targets (`type = "response"`,
  `"component"`, `"variance"`, `"prob"`, `"link"`; `conditioning =
  "prior"`, `"history"` for panel forecasts from earlier outcomes only, or
  the retrospective `"posterior"`), `fitted()`, `residuals()` (response,
  Pearson with the full mixture variance, and randomized quantile residuals
  from the mixture distribution function), `simulate()` (one class per row
  or panel), `tidy()`, `glance()` and `augment()`. New
  `component_probs()` / `class_probs()` return prior, history-conditioned
  or posterior class probabilities.
* `rootogram()` and `cureplot()` accept mixture models (prior mixture
  probabilities and means), and `margEffTable()` includes variables that
  enter only the membership model (delta-method standard errors with the
  full covariance matrix).
* Not available yet, with explanatory errors: confidence, prediction and
  tolerance intervals for mixture predictions; bootstrap standard errors
  of marginal effects; likelihood-ratio tests involving mixtures
  (`regCompTest()` refuses them, because the chi-square reference
  distribution does not hold when the number of components changes).
  Zero-inflated, hurdle and underreporting components inside mixtures and
  component-specific formulas are not implemented. Mixtures of families
  other than Poisson, NB1 and NB2 are labeled experimental: their
  component distributions are those of `countreg()` (one-component fits
  agree with it), but recovery of such mixtures has not been studied.
* New vignettes: "Finite mixtures of count distributions" and
  "Latent-class and finite-mixture count regression", and a recovery and
  coverage study, `inst/validation/finite-mixtures/` (its `CORES` setting
  runs the replicates in parallel with the same results). In 200
  replicates of each of seven designs (well and weakly separated Poisson
  mixtures, a latent-class NB2 regression with membership covariates and a
  shared dispersion, a panel latent-class model, and a Poisson + NB2
  mixture), every fit converged; Wald intervals covered 0.89 to 0.985
  (Monte Carlo standard errors about 0.015), and the label-invariant
  mixture mean was nearly unbiased in every design. With weakly separated
  components and n = 300 the low component's intercept was biased
  (-0.18) and its standard errors unstable, and 40% of the fits were
  reached by only one of five starts: use more starts and treat
  component-specific inference with caution in such data.
* A warning (and a note in `summary()`) flags components of the same
  family with nearly identical parameters: the fit is then effectively
  one component fewer and its shares are arbitrary.

## New count distributions

* Eight new count distributions, each with density, distribution, quantile
  and random-generation functions following base R conventions (recycling,
  `NA`/`NaN`, `NaN` with a warning for invalid parameters, `lower.tail`,
  `log`/`log.p`):
  * discrete Weibull (`ddweib()`, `pdweib()`, `qdweib()`, `rdweib()`);
  * double Poisson, exactly normalized (`ddpo()`, ...);
  * Delaporte in a mean parameterization (`ddelaporte()`, ...);
  * generalized Hermite of order `m` (`dghermite()`, ...);
  * hyper-Poisson (`dhpois()`, ...);
  * gamma renewal count (`dgcount()`, ...);
  * Polya-Aeppli, the Poisson-Tweedie distribution with power 1.5
    (`dpaeppli()`, ...).
* Probabilities are computed on the log scale, and upper tails
  (`lower.tail = FALSE`) are computed directly rather than as one minus the
  distribution function, so tiny tail probabilities keep their relative
  accuracy. Quantiles are the smallest count whose distribution function
  reaches the requested probability.
* Negative binomial-generalized exponential mixture (`dnbge()`, `pnbge()`,
  `qnbge()`, `rnbge()`): NB2 with a mean-one generalized exponential random
  effect, computed by the adaptive quadrature used for PGE.

## New regression families

* `countreg()` and `countreg.rp()` gain the families `"DW"`, `"DPO"`,
  `"DELAPORTE"`, `"HERMITE"` (order 2), `"HERMITE3"`, `"HERMITE4"`,
  `"HYPERPOIS"`, `"GAMMACOUNT"`, `"PA"` (Poisson-Tweedie with power 1.5) and
  `"NBGE"`. All are parameterized by the mean, so coefficients are effects on
  the expected count. For DW, DPO, HYPERPOIS and GAMMACOUNT, whose means have
  no closed form, the natural parameter of each observation is found by
  solving the mean equation in compiled code; these families are slower to
  fit and have a numerical domain (shape in [1e-6, 1e6], mean up to 1e9).
* Shares in [0, 1] (Delaporte and Hermite `rho`) use a logit link, labeled
  `logit(rho)`.
* Prediction, Pearson and quantile residuals, simulation, rootograms and
  summaries work for all new families. Likelihood-ratio tests against
  Poisson are reported for DPO, HYPERPOIS and GAMMACOUNT (interior
  restriction). Underreporting (binomial thinning) is available for
  DELAPORTE, NBGE, HERMITE and PA.
* `countreg.rp()` uses family-specific starting values for the distribution
  parameters when the registry provides them (the new families start at
  their Poisson or geometric limit rather than at log-scale -2).

## Confidence, prediction and tolerance intervals

* `predict()` gains `interval = "confidence"` (for the expected count) and
  `interval = "prediction"` (for one new count), with `level`, for
  `countreg()` models. Estimation uncertainty is represented by parameter
  vectors from the asymptotic normal distribution of the estimates,
  evaluated at deterministic Halton points (reproducible), or by a
  parametric bootstrap (`uncertainty = "bootstrap"`), which is also used
  when the covariance matrix is not positive definite. Prediction intervals
  are integer quantiles of the predictive distribution (the count PMF
  averaged over the parameter vectors).
* New `tolerance_interval()`: content tolerance intervals (content and
  confidence levels; two-sided, upper or lower; pointwise or simultaneous)
  for `countreg()` models. The default joint method returns the narrowest
  integer interval whose content is at least the requested fraction at
  every evaluated point of the confidence set; the envelope method uses
  per-tail worst cases; intercept-only Poisson models use the exact
  common-rate construction. Zero-inflated, hurdle and underreporting
  components are included. A reproducible coverage study is in
  `inst/validation/intervals/`.

## Starting values for the dispersion

* `countreg()` now starts the dispersion parameter of NB1, NB2, PIG1, PIG2
  and PIG at a moment estimate from the Poisson starting fit (the
  auxiliary regression of Cameron and Trivedi) instead of `alpha = 1`. With
  the old start, quasi-Newton optimizers could step onto the flat region
  `alpha -> 0` and stop there with a falsely "converged" Poisson fit: NB2
  with `method = "BFGS"` gave log-likelihood -228.48 instead of -195.56 on
  `warpbreaks`, and -182.29 instead of -180.11 on `InsectSprays`. Results
  now match `MASS::glm.nb()`. Fits that already reached the maximum are
  unchanged.

## Standard errors at flat likelihood directions

* When the Hessian is singular, standard errors still use a generalized
  inverse, but parameters involved in a direction where the likelihood is
  flat now get `NA` variances and standard errors instead of zero (a zero
  standard error read as perfect precision, e.g. the random-effect
  parameter of a Poisson-Lindley random-effects model at its boundary). The
  warning names those parameters.
* When the Hessian is invertible but not negative definite (for example an
  NB2 dispersion drifting toward zero on Poisson-like data), the negative
  variances already gave `NA` standard errors; a warning now names those
  parameters instead of leaving them silently `NA`.

## Vignettes

* Seventeen vignettes, grouped on the package website as Start here,
  Diagnose and interpret, Model structure, Finite mixtures and latent
  classes, Predict and validate, and Advanced. New (besides the two
  finite-mixture vignettes above): getting started; data, exposure and weights; diagnosing
  count models; predictions, marginal effects and elasticities; modeling
  distribution parameters; zeros, hurdles and sampling; underreporting and
  identification; predictive validation; prediction for panel data;
  numerical reliability; simulation studies; and adding a count
  distribution (for contributors). The three earlier vignettes are revised.
* `inst/validation/simulation-studies/count-regression-study.R` is a
  parameterized simulation study (Poisson, NB2 and a negative multinomial
  panel) with one random-number stream per replicate; its summary and
  settings are saved alongside it. The per-replicate results of this study
  and of the interval coverage study are kept in the source repository but
  left out of the package build.
* The vignettes use `countreg_rp()`.
* A control character in the negative-multinomial paragraph of the
  simulation-studies vignette is replaced by the intended `\alpha`.

## New data

* `abalone_rings`: shell-ring counts and measurements of 4,177 abalone from
  the UCI Machine Learning Repository (CC BY 4.0), a non-transportation
  example with conditionally underdispersed counts. The preparation script
  and source hashes are in `data-raw/` of the source repository; data terms
  are in the installed `DATA-LICENSES.md`.
* `stream_salamanders`: repeated counts of seven stream salamander groups at
  23 sites (Price et al., 2016; Dryad, CC0), prepared from the original
  deposit. Zero-heavy, overdispersed and grouped by site.

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

* `margEffTable()` no longer reports an effect for a variable that enters
  the model only through `offset()` in the formula (for example the
  exposure in `offset(log(exposure))`); such variables were treated as
  predictors, while offsets given with the `offset` argument were already
  excluded. A variable that is also a predictor is still evaluated, and
  offset variables can still be named in `vars`.
* `cureplot()` no longer prints the cureplots messages "Covariate:
  indvar_values" and "CURE data frame was provided...", which named the
  wrapper's internal objects; other messages are kept. The documentation
  of `n_resamples` now says that the red reference limits are drawn
  whatever its value (they were described as "no bands" at 0).

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

## Corrections from the follow-up review

Changed results are marked **(changed)**.

* The underreporting identification check also warns when standard errors
  are not available because the Hessian is singular (a flat likelihood
  ridge). Before, it returned silently in that case, so whether a weakly
  identified model warned depended on the platform's floating point.
* **Weights for every estimator.** `renb()`, `repois()`, `poisLind.re()` and
  `negmultinom()` accept `weights`, as `countreg()` and `countreg.rp()` do.
  In panel models these are frequency weights of the panels (constant
  within a panel): a panel with weight w counts as w independent panels in
  the log-likelihood, the robust covariance, BIC and the bootstrap. Fits
  with weights equal fits to the expanded panels.
* **Frequency-weight bootstrap (changed).** With `weights`, bootstrap
  standard errors (coefficients and `margEffTable()`) now draw new
  frequencies from Multinomial(sum(w), w / sum(w)) and refit with them in
  place of the weights. This equals the bootstrap of the expanded
  observations (or panels). Before, rows were resampled with their original
  weights, which is a different and much more variable sampling scheme.
  Bootstrap standard errors need whole-number weights; use
  `stderr = "robust"` for other weights.
* **Bootstrap acceptance.** A replicate counts only if its refit reports
  convergence (unknown convergence is a failure), has exactly the expected
  coefficients and finite estimates. `negmultinom()` uses the shared
  bootstrap. Failures are stored by reason in `model$boot_failures`.
* **Marginal-effect bootstraps** refit with the frozen, expanded formulas
  (a `y ~ .` formula no longer picks up helper columns), stored argument
  values (not objects that changed since the fit), a collision-safe panel
  column and panel keys of (sampled cluster occurrence, original panel), so
  a `cluster_var` cluster drawn twice keeps its panels distinct. They
  resample the analysis rows with the stored weights and panels.
  `cluster_var` with weights, and panels split across clusters, are
  rejected. Without `data`, effects average over the rows used in
  estimation, counting each with its frequency weight **(changed** for
  weighted models**)**.
* **Likelihood-ratio tests.** Fits store evidence of their analysis sample
  (rows, response, weights, offsets, panels and the variables they use).
  An LR test between two models requires the same sample (shared variables
  must agree row by row), reported convergence for both, nested designs
  (checked by column space, so equivalent bases such as `x` and `2 * x` are
  accepted) and matching random-parameter distributions and correlation. A
  larger model whose log-likelihood is materially below the restricted
  one's (LR below -max(0.01, 2e-6 |logLik|)) gives no test, with a reason
  in `LR_reason`, instead of LR = 0. The intercept-only base model of
  `regCompTest()` is not nested in a model without an intercept. Objects
  from older versions get no LR test between fitted models (refit them).
* **Mixture families (changed).** The PLN, PGE, PW and PLL probabilities
  are computed on the log scale by deterministic adaptive quadrature over the
  log of the mixing variable (`engine = "quadrature"`, the new default of
  the d/p/q functions, `countreg()` and `countreg.rp()`). Averaging over
  1500 Halton draws was badly biased in the tails: for PGE with shape 1
  (geometric with mean 1) and y = 1000 it gave log P = -3847 instead of
  -693.8. The new values agree with exact special cases (PGE with shape 1
  and PW with alpha 1 are geometric; P(Y = 0) of PGE is a beta function) and
  with independent quadrature to about 1e-9, for counts up to 1e5. Fits are
  faster too. `engine = "halton"` keeps the earlier method, and objects
  without a stored engine predict as before.
* **GP domains.** Generalized Poisson moments and predictions are `NaN`, with
  a warning, outside the parameter space (for GP2, 1 + a mu > 0), instead
  of textbook moments. In `countreg.rp()` the GP2 dispersion must be
  nonnegative, since random parameters with unbounded support violate the
  bound for some draws; use GP1 for underdispersed random-parameter models.
* **Moment existence.** Random-parameter predictions decide whether
  E[exp(x beta)] exists from the parameters alone, so `Simulated`
  predictions with lognormal parameters no longer need lamW. `Exact` uses
  log-scale MGFs and reports an overflow of a finite expectation as such.
  For lognormal parameters with x < 0, `Exact` uses a Laplace approximation
  (documented as approximate).
* **Quantile residuals** evaluate only F(y - 1) and F(y): directly for
  Poisson and the negative binomial families (no cost growth with the
  counts), and by one streaming pass with a work budget (an error before
  any work when it would be exceeded) for the others. Zero-inflated and
  hurdle CDFs are computed stably.
* **Rootograms (changed).** Expected frequencies sum the fitted
  probabilities of every row, including rows whose outcome is above
  `max_count`; before, those rows were dropped. New `overflow` argument for a
  `Y > max_count` bin; the tail is also returned in attribute `"tail"`. Rows
  with non-finite fitted means are left out of both frequencies, with a
  warning.
* **Distribution contracts.** Invalid parameters give `NaN` for that
  element, with one warning (base R's convention), and missing parameters
  `NA`, at every `q` including `q < 0` and `q = Inf` (before, those endpoints
  returned 0 or 1). `dgwar()`/`pgwar()` no longer stop for an invalid
  element. `plindley()` computes both tails directly. Zero-truncated
  simulation reports invalid probabilities and its count limit as `NA`
  draws with a warning, instead of returning a count.
* **COM controls.** `find_lambda_cpp()` honors its tolerance (relative
  error of the mean) and iteration limit and reports non-convergence
  instead of accepting a collapsed bracket; `find_lambda_diag_cpp()` gives
  the residual, iterations and status. Unused controls of the internal
  moment functions were removed.
* **`renb()`**: the log-likelihood is computed as three stable differences
  of log-gamma functions. For large `a` and `b`, terms of size ~1e200
  canceled and the optimizer could report a positive log-likelihood of
  ~1e203 **(changed**, fits that hit this region**)**.
* **`renb()` standard errors (changed).** `renb()` now has analytic
  scores, so its Hessian is the derivative of an exact gradient. Before, it
  was a second finite difference of the log-likelihood, which was
  inaccurate (it could be indefinite and give `NaN` standard errors or
  delta-method marginal-effect standard errors, notably for `ln(a-1)` and
  `ln(b)`). Estimates are unchanged; standard errors can change.
* `rpardists` values are checked before a single value is expanded.
* Documentation: the PGE equations use one rate convention and explain why
  the rate cancels; fitted objects carry a schema version
  (`model$fcr_schema`).

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
* Auxiliary formulas (distribution parameters, underreporting, zero model,
  random parameters, heterogeneity) keep their trained terms, so
  data-dependent bases such as `poly()`, `splines::ns()` and `scale()` are
  evaluated with the training values at prediction time (a row is predicted
  the same alone or in any batch).
* `negmultinom()` uses the shared data preparation and stores the training
  blueprint: predictions use the training factor levels, contrasts and
  bases and coefficients by name; a missing offset column is an error
  instead of silently dropping the offset; giving both a formula offset and
  `offset` is an error; EB predictions handle transformed responses; named
  `start.vals` are matched by name; convergence uses the optimizer's own
  return codes (a BFGS fit stopped at the iteration limit was reported as
  converged); and the stored BIC uses the number of panels, as `BIC()` does.
* `renb()`, `poisLind.re()` and `repois()` predictions use the estimates by
  name rather than a cached copy, so their delta-method marginal-effect
  standard errors are no longer zero.
* The `margEffTable()` bootstrap resamples the model's independent units by
  default (panels for panel models), gives every sampled cluster a fresh id
  (a cluster drawn twice is two panels), resamples a numeric weights vector
  with its rows, refits without an inner bootstrap, counts only converged
  refits, and errors for a `cluster_var` that is not in the data.
* Bootstrapped standard errors use only replicates whose refit converged
  and gave finite estimates for every parameter (failed optimizations were
  counted as successes before). The warning about failed replicates gives
  the reasons, the counts are stored (`boot_failures`), and the covariance
  uses the accepted replicates only (not pairwise deletion).
  `negmultinom()` applies the same rule.
* Frequency weights now give the same results as expanding the data: the
  robust (sandwich) covariance uses sum(w s s') instead of sum(w^2 s s'),
  and BIC uses the sum of the weights (of rows, or of panels in panel
  models) as the sample size. `countreg.rp()` weights that vary within a
  panel are an error (they were silently replaced by the first weight of
  the panel). The parsnip engine warns that importance weights are treated
  as frequencies.
* Hurdle models: `simulate()` draws membership first and positive counts
  from the exact zero-truncated distribution (by inversion), instead of
  redrawing zeros at most 1,000 times and keeping the remaining zeros (about
  37% of the "positive" draws were zero for a count mean of 0.001).
  Predictions, Pearson variances and probability matrices use
  1 - P(Y = 0) = -expm1(log P(Y = 0)), which keeps its precision when
  P(Y = 0) is close to 1.
* All count p-functions (`pplind()`, `ppinvgaus()`, `ppinvgamma()`,
  `ppLnorm()`, `pplindLnorm()`, `pplindGamma()`, `ppge()`, `ppoisweibull()`,
  `psichel()`, `pgwar()`, `pcom()`, `pgp1()`, `pgp2()`) follow one contract:
  NA gives NA (it gave 0), q < 0 and -Inf give 0, Inf gives 1 (several gave
  0 or failed with "result would be too long a vector"), the upper tail and
  log scale follow, and zero-length input gives `numeric(0)`. Small upper
  tails are summed directly on the log scale instead of taken as 1 - CDF.
* `cureplot()` uses the rows used in estimation by default, builds the
  outcome, predictions and plotting variable with one mask of complete rows
  (they could be misaligned after missing values), errors for an unknown
  `indvar`, and passes `method` only to random-parameter and negative
  multinomial models (the default `NULL` uses the estimator's default).
  `augment(data = )` predicts the rows of `data` in the order given
  (before, data of the training size were matched to the cached fitted
  values by position).
* Input checks before optimization, shared by all estimators: the response
  must be finite, non-negative whole numbers (non-integer counts are now an
  error, not a warning); design matrices must be finite after formula
  transformations (e.g. `log()` of a non-positive value) and of full column
  rank (aliased columns are named); weights must be finite with a positive
  total; `start.vals` must be finite; `max.iters` and `ndraws` must be
  whole numbers.
* Random-parameter predictions: Simulated predictions now report a
  nonexistent expected value (lognormal coefficient with positive x, gamma
  coefficient beyond its MGF domain) like Exact predictions, instead of a
  finite average of draws. `rp_mgf()` handles zero-scale lognormal and gamma
  coefficients exactly and computes uniform and triangular MGFs on the log
  scale.
* `countreg.rp()` has a `mix_draws` argument for the inner integration of
  the PLN, PGE, PLL and PW families (default `min(ndraws, 200)`), so the
  cost is proportional to `ndraws * mix_draws` instead of `ndraws^2`.
  Scrambled Halton draws are generated from a seed stored in the fit, so
  predictions reuse the estimation draws and leave the caller's random
  numbers unchanged.
* `repois()` and `negmultinom()` compute their panel log-likelihoods and
  scores without the cancellation in lgamma(r + Y) - lgamma(r) for large
  r = 1 / alpha, so they tend continuously to the Poisson log-likelihood.
* When fewer than 2 bootstrap replicates succeed, the fit is returned with
  NA standard errors and a warning, instead of an error.
* Documentation now matches the implementation: the NB2 PMF uses size
  1/alpha (as the variance and the code do); the PLL regression mean is
  exp(X beta + sigma^2 / 2) with the corresponding variance; the PLL integral
  uses the Poisson-Lindley scale lambda, not the mean; the PL regression
  section gives the base PMF evaluated at lambda (the expanded formula had
  the wrong power of theta + 2) and no longer describes an LR test against
  Poisson; the PGE integral is over the generalized-exponential
  distribution. Version references in code and messages refer to 0.1.1
  (the last release) instead of a nonexistent 0.1.3, and fitted objects
  record the package version that created them (`model$fcr_version`).
* Diagnostics: randomized quantile residuals are computed in chunks of rows
  sorted by the outcome (bounded memory; one large count no longer widens
  every row's probability grid); `predict(type = "prob")` validates
  `max_count` and adds the tail mass P(Y > max_count) with
  `overflow = TRUE`; rootograms count frequency weights; and
  `simulate(seed = )` restores the caller's random-number state.
* Packaging: the parsnip engine is registered even when another package has
  already registered a `count_reg` model; `DESCRIPTION` has `BugReports`
  and the GitHub URL; a `CITATION` file and a contributor guide with a
  release checklist were added; CI also checks the minimum R version; the
  slowest stochastic tests are skipped on CRAN.
* The minimum R version is now 4.5.0, because the current CRAN release of
  the gsl package (an import) requires it.
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

* The package description (DESCRIPTION) now summarizes the current
  scope: the 28 count families and their distribution functions,
  distribution-parameter, zero, underreporting, random-parameter,
  random-effects and finite-mixture models, the prediction, interval,
  effect and diagnostic tools, and references for latent-class
  regression, randomized quantile residuals and rootograms.
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

* `countreg.rp()` is deprecated; use `countreg_rp()`, now the documented
  name (`?countreg_rp`). `countreg.rp()` still works, with a
  once-per-session warning, and returns the same fit.
* `predict(object, data = ...)` is deprecated; use `newdata`.
* `print.level` in `renb()`, `poisLind.re()` and `repois()` is deprecated; use
  `verbose`. `panel_id` in `countreg_rp()` is deprecated; use `group_var`.
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
