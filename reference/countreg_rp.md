# Random Parameters Count Regression Models

Random Parameters Count Regression Models

## Usage

``` r
countreg_rp(
  formula,
  rpar_formula,
  data,
  family = "NB2",
  rpardists = NULL,
  dis_param_formula_1 = NULL,
  dis_param_formula_2 = NULL,
  het_mean_formula = NULL,
  het_var_formula = NULL,
  ndraws = 500,
  scrambled = FALSE,
  correlated = FALSE,
  group_var = NULL,
  weights = NULL,
  offset = NULL,
  method = "BFGS",
  max.iters = 1000,
  start.vals = NULL,
  verbose = FALSE,
  stderr = "normal",
  bootstraps = NULL,
  engine = "quadrature",
  keep_data = TRUE,
  panel_id = NULL,
  mix_draws = NULL,
  weight_type = c("frequency", "sampling"),
  draws = c("per_unit", "shared"),
  cluster = NULL
)
```

## Arguments

- formula:

  an R formula. This formula should specify the outcome and the
  independent variables that have fixed parameters.

- rpar_formula:

  a symbolic description of the model related specifically to the random
  parameters. This should not include an outcome variable. If the
  intercept is random, include it in this formula. If the intercept is
  fixed, include it in `formula` but not in `rpar_formula`.

- data:

  a dataframe that has all of the variables in the `formula` and
  `rpar_formula`.

- family:

  the name of the distribution/model type to estimate. Default is "NB2".
  Options include "Poisson", "NB1", "NB2", "NBP", "PIG", "Sichel", etc.
  (see
  [`count_families`](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md)
  for the full list). DW, DPO, HYPERPOIS and GAMMACOUNT solve their mean
  equation for every observation and draw, so they are much slower here
  than closed-form families; consider fewer `ndraws` while exploring.
  With `"GP2"`, the dispersion must be nonnegative: the underdispersed
  GP2 requires \\1 + a\mu \> 0\\, which random parameters with unbounded
  support violate for some draws of every observation. Use `"GP1"`,
  whose parameter space does not depend on \\\mu\\, for underdispersed
  random-parameter models.

- rpardists:

  an optional named vector whose names are the random parameters and
  values the distribution. The distribution options include normal
  ("n"), lognormal ("ln"), triangular ("t"), uniform ("u"), and gamma
  ("g"). Values are matched to the random parameters by name; every
  random parameter must be named. A single unnamed value applies to all
  random parameters. If not provided, normal is used. See Details for
  how each distribution's `St.Dev` parameter is defined.

- dis_param_formula_1:

  a symbolic description of the model for the first parameter of the
  count distribution (e.g., ln(alpha) for NB2).

- dis_param_formula_2:

  a symbolic description of the model for the second parameter of the
  count distribution (if applicable).

- het_mean_formula:

  an optional one-sided formula for heterogeneity in the means of the
  random parameters. Each random parameter gets its own coefficient for
  each variable (named `HetMean:<parameter>:<variable>`), and the effect
  is additive: \\\beta\_{ik} = \mu_k + z_i'\delta_k\\ (Wood and Gayah,
  2025, Eq. 2). For lognormal parameters the shift is on the log scale,
  and for gamma parameters it is multiplicative, \\\mu_k
  \exp(z_i'\delta_k)\\, so the mean stays positive. The intercept is
  dropped.

- het_var_formula:

  an optional one-sided formula for heterogeneity in the scales of the
  random parameters: \\\sigma\_{ik} = \sigma_k \exp(s_i'\phi)\\, with
  \\\phi\\ shared by all random parameters. The intercept is dropped.

- ndraws:

  the number of Halton draws to use for estimating the random
  parameters.

- scrambled:

  if the Halton draws should be scrambled. Scrambled draws are generated
  from a random seed that is stored in the fit, so predictions reuse the
  same draws, [`set.seed()`](https://rdrr.io/r/base/Random.html) before
  fitting makes the fit reproducible, and different seeds give different
  draws.

- correlated:

  if the random parameters should be correlated. If TRUE, only normal
  distributions are used.

- group_var:

  an optional variable name (string), or a vector with one value per
  row, defining the panel structure (repeated measures). If provided,
  each panel shares one draw of the random parameters, and the
  likelihood and standard errors account for the panel structure.

- weights:

  the name of a variable, or a numeric vector, of weights (see
  `weight_type`). With `group_var`, weights apply to whole panels (the
  panel log-likelihood is weighted), so they should be constant within a
  panel.

- offset:

  the name of a variable, or vector of variable names, to be used as an
  offset (summed when several are given). An offset can instead be given
  with [`offset()`](https://rdrr.io/r/stats/offset.html) in `formula`;
  use only one of the two.

- method:

  optimization method passed to
  [`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html): `"BFGS"`
  (default), `"NR"`, `"BHHH"`, `"BFGSR"`, `"CG"`, `"NM"` or `"SANN"`.
  BFGS uses the analytic gradient and is the most reliable choice. BHHH
  is often faster on simple models, but on harder ones it can crawl, or
  stop where the log-likelihood has stopped changing well short of the
  maximum; such stops are flagged (see Details). NR is reliable near the
  maximum but slow, because it differentiates the gradient numerically
  for its Hessian.

- max.iters:

  maximum number of iterations.

- start.vals:

  optional vector of starting values, one per parameter. Named values
  are matched to the parameters by name.

- verbose:

  logical.

- stderr:

  the type of standard errors: `"normal"` (Hessian-based, the default),
  `"robust"` (sandwich; clustered by panel when `group_var` is given),
  or `"boot"` (bootstrap that resamples whole panels and refits the
  model), or `"cluster"` (cluster-robust, with `cluster`; clusters must
  contain whole panels, so they may be coarser than the panels, e.g.
  sites within corridors).

- bootstraps:

  the number of bootstrap samples when `stderr = "boot"` (default 200).

- engine:

  the numerical method for the mixed Poisson families PLN, PGE, PW and
  PLL: `"quadrature"` (default; adaptive quadrature on the log scale,
  accurate for every count), `"halton"` (the average over Halton draws
  used by flexCountReg 0.1.2 and earlier, kept to reproduce earlier
  fits; it is badly biased for counts far in the tail) or, for PLN only,
  `"poilog"` (the poilog package).

- keep_data:

  logical; if `FALSE`, the training data are not stored in the fitted
  object. [`predict()`](https://rdrr.io/r/stats/predict.html) then needs
  `newdata`.

- panel_id:

  deprecated; use `group_var`.

- mix_draws:

  with `engine = "halton"`, the number of Halton draws for integrating
  over the mixing distribution of the PLN, PGE, PLL and PW families
  (inside each random-parameter draw). The default is
  `min(ndraws, 200)`. Not used by the default `engine = "quadrature"`,
  which integrates adaptively.

- weight_type:

  `"frequency"` (default): each weight is the number of independent
  observations a row stands for, so the log-likelihood, LR tests, AIC
  and BIC are those of the expanded data. `"sampling"`: survey or
  inverse-probability weights, which only set how much each unit counts.
  The weighted log-likelihood is then a pseudo-likelihood: standard
  errors default to the robust sandwich with meat sum((w s)(w s)') and
  the small-sample factor G / (G - 1) for G units, so they do not change
  when all weights are rescaled; the number of units is the number of
  panels (rows without `group_var`) with a positive weight; bootstrap
  standard errors are not available; and LR tests, AIC and BIC are not
  valid (they are reported, with a note, by
  [`summary()`](https://rdrr.io/r/base/summary.html) and are refused by
  [`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md)).

- draws:

  `"per_unit"` (default): each unit (panel, or row when there is no
  `group_var`) has its own `ndraws` draws, a consecutive block of one
  Halton sequence after its first 10 points (Train 2009, Sec. 9.3.3), so
  the simulation errors of different units are independent. `"shared"`:
  the same draws for every unit, as in earlier versions. Empirical Bayes
  predictions (`method = "Individual"`) regenerate each estimation
  unit's own draws, matched by panel id (by row position without
  panels); units not in the estimation data get new draws.

- cluster:

  clustering variable(s) for `stderr = "cluster"`: a one-sided formula
  of columns of `data` (`~ corridor`; two or more, `~ corridor + year`,
  give multiway clustering), their names, or a vector with one cluster
  id per row of `data`. Given without `stderr`, it sets
  `stderr = "cluster"`. See
  [`vcov.flexCountReg`](https://jwood-iastate.github.io/flexCountReg/reference/vcov.flexCountReg.md),
  which also gives cluster-robust covariances after the fit.

## Value

An object of class `flexCountReg` which is a list with the following
components:

- model: the fitted model object.

- data: the data frame used to fit the model.

- call: the matched call.

- formula: the formula used to fit the model.

## Details

`countreg.rp()` is the former name of `countreg_rp()`. It is deprecated
(a name with a dot is easily mistaken for an S3 method): it still works,
with a once-per-session warning, and calls `countreg_rp()`.

**Random parameter distributions.** Each random parameter is \\\beta =
\mu + \sigma \upsilon\\, where the `:Mean` and `:St.Dev` estimates are
\\\mu\\ and \\\sigma\\:

- normal ("n"): \\\upsilon \sim N(0, 1)\\, so \\\sigma\\ is the standard
  deviation.

- uniform ("u") and triangular ("t"): \\\upsilon\\ is uniform or
  symmetric triangular on \\\[-1, 1\]\\, so \\\beta\\ lies in \\\[\mu -
  \sigma, \mu + \sigma\]\\. Here \\\sigma\\ is a half-width; the
  standard deviation is \\\sigma/\sqrt{3}\\ (uniform) or
  \\\sigma/\sqrt{6}\\ (triangular). Models fitted with version 0.1.1 or
  earlier used \\\[\mu - \sigma/2, \mu + \sigma/2\]\\ for the uniform;
  their \\\sigma\\ is converted automatically in
  [`predict()`](https://rdrr.io/r/stats/predict.html).

- lognormal ("ln"): \\\beta = \exp(\mu + \sigma \upsilon)\\, always
  positive. The expected value of \\\exp(x\beta)\\ is infinite for \\x
  \> 0\\, so marginal (Exact or Simulated) predictions do not exist for
  observations with a positive value of a lognormal random parameter's
  variable. Use the variable's negative (so \\x \le 0\\), a gamma random
  parameter, or `method = "Individual"` predictions instead.

- gamma ("g"): mean \\\mu \> 0\\ and standard deviation \\\sigma\\, so
  the coefficient is always positive.

The mixing distribution of the PLN, PGE, PLL and PW families uses its
own Halton dimension, independent of the random parameters' draws.

**Gradient.** The optimizer uses the analytic gradient of the simulated
log-likelihood: each panel's score is the average, over draws weighted
by their posterior probabilities, of the family's score in the log mean
times the derivative of the log mean in each parameter. The Poisson,
NB1, NB2 and NBP families have closed-form scores; the other families
use central differences in the log mean and in each distribution
parameter, two density evaluations per parameter of the distribution
whatever the number of regression parameters. A numeric gradient needs
two likelihood evaluations per parameter, so fits are several times
faster. `options(flexCountReg.rp_gradient = "numeric")` restores the
numeric gradient, to check a fit.

**Convergence.** BHHH, NR and BFGSR stop when the log-likelihood changes
by less than a tolerance between iterations, which a slow optimizer can
also meet far from the maximum. After the fit, `countreg_rp()` computes
the Newton decrement \\g'(-H)^{-1}g\\ at the estimates (twice the
further gain a Newton step predicts; stored as
`model$newton_decrement`). Above 0.1, the fit is marked as not
converged, with a warning.

## References

Train, K. (2009). *Discrete Choice Methods with Simulation*, 2nd ed.
Cambridge University Press.

Wood, J.S., Gayah, V. (2025). Out-of-sample prediction and
interpretation for random parameter generalized linear models. *Accident
Analysis and Prevention*, 220, 108147.

## Examples

``` r
# \donttest{
# Load data
data("washington_roads")
washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)

# 1. Basic Random Parameters Negative Binomial (NB2)
rp_nb2 <- countreg_rp(Total_crashes ~ lnaadt + lnlength,
                      rpar_formula = ~ -1 + speed50,
                      data = washington_roads,
                      family = "NB2",
                      rpardists = c(speed50 = "n"),
                      ndraws = 100)
summary(rp_nb2)
#> Call:
#> countreg_rp(formula = Total_crashes ~ lnaadt + lnlength, rpar_formula = ~-1 + 
#>     speed50, data = washington_roads, family = "NB2", rpardists = c(speed50 = "n"), 
#>     ndraws = 100)
#> 
#> Model: countreg_rp | family: NB2 | optimizer: BFGS
#> Iterations: 41 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>                Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)      -8.760      0.442  -19.82  < 2e-16 ***
#> lnaadt            1.080      0.052   20.79  < 2e-16 ***
#> lnlength          0.756      0.069   10.95  < 2e-16 ***
#> speed50:Mean     -0.745      0.167   -4.47  7.7e-06 ***
#> speed50:St.Dev    0.640      0.201    3.18   0.0015 ** 
#> ln(alpha)        -1.297      0.316   -4.10  4.1e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1083.165 | AIC: 2178.33 | BIC: 2210.213 
#> N: 1501

# 2. Random Parameters with Panel Structure (if 'site_id' exists)
# rp_panel <- countreg_rp(Total_crashes ~ -1 + lnaadt,
#                        rpar_formula = ~ speed50,
#                        data = washington_roads,
#                        group_var = "site_id",
#                        family = "NB2",
#                        ndraws = 100)

# 3. Generalized Random Parameters Model with Heterogeneity
rp_gen <- countreg_rp(Total_crashes ~ lnaadt,
                      rpar_formula = ~ -1 + speed50,
                      dis_param_formula_1 = ~ lnlength, 
                      het_mean_formula = ~ AADT10kplus,
                      data = washington_roads,
                      family = "NB2",
                      rpardists = c(speed50 = "n"),
                      ndraws = 100)
summary(rp_gen)
#> Call:
#> countreg_rp(formula = Total_crashes ~ lnaadt, rpar_formula = ~-1 + 
#>     speed50, data = washington_roads, family = "NB2", rpardists = c(speed50 = "n"), 
#>     dis_param_formula_1 = ~lnlength, het_mean_formula = ~AADT10kplus, 
#>     ndraws = 100)
#> 
#> Model: countreg_rp | family: NB2 | optimizer: BFGS
#> Iterations: 45 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>                             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)                  -8.4867     0.4521  -18.77  < 2e-16 ***
#> lnaadt                        0.9575     0.0531   18.04  < 2e-16 ***
#> speed50:Mean                 -0.8133     0.1894   -4.29  1.8e-05 ***
#> speed50:St.Dev                0.7099     0.2279    3.12   0.0018 ** 
#> HetMean:speed50:AADT10kplus   1.5787     0.6084    2.60   0.0095 ** 
#> ln(alpha):(Intercept)        -1.2101     0.5325   -2.27   0.0231 *  
#> ln(alpha):lnlength           -0.5031     0.4293   -1.17   0.2413    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1139.169 | AIC: 2292.337 | BIC: 2329.534 
#> N: 1501

# 4. Random Parameters Poisson Model with panel specification
rp_poisson <- countreg_rp(Total_crashes ~ lnaadt,
                      rpar_formula = ~ -1 + speed50,
                      het_mean_formula = ~ AADT10kplus,
                      data = washington_roads,
                      family = "POISSON",
                      rpardists = c(speed50 = "n"),
                      ndraws = 100,
                      group_var = "ID")
summary(rp_poisson)
#> Call:
#> countreg_rp(formula = Total_crashes ~ lnaadt, rpar_formula = ~-1 + 
#>     speed50, data = washington_roads, family = "POISSON", rpardists = c(speed50 = "n"), 
#>     het_mean_formula = ~AADT10kplus, ndraws = 100, group_var = "ID")
#> 
#> Model: countreg_rp | family: POISSON | optimizer: BFGS
#> Iterations: 34 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>                             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)                  -8.6547     0.4108  -21.07  < 2e-16 ***
#> lnaadt                        0.9730     0.0467   20.84  < 2e-16 ***
#> speed50:Mean                 -0.8427     0.1611   -5.23  1.7e-07 ***
#> speed50:St.Dev                0.8015     0.1565    5.12  3.0e-07 ***
#> HetMean:speed50:AADT10kplus   1.1471     0.5484    2.09    0.036 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1161.328 | AIC: 2332.656 | BIC: 2353.798 
#> N: 1501 | panels: 507
# }
```
