# Random Parameters Count Regression Models

Random Parameters Count Regression Models

## Usage

``` r
countreg.rp(
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
  method = "BHHH",
  max.iters = 1000,
  start.vals = NULL,
  verbose = FALSE,
  stderr = "normal",
  bootstraps = NULL,
  engine = "halton",
  keep_data = TRUE,
  panel_id = NULL
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
  (See
  [`countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  for full list).

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

  if the Halton draws should be scrambled.

- correlated:

  if the random parameters should be correlated. If TRUE, only normal
  distributions are used.

- group_var:

  an optional variable name (string), or a vector with one value per
  row, defining the panel structure (repeated measures). If provided,
  each panel shares one draw of the random parameters, and the
  likelihood and standard errors account for the panel structure.

- weights:

  the name of a variable, or a numeric vector, of frequency weights.
  With `panel_id`, weights apply to whole panels (the panel
  log-likelihood is weighted), so they should be constant within a
  panel.

- offset:

  the name of a variable, or vector of variable names, to be used as an
  offset (summed when several are given). An offset can instead be given
  with [`offset()`](https://rdrr.io/r/stats/offset.html) in `formula`;
  use only one of the two.

- method:

  optimization method (e.g., "BHHH", "BFGS", "NM").

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
  model).

- bootstraps:

  the number of bootstrap samples when `stderr = "boot"` (default 200).

- engine:

  the engine to use for the integration if the Poisson Lognormal (PLN)
  distribution is used. Options include "halton" which uses Halton draws
  or "poilog" which uses the poilog package. The default is "halton".

- keep_data:

  logical; if `FALSE`, the training data are not stored in the fitted
  object. [`predict()`](https://rdrr.io/r/stats/predict.html) then needs
  `newdata`.

- panel_id:

  deprecated; use `group_var`.

## Value

An object of class `flexCountReg` which is a list with the following
components:

- model: the fitted model object.

- data: the data frame used to fit the model.

- call: the matched call.

- formula: the formula used to fit the model.

## Details

`countreg_rp()` is the same function as `countreg.rp()`, under a name
that tools do not mistake for an S3 method.

**Random parameter distributions.** Each random parameter is \\\beta =
\mu + \sigma \upsilon\\, where the `:Mean` and `:St.Dev` estimates are
\\\mu\\ and \\\sigma\\:

- normal ("n"): \\\upsilon \sim N(0, 1)\\, so \\\sigma\\ is the standard
  deviation.

- uniform ("u") and triangular ("t"): \\\upsilon\\ is uniform or
  symmetric triangular on \\\[-1, 1\]\\, so \\\beta\\ lies in \\\[\mu -
  \sigma, \mu + \sigma\]\\. Here \\\sigma\\ is a half-width; the
  standard deviation is \\\sigma/\sqrt{3}\\ (uniform) or
  \\\sigma/\sqrt{6}\\ (triangular). Models fitted before version 0.1.3
  used \\\[\mu - \sigma/2, \mu + \sigma/2\]\\ for the uniform; their
  \\\sigma\\ is converted automatically in
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

## References

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
rp_nb2 <- countreg.rp(Total_crashes ~ lnaadt + lnlength,
                      rpar_formula = ~ -1 + speed50,
                      data = washington_roads,
                      family = "NB2",
                      rpardists = c(speed50 = "n"),
                      ndraws = 100,
                      method = "BHHH")
summary(rp_nb2)
#> Call:
#> countreg.rp(formula = Total_crashes ~ lnaadt + lnlength, rpar_formula = ~-1 + 
#>     speed50, data = washington_roads, family = "NB2", rpardists = c(speed50 = "n"), 
#>     ndraws = 100, method = "BHHH")
#> 
#> Model: countreg_rp | family: NB2 | optimizer: BHHH
#> Iterations: 29 
#> Convergence: successive function values within relative tolerance limit (reltol) 
#> 
#> Coefficients (standard errors: Hessian-based):
#>                Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)     -8.7552     0.4133  -21.18  < 2e-16 ***
#> lnaadt           1.0799     0.0500   21.58  < 2e-16 ***
#> lnlength         0.7564     0.0693   10.92  < 2e-16 ***
#> speed50:Mean    -0.7445     0.1687   -4.41  1.0e-05 ***
#> speed50:St.Dev   0.7170     0.2705    2.65    0.008 ** 
#> ln(alpha)       -1.3182     0.3052   -4.32  1.6e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1082.98 | AIC: 2177.961 | BIC: 2209.844 
#> N: 1501

# 2. Random Parameters with Panel Structure (if 'site_id' exists)
# rp_panel <- countreg.rp(Total_crashes ~ -1 + lnaadt,
#                        rpar_formula = ~ speed50,
#                        data = washington_roads,
#                        group_var = "site_id",
#                        family = "NB2",
#                        ndraws = 100)

# 3. Generalized Random Parameters Model with Heterogeneity
rp_gen <- countreg.rp(Total_crashes ~ lnaadt,
                      rpar_formula = ~ -1 + speed50,
                      dis_param_formula_1 = ~ lnlength, 
                      het_mean_formula = ~ AADT10kplus,
                      data = washington_roads,
                      family = "NB2",
                      rpardists = c(speed50 = "n"),
                      ndraws = 100)
summary(rp_gen)
#> Call:
#> countreg.rp(formula = Total_crashes ~ lnaadt, rpar_formula = ~-1 + 
#>     speed50, data = washington_roads, family = "NB2", rpardists = c(speed50 = "n"), 
#>     dis_param_formula_1 = ~lnlength, het_mean_formula = ~AADT10kplus, 
#>     ndraws = 100)
#> 
#> Model: countreg_rp | family: NB2 | optimizer: BHHH
#> Iterations: 18 
#> Convergence: successive function values within relative tolerance limit (reltol) 
#> 
#> Coefficients (standard errors: Hessian-based):
#>                             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)                  -8.4758     0.4390  -19.31   <2e-16 ***
#> lnaadt                        0.9559     0.0518   18.46   <2e-16 ***
#> speed50:Mean                 -0.7810     0.1872   -4.17    3e-05 ***
#> speed50:St.Dev                0.7417     0.2514    2.95   0.0032 ** 
#> HetMean:speed50:AADT10kplus   1.5599     1.2085    1.29   0.1968    
#> ln(alpha):(Intercept)        -1.1705     0.3638   -3.22   0.0013 ** 
#> ln(alpha):lnlength           -0.4738     0.2842   -1.67   0.0955 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1139.336 | AIC: 2292.673 | BIC: 2329.87 
#> N: 1501

# 4. Random Parameters Poisson Model with panel specification
rp_poisson <- countreg.rp(Total_crashes ~ lnaadt,
                      rpar_formula = ~ -1 + speed50,
                      het_mean_formula = ~ AADT10kplus,
                      data = washington_roads,
                      family = "POISSON",
                      rpardists = c(speed50 = "n"),
                      ndraws = 100,
                      group_var = "ID")
summary(rp_poisson)
#> Call:
#> countreg.rp(formula = Total_crashes ~ lnaadt, rpar_formula = ~-1 + 
#>     speed50, data = washington_roads, family = "POISSON", rpardists = c(speed50 = "n"), 
#>     het_mean_formula = ~AADT10kplus, ndraws = 100, group_var = "ID")
#> 
#> Model: countreg_rp | family: POISSON | optimizer: BHHH
#> Iterations: 11 
#> Convergence: successive function values within relative tolerance limit (reltol) 
#> 
#> Coefficients (standard errors: Hessian-based):
#>                             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)                  -8.6482     0.2933  -29.48  < 2e-16 ***
#> lnaadt                        0.9723     0.0328   29.64  < 2e-16 ***
#> speed50:Mean                 -0.8228     0.1555   -5.29  1.2e-07 ***
#> speed50:St.Dev                0.8465     0.1844    4.59  4.4e-06 ***
#> HetMean:speed50:AADT10kplus   1.0877     0.5333    2.04    0.041 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1161.033 | AIC: 2332.066 | BIC: 2353.208 
#> N: 1501 | panels: 507
# }
```
