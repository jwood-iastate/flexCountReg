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
  panel_id = NULL,
  weights = NULL,
  offset = NULL,
  method = "BHHH",
  max.iters = 1000,
  start.vals = NULL,
  verbose = FALSE
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
  ("g"). If not provided, normal is used.

- dis_param_formula_1:

  a symbolic description of the model for the first parameter of the
  count distribution (e.g., ln(alpha) for NB2).

- dis_param_formula_2:

  a symbolic description of the model for the second parameter of the
  count distribution (if applicable).

- het_mean_formula:

  an optional symbolic description of the model for heterogeneity in the
  means of the random parameters.

- het_var_formula:

  an optional symbolic description of the model for heterogeneity in the
  variances of the random parameters.

- ndraws:

  the number of Halton draws to use for estimating the random
  parameters.

- scrambled:

  if the Halton draws should be scrambled.

- correlated:

  if the random parameters should be correlated. If TRUE, only normal
  distributions are used.

- panel_id:

  an optional variable name (string) or vector defining the panel
  structure (repeated measures). If provided, the standard errors and
  likelihood are estimated accounting for the panel structure.

- weights:

  variable name to be used as frequency weights.

- offset:

  variable name to be used as an offset.

- method:

  optimization method (e.g., "BHHH", "BFGS", "NM").

- max.iters:

  maximum number of iterations.

- start.vals:

  optional vector of starting values.

- verbose:

  logical.

## Value

An object of class \`countreg\` which is a list with the following
components:

- model: the fitted model object.

- data: the data frame used to fit the model.

- call: the matched call.

- formula: the formula used to fit the model.

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
#>  Total_crashes ~ lnaadt + lnlength 
#> 
#>  Method:  countreg.rp 
#> Iterations:  28 
#> Convergence:  successive function values within relative tolerance limit (reltol) 
#> Log-likelihood:  -1082.98 
#> 
#> Parameter Estimates:
#> # A tibble: 6 × 7
#>   parameter       coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>           <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)    -8.76        0.044  -199.           0     -8.84      -8.67 
#> 2 lnaadt          1.08        0.005   213.           0      1.07       1.09 
#> 3 lnlength        0.756       0.035    21.7          0      0.688      0.825
#> 4 speed50:Mean   -0.744       0.106    -7.04         0     -0.951     -0.537
#> 5 speed50:St.Dev  0.715       0.177     4.04         0      0.368      1.06 
#> 6 ln(alpha)      -1.32        0.293    -4.5          0     -1.89      -0.744

# 2. Random Parameters with Panel Structure (if 'site_id' exists)
# rp_panel <- countreg.rp(Total_crashes ~ -1 + lnaadt,
#                        rpar_formula = ~ speed50,
#                        data = washington_roads,
#                        panel_id = "site_id",
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
#>  Total_crashes ~ lnaadt 
#> 
#>  Method:  countreg.rp 
#> Iterations:  11 
#> Convergence:  successive function values within relative tolerance limit (reltol) 
#> Log-likelihood:  -1140.142 
#> 
#> Parameter Estimates:
#> # A tibble: 7 × 7
#>   parameter           coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>               <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)        -8.54      4.7 e-2  -182.       0       -8.64e+0   -8.45e+0
#> 2 lnaadt              0.965     6   e-3   174.       0        9.54e-1    9.76e-1
#> 3 speed50:Mean       -0.834     1.15e-1    -7.28     0       -1.06e+0   -6.1 e-1
#> 4 speed50:St.Dev      0.83      1.52e-1     5.46     0        5.32e-1    1.13e+0
#> 5 HetMean:AADT10kp… -18.4       7.00e+7     0        1       -1.37e+8    1.37e+8
#> 6 ln(alpha):(Inter…  -1.21      1.92e-1    -6.29     0       -1.58e+0   -8.31e-1
#> 7 ln(alpha):lnleng…  -0.49      1.45e-1    -3.38     0.001   -7.74e-1   -2.06e-1

# 4. Random Parameters Poisson Model with panel specification
rp_poisson <- countreg.rp(Total_crashes ~ lnaadt,
                      rpar_formula = ~ -1 + speed50,
                      dis_param_formula_1 = ~ lnlength, 
                      het_mean_formula = ~ AADT10kplus,
                      data = washington_roads,
                      family = "POISSON",
                      rpardists = c(speed50 = "n"),
                      ndraws = 100,
                      panel_id = "ID")
summary(rp_poisson)
#> Call:
#>  Total_crashes ~ lnaadt 
#> 
#>  Method:  countreg.rp 
#> Iterations:  14 
#> Convergence:  successive function values within relative tolerance limit (reltol) 
#> Log-likelihood:  -1161.137 
#> 
#> Parameter Estimates:
#> # A tibble: 5 × 7
#>   parameter           coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>               <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)        -8.66        0.025  -342.           0     -8.71      -8.61 
#> 2 lnaadt              0.974       0.003   342.           0      0.968      0.98 
#> 3 speed50:Mean       -0.838       0.124    -6.74         0     -1.08      -0.594
#> 4 speed50:St.Dev      0.869       0.151     5.76         0      0.574      1.16 
#> 5 HetMean:AADT10k… -110.         NA        NA           NA     NA         NA    
# }
```
