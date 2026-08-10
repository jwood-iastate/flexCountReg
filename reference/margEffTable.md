# Marginal Effects, Elasticities, and Pseudo-Elasticities for flexCountReg Models

Compute average marginal effects (AMEs), elasticities, or
pseudo-elasticities for fitted `flexCountReg` objects. Standard errors
are computed using either the delta method or a bootstrap, and the
results can be returned as a tibble, a gt table, or a LaTeX table.

## Usage

``` r
margEffTable(
  object,
  data = NULL,
  vars = NULL,
  measure = c("auto", "ame", "elasticity"),
  indicator = c("pseudo_elasticity", "discrete_change"),
  se = c("delta", "bootstrap"),
  bootstraps = 200,
  pred_method = NULL,
  confint_level = 0.95,
  tableType = c("tibble", "gt", "latex"),
  digits = 3,
  cluster_var = NULL,
  ...
)
```

## Arguments

- object:

  A fitted `flexCountReg` object.

- data:

  Optional data frame. Defaults to the data used to fit the model.

- vars:

  Optional character vector of variable names to evaluate. If `NULL`,
  the function uses the non-response variables appearing in the model
  formula that are present in `data`.

- measure:

  Character scalar. One of:

  "auto"

  :   Continuous variables receive AMEs; indicator variables receive
      pseudo-elasticities (or discrete changes if
      `indicator = "discrete_change"`).

  "ame"

  :   Compute average marginal effects for continuous variables.

  "elasticity"

  :   Compute elasticities for continuous variables. Indicator variables
      receive pseudo-elasticities (or discrete changes if
      `indicator = "discrete_change"`).

- indicator:

  Character scalar. For indicator variables, return
  `"pseudo_elasticity"` (default) or `"discrete_change"`.

- se:

  Character scalar. Standard error method: `"delta"` (default) or
  `"bootstrap"`.

- bootstraps:

  Integer. Number of bootstrap replications when `se = "bootstrap"`.

- pred_method:

  Optional prediction method forwarded to
  [`predict.flexCountReg()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md).
  For random-parameters models, this is commonly `"Exact"` or
  `"Simulated"`. If `NULL`, the function uses `"Exact"` for
  random-parameters models when possible.

- confint_level:

  Numeric scalar between 0 and 1. Confidence level for the confidence
  interval. Default is 0.95.

- tableType:

  Character scalar. One of `"tibble"`, `"gt"`, or `"latex"`.

- digits:

  Integer. Number of digits to round in the returned table.

- cluster_var:

  Optional character vector giving a clustering variable for bootstrap
  resampling. If supplied, bootstrap samples are drawn at the cluster
  level. If omitted, the function attempts row-level bootstrap
  resampling.

- ...:

  Additional arguments passed to
  [`predict.flexCountReg()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md).

## Value

A table object of the type requested by `tableType`. For `"tibble"`, the
returned object also carries an `effect_info` attribute with metadata
about the calculation.

## Details

Continuous variables are handled by numerical differentiation on the
expected-response scale. Indicator variables are handled by discrete
changes (or pseudo-elasticities, depending on `indicator`).

The function works on the expected-response scale and is intended to be
compatible with all fitted `flexCountReg` model types, provided that
[`predict.flexCountReg()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md)
can generate predictions for the model.

For continuous regressors, the average marginal effect is computed using
a central finite difference: \$\$ \text{AME}\_j =
\frac{1}{n}\sum\_{i=1}^n \frac{\mu_i(x\_{ij}+h)-\mu_i(x\_{ij}-h)}{2h}
\$\$ where \\\mu_i\\ is the predicted mean response.

For elasticities, the effect is computed as: \$\$ E_j =
\frac{1}{n}\sum\_{i=1}^n \left(\frac{x\_{ij}}{\mu_i}\frac{\partial
\mu_i}{\partial x\_{ij}}\right) \$\$

For indicator variables, the function computes a discrete change between
the observed baseline and a switched value. If
`indicator = "pseudo_elasticity"`, the result is reported as a percent
change relative to the baseline.

## Examples

``` r
# \donttest{
data("washington_roads")
washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)

nb2 <- countreg(
  Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
  data = washington_roads,
  family = "NB2"
)

margEffTable(nb2, tableType = "tibble")
#> Marginal Effects / Elasticities
#> Method: delta 
#> 
#>         term variable_type           effect_metric estimate std_error t_value
#>       lnaadt    continuous Average marginal effect    0.423     0.034  12.551
#>     lnlength    continuous Average marginal effect    0.394     0.039  10.203
#>      speed50     indicator   Pseudo-elasticity (%)  -36.701     6.759  -5.430
#>  AADT10kplus     indicator   Pseudo-elasticity (%)  124.032    29.632   4.186
#>  p_value sig lower_ci upper_ci n_obs successful_bootstraps
#>        0 ***    0.357    0.490  1501                    NA
#>        0 ***    0.318    0.469  1501                    NA
#>        0 ***  -49.949  -23.454  1501                    NA
#>        0 ***   65.955  182.109  1501                    NA
margEffTable(nb2, tableType = "gt")


  


Marginal Effects / Elasticities
```
