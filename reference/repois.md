# Random-Effects Poisson Panel Model

Estimate a random-effects Poisson panel model for clustered event
counts. This function follows the gamma random-effects representation
discussed by Guo (1996), which yields the same likelihood form used for
negative multinomial clustered-count models.

## Usage

``` r
repois(
  formula,
  group_var,
  data,
  method = "NM",
  max.iters = 1000,
  print.level = 0,
  bootstraps = NULL,
  offset = NULL
)
```

## Arguments

- formula:

  An `R` formula.

- group_var:

  The grouping variable(s) for the random effects (for example, an
  individual ID, corridor ID, or other panel identifier).

- data:

  A data frame containing all variables used in `formula`.

- method:

  A method used for maximum-likelihood optimization. The default is
  `"NM"`.

- max.iters:

  The maximum number of optimizer iterations. Default is `1000`.

- print.level:

  Integer controlling optimizer verbosity. Default is `0`.

- bootstraps:

  Optional integer specifying the number of bootstrap replications for
  standard errors. If `NULL`, no bootstrap is used.

- offset:

  Optional offset term provided as a string.

## Value

An object of class `countreg` with components including:

- `model`:

  The fitted model object.

- `data`:

  The data frame used to fit the model.

- `call`:

  The matched call.

- `formula`:

  The formula used to fit the model.

## Details

The function is intended for panel or clustered count data where
repeated observations share an unobserved cluster-specific effect.

This is a thin user-facing wrapper around the existing random-effects
count regression engine used in
[`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md).
The parameterization is presented as a Poisson panel model with a
cluster-level random effect, following the derivation in Guo (1996). The
fitted object is returned as class `countreg`, so it remains compatible
with existing methods such as
[`summary()`](https://rdrr.io/r/base/summary.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`regCompTable()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTable.md),
and
[`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md).

## References

Guo, G. (1996). Negative multinomial regression models for clustered
event counts. *Sociological Methodology*, 26, 113–132.

## See also

[`renb`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
[`predict.flexCountReg`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md),
[`summary.flexCountReg`](https://jwood-iastate.github.io/flexCountReg/reference/summary.flexCountReg.md),
[`regCompTable`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTable.md),
[`regCompTest`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md)

## Examples

``` r
# \donttest{
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT > 10000, 1, 0)

repois.mod <- repois(
  Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k,
  data = washington_roads,
  offset = "lnlength",
  group_var = "ID",
  method = "NM",
  max.iters = 1000
)

summary(repois.mod)
#> Call:
#>  Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k 
#> 
#>  Method:  RENB 
#> Iterations:  746 
#> Convergence:  successful convergence  
#> Log-likelihood:  -263.7988 
#> 
#> Parameter Estimates:
#> # A tibble: 8 × 7
#>   parameter           coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>               <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)        -7.05        0.127   -55.7      0         -7.30      -6.80 
#> 2 lnaadt              0.972       0.015    64.6      0          0.942      1.00 
#> 3 speed50            -0.985       0.301    -3.27     0.001     -1.58      -0.395
#> 4 ShouldWidth04      -0.414       0.212    -1.96     0.051     -0.828      0.001
#> 5 AADTover10k        -0.861       0.447    -1.92     0.054     -1.74       0.015
#> 6 ln(a)               3.01        0.115    26.2      0          2.78       3.23 
#> 7 ln(b)               0.604       0.122     4.97     0          0.366      0.843
#> 8 lnlength (Offset …  1          NA        NA       NA         NA         NA    
# }
```
