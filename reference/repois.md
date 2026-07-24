# Estimate a Random Effects Poisson Panel Model

Estimate a random effects Poisson panel model for clustered count data.
The model is also known as a negative multinomial model.

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

  an R formula.

- group_var:

  the grouping variable(s) for the random effects (e.g., individual ID
  or other panel ID variables).

- data:

  a dataframe that has all of the variables in the `formula`.

- method:

  a method to use for optimization in the maximum likelihood estimation.
  For options, see
  [`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html).

- max.iters:

  the maximum number of iterations to allow the optimization method to
  perform.

- print.level:

  integer specifying the verbosity of output during optimization.

- bootstraps:

  optional integer specifying the number of bootstrap samples to be used
  for estimating standard errors. If not specified, no bootstrapping is
  performed.

- offset:

  an optional offset term provided as a string.

## Value

An object of class `countreg` which is a list with the following
components:

- model: the fitted model object.

- data: the data frame used to fit the model.

- call: the matched call.

- formula: the formula used to fit the model.

## Details

This function estimates a random effects Poisson panel model. The
likelihood is integrated analytically over a gamma-distributed cluster
effect with mean 1 and variance `alpha`. The model estimates the
regression coefficients and one heterogeneity parameter `ln(alpha)`.

## Examples

``` r
# \donttest{
## RE Poisson Panel Model
data("washington_roads")
washington_roads$AADTover10k <-
  ifelse(washington_roads$AADT > 10000, 1, 0)

repois.mod <- repois(
  Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k,
  data = washington_roads,
  offset = "lnlength",
  group_var = "ID",
  method = "BFGS",
  max.iters = 1000
)

summary(repois.mod)
#> Call:
#>  Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k 
#> 
#>  Method:  REPoisson 
#> Iterations:  6 
#> Convergence:  successful convergence  
#> Log-likelihood:  -271.7189 
#> 
#> Parameter Estimates:
#> # A tibble: 7 × 7
#>   parameter           coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>               <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)        -9.82        1.26     -7.80     0        -12.3       -7.35 
#> 2 lnaadt              1.03        0.144     7.14     0          0.745      1.31 
#> 3 speed50            -0.919       0.291    -3.16     0.002     -1.49      -0.349
#> 4 ShouldWidth04      -0.486       0.252    -1.93     0.054     -0.981      0.009
#> 5 AADTover10k        -0.87        0.467    -1.86     0.062     -1.78       0.045
#> 6 lnalpha           -13.8         0.003 -5164.       0        -13.8      -13.8  
#> 7 lnlength (Offset…   1          NA        NA       NA         NA         NA    
# }
```
