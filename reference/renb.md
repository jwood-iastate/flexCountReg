# Estimate a Random Effects Negative Binomial regression model

Estimate a Random Effects Negative Binomial regression model

## Usage

``` r
renb(
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
  [`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html). Note that
  "BHHH" is not available for this function due to the implementation
  for the random effects.

- max.iters:

  the maximum number of iterations to allow the optimization method to
  perform.

- print.level:

  Integer specifying the verbosity of output during optimization.

- bootstraps:

  Optional integer specifying the number of bootstrap samples to be used
  for estimating standard errors. If not specified, no bootstrapping is
  performed.

- offset:

  an optional offset term provided as a string.

## Details

This function estimates a random effects negative binomial (RENB)
regression model. This model is based on the NB-1 model. The PDF for the
RENB is: \$\$f(y\_{it}\|\mu\_{it}, a, b) = \frac{\Gamma(a+b) +
\Gamma(a + \sum\_{t = 1}^{n_i} \mu\_{it}) + \Gamma(b +
\sum\_{t=1}^{n_i}y\_{it})} {\Gamma(a) \Gamma(b) \Gamma(a + b +
\sum\_{t=1}^{n_i}\mu\_{it} + \sum\_{t=1}^{n_i}y\_{it})}
\prod\_{t=1}^{n_i}
\frac{\Gamma(\mu\_{it}+y\_{it})}{\Gamma(\mu\_{it})\Gamma(y\_{it})}\$\$

## Examples

``` r
# \donttest{
## RENB Model
data("washington_roads")
washington_roads$AADTover10k <- 
  ifelse(washington_roads$AADT > 10000, 1, 0) # create a dummy variable
renb.mod <- renb(Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k,
                                data=washington_roads,
                                offset = "lnlength",
                                group_var="ID",
                                method="nm",
                                max.iters = 1000)
summary(renb.mod)
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
