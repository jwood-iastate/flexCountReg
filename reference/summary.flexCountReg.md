# Custom summary method for flexCountReg models

Custom summary method for flexCountReg models

## Usage

``` r
# S3 method for class 'flexCountReg'
summary(object, ...)
```

## Arguments

- object:

  A flexCountReg model object.

- ...:

  Optional parameters that include `confint_level` and `digits`.

## Value

Prints the model formula, method used for estimation, number of
iterations used, if the model converged, and the log-likelihood. Then,
it prints a table containing parameter estimates, standard errors,
t-statistics, p-values, and confidence intervals. Also quietly returns a
tibble with these values.

## Details

This summary method accounts for bootstrapped or robust standard errors
(when used).

## Note

Optional parameter `confint_level`: A numeric value between 0 and 1
indicating the confidence level for confidence intervals. Default is
0.95.

Optional parameter `digits`: Number of digits (decimal places) to round
to. Default is 3.

## Examples

``` r
# \donttest{
# NB2 Model
data("washington_roads")
washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
nb2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                data = washington_roads, family = "NB2",
                dis_param_formula_1 = ~ speed50, method='BFGS')
summary(nb2)
#> Call:
#>  Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus 
#> 
#>  Method:  countreg 
#> Iterations:  44 
#> Convergence:  successful convergence  
#> Log-likelihood:  -1064.876 
#> 
#> Parameter Estimates:
#> # A tibble: 7 × 7
#>   parameter           coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>               <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)        -7.40        0.476   -15.5      0         -8.34      -6.47 
#> 2 lnaadt              0.912       0.057    16.0      0          0.8        1.02 
#> 3 lnlength            0.843       0.07     12.0      0          0.706      0.98 
#> 4 speed50            -0.47        0.115    -4.08     0         -0.695     -0.244
#> 5 AADT10kplus         0.77        0.13      5.92     0          0.515      1.02 
#> 6 ln(alpha):(Interc… -1.62        0.369    -4.39     0         -2.34      -0.896
#> 7 ln(alpha):speed50   1.31        0.575     2.27     0.023      0.178      2.43 
# }
```
