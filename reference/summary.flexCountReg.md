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

  Optional parameters that include \`confint_level\` and \`digits\`.

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

Optional parameter \`confint_level\`: A numeric value between 0 and 1
indicating the confidence level for confidence intervals. Default is
0.95.

Optional parameter \`digits\`: Number of digits (decimal places) to
round to. Default is 3.

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
#> (Using bootstrapped standard errors)
#> # A tibble: 7 × 7
#>   parameter           coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>               <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)        -7.40        0.043  -172.       0         -7.49      -7.32 
#> 2 lnaadt              0.912       0.005   182.       0          0.902      0.921
#> 3 lnlength            0.843       0.037    22.9      0          0.771      0.915
#> 4 speed50            -0.47        0.102    -4.61     0         -0.67      -0.27 
#> 5 AADT10kplus         0.77        0.09      8.59     0          0.594      0.945
#> 6 ln(alpha):(Interc… -1.62        0.288    -5.62     0         -2.18      -1.06 
#> 7 ln(alpha):speed50   1.31        0.458     2.85     0.004      0.409      2.20 
# }
```
