# Summarize a flexCountReg model

Summarize a flexCountReg model

## Usage

``` r
# S3 method for class 'flexCountReg'
summary(object, level = 0.95, ...)

# S3 method for class 'summary.flexCountReg'
print(x, digits = x$digits, ...)

# S3 method for class 'summary.flexCountReg'
as_tibble(x, ...)
```

## Arguments

- object:

  A flexCountReg model object.

- level:

  the confidence level for the confidence intervals. Default is 0.95.

- ...:

  Optional arguments. The arguments `confint_level` and `digits` of
  earlier versions are still accepted.

- x:

  a `summary.flexCountReg` object.

- digits:

  the number of digits (decimal places) to round to when printing.
  Default is 3.

## Value

An object of class `summary.flexCountReg`, a list with `coefficients` (a
matrix with the estimate, standard error, z value, p-value and
confidence limits of every parameter), `vcov_type`, `logLik`, `AIC`,
`BIC`, `nobs`, `n_groups`, `n_dropped`, `converged`, `offset` and
`lr_test`. Printing it shows the results. Use
[`tibble::as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
on it for the coefficient table as a tibble (the value returned by
[`summary()`](https://rdrr.io/r/base/summary.html) in earlier versions),
or
[`tidy.flexCountReg`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-tidiers.md).

## Details

The standard errors are those chosen when the model was fit:
Hessian-based, robust (sandwich), or bootstrapped. The z statistics and
p-values use the normal distribution.

For
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
models of families that nest the Poisson model with a regular
restriction, a likelihood-ratio test against a Poisson model with the
same mean specification, rows, offset and weights is reported. For NB1,
NB2, PLN, PIG1, PIG2 and PIG with a scalar dispersion parameter the
Poisson model is on the boundary of the parameter space, and the p-value
uses a 50:50 mixture of chi-square distributions with 0 and 1 degrees of
freedom; for GP1, GP2 and COM the restriction is interior. No test is
reported for PL, PLG and PLL (which do not reduce to Poisson), for NBP,
SI, PW, PGE and GW (nonregular restrictions), for boundary families with
a dispersion formula, or for models with zero or underreporting
components.

## Examples

``` r
# \donttest{
data("washington_roads")
nb2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50,
                data = washington_roads, family = "NB2",
                dis_param_formula_1 = ~ speed50, method = "BFGS")
s <- summary(nb2)
s
#> Call:
#> countreg(formula = Total_crashes ~ lnaadt + lnlength + speed50, 
#>     data = washington_roads, family = "NB2", dis_param_formula_1 = ~speed50, 
#>     method = "BFGS")
#> 
#> Model: countreg | family: NB2 | optimizer: BFGS
#> Iterations: 47 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>                       Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)            -8.7564     0.4404  -19.88  < 2e-16 ***
#> lnaadt                  1.0803     0.0518   20.86  < 2e-16 ***
#> lnlength                0.7588     0.0691   10.98  < 2e-16 ***
#> speed50                -0.5407     0.1173   -4.61  4.0e-06 ***
#> ln(alpha):(Intercept)  -1.3233     0.3150   -4.20  2.7e-05 ***
#> ln(alpha):speed50       1.2670     0.4849    2.61    0.009 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1082.176 | AIC: 2176.352 | BIC: 2208.235 
#> N: 1501
s$coefficients
#>                         Estimate Std. Error    z value     Pr(>|z|)      2.5 %
#> (Intercept)           -8.7563978 0.44040867 -19.882437 5.775863e-88 -9.6195830
#> lnaadt                 1.0803267 0.05178872  20.860270 1.229765e-96  0.9788227
#> lnlength               0.7588171 0.06909183  10.982731 4.627175e-28  0.6233996
#> speed50               -0.5406739 0.11728001  -4.610112 4.024529e-06 -0.7705385
#> ln(alpha):(Intercept) -1.3232832 0.31500069  -4.200890 2.658677e-05 -1.9406732
#> ln(alpha):speed50      1.2669833 0.48492203   2.612757 8.981519e-03  0.3165536
#>                           97.5 %
#> (Intercept)           -7.8932127
#> lnaadt                 1.1818308
#> lnlength               0.8942346
#> speed50               -0.3108093
#> ln(alpha):(Intercept) -0.7058932
#> ln(alpha):speed50      2.2174130
tibble::as_tibble(s)
#> # A tibble: 6 × 7
#>   parameter          coeff `Std. Err.` `z value` `p-value` `lower CI` `upper CI`
#>   <chr>              <dbl>       <dbl>     <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)       -8.76        0.44     -19.9      0         -9.62      -7.89 
#> 2 lnaadt             1.08        0.052     20.9      0          0.979      1.18 
#> 3 lnlength           0.759       0.069     11.0      0          0.623      0.894
#> 4 speed50           -0.541       0.117     -4.61     0         -0.771     -0.311
#> 5 ln(alpha):(Inter… -1.32        0.315     -4.20     0         -1.94      -0.706
#> 6 ln(alpha):speed50  1.27        0.485      2.61     0.009      0.317      2.22 
# }
```
