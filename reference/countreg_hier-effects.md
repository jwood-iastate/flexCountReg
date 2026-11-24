# Fixed effects, random effects and variance components

[`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html) returns the
fixed coefficients of the mean,
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) the
conditional modes of the random effects (one data frame per grouping
factor, with their posterior covariances in the attribute `"postVar"`),
and [`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) the
estimated covariance matrices of the random effects.

## Usage

``` r
# S3 method for class 'countreg_hier'
fixef(object, ...)

# S3 method for class 'countreg_hier'
ranef(object, condVar = TRUE, ...)

# S3 method for class 'countreg_hier'
VarCorr(x, ...)

# S3 method for class 'VarCorr.countreg_hier'
print(x, digits = 4, ...)
```

## Arguments

- object:

  a
  [`countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md)
  fit.

- ...:

  unused.

- condVar:

  include the posterior covariances.

- x:

  a [`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) result.

- digits:

  number of significant digits to print.

## Value

[`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html): a named
vector. [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html): a
named list of data frames.
[`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html): a named list
of covariance matrices with attributes `"stddev"` and `"correlation"`.

## Examples

``` r
# \donttest{
data("washington_roads")
h <- countreg_hier(Total_crashes ~ lnaadt + offset(lnlength) + (1 | ID),
                   data = washington_roads, family = "POISSON")
fixef(h)
#> (Intercept)      lnaadt 
#>   -9.445363    1.148339 
head(ranef(h)$ID)
#>   (Intercept)
#> 1  -0.4459850
#> 2   0.4177328
#> 3  -0.4282198
#> 4   0.5740780
#> 5  -0.4678016
#> 6   0.1352786
VarCorr(h)
#>  Groups Name        Variance Std.Dev.
#>  ID     (Intercept) 0.4658   0.6825  
# }
```
