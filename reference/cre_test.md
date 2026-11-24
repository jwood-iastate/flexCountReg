# Test of correlated random effects (random against fixed effects)

Wald test that the coefficients of the group means added by `cre` are
all zero. Under the null hypothesis the random effects are uncorrelated
with the time-varying regressors, so the random-effects estimates are
consistent; a rejection favors fixed effects (or keeping the means in
the model). With robust or cluster-robust standard errors the test is a
robust version of the Hausman test (Wooldridge 2019).

## Usage

``` r
cre_test(object, vcov = NULL)
```

## Arguments

- object:

  a model fitted with `cre`
  ([`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md),
  [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)
  or
  [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md)).

- vcov:

  optional covariance matrix of the estimates (default: the fit's,
  `vcov(object)`).

## Value

An object of class `"htest"`.

## References

Mundlak, Y. (1978). On the pooling of time series and cross section
data. *Econometrica*, 46(1), 69-85.

Wooldridge, J.M. (2019). Correlated random effects models with
unbalanced panels. *Journal of Econometrics*, 211(1), 137-150.

## Examples

``` r
# \donttest{
data("washington_roads")
m <- repois(Total_crashes ~ lnaadt + offset(lnlength), group_var = "ID",
            data = washington_roads, cre = ~ lnaadt, stderr = "robust")
cre_test(m)
#> 
#>  Wald test of correlated random effects (group means)
#> 
#> data:  repois(formula = Total_crashes ~ lnaadt + offset(lnlength), group_var = "ID",      data = washington_roads, stderr = "robust", cre = ~lnaadt)
#> Wald chi-squared = 7.8392, df = 1, p-value = 0.005113
#> sample estimates:
#> mean_lnaadt 
#>    1.861682 
#> 
# }
```
