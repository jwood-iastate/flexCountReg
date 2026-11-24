# Check the basis dimension of smooth terms

For each smooth, compares the residuals of neighboring observations
(ordered by the smooth's covariates) with the residual variance. If the
basis dimension `k` is too small, the residuals keep a pattern in the
covariate, so neighbors are more alike than the overall variance
suggests and the k-index falls below 1.

## Usage

``` r
smooth_check(object, term = NULL, nperm = 400, seed = 1)
```

## Arguments

- object:

  a fitted model with smooth terms.

- term:

  optional term labels; by default every smooth.

- nperm:

  number of random permutations of the residuals for the p-value.

- seed:

  seed for the residual randomization and the permutations (the caller's
  random-number stream is left unchanged).

## Value

A data frame with `parameter`, `term`, `k` (the basis dimension after
the identifiability constraint), `edf`, `k_index`, `p_value` and the
residual type.

## Details

The k-index is \\\frac{1}{2}\overline{(r_i - r\_{(i)})^2} /
\mathrm{var}(r)\\, where \\r\_{(i)}\\ is the residual of the nearest
neighbor of observation \\i\\ in the smooth's (scaled) covariates; for
one covariate, the next observation in its order. The residuals are the
randomized quantile residuals of the model, which are N(0, 1) for a
correct count model whatever the family (Pearson residuals where they
are not available). The p-value is the share of residual permutations
with a smaller k-index: a small p-value with a k-index below 1 suggests
increasing `k` (refit and compare the EDF and the curve). Smooths of
factors only (random effects, Markov random fields) are not checked.
This follows
[`mgcv::k.check()`](https://rdrr.io/pkg/mgcv/man/k.check.html), adapted
to count residuals.

## Examples

``` r
# \donttest{
data("washington_roads")
m <- countreg(Total_crashes ~ s(lnaadt, k = 5) + offset(lnlength),
              data = washington_roads, family = "NB2")
smooth_check(m)
#>   parameter      term k      edf   k_index p_value residuals
#> 1      mean s(lnaadt) 4 3.520246 0.9383687  0.0075  quantile
# }
```
