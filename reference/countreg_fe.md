# Fixed-effects count regression

Estimates count regression models with fixed effects for one or more
grouping variables (sites, years, ...). Poisson models use the
conditional likelihood (one fixed-effect dimension) or IRLS with
alternating projections (several dimensions); other families, including
families made with
[`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md),
use the dummy-variable likelihood with the group effects concentrated
out.

## Usage

``` r
countreg_fe(
  formula,
  fixef,
  data,
  family = "POISSON",
  offset = NULL,
  weights = NULL,
  method = c("auto", "conditional", "demean", "concentrated"),
  stderr = c("cluster", "hessian"),
  cluster = NULL,
  jackknife = FALSE,
  time = NULL,
  max.iters = 1000,
  tol = 1e-10,
  keep_data = TRUE
)
```

## Arguments

- formula:

  model formula of the mean, `y ~ x1 + x2`, with an optional
  [`offset()`](https://rdrr.io/r/stats/offset.html). The intercept is
  absorbed by the fixed effects.

- fixef:

  a one-sided formula of the fixed-effect variables, e.g. `~ site` or
  `~ site + year`.

- data:

  a data frame.

- family:

  the count family: `"POISSON"` (default) or any other family of
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  (one fixed-effect dimension; enter further dimensions with few levels,
  such as years, as factors in `formula`), or a family made with
  [`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md).

- offset:

  optional name of an offset column (or use
  [`offset()`](https://rdrr.io/r/stats/offset.html) in the formula).

- weights:

  optional name of a column of frequency weights.

- method:

  `"auto"` (default: `"conditional"` for Poisson with one dimension,
  `"demean"` for Poisson with several, `"concentrated"` otherwise), or
  one of those.

- stderr:

  `"cluster"` (default; cluster-robust by `cluster`) or `"hessian"`.

- cluster:

  the clustering variable(s) for `stderr = "cluster"`: a one-sided
  formula, column name(s) or a vector (see
  [`vcov.flexCountReg`](https://jwood-iastate.github.io/flexCountReg/reference/vcov.flexCountReg.md)).
  The default is the first fixed-effect variable. Clusters must contain
  whole groups of that variable.

- jackknife:

  logical: apply the split-panel jackknife bias correction of Dhaene and
  Jochmans (2015) to the slopes (for families other than Poisson in
  short panels). Each group's observations, ordered by `time`, are split
  into a first and second half (sharing the middle period when the
  number of periods is odd); the corrected estimate is \\2\hat\theta -
  (\hat\theta_1 + \hat\theta_2)/2\\, for the slopes and the distribution
  parameters.

- time:

  optional column that orders the observations within groups for
  `jackknife` (default: the row order).

- max.iters:

  maximum number of iterations.

- tol:

  convergence tolerance (relative change of the log-likelihood or of the
  deviance).

- keep_data:

  logical; store the data in the fit (needed for predictions on the
  training rows without `newdata`).

## Value

An object of class `c("countreg_fe", "flexCountReg")`. `model$fixef`
holds the estimated effects of each dimension, `model$dropped` the
groups and rows that were dropped, and `model$not_identified` the
regressors without variation within groups.

## Details

Poisson fixed-effects estimates of the slopes are consistent with fixed
T and need only a correct conditional mean (Wooldridge 1999). For other
families the concentrated (dummy-variable) estimates have an
incidental-parameter bias of order 1/T. For NB2 the bias of the slopes
is small (Allison and Waterman 2002), but the dispersion is strongly
underestimated in short panels. In simulations with 200 groups (NB2,
size 2), the bias of ln(alpha) was -1.08, -0.55 and -0.24 with T = 3, 5
and 10 periods, and slope biases were 1-2\\ `jackknife = TRUE` reduced
the dispersion bias to -0.03 (T = 5) and 0.06 (T = 10) but did not
improve the slopes, and it failed with T = 3 (its half-panels are too
short; a warning is given below 4 periods). The conditional negative
binomial of Hausman, Hall and Griliches is a different model (its fixed
effect enters the dispersion) and is not provided.

Standard errors are cluster-robust by the first fixed-effect variable by
default, with the factor \\G/(G-1)\\ for \\G\\ clusters. They equal
those of
[`fixest::fepois()`](https://lrberge.github.io/fixest/reference/feglm.html)
and
[`fixest::fenegbin()`](https://lrberge.github.io/fixest/reference/femlm.html)
with `ssc(adj = FALSE)`; by default fixest also multiplies the variance
by \\(n-1)/(n-K)\\. The log-likelihood is that of the dummy-variable
model, and its degrees of freedom count the fixed effects, so AIC and
BIC compare with models that estimate them.

[`predict()`](https://rdrr.io/r/stats/predict.html) uses each known
group's estimated effect. A group not in the estimation data has no
fixed effect: its prediction is `NA`, unless `newgroup = "mean_fe"` is
given to use the weighted mean of the estimated effects.

## References

Allison, P.D. and Waterman, R.P. (2002). Fixed-effects negative binomial
regression models. *Sociological Methodology*, 32, 247-265.

Correia, S., Guimaraes, P. and Zylkin, T. (2020). Fast Poisson
estimation with high-dimensional fixed effects. *The Stata Journal*,
20(1), 95-115.

Dhaene, G. and Jochmans, K. (2015). Split-panel jackknife estimation of
fixed-effect models. *Review of Economic Studies*, 82(3), 991-1030.

Wooldridge, J.M. (1999). Distribution-free estimation of some nonlinear
panel data models. *Journal of Econometrics*, 90(1), 77-97.

## Examples

``` r
# \donttest{
data("washington_roads")
fe <- countreg_fe(Total_crashes ~ lnaadt + offset(lnlength),
                  fixef = ~ ID, data = washington_roads)
#> 804 row(s) in 273 group(s) with all-zero counts or a single observation were dropped: they carry no information on the slopes.
summary(fe)
#> Call:
#> countreg_fe(formula = Total_crashes ~ lnaadt + offset(lnlength), 
#>     fixef = ~ID, data = washington_roads)
#> 
#> Model: countreg_fe | family: POISSON | optimizer: conditional
#> Poisson fixed effects (conditional likelihood) 
#> Groups dropped (all-zero counts or one observation): 273 
#> Iterations: 4 
#> Convergence: converged 
#> 
#> Coefficients (standard errors: cluster-robust (234 clusters of ID)):
#>        Estimate Std. Error z value Pr(>|z|)
#> lnaadt   -0.543      0.491    -1.1     0.27
#> Offset: offset(lnlength) (coefficient fixed at 1)
#> 
#> Log-likelihood: -735.347 | AIC: 1940.695 | BIC: 3009.189 
#> N: 697 | panels: 234
# }
```
