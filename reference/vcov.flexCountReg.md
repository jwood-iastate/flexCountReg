# Covariance matrix of the estimates

[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the covariance
matrix chosen at fit time, or recomputes a Hessian-based, robust or
cluster-robust one from the stored scores, without refitting.

## Usage

``` r
# S3 method for class 'flexCountReg'
vcov(object, type = NULL, cluster = NULL, ...)
```

## Arguments

- object:

  a fitted `flexCountReg` model.

- type:

  `NULL` or `"fit"` (default): the covariance chosen at fit time
  (`stderr`). `"hessian"`: the inverse of the negative Hessian.
  `"robust"`: the sandwich with one independent unit per row, or per
  panel for panel models. For
  [`countreg_fe()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_fe.md)
  fits by the conditional or concentrated likelihood the units are the
  groups of the fixed effect (their scores do not split into rows), so
  `"robust"` is clustered by group, without the factor \\G/(G-1)\\; with
  the `"demean"` method the units are the rows (as `vcov = "hetero"` in
  fixest). `"cluster"`: cluster-robust with the clusters in `cluster`.

- cluster:

  for `type = "cluster"`: a one-sided formula of the clustering
  variable(s) in the model's data (`~ corridor`; two or more variables,
  `~ corridor + year`, give multiway clustering), their names, or a
  vector with one cluster id per row used in estimation.

- ...:

  not used.

## Value

A square matrix with the parameter names as dimnames.

## Details

Cluster-robust covariances use the scores of the model's independent
units (rows, or panels for
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
with `group_var`,
[`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
[`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md),
[`poislind_re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)),
summed within each cluster, with the small-sample factor \\G/(G-1)\\ for
\\G\\ clusters. Clusters may be coarser than the panels (for example
sites within corridors) but may not split a panel. Multiway clustering
follows Cameron, Gelbach and Miller (2011). With frequency weights, a
row with weight \\w\\ stands for \\w\\ observations in the same cluster.

The number of clusters is attached as the attribute `"n_clusters"`. The
scores and the bread are also available to the sandwich package through
[`estfun.flexCountReg`](https://jwood-iastate.github.io/flexCountReg/reference/estfun.flexCountReg.md).

## References

Cameron, A.C., Gelbach, J.B. and Miller, D.L. (2011). Robust inference
with multiway clustering. *Journal of Business and Economic Statistics*,
29(2), 238-249.

## Examples

``` r
# \donttest{
data("washington_roads")
m <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads,
              family = "NB2", method = "BFGS")
sqrt(diag(vcov(m)))
#> (Intercept)      lnaadt    lnlength   ln(alpha) 
#>  0.44451075  0.05291691  0.06960363  0.23366292 
sqrt(diag(vcov(m, type = "cluster", cluster = ~ ID)))
#> (Intercept)      lnaadt    lnlength   ln(alpha) 
#>  0.58639184  0.06790319  0.08712895  0.24190144 
# }
```
