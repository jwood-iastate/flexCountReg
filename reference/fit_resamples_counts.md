# Resampling count models, including models with several outcomes

Fits a workflow or model specification to the analysis set of each
resample and evaluates it on the assessment set, for every outcome.
tune's `fit_resamples()` computes metrics for one outcome only, so
models with several
([`countreg_mv`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mv.md),
[`countreg_alloc`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_alloc.md))
need this function. It also computes the out-of-sample log score of the
count distribution, which compares distributions better than errors of
the mean, for any engine.

## Usage

``` r
fit_resamples_counts(
  object,
  resamples,
  formula = NULL,
  metrics = NULL,
  log_score = FALSE,
  summarize = TRUE
)
```

## Arguments

- object:

  a workflows workflow or a
  [`count_reg`](https://jwood-iastate.github.io/flexCountReg/reference/count_reg.md)
  specification.

- resamples:

  an rsample resampling set, such as
  [`rsample::group_vfold_cv()`](https://rsample.tidymodels.org/reference/group_vfold_cv.html).

- formula:

  the model formula, for a model specification.

- metrics:

  a yardstick metric set of numeric metrics (default RMSE and MAE).

- log_score:

  logical; also compute the mean log probability of the held-out counts
  (`"log_score"`, larger is better). For models with several outcomes it
  is the score of each outcome's marginal distribution.

- summarize:

  logical; average over the resamples (`TRUE`) or return the values of
  each?

## Value

A tibble with columns `outcome`, `.metric`, `mean`, `n` (resamples with
a value) and `std_err`, or with `summarize = FALSE` `id`, `outcome`,
`.metric` and `.estimate`.

## Details

Rows without a prediction (such as new groups of a fixed-effects model)
are left out of the metrics. The options `marginal`, `re.form`,
`allow.new.levels` and `newgroup` of the engine argument `predict_args`
also apply to the log score.

## Examples

``` r
# \donttest{
if (requireNamespace("parsnip", quietly = TRUE) &&
    requireNamespace("rsample", quietly = TRUE) &&
    requireNamespace("yardstick", quietly = TRUE)) {
  spec <- count_reg(family = "NB2") |>
    parsnip::set_engine("flexCountReg", method = "BFGS")
  set.seed(1)
  folds <- rsample::group_vfold_cv(washington_roads, group = ID, v = 3)
  fit_resamples_counts(spec, folds, Total_crashes ~ lnaadt + lnlength,
                       log_score = TRUE)
}
#> # A tibble: 3 × 5
#>   outcome       .metric     mean     n std_err
#>   <chr>         <chr>      <dbl> <int>   <dbl>
#> 1 Total_crashes rmse       0.823     3  0.0581
#> 2 Total_crashes mae        0.489     3  0.0171
#> 3 Total_crashes log_score -0.742     3  0.0143
# }
```
