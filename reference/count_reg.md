# Count regression model specification for tidymodels

`count_reg()` defines a count regression model for the parsnip package,
fitted with the `"flexCountReg"` engine
([`countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)).
The count distribution is a main argument, so it can be tuned with tune.

## Usage

``` r
count_reg(mode = "regression", family = NULL, engine = "flexCountReg")
```

## Arguments

- mode:

  the model mode; only `"regression"` is available.

- family:

  the count distribution, one of the `family` values of
  [`countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  (for example `"NB2"`, `"PIG2"` or `"PLN"`). Default `"NB2"`.

- engine:

  the fitting engine; only `"flexCountReg"`.

## Value

A `count_reg` model specification.

## Details

Other
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
arguments are engine arguments, set with
`parsnip::set_engine("flexCountReg", ...)`: for example `method`,
`dis_param_formula_1`, `underreport_formula`, `offset` or `stderr`.
Variables used only in those auxiliary formulas, and offset columns,
must be kept in the data by any recipe (in a workflow, add them with
`workflows::add_variables()` or give them a non-predictor role). A
formula [`offset()`](https://rdrr.io/r/stats/offset.html) term can be
used with
`workflows::add_model(spec, formula = y ~ x + offset(log_exposure))`.
Case weights from tidymodels are passed to
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
as frequency weights.

Predictions (`type = "numeric"`) are the expected counts from
[`predict.flexCountReg`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md)
and have one row per row of `new_data`, with `NA` for rows with missing
values.

Requires the parsnip package. The engine is registered when flexCountReg
or parsnip is loaded.

## Examples

``` r
# \donttest{
if (requireNamespace("parsnip", quietly = TRUE)) {
  spec <- count_reg(family = "NB2") |>
    parsnip::set_engine("flexCountReg", method = "BFGS")
  fit <- parsnip::fit(spec, Total_crashes ~ lnaadt + lnlength,
                      data = washington_roads)
  predict(fit, new_data = washington_roads[1:5, ])
}
#> # A tibble: 5 × 1
#>   .pred
#>   <dbl>
#> 1 1.18 
#> 2 1.07 
#> 3 1.56 
#> 4 0.511
#> 5 1.22 
# }
```
