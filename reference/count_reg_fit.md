# Fit and predict wrappers used by the parsnip engines

Internal functions called by parsnip; use
[`count_reg`](https://jwood-iastate.github.io/flexCountReg/reference/count_reg.md)
or the estimators instead.

## Usage

``` r
count_reg_fit(
  formula,
  data,
  family = NULL,
  weights = NULL,
  ...,
  model = "countreg",
  predict_args = NULL,
  fit_engine = NULL
)

count_reg_predict(
  object,
  new_data,
  ...,
  format = c("numeric", "interval", "raw")
)
```

## Arguments

- formula, data, family, weights, ...:

  passed to the estimator. `family = NULL` uses the estimator's default.

- model:

  the name of the estimator.

- predict_args:

  a list of arguments for the model's
  [`predict()`](https://rdrr.io/r/stats/predict.html) method.

- fit_engine:

  the estimator's `engine` argument.

- object:

  a fitted model.

- new_data:

  data to predict.

- format:

  `"numeric"`, `"interval"` or `"raw"`.

## Value

`count_reg_fit()`: the estimator's model. `count_reg_predict()`: its
predictions.
