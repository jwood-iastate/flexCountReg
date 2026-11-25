# Calculate Mean Absolute Error (MAE)

This function calculates the Mean Absolute Error (MAE) between observed
and predicted values.

## Usage

``` r
mae(y, mu, ...)
```

## Arguments

- y:

  Numeric vector representing the observed values, or a data frame (see
  the section on yardstick).

- mu:

  Numeric vector representing the predicted values (with a data frame,
  the column of observed values).

- ...:

  with a data frame, further arguments of the yardstick metric (the
  column of predictions first).

## Value

Numeric value representing the MAE (with a data frame, the yardstick
tibble).

## Details

The MAE is calculated using the formula: \$\$MAE = \frac{1}{n}
\sum\_{i=1}^{n} \|y_i - \mu_i\|\$\$ Where \\y\\ is the vector of
observed values and \\\mu\\ is the vector of predicted values.

## yardstick

With a data frame, `mae(data, truth, estimate)` calls the yardstick
metric of the same name, and the function can be used in
[`yardstick::metric_set()`](https://yardstick.tidymodels.org/reference/metric_set.html).

## Examples

``` r
y <- c(1, 2, 3)
mu <- c(1.1, 1.9, 3.2)
mae(y, mu)
#> [1] 0.1333333
```
