# Calculate Root Mean Squared Error (RMSE)

This function calculates the Root Mean Squared Error (RMSE) between
observed and predicted values.

## Usage

``` r
rmse(y, mu, ...)
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

Numeric value representing the RMSE (with a data frame, the yardstick
tibble).

## Details

The RMSE is calculated using the formula: \$\$RMSE = \sqrt{\frac{1}{n}
\sum\_{i=1}^{n} (y_i - \mu_i)^2}\$\$ Where \\y\\ is the vector of
observed values and \\\mu\\ is the vector of predicted values.

## yardstick

With a data frame, `rmse(data, truth, estimate)` calls the yardstick
metric of the same name, and the function can be used in
[`yardstick::metric_set()`](https://yardstick.tidymodels.org/reference/metric_set.html).
Loading flexCountReg after tidymodels (which masks yardstick's function)
therefore does not change tidymodels code.

## Examples

``` r
y <- c(1, 2, 3)
mu <- c(1.1, 1.9, 3.2)
rmse(y, mu)
#> [1] 0.1414214
```
