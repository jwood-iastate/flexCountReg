# Calculate Root Mean Squared Error (RMSE)

This function calculates the Root Mean Squared Error (RMSE) between
observed and predicted values.

## Usage

``` r
rmse(y, mu)
```

## Arguments

- y:

  Numeric vector representing the observed values.

- mu:

  Numeric vector representing the predicted values.

## Value

Numeric value representing the RMSE.

## Details

The RMSE is calculated using the formula: \$\$RMSE = \sqrt{\frac{1}{n}
\sum\_{i=1}^{n} (y_i - \mu_i)^2}\$\$ Where \\y\\ is the vector of
observed values and \\\mu\\ is the vector of predicted values.

## Examples

``` r
y <- c(1, 2, 3)
mu <- c(1.1, 1.9, 3.2)
rmse(y, mu)
#> [1] 0.1414214
```
