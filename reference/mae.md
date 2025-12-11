# Calculate Mean Absolute Error (MAE)

This function calculates the Mean Absolute Error (MAE) between observed
and predicted values.

## Usage

``` r
mae(y, mu)
```

## Arguments

- y:

  Numeric vector representing the observed values.

- mu:

  Numeric vector representing the predicted values.

## Value

Numeric value representing the MAE.

## Details

The MAE is calculated using the formula: \$\$MAE = \frac{1}{n}
\sum\_{i=1}^{n} \|y_i - \mu_i\|\$\$ Where \\y\\ is the vector of
observed values and \\\mu\\ is the vector of predicted values.

## Examples

``` r
y <- c(1, 2, 3)
mu <- c(1.1, 1.9, 3.2)
mae(y, mu)
#> [1] 0.1333333
```
