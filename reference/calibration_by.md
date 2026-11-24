# Calibration of a count model by a covariate

Observed and predicted mean counts, probabilities of zero and variances
within groups (quantiles) of a covariate. A smooth can fit the mean well
and still miss the zeros or the dispersion; this table shows where.

## Usage

``` r
calibration_by(object, var, groups = 10)
```

## Arguments

- object:

  a fitted model.

- var:

  the name of a numeric covariate (a column of the data).

- groups:

  number of quantile groups.

## Value

A data frame with one row per group: the covariate's range, `n`,
observed and predicted mean, observed and predicted \\P(Y = 0)\\ and
observed and predicted variance. Predicted probabilities and variances
are averaged over the rows of each group (the predicted variance
includes the spread of the means within the group) and are `NA` for
models without probability predictions.
