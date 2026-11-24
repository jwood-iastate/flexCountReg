# Concurvity of smooth terms

Measures how much each smooth can be reproduced by the other terms of
its linear predictor: the nonlinear analogue of collinearity. High
concurvity makes a smooth's shape poorly determined and its standard
errors large.

## Usage

``` r
concurvity_fcr(object)
```

## Arguments

- object:

  a fitted model with smooth terms.

## Value

A data frame with `parameter`, `term`, `observed` and `worst`.

## Details

For a smooth \\f_j = X_j b_j\\ at the data rows, `observed` is the share
of \\\\f_j\\^2\\ explained by the other columns of the same linear
predictor (the parametric terms and the other smooths), and `worst` is
the largest squared canonical correlation between the smooth's basis and
those columns, the worst case over all possible estimates. Both are
between 0 and 1; values near 1 indicate a smooth that the rest of the
model can mimic. This follows
[`mgcv::concurvity()`](https://rdrr.io/pkg/mgcv/man/concurvity.html) in
spirit.
