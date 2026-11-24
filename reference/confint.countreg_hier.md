# Confidence intervals for a hierarchical count model

Wald intervals for every parameter, or profile-likelihood intervals for
the fixed coefficients, the distribution parameters and the standard
deviation of the first effect of each random-effect term.

## Usage

``` r
# S3 method for class 'countreg_hier'
confint(object, parm, level = 0.95, method = c("wald", "profile"), ...)
```

## Arguments

- object:

  a
  [`countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md)
  fit.

- parm:

  parameter names (by default all for `"wald"`, and the random-effect
  standard deviations for `"profile"`).

- level:

  confidence level.

- method:

  `"wald"` (estimate plus or minus z standard errors) or `"profile"`.

- ...:

  unused.

## Value

A matrix with the lower and upper limits.

## Details

A profile-likelihood interval contains the values whose profile
log-likelihood (maximized over the other parameters) is within
\\\chi^2_1(level)/2\\ of the maximum. It is not symmetric, and for a
standard deviation its lower limit can be 0, where the Wald interval is
unreliable. Each profile point refits the model, so profiling takes
several fits per limit. Correlations, and standard deviations of second
and later effects of a term, are not profiled (their Cholesky parameters
mix several variances).
