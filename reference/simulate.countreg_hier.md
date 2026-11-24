# Simulate responses from a hierarchical count model

Simulate responses from a hierarchical count model

## Usage

``` r
# S3 method for class 'countreg_hier'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  newdata = NULL,
  use_re = c("conditional", "new"),
  ...
)
```

## Arguments

- object:

  a
  [`countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md)
  fit.

- nsim:

  number of simulated response vectors.

- seed:

  optional random seed.

- newdata:

  optional data frame (default: the estimation rows).

- use_re:

  `"conditional"` keeps the estimated effects of the groups (conditional
  modes); `"new"` draws new effects for every group from their estimated
  distribution.

- ...:

  unused.

## Value

A data frame with `nsim` columns.
