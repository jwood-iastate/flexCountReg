# Update and refit a flexCountReg model

As [`update`](https://rdrr.io/r/stats/update.html), also for models
fitted with tidymodels (see
[`count_reg`](https://jwood-iastate.github.io/flexCountReg/reference/count_reg.md)),
which are refitted with their stored training data.

## Usage

``` r
# S3 method for class 'flexCountReg'
update(object, formula., ..., evaluate = TRUE)

# S3 method for class 'countreg_alloc'
update(object, formula., ..., evaluate = TRUE)

# S3 method for class 'countreg_mv'
update(object, formula., ..., evaluate = TRUE)
```

## Arguments

- object:

  a fitted model.

- formula.:

  changes to the formula (see
  [`update.formula`](https://rdrr.io/r/stats/update.formula.html)).

- ...:

  arguments to change or add.

- evaluate:

  if `FALSE`, return the call instead of the model.

## Value

The refitted model, or the call.
