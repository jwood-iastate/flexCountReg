# A distribution parameter of a user-defined count family

A distribution parameter of a user-defined count family

## Usage

``` r
count_param(link = c("log", "identity", "logit"), start = 0)
```

## Arguments

- link:

  the link of the parameter's linear predictor: `"log"` (default; a
  positive parameter), `"identity"` (any real value) or `"logit"` (a
  parameter in (0, 1)).

- start:

  the starting value of the linear predictor (on the link scale).

## Value

An object of class `"fcr_count_param"`.

## See also

[`new_count_family`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md)
