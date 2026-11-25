# Component (class) membership probabilities

Prior or posterior probabilities of the mixture components (latent
classes). `class_probs()` is the same function under the latent-class
name.

## Usage

``` r
component_probs(
  object,
  newdata = NULL,
  conditioning = c("prior", "posterior", "history"),
  history = NULL
)

class_probs(object, ...)
```

## Arguments

- object:

  a model from
  [`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md)
  or
  [`countreg_lc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md).

- newdata:

  optional data frame. Without it, the probabilities of the independent
  units of the training data are returned (one row per row, or per panel
  for panel models).

- conditioning:

  `"prior"` (default): membership probabilities from the covariates;
  `"posterior"`: conditional on the outcomes (each unit's own outcomes;
  a retrospective classification, not a prediction); `"history"`: for
  panel models, conditional on earlier outcomes in `history` (see
  [`predict.countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.countreg_mix.md)).

- history:

  earlier rows of the panels, for `conditioning = "history"`.

- ...:

  passed from `class_probs()`.

## Value

A matrix with one column per component (named by the component labels)
and one row per unit (training data) or per row of `newdata`. These are
probabilities of the components, not of counts (for count probabilities
use `predict(type = "prob")`).

## Details

The posterior probabilities give a descriptive classification (e.g. the
most probable class of each unit), which is uncertain: units with
probabilities near \\1/K\\ are not well classified. The entropy in
[`summary()`](https://rdrr.io/r/base/summary.html) summarizes that
uncertainty.

With random parameters within the classes, the posterior of a unit
integrates over its random parameters by simulation. Panels of the
training data passed in `newdata` (matched by `group_var`) use their
estimation draws and so reproduce the stored probabilities; new panels,
and the rows of a cross-sectional model passed as `newdata`, use a
common set of Halton draws (see
[`predict.countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.countreg_mix.md)).

## Examples

``` r
set.seed(1)
d <- data.frame(y = c(rpois(300, 1), rpois(200, 6)))
m <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 3,
                  seed = 1)
head(component_probs(m, conditioning = "posterior"))
#>      Component1  Component2
#> [1,]  0.9960297 0.003970273
#> [2,]  0.9754798 0.024520217
#> [3,]  0.9754798 0.024520217
#> [4,]  0.8631757 0.136824277
#> [5,]  0.9960297 0.003970273
#> [6,]  0.8631757 0.136824277
```
