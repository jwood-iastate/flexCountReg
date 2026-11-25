# Predictions for finite-mixture and latent-class models

Predictions for finite-mixture and latent-class models

## Usage

``` r
# S3 method for class 'countreg_mix'
predict(
  object,
  newdata = NULL,
  ...,
  type = c("response", "component", "variance", "prob", "link"),
  conditioning = c("prior", "history", "posterior"),
  history = NULL,
  new_data = NULL,
  interval = c("none", "confidence", "prediction"),
  level = 0.95
)
```

## Arguments

- object:

  a model from
  [`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md)
  or
  [`countreg_lc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md).

- newdata, new_data:

  optional data frame (give one or the other). When omitted, the data
  used to fit the model (rows removed for missing values give `NA`).

- ...:

  `counts`, `max_count`, `overflow` and `underflow` for `type = "prob"`;
  `ndraws`, `seed` and `max_count` for intervals.

- type:

  the prediction:

  `"response"`

  :   the expected count of the mixture, \\\sum_k \pi_k m_k\\, where
      \\m_k\\ is the actual expected count of component \\k\\ (for PLN
      and PLL it includes the lognormal factor).

  `"component"`

  :   a matrix of the component expected counts \\m_k\\, one column per
      component.

  `"variance"`

  :   the variance of the mixture, \\\sum_k \pi_k (V_k + m_k^2) - m^2\\,
      including the variation between components (`Inf` when a component
      with positive probability has infinite variance).

  `"prob"`

  :   a matrix of count probabilities \\P(Y = y) = \sum_k \pi_k f_k(y)\\
      for the counts \\y\\ in `counts`, or \\y = 0, \dots,\\ `max_count`
      (default the largest training count). `overflow = TRUE` and
      `underflow = TRUE` add columns with the probabilities above and
      below them.

  `"link"`

  :   a matrix of the component linear predictors \\x'\beta_k\\ plus
      offset. There is no single linear predictor of the mixture. With
      random parameters (`rpar_formula`), it is the linear predictor at
      the class's coefficients of the random variables: the means of
      normal, uniform, triangular and gamma random parameters, but for
      lognormal ones the location parameter \\\mu\\ (`Mean:`), whose
      coefficient is \\e^{\mu + \sigma Z}\\: the column then is not on
      the scale of the coefficient at all. It is not the log of the
      component's expected count either.

- conditioning:

  which class probabilities \\\pi_k\\ are used:

  `"prior"` (default)

  :   membership probabilities from the covariates only. The outcomes in
      `newdata` are never used, so these are predictions for new units
      and the basis of fitted values, residuals and diagnostics.

  `"history"`

  :   for panel models: the prior updated once with the likelihood of
      each panel's earlier outcomes in `history`, \\\pi^H\_{gk} \propto
      \pi\_{gk} \prod\_{i \in H_g} f_k(y_i)\\. Rows of panels without
      history use the prior. The outcomes of `newdata` are not used.

  `"posterior"`

  :   probabilities conditional on the outcomes in `newdata` itself
      (whole panels for panel models). This is a retrospective summary,
      not an out-of-sample prediction: the outcomes being "predicted"
      were used.

  With random parameters within the classes, `"history"` and
  `"posterior"` update only the class probabilities: within each class
  the random parameters keep their population distribution (they are not
  updated to the unit's own coefficients), so the component means
  \\m_k\\ are the same for every conditioning. The class likelihoods of
  the update integrate over the random parameters by simulation. For
  panels of the training data (identified by `group_var`) the
  `"posterior"` update uses each panel's estimation draws, so the
  training data reproduce the fitted posterior probabilities
  ([`component_probs()`](https://jwood-iastate.github.io/flexCountReg/reference/component_probs.md));
  other panels, `history` rows, and the rows of a cross-sectional model
  passed as `newdata` use one common set of `max(rpar_draws, 500)`
  Halton draws, so their probabilities differ from the fitted ones by
  simulation error. Without `newdata`, the training rows use their
  estimation draws.

- history:

  for `conditioning = "history"`: a data frame of earlier rows of the
  panels in `newdata`, with outcomes and the panel identifier. History
  rows with a missing outcome or covariate are dropped with a warning; a
  panel without usable history rows keeps its prior. For panel models
  the membership variables must be constant within each panel (in
  `newdata`, in `history` and between the two), as in fitting, so the
  class probabilities never depend on the row order.

- interval:

  `"none"` (default), `"confidence"` (for the expected count) or
  `"prediction"` (for a new count), for `type = "response"` with
  `conditioning = "prior"`; see the section *Intervals of other
  estimators* of
  [`predict.flexCountReg()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md).
  They use the asymptotic normal distribution of the estimates, which is
  poor when components are weakly separated or a mixing proportion is
  near 0.

- level:

  the level of the interval.

## Value

A numeric vector (`"response"`, `"variance"`) with one value per row of
the data, or a matrix (`"component"`, `"prob"`, `"link"`).

## See also

[`component_probs()`](https://jwood-iastate.github.io/flexCountReg/reference/component_probs.md)
for the class probabilities themselves.

## Examples

``` r
set.seed(1)
d <- data.frame(y = c(rpois(300, 1), rpois(200, 6)))
m <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 3,
                  seed = 1)
predict(m, newdata = d[1:3, , drop = FALSE])
#> [1] 3.01 3.01 3.01
predict(m, newdata = d[1:3, , drop = FALSE], type = "component")
#>      Component1 Component2
#> [1,]  0.9645638   6.082595
#> [2,]  0.9645638   6.082595
#> [3,]  0.9645638   6.082595
predict(m, newdata = d[1:3, , drop = FALSE], type = "prob",
        max_count = 5, overflow = TRUE)
#>             0         1         2         3          4          5        >5
#> [1,] 0.229734 0.2262613 0.1233191 0.0684353 0.06027532 0.06487828 0.2270966
#> [2,] 0.229734 0.2262613 0.1233191 0.0684353 0.06027532 0.06487828 0.2270966
#> [3,] 0.229734 0.2262613 0.1233191 0.0684353 0.06027532 0.06487828 0.2270966
```
