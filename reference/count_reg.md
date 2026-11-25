# Count regression model specification for tidymodels

`count_reg()` defines a count regression model for the parsnip package.
Every estimator of flexCountReg is an engine, so all of them can be used
in workflows, resampling, tuning and workflow sets. The count
distribution (`family`) is a main argument.

## Usage

``` r
count_reg(mode = "regression", family = NULL, engine = "flexCountReg")

# S3 method for class 'count_reg'
update(object, parameters = NULL, family = NULL, fresh = FALSE, ...)
```

## Arguments

- mode:

  the model mode; only `"regression"` is available.

- family:

  the count distribution, one of the `family` values of
  [`countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  (for example `"NB2"`, `"PIG2"` or `"PLN"`) or a family from
  [`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md).
  `NULL` uses the default of the engine's estimator (`"NB2"` for
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)).
  The engines `"renb"`, `"repois"` and `"poislind_re"` have a fixed
  distribution and do not accept it; `"negmultinom"` and `"countreg_fe"`
  accept the families of their estimators; `"countreg_mv"` accepts one
  family per outcome.

- engine:

  the fitting engine; see the section *Engines*.

- object:

  a `count_reg()` specification.

- parameters:

  a one-row tibble or named list of main arguments to update.

- fresh:

  logical; replace the arguments (`TRUE`) or modify them?

- ...:

  engine arguments to update.

## Value

A `count_reg` model specification.

## Engines

|  |  |  |
|----|----|----|
| **Engine** | **Estimator** | **Model** |
| `"flexCountReg"` | [`countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md) | any family, with zero-inflated, hurdle and underreporting components |
| `"countreg_rp"` | [`countreg_rp`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md) | random parameters |
| `"countreg_mix"` | [`countreg_mix`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md) | finite mixtures |
| `"countreg_lc"` | [`countreg_lc`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md) | latent classes |
| `"countreg_fe"` | [`countreg_fe`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_fe.md) | fixed effects |
| `"countreg_hier"` | [`countreg_hier`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md) | hierarchical models |
| `"countreg_alloc"` | [`countreg_alloc`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_alloc.md) | total plus allocation (several outcome columns) |
| `"countreg_mv"` | [`countreg_mv`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mv.md) | multivariate (several outcomes) |
| `"renb"` | [`renb`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md) | random-effects negative binomial |
| `"repois"` | [`repois`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md) | random-effects Poisson |
| `"poislind_re"` | [`poislind_re`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md) | random-effects Poisson-Lindley |
| `"negmultinom"` | [`negmultinom`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md) | negative multinomial |

All other arguments of the estimator are engine arguments, set with
[`parsnip::set_engine()`](https://parsnip.tidymodels.org/reference/set_engine.html),
for example
`set_engine("countreg_rp", rpar_formula = ~ -1 + x, ndraws = 200)`.
Arguments with a finite set of choices have tuning parameters (see
[`count_family`](https://jwood-iastate.github.io/flexCountReg/reference/count_family.md)).

The estimators' own `engine` argument (the numerical method) is called
`fit_engine` here, for example
`set_engine("countreg_hier", fit_engine = "laplace")`.

The engine argument `predict_args` (a list) holds arguments passed to
the model's [`predict()`](https://rdrr.io/r/stats/predict.html) method
for numeric predictions and intervals, for example
`predict_args = list(method = "Exact")` for `"countreg_rp"`,
`list(newgroup = "mean_fe")` for `"countreg_fe"` or
`list(marginal = TRUE)` for `"countreg_hier"`.

## Formulas, recipes and extra columns

The engines use the formula interface. In a workflow, the model receives
only the columns the preprocessor declares as outcomes and predictors
(not those with other roles), so every column the model needs must be a
predictor: offsets, grouping variables, and variables of auxiliary
formulas (`dis_param_formula_1`, `rpar_formula`, ...). Two ways to say
which of them enter the mean model:

- Give the model formula, as in
  `workflows::add_model(spec, formula = y ~ x + offset(log_exp))`, with
  [`workflows::add_variables()`](https://workflows.tidymodels.org/reference/add_variables.html)
  or a recipe as the preprocessor. This works for every model, including
  smooths (`s(x)`), formula offsets, random effects (`(1 | site)`),
  [`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
  terms and several outcomes (`cbind(y1, y2) ~ x`).

- Without a model formula, all predictors enter the mean model
  (`y ~ .`), except that columns named in `offset`, `group_var`,
  `panel_id`, `cluster`, `fixef`, `cre_group` or `time`, and the
  variables of `rpar_formula` for `"countreg_rp"` (which are random
  instead of fixed), are left out of the `.`. The outcome is then named
  `..y` in the fitted model.

With a recipe, do not transform the grouping columns (for example with
`step_dummy()` or `step_normalize()`) or an offset.

## Case weights

Case weights are passed to the estimator's `weights`. By default they
are frequency weights. parsnip passes the engine only their values, so
for importance weights with `"flexCountReg"` or `"countreg_rp"` set the
engine argument `weight_type = "sampling"`; a warning is given for
weights that are not whole numbers without it.

## Predictions

`type = "numeric"` gives the expected counts of the estimator's
[`predict()`](https://rdrr.io/r/stats/predict.html) method, with one row
per row of `new_data` and `NA` for rows with missing values. For models
with several outcomes there is one column per outcome or category
(`.pred_<name>`); `predict_args = list(type = "total")` gives the total
of `"countreg_alloc"` as `.pred`. The `"flexCountReg"` engine also has
`type = "conf_int"` and `type = "pred_int"` (see
[`predict.flexCountReg`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md)).
Everything else is available with `type = "raw"`, whose `opts` are
passed to [`predict()`](https://rdrr.io/r/stats/predict.html): for
example
`predict(fit, new_data, type = "raw", opts = list(type = "prob", max_count = 5))`
for count probabilities.

Predictions that use the observed outcome (Empirical Bayes predictions
of random-parameter models, posterior class probabilities) need it in
`new_data`; a workflow removes the outcome, so use the model itself:
`predict(hardhat::extract_fit_engine(wf), newdata = ...)`. All other
functions of the package
([`margEffTable`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md),
[`cureplot`](https://jwood-iastate.github.io/flexCountReg/reference/cureplot.md),
[`rootogram`](https://jwood-iastate.github.io/flexCountReg/reference/rootogram.md),
[`vcov()`](https://rdrr.io/r/stats/vcov.html), ...) also take the
extracted model.

## Several outcomes

The engines `"countreg_alloc"` and `"countreg_mv"` can be fitted and
predicted in workflows, but the resampling functions of tune compute
metrics for one outcome, so use `rsample` splits directly to validate
them.

## Loading

Requires the parsnip package. The engines are registered when
flexCountReg is loaded; load it before using a saved model or workflow.
Parallel workers of tune load it themselves.

## Examples

``` r
# \donttest{
if (requireNamespace("parsnip", quietly = TRUE)) {
  spec <- count_reg(family = "NB2") |>
    parsnip::set_engine("flexCountReg", method = "BFGS")
  fit <- parsnip::fit(spec, Total_crashes ~ lnaadt + lnlength,
                      data = washington_roads)
  predict(fit, new_data = washington_roads[1:5, ])
  predict(fit, new_data = washington_roads[1:5, ], type = "pred_int")

  # Another estimator
  re <- count_reg() |>
    parsnip::set_engine("repois", group_var = "ID", method = "BFGS")
  fit_re <- parsnip::fit(re, Total_crashes ~ lnaadt + lnlength,
                         data = washington_roads)
  tidy(fit_re)
}
#> # A tibble: 4 × 6
#>   term        component    estimate std.error statistic  p.value
#>   <chr>       <chr>           <dbl>     <dbl>     <dbl>    <dbl>
#> 1 (Intercept) mean           -9.00     0.489     -18.4  1.05e-75
#> 2 lnaadt      mean            1.09     0.0598     18.3  8.64e-75
#> 3 lnlength    mean            0.766    0.0847      9.04 1.51e-19
#> 4 lnalpha     distribution   -0.821    0.201      -4.10 4.19e- 5
# }
```
