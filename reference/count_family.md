# Tuning parameters for count regression models

dials parameters for
[`count_reg`](https://jwood-iastate.github.io/flexCountReg/reference/count_reg.md):
the count distribution (the main argument `family`) and engine
arguments. Arguments of
[`count_reg()`](https://jwood-iastate.github.io/flexCountReg/reference/count_reg.md)
marked with
[`tune()`](https://hardhat.tidymodels.org/reference/tune.html) use them
automatically.

## Usage

``` r
count_family(
  values = c("POISSON", "NB1", "NB2", "NBP", "GP1", "GP2", "PIG1", "PIG2", "PL", "PLN",
    "COM")
)

count_components(range = c(2L, 4L), trans = NULL)

count_zero_type(values = c("inflated", "hurdle"))

count_correlated(values = c(FALSE, TRUE))

count_re_dist(values = c("normal", "uniform", "triangular", "gamma"))
```

## Arguments

- values:

  the candidate values.

- range:

  the smallest and largest number of components or classes.

- trans:

  a transformation (see
  [`dials::new_quant_param()`](https://dials.tidymodels.org/reference/new-param.html)).

## Value

A dials parameter.

## Details

|  |  |  |
|----|----|----|
| **Parameter** | **Argument** | **Engine** |
| `count_family()` | `family` | all with a family |
| `count_components()` | `components`, `classes` | `"countreg_mix"`, `"countreg_lc"` |
| `count_zero_type()` | `zero_type` | `"flexCountReg"` |
| `count_correlated()` | `correlated` | `"countreg_rp"`, `"countreg_mix"`, `"countreg_lc"` |
| `count_re_dist()` | `re_dist` | `"countreg_hier"` (one term) |

The families supported by `"countreg_fe"` (`"POISSON"` and `"NB2"` by
default here) and `"negmultinom"` (`"NM1"`, `"NM2"` and `"NMP"`) are
their own default values. Any other engine argument can be tuned with an
explicit grid, or with [`update()`](https://rdrr.io/r/stats/update.html)
on the parameter set.

## Examples

``` r
if (requireNamespace("dials", quietly = TRUE)) {
  count_family()
  count_components()
}
#> # Mixture Components (quantitative)
#> Range: [2, 4]
```
