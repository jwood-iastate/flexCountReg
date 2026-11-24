# Time-varying effects in model formulas

`tv()` lets the effect of a variable change over time (temporal
instability). Used in a model formula, `tv(x, time = year)` replaces the
single coefficient of `x` by one effect for each value of `year`,
parameterized by a basis chosen for the question being asked. It can be
used in the mean formula and in the formulas of distribution parameters,
zero models and the other auxiliary formulas of every estimator, and
`tv(1, time = year)` gives a time-varying intercept.

## Usage

``` r
tv(
  x,
  time,
  basis = c("step", "level", "hinge", "smooth"),
  levels = NULL,
  new_time = c("error", "nearest"),
  ...
)
```

## Arguments

- x:

  the numeric variable whose effect varies, or `1` for the intercept.

- time:

  the time variable (numeric, such as a year, or a factor). Its distinct
  values in the estimation data define the periods, in increasing order
  (factor level order for a factor).

- basis:

  how the period effects \\\beta_t\\ are parameterized:

  `"step"` (default)

  :   the effect in the first period and the change at each later
      period: \\\beta_t = \beta_1 + \sum\_{k=2}^t \Delta_k\\. A zero
      change \\\Delta_k\\ means equal effects in periods \\k-1\\ and
      \\k\\. Columns `x:time<first>` and `x:time_step<k>`.

  `"level"`

  :   a separate effect for each period (columns `x:time<t>`). A zero
      coefficient means no effect in that period.

  `"hinge"`

  :   the cumulative "spline indicators" of Shabab et al. (2026):
      \\\beta_t = \sum_k \max(t - k + 1, 0)\\\delta_k\\, with \\t\\ the
      period's position (1, 2, ...). A zero \\\delta_k\\ means that the
      trend of the earlier periods continues at period \\k\\, not that
      the effect stays the same: dropping \\\delta_2\\ imposes \\\beta_2
      = 2\beta_1\\. Columns `x:time_hinge<k>`.

  `"smooth"`

  :   a penalized spline in time, \\\beta(t)\\, the varying-coefficient
      smooth `s(time, by = x)` of mgcv, so the data choose between a
      constant, a linear trend and a flexible path. `time` must be
      numeric. Arguments in `...` (such as `k` and `bs`) are passed to
      [`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html); `k` defaults
      to the number of periods, at most 10.

  The step, level and hinge bases give the same fit (the same
  log-likelihood) and the same period effects; they differ in what the
  individual coefficients and their tests mean.

- levels:

  the periods. Set from the estimation data; give it only to fix the
  periods in advance.

- new_time:

  what to do at prediction for a time value that was not in the
  estimation data: `"error"` (default), or `"nearest"` to use the effect
  of the nearest period (numeric time only), for example the last
  period's effect for a later year.

- ...:

  arguments passed to [`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html)
  with `basis = "smooth"`.

## Value

A numeric matrix with one column per basis function (it is meant to be
used inside a model formula).

## Details

With an intercept in the model, `tv(1, time = year)` leaves out its
first column: the intercept is the first period's level, and the
remaining columns are the changes (`"step"`), the differences from the
first period (`"level"`) or the hinge terms from the second period
(`"hinge"`). A variable should appear either alone or in `tv()`, not
both.

Use
[`temporal_effects()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_effects.md)
for the effect in each period with its standard error, and
[`temporal_tests()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_tests.md)
to test whether and when the effects change.

## References

Shabab, K. R., Bhowmik, T., Zaki, M. H. and Eluru, N. (2026). A joint
count and generalized ordered fractional split approach for addressing
temporal instability in road safety analysis. *Analytic Methods in
Accident Research*, 51, 100440.

## See also

[`temporal_effects()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_effects.md),
[`temporal_tests()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_tests.md)

## Examples

``` r
data("washington_roads")
m <- countreg(Total_crashes ~ tv(lnaadt, time = Year) + speed50 +
                offset(lnlength), data = washington_roads,
              family = "NB2", method = "BFGS")
temporal_effects(m)
#>   part   term time basis period estimate  std.error conf.low conf.high
#> 1 mean lnaadt Year  step   2016 1.130270 0.05200214 1.028348  1.232193
#> 2 mean lnaadt Year  step   2017 1.123700 0.05201507 1.021752  1.225648
#> 3 mean lnaadt Year  step   2018 1.120515 0.05185510 1.018881  1.222149
temporal_tests(m)
#>   part   term         test    estimate  std.error  statistic df   p.value
#> 1 mean lnaadt    all equal          NA         NA  0.5951265  2 0.7426256
#> 2 mean lnaadt 2017 vs 2016 -0.00657049 0.01297047 -0.5065731  1 0.6124544
#> 3 mean lnaadt 2018 vs 2017 -0.00318521 0.01308643 -0.2433979  1 0.8076972
```
