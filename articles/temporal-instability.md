# Time-varying effects: temporal instability

With several years of data, the effect of a variable may change over
time: vehicles, enforcement and driver behavior change while the roads
stay the same. This is **temporal instability**. `tv(x, time = year)` in
a model formula gives `x` one effect per year, and
[`temporal_effects()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_effects.md)
and
[`temporal_tests()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_tests.md)
report and test those effects.

``` r

library(flexCountReg)
data("washington_roads")
w <- washington_roads
table(w$Year)
#> 
#> 2016 2017 2018 
#>  501  500  500
```

## 1. One effect per year

``` r

m <- countreg(Total_crashes ~ tv(lnaadt, time = Year) + speed50 +
                offset(lnlength), data = w, family = "NB2",
              method = "BFGS")
coef(m)
#>          (Intercept)      lnaadt:Year2016 lnaadt:Year_step2017 
#>          -8.89992454           1.13027036          -0.00657049 
#> lnaadt:Year_step2018              speed50            ln(alpha) 
#>          -0.00318521          -0.56770930          -0.92122379
```

The default **step** basis has two kinds of coefficients:

- `lnaadt:Year2016` is the effect in the first year;
- `lnaadt:Year_step2017` and `lnaadt:Year_step2018` are the changes from
  the year before.

A change near zero means the effect is the same in two adjacent years.
[`temporal_effects()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_effects.md)
gives the effect in each year, with its standard error and confidence
interval:

``` r

temporal_effects(m)
#>   part   term time basis period estimate  std.error conf.low conf.high
#> 1 mean lnaadt Year  step   2016 1.130270 0.05200214 1.028348  1.232193
#> 2 mean lnaadt Year  step   2017 1.123700 0.05201507 1.021752  1.225648
#> 3 mean lnaadt Year  step   2018 1.120515 0.05185510 1.018881  1.222149
```

[`temporal_tests()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_tests.md)
tests whether the effects are all equal, and whether each adjacent pair
of years differs:

``` r

temporal_tests(m)
#>   part   term         test    estimate  std.error  statistic df   p.value
#> 1 mean lnaadt    all equal          NA         NA  0.5951265  2 0.7426256
#> 2 mean lnaadt 2017 vs 2016 -0.00657049 0.01297047 -0.5065731  1 0.6124544
#> 3 mean lnaadt 2018 vs 2017 -0.00318521 0.01308643 -0.2433979  1 0.8076972
```

Here the effect of traffic is stable over the three years.

## 2. Choosing the basis

The basis decides what each coefficient means. The fit and the year
effects are the same for every basis.

| `basis` | effect in year \\t\\ | a zero coefficient means |
|----|----|----|
| `"step"` (default) | \\\beta_1 + \sum\_{k=2}^t \Delta_k\\ | no change between years \\k-1\\ and \\k\\ |
| `"level"` | \\\beta_t\\ | no effect in year \\t\\ |
| `"hinge"` | \\\sum_k \max(t-k+1, 0)\\\delta_k\\ | the earlier trend continues at year \\k\\ |
| `"smooth"` | \\\beta(t)\\, a penalized spline | (the smoothness is estimated) |

``` r

fit_basis <- function(b) {
  countreg(stats::as.formula(paste0(
    "Total_crashes ~ tv(lnaadt, time = Year, basis = '", b,
    "') + speed50 + offset(lnlength)")),
    data = w, family = "NB2", method = "BFGS")
}
mh <- fit_basis("hinge")
c(step = logLik(m), hinge = logLik(mh))
#>      step     hinge 
#> -1090.263 -1090.263
cbind(step = temporal_effects(m)$estimate,
      hinge = temporal_effects(mh)$estimate)
#>          step    hinge
#> [1,] 1.130270 1.130270
#> [2,] 1.123700 1.123699
#> [3,] 1.120515 1.120514
```

### A caution about the hinge basis

The **hinge** basis is the “spline indicator” model of Shabab et al.
(2026). Its columns are cumulative: the indicator for year \\k\\ is
\\\max(t - k + 1, 0)\\. As a result, \\\delta_2 = \beta_2 - 2\beta_1\\.
It is not \\\beta_2 - \beta_1\\.

Dropping a hinge coefficient that is not significant therefore does not
make two years equal. It continues the trend of the earlier years. With
the paper’s own AADT estimates:

| year | separate effects | dropping \\\delta_2\\ imposes |
|------|------------------|-------------------------------|
| 1    | 0.122            | 0.122                         |
| 2    | 0.215            | \\2 \times 0.122 = 0.244\\    |

Equal effects in years 1 and 2 would be 0.122 in both. The basis also
numbers the years from 1, so the restriction depends on which year comes
first. To find years with equal effects, use the step basis, where a
zero coefficient means exactly that, or
[`temporal_tests()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_tests.md),
which tests the year effects themselves.

### A smooth path over time

With many years, `basis = "smooth"` fits the effect as a penalized
spline in time, the varying-coefficient smooth `s(Year, by = lnaadt)` of
mgcv. The data then choose between a constant effect, a linear trend and
a flexible path:

``` r

msm <- fit_basis("smooth")
temporal_effects(msm)
#>   part   term time  basis period estimate  std.error conf.low conf.high
#> 1 mean lnaadt Year smooth   2016 1.129769 0.05189659 1.028053  1.231484
#> 2 mean lnaadt Year smooth   2017 1.124881 0.05141355 1.024112  1.225650
#> 3 mean lnaadt Year smooth   2018 1.119993 0.05174021 1.018584  1.221402
```

Any mgcv basis can be used (for example `bs = "cr"`); the effects are
reported at the years of the estimation data. Smooths of time written
directly with `s()` are ordinary smooth terms, not time-varying effects,
and
[`temporal_effects()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_effects.md)
leaves them out.

## 3. Time-varying intercepts and distribution parameters

`tv(1, time = Year)` lets the intercept vary by year. With an intercept
in the model, the intercept is the first year’s level and the
[`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
columns are the changes after it.
[`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
also works in the formulas of the distribution parameters, here a
separate NB2 overdispersion for each year:

``` r

mi <- countreg(Total_crashes ~ tv(1, time = Year) + lnaadt +
                 offset(lnlength), data = w, family = "NB2",
               dis_param_formula_1 = ~ tv(1, time = Year, basis = "level"),
               method = "BFGS")
temporal_effects(mi)
#>        part        term time basis period   estimate std.error   conf.low
#> 1      mean (Intercept) Year  step   2016 -9.3636110 0.4574756 -10.260247
#> 2      mean (Intercept) Year  step   2017 -9.4307156 0.4588819 -10.330108
#> 3      mean (Intercept) Year  step   2018 -9.4248184 0.4570932 -10.320705
#> 4 ln(alpha) (Intercept) Year level   2016 -0.8703054 0.3746837  -1.604672
#> 5 ln(alpha) (Intercept) Year level   2017 -1.1701419 0.4819031  -2.114655
#> 6 ln(alpha) (Intercept) Year level   2018 -0.4144498 0.3105474  -1.023112
#>    conf.high
#> 1 -8.4669754
#> 2 -8.5313236
#> 3 -8.5289323
#> 4 -0.1359388
#> 5 -0.2256291
#> 6  0.1942120
```

## 4. Prediction for new years

[`predict()`](https://rdrr.io/r/stats/predict.html) uses the years of
the estimation data. A year that was not estimated is an error by
default, because no effect was estimated for it.
`tv(x, time = Year, new_time = "nearest")` uses the nearest year’s
effect instead, for example the last year’s effect when forecasting:

``` r

mn <- countreg(Total_crashes ~ tv(lnaadt, time = Year, new_time = "nearest") +
                 speed50 + offset(lnlength), data = w, family = "NB2",
               method = "BFGS")
nd <- w[1:3, ]
nd$Year <- 2019
predict(mn, newdata = nd)
#> [1] 0.7657030 0.6766678 1.1218439
```

## 5. Other models

[`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
works in the formulas of every estimator:
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md),
[`countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md),
[`countreg_fe()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_fe.md),
[`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
[`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md),
[`poislind_re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)
and
[`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md).
Here the effect of traffic varies by year in a random-parameters model:

``` r

r <- countreg_rp(Total_crashes ~ tv(lnaadt, time = Year) + offset(lnlength),
                 rpar_formula = ~ -1 + speed50, data = w, family = "NB2",
                 ndraws = 100, method = "BFGS", group_var = "ID")
temporal_effects(r)
#>   part   term time basis period estimate  std.error conf.low conf.high
#> 1 mean lnaadt Year  step   2016 1.130707 0.05300066 1.026827  1.234586
#> 2 mean lnaadt Year  step   2017 1.123370 0.05302431 1.019444  1.227296
#> 3 mean lnaadt Year  step   2018 1.121112 0.05277747 1.017670  1.224554
```

[`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
can also appear in `rpar_formula`.
[`temporal_effects()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_effects.md)
then reports the year effects of the random parameters’ means, and each
basis column has its own standard deviation.

In the multivariate models of
[`countreg_mv()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mv.md),
[`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
terms can differ by outcome, and
[`temporal_effects()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_effects.md)
and
[`temporal_tests()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_tests.md)
add a column `outcome`.

A [`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
term inside an interaction, such as `tv(lnaadt, time = Year):speed50`,
is not summarized by
[`temporal_effects()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_effects.md)
(it warns). Put the product inside
[`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
instead: `tv(lnaadt * speed50, time = Year)`.

## Reference

Shabab, K. R., Bhowmik, T., Zaki, M. H. and Eluru, N. (2026). A joint
count and generalized ordered fractional split approach for addressing
temporal instability in road safety analysis. *Analytic Methods in
Accident Research*, 51, 100440.
