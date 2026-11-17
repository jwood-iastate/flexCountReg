# Random parameters and panel models

Road segments, intersections or individuals observed repeatedly are not
independent, and the effect of a variable often differs from one unit to
another. flexCountReg handles both with random-parameter models
([`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md),
also available as
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md))
and random-effects models
([`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
[`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md),
[`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md)
and
[`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md)).

``` r

library(flexCountReg)
data("washington_roads")
length(unique(washington_roads$ID))   # road segments
#> [1] 507
nrow(washington_roads)                # segment-years
#> [1] 1501
```

## Random-parameter models

In a random-parameter model, a coefficient varies across units: . With
`group_var`, all observations of a segment share one draw of the random
parameter, which accounts for the panel structure. The model is
estimated by maximum simulated likelihood with Halton draws; real
analyses should use more draws (`ndraws`) than this example, and check
that the estimates are stable as the number grows.

``` r

rp <- countreg_rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                  data = washington_roads, family = "NB2",
                  rpardists = c(speed50 = "n"), offset = "lnlength",
                  group_var = "ID", ndraws = 100, method = "BFGS")
summary(rp)
#> Call:
#> countreg_rp(formula = Total_crashes ~ lnaadt, rpar_formula = ~-1 + 
#>     speed50, data = washington_roads, family = "NB2", rpardists = c(speed50 = "n"), 
#>     ndraws = 100, group_var = "ID", offset = "lnlength", method = "BFGS")
#> 
#> Model: countreg_rp | family: NB2 | optimizer: BFGS
#> Iterations: 35 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>                Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)     -8.9052     0.4516  -19.72  < 2e-16 ***
#> lnaadt           1.1251     0.0519   21.68  < 2e-16 ***
#> speed50:Mean    -0.6991     0.1467   -4.77  1.9e-06 ***
#> speed50:St.Dev   0.6413     0.1700    3.77  0.00016 ***
#> ln(alpha)       -1.1395     0.2690   -4.24  2.3e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Offset: lnlength (coefficient fixed at 1)
#> 
#> Log-likelihood: -1087.279 | AIC: 2184.559 | BIC: 2205.701 
#> N: 1501 | panels: 507
```

`speed50:Mean` and `speed50:St.Dev` are the mean and standard deviation
of the speed-limit effect across segments. Distributions other than the
normal are available through `rpardists`: lognormal (`"ln"`), gamma
(`"g"`), triangular (`"t"`) and uniform (`"u"`). For the triangular and
uniform, the `St.Dev` parameter is a half-width: the coefficient lies in
.

### Heterogeneity in means

The mean of a random parameter can depend on covariates, following Wood
and Gayah (2025, Eq. 2): . Each random parameter gets its own
coefficient for each heterogeneity variable:

``` r

rp_het <- countreg_rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                      het_mean_formula = ~ ShouldWidth04,
                      data = washington_roads, offset = "lnlength",
                      group_var = "ID", ndraws = 100, method = "BFGS")
tidy(rp_het)
#> # A tibble: 6 × 6
#>   term                          component estimate std.error statistic   p.value
#>   <chr>                         <chr>        <dbl>     <dbl>     <dbl>     <dbl>
#> 1 (Intercept)                   mean        -8.98     0.463     -19.4  7.74e- 84
#> 2 lnaadt                        mean         1.13     0.0532     21.3  6.10e-101
#> 3 speed50:Mean                  random_m…   -0.859    0.159      -5.42 6.08e-  8
#> 4 speed50:St.Dev                random_s…    0.551    0.170       3.24 1.17e-  3
#> 5 HetMean:speed50:ShouldWidth04 heteroge…    0.700    0.240       2.91 3.60e-  3
#> 6 ln(alpha)                     distribu…   -1.13     0.280      -4.05 5.17e-  5
```

### Predictions

There are three prediction methods for random-parameter models (Wood and
Gayah, 2025):

- `"Exact"` integrates over the random-parameter distributions
  analytically (with their moment generating functions), so it has no
  simulation error.
- `"Simulated"` averages over Halton draws.
- `"Individual"` (Empirical Bayes) uses each segment’s observed crashes
  to estimate its own coefficients. It needs the outcome, so it is for
  data with observed outcomes.

``` r

nd <- washington_roads[washington_roads$ID %in% unique(washington_roads$ID)[1:2], ]
cbind(observed = nd$Total_crashes,
      exact = predict(rp, newdata = nd, method = "Exact"),
      simulated = predict(rp, newdata = nd, method = "Simulated"),
      individual = predict(rp, newdata = nd, method = "Individual"))
#>      observed     exact simulated individual
#> [1,]        0 0.8546141 0.8422646  0.5316823
#> [2,]        2 0.7552404 0.7443269  1.0291883
#> [3,]        0 0.8495740 0.8372973  0.5285467
#> [4,]        0 0.7507863 0.7399372  1.0231186
#> [5,]        1 0.8957947 0.8828501  0.5573020
#> [6,]        3 0.7916325 0.7801931  1.0787808
```

The segment-specific coefficients from the Empirical Bayes method are
returned as an attribute:

``` r

eb <- predict(rp, newdata = washington_roads, method = "Individual")
summary(attr(eb, "eb_coefficients")[, "speed50"])
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -1.3126 -0.7051 -0.7051 -0.7080 -0.7051  0.1780
```

A caution for lognormal random parameters: the expected value of does
not exist when is lognormal and , so `"Exact"` returns `Inf` with a
warning there, and `"Simulated"` predictions depend on the number of
draws. Use the negative of the variable, a gamma random parameter, or
the Empirical Bayes predictions.

## Random-effects models

Random-effects models add a panel-level effect to every observation of a
unit. The negative binomial random-effects model of Hausman, Hall and
Griliches
([`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)),
the Poisson-Lindley random-effects model
([`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md))
and the negative multinomial models
([`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md))
are available:

``` r

re_models <- list(
  RENB = renb(Total_crashes ~ lnaadt + speed50, group_var = "ID",
              data = washington_roads, offset = "lnlength", method = "BFGS"),
  PoisLindRE = poisLind.re(Total_crashes ~ lnaadt + speed50, group_var = "ID",
                           data = washington_roads, offset = "lnlength",
                           method = "BFGS"),
  NM2 = negmultinom(Total_crashes ~ lnaadt + speed50, group_var = "ID",
                    data = washington_roads, offset = "lnlength")
)
do.call(rbind, lapply(names(re_models), function(nm) {
  cbind(model = nm,
        glance(re_models[[nm]])[, c("logLik", "AIC", "BIC", "nobs",
                                    "n_groups")])
}))
#>        model    logLik      AIC      BIC nobs n_groups
#> 1       RENB -1070.458 2150.916 2172.059 1501      507
#> 2 PoisLindRE -1070.830 2149.660 2166.574 1501      507
#> 3        NM2 -1070.454 2148.908 2165.822 1501      507
```

For panel models, [`BIC()`](https://rdrr.io/r/stats/AIC.html) uses the
number of panels (`n_groups`) as the sample size, because the panels are
the independent units.

### Standard errors

The likelihood of a panel model has one contribution per panel, so
robust (sandwich) standard errors are clustered by panel. Bootstrap
standard errors resample whole panels:

``` r

re_robust <- poisLind.re(Total_crashes ~ lnaadt + speed50, group_var = "ID",
                         data = washington_roads, offset = "lnlength",
                         method = "BFGS", stderr = "robust")
tidy(re_robust)
#> # A tibble: 4 × 6
#>   term        component    estimate std.error statistic  p.value
#>   <chr>       <chr>           <dbl>     <dbl>     <dbl>    <dbl>
#> 1 (Intercept) mean           -8.77     0.599     -14.6  1.36e-48
#> 2 lnaadt      mean            1.11     0.0702     15.9  1.41e-56
#> 3 speed50     mean           -0.561    0.131      -4.29 1.79e- 5
#> 4 ln(theta)   distribution   -7.13     0.0642   -111.   0    
```

Cross-sectional models can also be fitted with `stderr = "robust"` or
`"boot"`; for repeated observations of the same units, prefer a panel
model or clustered resampling.

## Cross-validation

When tuning or comparing these models by cross-validation (for example
with the tidymodels
[`count_reg()`](https://jwood-iastate.github.io/flexCountReg/reference/count_reg.md)
specification), resample whole panels, for example with
`rsample::group_vfold_cv(data, group = ID)`. Splitting the observations
of one segment between the training and test sets leaks information and
makes the models look better than they are.

## References

Wood, J.S., & Gayah, V. (2025). Out-of-sample prediction and
interpretation for random parameter generalized linear models. *Accident
Analysis and Prevention*, 220, 108147.

Hausman, J. A., Hall, B. H., & Griliches, Z. (1984). Econometric models
for count data with an application to the patents-R&D relationship.
*Econometrica*, 52(4), 909-938.
