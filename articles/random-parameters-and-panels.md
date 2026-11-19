# Random parameters and panel models

Road segments, intersections or individuals observed repeatedly are not
independent, and the effect of a variable often differs from one unit to
another. flexCountReg handles both with random-parameter models
([`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md))
and random-effects models
([`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
[`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md),
[`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md)
and
[`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md)).

``` r

library(flexCountReg)
```

## Why the panel structure matters: repeated seizure counts

[`MASS::epil`](https://rdrr.io/pkg/MASS/man/epil.html) records epileptic
seizure counts for 59 patients in a clinical trial of progabide (Thall
and Vail, 1990). Each patient was observed in four consecutive two-week
periods, so there are four rows per patient. The independent unit is the
patient, not the row.

``` r

epil <- MASS::epil
c(rows = nrow(epil), patients = length(unique(epil$subject)))
#>     rows patients 
#>      236       59
table(table(epil$subject))      # rows per patient
#> 
#>  4 
#> 59
# treatment does not change within a patient
all(tapply(epil$trt, epil$subject, function(v) length(unique(v))) == 1)
#> [1] TRUE
```

The response `y` is the count in a two-week follow-up period. The count
in the eight-week *baseline* period before treatment (`base`) is a
predictor, not a response: MASS provides it as `lbase`, the centred
logarithm of the baseline count divided by 4 (a two-week equivalent),
and the age as `lage`, the centred log age. `V4` marks the fourth visit.
(The data in MASS 7.3-65 and later include a documented correction to
one count.)

A pooled model treats the 236 rows as independent; a random-effects
model
([`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
the negative binomial model of Hausman, Hall and Griliches, 1984) gives
each patient a shared effect:

``` r

f_epil <- y ~ trt + lbase + lage + V4
pooled <- countreg(f_epil, data = epil, family = "NB2", method = "BFGS")
pooled_robust <- countreg(f_epil, data = epil, family = "NB2",
                          method = "BFGS", stderr = "robust")
panel <- renb(f_epil, group_var = "subject", data = epil, method = "BFGS")
se_trt <- function(m) sqrt(diag(vcov(m)))[["trtprogabide"]]
data.frame(
  model = c("pooled NB2", "pooled NB2, robust SE", "RENB (panel)"),
  treatment_effect = c(coef(pooled)[["trtprogabide"]],
                       coef(pooled_robust)[["trtprogabide"]],
                       coef(panel)[["trtprogabide"]]),
  std_error = c(se_trt(pooled), se_trt(pooled_robust), se_trt(panel)))
#>                   model treatment_effect std_error
#> 1            pooled NB2       -0.2342385 0.1008920
#> 2 pooled NB2, robust SE       -0.2342385 0.1059424
#> 3          RENB (panel)       -0.2324608 0.1431267
```

The pooled model’s standard error for the treatment effect is 0.101,
against 0.143 in the panel model. Rows from the same patient share that
patient’s seizure tendency, so the pooled model counts four rows as four
independent pieces of evidence and overstates its precision. Robust
standard errors computed row by row (0.106) do not fix this: they do not
group the rows by patient.

``` r

do.call(rbind, lapply(list(pooled = pooled, RENB = panel), function(m)
  glance(m)[, c("logLik", "AIC", "BIC", "nobs", "n_groups")]))
#> # A tibble: 2 × 5
#>   logLik   AIC   BIC  nobs n_groups
#> *  <dbl> <dbl> <dbl> <int>    <int>
#> 1  -651. 1314. 1335.   236       NA
#> 2  -632. 1277. 1292.   236       59
```

For panel models, [`BIC()`](https://rdrr.io/r/stats/AIC.html) uses the
number of panels (`n_groups`, 59 patients) as the sample size, because
the panels are the independent units. The same unit applies to every
form of inference: the likelihood has one contribution per panel, robust
standard errors are clustered by panel, bootstrap standard errors
resample whole panels, and validation splits should keep each panel in
one fold.

The rest of this vignette uses the `washington_roads` data:

``` r

data("washington_roads")
length(unique(washington_roads$ID))   # road segments
#> [1] 507
nrow(washington_roads)                # segment-years
#> [1] 1501
```

## Random-parameter models

In a random-parameter model, a coefficient varies across units:
$`\beta_i = \mu + \sigma \upsilon_i`$. With `group_var`, all
observations of a segment share one draw of the random parameter, which
accounts for the panel structure. The model is estimated by maximum
simulated likelihood with Halton draws; real analyses should use more
draws (`ndraws`) than this example, and check that the estimates are
stable as the number grows.

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
normal are available through `rpardists`. The two estimated parameters,
labelled `Mean` and `St.Dev`, do not always equal the actual mean and
standard deviation of the coefficient:

| `rpardists` | Coefficient | Actual mean | Actual SD |
|----|----|----|----|
| `"n"` normal | $`\beta = \mu + \sigma z`$ | $`\mu`$ | $`\sigma`$ |
| `"ln"` lognormal | $`\beta = \exp(\mu + \sigma z)`$ | $`\exp(\mu + \sigma^2/2)`$ | $`\exp(\mu + \sigma^2/2)\sqrt{\exp(\sigma^2) - 1}`$ |
| `"g"` gamma | gamma with mean $`\mu`$ and SD $`\sigma`$ | $`\mu`$ | $`\sigma`$ |
| `"t"` triangular | $`\mu + \sigma T`$, $`T`$ triangular on $`[-1, 1]`$ | $`\mu`$ | $`\sigma/\sqrt{6}`$ |
| `"u"` uniform | $`\mu + \sigma U`$, $`U`$ uniform on $`[-1, 1]`$ | $`\mu`$ | $`\sigma/\sqrt{3}`$ |

For the lognormal, $`\mu`$ and $`\sigma`$ are the mean and SD of
$`\log\beta`$, and the coefficient is positive; for a negative effect,
enter the negative of the variable. For the triangular and uniform,
`St.Dev` is a half-width: the coefficient lies in
$`[\mu - \sigma, \mu + \sigma]`$.

### Heterogeneity in means

The mean of a random parameter can depend on covariates, following Wood
and Gayah (2025, Eq. 2):
$`\beta_i = \mu + z_i'\delta + \sigma\upsilon_i`$. Each random parameter
gets its own coefficient for each heterogeneity variable:

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

A caution about the expected count: with a random parameter the expected
value of $`\exp(x\beta)`$ must exist. It is infinite when $`\beta`$ is
lognormal and $`x > 0`$, and when $`\beta`$ is gamma and
$`x\sigma^2/\mu \ge 1`$. In those cases `"Exact"` and `"Simulated"`
predictions return `Inf` with a warning (a simulated average would
otherwise depend on the number of draws without converging). This is a
property of the model, not a numerical failure: changing the
distribution only helps if the new distribution’s expectation exists for
the covariate values used. Normal, triangular and uniform parameters
always have a finite expectation.

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
#> Warning: The Hessian is not negative definite at the estimates, so the standard
#> errors of ln(theta) are NA. The estimate may be on the boundary of the
#> parameter space or not at a maximum.
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

The warning above comes from the Poisson-Lindley random-effects model:
its random-effect parameter is driven toward the edge of its parameter
space, where the likelihood is flat, so its standard error is not
identified and is reported as `NA`:

``` r

tidy(re_models$PoisLindRE)
#> # A tibble: 4 × 6
#>   term        component    estimate std.error statistic   p.value
#>   <chr>       <chr>           <dbl>     <dbl>     <dbl>     <dbl>
#> 1 (Intercept) mean           -8.77     0.512     -17.1   1.22e-65
#> 2 lnaadt      mean            1.11     0.0603     18.5   4.18e-76
#> 3 speed50     mean           -0.561    0.132      -4.24  2.19e- 5
#> 4 ln(theta)   distribution   -7.13    NA          NA    NA    
```

Its fit statistics are close to those of the other two models, but with
one parameter not identified its reported uncertainty is incomplete.
When models fit about equally well, prefer one whose parameters are all
identified, such as RENB or the negative multinomial model here.

### Standard errors

The likelihood of a panel model has one contribution per panel, so
robust (sandwich) standard errors are clustered by panel. Bootstrap
standard errors resample whole panels:

``` r

re_robust <- poisLind.re(Total_crashes ~ lnaadt + speed50, group_var = "ID",
                         data = washington_roads, offset = "lnlength",
                         method = "BFGS", stderr = "robust")
#> Warning: The Hessian is not negative definite at the estimates, so the standard
#> errors of ln(theta) are NA. The estimate may be on the boundary of the
#> parameter space or not at a maximum.
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

## Next steps

- [`vignette("prediction-for-panel-data", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/prediction-for-panel-data.md)
  separates predictions for new patients from those for known patients.
- [`vignette("predictive-validation", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/predictive-validation.md)
  builds grouped cross-validation step by step.
- [`vignette("numerical-reliability", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/numerical-reliability.md)
  shows how to check the number of draws and other numerical settings.

## References

Thall, P. F. and Vail, S. C. (1990). Some covariance models for
longitudinal count data with overdispersion. *Biometrics*, 46(3),
657-671.

Wood, J.S., & Gayah, V. (2025). Out-of-sample prediction and
interpretation for random parameter generalized linear models. *Accident
Analysis and Prevention*, 220, 108147.

Hausman, J. A., Hall, B. H., & Griliches, Z. (1984). Econometric models
for count data with an application to the patents-R&D relationship.
*Econometrica*, 52(4), 909-938.
