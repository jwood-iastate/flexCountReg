# Inference, prediction and model comparison

This vignette covers what to do with a fitted model: which standard
errors to use, how to predict and simulate, which residuals and
diagnostics are available, and how to compare models, including when a
likelihood-ratio (LR) test is valid and when it is not.

``` r

library(flexCountReg)
data("washington_roads")
d <- washington_roads
```

## Standard errors

Every estimator has a `stderr` argument with three choices. The choice
is made when the model is fit, and
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and the
delta-method standard errors of
[`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
all use it.

- `"normal"` (default): the inverse of the negative Hessian of the
  log-likelihood. It is valid when the model is correctly specified,
  including the count distribution.
- `"robust"`: the sandwich estimator. It remains valid for the mean
  coefficients when the distribution is misspecified but the mean model
  is right (quasi-likelihood). It does not rescue estimates from a wrong
  mean model (they are then estimates of something else, whatever their
  standard errors), and with sampling weights it is not a full
  survey-design analysis (strata, clusters and finite-population
  corrections are not modelled). For panel models
  ([`countreg.rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md)
  with `group_var`,
  [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
  [`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md),
  [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md),
  [`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md))
  the log-likelihood has one contribution per panel, so these standard
  errors are cluster-robust; for
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  they are computed row by row and do not account for repeated
  observations of the same unit.
- `"boot"`: a nonparametric bootstrap with `bootstraps` replicates. Rows
  are resampled for
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
  and whole panels for panel models. Use it when the sample is small,
  when the Hessian is close to singular, or as a check of the other two.

``` r

f <- Total_crashes ~ lnaadt + lnlength + speed50
nb_h <- countreg(f, d, family = "NB2", method = "BFGS")
nb_r <- countreg(f, d, family = "NB2", method = "BFGS", stderr = "robust")
set.seed(1)
nb_b <- countreg(f, d, family = "NB2", method = "BFGS", stderr = "boot",
                 bootstraps = 50)
round(cbind(Hessian = sqrt(diag(vcov(nb_h))),
            robust = sqrt(diag(vcov(nb_r))),
            bootstrap = sqrt(diag(vcov(nb_b)))), 4)
#>             Hessian robust bootstrap
#> (Intercept)  0.4393 0.4855    0.4754
#> lnaadt       0.0518 0.0558    0.0541
#> lnlength     0.0694 0.0695    0.0674
#> speed50      0.1078 0.1142    0.1061
#> ln(alpha)    0.2486 0.2483    0.2304
```

Similar standard errors, as here, are reassuring. A robust standard
error that is much larger than the Hessian one points to a misspecified
distribution or dependence between observations.

The bootstrap uses only replicates whose refit converged and gave finite
estimates. The number used is stored, with the reasons for any failures:

``` r

nb_b$model$n_boot_ok
#> [1] 50
nb_b$model$boot_failures
#> 
#>                error        not converged   wrong coefficients 
#>                    0                    0                    0 
#> non-finite estimates 
#>                    0
```

If many replicates fail, increase `max.iters`, try another `method`, or
simplify the model. If fewer than two succeed, the standard errors are
`NA`.

**Weights.** Every estimator takes `weights`, which are frequency
weights: a row with weight 3 counts as three independent observations.
The estimates, log-likelihood, robust covariance and BIC are those of
the expanded data. In panel models the weights belong to the panels (a
panel with weight 3 counts as three independent panels), so they must be
constant within a panel. The bootstrap resamples the expanded
observations: it draws new frequencies from a multinomial distribution
and refits with them in place of the weights, so it needs whole-number
weights (use `stderr = "robust"` for other weights).

## Marginal effects

[`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
reports average marginal effects (or elasticities) with delta-method
standard errors, or bootstrap standard errors with `se = "bootstrap"`.
For an indicator variable the default is the pseudo-elasticity, the
percentage change in the expected count.

``` r

margEffTable(nb_h, digits = 4)
#> Marginal Effects / Elasticities
#> Method: delta 
#> 
#>      term variable_type           effect_metric estimate std_error t_value
#>    lnaadt    continuous Average marginal effect   0.4989    0.0349 14.2766
#>  lnlength    continuous Average marginal effect   0.3511    0.0380  9.2511
#>   speed50     indicator   Pseudo-elasticity (%) -41.5527    6.3016 -6.5940
#>  p_value sig lower_ci upper_ci n_obs successful_bootstraps
#>        0 ***   0.4304   0.5674  1501                    NA
#>        0 ***   0.2767   0.4255  1501                    NA
#>        0 *** -53.9035 -29.2018  1501                    NA
```

The bootstrap of
[`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
resamples the model’s independent units (rows, or panels for panel
models, or the clusters in `cluster_var`) and refits the model without
an inner bootstrap.

## Prediction

[`predict()`](https://rdrr.io/r/stats/predict.html) returns one value
per row of `newdata`, in order. Rows with missing values give `NA`;
factor levels, contrasts and data-dependent bases such as
[`poly()`](https://rdrr.io/r/stats/poly.html) and
[`splines::ns()`](https://rdrr.io/r/splines/ns.html) are those of the
training data; offset columns must be present.

For
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
models, `type` chooses the quantity:

| `type` | Quantity |
|----|----|
| `"response"` (default) | expected observed count |
| `"true"` | expected true count (with underreporting) |
| `"prob_reported"` | probability that a count is reported |
| `"link"` | linear predictor $`X\beta`$ + offset |
| `"prob"` | matrix of $`P(Y = k)`$, $`k = 0, \dots,`$`max_count` |
| `"zero"` | zero-model probability (zero-inflated or hurdle models) |

``` r

nd <- d[1:3, ]
predict(nb_h, newdata = nd)
#> [1] 0.7798765 0.7097941 1.0432200
round(predict(nb_h, newdata = nd, type = "prob", max_count = 4,
              overflow = TRUE), 3)
#>          0     1     2     3     4    >4
#> [1,] 0.502 0.307 0.127 0.044 0.014 0.006
#> [2,] 0.531 0.301 0.116 0.037 0.011 0.004
#> [3,] 0.411 0.314 0.162 0.070 0.028 0.015
```

With `overflow = TRUE` the last column holds $`P(Y >`$`max_count`$`)`$,
so each row sums to one.

For random-parameter models, `method` chooses how the random
coefficients are handled:

- `"Simulated"` (default): the average over Halton draws.
- `"Exact"`: the closed-form expectation (Wood and Gayah, 2025). It is
  not available, and the prediction is infinite, when the expected value
  does not exist (a lognormal random parameter with a positive variable,
  or a gamma random parameter beyond its moment-generating function’s
  domain); `"Simulated"` reports the same cases.
- `"Individual"` (or `"EB"`): Empirical Bayes predictions that condition
  on each panel’s observed outcomes, so `newdata` must contain the
  outcome and the panel id. These are not forecasts for new units.

Negative multinomial models have `method = "Standard"` and `"EB"` in the
same sense.

## Three different intervals

Three intervals answer three different questions, and are easily
confused:

| Interval | Question | How |
|----|----|----|
| Confidence interval for a coefficient | How uncertain is an estimated effect? | [`confint()`](https://rdrr.io/r/stats/confint.html), `tidy(conf.int = TRUE)` |
| Confidence interval for the expected count | How uncertain is the estimated mean at given covariates? | `predict(interval = "confidence")` |
| Prediction interval for a new count | Where is one new observation at given covariates likely to fall? | `predict(interval = "prediction")` |

``` r

confint(nb_h)["speed50", ]
#>      2.5 %     97.5 % 
#> -0.7483606 -0.3257283
nd2 <- data.frame(lnaadt = log(c(3000, 15000)), lnlength = 0, speed50 = 1)
predict(nb_h, newdata = nd2, interval = "confidence")
#>         fit      lwr       upr
#> 1 0.5258891 0.428403 0.6438896
#> 2 3.0025638 2.313604 3.9167160
predict(nb_h, newdata = nd2, interval = "prediction", level = 0.9)
#>         fit lwr upr
#> 1 0.5258891   0   2
#> 2 3.0025638   0   8
```

The confidence interval for the mean narrows as the sample grows; the
prediction interval does not, because a new count varies around its mean
whatever the sample size. Prediction-interval limits are whole counts
and, because counts are discrete, usually cover a little more than the
stated level. Both intervals for
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
models account for estimation uncertainty using the asymptotic normal
distribution of the estimates, evaluated at deterministic Halton points
(or a parametric bootstrap with `uncertainty = "bootstrap"`); see
[`?predict.flexCountReg`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md).

A fourth kind, a *tolerance interval*, is meant to contain at least a
stated share of all counts at given covariates (for example 90% of road
segments with these characteristics) with a stated confidence; see
[`?tolerance_interval`](https://jwood-iastate.github.io/flexCountReg/reference/tolerance_interval.md).

## Fitted values, residuals and simulation

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) gives the
predictions for the rows used in estimation, and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) the residuals of
three types:

- `"response"`: observed minus fitted.
- `"pearson"`: divided by the standard deviation implied by the fitted
  distribution
  ([`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  models).
- `"quantile"`: randomized quantile residuals (Dunn and Smyth, 1996),
  which are standard normal when the model is correct; use them for
  residual plots and normality checks.

``` r

r <- residuals(nb_h, type = "quantile")
c(mean = mean(r), sd = sd(r))
#>        mean          sd 
#> -0.00566378  0.97827719
```

[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws new
responses from the fitted distribution for the same rows (with `seed`,
reproducible and without changing the random-number stream of the
session). The parameters are held at their estimates (a plug-in
simulation), so the spread of simulated statistics reflects the
randomness of the counts, not the uncertainty of the estimates. It is
the basis of parametric-bootstrap checks:

``` r

sims <- simulate(nb_h, nsim = 200, seed = 1)
# Observed share of zeros against its simulated distribution
c(observed = mean(d$Total_crashes == 0),
  simulated = mean(colMeans(sims == 0)))
#>  observed simulated 
#> 0.7335110 0.7307528
```

[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html) and
[`augment()`](https://generics.r-lib.org/reference/augment.html) return
tidy data frames, and
[`rootogram()`](https://jwood-iastate.github.io/flexCountReg/reference/rootogram.md)
and
[`cureplot()`](https://jwood-iastate.github.io/flexCountReg/reference/cureplot.md)
compare the fitted distribution and mean with the data (see the vignette
on choosing a count distribution).

## Comparing models

### AIC and BIC

AIC and BIC compare any models fit to the same observations, nested or
not, and are always reported. BIC uses the number of independent units:
rows, panels for panel models, and the sum of the frequency weights when
weights are used.

``` r

po <- countreg(f, d, family = "POISSON", method = "BFGS")
nb1 <- countreg(f, d, family = "NB1", method = "BFGS")
pl <- countreg(f, d, family = "PL", method = "BFGS")
do.call(rbind, lapply(list(Poisson = po, NB1 = nb1, NB2 = nb_h, PL = pl),
                      function(m) glance(m)[, c("logLik", "AIC", "BIC")]))
#> # A tibble: 4 × 3
#>   logLik   AIC   BIC
#> *  <dbl> <dbl> <dbl>
#> 1 -1101. 2209. 2230.
#> 2 -1089. 2188. 2215.
#> 3 -1085. 2180. 2206.
#> 4 -1086. 2182. 2209.
```

### Likelihood-ratio tests

An LR test compares a restricted model with a larger one. The statistic
$`2(\ell_\text{full} - \ell_\text{restricted})`$ has its usual
chi-square distribution only when:

1.  **The models are nested.** The restricted model is the larger model
    with some parameters fixed: the same estimator and model components,
    the restricted model’s parameters a subset of the larger model’s,
    and the same family or a nested family.
2.  **The restriction is regular.** The fixed values are inside the
    parameter space, every remaining parameter is identified under the
    restriction, and the restriction is not reached only in a limit.
3.  **The data are the same.** The same observations, outcome, offsets,
    weights and panels, and the same values of the variables both models
    use.
4.  **Both fits reached their maximum.** Both report convergence, and
    the larger model’s log-likelihood is not below the restricted
    model’s.

The families nest as follows:

| Restricted | Larger | Restriction | Null distribution |
|----|----|----|----|
| Poisson | NB1, NB2, PLN, PIG1, PIG2, PIG | dispersion = 0 | 50:50 mixture of $`\chi^2_{df-1}`$ and $`\chi^2_{df}`$ |
| Poisson | GP1, GP2 | dispersion = 0 | $`\chi^2_{df}`$ |
| Poisson | COM | $`\nu = 1`$ | $`\chi^2_{df}`$ |
| Poisson | DPO, HYPERPOIS, GAMMACOUNT | $`\phi = 1`$, $`\gamma = 1`$, shape $`= 1`$ | $`\chi^2_{df}`$ |
| NB1, NB2 | NBP | $`P = 1`$ or $`P = 2`$ | $`\chi^2_{df}`$ |
| PL | PLG, PLL | $`\alpha = 0`$, $`\sigma = 0`$ | 50:50 mixture |

For NB1, NB2 and the other families in the first row, a dispersion of 0
is on the boundary of the parameter space (a dispersion cannot be
negative). The statistic then has a 50:50 mixture distribution, and the
ordinary chi-square p-value would be twice too large.

These comparisons are **not** valid, and no p-value is reported:

- **Poisson vs PL, PLG or PLL.** The Poisson-Lindley family does not
  contain the Poisson distribution: its variance is at least
  $`\mu + \mu^2/2`$ for every $`\theta`$.
- **Poisson vs NBP or SI.** At the Poisson limit, NBP’s power $`P`$ and
  the Sichel $`\gamma`$ are not identified, so the usual asymptotics
  fail.
- **Poisson vs PW, PGE or GW.** Their distributions approach Poisson
  only as a parameter goes to infinity.
- **Poisson vs DELAPORTE, HERMITE, PA or NBGE.** Poisson is reached only
  at the edge of the parameter space (a share of 0 on the logit scale,
  or a dispersion going to zero), where other parameters are not
  identified.
- **A boundary family with a dispersion formula.** The restriction sends
  an intercept to minus infinity while the slopes become unidentified.
- **Models with and without a zero model or underreporting.** No zero
  inflation, or full reporting, is reached only in a limit (an intercept
  of minus or plus infinity on the logit scale).
- **Models with different random parameters.** A random-parameter
  variance of zero is on the boundary, often with other parameters
  unidentified.
- **Non-nested models**, such as NB1 vs NB2.

For these, compare models with AIC and BIC, or with a parametric
bootstrap of the LR statistic using
[`simulate()`](https://rdrr.io/r/stats/simulate.html).

[`summary()`](https://rdrr.io/r/base/summary.html) reports a test
against Poisson when it is valid:

``` r

summary(nb_h)$lr_test
#> $statistic
#> [1] 31.31106
#> 
#> $df
#> [1] 1
#> 
#> $p.value
#> [1] 1.099119e-08
#> 
#> $logLik_poisson
#> [1] -1100.597
#> 
#> $boundary
#> [1] TRUE
is.null(summary(pl)$lr_test)   # no test: PL does not nest Poisson
#> [1] TRUE
```

When a condition fails,
[`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md)
reports no LR test and gives the reason in `LR_reason`. Nesting of the
mean model is checked on the design matrices, so `y ~ x` is nested in
`y ~ I(2 * x) + z`.

[`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md)
applies the same rules to two fitted models (the one with fewer
parameters is the restricted model) or to a model and a base model it
fits on the same rows:

``` r

r <- regCompTest(po, model2 = nb_h)
r[c("LR", "LRdof", "LR_type", "LR_pvalue")]
#> $LR
#> [1] 31.31106
#> 
#> $LRdof
#> [1] 1
#> 
#> $LR_type
#> [1] "boundary"
#> 
#> $LR_pvalue
#> [1] 1.099119e-08
r <- regCompTest(nb1, model2 = countreg(f, d, family = "NBP",
                                        method = "BFGS"))
r[c("LR", "LRdof", "LR_type", "LR_pvalue")]
#> $LR
#> [1] 9.828597
#> 
#> $LRdof
#> [1] 1
#> 
#> $LR_type
#> [1] "interior"
#> 
#> $LR_pvalue
#> [1] 0.001718194
```

When the larger model has the lower log-likelihood, the larger model did
not reach its maximum (check its convergence), and no p-value is given.

### McFadden’s pseudo-R²

[`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md)
and
[`regCompTable()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTable.md)
report McFadden’s pseudo-R², $`1 - \ell / \ell_0`$, against an
intercept-only Poisson model on the same rows by default, so it is
comparable across families. `basemodel` and `variables` choose another
baseline.

``` r

regCompTable(list(Poisson = po, NB2 = nb_h))
#> # A tibble: 10 × 3
#>    Parameter    Poisson           NB2              
#>    <chr>        <chr>             <chr>            
#>  1 (Intercept)  -9.011 (0.417)*** -8.772 (0.439)***
#>  2 lnaadt       1.108 (0.048)***  1.082 (0.052)*** 
#>  3 lnlength     0.739 (0.059)***  0.762 (0.069)*** 
#>  4 speed50      -0.516 (0.097)*** -0.537 (0.108)***
#>  5 ln(alpha)    ---               -1.045 (0.249)***
#>  6 N Obs.       1501              1501             
#>  7 LL           -1100.597         -1084.942        
#>  8 AIC          2209.195          2179.884         
#>  9 BIC          2230.45           2206.453         
#> 10 Pseudo-R-Sq. 0.278             0.288
```

## Next steps

- [`vignette("getting-started", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/getting-started.md)
  applies these tools in a complete analysis.
- [`vignette("diagnosing-count-models", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/diagnosing-count-models.md)
  covers residuals, rootograms and held-out comparisons.
- [`vignette("random-parameters-and-panels", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/random-parameters-and-panels.md)
  explains the independent unit for repeated observations: the
  `washington_roads` rows used here are segment-years, so the row-level
  standard errors above treat repeated years of a segment as
  independent.
- [`vignette("predictions-and-effects", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/predictions-and-effects.md)
  interprets predictions and marginal effects in detail.
- [`vignette("prediction-for-panel-data", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/prediction-for-panel-data.md)
  separates predictions for new and for known units in panels.

## References

Dunn, P. K., & Smyth, G. K. (1996). Randomized quantile residuals.
*Journal of Computational and Graphical Statistics*, 5(3), 236-244.

Self, S. G., & Liang, K.-Y. (1987). Asymptotic properties of maximum
likelihood estimators and likelihood ratio tests under nonstandard
conditions. *Journal of the American Statistical Association*, 82(398),
605-610.

Wood, J. S., & Gayah, V. (2025). Out-of-sample prediction and
interpretation for random parameter generalized linear models. *Accident
Analysis and Prevention*, 220, 108147.
