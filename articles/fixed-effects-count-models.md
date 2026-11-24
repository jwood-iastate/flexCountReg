# Fixed-effects count models

Panel count data, such as crashes on the same road segments over several
years, contain differences between units that the regressors do not
capture: driver populations, geometry, enforcement. Random-effects
models
([`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
[`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md),
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
with a random intercept) treat these differences as draws from a
distribution that is independent of the regressors. When that is
doubtful (busier roads may also be different roads) a **fixed-effects**
model gives each unit its own intercept, so the slopes are estimated
only from changes *within* units.

[`countreg_fe()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_fe.md)
estimates fixed-effects count models. This article shows when to use it,
what it can and cannot estimate, and how its results compare with pooled
and random-effects models.

``` r

library(flexCountReg)
data("washington_roads")
```

## 1. The model

For unit \\i\\ (here a road segment) in period \\t\\,

\\ E\[y\_{it} \mid x\_{it}, \alpha_i\] = \exp(\alpha_i + x\_{it}'\beta +
o\_{it}), \\

with a fixed effect \\\alpha_i\\ for each unit and an optional offset
\\o\_{it}\\. Nothing is assumed about how \\\alpha_i\\ relates to
\\x\_{it}\\.

[`countreg_fe()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_fe.md)
chooses the estimator from the family and the number of fixed-effect
dimensions:

| Family | Fixed effects | Method | Notes |
|:---|:---|:---|:---|
| Poisson | one (`~ site`) | conditional likelihood | the effects are conditioned out; consistent for fixed T with only a correct mean (Wooldridge 1999) |
| Poisson | several (`~ site + year`) | IRLS with alternating projections | the method of [`fixest::fepois()`](https://lrberge.github.io/fixest/reference/feglm.html) |
| any other | one | concentrated likelihood | each unit’s intercept is solved for given slopes; small-T bias in the dispersion |

Poisson fixed effects are the usual first choice: the slopes are
consistent even when the counts are overdispersed or not Poisson at all,
as long as the mean is right, and the default cluster-robust standard
errors account for the overdispersion.

## 2. Site fixed effects

Here `washington_roads` has 507 segments, most observed in three years.
Traffic volume (`lnaadt`) changes from year to year within a segment,
but the speed indicator `speed50` does not:

``` r

fe <- countreg_fe(Total_crashes ~ lnaadt + speed50 + offset(lnlength),
                  fixef = ~ ID, data = washington_roads)
#> 804 row(s) in 273 group(s) with all-zero counts or a single observation were dropped: they carry no information on the slopes.
#> Warning: Regressor(s) without variation within the fixed-effect groups are not
#> identified and were removed: speed50.
summary(fe)
#> Call:
#> countreg_fe(formula = Total_crashes ~ lnaadt + speed50 + offset(lnlength), 
#>     fixef = ~ID, data = washington_roads)
#> 
#> Model: countreg_fe | family: POISSON | optimizer: conditional
#> Poisson fixed effects (conditional likelihood) 
#> Groups dropped (all-zero counts or one observation): 273 
#> Iterations: 4 
#> Convergence: converged 
#> 
#> Coefficients (standard errors: cluster-robust (234 clusters of ID)):
#>        Estimate Std. Error z value Pr(>|z|)
#> lnaadt   -0.543      0.491    -1.1     0.27
#> Offset: offset(lnlength) (coefficient fixed at 1)
#> 
#> Log-likelihood: -735.347 | AIC: 1940.695 | BIC: 3009.189 
#> N: 697 | panels: 234
```

Two things happened before estimation, and both are reported:

- **Segments with no crashes in any year were dropped.** Their fixed
  effect would be \\-\infty\\; they carry no information on the slopes.
  So were segments observed only once.
- **`speed50` was removed.** A regressor that is constant within every
  unit is absorbed by the unit’s fixed effect, so its coefficient is not
  identified. Time-invariant characteristics can only be studied with
  random-effects or pooled models.

``` r

length(fe$model$dropped$ID)     # segments dropped
#> [1] 273
fe$model$not_identified         # regressors removed
#> [1] "speed50"
```

The standard errors are cluster-robust by segment (`stderr = "cluster"`,
the default). Other covariance matrices are available after the fit:

``` r

sqrt(diag(vcov(fe)))                      # clustered by segment
#>    lnaadt 
#> 0.4912069
sqrt(diag(vcov(fe, type = "hessian")))    # Poisson model-based
#>    lnaadt 
#> 0.6052379
```

The model-based standard error is much smaller: it assumes the Poisson
variance, which these counts exceed.

## 3. Site and year fixed effects

Crash counts often share year-to-year changes (weather, reporting, the
economy). A second fixed-effect dimension absorbs them:

``` r

fe2 <- countreg_fe(Total_crashes ~ lnaadt + offset(lnlength),
                   fixef = ~ ID + Year, data = washington_roads)
#> 804 row(s) in 273 group(s) with all-zero counts or a single observation were dropped: they carry no information on the slopes.
summary(fe2)
#> Call:
#> countreg_fe(formula = Total_crashes ~ lnaadt + offset(lnlength), 
#>     fixef = ~ID + Year, data = washington_roads)
#> 
#> Model: countreg_fe | family: POISSON | optimizer: demean
#> Poisson fixed effects (IRLS with alternating projections) 
#> Groups dropped (all-zero counts or one observation): 273 
#> Iterations: 5 
#> Convergence: converged 
#> 
#> Coefficients (standard errors: cluster-robust (234 clusters of ID)):
#>        Estimate Std. Error z value Pr(>|z|)
#> lnaadt   -0.576      0.483   -1.19     0.23
#> Offset: offset(lnlength) (coefficient fixed at 1)
#> 
#> Log-likelihood: -735.118 | AIC: 1944.236 | BIC: 3021.824 
#> N: 697 | panels: 234
fe2$model$fixef$fe2                       # year effects
#>        2016        2017        2018 
#>  0.01735905 -0.03685296  0.01987559
```

With few levels, a dimension like `Year` can also enter the formula as a
factor; with several high-dimensional effects (sites and corridor-years,
say), the alternating projections keep the computation fast.

## 4. Comparison with pooled and random-effects models

The elasticity of crashes with respect to traffic volume from the
pooled, random-effects and fixed-effects Poisson models:

``` r

f <- Total_crashes ~ lnaadt + offset(lnlength)
pooled <- countreg(f, data = washington_roads, family = "POISSON",
                   method = "BFGS", cluster = ~ ID)
re <- repois(f, data = washington_roads, group_var = "ID")
est <- function(m) c(estimate = coef(m)[["lnaadt"]],
                     se = sqrt(diag(vcov(m)))[["lnaadt"]])
round(rbind(pooled = est(pooled), random_effects = est(re),
            fixed_effects = est(fe2)), 3)
#>                estimate    se
#> pooled            1.196 0.078
#> random_effects    1.144 0.061
#> fixed_effects    -0.576 0.483
```

The pooled and random-effects estimates use the differences *between*
segments in traffic volume, which are large; the fixed-effects estimate
uses only the changes *within* segments over three years, which are
small, so its standard error is much larger. A big gap between the
random-effects and fixed-effects estimates suggests that the segment
effects are related to traffic volume, which would bias the
random-effects estimate.

## 5. Other families

Any family of
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
(and any family made with
[`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md))
can be used with one fixed-effect dimension. The negative binomial (NB2)
fixed-effects model is the dummy-variable NB2, with the unit intercepts
concentrated out.

For the road segments, NB2 fixed effects reduce to the Poisson model:
once each segment has its own intercept, no overdispersion remains in
three years of counts, and the estimate of \\\alpha\\ goes to 0 (the fit
warns that the likelihood is flat in `ln(alpha)`). The Poisson
fixed-effects model of Section 2 is then the right one.

To show the NB2 model where it matters, simulate 300 units observed for
5 periods, with unit effects correlated with `x1` and overdispersed
counts (NB2 with \\\alpha = 0.5\\):

``` r

set.seed(11)
G <- 300; n_t <- 5
sim <- data.frame(id = rep(seq_len(G), each = n_t),
                  t = rep(seq_len(n_t), G))
a_i <- rnorm(G, 0, 0.7)
sim$x1 <- rnorm(nrow(sim)) + 0.5 * a_i[sim$id]
sim$x2 <- rnorm(nrow(sim))
sim$y <- rnbinom(nrow(sim), size = 2,
                 mu = exp(0.5 + a_i[sim$id] + 0.6 * sim$x1 - 0.4 * sim$x2))

fe_nb <- countreg_fe(y ~ x1 + x2, fixef = ~ id, data = sim, family = "NB2")
#> 55 row(s) in 11 group(s) with all-zero counts or a single observation were dropped: they carry no information on the slopes.
fe_nb_jk <- countreg_fe(y ~ x1 + x2, fixef = ~ id, data = sim,
                        family = "NB2", jackknife = TRUE, time = "t")
#> 55 row(s) in 11 group(s) with all-zero counts or a single observation were dropped: they carry no information on the slopes.
round(rbind(truth = c(0.6, -0.4, log(0.5)),
            NB2_FE = coef(fe_nb),
            jackknife = coef(fe_nb_jk)), 3)
#>              x1     x2 ln(alpha)
#> truth     0.600 -0.400    -0.693
#> NB2_FE    0.555 -0.416    -1.238
#> jackknife 0.550 -0.420    -0.707
```

The slopes are close to the truth, but the dispersion \\\ln\alpha\\ is
underestimated, because the fixed effects absorb much of the variation
between units. Two cautions apply to non-Poisson fixed effects in short
panels:

- The slopes have a small bias of order \\1/T\\ (Allison and Waterman
  2002); the dispersion has a large one.
- `jackknife = TRUE` applies the split-panel jackknife of Dhaene and
  Jochmans (2015). In simulations like this one, with 200 units, it
  removed most of the dispersion bias with 5 or 10 periods but did not
  improve the slopes, and with 3 periods its half-panels are too short
  (a warning is given below 4 periods).

The Hausman, Hall and Griliches “conditional negative binomial” is a
different model, whose fixed effect enters the dispersion, not the mean;
it is not provided.

## 6. Correlated random effects: a test of random against fixed effects

Random-effects models can be kept, and the question of section 4 tested,
with the device of Mundlak (1978): add each unit’s mean of the
time-varying regressors to the model. The random effect may then depend
on the regressors through their means, and the coefficient of the
regressor itself is a within-unit effect, as in the fixed-effects model.
`cre =` does this in
[`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md),
[`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)
and
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md):

``` r

cre <- repois(Total_crashes ~ lnaadt + offset(lnlength), group_var = "ID",
              data = washington_roads, cre = ~ lnaadt, method = "BFGS",
              stderr = "robust")
round(coef(cre)[c("lnaadt", "mean_lnaadt")], 3)
#>      lnaadt mean_lnaadt 
#>      -0.706       1.861
cre_test(cre)
#> 
#>  Wald test of correlated random effects (group means)
#> 
#> data:  repois(formula = Total_crashes ~ lnaadt + offset(lnlength), group_var = "ID",      data = washington_roads, method = "BFGS", stderr = "robust",      cre = ~lnaadt)
#> Wald chi-squared = 7.9557, df = 1, p-value = 0.004794
#> sample estimates:
#> mean_lnaadt 
#>    1.860804
```

[`cre_test()`](https://jwood-iastate.github.io/flexCountReg/reference/cre_test.md)
is a Wald test that the coefficients of the means are zero. With robust
standard errors it is a robust version of the Hausman test (Wooldridge
2019). A small p-value says that the segment effects are related to
traffic volume: the plain random-effects estimate of section 4 is then
biased, and the within-segment effect (the coefficient of `lnaadt` here,
or the fixed-effects estimate) is the one to report.

In simulations with 300 units and 5 periods, when the unit effects were
correlated with the regressor the plain random-effects slope was biased
(0.63 for a true 0.5), while the model with `cre` gave 0.505, equal to
the fixed-effects estimate, and
[`cre_test()`](https://jwood-iastate.github.io/flexCountReg/reference/cre_test.md)
rejected in every replication. With uncorrelated effects all three were
unbiased and the test rejected in 7% of replications at the 5% level.

The group means are computed on the estimation rows and stored with the
fit, so predictions for a known unit use its training means: rows of
`newdata` never change them.

## 7. Predictions

Predictions for a known segment use its estimated fixed effect. A
segment that was not in the estimation data has no fixed effect, so its
prediction is `NA` unless you ask for the mean effect explicitly:

``` r

nd <- washington_roads[washington_roads$ID %in%
                         names(fe2$model$fixef$fe1)[1:2], ]
predict(fe2, newdata = nd)
#> [1] 0.3414432 0.3402894 0.3244072 0.3233158 0.3341496 0.3363948

new_site <- nd[1, ]
new_site$ID <- "not-in-data"
predict(fe2, newdata = new_site)
#> 1 row(s) belong to groups not in the estimation data, which have no fixed effect: their predictions are NA (use newgroup = "mean_fe" to use the mean effect).
#> [1] NA
predict(fe2, newdata = new_site, newgroup = "mean_fe")
#> [1] 0.5565084
```

## 8. Agreement with other software

The Poisson estimates, log-likelihoods and dropped groups equal those of
[`fixest::fepois()`](https://lrberge.github.io/fixest/reference/feglm.html),
and the NB2 estimates those of
[`fixest::fenegbin()`](https://lrberge.github.io/fixest/reference/femlm.html)
and of [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) with
segment dummies. The cluster-robust standard errors use the factor
\\G/(G-1)\\ for \\G\\ clusters; they equal fixest’s with
`ssc(adj = FALSE)`, as fixest by default also multiplies the variance by
\\(n-1)/(n-K)\\.

## References

Allison, P.D. and Waterman, R.P. (2002). Fixed-effects negative binomial
regression models. *Sociological Methodology*, 32, 247-265.

Dhaene, G. and Jochmans, K. (2015). Split-panel jackknife estimation of
fixed-effect models. *Review of Economic Studies*, 82(3), 991-1030.

Mundlak, Y. (1978). On the pooling of time series and cross section
data. *Econometrica*, 46(1), 69-85.

Wooldridge, J.M. (1999). Distribution-free estimation of some nonlinear
panel data models. *Journal of Econometrics*, 90(1), 77-97.

Wooldridge, J.M. (2019). Correlated random effects models with
unbalanced panels. *Journal of Econometrics*, 211(1), 137-150.
