# Hierarchical count models

Counts are often grouped: crash counts by road segment and year,
segments within counties, animal counts by site and survey. Observations
in the same group share unmeasured conditions, so they are correlated,
and the groups themselves are of interest.
[`countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md)
fits **hierarchical** (mixed-effects, multilevel) count models with
normally distributed random effects at any number of grouping levels,
nested or crossed, for any count family of the package.

``` r

library(flexCountReg)
data("stream_salamanders")
s <- stream_salamanders
s$sample <- factor(s$sample)
```

The salamander data have 644 counts from 23 stream sites, each surveyed
four times (`sample`), for seven species and life stages (`spp`), at
sites affected or not by mountaintop mining (`mined`).

## 1. A random intercept

Random-effect terms are written as in lme4: `(1 | site)` gives each site
its own intercept, \\b_j \sim N(0, \sigma^2)\\:

``` r

h <- countreg_hier(count ~ spp + mined + (1 | site), data = s,
                   family = "NB2")
summary(h)
#> Hierarchical count model: NB2
#> Integration: adaptive Gauss-Hermite quadrature, 9 points
#> Call: countreg_hier(formula = count ~ spp + mined + (1 | site), data = s, family = "NB2")
#> 
#> Random effects:
#>  Groups Name        Variance Std.Dev.
#>  site   (Intercept) 0.2968   0.5447  
#> Number of obs: 644; groups: site, 23
#> 
#> Fixed coefficients and distribution parameters:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.58008    0.22595   2.567 0.010250 *  
#> sppPR       -1.31981    0.28752  -4.590 4.43e-06 ***
#> sppDM        0.36856    0.22349   1.649 0.099120 .  
#> sppEC-A     -0.70982    0.25300  -2.806 0.005022 ** 
#> sppEC-L      0.57134    0.21912   2.607 0.009124 ** 
#> sppDES-L     0.79285    0.21668   3.659 0.000253 ***
#> sppDF        0.31202    0.23294   1.340 0.180401    
#> minedyes    -2.26385    0.28455  -7.956 1.78e-15 ***
#> ln(alpha)    0.05965    0.13501   0.442 0.658617    
#> 
#> Random-effect standard deviations and correlations (Wald SEs):
#>                     Estimate Std. Error
#> sd_(Intercept)|site   0.5447     0.1329
#> 
#> Log-likelihood: -826.1554   AIC: 1672.311   BIC: 1716.988
```

The model is

\\ \log \mu\_{ij} = x\_{ij}'\beta + b_j, \qquad y\_{ij} \mid b_j \sim
\text{NB2}(\mu\_{ij}, \alpha), \qquad b_j \sim N(0, \sigma^2), \\

and the likelihood integrates each site’s \\b_j\\ out. With one grouping
factor the integral is computed by **adaptive Gauss-Hermite quadrature**
(9 points by default), which is accurate even for small groups.

The standard deviation 0.54 means that sites with the same covariates
differ in their expected counts by a factor of about \\\exp(2 \times
0.54) \approx 3\\ between the 16th and 84th percentiles.
[`re_lrtest()`](https://jwood-iastate.github.io/flexCountReg/reference/re_lrtest.md)
tests whether the variance is zero. The null value is on the boundary of
the parameter space, so the p-value is half that of the usual
\\\chi^2_1\\ test:

``` r

re_lrtest(h)
#>         term logLik_without logLik_with statistic df      p.value
#> 1 (1 | site)       -837.653   -826.1554  22.99518  1 8.120391e-07
```

## 2. Effects of the groups

[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) returns
the conditional modes of the site effects, the most likely value of each
site’s effect given its counts. They are shrunk towards zero, more so
for sites with little information.
[`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) returns the
estimated variance components and
[`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html) the fixed
coefficients:

``` r

head(ranef(h)$site)
#>     (Intercept)
#> R-1  0.20312721
#> R-2  0.08408975
#> R-3 -0.05997694
#> R-4  0.47914042
#> R-5  0.09240978
#> R-6  0.38168445
VarCorr(h)
#>  Groups Name        Variance Std.Dev.
#>  site   (Intercept) 0.2968   0.5447
```

## 3. Crossed and nested effects

Several terms can be combined. Here the site effects are **crossed**
with survey occasions shared by all sites:

``` r

h2 <- countreg_hier(count ~ spp + mined + (1 | site) + (1 | sample),
                    data = s, family = "POISSON")
VarCorr(h2)
#>  Groups Name        Variance Std.Dev.
#>  site   (Intercept) 0.33500  0.5788  
#>  sample (Intercept) 0.01012  0.1006
```

**Nested** effects, such as surveys within sites, are written
`(1 | site/sample)`, which is short for
`(1 | site) + (1 | site:sample)`. With more than one term the model uses
the **Laplace approximation** with sparse matrices, as
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) and `glmmTMB`
do. It handles thousands of groups at several levels.

## 4. Random slopes

`(1 + cover | site)` gives each site its own intercept and its own
effect of `cover`, with a correlation between them. `(cover || site)`
makes them uncorrelated:

``` r

h3 <- countreg_hier(count ~ cover + mined + (cover || site), data = s,
                    family = "NB2")
#> Singular fit: a random-effect variance is estimated at or near zero (see summary()).
VarCorr(h3)
#>  Groups Name        Variance  Std.Dev. 
#>  site   (Intercept) 2.626e-01 5.125e-01
#>  site   cover       1.913e-16 1.383e-08
```

The slope’s standard deviation is estimated at zero, so `cover` has the
same effect at every site. A variance on the boundary is a **singular
fit**: [`summary()`](https://rdrr.io/r/base/summary.html) says so, and
the Wald standard errors of that term’s variance components are not
reported, since they are not meaningful there. Remove the term, or test
it with
[`re_lrtest()`](https://jwood-iastate.github.io/flexCountReg/reference/re_lrtest.md).

Group-level predictors (such as `mined`, constant within sites) are
ordinary columns of the formula, and cross-level interactions such as
`mined:cover` are ordinary interactions.

## 5. Predictions

[`predict()`](https://rdrr.io/r/stats/predict.html) conditions on the
estimated effects of the groups by default. `re.form` chooses which
levels to condition on (`NA` for none), and groups not seen in
estimation get no effect:

``` r

nd <- s[c(1, 30), ]
predict(h2, newdata = nd)                         # site and sample
#> [1] 0.2103916 3.1660525
predict(h2, newdata = nd, re.form = ~ (1 | site)) # site only
#> [1] 0.2127103 3.4092077
predict(h2, newdata = nd, re.form = NA)           # population level
#> [1] 0.1958288 1.8879303
```

A population-level prediction with the effects at zero is the expected
count of a *median* site. The **mean** over sites is larger: for the log
link it is \\\exp(\eta + \sigma^2/2)\\. `marginal = TRUE` integrates the
effects that are not conditioned on out of the mean:

``` r

predict(h2, newdata = nd, re.form = NA, marginal = TRUE)
#> [1] 0.2327066 2.2434593
```

[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws new counts,
either for the same groups (`use_re = "conditional"`) or for new groups
with new effects (`use_re = "new"`).

## 6. Predictions from a separate history

To forecast a new period, the effects of the groups should come from
their past counts, not from the counts being forecast. `history` gives
the past data: the effects of its groups are computed from those counts
alone, at the fitted parameters, and groups without history are treated
as new groups.

``` r

data("washington_roads")
w <- washington_roads
past <- w[w$Year < 2018, ]
h_w <- countreg_hier(Total_crashes ~ lnaadt + speed50 + offset(lnlength) +
                       (1 | ID), data = past, family = "POISSON")
next_year <- w[w$Year == 2018, ]
last_year <- past[past$Year == 2017, ]
head(cbind(
  all_years = predict(h_w, newdata = next_year),
  last_year = predict(h_w, newdata = next_year, history = last_year),
  observed  = next_year$Total_crashes))
#>      all_years last_year observed
#> [1,] 0.5105095 0.5872832        1
#> [2,] 0.7575946 0.5288174        3
#> [3,] 1.0443883 0.8029405        0
#> [4,] 0.3603831 0.3997441        1
#> [5,] 0.6745255 0.8043373        0
#> [6,] 0.3803613 0.3131675        1
```

## 7. Intervals for the variance components

Wald intervals for a standard deviation are symmetric and can cross
zero. Profile-likelihood intervals follow the likelihood, so they are
asymmetric and stop at zero:

``` r

confint(h, "sd_(Intercept)|site")                      # Wald
#>                         2.5 %    97.5 %
#> sd_(Intercept)|site 0.2842245 0.8052748
confint(h, "sd_(Intercept)|site", method = "profile")  # profile
#>                         2.5 %    97.5 %
#> sd_(Intercept)|site 0.3273034 0.8768668
```

Each profile point refits the model, so a profile interval takes several
fits per limit.

## 8. Smooth terms

The fixed part can contain penalized smooths, `s()` and `t2()`. They are
fitted in their mixed-model form (as in gamm4): the wiggly part of each
smooth is a set of normal random coefficients whose variance is
\\1/\lambda\\, so the smoothness is estimated by the marginal likelihood
together with the other variance components.

``` r

hs <- countreg_hier(count ~ mined + s(cover) + (1 | site), data = s,
                    family = "POISSON")
hs$model$smooth_table[, c("term", "edf", "p.value")]
#>       term edf   p.value
#> 1 s(cover)   1 0.4639229
```

[`plot_smooths()`](https://jwood-iastate.github.io/flexCountReg/reference/plot_smooths.md),
[`smooth_estimates()`](https://jwood-iastate.github.io/flexCountReg/reference/smooth_estimates.md)
and the other smooth tools of
[`vignette("smooth-effects")`](https://jwood-iastate.github.io/flexCountReg/articles/smooth-effects.md)
work as for
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md).
Here the effect of `cover` is linear (EDF 1).

## 9. Non-normal random effects

With one grouping term, the effects can have other distributions,
through `re_dist`, and the likelihood is computed by simulation
(`engine = "msl"`, `re_draws` draws per group):

- `"uniform"` and `"triangular"`: effects with mean zero and the
  reported standard deviation, with bounded support;
- `"gamma"`: the group’s multiplicative effect \\\exp(b_j)\\ is gamma
  with mean 1, as in the classic negative binomial random-effects model.
  With a Poisson family this is
  [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md),
  computed by simulation.

``` r

hg <- countreg_hier(count ~ spp + mined + (1 | site), data = s,
                    family = "POISSON", re_dist = "gamma")
rp <- repois(count ~ spp + mined, group_var = "site", data = s,
             method = "BFGS")
c(simulated = logLik(hg), closed_form = logLik(rp))
#>   simulated closed_form 
#>   -972.6204   -972.5913
```

Group effects from the simulation engine are posterior means. The gamma
distribution needs a gamma quantile for every draw, so it is slower than
the others.

## 10. Correlated random effects

Random effects are assumed independent of the covariates. `cre` adds the
group means of time-varying covariates (Mundlak’s device), which absorbs
that correlation, and
[`cre_test()`](https://jwood-iastate.github.io/flexCountReg/reference/cre_test.md)
tests it:

``` r

hc <- countreg_hier(Total_crashes ~ lnaadt + speed50 + offset(lnlength) +
                      (1 | ID), data = w, family = "POISSON", cre = ~ lnaadt)
cre_test(hc)
#> 
#>  Wald test of correlated random effects (group means)
#> 
#> data:  countreg_hier(formula = Total_crashes ~ lnaadt + speed50 + offset(lnlength) +      (1 | ID), data = w, family = "POISSON", cre = ~lnaadt)
#> Wald chi-squared = 7.1607, df = 1, p-value = 0.007452
#> sample estimates:
#> mean_lnaadt 
#>    1.855634
```

## 11. Families and engines

Every count family can be used, including user-defined families
([`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md)),
and the distribution parameters can depend on covariates
(`dis_param_formula_1`, `dis_param_formula_2`). For the Poisson, NB1 and
NB2 families the derivatives of the log-density are analytic; for the
others they are computed numerically.

`engine` chooses the integration method:

| engine | when | method |
|----|----|----|
| `"aghq"` | one grouping term | adaptive Gauss-Hermite quadrature, `nAGQ` points per effect |
| `"laplace"` | any terms, nested or crossed; smooth terms | Laplace approximation |
| `"msl"` | one grouping term, any `re_dist` | maximum simulated likelihood |
| `"auto"` (default) |  | `"msl"` for non-normal effects, `"aghq"` for one term with at most two effects, otherwise `"laplace"` |

With `nAGQ = 1` the quadrature is the Laplace approximation. The
estimates equal those of
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) (Poisson) and
`glmmTMB` (Poisson, NB1 and NB2) with the same approximation.

## 12. Related models

- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  fits random parameters at one level, with free means and non-normal
  distributions (lognormal, triangular, gamma, …) and heterogeneity in
  their means and variances, by simulation.
- [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
  [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md)
  and
  [`poislind_re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)
  fit random effects with conjugate distributions in closed form.
- [`countreg_fe()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_fe.md)
  fits fixed-effects models, which make no assumption about how the
  group effects relate to the covariates.

## 13. Not yet available

Non-normal random effects at more than one level (nested simulation),
smooth terms in the distribution-parameter formulas, and `te()` and
`ti()` smooths (use `t2()`) are not yet available in
[`countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md).
