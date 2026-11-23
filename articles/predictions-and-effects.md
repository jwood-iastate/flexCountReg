# Interpreting predictions, marginal effects and elasticities

A count model’s coefficients are on the log scale, but questions are
usually about counts: *how many more absences would we expect if …?*
This vignette shows how to turn a fitted model into answers on the count
scale, for a single profile or averaged over a population, and how to
get their uncertainty.

``` r

library(flexCountReg)
```

## Four quantities

For a fitted mean $`\mu(x)`$ and a covariate $`x_j`$:

| Quantity | Definition | Units |
|----|----|----|
| Prediction | $`\mu(x)`$ at a stated profile $`x`$ | counts |
| Contrast | $`\mu(x^{(1)}) - \mu(x^{(0)})`$ for two profiles (e.g. two levels of a factor) | counts |
| Marginal effect | $`\partial\mu/\partial x_j`$ for a continuous $`x_j`$ | counts per unit of $`x_j`$ |
| Elasticity | $`(\partial\mu/\partial x_j)(x_j/\mu)`$ | % change in the count per % change in $`x_j`$ |

Each can be evaluated at one *profile* (a stated set of covariate
values) or *averaged* over the units of a population. The two generally
differ: with a log link, the average of the predictions is not the
prediction at average covariate values. State which one is reported.

## Categorical predictors: school absences

[`MASS::quine`](https://rdrr.io/pkg/MASS/man/quine.html) records the
number of days absent from school (`Days`) in a school year for 146
children in an Australian town, with four categorical predictors: ethnic
background `Eth` (A = Aboriginal, N = not), `Sex`, age group `Age` (F0
primary, F1 to F3 forms 1 to 3) and learner status `Lrn` (AL average, SL
slow) (Quine, 1975). The categories are those of the original study. The
associations below describe these data; they are not causal effects.

``` r

q <- MASS::quine
nb <- countreg(Days ~ Eth + Sex + Age + Lrn, data = q, family = "NB2",
               method = "BFGS")
summary(nb)$coefficients[, 1:2]
#>                Estimate Std. Error
#> (Intercept)  2.89457803  0.2279260
#> EthN        -0.56937207  0.1576087
#> SexM         0.08232043  0.1646847
#> AgeF1       -0.44842604  0.2376017
#> AgeF2        0.08808263  0.2415475
#> AgeF3        0.35690340  0.2466199
#> LrnSL        0.29210918  0.1829368
#> ln(alpha)   -0.24286194  0.1263215
```

The reference levels are the first level of each factor (A, F, F0, AL).

### A profile prediction and an average prediction

``` r

ref <- data.frame(Eth = "A", Sex = "F", Age = "F0", Lrn = "AL")
ref[] <- Map(function(v, nm) factor(v, levels = levels(q[[nm]])), ref,
             names(ref))
c(reference_profile = predict(nb, newdata = ref),
  average_over_children = mean(predict(nb)))
#>     reference_profile average_over_children 
#>              18.07587              16.47124
```

The first number is the expected absence of one hypothetical child in
the reference categories; the second averages the predictions over the
146 children in the data. They answer different questions.

### Rate ratios and differences

With a log link, $`\exp(\beta)`$ is a rate ratio, the same for every
profile. For `Eth`, $`\exp(\beta)`$ = 0.566: children of non-Aboriginal
background are expected to have 43.4% fewer absent days, other
categories held fixed. The *difference* in days depends on the profile:

``` r

d_A <- d_N <- q
d_A$Eth <- factor("A", levels = levels(q$Eth))
d_N$Eth <- factor("N", levels = levels(q$Eth))
diff_by_child <- predict(nb, newdata = d_N) - predict(nb, newdata = d_A)
c(average_difference = mean(diff_by_child),
  smallest = min(diff_by_child), largest = max(diff_by_child))
#> average_difference           smallest            largest 
#>          -9.246071         -12.461632          -5.011404
```

[`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
computes such average discrete changes (and their delta-method standard
errors) for two-level factors and continuous variables. Its result
agrees with the hand calculation:

``` r

me <- margEffTable(nb, vars = c("Eth", "Sex", "Lrn"),
                   indicator = "discrete_change", tableType = "tibble")
me[, c("term", "effect_metric", "estimate", "std_error", "lower_ci",
       "upper_ci")]
#> Marginal Effects / Elasticities
#> Method: delta 
#> 
#>  term   effect_metric estimate std_error lower_ci upper_ci
#>   Eth Discrete change   -9.246     2.751  -14.638   -3.854
#>   Sex Discrete change    1.360     2.729   -3.989    6.709
#>   Lrn Discrete change    4.976     3.278   -1.448   11.400
all.equal(me$estimate[me$term == "Eth"], round(mean(diff_by_child), 3))
#> [1] TRUE
```

With the default `indicator = "pseudo_elasticity"` the same changes are
reported as percentages; for a model without other components the
pseudo-elasticity of a dummy variable is exactly the rate ratio minus
one.

### A factor with more than two levels

A derivative makes no sense for a factor, and a single “effect” of a
four-level factor needs a choice of contrast.
[`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
therefore refuses `Age` with an explanation:

``` r

tryCatch(margEffTable(nb, vars = "Age"),
         error = function(e) conditionMessage(e))
#> [1] "Variable `Age` is a factor with more than two levels. Create binary indicators or supply a specific contrast before calling `margEffTable()`."
```

Construct the contrasts explicitly instead: here, the expected absences
if every child in the data were in each age group, all else as observed:

``` r

by_age <- sapply(levels(q$Age), function(a) {
  d <- q
  d$Age <- factor(a, levels = levels(q$Age))
  mean(predict(nb, newdata = d))
})
round(by_age, 2)
#>    F0    F1    F2    F3 
#> 16.57 10.58 18.09 23.67
round(by_age - by_age[["F0"]], 2)   # differences from F0, in days
#>    F0    F1    F2    F3 
#>  0.00 -5.99  1.53  7.11
```

## An interaction: yarn breaks

`warpbreaks` records the number of breaks in 54 looms of yarn (each a
fixed length), for two wools and three tension levels, 9 looms per
combination. With an interaction, the effect of wool depends on the
tension, so a single wool effect is not meaningful:

``` r

w <- datasets::warpbreaks
wb <- countreg(breaks ~ wool * tension, data = w, family = "NB2",
               method = "BFGS")
grid <- expand.grid(wool = levels(w$wool), tension = levels(w$tension))
cbind(grid, round(predict(wb, newdata = grid, interval = "confidence"), 2))
#>   wool tension   fit   lwr   upr
#> 1    A       L 44.55 36.13 54.76
#> 2    B       L 28.22 22.56 35.13
#> 3    A       M 24.00 19.05 30.05
#> 4    B       M 28.78 22.83 35.61
#> 5    A       H 24.55 19.46 30.58
#> 6    B       H 18.78 14.74 23.69
```

``` r

ci <- predict(wb, newdata = grid, interval = "confidence")
x <- as.integer(grid$tension) + ifelse(grid$wool == "A", -0.08, 0.08)
plot(x, ci$fit, pch = ifelse(grid$wool == "A", 1, 19), xaxt = "n",
     ylim = range(ci$lwr, ci$upr), xlab = "Tension",
     ylab = "Expected breaks")
axis(1, at = 1:3, labels = levels(w$tension))
segments(x, ci$lwr, x, ci$upr)
legend("topright", c("wool A", "wool B"), pch = c(1, 19), bty = "n")
```

![Expected breaks by tension for each wool, with 95% confidence
intervals.](predictions-and-effects_files/figure-html/warp-plot-1.png)

Expected breaks by tension for each wool, with 95% confidence intervals.

The difference between wools within each tension, with delta-method
standard errors, comes from
[`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
applied to the rows of that tension (the model is the same; only the
population over which the change is averaged differs):

``` r

wool_eff <- do.call(rbind, lapply(levels(w$tension), function(tn) {
  r <- margEffTable(wb, data = w[w$tension == tn, ], vars = "wool",
                    indicator = "discrete_change", tableType = "tibble")
  data.frame(tension = tn, wool_B_minus_A = r$estimate,
             std_error = r$std_error)
}))
wool_eff
#>   tension wool_B_minus_A std_error
#> 1       L        -16.333     5.802
#> 2       M          4.777     4.333
#> 3       H         -5.778     3.688
clear <- abs(wool_eff$wool_B_minus_A) > 1.96 * wool_eff$std_error
```

The difference between the wools is clearly different from zero (more
than 1.96 standard errors) only at L tension, where wool B has fewer
breaks than wool A; at the other tension levels the estimates are small
relative to their standard errors.

## Continuous predictors and elasticities: doctor visits

[`AER::DoctorVisits`](https://rdrr.io/pkg/AER/man/DoctorVisits.html)
(Cameron and Trivedi, 1986) records the number of doctor visits in the
past two weeks for 5,190 adults in the 1977-78 Australian Health Survey.
`age` is age in years divided by 100 and `income` is annual income in
tens of thousands of dollars, so a change of 0.01 in `age` is one year.

``` r

e <- new.env()
utils::data("DoctorVisits", package = "AER", envir = e)
dv <- e$DoctorVisits
dvm <- countreg(visits ~ gender + age + I(age^2) + income + illness +
                  health + private, data = dv, family = "NB2",
                method = "BFGS")
```

Age enters through two terms, so its effect is not a single coefficient.
The average marginal effect differentiates the expected count with
respect to `age` (changing both terms) for each person and averages. A
central finite difference by hand gives the same value:

``` r

me_dv <- margEffTable(dvm, vars = c("age", "income"), tableType = "tibble")
me_dv[, c("term", "effect_metric", "estimate", "std_error")]
#> Marginal Effects / Elasticities
#> Method: delta 
#> 
#>    term           effect_metric estimate std_error
#>     age Average marginal effect    0.264     0.057
#>  income Average marginal effect   -0.040     0.034
h <- 1e-4
d_up <- d_dn <- dv
d_up$age <- dv$age + h
d_dn$age <- dv$age - h
mean((predict(dvm, newdata = d_up) - predict(dvm, newdata = d_dn)) / (2 * h))
#> [1] 0.2641004
```

The effect is per unit of `age`, i.e. per 100 years: per year of age the
expected number of visits changes by 0.0026 on average. An elasticity is
free of units:

``` r

el <- margEffTable(dvm, vars = "income", measure = "elasticity",
                   tableType = "tibble")
el[, c("term", "effect_metric", "estimate", "std_error", "lower_ci",
       "upper_ci")]
#> Marginal Effects / Elasticities
#> Method: delta 
#> 
#>    term effect_metric estimate std_error lower_ci upper_ci
#>  income    Elasticity   -0.076     0.064   -0.201    0.049
```

A 1% higher income is associated with an expected number of visits about
0.08% lower; the 95% confidence interval includes zero.

For a nonlinear term, a curve is clearer than a single number. Here the
expected visits are averaged over the 5,190 people with everyone’s age
set to each value in turn:

``` r

ages <- seq(0.2, 0.7, by = 0.05)
avg_visits <- sapply(ages, function(a) {
  d <- dv
  d$age <- a
  mean(predict(dvm, newdata = d))
})
plot(100 * ages, avg_visits, type = "b", xlab = "Age (years)",
     ylab = "Average expected visits")
```

![Average expected doctor visits in two weeks by
age.](predictions-and-effects_files/figure-html/dv-curve-1.png)

Average expected doctor visits in two weeks by age.

## Uncertainty of effects

- **Delta method** (`se = "delta"`, the default): uses the model’s
  covariance matrix, the one chosen with `stderr` when it was fit
  (Hessian-based, robust or bootstrap). It is fast and is valid when the
  effect is a smooth function of the estimates and the sample is not
  small.
- **Bootstrap** (`se = "bootstrap"`): refits the model to resampled
  units (rows, panels, or the clusters in `cluster_var`) and recomputes
  the effect. Use it with a seed and enough replicates; the number of
  successful replicates is reported.

``` r

set.seed(11)
margEffTable(nb, vars = "Eth", indicator = "discrete_change",
             se = "bootstrap", bootstraps = 100,
             tableType = "tibble")[, c("estimate", "std_error",
                                       "successful_bootstraps")]
#> Marginal Effects / Elasticities
#> Method: bootstrap 
#> 
#>  estimate std_error successful_bootstraps
#>    -9.246     2.891                   100
```

### Averaging over a different population

Effects can be averaged over a population other than the estimation
sample by passing it as `data`. With the delta method the model is the
one already estimated, and only the averaging population changes; here,
the effects among the slow learners only:

``` r

slow <- q[q$Lrn == "SL", ]
margEffTable(nb, data = slow, vars = c("Eth", "Sex"),
             indicator = "discrete_change",
             tableType = "tibble")[, c("term", "estimate", "std_error")]
#> Marginal Effects / Elasticities
#> Method: delta 
#> 
#>  term estimate std_error
#>   Eth   -9.296     2.796
#>   Sex    1.371     2.770
```

The bootstrap does not yet separate the population used to re-estimate
the model from the population over which the effect is averaged: with
`data` given, it resamples those rows as if they were the estimation
sample. Use delta-method standard errors for effects averaged over a
different population.

## Pitfalls

- **A derivative for a factor.**
  [`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
  uses discrete changes for two-level factors and refuses factors with
  more levels; construct explicit contrasts for those.
- **Coefficients of model components.** In a zero-inflated, hurdle or
  underreporting model, the count-model coefficients describe one
  component. The change in the *overall* expected count also involves
  the zero or reporting model, so compute effects from predictions (as
  [`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
  does), not from the coefficients alone.
- **Associations, not causal effects.** All effects here are comparisons
  of model-based expected counts between profiles. A causal
  interpretation needs a separate argument about how the data were
  generated.

## Next steps

- [`vignette("getting-started", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/getting-started.md)
  shows predictions with confidence and prediction intervals.
- [`vignette("inference-prediction-and-comparison", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/inference-prediction-and-comparison.md)
  covers the choice of standard errors.
- [`vignette("distribution-parameter-regression", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/distribution-parameter-regression.md)
  shows covariates that change the spread, not only the mean.

## References

Cameron, A. C. and Trivedi, P. K. (1986). Econometric models based on
count data: comparisons and applications of some estimators and tests.
*Journal of Applied Econometrics*, 1(1), 29-53.

Quine, S. (1975). *Achievement Orientation of Aboriginal and White
Adolescents*. Doctoral dissertation, Australian National University.
