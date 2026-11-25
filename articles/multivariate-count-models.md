# Multivariate count models

Crash data often hold several counts for each site: crashes of each
severity, or of each type.
[`countreg_mv()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mv.md)
fits one count regression for each outcome in one model. Each outcome
has its own family, mean formula and distribution-parameter formulas.

This is the direct way to model counts by category. A total and its
split across categories is the other way; see
[`vignette("total-allocation-models")`](https://jwood-iastate.github.io/flexCountReg/articles/total-allocation-models.md).

``` r

library(flexCountReg)
data("washington_roads")
w <- washington_roads
w$pdo <- w$Total_crashes - w$Injury_crashes
w$injury <- w$Injury_crashes
colSums(w[, c("pdo", "injury")])
#>    pdo injury 
#>    638     57
```

## 1. One family for each outcome

Outcomes can behave differently. Property damage only (PDO) crashes are
frequent; injury crashes are rarer and can be more skewed. So each
outcome can have its own family:

``` r

m <- countreg_mv(cbind(pdo, injury) ~ lnaadt + speed50 + offset(lnlength),
                 data = w, family = c(pdo = "NB2", injury = "PIG"))
summary(m)
#> Multivariate count model
#> 
#> Rows:        1501 
#> Std. errors: normal 
#> Equations:   independent given the covariates
#> 
#> pdo (NB2):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -9.42403    0.47838 -19.700  < 2e-16 ***
#> lnaadt       1.17438    0.05483  21.420  < 2e-16 ***
#> speed50     -0.51500    0.11414  -4.512 6.42e-06 ***
#> ln(alpha)   -0.86569    0.23122  -3.744 0.000181 ***
#> 
#> injury (PIG):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -7.4115     1.1903  -6.227 4.77e-10 ***
#> lnaadt        0.6661     0.1418   4.698 2.63e-06 ***
#> speed50      -1.2201     0.4146  -2.943  0.00325 ** 
#> ln(eta)       0.7900     1.9900   0.397  0.69137    
#> 
#> Log-likelihood: -1247.0544  AIC: 2510.1088  BIC: 2552.6199
```

The coefficients are named `<outcome>:<term>`. `family` takes one code
for every outcome, or one per outcome. Different right-hand sides are
given as a list of formulas, one per outcome:

``` r

m2 <- countreg_mv(list(pdo = pdo ~ lnaadt + speed50 + offset(lnlength),
                       injury = injury ~ lnaadt + offset(lnlength)),
                  data = w, family = c(pdo = "NB2", injury = "PIG"))
AIC(m, m2)
#>    df      AIC
#> m   8 2510.109
#> m2  7 2519.182
```

`dis_param_formula_1` and `dis_param_formula_2` take one formula for
every outcome or a list named by outcome, so the dispersion of one
outcome can depend on covariates while the others stay constant.

## 2. Choosing the families

The equations are fitted together, but in this version they are
independent given the covariates: the log-likelihood is the sum of the
equations’ log-likelihoods, and each equation’s estimates equal those of
a separate
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
fit. So the family of each outcome can be chosen by comparing fits that
change only that outcome’s family:

``` r

fams <- list(c("NB2", "NB2"), c("NB2", "PIG"), c("NB2", "PLN"),
             c("PIG", "PIG"))
sapply(fams, function(f) {
  fit <- countreg_mv(cbind(pdo, injury) ~ lnaadt + speed50 +
                       offset(lnlength), data = w, family = f)
  c(pdo = f[1], injury = f[2], AIC = round(AIC(fit), 1))
})
#>        [,1]     [,2]     [,3]   [,4]    
#> pdo    "NB2"    "NB2"    "NB2"  "PIG"   
#> injury "NB2"    "PIG"    "PLN"  "PIG"   
#> AIC    "2510.2" "2510.1" "2510" "2509.9"
```

Here the choices fit almost equally well (the AICs are within 0.3), so
the simpler NB2 is a reasonable choice for both outcomes. With rarer or
more skewed outcomes the differences between families are larger.

## 3. Tests across equations

One fit gives the covariance of all the estimates, so hypotheses that
involve several equations can be tested. Is the effect of traffic the
same for PDO and injury crashes?

``` r

b <- coef(m)
V <- vcov(m)
a <- c("pdo:lnaadt", "injury:lnaadt")
diff <- b[a[1]] - b[a[2]]
se <- sqrt(V[a[1], a[1]] + V[a[2], a[2]] - 2 * V[a[1], a[2]])
c(difference = unname(diff), z = unname(diff / se),
  p = 2 * pnorm(-abs(unname(diff / se))))
#>   difference            z            p 
#> 0.5083084901 3.3436907879 0.0008267183
```

Traffic increases PDO crashes clearly more than injury crashes (p \<
0.001), so the share of injury crashes falls as traffic grows. The
covariance between the two estimates is zero here, because the equations
are independent; with correlated equations it would not be.

## 4. Predictions

[`predict()`](https://rdrr.io/r/stats/predict.html) gives a matrix of
expected counts, one column per outcome. `type = "prob"` gives, for each
outcome, the probabilities of the counts 0, 1, 2, …, and
[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws all outcomes
at once:

``` r

nd <- w[c(1, 50, 200), ]
predict(m, newdata = nd)
#>             pdo     injury
#> [1,] 0.77454644 0.03005740
#> [2,] 0.09928381 0.01309017
#> [3,] 0.78249900 0.04238618
lapply(predict(m, newdata = nd, type = "prob", max_count = 3), round, 3)
#> $pdo
#>          0     1     2     3
#> [1,] 0.511 0.299 0.124 0.044
#> [2,] 0.907 0.086 0.006 0.000
#> [3,] 0.508 0.299 0.125 0.045
#> 
#> $injury
#>          0     1     2 3
#> [1,] 0.971 0.028 0.001 0
#> [2,] 0.987 0.013 0.000 0
#> [3,] 0.960 0.038 0.002 0
```

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) has an
`outcome` column, and
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`confint()`](https://rdrr.io/r/stats/confint.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) work as
for other models.

## 5. Correlated outcomes

The counts of different severity levels at a site are usually
correlated: unobserved conditions that raise PDO crashes also raise
injury crashes. This version treats the equations as independent given
the covariates, so that correlation is not modeled. Multivariate models
with shared random effects (correlated random parameters across
equations) and copula models are planned.
