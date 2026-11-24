# Total-plus-allocation models: counts by severity

Crash counts often come by severity: property damage only (PDO), injury
and fatal crashes at each site. One way to model them is to describe two
things:

- how many crashes happen (the **total**);
- how they **split** across the severity levels.

[`countreg_alloc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_alloc.md)
fits both parts jointly. The total follows any count family of the
package. The split follows an ordered probit, ordered logit or
multinomial logit model, and each part has its own covariates.

``` r

library(flexCountReg)
data("washington_roads")
w <- washington_roads
w$pdo <- w$Total_crashes - w$Injury_crashes - w$Fatal_crashes
w$injury <- w$Injury_crashes
w$fatal <- w$Fatal_crashes
colSums(w[, c("pdo", "injury", "fatal")])
#>    pdo injury  fatal 
#>    633     57      5
```

## 1. The model

For a site with category counts \\y_1, \dots, y_J\\ (from the least to
the most severe) and total \\c = \sum_j y_j\\:

\\ P(y_1, \dots, y_J) = P(c) \times \frac{c!}{\prod_j y_j!} \prod_j
\pi_j^{y_j}. \\

\\P(c)\\ is the count model of the total, here NB2 with mean \\\mu =
\exp(x'\beta)\\. The \\\pi_j\\ are the allocation probabilities. With an
ordered probit they come from a latent severity propensity \\z'\gamma\\
and thresholds \\\psi_k\\:

\\ \pi_j = \Phi(\psi_j - z'\gamma) - \Phi(\psi\_{j-1} - z'\gamma). \\

``` r

m <- countreg_alloc(cbind(pdo, injury, fatal) ~ lnaadt + speed50 +
                      offset(lnlength),
                    alloc_formula = ~ lnaadt + speed50,
                    data = w, family = "NB2")
summary(m)
#> Total-plus-allocation count model
#> 
#> Total:       NB2 
#> Allocation:  ordered_probit with multinomial likelihood
#> Categories:  pdo < injury < fatal 
#> Units:      1501 (1101 with a zero total)
#> Std. errors: normal 
#> 
#> Total:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -8.89586    0.44648 -19.924  < 2e-16 ***
#> lnaadt       1.12442    0.05141  21.870  < 2e-16 ***
#> speed50     -0.56771    0.11020  -5.152 2.58e-07 ***
#> 
#> Distribution:
#>           Estimate Std. Error z value Pr(>|z|)    
#> ln(alpha)  -0.9126     0.2279  -4.005 6.21e-05 ***
#> 
#> Allocation:
#>                    Estimate Std. Error z value Pr(>|z|)   
#> Allocation:lnaadt  -0.21819    0.06993  -3.120  0.00181 **
#> Allocation:speed50 -0.42783    0.19748  -2.166  0.03028 * 
#> 
#> Thresholds:
#>                        Estimate Std. Error z value Pr(>|z|)
#> Threshold1:(Intercept)  -0.5968     0.6035  -0.989    0.323
#> Threshold2:(Intercept)   0.1197     0.1410   0.849    0.396
#> 
#> Log-likelihood: -1275.7071  AIC: 2567.4143  BIC: 2609.9254
```

A positive allocation coefficient shifts crashes toward the more severe
categories.

The thresholds are \\\psi_1\\ (`Threshold1`) and the log of the gap
\\\psi_2 - \psi_1\\ (`Threshold2`). This keeps them in order.

## 2. Zero totals are part of the model

Most segment-years have no crashes. Such a segment-year contributes
\\P(c = 0)\\ to the likelihood and nothing about the split, since there
are no crashes to split. The total is **not truncated** at one: zeros
are ordinary observations of the count model.

Data sometimes hold only units with at least one crash by construction,
for example when they are built from crash records. There the total
should follow the zero-truncated distribution \\P(c) / (1 - P(0))\\,
which `truncated = TRUE` gives. Fitting such data without it biases the
count model, and
[`countreg_alloc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_alloc.md)
says so when the data have no zero totals:

``` r

pos <- w[w$Total_crashes > 0, ]
f <- cbind(pdo, injury, fatal) ~ lnaadt + speed50 + offset(lnlength)
mt <- countreg_alloc(f, alloc_formula = ~ lnaadt + speed50, data = pos,
                     family = "NB2", truncated = TRUE)
mn <- countreg_alloc(f, alloc_formula = ~ lnaadt + speed50, data = pos,
                     family = "NB2")
#> No unit has a zero total. If the data hold only units with at least one event by construction, use `truncated = TRUE`; otherwise the untruncated total is right.
se <- function(x) sqrt(diag(vcov(x)))[1:4]
round(rbind(all_rows = coef(m)[1:4], se = se(m),
            positive_truncated = coef(mt)[1:4], se = se(mt),
            positive_not_truncated = coef(mn)[1:4]), 3)
#>                        (Intercept) lnaadt speed50 ln(alpha)
#> all_rows                    -8.896  1.124  -0.568    -0.913
#> se                           0.446  0.051   0.110     0.228
#> positive_truncated         -10.945  1.342  -0.162    -0.946
#> se                           1.092  0.120   0.196     0.462
#> positive_not_truncated      -3.491  0.569  -0.072    -3.509
```

On the 400 segment-years with crashes, the truncated fit estimates the
same count model as the fit on all rows, but from much less information:
its standard errors are more than twice as large, because the zeros
carried information too, and its estimates differ from the full fit by
up to about two of those standard errors. The untruncated fit on the
same rows is badly biased: it mistakes the missing zeros for a much
higher crash rate and almost no overdispersion.

The expected totals of a truncated fit are those of the untruncated
distribution. `type = "total_positive"` (or `conditional = TRUE` for the
category counts) gives the expected count *given* at least one crash.

## 3. Predictions

[`predict()`](https://rdrr.io/r/stats/predict.html) returns the expected
count in each category, \\E\[c\]\\\pi_j\\, by default. `type = "prob"`
gives the probabilities \\\pi_j\\, and `type = "total"` gives the
expected total:

``` r

nd <- w[c(1, 50, 200), ]
predict(m, newdata = nd)
#>           pdo     injury        fatal
#> 1   0.7667519 0.02801566 0.0014201332
#> 50  0.1031854 0.01093792 0.0009719956
#> 200 0.7648239 0.04933585 0.0033316044
predict(m, newdata = nd, type = "prob")
#>           pdo     injury       fatal
#> 1   0.9630291 0.03518726 0.001783666
#> 50  0.8965213 0.09503359 0.008445140
#> 200 0.9355743 0.06035031 0.004075400
```

[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws new category
counts: a total from the count model, then its split from the
multinomial distribution.

## 4. Generalized thresholds and the multinomial logit

`threshold_formula` lets the thresholds depend on covariates, so a
variable can change the shares of the middle categories differently from
the extremes:

``` r

mg <- countreg_alloc(cbind(pdo, injury, fatal) ~ lnaadt + speed50 +
                       offset(lnlength), alloc_formula = ~ lnaadt,
                     threshold_formula = ~ speed50, data = w,
                     family = "NB2")
AIC(m, mg)
#>    df      AIC
#> m   8 2567.414
#> mg  9 2568.491
```

For categories without an order, such as crash types, use
`link = "mnl"`. Each category except the first then has its own
coefficients.

## 5. Multinomial or fractional likelihood

The default likelihood raises each \\\pi_j\\ to the category **count**
\\y_j\\. That makes it a proper likelihood of the category counts, so
likelihood-ratio tests, AIC and BIC are valid.

Shabab et al. (2026) instead raise \\\pi_j\\ to the observed **share**
\\y_j / c\\. Then a site-year with one crash carries the same weight in
the split as one with a hundred. `likelihood = "fractional"` gives that
version. It is a quasi-likelihood, so:

- standard errors come from the sandwich, clustered by `group_var` when
  it is given;
- AIC, BIC and likelihood-ratio tests are not available.

``` r

mf <- countreg_alloc(cbind(pdo, injury, fatal) ~ lnaadt + speed50 +
                       offset(lnlength), alloc_formula = ~ lnaadt + speed50,
                     data = w, family = "NB2", likelihood = "fractional",
                     group_var = "ID")
cbind(multinomial = coef(m), fractional = coef(mf))
#>                        multinomial  fractional
#> (Intercept)             -8.8958570 -8.89587118
#> lnaadt                   1.1244164  1.12441836
#> speed50                 -0.5677125 -0.56771887
#> ln(alpha)               -0.9125921 -0.91257179
#> Allocation:lnaadt       -0.2181929 -0.13393479
#> Allocation:speed50      -0.4278347 -0.35836783
#> Threshold1:(Intercept)  -0.5968112  0.14842063
#> Threshold2:(Intercept)   0.1196582  0.08127466
```

## 6. Effects that change over time

[`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
works in every formula of
[`countreg_alloc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_alloc.md),
including the allocation and threshold formulas. See
[`vignette("temporal-instability")`](https://jwood-iastate.github.io/flexCountReg/articles/temporal-instability.md).

## 7. Random parameters

Sites differ in ways the covariates do not capture: some have more
crashes than their traffic and speed limit predict, and some have more
severe crashes. `rpar_total` and `rpar_alloc` add random parameters to
the total and to the severity propensity. With `group_var = "ID"`, each
road segment gets one draw of its random parameters, shared by its three
years. `correlation = "across"` lets the two be correlated:

``` r

mr <- countreg_alloc(cbind(pdo, injury, fatal) ~ lnaadt + speed50 +
                       offset(lnlength), alloc_formula = ~ lnaadt + speed50,
                     data = w, family = "POISSON", group_var = "ID",
                     rpar_total = ~ 1, rpar_alloc = ~ 1,
                     correlation = "across", rpar_draws = 200, seed = 1)
summary(mr)
#> Total-plus-allocation count model
#> 
#> Total:       POISSON 
#> Allocation:  ordered_probit with multinomial likelihood
#> Categories:  pdo < injury < fatal 
#> Units:      507 (1501 rows, 1101 with a zero total)
#> Std. errors: normal 
#> 
#> Total:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -8.93059    0.50334 -17.743  < 2e-16 ***
#> lnaadt       1.10937    0.05881  18.864  < 2e-16 ***
#> speed50     -0.59187    0.13227  -4.475 7.65e-06 ***
#> 
#> Allocation:
#>                    Estimate Std. Error z value Pr(>|z|)   
#> Allocation:lnaadt  -0.23373    0.07423  -3.149  0.00164 **
#> Allocation:speed50 -0.41417    0.20392  -2.031  0.04225 * 
#> 
#> Thresholds:
#>                        Estimate Std. Error z value Pr(>|z|)
#> Threshold1:(Intercept)  -0.7951     0.6437  -1.235    0.217
#> Threshold2:(Intercept)   0.1443     0.1442   1.001    0.317
#> 
#> Random parameters (Cholesky factor):
#>                                               Estimate Std. Error z value
#> sd:Total:(Intercept)                           0.64050    0.06507   9.843
#> chol:Allocation:(Intercept)|Total:(Intercept) -0.19673    0.10527  -1.869
#> sd:Allocation:(Intercept)                     -0.11052    0.27038  -0.409
#>                                               Pr(>|z|)    
#> sd:Total:(Intercept)                            <2e-16 ***
#> chol:Allocation:(Intercept)|Total:(Intercept)   0.0617 .  
#> sd:Allocation:(Intercept)                       0.6827    
#> 
#> Random parameters: standard deviations and correlations
#>                           SD Total:(Intercept)
#> Total:(Intercept)      0.640                  
#> Allocation:(Intercept) 0.226            -0.872
#> (200 draws per unit of ID; 507 units)
#> 
#> Log-likelihood: -1252.9034  AIC: 2525.8068  BIC: 2578.9457
AIC(m, mr)
#>    df      AIC
#> m   8 2567.414
#> mr 10 2525.807
```

The random intercepts fit much better than the NB2 model without them.
The segment effect on the total (a Poisson-lognormal total) takes over
the role of the NB2 overdispersion. The negative correlation says that
segments with more crashes than their covariates predict tend to have
relatively less severe crashes.

The likelihood is simulated: each segment’s likelihood is averaged over
`rpar_draws` draws of its random parameters (randomly shifted Halton
points), with analytic gradients. `rpar_threshold` adds a random term to
the gaps between thresholds, which changes the share of the middle
categories. Random slopes are more terms in the same formulas, such as
`rpar_total = ~ 1 + lnaadt`.

Predictions for a new segment average over the random parameters
(`random = "marginal"`, the default). For a segment in the estimation
data, `random = "posterior"` uses its own crash history (Empirical
Bayes):

``` r

nd <- w[w$ID == w$ID[1], ]
round(cbind(predict(mr, newdata = nd),
            predict(mr, newdata = nd, random = "posterior")), 3)
#>        pdo injury fatal   pdo injury fatal
#> 1    0.754  0.028 0.001 0.535  0.025 0.001
#> 502  0.750  0.028 0.001 0.532  0.025 0.001
#> 1002 0.791  0.029 0.001 0.561  0.026 0.001
```

## Reference

Shabab, K. R., Bhowmik, T., Zaki, M. H. and Eluru, N. (2026). A joint
count and generalized ordered fractional split approach for addressing
temporal instability in road safety analysis. *Analytic Methods in
Accident Research*, 51, 100440.
