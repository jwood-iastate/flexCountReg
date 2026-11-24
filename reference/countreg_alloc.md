# Total-plus-allocation count models

Fits a joint model of a total count and its allocation across
categories, such as the number of crashes at a site and their split by
severity. The total follows any count family of the package (with
covariates in its mean and distribution parameters), and the split
follows an ordered probit, ordered logit or multinomial logit model,
with covariates in the propensity and, for the ordered links, in the
thresholds (generalized thresholds).

## Usage

``` r
countreg_alloc(
  formula,
  data,
  family = "NB2",
  alloc_formula = NULL,
  threshold_formula = ~1,
  link = c("ordered_probit", "ordered_logit", "mnl"),
  likelihood = c("multinomial", "fractional"),
  dis_param_formula_1 = NULL,
  dis_param_formula_2 = NULL,
  offset = NULL,
  weights = NULL,
  truncated = FALSE,
  group_var = NULL,
  stderr = NULL,
  method = "BFGS",
  max.iters = 1000,
  start.vals = NULL,
  ndraws = 1500,
  engine = "quadrature",
  verbose = FALSE
)
```

## Arguments

- formula:

  a formula with the category counts on the left, as
  `cbind(cat_1, ..., cat_J)` from the lowest to the highest category
  (for the ordered links), and the total's mean model on the right. The
  total is the sum of the categories. An
  [`offset()`](https://rdrr.io/r/stats/offset.html) term (an exposure)
  belongs here.

- data:

  a data frame.

- family:

  the count family of the total (see
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)),
  or a family made with
  [`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md).

- alloc_formula:

  a one-sided formula for the allocation propensity. By default the
  right-hand side of `formula` without its offsets. For the ordered
  links its intercept is dropped (the thresholds absorb it).

- threshold_formula:

  a one-sided formula for the thresholds of the ordered links
  (generalized thresholds). Default `~ 1`, fixed thresholds. It must
  have an intercept.

- link:

  `"ordered_probit"` (default), `"ordered_logit"` or `"mnl"`
  (multinomial logit, for unordered categories; the first category is
  the base).

- likelihood:

  `"multinomial"` (default) or `"fractional"`; see Details.

- dis_param_formula_1, dis_param_formula_2:

  formulas for the total's distribution parameters, as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md).

- offset:

  optional names of offset columns (log scale), instead of
  [`offset()`](https://rdrr.io/r/stats/offset.html) in `formula`.

- weights:

  optional frequency weights: a column name or a numeric vector.

- truncated:

  `TRUE` when the data contain only units with a positive total (for
  example data built from crash records): the total then follows the
  zero-truncated distribution \\P(c) / (1 - P(0))\\.

- group_var:

  optional column naming the units for cluster-robust standard errors
  (such as a site observed over several years).

- stderr:

  `"normal"` (from the Hessian), `"robust"` (sandwich) or `"cluster"`
  (sandwich by `group_var`). The default is `"normal"` for the
  multinomial likelihood and `"cluster"` (or `"robust"` without
  `group_var`) for the fractional one, whose Hessian alone does not give
  valid standard errors.

- method, max.iters, start.vals, verbose:

  optimizer settings, as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md).

- ndraws, engine:

  integration settings of families without a closed form, as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md).

## Value

An object of class `"countreg_alloc"`, with methods
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html).

## Details

**The two parts.** With \\c = \sum_j y_j\\ the total, the multinomial
likelihood of a unit is \$\$P(c) \frac{c!}{\prod_j y_j!} \prod_j
\pi_j^{y_j},\$\$ the probability of the total times the multinomial
probability of its split. This is a proper joint likelihood of the
category counts, so likelihood-ratio tests, AIC and BIC are valid, and
so are simulation and prediction intervals. With a Poisson total the
categories are independent Poisson counts with means \\\mu \pi_j\\; with
an NB2 total they are negative multinomial.

**Zero totals are part of the model.** The total's distribution includes
zero. A unit with no events contributes \\P(c = 0)\\ and nothing about
the split, because a zero total has no shares. The total is *not*
truncated at one: the indicator \\\omega\\ of Shabab et al. (2026),
which switches off the allocation term for zero totals, does the same
thing. Only when the data hold units with at least one event by
construction should `truncated = TRUE` be used; fitting such data with
the untruncated total biases the count model, and `countreg_alloc()`
says so when it sees no zero totals.

**Fractional likelihood.** `likelihood = "fractional"` raises \\\pi_j\\
to the observed share \\y_j / c\\ instead of the count, as in Shabab et
al. (2026). A unit with one event then carries the same weight in the
split as a unit with a thousand. This is a quasi-likelihood: the
estimates are consistent, but standard errors must come from the
sandwich (the default), and likelihood-ratio tests, AIC and BIC are not
available.

**Thresholds.** For the ordered links the thresholds are \\\psi_1 =
v'\theta_1\\ and \\\psi_k = \psi\_{k-1} + \exp(v'\theta_k)\\, which
keeps them ordered for every \\v\\. The coefficients `Threshold1:` are
on the scale of the latent propensity and `Threshold<k>:` (k \>= 2) on
the log scale of the gap between thresholds.

Time-varying effects
([`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md))
can be used in every formula.

## References

Shabab, K. R., Bhowmik, T., Zaki, M. H. and Eluru, N. (2026). A joint
count and generalized ordered fractional split approach for addressing
temporal instability in road safety analysis. *Analytic Methods in
Accident Research*, 51, 100440.

## See also

[`predict.countreg_alloc()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.countreg_alloc.md),
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
[`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md)

## Examples

``` r
# \donttest{
data("washington_roads")
w <- washington_roads
# Illustrative split of the total crashes into three categories
set.seed(1)
sp <- t(vapply(w$Total_crashes, function(n)
  stats::rmultinom(1, n, c(0.6, 0.3, 0.1)), numeric(3)))
w$pdo <- sp[, 1]; w$minor <- sp[, 2]; w$severe <- sp[, 3]
m <- countreg_alloc(cbind(pdo, minor, severe) ~ lnaadt + offset(lnlength),
                    alloc_formula = ~ speed50, data = w, family = "NB2")
summary(m)
#> Total-plus-allocation count model
#> 
#> Total:       NB2 
#> Allocation:  ordered_probit with multinomial likelihood
#> Categories:  pdo < minor < severe 
#> Units:      1501 (1101 with a zero total)
#> Std. errors: normal 
#> 
#> Total:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -9.38253    0.45195  -20.76   <2e-16 ***
#> lnaadt       1.16465    0.05252   22.18   <2e-16 ***
#> 
#> Distribution:
#>           Estimate Std. Error z value Pr(>|z|)    
#> ln(alpha)  -0.7771     0.2133  -3.644 0.000269 ***
#> 
#> Allocation:
#>                    Estimate Std. Error z value Pr(>|z|)
#> Allocation:speed50   0.1289     0.1129   1.141    0.254
#> 
#> Thresholds:
#>                        Estimate Std. Error z value Pr(>|z|)    
#> Threshold1:(Intercept)  0.28603    0.05307   5.390 7.05e-08 ***
#> Threshold2:(Intercept)  0.09164    0.06049   1.515     0.13    
#> 
#> Log-likelihood: -1585.6988  AIC: 3183.3975  BIC: 3215.2809 
head(predict(m))
#>         pdo     minor     severe
#> 1 0.6964770 0.4117047 0.13011486
#> 2 0.6154913 0.3638321 0.11498522
#> 3 1.0204198 0.6031953 0.19063340
#> 4 0.2267600 0.1340434 0.04236298
#> 5 0.7288713 0.4308538 0.13616671
#> 6 0.3401399 0.2010651 0.06354447
# }
```
