# Total-plus-allocation count models

Fits a joint model of a total count and its allocation across
categories, such as the number of crashes at a site and their split by
severity (ordered categories) or by crash type (unordered categories).
The total follows any count family of the package (with covariates in
its mean and distribution parameters). The split follows an ordered
probit or ordered logit model for ordered categories, with covariates in
the propensity and in the thresholds (generalized thresholds), or a
multinomial logit model (`link = "mnl"`) for unordered categories, with
covariates whose effects differ by category.

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
  verbose = FALSE,
  rpar_total = NULL,
  rpar_alloc = NULL,
  rpar_threshold = NULL,
  correlation = c("none", "within", "across"),
  rpar_draws = 500,
  seed = NULL,
  shared_formula = NULL,
  shared_in = c("total", "alloc"),
  shared_nodes = 20,
  .setup_only = FALSE
)
```

## Arguments

- formula:

  a formula with the category counts on the left, as
  `cbind(cat_1, ..., cat_J)`, and the total's mean model on the right.
  For the ordered links the columns go from the lowest to the highest
  category; for `"mnl"` the first column is the base category and the
  order of the others does not matter. The total is the sum of the
  categories. An [`offset()`](https://rdrr.io/r/stats/offset.html) term
  (an exposure) belongs here.

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

  `"ordered_probit"` (default) or `"ordered_logit"` for ordered
  categories such as severity levels; `"mnl"` (multinomial logit) for
  unordered categories such as crash types, collision types or vehicle
  types. See Details.

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

- rpar_total, rpar_alloc, rpar_threshold:

  one-sided formulas of random parameters in the total's mean, the
  allocation propensity and the gaps between thresholds. Each term adds
  a normal random deviation, with mean zero, to the corresponding
  coefficient, so the coefficient in the equation's own formula is the
  mean (`~ 1` in `rpar_alloc` is a random intercept of the propensity,
  even though its fixed intercept is absorbed by the thresholds).
  `rpar_threshold` shifts the log gaps of all thresholds above the first
  by one common random amount, so it changes the share of the middle
  categories (a common shift of all thresholds would be the same as a
  random propensity intercept). For `link = "mnl"`, `rpar_alloc` gives
  each category but the first its own random coefficients. The random
  coefficients are drawn once per unit of `group_var` (shared by all its
  rows, such as the years of a site), or once per row without
  `group_var`.

- correlation:

  correlation of the random parameters: `"none"` (independent),
  `"within"` (correlated within each equation) or `"across"` (correlated
  within and across the total and the allocation, which lets a unit's
  unobserved crash frequency and its severity propensity be correlated).

- rpar_draws:

  the number of simulation draws per unit (randomly shifted Halton
  points, a separate block and shift for each unit).

- seed:

  seed of the random shifts of the draws (stored in the fit, so
  predictions use the same draws); random by default.

- shared_formula:

  optional one-sided formula of a shared site-year (row-level) random
  term: one standard normal draw \\u\\ per row that enters each equation
  in `shared_in` with a loading \\s'\kappa_e\\, where \\s\\ are the
  variables of `shared_formula` (`~ 1` for a constant loading) and
  \\\kappa_e\\ its coefficients in equation \\e\\. It correlates a row's
  total with its split (the term \\\eta\_{it}\\ of Shabab et al. 2026);
  see Details.

- shared_in:

  the equations the shared term enters: any of `"total"` (the log mean),
  `"alloc"` (the propensity; for `"mnl"`, each category but the first)
  and `"threshold"` (the log gaps between thresholds). Default
  `c("total", "alloc")`.

- shared_nodes:

  the number of Gauss-Hermite nodes for the integral over the shared
  term of each row. After the fit, the log-likelihood is recomputed with
  twice the nodes, with a warning if it changes by more than 0.05 (large
  loadings, above about 1 in the total, can need 40 or more nodes).

- .setup_only:

  internal use: return the model setup without fitting.

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

**Unordered categories.** With `link = "mnl"` the categories are types
with no order (crash types, collision types, vehicle types) and \\\pi_j
= \exp(z'\delta_j) / \sum_k \exp(z'\delta_k)\\, with \\\delta_1 = 0\\
for the first category (the base). Each other category has its own
coefficients, named `Allocation[<category>]:`, the effects on the log
odds of that category against the base; a different base gives the same
fit. With a Poisson total and the count model's covariates in
`alloc_formula`, the model is the same as independent Poisson
regressions of the categories; with an overdispersed total (such as NB2)
the categories share the total's unobserved heterogeneity. The split
assumes independence of irrelevant alternatives; random coefficients in
`rpar_alloc` (one set for each category but the first) relax it. The
threshold options (`threshold_formula`, `rpar_threshold`,
`shared_in = "threshold"`) are for the ordered links only.

**Thresholds.** For the ordered links the thresholds are \\\psi_1 =
v'\theta_1\\ and \\\psi_k = \psi\_{k-1} + \exp(v'\theta_k)\\, which
keeps them ordered for every \\v\\. The coefficients `Threshold1:` are
on the scale of the latent propensity and `Threshold<k>:` (k \>= 2) on
the log scale of the gap between thresholds.

**Random parameters.** With `rpar_total`, `rpar_alloc` or
`rpar_threshold`, the likelihood of a unit is the average over
`rpar_draws` draws of its coefficients of the product of its rows'
likelihoods (maximum simulated likelihood), with analytic gradients. The
random coefficients are \\b = L u\\, with \\u\\ standard normal and
\\L\\ lower triangular; its elements are reported as `sd:` (the
diagonal) and `chol:` (below the diagonal), and
[`summary()`](https://rdrr.io/r/base/summary.html) shows the implied
standard deviations and correlations. With the fractional likelihood the
integrand is the product of the rows' quasi-likelihoods, as in Shabab et
al. (2026).

**Shared site-year term.** Random parameters drawn per unit of
`group_var` are shared by all the rows (years) of the unit. A row's own
total and split can also share an unobserved factor: `shared_formula`
adds \\\lambda\_{e,it} u\_{it}\\ to each equation \\e\\ in `shared_in`,
with \\u\_{it}\\ standard normal and drawn anew for each row, and
loading \\\lambda\_{e,it} = s\_{it}'\kappa_e\\. The covariance of the
term in two equations is \\\lambda\_{e}\lambda\_{e'}\\: equal signs of
the loadings mean that a site-year with more events than expected also
has a more severe split (with the ordered links), opposite signs a less
severe one. Changing the signs of all loadings gives the same model, so
only their relative signs are identified. Shabab et al. (2026) use one
loading for all equations, \\\eta\_{it} = k's\_{it}\\, so their
correlation is always positive; separate loadings include that model and
allow a negative one. The integral over \\u\_{it}\\ is computed by
Gauss-Hermite quadrature within each draw of the unit's random
parameters, adapted to each row (the nodes are centred and scaled at the
row's posterior mean and standard deviation of \\u\_{it}\\). A shared
term in the total also adds overdispersion to it, so with a family that
has its own dispersion parameter (such as NB2) the two compete and can
be weakly identified. With `truncated = TRUE` the row's likelihood is
divided by its marginal probability of a positive total, \\\int (1 -
P_0(u))\\\phi(u)\\du\\, as the selection on a positive total happens
after \\u\_{it}\\ is drawn.

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
# Ordered categories: crashes by severity
w$pdo <- w$Total_crashes - w$Injury_crashes - w$Fatal_crashes
m <- countreg_alloc(cbind(pdo, Injury_crashes, Fatal_crashes) ~ lnaadt +
                      offset(lnlength), alloc_formula = ~ speed50,
                    data = w, family = "NB2")
summary(m)
#> Total-plus-allocation count model
#> 
#> Total:       NB2 
#> Allocation:  ordered_probit with multinomial likelihood
#> Categories:  pdo < Injury_crashes < Fatal_crashes 
#> Units:      1501 (1101 with a zero total)
#> Std. errors: normal 
#> 
#> Total:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -9.38253    0.45195  -20.76   <2e-16 ***
#> lnaadt       1.16464    0.05252   22.18   <2e-16 ***
#> 
#> Distribution:
#>           Estimate Std. Error z value Pr(>|z|)    
#> ln(alpha)  -0.7771     0.2133  -3.644 0.000269 ***
#> 
#> Allocation:
#>                    Estimate Std. Error z value Pr(>|z|)  
#> Allocation:speed50  -0.3610     0.1932  -1.868   0.0617 .
#> 
#> Thresholds:
#>                        Estimate Std. Error z value Pr(>|z|)    
#> Threshold1:(Intercept)  1.28728    0.07255  17.743   <2e-16 ***
#> Threshold2:(Intercept)  0.10665    0.14192   0.751    0.452    
#> 
#> Log-likelihood: -1294.2632  AIC: 2600.5263  BIC: 2632.4097 
head(predict(m))
#>         pdo Injury_crashes Fatal_crashes
#> 1 1.1768124     0.05791328   0.003570153
#> 2 1.0399738     0.05117918   0.003155019
#> 3 1.7241671     0.08484969   0.005230689
#> 4 0.3831482     0.01885549   0.001162375
#> 5 1.2315479     0.06060692   0.003736206
#> 6 0.5747224     0.02828323   0.001743563
# Unordered categories: crashes by type, other crashes as the base
w$other <- w$Total_crashes - w$Animal - w$Rollover
mt <- countreg_alloc(cbind(other, Animal, Rollover) ~ lnaadt +
                       offset(lnlength), alloc_formula = ~ lnaadt,
                     link = "mnl", data = w, family = "NB2")
summary(mt)
#> Total-plus-allocation count model
#> 
#> Total:       NB2 
#> Allocation:  mnl with multinomial likelihood
#> Categories:  other, Animal, Rollover (unordered; base: other )
#> Units:      1501 (1101 with a zero total)
#> Std. errors: normal 
#> 
#> Total:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -9.38260    0.45195  -20.76   <2e-16 ***
#> lnaadt       1.16465    0.05252   22.18   <2e-16 ***
#> 
#> Distribution:
#>           Estimate Std. Error z value Pr(>|z|)    
#> ln(alpha)  -0.7772     0.2133  -3.644 0.000269 ***
#> 
#> Allocation:
#>                                   Estimate Std. Error z value Pr(>|z|)   
#> Allocation[Animal]:(Intercept)   -0.003171   1.046819  -0.003  0.99758   
#> Allocation[Animal]:lnaadt        -0.222864   0.121159  -1.839  0.06585 . 
#> Allocation[Rollover]:(Intercept)  1.607992   1.647893   0.976  0.32917   
#> Allocation[Rollover]:lnaadt      -0.571248   0.198816  -2.873  0.00406 **
#> 
#> Log-likelihood: -1408.7772  AIC: 2831.5545  BIC: 2868.7517 
predict(mt, newdata = w[c(1, 50, 200), ], type = "prob")
#>          other    Animal   Rollover
#> [1,] 0.8583601 0.1160507 0.02558911
#> [2,] 0.7413005 0.1713395 0.08736004
#> [3,] 0.8814390 0.1012544 0.01730666
# }
```
