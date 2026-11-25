# Hierarchical (mixed-effects) count regression

Fits count models with normally distributed random effects at one or
more grouping levels, nested or crossed, written as in lme4:
`(1 | site)` for a random intercept, `(1 + lnaadt | county)` for
correlated random intercepts and slopes, `(1 | county/segment)` for
segments nested in counties and `(1 | site) + (1 | year)` for crossed
effects. Any count family of
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
can be used, including user-defined families, and the distribution
parameters can depend on covariates.

## Usage

``` r
countreg_hier(
  formula,
  data,
  family = "NB2",
  dis_param_formula_1 = NULL,
  dis_param_formula_2 = NULL,
  offset = NULL,
  weights = NULL,
  engine = c("auto", "aghq", "laplace", "msl"),
  nAGQ = NULL,
  re_dist = "normal",
  re_draws = 500,
  cre = NULL,
  cre_group = NULL,
  method = "BFGS",
  max.iters = 1000,
  start.vals = NULL,
  ndraws = 1500,
  seed = 123,
  verbose = FALSE,
  keep_data = TRUE
)
```

## Arguments

- formula:

  a model formula with random-effect terms. The fixed part is as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  (including [`offset()`](https://rdrr.io/r/stats/offset.html)), and can
  have penalized smooth terms `s()` and `t2()` of mgcv. Group-level
  predictors are ordinary columns (constant within groups), and
  cross-level interactions are ordinary interactions.

- data:

  a data frame.

- family:

  a count family: a name from
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)'s
  list, or a family made with
  [`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md).

- dis_param_formula_1, dis_param_formula_2:

  optional one-sided formulas for the family's distribution parameters.

- offset:

  optional names of offset columns (or use
  [`offset()`](https://rdrr.io/r/stats/offset.html) in the formula).

- weights:

  optional frequency weights of the rows (a column name or a numeric
  vector).

- engine:

  the approximation of the integral over the random effects: `"aghq"`
  (adaptive Gauss-Hermite quadrature; one grouping term only),
  `"laplace"` (any number of terms, nested or crossed, and smooth
  terms), `"msl"` (maximum simulated likelihood; one grouping term, any
  `re_dist`), or `"auto"` (simulation for non-normal effects, quadrature
  for a single grouping term with at most two effects, otherwise
  Laplace).

- nAGQ:

  number of quadrature points per effect for `"aghq"` (default 9 for one
  effect, 5 for two or three). `nAGQ = 1` is the Laplace approximation.

- re_dist:

  distribution of the random effects of the grouping term: `"normal"`
  (default), `"uniform"`, `"triangular"` (symmetric) or `"gamma"`; one
  value, or one per effect. Effects other than normal need a single
  grouping term and use `engine = "msl"`. With `"gamma"` (one effect)
  the multiplicative group effect \\\exp(b_j)\\ is gamma with mean 1 and
  standard deviation `sd_...`, as in the negative binomial
  random-effects model; with the others \\b_j\\ has mean 0 and the
  reported standard deviation.

- re_draws:

  number of (randomly shifted) Halton draws per group for
  `engine = "msl"`.

- cre:

  optional one-sided formula of time-varying regressors whose group
  means are added to the mean model (correlated random effects; see
  [`cre_test()`](https://jwood-iastate.github.io/flexCountReg/reference/cre_test.md)).

- cre_group:

  the grouping column(s) for the means of `cre`; by default those of the
  random-effect terms, when all terms share them.

- method:

  optimizer for [`optim()`](https://rdrr.io/r/stats/optim.html):
  `"BFGS"` (default) or `"Nelder-Mead"`.

- max.iters:

  maximum number of optimizer iterations.

- start.vals:

  optional starting values for the fixed coefficients and distribution
  parameters (as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)),
  optionally followed by the random-effect standard deviations and
  correlations.

- ndraws:

  number of Halton draws for families whose density is simulated (as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)).

- seed:

  seed of the random shifts of the Halton draws of `engine = "msl"`.

- verbose:

  print the optimizer's progress.

- keep_data:

  keep the training data in the fit (needed for predictions without
  `newdata`).

## Value

An object of class `c("countreg_hier", "flexCountReg")`. Methods:
[`summary()`](https://rdrr.io/r/base/summary.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html),
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html),
[`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html) and
[`re_lrtest()`](https://jwood-iastate.github.io/flexCountReg/reference/re_lrtest.md).

## Details

The random effects of grouping term \\l\\ are \\b\_{lj} \sim N(0,
\Sigma_l)\\ for each of its groups \\j\\, and enter the log of the mean:
\$\$\log \mu_i = x_i'\beta + o_i + \sum_l z\_{li}' b\_{l, g_l(i)}.\$\$
Each \\\Sigma_l = L_l L_l'\\ is estimated through its Cholesky factor
with an unconstrained diagonal, so a variance can be estimated as zero
(a singular fit, reported by
[`summary()`](https://rdrr.io/r/base/summary.html)).

The marginal likelihood integrates the random effects out. With
`engine = "aghq"` each group's integral uses Gauss-Hermite nodes
centered at the group's mode and scaled by its curvature; this is
accurate for small groups and is the default for one grouping factor.
`engine = "laplace"` uses the Laplace approximation, as
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) (with
`nAGQ = 1`) and `glmmTMB` do, which handles crossed and nested terms
with sparse matrices. The derivatives of the log-density are analytic
for the Poisson, NB1 and NB2 families and numerical for the others.

Smooth terms are fitted in their mixed-model form: the penalized part of
each smooth is a block of independent normal coefficients with variance
\\1/\lambda\\, so the smoothing parameters are estimated by the
(Laplace) marginal likelihood, as in gamm4. Their coefficients are
reported on mgcv's basis, with the covariance of the coefficients and
effects given the variance parameters, so
[`plot_smooths()`](https://jwood-iastate.github.io/flexCountReg/reference/plot_smooths.md)
and the other smooth tools work as for
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md).

With `engine = "msl"` the group likelihoods are averages over `re_draws`
draws of each group's effects, so effects can have any of the `re_dist`
distributions; the effects of the groups
([`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html)) are then
posterior means.

Standard errors come from the numerical Hessian of the marginal
log-likelihood. The variance components are reported as standard
deviations and correlations (`sd_(Intercept)|site`); their Wald standard
errors are rough near zero, where the likelihood-ratio test of
[`re_lrtest()`](https://jwood-iastate.github.io/flexCountReg/reference/re_lrtest.md)
is preferable.

[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
fits random parameters with a free mean and non-normal distributions at
one level;
[`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)
and
[`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md)
fit conjugate random effects in closed form.

## See also

[`predict.countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.countreg_hier.md),
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html)

## Examples

``` r
# \donttest{
data("washington_roads")
h <- countreg_hier(Total_crashes ~ lnaadt + speed50 + offset(lnlength) +
                     (1 | ID), data = washington_roads, family = "NB2")
summary(h)
#> Hierarchical count model: NB2
#> Integration: adaptive Gauss-Hermite quadrature, 9 points
#> Call: countreg_hier(formula = Total_crashes ~ lnaadt + speed50 + offset(lnlength) + (1 | ID), data = washington_roads, family = "NB2")
#> 
#> Random effects:
#>  Groups Name        Variance Std.Dev.
#>  ID     (Intercept) 0.408    0.6388  
#> Number of obs: 1501; groups: ID, 507
#> 
#> Fixed coefficients and distribution parameters:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -9.00121    0.50404 -17.858  < 2e-16 ***
#> lnaadt       1.11659    0.05893  18.948  < 2e-16 ***
#> speed50     -0.57519    0.13021  -4.417 9.99e-06 ***
#> ln(alpha)   -6.80676   47.43429  -0.143    0.886    
#> 
#> Random-effect standard deviations and correlations (Wald SEs):
#>                   Estimate Std. Error
#> sd_(Intercept)|ID   0.6388    0.06619
#> 
#> Log-likelihood: -1069.365   AIC: 2148.73   BIC: 2175.30
VarCorr(h)
#>  Groups Name        Variance Std.Dev.
#>  ID     (Intercept) 0.408    0.6388  
# }
```
