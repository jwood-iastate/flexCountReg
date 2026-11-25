# Multivariate count regression

Fits several count outcomes observed on the same units, such as crashes
by severity or by type at each site, each with its own count family,
mean formula and distribution-parameter formulas.

## Usage

``` r
countreg_mv(
  formula,
  data,
  family = "NB2",
  dis_param_formula_1 = NULL,
  dis_param_formula_2 = NULL,
  offset = NULL,
  weights = NULL,
  group_var = NULL,
  stderr = "normal",
  method = "BFGS",
  max.iters = 1000,
  start.vals = NULL,
  ndraws = 1500,
  engine = "quadrature",
  verbose = FALSE,
  rpar_formula = NULL,
  correlation = c("across", "within", "none"),
  rpar_draws = 500,
  seed = NULL,
  copula = c("none", "gaussian", "frank", "clayton", "gumbel"),
  .setup_only = FALSE
)
```

## Arguments

- formula:

  either a formula with the outcomes on the left as
  `cbind(y_1, ..., y_J)` and one right-hand side for every equation, or
  a list of two-sided formulas, one per outcome (named by the list
  names, or by their responses). An
  [`offset()`](https://rdrr.io/r/stats/offset.html) term (an exposure)
  belongs in the formulas.

- data:

  a data frame.

- family:

  the count family of each outcome (see
  [count_families](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md)):
  one code for every outcome, or one per outcome (named by outcome or in
  their order). Families made with
  [`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md)
  are accepted.

- dis_param_formula_1, dis_param_formula_2:

  formulas for the outcomes' distribution parameters, as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md):
  `NULL` (constant), one one-sided formula for every outcome that has
  the parameter, or a list named by outcome (outcomes not in the list
  keep a constant parameter).

- offset:

  optional names of offset columns (log scale) used by every equation,
  instead of [`offset()`](https://rdrr.io/r/stats/offset.html) in the
  formulas.

- weights:

  optional frequency weights: a column name or a numeric vector.

- group_var:

  optional column naming the units for cluster-robust standard errors
  (such as a site observed over several years).

- stderr:

  `"normal"` (from the Hessian, the default), `"robust"` (sandwich) or
  `"cluster"` (sandwich by `group_var`).

- method, max.iters, start.vals, verbose:

  optimizer settings, as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md);
  the default method is `"BFGS"`.

- ndraws, engine:

  integration settings of families without a closed form, as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md).

- rpar_formula:

  optional one-sided formula of random parameters: one formula for every
  outcome, or a list named by outcome. Each term adds a normal random
  deviation, with mean zero, to the coefficient of that variable in the
  outcome's equation (a variable not in the equation's formula gets a
  random coefficient with mean zero). The intercept is random only when
  the formula has the term `1` itself (`~ 1`, `~ 1 + x`). The random
  parameters are drawn once per unit of `group_var` (shared by all its
  rows, such as the years of a site), or once per row without
  `group_var`.

- correlation:

  correlation of the random parameters: `"across"` (default; correlated
  within and across equations, so the outcomes of a unit share
  unobserved factors), `"within"` (correlated within each equation only)
  or `"none"`.

- rpar_draws:

  the number of simulation draws per unit (randomly shifted Halton
  points, a separate block and shift for each unit).

- seed:

  seed of the random shifts of the draws (stored in the fit, so
  predictions use the same draws); random by default.

- copula:

  `"none"` (default), or a copula that joins the outcomes'
  distributions: `"gaussian"` (two outcomes), `"frank"`, `"clayton"` or
  `"gumbel"`. See Details.

- .setup_only:

  internal use: return the model setup without fitting.

## Value

An object of class `"countreg_mv"`, with methods
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html).

## Details

**Independent equations.** Without `rpar_formula`, the outcomes are
independent given the covariates, so the log-likelihood is the sum of
the equations' log-likelihoods and every equation's estimates equal
those of a separate
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
fit. Fitting them together gives one object with joint predictions,
simulations and tests across equations (for example, a Wald test that a
coefficient is the same in two equations). The rows used are those
complete for every equation.

**Correlated outcomes.** Counts of different categories at a site are
usually correlated: unobserved conditions that raise one also raise the
others. With `rpar_formula`, each equation's log mean gets random terms
\\x_r' b\_{uj}\\, with \\b_u = L z_u\\ drawn per unit and correlated
across equations (`correlation = "across"`).

**The multivariate Poisson-lognormal model.** Poisson outcomes with
correlated normal random intercepts (`family = "POISSON"`,
`rpar_formula = ~ 1`, `correlation = "across"`) give the multivariate
Poisson-lognormal model (MVPLN) of, for example, Ma et al. (2008). Other
families keep their own overdispersion within each outcome, with the
random intercepts adding the correlation between outcomes, and every
equation can have a different family. A positive correlation of the
random intercepts means that a site with more crashes of one category
than its covariates predict also tends to have more of the other. The
correlation can be negative too, unlike in models with a shared gamma
effect such as
[`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md).
The likelihood is simulated (`rpar_draws` draws per unit) with analytic
gradients; the elements of \\L\\ are reported as `sd:` and `chol:`, and
[`summary()`](https://rdrr.io/r/base/summary.html) shows the implied
standard deviations and correlations.
[`predict()`](https://rdrr.io/r/stats/predict.html) gives joint and
conditional probabilities of the outcomes, and Empirical Bayes
predictions for units in the estimation data.

**Copulas.** A copula joins the outcomes' own distributions, each with
its family, into one joint distribution: \\P(Y_1 \le y_1, \dots, Y_J \le
y_J) = C(F_1(y_1), \dots, F_J(y_J); \theta)\\. For counts, the
probability of the observed combination is the exact volume of the
copula over the rectangle \\\[F_j(y_j - 1), F_j(y_j)\]\\ (the sum over
its \\2^J\\ corners; Nikoloulopoulos and Karlis 2009), not a continuous
approximation such as the midpoint or jittering. The margins are those
of the equations, so expected counts and each outcome's probabilities do
not change; the copula changes the joint and conditional probabilities.
The families differ in the shape of the dependence:

- `"gaussian"`: symmetric, without tail dependence, any sign
  (`copula:atanh(rho)`); two outcomes.

- `"frank"`: symmetric, without tail dependence; any sign for two
  outcomes (`copula:theta`), positive for more (`copula:ln(theta)`).

- `"clayton"`: positive, stronger in the lower tail (joint low counts;
  `copula:ln(theta)`).

- `"gumbel"`: positive, stronger in the upper tail (joint high counts;
  `copula:ln(theta - 1)`).

All parameters are estimated together by maximum likelihood, started
from the separate equations and the copula parameter that is best for
them (inference functions for margins); the gradient is numerical.
[`summary()`](https://rdrr.io/r/base/summary.html) reports the copula
parameter and its Kendall's tau. With counts, the copula is not unique
outside the range of the margins (Genest and Neslehova 2007): its
parameter describes this model's dependence, and the tau is that of the
copula, not the rank correlation of the counts, which is weaker for
small counts. Copulas cannot yet be combined with `rpar_formula`.

Different outcomes often need different families: frequent property
damage crashes may be close to NB2, while rare fatal crashes are more
skewed and zero-heavy. Each equation can use any family of the package.

## References

Genest, C., and Neslehova, J. (2007). A primer on copulas for count
data. *ASTIN Bulletin*, 37(2), 475-515.

Ma, J., Kockelman, K. M., and Damien, P. (2008). A multivariate
Poisson-lognormal regression model for prediction of crash counts by
severity, using Bayesian methods. *Accident Analysis and Prevention*,
40(3), 964-975.

Nikoloulopoulos, A. K., and Karlis, D. (2009). Modeling multivariate
count data using copulas. *Communications in Statistics - Simulation and
Computation*, 39(1), 172-187.

## See also

[`predict.countreg_mv()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.countreg_mv.md),
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
[`countreg_alloc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_alloc.md)
for a total and its split across categories,
[`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md)

## Examples

``` r
# \donttest{
data("washington_roads")
w <- washington_roads
w$pdo <- w$Total_crashes - w$Injury_crashes
m <- countreg_mv(cbind(pdo, Injury_crashes) ~ lnaadt + speed50 +
                   offset(lnlength), data = w,
                 family = c(pdo = "NB2", Injury_crashes = "PIG"))
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
#> Injury_crashes (PIG):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -7.4115     1.1903  -6.227 4.77e-10 ***
#> lnaadt        0.6661     0.1418   4.698 2.63e-06 ***
#> speed50      -1.2201     0.4146  -2.943  0.00325 ** 
#> ln(eta)       0.7900     1.9897   0.397  0.69132    
#> 
#> Log-likelihood: -1247.0544  AIC: 2510.1088  BIC: 2552.6199 
head(predict(m))
#>            pdo Injury_crashes
#> [1,] 0.7745464     0.03005740
#> [2,] 0.6844829     0.02656235
#> [3,] 1.1348006     0.04403759
#> [4,] 0.2521779     0.00978613
#> [5,] 0.8105719     0.03145542
#> [6,] 0.3782669     0.01467920
# }
```
