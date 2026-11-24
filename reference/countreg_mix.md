# Finite-mixture and latent-class count regression

`countreg_mix()` estimates finite mixtures of count distributions and
finite-mixture (latent-class) count regressions by maximum likelihood.
`countreg_lc()` is the same model under the latent-class vocabulary
(`classes`, `membership_formula`): both functions normalize their
arguments into one specification and call one estimation engine, so
equivalent calls give the same fit.

For unit \\u\\ (a row, or a panel when `group_var` is given) the
probability of the observed counts is \$\$L_u = \sum\_{k=1}^K \pi\_{uk}
\prod\_{i \in u} f_k(y_i \mid x_i, \theta_k),\$\$ where \\f_k\\ is the
probability mass function of component \\k\\ (any family of
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md);
see
[count_families](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md))
and \\\pi\_{uk}\\ are the mixing (class-membership) probabilities. The
log-likelihood is \\\sum_u w_u \log L_u\\, with frequency weights
\\w_u\\.

## Usage

``` r
countreg_mix(
  formula,
  data,
  family = "POISSON",
  components = 2L,
  mixing_formula = ~1,
  fixed_mixing = NULL,
  shared = NULL,
  dis_param_formula_1 = NULL,
  dis_param_formula_2 = NULL,
  offset = NULL,
  weights = NULL,
  group_var = NULL,
  method = "BFGS",
  max.iters = 1000,
  start.vals = NULL,
  n_starts = 10L,
  seed = NULL,
  stderr = c("normal", "robust", "boot"),
  bootstraps = NULL,
  ndraws = 1500,
  engine = "quadrature",
  keep_data = TRUE,
  verbose = FALSE,
  control = mixture_control()
)

countreg_lc(
  formula,
  data,
  family = "POISSON",
  classes = 2L,
  membership_formula = ~1,
  fixed_mixing = NULL,
  shared = NULL,
  dis_param_formula_1 = NULL,
  dis_param_formula_2 = NULL,
  offset = NULL,
  weights = NULL,
  group_var = NULL,
  method = "BFGS",
  max.iters = 1000,
  start.vals = NULL,
  n_starts = 10L,
  seed = NULL,
  stderr = c("normal", "robust", "boot"),
  bootstraps = NULL,
  ndraws = 1500,
  engine = "quadrature",
  keep_data = TRUE,
  verbose = FALSE,
  control = lc_control()
)
```

## Arguments

- formula:

  a model formula for the counts, as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md).
  Use `y ~ 1` for a mixture of distributions without covariates. Offsets
  can be given with [`offset()`](https://rdrr.io/r/stats/offset.html) or
  the `offset` argument; the same offset applies to every component.

- data:

  a data frame.

- family:

  a family code (see
  [count_families](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md))
  used by every component, or one code per component (e.g.
  `c(base = "POISSON", tail = "NB2")`): a jointly estimated mixture of
  different families. Names become display labels of the components.

- components, classes:

  the number of components (classes) \\K\\. When `family` has several
  entries, it can be omitted.

- mixing_formula, membership_formula:

  a one-sided formula for the mixing (class-membership) probabilities, a
  multinomial logit with the last component as reference: `~ 1`
  (default) gives constant proportions, `~ z` lets them depend on `z` (a
  mixture of experts). For panel models, its variables must be constant
  within each panel.

- fixed_mixing:

  optional known mixing probabilities (one per component, nonnegative,
  summing to 1). They are not estimated and add no parameters. These are
  *not* frequency weights or starting values. Zero probabilities are not
  allowed (the component would be unidentified).

- shared:

  optional list saying which parameters are the same in all components:
  `mean`, a character vector of formula terms (use `"(Intercept)"` for
  the intercept; a factor is shared as a whole term), and `dis1`, `dis2`
  (`TRUE`/`FALSE`) for the blocks of the first and second distribution
  parameter. Distribution parameters can only be shared between
  components of the same family. By default nothing is shared.

- dis_param_formula_1, dis_param_formula_2:

  optional one-sided formulas for the first and second distribution
  parameter of the components, as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md).
  They apply to every component whose family has the parameter
  (coefficients are component-specific unless shared).

- offset, weights:

  an offset (column names, log scale) and frequency weights, as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md).
  With `group_var`, weights must be constant within each panel: a weight
  stands for that many independent copies of the panel.

- group_var:

  optional column name(s) identifying panels. All rows of a panel belong
  to the same (unknown) component, so the panel, not the row, is the
  independent unit; rows are independent given the component.

- method:

  optimizer for the direct maximization of the observed log-likelihood
  (see
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md));
  default `"BFGS"`.

- max.iters:

  maximum number of iterations of the direct maximization.

- start.vals:

  optional starting values: a vector with one value per parameter (named
  as in [`coef()`](https://rdrr.io/r/stats/coef.html), or in that
  order). It is used as the first start; the others are generated.

- n_starts:

  number of starts. Starts are reproducible from `seed` and are recorded
  in the fitted object (`fit$model$starts`).

- seed:

  optional random seed for the generated starts (and the bootstrap). The
  caller's random-number stream is not changed. When `NULL`, a seed is
  drawn from the current stream and stored.

- stderr:

  `"normal"` (from the Hessian of the observed-data log-likelihood),
  `"robust"` (sandwich with unit-level scores) or `"boot"`
  (nonparametric bootstrap of units, refitting the full mixture with
  several starts and aligning component labels).

- bootstraps:

  number of bootstrap replicates for `stderr = "boot"` (default 200).

- ndraws, engine:

  numerical settings of the mixed-Poisson families (PLN, PGE, PW, PLL)
  within each component, as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md).
  They are unrelated to the number of components.

- keep_data:

  logical; store the training data in the fitted object.

- verbose:

  logical; report the progress of each start.

- control:

  a
  [`mixture_control()`](https://jwood-iastate.github.io/flexCountReg/reference/mixture_control.md)
  object.

## Value

An object of class `c("countreg_mix", "flexCountReg")`
(`countreg_mix()`) or `c("countreg_lc", "countreg_mix", "flexCountReg")`
(`countreg_lc()`). The methods
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`predict()`](https://rdrr.io/r/stats/predict.html)
([`predict.countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.countreg_mix.md)),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html) and
[`augment()`](https://generics.r-lib.org/reference/augment.html) work as
for other models, and
[`component_probs()`](https://jwood-iastate.github.io/flexCountReg/reference/component_probs.md)
returns prior or posterior class probabilities.
[`rootogram()`](https://jwood-iastate.github.io/flexCountReg/reference/rootogram.md),
[`cureplot()`](https://jwood-iastate.github.io/flexCountReg/reference/cureplot.md)
and
[`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
(delta-method standard errors) accept mixture models.

## Details

### Estimation

Each start is a response-based or random preliminary partition of the
units, turned into parameters by one M-step. Starts are improved by the
EM algorithm (posterior class probabilities in the E-step; weighted
updates of the component and mixing parameters in the M-step, with
shared parameters updated jointly), and the best distinct solutions are
polished by maximizing the observed-data log-likelihood directly with
`maxLik`. The selected solution is the best polished one that converged.
The start table (`fit$model$starts`) records every start's strategy,
log-likelihood, EM iterations and status, and whether the best solution
was reached by several starts; mixtures often have several local maxima,
so a solution found by one start only deserves caution.

Components of the same family are exchangeable: their numbering is
arbitrary. After estimation they are numbered in increasing order of
their average expected count over the training rows (a display
convention). Mixtures of different families keep the order given in
`family`.

### Parameters

The coefficients are named `Component<k>:Mean:<term>`,
`Component<k>:<parameter>` (e.g. `Component2:ln(alpha)`), `Shared:...`
for shared parameters, and `Mixing:Component<k>:<term>` for the
membership logit of component \\k\\ against the last component.
Component expected counts use each family's own mean (for PLN and PLL it
includes the lognormal factor), so shared link-scale coefficients do not
imply equal component means.

### Inference and model choice

The covariance matrix is the inverse of the observed information of the
whole free parameter vector (not of the separate M-steps). Near a
boundary (an almost empty component, or two almost identical components)
it is unreliable; use the bootstrap and the start table. The usual
chi-square likelihood-ratio test does **not** apply to comparing numbers
of components (K versus K + 1); compare such models with AIC/BIC on the
same sample, start stability and held-out log scores.
[`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md)
refuses mixture models for this reason.

### Interpretation

Classes are a modeling device. A good fit with two components does not
show that two real subpopulations exist; the same mixture can describe a
skewed or heavy-tailed distribution.

Zero-inflated, hurdle and underreporting components, component-specific
formulas, and confidence/prediction/tolerance intervals for mixture
predictions are not available yet.

## References

McLachlan, G. J., & Peel, D. (2000). *Finite Mixture Models*. Wiley.

Grün, B., & Leisch, F. (2008). FlexMix version 2: Finite mixtures with
concomitant variables and varying and constant parameters. *Journal of
Statistical Software*, 28(4), 1-35.

Wedel, M., DeSarbo, W. S., Bult, J. R., & Ramaswamy, V. (1993). A latent
class Poisson regression model for heterogeneous count data. *Journal of
Applied Econometrics*, 8(4), 397-411.

## See also

[`mixture_control()`](https://jwood-iastate.github.io/flexCountReg/reference/mixture_control.md),
[`component_probs()`](https://jwood-iastate.github.io/flexCountReg/reference/component_probs.md),
[`predict.countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.countreg_mix.md),
[count_families](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md)

## Examples

``` r
# A two-component Poisson mixture of simulated counts
set.seed(1)
d <- data.frame(y = c(rpois(300, 1), rpois(200, 6)))
m <- countreg_mix(y ~ 1, data = d, family = "POISSON", components = 2,
                  n_starts = 5, seed = 42)
summary(m)
#> Call:
#> countreg_mix(formula = y ~ 1, data = d, family = "POISSON", components = 2, 
#>     n_starts = 5, seed = 42)
#> 
#> Finite-mixture count model with 2 components
#>   Component 1: POISSON
#>   Component 2: POISSON
#> Mixing: constant proportions (reference: Component 2)
#> Optimizer: BFGS + Newton-Raphson refinement | gradient close to zero (gradtol) 
#> 
#> Coefficients (standard errors: observed information (Hessian)):
#> 
#> Component 1 (POISSON):
#>                             Estimate Std. Error z value Pr(>|z|)
#> Component1:Mean:(Intercept)  -0.0361     0.0828   -0.44     0.66
#> 
#> Component 2 (POISSON):
#>                             Estimate Std. Error z value Pr(>|z|)    
#> Component2:Mean:(Intercept)   1.8054     0.0374    48.2   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Mixing (multinomial logit):
#>                               Estimate Std. Error z value Pr(>|z|)    
#> Mixing:Component1:(Intercept)    0.407      0.117    3.49  0.00048 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Component sizes:
#>       label  family  prior posterior effective_n modal_n mean_count
#>  Component1 POISSON 0.6003    0.6003       300.2     325     0.9646
#>  Component2 POISSON 0.3997    0.3997       199.8     175     6.0830
#> (prior: mean membership probability; posterior: mean posterior probability;
#>  effective_n: sum of posterior probabilities; modal_n: units whose most
#>  probable component it is; mean_count: average expected count)
#> 
#> Log-likelihood: -1107.547 | AIC: 2221.094 | BIC: 2233.738 | ICL: 2398.952 
#> N: 500
#> Entropy: 82.607 | relative entropy: 0.762 
#> Starts: 5 | 5 reached the selected log-likelihood
#> Observed information: positive definite (condition number 14.7 ) | max |gradient|: 4.82e-08 
#> No likelihood-ratio test of the number of components is reported (non-standard distribution).
component_probs(m, conditioning = "prior")[1, ]
#> Component1 Component2 
#>   0.600347   0.399653 

# \donttest{
# The same model through the latent-class interface
lc <- countreg_lc(y ~ 1, data = d, family = "POISSON", classes = 2,
                  n_starts = 5, seed = 42)
all.equal(logLik(m), logLik(lc))
#> [1] TRUE

# Latent-class NB2 regression with a membership covariate
data("abalone_rings")
fit <- countreg_lc(rings ~ shell_weight, data = abalone_rings,
                   family = "NB2", classes = 2,
                   membership_formula = ~ sex, n_starts = 5, seed = 1)
summary(fit)
#> Call:
#> countreg_lc(formula = rings ~ shell_weight, data = abalone_rings, 
#>     family = "NB2", classes = 2, membership_formula = ~sex, n_starts = 5, 
#>     seed = 1)
#> 
#> Latent-class count model with 2 classes
#>   Class 1: NB2
#>   Class 2: NB2
#> Mixing: multinomial logit on covariates (reference: Class 2)
#> Optimizer: BFGS + Newton-Raphson refinement | successive function values within relative tolerance limit (reltol) 
#> 
#> Coefficients (standard errors: observed information (Hessian)):
#> 
#> Class 1 (NB2):
#>                              Estimate Std. Error z value Pr(>|z|)    
#> Component1:Mean:(Intercept)    2.0974     0.0143   146.4   <2e-16 ***
#> Component1:Mean:shell_weight   0.9677     0.0430    22.5   <2e-16 ***
#> Component1:ln(alpha)         -20.6986   209.9017    -0.1     0.92    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Class 2 (NB2):
#>                              Estimate Std. Error z value Pr(>|z|)    
#> Component2:Mean:(Intercept)    1.7141     0.0189   90.66   <2e-16 ***
#> Component2:Mean:shell_weight   2.5623     0.1110   23.09   <2e-16 ***
#> Component2:ln(alpha)         -20.2632   380.0516   -0.05     0.96    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Mixing (multinomial logit):
#>                               Estimate Std. Error z value Pr(>|z|)    
#> Mixing:Component1:(Intercept)    4.323      0.772    5.60  2.1e-08 ***
#> Mixing:Component1:sexI          -9.328      1.565   -5.96  2.5e-09 ***
#> Mixing:Component1:sexM           8.407    120.117    0.07     0.94    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Class sizes:
#>   label family  prior posterior effective_n modal_n mean_count
#>  Class1    NB2 0.6768    0.6768      2826.8    2834      10.36
#>  Class2    NB2 0.3232    0.3232      1350.2    1343      10.97
#> (prior: mean membership probability; posterior: mean posterior probability;
#>  effective_n: sum of posterior probabilities; modal_n: units whose most
#>  probable class it is; mean_count: average expected count)
#> 
#> Log-likelihood: -9681.481 | AIC: 19380.96 | BIC: 19438 | ICL: 19701.21 
#> N: 4177
#> Entropy: 131.607 | relative entropy: 0.955 
#> Starts: 5 | 1 reached the selected log-likelihood
#> Note: only one start reached this solution; try more starts.
#> Observed information: singular | max |gradient|: 0.0019 
#> No likelihood-ratio test of the number of classes is reported (non-standard distribution).
# }
```
