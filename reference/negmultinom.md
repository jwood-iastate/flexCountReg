# Negative Multinomial Family Panel Regression (NM1, NM2, NM-P)

Estimate a negative multinomial (NM) regression model for clustered or
panel count data. Within each cluster (e.g., a site observed over
several years), the counts are conditionally independent Poisson
variables that share a common gamma-distributed multiplicative frailty
with mean 1. Integrating over the frailty gives the negative multinomial
likelihood in closed form.

Three members of the family are available, mirroring the NB-1, NB-2, and
NB-P forms of the negative binomial in
[`countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md):

- `"NM2"` (default): the conventional negative multinomial. Every
  observation, and every within-cluster aggregate, has an NB-2 marginal
  with the same \\\alpha\\. This is the model estimated by
  [`repois`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md).

- `"NM1"`: the frailty variance is calibrated so that the cluster total
  follows an NB-1 distribution.

- `"NMP"`: the frailty variance is calibrated so that the cluster total
  follows an NB-P distribution with power parameter \\p\\.

Both \\\ln(\alpha)\\ and \\p\\ can be specified as functions of
cluster-level variables.

## Usage

``` r
negmultinom(
  formula,
  group_var,
  data,
  family = "NM2",
  offset = NULL,
  dis_param_formula_1 = NULL,
  dis_param_formula_2 = NULL,
  method = "BFGS",
  max.iters = 1000,
  start.vals = NULL,
  stderr = "normal",
  bootstraps = NULL,
  verbose = FALSE,
  keep_data = TRUE,
  weights = NULL,
  sp = NULL
)
```

## Arguments

- formula:

  an R formula for the mean, e.g. `y ~ x1 + x2`. An
  [`offset()`](https://rdrr.io/r/stats/offset.html) term may be included
  in the formula.

- group_var:

  character vector naming the variable(s) that identify the clusters
  (e.g., a site or panel ID). If more than one variable is given,
  clusters are defined by their interaction.

- data:

  a data frame containing all variables used in the model.

- family:

  the member of the negative multinomial family to estimate: `"NM2"`
  (default), `"NM1"`, or `"NMP"`. Case and non-alphanumeric characters
  are ignored, so `"NM-P"` and `"nm2"` are accepted.

- offset:

  optional character vector naming offset variable(s) in `data`.
  Multiple offsets are summed. Do not also use
  [`offset()`](https://rdrr.io/r/stats/offset.html) in `formula` for the
  same variable.

- dis_param_formula_1:

  optional one-sided formula for \\\ln(\alpha)\\, e.g. `~ x3`. If
  `NULL`, a single constant \\\ln(\alpha)\\ is estimated. See Details on
  within-cluster variation.

- dis_param_formula_2:

  optional one-sided formula for \\p\\ (only used when
  `family = "NMP"`). The identity link is used, so \\p = W\delta\\. If
  `NULL`, a single constant \\p\\ is estimated.

- method:

  optimization method passed to
  [`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html). One of
  `"BFGS"` (default), `"NR"`, `"BHHH"`, `"NM"`, `"CG"`, `"BFGSR"`, or
  `"SANN"`. An analytic gradient is supplied to the gradient-based
  methods.

- max.iters:

  maximum number of iterations for the optimizer.

- start.vals:

  optional numeric vector of starting values, ordered as the mean
  coefficients, then the \\\ln(\alpha)\\ coefficient(s), then (for
  `"NMP"`) the \\p\\ coefficient(s).

- stderr:

  type of standard errors: `"normal"` (default; inverse of the negative
  Hessian), `"robust"` (sandwich estimator using cluster-level scores,
  i.e., cluster-robust), or `"boot"` (cluster bootstrap; requires
  `bootstraps`).

- bootstraps:

  integer number of cluster bootstrap replications used when
  `stderr = "boot"`.

- verbose:

  logical. If `TRUE`, optimizer progress is printed.

- keep_data:

  logical; if `FALSE`, the training data are not stored in the fitted
  object (saving memory).
  [`predict()`](https://rdrr.io/r/stats/predict.html) then needs
  `newdata`, and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)/[`residuals()`](https://rdrr.io/r/stats/residuals.html)
  are unavailable.

- weights:

  the name of a column of `data`, or a numeric vector with one value per
  row of `data`, of frequency weights of the panels: a panel with weight
  w counts as w independent panels in the log-likelihood, the robust
  covariance, BIC and the bootstrap. Weights must be constant within
  each panel. Bootstrap standard errors need whole-number weights.

- sp:

  optional fixed smoothing parameters for smooth terms (`s()`, `te()`,
  ...) in the formula(s); by default they are selected automatically.
  Smooth terms are fitted by penalized maximum likelihood with Bayesian
  standard errors (see the "Smooth terms" section of
  [`countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md));
  bootstrap standard errors are not available with them.

## Value

An object of class `flexCountReg`, a list with components:

- `model`: the fitted model. Along with the `maxLik` output, this
  includes `estimate`, `se`, `vcov`, `family`, `alpha` and `p`
  (cluster-level values), `fitted.values` (marginal means
  \\\mu\_{it}\\), `eb.fitted` (EB expected frequencies),
  `frailty.posterior` (cluster posterior mean of \\U_i\\), `n_groups`,
  `nobs`, `aic`, and `bic`.

- `data`: the data frame used to fit the model.

- `call`: the matched call.

- `formula`: the formula used to fit the model.

## Details

### Model structure

Let \\i\\ index clusters and \\t\\ index the observations within a
cluster. The mean is modeled as \$\$\mu\_{it} = \exp(X\_{it}\beta +
\text{offset}\_{it}).\$\$ Conditional on a cluster frailty \\U_i\\,
\$\$Y\_{it} \mid U_i \sim \text{Poisson}(\mu\_{it} U_i),\$\$
independently over \\t\\, and \\U_i \sim \text{Gamma}(r_i, 1/r_i)\\ so
that \\E(U_i) = 1\\ and \\\text{Var}(U_i) = 1/r_i\\.

Define the cluster totals over the reference collection \\R\\ (all of a
cluster's observations used in estimation): \$\$y\_{iR} = \sum_t
y\_{it}, \qquad \mu\_{iR} = \sum_t \mu\_{it}.\$\$

The three family members differ only in the gamma shape parameter
\\r_i\\: \$\$\text{NM2: } r_i = \frac{1}{\alpha_i}, \qquad \text{NM1: }
r_i = \frac{\mu\_{iR}}{\alpha_i}, \qquad \text{NM-P: } r_i =
\frac{\mu\_{iR}^{2-p_i}}{\alpha_i}.\$\$ NM2 and NM1 are the special
cases \\p = 2\\ and \\p = 1\\ of NM-P.

### Probability mass function

Integrating over \\U_i\\ gives the negative multinomial PMF:
\$\$P(\mathbf{y}\_i \mid \boldsymbol{\mu}\_i) = \frac{\Gamma(r_i +
y\_{iR})}{\Gamma(r_i)\prod_t y\_{it}!} \left(\frac{r_i}{r_i +
\mu\_{iR}}\right)^{r_i} \prod_t \left(\frac{\mu\_{it}}{r_i +
\mu\_{iR}}\right)^{y\_{it}}.\$\$

The log-likelihood for cluster \\i\\ is: \$\$LL_i = \ln\Gamma(r_i +
y\_{iR}) - \ln\Gamma(r_i) + r_i\ln(r_i) - (r_i + y\_{iR})\ln(r_i +
\mu\_{iR}) + \sum_t y\_{it}\ln(\mu\_{it}) - \sum_t \ln(y\_{it}!).\$\$

### Mean and variance

For every family member, \\E(Y\_{it}) = \mu\_{it}\\. The variance and
covariance are: \$\$\text{Var}(Y\_{it}) = \mu\_{it} +
\frac{\mu\_{it}^2}{r_i}, \qquad \text{Cov}(Y\_{is}, Y\_{it}) =
\frac{\mu\_{is}\mu\_{it}}{r_i}.\$\$ The cluster total therefore
satisfies \$\$\text{Var}(Y\_{iR}) = \mu\_{iR} +
\alpha_i\mu\_{iR}^{p_i},\$\$ i.e., it is NB-2, NB-1, or NB-P for
`"NM2"`, `"NM1"`, and `"NMP"`, respectively.

### Interpretation of NM1 and NM-P

For NM2, \\r_i = 1/\alpha_i\\ does not depend on the mean, so each
individual observation (and any within-cluster subtotal) is NB-2 with
the same \\\alpha\\. The model is coherent across temporal aggregation.

For \\p \neq 2\\, the frailty distribution depends on \\\mu\_{iR}\\.
Only the cluster total over the reference collection has the NB-P (or
NB-1) variance function; individual observations have
\\\text{Var}(Y\_{it}) = \mu\_{it} + \alpha_i\mu\_{iR}^{p_i -
2}\mu\_{it}^2\\, which is not NB-P with the same \\\alpha\\ and \\p\\.
The reference collection used here is all rows of a cluster that remain
after removing missing values. Changing the set of observations in a
cluster (e.g., adding or dropping years) changes the implied frailty
distribution, so NM1 and NM-P should be interpreted as
reference-window-specific joint models.

### Parameters as functions of variables

When `dis_param_formula_1` or `dis_param_formula_2` is supplied,
\$\$\ln(\alpha_i) = \bar{Z}\_i\gamma, \qquad p_i = \bar{W}\_i\delta,\$\$
where \\\bar{Z}\_i\\ and \\\bar{W}\_i\\ are the within-cluster means of
the corresponding design matrices. Because \\\alpha\\ and \\p\\ describe
the shared cluster frailty, they must be constant within a cluster.
Cluster-level variables (constant within a cluster) are therefore
recommended; if a variable varies within a cluster, a warning is issued
and its cluster mean is used.

The \\\ln(\alpha)\\ parameters use a log link so that \\\alpha \> 0\\.
The \\p\\ parameters use the identity link, so \\p\\ is unrestricted.

### Empirical Bayes estimates

Given the cluster history, the posterior of the frailty is \\U_i \mid
\mathbf{y}\_i \sim \text{Gamma}(r_i + y\_{iR}, 1/(r_i + \mu\_{iR}))\\.
The returned model includes the empirical Bayes (EB) expected frequency
for each observation: \$\$E(\Lambda\_{it} \mid \mathbf{y}\_i) =
\mu\_{it}\frac{r_i + y\_{iR}}{r_i + \mu\_{iR}}.\$\$ For NM2 this is
\\\mu\_{it}(1 + \alpha y\_{iR})/(1 + \alpha\mu\_{iR})\\. EB estimates
for new data can be obtained with
`predict(model, newdata, method = "EB")` when `newdata` contains the
outcome and grouping variables.

### Standard errors

The log-likelihood is evaluated cluster by cluster, so the `"robust"`
sandwich estimator and the `"BHHH"` optimizer both use cluster-level
scores. The bootstrap resamples whole clusters.

## References

Cameron, A. C., & Trivedi, P. K. (2013). Regression Analysis of Count
Data (2nd ed.). Cambridge University Press.

Greene, W. (2008). Functional forms for the negative binomial model for
count data. Economics Letters, 99(3), 585-590.

Hausman, J., Hall, B. H., & Griliches, Z. (1984). Econometric models for
count data with an application to the patents-R&D relationship.
Econometrica, 52(4), 909-938.

## See also

[`repois`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md),
[`countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
[`predict.flexCountReg`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md)

## Examples

``` r
# \donttest{
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT > 10000, 1, 0)

# Conventional negative multinomial (NM2)
nm2 <- negmultinom(Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k,
                   group_var = "ID", data = washington_roads,
                   offset = "lnlength")
summary(nm2)
#> Call:
#> negmultinom(formula = Animal ~ lnaadt + speed50 + ShouldWidth04 + 
#>     AADTover10k, group_var = "ID", data = washington_roads, offset = "lnlength")
#> 
#> Model: negmultinom | family: NM2 | optimizer: BFGS
#> Iterations: 34 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)    -9.7059     1.3624   -7.12  1.0e-12 ***
#> lnaadt          1.0035     0.1589    6.32  2.7e-10 ***
#> speed50        -0.8214     0.3274   -2.51    0.012 *  
#> ShouldWidth04  -0.4252     0.2813   -1.51    0.131    
#> AADTover10k    -0.9110     0.5132   -1.78    0.076 .  
#> ln(alpha)      -0.0226     0.4723   -0.05    0.962    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Offset: lnlength (coefficient fixed at 1)
#> 
#> Log-likelihood: -265.124 | AIC: 542.248 | BIC: 567.619 
#> N: 1501 | panels: 507

# NM1: cluster totals follow an NB-1 distribution
nm1 <- negmultinom(Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k,
                   group_var = "ID", data = washington_roads,
                   family = "NM1", offset = "lnlength")
summary(nm1)
#> Call:
#> negmultinom(formula = Animal ~ lnaadt + speed50 + ShouldWidth04 + 
#>     AADTover10k, group_var = "ID", data = washington_roads, family = "NM1", 
#>     offset = "lnlength")
#> 
#> Model: negmultinom | family: NM1 | optimizer: BFGS
#> Iterations: 30 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)     -9.263      1.374   -6.74  1.6e-11 ***
#> lnaadt           0.959      0.159    6.03  1.7e-09 ***
#> speed50         -0.919      0.334   -2.75   0.0059 ** 
#> ShouldWidth04   -0.436      0.279   -1.56   0.1179    
#> AADTover10k     -0.896      0.549   -1.63   0.1028    
#> ln(alpha)       -0.722      0.438   -1.65   0.0988 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Offset: lnlength (coefficient fixed at 1)
#> 
#> Log-likelihood: -264.762 | AIC: 541.524 | BIC: 566.895 
#> N: 1501 | panels: 507

# NM-P with ln(alpha) as a function of a site-level variable and
# cluster-robust standard errors
nmp <- negmultinom(Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k,
                   group_var = "ID", data = washington_roads,
                   family = "NMP", offset = "lnlength",
                   dis_param_formula_1 = ~ speed50,
                   stderr = "robust")
summary(nmp)
#> Call:
#> negmultinom(formula = Animal ~ lnaadt + speed50 + ShouldWidth04 + 
#>     AADTover10k, group_var = "ID", data = washington_roads, family = "NMP", 
#>     offset = "lnlength", dis_param_formula_1 = ~speed50, stderr = "robust")
#> 
#> Model: negmultinom | family: NMP | optimizer: BFGS
#> Iterations: 41 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: robust (sandwich)):
#>                       Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)             -9.940      1.408   -7.06  1.7e-12 ***
#> lnaadt                   1.032      0.163    6.33  2.5e-10 ***
#> speed50                 -0.748      0.465   -1.61     0.11    
#> ShouldWidth04           -0.429      0.284   -1.51     0.13    
#> AADTover10k             -0.888      0.557   -1.59     0.11    
#> ln(alpha):(Intercept)   -0.354      0.817   -0.43     0.66    
#> ln(alpha):speed50        2.330      5.176    0.45     0.65    
#> p                        2.226      2.538    0.88     0.38    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Offset: lnlength (coefficient fixed at 1)
#> 
#> Log-likelihood: -263.556 | AIC: 543.111 | BIC: 576.94 
#> N: 1501 | panels: 507

# Marginal means and empirical Bayes expected frequencies
head(predict(nm2))
#> [1] 0.09294448 0.08213698 0.13617446 0.03026099 0.09726747 0.04539149
head(predict(nm2, method = "EB"))
#> [1] 0.07283897 0.06603019 0.28655720 0.02776571 0.07546753 0.03999942
# }
```
