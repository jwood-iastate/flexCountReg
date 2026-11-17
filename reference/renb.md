# Estimate a Random Effects Negative Binomial regression model

Estimate a Random Effects Negative Binomial regression model

## Usage

``` r
renb(
  formula,
  group_var,
  data,
  method = "NM",
  max.iters = 1000,
  verbose = FALSE,
  stderr = "normal",
  bootstraps = NULL,
  offset = NULL,
  weights = NULL,
  keep_data = TRUE,
  print.level = NULL
)
```

## Arguments

- formula:

  an R formula.

- group_var:

  the grouping variable(s) for the random effects (e.g., individual ID
  or other panel ID variables). Required. Several columns are combined
  into one panel identifier.

- data:

  a dataframe that has all of the variables in the `formula`.

- method:

  a method to use for optimization in the maximum likelihood estimation.
  For options, see
  [`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html). Note that
  "BHHH" is not available for this function due to the implementation
  for the random effects.

- max.iters:

  the maximum number of iterations to allow the optimization method to
  perform.

- verbose:

  logical; if `TRUE`, print the progress of the optimization.

- stderr:

  the type of standard errors: `"normal"` (Hessian-based, the default),
  `"robust"` (sandwich; clustered by panel because the likelihood has
  one contribution per panel), or `"boot"` (bootstrap that resamples
  whole panels). Giving `bootstraps` without `stderr` also requests
  bootstrapped standard errors.

- bootstraps:

  Optional integer specifying the number of bootstrap samples to be used
  for estimating standard errors. Whole panels (groups) are resampled,
  so the standard errors account for the correlation within panels. If
  not specified, no bootstrapping is performed.

- offset:

  an optional character vector of offset column names (summed when
  several are given; on the log scale). An offset can instead be given
  with [`offset()`](https://rdrr.io/r/stats/offset.html) in `formula`;
  use only one of the two.

- weights:

  the name of a column of `data`, or a numeric vector with one value per
  row of `data`, of frequency weights of the panels: a panel with weight
  w counts as w independent panels in the log-likelihood, the robust
  covariance, BIC and the bootstrap. Weights must be constant within
  each panel. Bootstrap standard errors need whole-number weights.

- keep_data:

  logical; if `FALSE`, the training data are not stored in the fitted
  object (saving memory).
  [`predict()`](https://rdrr.io/r/stats/predict.html) then needs
  `newdata`, and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)/[`residuals()`](https://rdrr.io/r/stats/residuals.html)
  are unavailable.

- print.level:

  deprecated; use `verbose`.

## Value

An object of class `flexCountReg` which is a list with the following
components:

- model: the fitted model object.

- data: the data frame used to fit the model.

- call: the matched call.

- formula: the formula used to fit the model.

## Details

This function estimates a random effects negative binomial (RENB)
regression model. This model is based on the NB-1 model. The PDF for the
RENB is: \$\$f(y\_{it}\|\lambda\_{it}, a, b) = \frac{\Gamma(a+b)
\Gamma(a + \sum\_{t = 1}^{n_i} \lambda\_{it}) \Gamma(b +
\sum\_{t=1}^{n_i}y\_{it})} {\Gamma(a) \Gamma(b) \Gamma(a + b +
\sum\_{t=1}^{n_i}\lambda\_{it} + \sum\_{t=1}^{n_i}y\_{it})}
\prod\_{t=1}^{n_i}
\frac{\Gamma(\lambda\_{it}+y\_{it})}{\Gamma(\lambda\_{it})\Gamma(y\_{it}+1)}\$\$

Where \\y\_{it}\\ is the count outcome for individual \\i\\ at time
\\t\\, and \\\lambda\_{it}\\ is the latent Poisson mean parameter for
individual \\i\\ at time \\t\\. The parameters \\a\\ and \\b\\ are the
shape parameters for the beta distribution that is used to model the
random effects. The RENB model allows for overdispersion in the count
data and accounts for unobserved heterogeneity across individuals by
including random effects in the model. This formulation follows the
approach described in the paper by Hausman, Hall, and Griliches (1984)
for modeling panel data with random effects.

The marginal mean and marginal variance of the RENB model are given by:
\$\$E\[y\_{it}\] =
\lambda\_{it}\frac{b}{a-1}=\mu_it=exp(X\_{it}\beta)\$\$

\$\$Var\[y\_{it}\] = \frac{a+b-1}{a-2}\mu\_{it} + \frac{a+b-1}{b(a-2)}
\mu\_{it}^2\$\$

Thus, the formulation of the model estimated here allows the use of the
estimated coefficients to directly compute the marginal mean.

Note that the RENB model is a panel data model, and the `group_var`
argument must be specified to indicate the grouping variable(s) for the
random effects. The model is estimated using maximum likelihood
estimation, and the optimization is performed using the
[`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html) package. The user
can specify the optimization method and maximum iterations.

## References

Hausman, Jerry A., Bronwyn H. Hall, and Zvi Griliches. "Econometric
models for count data with an application to the patents-R&D
relationship." Econometrica: Journal of the Econometric Society (1984):
909-938.

## Examples

``` r
# \donttest{
## RENB Model
data("washington_roads")
washington_roads$AADTover10k <- 
  ifelse(washington_roads$AADT > 10000, 1, 0) # create a dummy variable
renb.mod <- renb(Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k,
                                data=washington_roads,
                                offset = "lnlength",
                                group_var="ID",
                                method="nm",
                                max.iters = 1000)
summary(renb.mod)
#> Call:
#> renb(formula = Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k, 
#>     group_var = "ID", data = washington_roads, method = "nm", 
#>     max.iters = 1000, offset = "lnlength")
#> 
#> Model: renb | optimizer: NM
#> Iterations: 680 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)     -9.466      1.435   -6.60  4.2e-11 ***
#> lnaadt           0.978      0.167    5.85  4.9e-09 ***
#> speed50         -0.966      0.338   -2.86   0.0042 ** 
#> ShouldWidth04   -0.363      0.277   -1.31   0.1899    
#> AADTover10k     -0.946      0.546   -1.73   0.0831 .  
#> ln(a-1)          2.977      0.273   10.92  < 2e-16 ***
#> ln(b)            0.573      0.331    1.73   0.0836 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Offset: lnlength (coefficient fixed at 1)
#> 
#> Log-likelihood: -263.831 | AIC: 541.661 | BIC: 571.261 
#> N: 1501 | panels: 507
# }
```
