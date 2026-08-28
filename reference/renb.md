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
  print.level = 0,
  bootstraps = NULL,
  offset = NULL
)
```

## Arguments

- formula:

  an R formula.

- group_var:

  the grouping variable(s) for the random effects (e.g., individual ID
  or other panel ID variables).

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

- print.level:

  Integer specifying the verbosity of output during optimization.

- bootstraps:

  Optional integer specifying the number of bootstrap samples to be used
  for estimating standard errors. If not specified, no bootstrapping is
  performed.

- offset:

  an optional offset term provided as a string.

## Value

An object of class `countreg` which is a list with the following
components:

- model: the fitted model object.

- data: the data frame used to fit the model.

- call: the matched call.

- formula: the formula used to fit the model.

## Details

This function estimates a random effects negative binomial (RENB)
regression model. This model is based on the NB-1 model. The PDF for the
RENB is: \$\$f(y\_{it}\|\lambda\_{it}, a, b) = \frac{\Gamma(a+b)
\Gamma(a + \sum\_{t = 1}^{n_i} \\lambda\_{it}) \Gamma(b +
\sum\_{t=1}^{n_i}y\_{it})} {\Gamma(a) \Gamma(b) \Gamma(a + b +
\sum\_{t=1}^{n_i}\lambda\_{it} + \sum\_{t=1}^{n_i}y\_{it})}
\prod\_{t=1}^{n_i}
\frac{\Gamma(\lambda\_{it}+y\_{it})}{\Gamma(\lambda\_{it})\Gamma(y\_{it})}\$\$

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
models for count data with an application to the patents–R&D
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
#> Error in exp(unlist(fit$estimate[(length(fit$estimate) - 1)])): non-numeric argument to mathematical function
summary(renb.mod)
#> Error: object 'renb.mod' not found
# }
```
