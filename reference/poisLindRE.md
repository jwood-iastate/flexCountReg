# Function for estimating a Random Effects Poisson-Lindley regression model

Function for estimating a Random Effects Poisson-Lindley regression
model

## Usage

``` r
poisLind.re(
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

  the grouping variable(s) indicating random effects (e.g., individual
  ID).

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

The function `poisLindRE` is similar to the `poisLind` function, but it
includes an additional argument `group_var` that specifies the grouping
variable for the random effects. The function estimates a Random Effects
Poisson-Lindley regression model using maximum likelihood. It is similar
to `poisLind`, but includes additional terms to account for the random
effects.

The Random Effects Poisson-Lindley model is useful for panel data and
assumes that the random effects follow a gamma distribution. The PDF is
\$\$ f(y\_{it}\|\mu\_{it},\theta)=\frac{\theta^2}{\theta+1}
\prod\_{t=1}^{n_i}\frac{\left(\mu\_{it}\frac{\theta(\theta+1)}
{\theta+2}\right)^{y\_{it}}}{y\_{it}!} \cdot \frac{
\left(\sum\_{t=1}^{n_i}y\_{it}\right)!
\left(\sum\_{t=1}^{n_i}\mu\_{it}\frac{\theta(\theta+1)}{\theta+2} +
\theta + \sum\_{t=1}^{n_i}y\_{it} + 1\right) }{
\left(\sum\_{t=1}^{n_i}\mu\_{it}\frac{\theta(\theta+1)}{\theta+2} +
\theta\right)^{\sum\_{t=1}^{n_i}y\_{it}+2} } \$\$

The log-likelihood function is: \$\$ LL = 2\log(\theta) -
\log(\theta+1) + \sum\_{t=1}^{n_i} y\_{it}\log(\mu\_{it}) +
\sum\_{t=1}^{n_i} y\_{it}\log\\\left( \frac{\theta(\theta+1)}{\theta+2}
\right) - \sum\_{t=1}^{n_i}\log(y\_{it}!) + \log\\\left(
\left(\sum\_{t=1}^{n_i}y\_{it}\right)! \right) + \log\\\left(
\sum\_{t=1}^{n_i}\mu\_{it}\frac{\theta(\theta+1)}{\theta+2} + \theta +
\sum\_{t=1}^{n_i}y\_{it} + 1 \right) - \left(\sum\_{t=1}^{n_i}y\_{it} +
2\right) \log\\\left(
\sum\_{t=1}^{n_i}\mu\_{it}\frac{\theta(\theta+1)}{\theta+2} + \theta
\right) \$\$

The mean and variance are: \$\$\mu\_{it}=\exp(X\_{it} \beta)\$\$ \$\$
V(\mu\_{it})=\mu\_{it}+ \left(1-\frac{2}{(\theta+2)^2}\right)\mu\_{it}^2
\$\$

## Examples

``` r
# \donttest{
data("washington_roads")
washington_roads$AADTover10k <-
  ifelse(washington_roads$AADT > 10000, 1, 0)

poislind.mod <- poisLind.re(
  Animal ~ lnaadt + lnlength + speed50 +
    ShouldWidth04 + AADTover10k,
  data      = washington_roads,
  group_var = "ID",
  method    = "NM",
  max.iters = 1000
)
summary(poislind.mod)
#> Call:
#>  Animal ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k 
#> 
#>  Method:  poisLindRE 
#> Iterations:  712 
#> Convergence:  successful convergence  
#> Log-likelihood:  -263.098 
#> 
#> Parameter Estimates:
#> # A tibble: 8 × 7
#>   parameter           coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>               <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)        -9.61        1.43     -6.71     0        -12.4       -6.8  
#> 2 lnaadt              1.03        0.168     6.15     0          0.705      1.36 
#> 3 lnlength            1.46        0.236     6.19     0          0.996      1.92 
#> 4 speed50            -0.856       0.337    -2.54     0.011     -1.52      -0.195
#> 5 ShouldWidth04      -0.434       0.292    -1.49     0.137     -1.01       0.138
#> 6 AADTover10k        -0.758       0.526    -1.44     0.15      -1.79       0.273
#> 7 ln(theta)           3.50        2.13      1.64     0.1       -0.671      7.68 
#> 8 Offset (Offset va…  1          NA        NA       NA         NA         NA    
# }
```
