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
#> Warning: NaNs produced
summary(poislind.mod)
#> Call:
#>  Animal ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k 
#> 
#>  Method:  poisLindRE 
#> Iterations:  1002 
#> Convergence:  iteration limit exceeded  
#> Log-likelihood:  91004.96 
#> 
#> Parameter Estimates:
#> # A tibble: 7 × 7
#>   parameter       coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>           <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)    -5.55        0.151   -36.6          0     -5.84      -5.25 
#> 2 lnaadt          0.175       0.018     9.44         0      0.138      0.211
#> 3 lnlength        8.67        0.262    33.1          0      8.16       9.19 
#> 4 speed50         4.26        0.131    32.5          0      4.00       4.52 
#> 5 ShouldWidth04  -2.92        0.262   -11.1          0     -3.43      -2.40 
#> 6 AADTover10k    11.8         0.262    45.2          0     11.3       12.4  
#> 7 ln(theta)     -46.6        NA        NA           NA     NA         NA    
# }
```
