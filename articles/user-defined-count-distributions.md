# Fitting your own count distribution

flexCountReg has many count distributions built in (see
[`?count_families`](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md)).
When the one you need is not among them, you can define it yourself with
[`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md)
and use it with
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
and
[`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md),
exactly like a built-in family. You write one function, the log
probability mass function; everything else (fitting, predictions,
residuals, rootograms, simulation) works from it.

This article builds a distribution that is not in the package, checks
it, fits it, and compares it with the built-in negative binomial.

``` r

library(flexCountReg)
data("washington_roads")
```

## 1. The distribution

As an example, take a **one-inflated negative binomial**: a negative
binomial (NB2) with an extra probability $`\omega`$ of a count of
exactly one,

``` math
P(Y = y) = (1 - \omega)\, f_{NB}(y \mid \mu, \alpha) + \omega\, [y = 1],
```

where $`f_{NB}`$ has mean $`\mu`$ and variance $`\mu + \alpha\mu^2`$.
Excess ones appear in data where a first event changes what follows (for
example, a site that is treated after its first crash, or repeat-offense
counts). The distribution has two parameters besides $`\mu`$: the
dispersion $`\alpha > 0`$ and the probability $`\omega \in (0, 1)`$. Its
mean is $`(1 - \omega)\mu + \omega`$, so $`\mu`$ is *not* the mean.

## 2. Define the family

A family needs a name, its parameters and the log PMF. Each parameter
has a link, which sets the scale of its linear predictor: `"log"` for a
positive parameter, `"logit"` for a probability and `"identity"` for any
real value. `start` is the starting value on that scale.

``` r

oinb <- new_count_family(
  name = "OINB",
  params = list(alpha = count_param("log", start = -1),
                omega = count_param("logit", start = -3)),
  logpmf = function(y, mu, par) {
    lnb <- dnbinom(y, size = 1 / par$alpha, mu = mu, log = TRUE)
    ifelse(y == 1,
           log((1 - par$omega) * exp(lnb) + par$omega),
           log1p(-par$omega) + lnb)
  },
  mean_is_mu = FALSE,
  mean = function(mu, par) (1 - par$omega) * mu + par$omega,
  variance = function(mu, par) {
    m <- (1 - par$omega) * mu + par$omega
    ey2 <- (1 - par$omega) * (mu + par$alpha * mu^2 + mu^2) + par$omega
    ey2 - m^2
  })
oinb
#> Count family "OINB" with parameter(s) alpha (log link), omega (logit link)
#>   mean: user function
#>   supplied: logpmf, mean, variance
```

`logpmf` receives vectors `y` and `mu` and a list `par` with one vector
per parameter, all of the same length, and must return one
log-probability per element. Only `logpmf` is required:

- without `mean` (and with `mean_is_mu = FALSE`) and without `variance`,
  the moments are computed by summing the PMF;
- without `cdf`, the CDF is the cumulative sum of the PMF;
- without `rand`, random draws are made by inverting that CDF.

Giving them is faster, and
[`check_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/check_count_family.md)
checks them against the PMF. Parameter values outside their space should
give `-Inf` or `NaN` (errors are caught too);
`valid = function(mu, par)` can state the parameter space explicitly.

## 3. Check it before fitting

[`check_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/check_count_family.md)
evaluates the family over a grid of means and parameter values: the
probabilities must sum to one, the mean and variance functions must
agree with the PMF, the function must be vectorized, and missing or
invalid inputs must not cause errors.

``` r

chk <- check_count_family(oinb, mu = c(0.2, 2, 20),
                          par = list(alpha = c(0.2, 2),
                                     omega = c(0.05, 0.5)))
#>                        check   ok                                    detail
#>                PMF sums to 1 TRUE               largest |sum - 1| = 2.1e-13
#>       mean() matches the PMF TRUE          largest relative error = 9.6e-11
#>   variance() matches the PMF TRUE          largest relative error = 2.1e-09
#>            vectorized over y TRUE                                        ok
#>              NA mean handled TRUE                                  no error
#>  invalid mean gives no error TRUE no error (the row is given probability 0)
```

A family that fails a check gives a warning naming the failed checks.
For example, a variance that scales the negative binomial variance by
$`1 - \omega`$ and ignores the extra ones:

``` r

wrong <- oinb
wrong$variance <- function(mu, par) (1 - par$omega) * (mu + par$alpha * mu^2)
chk <- check_count_family(wrong, mu = c(0.2, 2, 20),
                          par = list(alpha = c(0.2, 2),
                                     omega = c(0.05, 0.5)))
#>                        check    ok                                    detail
#>                PMF sums to 1  TRUE               largest |sum - 1| = 2.1e-13
#>       mean() matches the PMF  TRUE          largest relative error = 9.6e-11
#>   variance() matches the PMF FALSE             largest relative error = 0.64
#>            vectorized over y  TRUE                                        ok
#>              NA mean handled  TRUE                                  no error
#>  invalid mean gives no error  TRUE no error (the row is given probability 0)
#> Warning: Failed checks: variance() matches the PMF
```

## 4. Fit it

Pass the family object as `family`. The coefficients of the distribution
parameters are labelled by their link: `ln(alpha)` and `logit(omega)`.

``` r

f <- Total_crashes ~ lnaadt + lnlength + speed50
m_oinb <- countreg(f, data = washington_roads, family = oinb,
                   method = "BFGS")
summary(m_oinb)
#> Call:
#> countreg(formula = f, data = washington_roads, family = oinb, 
#>     method = "BFGS")
#> 
#> Model: countreg | family: OINB | optimizer: BFGS
#> Iterations: 40 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   -8.8891     0.5746  -15.47  < 2e-16 ***
#> lnaadt         1.0964     0.0680   16.12  < 2e-16 ***
#> lnlength       0.7704     0.0748   10.30  < 2e-16 ***
#> speed50       -0.5401     0.1095   -4.93  8.1e-07 ***
#> ln(alpha)     -1.0264     0.2534   -4.05  5.1e-05 ***
#> logit(omega)  -5.5742     2.8808   -1.93    0.053 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Log-likelihood: -1084.923 | AIC: 2181.847 | BIC: 2213.73 
#> N: 1501
```

Compare it with the built-in negative binomial, which is the special
case $`\omega = 0`$:

``` r

m_nb2 <- countreg(f, data = washington_roads, family = "NB2",
                  method = "BFGS")
c(NB2 = AIC(m_nb2), OINB = AIC(m_oinb))
#>      NB2     OINB 
#> 2179.884 2181.847
```

These road segments have no excess ones. The estimate of $`\omega`$ is
0.0038, at the boundary of its space, so its standard error is not
available (the warning above), and the AIC prefers the negative
binomial. The frequency table in Section 5 shows why: there are fewer
ones than the negative binomial predicts, not more.

When the extra ones are there, the family finds them. Here 20% of the
counts of a negative binomial sample are replaced by ones:

``` r

set.seed(42)
n <- 2000
sim <- data.frame(x = rnorm(n))
sim$y <- rnbinom(n, size = 2, mu = exp(0.8 + 0.5 * sim$x))
sim$y[runif(n) < 0.2] <- 1
s_oinb <- countreg(y ~ x, data = sim, family = oinb, method = "BFGS")
s_nb2 <- countreg(y ~ x, data = sim, family = "NB2", method = "BFGS")
c(NB2 = AIC(s_nb2), OINB = AIC(s_oinb))
#>      NB2     OINB 
#> 7538.036 7297.448
# omega (true 0.2) and alpha (true 0.5)
c(omega = plogis(coef(s_oinb)[["logit(omega)"]]),
  alpha = exp(coef(s_oinb)[["ln(alpha)"]]))
#>     omega     alpha 
#> 0.2066053 0.5075688
```

As with the built-in families, each distribution parameter can have its
own regression: `dis_param_formula_1` models the first parameter
($`\alpha`$) and `dis_param_formula_2` the second ($`\omega`$), on their
link scales.

``` r

m_oinb2 <- countreg(f, data = washington_roads, family = oinb,
                    dis_param_formula_2 = ~ speed50, method = "BFGS")
coef(m_oinb2)[c("logit(omega):(Intercept)", "logit(omega):speed50")]
#> logit(omega):(Intercept)     logit(omega):speed50 
#>                -4.749657                -3.187266
```

## 5. Everything else works from the PMF

Predictions of the mean use the family’s `mean` function; probabilities,
quantile residuals, rootograms and simulation use the PMF.

``` r

nd <- washington_roads[1:3, ]
predict(m_oinb, newdata = nd)                      # E[Y]
#> [1] 0.7789732 0.7085595 1.0441579
round(predict(m_oinb, newdata = nd, type = "prob")[, 1:5], 3)
#>          0     1     2     3     4
#> [1,] 0.501 0.309 0.126 0.044 0.014
#> [2,] 0.530 0.303 0.115 0.037 0.011
#> [3,] 0.410 0.316 0.161 0.070 0.028
head(residuals(m_oinb, type = "quantile"))
#> [1] -1.2686338  1.2134864  1.1183553 -0.1828045 -1.3905276  0.5135716
sims <- simulate(m_oinb, nsim = 2, seed = 1)
```

Observed and expected frequencies of the small counts show whether extra
ones are needed. For the road data both models predict more ones than
were observed:

``` r

expected <- function(m) {
  colSums(predict(m, newdata = washington_roads, type = "prob")[, 1:5])
}
rbind(observed = tabulate(washington_roads$Total_crashes + 1, 5),
      NB2 = round(expected(m_nb2), 1),
      OINB = round(expected(m_oinb), 1))
#>               0     1    2    3    4
#> observed 1101.0 242.0 91.0 30.0 23.0
#> NB2      1095.1 254.2 84.0 35.2 16.3
#> OINB     1094.5 256.0 82.9 34.9 16.2
```

The family also works in random-parameter models and finite mixtures:

``` r

m_rp <- countreg_rp(Total_crashes ~ lnaadt + lnlength,
                    rpar_formula = ~ -1 + speed50,
                    data = washington_roads, family = oinb,
                    group_var = "ID", ndraws = 50)
coef(m_rp)
#>    (Intercept)         lnaadt       lnlength   speed50:Mean speed50:St.Dev 
#>     -8.8709961      1.0940735      0.7680520     -0.7064730      0.5809304 
#>      ln(alpha)   logit(omega) 
#>     -1.2604041     -5.5365437
```

## 6. Saving and sharing a fitted model

The fitted model stores the family, so it can be saved and used in
another R session, on another computer, or in parallel workers, without
defining the family again:

``` r

path <- tempfile(fileext = ".rds")
saveRDS(m_oinb, path)
m_again <- readRDS(path)
predict(m_again, newdata = nd)
#> [1] 0.7789732 0.7085595 1.0441579
```

Two families with the same `name` but different definitions are kept
apart: each is identified by its name and a hash of its definition.

## 7. What a user family can and cannot do

| Feature | Available |
|:---|:---|
| [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md), with distribution-parameter formulas | yes |
| zero-inflated and hurdle models (`zero_formula`) | yes |
| [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md) (random parameters, panels) | yes |
| [`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md) and [`countreg_lc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md) | yes, slower than built-in NB and Poisson mixtures |
| predictions, probabilities, residuals, rootograms, [`simulate()`](https://rdrr.io/r/stats/simulate.html) | yes |
| robust and cluster-robust standard errors, `vcov(type = "cluster")` | yes |
| underreporting (`underreport_formula`) | no: a user family is not known to be closed under binomial thinning |
| more than two distribution parameters | not yet |

The gradient for a user family is computed by finite differences of the
log PMF in the mean and in each distribution parameter. That is accurate
and needs only two extra PMF evaluations per parameter, but a
closed-form built-in family is faster. If a distribution is widely
useful, it can be added to the package as a built-in family;
[`vignette("extending-count-distributions")`](https://jwood-iastate.github.io/flexCountReg/articles/extending-count-distributions.md)
explains how.
