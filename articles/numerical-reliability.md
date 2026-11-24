# Numerical reliability: convergence, integration and difficult fits

An optimizer that stops and reports success has found a point where its
stopping rule was met. That is not the same as a reliable maximum
likelihood estimate with trustworthy standard errors. This vignette
shows what to inspect after a fit, how to recognize boundary estimates,
weak identification and scaling problems, which controls govern the
numerical integrals of the more complex models, and how to tell apart
three different reasons why a prediction cannot be computed.

``` r

library(flexCountReg)
data("washington_roads")
```

## What a fit records

The `model` element of a fitted
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
or
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
object is the result of the optimizer (from the maxLik package) with
some additions: `converged`, the return `code` and `message`, the number
of `iterations`, the log-likelihood `maximum`, the `gradient` and the
`hessian` at the estimates. A small helper collects them:

``` r

fit_status <- function(m) {
  f <- m$model
  ev <- eigen(-f$hessian, symmetric = TRUE, only.values = TRUE)$values
  data.frame(converged = isTRUE(f$converged), code = f$code,
             iterations = unname(f$iterations[1]),
             logLik = round(f$maximum, 3),
             max_abs_gradient = signif(max(abs(f$gradient)), 3),
             min_eigen = signif(min(ev), 3),
             condition = signif(max(ev) / min(ev), 3),
             message = trimws(f$message))
}
```

- `max_abs_gradient` should be close to 0 at a maximum. What counts as
  close depends on the scale of the parameters.
- `min_eigen` is the smallest eigenvalue of the negative Hessian. It
  must be positive at a maximum; a value near 0 means a direction in
  which the likelihood is nearly flat.
- `condition` is the ratio of the largest to the smallest eigenvalue. A
  very large value means that the parameters are on very different
  scales or nearly redundant.

## A fit-status table

The same NB2 model is fitted with four optimizers, and then with the
traffic volume entered in vehicles per day instead of on the log scale,
both as it is and divided by 10,000:

``` r

f <- Total_crashes ~ lnaadt + lnlength + speed50
fits <- lapply(c(NM = "NM", BFGS = "BFGS", NR = "NR", BHHH = "BHHH"),
               function(meth) countreg(f, data = washington_roads,
                                       family = "NB2", method = meth))
washington_roads$aadt <- exp(washington_roads$lnaadt)
washington_roads$aadt_10k <- washington_roads$aadt / 1e4
fits$NM_raw_aadt <- countreg(Total_crashes ~ aadt + speed50,
                             data = washington_roads, family = "NB2")
#> Warning: The numerical Hessian is uncertain: standard errors change by up to
#> 5.4% between differencing steps. Some parameters may be weakly identified (a
#> nearly flat log-likelihood); treat their standard errors with care.
fits$NM_aadt_10k <- countreg(Total_crashes ~ aadt_10k + speed50,
                             data = washington_roads, family = "NB2")
status_table <- do.call(rbind, lapply(fits, fit_status))
status_table
#>             converged code iterations    logLik max_abs_gradient min_eigen
#> NM               TRUE    0        446 -1084.942         7.50e-02      5.09
#> BFGS             TRUE    0         34 -1084.942         3.48e-04      5.09
#> NR               TRUE    1          3 -1084.942         1.06e-08      5.09
#> BHHH             TRUE    8          5 -1084.942         3.42e-02      5.09
#> NM_raw_aadt      TRUE    0        217 -1159.276         7.06e+02     34.70
#> NM_aadt_10k      TRUE    0        119 -1159.277         8.15e-02     33.40
#>             condition
#> NM           7.64e+03
#> BFGS         7.64e+03
#> NR           7.64e+03
#> BHHH         7.64e+03
#> NM_raw_aadt  7.41e+08
#> NM_aadt_10k  1.94e+01
#>                                                                         message
#> NM                                                       successful convergence
#> BFGS                                                     successful convergence
#> NR                                             gradient close to zero (gradtol)
#> BHHH        successive function values within relative tolerance limit (reltol)
#> NM_raw_aadt                                              successful convergence
#> NM_aadt_10k                                              successful convergence
```

The four optimizers agree on the log-likelihood and the coefficients
(`round(sapply(fits[1:4], coef), 4)` shows them), but they stop for
different reasons. Nelder-Mead (`"NM"`, the default) does not use
derivatives and stops when the simplex is small: its gradient is the
largest of the four. Newton-Raphson (`"NR"`) stops at a gradient of
essentially zero. A `code` of 0 means success for NM and BFGS, while NR
and BHHH report success with codes 1, 2 or 8 (`converged` already
accounts for this).

With volume in vehicles per day, Nelder-Mead reports convergence with a
largest gradient of 706: a coefficient per vehicle is tiny, and the
condition number of 7.41e+08 shows the poor scaling. Dividing the volume
by 10,000 brings the condition number down to 19.4. Rescaling a variable
changes only the units of its coefficient, not the model. More
iterations would not have fixed the scaling problem.

## Boundaries, weak identification and rank deficiency

### A parameter at the edge of its space

When the data are close to Poisson, the NB2 dispersion $`\alpha`$ is
close to 0, the edge of its space, and $`\ln\alpha`$ drifts towards
$`-\infty`$. Twenty samples of Poisson data, each fitted with NB2:

``` r

boundary_fit <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = rnorm(300))
  d$y <- rpois(300, exp(0.5 + 0.4 * d$x))
  warned <- FALSE
  m <- withCallingHandlers(
    countreg(y ~ x, data = d, family = "NB2", method = "BFGS"),
    warning = function(w) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    })
  co <- summary(m)$coefficients
  c(ln_alpha = co["ln(alpha)", 1], se = co["ln(alpha)", 2],
    converged = m$model$converged, warned = warned)
}
boundary <- as.data.frame(t(sapply(1:20, boundary_fit)))
boundary$at_start <- abs(boundary$ln_alpha - log(1e-3)) < 0.01
round(boundary, 3)
#>    ln_alpha     se converged warned at_start
#> 1    -3.451  1.454         1      0        0
#> 2    -3.705  1.809         1      0        0
#> 3    -2.680  0.724         1      0        0
#> 4   -11.233 60.270         1      0        0
#> 5    -6.181 20.512         1      0        0
#> 6    -6.908  6.702         1      0        1
#> 7    -3.722  1.703         1      0        0
#> 8    -6.908 17.422         1      0        1
#> 9    -3.977  2.366         1      0        0
#> 10   -3.515  1.705         1      0        0
#> 11   -3.368  1.295         1      0        0
#> 12  -10.750 53.640         1      0        0
#> 13   -6.908  7.080         1      0        1
#> 14   -2.627  0.754         1      0        0
#> 15  -11.319 48.727         1      0        0
#> 16   -6.908  7.839         1      0        1
#> 17   -3.164  1.183         1      0        0
#> 18   -6.908  6.052         1      0        1
#> 19   -3.447  1.462         1      0        0
#> 20   -6.908 11.390         1      0        1
```

Every fit reports convergence. The dispersion starts from a moment
estimate clipped at $`\alpha = 10^{-3}`$ ($`\ln\alpha = -6.908`$), and 6
of the 20 fits stopped within 0.01 of it: at the edge of the space the
likelihood is so flat that the gradient is below the optimizer’s
tolerance. In 0 fits the Hessian is not negative definite, and
flexCountReg warns that the standard error of `ln(alpha)` is `NA`. The
remaining standard errors are large. In every case the right conclusion
is that $`\alpha`$ is near 0: Poisson is adequate. Neither the estimate
of $`\ln\alpha`$ nor its Wald standard error is meaningful there, and a
likelihood-ratio test against Poisson needs the boundary-corrected
reference distribution that
[`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md)
and [`summary()`](https://rdrr.io/r/base/summary.html) use.

### Weak identification

A model can be identified in principle but not by the data at hand: the
likelihood has a long, nearly flat ridge. Standard errors computed from
the curvature can then be finite and even small while the estimates are
unstable. Small standard errors are not proof of identification.
flexCountReg warns in the cases it can detect (for example a constant
zero-inflation probability driven to 0, or an underreporting model whose
reporting and count coefficients are nearly collinear). The simulations
in
[`vignette("underreporting-and-identification", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/underreporting-and-identification.md)
show that refitting replicate samples is the most reliable check.

### Rank deficiency

A variable that is an exact linear combination of others is not
estimable. flexCountReg stops instead of reporting arbitrary
coefficients:

``` r

d <- data.frame(x = rnorm(100))
d$x2 <- 2 * d$x
d$y <- rpois(100, exp(0.3 + 0.2 * d$x))
res <- tryCatch(countreg(y ~ x + x2, data = d, family = "POISSON"),
                error = function(e) conditionMessage(e))
res
#> [1] "The design matrix of `formula` is rank deficient: x2 is a linear combination of other columns (aliased). Remove it."
```

## Starting values and optimizers

When different starts or optimizers end at materially different
log-likelihoods, the lower ones are not maxima, or the likelihood has
several local maxima. Refit from a deliberately poor start with
`start.vals` (one value per parameter, in the order of
[`coef()`](https://rdrr.io/r/stats/coef.html)):

``` r

bad_start <- c(0, 0, 0, 0, 3)
from_bad <- countreg(f, data = washington_roads, family = "NB2",
                     method = "BFGS", start.vals = bad_start)
rbind(default_start = fit_status(fits$BFGS)[, c("converged", "logLik")],
      poor_start = fit_status(from_bad)[, c("converged", "logLik")])
#>               converged    logLik
#> default_start      TRUE -1084.942
#> poor_start         TRUE -1084.942
```

Here the poor start reaches the same maximum. When it does not, report
the best fit, and treat the disagreement as a warning about the model or
the data rather than a nuisance to suppress. Finite-mixture and
latent-class models
([`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md),
[`countreg_lc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md))
have several local maxima as a rule; they run several starts themselves
and record them in `fit$model$starts`, together with the definiteness
and conditioning of the observed information (`fit$model$numerical`);
see
[`vignette("latent-class-regression", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/latent-class-regression.md).
Raising `max.iters` is only useful when the status shows that the
iteration limit stopped the optimizer.

## Numerical integration: outer draws and inner integrals

Some models need integrals that have no closed form. There are two
different kinds, controlled by different arguments:

| Model | Integral | Control |
|----|----|----|
| [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md) with PLN, PGE, PW or PLL | the mixing distribution of the family (inner) | `engine`: `"quadrature"` (default, adaptive; `ndraws` is not used) or `"halton"` (`ndraws` draws) |
| [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md) | the random parameters (outer) | `ndraws` Halton draws for each panel (`draws = "per_unit"`, default) or the same draws for every panel (`draws = "shared"`); `scrambled` for scrambled draws |
| [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md) with PLN, PGE, PW or PLL | the mixing distribution, inside each outer draw | `engine`; `mix_draws` only with `engine = "halton"` |

### Inner integrals: the Poisson-lognormal family

On a subsample of 400 road segments, the Poisson-lognormal model is
fitted with adaptive quadrature and with Halton draws of increasing
number:

``` r

set.seed(1)
sub <- washington_roads[sample(nrow(washington_roads), 400), ]
pln <- function(engine, ndraws) {
  m <- countreg(f, data = sub, family = "PLN", method = "BFGS",
                engine = engine, ndraws = ndraws)
  data.frame(engine = engine, ndraws = ndraws,
             logLik = m$model$maximum,
             ln_sigma = coef(m)[["ln(sigma)"]],
             se_ln_sigma = summary(m)$coefficients["ln(sigma)", 2])
}
quad <- rbind(pln("quadrature", 50), pln("quadrature", 1500))
halton <- do.call(rbind, lapply(c(25, 50, 100, 200, 400),
                                function(r) pln("halton", r)))
rbind(quad, halton)
#>       engine ndraws    logLik   ln_sigma se_ln_sigma
#> 1 quadrature     50 -294.4796 -0.5551128   0.2282893
#> 2 quadrature   1500 -294.4796 -0.5551128   0.2282893
#> 3     halton     25 -294.4376 -0.3531253   0.2340590
#> 4     halton     50 -294.4294 -0.4321822   0.2294705
#> 5     halton    100 -294.4224 -0.4745107   0.2281815
#> 6     halton    200 -294.4255 -0.5017463   0.2273157
#> 7     halton    400 -294.4333 -0.5200564   0.2268872
```

The two quadrature fits are identical: quadrature ignores `ndraws`, and
increasing it there only suggests a check that was not made. With Halton
draws, $`\ln\sigma`$ moves steadily toward the quadrature value as the
draws increase:

``` r

plot(halton$ndraws, halton$ln_sigma, type = "b", log = "x", pch = 19,
     ylim = range(halton$ln_sigma, quad$ln_sigma),
     xlab = "Halton draws (log scale)", ylab = "Estimate of ln(sigma)",
     main = "Poisson-lognormal: Halton draws against quadrature")
abline(h = quad$ln_sigma[1], lty = 2)
legend("bottomright", c("Halton", "adaptive quadrature"), lty = c(1, 2),
       pch = c(19, NA), bty = "n")
```

![](numerical-reliability_files/figure-html/pln-plot-1.png)

A simulated likelihood with too few draws is biased, here toward a
smaller $`\sigma`$, and its log-likelihood is not comparable with the
quadrature one. The Halton engine is kept to reproduce fits from earlier
versions.

### Outer draws: random parameters

For a random-parameter model, compare the estimates as the number of
outer draws increases (on the same subsample):

``` r

rp_fit <- function(r) {
  m <- countreg_rp(Total_crashes ~ lnaadt + lnlength,
                   rpar_formula = ~ -1 + speed50, data = sub,
                   family = "NB2", ndraws = r, method = "BFGS")
  co <- summary(m)$coefficients
  data.frame(ndraws = r, converged = m$model$converged,
             logLik = round(m$model$maximum, 3),
             mean = co["speed50:Mean", 1], sd = co["speed50:St.Dev", 1],
             se_sd = co["speed50:St.Dev", 2])
}
rp_draws <- do.call(rbind, lapply(c(50, 100, 200), rp_fit))
rp_draws
#>   ndraws converged   logLik       mean        sd     se_sd
#> 1     50      TRUE -294.177 -0.5708183 0.5124562 0.4916367
#> 2    100      TRUE -294.278 -0.5321927 0.4309691 0.5723279
#> 3    200      TRUE -294.293 -0.5278457 0.4168919 0.5763524
```

The estimated standard deviation of the random parameter changes with
the draws by an amount that is small relative to its standard error
here, but it has not settled. Choose the number of draws by this kind of
comparison for the model at hand, and report it. Larger models and
studies belong in scripts outside of a vignette build.

## Three reasons a prediction is not a number

``` r

rp_ln <- countreg_rp(Total_crashes ~ lnaadt + lnlength,
                     rpar_formula = ~ -1 + speed50,
                     rpardists = c(speed50 = "ln"), data = sub,
                     family = "NB2", ndraws = 100, method = "BFGS")
nd <- data.frame(lnaadt = log(5000), lnlength = 0, speed50 = c(0, 1))
predict(rp_ln, newdata = nd, method = "Exact")
#> Warning: 1 prediction(s) are infinite: the expected value does not exist for a
#> lognormal random parameter with a positive x, or a gamma random parameter with
#> x * sigma^2 / mu >= 1. See ?countreg_rp.
#> [1] 1.574824      Inf
```

1.  **The expectation does not exist.** With a lognormal random
    parameter $`\beta`$, the expected count $`E[\exp(x\beta)]`$ is
    infinite for every $`x > 0`$: the lognormal distribution has no
    moment-generating function for positive arguments. flexCountReg
    returns `Inf` with a warning, for both `"Exact"` and `"Simulated"`
    predictions. This is a property of the model, not a numerical
    failure. (Here the fit also put the lognormal parameter at the edge
    of its space, with `NA` standard errors: a lognormal effect is
    positive, while the effect of `speed50` is negative.)

2.  **The expectation exists but exceeds the floating-point range.**
    With a normal random parameter,
    $`E[\exp(x\beta)] = \exp(x\mu + x^2\sigma^2/2)`$ is finite for every
    $`x`$, but it overflows for extreme $`x`$. The value below is far
    outside the data (`speed50` is 0 or 1) and is used only to show the
    message:

    ``` r

    rp_n <- countreg_rp(Total_crashes ~ lnaadt + lnlength,
                        rpar_formula = ~ -1 + speed50, data = sub,
                        family = "NB2", ndraws = 100, method = "BFGS")
    predict(rp_n, newdata = data.frame(lnaadt = log(5000), lnlength = 0,
                                       speed50 = c(1, 100)),
            method = "Exact")
    #> Warning: 1 Exact prediction(s) exist but exceed the largest representable
    #> number (overflow); they are Inf.
    #> [1] 1.140991      Inf
    ```

3.  **A valid quantity was not approximated.** Tolerance intervals need
    the predictive distribution’s upper tail. With a deliberately small
    `max_count`, the tail is not resolved, and the interval is reported
    as unavailable with a status rather than truncated silently:

    ``` r

    nb <- fits$BFGS
    at <- data.frame(lnaadt = log(c(2000, 20000)), lnlength = log(2),
                     speed50 = 0)
    tolerance_interval(nb, newdata = at, max_count = 15)
    #> Tolerance intervals: content 0.95, confidence 0.95, two-sided, pointwise
    #> Method: joint (normal)
    #>          fit lower upper min_content
    #> 1  0.9836129    NA    NA          NA
    #> 2 11.8924520    NA    NA          NA
    #>                                          status
    #> 1 upper tail not resolved within max_count = 15
    #> 2 upper tail not resolved within max_count = 15
    tolerance_interval(nb, newdata = at)
    #> Tolerance intervals: content 0.95, confidence 0.95, two-sided, pointwise
    #> Method: joint (normal)
    #>          fit lower upper min_content status
    #> 1  0.9836129     0     4   0.9620973     ok
    #> 2 11.8924520     0    40   0.9512015     ok
    ```

## Reporting an unresolved problem

When a fit remains doubtful after these checks, a useful report
contains: the code and data (or a simulated dataset with the same
problem), the seed, the fit-status table, the warnings, and the
versions:

``` r

c(R = R.version.string,
  flexCountReg = as.character(utils::packageVersion("flexCountReg")),
  maxLik = as.character(utils::packageVersion("maxLik")))
#>                              R                   flexCountReg 
#> "R version 4.6.2 (2026-10-29)"                        "0.1.2" 
#>                         maxLik 
#>                       "1.6.10"
```

Issues can be reported at
<https://github.com/jwood-iastate/flexCountReg/issues>.

## Summary

- Check the gradient, the curvature and the condition of the Hessian,
  not only `converged`.
- Rescale variables that make the condition number very large.
- Treat estimates that do not move from their starting values, `NA`
  standard errors and nearly flat directions as information about the
  model, not as noise to suppress.
- Know which control governs which integral, and compare results as it
  changes.
- An infinite prediction can be a property of the model, an overflow, or
  an approximation that was not achieved; the warnings and statuses say
  which.
