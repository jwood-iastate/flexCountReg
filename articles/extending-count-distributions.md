# Adding a count distribution to flexCountReg

This article is for contributors. It explains what a new count
distribution has to provide for it to work throughout flexCountReg, by
tracing the existing NB2 family through the source code. Code that uses
`flexCountReg:::` reaches *internal* objects: they are shown here to
explain the implementation, may change between versions, and are not
part of the user-facing interface. Users only see the `family` argument
of
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
and
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
and the exported distribution functions
([`ddpo()`](https://jwood-iastate.github.io/flexCountReg/reference/DoublePoisson.md),
[`pdpo()`](https://jwood-iastate.github.io/flexCountReg/reference/DoublePoisson.md)
and so on).

``` r

library(flexCountReg)
```

## 1. Define the distribution first

Before writing code, state for the new family:

- its support (here, the counts $`0, 1, 2, \dots`$);
- its natural parameters and their spaces;
- how the regression enters: in flexCountReg the mean model is
  $`\mu = \exp(X\beta)`$ (times an offset), and each distribution
  parameter is $`\exp(W\delta)`$ of its own linear predictor, or a
  logistic function of it for a probability (labels such as `ln(alpha)`
  and `logit(rho)` show the scale);
- whether $`\mu`$ is the actual mean. It is for NB2. For the
  Poisson-lognormal family the mean is $`\mu \exp(\sigma^2/2)`$; for the
  generalized Poisson with negative dispersion $`\mu`$ is a location
  parameter; for the double Poisson, hyper-Poisson, gamma count and
  discrete Weibull families flexCountReg solves numerically for the
  natural parameter that makes the mean equal to $`\mu`$.

These choices decide how coefficients are interpreted, so they belong in
the documentation
([`?count_families`](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md))
as well as in the code.

## 2. Source-file map

| Piece | Where | NB2 | A family added in 0.1.2 (DPO) |
|----|----|----|----|
| Public d/p/q/r functions | `R/<name>.R` | base R [`dnbinom()`](https://rdrr.io/r/stats/NegBinomial.html) etc. | `R/dpo.R`: [`ddpo()`](https://jwood-iastate.github.io/flexCountReg/reference/DoublePoisson.md), [`pdpo()`](https://jwood-iastate.github.io/flexCountReg/reference/DoublePoisson.md), [`qdpo()`](https://jwood-iastate.github.io/flexCountReg/reference/DoublePoisson.md), [`rdpo()`](https://jwood-iastate.github.io/flexCountReg/reference/DoublePoisson.md) |
| Shared argument handling | `R/count_dist_utils.R` | – | `.fcr_d()`, `.fcr_p()`, `.fcr_q()`, `.fcr_r()` |
| Compiled code (if any) | `src/*.cpp`, then [`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html) | – | `src/meanParam.cpp` (mean equation) |
| Regression family entry | `R/families.R` (`.fcr_family_registry`) | `NB2` entry | `DPO` entry |
| Mean-parameterization helpers | `R/families_mp.R` | – | `.tab_mp()`, `.dpo_mp_rand()` |
| Underreporting (opt-in) | `R/families.R`: `.fcr_thinning_families`, `.fcr_thin_params()` | listed | not closed under thinning |
| Likelihood-ratio nesting (opt-in) | `R/nesting.R`: `.fcr_family_nesting` | nests Poisson (boundary) | nests Poisson (interior) |
| Overview documentation | `R/count_families.R` ([`?count_families`](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md)) | row in each table | row in each table |
| Reference index | `_pkgdown.yml` “Count distributions” (also lists them in the README) | – | `DoublePoisson` |
| Tests | `tests/testthat/` | several | `test-new-distributions.R`, `test-new-families-regression.R`, `test-count-families-doc.R` |
| News | `NEWS.md` | – | entry for 0.1.2 |

## 3. The regression family entry

Every family that
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
and
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
can use has one entry in the internal list `.fcr_family_registry`, built
with `.fcr_family()`. The NB2 entry:

``` r

nb2 <- flexCountReg:::.fcr_family_registry$NB2
names(nb2)
#>  [1] "params"   "dens"     "variance" "rand"     "mean_adj" "start"   
#>  [7] "mean"     "cdf"      "start_rp" "valid"
nb2$params
#> [[1]]
#> [1] "ln(alpha)"
#> 
#> [[2]]
#> NULL
names(formals(nb2$dens))
#> [1] "y"         "predicted" "alpha"     "sigma"     "..."       "log"
body(nb2$dens)
#> .fcr_quiet_dens(stats::dnbinom(y, size = 1/alpha, mu = predicted, 
#>     log = log), log)
```

The fields (documented at the top of `R/families.R`):

| Field | Meaning |
|----|----|
| `params` | labels of the (up to two) distribution parameters, which are also the coefficient names; `NULL` for none |
| `dens` | `function(y, predicted, alpha, sigma, ..., log)`: $`P(Y = y)`$ given the mean parameter `predicted` and the distribution parameters on their natural scale; returns probability 0 (or `-Inf`) outside the parameter space, without warnings, because the optimizer tries such values |
| `variance` | variance at mean parameter $`\mu`$ (Pearson residuals, documentation tests) |
| `rand` | random draws ([`simulate()`](https://rdrr.io/r/stats/simulate.html), parametric bootstrap) |
| `cdf` | optional log CDF or log survival function, computed directly (quantile residuals, intervals); without it the PMF is summed |
| `mean`, `mean_adj` | the mean when it is not $`\mu`$ |
| `start`, `start_rp` | starting values of the distribution parameters (on the linear-predictor scale) for [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md) and [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md) |

The registry density must agree with an independent calculation:

``` r

y <- 0:6
all.equal(nb2$dens(y, predicted = 2.5, alpha = 0.4, sigma = NULL),
          dnbinom(y, size = 1 / 0.4, mu = 2.5))
#> [1] TRUE
```

The hooks present for each family, read from the registry, together with
the opt-in features:

``` r

reg <- flexCountReg:::.fcr_family_registry
nest <- flexCountReg:::.fcr_family_nesting
label <- function(x) if (is.null(x)) "--" else x
capabilities <- data.frame(
  family = names(reg),
  param_1 = vapply(reg, function(f) label(f$params[[1]]), ""),
  param_2 = vapply(reg, function(f) label(f$params[[2]]), ""),
  variance = vapply(reg, function(f) !is.null(f$variance), TRUE),
  cdf = vapply(reg, function(f) !is.null(f$cdf), TRUE),
  own_mean = vapply(reg, function(f) !is.null(f$mean), TRUE),
  thinning = names(reg) %in% flexCountReg:::.fcr_thinning_families,
  lr_vs_poisson = names(reg) %in% nest$full[nest$restricted == "POISSON"],
  row.names = NULL)
capabilities
#>        family     param_1   param_2 variance   cdf own_mean thinning
#> 1     POISSON          --        --     TRUE  TRUE    FALSE     TRUE
#> 2         GP1   ln(phi+1)        --     TRUE FALSE     TRUE    FALSE
#> 3         GP2 ln(alpha+1)        --     TRUE FALSE     TRUE    FALSE
#> 4         NB1   ln(alpha)        --     TRUE  TRUE    FALSE     TRUE
#> 5         NB2   ln(alpha)        --     TRUE  TRUE    FALSE     TRUE
#> 6         NBP   ln(alpha)     ln(p)     TRUE  TRUE    FALSE     TRUE
#> 7         PLN   ln(sigma)        --     TRUE FALSE    FALSE     TRUE
#> 8         PGE   ln(shape)        --     TRUE FALSE    FALSE     TRUE
#> 9        PIG1     ln(eta)        --     TRUE FALSE    FALSE     TRUE
#> 10       PIG2     ln(eta)        --     TRUE FALSE    FALSE     TRUE
#> 11        PIG     ln(eta)        --     TRUE FALSE    FALSE     TRUE
#> 12         PL   ln(theta)        --     TRUE FALSE    FALSE     TRUE
#> 13        PLG   ln(theta) ln(alpha)     TRUE FALSE    FALSE     TRUE
#> 14        PLL   ln(theta) ln(sigma)     TRUE FALSE    FALSE     TRUE
#> 15         PW   ln(alpha)        --     TRUE FALSE    FALSE     TRUE
#> 16         SI       gamma ln(sigma)     TRUE FALSE    FALSE     TRUE
#> 17         GW       ln(k) ln(rho-1)     TRUE FALSE    FALSE    FALSE
#> 18        COM      ln(nu)        --     TRUE FALSE    FALSE    FALSE
#> 19         DW   ln(shape)        --     TRUE  TRUE    FALSE    FALSE
#> 20        DPO     ln(phi)        --     TRUE  TRUE    FALSE    FALSE
#> 21  DELAPORTE  logit(rho)     ln(k)     TRUE  TRUE    FALSE     TRUE
#> 22    HERMITE  logit(rho)        --     TRUE  TRUE    FALSE     TRUE
#> 23   HERMITE3  logit(rho)        --     TRUE  TRUE    FALSE    FALSE
#> 24   HERMITE4  logit(rho)        --     TRUE  TRUE    FALSE    FALSE
#> 25  HYPERPOIS   ln(gamma)        --     TRUE  TRUE    FALSE    FALSE
#> 26 GAMMACOUNT   ln(shape)        --     TRUE  TRUE    FALSE    FALSE
#> 27         PA     ln(phi)        --     TRUE  TRUE    FALSE     TRUE
#> 28       NBGE   ln(shape) ln(alpha)     TRUE  TRUE    FALSE     TRUE
#>    lr_vs_poisson
#> 1          FALSE
#> 2           TRUE
#> 3           TRUE
#> 4           TRUE
#> 5           TRUE
#> 6          FALSE
#> 7           TRUE
#> 8          FALSE
#> 9           TRUE
#> 10          TRUE
#> 11          TRUE
#> 12         FALSE
#> 13         FALSE
#> 14         FALSE
#> 15         FALSE
#> 16         FALSE
#> 17         FALSE
#> 18          TRUE
#> 19         FALSE
#> 20          TRUE
#> 21         FALSE
#> 22         FALSE
#> 23         FALSE
#> 24         FALSE
#> 25          TRUE
#> 26          TRUE
#> 27         FALSE
#> 28         FALSE
```

## 4. Public distribution functions

A family is more than a likelihood: users need its d/p/q/r functions.
The families added in 0.1.2 implement four small internal functions and
pass them to shared wrappers in `R/count_dist_utils.R`, which give every
distribution the argument conventions of base R: recycling, zero-length
inputs, `NA` and `NaN`, `NaN` with a warning for invalid parameters,
probability 0 with a warning for non-integer `x`, `lower.tail` and
`log.p`. For the double Poisson:

``` r

body(ddpo)
#> {
#>     .fcr_d(x, list(lambda = lambda, phi = phi), .dpo_valid, .dpo_logpmf, 
#>         log)
#> }
body(qdpo)
#> {
#>     .fcr_q(p, list(lambda = lambda, phi = phi), .dpo_valid, .dpo_logcdf, 
#>         .dpo_qstart, lower.tail, log.p)
#> }
```

- `.dpo_valid()` says which parameter values are in the parameter space;
- `.dpo_logpmf()` returns $`\log P(Y = y)`$, computed on the log scale;
- `.dpo_logcdf()` returns $`\log P(Y \le y)`$ or $`\log P(Y > y)`$, each
  computed directly, so that tiny upper-tail probabilities keep their
  relative precision (never as $`1 - P(Y \le y)`$);
- `.dpo_qstart()` gives a rough starting point for the quantile search.

The contract can be checked numerically. For the Delaporte distribution:

``` r

mu <- 3; rho <- 0.4; k <- 2
x <- 0:200
p <- ddelaporte(x, mu = mu, rho = rho, k = k)
c(total = sum(p),
  mean = sum(x * p),
  cdf_matches_pmf = isTRUE(all.equal(cumsum(p)[1:20],
                                     pdelaporte(0:19, mu, rho, k))),
  log_tail_direct = pdelaporte(150, mu, rho, k, lower.tail = FALSE,
                               log.p = TRUE))
#>           total            mean cdf_matches_pmf log_tail_direct 
#>          1.0000          3.0000          1.0000       -140.5794
# a quantile is the smallest count whose CDF reaches the probability
q <- qdelaporte(0.9, mu, rho, k)
c(q = q, cdf_q = pdelaporte(q, mu, rho, k),
  cdf_q_minus_1 = pdelaporte(q - 1, mu, rho, k))
#>             q         cdf_q cdf_q_minus_1 
#>     6.0000000     0.9493058     0.8967050
# invalid parameters and non-integer counts
suppressWarnings(c(ddelaporte(2, mu = -1, rho = rho, k = k),
                   ddelaporte(2.5, mu = mu, rho = rho, k = k)))
#> [1] NaN   0
```

The log upper tail at 150 is finite, where $`1 - P(Y \le 150)`$ would be
0 in floating point.

## 5. Through a fitted model

Once registered, the family is used by every part of
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md):

``` r

set.seed(4)
d <- data.frame(x = rnorm(300))
d$y <- rnbinom(300, size = 2, mu = exp(0.5 + 0.4 * d$x))
fit <- countreg(y ~ x, data = d, family = "NB2", method = "BFGS")
mu_hat <- predict(fit)
alpha_hat <- exp(coef(fit)[["ln(alpha)"]])
# predict(type = "prob") uses the registry density
all.equal(unname(predict(fit, type = "prob", max_count = 3)[1, ]),
          nb2$dens(0:3, predicted = mu_hat[1], alpha = alpha_hat,
                   sigma = NULL))
#> [1] TRUE
# simulate() uses the registry's random draws
sims <- simulate(fit, nsim = 4000, seed = 1)
c(mean_ratio = mean(rowMeans(sims) / mu_hat),
  var_ratio = mean(apply(sims, 1, var) /
                     nb2$variance(mu_hat, alpha_hat, NULL)))
#> mean_ratio  var_ratio 
#>  0.9999219  0.9997107
# quantile residuals use the cdf hook
summary(residuals(fit, type = "quantile"))
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> -3.314260 -0.650915 -0.022552  0.006811  0.665310  2.549137
```

Offsets, weights, the rows used in the fit, and the alignment of
predictions with `newdata` are handled by the shared data preparation,
not by the family. A new family needs no code for them, but its tests
should include a model with an offset.

## 6. Capabilities are opt-in

A registry entry makes a family available in
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
(including zero models, which are built from `dens`) and in
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md).
Other features are enabled explicitly, because each depends on a
mathematical property that must be established for the family:

- **Underreporting** requires that binomial thinning of the family is
  the same family, with known parameters (`.fcr_thinning_families`,
  `.fcr_thin_params()`); the reported-count variance law is tested for
  each listed family.
- **Likelihood-ratio tests** against Poisson or another family require a
  nesting and its type: an interior restriction (ordinary chi-square) or
  a boundary one (a mixture of chi-squares). Families are listed in
  `.fcr_family_nesting` only when the restriction is regular.
- **Analytic gradients** exist only for Poisson, NB1 and NB2; other
  families use numerical derivatives.
- **Moment starting values** for the dispersion exist for NB1, NB2,
  PIG1, PIG2 and PIG (`.fcr_moment_disp_start()`); other families use
  `start`.
- **Exact random-parameter predictions** use the moment-generating
  function of the random parameter and the family’s mean adjustment, so
  `mean_adj` and `mean` must be right.

Diagnostics follow the same contract: Pearson residuals need `variance`,
quantile residuals and rootograms need `dens` and preferably `cdf`, and
marginal effects need the mean.

**Finite mixtures and latent classes**
([`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md),
[`countreg_lc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md))
use the same registry entry for each component: the log-scale `dens` for
the likelihood and posterior class probabilities, `mean` (or `mean_adj`)
for component and mixture expected counts, `variance` for the mixture
variance and Pearson residuals, `cdf` (or a summed PMF) for the mixture
distribution function and quantile residuals, and `rand` for
[`simulate()`](https://rdrr.io/r/stats/simulate.html). A new family
therefore works as a mixture component as soon as its entry is complete,
but it is reported as *experimental* until it is added to
`.fcr_fm_validated` together with tests: one-component agreement with
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
an independently coded two-component likelihood, and a recovery study
such as `inst/validation/finite-mixtures/mixture-recovery-study.R`.
Because a registry `dens` returns probability 0 (`-Inf` on the log
scale) for an invalid parameter value, as
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
needs, the mixture engine cannot tell such a value from a count with
probability zero by the density alone. It therefore checks the component
parameters itself and rejects a candidate whose parameters are not
finite and positive, or whose density raises an error, before it enters
the likelihood.

## 7. Tests

At least four kinds of tests are needed:

1.  **Special cases** against independent implementations, for example
    NB2 with $`\alpha \to 0`$ against
    [`dpois()`](https://rdrr.io/r/stats/Poisson.html), or a parameter
    value that gives a geometric distribution.
2.  **The distribution contract**: probabilities sum to 1, the CDF
    equals the summed PMF, quantiles are the smallest count reaching the
    probability, random draws match the probabilities, base-R argument
    handling, large means and extreme parameters stay finite.
3.  **Registry agreement**: the registry density and CDF equal the
    public functions at the regression parameterization, and the
    variance equals the summed second moment.
4.  **Estimation and prediction**:
    [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
    recovers the parameters of simulated data, and
    [`predict()`](https://rdrr.io/r/stats/predict.html),
    [`simulate()`](https://rdrr.io/r/stats/simulate.html) and the
    residuals agree with the registry.

A small recovery test, in the form used in `tests/testthat/`:

``` r

testthat::test_that("countreg() recovers NB2 parameters", {
  set.seed(11)
  n <- 4000
  d <- data.frame(x = rnorm(n), t = runif(n, 0.5, 2))
  d$y <- rnbinom(n, size = 1 / 0.5, mu = d$t * exp(0.5 + 0.4 * d$x))
  m <- countreg(y ~ x + offset(log(t)), data = d, family = "NB2",
                method = "BFGS")
  testthat::expect_true(m$model$converged)
  est <- coef(m)
  se <- summary(m)$coefficients[, 2]
  truth <- c(0.5, 0.4, log(0.5))
  # within four standard errors: fails rarely by chance, often if wrong
  testthat::expect_true(all(abs(est - truth) < 4 * se))
})
#> Test passed with 2 successes 🥇.
```

Tests involving random numbers set a seed and use tolerances justified
by the standard errors, not tolerances tuned until the test passes.

## 8. Contributor checklist

Support, parameters, links and the meaning of $`\mu`$ written down

Public `d`/`p`/`q`/`r` functions with roxygen documentation and
examples, exported (roxygen writes `NAMESPACE`)

Log-scale PMF and directly computed tails; vectorization, `NA`, invalid
parameters

Compiled code, if any, with
[`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html)
run

Registry entry: `params`, `dens`, `variance`, `rand`, `cdf`, `mean` or
`mean_adj`, `start`, `start_rp`

Opt-in features only where established: thinning, nesting

A row in each table of
[`?count_families`](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md),
and the help topic in the “Count distributions” section of
`_pkgdown.yml`

Tests of all four kinds above, and the family added to the test lists
that cover every family

A `NEWS.md` entry

`R CMD check` passes, including examples and vignettes
