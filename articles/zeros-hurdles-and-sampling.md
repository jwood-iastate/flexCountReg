# Understanding zeros, hurdles and the sampling process

Many count datasets have a high density of zeros. That alone says
little: a count distribution with a small mean or a large variance
produces many zeros by itself. This vignette shows how to tell whether a
model needs a separate zero process, what zero-inflated and hurdle
models assume, how flexCountReg labels their components, and why a
zero-heavy sample is not the same thing as a zero-truncated sampling
design. For discussion on whether a zero-inflated or hurdle model is
appropriate, see Zeileis et al. (2008) and Lord et al.  (2007). Note
that there should be a substantive reason for believing that a separate
zero process exists, and that the zero-inflated and hurdle models are
not interchangeable.

``` r

library(flexCountReg)
```

## Three ways to produce zeros

| Process | Zeros come from | Model in flexCountReg |
|----|----|----|
| Ordinary count | the count distribution only | any family |
| Zero-inflated | a separate “always zero” group *and* the count distribution | `zero_formula`, `zero_type = "inflated"` |
| Hurdle | a binary process for zero vs positive; positive counts from a zero-truncated count distribution | `zero_formula`, `zero_type = "hurdle"` |

The component labels differ between the two zero models, and so do the
meanings of their coefficients:

|  | Zero-inflated | Hurdle |
|----|----|----|
| `Zero:` coefficients model (logit) | $`\pi`$ = probability of a structural zero | $`q`$ = probability of a *positive* count |
| A positive `Zero:` coefficient means | more structural zeros | fewer zeros |
| $`P(Y = 0)`$ | $`\pi + (1 - \pi) f(0)`$ | $`1 - q`$ |
| `predict(type = "zero")` returns | $`\pi`$ | $`P(Y = 0) = 1 - q`$ |

These are the conventions of the pscl package (Zeileis et al., 2008).

## A simulation: what each model finds

Three datasets of 1,000 counts with the same covariate `x` are
simulated: an ordinary NB2 count, a zero-inflated NB2 count (25%
structural zeros) and a hurdle count. All three have many zeros:

``` r

set.seed(12)
n <- 1000
x <- rnorm(n)
mu <- exp(0.3 + 0.5 * x)
y_nb <- rnbinom(n, size = 1, mu = mu)
y_zi <- ifelse(runif(n) < 0.25, 0, rnbinom(n, size = 1, mu = mu))
# hurdle: positive with probability 0.55, then a zero-truncated NB2 count
draw_pos <- function(m) { repeat { v <- rnbinom(1, size = 1, mu = m)
  if (v > 0) return(v) } }
y_h <- ifelse(runif(n) < 0.55, vapply(mu, draw_pos, 0), 0)
sims <- data.frame(x, y_nb, y_zi, y_h)
sapply(sims[-1], function(y) mean(y == 0))
#>  y_nb  y_zi   y_h 
#> 0.417 0.585 0.445
```

For each dataset, fit NB2, a zero-inflated NB2 and a hurdle NB2 with a
constant zero model, and compare AIC:

``` r

fit_three <- function(resp) {
  f <- reformulate("x", resp)
  list(NB2 = countreg(f, data = sims, family = "NB2", method = "BFGS"),
       ZINB = suppressWarnings(
         countreg(f, data = sims, family = "NB2", zero_formula = ~ 1,
                  method = "BFGS")),
       hurdle = countreg(f, data = sims, family = "NB2", zero_formula = ~ 1,
                         zero_type = "hurdle", method = "BFGS"))
}
sim_fits <- lapply(c(ordinary = "y_nb", zero_inflated = "y_zi",
                     hurdle = "y_h"), fit_three)
aic <- t(sapply(sim_fits, function(m) sapply(m, AIC)))
round(aic, 1)
#>                  NB2   ZINB hurdle
#> ordinary      3303.8 3305.7 3352.6
#> zero_inflated 2704.1 2700.5 2727.1
#> hurdle        3205.6 3205.0 3192.8
```

The lowest AIC in each row belongs to NB2, ZINB, hurdle (ordinary,
zero-inflated and hurdle data in turn). For the ordinary count, the
estimated structural-zero probability is small (0.032) and the extra
parameter only costs AIC. For the zero-inflated data, with 25%
structural zeros and 1,000 rows, the zero-inflated model improves the
AIC of NB2 by only 3.7: an NB2 distribution with a larger dispersion can
reproduce much of the same zero share. Real zero processes are rarely as
clean as these, so expect small differences.

## Real data: articles by PhD biochemists

[`pscl::bioChemists`](https://rdrr.io/pkg/pscl/man/bioChemists.html)
(Long, 1990) records the number of articles `art` published in the last
three years of the PhD by 915 biochemists, with gender `fem`, marital
status `mar`, the number of children under six `kid5`, the prestige of
the PhD department `phd` and the number of articles of the mentor `ment`
(not of the student).

``` r

e <- new.env()
utils::data("bioChemists", package = "pscl", envir = e)
b <- e$bioChemists
c(n = nrow(b), share_zero = mean(b$art == 0), mean = mean(b$art),
  variance = var(b$art))
#>           n  share_zero        mean    variance 
#> 915.0000000   0.3005464   1.6928962   3.7097416
```

### Baselines: how many zeros does an ordinary model expect?

``` r

f <- art ~ fem + mar + kid5 + phd + ment
pois <- countreg(f, data = b, family = "POISSON", method = "BFGS")
nb2 <- countreg(f, data = b, family = "NB2", method = "BFGS")
zeros <- c(observed = mean(b$art == 0),
           poisson = mean(predict(pois, type = "prob", max_count = 0)[, 1]),
           nb2 = mean(predict(nb2, type = "prob", max_count = 0)[, 1]))
round(zeros, 3)
#> observed  poisson      nb2 
#>    0.301    0.209    0.304
```

The Poisson model expects too few zeros, but the NB2 model, which allows
for overdispersion, already expects about as many zeros as observed. The
“excess” zeros relative to Poisson are explained by overdispersion; they
do not by themselves call for a separate zero process.

### Zero-inflated and hurdle models

``` r

zinb_const <- countreg(f, data = b, family = "NB2", zero_formula = ~ 1,
                       method = "BFGS")
zinb <- countreg(f, data = b, family = "NB2",
                 zero_formula = ~ ment + kid5, method = "BFGS")
hnb <- countreg(f, data = b, family = "NB2",
                zero_formula = ~ fem + mar + kid5 + phd + ment,
                zero_type = "hurdle", method = "BFGS")
```

With a constant zero model, the structural-zero probability is driven to
the boundary: the warning above says that its likelihood is flat, and
its standard error is `NA`:

``` r

summary(zinb_const)$coefficients["Zero:(Intercept)", 1:2]
#>   Estimate Std. Error 
#>  -12.32427   67.71443
plogis(coef(zinb_const)[["Zero:(Intercept)"]])
#> [1] 4.442605e-06
```

With covariates, the zero model can describe *which* students have more
zeros. In the zero-inflated model a negative `Zero:ment` coefficient
means fewer structural zeros among students of productive mentors; in
the hurdle model the `Zero:` coefficients describe the probability of
publishing at least one article:

``` r

zc <- function(m) {
  co <- summary(m)$coefficients
  co[grep("^Zero:", rownames(co)), 1:2]
}
zc(zinb)
#>                    Estimate Std. Error
#> Zero:(Intercept) -0.8786574  0.4036952
#> Zero:ment        -0.6113860  0.2458692
#> Zero:kid5         0.1564024  0.3617611
zc(hnb)
#>                     Estimate Std. Error
#> Zero:(Intercept)  0.23680758 0.29551909
#> Zero:femWomen    -0.25115052 0.15910530
#> Zero:marMarried   0.32623628 0.18081834
#> Zero:kid5        -0.28525098 0.11113046
#> Zero:phd          0.02221526 0.07955718
#> Zero:ment         0.08012173 0.01301809
```

flexCountReg’s likelihoods and zero-model coefficients agree with those
of pscl:

``` r

z_ref <- pscl::zeroinfl(art ~ fem + mar + kid5 + phd + ment | ment + kid5,
                        data = b, dist = "negbin")
h_ref <- pscl::hurdle(art ~ fem + mar + kid5 + phd + ment, data = b,
                      dist = "negbin")
c(zinb = as.numeric(logLik(zinb)), pscl_zinb = as.numeric(logLik(z_ref)),
  hurdle = as.numeric(logLik(hnb)), pscl_hurdle = as.numeric(logLik(h_ref)))
#>        zinb   pscl_zinb      hurdle pscl_hurdle 
#>   -1553.194   -1553.194   -1552.597   -1552.597
```

### Probability accounting

Each model’s fitted probabilities can be compared with the observed
frequencies, count by count, on the same rows. The overall expected
count (including the zero component) is the mean prediction:

``` r

models <- list(Poisson = pois, NB2 = nb2, ZINB = zinb, hurdle = hnb)
probs <- sapply(models, function(m)
  colMeans(predict(m, type = "prob", max_count = 5)))
acct <- rbind(round(cbind(observed = tabulate(b$art + 1, 6) / nrow(b),
                          probs), 3),
              mean_count = round(c(mean(b$art),
                                   sapply(models, function(m)
                                     mean(predict(m)))), 3))
acct
#>            observed Poisson   NB2  ZINB hurdle
#> 0             0.301   0.209 0.304 0.313  0.301
#> 1             0.269   0.310 0.272 0.255  0.277
#> 2             0.195   0.242 0.180 0.181  0.178
#> 3             0.092   0.135 0.106 0.111  0.105
#> 4             0.073   0.061 0.060 0.063  0.060
#> 5             0.030   0.025 0.033 0.035  0.034
#> mean_count    1.693   1.693 1.705 1.696  1.698
```

The mean predictions of all four models are close to the observed mean,
so the mean alone cannot separate them. The NB2, zero-inflated and
hurdle models reproduce the zeros and small counts about equally well;
the Poisson model does not.

``` r

data.frame(model = names(models),
           logLik = sapply(models, function(m) as.numeric(logLik(m))),
           AIC = sapply(models, AIC), row.names = NULL)
#>     model    logLik      AIC
#> 1 Poisson -1651.056 3314.113
#> 2     NB2 -1560.958 3135.917
#> 3    ZINB -1553.194 3126.388
#> 4  hurdle -1552.597 3131.193
```

The zero-inflated and hurdle models with covariates have lower AIC than
NB2. That shows that the probability of a zero varies with the
covariates in a way that the NB2 mean model does not capture; it does
not establish that some students are “structural non-publishers”. A
large share of zeros, or a better fit of a zero model, is evidence about
the fit, not proof of a particular mechanism.

### Which comparisons are valid

The NB2 model is the zero-inflated NB2 model with $`\pi = 0`$, which is
on the boundary of the parameter space (an intercept of $`-\infty`$ on
the logit scale), so the usual likelihood-ratio test does not apply and
flexCountReg does not report one. Compare such models with AIC or BIC on
the same rows, or on held-out data
([`vignette("predictive-validation", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/predictive-validation.md)).
When a zero-model probability is estimated near 0 or 1, its Wald
standard errors are unreliable, as the `NA` above shows.

## A zero-truncated sample is a different problem

A hurdle model says that zeros happen and come from a separate process.
A *zero-truncated sampling design* says that zeros cannot be observed at
all: for example, counts of visits recorded only for people who visited
at least once. Having no zeros in the data, or few zeros by chance, is
not the same as truncation. Under truncation the likelihood of each
observation is $`f(y) / (1 - f(0))`$ for $`y \ge 1`$; this is the
positive part of a hurdle model, but fitting a hurdle model to truncated
data (where no zeros exist to estimate the binary part) is not a
substitute. flexCountReg does not yet provide a zero-truncated
estimator.

## Reproducibility and next steps

The simulation uses `set.seed(12)`; the bioChemists data come from the
pscl package.

- [`vignette("diagnosing-count-models", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/diagnosing-count-models.md)
  shows calibration checks for zeros and large counts.
- [`vignette("choosing-a-count-distribution", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/choosing-a-count-distribution.md)
  compares observed and expected zeros on road-crash data.

## References

Long, J. S. (1990). The origins of sex differences in science. *Social
Forces*, 68(4), 1297-1316.

Lord, Dominique, Simon Washington, and John N. Ivan. (2007). Further
notes on  
the application of zero-inflated models in highway safety. *Accident
Analysis & Prevention* 39(1), 53-57.

Zeileis, A., Kleiber, C. and Jackman, S. (2008). Regression models for
count data in R. *Journal of Statistical Software*, 27(8), 1-25.
