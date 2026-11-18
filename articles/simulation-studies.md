# Simulation studies for count models

When the data-generating process is known, a simulation study shows how
an estimation procedure behaves: whether its estimates are biased,
whether its standard errors match the actual variability, and whether
its intervals cover the true value as often as they claim. This vignette
builds a small study step by step. Its replicate counts are chosen to
keep the vignette quick; they illustrate the workflow and do not
establish the properties of any method. A parameterized script for a
full-size study is included with the package (see the end).

``` r

library(flexCountReg)
```

## 1. Design the study before generating data

The estimand is the coefficient $`\beta_x`$ of a covariate $`x`$ in a
log-linear count model with an exposure offset:

``` math
E[Y_i] = t_i \exp(\beta_0 + \beta_x x_i), \quad
\beta_0 = 0.5, \; \beta_x = 0.4, \;
x_i \sim N(0, 1), \; t_i \sim U(0.5, 2).
```

The counts are Poisson or NB2 with $`\alpha = 0.5`$, with $`n = 100`$ or
$`400`$ independent rows. Two procedures are compared on every dataset:
a Poisson model and an NB2 model, both with Hessian-based standard
errors and 95% Wald intervals.

``` r

scenarios <- expand.grid(family = c("POISSON", "NB2"), n = c(100, 400),
                         stringsAsFactors = FALSE)
scenarios$beta0 <- 0.5
scenarios$beta_x <- 0.4
scenarios$alpha <- ifelse(scenarios$family == "NB2", 0.5, NA)
scenarios$id <- seq_len(nrow(scenarios))
scenarios
#>    family   n beta0 beta_x alpha id
#> 1 POISSON 100   0.5    0.4    NA  1
#> 2     NB2 100   0.5    0.4   0.5  2
#> 3 POISSON 400   0.5    0.4    NA  3
#> 4     NB2 400   0.5    0.4   0.5  4
```

## 2. A generator that can be checked

The data are generated with base R’s
[`rpois()`](https://rdrr.io/r/stats/Poisson.html) and
[`rnbinom()`](https://rdrr.io/r/stats/NegBinomial.html), so the
generator is independent of the code being evaluated:

``` r

generate <- function(sc) {
  x <- rnorm(sc$n)
  t <- runif(sc$n, 0.5, 2)
  mu <- t * exp(sc$beta0 + sc$beta_x * x)
  y <- if (sc$family == "POISSON") rpois(sc$n, mu) else
    rnbinom(sc$n, size = 1 / sc$alpha, mu = mu)
  data.frame(y, x, log_t = log(t))
}
```

Two checks against independent results. First, with $`\alpha = 1`$ the
NB2 distribution is geometric, so the generator’s frequencies must match
[`dgeom()`](https://rdrr.io/r/stats/Geometric.html):

``` r

set.seed(1)
draws <- rnbinom(1e5, size = 1, mu = 2)
round(rbind(simulated = tabulate(draws + 1, 6) / 1e5,
            geometric = dgeom(0:5, prob = 1 / 3)), 4)
#>             [,1]   [,2]   [,3]   [,4]   [,5]   [,6]
#> simulated 0.3318 0.2232 0.1480 0.0987 0.0656 0.0436
#> geometric 0.3333 0.2222 0.1481 0.0988 0.0658 0.0439
```

Second, the package’s own
[`simulate()`](https://rdrr.io/r/stats/simulate.html) method, used later
for fitted models, must reproduce the NB2 mean and variance at the
fitted parameters:

``` r

set.seed(2)
d <- generate(scenarios[scenarios$family == "NB2" & scenarios$n == 400, ])
fit <- countreg(y ~ x + offset(log_t), data = d, family = "NB2",
                method = "BFGS")
sims <- simulate(fit, nsim = 2000, seed = 3)
mu <- predict(fit)
alpha <- exp(coef(fit)[["ln(alpha)"]])
summary(data.frame(mean_ratio = rowMeans(sims) / mu,
                   var_ratio = apply(sims, 1, var) / (mu + alpha * mu^2)))
#>    mean_ratio       var_ratio     
#>  Min.   :0.9178   Min.   :0.8341  
#>  1st Qu.:0.9855   1st Qu.:0.9678  
#>  Median :1.0014   Median :1.0038  
#>  Mean   :1.0008   Mean   :1.0048  
#>  3rd Qu.:1.0173   3rd Qu.:1.0403  
#>  Max.   :1.0688   Max.   :1.1824
```

Both ratios are centred on 1, with the spread expected from 2,000 draws.

## 3. One replicate, with its status

A replicate generates one dataset, applies each procedure, and returns
one row per procedure. Failures and warnings are recorded, not dropped:

``` r

fit_one <- function(d, family) {
  warn <- character(0)
  m <- tryCatch(
    withCallingHandlers(
      countreg(y ~ x + offset(log_t), data = d, family = family,
               method = "BFGS"),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = function(e) e)
  if (inherits(m, "error")) {
    return(data.frame(estimate = NA, se = NA, status = "error",
                      note = conditionMessage(m)))
  }
  co <- summary(m)$coefficients
  status <- if (!isTRUE(m$model$converged)) "not converged" else
    if (length(warn)) "warning" else "ok"
  data.frame(estimate = co["x", 1], se = co["x", 2], status = status,
             note = paste(warn, collapse = "; "))
}

one_replicate <- function(sc, rep) {
  d <- generate(sc)
  out <- do.call(rbind, lapply(c("POISSON", "NB2"), function(proc)
    cbind(scenario = sc$id, rep = rep, procedure = proc, fit_one(d, proc))))
  out$truth <- sc$beta_x
  out
}
```

## 4. Reproducible random numbers

Each replicate gets its own stream of the L’Ecuyer-CMRG generator,
derived from one seed. The results do not depend on the order in which
replicates are run, so the same streams can be used by parallel workers
(for example with
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) or
`future.apply`), and any single replicate can be rerun on its own:

``` r

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
make_streams <- function(k) {
  s <- .Random.seed
  lapply(seq_len(k), function(i) s <<- parallel::nextRNGStream(s))
}
R <- 100
streams <- make_streams(nrow(scenarios) * R)
run_scenario <- function(sc) {
  do.call(rbind, lapply(seq_len(R), function(r) {
    assign(".Random.seed", streams[[(sc$id - 1) * R + r]],
           envir = globalenv())
    one_replicate(sc, r)
  }))
}
results <- do.call(rbind, lapply(split(scenarios, scenarios$id),
                                 run_scenario))
results <- merge(scenarios[c("id", "family", "n")], results,
                 by.x = "id", by.y = "scenario")
head(results)
#>   id  family   n rep procedure  estimate         se  status
#> 1  1 POISSON 100   1   POISSON 0.4112948 0.07481077      ok
#> 2  1 POISSON 100   1       NB2 0.4113006 0.07481521      ok
#> 3  1 POISSON 100   2   POISSON 0.5252295 0.07525128      ok
#> 4  1 POISSON 100   2       NB2 0.5252188 0.07525494 warning
#> 5  1 POISSON 100   3   POISSON 0.5619560 0.06897158      ok
#> 6  1 POISSON 100   3       NB2 0.5619552 0.06897677      ok
#>                                                                                                                                                                                 note
#> 1                                                                                                                                                                                   
#> 2                                                                                                                                                                                   
#> 3                                                                                                                                                                                   
#> 4 The Hessian is not negative definite at the estimates, so the standard errors of ln(alpha) are NA. The estimate may be on the boundary of the parameter space or not at a maximum.
#> 5                                                                                                                                                                                   
#> 6                                                                                                                                                                                   
#>   truth
#> 1   0.4
#> 2   0.4
#> 3   0.4
#> 4   0.4
#> 5   0.4
#> 6   0.4
table(results$procedure, results$status)
#>          
#>            ok warning
#>   NB2     388      12
#>   POISSON 400       0
```

The result is a tidy table with one row per scenario, replicate and
procedure. Notes on warnings are kept with the row they belong to.

## 5. Summaries with Monte Carlo standard errors

For each scenario and procedure: bias, empirical standard deviation of
the estimates, mean reported standard error, root mean squared error and
the coverage of the 95% Wald interval, each with its Monte Carlo
standard error (MCSE). Summaries use the successful fits; the attempted
and successful counts are shown so that the conditioning is visible.

``` r

summarize_sim <- function(res) {
  do.call(rbind, lapply(split(res, list(res$id, res$procedure)),
                        function(d) {
    ok <- d$status %in% c("ok", "warning") & is.finite(d$se)
    e <- d$estimate[ok]
    se <- d$se[ok]
    k <- sum(ok)
    cover <- abs(e - d$truth[ok]) <= qnorm(0.975) * se
    data.frame(data = d$family[1], n = d$n[1], procedure = d$procedure[1],
               attempted = nrow(d), used = k,
               bias = mean(e - d$truth[ok]), mcse_bias = sd(e) / sqrt(k),
               emp_sd = sd(e), mean_se = mean(se),
               rmse = sqrt(mean((e - d$truth[ok])^2)),
               coverage = mean(cover),
               mcse_coverage = sqrt(mean(cover) * (1 - mean(cover)) / k))
  }))
}
summary_table <- summarize_sim(results)
summary_table[order(summary_table$data, summary_table$n), -c(4, 5)] |>
  format(digits = 3)
#>              data   n procedure      bias mcse_bias emp_sd mean_se   rmse
#> 2.NB2         NB2 100       NB2  0.002615   0.01033 0.1033  0.1039 0.1028
#> 2.POISSON     NB2 100   POISSON  0.006447   0.01077 0.1077  0.0697 0.1074
#> 4.NB2         NB2 400       NB2 -0.003456   0.00541 0.0541  0.0516 0.0540
#> 4.POISSON     NB2 400   POISSON -0.001721   0.00586 0.0586  0.0338 0.0584
#> 1.NB2     POISSON 100       NB2  0.009175   0.00607 0.0607  0.0714 0.0611
#> 1.POISSON POISSON 100   POISSON  0.008750   0.00606 0.0606  0.0700 0.0609
#> 3.NB2     POISSON 400       NB2  0.000857   0.00321 0.0321  0.0341 0.0320
#> 3.POISSON POISSON 400   POISSON  0.000782   0.00320 0.0320  0.0335 0.0318
#>           coverage mcse_coverage
#> 2.NB2         0.94        0.0237
#> 2.POISSON     0.82        0.0384
#> 4.NB2         0.92        0.0271
#> 4.POISSON     0.76        0.0427
#> 1.NB2         0.97        0.0171
#> 1.POISSON     0.97        0.0171
#> 3.NB2         0.95        0.0218
#> 3.POISSON     0.95        0.0218
summary_table[, c("data", "n", "procedure", "attempted", "used")]
#>              data   n procedure attempted used
#> 1.NB2     POISSON 100       NB2       100  100
#> 2.NB2         NB2 100       NB2       100  100
#> 3.NB2     POISSON 400       NB2       100  100
#> 4.NB2         NB2 400       NB2       100  100
#> 1.POISSON POISSON 100   POISSON       100  100
#> 2.POISSON     NB2 100   POISSON       100  100
#> 3.POISSON POISSON 400   POISSON       100  100
#> 4.POISSON     NB2 400   POISSON       100  100
```

On NB2 data, the Poisson procedure’s reported standard errors are much
smaller than the actual spread of its estimates (for $`n = 400`$, a mean
standard error of 0.034 against an empirical standard deviation of
0.059), so its intervals cover the true value only 76% of the time (MCSE
4.3 points). The NB2 procedure’s coverage is 92%. Both estimators are
nearly unbiased: the Poisson estimator is consistent for $`\beta_x`$
under NB2 data, but its standard errors are not.

With 100 replicates, a coverage estimate near 95% has an MCSE of about
2.2 percentage points, too large to distinguish 93% from 95%.
Differences smaller than about two MCSEs are not evidence of anything.

12 of the 800 fits gave a warning. These are NB2 fits to Poisson data,
where the dispersion approaches its boundary (see
[`vignette("numerical-reliability", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/numerical-reliability.md)).
They are kept in the summaries, because the coefficient of $`x`$ is
still estimated well; a study whose target was the dispersion would have
to treat them separately.

## 6. Panel data: keep the independent unit

With repeated measurements, the independent unit is the group, and the
generator must give all rows of a group one shared effect. Here, 100
units are each observed 5 times, with a gamma-distributed unit effect
(mean 1, variance 0.5), a unit-level covariate $`z`$ and a row-level
covariate $`x`$:

``` r

generate_panel <- function(n_units = 100, T = 5, theta = 0.5) {
  id <- rep(seq_len(n_units), each = T)
  u <- rgamma(n_units, shape = 1 / theta, rate = 1 / theta)[id]
  z <- rnorm(n_units)[id]
  x <- rnorm(n_units * T)
  data.frame(id, z, x, y = rpois(n_units * T, u * exp(0.2 + 0.4 * z +
                                                        0.3 * x)))
}
```

The negative multinomial model of
[`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md)
(family `"NM2"`, the default) matches this process: Poisson counts given
a gamma unit effect with mean 1 and variance \$lpha\$ (here 0.5), shared
by the rows of a unit. The pooled Poisson model ignores the grouping:

``` r

panel_rep <- function(r) {
  assign(".Random.seed", panel_streams[[r]], envir = globalenv())
  d <- generate_panel()
  est <- function(m, nm) {
    co <- summary(m)$coefficients
    c(estimate = co[nm, 1], se = co[nm, 2])
  }
  pooled <- countreg(y ~ z + x, data = d, family = "POISSON",
                     method = "BFGS")
  nm <- negmultinom(y ~ z + x, group_var = "id", data = d)
  rbind(data.frame(rep = r, procedure = "pooled Poisson",
                   term = c("z", "x"),
                   rbind(est(pooled, "z"), est(pooled, "x")),
                   converged = pooled$model$converged),
        data.frame(rep = r, procedure = "negative multinomial",
                   term = c("z", "x"), rbind(est(nm, "z"), est(nm, "x")),
                   converged = nm$model$converged))
}
R_panel <- 100
panel_streams <- make_streams(R_panel)
panel <- do.call(rbind, lapply(seq_len(R_panel), panel_rep))
panel$truth <- ifelse(panel$term == "z", 0.4, 0.3)
panel$cover <- abs(panel$estimate - panel$truth) <=
  qnorm(0.975) * panel$se
aggregate(cbind(estimate, se, cover, converged) ~ procedure + term,
          data = panel, FUN = mean)
#>              procedure term  estimate         se cover converged
#> 1 negative multinomial    x 0.3026586 0.04191190  0.94         1
#> 2       pooled Poisson    x 0.3004392 0.03796208  0.88         1
#> 3 negative multinomial    z 0.4021180 0.08419030  0.98         1
#> 4       pooled Poisson    z 0.3933019 0.03882494  0.72         1
aggregate(estimate ~ procedure + term, data = panel, FUN = sd)
#>              procedure term   estimate
#> 1 negative multinomial    x 0.04688169
#> 2       pooled Poisson    x 0.05319242
#> 3 negative multinomial    z 0.07529461
#> 4       pooled Poisson    z 0.07446998
```

The pooled model’s standard error for the unit-level covariate $`z`$ is
far below the spread of its estimates, and its intervals undercover; the
negative multinomial model accounts for the shared effect. Bootstrap or
cross-validation in such a study must also resample whole units, never
rows (`stderr = "boot"` in the panel models does this).

## 7. A full-size study

The package includes a parameterized script for a larger version of this
study:

``` r

script <- system.file("validation", "simulation-studies",
                      "count-regression-study.R", package = "flexCountReg")
file.exists(script)
#> [1] TRUE
```

(`script` holds its location on your computer.) It runs the scenarios
above with a chosen number of replicates (the environment variable
`REPS`, default 1000) and seed (`SEED`), using one random-number stream
per replicate, and writes the replicate results, the summary table, and
a record of the settings, R and package versions and the random-number
generator. Run it from a terminal, for example:

    REPS=1000 SEED=2026 Rscript count-regression-study.R

## Checklist

- Write down the estimand, the scenarios and the procedures first.
- Check the generator against independent special cases.
- Return a status for every fit; never drop failures silently.
- Report MCSEs and the number of attempted and used replicates.
- Use independent, reproducible random-number streams per replicate.
- In panel designs, generate and resample whole units.
- Record the seed, settings and software versions with the results.
