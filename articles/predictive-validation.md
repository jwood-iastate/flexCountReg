# Evaluating predictive performance without information leakage

How well will a model predict counts it has not seen? A held-out
evaluation answers this, but only if the held-out data match the
prediction target and nothing about them leaks into the fitted model.
This vignette builds the evaluation by hand, so every step is visible.

``` r

library(flexCountReg)
```

## Define the target first

| Target | Example | Hold out |
|----|----|----|
| A new independent unit | a new child, a new road segment | rows (for independent rows) |
| A new unit from a population with repeated measurements | a new patient | whole units (all their rows) |
| A later observation of a known unit | next year’s crashes on a known segment | later periods, training on earlier ones only |

The folds, the model fitting and the scores all follow from the target.

## New independent units: school absences

[`MASS::quine`](https://rdrr.io/pkg/MASS/man/quine.html) has one row per
child (146 children), so rows are independent and a random split of the
rows matches the target “a new child from the same population”. Every
row gets a fold, recorded in a table with the row identifier:

``` r

q <- MASS::quine
q$row_id <- seq_len(nrow(q))
set.seed(2026)
K <- 5
q$fold <- sample(rep(seq_len(K), length.out = nrow(q)))
table(q$fold)
#> 
#>  1  2  3  4  5 
#> 30 29 29 29 29
```

Each model is refit on the training folds only, including everything
that is estimated from data (here the factor levels are fixed, and every
level of every factor occurs in every training set):

``` r

all(sapply(seq_len(K), function(k) {
  tr <- q[q$fold != k, ]
  all(sapply(c("Eth", "Sex", "Age", "Lrn"), function(v)
    all(table(tr[[v]]) > 0)))
}))
#> [1] TRUE
```

(If a level were missing from a training set,
[`predict()`](https://rdrr.io/r/stats/predict.html) would stop rather
than invent a coefficient for it; such a split has to be redone or the
level merged.)

The candidates use the same folds and the same held-out rows. For each
held-out child, the out-of-fold table records the prediction, the log
probability of the observed count and the predicted probability of zero.
Failed or non-converged fits are recorded, not dropped:

``` r

f <- Days ~ Eth + Sex + Age + Lrn
families <- c("POISSON", "NB2")
oof <- do.call(rbind, lapply(seq_len(K), function(k) {
  train <- q[q$fold != k, ]
  test <- q[q$fold == k, ]
  do.call(rbind, lapply(families, function(fam) {
    fit <- tryCatch(countreg(f, data = train, family = fam,
                             method = "BFGS"),
                    error = function(e) NULL)
    ok <- !is.null(fit) && isTRUE(fit$model$converged)
    if (!ok) {
      return(data.frame(row_id = test$row_id, fold = k, family = fam,
                        y = test$Days, mu = NA, log_p = NA, p_zero = NA,
                        status = "fit failed"))
    }
    P <- predict(fit, newdata = test, type = "prob",
                 max_count = max(test$Days))
    data.frame(row_id = test$row_id, fold = k, family = fam,
               y = test$Days, mu = predict(fit, newdata = test),
               log_p = log(P[cbind(seq_len(nrow(test)), test$Days + 1)]),
               p_zero = P[, 1], status = "ok")
  }))
}))
head(oof)
#>   row_id fold  family  y       mu      log_p       p_zero status
#> 1      1    1 POISSON  2 29.54024 -23.461884 1.481955e-13     ok
#> 2      4    1 POISSON  5 21.55623 -10.990394 4.347608e-10     ok
#> 3      5    1 POISSON  5 21.55623 -10.990394 4.347608e-10     ok
#> 4      8    1 POISSON 22 21.55623  -2.472784 4.347608e-10     ok
#> 5     17    1 POISSON 57 31.11360 -11.563683 3.072810e-14     ok
#> 6     20    1 POISSON 16 22.70434  -3.415303 1.379210e-10     ok
table(oof$family, oof$status)
#>          
#>            ok
#>   NB2     146
#>   POISSON 146
```

The log score uses the probability of the *actual* held-out count, read
from the column for that count (the probabilities run up to the largest
held-out count, so no count falls outside them).

### Scores

``` r

scores <- do.call(rbind, lapply(split(oof, oof$family), function(d) {
  data.frame(family = d$family[1], MAE = mean(abs(d$y - d$mu)),
             RMSE = sqrt(mean((d$y - d$mu)^2)),
             log_score = mean(d$log_p),
             expected_zeros = mean(d$p_zero),
             observed_zeros = mean(d$y == 0))
}))
scores
#>          family      MAE     RMSE log_score expected_zeros observed_zeros
#> NB2         NB2 11.57235 15.49119 -3.791640   0.0416381979     0.06164384
#> POISSON POISSON 11.47100 15.29440 -8.325457   0.0001116341     0.06164384
```

MAE and RMSE judge only the predicted *mean*, which is nearly the same
for both models (they share the mean model), so they cannot tell the
models apart. The log score judges the whole predicted distribution: the
NB2 model assigns far more probability to the absences actually observed
(a higher, i.e. less negative, log score). The Poisson model expects
almost no children with zero absences (0.01%), against 6.2% observed.

``` r

by_fold <- tapply(oof$log_p, list(oof$fold, oof$family), mean)
round(by_fold, 2)
#>     NB2 POISSON
#> 1 -3.50   -6.04
#> 2 -3.94   -8.70
#> 3 -3.77   -8.16
#> 4 -3.73   -9.60
#> 5 -4.02   -9.21
```

NB2 has the better log score in 5 of the 5 folds. The folds are not
independent replicates (their training sets overlap), so this is a
description of consistency, not a significance test.

## New patients: group the folds

[`MASS::epil`](https://rdrr.io/pkg/MASS/man/epil.html) has four rows per
patient. For the target “a new patient”, all rows of a patient must be
in the same fold: otherwise the training set contains the patient’s
other visits, and the evaluation is no longer about new patients.

``` r

e <- MASS::epil
ids <- unique(e$subject)
set.seed(7)
patient_fold <- setNames(sample(rep(1:5, length.out = length(ids))), ids)
e$fold <- patient_fold[as.character(e$subject)]
# no patient is split between folds
all(tapply(e$fold, e$subject, function(v) length(unique(v))) == 1)
#> [1] TRUE
```

Two models are compared: a pooled NB2 model and the random-effects model
[`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md).
For a patient not in the training data,
[`predict()`](https://rdrr.io/r/stats/predict.html) on the
random-effects model gives the population prediction (the patient’s own
effect is unknown):

``` r

fe <- y ~ trt + lbase + lage + V4
held_out <- do.call(rbind, lapply(1:5, function(k) {
  train <- e[e$fold != k, ]
  test <- e[e$fold == k, ]
  pooled <- countreg(fe, data = train, family = "NB2", method = "BFGS")
  panel <- renb(fe, group_var = "subject", data = train, method = "BFGS")
  data.frame(fold = k, subject = test$subject, y = test$y,
             pooled = predict(pooled, newdata = test),
             renb = predict(panel, newdata = test))
}))
#> Warning: The Hessian is singular; standard errors use a generalized inverse and
#> may be unreliable. The model may not be identified.
sapply(held_out[c("pooled", "renb")], function(p)
  c(MAE = mean(abs(held_out$y - p)), RMSE = sqrt(mean((held_out$y - p)^2))))
#>        pooled      renb
#> MAE  4.693061  309307.9
#> RMSE 9.193978 4735626.3
```

Only mean-based scores are shown for the random-effects model: a log
score would need the predictive distribution of a new patient’s count,
integrating over the random effect, and `predict(type = "prob")` is
available for
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
models only.

For comparison, the same evaluation with the rows split at random:

``` r

set.seed(7)
e$row_fold <- sample(rep(1:5, length.out = nrow(e)))
row_split <- do.call(rbind, lapply(1:5, function(k) {
  panel <- renb(fe, group_var = "subject", data = e[e$row_fold != k, ],
                method = "BFGS")
  test <- e[e$row_fold == k, ]
  data.frame(y = test$y, renb = predict(panel, newdata = test))
}))
c(MAE_row_split = mean(abs(row_split$y - row_split$renb)))
#> MAE_row_split 
#>      4.374853
```

The row split gives a different error (4.37 against 3.0930787^{5} with
patient-grouped folds), but it does not estimate the error for new
patients: each test row’s patient also contributed rows to the training
set. The patient-grouped split is the one that matches the target.

## Later visits of known patients

Forecasting a known patient’s next visit is a different target: train on
earlier periods only and condition on that patient’s history.
flexCountReg does not yet provide a forecast that conditions on a unit’s
past outcomes; its Empirical Bayes predictions (`method = "Individual"`
for random-parameter models) condition on the outcomes of the rows being
predicted, so using them to score those same rows would leak the answer.
This target is therefore not demonstrated here.

## Summary

- Choose the folds from the target: independent rows, whole units, or
  time order.
- Refit everything inside each training set and evaluate every candidate
  on the same held-out rows.
- Report failed fits instead of comparing models on whichever rows
  succeeded.
- Use the log score of the actual held-out count to compare predicted
  distributions; MAE and RMSE compare only means.

## Next steps

- [`vignette("diagnosing-count-models", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/diagnosing-count-models.md)
  covers in-sample diagnostics and a single held-out split.
- [`vignette("random-parameters-and-panels", package = "flexCountReg")`](https://jwood-iastate.github.io/flexCountReg/articles/random-parameters-and-panels.md)
  explains why the independent unit matters.
