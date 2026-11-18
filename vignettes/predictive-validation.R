## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4)


## ----setup--------------------------------------------------------------------
library(flexCountReg)


## ----folds--------------------------------------------------------------------
q <- MASS::quine
q$row_id <- seq_len(nrow(q))
set.seed(2026)
K <- 5
q$fold <- sample(rep(seq_len(K), length.out = nrow(q)))
table(q$fold)


## ----levels-------------------------------------------------------------------
all(sapply(seq_len(K), function(k) {
  tr <- q[q$fold != k, ]
  all(sapply(c("Eth", "Sex", "Age", "Lrn"), function(v)
    all(table(tr[[v]]) > 0)))
}))


## ----cv-----------------------------------------------------------------------
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
table(oof$family, oof$status)


## ----scores-------------------------------------------------------------------
scores <- do.call(rbind, lapply(split(oof, oof$family), function(d) {
  data.frame(family = d$family[1], MAE = mean(abs(d$y - d$mu)),
             RMSE = sqrt(mean((d$y - d$mu)^2)),
             log_score = mean(d$log_p),
             expected_zeros = mean(d$p_zero),
             observed_zeros = mean(d$y == 0))
}))
scores


## ----score-folds--------------------------------------------------------------
by_fold <- tapply(oof$log_p, list(oof$fold, oof$family), mean)
round(by_fold, 2)


## ----epil-folds---------------------------------------------------------------
e <- MASS::epil
ids <- unique(e$subject)
set.seed(7)
patient_fold <- setNames(sample(rep(1:5, length.out = length(ids))), ids)
e$fold <- patient_fold[as.character(e$subject)]
# no patient is split between folds
all(tapply(e$fold, e$subject, function(v) length(unique(v))) == 1)


## ----epil-cv------------------------------------------------------------------
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
sapply(held_out[c("pooled", "renb")], function(p)
  c(MAE = mean(abs(held_out$y - p)), RMSE = sqrt(mean((held_out$y - p)^2))))


## ----epil-rows----------------------------------------------------------------
set.seed(7)
e$row_fold <- sample(rep(1:5, length.out = nrow(e)))
row_split <- do.call(rbind, lapply(1:5, function(k) {
  panel <- renb(fe, group_var = "subject", data = e[e$row_fold != k, ],
                method = "BFGS")
  test <- e[e$row_fold == k, ]
  data.frame(y = test$y, renb = predict(panel, newdata = test))
}))
c(MAE_row_split = mean(abs(row_split$y - row_split$renb)))


## ----row-vs-group, include = FALSE--------------------------------------------
grouped_mae <- mean(abs(held_out$y - held_out$renb))
row_mae <- mean(abs(row_split$y - row_split$renb))

