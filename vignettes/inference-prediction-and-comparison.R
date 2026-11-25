## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(flexCountReg)
data("washington_roads")
d <- washington_roads

## ----stderr-------------------------------------------------------------------
f <- Total_crashes ~ lnaadt + lnlength + speed50
nb_h <- countreg(f, d, family = "NB2", method = "BFGS")
nb_r <- countreg(f, d, family = "NB2", method = "BFGS", stderr = "robust")
set.seed(1)
nb_b <- countreg(f, d, family = "NB2", method = "BFGS", stderr = "boot",
                 bootstraps = 50)
round(cbind(Hessian = sqrt(diag(vcov(nb_h))),
            robust = sqrt(diag(vcov(nb_r))),
            bootstrap = sqrt(diag(vcov(nb_b)))), 4)

## ----boot-info----------------------------------------------------------------
nb_b$model$n_boot_ok
nb_b$model$boot_failures

## ----margeff------------------------------------------------------------------
margEffTable(nb_h, digits = 4)

## ----predict------------------------------------------------------------------
nd <- d[1:3, ]
predict(nb_h, newdata = nd)
round(predict(nb_h, newdata = nd, type = "prob", max_count = 4,
              overflow = TRUE), 3)

## ----predict-counts-----------------------------------------------------------
round(predict(nb_h, newdata = nd, type = "prob", counts = 3:5,
              overflow = TRUE, underflow = TRUE), 4)

## ----intervals----------------------------------------------------------------
confint(nb_h)["speed50", ]
nd2 <- data.frame(lnaadt = log(c(3000, 15000)), lnlength = 0, speed50 = 1)
predict(nb_h, newdata = nd2, interval = "confidence")
predict(nb_h, newdata = nd2, interval = "prediction", level = 0.9)

## ----residuals----------------------------------------------------------------
r <- residuals(nb_h, type = "quantile")
c(mean = mean(r), sd = sd(r))

## ----simulate-----------------------------------------------------------------
sims <- simulate(nb_h, nsim = 200, seed = 1)
# Observed share of zeros against its simulated distribution
c(observed = mean(d$Total_crashes == 0),
  simulated = mean(colMeans(sims == 0)))

## ----glance-------------------------------------------------------------------
po <- countreg(f, d, family = "POISSON", method = "BFGS")
nb1 <- countreg(f, d, family = "NB1", method = "BFGS")
pl <- countreg(f, d, family = "PL", method = "BFGS")
do.call(rbind, lapply(list(Poisson = po, NB1 = nb1, NB2 = nb_h, PL = pl),
                      function(m) glance(m)[, c("logLik", "AIC", "BIC")]))

## ----summary-lr---------------------------------------------------------------
summary(nb_h)$lr_test
is.null(summary(pl)$lr_test)   # no test: PL does not nest Poisson

## ----regcomptest--------------------------------------------------------------
r <- regCompTest(po, model2 = nb_h)
r[c("LR", "LRdof", "LR_type", "LR_pvalue")]
r <- regCompTest(nb1, model2 = countreg(f, d, family = "NBP",
                                        method = "BFGS"))
r[c("LR", "LRdof", "LR_type", "LR_pvalue")]

## ----table--------------------------------------------------------------------
regCompTable(list(Poisson = po, NB2 = nb_h))

