## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 4)

## ----setup--------------------------------------------------------------------
library(flexCountReg)
data("washington_roads")

## ----fit----------------------------------------------------------------------
m <- countreg(Total_crashes ~ s(lnaadt) + speed50 + offset(lnlength),
              data = washington_roads, family = "NB2")
summary(m)

## ----plot---------------------------------------------------------------------
plot_smooths(m)

## ----estimates----------------------------------------------------------------
head(smooth_estimates(m, n = 50, seed = 1))

## ----sp-----------------------------------------------------------------------
m$model$smooth_lambda
stiff <- countreg(Total_crashes ~ s(lnaadt) + speed50 + offset(lnlength),
                  data = washington_roads, family = "NB2", sp = 100)
stiff$model$smooth_table[, c("term", "edf", "lambda")]

## ----aic----------------------------------------------------------------------
lin <- countreg(Total_crashes ~ lnaadt + speed50 + offset(lnlength),
                data = washington_roads, family = "NB2", method = "BFGS")
c(linear = AIC(lin), smooth = AIC(m))

## ----dispersion---------------------------------------------------------------
m2 <- countreg(Total_crashes ~ te(lnaadt, lnlength) + speed50,
               data = washington_roads, family = "NB2",
               dis_param_formula_1 = ~ s(lnlength, k = 5))
m2$model$smooth_table[, c("parameter", "term", "edf", "p.value")]

## ----plot2, fig.width = 7, fig.height = 3.5-----------------------------------
op <- par(mfrow = c(1, 2))
plot_smooths(m2)
par(op)

## ----rp-----------------------------------------------------------------------
rp <- countreg_rp(Total_crashes ~ s(lnaadt, k = 8) + offset(lnlength),
                  rpar_formula = ~ -1 + speed50, data = washington_roads,
                  family = "NB2", group_var = "ID", ndraws = 50)
summary(rp)

## ----panel--------------------------------------------------------------------
re <- repois(Total_crashes ~ s(lnaadt, k = 8) + speed50 + offset(lnlength),
             group_var = "ID", data = washington_roads)
re$model$smooth_table[, c("term", "edf", "lambda", "p.value")]

## ----mixture------------------------------------------------------------------
mx <- countreg_mix(Total_crashes ~ s(lnaadt, k = 8) + offset(lnlength),
                   data = washington_roads, family = "POISSON",
                   components = 2, n_starts = 3, seed = 1)
mx$model$smooth_table[, c("parameter", "term", "edf", "lambda")]

