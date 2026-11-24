## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(flexCountReg)
data("washington_roads")
w <- washington_roads
table(w$Year)

## ----step---------------------------------------------------------------------
m <- countreg(Total_crashes ~ tv(lnaadt, time = Year) + speed50 +
                offset(lnlength), data = w, family = "NB2",
              method = "BFGS")
coef(m)

## ----effects------------------------------------------------------------------
temporal_effects(m)

## ----tests--------------------------------------------------------------------
temporal_tests(m)

## ----bases--------------------------------------------------------------------
fit_basis <- function(b) {
  countreg(stats::as.formula(paste0(
    "Total_crashes ~ tv(lnaadt, time = Year, basis = '", b,
    "') + speed50 + offset(lnlength)")),
    data = w, family = "NB2", method = "BFGS")
}
mh <- fit_basis("hinge")
c(step = logLik(m), hinge = logLik(mh))
cbind(step = temporal_effects(m)$estimate,
      hinge = temporal_effects(mh)$estimate)

## ----smooth-------------------------------------------------------------------
msm <- fit_basis("smooth")
temporal_effects(msm)

## ----intercept----------------------------------------------------------------
mi <- countreg(Total_crashes ~ tv(1, time = Year) + lnaadt +
                 offset(lnlength), data = w, family = "NB2",
               dis_param_formula_1 = ~ tv(1, time = Year, basis = "level"),
               method = "BFGS")
temporal_effects(mi)

## ----predict------------------------------------------------------------------
mn <- countreg(Total_crashes ~ tv(lnaadt, time = Year, new_time = "nearest") +
                 speed50 + offset(lnlength), data = w, family = "NB2",
               method = "BFGS")
nd <- w[1:3, ]
nd$Year <- 2019
predict(mn, newdata = nd)

## ----rp-----------------------------------------------------------------------
r <- countreg_rp(Total_crashes ~ tv(lnaadt, time = Year) + offset(lnlength),
                 rpar_formula = ~ -1 + speed50, data = w, family = "NB2",
                 ndraws = 100, method = "BFGS", group_var = "ID")
temporal_effects(r)

