## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(flexCountReg)
data("washington_roads")
w <- washington_roads
w$pdo <- w$Total_crashes - w$Injury_crashes - w$Fatal_crashes
w$injury <- w$Injury_crashes
w$fatal <- w$Fatal_crashes
colSums(w[, c("pdo", "injury", "fatal")])

## ----fit----------------------------------------------------------------------
m <- countreg_alloc(cbind(pdo, injury, fatal) ~ lnaadt + speed50 +
                      offset(lnlength),
                    alloc_formula = ~ lnaadt + speed50,
                    data = w, family = "NB2")
summary(m)

## ----truncated----------------------------------------------------------------
pos <- w[w$Total_crashes > 0, ]
f <- cbind(pdo, injury, fatal) ~ lnaadt + speed50 + offset(lnlength)
mt <- countreg_alloc(f, alloc_formula = ~ lnaadt + speed50, data = pos,
                     family = "NB2", truncated = TRUE)
mn <- countreg_alloc(f, alloc_formula = ~ lnaadt + speed50, data = pos,
                     family = "NB2")
se <- function(x) sqrt(diag(vcov(x)))[1:4]
round(rbind(all_rows = coef(m)[1:4], se = se(m),
            positive_truncated = coef(mt)[1:4], se = se(mt),
            positive_not_truncated = coef(mn)[1:4]), 3)

## ----predict------------------------------------------------------------------
nd <- w[c(1, 50, 200), ]
predict(m, newdata = nd)
predict(m, newdata = nd, type = "prob")

## ----thresholds---------------------------------------------------------------
mg <- countreg_alloc(cbind(pdo, injury, fatal) ~ lnaadt + speed50 +
                       offset(lnlength), alloc_formula = ~ lnaadt,
                     threshold_formula = ~ speed50, data = w,
                     family = "NB2")
AIC(m, mg)

## ----types--------------------------------------------------------------------
w$other <- w$Total_crashes - w$Animal - w$Rollover
colSums(w[, c("other", "Animal", "Rollover")])
mty <- countreg_alloc(cbind(other, Animal, Rollover) ~ lnaadt + speed50 +
                        offset(lnlength), alloc_formula = ~ lnaadt + speed50,
                      link = "mnl", data = w, family = "NB2")
summary(mty)

## ----types-predict------------------------------------------------------------
predict(mty, newdata = nd, type = "prob")
predict(mty, newdata = nd)

## ----fractional---------------------------------------------------------------
mf <- countreg_alloc(cbind(pdo, injury, fatal) ~ lnaadt + speed50 +
                       offset(lnlength), alloc_formula = ~ lnaadt + speed50,
                     data = w, family = "NB2", likelihood = "fractional",
                     group_var = "ID")
cbind(multinomial = coef(m), fractional = coef(mf))

## ----random-------------------------------------------------------------------
mr <- countreg_alloc(cbind(pdo, injury, fatal) ~ lnaadt + speed50 +
                       offset(lnlength), alloc_formula = ~ lnaadt + speed50,
                     data = w, family = "POISSON", group_var = "ID",
                     rpar_total = ~ 1, rpar_alloc = ~ 1,
                     correlation = "across", rpar_draws = 200, seed = 1)
summary(mr)
AIC(m, mr)

## ----posterior----------------------------------------------------------------
nd <- w[w$ID == w$ID[1], ]
round(cbind(predict(mr, newdata = nd),
            predict(mr, newdata = nd, random = "posterior")), 3)

## ----shared-------------------------------------------------------------------
ms <- countreg_alloc(cbind(pdo, injury, fatal) ~ lnaadt + speed50 +
                       offset(lnlength), alloc_formula = ~ lnaadt + speed50,
                     data = w, family = "POISSON", shared_formula = ~ 1)
summary(ms)
AIC(m, ms)

