## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 3.5)

## ----setup--------------------------------------------------------------------
library(flexCountReg)
data("washington_roads")

## ----fits---------------------------------------------------------------------
f <- Total_crashes ~ lnaadt + speed50
fit <- function(family) {
  countreg(f, data = washington_roads, family = family, offset = "lnlength",
           method = "BFGS")
}
models <- list(
  Poisson = fit("Poisson"),
  NB1 = fit("NB1"),
  NB2 = fit("NB2"),
  PIG1 = fit("PIG1"),
  GP1 = fit("GP1")
)

## ----glance-------------------------------------------------------------------
do.call(rbind, lapply(names(models), function(nm) {
  cbind(model = nm, glance(models[[nm]])[, c("logLik", "AIC", "BIC", "df")])
}))

## ----summary------------------------------------------------------------------
summary(models$NB2)

## ----lrtest-------------------------------------------------------------------
regCompTest(models$Poisson, model2 = models$NB2)$statistics

## ----table--------------------------------------------------------------------
regCompTable(models[c("Poisson", "NB1", "NB2")])

## ----rootograms---------------------------------------------------------------
oldpar <- par(mfrow = c(1, 2))
rootogram(models$Poisson, main = "Poisson")
rootogram(models$NB2, main = "NB2")
par(oldpar)

## ----qq-----------------------------------------------------------------------
set.seed(1)
oldpar <- par(mfrow = c(1, 2))
for (nm in c("Poisson", "NB2")) {
  r <- residuals(models[[nm]], type = "quantile")
  qqnorm(r, main = nm, pch = 20, col = "grey40")
  abline(0, 1, col = "firebrick")
}
par(oldpar)

## ----cure, eval = FALSE-------------------------------------------------------
# cureplot(models$NB2, indvar = "lnaadt", n_resamples = 20)

## ----generalized--------------------------------------------------------------
nb2_gen <- countreg(f, data = washington_roads, family = "NB2",
                    offset = "lnlength",
                    dis_param_formula_1 = ~ speed50, method = "BFGS")
tidy(nb2_gen, conf.int = TRUE)

## ----zinb---------------------------------------------------------------------
zinb <- countreg(f, data = washington_roads, family = "NB2",
                 offset = "lnlength", zero_formula = ~ 1, method = "NR")
c(NB2 = as.numeric(logLik(models$NB2)), ZINB = as.numeric(logLik(zinb)))
# Estimated probability of a structural zero
unname(plogis(coef(zinb)["Zero:(Intercept)"]))

## ----predict------------------------------------------------------------------
nd <- washington_roads[1:3, ]
predict(models$NB2, newdata = nd)
round(predict(models$NB2, newdata = nd, type = "prob", max_count = 4), 3)

## ----margeff------------------------------------------------------------------
margEffTable(models$NB2)

