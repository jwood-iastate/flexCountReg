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
  GP1 = fit("GP1"),
  PA = fit("PA"),
  DELAPORTE = fit("DELAPORTE")
)
sapply(models, nobs)

## ----glance-------------------------------------------------------------------
stats <- do.call(rbind, lapply(names(models), function(nm) {
  cbind(model = nm, glance(models[[nm]])[, c("logLik", "AIC", "BIC", "df")])
}))
stats$dAIC <- stats$AIC - min(stats$AIC)
stats[order(stats$AIC), ]

## ----best, include = FALSE----------------------------------------------------
ord <- stats[order(stats$AIC), ]
close <- ord$model[ord$dAIC < 2]
nb1_gap <- min(stats$dAIC[stats$model %in% c("NB1", "PIG1", "GP1")])

## ----delaporte----------------------------------------------------------------
summary(models$DELAPORTE)$coefficients[c("logit(rho)", "ln(k)"), 1:2]

## ----summary------------------------------------------------------------------
summary(models$NB2)

## ----lrtest-------------------------------------------------------------------
regCompTest(models$Poisson, model2 = models$NB2)$statistics

## ----not-nested---------------------------------------------------------------
regCompTest(models$NB1, model2 = models$NB2)$statistics

## ----table--------------------------------------------------------------------
regCompTable(models[c("Poisson", "NB1", "NB2")])

## ----shortcuts----------------------------------------------------------------
y <- washington_roads$Total_crashes
pearson <- function(m) {
  sum(residuals(m, type = "pearson")^2) / (nobs(m) - length(coef(m)))
}
expected_zeros <- sapply(models[c("Poisson", "NB2")], function(m)
  mean(predict(m, type = "prob", max_count = 0)[, 1]))
c(variance_to_mean = var(y) / mean(y),
  pearson_poisson = pearson(models$Poisson),
  observed_zeros = mean(y == 0),
  expected_zeros_poisson = expected_zeros[["Poisson"]],
  expected_zeros_nb2 = expected_zeros[["NB2"]])

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
  abline(0, 1, lty = 2)
}
par(oldpar)

## ----cure, message = FALSE, fig.height = 3------------------------------------
cureplot(models$NB2, indvar = "lnaadt")

## ----generalized--------------------------------------------------------------
nb2_gen <- countreg(f, data = washington_roads, family = "NB2",
                    offset = "lnlength",
                    dis_param_formula_1 = ~ speed50, method = "BFGS")
tidy(nb2_gen, conf.int = TRUE)

## ----zinb---------------------------------------------------------------------
zinb <- countreg(f, data = washington_roads, family = "NB2",
                 offset = "lnlength", zero_formula = ~ 1, method = "NR")
pi_zero <- unname(plogis(coef(zinb)["Zero:(Intercept)"]))
c(NB2 = as.numeric(logLik(models$NB2)), ZINB = as.numeric(logLik(zinb)),
  zero_probability = pi_zero)

## ----predict------------------------------------------------------------------
nd <- washington_roads[1:3, ]
predict(models$NB2, newdata = nd)
round(predict(models$NB2, newdata = nd, type = "prob", max_count = 4), 3)
predict(models$NB2, newdata = nd, interval = "prediction", level = 0.9)

## ----margeff------------------------------------------------------------------
margEffTable(models$NB2)

