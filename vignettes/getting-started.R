## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4)


## ----setup--------------------------------------------------------------------
library(flexCountReg)
d <- datasets::InsectSprays
str(d)


## ----inspect------------------------------------------------------------------
anyNA(d)
table(d$spray)
range(d$count)


## ----describe-plot, fig.cap = "Insect counts by spray (each point is one unit; the bar marks the mean)."----
set.seed(1)   # the horizontal jitter is random
stripchart(count ~ spray, data = d, vertical = TRUE, method = "jitter",
           pch = 1, xlab = "Spray", ylab = "Insect count")
means <- tapply(d$count, d$spray, mean)
segments(seq_along(means) - 0.25, means, seq_along(means) + 0.25, means,
         lwd = 3)


## ----describe-table-----------------------------------------------------------
grp <- data.frame(mean = tapply(d$count, d$spray, mean),
                  variance = tapply(d$count, d$spray, var))
grp$ratio <- grp$variance / grp$mean
round(grp, 2)


## ----poisson------------------------------------------------------------------
pois <- countreg(count ~ spray, data = d, family = "POISSON",
                 method = "BFGS")
summary(pois)


## ----poisson-tidy-------------------------------------------------------------
tidy(pois, conf.int = TRUE, exponentiate = TRUE)


## ----glm-check----------------------------------------------------------------
ref <- glm(count ~ spray, data = d, family = poisson())
all.equal(unname(coef(pois)), unname(coef(ref)), tolerance = 1e-6)
all.equal(as.numeric(predict(pois)), as.numeric(fitted(ref)),
          tolerance = 1e-6)


## ----pearson------------------------------------------------------------------
pearson <- function(fit, k) sum(residuals(fit, type = "pearson")^2) /
  (nobs(fit) - k)
pearson(pois, length(coef(pois)))


## ----rootogram-pois, fig.cap = "Hanging rootogram of the Poisson model."------
rootogram(pois, main = "Poisson")


## ----qres-pois, fig.cap = "Randomized quantile residuals of the Poisson model against the fitted means."----
set.seed(2026)
qr_pois <- residuals(pois, type = "quantile")
plot(predict(pois), qr_pois, xlab = "Fitted mean",
     ylab = "Quantile residual", main = "Poisson")
abline(h = c(-2, 0, 2), lty = c(2, 1, 2))


## ----alternatives-------------------------------------------------------------
nb2 <- countreg(count ~ spray, data = d, family = "NB2", method = "BFGS")
nb1 <- countreg(count ~ spray, data = d, family = "NB1", method = "BFGS")
c(poisson = nobs(pois), nb2 = nobs(nb2), nb1 = nobs(nb1))


## ----nb2-summary--------------------------------------------------------------
summary(nb2)


## ----compare------------------------------------------------------------------
cmp <- data.frame(
  model = c("Poisson", "NB2", "NB1"),
  logLik = c(logLik(pois), logLik(nb2), logLik(nb1)),
  parameters = c(length(coef(pois)), length(coef(nb2)), length(coef(nb1))),
  AIC = c(AIC(pois), AIC(nb2), AIC(nb1)),
  BIC = c(BIC(pois), BIC(nb2), BIC(nb1)),
  LR_p_vs_Poisson = c(NA, summary(nb2)$lr_test$p.value,
                      summary(nb1)$lr_test$p.value))
cmp[, -1] <- round(cmp[, -1], 3)
cmp


## ----not-nested---------------------------------------------------------------
regCompTest(nb1, model2 = nb2)$statistics


## ----effects------------------------------------------------------------------
round(cbind(Poisson = coef(pois)[1:6], NB2 = coef(nb2)[1:6],
            NB1 = coef(nb1)[1:6]), 3)


## ----rootogram-nb2, fig.cap = "Hanging rootogram of the NB2 model."-----------
rootogram(nb2, main = "NB2")


## ----predict------------------------------------------------------------------
new_spray <- data.frame(spray = factor(levels(d$spray),
                                       levels = levels(d$spray)))
ci <- predict(nb2, newdata = new_spray, interval = "confidence")
pr <- predict(nb2, newdata = new_spray, interval = "prediction",
              level = 0.90)
pred <- data.frame(spray = new_spray$spray, expected = ci$fit,
                   ci_lower = ci$lwr, ci_upper = ci$upr,
                   pi90_lower = pr$lwr, pi90_upper = pr$upr)
pred[, 2:4] <- round(pred[, 2:4], 2)
pred
# Sprays below the overall mean, and whether their confidence intervals
# are separated from those of the other sprays
low <- ci$fit < mean(d$count)
separated <- max(ci$upr[low]) < min(ci$lwr[!low])


## ----predict-plot, fig.cap = "Expected counts with 95% confidence intervals (thick) and 90% prediction intervals for a new unit (thin)."----
x <- seq_len(nrow(pred))
plot(x, pred$expected, ylim = range(0, pr$upr), xaxt = "n", pch = 19,
     xlab = "Spray", ylab = "Insect count")
axis(1, at = x, labels = pred$spray)
segments(x, pr$lwr, x, pr$upr, lwd = 1)
segments(x, ci$lwr, x, ci$upr, lwd = 4)


## ----new-level----------------------------------------------------------------
tryCatch(predict(nb2, newdata = data.frame(spray = "G")),
         error = function(e) conditionMessage(e))


## ----script, eval = FALSE-----------------------------------------------------
# library(flexCountReg)
# d <- datasets::InsectSprays
# pois <- countreg(count ~ spray, data = d, family = "POISSON",
#                  method = "BFGS")
# nb2 <- countreg(count ~ spray, data = d, family = "NB2", method = "BFGS")
# summary(nb2)
# rootogram(nb2)
# new_spray <- data.frame(spray = factor(levels(d$spray),
#                                        levels = levels(d$spray)))
# predict(nb2, newdata = new_spray, interval = "confidence")
# predict(nb2, newdata = new_spray, interval = "prediction", level = 0.90)

