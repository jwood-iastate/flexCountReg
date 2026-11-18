## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4)

## ----setup--------------------------------------------------------------------
library(flexCountReg)

## ----simulate-----------------------------------------------------------------
set.seed(101)
n <- 600
x <- rnorm(n)
z <- rbinom(n, 1, 0.5)
true_class <- ifelse(runif(n) < plogis(-0.8 + 1.6 * z), 1, 2)
mu <- ifelse(true_class == 1, exp(0.3 + 0.2 * x), exp(1.6 + 0.6 * x))
y <- rnbinom(n, size = 1 / 0.25, mu = mu)
d <- data.frame(y, x, z)

## ----interfaces---------------------------------------------------------------
lc <- countreg_lc(y ~ x, data = d, family = "NB2", classes = 2,
                  membership_formula = ~ z, n_starts = 4, seed = 1)
mx <- countreg_mix(y ~ x, data = d, family = "NB2", components = 2,
                   mixing_formula = ~ z, n_starts = 4, seed = 1)
c(latent_class = logLik(lc), mixture = logLik(mx))
all.equal(coef(lc), coef(mx))
class(lc)

## ----shared-------------------------------------------------------------------
lc_sd <- countreg_lc(y ~ x, data = d, family = "NB2", classes = 2,
                     membership_formula = ~ z,
                     shared = list(dis1 = TRUE), n_starts = 4, seed = 1)
lc_sx <- countreg_lc(y ~ x, data = d, family = "NB2", classes = 2,
                     membership_formula = ~ z,
                     shared = list(dis1 = TRUE, mean = "x"),
                     n_starts = 4, seed = 1)
data.frame(model = c("class-specific dispersion", "shared dispersion",
                     "shared dispersion and slope"),
           logLik = sapply(list(lc, lc_sd, lc_sx), logLik),
           df = sapply(list(lc, lc_sd, lc_sx), function(m) length(coef(m))),
           AIC = sapply(list(lc, lc_sd, lc_sx), AIC),
           BIC = sapply(list(lc, lc_sd, lc_sx), BIC))

## ----summary------------------------------------------------------------------
summary(lc_sd)

## ----choose-k-----------------------------------------------------------------
lc_1 <- countreg_lc(y ~ x, data = d, family = "NB2", classes = 1)
lc_3 <- countreg_lc(y ~ x, data = d, family = "NB2", classes = 3,
                    membership_formula = ~ z,
                    shared = list(dis1 = TRUE), n_starts = 4, seed = 1)
fits <- list(`K = 1` = lc_1, `K = 2` = lc_sd, `K = 3` = lc_3)
data.frame(K = 1:3,
           logLik = sapply(fits, logLik),
           df = sapply(fits, function(m) length(coef(m))),
           AIC = sapply(fits, AIC), BIC = sapply(fits, BIC),
           ICL = sapply(fits, function(m) m$model$ICL),
           starts_reaching_best = sapply(fits, function(m)
             sum(m$model$starts$reached_best)))

## ----probs--------------------------------------------------------------------
prior <- class_probs(lc_sd, conditioning = "prior")
post <- class_probs(lc_sd, conditioning = "posterior")
colnames(prior) <- paste0("prior_", colnames(prior))
colnames(post) <- paste0("posterior_", colnames(post))
head(cbind(d[, c("y", "x", "z")], round(prior, 3), round(post, 3)))
aggregate(cbind(prior_class1 = prior[, 1]) ~ z, data = d, FUN = mean)
table(estimated = max.col(post), generating = true_class)

## ----predict------------------------------------------------------------------
nd <- data.frame(x = c(-1, 0, 1, 0), z = c(0, 0, 0, 1), y = NA)
cbind(nd[, c("x", "z")],
      predict(lc_sd, newdata = nd, type = "component"),
      prior_class1 = class_probs(lc_sd, newdata = nd)[, 1],
      mixture = predict(lc_sd, newdata = nd, type = "response"))

## ----margeff------------------------------------------------------------------
margEffTable(lc_sd)

## ----mixed-families-----------------------------------------------------------
mf <- countreg_mix(y ~ x, data = d, family = c(low = "POISSON", high = "NB2"),
                   mixing_formula = ~ z, n_starts = 4, seed = 1)
coef(mf)
mf$model$classes[, c("label", "family", "prior_share", "mean_count")]
c(AIC_mixed = AIC(mf), AIC_NB2_NB2 = AIC(lc_sd))

## ----panel-data---------------------------------------------------------------
set.seed(202)
G <- 120; Tt <- 6
pd <- data.frame(id = rep(seq_len(G), each = Tt), t = rep(seq_len(Tt), G))
w_panel <- rbinom(G, 1, 0.5)
cls_panel <- ifelse(runif(G) < plogis(-0.5 + 1.5 * w_panel), 1, 2)
pd$w <- w_panel[pd$id]
pd$x <- rnorm(nrow(pd))
pd$y <- rpois(nrow(pd), exp(c(0, 1.2)[cls_panel[pd$id]] + 0.4 * pd$x))
pfit <- countreg_lc(y ~ x, data = pd, classes = 2, membership_formula = ~ w,
                    group_var = "id", shared = list(mean = "x"),
                    n_starts = 4, seed = 1)
coef(pfit)
c(rows = nobs(pfit), panels = attr(logLik(pfit), "nobs"))

## ----panel-history------------------------------------------------------------
history <- subset(pd, t <= 4)
target <- subset(pd, t >= 5)
head(history[, c("id", "t", "w", "x", "y")], 4)
head(target[, c("id", "t", "w", "x", "y")], 4)

pred_new <- predict(pfit, newdata = target)                      # prior
pred_hist <- predict(pfit, newdata = target, conditioning = "history",
                     history = history)
head(cbind(target[, c("id", "t", "y")], new_panel = pred_new,
           with_history = pred_hist))
c(MAE_new_panel = mean(abs(target$y - pred_new)),
  MAE_with_history = mean(abs(target$y - pred_hist)))

## ----panel-leakage------------------------------------------------------------
target_changed <- target
target_changed$y <- 999
all.equal(pred_hist,
          predict(pfit, newdata = target_changed, conditioning = "history",
                  history = history))

## ----diagnostics, fig.height = 3.5--------------------------------------------
oldpar <- par(mfrow = c(1, 2))
rootogram(lc_sd, main = "Two-class NB2")
qr <- residuals(lc_sd, type = "quantile")
qqnorm(qr, main = "Quantile residuals"); qqline(qr)
par(oldpar)

## ----cure, message = FALSE----------------------------------------------------
cureplot(lc_sd, indvar = "x")

## ----abalone------------------------------------------------------------------
data("abalone_rings")
set.seed(3)
ab <- abalone_rings[sample(nrow(abalone_rings), 1000), ]
ab$log_shell <- log(ab$shell_weight)
ab_1 <- countreg(rings ~ log_shell, data = ab, family = "POISSON",
                 method = "BFGS")
ab_2 <- countreg_lc(rings ~ log_shell, data = ab, family = "POISSON",
                    classes = 2, n_starts = 4, seed = 1)
coef(ab_2)
ab_2$model$classes[, c("label", "prior_share", "mean_count")]
ab_2$model$numerical$information
data.frame(model = c("Poisson", "2-class Poisson"),
           logLik = c(logLik(ab_1), logLik(ab_2)),
           AIC = c(AIC(ab_1), AIC(ab_2)), BIC = c(BIC(ab_1), BIC(ab_2)))

## ----inference, eval = FALSE--------------------------------------------------
# # Sandwich standard errors from unit-level (row or panel) scores
# fit_rob <- update(lc_sd, stderr = "robust")
# # Nonparametric bootstrap of units: each replicate refits the whole mixture
# # from the original estimates plus `boot_starts` new starts, and its classes
# # are aligned with the original ones before computing the covariance
# fit_boot <- update(lc_sd, stderr = "boot", bootstraps = 200,
#                    control = lc_control(boot_starts = 3))
# fit_boot$model$n_boot_ok        # replicates used
# fit_boot$model$n_boot_ambiguous # replicates with ambiguous class matching

