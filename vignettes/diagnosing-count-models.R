## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 3.6)

## ----setup--------------------------------------------------------------------
library(flexCountReg)

## ----simulate-----------------------------------------------------------------
set.seed(101)
n <- 600
x <- runif(n, -1.5, 1.5)
sim <- data.frame(
  x = x,
  y_ok = rpois(n, exp(0.5 + 0.8 * x)),
  y_curve = rpois(n, exp(0.3 + 0.8 * x + 0.6 * x^2)),
  y_disp = rnbinom(n, size = 1 / 0.8, mu = exp(0.5 + 0.8 * x)))

ok <- countreg(y_ok ~ x, data = sim, family = "POISSON", method = "BFGS")
curve_lin <- countreg(y_curve ~ x, data = sim, family = "POISSON",
                      method = "BFGS")
disp_pois <- countreg(y_disp ~ x, data = sim, family = "POISSON",
                      method = "BFGS")

## ----pearson------------------------------------------------------------------
pearson <- function(fit) {
  sum(residuals(fit, type = "pearson")^2) / (nobs(fit) - length(coef(fit)))
}
round(c(ok = pearson(ok), curve = pearson(curve_lin),
        dispersion = pearson(disp_pois)), 2)

## ----resid-x, fig.cap = "Response residuals against x for the curved mean, before (left) and after (right) adding x^2."----
curve_quad <- countreg(y_curve ~ x + I(x^2), data = sim, family = "POISSON",
                       method = "BFGS")
op <- par(mfrow = c(1, 2))
for (f in list(curve_lin, curve_quad)) {
  r <- residuals(f, type = "response")
  plot(sim$x, r, col = "grey50", xlab = "x", ylab = "Response residual",
       main = deparse(formula(f)[[3]]))
  lines(lowess(sim$x, r), lwd = 3)
  abline(h = 0, lty = 2)
}
par(op)
round(c(before = pearson(curve_lin), after = pearson(curve_quad)), 2)

## ----rootogram-disp, fig.cap = "Hanging rootograms for the overdispersed count: Poisson (left) and NB2 (right)."----
disp_nb2 <- countreg(y_disp ~ x, data = sim, family = "NB2", method = "BFGS")
op <- par(mfrow = c(1, 2))
rootogram(disp_pois, max_count = 12, main = "Poisson")
rootogram(disp_nb2, max_count = 12, main = "NB2")
par(op)

## ----rootogram-cutoff---------------------------------------------------------
full <- rootogram(disp_pois, plot = FALSE)
cut8 <- rootogram(disp_pois, max_count = 8, plot = FALSE)
all.equal(full$expected[1:9], cut8$expected)
attr(cut8, "tail")    # observed and expected counts above 8

## ----qres-seeds, fig.cap = "Normal Q-Q plots of quantile residuals for the correctly specified model, with two different seeds."----
set.seed(1)
q1 <- residuals(ok, type = "quantile")
set.seed(2)
q2 <- residuals(ok, type = "quantile")
op <- par(mfrow = c(1, 2))
qqnorm(q1, main = "seed 1"); qqline(q1)
qqnorm(q2, main = "seed 2"); qqline(q2)
par(op)
cor(q1, q2)

## ----qres-disp, fig.cap = "Quantile residuals against fitted means for the overdispersed count: Poisson (left) and NB2 (right)."----
set.seed(3)
op <- par(mfrow = c(1, 2))
for (f in list(disp_pois, disp_nb2)) {
  plot(predict(f), residuals(f, type = "quantile"), xlab = "Fitted mean",
       ylab = "Quantile residual", main = f$model$family)
  abline(h = c(-2, 0, 2), lty = c(2, 1, 2))
}
par(op)

## ----calibration--------------------------------------------------------------
calib <- function(fit, y, k = 8) {
  P <- predict(fit, type = "prob", max_count = 30, overflow = TRUE)
  c(zeros_expected = mean(P[, 1]),
    above_k_expected = mean(rowSums(P[, -(1:(k + 1))])))
}
rbind(observed = c(mean(sim$y_disp == 0), mean(sim$y_disp > 8)),
      poisson = calib(disp_pois, sim$y_disp),
      nb2 = calib(disp_nb2, sim$y_disp))

## ----cure, fig.cap = "CURE plots for the curved mean: linear (top) and quadratic (bottom) models.", fig.height = 3.2, message = FALSE----
cureplot(curve_lin, indvar = "x")
cureplot(curve_quad, indvar = "x")

## ----step-sim-----------------------------------------------------------------
set.seed(102)
step_dat <- data.frame(x = runif(n, -1.5, 1.5))
step_dat$above <- as.numeric(step_dat$x > 0.5)
step_dat$y <- rpois(n, exp(0.4 + 0.5 * step_dat$x + 0.7 * step_dat$above))
no_ind <- countreg(y ~ x, data = step_dat, family = "POISSON",
                   method = "BFGS")
with_ind <- countreg(y ~ x + above, data = step_dat, family = "POISSON",
                     method = "BFGS")

## ----cure-step, fig.cap = "CURE plots along x without (top) and with (bottom) the indicator for x > 0.5.", fig.height = 3.2, message = FALSE----
cureplot(no_ind, indvar = "x")
cureplot(with_ind, indvar = "x")

## ----step-slopes--------------------------------------------------------------
seg_mean <- function(fit, data, v, breaks) {
  r <- data$y - predict(fit, newdata = data)
  round(tapply(r, cut(data[[v]], breaks, include.lowest = TRUE), mean), 3)
}
rbind(without_indicator = seg_mean(no_ind, step_dat, "x", c(-1.5, 0.5, 1.5)),
      with_indicator = seg_mean(with_ind, step_dat, "x", c(-1.5, 0.5, 1.5)))

## ----wa-base------------------------------------------------------------------
data("washington_roads")
wa <- washington_roads
wa$AADT10kplus <- as.numeric(wa$AADT > 10000)
nb_base <- countreg(Total_crashes ~ lnaadt + lnlength + speed50, data = wa,
                    family = "NB2", method = "BFGS")

## ----cure-wa-base, fig.cap = "CURE plot along lnaadt for the model without the high-volume indicator.", fig.height = 3.2, message = FALSE----
cureplot(nb_base, indvar = "lnaadt")

## ----wa-updated---------------------------------------------------------------
nb_ind <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
                     AADT10kplus, data = wa, family = "NB2",
                   method = "BFGS")
nb_int <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
                     AADT10kplus + AADT10kplus:lnaadt, data = wa,
                   family = "NB2", method = "BFGS")

## ----cure-wa-ind, fig.cap = "CURE plot along lnaadt with the high-volume indicator.", fig.height = 3.2, message = FALSE----
cureplot(nb_ind, indvar = "lnaadt")

## ----wa-compare---------------------------------------------------------------
wa_seg <- function(fit) {
  r <- wa$Total_crashes - predict(fit, newdata = wa)
  round(tapply(r, cut(wa$lnaadt, c(-Inf, 8.5, log(10000), Inf)), mean), 3)
}
cbind(rbind(base = wa_seg(nb_base), indicator = wa_seg(nb_ind),
            interaction = wa_seg(nb_int)),
      AIC = round(c(AIC(nb_base), AIC(nb_ind), AIC(nb_int)), 1))

## ----holdout------------------------------------------------------------------
set.seed(7)
train_rows <- sample(n, 400)
train <- sim[train_rows, ]
test <- sim[-train_rows, ]
log_score <- function(fit, test) {
  P <- predict(fit, newdata = test, type = "prob",
               max_count = max(test$y_disp))
  mean(log(P[cbind(seq_len(nrow(test)), test$y_disp + 1)]))
}
fit_p <- countreg(y_disp ~ x, data = train, family = "POISSON",
                  method = "BFGS")
fit_nb <- countreg(y_disp ~ x, data = train, family = "NB2",
                   method = "BFGS")
round(c(poisson = log_score(fit_p, test), nb2 = log_score(fit_nb, test)), 3)

## ----insects, fig.cap = "Quantile residuals of the NB1 model for InsectSprays."----
insects <- countreg(count ~ spray, data = datasets::InsectSprays,
                    family = "NB1", method = "BFGS")
set.seed(4)
plot(predict(insects), residuals(insects, type = "quantile"),
     xlab = "Fitted mean", ylab = "Quantile residual")
abline(h = c(-2, 0, 2), lty = c(2, 1, 2))

