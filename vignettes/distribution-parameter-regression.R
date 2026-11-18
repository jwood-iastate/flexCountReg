## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4)


## ----setup--------------------------------------------------------------------
library(flexCountReg)


## ----simulate-----------------------------------------------------------------
set.seed(5)
n <- 800
sim <- data.frame(x = runif(n, -1, 1), z = rbinom(n, 1, 0.5))
sim$y <- rnbinom(n, size = 1 / exp(-1 + 1.2 * sim$z),
                 mu = exp(1 + 0.5 * sim$x))


## ----fits---------------------------------------------------------------------
nb_const <- countreg(y ~ x, data = sim, family = "NB2", method = "BFGS")
nb_disp <- countreg(y ~ x, data = sim, family = "NB2",
                    dis_param_formula_1 = ~ z, method = "BFGS")
summary(nb_disp)$coefficients[, 1:2]


## ----alpha-natural------------------------------------------------------------
d <- coef(nb_disp)
alpha <- exp(c(z0 = d[["ln(alpha):(Intercept)"]],
               z1 = d[["ln(alpha):(Intercept)"]] + d[["ln(alpha):z"]]))
round(rbind(estimate = alpha, truth = exp(c(-1, -1 + 1.2))), 3)


## ----lr-----------------------------------------------------------------------
regCompTest(nb_disp, model2 = nb_const)$statistics


## ----distribution-------------------------------------------------------------
profiles <- data.frame(x = 0, z = c(0, 1))
P <- predict(nb_disp, newdata = profiles, type = "prob", max_count = 200,
             overflow = TRUE)
k <- 0:200
mean_k <- as.vector(P[, 1:201] %*% k)
var_k <- as.vector(P[, 1:201] %*% k^2) - mean_k^2
cbind(profiles, mean = mean_k, variance = var_k, P_zero = P[, 1],
      P_above_10 = rowSums(P[, -(1:11)]))
all.equal(mean_k, as.numeric(predict(nb_disp, newdata = profiles)))


## ----pmf-plot, fig.cap = "Fitted probabilities for the two profiles: same mean, different dispersion."----
barplot(P[, 1:13], beside = TRUE, names.arg = 0:12,
        legend.text = c("z = 0", "z = 1"), xlab = "Count",
        ylab = "Probability", args.legend = list(bty = "n"))


## ----misspecified-------------------------------------------------------------
set.seed(4)
n2 <- 800
sim2 <- data.frame(x = runif(n2, -1, 1), g = rbinom(n2, 1, 0.5))
sim2$y <- rnbinom(n2, size = 1 / exp(-1),
                  mu = exp(1 + 0.5 * sim2$x + 0.8 * sim2$g))
wrong_mean <- countreg(y ~ x, data = sim2, family = "NB2",
                       dis_param_formula_1 = ~ g, method = "BFGS")
right_mean <- countreg(y ~ x + g, data = sim2, family = "NB2",
                       dis_param_formula_1 = ~ g, method = "BFGS")
disp_rows <- c("ln(alpha):(Intercept)", "ln(alpha):g")
cw <- summary(wrong_mean)$coefficients[disp_rows, 1:2]
cr <- summary(right_mean)$coefficients[disp_rows, 1:2]
round(cbind(wrong_mean = cw, right_mean = cr), 3)


## ----abalone------------------------------------------------------------------
data("abalone_rings")
set.seed(3)
ab <- abalone_rings[sample(nrow(abalone_rings), 600), ]
f_ab <- rings ~ sex + length + shell_weight
ab_pois <- countreg(f_ab, data = ab, family = "POISSON", method = "BFGS")
sum(residuals(ab_pois, type = "pearson")^2) / (nobs(ab_pois) - 5)


## ----abalone-fits-------------------------------------------------------------
ab_dpo <- countreg(f_ab, data = ab, family = "DPO", method = "BFGS")
ab_gc <- countreg(f_ab, data = ab, family = "GAMMACOUNT", method = "BFGS")
data.frame(
  model = c("Poisson", "DPO", "GAMMACOUNT"),
  logLik = c(logLik(ab_pois), logLik(ab_dpo), logLik(ab_gc)),
  AIC = c(AIC(ab_pois), AIC(ab_dpo), AIC(ab_gc)),
  shape = c(NA, exp(coef(ab_dpo)[["ln(phi)"]]),
            exp(coef(ab_gc)[["ln(shape)"]])))


## ----abalone-lr---------------------------------------------------------------
summary(ab_dpo)$lr_test[c("statistic", "df", "p.value")]


## ----abalone-moments----------------------------------------------------------
one <- ab[1, ]
P1 <- predict(ab_dpo, newdata = one, type = "prob", max_count = 60,
              overflow = TRUE)
k1 <- 0:60
m1 <- sum(P1[1, 1:61] * k1)
c(predicted_mean = predict(ab_dpo, newdata = one), mean_from_pmf = m1,
  variance_from_pmf = sum(P1[1, 1:61] * k1^2) - m1^2)

