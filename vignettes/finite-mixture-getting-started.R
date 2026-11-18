## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4)

## ----setup--------------------------------------------------------------------
library(flexCountReg)

## ----simulate-----------------------------------------------------------------
set.seed(2026)
n <- 500
true_class <- ifelse(runif(n) < 0.6, 1, 2)
y <- rpois(n, lambda = c(1, 5)[true_class])
d <- data.frame(y = y)
table(y)
c(mean = mean(y), variance = var(y))

## ----oracle-------------------------------------------------------------------
p <- c(0.6, 0.4)
m_k <- c(1, 5)
m <- sum(p * m_k)
V <- sum(p * (m_k + m_k^2)) - m^2
c(mean = m, variance = V)

## ----fit----------------------------------------------------------------------
fit2 <- countreg_mix(y ~ 1, data = d, family = "POISSON", components = 2,
                     n_starts = 5, seed = 42)
summary(fit2)

## ----estimates----------------------------------------------------------------
cf <- coef(fit2)
exp(cf[c("Component1:Mean:(Intercept)", "Component2:Mean:(Intercept)")])
plogis(cf["Mixing:Component1:(Intercept)"])
fit2$model$classes[, c("label", "prior_share", "effective_size",
                       "mean_count")]

## ----fitted-moments-----------------------------------------------------------
nd <- data.frame(y = 0)
c(mean = unname(predict(fit2, newdata = nd, type = "response")),
  variance = unname(predict(fit2, newdata = nd, type = "variance")))
p_hat <- fit2$model$classes$prior_share
m_hat <- predict(fit2, newdata = nd, type = "component")[1, ]
c(mean = sum(p_hat * m_hat),
  variance = sum(p_hat * (m_hat + m_hat^2)) - sum(p_hat * m_hat)^2)

## ----starts-------------------------------------------------------------------
fit2$model$starts[, c("start", "strategy", "em_iterations", "em_status",
                      "loglik_start", "polished", "loglik", "reached_best",
                      "selected")]

## ----freq-table---------------------------------------------------------------
tab <- as.data.frame(table(y), stringsAsFactors = FALSE)
tab$y <- as.numeric(tab$y)
head(tab)
fit2_tab <- countreg_mix(y ~ 1, data = tab, weights = "Freq",
                         family = "POISSON", components = 2,
                         n_starts = 5, seed = 42)
c(raw = logLik(fit2), table = logLik(fit2_tab))
all.equal(coef(fit2), coef(fit2_tab), tolerance = 1e-5)
c(raw = BIC(fit2), table = BIC(fit2_tab))

## ----pmf-plot-----------------------------------------------------------------
max_count <- 12
P <- predict(fit2, newdata = nd, type = "prob", max_count = max_count,
             overflow = TRUE)
round(P, 4)
grid <- 0:max_count
contrib <- sapply(1:2, function(k) p_hat[k] * dpois(grid, m_hat[k]))
all.equal(rowSums(contrib), unname(P[1, seq_along(grid)]))

obs <- as.vector(table(factor(pmin(y, max_count + 1),
                              levels = 0:(max_count + 1)))) / n
plot(c(grid, max_count + 1), P[1, ], type = "h", lwd = 6, col = "grey70",
     xaxt = "n", xlab = "Count", ylab = "Probability",
     ylim = c(0, max(P, obs, contrib)),
     main = "Fitted mixture (bars), components (lines), data (points)")
axis(1, at = c(grid, max_count + 1),
     labels = c(grid, paste0(">", max_count)))
lines(grid, contrib[, 1], type = "b", pch = 1, lty = 1)
lines(grid, contrib[, 2], type = "b", pch = 2, lty = 2)
points(c(grid, max_count + 1), obs, pch = 19, col = "firebrick")
legend("topright", c("mixture", "pi_1 f_1", "pi_2 f_2", "observed"),
       lty = c(1, 1, 2, NA), lwd = c(6, 1, 1, NA), pch = c(NA, 1, 2, 19),
       col = c("grey70", "black", "black", "firebrick"), bty = "n")

## ----posterior----------------------------------------------------------------
post <- component_probs(fit2, conditioning = "posterior")
round(cbind(y = y, post)[order(y)[c(1, 250, 350, 450, 500)], ], 3)

## ----rootogram----------------------------------------------------------------
rootogram(fit2, main = "Two-component Poisson mixture")

## ----compare------------------------------------------------------------------
fit1 <- countreg_mix(y ~ 1, data = d, family = "POISSON", components = 1)
nb2 <- countreg(y ~ 1, data = d, family = "NB2", method = "BFGS")
fit3 <- countreg_mix(y ~ 1, data = d, family = "POISSON", components = 3,
                     n_starts = 5, seed = 42)
data.frame(model = c("Poisson", "NB2", "2 Poisson components",
                     "3 Poisson components"),
           logLik = sapply(list(fit1, nb2, fit2, fit3), logLik),
           df = sapply(list(fit1, nb2, fit2, fit3),
                       function(m) length(coef(m))),
           AIC = sapply(list(fit1, nb2, fit2, fit3), AIC),
           BIC = sapply(list(fit1, nb2, fit2, fit3), BIC))
fit3$model$classes[, c("label", "prior_share", "mean_count")]

## ----fixed--------------------------------------------------------------------
fit_fixed <- countreg_mix(y ~ 1, data = d, family = "POISSON",
                          components = 2, fixed_mixing = c(0.6, 0.4),
                          n_starts = 3, seed = 1)
coef(fit_fixed)
fit_fixed$model$classes[, c("label", "prior_share", "posterior_share")]
c(df_estimated = length(coef(fit2)), df_fixed = length(coef(fit_fixed)))

