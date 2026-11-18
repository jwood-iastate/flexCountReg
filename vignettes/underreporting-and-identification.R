## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4)


## ----setup--------------------------------------------------------------------
library(flexCountReg)


## ----intercept-only-----------------------------------------------------------
set.seed(1)
toy <- data.frame(x = rnorm(200))
toy$y <- rpois(200, exp(0.5 + 0.3 * toy$x))
res <- tryCatch(
  countreg(y ~ x, data = toy, family = "NB2", underreport_formula = ~ 1),
  error = function(e) conditionMessage(e))
cat(strwrap(res), sep = "\n")


## ----sim----------------------------------------------------------------------
simulate_reports <- function(n, seed, gamma) {
  set.seed(seed)
  x <- rnorm(n)
  w <- runif(n)
  mu <- exp(1 + 0.5 * x)
  p <- plogis(gamma[1] + gamma[2] * w)
  y_true <- rnbinom(n, size = 2, mu = mu)
  data.frame(x, w, y_true, y = rbinom(n, y_true, p), mu, p)
}
designs <- list(weak = c(1.5, -2), strong = c(3, -5))
sapply(designs, function(g) round(plogis(g[1] + g[2] * c(w0 = 0, w1 = 1)), 2))


## ----fit-one------------------------------------------------------------------
fits <- lapply(designs, function(g) {
  d <- simulate_reports(2000, seed = 8, gamma = g)
  countreg(y ~ x, data = d, family = "NB2", underreport_formula = ~ w,
           method = "BFGS")
})
lapply(fits, function(m) round(summary(m)$coefficients[, 1:2], 3))


## ----replicates---------------------------------------------------------------
recover_one <- function(seed, gamma) {
  d <- simulate_reports(2000, seed, gamma)
  warned <- FALSE
  m <- tryCatch(withCallingHandlers(
    countreg(y ~ x, data = d, family = "NB2", underreport_formula = ~ w,
             method = "BFGS"),
    warning = function(w) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    }), error = function(e) NULL)
  if (is.null(m) || !isTRUE(m$model$converged)) {
    return(c(NA, NA, NA, NA, converged = 0, warned = warned))
  }
  co <- coef(m)[c("(Intercept)", "x", "Underreporting:(Intercept)",
                  "Underreporting:w")]
  c(co, converged = 1, warned = warned)
}
recovery <- lapply(designs, function(g)
  t(sapply(101:120, recover_one, gamma = g)))
truth <- function(g) c(1, 0.5, g)
recovery_table <- do.call(rbind, lapply(names(designs), function(k) {
  r <- recovery[[k]]
  ok <- r[, "converged"] == 1
  data.frame(design = k,
             parameter = c("(Intercept)", "x", "Rep:(Intercept)", "Rep:w"),
             truth = truth(designs[[k]]),
             mean = colMeans(r[ok, 1:4]),
             median = apply(r[ok, 1:4], 2, median),
             sd = apply(r[ok, 1:4], 2, sd),
             converged = sum(ok), warned = sum(r[, "warned"]),
             row.names = NULL)
}))
recovery_table[, 4:6] <- round(recovery_table[, 4:6], 3)
recovery_table


## ----recovery-summary, include = FALSE----------------------------------------
warned <- sapply(recovery, function(r) sum(r[, "warned"]))
sd_rep <- sapply(recovery, function(r) sd(r[, 3], na.rm = TRUE))


## ----predictions, fig.height = 4.5--------------------------------------------
m <- fits$strong
grid <- data.frame(x = 0, w = seq(0, 1, length.out = 50))
pred <- data.frame(grid,
                   reported = predict(m, newdata = grid),
                   true = predict(m, newdata = grid, type = "true"),
                   p = predict(m, newdata = grid, type = "prob_reported"))
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(pred$w, pred$true, type = "l", lwd = 2, ylim = c(0, 3.2),
     xlab = "w", ylab = "Expected count at x = 0", main = "Counts")
lines(pred$w, pred$reported, lwd = 2, col = 2)
abline(h = exp(1), lty = 2)
lines(pred$w, exp(1) * plogis(3 - 5 * pred$w), lty = 2, col = 2)
legend("bottomleft", c("true", "reported"), col = 1:2, lwd = 2, bty = "n")
plot(pred$w, pred$p, type = "l", lwd = 2, col = 4, ylim = c(0, 1),
     xlab = "w", ylab = "Probability", main = "Reporting probability")
lines(pred$w, plogis(3 - 5 * pred$w), lty = 2, col = 4)
par(op)


## ----sensitivity--------------------------------------------------------------
d <- simulate_reports(2000, seed = 8, gamma = designs$strong)
d$log_p <- log(d$p)
models <- list(
  logit = fits$strong,
  probit = countreg(y ~ x, data = d, family = "NB2",
                    underreport_formula = ~ w, underreport_family = "probit",
                    method = "BFGS"),
  ignore = countreg(y ~ x + w, data = d, family = "NB2", method = "BFGS"),
  known_p = countreg(y ~ x + offset(log_p), data = d, family = "NB2",
                     method = "BFGS"))
at <- data.frame(x = 0, w = c(0, 1), log_p = 0)
sens <- t(sapply(models, function(m) {
  has_rep <- any(grepl("^Underreporting:", names(coef(m))))
  c(beta_x = coef(m)[["x"]],
    true_mean_x0 = if (has_rep) predict(m, newdata = at[1, ], type = "true")
                   else NA,
    p_w0 = if (has_rep) predict(m, newdata = at[1, ],
                                type = "prob_reported") else NA,
    p_w1 = if (has_rep) predict(m, newdata = at[2, ],
                                type = "prob_reported") else NA,
    AIC = AIC(m))
}))
round(sens, 3)


## ----both---------------------------------------------------------------------
both <- withCallingHandlers(
  countreg(y ~ x + w, data = d, family = "NB2", underreport_formula = ~ w,
           method = "BFGS"),
  warning = function(w) {
    message("Warning: ", conditionMessage(w))
    invokeRestart("muffleWarning")
  })
round(summary(both)$coefficients[, 1:2], 3)
round(cov2cor(vcov(both))[c("(Intercept)", "w"),
                          c("Underreporting:(Intercept)",
                            "Underreporting:w")], 3)

