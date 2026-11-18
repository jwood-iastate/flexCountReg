## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4)


## ----setup--------------------------------------------------------------------
library(flexCountReg)


## ----data---------------------------------------------------------------------
epil <- MASS::epil
table(table(epil$subject))     # rows per patient
# the periods are 1 to 4 in order for every patient
all(tapply(epil$period, epil$subject, function(p) identical(p, 1:4)))
epil$progabide <- as.numeric(epil$trt == "progabide")


## ----fit----------------------------------------------------------------------
fe <- y ~ -1 + progabide + lbase + lage + V4
rp <- countreg_rp(fe, rpar_formula = ~ 1, data = epil, family = "NB2",
                  group_var = "subject", ndraws = 200, method = "BFGS")
rp$model$converged
round(summary(rp)$coefficients[, 1:2], 3)


## ----one-patient--------------------------------------------------------------
p1 <- epil[epil$subject == 1, ]
data.frame(period = p1$period, y = p1$y,
           new_patient = predict(rp, newdata = p1, method = "Exact"),
           simulated = predict(rp, newdata = p1, method = "Simulated"),
           empirical_bayes = predict(rp, newdata = p1, method = "Individual"))


## ----unseen-------------------------------------------------------------------
p_new <- p1
p_new$subject <- 999
all.equal(predict(rp, newdata = p_new, method = "Exact"),
          predict(rp, newdata = p1, method = "Exact"))


## ----profiles, fig.height = 5-------------------------------------------------
show <- c(2, 10, 25, 49)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (s in show) {
  d <- epil[epil$subject == s, ]
  pop <- predict(rp, newdata = d, method = "Exact")
  eb <- predict(rp, newdata = d, method = "Individual")
  plot(d$period, d$y, pch = 19, ylim = range(0, d$y, pop, eb),
       xlab = "Period", ylab = "Seizures", main = paste("Patient", s),
       xaxt = "n")
  axis(1, at = 1:4)
  lines(d$period, pop, lty = 2)
  lines(d$period, eb, lwd = 2, col = 4)
}
par(op)


## ----grouped-cv---------------------------------------------------------------
ids <- unique(epil$subject)
set.seed(7)
fold_of <- setNames(sample(rep(1:5, length.out = length(ids))), ids)
epil$fold <- fold_of[as.character(epil$subject)]
held_out <- do.call(rbind, lapply(1:5, function(k) {
  train <- epil[epil$fold != k, ]
  test <- epil[epil$fold == k, ]
  m_rp <- countreg_rp(fe, rpar_formula = ~ 1, data = train, family = "NB2",
                      group_var = "subject", ndraws = 200, method = "BFGS")
  m_pool <- countreg(y ~ progabide + lbase + lage + V4, data = train,
                     family = "NB2", method = "BFGS")
  data.frame(fold = k, subject = test$subject, y = test$y,
             random_intercept = predict(m_rp, newdata = test,
                                        method = "Exact"),
             pooled = predict(m_pool, newdata = test),
             converged = m_rp$model$converged)
}))
all(held_out$converged)
sapply(held_out[c("random_intercept", "pooled")], function(p)
  c(MAE = mean(abs(held_out$y - p)), RMSE = sqrt(mean((held_out$y - p)^2))))


## ----in-sample----------------------------------------------------------------
c(empirical_bayes_in_sample = mean(abs(epil$y - predict(rp, newdata = epil,
                                                         method = "Individual"))),
  population_in_sample = mean(abs(epil$y - predict(rp, newdata = epil,
                                                    method = "Exact"))))


## ----forecast-limit-----------------------------------------------------------
p1_hist <- p1
p1_hist$y[4] <- NA
predict(rp, newdata = p1_hist, method = "Individual")


## ----totals-------------------------------------------------------------------
co <- coef(rp)
patient <- epil[epil$subject == 1, ]
x_fixed <- as.vector(model.matrix(~ -1 + progabide + lbase + lage + V4,
                                  patient) %*%
                       co[c("progabide", "lbase", "lage", "V4")])
size <- 1 / exp(co[["ln(alpha)"]])           # NB2: Var = mu + alpha mu^2
mu_b <- co[["(Intercept):Mean"]]
sd_b <- abs(co[["(Intercept):St.Dev"]])
set.seed(42)
R <- 20000
# joint: one intercept per simulated patient, shared by the four periods
b <- rnorm(R, mu_b, sd_b)
joint <- rowSums(sapply(1:4, function(t)
  rnbinom(R, size = size, mu = exp(b + x_fixed[t]))))
# independent: a new intercept for every period (wrong for one patient)
indep <- rowSums(sapply(1:4, function(t)
  rnbinom(R, size = size, mu = exp(rnorm(R, mu_b, sd_b) + x_fixed[t]))))
rbind(joint = c(mean = mean(joint), sd = sd(joint),
                quantile(joint, c(0.05, 0.95)), P_over_40 = mean(joint > 40)),
      independent = c(mean = mean(indep), sd = sd(indep),
                      quantile(indep, c(0.05, 0.95)),
                      P_over_40 = mean(indep > 40)))
sum(predict(rp, newdata = patient, method = "Exact"))

