## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(flexCountReg)
data("washington_roads")

## ----define-------------------------------------------------------------------
oinb <- new_count_family(
  name = "OINB",
  params = list(alpha = count_param("log", start = -1),
                omega = count_param("logit", start = -3)),
  logpmf = function(y, mu, par) {
    lnb <- dnbinom(y, size = 1 / par$alpha, mu = mu, log = TRUE)
    ifelse(y == 1,
           log((1 - par$omega) * exp(lnb) + par$omega),
           log1p(-par$omega) + lnb)
  },
  mean_is_mu = FALSE,
  mean = function(mu, par) (1 - par$omega) * mu + par$omega,
  variance = function(mu, par) {
    m <- (1 - par$omega) * mu + par$omega
    ey2 <- (1 - par$omega) * (mu + par$alpha * mu^2 + mu^2) + par$omega
    ey2 - m^2
  })
oinb

## ----check--------------------------------------------------------------------
chk <- check_count_family(oinb, mu = c(0.2, 2, 20),
                          par = list(alpha = c(0.2, 2),
                                     omega = c(0.05, 0.5)))

## ----check-bad----------------------------------------------------------------
wrong <- oinb
wrong$variance <- function(mu, par) (1 - par$omega) * (mu + par$alpha * mu^2)
chk <- check_count_family(wrong, mu = c(0.2, 2, 20),
                          par = list(alpha = c(0.2, 2),
                                     omega = c(0.05, 0.5)))

## ----fit----------------------------------------------------------------------
f <- Total_crashes ~ lnaadt + lnlength + speed50
m_oinb <- countreg(f, data = washington_roads, family = oinb,
                   method = "BFGS")
summary(m_oinb)

## ----compare------------------------------------------------------------------
m_nb2 <- countreg(f, data = washington_roads, family = "NB2",
                  method = "BFGS")
c(NB2 = AIC(m_nb2), OINB = AIC(m_oinb))

## ----simulated----------------------------------------------------------------
set.seed(42)
n <- 2000
sim <- data.frame(x = rnorm(n))
sim$y <- rnbinom(n, size = 2, mu = exp(0.8 + 0.5 * sim$x))
sim$y[runif(n) < 0.2] <- 1
s_oinb <- countreg(y ~ x, data = sim, family = oinb, method = "BFGS")
s_nb2 <- countreg(y ~ x, data = sim, family = "NB2", method = "BFGS")
c(NB2 = AIC(s_nb2), OINB = AIC(s_oinb))
# omega (true 0.2) and alpha (true 0.5)
c(omega = plogis(coef(s_oinb)[["logit(omega)"]]),
  alpha = exp(coef(s_oinb)[["ln(alpha)"]]))

## ----dis-formula--------------------------------------------------------------
m_oinb2 <- countreg(f, data = washington_roads, family = oinb,
                    dis_param_formula_2 = ~ speed50, method = "BFGS")
coef(m_oinb2)[c("logit(omega):(Intercept)", "logit(omega):speed50")]

## ----methods------------------------------------------------------------------
nd <- washington_roads[1:3, ]
predict(m_oinb, newdata = nd)                      # E[Y]
round(predict(m_oinb, newdata = nd, type = "prob")[, 1:5], 3)
head(residuals(m_oinb, type = "quantile"))
sims <- simulate(m_oinb, nsim = 2, seed = 1)

## ----frequencies--------------------------------------------------------------
expected <- function(m) {
  colSums(predict(m, newdata = washington_roads, type = "prob")[, 1:5])
}
rbind(observed = tabulate(washington_roads$Total_crashes + 1, 5),
      NB2 = round(expected(m_nb2), 1),
      OINB = round(expected(m_oinb), 1))

## ----rp-----------------------------------------------------------------------
m_rp <- countreg_rp(Total_crashes ~ lnaadt + lnlength,
                    rpar_formula = ~ -1 + speed50,
                    data = washington_roads, family = oinb,
                    group_var = "ID", ndraws = 50)
coef(m_rp)

## ----save---------------------------------------------------------------------
path <- tempfile(fileext = ".rds")
saveRDS(m_oinb, path)
m_again <- readRDS(path)
predict(m_again, newdata = nd)

