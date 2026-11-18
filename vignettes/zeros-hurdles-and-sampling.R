## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4)
has_pscl <- requireNamespace("pscl", quietly = TRUE)


## ----setup--------------------------------------------------------------------
library(flexCountReg)


## ----simulate-----------------------------------------------------------------
set.seed(12)
n <- 1000
x <- rnorm(n)
mu <- exp(0.3 + 0.5 * x)
y_nb <- rnbinom(n, size = 1, mu = mu)
y_zi <- ifelse(runif(n) < 0.25, 0, rnbinom(n, size = 1, mu = mu))
# hurdle: positive with probability 0.55, then a zero-truncated NB2 count
draw_pos <- function(m) { repeat { v <- rnbinom(1, size = 1, mu = m)
  if (v > 0) return(v) } }
y_h <- ifelse(runif(n) < 0.55, vapply(mu, draw_pos, 0), 0)
sims <- data.frame(x, y_nb, y_zi, y_h)
sapply(sims[-1], function(y) mean(y == 0))


## ----sim-fits-----------------------------------------------------------------
fit_three <- function(resp) {
  f <- reformulate("x", resp)
  list(NB2 = countreg(f, data = sims, family = "NB2", method = "BFGS"),
       ZINB = suppressWarnings(
         countreg(f, data = sims, family = "NB2", zero_formula = ~ 1,
                  method = "BFGS")),
       hurdle = countreg(f, data = sims, family = "NB2", zero_formula = ~ 1,
                         zero_type = "hurdle", method = "BFGS"))
}
sim_fits <- lapply(c(ordinary = "y_nb", zero_inflated = "y_zi",
                     hurdle = "y_h"), fit_three)
aic <- t(sapply(sim_fits, function(m) sapply(m, AIC)))
round(aic, 1)


## ----sim-best, include = FALSE------------------------------------------------
best <- colnames(aic)[apply(aic, 1, which.min)]
gap_zi <- aic["zero_inflated", "NB2"] - aic["zero_inflated", "ZINB"]


## ----biochem, eval = has_pscl-------------------------------------------------
e <- new.env()
utils::data("bioChemists", package = "pscl", envir = e)
b <- e$bioChemists
c(n = nrow(b), share_zero = mean(b$art == 0), mean = mean(b$art),
  variance = var(b$art))


## ----biochem-skip, echo = FALSE, results = "asis", eval = !has_pscl-----------
# cat("*The pscl package is not installed, so this section is not run.*")


## ----baselines, eval = has_pscl-----------------------------------------------
f <- art ~ fem + mar + kid5 + phd + ment
pois <- countreg(f, data = b, family = "POISSON", method = "BFGS")
nb2 <- countreg(f, data = b, family = "NB2", method = "BFGS")
zeros <- c(observed = mean(b$art == 0),
           poisson = mean(predict(pois, type = "prob", max_count = 0)[, 1]),
           nb2 = mean(predict(nb2, type = "prob", max_count = 0)[, 1]))
round(zeros, 3)


## ----zero-models, eval = has_pscl---------------------------------------------
zinb_const <- countreg(f, data = b, family = "NB2", zero_formula = ~ 1,
                       method = "BFGS")
zinb <- countreg(f, data = b, family = "NB2",
                 zero_formula = ~ ment + kid5, method = "BFGS")
hnb <- countreg(f, data = b, family = "NB2",
                zero_formula = ~ fem + mar + kid5 + phd + ment,
                zero_type = "hurdle", method = "BFGS")


## ----zinb-const, eval = has_pscl----------------------------------------------
summary(zinb_const)$coefficients["Zero:(Intercept)", 1:2]
plogis(coef(zinb_const)[["Zero:(Intercept)"]])


## ----zero-coefs, eval = has_pscl----------------------------------------------
zc <- function(m) {
  co <- summary(m)$coefficients
  co[grep("^Zero:", rownames(co)), 1:2]
}
zc(zinb)
zc(hnb)


## ----pscl-check, eval = has_pscl----------------------------------------------
z_ref <- pscl::zeroinfl(art ~ fem + mar + kid5 + phd + ment | ment + kid5,
                        data = b, dist = "negbin")
h_ref <- pscl::hurdle(art ~ fem + mar + kid5 + phd + ment, data = b,
                      dist = "negbin")
c(zinb = as.numeric(logLik(zinb)), pscl_zinb = as.numeric(logLik(z_ref)),
  hurdle = as.numeric(logLik(hnb)), pscl_hurdle = as.numeric(logLik(h_ref)))


## ----accounting, eval = has_pscl----------------------------------------------
models <- list(Poisson = pois, NB2 = nb2, ZINB = zinb, hurdle = hnb)
probs <- sapply(models, function(m)
  colMeans(predict(m, type = "prob", max_count = 5)))
acct <- rbind(round(cbind(observed = tabulate(b$art + 1, 6) / nrow(b),
                          probs), 3),
              mean_count = round(c(mean(b$art),
                                   sapply(models, function(m)
                                     mean(predict(m)))), 3))
acct


## ----compare, eval = has_pscl-------------------------------------------------
data.frame(model = names(models),
           logLik = sapply(models, function(m) as.numeric(logLik(m))),
           AIC = sapply(models, AIC), row.names = NULL)

