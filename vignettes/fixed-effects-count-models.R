## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(flexCountReg)
data("washington_roads")

## ----fe-----------------------------------------------------------------------
fe <- countreg_fe(Total_crashes ~ lnaadt + speed50 + offset(lnlength),
                  fixef = ~ ID, data = washington_roads)
summary(fe)

## ----dropped------------------------------------------------------------------
length(fe$model$dropped$ID)     # segments dropped
fe$model$not_identified         # regressors removed

## ----vcov---------------------------------------------------------------------
sqrt(diag(vcov(fe)))                      # clustered by segment
sqrt(diag(vcov(fe, type = "hessian")))    # Poisson model-based

## ----twoway-------------------------------------------------------------------
fe2 <- countreg_fe(Total_crashes ~ lnaadt + offset(lnlength),
                   fixef = ~ ID + Year, data = washington_roads)
summary(fe2)
fe2$model$fixef$fe2                       # year effects

## ----compare------------------------------------------------------------------
f <- Total_crashes ~ lnaadt + offset(lnlength)
pooled <- countreg(f, data = washington_roads, family = "POISSON",
                   method = "BFGS", cluster = ~ ID)
re <- repois(f, data = washington_roads, group_var = "ID")
est <- function(m) c(estimate = coef(m)[["lnaadt"]],
                     se = sqrt(diag(vcov(m)))[["lnaadt"]])
round(rbind(pooled = est(pooled), random_effects = est(re),
            fixed_effects = est(fe2)), 3)

## ----nb2----------------------------------------------------------------------
set.seed(11)
G <- 300; n_t <- 5
sim <- data.frame(id = rep(seq_len(G), each = n_t),
                  t = rep(seq_len(n_t), G))
a_i <- rnorm(G, 0, 0.7)
sim$x1 <- rnorm(nrow(sim)) + 0.5 * a_i[sim$id]
sim$x2 <- rnorm(nrow(sim))
sim$y <- rnbinom(nrow(sim), size = 2,
                 mu = exp(0.5 + a_i[sim$id] + 0.6 * sim$x1 - 0.4 * sim$x2))

fe_nb <- countreg_fe(y ~ x1 + x2, fixef = ~ id, data = sim, family = "NB2")
fe_nb_jk <- countreg_fe(y ~ x1 + x2, fixef = ~ id, data = sim,
                        family = "NB2", jackknife = TRUE, time = "t")
round(rbind(truth = c(0.6, -0.4, log(0.5)),
            NB2_FE = coef(fe_nb),
            jackknife = coef(fe_nb_jk)), 3)

## ----predict------------------------------------------------------------------
nd <- washington_roads[washington_roads$ID %in%
                         names(fe2$model$fixef$fe1)[1:2], ]
predict(fe2, newdata = nd)

new_site <- nd[1, ]
new_site$ID <- "not-in-data"
predict(fe2, newdata = new_site)
predict(fe2, newdata = new_site, newgroup = "mean_fe")

