## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(flexCountReg)
data("stream_salamanders")
s <- stream_salamanders
s$sample <- factor(s$sample)

## ----fit----------------------------------------------------------------------
h <- countreg_hier(count ~ spp + mined + (1 | site), data = s,
                   family = "NB2")
summary(h)

## ----lrtest-------------------------------------------------------------------
re_lrtest(h)

## ----ranef--------------------------------------------------------------------
head(ranef(h)$site)
VarCorr(h)

## ----crossed------------------------------------------------------------------
h2 <- countreg_hier(count ~ spp + mined + (1 | site) + (1 | sample),
                    data = s, family = "POISSON")
VarCorr(h2)

## ----slopes-------------------------------------------------------------------
h3 <- countreg_hier(count ~ cover + mined + (cover || site), data = s,
                    family = "NB2")
VarCorr(h3)

## ----predict------------------------------------------------------------------
nd <- s[c(1, 30), ]
predict(h2, newdata = nd)                         # site and sample
predict(h2, newdata = nd, re.form = ~ (1 | site)) # site only
predict(h2, newdata = nd, re.form = NA)           # population level

## ----marginal-----------------------------------------------------------------
predict(h2, newdata = nd, re.form = NA, marginal = TRUE)

## ----history------------------------------------------------------------------
data("washington_roads")
w <- washington_roads
past <- w[w$Year < 2018, ]
h_w <- countreg_hier(Total_crashes ~ lnaadt + speed50 + offset(lnlength) +
                       (1 | ID), data = past, family = "POISSON")
next_year <- w[w$Year == 2018, ]
last_year <- past[past$Year == 2017, ]
head(cbind(
  all_years = predict(h_w, newdata = next_year),
  last_year = predict(h_w, newdata = next_year, history = last_year),
  observed  = next_year$Total_crashes))

## ----profile------------------------------------------------------------------
confint(h, "sd_(Intercept)|site")                      # Wald
confint(h, "sd_(Intercept)|site", method = "profile")  # profile

## ----smooth-------------------------------------------------------------------
hs <- countreg_hier(count ~ mined + s(cover) + (1 | site), data = s,
                    family = "POISSON")
hs$model$smooth_table[, c("term", "edf", "p.value")]

## ----gamma--------------------------------------------------------------------
hg <- countreg_hier(count ~ spp + mined + (1 | site), data = s,
                    family = "POISSON", re_dist = "gamma")
rp <- repois(count ~ spp + mined, group_var = "site", data = s,
             method = "BFGS")
c(simulated = logLik(hg), closed_form = logLik(rp))

## ----cre----------------------------------------------------------------------
hc <- countreg_hier(Total_crashes ~ lnaadt + speed50 + offset(lnlength) +
                      (1 | ID), data = w, family = "POISSON", cre = ~ lnaadt)
cre_test(hc)

