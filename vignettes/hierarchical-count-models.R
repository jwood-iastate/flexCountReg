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

