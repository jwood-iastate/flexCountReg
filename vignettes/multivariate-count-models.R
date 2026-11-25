## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(flexCountReg)
data("washington_roads")
w <- washington_roads
w$pdo <- w$Total_crashes - w$Injury_crashes
w$injury <- w$Injury_crashes
colSums(w[, c("pdo", "injury")])

## ----fit----------------------------------------------------------------------
m <- countreg_mv(cbind(pdo, injury) ~ lnaadt + speed50 + offset(lnlength),
                 data = w, family = c(pdo = "NB2", injury = "PIG"))
summary(m)

## ----formulas-----------------------------------------------------------------
m2 <- countreg_mv(list(pdo = pdo ~ lnaadt + speed50 + offset(lnlength),
                       injury = injury ~ lnaadt + offset(lnlength)),
                  data = w, family = c(pdo = "NB2", injury = "PIG"))
AIC(m, m2)

## ----families-----------------------------------------------------------------
fams <- list(c("NB2", "NB2"), c("NB2", "PIG"), c("NB2", "PLN"),
             c("PIG", "PIG"))
sapply(fams, function(f) {
  fit <- countreg_mv(cbind(pdo, injury) ~ lnaadt + speed50 +
                       offset(lnlength), data = w, family = f)
  c(pdo = f[1], injury = f[2], AIC = round(AIC(fit), 1))
})

## ----wald---------------------------------------------------------------------
b <- coef(m)
V <- vcov(m)
a <- c("pdo:lnaadt", "injury:lnaadt")
diff <- b[a[1]] - b[a[2]]
se <- sqrt(V[a[1], a[1]] + V[a[2], a[2]] - 2 * V[a[1], a[2]])
c(difference = unname(diff), z = unname(diff / se),
  p = 2 * pnorm(-abs(unname(diff / se))))

## ----predict------------------------------------------------------------------
nd <- w[c(1, 50, 200), ]
predict(m, newdata = nd)
lapply(predict(m, newdata = nd, type = "prob", max_count = 3), round, 3)

## ----rootogram, fig.width = 7, fig.height = 3.8-------------------------------
rootogram(m)

## ----mvpln--------------------------------------------------------------------
mr <- countreg_mv(cbind(pdo, injury) ~ lnaadt + speed50 + offset(lnlength),
                  data = w, family = "POISSON", rpar_formula = ~ 1,
                  group_var = "ID", rpar_draws = 200, seed = 1)
summary(mr)
AIC(m, mr)

## ----conditional--------------------------------------------------------------
nd <- w[c(1, 50, 200), ]
predict(mr, newdata = nd, type = "joint", y = c(pdo = 1, injury = 1))
round(predict(mr, newdata = nd, type = "prob", max_count = 2)$injury, 4)
round(predict(mr, newdata = nd, type = "conditional", outcome = "injury",
              given = c(pdo = 3), max_count = 2), 4)
round(predict(mr, newdata = nd, type = "conditional", outcome = "injury",
              given = c(pdo = 0), max_count = 2), 4)

## ----posterior----------------------------------------------------------------
seg <- w[w$ID == w$ID[1], ]
round(cbind(predict(mr, newdata = seg),
            predict(mr, newdata = seg, random = "posterior")), 3)

## ----copula-------------------------------------------------------------------
mc <- countreg_mv(cbind(pdo, injury) ~ lnaadt + speed50 + offset(lnlength),
                  data = w, family = c(pdo = "NB2", injury = "PIG"),
                  copula = "clayton")
summary(mc)

## ----copula-compare-----------------------------------------------------------
fams <- c(pdo = "NB2", injury = "PIG")
f <- cbind(pdo, injury) ~ lnaadt + speed50 + offset(lnlength)
cop <- lapply(c("gaussian", "frank", "clayton", "gumbel"), function(cf)
  countreg_mv(f, data = w, family = fams, copula = cf))
round(c(independent = AIC(m), gaussian = AIC(cop[[1]]),
        frank = AIC(cop[[2]]), clayton = AIC(cop[[3]]),
        gumbel = AIC(cop[[4]]), MVPLN = AIC(mr)))

## ----both, eval = FALSE-------------------------------------------------------
# mb <- countreg_mv(cbind(pdo, injury) ~ lnaadt + speed50 + offset(lnlength),
#                   data = w, family = "POISSON", rpar_formula = ~ 1,
#                   group_var = "ID", rpar_draws = 200, seed = 1,
#                   copula = "frank")

