## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(flexCountReg)
data("washington_roads")
length(unique(washington_roads$ID))   # road segments
nrow(washington_roads)                # segment-years

## ----rp-----------------------------------------------------------------------
rp <- countreg_rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                  data = washington_roads, family = "NB2",
                  rpardists = c(speed50 = "n"), offset = "lnlength",
                  group_var = "ID", ndraws = 100, method = "BFGS")
summary(rp)

## ----het----------------------------------------------------------------------
rp_het <- countreg_rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                      het_mean_formula = ~ ShouldWidth04,
                      data = washington_roads, offset = "lnlength",
                      group_var = "ID", ndraws = 100, method = "BFGS")
tidy(rp_het)

## ----predict------------------------------------------------------------------
nd <- washington_roads[washington_roads$ID %in% unique(washington_roads$ID)[1:2], ]
cbind(observed = nd$Total_crashes,
      exact = predict(rp, newdata = nd, method = "Exact"),
      simulated = predict(rp, newdata = nd, method = "Simulated"),
      individual = predict(rp, newdata = nd, method = "Individual"))

## ----eb-----------------------------------------------------------------------
eb <- predict(rp, newdata = washington_roads, method = "Individual")
summary(attr(eb, "eb_coefficients")[, "speed50"])

## ----re-----------------------------------------------------------------------
re_models <- list(
  RENB = renb(Total_crashes ~ lnaadt + speed50, group_var = "ID",
              data = washington_roads, offset = "lnlength", method = "BFGS"),
  PoisLindRE = poisLind.re(Total_crashes ~ lnaadt + speed50, group_var = "ID",
                           data = washington_roads, offset = "lnlength",
                           method = "BFGS"),
  NM2 = negmultinom(Total_crashes ~ lnaadt + speed50, group_var = "ID",
                    data = washington_roads, offset = "lnlength")
)
do.call(rbind, lapply(names(re_models), function(nm) {
  cbind(model = nm,
        glance(re_models[[nm]])[, c("logLik", "AIC", "BIC", "nobs",
                                    "n_groups")])
}))

## ----se-----------------------------------------------------------------------
re_robust <- poisLind.re(Total_crashes ~ lnaadt + speed50, group_var = "ID",
                         data = washington_roads, offset = "lnlength",
                         method = "BFGS", stderr = "robust")
tidy(re_robust)

