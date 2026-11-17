## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")


## ----setup--------------------------------------------------------------------
library(flexCountReg)


## ----epil---------------------------------------------------------------------
epil <- MASS::epil
c(rows = nrow(epil), patients = length(unique(epil$subject)))
table(table(epil$subject))      # rows per patient
# treatment does not change within a patient
all(tapply(epil$trt, epil$subject, function(v) length(unique(v))) == 1)


## ----epil-fits----------------------------------------------------------------
f_epil <- y ~ trt + lbase + lage + V4
pooled <- countreg(f_epil, data = epil, family = "NB2", method = "BFGS")
pooled_robust <- countreg(f_epil, data = epil, family = "NB2",
                          method = "BFGS", stderr = "robust")
panel <- renb(f_epil, group_var = "subject", data = epil, method = "BFGS")
se_trt <- function(m) sqrt(diag(vcov(m)))[["trtprogabide"]]
data.frame(
  model = c("pooled NB2", "pooled NB2, robust SE", "RENB (panel)"),
  treatment_effect = c(coef(pooled)[["trtprogabide"]],
                       coef(pooled_robust)[["trtprogabide"]],
                       coef(panel)[["trtprogabide"]]),
  std_error = c(se_trt(pooled), se_trt(pooled_robust), se_trt(panel)))


## ----epil-glance--------------------------------------------------------------
do.call(rbind, lapply(list(pooled = pooled, RENB = panel), function(m)
  glance(m)[, c("logLik", "AIC", "BIC", "nobs", "n_groups")]))


## ----roads--------------------------------------------------------------------
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


## ----re-flat------------------------------------------------------------------
tidy(re_models$PoisLindRE)


## ----se-----------------------------------------------------------------------
re_robust <- poisLind.re(Total_crashes ~ lnaadt + speed50, group_var = "ID",
                         data = washington_roads, offset = "lnlength",
                         method = "BFGS", stderr = "robust")
tidy(re_robust)

