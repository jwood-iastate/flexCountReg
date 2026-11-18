## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4)
has_aer <- requireNamespace("AER", quietly = TRUE)


## ----setup--------------------------------------------------------------------
library(flexCountReg)


## ----quine--------------------------------------------------------------------
q <- MASS::quine
nb <- countreg(Days ~ Eth + Sex + Age + Lrn, data = q, family = "NB2",
               method = "BFGS")
summary(nb)$coefficients[, 1:2]


## ----quine-profile------------------------------------------------------------
ref <- data.frame(Eth = "A", Sex = "F", Age = "F0", Lrn = "AL")
ref[] <- Map(function(v, nm) factor(v, levels = levels(q[[nm]])), ref,
             names(ref))
c(reference_profile = predict(nb, newdata = ref),
  average_over_children = mean(predict(nb)))


## ----eth-contrast-------------------------------------------------------------
d_A <- d_N <- q
d_A$Eth <- factor("A", levels = levels(q$Eth))
d_N$Eth <- factor("N", levels = levels(q$Eth))
diff_by_child <- predict(nb, newdata = d_N) - predict(nb, newdata = d_A)
c(average_difference = mean(diff_by_child),
  smallest = min(diff_by_child), largest = max(diff_by_child))


## ----margeff-quine------------------------------------------------------------
me <- margEffTable(nb, vars = c("Eth", "Sex", "Lrn"),
                   indicator = "discrete_change", tableType = "tibble")
me[, c("term", "effect_metric", "estimate", "std_error", "lower_ci",
       "upper_ci")]
all.equal(me$estimate[me$term == "Eth"], round(mean(diff_by_child), 3))


## ----age-refused--------------------------------------------------------------
tryCatch(margEffTable(nb, vars = "Age"),
         error = function(e) conditionMessage(e))


## ----age-contrast-------------------------------------------------------------
by_age <- sapply(levels(q$Age), function(a) {
  d <- q
  d$Age <- factor(a, levels = levels(q$Age))
  mean(predict(nb, newdata = d))
})
round(by_age, 2)
round(by_age - by_age[["F0"]], 2)   # differences from F0, in days


## ----warp---------------------------------------------------------------------
w <- datasets::warpbreaks
wb <- countreg(breaks ~ wool * tension, data = w, family = "NB2",
               method = "BFGS")
grid <- expand.grid(wool = levels(w$wool), tension = levels(w$tension))
cbind(grid, round(predict(wb, newdata = grid, interval = "confidence"), 2))


## ----warp-plot, fig.cap = "Expected breaks by tension for each wool, with 95% confidence intervals."----
ci <- predict(wb, newdata = grid, interval = "confidence")
x <- as.integer(grid$tension) + ifelse(grid$wool == "A", -0.08, 0.08)
plot(x, ci$fit, pch = ifelse(grid$wool == "A", 1, 19), xaxt = "n",
     ylim = range(ci$lwr, ci$upr), xlab = "Tension",
     ylab = "Expected breaks")
axis(1, at = 1:3, labels = levels(w$tension))
segments(x, ci$lwr, x, ci$upr)
legend("topright", c("wool A", "wool B"), pch = c(1, 19), bty = "n")


## ----warp-effects-------------------------------------------------------------
wool_eff <- do.call(rbind, lapply(levels(w$tension), function(tn) {
  r <- margEffTable(wb, data = w[w$tension == tn, ], vars = "wool",
                    indicator = "discrete_change", tableType = "tibble")
  data.frame(tension = tn, wool_B_minus_A = r$estimate,
             std_error = r$std_error)
}))
wool_eff
clear <- abs(wool_eff$wool_B_minus_A) > 1.96 * wool_eff$std_error


## ----dv-available, echo = FALSE, results = "asis"-----------------------------
if (!has_aer) {
  cat("*This section uses the `DoctorVisits` data of the AER package, which",
      "is not installed, so its code is not run here.*")
}


## ----dv, eval = has_aer-------------------------------------------------------
e <- new.env()
utils::data("DoctorVisits", package = "AER", envir = e)
dv <- e$DoctorVisits
dvm <- countreg(visits ~ gender + age + I(age^2) + income + illness +
                  health + private, data = dv, family = "NB2",
                method = "BFGS")


## ----dv-ame, eval = has_aer---------------------------------------------------
me_dv <- margEffTable(dvm, vars = c("age", "income"), tableType = "tibble")
me_dv[, c("term", "effect_metric", "estimate", "std_error")]
h <- 1e-4
d_up <- d_dn <- dv
d_up$age <- dv$age + h
d_dn$age <- dv$age - h
mean((predict(dvm, newdata = d_up) - predict(dvm, newdata = d_dn)) / (2 * h))


## ----dv-elasticity, eval = has_aer--------------------------------------------
el <- margEffTable(dvm, vars = "income", measure = "elasticity",
                   tableType = "tibble")
el[, c("term", "effect_metric", "estimate", "std_error", "lower_ci",
       "upper_ci")]


## ----dv-elasticity-text, include = FALSE, eval = has_aer----------------------
el_txt <- paste0(abs(round(el$estimate, 2)), "% ",
                 if (el$estimate < 0) "lower" else "higher")
el_ci <- if (el$lower_ci <= 0 && el$upper_ci >= 0) "includes" else
  "excludes"


## ----dv-curve, eval = has_aer, fig.cap = "Average expected doctor visits in two weeks by age."----
ages <- seq(0.2, 0.7, by = 0.05)
avg_visits <- sapply(ages, function(a) {
  d <- dv
  d$age <- a
  mean(predict(dvm, newdata = d))
})
plot(100 * ages, avg_visits, type = "b", xlab = "Age (years)",
     ylab = "Average expected visits")


## ----bootstrap----------------------------------------------------------------
set.seed(11)
margEffTable(nb, vars = "Eth", indicator = "discrete_change",
             se = "bootstrap", bootstraps = 100,
             tableType = "tibble")[, c("estimate", "std_error",
                                       "successful_bootstraps")]


## ----new-population-----------------------------------------------------------
slow <- q[q$Lrn == "SL", ]
margEffTable(nb, data = slow, vars = c("Eth", "Sex"),
             indicator = "discrete_change",
             tableType = "tibble")[, c("term", "estimate", "std_error")]

