## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")


## ----setup--------------------------------------------------------------------
library(flexCountReg)


## ----simulate-----------------------------------------------------------------
set.seed(42)
n <- 400
d <- data.frame(id = seq_len(n), x = rnorm(n),
                exposure = round(runif(n, 0.2, 5), 2))
d$y <- rpois(n, d$exposure * exp(-0.5 + 0.6 * d$x))
d$log_exposure <- log(d$exposure)
head(d)


## ----offsets------------------------------------------------------------------
# 1. offset() in the formula, with the exposure on its original scale
m_formula <- countreg(y ~ x + offset(log(exposure)), data = d,
                      family = "POISSON", method = "BFGS")
# 2. the `offset` argument, naming a column that is already logged
m_column <- countreg(y ~ x, data = d, family = "POISSON",
                     offset = "log_exposure", method = "BFGS")
rbind(formula = coef(m_formula), column = coef(m_column))
all.equal(coef(m_formula), coef(m_column))


## ----offset-glm---------------------------------------------------------------
coef(glm(y ~ x + offset(log(exposure)), data = d, family = poisson()))


## ----both-offsets-------------------------------------------------------------
tryCatch(countreg(y ~ x + offset(log(exposure)), data = d,
                  family = "POISSON", offset = "log_exposure"),
         error = function(e) conditionMessage(e))


## ----estimate-exposure--------------------------------------------------------
m_free <- countreg(y ~ x + log_exposure, data = d, family = "POISSON",
                   method = "BFGS")
summary(m_free)$coefficients[, 1:2]


## ----exposure-check, include = FALSE------------------------------------------
b_exp <- coef(m_free)[["log_exposure"]]
se_exp <- summary(m_free)$coefficients["log_exposure", 2]
near_one <- abs(b_exp - 1) < 2 * se_exp


## ----rate-response------------------------------------------------------------
tryCatch(countreg(I(y / exposure) ~ x, data = d, family = "POISSON",
                  method = "BFGS"),
         error = function(e) conditionMessage(e))


## ----log-twice----------------------------------------------------------------
tryCatch(suppressWarnings(
  countreg(y ~ x + offset(log(log_exposure)), data = d,
           family = "POISSON", method = "BFGS")),
  error = function(e) conditionMessage(e))


## ----bad-exposure-------------------------------------------------------------
d_bad <- d
d_bad$exposure[3] <- 0
tryCatch(countreg(y ~ x + offset(log(exposure)), data = d_bad,
                  family = "POISSON", method = "BFGS"),
         error = function(e) conditionMessage(e))


## ----predict-exposure---------------------------------------------------------
predict(m_formula, newdata = data.frame(x = 0, exposure = c(1, 2, 10)))


## ----weights-data-------------------------------------------------------------
compact <- data.frame(x = rep(c(0, 1), each = 4), y = rep(0:3, 2),
                      freq = c(5, 8, 4, 1, 2, 6, 6, 3))
compact
expanded <- compact[rep(seq_len(nrow(compact)), compact$freq),
                    c("x", "y")]
nrow(expanded)


## ----weights-fit--------------------------------------------------------------
f_compact <- countreg(y ~ x, data = compact, family = "POISSON",
                      weights = "freq", method = "BFGS")
f_expanded <- countreg(y ~ x, data = expanded, family = "POISSON",
                       method = "BFGS")
rbind(compact = c(coef(f_compact), logLik = logLik(f_compact),
                  BIC = BIC(f_compact)),
      expanded = c(coef(f_expanded), logLik = logLik(f_expanded),
                   BIC = BIC(f_expanded)))
rbind(compact = sqrt(diag(vcov(f_compact))),
      expanded = sqrt(diag(vcov(f_expanded))))


## ----fractional-weights-------------------------------------------------------
tryCatch(countreg(y ~ x, data = transform(compact, freq = freq / 2),
                  family = "POISSON", weights = "freq", method = "BFGS",
                  stderr = "boot", bootstraps = 20),
         error = function(e) conditionMessage(e))


## ----missing-data-------------------------------------------------------------
dm <- d
dm$y_nb <- rnbinom(n, size = 2, mu = dm$exposure * exp(-0.5 + 0.6 * dm$x))
dm$z <- rnorm(n)
dm$x[c(4, 9)] <- NA
dm$z[c(9, 15, 22)] <- NA
m_mean <- countreg(y_nb ~ x + offset(log(exposure)), data = dm, family = "NB2",
                   method = "BFGS")
m_disp <- countreg(y_nb ~ x + offset(log(exposure)), data = dm, family = "NB2",
                   dis_param_formula_1 = ~ z, method = "BFGS")
c(mean_model = nobs(m_mean), dispersion_model = nobs(m_disp))


## ----dropped-rows-------------------------------------------------------------
dm$id[setdiff(seq_len(nrow(dm)), m_mean$model$rows_used)]
dm$id[setdiff(seq_len(nrow(dm)), m_disp$model$rows_used)]


## ----compare-different--------------------------------------------------------
tryCatch(regCompTest(m_disp, model2 = m_mean),
         error = function(e) conditionMessage(e))


## ----compare-common-----------------------------------------------------------
vars <- c("y_nb", "x", "z", "exposure")
common <- dm[complete.cases(dm[vars]), ]
m_mean2 <- countreg(y_nb ~ x + offset(log(exposure)), data = common,
                    family = "NB2", method = "BFGS")
m_disp2 <- countreg(y_nb ~ x + offset(log(exposure)), data = common,
                    family = "NB2", dis_param_formula_1 = ~ z,
                    method = "BFGS")
c(nobs(m_mean2), nobs(m_disp2))
regCompTest(m_disp2, model2 = m_mean2)$statistics

