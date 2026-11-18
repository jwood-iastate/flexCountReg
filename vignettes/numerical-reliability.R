## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4)


## ----setup--------------------------------------------------------------------
library(flexCountReg)
data("washington_roads")


## ----status-helper------------------------------------------------------------
fit_status <- function(m) {
  f <- m$model
  ev <- eigen(-f$hessian, symmetric = TRUE, only.values = TRUE)$values
  data.frame(converged = isTRUE(f$converged), code = f$code,
             iterations = unname(f$iterations[1]),
             logLik = round(f$maximum, 3),
             max_abs_gradient = signif(max(abs(f$gradient)), 3),
             min_eigen = signif(min(ev), 3),
             condition = signif(max(ev) / min(ev), 3),
             message = trimws(f$message))
}


## ----status-table-------------------------------------------------------------
f <- Total_crashes ~ lnaadt + lnlength + speed50
fits <- lapply(c(NM = "NM", BFGS = "BFGS", NR = "NR", BHHH = "BHHH"),
               function(meth) countreg(f, data = washington_roads,
                                       family = "NB2", method = meth))
washington_roads$aadt <- exp(washington_roads$lnaadt)
washington_roads$aadt_10k <- washington_roads$aadt / 1e4
fits$NM_raw_aadt <- countreg(Total_crashes ~ aadt + speed50,
                             data = washington_roads, family = "NB2")
fits$NM_aadt_10k <- countreg(Total_crashes ~ aadt_10k + speed50,
                             data = washington_roads, family = "NB2")
status_table <- do.call(rbind, lapply(fits, fit_status))
status_table


## ----status-text, include = FALSE---------------------------------------------
g_raw <- status_table["NM_raw_aadt", "max_abs_gradient"]
c_raw <- status_table["NM_raw_aadt", "condition"]
c_10k <- status_table["NM_aadt_10k", "condition"]


## ----boundary-----------------------------------------------------------------
boundary_fit <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = rnorm(300))
  d$y <- rpois(300, exp(0.5 + 0.4 * d$x))
  warned <- FALSE
  m <- withCallingHandlers(
    countreg(y ~ x, data = d, family = "NB2", method = "BFGS"),
    warning = function(w) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    })
  co <- summary(m)$coefficients
  c(ln_alpha = co["ln(alpha)", 1], se = co["ln(alpha)", 2],
    converged = m$model$converged, warned = warned)
}
boundary <- as.data.frame(t(sapply(1:20, boundary_fit)))
boundary$at_start <- abs(boundary$ln_alpha - log(1e-3)) < 0.01
round(boundary, 3)


## ----boundary-text, include = FALSE-------------------------------------------
n_start <- sum(boundary$at_start)
n_na <- sum(is.na(boundary$se))


## ----rank-deficient-----------------------------------------------------------
d <- data.frame(x = rnorm(100))
d$x2 <- 2 * d$x
d$y <- rpois(100, exp(0.3 + 0.2 * d$x))
res <- tryCatch(countreg(y ~ x + x2, data = d, family = "POISSON"),
                error = function(e) conditionMessage(e))
res


## ----starts-------------------------------------------------------------------
bad_start <- c(0, 0, 0, 0, 3)
from_bad <- countreg(f, data = washington_roads, family = "NB2",
                     method = "BFGS", start.vals = bad_start)
rbind(default_start = fit_status(fits$BFGS)[, c("converged", "logLik")],
      poor_start = fit_status(from_bad)[, c("converged", "logLik")])


## ----pln-engines--------------------------------------------------------------
set.seed(1)
sub <- washington_roads[sample(nrow(washington_roads), 400), ]
pln <- function(engine, ndraws) {
  m <- countreg(f, data = sub, family = "PLN", method = "BFGS",
                engine = engine, ndraws = ndraws)
  data.frame(engine = engine, ndraws = ndraws,
             logLik = m$model$maximum,
             ln_sigma = coef(m)[["ln(sigma)"]],
             se_ln_sigma = summary(m)$coefficients["ln(sigma)", 2])
}
quad <- rbind(pln("quadrature", 50), pln("quadrature", 1500))
halton <- do.call(rbind, lapply(c(25, 50, 100, 200, 400),
                                function(r) pln("halton", r)))
rbind(quad, halton)


## ----pln-plot-----------------------------------------------------------------
plot(halton$ndraws, halton$ln_sigma, type = "b", log = "x", pch = 19,
     ylim = range(halton$ln_sigma, quad$ln_sigma),
     xlab = "Halton draws (log scale)", ylab = "Estimate of ln(sigma)",
     main = "Poisson-lognormal: Halton draws against quadrature")
abline(h = quad$ln_sigma[1], lty = 2)
legend("bottomright", c("Halton", "adaptive quadrature"), lty = c(1, 2),
       pch = c(19, NA), bty = "n")


## ----rp-draws-----------------------------------------------------------------
rp_fit <- function(r) {
  m <- countreg_rp(Total_crashes ~ lnaadt + lnlength,
                   rpar_formula = ~ -1 + speed50, data = sub,
                   family = "NB2", ndraws = r, method = "BFGS")
  co <- summary(m)$coefficients
  data.frame(ndraws = r, converged = m$model$converged,
             logLik = round(m$model$maximum, 3),
             mean = co["speed50:Mean", 1], sd = co["speed50:St.Dev", 1],
             se_sd = co["speed50:St.Dev", 2])
}
rp_draws <- do.call(rbind, lapply(c(50, 100, 200), rp_fit))
rp_draws


## ----lognormal----------------------------------------------------------------
rp_ln <- countreg_rp(Total_crashes ~ lnaadt + lnlength,
                     rpar_formula = ~ -1 + speed50,
                     rpardists = c(speed50 = "ln"), data = sub,
                     family = "NB2", ndraws = 100, method = "BFGS")
nd <- data.frame(lnaadt = log(5000), lnlength = 0, speed50 = c(0, 1))
predict(rp_ln, newdata = nd, method = "Exact")


## ----overflow-----------------------------------------------------------------
rp_n <- countreg_rp(Total_crashes ~ lnaadt + lnlength,
                    rpar_formula = ~ -1 + speed50, data = sub,
                    family = "NB2", ndraws = 100, method = "BFGS")
predict(rp_n, newdata = data.frame(lnaadt = log(5000), lnlength = 0,
                                   speed50 = c(1, 100)),
        method = "Exact")


## ----unresolved---------------------------------------------------------------
nb <- fits$BFGS
at <- data.frame(lnaadt = log(c(2000, 20000)), lnlength = log(2),
                 speed50 = 0)
tolerance_interval(nb, newdata = at, max_count = 15)
tolerance_interval(nb, newdata = at)


## ----versions-----------------------------------------------------------------
c(R = R.version.string,
  flexCountReg = as.character(utils::packageVersion("flexCountReg")),
  maxLik = as.character(utils::packageVersion("maxLik")))

