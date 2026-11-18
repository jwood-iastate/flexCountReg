## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4)
old_rng <- RNGkind()

## ----setup--------------------------------------------------------------------
library(flexCountReg)

## ----scenarios----------------------------------------------------------------
scenarios <- expand.grid(family = c("POISSON", "NB2"), n = c(100, 400),
                         stringsAsFactors = FALSE)
scenarios$beta0 <- 0.5
scenarios$beta_x <- 0.4
scenarios$alpha <- ifelse(scenarios$family == "NB2", 0.5, NA)
scenarios$id <- seq_len(nrow(scenarios))
scenarios

## ----generator----------------------------------------------------------------
generate <- function(sc) {
  x <- rnorm(sc$n)
  t <- runif(sc$n, 0.5, 2)
  mu <- t * exp(sc$beta0 + sc$beta_x * x)
  y <- if (sc$family == "POISSON") rpois(sc$n, mu) else
    rnbinom(sc$n, size = 1 / sc$alpha, mu = mu)
  data.frame(y, x, log_t = log(t))
}

## ----check-geometric----------------------------------------------------------
set.seed(1)
draws <- rnbinom(1e5, size = 1, mu = 2)
round(rbind(simulated = tabulate(draws + 1, 6) / 1e5,
            geometric = dgeom(0:5, prob = 1 / 3)), 4)

## ----check-simulate-----------------------------------------------------------
set.seed(2)
d <- generate(scenarios[scenarios$family == "NB2" & scenarios$n == 400, ])
fit <- countreg(y ~ x + offset(log_t), data = d, family = "NB2",
                method = "BFGS")
sims <- simulate(fit, nsim = 2000, seed = 3)
mu <- predict(fit)
alpha <- exp(coef(fit)[["ln(alpha)"]])
summary(data.frame(mean_ratio = rowMeans(sims) / mu,
                   var_ratio = apply(sims, 1, var) / (mu + alpha * mu^2)))

## ----replicate----------------------------------------------------------------
fit_one <- function(d, family) {
  warn <- character(0)
  m <- tryCatch(
    withCallingHandlers(
      countreg(y ~ x + offset(log_t), data = d, family = family,
               method = "BFGS"),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = function(e) e)
  if (inherits(m, "error")) {
    return(data.frame(estimate = NA, se = NA, status = "error",
                      note = conditionMessage(m)))
  }
  co <- summary(m)$coefficients
  status <- if (!isTRUE(m$model$converged)) "not converged" else
    if (length(warn)) "warning" else "ok"
  data.frame(estimate = co["x", 1], se = co["x", 2], status = status,
             note = paste(warn, collapse = "; "))
}

one_replicate <- function(sc, rep) {
  d <- generate(sc)
  out <- do.call(rbind, lapply(c("POISSON", "NB2"), function(proc)
    cbind(scenario = sc$id, rep = rep, procedure = proc, fit_one(d, proc))))
  out$truth <- sc$beta_x
  out
}

## ----streams------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
make_streams <- function(k) {
  s <- .Random.seed
  lapply(seq_len(k), function(i) s <<- parallel::nextRNGStream(s))
}
R <- 100
streams <- make_streams(nrow(scenarios) * R)
run_scenario <- function(sc) {
  do.call(rbind, lapply(seq_len(R), function(r) {
    assign(".Random.seed", streams[[(sc$id - 1) * R + r]],
           envir = globalenv())
    one_replicate(sc, r)
  }))
}
results <- do.call(rbind, lapply(split(scenarios, scenarios$id),
                                 run_scenario))
results <- merge(scenarios[c("id", "family", "n")], results,
                 by.x = "id", by.y = "scenario")
head(results)
table(results$procedure, results$status)

## ----summarize----------------------------------------------------------------
summarize_sim <- function(res) {
  do.call(rbind, lapply(split(res, list(res$id, res$procedure)),
                        function(d) {
    ok <- d$status %in% c("ok", "warning") & is.finite(d$se)
    e <- d$estimate[ok]
    se <- d$se[ok]
    k <- sum(ok)
    cover <- abs(e - d$truth[ok]) <= qnorm(0.975) * se
    data.frame(data = d$family[1], n = d$n[1], procedure = d$procedure[1],
               attempted = nrow(d), used = k,
               bias = mean(e - d$truth[ok]), mcse_bias = sd(e) / sqrt(k),
               emp_sd = sd(e), mean_se = mean(se),
               rmse = sqrt(mean((e - d$truth[ok])^2)),
               coverage = mean(cover),
               mcse_coverage = sqrt(mean(cover) * (1 - mean(cover)) / k))
  }))
}
summary_table <- summarize_sim(results)
summary_table[order(summary_table$data, summary_table$n), -c(4, 5)] |>
  format(digits = 3)
summary_table[, c("data", "n", "procedure", "attempted", "used")]

## ----summary-text, include = FALSE--------------------------------------------
st <- summary_table
pois_on_nb <- st[st$data == "NB2" & st$procedure == "POISSON", ]
nb_on_nb <- st[st$data == "NB2" & st$procedure == "NB2", ]
n_warn <- sum(results$status == "warning")

## ----panel-generator----------------------------------------------------------
generate_panel <- function(n_units = 100, T = 5, theta = 0.5) {
  id <- rep(seq_len(n_units), each = T)
  u <- rgamma(n_units, shape = 1 / theta, rate = 1 / theta)[id]
  z <- rnorm(n_units)[id]
  x <- rnorm(n_units * T)
  data.frame(id, z, x, y = rpois(n_units * T, u * exp(0.2 + 0.4 * z +
                                                        0.3 * x)))
}

## ----panel-study--------------------------------------------------------------
panel_rep <- function(r) {
  assign(".Random.seed", panel_streams[[r]], envir = globalenv())
  d <- generate_panel()
  est <- function(m, nm) {
    co <- summary(m)$coefficients
    c(estimate = co[nm, 1], se = co[nm, 2])
  }
  pooled <- countreg(y ~ z + x, data = d, family = "POISSON",
                     method = "BFGS")
  nm <- negmultinom(y ~ z + x, group_var = "id", data = d)
  rbind(data.frame(rep = r, procedure = "pooled Poisson",
                   term = c("z", "x"),
                   rbind(est(pooled, "z"), est(pooled, "x")),
                   converged = pooled$model$converged),
        data.frame(rep = r, procedure = "negative multinomial",
                   term = c("z", "x"), rbind(est(nm, "z"), est(nm, "x")),
                   converged = nm$model$converged))
}
R_panel <- 100
panel_streams <- make_streams(R_panel)
panel <- do.call(rbind, lapply(seq_len(R_panel), panel_rep))
panel$truth <- ifelse(panel$term == "z", 0.4, 0.3)
panel$cover <- abs(panel$estimate - panel$truth) <=
  qnorm(0.975) * panel$se
aggregate(cbind(estimate, se, cover, converged) ~ procedure + term,
          data = panel, FUN = mean)
aggregate(estimate ~ procedure + term, data = panel, FUN = sd)

## ----script-path--------------------------------------------------------------
script <- system.file("validation", "simulation-studies",
                      "count-regression-study.R", package = "flexCountReg")
file.exists(script)

## ----restore-rng, include = FALSE---------------------------------------------
do.call(RNGkind, as.list(old_rng))

