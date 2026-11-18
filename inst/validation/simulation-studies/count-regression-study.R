# Simulation study of Poisson and NB2 estimators of a regression coefficient,
# the full-size version of vignette("simulation-studies").
#
# Part 1 (independent rows): E[Y] = t exp(0.5 + 0.4 x), x ~ N(0, 1),
# t ~ U(0.5, 2), with Poisson or NB2 (alpha = 0.5) counts and n = 100 or
# 400. Each dataset is fitted with Poisson and NB2 models (BFGS, Hessian
# standard errors); the estimand is the coefficient of x (0.4).
# Part 2 (panel): 100 units observed 5 times, Y | u ~ Poisson(u exp(0.2 +
# 0.4 z + 0.3 x)), u ~ Gamma(mean 1, variance 0.5) shared within a unit, z
# a unit-level and x a row-level covariate. Each dataset is fitted with a
# pooled Poisson model and negmultinom() (NM2); the estimands are the
# coefficients of z (0.4) and x (0.3).
#
# Every replicate has its own L'Ecuyer-CMRG random-number stream derived
# from SEED, so results do not depend on the order of the replicates and
# the loop can be run in parallel with the same streams. Every fit gets a
# status ("ok", "warning", "not converged" or "error"); summaries use the
# fits with a finite standard error and report how many were attempted and
# used, with Monte Carlo standard errors.
#
# Run from the package root (development version) or anywhere with the
# package installed:
#   Rscript inst/validation/simulation-studies/count-regression-study.R
# Environment variables: REPS (default 1000), SEED (default 2026), OUT
# (results directory; default this script's directory in the package
# source).

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(flexCountReg)
}

reps <- as.integer(Sys.getenv("REPS", "1000"))
seed <- as.integer(Sys.getenv("SEED", "2026"))
out_dir <- Sys.getenv("OUT", file.path("inst", "validation",
                                       "simulation-studies"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

RNGkind("L'Ecuyer-CMRG")
set.seed(seed)
make_streams <- function(k) {
  s <- .Random.seed
  lapply(seq_len(k), function(i) s <<- parallel::nextRNGStream(s))
}

# ---- Part 1: independent rows ----------------------------------------------

scenarios <- expand.grid(family = c("POISSON", "NB2"), n = c(100, 400),
                         stringsAsFactors = FALSE)
scenarios$beta0 <- 0.5
scenarios$beta_x <- 0.4
scenarios$alpha <- ifelse(scenarios$family == "NB2", 0.5, NA)
scenarios$id <- seq_len(nrow(scenarios))

generate <- function(sc) {
  x <- stats::rnorm(sc$n)
  t <- stats::runif(sc$n, 0.5, 2)
  mu <- t * exp(sc$beta0 + sc$beta_x * x)
  y <- if (sc$family == "POISSON") stats::rpois(sc$n, mu) else
    stats::rnbinom(sc$n, size = 1 / sc$alpha, mu = mu)
  data.frame(y, x, log_t = log(t))
}

# Fit with a status; `terms` are the coefficients to record
fit_status <- function(expr, terms) {
  warn <- character(0)
  m <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warn <<- c(warn, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  if (inherits(m, "error")) {
    return(data.frame(term = terms, estimate = NA_real_, se = NA_real_,
                      status = "error", note = conditionMessage(m)))
  }
  co <- summary(m)$coefficients
  status <- if (!isTRUE(m$model$converged)) "not converged" else
    if (length(warn)) "warning" else "ok"
  data.frame(term = terms, estimate = co[terms, 1], se = co[terms, 2],
             status = status, note = paste(unique(warn), collapse = "; "),
             row.names = NULL)
}

streams <- make_streams(nrow(scenarios) * reps)
t0 <- proc.time()[3]
rows <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  sc <- scenarios[i, ]
  do.call(rbind, lapply(seq_len(reps), function(r) {
    assign(".Random.seed", streams[[(i - 1) * reps + r]], envir = globalenv())
    d <- generate(sc)
    do.call(rbind, lapply(c("POISSON", "NB2"), function(proc) {
      cbind(part = "rows", family = sc$family, n = sc$n, rep = r,
            procedure = proc,
            fit_status(countreg(y ~ x + offset(log_t), data = d,
                                family = proc, method = "BFGS"), "x"),
            truth = sc$beta_x)
    }))
  }))
}))
time_rows <- proc.time()[3] - t0

# ---- Part 2: panel ----------------------------------------------------------

generate_panel <- function(n_units = 100, T = 5, theta = 0.5) {
  id <- rep(seq_len(n_units), each = T)
  u <- stats::rgamma(n_units, shape = 1 / theta, rate = 1 / theta)[id]
  z <- stats::rnorm(n_units)[id]
  x <- stats::rnorm(n_units * T)
  data.frame(id, z, x,
             y = stats::rpois(n_units * T, u * exp(0.2 + 0.4 * z + 0.3 * x)))
}

panel_streams <- make_streams(reps)
t0 <- proc.time()[3]
panel <- do.call(rbind, lapply(seq_len(reps), function(r) {
  assign(".Random.seed", panel_streams[[r]], envir = globalenv())
  d <- generate_panel()
  rbind(cbind(part = "panel", family = "NM2", n = 100, rep = r,
              procedure = "pooled Poisson",
              fit_status(countreg(y ~ z + x, data = d, family = "POISSON",
                                  method = "BFGS"), c("z", "x")),
              truth = c(0.4, 0.3)),
        cbind(part = "panel", family = "NM2", n = 100, rep = r,
              procedure = "negmultinom",
              fit_status(negmultinom(y ~ z + x, group_var = "id", data = d),
                         c("z", "x")),
              truth = c(0.4, 0.3)))
}))
time_panel <- proc.time()[3] - t0

# ---- Summaries --------------------------------------------------------------

results <- rbind(rows, panel)
summarize <- function(d) {
  ok <- is.finite(d$se) & d$status %in% c("ok", "warning")
  e <- d$estimate[ok]
  k <- sum(ok)
  cover <- abs(e - d$truth[ok]) <= stats::qnorm(0.975) * d$se[ok]
  data.frame(part = d$part[1], data = d$family[1], n = d$n[1],
             procedure = d$procedure[1], term = d$term[1],
             attempted = nrow(d), used = k,
             warnings = sum(d$status == "warning"),
             failures = sum(!ok),
             bias = mean(e - d$truth[ok]), mcse_bias = stats::sd(e) / sqrt(k),
             emp_sd = stats::sd(e), mean_se = mean(d$se[ok]),
             rmse = sqrt(mean((e - d$truth[ok])^2)),
             coverage = mean(cover),
             mcse_coverage = sqrt(mean(cover) * (1 - mean(cover)) / k))
}
groups <- split(results, list(results$part, results$family, results$n,
                              results$procedure, results$term), drop = TRUE)
summary_table <- do.call(rbind, lapply(groups, summarize))
rownames(summary_table) <- NULL

utils::write.csv(results, file.path(out_dir, "study-replicates.csv"),
                 row.names = FALSE)
utils::write.csv(summary_table, file.path(out_dir, "study-summary.csv"),
                 row.names = FALSE)
info <- c(
  paste("flexCountReg", as.character(utils::packageVersion("flexCountReg"))),
  R.version.string,
  paste("RNG:", paste(RNGkind(), collapse = ", ")),
  paste("seed:", seed),
  paste("replications per scenario:", reps),
  paste("date:", format(Sys.Date())),
  sprintf("time: rows %.0f s, panel %.0f s", time_rows, time_panel),
  "statuses:",
  utils::capture.output(print(table(paste(results$part, results$procedure),
                                    results$status))))
writeLines(info, file.path(out_dir, "study-info.txt"))
print(summary_table, digits = 3)
