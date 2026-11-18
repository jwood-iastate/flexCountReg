# Recovery and Wald-interval coverage of finite-mixture and latent-class
# count models fitted with countreg_mix() / countreg_lc().
#
# Scenarios (true values in `truth` below):
#   pois2_sep    y ~ 1, Poisson means 1 and 6, proportions 0.6/0.4
#                (well separated), n = 300 and 1000.
#   pois2_weak   y ~ 1, Poisson means 2 and 4, proportions 0.5/0.5
#                (weakly separated), n = 300 and 1000.
#   nb2_lc       latent-class NB2 regression: class means
#                exp(0.2 + 0.3 x) and exp(1.5 - 0.4 x), shared alpha = 0.3,
#                membership logit -0.5 + 1.2 z (z binary, only in the
#                membership model), n = 600.
#   panel_lc     panel latent-class Poisson: 120 panels of 3 to 7 rows, one
#                class per panel, class intercepts 0 and 1.2, shared slope
#                0.5 for a row-level x, membership logit -0.3 + 1.0 w with a
#                panel-level w.
#   pois_nb2     mixed families: Poisson component exp(0.2 + 0.6 x) and NB2
#                component exp(1.5 - 0.3 x) with alpha = 0.5, proportions
#                0.3/0.7, n = 800.
#
# Estimands and alignment. Components of the same family are numbered by
# increasing average expected count after every fit (the same convention
# as the true values), so component-specific estimands are aligned by that
# rule; mixed-family components are identified by their family. The
# label-invariant estimand "mean_at_ref" is the mixture mean at a reference
# covariate value (x = 0, z = 0 / w = 0), compared with its true value.
#
# Each replicate has its own L'Ecuyer-CMRG stream derived from SEED, and the
# starts of each fit use a seed drawn from that stream, so results do not
# depend on the order in which replicates are run. Every fit gets a status:
# "ok", "warning" (estimates usable; e.g. a standard error unavailable),
# "not converged" or "error"; summaries use converged fits with a finite
# standard error and report attempted and used counts with Monte Carlo
# standard errors (MCSE).
#
# Run from the package root (development version) or anywhere with the
# package installed:
#   Rscript inst/validation/finite-mixtures/mixture-recovery-study.R
# Environment variables: REPS (default 200), SEED (default 2027), STARTS
# (starts per fit, default 5), OUT (results directory; default this
# script's directory in the package source).

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(flexCountReg)
}

reps <- as.integer(Sys.getenv("REPS", "200"))
seed <- as.integer(Sys.getenv("SEED", "2027"))
n_starts <- as.integer(Sys.getenv("STARTS", "5"))
out_dir <- Sys.getenv("OUT", file.path("inst", "validation",
                                       "finite-mixtures"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

RNGkind("L'Ecuyer-CMRG")
set.seed(seed)
make_streams <- function(k) {
  s <- .Random.seed
  lapply(seq_len(k), function(i) s <<- parallel::nextRNGStream(s))
}

# ---- Scenarios ------------------------------------------------------------------

scenarios <- list(
  list(id = "pois2_sep", n = 300), list(id = "pois2_sep", n = 1000),
  list(id = "pois2_weak", n = 300), list(id = "pois2_weak", n = 1000),
  list(id = "nb2_lc", n = 600), list(id = "panel_lc", n = 120),
  list(id = "pois_nb2", n = 800))

truth <- list(
  pois2_sep = c("Component1:Mean:(Intercept)" = log(1),
                "Component2:Mean:(Intercept)" = log(6),
                "Mixing:Component1:(Intercept)" = log(0.6 / 0.4),
                mean_at_ref = 0.6 * 1 + 0.4 * 6),
  pois2_weak = c("Component1:Mean:(Intercept)" = log(2),
                 "Component2:Mean:(Intercept)" = log(4),
                 "Mixing:Component1:(Intercept)" = 0,
                 mean_at_ref = 3),
  nb2_lc = c("Component1:Mean:(Intercept)" = 0.2,
             "Component1:Mean:x" = 0.3,
             "Component2:Mean:(Intercept)" = 1.5,
             "Component2:Mean:x" = -0.4,
             "Shared:ln(alpha)" = log(0.3),
             "Mixing:Component1:(Intercept)" = -0.5,
             "Mixing:Component1:z" = 1.2,
             mean_at_ref = stats::plogis(-0.5) * exp(0.2) +
               (1 - stats::plogis(-0.5)) * exp(1.5)),
  panel_lc = c("Component1:Mean:(Intercept)" = 0,
               "Component2:Mean:(Intercept)" = 1.2,
               "Shared:Mean:x" = 0.5,
               "Mixing:Component1:(Intercept)" = -0.3,
               "Mixing:Component1:w" = 1.0,
               mean_at_ref = stats::plogis(-0.3) * 1 +
                 (1 - stats::plogis(-0.3)) * exp(1.2)),
  pois_nb2 = c("Component1:Mean:(Intercept)" = 0.2,
               "Component1:Mean:x" = 0.6,
               "Component2:Mean:(Intercept)" = 1.5,
               "Component2:Mean:x" = -0.3,
               "Component2:ln(alpha)" = log(0.5),
               "Mixing:Component1:(Intercept)" = log(0.3 / 0.7),
               mean_at_ref = 0.3 * exp(0.2) + 0.7 * exp(1.5)))

generate <- function(sc) {
  n <- sc$n
  switch(sc$id,
    pois2_sep = {
      cl <- stats::runif(n) < 0.6
      data.frame(y = stats::rpois(n, ifelse(cl, 1, 6)))
    },
    pois2_weak = {
      cl <- stats::runif(n) < 0.5
      data.frame(y = stats::rpois(n, ifelse(cl, 2, 4)))
    },
    nb2_lc = {
      x <- stats::rnorm(n)
      z <- stats::rbinom(n, 1, 0.5)
      cl <- stats::runif(n) < stats::plogis(-0.5 + 1.2 * z)
      mu <- ifelse(cl, exp(0.2 + 0.3 * x), exp(1.5 - 0.4 * x))
      data.frame(y = stats::rnbinom(n, size = 1 / 0.3, mu = mu), x, z)
    },
    panel_lc = {
      len <- sample(3:7, n, replace = TRUE)
      id <- rep(seq_len(n), len)
      w <- stats::rbinom(n, 1, 0.5)
      cl <- (stats::runif(n) < stats::plogis(-0.3 + w))[id]
      x <- stats::rnorm(length(id))
      data.frame(id, w = w[id], x,
                 y = stats::rpois(length(id), exp(ifelse(cl, 0, 1.2) +
                                                    0.5 * x)))
    },
    pois_nb2 = {
      x <- stats::rnorm(n)
      cl <- stats::runif(n) < 0.3
      y <- ifelse(cl, stats::rpois(n, exp(0.2 + 0.6 * x)),
                  stats::rnbinom(n, size = 2, mu = exp(1.5 - 0.3 * x)))
      data.frame(y, x)
    })
}

fit_model <- function(sc, d, s) {
  switch(sc$id,
    pois2_sep = , pois2_weak =
      countreg_mix(y ~ 1, data = d, components = 2, n_starts = n_starts,
                   seed = s),
    nb2_lc = countreg_lc(y ~ x, data = d, family = "NB2", classes = 2,
                         membership_formula = ~ z,
                         shared = list(dis1 = TRUE), n_starts = n_starts,
                         seed = s),
    panel_lc = countreg_lc(y ~ x, data = d, classes = 2,
                           membership_formula = ~ w, group_var = "id",
                           shared = list(mean = "x"), n_starts = n_starts,
                           seed = s),
    pois_nb2 = countreg_mix(y ~ x, data = d,
                            family = c("POISSON", "NB2"),
                            n_starts = n_starts, seed = s))
}

ref_data <- function(sc) {
  switch(sc$id, pois2_sep = , pois2_weak = data.frame(y = 0),
         nb2_lc = data.frame(x = 0, z = 0),
         panel_lc = data.frame(x = 0, w = 0, id = 1),
         pois_nb2 = data.frame(x = 0))
}

# Mean at the reference point and its delta-method standard error
mean_at_ref <- function(m, sc) {
  nd <- ref_data(sc)
  f <- function(th) {
    mm <- m
    mm$model$estimate[] <- th
    stats::predict(mm, newdata = nd)
  }
  th <- stats::coef(m)
  g <- vapply(seq_along(th), function(j) {
    h <- 1e-6 * max(1, abs(th[j]))
    a <- th
    b <- th
    a[j] <- a[j] + h
    b[j] <- b[j] - h
    (f(a) - f(b)) / (2 * h)
  }, numeric(1))
  c(estimate = f(th), se = sqrt(drop(t(g) %*% stats::vcov(m) %*% g)))
}

one_fit <- function(sc, r, stream) {
  assign(".Random.seed", stream, envir = globalenv())
  d <- generate(sc)
  s <- sample.int(.Machine$integer.max, 1L)
  tr <- truth[[sc$id]]
  terms <- names(tr)
  warn <- character(0)
  t0 <- proc.time()[3]
  m <- tryCatch(withCallingHandlers(suppressMessages(fit_model(sc, d, s)),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  secs <- proc.time()[3] - t0
  base <- data.frame(scenario = sc$id, n = sc$n, rep = r, term = terms,
                     truth = unname(tr), secs = secs)
  if (inherits(m, "error")) {
    return(cbind(base, estimate = NA_real_, se = NA_real_, status = "error",
                 reached_best = NA_integer_, note = conditionMessage(m)))
  }
  cf <- stats::coef(m)
  se <- m$model$se
  ref <- mean_at_ref(m, sc)
  est <- c(cf, mean_at_ref = ref[["estimate"]])[terms]
  ses <- c(se, mean_at_ref = ref[["se"]])[terms]
  status <- if (!isTRUE(m$model$converged)) "not converged" else
    if (length(warn)) "warning" else "ok"
  cbind(base, estimate = unname(est), se = unname(ses), status = status,
        reached_best = sum(m$model$starts$reached_best),
        note = paste(unique(warn), collapse = "; "))
}

streams <- make_streams(length(scenarios) * reps)
t_start <- proc.time()[3]
results <- do.call(rbind, lapply(seq_along(scenarios), function(i) {
  sc <- scenarios[[i]]
  message("Scenario ", sc$id, " (n = ", sc$n, ")")
  do.call(rbind, lapply(seq_len(reps), function(r) {
    one_fit(sc, r, streams[[(i - 1) * reps + r]])
  }))
}))
elapsed <- proc.time()[3] - t_start

# ---- Summaries ----------------------------------------------------------------

summarize <- function(d) {
  ok <- is.finite(d$se) & d$status %in% c("ok", "warning")
  e <- d$estimate[ok]
  k <- sum(ok)
  cover <- abs(e - d$truth[ok]) <= stats::qnorm(0.975) * d$se[ok]
  data.frame(scenario = d$scenario[1], n = d$n[1], term = d$term[1],
             truth = d$truth[1], attempted = nrow(d),
             converged = sum(d$status %in% c("ok", "warning")), used = k,
             bias = mean(e - d$truth[ok]),
             mcse_bias = stats::sd(e) / sqrt(k),
             emp_sd = stats::sd(e), mean_se = mean(d$se[ok]),
             rmse = sqrt(mean((e - d$truth[ok])^2)),
             coverage = mean(cover),
             mcse_coverage = sqrt(mean(cover) * (1 - mean(cover)) / k),
             single_start_solutions = sum(d$reached_best[ok] < 2,
                                          na.rm = TRUE),
             mean_secs = mean(d$secs))
}
groups <- split(results, list(results$scenario, results$n, results$term),
                drop = TRUE)
summary_table <- do.call(rbind, lapply(groups, summarize))
summary_table <- summary_table[order(summary_table$scenario,
                                     summary_table$n), ]
rownames(summary_table) <- NULL

utils::write.csv(results, file.path(out_dir, "mixture-replicates.csv"),
                 row.names = FALSE)
utils::write.csv(summary_table, file.path(out_dir, "mixture-summary.csv"),
                 row.names = FALSE)
info <- c(
  paste("flexCountReg", as.character(utils::packageVersion("flexCountReg"))),
  R.version.string,
  paste("RNG:", paste(RNGkind(), collapse = ", ")),
  paste("seed:", seed),
  paste("replications per scenario:", reps),
  paste("starts per fit:", n_starts),
  paste("date:", format(Sys.Date())),
  sprintf("time: %.0f s", elapsed),
  "statuses (one row per fit and term):",
  utils::capture.output(print(table(paste(results$scenario, results$n),
                                    results$status))))
writeLines(info, file.path(out_dir, "mixture-info.txt"))
print(summary_table, digits = 3)
