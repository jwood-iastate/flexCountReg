# Coverage study for the interval methods of flexCountReg (countreg models).
#
# For each setting (family, sample size), estimation samples are simulated
# from known parameters with a fixed design, the model is fit, and for each
# scenario x0 the following are computed from that sample alone:
#   * tolerance_interval(): joint and envelope methods (normal uncertainty),
#     content 0.90, confidence 0.95 -> record whether the TRUE content
#     F(U | x0) - F(L - 1 | x0) is >= 0.90 (target: fraction >= 0.95);
#   * predict(interval = "prediction", level = 0.90) -> record the TRUE
#     probability that a new count falls in the interval (target: average
#     >= 0.90);
#   * predict(interval = "confidence", level = 0.95) -> record whether the
#     interval contains the true mean (target: 0.95).
# Monte Carlo standard errors are reported. Run from the package root with
# the development version installed or loaded:
#   Rscript inst/validation/intervals/coverage-study.R
# Environment variables: REPS (default 500), OUT (results directory),
# SETTINGS (optional subset, e.g. "NB2:50,POISSON:200").
# Per-replicate results are written to coverage-replicates.csv; intervals
# that are unavailable (status other than "ok") are counted separately and
# excluded from the coverage fractions, which are therefore conditional on
# availability.

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(flexCountReg)
}

reps <- as.integer(Sys.getenv("REPS", "500"))
out_dir <- Sys.getenv("OUT", file.path("inst", "validation", "intervals"))
content <- 0.90
confidence <- 0.95
beta <- c(0.5, 0.5)
alpha <- 0.5
x0s <- c(0, 2)

true_cdf <- function(fam, q, mu) {
  if (fam == "POISSON") stats::ppois(q, mu) else
    stats::pnbinom(q, size = 1 / alpha, mu = mu)
}

one_rep <- function(fam, n, x) {
  mu <- exp(beta[1] + beta[2] * x)
  y <- if (fam == "POISSON") stats::rpois(n, mu) else
    stats::rnbinom(n, size = 1 / alpha, mu = mu)
  fit <- tryCatch(suppressWarnings(countreg(
    y ~ x, data = data.frame(y = y, x = x), family = fam, method = "BFGS")),
    error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$model$converged)) {
    return(data.frame(status = "fit failed"))
  }
  scen <- data.frame(x = x0s)
  mu0 <- exp(beta[1] + beta[2] * x0s)
  res <- tryCatch(suppressMessages({
    tj <- tolerance_interval(fit, scen, content, confidence, method = "joint")
    te <- tolerance_interval(fit, scen, content, confidence,
                             method = "envelope")
    pr <- predict(fit, scen, interval = "prediction", level = 0.90)
    ci <- predict(fit, scen, interval = "confidence", level = 0.95)
    cont <- function(L, U) true_cdf(fam, U, mu0) - true_cdf(fam, L - 1, mu0)
    data.frame(
      x0 = x0s, status = "ok",
      uncertainty = attr(tj, "tolerance")$uncertainty,
      ti_status = tj$status,
      ti_joint_ok = cont(tj$lower, tj$upper) >= content,
      ti_joint_width = tj$upper - tj$lower,
      ti_env_ok = cont(te$lower, te$upper) >= content,
      ti_env_width = te$upper - te$lower,
      pi_cover = cont(pr$lwr, pr$upr),
      pi_width = pr$upr - pr$lwr,
      ci_cover = ci$lwr <= mu0 & mu0 <= ci$upr)
  }), error = function(e) data.frame(status = conditionMessage(e)))
  res
}

set.seed(20261117)
settings <- expand.grid(fam = c("POISSON", "NB2"), n = c(50, 200),
                        stringsAsFactors = FALSE)
sel <- Sys.getenv("SETTINGS", "")
if (nzchar(sel)) {
  keep <- paste(settings$fam, settings$n, sep = ":") %in%
    strsplit(sel, ",")[[1]]
  settings <- settings[keep, , drop = FALSE]
}
all <- list()
for (i in seq_len(nrow(settings))) {
  s <- settings[i, ]
  x <- stats::rnorm(s$n)   # fixed design within a setting
  t0 <- proc.time()[3]
  r <- lapply(seq_len(reps), function(b) {
    out <- one_rep(s$fam, s$n, x)
    out$rep <- b
    out
  })
  r <- do.call(rbind, lapply(r, function(d) {
    cols <- c("x0", "status", "uncertainty", "ti_status", "ti_joint_ok",
              "ti_joint_width",
              "ti_env_ok", "ti_env_width", "pi_cover", "pi_width",
              "ci_cover", "rep")
    d[setdiff(cols, names(d))] <- NA
    d[cols]
  }))
  r$fam <- s$fam
  r$n <- s$n
  all[[i]] <- r
  message(sprintf("%s n=%d: %.0f s", s$fam, s$n, proc.time()[3] - t0))
}
res <- do.call(rbind, all)
ok <- res[res$status == "ok", ]
mcse <- function(p, m) sqrt(p * (1 - p) / m)
summ <- do.call(rbind, lapply(split(ok, list(ok$fam, ok$n, ok$x0),
                                    drop = TRUE), function(d) {
  m <- nrow(d)
  av <- d$ti_status == "ok"
  ma <- sum(av)
  data.frame(
    family = d$fam[1], n = d$n[1], x0 = d$x0[1], samples = m,
    bootstrap_used = sum(d$uncertainty == "bootstrap"),
    ti_unavailable = m - ma,
    ti_joint_conf = mean(d$ti_joint_ok[av]),
    ti_joint_se = mcse(mean(d$ti_joint_ok[av]), ma),
    ti_joint_width = mean(d$ti_joint_width[av]),
    ti_env_conf = mean(d$ti_env_ok[av], na.rm = TRUE),
    ti_env_width = mean(d$ti_env_width[av], na.rm = TRUE),
    pi_coverage = mean(d$pi_cover), pi_width = mean(d$pi_width),
    ci_coverage = mean(d$ci_cover),
    ci_se = mcse(mean(d$ci_cover), m))
}))
summ <- summ[order(summ$family, summ$n, summ$x0), ]
failures <- table(res$status[res$status != "ok"])
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tag <- if (nzchar(sel)) paste0("-", gsub("[:,]", "_", sel)) else ""
utils::write.csv(summ, file.path(out_dir, paste0("coverage-summary", tag,
                                                 ".csv")),
                 row.names = FALSE)
utils::write.csv(res, file.path(out_dir, paste0("coverage-replicates", tag,
                                                ".csv")),
                 row.names = FALSE)
ti_bad <- table(ok$ti_status[ok$ti_status != "ok"])
writeLines(c(
  paste("flexCountReg", as.character(utils::packageVersion("flexCountReg"))),
  paste(R.version.string), paste("replications per setting:", reps),
  paste("date:", format(Sys.Date())),
  "fit failures:", if (length(failures)) paste(names(failures), failures)
  else "none",
  "unavailable tolerance intervals:",
  if (length(ti_bad)) paste(names(ti_bad), ti_bad) else "none"),
  file.path(out_dir, paste0("coverage-info", tag, ".txt")))
print(summ, row.names = FALSE, digits = 3)
print(failures)
print(ti_bad)
