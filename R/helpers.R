# Helper functions for the flexCountReg package

#' @importFrom stats model.frame model.matrix model.response
#' @importFrom stats dnbinom dpois qnorm qgamma dnorm
#' @importFrom stats na.omit optimize pchisq rgamma rnorm sd
#' @importFrom randtoolbox halton
NULL

# ---- Internal helpers for estimators ----------------------------------------

# Validate and normalize the optimizer name. An unknown method is an error.
.fcr_check_method <- function(method) {
  valid <- c("SANN", "NM", "BFGS", "BFGSR", "CG", "NR", "BHHH")
  m <- toupper(gsub("[^[:alnum:]]", "", method))
  m[m == "NELDERMEAD"] <- "NM"
  if (length(m) != 1L || !m %in% valid) {
    stop("`method` must be one of: ", toString(paste0('"', valid, '"')), ".",
         call. = FALSE)
  }
  m
}

# Validate user starting values: right length, matched by name when named.
.fcr_check_start <- function(start.vals, x_names) {
  start.vals <- unlist(start.vals)
  if (!is.numeric(start.vals)) {
    stop("`start.vals` must be numeric.", call. = FALSE)
  }
  if (any(!is.finite(start.vals))) {
    stop("`start.vals` must be finite (no NA, NaN or Inf).", call. = FALSE)
  }
  if (length(start.vals) != length(x_names)) {
    stop("`start.vals` has ", length(start.vals), " value(s), but the model ",
         "has ", length(x_names), " parameter(s): ", toString(x_names), ".",
         call. = FALSE)
  }
  nm <- names(start.vals)
  if (!is.null(nm) && all(nzchar(nm))) {
    if (!setequal(nm, x_names)) {
      stop("Names of `start.vals` do not match the model parameters. ",
           "Expected: ", toString(x_names), ".", call. = FALSE)
    }
    start.vals <- start.vals[x_names]
  }
  names(start.vals) <- x_names
  start.vals
}

# Probability of an event being reported, given the underreporting linear
# predictor. Both links use the same sign convention: a positive coefficient
# means the event is more likely to be reported.
.report_prob <- function(lin, family = c("logit", "probit")) {
  switch(match.arg(family),
         logit = stats::plogis(lin),
         probit = stats::pnorm(lin))
}

# Did maxLik converge? Return codes differ by optimizer.
.fcr_converged <- function(fit, method) {
  code <- fit$code
  if (is.null(code)) return(NA)
  if (method %in% c("NR", "BHHH", "BFGSR")) {
    code %in% c(1L, 2L, 8L)
  } else {
    code == 0L
  }
}

# Moment-based starting value of the log dispersion from a Poisson fit with
# means mu, by the auxiliary regression of Cameron and Trivedi (2013, Sec.
# 3.4) through the origin: for Var(Y) = mu + a mu^2 (NB2, PIG2, PIG),
# a = sum((y - mu)^2 - y) / sum(mu^2); for Var(Y) = mu (1 + a) (NB1, PIG1),
# a = sum((y - mu)^2 - y) / sum(mu). (A mean of the ratios is unstable when
# many fitted means are small.) Clipped to [1e-3, 100]; NULL for other
# families.
.fcr_moment_disp_start <- function(family, y, mu, w = NULL) {
  if (!family %in% c("NB2", "PIG2", "PIG", "NB1", "PIG1")) return(NULL)
  w <- if (is.null(w)) rep(1, length(y)) else w
  ok <- is.finite(mu) & mu > 0 & w > 0
  if (!any(ok)) return(NULL)
  r <- (y[ok] - mu[ok])^2 - y[ok]
  den <- if (family %in% c("NB1", "PIG1")) mu[ok] else mu[ok]^2
  a <- sum(w[ok] * r) / sum(w[ok] * den)
  if (!is.finite(a)) return(NULL)
  log(min(max(a, 1e-3), 100))
}

# Covariance from the Hessian, with a generalized-inverse fallback. The
# generalized inverse gives zero variance in directions where the likelihood
# is flat (the null space of the Hessian), which would read as perfect
# precision; the variances and covariances of the parameters involved in
# those directions are set to NA instead.
.fcr_hessian_vcov <- function(fit) {
  H <- fit$hessian
  V <- tryCatch(solve(-H), error = function(e) NULL)
  if (is.null(V)) {
    V <- MASS::ginv(-H)
    flat <- .fcr_flat_params(-H)
    if (any(flat)) {
      V[flat, ] <- NA_real_
      V[, flat] <- NA_real_
    }
    nm <- colnames(H) %||% names(fit$estimate)
    warning("The Hessian is singular; standard errors use a generalized ",
            "inverse and may be unreliable. The model may not be identified.",
            if (any(flat)) paste0(
              " The likelihood is flat in a direction involving ",
              toString(if (!is.null(nm)) nm[flat] else which(flat)),
              "; their standard errors are NA."),
            call. = FALSE)
  } else {
    # Invertible, but not a maximum in every direction (for example a
    # parameter drifting to the boundary of its space): negative variances
    bad <- !is.finite(diag(V)) | diag(V) < 0
    if (any(bad)) {
      nm <- colnames(H) %||% names(fit$estimate)
      warning("The Hessian is not negative definite at the estimates, so ",
              "the standard errors of ",
              toString(if (!is.null(nm)) nm[bad] else which(bad)),
              " are NA. The estimate may be on the boundary of the ",
              "parameter space or not at a maximum.", call. = FALSE)
    }
  }
  V
}

# Parameters with a non-negligible loading on the null space of a symmetric
# matrix (eigenvalues below the tolerance used by MASS::ginv)
.fcr_flat_params <- function(M) {
  e <- tryCatch(eigen((M + t(M)) / 2, symmetric = TRUE),
                error = function(err) NULL)
  if (is.null(e) || any(!is.finite(e$values))) return(rep(FALSE, nrow(M)))
  tol <- sqrt(.Machine$double.eps) * max(abs(e$values))
  null <- abs(e$values) <= tol
  if (!any(null)) return(rep(FALSE, nrow(M)))
  rowSums(abs(e$vectors[, null, drop = FALSE])) > 1e-6
}

# Nonparametric bootstrap of a fit.
#
# Without `freq` the n independent units (rows or panels) are resampled with
# replacement and the fit is redone with refit(idx), where a unit drawn
# twice appears twice in `idx`.
#
# With integer frequency weights `freq` (one per unit), a unit with weight w
# stands for w independent units. The bootstrap of those sum(freq) units is
# done on the compressed data: new frequencies are drawn from
# Multinomial(sum(freq), freq / sum(freq)) and the fit is redone with
# refit(idx, w), where `idx` are the units drawn at least once and `w` their
# new frequencies, which replace (not multiply) the original weights.
#
# A replicate is accepted only when the refit returned a result, reported
# convergence (TRUE; FALSE or unknown is a failure), has exactly the expected
# unique coefficient names, and finite estimates. Returns a B x p matrix with
# NA rows for failed replicates and attribute "failures" (counts by reason).
.fcr_bootstrap <- function(n, B, coef_names, refit, freq = NULL) {
  if (!is.null(freq)) .fcr_check_boot_freq(freq, n)
  draws <- matrix(NA_real_, nrow = B, ncol = length(coef_names),
                  dimnames = list(NULL, coef_names))
  why <- character(B)
  for (b in seq_len(B)) {
    est <- tryCatch(suppressWarnings(suppressMessages(
      if (is.null(freq)) {
        refit(sample.int(n, n, replace = TRUE))
      } else {
        counts <- .fcr_boot_frequencies(freq)
        pos <- which(counts > 0)
        refit(pos, counts[pos])
      })), error = function(e) NULL)
    why[b] <- .fcr_boot_check_replicate(est, coef_names)
    if (!nzchar(why[b])) draws[b, ] <- est[coef_names]
  }
  failures <- table(factor(why[nzchar(why)], levels = .fcr_boot_reasons))
  n_ok <- sum(!nzchar(why))
  reasons <- failures[failures > 0]
  if (n_ok < 2) {
    # The estimates are still usable; their standard errors are NA
    warning("Fewer than 2 bootstrap replicates succeeded (",
            paste(names(reasons), reasons, sep = ": ", collapse = ", "),
            "); bootstrapped standard errors are NA.", call. = FALSE)
  } else if (n_ok < 0.9 * B) {
    warning(B - n_ok, " of ", B, " bootstrap replicates failed (",
            paste(names(reasons), reasons, sep = ": ", collapse = ", "),
            "); the bootstrapped standard errors use the ", n_ok,
            " that succeeded.", call. = FALSE)
  }
  attr(draws, "failures") <- failures
  draws
}

.fcr_boot_reasons <- c("error", "not converged", "wrong coefficients",
                       "non-finite estimates")

# "" when a refit result is an acceptable bootstrap replicate, otherwise the
# reason it is not.
.fcr_boot_check_replicate <- function(est, coef_names) {
  if (is.null(est)) return("error")
  if (!isTRUE(attr(est, "converged"))) return("not converged")
  nm <- names(est)
  if (!is.numeric(est) || is.null(nm) || anyNA(nm) || anyDuplicated(nm) ||
      length(est) != length(coef_names) || !setequal(nm, coef_names)) {
    return("wrong coefficients")
  }
  if (!all(is.finite(est))) return("non-finite estimates")
  ""
}

# Bootstrap frequency weights must be nonnegative whole numbers with a
# positive total: other weights (importance or survey weights) do not
# describe a sample of independent units that can be resampled.
.fcr_check_boot_freq <- function(freq, n) {
  if (length(freq) != n || !is.numeric(freq) || anyNA(freq) ||
      any(!is.finite(freq)) || any(freq < 0) ||
      any(abs(freq - round(freq)) > 1e-8) || sum(freq) < 1) {
    stop("Bootstrap standard errors with `weights` need frequency weights: ",
         "nonnegative whole numbers (each row stands for that many ",
         "independent observations) with a positive total. Use ",
         "stderr = \"robust\" for other weights.", call. = FALSE)
  }
  invisible(TRUE)
}

# A column name based on `base` that is not in `existing` (for internal
# helper columns added to user data).
.fcr_unused_name <- function(existing, base) {
  nm <- base
  k <- 0L
  while (nm %in% existing) {
    k <- k + 1L
    nm <- paste0(base, "_", k)
  }
  nm
}

# One bootstrap draw of frequencies for units with frequencies `freq`: the
# number of times each unit is drawn when sum(freq) units are resampled.
.fcr_boot_frequencies <- function(freq) {
  freq <- round(freq)
  N <- sum(freq)
  as.vector(stats::rmultinom(1L, N, freq / N))
}

# Covariance of the accepted bootstrap replicates (complete rows only, so
# every entry uses the same replicates).
.fcr_boot_cov <- function(draws) {
  ok <- stats::complete.cases(draws)
  if (sum(ok) < 2) {
    return(matrix(NA_real_, ncol(draws), ncol(draws),
                  dimnames = list(colnames(draws), colnames(draws))))
  }
  V <- stats::cov(draws[ok, , drop = FALSE])
  dimnames(V) <- list(colnames(draws), colnames(draws))
  V
}

# Estimates of a refit, with its convergence status attached.
.fcr_refit_result <- function(estimate, converged) {
  attr(estimate, "converged") <- length(converged) == 1L &&
    isTRUE(as.logical(converged))
  estimate
}

# --- Internal Likelihood Functions ---

# Random-parameter draws for each observation and draw.
#
#   hdraws             ndraws x K matrix of uniform (Halton) draws
#   random_coefs_means length-K means (log-scale location for "ln")
#   rand_var_params    length-K scales, or the Cholesky elements if correlated
#   rpardists          length-K distribution codes ("n", "ln", "t", "u", "g")
#   X_rand             N x K design matrix of the random parameters
#   het_mean_coefs     heterogeneity in means. For het_form = "additive" a
#                      (n_het x K) matrix: beta_ik = mu_k + z_i' delta_k
#                      ("n", "t", "u"; on the log scale for "ln"), and
#                      beta_ik = mu_k * exp(z_i' delta_k) for "g". For
#                      het_form = "legacy" (fits from flexCountReg 0.1.1 and earlier) a
#                      vector shared by all parameters: mu_ik = mu_k exp(z_i'd).
#   het_var_coefs      heterogeneity in scales: sigma_ik = sigma_k exp(s_i' phi)
#
# Scales: normal and lognormal use the standard deviation; triangular and
# uniform use the half-width, so beta lies in [mu - sigma, mu + sigma]
# (Wood and Gayah, 2025, Appendix A); gamma uses mean mu and SD sigma.
#
# Returns list(xb_rand_mat = N x ndraws matrix of sum_k x_ik beta_ikr).
generate_random_draws <- function(hdraws, random_coefs_means, rand_var_params,
                                  rpardists, rpar, X_rand,
                                  het_mean_coefs = NULL, X_het_mean = NULL,
                                  het_var_coefs = NULL, X_het_var = NULL,
                                  correlated = FALSE,
                                  het_form = c("additive", "legacy")) {
  het_form <- match.arg(het_form)
  n_obs <- nrow(X_rand)
  n_draws <- nrow(hdraws)
  n_params <- length(rpar)
  if (is.null(rpardists)) rpardists <- rep("n", n_params)
  
  loc <- .fcr_rp_locations(random_coefs_means, rpardists, n_obs,
                           het_mean_coefs, X_het_mean, het_var_coefs,
                           X_het_var, het_form)
  mu_mat <- loc$mu_mat
  var_multiplier <- loc$var_multiplier
  
  if (correlated) {
    # --- Correlated (normal only) ---
    L <- .fcr_chol_unpack(rand_var_params, n_params)
    corr_draws <- stats::qnorm(hdraws[, seq_len(n_params), drop = FALSE]) %*%
      t(L)
    xb_deterministic <- rowSums(X_rand * mu_mat)
    xb_stochastic <- (X_rand * var_multiplier) %*% t(corr_draws)
    return(list(xb_rand_mat = xb_stochastic + xb_deterministic))
  }
  
  # --- Uncorrelated ---
  xb_rand_mat <- matrix(0, nrow = n_obs, ncol = n_draws)
  for (k in seq_len(n_params)) {
    mu_i <- mu_mat[, k]
    sd_i <- abs(rand_var_params[k]) * var_multiplier
    u <- hdraws[, k]
    dist_type <- rpardists[k]
    
    if (dist_type == "g") {
      shape_vec <- (mu_i^2) / (sd_i^2)
      rate_vec <- mu_i / (sd_i^2)
      param_draws <- matrix(0, nrow = n_obs, ncol = n_draws)
      for (d in seq_len(n_draws)) {
        param_draws[, d] <- suppressWarnings(
          stats::qgamma(u[d], shape = shape_vec, rate = rate_vec))
      }
    } else {
      std_draws <- switch(dist_type,
                          "n" = stats::qnorm(u),
                          "ln" = stats::qnorm(u),
                          "t" = qtri(u, 0, 1),
                          "u" = 2 * u - 1,
                          stop("Unknown random-parameter distribution: ",
                               dist_type, call. = FALSE))
      param_draws <- outer(mu_i, rep(1, n_draws)) + outer(sd_i, std_draws)
      if (dist_type == "ln") param_draws <- exp(param_draws)
    }
    xb_rand_mat <- xb_rand_mat + X_rand[, k] * param_draws
  }
  list(xb_rand_mat = xb_rand_mat)
}

# Per-observation locations (N x K matrix; log scale for "ln") and scale
# multipliers (length N) of the random parameters. See generate_random_draws().
.fcr_rp_locations <- function(random_coefs_means, rpardists, n_obs,
                              het_mean_coefs = NULL, X_het_mean = NULL,
                              het_var_coefs = NULL, X_het_var = NULL,
                              het_form = c("additive", "legacy")) {
  het_form <- match.arg(het_form)
  n_params <- length(random_coefs_means)
  if (is.null(rpardists)) rpardists <- rep("n", n_params)
  mu_mat <- matrix(random_coefs_means, nrow = n_obs, ncol = n_params,
                   byrow = TRUE)
  log_mult <- matrix(0, n_obs, n_params)
  if (!is.null(X_het_mean) && !is.null(het_mean_coefs)) {
    if (het_form == "legacy") {
      log_mult[] <- as.vector(X_het_mean %*% het_mean_coefs)
    } else {
      shift <- X_het_mean %*% matrix(het_mean_coefs, ncol = n_params)
      add <- rpardists != "g"
      mu_mat[, add] <- mu_mat[, add] + shift[, add]
      log_mult[, !add] <- shift[, !add]
    }
  }
  var_multiplier <- rep(1, n_obs)
  if (!is.null(X_het_var) && !is.null(het_var_coefs)) {
    var_multiplier <- exp(as.vector(X_het_var %*% het_var_coefs))
  }
  list(mu_mat = mu_mat * exp(log_mult), var_multiplier = var_multiplier)
}

# E[exp(t * beta)] for one random parameter (the moment generating function),
# with per-observation t, mu and sigma, using the conventions of
# generate_random_draws(). Wood and Gayah (2025), Appendix A.
#   n: exp(t mu + t^2 sigma^2 / 2)
#   u: exp(t mu) sinh(t sigma) / (t sigma)        (beta on [mu - s, mu + s])
#   t: exp(t (mu - s)) (exp(t s) - 1)^2 / (t s)^2 (symmetric triangular)
#   g: (1 - t sigma^2 / mu)^(-mu^2 / sigma^2), Inf when t sigma^2 / mu >= 1
#   ln: Laplace (saddlepoint) approximation for t < 0; Inf for t > 0, where
#       the lognormal MGF does not exist.
rp_mgf <- function(dist, t, mu, sigma) {
  exp(rp_log_mgf(dist, t, mu, sigma))
}

# Does E[exp(t * beta)] exist (is it finite) for a random coefficient beta?
# Decided from the parameters alone, without evaluating the MGF:
#   normal, uniform, triangular: always, for finite parameters;
#   lognormal: for t <= 0, or when sigma = 0 (a fixed coefficient exp(mu));
#   gamma (mean mu, SD sigma): when mu > 0 and t sigma^2 / mu < 1, or when
#   sigma = 0 (a fixed coefficient mu).
# NA for missing inputs. Overflow of a finite value is a different matter.
.fcr_rp_mean_exists <- function(dist, t, mu, sigma) {
  n <- length(t)
  mu <- rep_len(mu, n)
  sigma <- rep_len(abs(sigma), n)
  fin <- is.finite(t) & is.finite(mu) & is.finite(sigma)
  out <- switch(dist,
    n = , u = , t = fin,
    ln = fin & (sigma == 0 | t <= 0),
    g = fin & mu > 0 & (sigma == 0 | t * sigma^2 / mu < 1),
    stop("Unknown random-parameter distribution: ", dist, call. = FALSE))
  out[is.na(t) | is.na(mu) | is.na(sigma)] <- NA
  out
}

# log E[exp(t * beta)], with Inf where the expected value does not exist.
# Working on the log scale avoids overflow when the values are combined
# across random coefficients. For a lognormal coefficient and t < 0 the
# value is the Laplace approximation of Asmussen, Jensen and
# Rojas-Nandayapa (2016), which is accurate but not exact (the lognormal
# MGF has no closed form); it needs the lamW package.
rp_log_mgf <- function(dist, t, mu, sigma) {
  n <- length(t)
  mu <- rep_len(mu, n)
  sigma <- rep_len(abs(sigma), n)
  ts <- t * sigma
  small <- !is.na(ts) & abs(ts) < 1e-8
  zero <- !is.na(sigma) & sigma == 0
  exists <- .fcr_rp_mean_exists(dist, t, mu, sigma)
  out <- switch(dist,
    n = t * mu + 0.5 * t^2 * sigma^2,
    # Uniform on [mu - sigma, mu + sigma]: exp(t mu) sinh(ts) / ts
    u = t * mu + ifelse(small, 0, .fcr_log_sinhc(ts)),
    # Triangular on [mu - sigma, mu + sigma] with mode mu:
    # exp(t (mu - sigma)) (exp(ts) - 1)^2 / ts^2
    t = ifelse(small, t * mu,
               t * (mu - sigma) + 2 * (.fcr_log_abs_expm1(ts) -
                                         log(abs(ts)))),
    g = {
      b <- 1 - t * sigma^2 / mu
      o <- rep(Inf, n)
      ok <- !is.na(exists) & exists & !zero
      o[ok] <- -(mu[ok]^2) / sigma[ok]^2 * log(b[ok])
      o[zero] <- t[zero] * mu[zero]   # degenerate at mu
      o[!is.na(mu) & mu <= 0 & !zero] <- NaN
      o
    },
    ln = {
      o <- rep(Inf, n)
      o[!is.na(t) & t == 0] <- 0
      o[zero] <- t[zero] * exp(mu[zero])   # degenerate at exp(mu)
      neg <- !is.na(t) & t < 0 & !zero
      if (any(neg)) {
        a <- -exp(mu[neg]) * t[neg] * sigma[neg]^2
        rlang::check_installed(
          "lamW", reason = "for Exact lognormal random-parameter predictions.")
        w <- lamW::lambertW0(a)
        o[neg] <- -(w^2 + 2 * w) / (2 * sigma[neg]^2) - 0.5 * log1p(w)
      }
      o
    },
    stop("Unknown random-parameter distribution: ", dist, call. = FALSE))
  out[is.na(t)] <- NA
  out
}

# log(sinh(x) / x) and log|exp(x) - 1|, stable for large |x|
.fcr_log_sinhc <- function(x) {
  ax <- abs(x)
  ax + log1p(-exp(-2 * ax)) - log(2) - log(ax)
}
.fcr_log_abs_expm1 <- function(x) {
  ifelse(x > 30, x + log1p(-exp(-x)), log(abs(expm1(x))))
}

# Rows of a prepare_model_data() result for a resample of clusters. `idx`
# indexes the cluster levels; clusters drawn more than once become separate
# clusters. `w` (frequency bootstrap) gives the new frequency of each drawn
# cluster; `panel_w` is then the panel weight vector of the subset (NULL
# for unweighted models).
.fcr_cluster_subset <- function(prep, idx, w = NULL) {
  rows <- split(seq_along(prep$y), prep$group)[idx]
  r <- unlist(rows, use.names = FALSE)
  list(y = prep$y[r], X = prep$X[r, , drop = FALSE], offset = prep$offset[r],
       weights = if (!is.null(w)) rep(w, lengths(rows)) else prep$weights[r],
       panel_w = w,
       group = factor(rep(seq_along(rows), lengths(rows))))
}

# Cluster bootstrap: resample whole clusters (with frequency weights `freq`,
# draw new cluster frequencies), refit with refit(subset) -> named
# estimates, and return the covariance of the draws.
.fcr_cluster_bootstrap <- function(prep, B, coef_names, refit, freq = NULL) {
  draws <- .fcr_bootstrap(nlevels(prep$group), B, coef_names,
                          function(idx, w = NULL) {
    refit(.fcr_cluster_subset(prep, idx, w))
  }, freq = freq)
  list(vcov = .fcr_boot_cov(draws), n_ok = sum(stats::complete.cases(draws)),
       failures = attr(draws, "failures"))
}

# Panel frequency weights of a panel model: one weight per level of
# prep$group (the independent units), or NULL without `weights`. A panel
# with weight w counts as w independent panels, so the weight must be
# constant within each panel.
.fcr_panel_weights <- function(prep, weights) {
  if (is.null(weights)) return(NULL)
  w <- prep$weights
  g <- prep$group
  rng <- tapply(w, g, function(v) diff(range(v)))
  if (any(rng > 0, na.rm = TRUE)) {
    stop("`weights` vary within panels. Weights are frequency weights of ",
         "the panels (the independent units), so they must be constant ",
         "within each panel.", call. = FALSE)
  }
  as.numeric(tapply(w, g, `[`, 1)[levels(g)])
}

# Per-panel log-likelihood contributions times the panel frequency weights
# (a weight of 0 removes the panel, even when its contribution is -Inf).
.fcr_weight_ll <- function(ll, w) {
  if (is.null(w)) return(ll)
  out <- ll * w
  out[w == 0] <- 0
  out
}

# Common post-fit bookkeeping for the random-effects estimators.
.fcr_finish_re <- function(fit, prep, x_names, stderr, bootstraps, refit,
                           w_panel = NULL, data = NULL, group_var = NULL) {
  names(fit$estimate) <- x_names
  if (stderr == "boot") {
    if (is.null(bootstraps)) {
      bootstraps <- 200L
      message("`bootstraps` not given; using 200 bootstrap samples.")
    }
    if (!is.numeric(bootstraps) || length(bootstraps) != 1L ||
        bootstraps < 2) {
      stop("`bootstraps` must be a single number of at least 2.",
           call. = FALSE)
    }
    bs <- .fcr_cluster_bootstrap(prep, as.integer(bootstraps), x_names, refit,
                                 freq = w_panel)
    V <- bs$vcov
    fit$boot_failures <- bs$failures
    fit$n_boot_ok <- bs$n_ok
    fit$successful_bootstraps <- bs$n_ok
    fit$vcov_type <- "bootstrap"
    fit$bootstraps <- as.integer(bootstraps)
  } else if (stderr == "robust") {
    # One log-likelihood contribution per panel, so this is cluster-robust
    V <- .fcr_robust_vcov(fit, w_panel)
    fit$vcov_type <- "robust"
  } else {
    V <- .fcr_hessian_vcov(fit)
    fit$vcov_type <- "hessian"
  }
  fit$stderr <- stderr
  dimnames(V) <- list(x_names, x_names)
  fit$vcov <- V
  fit$se <- .fcr_se(V)
  names(fit$se) <- x_names
  fit$bootstrapped_se <- if (stderr == "boot") fit$se else NULL
  fit$blueprint <- prep$blueprint
  fit$nobs <- length(prep$y)
  fit$n_groups <- nlevels(prep$group)
  fit$panel_index <- as.integer(prep$group)
  if (!is.null(data)) fit$lr_sample <- .fcr_lr_sample(prep, data, group_var)
  fit$weights <- prep$weights
  fit$panel_weights <- w_panel
  # Independent units for BIC: the panels, counted with their weights
  fit$n_bic <- if (is.null(w_panel)) nlevels(prep$group) else sum(w_panel)
  fit$rows_used <- prep$rows_used
  fit$n_dropped <- prep$n_dropped
  fit$converged <- .fcr_converged(fit, fit$fcr_method)
  fit
}


# Once-per-session deprecation warning for a renamed argument.
.fcr_deprecated_arg <- function(fun, old, new) {
  rlang::warn(
    paste0("The `", old, "` argument of ", fun, "() is deprecated; use `",
           new, "` instead."),
    class = "lifecycle_warning_deprecated",
    .frequency = "once", .frequency_id = paste0("fcr_", fun, "_", old))
}

# Resolve print.level/verbose to a maxLik printLevel.
.fcr_print_level <- function(fun, verbose, print.level) {
  if (!is.null(print.level)) {
    .fcr_deprecated_arg(fun, "print.level", "verbose")
    return(as.integer(print.level))
  }
  if (isTRUE(verbose)) 2L else 0L
}

# Validate stderr; `bootstraps` given without `stderr` means bootstrap, as in
# earlier versions.
.fcr_resolve_stderr <- function(stderr, stderr_missing, bootstraps) {
  if (!is.null(bootstraps) &&
      (!is.numeric(bootstraps) || length(bootstraps) != 1L ||
       is.na(bootstraps) || bootstraps < 2)) {
    stop("`bootstraps` must be a single number of at least 2.", call. = FALSE)
  }
  if (stderr_missing && !is.null(bootstraps)) return("boot")
  stderr <- match.arg(tolower(stderr), c("normal", "robust", "boot"))
  stderr
}

# Standard errors from a covariance matrix: NA (not 0) where the variance is
# negative or not finite, which happens for weakly identified parameters.
# With k given, a V that is not a k x k matrix gives k NAs.
.fcr_se <- function(V, k = NULL) {
  if (!is.null(k) && !(is.matrix(V) && all(dim(V) == k))) {
    return(rep(NA_real_, k))
  }
  v <- diag(V)
  v[!is.finite(v) | v < 0] <- NA_real_
  sqrt(v)
}

# Zero-inflated or hurdle probabilities from the count family's P(Y = y)
# (probs), P(Y = 0) (f0) and the zero-model probability zp: the structural
# zero probability (inflated) or the probability of a positive count
# (hurdle).
# The same on the log scale: log-probabilities lp and log P(Y = 0) lf0 of the
# count family, and the zero model's linear predictor eta (zp = plogis(eta)).
.fcr_zero_adjust_log <- function(y, lp, lf0, eta, zero_type) {
  l_zp <- stats::plogis(eta, log.p = TRUE)                       # log zp
  l_1mzp <- stats::plogis(eta, lower.tail = FALSE, log.p = TRUE) # log(1-zp)
  if (zero_type == "inflated") {
    ifelse(y == 0, .fcr_logadd(l_zp, l_1mzp + lf0), l_1mzp + lp)
  } else {
    ifelse(y == 0, l_1mzp, l_zp + lp - .fcr_log1mexp(lf0))
  }
}

# Elementwise log(exp(a) + exp(b)) and log(1 - exp(a)) for a <= 0.
.fcr_logadd <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}
.fcr_log1mexp <- function(a) {
  ifelse(a > -log(2), log(-expm1(a)), log1p(-exp(a)))
}

.fcr_zero_adjust <- function(y, probs, f0, zp, zero_type) {
  if (zero_type == "inflated") {
    ifelse(y == 0, zp + (1 - zp) * f0, (1 - zp) * probs)
  } else {
    ifelse(y == 0, 1 - zp, zp * probs / (1 - f0))
  }
}

# Analytic per-observation gradient (scores) of the countreg() log-likelihood
# for the Poisson, NB1 and NB2 families, with an optional dispersion formula
# (no underreporting or zero model). Rows are multiplied by the weights.
# p = c(beta, delta); eta = X beta + offset; mu = exp(eta);
# alpha = exp(W delta). Returns an n x length(p) matrix.
.fcr_analytic_grad <- function(family, p, y, X, offset, W_alpha, w) {
  kx <- ncol(X)
  beta <- p[seq_len(kx)]
  mu <- exp(as.vector(X %*% beta) + offset)
  if (family == "POISSON") {
    G <- (y - mu) * X
  } else {
    delta <- p[kx + seq_len(ncol(W_alpha))]
    alpha <- exp(as.vector(W_alpha %*% delta))
    if (family == "NB2") {
      r <- 1 / alpha
      d_eta <- (y - mu) / (1 + alpha * mu)
      d_r <- digamma(y + r) - digamma(r) + log(r / (r + mu)) +
        (mu - y) / (r + mu)
      d_lnalpha <- -r * d_r                      # dr / dln(alpha) = -r
    } else {                                     # NB1: size r = mu / alpha
      r <- mu / alpha
      psi <- digamma(y + r) - digamma(r)
      d_eta <- r * (psi - log1p(alpha))
      d_lnalpha <- -r * (psi - log1p(alpha)) - r * alpha / (1 + alpha) +
        y / (1 + alpha)
    }
    G <- cbind(d_eta * X, d_lnalpha * W_alpha)
  }
  G <- G * w
  G[w == 0, ] <- 0
  G
}

# ---- Cholesky factor of correlated random parameters -----------------------
# The lower-triangular factor L (Var = L L') is packed column by column, the
# order of L[lower.tri(L, diag = TRUE)]. Labels, starting values, the
# likelihood and predictions all use these helpers.
.fcr_chol_index <- function(K) {
  which(lower.tri(diag(K), diag = TRUE), arr.ind = TRUE)
}

.fcr_chol_names <- function(rpar) {
  idx <- .fcr_chol_index(length(rpar))
  paste0("Chol:", rpar[idx[, 1]], "_", rpar[idx[, 2]])
}

.fcr_chol_pack <- function(L) {
  L[lower.tri(L, diag = TRUE)]
}

.fcr_chol_unpack <- function(v, K) {
  if (length(v) != K * (K + 1) / 2) {
    stop("A Cholesky factor for ", K, " random parameters needs ",
         K * (K + 1) / 2, " values, not ", length(v), ".", call. = FALSE)
  }
  L <- matrix(0, K, K)
  L[lower.tri(L, diag = TRUE)] <- v
  L
}

# Robust (sandwich) covariance for a maxLik fit whose log-likelihood
# contributions (rows or panels) are multiplied by frequency weights w. The
# scores are then w_i s_i, and the meat must be sum_i w_i s_i s_i' (as for
# the expanded data), not sum_i w_i^2 s_i s_i'. With w = 1 this is
# sandwich::sandwich(fit).
.fcr_robust_vcov <- function(fit, w = NULL) {
  G <- sandwich::estfun(fit)
  if (is.null(w)) w <- rep(1, nrow(G))
  keep <- w > 0
  meat <- crossprod(G[keep, , drop = FALSE] / sqrt(w[keep]))
  bread <- .fcr_hessian_vcov(fit)
  V <- bread %*% meat %*% bread
  dimnames(V) <- dimnames(bread)
  V
}

# Scalar optimizer and integration controls: whole numbers, at least 0 for
# max.iters (0 evaluates the model at the starting values) and at least 1
# for ndraws.
.fcr_check_controls <- function(max.iters = NULL, ndraws = NULL) {
  for (nm in c("max.iters", "ndraws")) {
    v <- get(nm)
    if (is.null(v)) next
    lo <- if (nm == "ndraws") 1 else 0
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < lo ||
        v != round(v)) {
      stop("`", nm, "` must be a single whole number of at least ", lo, ".",
           call. = FALSE)
    }
  }
  invisible(NULL)
}

# Halton draws (n x dim). Scrambled draws are generated from `seed`, without
# changing the caller's random-number stream, so the same draws can be
# generated again for prediction.
.fcr_halton <- function(n, dim, scrambled = FALSE, seed = NULL) {
  if (!scrambled || is.null(seed)) {
    return(as.matrix(randtoolbox::halton(n, dim, mixed = scrambled)))
  }
  old <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  as.matrix(randtoolbox::halton(n, dim, mixed = TRUE))
}

