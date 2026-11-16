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

# Covariance from the Hessian, with a generalized-inverse fallback.
.fcr_hessian_vcov <- function(fit) {
  H <- fit$hessian
  V <- tryCatch(solve(-H), error = function(e) NULL)
  if (is.null(V)) {
    warning("The Hessian is singular; standard errors use a generalized ",
            "inverse and may be unreliable. The model may not be identified.",
            call. = FALSE)
    V <- MASS::ginv(-H)
  }
  V
}

# Nonparametric bootstrap of a set of estimates.
#   n          number of rows (or clusters) to resample
#   B          number of replicates
#   coef_names names of the estimates, in order
#   refit      function(idx) returning a named numeric vector of estimates
# Returns a B x length(coef_names) matrix; failed replicates are NA rows.
.fcr_bootstrap <- function(n, B, coef_names, refit) {
  draws <- matrix(NA_real_, nrow = B, ncol = length(coef_names),
                  dimnames = list(NULL, coef_names))
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    est <- tryCatch(suppressWarnings(suppressMessages(refit(idx))),
                    error = function(e) NULL)
    if (!is.null(est) && all(coef_names %in% names(est))) {
      draws[b, ] <- est[coef_names]
    }
  }
  n_ok <- sum(stats::complete.cases(draws))
  if (n_ok < 2) {
    stop("Fewer than 2 bootstrap replicates succeeded; bootstrapped ",
         "standard errors cannot be computed.", call. = FALSE)
  }
  if (n_ok < 0.9 * B) {
    warning(B - n_ok, " of ", B, " bootstrap replicates failed; the ",
            "bootstrapped standard errors use the ", n_ok, " that succeeded.",
            call. = FALSE)
  }
  draws
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
#                      het_form = "legacy" (fits before flexCountReg 0.1.3) a
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
    L <- matrix(0, n_params, n_params)
    L[lower.tri(L, diag = TRUE)] <- rand_var_params
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
  n <- length(t)
  mu <- rep_len(mu, n)
  sigma <- rep_len(abs(sigma), n)
  ts <- t * sigma
  small <- abs(ts) < 1e-8
  switch(dist,
    n = exp(t * mu + 0.5 * t^2 * sigma^2),
    u = ifelse(small, exp(t * mu), exp(t * mu) * sinh(ts) / ts),
    t = ifelse(small, exp(t * mu),
               exp(t * (mu - sigma)) * (exp(ts) - 1)^2 / ts^2),
    g = {
      b <- 1 - t * sigma^2 / mu
      out <- ifelse(b > 0, b^(-(mu^2) / sigma^2), Inf)
      out[mu <= 0] <- NaN
      out
    },
    ln = {
      out <- rep(Inf, n)
      out[t == 0] <- 1
      neg <- !is.na(t) & t < 0
      if (any(neg)) {
        a <- -exp(mu[neg]) * t[neg] * sigma[neg]^2
        rlang::check_installed(
          "lamW", reason = "for lognormal random-parameter predictions.")
        w <- lamW::lambertW0(a)
        out[neg] <- exp(-(w^2 + 2 * w) / (2 * sigma[neg]^2)) / sqrt(1 + w)
      }
      out[is.na(t)] <- NA
      out
    },
    stop("Unknown random-parameter distribution: ", dist, call. = FALSE))
}

# Rows of a prepare_model_data() result for a resample of clusters. `idx`
# indexes the cluster levels; clusters drawn more than once become separate
# clusters.
.fcr_cluster_subset <- function(prep, idx) {
  rows <- split(seq_along(prep$y), prep$group)[idx]
  r <- unlist(rows, use.names = FALSE)
  list(y = prep$y[r], X = prep$X[r, , drop = FALSE], offset = prep$offset[r],
       weights = prep$weights[r],
       group = factor(rep(seq_along(rows), lengths(rows))))
}

# Cluster bootstrap: resample whole clusters, refit with
# refit(subset) -> named estimates, and return the covariance of the draws.
.fcr_cluster_bootstrap <- function(prep, B, coef_names, refit) {
  draws <- .fcr_bootstrap(nlevels(prep$group), B, coef_names, function(idx) {
    refit(.fcr_cluster_subset(prep, idx))
  })
  V <- stats::cov(draws, use = "pairwise.complete.obs")
  dimnames(V) <- list(coef_names, coef_names)
  list(vcov = V, n_ok = sum(stats::complete.cases(draws)))
}

# Common post-fit bookkeeping for the random-effects estimators.
.fcr_finish_re <- function(fit, prep, x_names, stderr, bootstraps, refit) {
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
    bs <- .fcr_cluster_bootstrap(prep, as.integer(bootstraps), x_names, refit)
    V <- bs$vcov
    fit$n_boot_ok <- bs$n_ok
    fit$successful_bootstraps <- bs$n_ok
    fit$vcov_type <- "bootstrap"
    fit$bootstraps <- as.integer(bootstraps)
  } else if (stderr == "robust") {
    # One log-likelihood contribution per panel, so this is cluster-robust
    V <- sandwich::sandwich(fit)
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
  if (stderr_missing && !is.null(bootstraps)) return("boot")
  stderr <- match.arg(tolower(stderr), c("normal", "robust", "boot"))
  stderr
}

# Standard errors from a covariance matrix: NA (not 0) where the variance is
# negative or not finite, which happens for weakly identified parameters.
.fcr_se <- function(V) {
  v <- diag(V)
  v[!is.finite(v) | v < 0] <- NA_real_
  sqrt(v)
}

# Zero-inflated or hurdle probabilities from the count family's P(Y = y)
# (probs), P(Y = 0) (f0) and the zero-model probability zp: the structural
# zero probability (inflated) or the probability of a positive count
# (hurdle).
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
