# Log-likelihood of countreg_rp() and its analytic gradient.
#
# The model: for panel g with rows t and simulation draws r = 1..R,
#   L_g = (1 / R) sum_r prod_{t in g} f(y_t | mu_tr, alpha_t, sigma_t),
#   log(mu_tr) = eta_tr = x_t' beta + offset_t + sum_k x_tk beta_tkr,
# with beta_tkr the r-th draw of the k-th random parameter. With the
# posterior draw weights w_gr = L_gr / sum_r L_gr (the weights the
# Empirical Bayes predictor uses), the panel score is
#   d log L_g / d theta = sum_r w_gr sum_t [ s_tr d eta_tr / d theta
#                                            + d log f_tr / d theta_dist ],
# where s_tr = d log f / d eta is the family's score in the log mean
# (.fcr_family_score()). The derivative of eta_tr is simple for every
# parameter type:
#   fixed beta_j                      x_tj
#   mean of an n, u or t parameter    x_tk
#   lognormal location mu_k           x_tk beta_tkr
#   scale sigma_k                     x_tk z_kr (times beta_tkr, lognormal)
#   Cholesky element L_kl             x_tk z_lr
#   heterogeneity in means delta_k    x_tk z_t (times beta_tkr, lognormal)
#   heterogeneity in scales phi       sum_k x_tk sigma_tk z_kr s_t
# Gamma random parameters (quantile draws) use the derivative of the gamma
# quantile in its shape, by central differences. So one likelihood pass and
# one score pass give the gradient, instead of two likelihood passes per
# parameter for a numeric gradient.
#
# `d` is the fixed part of the model (data, draws and sizes), built once by
# countreg_rp().

# Split the parameter vector into its parts.
.fcr_rp_unpack <- function(p, d) {
  idx <- 0
  take <- function(k) {
    if (k == 0) return(NULL)
    out <- p[idx + seq_len(k)]
    idx <<- idx + k
    out
  }
  out <- list(fixed = take(d$N_fixed), rand_means = take(d$N_rand),
              rand_vars = take(d$n_var_params),
              het_mean = take(d$N_het_mean), het_var = take(d$N_het_var),
              dis_1 = take(d$N_dis_1), dis_2 = take(d$N_dis_2))
  if (!is.null(out$het_mean)) {
    out$het_mean <- matrix(out$het_mean, ncol = d$N_rand)
  }
  out
}

# Evaluate the log-likelihood at p. Returns the per-panel contributions
# (`ll`, multiplied by the panel weights) and the pieces the score needs.
.fcr_rp_eval <- function(p, d) {
  pp <- .fcr_rp_unpack(p, d)
  N <- d$N_obs
  vec_1 <- if (d$N_dis_1 == 0) NULL else if (is.null(d$X_dis_1))
    rep(exp(pp$dis_1), N) else exp(as.vector(d$X_dis_1 %*% pp$dis_1))
  vec_2 <- if (d$N_dis_2 == 0) NULL else if (is.null(d$X_dis_2))
    rep(exp(pp$dis_2), N) else exp(as.vector(d$X_dis_2 %*% pp$dis_2))
  G <- length(d$w_panel)

  # GP2 with negative dispersion a requires 1 + a * mu > 0, which
  # unbounded random parameters violate for some draws of every
  # observation: the model is only defined for a >= 0 (alpha = a + 1)
  if (d$family == "GP2" && any(vec_1 < 1)) {
    return(list(p = p, ll = rep(-Inf, G), valid = FALSE))
  }
  mu_fixed <- exp(as.vector(d$X_Fixed %*% pp$fixed) + d$X_offset)

  draws_info <- generate_random_draws(hdraws = d$hdraws,
                                      random_coefs_means = pp$rand_means,
                                      rand_var_params = pp$rand_vars,
                                      rpardists = d$rpardists,
                                      rpar = d$rpar,
                                      X_rand = d$X_rand,
                                      het_mean_coefs = pp$het_mean,
                                      X_het_mean = d$X_het_mean,
                                      het_var_coefs = pp$het_var,
                                      X_het_var = d$X_het_var,
                                      correlated = d$correlated)
  mu_mat <- exp(draws_info$xb_rand_mat) * mu_fixed
  D <- ncol(mu_mat)

  probs_flat <- tryCatch(
    d$probFunc(y = rep(d$y, times = D),
               predicted = as.vector(mu_mat),
               alpha = if (!is.null(vec_1)) rep(vec_1, times = D),
               sigma = if (!is.null(vec_2)) rep(vec_2, times = D),
               haltons = d$dist_haltons,
               normed_haltons = d$normed_dist_haltons,
               engine = d$engine, log = TRUE),
    error = function(e) rep(NA_real_, N * D))

  # Panel log-likelihood per draw, averaged over draws on the log scale
  log_probs <- matrix(probs_flat, nrow = N, ncol = D)
  sum_log <- rowsum(log_probs, group = d$panel_group, reorder = TRUE)
  m <- sum_log[cbind(seq_len(nrow(sum_log)),
                     max.col(sum_log, ties.method = "first"))]
  ll_panel <- m + log(rowMeans(exp(sum_log - m)))
  ll_panel[!is.finite(m)] <- m[!is.finite(m)]
  # Panels with zero weight contribute exactly 0 (also when -Inf)
  ll_panel[d$w_panel == 0] <- 0

  list(p = p, ll = ll_panel * d$w_panel, valid = TRUE, pars = pp,
       vec_1 = vec_1, vec_2 = vec_2, mu_mat = mu_mat, log_probs = log_probs,
       sum_log = sum_log, m = m)
}

# Per-panel scores (G x P, multiplied by the panel weights) from an
# evaluated state `st` of .fcr_rp_eval().
.fcr_rp_score <- function(st, d) {
  P <- length(st$p)
  G <- length(d$w_panel)
  if (!isTRUE(st$valid)) return(matrix(NA_real_, G, P))
  N <- d$N_obs
  K <- d$N_rand
  R <- ncol(st$mu_mat)
  pp <- st$pars

  # Posterior draw weights of each panel, expanded to its rows
  Wg <- exp(st$sum_log - st$m)
  Wg <- Wg / rowSums(Wg)
  Wobs <- Wg[d$panel_idx, , drop = FALSE]
  weigh <- function(v) {
    A <- Wobs * matrix(v, N, R)
    A[!is.na(Wobs) & Wobs == 0] <- 0
    A
  }

  sc <- .fcr_family_score(
    d$family, rep(d$y, times = R), as.vector(st$mu_mat),
    alpha = if (!is.null(st$vec_1)) rep(st$vec_1, times = R),
    sigma = if (!is.null(st$vec_2)) rep(st$vec_2, times = R),
    center = as.vector(st$log_probs),
    dens_args = list(haltons = d$dist_haltons,
                     normed_haltons = d$normed_dist_haltons,
                     engine = d$engine))
  A <- weigh(sc$eta)
  rA <- rowSums(A)

  # Random parameters: locations, scale multipliers and scales
  loc <- .fcr_rp_locations(pp$rand_means, d$rpardists, N, pp$het_mean,
                           d$X_het_mean, pp$het_var, d$X_het_var)
  v <- loc$var_multiplier
  X_rand <- d$X_rand
  # M[, k] = sum_r A d beta_tkr / d loc_tk; dloc[, k] = d loc_tk / d mu_k
  M <- matrix(rA, N, K)
  dloc <- matrix(1, N, K)
  is_g <- d$rpardists == "g"
  if (any(is_g) && !is.null(pp$het_mean)) {
    shift <- d$X_het_mean %*% pp$het_mean
    dloc[, is_g] <- exp(shift[, is_g])
  }

  if (d$correlated) {
    Z <- stats::qnorm(d$hdraws[, seq_len(K), drop = FALSE])
    L <- .fcr_chol_unpack(pp$rand_vars, K)
    AZ <- A %*% Z
    ci <- .fcr_chol_index(K)
    var_cols <- X_rand[, ci[, 1], drop = FALSE] * v *
      AZ[, ci[, 2], drop = FALSE]
    # d eta / d phi = s_t times the stochastic part of eta
    hv <- rowSums((X_rand * v) * (AZ %*% t(L)))
  } else {
    C <- matrix(0, N, K)        # sum_r A d beta_tkr / d sd_tk
    s_k <- pp$rand_vars
    for (k in seq_len(K)) {
      u <- d$hdraws[, k]
      sd_i <- abs(s_k[k]) * v
      dist <- d$rpardists[k]
      if (dist == "g") {
        g <- .fcr_rp_gamma_derivs(u, loc$mu_mat[, k], sd_i)
        M[, k] <- rowSums(A * g$dmu)
        C[, k] <- rowSums(A * g$dsd)
      } else {
        z <- switch(dist,
                    "n" = stats::qnorm(u),
                    "ln" = stats::qnorm(u),
                    "t" = qtri(u, 0, 1),
                    "u" = 2 * u - 1)
        if (dist == "ln") {
          AB <- A * exp(outer(loc$mu_mat[, k], rep(1, R)) + outer(sd_i, z))
          M[, k] <- rowSums(AB)
          C[, k] <- as.vector(AB %*% z)
        } else {
          C[, k] <- as.vector(A %*% z)
        }
      }
    }
    var_cols <- sweep(X_rand * C * v, 2, sign(s_k), `*`)
    hv <- rowSums(X_rand * C * outer(v, abs(s_k)))
  }

  cols <- list(d$X_Fixed * rA, X_rand * M * dloc, var_cols)
  if (d$N_het_mean > 0) {
    # d loc_tk / d delta_k = z_t (times loc_tk for a gamma parameter)
    gm <- matrix(1, N, K)
    gm[, is_g] <- loc$mu_mat[, is_g]
    f <- X_rand * M * gm
    cols <- c(cols, list(do.call(cbind, lapply(seq_len(K), function(k)
      d$X_het_mean * f[, k]))))
  }
  if (d$N_het_var > 0) cols <- c(cols, list(d$X_het_var * hv))
  if (d$N_dis_1 > 0) {
    h1 <- rowSums(weigh(sc$alpha))
    cols <- c(cols, list(if (is.null(d$X_dis_1)) h1 else d$X_dis_1 * h1))
  }
  if (d$N_dis_2 > 0) {
    h2 <- rowSums(weigh(sc$sigma))
    cols <- c(cols, list(if (is.null(d$X_dis_2)) h2 else d$X_dis_2 * h2))
  }
  Grow <- do.call(cbind, cols)
  Gp <- rowsum(Grow, d$panel_group, reorder = TRUE) * d$w_panel
  Gp[d$w_panel == 0, ] <- 0
  Gp[!is.finite(st$m) & d$w_panel != 0, ] <- NA_real_
  dimnames(Gp) <- NULL
  Gp
}

# Derivatives of gamma random-parameter draws beta = qgamma(u, a, b), with
# shape a = mu^2 / sd^2 and rate b = mu / sd^2, in mu and sd (N x R each).
# Only d qgamma / d shape has no closed form; it uses central differences.
.fcr_rp_gamma_derivs <- function(u, mu, sd) {
  a <- mu^2 / sd^2
  b <- mu / sd^2
  h <- .Machine$double.eps^(1 / 3) * a
  q <- function(shape) suppressWarnings(
    vapply(u, function(ur) stats::qgamma(ur, shape = shape),
           numeric(length(mu))))
  Q <- q(a)
  dQ <- (q(a + h) - q(a - h)) / (2 * h)
  if (!is.matrix(Q)) {
    Q <- matrix(Q, nrow = length(mu))
    dQ <- matrix(dQ, nrow = length(mu))
  }
  beta <- Q / b
  list(dmu = 2 * dQ - beta / mu,
       dsd = 2 * (beta - mu * dQ) / sd)
}
