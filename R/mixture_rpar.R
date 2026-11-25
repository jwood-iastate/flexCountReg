# Random parameters within the components (classes) of finite mixtures:
# countreg_mix() and countreg_lc() with `rpar_formula`.
#
# For unit u (a row, or a panel with `group_var`) and component k,
#   L_uk = (1 / R) sum_r prod_{i in u} f_k(y_i | x_i, beta_k + b_r),
# with b_r the r-th draw of the random deviations of the unit (the same
# draws for every component), and
#   L_u = sum_k pi_uk L_uk.
# The coefficient of a random variable in component k is the mean (the
# location, for lognormal parameters) of its random parameter, so the
# components differ in the means of the random parameters, and in their
# scales unless shared$rpar = TRUE. Each component's L_uk is the simulated
# panel likelihood of countreg_rp(), evaluated by its kernel
# (.fcr_rp_eval() and .fcr_rp_score() in rp_gradient.R): with unit weights
# w_u tau_uk the same kernel gives the EM objective of component k and its
# gradient.

# Most EM iterations per start with random parameters (see .fcr_fm_em()).
.fcr_fm_rp_em_maxit <- 25L

# Validate the random-parameter arguments (NULL without random parameters).
.fcr_fm_rp_settings <- function(a) {
  f <- a$rpar_formula
  if (is.null(f)) {
    if (!is.null(a$rpardists) || isTRUE(a$correlated)) {
      stop("`rpardists` and `correlated` need `rpar_formula`.",
           call. = FALSE)
    }
    return(NULL)
  }
  f <- tryCatch(stats::as.formula(f), error = function(e) NULL)
  if (is.null(f) || length(f) != 2L) {
    stop("`rpar_formula` must be a one-sided formula such as ~ x.",
         call. = FALSE)
  }
  if (.fcr_has_smooth(f)) {
    stop("Smooth terms are not available in `rpar_formula`.", call. = FALSE)
  }
  if (.fcr_has_formula_offset(f)) {
    stop("`rpar_formula` cannot have an offset.", call. = FALSE)
  }
  if (!is.logical(a$correlated) || length(a$correlated) != 1L ||
      is.na(a$correlated)) {
    stop("`correlated` must be TRUE or FALSE.", call. = FALSE)
  }
  .fcr_check_controls(ndraws = a$rpar_draws)
  if (a$correlated && !is.null(a$rpardists) && any(unlist(a$rpardists) != "n")) {
    stop("Correlated random parameters are normal; drop `rpardists` or use ",
         "correlated = FALSE.", call. = FALSE)
  }
  list(formula = f, dists = a$rpardists, correlated = a$correlated,
       draws = as.integer(a$rpar_draws),
       intercept = .fcr_explicit_intercept(f))
}

# Sign changes (+1/-1 per parameter) that make the diagonal of a packed
# Cholesky factor positive (layout RS with pairs k, l and positions idx):
# changing the signs of column l of L leaves L z with the same
# distribution, so only the signs relative to the diagonal are identified.
.fcr_chol_signs <- function(theta, RS) {
  sg <- rep(1, length(theta))
  pr <- RS$pairs
  dg <- theta[RS$idx][pr[, "k"] == pr[, "l"]]
  col_sign <- stats::setNames(ifelse(dg < 0, -1, 1),
                              pr[pr[, "k"] == pr[, "l"], "l"])
  sg[RS$idx] <- col_sign[as.character(pr[, "l"])]
  sg
}

# The random-parameter design of a formula: its model matrix without the
# intercept column unless the formula has the term 1 itself (NULL stays
# NULL; a matrix with no columns left becomes NULL).
.fcr_rpar_matrix <- function(M, f) {
  if (is.null(M) || is.null(f)) return(NULL)
  if (!.fcr_explicit_intercept(f)) {
    M <- M[, colnames(M) != "(Intercept)", drop = FALSE]
  }
  if (!ncol(M)) NULL else M
}

# Does the right side of a one-sided formula contain the term 1 itself
# (~ 1, ~ 1 + x)? An implicit intercept (~ x) does not count: in
# `rpar_formula` it would make every intercept random by accident.
.fcr_explicit_intercept <- function(f) {
  has_one <- function(e) {
    if (is.numeric(e)) return(identical(as.numeric(e), 1))
    if (is.call(e) && identical(e[[1L]], as.name("+"))) {
      return(any(vapply(as.list(e)[-1L], has_one, logical(1))))
    }
    if (is.call(e) && identical(e[[1L]], as.name("("))) return(has_one(e[[2L]]))
    FALSE
  }
  has_one(f[[length(f)]])
}

# The random parameters of the training problem: their columns in the mean
# design X, distributions and per-unit draws (the same draws for every
# component).
.fcr_fm_rp_setup <- function(st, prep, X, unit) {
  rp <- st$rp
  if (is.null(rp)) return(NULL)
  if (.fcr_prep_has_smooth(prep)) {
    stop("Smooth terms are not yet available in mixtures with random ",
         "parameters.", call. = FALSE)
  }
  cols <- colnames(prep$aux$rpar)
  if (!rp$intercept) cols <- setdiff(cols, "(Intercept)")
  if (!length(cols)) {
    stop("`rpar_formula` has no random parameters. For a random intercept ",
         "write ~ 1 (or ~ 1 + x).", call. = FALSE)
  }
  miss <- setdiff(cols, colnames(X))
  if (length(miss)) {
    stop("Random parameters must also be in `formula`: the coefficient of ",
         "each in a component is the mean of its random parameter. Not in ",
         "`formula`: ", toString(miss), ".", call. = FALSE)
  }
  dists <- .fcr_match_rpardists(rp$dists, cols, rp$correlated)
  j <- match(cols, colnames(X))
  n <- nrow(X)
  grp <- if (is.null(unit)) seq_len(n) else unit
  h <- .fcr_unit_halton(grp, rp$draws, length(cols), scrambled = TRUE,
                        seed = rp$seed)
  out <- list(names = cols, cols = j, fixed = setdiff(seq_len(ncol(X)), j),
              dists = dists, correlated = rp$correlated, R = rp$draws,
              seed = rp$seed, group = grp,
              draws = .fcr_prepare_draws(h, dists, n, rp$correlated))
  # Class-specific draw signs of a fit (see .fcr_fm_rp_set_signs())
  if (!is.null(rp$signs)) out <- .fcr_fm_rp_set_signs(out, rp$signs)
  out
}

# Draw signs of each class: a list (one per class) of +1/-1 per random
# parameter, multiplying that class's standardized draws. All classes start
# with the same draws (signs 1). After a fit, the columns of each class's
# Cholesky factor are re-signed so that its diagonal is positive (see
# .fcr_fm_rp_signs()); changing the sign of column l of L and of the draws
# z_l together leaves L z, and so the simulated likelihood, unchanged.
# Different classes may need different column signs, so the signs belong to
# the class. Signs that are all 1 leave the shared draws.
.fcr_fm_rp_set_signs <- function(rp, signs) {
  rp$signs <- signs
  rp$class_draws <- lapply(signs, function(s) {
    if (is.null(s) || all(s == 1)) return(NULL)
    .fcr_fm_rp_sign_draws(rp$draws, s)
  })
  rp
}

# Prepared draws with the standardized draws of parameter l multiplied by
# s[l] (gamma parameters use their uniform draws and keep their signs: their
# scales are not sign-identified).
.fcr_fm_rp_sign_draws <- function(dr, s) {
  for (l in seq_along(s)) {
    if (s[l] != 1 && !is.null(dr$Z[[l]])) dr$Z[[l]] <- -dr$Z[[l]]
  }
  dr
}

# Column sign changes (+1/-1 per random parameter) of each class that make
# the diagonal of its Cholesky factor positive (all 1 without correlation:
# uncorrelated scales enter the likelihood through their absolute values).
.fcr_fm_rp_col_signs <- function(theta, spec) {
  rp <- spec$rp
  Kr <- length(rp$cols)
  lapply(seq_len(spec$K), function(k) {
    if (!rp$correlated) return(rep(1, Kr))
    idx <- spec$map$rv_idx[[k]]
    ci <- .fcr_chol_index(Kr)
    dg <- unname(theta[idx][ci[, 1] == ci[, 2]])
    ifelse(dg < 0, -1, 1)
  })
}

# The countreg_rp() kernel's design for component k, with unit weights w
# (1 by default: the unweighted unit log-likelihoods).
.fcr_fm_rp_design <- function(spec, k, w = NULL) {
  rp <- spec$rp
  des <- spec$des
  map <- spec$map
  fam <- map$families[k]
  n1 <- length(map$d1_idx[[k]])
  n2 <- length(map$d2_idx[[k]])
  list(family = fam, probFunc = get_probFunc(fam), engine = spec$engine,
       y = spec$y, N_obs = spec$n,
       X_Fixed = des$X[, rp$fixed, drop = FALSE], X_offset = des$offset,
       X_rand = des$X[, rp$cols, drop = FALSE],
       X_dis_1 = if (n1) des$D1, X_dis_2 = if (n2) des$D2,
       X_het_mean = NULL, X_het_var = NULL,
       N_fixed = length(rp$fixed), N_rand = length(rp$cols),
       n_var_params = length(map$rv_idx[[k]]), N_het_mean = 0, N_het_var = 0,
       N_dis_1 = n1, N_dis_2 = n2,
       draws = rp$class_draws[[k]] %||% rp$draws, rpardists = rp$dists,
       rpar = rp$names, correlated = rp$correlated,
       dist_haltons = spec$haltons,
       normed_dist_haltons = spec$normed_haltons,
       panel_group = rp$group, panel_idx = rp$group,
       w_panel = w %||% rep(1, spec$n_u))
}

# Component k's parameters in the kernel's order (fixed coefficients,
# random-parameter means, scales, distribution parameters) and their
# positions in theta.
.fcr_fm_rp_par <- function(theta, spec, k) {
  map <- spec$map
  mi <- map$mean_idx[[k]]
  idx <- c(mi[spec$rp$fixed], mi[spec$rp$cols], map$rv_idx[[k]],
           map$d1_idx[[k]], map$d2_idx[[k]])
  list(p = unname(theta[idx]), idx = idx)
}

# Unit log-likelihoods of every component (n_u x K) and the kernel states
# (for the gradient). NULL when a component is invalid.
.fcr_fm_rp_A <- function(theta, spec) {
  K <- spec$K
  A <- matrix(NA_real_, spec$n_u, K)
  states <- vector("list", K)
  for (k in seq_len(K)) {
    pk <- .fcr_fm_rp_par(theta, spec, k)
    st <- .fcr_rp_eval(pk$p, .fcr_fm_rp_design(spec, k))
    if (!isTRUE(st$valid)) return(NULL)
    A[, k] <- st$ll
    states[[k]] <- st
  }
  list(A = A, states = states)
}

# EM objective of the components `comps` (the weighted simulated
# log-likelihoods sum_u w_u tau_uk log L_uk) at theta[idx] = par, or its
# gradient in par.
.fcr_fm_rp_Q <- function(par, idx, theta, spec, tau, comps, grad = FALSE) {
  theta[idx] <- par
  total <- 0
  g <- numeric(spec$map$n_par)
  for (k in comps) {
    pk <- .fcr_fm_rp_par(theta, spec, k)
    d <- .fcr_fm_rp_design(spec, k, spec$w_unit * tau[, k])
    st <- .fcr_rp_eval(pk$p, d)
    if (!isTRUE(st$valid)) return(if (grad) rep(NA_real_, length(idx)) else
      NA_real_)
    v <- sum(st$ll)
    if (!is.finite(v)) return(if (grad) rep(NA_real_, length(idx)) else
      NA_real_)
    total <- total + v
    if (grad) {
      S <- colSums(.fcr_rp_score(st, d))
      for (j in seq_along(pk$idx)) g[pk$idx[j]] <- g[pk$idx[j]] + S[j]
    }
  }
  if (grad) g[idx] else total
}

# Per-unit gradient of the observed log-likelihood from the components'
# kernel scores weighted by w_u tau_uk (the mixing part is added by the
# caller).
.fcr_fm_rp_grad <- function(e, spec, G) {
  for (k in seq_len(spec$K)) {
    pk <- .fcr_fm_rp_par(e$theta, spec, k)
    d <- .fcr_fm_rp_design(spec, k, spec$w_unit * e$tau[, k])
    S <- .fcr_rp_score(e$rp_states[[k]], d)
    S[is.na(S)] <- 0
    for (j in seq_along(pk$idx)) G[, pk$idx[j]] <- G[, pk$idx[j]] + S[, j]
  }
  G
}

# Average expected count of each component over the training rows, with
# the random parameters integrated out: exactly (.fcr_fm_rp_moments()) when
# the family's mean is proportional to mu, otherwise over the estimation
# draws. Inf when the expected count of a row does not exist.
.fcr_fm_rp_avg_mean <- function(theta, spec) {
  pars <- .fcr_fm_unpack(theta, spec$map)
  w <- spec$w_row
  pos <- w > 0
  vapply(seq_len(spec$K), function(k) {
    fam <- spec$map$families[k]
    cv <- .fcr_fm_comp_values(spec$des, pars, k)
    mom <- .fcr_fm_rp_moments(spec$rp, spec$des, pars, k)
    if (.fcr_mean_is_linear(fam, cv$alpha)) {
      mu_f <- exp(as.vector(spec$des$X[, spec$rp$fixed, drop = FALSE] %*%
                              pars$beta[[k]][spec$rp$fixed]) +
                    spec$des$offset)
      m <- suppressWarnings(.fcr_family_mean(fam, mu_f * exp(mom$logf),
                                             cv$alpha, cv$sigma))
    } else {
      pk <- .fcr_fm_rp_par(theta, spec, k)
      st <- .fcr_rp_eval(pk$p, .fcr_fm_rp_design(spec, k))
      if (!isTRUE(st$valid)) return(NA_real_)
      R <- ncol(st$mu_mat)
      m <- suppressWarnings(.fcr_family_mean(
        fam, as.vector(st$mu_mat),
        if (!is.null(cv$alpha)) rep(cv$alpha, R),
        if (!is.null(cv$sigma)) rep(cv$sigma, R)))
      m <- rowMeans(matrix(m, spec$n, R))
    }
    m[!is.na(mom$mean_ok) & !mom$mean_ok] <- Inf
    sum(w[pos] * m[pos], na.rm = TRUE) / sum(w)
  }, numeric(1))
}

# Starting values of the random-parameter scales: 0.1 (a Cholesky factor
# 0.1 I when correlated); gamma parameters need a positive mean.
.fcr_fm_rp_start <- function(theta, spec) {
  rp <- spec$rp
  map <- spec$map
  Kr <- length(rp$cols)
  v <- if (rp$correlated) .fcr_chol_pack(diag(0.1, Kr)) else rep(0.1, Kr)
  for (k in seq_len(spec$K)) {
    theta[map$rv_idx[[k]]] <- v
    is_g <- rp$dists == "g"
    if (any(is_g)) {
      mi <- map$mean_idx[[k]][rp$cols[is_g]]
      theta[mi] <- pmax(abs(theta[mi]), 0.1)
      if (!rp$correlated) {
        theta[map$rv_idx[[k]][is_g]] <- pmax(0.5 * theta[mi], 0.05)
      }
    }
  }
  theta
}

# ---- Prediction ------------------------------------------------------------------

# Number of draws for predictions: at least this many (more than the
# estimation draws need: one set of draws serves every row).
.fcr_fm_rp_pred_draws <- 500L

# Expected counts of each row and component at R draws of the random
# parameters (n x R): plain Halton points, the same for every row, so a
# panel's rows share their draws; or, with `draws` (from
# .fcr_fm_rp_train_draws()), the estimation draws of the rows' training
# units.
.fcr_fm_rp_mu_draws <- function(rp, des, pars, k, draws = NULL) {
  cols <- match(rp$names, colnames(des$X))
  fixed <- setdiff(seq_len(ncol(des$X)), cols)
  b <- pars$beta[[k]]
  mu_f <- exp(as.vector(des$X[, fixed, drop = FALSE] %*% b[fixed]) +
                des$offset)
  # Without training draws: a Gauss-Hermite grid (.fcr_normal_nodes()), as
  # uniforms Phi(z) for generate_random_draws(), with its weights
  w <- NULL
  H <- draws[[k]]
  if (is.null(H)) {
    nd <- .fcr_normal_nodes(length(cols))
    H <- .fcr_fm_rp_node_draws(nd$Z, rp, nrow(des$X))
    w <- nd$w
  }
  xb <- generate_random_draws(H, b[cols], pars$rv[[k]], rp$dists, rp$names,
                              des$X[, cols, drop = FALSE],
                              correlated = rp$correlated)$xb_rand_mat
  list(mu_draws = mu_f * exp(xb), mu_fixed = mu_f, draw_w = w)
}

# The estimation draws (with each class's signs) of rows whose training
# units are `tu` (one unit index per row; rows of a panel share its draws):
# a list with the prepared draws of each of the K classes. They are
# regenerated from the fit's draw seed and number of training units, as in
# .fcr_fm_rp_setup().
.fcr_fm_rp_train_draws <- function(rp, tu, K) {
  h <- .fcr_unit_halton(seq_len(rp$n_units), rp$draws, length(rp$names),
                        scrambled = TRUE, seed = rp$seed)
  dr <- .fcr_prepare_draws(.fcr_draw_rows(h, tu), rp$dists, length(tu),
                           rp$correlated)
  lapply(seq_len(K), function(k) {
    s <- rp$signs[[k]]
    if (is.null(s) || all(s == 1)) dr else .fcr_fm_rp_sign_draws(dr, s)
  })
}

# Exact moments of component k over its random parameters. With
# mu_i(b) = mu_fixed_i exp(x_i' b), E[mu_i(b)] = mu_fixed_i exp(logf_i),
# where logf_i is the log moment generating function of x_i' b at 1: the
# sum of the parameters' log MGFs (rp_log_mgf(); Wood and Gayah, 2025), or
# x' m + ||x' L||^2 / 2 for correlated normal parameters. `mean_ok` says
# whether the expected count exists (not for a lognormal parameter with
# x > 0, nor a gamma parameter with x sd^2 / mean >= 1), `var_ok` whether
# the second moment E[mu_i(b)^2] (the MGF at 2 x) does.
.fcr_fm_rp_moments <- function(rp, des, pars, k) {
  cols <- match(rp$names, colnames(des$X))
  Xr <- des$X[, cols, drop = FALSE]
  b <- pars$beta[[k]][cols]
  v <- pars$rv[[k]]
  n <- nrow(Xr)
  if (rp$correlated) {
    L <- .fcr_chol_unpack(v, length(cols))
    logf <- as.vector(Xr %*% b) + 0.5 * rowSums((Xr %*% L)^2)
    ok <- !is.na(logf)
    return(list(logf = logf, mean_ok = ok, var_ok = ok))
  }
  logf <- numeric(n)
  mean_ok <- var_ok <- rep(TRUE, n)
  for (j in seq_along(cols)) {
    dist <- rp$dists[[j]]
    t <- Xr[, j]
    s <- abs(v[j])
    e1 <- .fcr_rp_mean_exists(dist, t, b[j], s)
    e2 <- .fcr_rp_mean_exists(dist, 2 * t, b[j], s)
    lj <- rep(Inf, n)
    ok <- !is.na(e1) & e1
    if (any(ok)) {
      lj[ok] <- if (dist == "ln") .fcr_fm_ln_log_mgf(t[ok], b[j], s) else
        rp_log_mgf(dist, t[ok], b[j], s)
    }
    lj[is.na(e1)] <- NA
    logf <- logf + lj
    mean_ok <- mean_ok & e1
    var_ok <- var_ok & e2
  }
  list(logf = logf, mean_ok = mean_ok, var_ok = var_ok)
}

# log E[exp(t exp(m + s Z))], Z standard normal, for a lognormal random
# parameter with t <= 0 (Inf for t > 0, where it does not exist). The
# integrand exp(t e^(m + s z)) phi(z) is integrated by Gauss-Hermite
# quadrature centered at the mode of its logarithm,
#   h(z) = t e^(m + s z) - z^2 / 2,
# and scaled by its curvature (adaptive Gauss-Hermite): with
# u = -s z_hat, u e^u = -t e^m s^2 and h''(z_hat) = -(1 + u), so the
# integrand is close to a normal density in the rescaled variable and 40
# nodes give the integral to near machine precision.
.fcr_fm_ln_log_mgf <- function(t, m, s, nodes = 40L) {
  n <- length(t)
  m <- rep_len(m, n)
  s <- rep_len(abs(s), n)
  out <- rep(Inf, n)
  out[is.na(t) | is.na(m) | is.na(s)] <- NA
  zero <- !is.na(out) & (t == 0 | s == 0)
  out[zero] <- t[zero] * exp(m[zero])
  neg <- !is.na(out) & !zero & t < 0
  if (!any(neg)) return(out)
  t <- t[neg]
  m <- m[neg]
  s <- s[neg]
  # u + log(u) = c, with c = log(-t e^m s^2), by Newton's method
  cc <- log(-t) + m + 2 * log(s)
  u <- log1p(exp(pmin(cc, 700)))
  u[cc > 700] <- cc[cc > 700]
  for (i in seq_len(100L)) {
    step <- (u + log(u) - cc) / (1 + 1 / u)
    u <- pmax(u - step, u / 10)
    if (all(abs(step) <= 1e-12 * pmax(1, u))) break
  }
  z_hat <- -u / s
  sig <- 1 / sqrt(1 + u)
  gh <- .fcr_gh_rule(nodes)
  Zq <- z_hat + outer(sig, gh$z)
  H <- t * exp(m + s * Zq) - Zq^2 / 2 +
    matrix(gh$z^2 / 2 + log(gh$w), length(t), nodes, byrow = TRUE)
  out[neg] <- log(sig) + .fcr_lse_rows(H)
  out
}

# Prepared draws (as .fcr_prepare_draws()) from standard normal nodes Z
# (nodes x K) for n rows: each distribution's standardized draw computed
# from z directly, with the upper tail for z > 0, so far nodes keep their
# precision (Phi(z) rounds to 1 beyond z of about 8.3). Gamma parameters
# use uniforms, kept just below 1.
.fcr_fm_rp_node_draws <- function(Z, rp, n) {
  K <- ncol(Z)
  R <- nrow(Z)
  dists <- if (rp$correlated) rep("n", K) else rp$dists
  rows <- function(v) matrix(v, n, R, byrow = TRUE)
  Zl <- vector("list", K)
  Ul <- vector("list", K)
  for (j in seq_len(K)) {
    z <- Z[, j]
    lo <- stats::pnorm(-abs(z))          # the smaller tail, accurately
    switch(dists[j],
      n = , ln = { Zl[[j]] <- rows(z) },
      u = { Zl[[j]] <- rows(sign(z) * (1 - 2 * lo)) },
      t = { Zl[[j]] <- rows(sign(z) * (1 - sqrt(2 * lo))) },
      g = { Ul[[j]] <- rows(pmin(pmax(ifelse(z < 0, lo, 1 - lo), 1e-300),
                                 1 - 1e-15)) })
  }
  structure(list(U = Ul, Z = Zl, R = R), class = "fcr_prepared_draws")
}

# log f_k(y) at each draw (length(rows) x R).
.fcr_fm_logpmf_draws <- function(ps, cv, family, y, rows) {
  R <- ncol(cv$mu_draws)
  sv <- list(mu = as.vector(cv$mu_draws[rows, , drop = FALSE]),
             alpha = if (!is.null(cv$alpha)) rep(cv$alpha[rows], R),
             sigma = if (!is.null(cv$sigma)) rep(cv$sigma[rows], R))
  lp <- .fcr_fm_logpmf(ps, sv, family, rep(y[rows], R))
  # An expected count that overflows (far nodes of the integration grid)
  # gives a finite count probability 0, not NA
  lp[is.na(lp) & is.infinite(sv$mu)] <- -Inf
  matrix(lp, length(rows), R)
}

# log of the row means of exp(M), or of their weighted means with weights w
# (one per column, summing to 1).
.fcr_lmean_rows <- function(M, w = NULL) {
  if (is.null(w)) return(.fcr_lse_rows(M) - log(ncol(M)))
  .fcr_lse_rows(M + matrix(log(w), nrow(M), ncol(M), byrow = TRUE))
}

# Row means of M, or weighted means with weights w.
.fcr_wmean_rows <- function(M, w = NULL) {
  if (is.null(w)) rowMeans(M) else as.vector(M %*% w)
}

# Log-likelihood of the units `gi` (rows with `use`) under each component
# (G x K, rows named by unit): the sum of the rows' log-probabilities, or
# with random parameters the log of the average over the draws of the
# product of the unit's row probabilities.
.fcr_fm_unit_A <- function(r, y, gi, use) {
  if (is.null(r$cvs[[1]]$mu_draws)) {
    LP <- .fcr_fm_row_logpmf(r, y)
    LP[!use, ] <- 0
    return(rowsum(LP, gi, reorder = TRUE))
  }
  rows <- which(use)
  A <- NULL
  for (k in seq_len(r$ps$K)) {
    cv <- r$cvs[[k]]
    LPd <- matrix(0, r$n, ncol(cv$mu_draws))
    if (length(rows)) {
      LPd[rows, ] <- .fcr_fm_logpmf_draws(r$ps, cv, r$ps$families[k], y, rows)
    }
    S <- rowsum(LPd, gi, reorder = TRUE)
    a <- .fcr_lmean_rows(S, cv$draw_w)
    A <- if (is.null(A)) matrix(a, ncol = 1, dimnames = list(rownames(S))) else
      cbind(A, a)
  }
  A
}

# Sign changes (+1/-1 per parameter) that make the random-parameter
# standard deviations, or the diagonal of each Cholesky factor, positive.
# For a Cholesky factor, changing the sign of column j changes L z into
# L D z with D = diag(+-1), and D z has the same distribution as z.
.fcr_fm_rp_signs <- function(theta, spec) {
  sg <- rep(1, length(theta))
  rp <- spec$rp
  Kr <- length(rp$cols)
  for (idx in unique(spec$map$rv_idx)) {
    if (!length(idx)) next
    if (!rp$correlated) {
      sg[idx] <- ifelse(theta[idx] < 0, -1, 1)
    } else {
      ci <- .fcr_chol_index(Kr)
      dg <- theta[idx][ci[, 1] == ci[, 2]]
      col_sign <- ifelse(dg < 0, -1, 1)
      sg[idx] <- col_sign[ci[, 2]]
    }
  }
  sg
}
