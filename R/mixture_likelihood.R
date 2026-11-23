# Finite mixtures: component distributions, the observed-data likelihood,
# posterior probabilities and scores. One evaluator serves estimation,
# prediction, residuals and simulation.
#
# Notation. Units u are rows (cross-sectional models) or panels (one
# component per panel). a_uk is the log-likelihood of unit u under component
# k (the sum of its rows' log-probabilities for a panel), log pi_uk the log
# prior (mixing) probability, and
#   l_u = log sum_k exp(log pi_uk + a_uk),
#   tau_uk = exp(log pi_uk + a_uk - l_u).
# The observed log-likelihood is sum_u w_u l_u, with w_u the frequency weight
# of the unit, outside the logarithm of the mixture.

# Row-wise log-sum-exp of a matrix (rows of -Inf give -Inf).
.fcr_lse_rows <- function(M) {
  if (ncol(M) == 1L) return(M[, 1])
  m <- do.call(pmax, lapply(seq_len(ncol(M)), function(j) M[, j]))
  fin <- is.finite(m)
  out <- m
  if (any(fin)) {
    out[fin] <- m[fin] + log(rowSums(exp(M[fin, , drop = FALSE] - m[fin])))
  }
  out
}

# Design of a set of rows: mean design X, offset, distribution-parameter
# designs D1/D2 (a column of ones when the parameter has no formula) and the
# mixing design Z (rows, not units). Built from forge_new_data() so new data
# use the training terms, factor levels and contrasts.
.fcr_fm_design <- function(spec, data) {
  fd <- forge_new_data(spec$blueprint, data)
  n <- nrow(fd$X)
  ones <- function() matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))
  list(X = fd$X, offset = fd$offset,
       D1 = if (any(spec$map$has1)) fd$aux$dis1 %||% ones(),
       D2 = if (any(spec$map$has2)) fd$aux$dis2 %||% ones(),
       Z = if (spec$mix_mode != "fixed" && spec$K > 1L)
         fd$aux$mixing %||% ones(),
       n = n)
}

# Membership (mixing) design of panels. A whole panel belongs to one class,
# so its rows must share one membership design row; otherwise the panel's
# prior would depend on which row is taken, and so on the row order. One
# validator serves fitting, prior and posterior prediction and history
# conditioning. Z: the row design (NULL for fixed or one-component mixing,
# which has nothing to check); gi: the panel index (1..G) of each row.
# Rows with missing membership values are not compared: a panel's design is
# the common value of its complete rows (NA when it has none). Returns the
# G-row panel design, or NULL when Z is NULL. `where` is added to the error
# message (e.g. " (in `history`)").
.fcr_fm_panel_design <- function(Z, gi, G = max(0L, gi), where = "") {
  if (is.null(Z)) return(NULL)
  Zu <- matrix(NA_real_, G, ncol(Z), dimnames = list(NULL, colnames(Z)))
  cc <- which(stats::complete.cases(Z))
  if (!length(cc)) return(Zu)
  ref <- cc[match(seq_len(G), gi[cc])]
  has <- !is.na(ref)
  Zu[has, ] <- Z[ref[has], ]
  Zc <- Z[cc, , drop = FALSE]
  dif <- abs(Zc - Zu[gi[cc], , drop = FALSE])
  bad <- colnames(Z)[colSums(dif > 1e-10 * pmax(1, abs(Zc))) > 0]
  if (length(bad)) {
    stop("Mixing (membership) variables must be constant within each ",
         "panel", where, ", because the whole panel belongs to one class. ",
         "Not constant: ", toString(bad), ". Use a panel-level summary ",
         "(e.g. the panel mean) as a new column if that is what you ",
         "intend.", call. = FALSE)
  }
  Zu
}

# Component k's parameters for the rows of a design: the mean parameter
# mu = exp(X beta_k + offset) (the `predicted` argument of the registry
# density) and the distribution parameters alpha = exp(D1 delta1_k), sigma =
# exp(D2 delta2_k) on the registry's scale.
.fcr_fm_comp_values <- function(des, pars, k) {
  mu <- exp(as.vector(des$X %*% pars$beta[[k]]) + des$offset)
  alpha <- if (length(pars$d1[[k]])) exp(as.vector(des$D1 %*% pars$d1[[k]]))
  sigma <- if (length(pars$d2[[k]])) exp(as.vector(des$D2 %*% pars$d2[[k]]))
  list(mu = mu, alpha = alpha, sigma = sigma)
}

# Is component `family` defined at its parameters `cv` on the rows `rows`
# (all rows when NULL)? Besides finite positive values this checks the
# family's coupled conditions (GP2: 1 + (alpha - 1) mu > 0) through
# .fcr_family_valid(). An undefined component is not probability zero: its
# registry density would quietly return 0 for every count, and another
# component with positive weight could then hide it (a nominal 50:50
# mixture of a nonexistent GP2 and a Poisson sums to 0.5). So the whole
# candidate is rejected. A legitimate zero, such as a count outside the
# finite support of a valid underdispersed GP2, is not affected: the
# component is valid and only that probability is 0.
.fcr_fm_valid <- function(cv, family, rows = NULL) {
  sub <- function(v) if (!is.null(v) && !is.null(rows)) v[rows] else v
  all(.fcr_family_valid(family, sub(cv$mu), sub(cv$alpha), sub(cv$sigma)))
}

# Log-probabilities log f_k(y) of component values `cv` (family `family`).
# Rows where the component is undefined (see .fcr_fm_valid()) and errors of
# the numerical code give NA, never a quiet probability 0: estimation then
# rejects the candidate, and predictions report NA for those rows.
.fcr_fm_logpmf <- function(spec, cv, family, y) {
  dens <- get_probFunc(family)
  n <- length(y)
  ok <- rep_len(.fcr_family_valid(family, cv$mu, cv$alpha, cv$sigma), n)
  out <- rep(NA_real_, n)
  if (!any(ok)) return(out)
  sub <- function(v) if (!is.null(v)) rep_len(v, n)[ok]
  out[ok] <- tryCatch(
    as.vector(dens(y = y[ok], predicted = sub(cv$mu), alpha = sub(cv$alpha),
                   sigma = sub(cv$sigma), haltons = spec$haltons,
                   normed_haltons = spec$normed_haltons,
                   engine = spec$engine, log = TRUE)),
    error = function(e) NA_real_)
  out
}

# log P(Y <= q) (lower.tail = TRUE) or log P(Y > q) of a component, from the
# family's own tail function when it has one, otherwise by summing the PMF
# over 0..q in one pass (with a work budget checked first).
.fcr_fm_logcdf <- function(spec, cv, family, q, lower.tail = TRUE,
                           max_work = 5e7) {
  n <- length(cv$mu)
  q <- rep_len(as.numeric(q), n)
  out <- rep(NA_real_, n)
  neg <- !is.na(q) & q < 0
  out[neg] <- if (lower.tail) -Inf else 0
  # Rows where the component is undefined stay NA (no distribution
  # function), like .fcr_fm_logpmf()
  defined <- .fcr_family_valid(family, cv$mu, cv$alpha, cv$sigma)
  out[neg & !defined] <- NA_real_
  todo <- which(!is.na(q) & q >= 0 & defined)
  if (!length(todo)) return(out)
  sub <- function(v, i) if (!is.null(v)) v[i]
  fam <- .fcr_family_or_null(family)
  if (!is.null(fam$cdf)) {
    out[todo] <- suppressWarnings(fam$cdf(q[todo], cv$mu[todo],
                                          sub(cv$alpha, todo),
                                          sub(cv$sigma, todo), lower.tail))
    return(out)
  }
  work <- sum(q[todo] + 1)
  if (work > max_work) {
    stop("The ", family, " component distribution function has no closed ",
         "form and would need ", format(work, big.mark = ","),
         " probability evaluations (the budget is ",
         format(max_work, big.mark = ","), ").", call. = FALSE)
  }
  acc <- rep(-Inf, n)
  active <- todo
  k <- 0
  while (length(active)) {
    sv <- list(mu = cv$mu[active], alpha = sub(cv$alpha, active),
               sigma = sub(cv$sigma, active))
    acc[active] <- .fcr_logadd(acc[active],
                               .fcr_fm_logpmf(spec, sv, family,
                                              rep(k, length(active))))
    active <- active[q[active] > k]
    k <- k + 1
  }
  low <- pmin(acc[todo], 0)
  out[todo] <- if (lower.tail) low else .fcr_log1mexp(low)
  out
}

# Log prior probabilities, one row per row of Z (or n rows for fixed or
# one-component models).
.fcr_fm_logprior <- function(spec, pars, Z, n) {
  K <- spec$K
  if (K == 1L) return(matrix(0, n, 1))
  if (spec$mix_mode == "fixed") {
    return(matrix(log(spec$fixed_mixing), n, K, byrow = TRUE))
  }
  eta <- Z %*% pars$gamma
  eta - .fcr_lse_rows(eta)
}

# Evaluate the mixture at theta for the training data in `spec`. Returns NULL
# when the candidate is invalid. Otherwise a list with the n x K row
# log-probabilities LP, the unit log-likelihoods A (n_u x K), log priors,
# unit log-likelihood l_u and posterior probabilities tau.
.fcr_fm_eval <- function(theta, spec, comps = seq_len(spec$K), LP = NULL) {
  pars <- .fcr_fm_unpack(theta, spec$map)
  K <- spec$K
  pos <- spec$w_row > 0
  if (is.null(LP)) LP <- matrix(NA_real_, spec$n, K)
  for (k in comps) {
    cv <- .fcr_fm_comp_values(spec$des, pars, k)
    if (!.fcr_fm_valid(cv, spec$map$families[k], pos)) return(NULL)
    LP[, k] <- .fcr_fm_logpmf(spec, cv, spec$map$families[k], spec$y)
  }
  if (anyNA(LP[pos, , drop = FALSE])) return(NULL)
  LP[!pos, ] <- 0
  A <- if (is.null(spec$unit)) LP else rowsum(LP, spec$unit, reorder = TRUE)
  logpi <- .fcr_fm_logprior(spec, pars, spec$Zu, spec$n_u)
  if (anyNA(logpi)) return(NULL)
  lc <- logpi + A
  ll <- .fcr_lse_rows(lc)
  tau <- exp(lc - ll)
  tau[!is.finite(ll), ] <- NA_real_
  list(LP = LP, A = A, logpi = logpi, ll = ll, tau = tau, pars = pars)
}

# Observed log-likelihood contributions (one per unit, times its weight).
.fcr_fm_loglik_units <- function(theta, spec) {
  e <- .fcr_fm_eval(theta, spec)
  if (is.null(e)) return(rep(NA_real_, spec$n_u))
  v <- spec$w_unit * e$ll
  v[spec$w_unit == 0] <- 0
  v
}

# Can the scores be computed analytically? (Poisson, NB1 and NB2 components.)
.fcr_fm_analytic <- function(families) {
  all(families %in% c("POISSON", "NB1", "NB2"))
}

# Per-row scores of component k's own parameters (mean, then dis1), times the
# row weights wk, from the analytic scores of countreg().
.fcr_fm_comp_scores <- function(spec, pars, k, wk) {
  fam <- spec$map$families[k]
  p <- c(pars$beta[[k]], pars$d1[[k]])
  .fcr_analytic_grad(fam, p, spec$y, spec$des$X, spec$des$offset,
                     if (fam == "POISSON") NULL else spec$des$D1, wk)
}

# Analytic per-unit gradient of the observed log-likelihood (n_u x n_par):
# d l_u / d theta_k = tau_uk * sum_{i in u} d log f_k(y_i) / d theta_k, and
# for the mixing coefficients (tau_uk - pi_uk) z_u.
.fcr_fm_grad_units <- function(theta, spec) {
  map <- spec$map
  e <- .fcr_fm_eval(theta, spec)
  G <- matrix(NA_real_, spec$n_u, map$n_par)
  if (is.null(e)) return(G)
  G[] <- 0
  tau_row <- if (is.null(spec$unit)) e$tau else e$tau[spec$unit, , drop = FALSE]
  for (k in seq_len(spec$K)) {
    S <- .fcr_fm_comp_scores(spec, e$pars, k, spec$w_row * tau_row[, k])
    if (!is.null(spec$unit)) S <- rowsum(S, spec$unit, reorder = TRUE)
    idx <- c(map$mean_idx[[k]], map$d1_idx[[k]])
    for (j in seq_along(idx)) G[, idx[j]] <- G[, idx[j]] + S[, j]
  }
  if (!is.null(map$mix_idx)) {
    pi <- exp(e$logpi)
    for (k in seq_len(spec$K - 1L)) {
      G[, map$mix_idx[k, ]] <- (spec$w_unit * (e$tau[, k] - pi[, k])) *
        spec$Zu
    }
  }
  G[spec$w_unit == 0, ] <- 0
  G
}

# ---- EM steps -------------------------------------------------------------------

# Expected complete-data log-likelihood of the component parameters at
# positions `idx` (the other entries of theta fixed), for responsibilities
# tau_row (n x K, at the rows) and components `comps`.
.fcr_fm_Q <- function(par, idx, theta, spec, tau_row, comps) {
  theta[idx] <- par
  pars <- .fcr_fm_unpack(theta, spec$map)
  pos <- spec$w_row > 0
  total <- 0
  for (k in comps) {
    cv <- .fcr_fm_comp_values(spec$des, pars, k)
    if (!.fcr_fm_valid(cv, spec$map$families[k], pos)) return(NA_real_)
    lp <- .fcr_fm_logpmf(spec, cv, spec$map$families[k], spec$y)
    wk <- spec$w_row * tau_row[, k]
    use <- wk > 0
    if (anyNA(lp[use])) return(NA_real_)
    total <- total + sum(wk[use] * lp[use])
  }
  total
}

.fcr_fm_Q_grad <- function(par, idx, theta, spec, tau_row, comps) {
  theta[idx] <- par
  pars <- .fcr_fm_unpack(theta, spec$map)
  g <- numeric(spec$map$n_par)
  for (k in comps) {
    S <- .fcr_fm_comp_scores(spec, pars, k, spec$w_row * tau_row[, k])
    ii <- c(spec$map$mean_idx[[k]], spec$map$d1_idx[[k]])
    cs <- colSums(S)
    for (j in seq_along(ii)) g[ii[j]] <- g[ii[j]] + cs[j]
  }
  g[idx]
}

# Groups of component parameters that can be updated separately: one per
# component when nothing is shared, otherwise all together.
.fcr_fm_mstep_groups <- function(map) {
  if (map$any_shared) {
    return(list(list(idx = map$comp_par, comps = seq_len(map$K))))
  }
  lapply(seq_len(map$K), function(k)
    list(idx = map$comp_specific[[k]], comps = k))
}

# Component M-step: increase the expected complete-data log-likelihood
# (a generalized EM step with at most `maxit` quasi-Newton iterations).
.fcr_fm_mstep_comp <- function(theta, spec, tau, maxit) {
  tau_row <- if (is.null(spec$unit)) tau else tau[spec$unit, , drop = FALSE]
  analytic <- spec$analytic
  for (g in .fcr_fm_mstep_groups(spec$map)) {
    if (!length(g$idx)) next
    q0 <- .fcr_fm_Q(theta[g$idx], g$idx, theta, spec, tau_row, g$comps)
    if (!is.finite(q0)) next
    fn <- function(par) {
      v <- .fcr_fm_Q(par, g$idx, theta, spec, tau_row, g$comps)
      if (is.finite(v)) -v else 1e300
    }
    gr <- if (analytic) function(par)
      -.fcr_fm_Q_grad(par, g$idx, theta, spec, tau_row, g$comps)
    opt <- tryCatch(stats::optim(theta[g$idx], fn, gr, method = "BFGS",
                                 control = list(maxit = maxit, reltol = 1e-10)),
                    error = function(e) NULL)
    if (!is.null(opt) && is.finite(opt$value) && -opt$value >= q0) {
      theta[g$idx] <- opt$par
    }
  }
  theta
}

# Mixing M-step: constant probabilities are the weighted mean
# responsibilities; covariate-dependent probabilities maximize the weighted
# multinomial-logit objective sum_u w_u sum_k tau_uk log pi_uk. Returns NULL
# when a component's estimated share is zero (an empty component).
.fcr_fm_mstep_mix <- function(theta, spec, tau, maxit) {
  map <- spec$map
  if (is.null(map$mix_idx)) return(theta)
  K <- spec$K
  w <- spec$w_unit
  if (spec$mix_mode == "constant") {
    p <- colSums(w * tau) / sum(w)
    if (any(!(p > 0))) return(NULL)
    theta[map$mix_idx[, 1]] <- log(p[-K]) - log(p[K])
    return(theta)
  }
  idx <- as.vector(map$mix_idx)
  Z <- spec$Zu
  obj <- function(par) {
    th <- theta
    th[idx] <- par
    G <- .fcr_fm_unpack(th, map)$gamma
    eta <- Z %*% G
    lp <- eta - .fcr_lse_rows(eta)
    -sum(w * rowSums(ifelse(tau > 0, tau * lp, 0)))
  }
  grad <- function(par) {
    th <- theta
    th[idx] <- par
    G <- .fcr_fm_unpack(th, map)$gamma
    eta <- Z %*% G
    pi <- exp(eta - .fcr_lse_rows(eta))
    g <- numeric(map$n_par)
    for (k in seq_len(K - 1L)) {
      g[map$mix_idx[k, ]] <- colSums((w * (tau[, k] - pi[, k])) * Z)
    }
    -g[idx]
  }
  v0 <- obj(theta[idx])
  opt <- tryCatch(stats::optim(theta[idx], obj, grad, method = "BFGS",
                               control = list(maxit = maxit, reltol = 1e-10)),
                  error = function(e) NULL)
  if (!is.null(opt) && is.finite(opt$value) && opt$value <= v0) {
    theta[idx] <- opt$par
  }
  theta
}
