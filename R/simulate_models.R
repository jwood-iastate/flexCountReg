# simulate() for the estimators other than countreg(): random parameters
# (countreg_rp), random effects (renb, repois, poislind_re, negmultinom)
# and fixed effects (countreg_fe). Each returns a function that draws one
# vector of counts for the training rows, drawing every unit's random terms
# anew (a parametric bootstrap draws the whole data-generating process).

# The panel (unit) of each training row: the panel index stored by the fit
# (renb, repois, poislind_re and others store only that), otherwise the
# levels of the grouping variable(s) in the training rows, otherwise each
# row its own unit.
.fcr_sim_units <- function(object, d) {
  pi <- object$model$panel_index
  if (!is.null(pi) && length(pi) == nrow(d)) {
    return(as.integer(factor(pi)))
  }
  gv <- object$model$group_var %||% object$model$panel_id
  if (is.null(gv) || !all(gv %in% names(d))) return(seq_len(nrow(d)))
  as.integer(factor(interaction(d[, gv, drop = FALSE], drop = TRUE,
                                sep = ":")))
}

# countreg_rp(): each unit's random coefficients from their fitted
# distributions (correlated, heterogeneity in means and variances), then the
# counts from the family.
.fcr_sim_rp <- function(object, d) {
  P <- .fcr_rp_parts(object, d)
  n <- nrow(d)
  unit <- .fcr_sim_units(object, d)
  G <- max(unit)
  function() {
    U <- matrix(stats::runif(G * P$N_rand), G, P$N_rand)
    hd <- structure(lapply(seq_len(P$N_rand), function(k)
      matrix(U[unit, k], n, 1)), class = "fcr_unit_draws")
    xb <- generate_random_draws(
      hd, P$random_coefs_means, P$rand_var_params, P$rpardists, P$rpar,
      P$X_rand, het_mean_coefs = P$het_mean_coefs,
      X_het_mean = P$X_het_mean, het_var_coefs = P$het_var_coefs,
      X_het_var = P$X_het_var, correlated = isTRUE(P$correlated),
      het_form = P$het_form)$xb_rand_mat[, 1]
    mu <- P$mu_fixed * exp(xb)
    .fcr_family_r(P$family, n, mu, P$vec_1, P$vec_2)
  }
}

# Mean of the training rows of a random-effects panel model (exp(X beta +
# offset), as predict()).
.fcr_sim_re_mu <- function(object, d) {
  as.vector(stats::predict(object, newdata = d))
}

# repois(): Poisson with a gamma panel effect of mean 1 and variance alpha.
.fcr_sim_repois <- function(object, d) {
  mu <- .fcr_sim_re_mu(object, d)
  unit <- .fcr_sim_units(object, d)
  est <- object$model$estimate
  alpha <- exp(est[[intersect(c("lnalpha", "ln(alpha)"), names(est))[1]]])
  function() {
    nu <- stats::rgamma(max(unit), shape = 1 / alpha, rate = 1 / alpha)
    stats::rpois(length(mu), mu * nu[unit])
  }
}

# renb() (Hausman, Hall and Griliches 1984): given the panel's p_i,
# y_it ~ NB(size = lambda_it, prob = p_i), with p_i ~ Beta(a, b) and
# lambda_it = mu_it (a - 1) / b.
.fcr_sim_renb <- function(object, d) {
  mu <- .fcr_sim_re_mu(object, d)
  unit <- .fcr_sim_units(object, d)
  est <- object$model$estimate
  a <- 1 + exp(est[["ln(a-1)"]])
  b <- exp(est[["ln(b)"]])
  lambda <- mu * (a - 1) / b
  function() {
    p <- stats::rbeta(max(unit), a, b)
    stats::rnbinom(length(mu), size = lambda, prob = p[unit])
  }
}

# poislind_re(): Poisson with a Lindley(theta) panel effect, scaled to mean
# 1: y_it ~ Poisson(mu_it nu_i theta (theta + 1) / (theta + 2)).
.fcr_sim_poislind <- function(object, d) {
  mu <- .fcr_sim_re_mu(object, d)
  unit <- .fcr_sim_units(object, d)
  theta <- exp(object$model$estimate[["ln(theta)"]])
  scale <- theta * (theta + 1) / (theta + 2)
  function() {
    nu <- .fcr_rlindley_mix(max(unit), theta)
    stats::rpois(length(mu), mu * nu[unit] * scale)
  }
}

# negmultinom(): Poisson with a gamma cluster frailty of mean 1 and shape
# r_i (NM2: 1 / alpha_i; NM1: mu_iR / alpha_i; NM-P: mu_iR^(2 - p_i) /
# alpha_i), from the model's own components.
.fcr_sim_negmultinom <- function(object, d) {
  model <- object$model
  bp <- model$blueprint %||% .fcr_blueprint(object)
  fd <- forge_new_data(bp, d)
  unit <- .fcr_sim_units(object, d)
  G <- max(unit)
  one_col <- matrix(1, G, 1, dimnames = list(NULL, "(Intercept)"))
  aux <- fd$aux
  Zc <- if (is.null(aux$dis1)) one_col else
    suppressWarnings(.nm_cluster_design(aux$dis1, unit, G,
                                        "dis_param_formula_1"))
  Wc <- if (model$family != "NMP") NULL else if (is.null(aux$dis2)) one_col
    else suppressWarnings(.nm_cluster_design(aux$dis2, unit, G,
                                             "dis_param_formula_2"))
  prep <- list(y = rep(0, nrow(d)), X = fd$X, offset = fd$offset, g = unit,
               G = G, Zc = Zc, Wc = Wc, n_beta = ncol(fd$X),
               n_alpha = ncol(Zc), n_p = if (is.null(Wc)) 0L else ncol(Wc),
               Y = rep(0, G))
  cp <- .nm_components(as.numeric(model$estimate), prep, model$family)
  function() {
    U <- stats::rgamma(G, shape = cp$r, rate = cp$r)
    stats::rpois(length(cp$mu), cp$mu * U[unit])
  }
}

# countreg_fe(): the counts given the estimated group effects. Groups left
# out of estimation have an effect of -Inf (all-zero groups: their counts
# are 0) or one fitted exactly (single-observation groups: the mean is the
# observed count).
.fcr_sim_fe <- function(object, d) {
  m <- object$model
  fd <- forge_new_data(.fcr_blueprint(object), d)
  b <- m$estimate
  xb <- names(b)[names(b) %in% colnames(fd$X)]
  eta <- as.vector(fd$X[, xb, drop = FALSE] %*% b[xb]) + fd$offset
  tot <- rep(0, nrow(d))
  known <- rep(TRUE, nrow(d))
  for (k in seq_along(m$fixef_vars)) {
    a <- m$fixef[[k]]
    val <- a[match(as.character(d[[m$fixef_vars[k]]]), names(a))]
    known <- known & !is.na(val)
    tot <- tot + ifelse(is.na(val), 0, val)
  }
  mu <- exp(eta + tot)
  y <- .fcr_response(object, d)
  mu[!known] <- y[!known]
  pars <- get_params(m$family)
  n <- length(mu)
  alpha <- if (!is.null(pars[[1]])) rep(exp(b[[pars[[1]]]]), n)
  sigma <- if (!is.null(pars[[2]])) rep(exp(b[[pars[[2]]]]), n)
  function() {
    out <- numeric(n)
    pos <- mu > 0
    out[pos] <- .fcr_family_r(m$family, sum(pos), mu[pos],
                              if (!is.null(alpha)) alpha[pos],
                              if (!is.null(sigma)) sigma[pos])
    out
  }
}
