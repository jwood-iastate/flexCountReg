#' Hierarchical (mixed-effects) count regression
#'
#' Fits count models with normally distributed random effects at one or more
#' grouping levels, nested or crossed, written as in lme4:
#' \code{(1 | site)} for a random intercept, \code{(1 + lnaadt | county)} for
#' correlated random intercepts and slopes, \code{(1 | county/segment)} for
#' segments nested in counties and \code{(1 | site) + (1 | year)} for crossed
#' effects. Any count family of \code{\link{countreg}()} can be used,
#' including user-defined families, and the distribution parameters can
#' depend on covariates.
#'
#' @param formula a model formula with random-effect terms. The fixed part
#'   is as in \code{\link{countreg}()} (including \code{offset()}).
#'   Group-level predictors are ordinary columns (constant within groups),
#'   and cross-level interactions are ordinary interactions.
#' @param data a data frame.
#' @param family a count family: a name from \code{\link{countreg}()}'s list,
#'   or a family made with \code{\link{new_count_family}()}.
#' @param dis_param_formula_1,dis_param_formula_2 optional one-sided formulas
#'   for the family's distribution parameters.
#' @param offset optional names of offset columns (or use \code{offset()} in
#'   the formula).
#' @param weights optional frequency weights of the rows (a column name or a
#'   numeric vector).
#' @param engine the approximation of the integral over the random effects:
#'   \code{"aghq"} (adaptive Gauss-Hermite quadrature; one grouping term
#'   only), \code{"laplace"} (any number of terms, nested or crossed), or
#'   \code{"auto"} (quadrature for a single grouping term with at most two
#'   effects, otherwise Laplace).
#' @param nAGQ number of quadrature points per effect for \code{"aghq"}
#'   (default 9 for one effect, 5 for two or three). \code{nAGQ = 1} is the
#'   Laplace approximation.
#' @param method optimizer for \code{\link[stats]{optim}()}: \code{"BFGS"}
#'   (default) or \code{"Nelder-Mead"}.
#' @param max.iters maximum number of optimizer iterations.
#' @param start.vals optional starting values for the fixed coefficients and
#'   distribution parameters (as in \code{countreg()}), optionally followed
#'   by the random-effect standard deviations and correlations.
#' @param ndraws number of Halton draws for families whose density is
#'   simulated (as in \code{countreg()}).
#' @param verbose print the optimizer's progress.
#' @param keep_data keep the training data in the fit (needed for
#'   predictions without \code{newdata}).
#'
#' @details The random effects of grouping term \eqn{l} are
#'   \eqn{b_{lj} \sim N(0, \Sigma_l)} for each of its groups \eqn{j}, and
#'   enter the log of the mean:
#'   \deqn{\log \mu_i = x_i'\beta + o_i + \sum_l z_{li}' b_{l, g_l(i)}.}
#'   Each \eqn{\Sigma_l = L_l L_l'} is estimated through its Cholesky factor
#'   with an unconstrained diagonal, so a variance can be estimated as zero
#'   (a singular fit, reported by \code{summary()}).
#'
#'   The marginal likelihood integrates the random effects out. With
#'   \code{engine = "aghq"} each group's integral uses Gauss-Hermite nodes
#'   centered at the group's mode and scaled by its curvature; this is
#'   accurate for small groups and is the default for one grouping factor.
#'   \code{engine = "laplace"} uses the Laplace approximation, as
#'   \code{lme4::glmer()} (with \code{nAGQ = 1}) and \code{glmmTMB} do, which
#'   handles crossed and nested terms with sparse matrices. The derivatives
#'   of the log-density are analytic for the Poisson, NB1 and NB2 families
#'   and numerical for the others.
#'
#'   Standard errors come from the numerical Hessian of the marginal
#'   log-likelihood. The variance components are reported as standard
#'   deviations and correlations (\code{sd_(Intercept)|site}); their Wald
#'   standard errors are rough near zero, where the likelihood-ratio test
#'   of \code{\link{re_lrtest}()} is preferable.
#'
#'   \code{\link{countreg_rp}()} fits random parameters with a free mean and
#'   non-normal distributions at one level; \code{renb()} and
#'   \code{repois()} fit conjugate random effects in closed form.
#'
#' @returns An object of class \code{c("countreg_hier", "flexCountReg")}.
#'   Methods: \code{summary()}, \code{predict()}, \code{fitted()},
#'   \code{residuals()}, \code{simulate()}, \code{fixef()}, \code{ranef()},
#'   \code{VarCorr()}, \code{coef()}, \code{vcov()}, \code{logLik()} and
#'   \code{\link{re_lrtest}()}.
#' @seealso \code{\link{predict.countreg_hier}()}, \code{\link{ranef}()}
#' @examples
#' \donttest{
#' data("washington_roads")
#' h <- countreg_hier(Total_crashes ~ lnaadt + speed50 + offset(lnlength) +
#'                      (1 | ID), data = washington_roads, family = "NB2")
#' summary(h)
#' VarCorr(h)
#' }
#' @export
countreg_hier <- function(formula, data, family = "NB2",
                          dis_param_formula_1 = NULL,
                          dis_param_formula_2 = NULL,
                          offset = NULL, weights = NULL,
                          engine = c("auto", "aghq", "laplace"),
                          nAGQ = NULL, method = "BFGS", max.iters = 1000,
                          start.vals = NULL, ndraws = 1500, verbose = FALSE,
                          keep_data = TRUE) {
  call <- match.call()
  engine <- match.arg(engine)
  method <- match.arg(method, c("BFGS", "Nelder-Mead"))
  .fcr_check_controls(max.iters, ndraws)
  rlang::check_installed("Matrix")
  family_spec <- if (inherits(family, "fcr_count_family")) family
  family <- .fcr_resolve_family(family)
  params <- get_params(family)
  .fcr_check_dis_formulas(family, params, dis_param_formula_1,
                          dis_param_formula_2)
  data <- as.data.frame(data)

  sp <- .fcr_hier_split(formula)
  if (.fcr_has_smooth(sp$fixed) || .fcr_has_smooth(dis_param_formula_1) ||
      .fcr_has_smooth(dis_param_formula_2)) {
    stop("Smooth terms are not yet available in countreg_hier().",
         call. = FALSE)
  }
  # Rows with missing values in the random-effect terms are removed first;
  # prepare_model_data() removes those missing in the fixed part
  re_vars <- unique(unlist(lapply(sp$bars, all.vars)))
  miss <- setdiff(re_vars, names(data))
  if (length(miss)) {
    stop("Variable(s) not found in `data`: ", toString(miss), call. = FALSE)
  }
  keep_re <- stats::complete.cases(data[re_vars])
  if (is.numeric(weights)) weights <- weights[keep_re]
  prep <- prepare_model_data(sp$fixed, data[keep_re, , drop = FALSE],
                             aux_formulas = list(dis1 = dis_param_formula_1,
                                                 dis2 = dis_param_formula_2),
                             offset = offset, weights = weights)
  prep$rows_used <- which(keep_re)[prep$rows_used]
  prep$n_dropped <- nrow(data) - length(prep$rows_used)
  if (prep$n_dropped > 0) {
    message(prep$n_dropped, " row(s) with missing values were removed.")
  }
  d <- data[prep$rows_used, , drop = FALSE]
  terms <- .fcr_hier_terms(sp$bars, d, sp$env)

  y <- as.numeric(prep$y)
  n <- length(y)
  X <- prep$X
  W1 <- prep$aux$dis1
  W2 <- prep$aux$dis2
  k_beta <- ncol(X)
  n_alpha <- if (is.null(params[[1]])) 0L else if (is.null(W1)) 1L else
    ncol(W1)
  n_sigma <- if (is.null(params[[2]])) 0L else if (is.null(W2)) 1L else
    ncol(W2)
  dist_names <- c(
    if (n_alpha) {
      if (is.null(W1)) params[[1]] else paste0(params[[1]], ":", colnames(W1))
    },
    if (n_sigma) {
      if (is.null(W2)) params[[2]] else paste0(params[[2]], ":", colnames(W2))
    })
  k_fix <- k_beta + n_alpha + n_sigma
  n_theta <- sum(vapply(terms, function(t) as.integer(t$ntheta), 0L))

  # Engine
  q_all <- vapply(terms, function(t) as.integer(t$q), 0L)
  if (engine == "auto") {
    engine <- if (length(terms) == 1L && q_all[1] <= 2L) "aghq" else
      "laplace"
  }
  if (engine == "aghq") {
    if (length(terms) > 1L) {
      stop("`engine = \"aghq\"` needs a single grouping term; this model ",
           "has ", length(terms), ". Use engine = \"laplace\".",
           call. = FALSE)
    }
    if (is.null(nAGQ)) nAGQ <- if (q_all[1] == 1L) 9L else 5L
    nAGQ <- as.integer(nAGQ)
    if (length(nAGQ) != 1L || is.na(nAGQ) || nAGQ < 1L) {
      stop("`nAGQ` must be a positive whole number.", call. = FALSE)
    }
    if (nAGQ^q_all[1] > 1000) {
      stop("`nAGQ`^", q_all[1], " = ", nAGQ^q_all[1], " quadrature nodes ",
           "per group is too many; use a smaller `nAGQ` or the Laplace ",
           "engine.", call. = FALSE)
    }
  } else {
    nAGQ <- 1L
  }
  rule <- .fcr_gh_grid(nAGQ, q_all[1])

  haltons <- randtoolbox::halton(ndraws)
  dens_args <- list(haltons = haltons, normed_haltons = stats::qnorm(haltons),
                    engine = "quadrature")
  logf <- .fcr_hier_logf(family, y, prep$weights, dens_args)
  Z <- .fcr_hier_Z(terms, n)
  tmpl <- .fcr_hier_lambda_template(terms)

  unpack <- function(p) {
    beta <- p[seq_len(k_beta)]
    pos <- k_beta
    alpha <- sigma <- NULL
    if (n_alpha) {
      a <- p[pos + seq_len(n_alpha)]
      alpha <- if (is.null(W1)) rep(exp(a), n) else exp(as.vector(W1 %*% a))
      pos <- pos + n_alpha
    }
    if (n_sigma) {
      s <- p[pos + seq_len(n_sigma)]
      sigma <- if (is.null(W2)) rep(exp(s), n) else exp(as.vector(W2 %*% s))
      pos <- pos + n_sigma
    }
    list(beta = beta, alpha = alpha, sigma = sigma,
         theta = p[pos + seq_len(n_theta)],
         eta0 = as.vector(X %*% beta) + prep$offset)
  }

  # Warm start of the modes across evaluations
  state <- new.env(parent = emptyenv())
  state$u <- rep(0, tmpl$Q)
  mode_at <- function(p) {
    up <- unpack(p)
    Lam <- .fcr_hier_lambda(tmpl, up$theta)
    A <- Z %*% Lam
    md <- .fcr_hier_mode(state$u, up$eta0, A, logf, up$alpha, up$sigma)
    if (is.null(md) || any(!is.finite(md$u))) {
      md <- .fcr_hier_mode(rep(0, tmpl$Q), up$eta0, A, logf, up$alpha,
                           up$sigma)
    }
    list(md = md, up = up, A = A, Lam = Lam)
  }
  loglik <- function(p) {
    r <- mode_at(p)
    if (is.null(r$md)) return(-Inf)
    ll <- if (engine == "aghq") {
      .fcr_hier_aghq(r$md, terms[[1]], r$up$eta0, r$A, logf, r$up$alpha,
                     r$up$sigma, rule)
    } else {
      .fcr_hier_laplace(r$md)
    }
    if (is.finite(ll) && all(is.finite(r$md$u)) &&
        max(abs(r$md$u)) < 50) {
      state$u <- r$md$u
    }
    ll
  }

  # ---- Starting values -------------------------------------------------------
  vc_names <- names(.fcr_hier_vc_vector(terms,
                                        .fcr_hier_theta_start(terms)))
  fix_names <- c(colnames(X), dist_names)
  theta0 <- .fcr_hier_theta_start(terms)
  if (!is.null(start.vals)) {
    if (length(start.vals) == k_fix) {
      start <- c(unname(start.vals), theta0)
    } else if (length(start.vals) == k_fix + length(vc_names)) {
      start <- c(unname(start.vals[seq_len(k_fix)]),
                 .fcr_hier_theta_from_vc(terms, start.vals[-seq_len(k_fix)]))
    } else {
      stop("`start.vals` must have ", k_fix, " values (fixed coefficients ",
           "and distribution parameters), or ", k_fix + length(vc_names),
           " with the random-effect standard deviations and correlations.",
           call. = FALSE)
    }
  } else {
    st <- tryCatch({
      m0 <- suppressMessages(suppressWarnings(countreg(
        sp$fixed, data = d, family = family_spec %||% family,
        dis_param_formula_1 = dis_param_formula_1,
        dis_param_formula_2 = dis_param_formula_2, offset = offset,
        weights = prep$weights, method = "BFGS", ndraws = ndraws,
        keep_data = FALSE)))
      unname(m0$model$estimate[seq_len(k_fix)])
    }, error = function(e) NULL)
    if (is.null(st) || any(!is.finite(st))) {
      g0 <- stats::glm.fit(X, y, weights = prep$weights, offset = prep$offset,
                           family = stats::poisson())
      b0 <- g0$coefficients
      b0[is.na(b0)] <- 0
      st <- c(b0, rep_len(.fcr_get_family(family)$start, n_alpha + n_sigma))
    }
    start <- c(st, theta0)
  }
  int_names <- c(fix_names, paste0("theta", seq_len(n_theta)))
  names(start) <- int_names

  # ---- Optimization ----------------------------------------------------------
  nll <- function(p) {
    v <- -loglik(p)
    if (!is.finite(v)) 1e100 else v
  }
  ngrad <- function(p) {
    h <- 1e-5 * pmax(1, abs(p))
    vapply(seq_along(p), function(i) {
      e <- numeric(length(p))
      e[i] <- h[i]
      (nll(p + e) - nll(p - e)) / (2 * h[i])
    }, 0)
  }
  if (!is.finite(nll(start)) || nll(start) >= 1e100) {
    stop("The log-likelihood is not finite at the starting values.",
         call. = FALSE)
  }
  opt <- stats::optim(start, nll, gr = if (method == "BFGS") ngrad,
                      method = method,
                      control = list(maxit = max.iters,
                                     trace = as.integer(isTRUE(verbose)),
                                     reltol = 1e-12))
  p_hat <- opt$par
  # Newton steps on the numerical Hessian finish the optimization
  H <- NULL
  for (it in 1:5) {
    H <- .fcr_final_hessian(loglik, p_hat)
    if (is.null(H)) break
    g <- -ngrad(p_hat)
    step <- tryCatch(-solve(H, g), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) break
    if (sum(g * step) < 1e-8) break
    p_new <- p_hat + step
    if (loglik(p_new) >= loglik(p_hat)) p_hat <- p_new else break
  }
  # Report with a non-negative diagonal of each L
  up <- unpack(p_hat)
  vc <- .fcr_hier_vc_vector(terms, up$theta)
  p_hat[k_fix + seq_len(n_theta)] <- .fcr_hier_theta_from_vc(terms, vc)
  names(p_hat) <- int_names
  ll_hat <- loglik(p_hat)
  H <- .fcr_final_hessian(loglik, p_hat)
  grad_hat <- -ngrad(p_hat)

  # ---- Covariance on the reported scale (sd, cor) -----------------------------
  to_report <- function(p) {
    c(p[seq_len(k_fix)],
      .fcr_hier_vc_vector(terms, p[k_fix + seq_len(n_theta)]))
  }
  est <- to_report(p_hat)
  names(est) <- c(fix_names, vc_names)
  V_int <- if (!is.null(H)) {
    fit_h <- list(hessian = H)
    tryCatch(.fcr_hessian_vcov(fit_h), error = function(e) NULL)
  }
  if (is.null(V_int)) V_int <- matrix(NA_real_, length(p_hat), length(p_hat))
  Jac <- matrix(0, length(est), length(p_hat))
  for (i in seq_along(p_hat)) {
    h <- 1e-6 * max(1, abs(p_hat[i]))
    e <- numeric(length(p_hat))
    e[i] <- h
    Jac[, i] <- (to_report(p_hat + e) - to_report(p_hat - e)) / (2 * h)
  }
  V <- Jac %*% V_int %*% t(Jac)
  dimnames(V) <- list(names(est), names(est))

  # ---- Conditional modes and fitted values -------------------------------------
  state$u <- rep(0, tmpl$Q)
  fin <- mode_at(p_hat)
  b <- as.vector(fin$Lam %*% fin$md$u)
  theta_diag <- unlist(lapply(terms, function(t) {
    L <- .fcr_hier_L(p_hat[k_fix + t$toff + seq_len(t$ntheta)], t$q)
    abs(diag(L))
  }))
  singular <- any(theta_diag < 1e-4)
  # At a zero variance or a correlation of +-1 the Wald standard errors of
  # that term's variance components are not meaningful
  pos <- k_fix
  for (t in terms) {
    nv <- t$q + t$q * (t$q - 1L) / 2L
    ix <- pos + seq_len(nv)
    v <- est[ix]
    if (any(v[seq_len(t$q)] < 1e-4) ||
        (t$q > 1L && any(abs(v[-seq_len(t$q)]) > 1 - 1e-4))) {
      V[ix, ] <- NA_real_
      V[, ix] <- NA_real_
    }
    pos <- pos + nv
  }
  eta <- fin$md$eta
  mu <- .fcr_family_mean(family, exp(eta), up$alpha, up$sigma)

  conv <- opt$convergence == 0L && all(is.finite(grad_hat)) &&
    max(abs(grad_hat)) < 1e-2
  if (!conv) {
    warning("The optimizer may not have converged (code ", opt$convergence,
            ", largest gradient ", signif(max(abs(grad_hat)), 3), ").",
            call. = FALSE)
  }
  if (singular) {
    message("Singular fit: a random-effect variance is estimated at or ",
            "near zero (see summary()).")
  }

  terms_store <- lapply(terms, function(t) {
    t$X <- NULL
    t$g <- NULL
    t
  })
  fit <- list(
    estimate = est, maximum = ll_hat, vcov = V, vcov_type = "hessian",
    stderr = "normal", se = .fcr_se(V), hessian_internal = H,
    par_internal = p_hat, gradient = grad_hat, iterations = opt$counts,
    code = opt$convergence, converged = conv, method = method,
    fcr_method = method,
    family = family, model = family, family_spec = family_spec,
    params = params, formula = formula, fixed_formula = prep$formula,
    dis_param_formula_1 = prep$aux_formulas$dis1,
    dis_param_formula_2 = prep$aux_formulas$dis2,
    offset = prep$blueprint$offset_label, blueprint = prep$blueprint,
    re_terms = terms_store, k_beta = k_beta, n_alpha = n_alpha,
    n_sigma = n_sigma, k_fix = k_fix, n_theta = n_theta,
    engine = engine, nAGQ = nAGQ, ndraws = ndraws,
    u = fin$md$u, b = b, singular = singular,
    n_groups = stats::setNames(vapply(terms, function(t) as.integer(t$J), 0L),
                               vapply(terms, `[[`, "", "group")),
    y = y, observed = y, nobs = n, n_bic = n, weights = prep$weights,
    rows_used = prep$rows_used, n_dropped = prep$n_dropped,
    eta = eta, predictions = mu, fitted.values = mu, residuals = y - mu,
    alpha = up$alpha, sigma = up$sigma, beta_pred = up$beta,
    n_beta = k_beta, modelType = "countreg_hier", condvar = NULL)
  # Posterior covariance of the effects (from the curvature at the mode)
  fit$condvar <- tryCatch(.fcr_hier_condvar(fin, terms), error = function(e)
    NULL)
  fit$lr_sample <- .fcr_lr_sample(prep, data, NULL)

  obj <- .createFlexCountReg(model = fit, data = data, keep_data = keep_data,
                             call = call, formula = formula)
  obj
}

# Posterior covariance of each group's effects b_j = L u_j:
# L (M^{-1})_jj L', with M the negative Hessian of h at the mode. Returns,
# per term, an array q x q x J (NULL for very large models).
.fcr_hier_condvar <- function(fin, terms) {
  Q <- ncol(fin$md$M)
  if (Q > 3000) return(NULL)
  Minv <- as.matrix(Matrix::solve(fin$md$ch, Matrix::Diagonal(Q)))
  Lam <- as.matrix(fin$Lam)
  lapply(terms, function(t) {
    q <- t$q
    out <- array(0, c(q, q, t$J))
    for (j in seq_len(t$J)) {
      ix <- t$uoff + (j - 1L) * q + seq_len(q)
      Lj <- Lam[ix, ix, drop = FALSE]
      out[, , j] <- Lj %*% Minv[ix, ix, drop = FALSE] %*% t(Lj)
    }
    out
  })
}
