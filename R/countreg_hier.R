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
#'   is as in \code{\link{countreg}()} (including \code{offset()}), and can
#'   have penalized smooth terms \code{s()} and \code{t2()} of mgcv.
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
#'   only), \code{"laplace"} (any number of terms, nested or crossed, and
#'   smooth terms), \code{"msl"} (maximum simulated likelihood; one grouping
#'   term, any \code{re_dist}), or \code{"auto"} (simulation for non-normal
#'   effects, quadrature for a single grouping term with at most two
#'   effects, otherwise Laplace).
#' @param nAGQ number of quadrature points per effect for \code{"aghq"}
#'   (default 9 for one effect, 5 for two or three). \code{nAGQ = 1} is the
#'   Laplace approximation.
#' @param re_dist distribution of the random effects of the grouping term:
#'   \code{"normal"} (default), \code{"uniform"}, \code{"triangular"}
#'   (symmetric) or \code{"gamma"}; one value, or one per effect. Effects
#'   other than normal need a single grouping term and use
#'   \code{engine = "msl"}. With \code{"gamma"} (one effect) the
#'   multiplicative group effect \eqn{\exp(b_j)} is gamma with mean 1 and
#'   standard deviation \code{sd_...}, as in the negative binomial
#'   random-effects model; with the others \eqn{b_j} has mean 0 and the
#'   reported standard deviation.
#' @param re_draws number of (randomly shifted) Halton draws per group for
#'   \code{engine = "msl"}.
#' @param cre optional one-sided formula of time-varying regressors whose
#'   group means are added to the mean model (correlated random effects; see
#'   \code{\link{cre_test}()}).
#' @param cre_group the grouping column(s) for the means of \code{cre};
#'   by default those of the random-effect terms, when all terms share them.
#' @param seed seed of the random shifts of the Halton draws of
#'   \code{engine = "msl"}.
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
#'   Smooth terms are fitted in their mixed-model form: the penalized part
#'   of each smooth is a block of independent normal coefficients with
#'   variance \eqn{1/\lambda}, so the smoothing parameters are estimated by
#'   the (Laplace) marginal likelihood, as in gamm4. Their coefficients are
#'   reported on mgcv's basis, with the covariance of the coefficients and
#'   effects given the variance parameters, so \code{plot_smooths()} and the
#'   other smooth tools work as for \code{countreg()}.
#'
#'   With \code{engine = "msl"} the group likelihoods are averages over
#'   \code{re_draws} draws of each group's effects, so effects can have any
#'   of the \code{re_dist} distributions; the effects of the groups
#'   (\code{ranef()}) are then posterior means.
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
                          engine = c("auto", "aghq", "laplace", "msl"),
                          nAGQ = NULL, re_dist = "normal", re_draws = 500,
                          cre = NULL, cre_group = NULL,
                          method = "BFGS", max.iters = 1000,
                          start.vals = NULL, ndraws = 1500, seed = 123,
                          verbose = FALSE, keep_data = TRUE) {
  call <- match.call()
  engine <- match.arg(engine)
  method <- match.arg(method, c("BFGS", "Nelder-Mead"))
  .fcr_check_controls(max.iters, ndraws)
  rlang::check_installed("Matrix")
  family_spec <- if (inherits(family, "fcr_count_family")) family
  family <- .fcr_resolve_family(family)
  data <- as.data.frame(data)

  # Correlated random effects: group means of time-varying regressors
  cre_info <- NULL
  if (!is.null(cre)) {
    grp <- cre_group %||% .fcr_hier_cre_group(formula)
    cp <- .fcr_cre_prepare(cre, formula, data, grp,
                           formulas = list(dis_param_formula_1,
                                           dis_param_formula_2),
                           extra_vars = c(if (is.character(offset)) offset,
                                          if (is.character(weights)) weights))
    data <- cp$data
    formula <- cp$formula
    cre_info <- cp$info
  }

  S <- .fcr_hier_setup(formula, data, family, family_spec,
                       dis_param_formula_1, dis_param_formula_2, offset,
                       weights, engine, nAGQ, re_dist, re_draws, ndraws,
                       seed)
  if (S$prep$n_dropped > 0) {
    message(S$prep$n_dropped, " row(s) with missing values were removed.")
  }

  # ---- Starting values -------------------------------------------------------
  start <- .fcr_hier_start(S, start.vals, ndraws)

  # ---- Optimization ----------------------------------------------------------
  opt <- .fcr_hier_optimize(S, start, method, max.iters, verbose)
  p_hat <- opt$par

  fit <- .fcr_hier_finish(S, p_hat, opt)
  fit$formula <- formula
  fit$family_spec <- family_spec
  fit$method <- method
  fit$fcr_method <- method
  fit$cre <- cre_info
  fit$spec <- list(dis_param_formula_1 = dis_param_formula_1,
                   dis_param_formula_2 = dis_param_formula_2,
                   offset = offset,
                   weights = if (!is.null(weights)) weights,
                   engine = S$engine, nAGQ = S$nAGQ, re_dist = S$re_dist,
                   re_draws = re_draws, ndraws = ndraws, seed = seed)
  fit$lr_sample <- .fcr_lr_sample(S$prep, data, NULL)
  if (isTRUE(fit$singular)) {
    message("Singular fit: a random-effect variance is estimated at or ",
            "near zero (see summary()).")
  }
  .createFlexCountReg(model = fit, data = data, keep_data = keep_data,
                      call = call, formula = formula)
}

# Grouping variables for `cre`: those of the random-effect terms, when all
# terms share them.
.fcr_hier_cre_group <- function(formula) {
  bars <- reformulas::findbars(stats::as.formula(formula))
  vs <- unique(lapply(bars, function(b) sort(all.vars(b[[3]]))))
  if (length(vs) != 1L) {
    stop("The random-effect terms have different grouping variables; give ",
         "the group of the means with `cre_group`.", call. = FALSE)
  }
  vs[[1]]
}

# Model setup: data, design, engine and the log-likelihood of the internal
# parameters p = c(beta (fixed columns), distribution parameters, theta).
.fcr_hier_setup <- function(formula, data, family, family_spec = NULL,
                            dis1 = NULL, dis2 = NULL, offset = NULL,
                            weights = NULL, engine = "auto", nAGQ = NULL,
                            re_dist = "normal", re_draws = 500,
                            ndraws = 1500, seed = 123) {
  params <- get_params(family)
  .fcr_check_dis_formulas(family, params, dis1, dis2)
  sp <- .fcr_hier_split(formula)
  if (.fcr_has_smooth(dis1) || .fcr_has_smooth(dis2)) {
    stop("Smooth terms in the distribution-parameter formulas are not yet ",
         "available in countreg_hier().", call. = FALSE)
  }
  # Rows with missing values in the random-effect terms are removed first;
  # prepare_model_data() removes those missing in the fixed part
  re_vars <- unique(unlist(lapply(sp$bars, all.vars)))
  miss <- setdiff(re_vars, names(data))
  if (length(miss)) {
    stop("Variable(s) not found in `data`: ", toString(miss), call. = FALSE)
  }
  keep_re <- if (length(re_vars)) stats::complete.cases(data[re_vars]) else
    rep(TRUE, nrow(data))
  if (is.numeric(weights)) weights <- weights[keep_re]
  prep <- prepare_model_data(sp$fixed, data[keep_re, , drop = FALSE],
                             aux_formulas = list(dis1 = dis1, dis2 = dis2),
                             offset = offset, weights = weights)
  prep$rows_used <- which(keep_re)[prep$rows_used]
  prep$n_dropped <- nrow(data) - length(prep$rows_used)
  d <- data[prep$rows_used, , drop = FALSE]
  bar_terms <- if (length(sp$bars)) .fcr_hier_terms(sp$bars, d, sp$env) else
    list()

  # Smooths of the mean: fixed null-space columns and random blocks
  ents <- prep$smooth$main
  smooth_cols <- unlist(lapply(ents, `[[`, "cols"))
  param_cols <- setdiff(seq_len(ncol(prep$X)), smooth_cols)
  X <- prep$X[, param_cols, drop = FALSE]
  smooth <- NULL
  sm_terms <- list()
  if (length(ents)) {
    smooth <- .fcr_hier_smooth_terms(ents, prep$X)
    if (!is.null(smooth$X) && ncol(smooth$X)) X <- cbind(X, smooth$X)
    sm_terms <- smooth$terms
  }
  terms <- .fcr_hier_offsets(c(bar_terms, sm_terms))
  n_bar <- length(bar_terms)

  y <- as.numeric(prep$y)
  n <- length(y)
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
  ntheta <- vapply(terms, function(t) as.integer(t$ntheta), 0L)
  n_theta <- sum(ntheta)
  n_theta_bar <- sum(ntheta[seq_len(n_bar)])

  # ---- Engine ----------------------------------------------------------------
  q_all <- vapply(bar_terms, function(t) as.integer(t$q), 0L)
  re_dist <- .fcr_hier_check_dist(re_dist, bar_terms)
  nonnormal <- any(re_dist != "normal")
  if (engine == "auto") {
    engine <- if (nonnormal) "msl" else if (n_bar == 1L &&
                 !length(sm_terms) && q_all[1] <= 2L) "aghq" else "laplace"
  }
  if (nonnormal && engine != "msl") {
    stop("Non-normal random effects (`re_dist`) need engine = \"msl\".",
         call. = FALSE)
  }
  if (engine %in% c("aghq", "msl")) {
    if (n_bar != 1L || length(sm_terms)) {
      stop("`engine = \"", engine, "\"` needs a single grouping term and ",
           "no smooth terms; use engine = \"laplace\".", call. = FALSE)
    }
  }
  if (engine == "aghq") {
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
  rule <- if (engine == "aghq") .fcr_gh_grid(nAGQ, q_all[1])
  msl <- NULL
  if (engine == "msl") {
    if (!is.numeric(re_draws) || length(re_draws) != 1L || re_draws < 10) {
      stop("`re_draws` must be a number of at least 10.", call. = FALSE)
    }
    t1 <- bar_terms[[1]]
    msl <- list(U = .fcr_hier_draws(t1$J, re_draws, t1$q, seed),
                X = t1$X, g = as.integer(t1$g))
  }

  haltons <- randtoolbox::halton(ndraws)
  dens_args <- list(haltons = haltons, normed_haltons = stats::qnorm(haltons),
                    engine = "quadrature")
  logf <- .fcr_hier_logf(family, y, prep$weights, dens_args)
  logf_for <- .fcr_hier_logf_stack(family, y, prep$weights, dens_args)
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
  msl_at <- function(p, posterior = FALSE) {
    up <- unpack(p)
    L <- .fcr_hier_L(up$theta, q_all[1])
    .fcr_hier_msl(up$eta0, msl$X, msl$g, msl$U, re_dist, L, logf_for,
                  up$alpha, up$sigma, posterior = posterior)
  }
  loglik <- function(p) {
    if (engine == "msl") return(msl_at(p)$ll)
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

  list(prep = prep, d = d, sp = sp, y = y, n = n, X = X, W1 = W1, W2 = W2,
       family = family, family_spec = family_spec, params = params,
       terms = terms, bar_terms = terms[seq_len(n_bar)], n_bar = n_bar,
       smooth = smooth, param_cols = param_cols, k_beta = k_beta,
       n_alpha = n_alpha, n_sigma = n_sigma, k_fix = k_fix,
       dist_names = dist_names, n_theta = n_theta, n_theta_bar = n_theta_bar,
       engine = engine, nAGQ = nAGQ, re_dist = re_dist, q_all = q_all,
       logf = logf, Z = Z, tmpl = tmpl, state = state, unpack = unpack,
       mode_at = mode_at, msl_at = msl_at, loglik = loglik,
       int_names = c(colnames(X), dist_names,
                     paste0("theta", seq_len(n_theta))))
}

# Effect distribution of the (single) term: one per effect.
.fcr_hier_check_dist <- function(re_dist, bar_terms) {
  if (is.null(re_dist)) re_dist <- "normal"
  if (!is.character(re_dist) || !all(re_dist %in% .fcr_re_dists)) {
    stop("`re_dist` must be one of ", toString(dQuote(.fcr_re_dists, FALSE)),
         " (one value, or one per effect of the grouping term).",
         call. = FALSE)
  }
  if (all(re_dist == "normal")) {
    return(rep("normal", max(1L, if (length(bar_terms))
      bar_terms[[1]]$q else 1L)))
  }
  if (length(bar_terms) != 1L) {
    stop("Non-normal random effects (`re_dist`) are available for a single ",
         "grouping term.", call. = FALSE)
  }
  q <- bar_terms[[1]]$q
  if (!length(re_dist) %in% c(1L, q)) {
    stop("`re_dist` must have one value or one per effect (", q, ").",
         call. = FALSE)
  }
  re_dist <- rep_len(re_dist, q)
  if (any(re_dist == "gamma") && q > 1L) {
    stop("`re_dist = \"gamma\"` is available for a term with one effect, ",
         "such as (1 | site).", call. = FALSE)
  }
  re_dist
}

# Starting values of the internal parameters.
.fcr_hier_start <- function(S, start.vals, ndraws) {
  theta0 <- c(.fcr_hier_theta_start(S$bar_terms),
              rep(1, S$n_theta - S$n_theta_bar))
  vc_names <- if (S$n_bar) names(.fcr_hier_vc_vector(
    S$bar_terms, .fcr_hier_theta_start(S$bar_terms))) else character(0)
  k_fix <- S$k_fix
  if (!is.null(start.vals)) {
    if (length(start.vals) == k_fix) {
      start <- c(unname(start.vals), theta0)
    } else if (length(start.vals) == k_fix + length(vc_names)) {
      start <- c(unname(start.vals[seq_len(k_fix)]),
                 .fcr_hier_theta_from_vc(S$bar_terms,
                                         start.vals[-seq_len(k_fix)]),
                 theta0[-seq_len(S$n_theta_bar)])
    } else {
      stop("`start.vals` must have ", k_fix, " values (fixed coefficients ",
           "and distribution parameters), or ", k_fix + length(vc_names),
           " with the random-effect standard deviations and correlations.",
           call. = FALSE)
    }
  } else {
    fam <- S$family_spec %||% S$family
    st <- tryCatch({
      m0 <- suppressMessages(suppressWarnings(countreg(
        S$sp$fixed, data = S$d, family = fam,
        dis_param_formula_1 = S$prep$aux_formulas$dis1,
        dis_param_formula_2 = S$prep$aux_formulas$dis2,
        offset = S$prep$blueprint$offset_cols, weights = S$prep$weights,
        method = "BFGS", ndraws = ndraws, keep_data = FALSE,
        sp = if (!is.null(S$smooth)) rep(1e3, length(unlist(lapply(
          S$prep$smooth$main, `[[`, "S")))))))
      est <- m0$model$estimate
      k0 <- ncol(S$prep$X)
      b0 <- est[seq_len(k0)]
      # Smooths start at their null-space part (least squares)
      bi <- b0[S$param_cols]
      if (!is.null(S$smooth)) {
        fx <- S$X[, -seq_along(S$param_cols), drop = FALSE]
        sm_fit <- S$prep$X[, -S$param_cols, drop = FALSE] %*%
          b0[-S$param_cols]
        bf <- if (ncol(fx)) stats::lm.fit(fx, sm_fit)$coefficients else
          numeric(0)
        bf[is.na(bf)] <- 0
        bi <- c(bi, bf)
      }
      unname(c(bi, est[k0 + seq_len(S$n_alpha + S$n_sigma)]))
    }, error = function(e) NULL)
    if (is.null(st) || length(st) != k_fix || any(!is.finite(st))) {
      g0 <- stats::glm.fit(S$X, S$y, weights = S$prep$weights,
                           offset = S$prep$offset, family = stats::poisson())
      b0 <- g0$coefficients
      b0[is.na(b0)] <- 0
      st <- c(b0, rep_len(.fcr_get_family(S$family)$start,
                          S$n_alpha + S$n_sigma))
    }
    start <- c(st, theta0)
  }
  names(start) <- S$int_names
  start
}

# Central-difference gradient of the negative log-likelihood.
.fcr_hier_ngrad <- function(nll) {
  function(p) {
    h <- 1e-5 * pmax(1, abs(p))
    vapply(seq_along(p), function(i) {
      e <- numeric(length(p))
      e[i] <- h[i]
      (nll(p + e) - nll(p - e)) / (2 * h[i])
    }, 0)
  }
}

# Maximize the log-likelihood (optionally with some parameters fixed at
# `start`): optim() followed by Newton steps on the numerical Hessian.
.fcr_hier_optimize <- function(S, start, method = "BFGS", max.iters = 1000,
                               verbose = FALSE, fixed = NULL, polish = TRUE) {
  free <- setdiff(seq_along(start), fixed)
  full <- function(pf) {
    p <- start
    p[free] <- pf
    p
  }
  ll_free <- function(pf) S$loglik(full(pf))
  nll <- function(pf) {
    v <- -ll_free(pf)
    if (!is.finite(v)) 1e100 else v
  }
  ngrad <- .fcr_hier_ngrad(nll)
  if (nll(start[free]) >= 1e100) {
    stop("The log-likelihood is not finite at the starting values.",
         call. = FALSE)
  }
  opt <- stats::optim(start[free], nll, gr = if (method == "BFGS") ngrad,
                      method = method,
                      control = list(maxit = max.iters,
                                     trace = as.integer(isTRUE(verbose)),
                                     reltol = 1e-12))
  pf <- opt$par
  # Newton steps on the numerical Hessian finish the optimization
  for (it in seq_len(if (polish) 5L else 0L)) {
    H <- .fcr_final_hessian(ll_free, pf)
    if (is.null(H)) break
    g <- -ngrad(pf)
    step <- tryCatch(-solve(H, g), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) break
    if (sum(g * step) < 1e-8) break
    p_new <- pf + step
    if (ll_free(p_new) >= ll_free(pf)) pf <- p_new else break
  }
  list(par = full(pf), value = ll_free(pf), counts = opt$counts,
       convergence = opt$convergence)
}

# The fitted model from the internal estimates: reported coefficients,
# covariance, effects, fitted values and smooth summaries.
.fcr_hier_finish <- function(S, p_hat, opt) {
  k_fix <- S$k_fix
  bar <- S$bar_terms
  ith <- k_fix + seq_len(S$n_theta)
  ith_bar <- k_fix + seq_len(S$n_theta_bar)
  ith_sm <- setdiff(ith, ith_bar)
  # Report with a non-negative diagonal of each L
  if (S$n_bar) {
    vc <- .fcr_hier_vc_vector(bar, p_hat[ith_bar])
    p_hat[ith_bar] <- .fcr_hier_theta_from_vc(bar, vc)
  }
  p_hat[ith_sm] <- abs(p_hat[ith_sm])
  names(p_hat) <- S$int_names
  ll_hat <- S$loglik(p_hat)
  H <- .fcr_final_hessian(S$loglik, p_hat)
  nll <- function(p) {
    v <- -S$loglik(p)
    if (!is.finite(v)) 1e100 else v
  }
  grad_hat <- -.fcr_hier_ngrad(nll)(p_hat)

  # Internal scale -> (beta, distribution parameters, sd and cor)
  to_report <- function(p) {
    c(p[seq_len(k_fix)],
      if (S$n_bar) .fcr_hier_vc_vector(bar, p[ith_bar]))
  }
  rep0 <- to_report(p_hat)
  V_int <- if (!is.null(H)) {
    tryCatch(.fcr_hessian_vcov(list(hessian = H)), error = function(e) NULL)
  }
  if (is.null(V_int)) V_int <- matrix(NA_real_, length(p_hat), length(p_hat))
  Jac <- matrix(0, length(rep0), length(p_hat))
  for (i in seq_along(p_hat)) {
    h <- 1e-6 * max(1, abs(p_hat[i]))
    e <- numeric(length(p_hat))
    e[i] <- h
    Jac[, i] <- (to_report(p_hat + e) - to_report(p_hat - e)) / (2 * h)
  }
  V0 <- Jac %*% V_int %*% t(Jac)
  vc_names <- if (S$n_bar) names(.fcr_hier_vc_vector(bar, p_hat[ith_bar]))
  up <- S$unpack(p_hat)

  # Effects at the estimates
  S$state$u <- rep(0, S$tmpl$Q)
  fin <- NULL
  condvar <- NULL
  if (S$engine == "msl") {
    post <- S$msl_at(p_hat, posterior = TRUE)
    t1 <- bar[[1]]
    u_b <- as.vector(t(post$b))   # group-major
    b <- u_b
    u <- NULL
    eta <- up$eta0 + rowSums(t1$X * post$b[as.integer(t1$g), , drop = FALSE])
    condvar <- list(post$condvar)
  } else {
    fin <- S$mode_at(p_hat)
    u <- fin$md$u
    bb <- as.vector(fin$Lam %*% u)
    Qbar <- sum(vapply(bar, function(t) t$J * t$q, 0))
    b <- bb[seq_len(Qbar)]
    eta <- fin$md$eta
    condvar <- tryCatch(.fcr_hier_condvar(fin, bar), error = function(e)
      NULL)
  }

  # Fixed coefficients on the original columns (smooths: original basis)
  kx <- ncol(S$prep$X)
  beta_full <- numeric(kx)
  names(beta_full) <- colnames(S$prep$X)
  np <- length(S$param_cols)
  beta_full[S$param_cols] <- p_hat[seq_len(np)]
  smooth_out <- NULL
  if (!is.null(S$smooth)) {
    sm <- .fcr_hier_smooth_cov(S, fin, p_hat)
    beta_full <- sm$beta
    V_fix <- sm$V
  }
  est <- c(beta_full, p_hat[S$k_beta + seq_len(S$n_alpha + S$n_sigma)],
           rep0[-seq_len(k_fix)])
  names(est) <- c(colnames(S$prep$X), S$dist_names, vc_names)
  # Covariance on the reported scale
  if (is.null(S$smooth)) {
    V <- V0
  } else {
    nr <- length(est)
    V <- matrix(0, nr, nr)
    other_rep <- kx + seq_len(nr - kx)
    other_0 <- S$k_beta + seq_len(nr - kx)
    V[seq_len(kx), seq_len(kx)] <- V_fix
    V[other_rep, other_rep] <- V0[other_0, other_0]
    V[S$param_cols, other_rep] <- V0[seq_len(np), other_0]
    V[other_rep, S$param_cols] <- V0[other_0, seq_len(np)]
  }
  dimnames(V) <- list(names(est), names(est))

  theta_diag <- unlist(lapply(bar, function(t) {
    L <- .fcr_hier_L(p_hat[k_fix + t$toff + seq_len(t$ntheta)], t$q)
    abs(diag(L))
  }))
  singular <- any(theta_diag < 1e-4)
  # At a zero variance or a correlation of +-1 the Wald standard errors of
  # that term's variance components are not meaningful
  pos <- kx + S$n_alpha + S$n_sigma
  for (t in bar) {
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
  mu <- .fcr_family_mean(S$family, exp(eta), up$alpha, up$sigma)

  conv <- opt$convergence == 0L && all(is.finite(grad_hat)) &&
    max(abs(grad_hat)) < 1e-2
  if (!conv) {
    warning("The optimizer may not have converged (code ", opt$convergence,
            ", largest gradient ", signif(max(abs(grad_hat)), 3), ").",
            call. = FALSE)
  }

  terms_store <- lapply(bar, function(t) {
    t$X <- NULL
    t$g <- NULL
    t
  })
  prep <- S$prep
  fit <- list(
    estimate = est, maximum = ll_hat, vcov = V,
    vcov_type = if (is.null(S$smooth)) "hessian" else "bayesian",
    stderr = "normal", se = .fcr_se(V), hessian_internal = H,
    par_internal = p_hat, gradient = grad_hat,
    iterations = unname(opt$counts[1]), counts = opt$counts,
    message = if (opt$convergence == 0L) "converged" else
      paste("optim() convergence code", opt$convergence),
    numerical = .fcr_hier_numerical(H, grad_hat),
    code = opt$convergence, converged = conv,
    family = S$family, model = S$family, params = S$params,
    fixed_formula = prep$formula,
    dis_param_formula_1 = prep$aux_formulas$dis1,
    dis_param_formula_2 = prep$aux_formulas$dis2,
    offset = prep$blueprint$offset_label, blueprint = prep$blueprint,
    re_terms = terms_store, k_beta = length(S$param_cols),
    k_beta_int = S$k_beta, n_alpha = S$n_alpha, n_sigma = S$n_sigma,
    k_fix = kx + S$n_alpha + S$n_sigma, k_fix_int = k_fix,
    n_theta = S$n_theta_bar, beta_full = beta_full,
    engine = S$engine, nAGQ = S$nAGQ, re_dist = S$re_dist,
    u = u, b = b, condvar = condvar, singular = singular,
    n_groups = stats::setNames(vapply(bar, function(t) as.integer(t$J), 0L),
                               vapply(bar, `[[`, "", "group")),
    y = S$y, observed = S$y, nobs = S$n, n_bic = S$n,
    weights = prep$weights, rows_used = prep$rows_used,
    n_dropped = prep$n_dropped, df_total = length(p_hat),
    eta = eta, predictions = mu, fitted.values = mu, residuals = S$y - mu,
    alpha = up$alpha, sigma = up$sigma, beta_pred = beta_full,
    n_beta = kx, modelType = "countreg_hier")
  if (!is.null(S$smooth)) {
    fit$smooth_table <- sm$table
    fit$smooth_info <- sm$info
    fit$smooth_lambda <- sm$lambda
    fit$edf <- sm$edf
  }
  fit
}

# Numerical checks of the optimum (for the app and diagnostics): the
# gradient and the eigenvalues of the information on the internal scale.
.fcr_hier_numerical <- function(H, grad) {
  ev <- if (is.matrix(H) && all(is.finite(H)))
    eigen(-(H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  list(max_abs_gradient = if (length(grad)) max(abs(grad)) else NA_real_,
       min_eigenvalue = if (length(ev)) min(ev) else NA_real_,
       information = if (!length(ev)) "unavailable" else
         if (min(ev) > sqrt(.Machine$double.eps) * max(abs(ev)))
           "positive definite" else "not positive definite",
       condition = if (length(ev) && min(ev) > 0) max(ev) / min(ev) else
         NA_real_)
}

# Smooth coefficients on their original basis, their Bayesian covariance
# and EDF, from the joint precision of (beta, u) at the mode,
#   C = [X A]' D [X A] + blockdiag(0, I),
# whose inverse is the covariance of the fixed coefficients and the
# effects given the variance parameters (as for lme4's fixed effects).
.fcr_hier_smooth_cov <- function(S, fin, p_hat) {
  kb <- S$k_beta
  np <- length(S$param_cols)
  XA <- cbind(Matrix::Matrix(S$X, sparse = TRUE), fin$A)
  D <- fin$md$D
  Q <- ncol(fin$A)
  C <- Matrix::crossprod(XA, XA * D) +
    Matrix::bdiag(Matrix::Matrix(0, kb, kb, sparse = TRUE),
                  Matrix::Diagonal(Q))
  sm_terms <- S$terms[setdiff(seq_along(S$terms), seq_len(S$n_bar))]
  u_idx <- unlist(lapply(sm_terms, function(t) t$uoff + seq_len(t$J)))
  need <- c(seq_len(kb), kb + u_idx)
  E <- Matrix::sparseMatrix(i = need, j = seq_along(need), x = 1,
                            dims = c(kb + Q, length(need)))
  ch <- Matrix::Cholesky(Matrix::forceSymmetric(C), perm = TRUE, LDL = FALSE)
  Ci <- as.matrix(Matrix::solve(ch, E))[need, , drop = FALSE]
  Ci <- (Ci + t(Ci)) / 2
  # Linear map from (beta_int, u_smooth) to the original fixed columns
  kx <- ncol(S$prep$X)
  Tm <- matrix(0, kx, length(need))
  Tm[cbind(S$param_cols, seq_len(np))] <- 1
  theta <- p_hat[S$k_fix + seq_len(S$n_theta)]
  fpos <- np
  edf <- numeric(length(S$smooth$info))
  lambda <- numeric(0)
  gterms <- list()
  for (si in seq_along(S$smooth$info)) {
    inf <- S$smooth$info[[si]]
    # transformed coefficients: [rand blocks..., fixed]
    cols_t <- integer(0)
    scale <- numeric(0)
    for (ti in inf$terms) {
      t <- S$terms[[S$n_bar + ti]]
      pos_u <- match(t$uoff + seq_len(t$J), u_idx)
      cols_t <- c(cols_t, kb + pos_u)
      scale <- c(scale, rep(theta[t$toff + 1L], t$J))
      lambda <- c(lambda, 1 / max(theta[t$toff + 1L]^2, 1e-300))
      gterms[[length(gterms) + 1L]] <- list(
        component = "main", smooth = inf$label, label = "mean",
        idx = inf$cols)
      edf[si] <- edf[si] + sum(1 - diag(Ci)[kb + pos_u])
    }
    fcols <- fpos + seq_len(inf$n_fixed)
    fpos <- fpos + inf$n_fixed
    cols_t <- c(cols_t, fcols)
    scale <- c(scale, rep(1, inf$n_fixed))
    edf[si] <- edf[si] + inf$n_fixed
    M <- inf$U %*% diag(inf$D * scale, length(scale))
    Tm[inf$cols, cols_t] <- M
  }
  v <- c(p_hat[seq_len(kb)], fin$md$u[u_idx])
  beta <- drop(Tm %*% v)
  names(beta) <- colnames(S$prep$X)
  V <- Tm %*% Ci %*% t(Tm)
  dimnames(V) <- list(names(beta), names(beta))
  # EDF per coefficient for the smooth table (each smooth's total on its
  # first coefficient)
  nr <- kx + S$n_alpha + S$n_sigma
  edf_vec <- numeric(nr)
  for (si in seq_along(S$smooth$info)) {
    edf_vec[S$smooth$info[[si]]$cols[1]] <- edf[si]
  }
  gam <- list(terms = gterms, edf = edf_vec, theta = beta, Vb = V,
              lambda = lambda)
  so <- .fcr_smooth_outputs(gam, list(main = S$prep$smooth$main),
                            list(main = S$prep$X))
  list(beta = beta, V = V, edf = edf_vec, table = so$table, info = so$info,
       lambda = lambda)
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
