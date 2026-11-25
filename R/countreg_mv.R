# Multivariate count models: countreg_mv().
#
# J count outcomes observed on the same units (crashes by severity or type
# at a site). Each outcome has its own family, mean formula and
# distribution-parameter formulas. Stage 1: the equations are independent
# given the covariates, so the log-likelihood is the sum of the equations'
# log-likelihoods, and the fit equals separate countreg() fits; the
# dependent forms (correlated random parameters, copulas) build on the same
# equation engine.
#
# Parameters, in order: for each equation j, its mean coefficients
# ("<outcome>:<term>") and distribution parameters ("<outcome>:ln(alpha)",
# or "<outcome>:ln(alpha):<term>" with a formula).

#' Multivariate count regression
#'
#' Fits several count outcomes observed on the same units, such as crashes
#' by severity or by type at each site, each with its own count family,
#' mean formula and distribution-parameter formulas.
#'
#' @param formula either a formula with the outcomes on the left as
#'   `cbind(y_1, ..., y_J)` and one right-hand side for every equation, or a
#'   list of two-sided formulas, one per outcome (named by the list names,
#'   or by their responses). An `offset()` term (an exposure) belongs in the
#'   formulas.
#' @param data a data frame.
#' @param family the count family of each outcome (see [count_families]):
#'   one code for every outcome, or one per outcome (named by outcome or in
#'   their order). Families made with [new_count_family()] are accepted.
#' @param dis_param_formula_1,dis_param_formula_2 formulas for the
#'   outcomes' distribution parameters, as in [countreg()]: `NULL`
#'   (constant), one one-sided formula for every outcome that has the
#'   parameter, or a list named by outcome (outcomes not in the list keep a
#'   constant parameter).
#' @param offset optional names of offset columns (log scale) used by every
#'   equation, instead of `offset()` in the formulas.
#' @param weights optional frequency weights: a column name or a numeric
#'   vector.
#' @param group_var optional column naming the units for cluster-robust
#'   standard errors (such as a site observed over several years).
#' @param stderr `"normal"` (from the Hessian, the default), `"robust"`
#'   (sandwich) or `"cluster"` (sandwich by `group_var`).
#' @param method,max.iters,start.vals,verbose optimizer settings, as in
#'   [countreg()]; the default method is `"BFGS"`.
#' @param ndraws,engine integration settings of families without a closed
#'   form, as in [countreg()].
#' @param rpar_formula optional one-sided formula of random parameters:
#'   one formula for every outcome, or a list named by outcome. Each term
#'   adds a normal random deviation, with mean zero, to the coefficient of
#'   that variable in the outcome's equation (a variable not in the
#'   equation's formula gets a random coefficient with mean zero). The
#'   intercept is random only when the formula has the term `1` itself
#'   (`~ 1`, `~ 1 + x`). The random parameters are drawn once per unit of
#'   `group_var` (shared by all its rows, such as the years of a site), or
#'   once per row without `group_var`.
#' @param correlation correlation of the random parameters: `"across"`
#'   (default; correlated within and across equations, so the outcomes of a
#'   unit share unobserved factors), `"within"` (correlated within each
#'   equation only) or `"none"`.
#' @param rpar_draws the number of simulation draws per unit (randomly
#'   shifted Halton points, a separate block and shift for each unit).
#' @param seed seed of the random shifts of the draws (stored in the fit,
#'   so predictions use the same draws); random by default.
#' @param copula `"none"` (default), or a copula that joins the outcomes'
#'   distributions: `"gaussian"` (two outcomes), `"frank"`, `"clayton"` or
#'   `"gumbel"`. See Details.
#' @param .setup_only internal use: return the model setup without fitting.
#'
#' @details
#' **Independent equations.** Without `rpar_formula`, the outcomes are
#' independent given the covariates, so the log-likelihood is the sum of the equations'
#' log-likelihoods and every equation's estimates equal those of a
#' separate [countreg()] fit. Fitting them together gives one object with
#' joint predictions, simulations and tests across equations (for
#' example, a Wald test that a coefficient is the same in two equations).
#' The rows used are those complete for every equation.
#'
#' **Correlated outcomes.** Counts of different categories at a site are
#' usually correlated: unobserved conditions that raise one also raise the
#' others. With `rpar_formula`, each equation's log mean gets random terms
#' \eqn{x_r' b_{uj}}, with \eqn{b_u = L z_u} drawn per unit and
#' correlated across equations (`correlation = "across"`).
#'
#' **The multivariate Poisson-lognormal model.** Poisson outcomes with
#' correlated normal random intercepts (`family = "POISSON"`,
#' `rpar_formula = ~ 1`, `correlation = "across"`) give the multivariate
#' Poisson-lognormal model (MVPLN) of, for example, Ma et al. (2008).
#' Other families keep their own overdispersion within each outcome, with
#' the random intercepts adding the correlation between outcomes, and every
#' equation can have a different family.
#' A positive correlation of the random intercepts means that a site with
#' more crashes of one category than its covariates predict also tends to
#' have more of the other. The correlation can be negative too, unlike in
#' models with a shared gamma effect such as [negmultinom()]. The
#' likelihood is simulated (`rpar_draws` draws per unit) with analytic
#' gradients; the elements of \eqn{L} are reported as `sd:` and `chol:`,
#' and `summary()` shows the implied standard deviations and correlations.
#' `predict()` gives joint and conditional probabilities of the outcomes,
#' and Empirical Bayes predictions for units in the estimation data.
#'
#' **Copulas.** A copula joins the outcomes' own distributions, each with
#' its family, into one joint distribution:
#' \eqn{P(Y_1 \le y_1, \dots, Y_J \le y_J) = C(F_1(y_1), \dots,
#' F_J(y_J); \theta)}. For counts, the probability of the observed
#' combination is the exact volume of the copula over the rectangle
#' \eqn{[F_j(y_j - 1), F_j(y_j)]} (the sum over its \eqn{2^J} corners;
#' Nikoloulopoulos and Karlis 2009), not a continuous approximation such as
#' the midpoint or jittering. The margins are those of the equations, so
#' expected counts and each outcome's probabilities do not change; the
#' copula changes the joint and conditional probabilities. The families
#' differ in the shape of the dependence:
#' * `"gaussian"`: symmetric, without tail dependence, any sign
#'   (`copula:atanh(rho)`); two outcomes.
#' * `"frank"`: symmetric, without tail dependence; any sign for two
#'   outcomes (`copula:theta`), positive for more (`copula:ln(theta)`).
#' * `"clayton"`: positive, stronger in the lower tail (joint low counts;
#'   `copula:ln(theta)`).
#' * `"gumbel"`: positive, stronger in the upper tail (joint high counts;
#'   `copula:ln(theta - 1)`).
#'
#' All parameters are estimated together by maximum likelihood, started
#' from the separate equations and the copula parameter that is best for
#' them (inference functions for margins); the gradient is numerical.
#' `summary()` reports the copula parameter and its Kendall's tau. With
#' counts, the copula is not unique outside the range of the margins
#' (Genest and Neslehova 2007): its parameter describes this model's
#' dependence, and the tau is that of the copula, not the rank correlation
#' of the counts, which is weaker for small counts. Copulas cannot yet be
#' combined with `rpar_formula`.
#'
#' Different outcomes often need different families: frequent property
#' damage crashes may be close to NB2, while rare fatal crashes are more
#' skewed and zero-heavy. Each equation can use any family of the package.
#'
#' @return An object of class `"countreg_mv"`, with methods `coef()`,
#'   `vcov()`, `logLik()`, `summary()`, `predict()`, `fitted()`,
#'   `residuals()`, `simulate()`, `tidy()`, `glance()` and `confint()`.
#'
#' @references
#' Genest, C., and Neslehova, J. (2007). A primer on copulas for count
#' data. *ASTIN Bulletin*, 37(2), 475-515.
#'
#' Ma, J., Kockelman, K. M., and Damien, P. (2008). A multivariate
#' Poisson-lognormal regression model for prediction of crash counts by
#' severity, using Bayesian methods. *Accident Analysis and Prevention*,
#' 40(3), 964-975.
#'
#' Nikoloulopoulos, A. K., and Karlis, D. (2009). Modeling multivariate
#' count data using copulas. *Communications in Statistics - Simulation and
#' Computation*, 39(1), 172-187.
#'
#' @seealso [predict.countreg_mv()], [countreg()], [countreg_alloc()] for a
#'   total and its split across categories, [negmultinom()]
#' @examples
#' \donttest{
#' data("washington_roads")
#' w <- washington_roads
#' w$pdo <- w$Total_crashes - w$Injury_crashes
#' m <- countreg_mv(cbind(pdo, Injury_crashes) ~ lnaadt + speed50 +
#'                    offset(lnlength), data = w,
#'                  family = c(pdo = "NB2", Injury_crashes = "PIG"))
#' summary(m)
#' head(predict(m))
#' }
#' @export
countreg_mv <- function(formula, data, family = "NB2",
                        dis_param_formula_1 = NULL,
                        dis_param_formula_2 = NULL, offset = NULL,
                        weights = NULL, group_var = NULL, stderr = "normal",
                        method = "BFGS", max.iters = 1000, start.vals = NULL,
                        ndraws = 1500, engine = "quadrature", verbose = FALSE,
                        rpar_formula = NULL,
                        correlation = c("across", "within", "none"),
                        rpar_draws = 500, seed = NULL,
                        copula = c("none", "gaussian", "frank", "clayton",
                                   "gumbel"),
                        .setup_only = FALSE) {
  cl <- match.call()
  correlation <- match.arg(correlation)
  copula <- match.arg(copula)
  .fcr_check_controls(ndraws = rpar_draws)
  data <- as.data.frame(data)
  method <- .fcr_check_method(method)
  .fcr_check_controls(max.iters, ndraws)
  engine <- match.arg(engine, c("quadrature", "halton", "poilog"))
  stderr <- match.arg(stderr, c("normal", "robust", "cluster"))
  if (stderr == "cluster" && is.null(group_var)) {
    stop("`stderr = \"cluster\"` needs `group_var`.", call. = FALSE)
  }

  # ---- Equations ------------------------------------------------------------
  forms <- .fcr_mv_formulas(formula, data)
  outcomes <- names(forms)
  J <- length(forms)
  fams <- .fcr_mv_per_outcome(family, outcomes, "family")
  family_specs <- lapply(fams, function(f)
    if (inherits(f, "fcr_count_family")) f)
  fams <- vapply(fams, .fcr_resolve_family, character(1))
  d1 <- .fcr_mv_dis_formulas(dis_param_formula_1, outcomes,
                             "dis_param_formula_1")
  d2 <- .fcr_mv_dis_formulas(dis_param_formula_2, outcomes,
                             "dis_param_formula_2")
  for (j in seq_len(J)) {
    params <- get_params(fams[j])
    # One formula given for every outcome only applies where the family has
    # the parameter (a formula named for an outcome must fit its family)
    if (is.null(params[[1]]) && inherits(dis_param_formula_1, "formula")) {
      d1[j] <- list(NULL)
    }
    if (is.null(params[[2]]) && inherits(dis_param_formula_2, "formula")) {
      d2[j] <- list(NULL)
    }
    .fcr_check_dis_formulas(fams[j], params, d1[[j]], d2[[j]])
  }
  rpf <- .fcr_mv_dis_formulas(rpar_formula, outcomes, "rpar_formula")
  for (f in Filter(Negate(is.null), rpf)) {
    if (.fcr_has_formula_offset(f)) {
      stop("`rpar_formula` cannot have an offset.", call. = FALSE)
    }
    if (.fcr_has_smooth(f)) {
      stop("Smooth terms are not available in `rpar_formula`.",
           call. = FALSE)
    }
  }
  if (engine == "poilog" && "PLN" %in% fams) {
    rlang::check_installed("poilog", reason = "for `engine = \"poilog\"`.")
  }

  # ---- Data: the rows complete for every equation ----------------------------
  prep_eq <- function(j, dat, w) {
    prepare_model_data(forms[[j]], dat,
                       aux_formulas = list(dis1 = d1[[j]], dis2 = d2[[j]],
                                           rpar = rpf[[j]]),
                       offset = offset, weights = w, group_var = group_var)
  }
  rows <- seq_len(nrow(data))
  for (j in seq_len(J)) {
    p1 <- prep_eq(j, data, weights)
    rows <- intersect(rows, p1$rows_used)
  }
  if (!length(rows)) stop("No row is complete for every equation.",
                          call. = FALSE)
  n_dropped <- nrow(data) - length(rows)
  if (n_dropped > 0) {
    message(n_dropped, " row(s) with missing values were removed.")
  }
  dsub <- data[rows, , drop = FALSE]
  wsub <- if (is.numeric(weights)) weights[rows] else weights
  preps <- lapply(seq_len(J), function(j) prep_eq(j, dsub, wsub))
  if (any(vapply(preps, .fcr_prep_has_smooth, logical(1)))) {
    stop("Smooth terms are not yet available in countreg_mv().",
         call. = FALSE)
  }
  w <- preps[[1]]$weights

  # ---- Parameter layout --------------------------------------------------------
  haltons <- randtoolbox::halton(ndraws)
  dens_args <- list(haltons = haltons, normed_haltons = stats::qnorm(haltons),
                    engine = engine)
  eqs <- vector("list", J)
  names_all <- character(0)
  for (j in seq_len(J)) {
    pr <- preps[[j]]
    S <- .fcr_alloc_layout(fams[j], get_params(fams[j]), pr$X, pr$aux$dis1,
                           pr$aux$dis2, matrix(0, nrow(pr$X), 0), NULL, 2L,
                           c("a", "b"), FALSE)
    off <- length(names_all)
    nm <- paste0(outcomes[j], ":", S$names[seq_len(S$n_count)])
    eqs[[j]] <- list(
      outcome = outcomes[j], family = fams[j], y = as.numeric(pr$y),
      X = pr$X, off = pr$offset, W1 = pr$aux$dis1, W2 = pr$aux$dis2,
      S = list(idx_beta = off + S$idx_beta, idx_a = off + S$idx_a,
               idx_s = off + S$idx_s, n1 = S$n1, n2 = S$n2, kx = S$kx),
      dens_args = dens_args,
      Xr = .fcr_rpar_matrix(pr$aux$rpar, rpf[[j]]))
    names_all <- c(names_all, nm)
  }
  E <- list(eqs = eqs, J = J, n = length(rows), w = w, names = names_all,
            outcomes = outcomes)

  # ---- Random parameters -----------------------------------------------------------
  has_rpar <- any(!vapply(rpf, is.null, logical(1)))
  has_cop <- copula != "none"
  if (has_cop && has_rpar) {
    stop("Copulas and random parameters cannot yet be combined in ",
         "countreg_mv(); use one or the other.", call. = FALSE)
  }
  if (has_cop) {
    if (copula == "gaussian" && J != 2L) {
      stop("The Gaussian copula is available for two outcomes; use ",
           "\"frank\", \"clayton\" or \"gumbel\" for more.", call. = FALSE)
    }
    E$cop <- list(family = copula, idx = length(names_all) + 1L)
    names_all <- c(names_all, .fcr_cop_par_name(copula, J))
    E$names <- names_all
  }
  draw_seed <- NULL
  if (has_rpar) {
    RS <- .fcr_mv_rlayout(eqs, correlation, length(names_all))
    if (is.null(RS)) {
      stop("`rpar_formula` has no terms. For a random intercept write ~ 1 ",
           "(or ~ 1 + x).", call. = FALSE)
    }
    names_all <- c(names_all, RS$names)
    unit <- if (is.null(group_var)) factor(seq_len(E$n)) else preps[[1]]$group
    draw_seed <- seed %||% sample.int(.Machine$integer.max, 1L)
    E[c("RS", "g", "G", "units", "names")] <- list(
      RS, as.integer(unit), nlevels(unit), levels(unit), names_all)
    E$UA <- .fcr_alloc_draws(E$G, RS$K, rpar_draws, draw_seed)
  }
  if (.setup_only) return(E)

  # ---- Fit ---------------------------------------------------------------------
  if (has_cop) {
    start <- if (is.null(start.vals)) .fcr_mv_cop_start(E) else
      .fcr_check_start(start.vals, names_all)
    fn <- function(p) .fcr_mv_cop_ll(p, E)
    # No analytic scores of the distribution functions: numerical gradient
    gr <- NULL
  } else {
    start <- if (is.null(start.vals)) .fcr_mv_start(E) else
      .fcr_check_start(start.vals, names_all)
    llf <- if (has_rpar) .fcr_mv_sim_ll else .fcr_mv_ll
    fn <- function(p) llf(p, E, grad = FALSE)$ll
    gr <- function(p) llf(p, E, grad = TRUE)$G
  }
  names(start) <- names_all
  fit <- maxLik::maxLik(fn, grad = gr, start = start, method = method,
                        finalHessian = method == "BHHH",
                        control = list(iterlim = max.iters,
                                       printLevel = if (verbose) 2 else 0))
  # A few Newton steps (numerical Hessian of the analytic gradient) tighten
  # the solution after a quasi-Newton optimizer; kept only if the
  # log-likelihood does not decrease
  if (method != "NR" && !is.null(gr)) {
    f2 <- tryCatch(maxLik::maxLik(fn, grad = gr, start = fit$estimate,
                                  method = "NR",
                                  control = list(iterlim = 20)),
                   error = function(e) NULL)
    if (!is.null(f2) && is.finite(f2$maximum) &&
        f2$maximum >= fit$maximum - 1e-8 &&
        isTRUE(.fcr_converged(f2, "NR"))) {
      fit <- f2
      fit$refined_from <- method
    }
  }
  fit <- .fcr_set_final_hessian(fit, fn, gr, method)
  est <- fit$estimate
  names(est) <- names_all
  if (has_rpar) {
    sg <- .fcr_chol_signs(est, E$RS)
    est <- est * sg
    fit$hessian <- fit$hessian * outer(sg, sg)
  }
  conv <- isTRUE(.fcr_converged(fit, if (is.null(fit$refined_from)) method
                                 else "NR"))
  if (!conv) {
    warning("The optimizer did not report convergence: ", fit$message,
            call. = FALSE)
  }
  U <- if (is.null(gr)) maxLik::numericGradient(fn, est) else gr(est)
  colnames(U) <- names_all
  Hv <- .fcr_hessian_vcov(list(hessian = fit$hessian, estimate = est))
  V <- switch(stderr,
    normal = Hv,
    robust = if (has_rpar) {
      # The units (groups) are the independent observations
      Hv %*% crossprod(U) %*% Hv
    } else {
      keep <- w > 0
      Hv %*% crossprod(U[keep, , drop = FALSE] / sqrt(w[keep])) %*% Hv
    },
    cluster = {
      Ug <- if (has_rpar) U else rowsum(U, preps[[1]]$group)
      G <- nrow(Ug)
      Hv %*% (crossprod(Ug) * G / max(G - 1, 1)) %*% Hv
    })
  dimnames(V) <- list(names_all, names_all)

  out <- list(
    call = cl, coefficients = est, vcov = V, loglik = sum(fn(est)),
    gradientObs = U, hessian = fit$hessian, converged = conv,
    iterations = unname(fit$iterations %||% NA), message = fit$message,
    outcomes = outcomes, family = stats::setNames(fams, outcomes),
    family_spec = family_specs, formulas = lapply(preps, `[[`, "formula"),
    dis_param_formula_1 = lapply(preps, function(p) p$aux_formulas$dis1),
    dis_param_formula_2 = lapply(preps, function(p) p$aux_formulas$dis2),
    layout = lapply(eqs, `[[`, "S"),
    blueprints = lapply(preps, `[[`, "blueprint"), stderr = stderr,
    group_var = group_var, y = vapply(eqs, `[[`, numeric(length(rows)), "y"),
    weights = w, nobs = sum(w > 0), n_obs_ic = sum(w),
    rows_used = rows, n_dropped = n_dropped,
    data_used = dsub, engine = engine, ndraws = ndraws,
    n_units = if (has_rpar) E$G else sum(w > 0),
    copula = if (has_cop) list(
      family = copula, eta = unname(est[E$cop$idx]),
      theta = .fcr_cop_theta(est[E$cop$idx], copula, J),
      idx = E$cop$idx),
    rpar = if (has_rpar) list(
      layout = E$RS, formulas = rpf, correlation = correlation,
      draws = rpar_draws, seed = draw_seed, units = E$units,
      unit_of_row = E$g, posterior = .fcr_mv_sim_ll(est, E, grad = TRUE)$W))
  colnames(out$y) <- outcomes
  names(out$formulas) <- outcomes
  class(out) <- "countreg_mv"
  out$fitted.values <- predict(out, newdata = dsub)
  out
}

# The equations' formulas, named by outcome: from cbind(...) ~ rhs or a list
# of two-sided formulas.
.fcr_mv_formulas <- function(formula, data) {
  if (is.list(formula) && !inherits(formula, "formula")) {
    forms <- lapply(formula, stats::as.formula)
    if (length(forms) < 2L) {
      stop("`formula` needs two or more outcomes.", call. = FALSE)
    }
    if (any(vapply(forms, length, 1L) != 3L)) {
      stop("Every formula in the list must have an outcome on the left.",
           call. = FALSE)
    }
    nm <- names(forms)
    resp <- vapply(forms, function(f) paste(deparse(f[[2L]]), collapse = ""),
                   character(1))
    if (is.null(nm)) nm <- resp
    nm[!nzchar(nm)] <- resp[!nzchar(nm)]
  } else {
    formula <- stats::as.formula(formula)
    lhs <- if (length(formula) == 3L) formula[[2L]]
    if (is.null(lhs) || !is.call(lhs) || !identical(lhs[[1L]], as.name("cbind")) ||
        length(lhs) < 3L) {
      stop("`formula` needs two or more outcomes, as cbind(y_1, y_2) ~ x, ",
           "or a list of formulas, one per outcome.", call. = FALSE)
    }
    parts <- as.list(lhs)[-1L]
    nm <- names(parts)
    resp <- vapply(parts, function(e) paste(deparse(e), collapse = ""),
                   character(1))
    if (is.null(nm)) nm <- resp
    nm[!nzchar(nm)] <- resp[!nzchar(nm)]
    forms <- lapply(parts, function(e) {
      f <- formula
      f[[2L]] <- e
      f
    })
  }
  if (anyDuplicated(nm)) {
    stop("The outcomes must have different names: ", toString(nm[duplicated(nm)]),
         ".", call. = FALSE)
  }
  stats::setNames(forms, nm)
}

# One value per outcome from a single value or a vector/list named by
# outcome (or in their order).
.fcr_mv_per_outcome <- function(x, outcomes, what) {
  J <- length(outcomes)
  if (inherits(x, "fcr_count_family")) x <- list(x)
  x <- as.list(x)
  if (length(x) == 1L && is.null(names(x))) return(rep(x, J))
  if (!is.null(names(x)) && all(nzchar(names(x)))) {
    bad <- setdiff(names(x), outcomes)
    miss <- setdiff(outcomes, names(x))
    if (length(bad) || length(miss)) {
      stop("Names of `", what, "` must be the outcomes (",
           toString(outcomes), ").",
           if (length(miss)) paste0(" Missing: ", toString(miss), "."),
           if (length(bad)) paste0(" Unknown: ", toString(bad), "."),
           call. = FALSE)
    }
    return(unname(x[outcomes]))
  }
  if (length(x) != J) {
    stop("`", what, "` needs one value, or one per outcome (", J, ").",
         call. = FALSE)
  }
  unname(x)
}

# Distribution-parameter formulas per outcome: NULL, one formula for all, or
# a list named by outcome.
.fcr_mv_dis_formulas <- function(f, outcomes, what) {
  J <- length(outcomes)
  if (is.null(f)) return(vector("list", J))
  if (inherits(f, "formula")) return(rep(list(f), J))
  if (!is.list(f) || is.null(names(f)) || any(!nzchar(names(f)))) {
    stop("`", what, "` must be a formula or a list of formulas named by ",
         "outcome.", call. = FALSE)
  }
  bad <- setdiff(names(f), outcomes)
  if (length(bad)) {
    stop("Unknown outcome(s) in `", what, "`: ", toString(bad), ".",
         call. = FALSE)
  }
  out <- vector("list", J)
  for (j in seq_len(J)) {
    if (!is.null(f[[outcomes[j]]])) out[[j]] <- stats::as.formula(f[[outcomes[j]]])
  }
  out
}

# Row log-likelihoods of equation `eq` (unweighted) and, with grad, the
# derivatives with respect to its log mean and log distribution parameters.
# `add_eta` is added to the log mean (a random term).
.fcr_mv_eq_parts <- function(p, eq, add_eta = NULL, grad = TRUE) {
  S <- eq$S
  n <- length(eq$y)
  fam <- .fcr_get_family(eq$family)
  lp <- as.vector(eq$X %*% p[S$idx_beta]) + eq$off
  if (!is.null(add_eta)) lp <- lp + add_eta
  mu <- exp(lp)
  dp <- .fcr_alloc_dispar(p, S, eq$W1, eq$W2, n)
  ll <- tryCatch(.fcr_quiet_dens(as.vector(do.call(fam$dens, c(
    list(y = eq$y, predicted = mu, alpha = dp$alpha, sigma = dp$sigma),
    eq$dens_args, list(log = TRUE)))), log = TRUE),
    error = function(e) rep(NA_real_, n))
  ll[!.fcr_family_valid(eq$family, mu, dp$alpha, dp$sigma)] <- -Inf
  if (!grad) return(list(ll = ll, mu = mu, alpha = dp$alpha, sigma = dp$sigma))
  sc <- .fcr_family_score(eq$family, eq$y, mu, dp$alpha, dp$sigma,
                          dens_args = eq$dens_args)
  list(ll = ll, mu = mu, alpha = dp$alpha, sigma = dp$sigma, s_eta = sc$eta,
       s_a = sc$alpha, s_s = sc$sigma)
}

# Gradient rows of equation `eq`'s own parameters, from its parts and the
# row weights w.
.fcr_mv_eq_grad <- function(pt, eq, w, G) {
  S <- eq$S
  G[, S$idx_beta] <- (w * pt$s_eta) * eq$X
  if (S$n1) G[, S$idx_a] <- (w * pt$s_a) * (if (is.null(eq$W1)) 1 else eq$W1)
  if (S$n2) G[, S$idx_s] <- (w * pt$s_s) * (if (is.null(eq$W2)) 1 else eq$W2)
  G
}

# Log-likelihood of independent equations: one contribution per row (the
# sum over equations, times the row weight), and the gradient rows.
.fcr_mv_ll <- function(p, E, grad = FALSE) {
  w <- E$w
  ll <- numeric(E$n)
  G <- if (grad) matrix(0, E$n, length(p))
  for (eq in E$eqs) {
    pt <- .fcr_mv_eq_parts(p, eq, grad = grad)
    ll <- ll + pt$ll
    if (grad) G <- .fcr_mv_eq_grad(pt, eq, w, G)
  }
  ll <- w * ll
  ll[w == 0] <- 0
  if (!grad) {
    ll[!is.finite(ll)] <- NA_real_
    return(list(ll = ll))
  }
  G[!is.finite(G)] <- 0
  list(ll = ll, G = G)
}

# Starting values: a Poisson fit of each equation and the family's
# dispersion start.
.fcr_mv_start <- function(E) {
  st <- numeric(length(E$names))
  for (eq in E$eqs) {
    S <- eq$S
    pm <- suppressWarnings(stats::glm.fit(eq$X, eq$y, weights = E$w,
                                          offset = eq$off,
                                          family = stats::poisson()))
    b <- pm$coefficients
    b[is.na(b)] <- 0
    st[S$idx_beta] <- b
    fs <- .fcr_get_family(eq$family)$start
    if (S$n1) {
      a0 <- .fcr_moment_disp_start(eq$family, eq$y, pm$fitted.values, E$w) %||%
        fs[1L]
      pos <- if (is.null(eq$W1)) 1L else match("(Intercept)", colnames(eq$W1))
      if (!is.na(pos)) st[S$idx_a[pos]] <- a0
    }
    if (S$n2) {
      pos <- if (is.null(eq$W2)) 1L else match("(Intercept)", colnames(eq$W2))
      if (!is.na(pos)) st[S$idx_s[pos]] <- rep_len(fs, 2L)[2L]
    }
  }
  if (!is.null(E$RS)) {
    dg <- E$RS$pairs[, "k"] == E$RS$pairs[, "l"]
    st[E$RS$idx] <- ifelse(dg, 0.1, 0)
  }
  st
}

# Starting values of a copula model (inference functions for margins): the
# independent equations' estimates, then the copula parameter that
# maximizes the likelihood with those margins.
.fcr_mv_cop_start <- function(E) {
  k <- E$cop$idx
  Ei <- E
  Ei$cop <- NULL
  st <- .fcr_mv_start(Ei)[-k]
  fi <- tryCatch(maxLik::maxLik(
    function(p) .fcr_mv_ll(p, Ei)$ll,
    grad = function(p) .fcr_mv_ll(p, Ei, grad = TRUE)$G, start = st,
    method = "BFGS", control = list(iterlim = 500)), error = function(e) NULL)
  pm <- if (is.null(fi)) st else fi$estimate
  rng <- switch(E$cop$family, gaussian = c(-3, 3),
                frank = if (E$J == 2L) c(-30, 30) else c(-5, 4),
                clayton = c(-6, 3), gumbel = c(-6, 3))
  f <- function(e) {
    v <- sum(.fcr_mv_cop_ll(c(pm, e), E))
    if (is.finite(v)) v else -1e300
  }
  eta <- stats::optimize(f, rng, maximum = TRUE)$maximum
  c(pm, eta)
}

# ---- Random parameters ------------------------------------------------------

# Layout of the random terms: one per column of each equation's random
# design (its equation and column), and the Cholesky factor's free elements
# (none: diagonal; within: blocks by equation; across: full).
.fcr_mv_rlayout <- function(eqs, correlation, offset) {
  eq <- integer(0)
  col <- integer(0)
  labels <- character(0)
  for (j in seq_along(eqs)) {
    Xr <- eqs[[j]]$Xr
    if (is.null(Xr)) next
    eq <- c(eq, rep(j, ncol(Xr)))
    col <- c(col, seq_len(ncol(Xr)))
    labels <- c(labels, paste0(eqs[[j]]$outcome, ":", colnames(Xr)))
  }
  K <- length(eq)
  if (!K) return(NULL)
  pairs <- do.call(rbind, lapply(seq_len(K), function(k) {
    l <- switch(correlation,
      none = k,
      within = which(eq == eq[k] & seq_len(K) <= k),
      across = seq_len(k))
    cbind(k = k, l = l)
  }))
  nm <- ifelse(pairs[, "k"] == pairs[, "l"],
               paste0("sd:", labels[pairs[, "k"]]),
               paste0("chol:", labels[pairs[, "k"]], "|",
                      labels[pairs[, "l"]]))
  list(K = K, eq = eq, col = col, labels = labels, pairs = pairs,
       names = nm, idx = offset + seq_len(nrow(pairs)))
}

# Each equation's random term in the log mean for random coefficients B
# (rows x K); NULL for equations without random parameters.
.fcr_mv_add <- function(B, E) {
  RS <- E$RS
  lapply(seq_along(E$eqs), function(j) {
    k <- which(RS$eq == j)
    if (!length(k)) return(NULL)
    rowSums(E$eqs[[j]]$Xr[, RS$col[k], drop = FALSE] * B[, k, drop = FALSE])
  })
}

# Simulated log-likelihood: one contribution per unit,
# log((1/R) sum_r prod_(rows of the unit) prod_j f_j(y_ij | b_r)), and its
# gradient through the posterior weights of the draws.
.fcr_mv_sim_ll <- function(p, E, grad = FALSE) {
  RS <- E$RS
  L <- .fcr_alloc_chol(p, RS)
  g <- E$g
  G <- E$G
  R <- dim(E$UA)[3L]
  w <- E$w
  row_ll <- function(adds) {
    ll <- numeric(E$n)
    for (j in seq_along(E$eqs)) {
      ll <- ll + .fcr_mv_eq_parts(p, E$eqs[[j]], adds[[j]], grad = FALSE)$ll
    }
    ll <- w * ll
    ll[w == 0] <- 0
    ll
  }
  LL <- matrix(NA_real_, G, R)
  for (r in seq_len(R)) {
    B <- (matrix(E$UA[, , r], G, RS$K) %*% t(L))[g, , drop = FALSE]
    LL[, r] <- rowsum(row_ll(.fcr_mv_add(B, E)), g, reorder = TRUE)[, 1L]
  }
  LL[is.nan(LL)] <- -Inf
  mx <- apply(LL, 1L, max)
  mx[!is.finite(mx)] <- 0
  sm <- rowSums(exp(LL - mx))
  lg <- mx + log(sm) - log(R)
  if (!grad) {
    lg[!is.finite(lg)] <- NA_real_
    return(list(ll = lg))
  }
  W <- exp(LL - mx) / sm
  W[!is.finite(W)] <- 0
  np <- length(p)
  Gg <- matrix(0, G, np)
  for (r in seq_len(R)) {
    if (!any(W[, r] > 1e-12)) next
    Ug <- matrix(E$UA[, , r], G, RS$K)
    B <- (Ug %*% t(L))[g, , drop = FALSE]
    adds <- .fcr_mv_add(B, E)
    Gr <- matrix(0, E$n, np)
    Db <- matrix(0, E$n, RS$K)
    for (j in seq_along(E$eqs)) {
      eq <- E$eqs[[j]]
      pt <- .fcr_mv_eq_parts(p, eq, adds[[j]], grad = TRUE)
      Gr <- .fcr_mv_eq_grad(pt, eq, w, Gr)
      for (k in which(RS$eq == j)) {
        Db[, k] <- w * pt$s_eta * eq$Xr[, RS$col[k]]
      }
    }
    Urow <- Ug[g, , drop = FALSE]
    Gr[, RS$idx] <- Db[, RS$pairs[, "k"], drop = FALSE] *
      Urow[, RS$pairs[, "l"], drop = FALSE]
    Gr[!is.finite(Gr)] <- 0
    Gg <- Gg + rowsum(Gr, g, reorder = TRUE) * W[, r]
  }
  lg[!is.finite(lg)] <- 0
  list(ll = lg, G = Gg, W = W)
}

# Number of draws for predictions: at least this many.
.fcr_mv_pred_draws <- 500L

# Average of FUN(adds) over the distribution of the random parameters for
# the rows of the design E: marginal (plain Halton points, the same for
# every row) or, for the rows of units in the estimation data, posterior
# (their draws weighted by their posterior weights). FUN returns an n x m
# matrix. Without random parameters, FUN(NULL).
.fcr_mv_integrate <- function(object, E, newdata, FUN,
                              random = "marginal") {
  rp <- object$rpar
  if (is.null(rp)) return(FUN(NULL))
  RS <- rp$layout
  L <- .fcr_alloc_chol(object$coefficients, RS)
  n <- E$n
  Rm <- max(rp$draws, .fcr_mv_pred_draws)
  Um <- matrix(stats::qnorm(.fcr_halton(Rm, RS$K)), Rm, RS$K)
  acc <- 0
  for (r in seq_len(Rm)) {
    B <- matrix(as.vector(Um[r, ] %*% t(L)), n, RS$K, byrow = TRUE)
    acc <- acc + FUN(.fcr_mv_add(B, E)) / Rm
  }
  if (random == "posterior") {
    gv <- object$group_var
    if (is.null(gv) || !all(gv %in% names(newdata))) {
      stop("Posterior predictions need the grouping column(s) ",
           toString(gv %||% "(group_var)"), " in `newdata` (and a fit with ",
           "`group_var`).", call. = FALSE)
    }
    key <- as.character(interaction(newdata[, gv, drop = FALSE],
                                    drop = TRUE, sep = ":"))
    ui <- match(key, rp$units)
    seen <- which(!is.na(ui))
    if (length(seen)) {
      UA <- .fcr_alloc_draws(length(rp$units), RS$K, rp$draws, rp$seed)
      post <- 0
      for (r in seq_len(rp$draws)) {
        B <- matrix(0, n, RS$K)
        B[seen, ] <- matrix(UA[ui[seen], , r], length(seen), RS$K) %*% t(L)
        post <- post + rp$posterior[ui, r] * FUN(.fcr_mv_add(B, E))
      }
      acc[seen, ] <- post[seen, ]
    }
  }
  acc
}

# Row log-probabilities of counts Y (n x J matrix, columns in outcome
# order; NA columns are left out) at random terms `adds`.
.fcr_mv_logp <- function(p, E, Y, adds = NULL, which = seq_along(E$eqs)) {
  out <- numeric(E$n)
  for (j in which) {
    eq <- E$eqs[[j]]
    eq$y <- Y[, j]
    out <- out + .fcr_mv_eq_parts(p, eq, adds[[j]], grad = FALSE)$ll
  }
  out
}

# The equations of a fit for new data.
.fcr_mv_design <- function(object, newdata) {
  haltons <- randtoolbox::halton(object$ndraws)
  dens_args <- list(haltons = haltons, normed_haltons = stats::qnorm(haltons),
                    engine = object$engine)
  eqs <- lapply(seq_along(object$outcomes), function(j) {
    fd <- forge_new_data(object$blueprints[[j]], newdata)
    list(outcome = object$outcomes[j], family = object$family[[j]],
         y = rep(0, nrow(fd$X)), X = fd$X, off = fd$offset,
         W1 = fd$aux$dis1, W2 = fd$aux$dis2, S = object$layout[[j]],
         dens_args = dens_args,
         Xr = .fcr_rpar_matrix(fd$aux$rpar, object$rpar$formulas[[j]]))
  })
  list(eqs = eqs, J = length(eqs), n = nrow(newdata),
       outcomes = object$outcomes, RS = object$rpar$layout)
}

# Probability of the counts Y (n x J; NA = outcome left out) under the
# fitted copula: the C-volume of the rectangle of each row.
.fcr_mv_cop_prob <- function(object, E, Y) {
  bd <- .fcr_mv_cop_bounds(object$coefficients, E, Y)
  if (is.null(bd)) return(rep(NA_real_, nrow(Y)))
  pmax(.fcr_cop_rect(bd$hi, bd$lo, object$copula$family,
                     object$copula$theta), 0)
}

# Probabilities of the counts 0..mc of equation j at random terms `add`
# (n x (mc + 1)).
.fcr_mv_pmf <- function(p, eq, add, mc) {
  P <- vapply(0:mc, function(k) {
    eq$y <- rep(k, length(eq$y))
    exp(.fcr_mv_eq_parts(p, eq, add, grad = FALSE)$ll)
  }, numeric(length(eq$y)))
  matrix(P, nrow = length(eq$y))
}

# Counts of the outcomes as an n x J matrix (NA for outcomes not given)
# from a named vector (the same for every row) or a matrix or data frame
# with columns named by outcome.
.fcr_mv_count_matrix <- function(y, outcomes, n, what) {
  if (is.data.frame(y)) y <- as.matrix(y)
  nm <- if (is.matrix(y)) colnames(y) else names(y)
  if (is.null(nm) || any(!nzchar(nm))) {
    stop("`", what, "` must be named by outcome (", toString(outcomes), ").",
         call. = FALSE)
  }
  bad <- setdiff(nm, outcomes)
  if (length(bad)) {
    stop("Unknown outcome(s) in `", what, "`: ", toString(bad), ".",
         call. = FALSE)
  }
  Y <- matrix(NA_real_, n, length(outcomes), dimnames = list(NULL, outcomes))
  if (is.matrix(y)) {
    if (nrow(y) != n && nrow(y) != 1L) {
      stop("`", what, "` needs one row, or one per row of `newdata`.",
           call. = FALSE)
    }
    for (c in nm) Y[, c] <- rep_len(y[, c], n)
  } else {
    for (c in nm) Y[, c] <- y[[c]]
  }
  if (any(!.fcr_is_count(Y[, nm][!is.na(Y[, nm])]))) {
    stop("`", what, "` must hold non-negative whole numbers.", call. = FALSE)
  }
  Y
}

# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------

#' Predictions from multivariate count models
#'
#' @param object a [countreg_mv()] fit.
#' @param newdata optional new data; the estimation data by default.
#' @param type the prediction:
#'   \describe{
#'   \item{`"response"` (default)}{a matrix of the expected counts, one
#'     column per outcome.}
#'   \item{`"link"`}{the linear predictors (log of the mean parameter, at
#'     the mean of the random parameters).}
#'   \item{`"variance"`}{the variances (including the variance from the
#'     random parameters).}
#'   \item{`"prob"`}{a list with one matrix per outcome of the
#'     probabilities of the counts `0:max_count`.}
#'   \item{`"joint"`}{the joint probability of the counts `y` of all
#'     outcomes, one value per row.}
#'   \item{`"conditional"`}{a matrix of the probabilities of the counts
#'     `0:max_count` of `outcome` given the counts `given` of other
#'     outcomes, such as the probability of a fatal crash at a site given
#'     its injury crashes.}
#'   }
#' @param max_count for `"prob"` and `"conditional"`: the largest count
#'   (default: the largest observed count of the outcome).
#' @param y for `"joint"`: the counts of every outcome, a vector named by
#'   outcome (the same for every row) or a matrix or data frame with one
#'   column per outcome.
#' @param outcome,given for `"conditional"`: the name of the outcome, and
#'   the counts of the outcomes conditioned on, named as `y` (they can be
#'   some of the other outcomes).
#' @param random for fits with random parameters: `"marginal"` (default)
#'   averages over their distribution (the prediction for a new unit);
#'   `"posterior"` averages over each unit's posterior distribution given
#'   its counts in the estimation data (Empirical Bayes; units not in the
#'   estimation data get the marginal prediction). `newdata` then needs the
#'   `group_var` column(s).
#' @param ... unused.
#' @return A matrix with one column per outcome (`"response"`, `"link"`,
#'   `"variance"`), a list of matrices (`"prob"`), a vector (`"joint"`) or
#'   a matrix with one column per count (`"conditional"`).
#' @details Without random parameters the outcomes are independent given
#'   the covariates, so joint probabilities are products of the outcomes'
#'   probabilities and conditional probabilities equal the marginal ones.
#'   With random parameters shared by the equations they are not: the
#'   probabilities are averages over the random parameters of products of
#'   the outcomes' probabilities.
#' @export
predict.countreg_mv <- function(object, newdata = NULL,
                                type = c("response", "link", "variance",
                                         "prob", "joint", "conditional"),
                                max_count = NULL, y = NULL, outcome = NULL,
                                given = NULL,
                                random = c("marginal", "posterior"), ...) {
  type <- match.arg(type)
  random <- match.arg(random)
  if (is.null(newdata)) newdata <- object$data_used
  newdata <- as.data.frame(newdata)
  E <- .fcr_mv_design(object, newdata)
  p <- object$coefficients
  outs <- object$outcomes
  J <- length(outs)
  n <- E$n
  vals <- function(adds, j) .fcr_mv_eq_parts(p, E$eqs[[j]], adds[[j]],
                                             grad = FALSE)
  integrate <- function(FUN) .fcr_mv_integrate(object, E, newdata, FUN,
                                               random)
  named <- function(M) {
    M <- matrix(M, nrow = n)
    colnames(M) <- outs
    M
  }
  if (type == "link") {
    return(named(vapply(seq_len(J), function(j) log(vals(NULL, j)$mu),
                        numeric(n))))
  }
  if (type %in% c("response", "variance")) {
    A <- integrate(function(adds) {
      do.call(cbind, lapply(seq_len(J), function(j) {
        v <- vals(adds, j)
        fam <- E$eqs[[j]]$family
        m <- .fcr_family_mean(fam, v$mu, v$alpha, v$sigma)
        if (type == "response") return(m)
        cbind(m, rep_len(.fcr_family_variance(fam, v$mu, v$alpha, v$sigma),
                         n), m^2)
      }))
    })
    if (type == "response") return(named(A))
    m <- A[, 3L * seq_len(J) - 2L, drop = FALSE]
    return(named(A[, 3L * seq_len(J) - 1L, drop = FALSE] +
                   A[, 3L * seq_len(J), drop = FALSE] - m^2))
  }
  if (type == "prob") {
    out <- lapply(seq_len(J), function(j) {
      mc <- max_count %||% max(object$y[, j])
      P <- integrate(function(adds) .fcr_mv_pmf(p, E$eqs[[j]], adds[[j]], mc))
      colnames(P) <- 0:mc
      P
    })
    return(stats::setNames(out, outs))
  }
  if (type == "joint") {
    if (is.null(y)) stop("`type = \"joint\"` needs `y`.", call. = FALSE)
    Y <- .fcr_mv_count_matrix(y, outs, n, "y")
    if (anyNA(Y)) {
      stop("`y` must give the counts of every outcome (", toString(outs),
           ").", call. = FALSE)
    }
    if (!is.null(object$copula)) return(.fcr_mv_cop_prob(object, E, Y))
    return(as.vector(integrate(function(adds)
      matrix(exp(.fcr_mv_logp(p, E, Y, adds)), ncol = 1L))))
  }
  # conditional
  if (is.null(outcome) || length(outcome) != 1L || !outcome %in% outs) {
    stop("`outcome` must be one of ", toString(outs), ".", call. = FALSE)
  }
  if (is.null(given)) stop("`type = \"conditional\"` needs `given`.",
                           call. = FALSE)
  Yg <- .fcr_mv_count_matrix(given, outs, n, "given")
  t <- match(outcome, outs)
  if (!all(is.na(Yg[, t]))) {
    stop("`given` cannot include the outcome itself.", call. = FALSE)
  }
  gj <- which(colSums(!is.na(Yg)) > 0)
  mc <- max_count %||% max(object$y[, t])
  if (!is.null(object$copula)) {
    den <- .fcr_mv_cop_prob(object, E, Yg)
    P <- vapply(0:mc, function(k) {
      Yk <- Yg
      Yk[, t] <- k
      .fcr_mv_cop_prob(object, E, Yk)
    }, numeric(n)) / den
    P <- matrix(P, nrow = n)
    colnames(P) <- 0:mc
    return(P)
  }
  A <- integrate(function(adds) {
    base <- exp(.fcr_mv_logp(p, E, Yg, adds, which = gj))
    cbind(.fcr_mv_pmf(p, E$eqs[[t]], adds[[t]], mc) * base, base)
  })
  P <- A[, seq_len(mc + 1L), drop = FALSE] / A[, mc + 2L]
  colnames(P) <- 0:mc
  P
}

#' @export
fitted.countreg_mv <- function(object, ...) object$fitted.values

#' @export
residuals.countreg_mv <- function(object, type = c("response", "pearson"),
                                  ...) {
  type <- match.arg(type)
  r <- object$y - object$fitted.values
  if (type == "pearson") {
    r <- r / sqrt(predict(object, newdata = object$data_used,
                          type = "variance"))
  }
  r
}

#' @export
coef.countreg_mv <- function(object, ...) object$coefficients

#' @export
vcov.countreg_mv <- function(object, ...) object$vcov

#' @export
nobs.countreg_mv <- function(object, ...) object$nobs

#' @export
logLik.countreg_mv <- function(object, ...) {
  structure(object$loglik, df = length(object$coefficients),
            nobs = object$n_obs_ic, class = "logLik")
}

#' @export
print.countreg_mv <- function(x, ...) {
  cat("Multivariate count model:", length(x$outcomes), "outcomes\n")
  cat("Families:", paste0(x$outcomes, " = ",
                          vapply(x$family, .fcr_family_label, ""),
                          collapse = ", "), "\n\n")
  print(round(x$coefficients, 4))
  invisible(x)
}

# The outcome (equation) of each coefficient.
.fcr_mv_coef_outcome <- function(object) {
  out <- character(length(object$coefficients))
  for (j in seq_along(object$outcomes)) {
    S <- object$layout[[j]]
    out[c(S$idx_beta, S$idx_a, S$idx_s)] <- object$outcomes[j]
  }
  if (!is.null(object$rpar)) {
    out[object$rpar$layout$idx] <- "Random parameters (Cholesky factor)"
  }
  if (!is.null(object$copula)) out[object$copula$idx] <- "Copula"
  out
}

#' @export
summary.countreg_mv <- function(object, ...) {
  est <- object$coefficients
  se <- sqrt(pmax(diag(object$vcov), 0))
  z <- est / se
  tab <- data.frame(Estimate = est, `Std. Error` = se, `z value` = z,
                    `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)),
                    check.names = FALSE)
  structure(list(table = tab, part = .fcr_mv_coef_outcome(object),
                 object = object), class = "summary.countreg_mv")
}

#' @export
print.summary.countreg_mv <- function(x, digits = 4, ...) {
  o <- x$object
  cat("Multivariate count model\n\n")
  cat("Rows:       ", o$nobs, "\n")
  cat("Std. errors:", o$stderr, "\n")
  rp <- o$rpar
  cop <- o$copula
  if (!is.null(cop)) {
    cat("Dependence:  ", .fcr_cop_label(cop$family), "copula\n\n")
  } else if (is.null(rp)) {
    cat("Equations:   independent given the covariates\n\n")
  } else {
    cat("Units:      ", o$n_units, "\n")
    cat("Random parameters:", rp$correlation, "correlation,", rp$draws,
        "draws per", if (is.null(o$group_var)) "row" else
          paste("unit of", toString(o$group_var)), "\n\n")
  }
  for (j in seq_along(o$outcomes)) {
    oc <- o$outcomes[j]
    cat(oc, " (", .fcr_family_label(o$family[[j]]), "):\n", sep = "")
    tb <- x$table[x$part == oc, , drop = FALSE]
    rownames(tb) <- sub(paste0("^", .fcr_regex_escape(oc), ":"), "",
                        rownames(tb))
    stats::printCoefmat(tb, digits = digits, signif.legend = FALSE)
    cat("\n")
  }
  if (!is.null(rp)) {
    cat("Random parameters (Cholesky factor):\n")
    stats::printCoefmat(x$table[rp$layout$idx, , drop = FALSE],
                        digits = digits, signif.legend = FALSE)
    vc <- .fcr_alloc_varcorr(o)
    sdv <- sqrt(diag(vc))
    tab <- cbind(SD = sdv)
    if (length(sdv) > 1L) {
      cr <- stats::cov2cor(vc)
      cr[upper.tri(cr, diag = TRUE)] <- NA
      tab <- cbind(tab, cr[, -ncol(cr), drop = FALSE])
    }
    cat("\nStandard deviations and correlations:\n")
    print(round(tab, 3), na.print = "")
    cat("\n")
  }
  if (!is.null(cop)) {
    cat("Copula:\n")
    stats::printCoefmat(x$table[cop$idx, , drop = FALSE], digits = digits,
                        signif.legend = FALSE)
    dep <- .fcr_mv_cop_dependence(o)
    cat("\n", .fcr_cop_label(cop$family), " copula: ", dep$label, " = ",
        format(dep$theta, digits = 4), ", Kendall's tau = ",
        format(dep$tau, digits = 3), " (SE ", format(dep$tau_se, digits = 2),
        ")\n\n", sep = "")
  }
  k <- length(o$coefficients)
  cat("Log-likelihood:", format(o$loglik, digits = 8), " AIC:",
      format(-2 * o$loglik + 2 * k, digits = 8), " BIC:",
      format(-2 * o$loglik + log(o$n_obs_ic) * k, digits = 8), "\n")
  if (!o$converged) cat("Warning: the optimizer did not converge.\n")
  invisible(x)
}

.fcr_regex_escape <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)

#' @export
tidy.countreg_mv <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  s <- summary(x)
  oc <- s$part
  term <- mapply(function(t, o) sub(paste0("^", .fcr_regex_escape(o), ":"),
                                    "", t), rownames(s$table), oc,
                 USE.NAMES = FALSE)
  comp <- rep("mean", length(term))
  for (j in seq_along(x$outcomes)) {
    S <- x$layout[[j]]
    comp[c(S$idx_a, S$idx_s)] <- "distribution"
  }
  if (!is.null(x$copula)) {
    comp[x$copula$idx] <- "copula"
    oc[x$copula$idx] <- NA_character_
    term[x$copula$idx] <- rownames(s$table)[x$copula$idx]
  }
  if (!is.null(x$rpar)) {
    RS <- x$rpar$layout
    comp[RS$idx] <- "random"
    oc[RS$idx] <- x$outcomes[RS$eq[RS$pairs[, "k"]]]
    term[RS$idx] <- RS$names
  }
  out <- data.frame(outcome = oc, term = term, component = comp,
                    estimate = s$table[, 1], std.error = s$table[, 2],
                    statistic = s$table[, 3], p.value = s$table[, 4],
                    stringsAsFactors = FALSE, row.names = NULL)
  if (conf.int) {
    zq <- stats::qnorm(1 - (1 - conf.level) / 2)
    out$conf.low <- out$estimate - zq * out$std.error
    out$conf.high <- out$estimate + zq * out$std.error
  }
  tibble::as_tibble(out)
}

#' @export
glance.countreg_mv <- function(x, ...) {
  k <- length(x$coefficients)
  tibble::tibble(logLik = x$loglik, AIC = -2 * x$loglik + 2 * k,
                 BIC = -2 * x$loglik + log(x$n_obs_ic) * k, nobs = x$nobs,
                 npar = k, outcomes = length(x$outcomes),
                 converged = x$converged)
}

#' @export
simulate.countreg_mv <- function(object, nsim = 1, seed = NULL,
                                 newdata = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(newdata)) newdata <- object$data_used
  newdata <- as.data.frame(newdata)
  E <- .fcr_mv_design(object, newdata)
  p <- object$coefficients
  n <- E$n
  rp <- object$rpar
  # With random parameters, each simulation draws new coefficients for
  # each unit (the group_var units of newdata, or each row)
  unit <- if (!is.null(rp) && !is.null(object$group_var) &&
              all(object$group_var %in% names(newdata))) {
    as.integer(factor(interaction(newdata[, object$group_var, drop = FALSE],
                                  drop = TRUE, sep = ":")))
  } else seq_len(n)
  cop <- object$copula
  one <- function() {
    if (!is.null(cop)) {
      # Uniforms from the copula, then each outcome's quantiles
      U <- .fcr_rcopula(n, length(E$eqs), cop$family, cop$theta)
      Y <- vapply(seq_along(E$eqs), function(j) {
        v <- .fcr_mv_eq_values(p, E$eqs[[j]])
        .fcr_mv_qcount(E$eqs[[j]], v, U[, j])
      }, numeric(n))
      Y <- matrix(Y, nrow = n)
      colnames(Y) <- object$outcomes
      return(Y)
    }
    adds <- NULL
    if (!is.null(rp)) {
      L <- .fcr_alloc_chol(p, rp$layout)
      Bu <- matrix(stats::rnorm(max(unit) * rp$layout$K), max(unit)) %*% t(L)
      adds <- .fcr_mv_add(Bu[unit, , drop = FALSE], E)
    }
    Y <- vapply(seq_along(E$eqs), function(j) {
      v <- .fcr_mv_eq_parts(p, E$eqs[[j]], adds[[j]], grad = FALSE)
      .fcr_family_r(E$eqs[[j]]$family, n, v$mu, v$alpha, v$sigma)
    }, numeric(n))
    Y <- matrix(Y, nrow = n)
    colnames(Y) <- object$outcomes
    Y
  }
  stats::setNames(replicate(nsim, one(), simplify = FALSE),
                  paste0("sim_", seq_len(nsim)))
}

.fcr_cop_label <- function(family) {
  c(gaussian = "Gaussian", frank = "Frank", clayton = "Clayton",
    gumbel = "Gumbel")[[family]]
}

# The copula parameter on its own scale and its Kendall's tau, with a
# delta-method standard error of tau.
.fcr_mv_cop_dependence <- function(object) {
  cop <- object$copula
  J <- length(object$outcomes)
  eta <- cop$eta
  se_eta <- sqrt(max(object$vcov[cop$idx, cop$idx], 0))
  tau_of <- function(e) .fcr_cop_tau(.fcr_cop_theta(e, cop$family, J),
                                     cop$family)
  h <- 1e-5 * max(1, abs(eta))
  dtau <- (tau_of(eta + h) - tau_of(eta - h)) / (2 * h)
  list(label = if (cop$family == "gaussian") "rho" else "theta",
       theta = cop$theta, tau = tau_of(eta), tau_se = abs(dtau) * se_eta)
}
