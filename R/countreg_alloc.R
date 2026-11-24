# Total-plus-allocation count models: countreg_alloc().
#
# The total count of each unit follows a count family (any family of the
# package, with distribution-parameter formulas), and its split across J
# ordered (or unordered) categories follows an allocation model:
#
#   P(y_1, ..., y_J) = P(c) * c! / prod(y_j!) * prod(pi_j^y_j)   (multinomial)
#
# with c = sum(y_j). The "fractional" likelihood of Shabab et al. (2026)
# raises pi_j to the observed share y_j / c instead; it is a
# quasi-likelihood. Rows with a zero total carry no information about the
# split (all exponents are 0) and contribute P(c = 0) only; the total is not
# truncated unless `truncated = TRUE`.
#
# Allocation links:
#   ordered_probit / ordered_logit: latent propensity z'gamma, thresholds
#     psi_1 = v'theta_1, psi_k = psi_(k-1) + exp(v'theta_k), and
#     pi_j = F(psi_j - z'gamma) - F(psi_(j-1) - z'gamma).
#   mnl: pi_j = exp(z'delta_j) / sum_k exp(z'delta_k), delta_1 = 0.
#
# Parameters, in order: the count model's mean coefficients and its
# distribution parameters (named as in countreg()), the allocation
# coefficients ("Allocation:<col>", or "Allocation[<category>]:<col>" for
# mnl) and the thresholds ("Threshold<k>:<col>", k >= 2 on the log-gap
# scale).

#' Total-plus-allocation count models
#'
#' Fits a joint model of a total count and its allocation across
#' categories, such as the number of crashes at a site and their split by
#' severity. The total follows any count family of the package (with
#' covariates in its mean and distribution parameters), and the split
#' follows an ordered probit, ordered logit or multinomial logit model,
#' with covariates in the propensity and, for the ordered links, in the
#' thresholds (generalized thresholds).
#'
#' @param formula a formula with the category counts on the left, as
#'   `cbind(cat_1, ..., cat_J)` from the lowest to the highest category (for
#'   the ordered links), and the total's mean model on the right. The total
#'   is the sum of the categories. An `offset()` term (an exposure) belongs
#'   here.
#' @param data a data frame.
#' @param family the count family of the total (see [countreg()]), or a
#'   family made with [new_count_family()].
#' @param alloc_formula a one-sided formula for the allocation propensity.
#'   By default the right-hand side of `formula` without its offsets. For
#'   the ordered links its intercept is dropped (the thresholds absorb it).
#' @param threshold_formula a one-sided formula for the thresholds of the
#'   ordered links (generalized thresholds). Default `~ 1`, fixed
#'   thresholds. It must have an intercept.
#' @param link `"ordered_probit"` (default), `"ordered_logit"` or `"mnl"`
#'   (multinomial logit, for unordered categories; the first category is
#'   the base).
#' @param likelihood `"multinomial"` (default) or `"fractional"`; see
#'   Details.
#' @param dis_param_formula_1,dis_param_formula_2 formulas for the total's
#'   distribution parameters, as in [countreg()].
#' @param offset optional names of offset columns (log scale), instead of
#'   `offset()` in `formula`.
#' @param weights optional frequency weights: a column name or a numeric
#'   vector.
#' @param truncated `TRUE` when the data contain only units with a positive
#'   total (for example data built from crash records): the total then
#'   follows the zero-truncated distribution \eqn{P(c) / (1 - P(0))}.
#' @param group_var optional column naming the units for cluster-robust
#'   standard errors (such as a site observed over several years).
#' @param stderr `"normal"` (from the Hessian), `"robust"` (sandwich) or
#'   `"cluster"` (sandwich by `group_var`). The default is `"normal"` for
#'   the multinomial likelihood and `"cluster"` (or `"robust"` without
#'   `group_var`) for the fractional one, whose Hessian alone does not give
#'   valid standard errors.
#' @param method,max.iters,start.vals,verbose optimizer settings, as in
#'   [countreg()].
#' @param ndraws,engine integration settings of families without a closed
#'   form, as in [countreg()].
#' @param rpar_total,rpar_alloc,rpar_threshold one-sided formulas of
#'   random parameters in the total's mean, the allocation propensity and
#'   the gaps between thresholds. Each term adds a normal random deviation,
#'   with mean zero, to the corresponding coefficient, so the coefficient
#'   in the equation's own formula is the mean (`~ 1` in `rpar_alloc` is a
#'   random intercept of the propensity, even though its fixed intercept is
#'   absorbed by the thresholds). `rpar_threshold` shifts the log gaps of
#'   all thresholds above the first by one common random amount, so it
#'   changes the share of the middle categories (a common shift of all
#'   thresholds would be the same as a random propensity intercept). For
#'   `link = "mnl"`, `rpar_alloc` gives each category but the first its own
#'   random coefficients. The random coefficients are drawn once per unit
#'   of `group_var` (shared by all its rows, such as the years of a site),
#'   or once per row without `group_var`.
#' @param correlation correlation of the random parameters: `"none"`
#'   (independent), `"within"` (correlated within each equation) or
#'   `"across"` (correlated within and across the total and the
#'   allocation, which lets a unit's unobserved crash frequency and its
#'   severity propensity be correlated).
#' @param rpar_draws the number of simulation draws per unit (randomly
#'   shifted Halton points, a separate block and shift for each unit).
#' @param seed seed of the random shifts of the draws (stored in the fit,
#'   so predictions use the same draws); random by default.
#' @param .setup_only internal use: return the model setup without
#'   fitting.
#'
#' @details
#' **The two parts.** With \eqn{c = \sum_j y_j} the total, the multinomial
#' likelihood of a unit is
#' \deqn{P(c) \frac{c!}{\prod_j y_j!} \prod_j \pi_j^{y_j},}
#' the probability of the total times the multinomial probability of its
#' split. This is a proper joint likelihood of the category counts, so
#' likelihood-ratio tests, AIC and BIC are valid, and so are simulation and
#' prediction intervals. With a Poisson total the categories are
#' independent Poisson counts with means \eqn{\mu \pi_j}; with an NB2
#' total they are negative multinomial.
#'
#' **Zero totals are part of the model.** The total's distribution includes
#' zero. A unit with no events contributes \eqn{P(c = 0)} and nothing about
#' the split, because a zero total has no shares. The total is *not*
#' truncated at one: the indicator \eqn{\omega} of Shabab et al. (2026),
#' which switches off the allocation term for zero totals, does the same
#' thing. Only when the data hold units with at least one event by
#' construction should `truncated = TRUE` be used; fitting such data with
#' the untruncated total biases the count model, and `countreg_alloc()`
#' says so when it sees no zero totals.
#'
#' **Fractional likelihood.** `likelihood = "fractional"` raises
#' \eqn{\pi_j} to the observed share \eqn{y_j / c} instead of the count, as
#' in Shabab et al. (2026). A unit with one event then carries the same
#' weight in the split as a unit with a thousand. This is a
#' quasi-likelihood: the estimates are consistent, but standard errors
#' must come from the sandwich (the default), and likelihood-ratio tests,
#' AIC and BIC are not available.
#'
#' **Thresholds.** For the ordered links the thresholds are
#' \eqn{\psi_1 = v'\theta_1} and \eqn{\psi_k = \psi_{k-1} +
#' \exp(v'\theta_k)}, which keeps them ordered for every \eqn{v}. The
#' coefficients `Threshold1:` are on the scale of the latent propensity and
#' `Threshold<k>:` (k >= 2) on the log scale of the gap between thresholds.
#'
#' **Random parameters.** With `rpar_total`, `rpar_alloc` or
#' `rpar_threshold`, the likelihood of a unit is the average over
#' `rpar_draws` draws of its coefficients of the product of its rows'
#' likelihoods (maximum simulated likelihood), with analytic gradients. The
#' random coefficients are \eqn{b = L u}, with \eqn{u} standard normal and
#' \eqn{L} lower triangular; its elements are reported as `sd:` (the
#' diagonal) and `chol:` (below the diagonal), and `summary()` shows the
#' implied standard deviations and correlations. With the fractional
#' likelihood the integrand is the product of the rows' quasi-likelihoods,
#' as in Shabab et al. (2026).
#'
#' Time-varying effects ([tv()]) can be used in every formula.
#'
#' @return An object of class `"countreg_alloc"`, with methods `coef()`,
#'   `vcov()`, `logLik()`, `summary()`, `predict()`, `fitted()`,
#'   `simulate()`, `tidy()`, `glance()` and `confint()`.
#'
#' @references
#' Shabab, K. R., Bhowmik, T., Zaki, M. H. and Eluru, N. (2026). A joint
#' count and generalized ordered fractional split approach for addressing
#' temporal instability in road safety analysis. *Analytic Methods in
#' Accident Research*, 51, 100440.
#'
#' @seealso [predict.countreg_alloc()], [countreg()], [negmultinom()]
#' @examples
#' \donttest{
#' data("washington_roads")
#' w <- washington_roads
#' # Illustrative split of the total crashes into three categories
#' set.seed(1)
#' sp <- t(vapply(w$Total_crashes, function(n)
#'   stats::rmultinom(1, n, c(0.6, 0.3, 0.1)), numeric(3)))
#' w$pdo <- sp[, 1]; w$minor <- sp[, 2]; w$severe <- sp[, 3]
#' m <- countreg_alloc(cbind(pdo, minor, severe) ~ lnaadt + offset(lnlength),
#'                     alloc_formula = ~ speed50, data = w, family = "NB2")
#' summary(m)
#' head(predict(m))
#' }
#' @export
countreg_alloc <- function(formula, data, family = "NB2",
                           alloc_formula = NULL, threshold_formula = ~ 1,
                           link = c("ordered_probit", "ordered_logit", "mnl"),
                           likelihood = c("multinomial", "fractional"),
                           dis_param_formula_1 = NULL,
                           dis_param_formula_2 = NULL, offset = NULL,
                           weights = NULL, truncated = FALSE,
                           group_var = NULL, stderr = NULL,
                           method = "BFGS", max.iters = 1000,
                           start.vals = NULL, ndraws = 1500,
                           engine = "quadrature", verbose = FALSE,
                           rpar_total = NULL, rpar_alloc = NULL,
                           rpar_threshold = NULL,
                           correlation = c("none", "within", "across"),
                           rpar_draws = 500, seed = NULL,
                           .setup_only = FALSE) {
  cl <- match.call()
  link <- match.arg(link)
  correlation <- match.arg(correlation)
  .fcr_check_controls(ndraws = rpar_draws)
  likelihood <- match.arg(likelihood)
  ordered <- link != "mnl"
  family_spec <- if (inherits(family, "fcr_count_family")) family
  family <- .fcr_resolve_family(family)
  params <- get_params(family)
  .fcr_check_dis_formulas(family, params, dis_param_formula_1,
                          dis_param_formula_2)
  method <- .fcr_check_method(method)
  .fcr_check_controls(max.iters, ndraws)
  engine <- match.arg(engine, c("quadrature", "halton", "poilog"))
  if (!is.logical(truncated) || length(truncated) != 1L || is.na(truncated)) {
    stop("`truncated` must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(stderr)) {
    stderr <- if (likelihood == "fractional") {
      if (is.null(group_var)) "robust" else "cluster"
    } else "normal"
  }
  stderr <- match.arg(stderr, c("normal", "robust", "cluster"))
  if (stderr == "cluster" && is.null(group_var)) {
    stop("`stderr = \"cluster\"` needs `group_var`.", call. = FALSE)
  }
  if (likelihood == "fractional" && stderr == "normal") {
    warning("Standard errors from the Hessian are not valid for the ",
            "fractional quasi-likelihood; use stderr = \"robust\" or ",
            "\"cluster\".", call. = FALSE)
  }
  if (!ordered && !.fcr_is_intercept_only(threshold_formula)) {
    stop("`threshold_formula` applies only to the ordered links.",
         call. = FALSE)
  }
  if (!ordered && !is.null(rpar_threshold)) {
    stop("`rpar_threshold` applies only to the ordered links.",
         call. = FALSE)
  }
  has_rpar <- !is.null(rpar_total) || !is.null(rpar_alloc) ||
    !is.null(rpar_threshold)
  for (nm in c("rpar_total", "rpar_alloc", "rpar_threshold")) {
    f <- get(nm)
    if (!is.null(f) && .fcr_has_formula_offset(f)) {
      stop("`", nm, "` cannot have an offset.", call. = FALSE)
    }
  }

  # ---- Categories -----------------------------------------------------------
  formula <- stats::as.formula(formula)
  if (length(formula) != 3L) {
    stop("`formula` needs the category counts on the left, as ",
         "cbind(cat_1, ..., cat_J).", call. = FALSE)
  }
  data <- as.data.frame(data)
  lhs <- formula[[2L]]
  Y_all <- eval(lhs, data, environment(formula))
  Y_all <- as.matrix(Y_all)
  if (!is.numeric(Y_all) || nrow(Y_all) != nrow(data) || ncol(Y_all) < 2L) {
    stop("The left side of `formula` must give two or more category ",
         "count columns, as cbind(cat_1, ..., cat_J).", call. = FALSE)
  }
  categories <- colnames(Y_all)
  if (is.null(categories) || any(!nzchar(categories))) {
    categories <- if (is.call(lhs)) {
      vapply(as.list(lhs)[-1L], function(e) paste(deparse(e), collapse = ""),
             character(1))
    } else paste0("cat", seq_len(ncol(Y_all)))
  }
  colnames(Y_all) <- categories
  J <- ncol(Y_all)

  # Count model: the total as the response
  tot_name <- ".fcr_total"
  total_all <- rowSums(Y_all)
  data[[tot_name]] <- total_all
  count_formula <- formula
  count_formula[[2L]] <- as.name(tot_name)
  rhs_noff <- .fcr_drop_offsets(stats::delete.response(stats::terms(formula)))
  if (is.null(alloc_formula)) alloc_formula <- rhs_noff
  alloc_formula <- stats::as.formula(alloc_formula)
  threshold_formula <- stats::as.formula(threshold_formula)
  if (ordered && !attr(stats::terms(threshold_formula), "intercept")) {
    stop("`threshold_formula` must have an intercept.", call. = FALSE)
  }

  prep <- prepare_model_data(
    count_formula, data,
    aux_formulas = list(dis1 = dis_param_formula_1,
                        dis2 = dis_param_formula_2,
                        alloc = alloc_formula,
                        thr = if (ordered) threshold_formula,
                        rp_total = rpar_total, rp_alloc = rpar_alloc,
                        rp_thr = rpar_threshold),
    offset = offset, weights = weights, group_var = group_var)
  if (prep$n_dropped > 0) {
    message(prep$n_dropped, " row(s) with missing values were removed.")
  }
  Y <- Y_all[prep$rows_used, , drop = FALSE]
  if (any(!.fcr_is_count(Y))) {
    stop("The category counts must be non-negative whole numbers.",
         call. = FALSE)
  }
  ctot <- prep$y
  w <- prep$weights
  if (truncated && any(ctot[w > 0] == 0)) {
    stop("`truncated = TRUE`, but ", sum(ctot[w > 0] == 0), " row(s) have ",
         "a zero total. A zero-truncated total is only for data that hold ",
         "units with at least one event.", call. = FALSE)
  }
  if (!truncated && all(ctot[w > 0] > 0)) {
    message("No unit has a zero total. If the data hold only units with at ",
            "least one event by construction, use `truncated = TRUE`; ",
            "otherwise the untruncated total is right.")
  }

  X <- prep$X
  W1 <- prep$aux$dis1
  W2 <- prep$aux$dis2
  Z <- prep$aux$alloc
  if (ordered) Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  V <- if (ordered) prep$aux$thr
  if (!is.null(rpar_threshold) && J < 3L) {
    stop("`rpar_threshold` needs three or more categories (it acts on the ",
         "gaps between thresholds).", call. = FALSE)
  }

  # ---- Parameter layout -----------------------------------------------------
  S <- .fcr_alloc_layout(family, params, X, W1, W2, Z, V, J, categories,
                         ordered)
  haltons <- randtoolbox::halton(ndraws)
  dens_args <- list(haltons = haltons, normed_haltons = stats::qnorm(haltons),
                    engine = engine)
  R <- if (likelihood == "multinomial") Y else {
    sh <- Y / ifelse(ctot > 0, ctot, 1)
    sh[ctot == 0, ] <- 0
    sh
  }
  const <- if (likelihood == "multinomial") {
    lgamma(ctot + 1) - rowSums(lgamma(Y + 1))
  } else rep(0, length(ctot))
  E <- list(family = family, link = link, y = ctot, R = R, const = const,
            w = w, X = X, off = prep$offset, W1 = W1, W2 = W2, Z = Z, V = V,
            J = J, S = S, truncated = truncated, dens_args = dens_args,
            categories = categories)

  # ---- Random parameters ----------------------------------------------------
  RS <- NULL
  draw_seed <- NULL
  if (has_rpar) {
    RS <- .fcr_alloc_rlayout(prep$aux$rp_total, prep$aux$rp_alloc,
                             prep$aux$rp_thr, categories, ordered,
                             correlation, length(S$names))
    S$names <- c(S$names, RS$names)
    unit <- if (is.null(group_var)) factor(seq_len(nrow(X))) else prep$group
    E$g <- as.integer(unit)
    E$G <- nlevels(unit)
    E$units <- levels(unit)
    draw_seed <- seed %||% sample.int(.Machine$integer.max, 1L)
    E$UA <- .fcr_alloc_draws(E$G, RS$K, rpar_draws, draw_seed)
    E$Xr <- prep$aux$rp_total
    E$Zr <- prep$aux$rp_alloc
    E$Vr <- prep$aux$rp_thr
    E$RS <- RS
    S$RS <- RS
    E$S <- S
  }

  if (.setup_only) return(E)

  # ---- Fit -------------------------------------------------------------------
  start <- if (is.null(start.vals)) .fcr_alloc_start(E) else
    .fcr_check_start(start.vals, S$names)
  names(start) <- S$names
  fn <- if (has_rpar) function(p) .fcr_alloc_sim_ll(p, E, grad = FALSE)$ll else
    function(p) .fcr_alloc_ll(p, E, grad = FALSE)$ll
  gr <- if (has_rpar) function(p) .fcr_alloc_sim_ll(p, E, grad = TRUE)$G else
    function(p) .fcr_alloc_ll(p, E, grad = TRUE)$G
  fit <- maxLik::maxLik(fn, grad = gr, start = start, method = method,
                        finalHessian = method == "BHHH",
                        control = list(iterlim = max.iters,
                                       printLevel = if (verbose) 2 else 0))
  fit <- .fcr_set_final_hessian(fit, fn, gr, method)
  est <- fit$estimate
  names(est) <- S$names
  conv <- isTRUE(fit$code %in% c(0L, 1L, 2L, 8L)) ||
    grepl("converge|successful", fit$message %||% "", ignore.case = TRUE)
  if (!conv) {
    warning("The optimizer did not report convergence: ", fit$message,
            call. = FALSE)
  }
  U <- gr(est)
  colnames(U) <- S$names
  Hv <- .fcr_hessian_vcov(list(hessian = fit$hessian, estimate = est))
  V_cov <- switch(stderr,
    normal = Hv,
    robust = if (has_rpar) {
      # The units (groups) are the independent observations
      Hv %*% crossprod(U) %*% Hv
    } else {
      keep <- w > 0
      meat <- crossprod(U[keep, , drop = FALSE] / sqrt(w[keep]))
      Hv %*% meat %*% Hv
    },
    cluster = {
      Ug <- if (has_rpar) U else rowsum(U, prep$group)
      G <- nrow(Ug)
      Hv %*% (crossprod(Ug) * G / max(G - 1, 1)) %*% Hv
    })
  dimnames(V_cov) <- list(S$names, S$names)

  out <- list(
    call = cl, coefficients = est, vcov = V_cov, loglik = sum(fn(est)),
    gradientObs = U, hessian = fit$hessian, converged = conv,
    iterations = unname(fit$iterations %||% NA), message = fit$message,
    family = family, family_spec = family_spec, link = link,
    likelihood = likelihood, truncated = truncated, stderr = stderr,
    categories = categories, layout = S, formula = formula,
    alloc_formula = alloc_formula,
    threshold_formula = if (ordered) threshold_formula,
    dis_param_formula_1 = prep$aux_formulas$dis1,
    dis_param_formula_2 = prep$aux_formulas$dis2,
    group_var = group_var, nobs = sum(w > 0),
    n_units = if (has_rpar) E$G else sum(w), n_obs_ic = sum(w),
    rpar = if (has_rpar) list(
      layout = RS, formulas = list(total = rpar_total, alloc = rpar_alloc,
                                   threshold = rpar_threshold),
      correlation = correlation, draws = rpar_draws, seed = draw_seed,
      units = E$units, unit_of_row = E$g,
      posterior = .fcr_alloc_sim_ll(est, E, grad = TRUE)$W),
    y = Y, total = ctot, weights = w, offset = prep$offset,
    rows_used = prep$rows_used, blueprint = prep$blueprint,
    data_used = data[prep$rows_used, setdiff(names(data), tot_name),
                     drop = FALSE],
    engine = engine, ndraws = ndraws,
    fitted.values = NULL)
  class(out) <- "countreg_alloc"
  out$fitted.values <- if (has_rpar) {
    predict(out, newdata = out$data_used)
  } else .fcr_alloc_predict_core(est, E)$counts
  out
}

.fcr_is_intercept_only <- function(f) {
  if (is.null(f)) return(TRUE)
  tt <- stats::terms(stats::as.formula(f))
  !length(attr(tt, "term.labels")) && attr(tt, "intercept") == 1L
}

# A one-sided formula of the terms of `tt` without offset() terms.
.fcr_drop_offsets <- function(tt) {
  labs <- attr(tt, "term.labels")
  env <- environment(tt) %||% parent.frame()
  if (!length(labs)) {
    f <- if (attr(tt, "intercept")) ~ 1 else ~ 0
  } else {
    f <- stats::reformulate(labs, intercept = attr(tt, "intercept") == 1L)
  }
  environment(f) <- env
  f
}

# Parameter blocks and names.
.fcr_alloc_layout <- function(family, params, X, W1, W2, Z, V, J,
                              categories, ordered) {
  kx <- ncol(X)
  n1 <- if (is.null(params[[1]])) 0L else if (is.null(W1)) 1L else ncol(W1)
  n2 <- if (is.null(params[[2]])) 0L else if (is.null(W2)) 1L else ncol(W2)
  kz <- ncol(Z)
  nm <- colnames(X)
  if (n1) nm <- c(nm, if (is.null(W1)) params[[1]] else
    paste0(params[[1]], ":", colnames(W1)))
  if (n2) nm <- c(nm, if (is.null(W2)) params[[2]] else
    paste0(params[[2]], ":", colnames(W2)))
  pos <- length(nm)
  if (ordered) {
    if (!kz) {
      alloc_names <- character(0)
    } else alloc_names <- paste0("Allocation:", colnames(Z))
    kv <- ncol(V)
    thr_names <- unlist(lapply(seq_len(J - 1L), function(k)
      paste0("Threshold", k, ":", colnames(V))))
  } else {
    alloc_names <- unlist(lapply(categories[-1L], function(cn)
      paste0("Allocation[", cn, "]:", colnames(Z))))
    kv <- 0L
    thr_names <- character(0)
  }
  idx_beta <- seq_len(kx)
  idx_a <- kx + seq_len(n1)
  idx_s <- kx + n1 + seq_len(n2)
  idx_g <- pos + seq_along(alloc_names)
  idx_t <- pos + length(alloc_names) + seq_along(thr_names)
  list(names = c(nm, alloc_names, thr_names), kx = kx, n1 = n1, n2 = n2,
       kz = kz, kv = kv, idx_beta = idx_beta, idx_a = idx_a, idx_s = idx_s,
       idx_g = idx_g, idx_t = idx_t, n_count = pos, ordered = ordered,
       J = J)
}

# Distribution parameters of the total at parameter vector p.
.fcr_alloc_dispar <- function(p, S, W1, W2, n) {
  alpha <- if (S$n1) {
    if (is.null(W1)) rep(exp(p[S$idx_a]), n) else
      exp(as.vector(W1 %*% p[S$idx_a]))
  }
  sigma <- if (S$n2) {
    if (is.null(W2)) rep(exp(p[S$idx_s]), n) else
      exp(as.vector(W2 %*% p[S$idx_s]))
  }
  list(alpha = alpha, sigma = sigma)
}

.fcr_alloc_cdf <- function(link) {
  switch(link,
    ordered_probit = list(F = stats::pnorm, f = stats::dnorm),
    ordered_logit = list(F = stats::plogis, f = stats::dlogis))
}

# Thresholds (n x (J - 1)) and the exponentiated gaps (n x (J - 1); the
# first column is 1). `add_gap` (a vector) is added to the linear
# predictor of every gap (a random threshold term).
.fcr_alloc_thresholds <- function(p, S, V, add_gap = NULL) {
  J <- S$J
  th <- matrix(p[S$idx_t], S$kv, J - 1L)
  lin <- V %*% th
  glin <- lin[, -1L, drop = FALSE]
  if (!is.null(add_gap) && ncol(glin)) glin <- glin + add_gap
  gap <- cbind(1, exp(glin))
  psi <- lin
  if (J > 2L) {
    for (k in 2:(J - 1L)) psi[, k] <- psi[, k - 1L] + gap[, k]
  }
  list(psi = psi, gap = gap)
}

# Allocation probabilities (n x J) and, for the ordered links, the pieces
# of their derivatives. `add$xb` is added to the propensity (ordered: a
# vector; mnl: an n x (J - 1) matrix), `add$gap` to the threshold gaps.
.fcr_alloc_probs <- function(p, E, add = NULL) {
  S <- E$S
  n <- nrow(E$X)
  if (!S$ordered) {
    eta <- matrix(0, n, S$J)
    if (S$kz) {
      D <- matrix(p[S$idx_g], S$kz, S$J - 1L)
      eta[, -1L] <- E$Z %*% D
    }
    if (!is.null(add$xb)) eta[, -1L] <- eta[, -1L] + add$xb
    eta <- eta - apply(eta, 1L, max)
    ex <- exp(eta)
    pr <- ex / rowSums(ex)
    return(list(pi = pr, logpi = log(pr)))
  }
  cf <- .fcr_alloc_cdf(E$link)
  xb <- if (S$kz) as.vector(E$Z %*% p[S$idx_g]) else rep(0, n)
  if (!is.null(add$xb)) xb <- xb + add$xb
  th <- .fcr_alloc_thresholds(p, S, E$V, add$gap)
  u <- th$psi - xb
  # P(Y = j) = F(u_j) - F(u_(j-1)), from the upper tails when u_(j-1) > 0
  Fl <- cbind(0, cf$F(u), 1)
  Fu <- cbind(1, cf$F(-u), 0)
  lo <- cbind(-Inf, u)
  pr <- ifelse(lo > 0, Fu[, -(S$J + 1L), drop = FALSE] - Fu[, -1L, drop = FALSE],
               Fl[, -1L, drop = FALSE] - Fl[, -(S$J + 1L), drop = FALSE])
  pr <- pmax(pr, 0)
  list(pi = pr, logpi = log(pr), u = u, f = cf$f(u), gap = th$gap)
}

# Row log-likelihood contributions (weighted) and the derivatives of each
# row's contribution with respect to the linear predictors: the count's
# log mean (s_eta) and log distribution parameters (s_a, s_s), the
# propensity (dxb; an n x (J - 1) matrix for mnl), the fixed threshold
# coefficients (thrG, n x length(idx_t)) and a common shift of the log
# gaps (dgap). `add` holds random terms: add$eta (the count's log mean),
# add$xb and add$gap (see .fcr_alloc_probs()).
.fcr_alloc_rowparts <- function(p, E, add = NULL, grad = TRUE) {
  S <- E$S
  n <- length(E$y)
  fam <- .fcr_get_family(E$family)
  lp <- as.vector(E$X %*% p[S$idx_beta]) + E$off
  if (!is.null(add$eta)) lp <- lp + add$eta
  mu <- exp(lp)
  dp <- .fcr_alloc_dispar(p, S, E$W1, E$W2, n)
  logf <- function(yv) {
    out <- tryCatch(.fcr_quiet_dens(as.vector(do.call(fam$dens, c(
      list(y = yv, predicted = mu, alpha = dp$alpha, sigma = dp$sigma),
      E$dens_args, list(log = TRUE)))), log = TRUE),
      error = function(e) rep(NA_real_, n))
    ok <- .fcr_family_valid(E$family, mu, dp$alpha, dp$sigma)
    out[!ok] <- -Inf
    out
  }
  lc <- logf(E$y)
  if (E$truncated) {
    l0 <- logf(rep(0, n))
    lc <- lc - .fcr_log1mexp(l0)
  }
  pr <- .fcr_alloc_probs(p, E, add)
  la <- rowSums(ifelse(E$R > 0, E$R * pr$logpi, 0))
  ll <- E$w * (lc + la + E$const)
  ll[E$w == 0] <- 0
  if (!grad) return(list(ll = ll))

  sc <- .fcr_family_score(E$family, E$y, mu, dp$alpha, dp$sigma,
                          dens_args = E$dens_args)
  s_eta <- sc$eta
  s_a <- sc$alpha
  s_s <- sc$sigma
  if (E$truncated) {
    s0 <- .fcr_family_score(E$family, rep(0, n), mu, dp$alpha, dp$sigma,
                            dens_args = E$dens_args)
    # d log(1 - P0) = -P0 / (1 - P0) d log P0
    r0 <- exp(l0 - .fcr_log1mexp(l0))
    s_eta <- s_eta + r0 * s0$eta
    if (!is.null(s_a)) s_a <- s_a + r0 * s0$alpha
    if (!is.null(s_s)) s_s <- s_s + r0 * s0$sigma
  }
  out <- list(ll = ll, s_eta = E$w * s_eta,
              s_a = if (!is.null(s_a)) E$w * s_a,
              s_s = if (!is.null(s_s)) E$w * s_s)
  if (!S$ordered) {
    rs <- rowSums(E$R)
    out$dxb <- E$w * (E$R[, -1L, drop = FALSE] -
                        rs * pr$pi[, -1L, drop = FALSE])
    return(out)
  }
  J <- S$J
  ratio <- ifelse(E$R > 0, E$R / pr$pi, 0)
  # D_k = f(u_k) (r_k / pi_k - r_(k+1) / pi_(k+1)), k = 1, ..., J - 1
  D <- pr$f * (ratio[, -J, drop = FALSE] - ratio[, -1L, drop = FALSE])
  out$dxb <- -E$w * rowSums(D)
  # psi_k depends on theta_m for m <= k
  tail_sum <- D[, rev(seq_len(J - 1L)), drop = FALSE]
  tail_sum <- t(apply(tail_sum, 1L, cumsum))
  if (J - 1L == 1L) tail_sum <- t(tail_sum)
  tail_sum <- tail_sum[, rev(seq_len(J - 1L)), drop = FALSE]
  tg <- tail_sum * pr$gap * E$w
  thrG <- matrix(0, n, length(S$idx_t))
  for (m in seq_len(J - 1L)) {
    thrG[, (m - 1L) * S$kv + seq_len(S$kv)] <- tg[, m] * E$V
  }
  out$thrG <- thrG
  out$dgap <- if (J > 2L) rowSums(tg[, -1L, drop = FALSE]) else rep(0, n)
  out
}

# Gradient rows of the fixed parameters from .fcr_alloc_rowparts().
.fcr_alloc_fixed_grad <- function(pt, E, np) {
  S <- E$S
  G <- matrix(0, length(E$y), np)
  G[, S$idx_beta] <- pt$s_eta * E$X
  if (S$n1) G[, S$idx_a] <- pt$s_a * (if (is.null(E$W1)) 1 else E$W1)
  if (S$n2) G[, S$idx_s] <- pt$s_s * (if (is.null(E$W2)) 1 else E$W2)
  if (S$kz) {
    if (S$ordered) {
      G[, S$idx_g] <- pt$dxb * E$Z
    } else {
      for (j in 2:S$J) {
        cols <- S$idx_g[(j - 2L) * S$kz + seq_len(S$kz)]
        G[, cols] <- pt$dxb[, j - 1L] * E$Z
      }
    }
  }
  if (S$ordered) G[, S$idx_t] <- pt$thrG
  G
}

# Log-likelihood without random parameters: one contribution per row.
.fcr_alloc_ll <- function(p, E, grad = FALSE) {
  pt <- .fcr_alloc_rowparts(p, E, grad = grad)
  ll <- pt$ll
  ll[!is.finite(ll)] <- if (grad) 0 else NA_real_
  if (!grad) return(list(ll = ll))
  G <- .fcr_alloc_fixed_grad(pt, E, length(p))
  G[!is.finite(G)] <- 0
  list(ll = ll, G = G)
}

# ---- Random parameters ------------------------------------------------------

# Layout of the random terms: one entry per random coefficient (its
# equation, its column in the random design, its category for mnl) and the
# Cholesky factor's free elements.
.fcr_alloc_rlayout <- function(Xr, Zr, Vr, categories, ordered, correlation,
                               offset) {
  terms <- list()
  add <- function(type, M, cat = NA_integer_, prefix) {
    if (is.null(M) || !ncol(M)) return()
    for (c in seq_len(ncol(M))) {
      terms[[length(terms) + 1L]] <<- list(
        type = type, col = c, cat = cat,
        label = paste0(prefix, colnames(M)[c]))
    }
  }
  add("total", Xr, prefix = "Total:")
  if (ordered) {
    add("alloc", Zr, prefix = "Allocation:")
  } else if (!is.null(Zr)) {
    for (j in 2:length(categories)) {
      add("alloc", Zr, cat = j,
          prefix = paste0("Allocation[", categories[j], "]:"))
    }
  }
  add("gap", Vr, prefix = "Threshold gaps:")
  K <- length(terms)
  if (!K) return(NULL)
  type <- vapply(terms, `[[`, "", "type")
  eq <- match(type, c("total", "alloc", "gap"))
  if (!ordered) {
    # Each mnl category is its own equation for "within"
    cats <- vapply(terms, `[[`, 0L, "cat")
    eq[type == "alloc"] <- 10L + cats[type == "alloc"]
  }
  pairs <- do.call(rbind, lapply(seq_len(K), function(k) {
    l <- switch(correlation,
      none = k,
      within = which(eq == eq[k] & seq_len(K) <= k),
      across = seq_len(k))
    cbind(k = k, l = l)
  }))
  labels <- vapply(terms, `[[`, "", "label")
  nm <- ifelse(pairs[, "k"] == pairs[, "l"],
               paste0("sd:", labels[pairs[, "k"]]),
               paste0("chol:", labels[pairs[, "k"]], "|",
                      labels[pairs[, "l"]]))
  list(K = K, terms = terms, type = type, labels = labels, pairs = pairs,
       names = nm, idx = offset + seq_len(nrow(pairs)),
       col = vapply(terms, `[[`, 0L, "col"),
       cat = vapply(terms, `[[`, 0L, "cat"))
}

.fcr_alloc_chol <- function(p, RS) {
  L <- matrix(0, RS$K, RS$K)
  L[RS$pairs] <- p[RS$idx]
  L
}

# Random terms of every row for coefficient values B (rows x K).
.fcr_alloc_add <- function(B, E) {
  RS <- E$RS
  n <- nrow(B)
  add <- list()
  te <- which(RS$type == "total")
  if (length(te)) add$eta <- rowSums(E$Xr[, RS$col[te], drop = FALSE] *
                                       B[, te, drop = FALSE])
  ta <- which(RS$type == "alloc")
  if (length(ta)) {
    if (E$S$ordered) {
      add$xb <- rowSums(E$Zr[, RS$col[ta], drop = FALSE] *
                          B[, ta, drop = FALSE])
    } else {
      add$xb <- matrix(0, n, E$S$J - 1L)
      for (j in 2:E$S$J) {
        tj <- ta[RS$cat[ta] == j]
        add$xb[, j - 1L] <- rowSums(E$Zr[, RS$col[tj], drop = FALSE] *
                                      B[, tj, drop = FALSE])
      }
    }
  }
  tg <- which(RS$type == "gap")
  if (length(tg)) add$gap <- rowSums(E$Vr[, RS$col[tg], drop = FALSE] *
                                       B[, tg, drop = FALSE])
  add
}

# Derivatives of the row contributions with respect to each random
# coefficient (rows x K).
.fcr_alloc_rand_deriv <- function(pt, E) {
  RS <- E$RS
  Db <- matrix(0, length(E$y), RS$K)
  for (k in seq_len(RS$K)) {
    c <- RS$col[k]
    Db[, k] <- switch(RS$type[k],
      total = pt$s_eta * E$Xr[, c],
      alloc = if (E$S$ordered) pt$dxb * E$Zr[, c] else
        pt$dxb[, RS$cat[k] - 1L] * E$Zr[, c],
      gap = pt$dgap * E$Vr[, c])
  }
  Db
}

# Simulated log-likelihood: one contribution per unit (group),
# log((1/R) sum_r prod_(rows of the unit) P(row | b_r)), and its gradient
# through the posterior weights of the draws.
.fcr_alloc_sim_ll <- function(p, E, grad = FALSE) {
  RS <- E$RS
  L <- .fcr_alloc_chol(p, RS)
  g <- E$g
  G <- E$G
  R <- dim(E$UA)[3L]
  LL <- matrix(NA_real_, G, R)
  for (r in seq_len(R)) {
    B <- (matrix(E$UA[, , r], G, RS$K) %*% t(L))[g, , drop = FALSE]
    ll <- .fcr_alloc_rowparts(p, E, .fcr_alloc_add(B, E), grad = FALSE)$ll
    LL[, r] <- rowsum(ll, g, reorder = TRUE)[, 1L]
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
    pt <- .fcr_alloc_rowparts(p, E, .fcr_alloc_add(B, E), grad = TRUE)
    Gr <- .fcr_alloc_fixed_grad(pt, E, np)
    Db <- .fcr_alloc_rand_deriv(pt, E)
    Urow <- Ug[g, , drop = FALSE]
    Gr[, RS$idx] <- Db[, RS$pairs[, "k"], drop = FALSE] *
      Urow[, RS$pairs[, "l"], drop = FALSE]
    Gr[!is.finite(Gr)] <- 0
    Gg <- Gg + rowsum(Gr, g, reorder = TRUE) * W[, r]
  }
  lg[!is.finite(lg)] <- 0
  list(ll = lg, G = Gg, W = W)
}

# Standard normal draws for each unit: G x K x R, from randomly shifted
# Halton points (each unit its own block and shift).
.fcr_alloc_draws <- function(G, K, R, seed) {
  h <- .fcr_unit_halton(seq_len(G), R, K, scrambled = TRUE, seed = seed)
  UA <- array(0, c(G, K, R))
  for (k in seq_len(K)) UA[, k, ] <- stats::qnorm(h[[k]])
  UA
}

# log(1 - exp(x)) for x <= 0.
.fcr_log1mexp <- function(x) {
  ifelse(x > -log(2), log(-expm1(x)), log1p(-exp(x)))
}

# Starting values: a Poisson fit of the total, the family's dispersion
# start, and the observed shares for the allocation.
.fcr_alloc_start <- function(E) {
  S <- E$S
  pm <- stats::glm.fit(E$X, E$y, weights = E$w, offset = E$off,
                       family = stats::poisson())
  b <- pm$coefficients
  b[is.na(b)] <- 0
  st <- numeric(length(S$names))
  st[S$idx_beta] <- b
  fs <- .fcr_get_family(E$family)$start
  if (S$n1) {
    st[S$idx_a] <- 0
    a0 <- .fcr_moment_disp_start(E$family, E$y, pm$fitted.values, E$w) %||%
      fs[1L]
    pos <- if (is.null(E$W1)) 1L else match("(Intercept)", colnames(E$W1))
    if (!is.na(pos)) st[S$idx_a[pos]] <- a0
  }
  if (S$n2) {
    st[S$idx_s] <- 0
    pos <- if (is.null(E$W2)) 1L else match("(Intercept)", colnames(E$W2))
    if (!is.na(pos)) st[S$idx_s[pos]] <- rep_len(fs, 2L)[2L]
  }
  share <- colSums(E$R * E$w)
  share <- pmax(share / sum(share), 1e-4)
  share <- share / sum(share)
  if (S$ordered) {
    Finv <- switch(E$link, ordered_probit = stats::qnorm,
                   ordered_logit = stats::qlogis)
    q <- Finv(pmin(cumsum(share)[-S$J], 1 - 1e-6))
    ic <- match("(Intercept)", colnames(E$V))
    if (!is.na(ic)) {
      for (k in seq_len(S$J - 1L)) {
        val <- if (k == 1L) q[1L] else log(max(q[k] - q[k - 1L], 1e-3))
        st[S$idx_t[(k - 1L) * S$kv + ic]] <- val
      }
    }
  } else {
    ic <- match("(Intercept)", colnames(E$Z))
    if (!is.na(ic)) {
      for (j in 2:S$J) {
        st[S$idx_g[(j - 2L) * S$kz + ic]] <- log(share[j] / share[1L])
      }
    }
  }
  if (!is.null(E$RS)) {
    dg <- E$RS$pairs[, "k"] == E$RS$pairs[, "l"]
    st[E$RS$idx] <- ifelse(dg, 0.1, 0)
  }
  st
}

# Predictions at parameters `p` for the design in `E`, with random terms
# `add` (see .fcr_alloc_add()).
.fcr_alloc_predict_core <- function(p, E, add = NULL) {
  S <- E$S
  n <- nrow(E$X)
  lp <- as.vector(E$X %*% p[S$idx_beta]) + E$off
  if (!is.null(add$eta)) lp <- lp + add$eta
  mu <- exp(lp)
  dp <- .fcr_alloc_dispar(p, S, E$W1, E$W2, n)
  mean_total <- .fcr_family_mean(E$family, mu, dp$alpha, dp$sigma)
  pr <- .fcr_alloc_probs(p, E, add)$pi
  colnames(pr) <- E$categories
  counts <- mean_total * pr
  list(mu = mu, alpha = dp$alpha, sigma = dp$sigma, total = mean_total,
       prob = pr, counts = counts)
}

# The design of a fit for new data (or the estimation data).
.fcr_alloc_design <- function(object, newdata = NULL) {
  S <- object$layout
  if (is.null(newdata)) {
    stop("internal: newdata needed", call. = FALSE)
  }
  fd <- forge_new_data(object$blueprint, newdata)
  Z <- fd$aux$alloc
  if (S$ordered) Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  haltons <- randtoolbox::halton(object$ndraws)
  list(family = object$family, link = object$link, X = fd$X,
       off = fd$offset, W1 = fd$aux$dis1, W2 = fd$aux$dis2, Z = Z,
       V = if (S$ordered) fd$aux$thr, S = S, truncated = object$truncated,
       categories = object$categories,
       dens_args = list(haltons = haltons,
                        normed_haltons = stats::qnorm(haltons),
                        engine = object$engine),
       y = rep(0, nrow(fd$X)), Xr = fd$aux$rp_total, Zr = fd$aux$rp_alloc,
       Vr = fd$aux$rp_thr, RS = object$rpar$layout)
}

# P(total = 0) at the predictions `pc` of .fcr_alloc_predict_core().
.fcr_alloc_p0 <- function(object, pc, E) {
  fam <- .fcr_get_family(object$family)
  exp(.fcr_quiet_dens(as.vector(do.call(fam$dens, c(
    list(y = rep(0, length(pc$mu)), predicted = pc$mu, alpha = pc$alpha,
         sigma = pc$sigma), E$dens_args, list(log = TRUE)))), log = TRUE))
}

# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------

#' Predictions from total-plus-allocation models
#'
#' @param object a [countreg_alloc()] fit.
#' @param newdata optional new data; the estimation data by default.
#' @param type `"counts"` (default): the expected count in each category,
#'   \eqn{E[c]\,\pi_j}; `"prob"`: the allocation probabilities
#'   \eqn{\pi_j}; `"total"`: the expected total \eqn{E[c]};
#'   `"total_positive"`: the expected total given at least one event,
#'   \eqn{E[c] / (1 - P(0))}.
#' @param conditional for `type = "counts"`: `TRUE` gives the expected
#'   counts given at least one event, \eqn{E[c \mid c > 0]\,\pi_j}. The
#'   default is `FALSE`, the unconditional expected counts, also for a fit
#'   with `truncated = TRUE` (whose estimated distribution is that of the
#'   untruncated total).
#' @param random for fits with random parameters: `"marginal"` (default)
#'   averages over the distribution of the random parameters (the
#'   prediction for a new unit); `"posterior"` averages over each unit's
#'   posterior distribution given its counts in the estimation data
#'   (Empirical Bayes; units not in the estimation data get the marginal
#'   prediction). `newdata` then needs the `group_var` column(s).
#' @param ... unused.
#' @return A matrix with one column per category (`"counts"`, `"prob"`), or
#'   a vector (`"total"`, `"total_positive"`). With random parameters,
#'   `"prob"` is the average allocation probability and `"counts"` the
#'   average expected count over the draws.
#' @export
predict.countreg_alloc <- function(object, newdata = NULL,
                                   type = c("counts", "prob", "total",
                                            "total_positive"),
                                   conditional = FALSE,
                                   random = c("marginal", "posterior"),
                                   ...) {
  type <- match.arg(type)
  random <- match.arg(random)
  if (is.null(newdata)) {
    newdata <- .fcr_alloc_rebuild_data(object)
  }
  E <- .fcr_alloc_design(object, newdata)
  p <- object$coefficients
  need_p0 <- type == "total_positive" || (type == "counts" && conditional)
  if (is.null(object$rpar)) {
    pc <- .fcr_alloc_predict_core(p, E)
    p0 <- if (need_p0) .fcr_alloc_p0(object, pc, E)
    return(switch(type,
      counts = if (conditional) pc$counts / (1 - p0) else pc$counts,
      prob = pc$prob,
      total = pc$total,
      total_positive = pc$total / (1 - p0)))
  }
  # Random parameters: average over the draws, with equal weights
  # (marginal) or with each unit's posterior weights (posterior; units not
  # in the estimation data get the marginal prediction)
  rp <- object$rpar
  RS <- rp$layout
  L <- .fcr_alloc_chol(p, RS)
  n <- nrow(E$X)
  R <- rp$draws
  Wrow <- matrix(1 / R, n, R)
  Urow <- NULL
  if (random == "posterior") {
    if (is.null(object$group_var) ||
        !all(object$group_var %in% names(newdata))) {
      stop("Posterior predictions need the grouping column(s) ",
           toString(object$group_var %||% "(group_var)"), " in `newdata` ",
           "(and a fit with `group_var`).", call. = FALSE)
    }
    key <- as.character(interaction(newdata[, object$group_var,
                                            drop = FALSE], drop = TRUE,
                                    sep = ":"))
    ui <- match(key, rp$units)
    seen <- !is.na(ui)
    if (any(seen)) {
      UA <- .fcr_alloc_draws(length(rp$units), RS$K, R, rp$seed)
      Wrow[seen, ] <- rp$posterior[ui[seen], , drop = FALSE]
    }
  }
  # Marginal draws: plain Halton points, the same for every row
  Um <- matrix(stats::qnorm(.fcr_halton(R, RS$K)), R, RS$K)
  acc <- list(counts = 0, prob = 0, total = 0, p0 = 0)
  for (r in seq_len(R)) {
    B <- matrix(Um[r, ] %*% t(L), n, RS$K, byrow = TRUE)
    if (random == "posterior" && any(seen)) {
      B[seen, ] <- (matrix(UA[ui[seen], , r], sum(seen), RS$K)) %*% t(L)
    }
    pc <- .fcr_alloc_predict_core(p, E, .fcr_alloc_add(B, E))
    wr <- Wrow[, r]
    acc$counts <- acc$counts + wr * pc$counts
    acc$prob <- acc$prob + wr * pc$prob
    acc$total <- acc$total + wr * pc$total
    if (need_p0) acc$p0 <- acc$p0 + wr * .fcr_alloc_p0(object, pc, E)
  }
  colnames(acc$counts) <- colnames(acc$prob) <- object$categories
  switch(type,
    counts = if (conditional) acc$counts / (1 - acc$p0) else acc$counts,
    prob = acc$prob,
    total = acc$total,
    total_positive = acc$total / (1 - acc$p0))
}

# The estimation rows of the data, for predictions without newdata.
.fcr_alloc_rebuild_data <- function(object) {
  d <- object$data_used
  if (!is.null(d)) return(d)
  d <- tryCatch(eval(object$call$data, environment(object$formula)),
                error = function(e) NULL)
  if (is.null(d)) {
    stop("The estimation data are not available; supply `newdata`.",
         call. = FALSE)
  }
  as.data.frame(d)[object$rows_used, , drop = FALSE]
}

#' @export
fitted.countreg_alloc <- function(object, ...) object$fitted.values

#' @export
coef.countreg_alloc <- function(object, ...) object$coefficients

#' @export
vcov.countreg_alloc <- function(object, ...) object$vcov

#' @export
nobs.countreg_alloc <- function(object, ...) object$nobs

#' @export
logLik.countreg_alloc <- function(object, ...) {
  if (object$likelihood == "fractional") {
    stop("The fractional likelihood is a quasi-likelihood: logLik(), AIC() ",
         "and BIC() are not available. Use likelihood = \"multinomial\" for ",
         "a proper likelihood.", call. = FALSE)
  }
  structure(object$loglik, df = length(object$coefficients),
            nobs = object$n_obs_ic %||% object$n_units, class = "logLik")
}

#' @export
print.countreg_alloc <- function(x, ...) {
  cat("Total-plus-allocation count model\n")
  cat("Total:", .fcr_family_label(x$family),
      if (x$truncated) "(zero-truncated)", " Allocation:", x$link,
      " Likelihood:", x$likelihood, "\n")
  cat("Categories:", paste(x$categories, collapse = ", "), "\n\n")
  print(round(x$coefficients, 4))
  invisible(x)
}

#' @export
summary.countreg_alloc <- function(object, ...) {
  est <- object$coefficients
  se <- sqrt(pmax(diag(object$vcov), 0))
  z <- est / se
  tab <- data.frame(Estimate = est, `Std. Error` = se, `z value` = z,
                    `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)),
                    check.names = FALSE)
  S <- object$layout
  part <- rep("Total", length(est))
  part[c(S$idx_a, S$idx_s)] <- "Distribution"
  part[S$idx_g] <- "Allocation"
  part[S$idx_t] <- "Thresholds"
  rp <- object$rpar
  if (!is.null(rp)) {
    part[rp$layout$idx] <- "Random parameters (Cholesky factor)"
  }
  structure(list(table = tab, part = part, object = object),
            class = "summary.countreg_alloc")
}

#' @export
print.summary.countreg_alloc <- function(x, digits = 4, ...) {
  o <- x$object
  cat("Total-plus-allocation count model\n\n")
  cat("Total:      ", .fcr_family_label(o$family),
      if (o$truncated) "(zero-truncated)", "\n")
  cat("Allocation: ", o$link, "with", o$likelihood, "likelihood\n")
  cat("Categories: ", paste(o$categories, collapse = " < "), "\n")
  nz <- sum(o$total == 0 & o$weights > 0)
  if (is.null(o$rpar)) {
    cat("Units:      ", o$n_units, " (", nz, " with a zero total)\n", sep = "")
  } else {
    cat("Units:      ", o$n_units, " (", o$nobs, " rows, ", nz,
        " with a zero total)\n", sep = "")
  }
  cat("Std. errors:", o$stderr, "\n\n")
  for (pt in unique(x$part)) {
    cat(pt, ":\n", sep = "")
    stats::printCoefmat(x$table[x$part == pt, , drop = FALSE],
                        digits = digits, signif.legend = FALSE)
    cat("\n")
  }
  if (!is.null(o$rpar)) {
    vc <- .fcr_alloc_varcorr(o)
    cat("Random parameters: standard deviations and correlations\n")
    sdv <- sqrt(diag(vc))
    out <- cbind(SD = sdv)
    if (length(sdv) > 1L) {
      cr <- stats::cov2cor(vc)
      cr[upper.tri(cr, diag = TRUE)] <- NA
      out <- cbind(out, cr[, -ncol(cr), drop = FALSE])
    }
    print(round(out, 3), na.print = "")
    cat("(", o$rpar$draws, " draws per ",
        if (is.null(o$group_var)) "row" else
          paste0("unit of ", toString(o$group_var)),
        "; ", o$n_units, " units)\n\n", sep = "")
  }
  if (o$likelihood == "multinomial") {
    k <- length(o$coefficients)
    cat("Log-likelihood:", format(o$loglik, digits = 8), " AIC:",
        format(-2 * o$loglik + 2 * k, digits = 8), " BIC:",
        format(-2 * o$loglik + log(o$n_obs_ic %||% o$n_units) * k,
               digits = 8), "\n")
  } else {
    cat("Fractional quasi-likelihood:", format(o$loglik, digits = 8),
        "(no AIC, BIC or likelihood-ratio tests)\n")
  }
  if (!o$converged) cat("Warning: the optimizer did not converge.\n")
  invisible(x)
}

# Covariance matrix of the random coefficients, LL'.
.fcr_alloc_varcorr <- function(object) {
  RS <- object$rpar$layout
  L <- .fcr_alloc_chol(object$coefficients, RS)
  vc <- L %*% t(L)
  dimnames(vc) <- list(RS$labels, RS$labels)
  vc
}

#' @export
tidy.countreg_alloc <- function(x, conf.int = FALSE, conf.level = 0.95,
                                ...) {
  s <- summary(x)
  out <- data.frame(term = rownames(s$table),
                    component = ifelse(startsWith(s$part, "Random"), "random",
                       tolower(s$part)),
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
glance.countreg_alloc <- function(x, ...) {
  k <- length(x$coefficients)
  multi <- x$likelihood == "multinomial"
  tibble::tibble(
    logLik = if (multi) x$loglik else NA_real_,
    AIC = if (multi) -2 * x$loglik + 2 * k else NA_real_,
    BIC = if (multi) -2 * x$loglik + log(x$n_obs_ic %||% x$n_units) * k else
      NA_real_,
    nobs = x$nobs, npar = k, converged = x$converged)
}

#' @export
simulate.countreg_alloc <- function(object, nsim = 1, seed = NULL,
                                    newdata = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(newdata)) newdata <- .fcr_alloc_rebuild_data(object)
  E <- .fcr_alloc_design(object, newdata)
  p <- object$coefficients
  n <- nrow(E$X)
  # With random parameters, each simulation draws new coefficients for
  # each unit (the group_var units of newdata, or each row)
  unit <- if (!is.null(object$rpar) && !is.null(object$group_var) &&
              all(object$group_var %in% names(newdata))) {
    as.integer(factor(interaction(newdata[, object$group_var, drop = FALSE],
                                  drop = TRUE, sep = ":")))
  } else seq_len(n)
  pc_draw <- function() {
    if (is.null(object$rpar)) return(.fcr_alloc_predict_core(p, E))
    RS <- object$rpar$layout
    L <- .fcr_alloc_chol(p, RS)
    Bu <- matrix(stats::rnorm(max(unit) * RS$K), max(unit), RS$K) %*% t(L)
    .fcr_alloc_predict_core(p, E, .fcr_alloc_add(Bu[unit, , drop = FALSE],
                                                 E))
  }
  pc <- pc_draw()
  one <- function() {
    if (!is.null(object$rpar)) pc <- pc_draw()
    tot <- .fcr_family_r(object$family, n, pc$mu, pc$alpha, pc$sigma)
    if (object$truncated) {
      for (i in seq_len(100L)) {
        z <- which(tot == 0)
        if (!length(z)) break
        tot[z] <- .fcr_family_r(object$family, length(z), pc$mu[z],
                                pc$alpha[z], pc$sigma[z])
      }
    }
    out <- t(vapply(seq_len(n), function(i)
      as.numeric(stats::rmultinom(1L, tot[i], pc$prob[i, ])),
      numeric(ncol(pc$prob))))
    colnames(out) <- object$categories
    out
  }
  stats::setNames(replicate(nsim, one(), simplify = FALSE),
                  paste0("sim_", seq_len(nsim)))
}
