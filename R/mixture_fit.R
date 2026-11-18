#' Finite-mixture and latent-class count regression
#'
#' @description
#' `countreg_mix()` estimates finite mixtures of count distributions and
#' finite-mixture (latent-class) count regressions by maximum likelihood.
#' `countreg_lc()` is the same model under the latent-class vocabulary
#' (`classes`, `membership_formula`): both functions normalize their
#' arguments into one specification and call one estimation engine, so
#' equivalent calls give the same fit.
#'
#' For unit \eqn{u} (a row, or a panel when `group_var` is given) the
#' probability of the observed counts is
#' \deqn{L_u = \sum_{k=1}^K \pi_{uk} \prod_{i \in u} f_k(y_i \mid x_i,
#' \theta_k),}
#' where \eqn{f_k} is the probability mass function of component \eqn{k}
#' (any family of [countreg()]; see [count_families]) and \eqn{\pi_{uk}} are
#' the mixing (class-membership) probabilities. The log-likelihood is
#' \eqn{\sum_u w_u \log L_u}, with frequency weights \eqn{w_u}.
#'
#' @param formula a model formula for the counts, as in [countreg()]. Use
#'   `y ~ 1` for a mixture of distributions without covariates. Offsets can
#'   be given with `offset()` or the `offset` argument; the same offset
#'   applies to every component.
#' @param data a data frame.
#' @param family a family code (see [count_families]) used by every
#'   component, or one code per component (e.g.
#'   `c(base = "POISSON", tail = "NB2")`): a jointly estimated mixture of
#'   different families. Names become display labels of the components.
#' @param components,classes the number of components (classes) \eqn{K}.
#'   When `family` has several entries, it can be omitted.
#' @param mixing_formula,membership_formula a one-sided formula for the
#'   mixing (class-membership) probabilities, a multinomial logit with the
#'   last component as reference: `~ 1` (default) gives constant
#'   proportions, `~ z` lets them depend on `z` (a mixture of experts). For
#'   panel models, its variables must be constant within each panel.
#' @param fixed_mixing optional known mixing probabilities (one per
#'   component, nonnegative, summing to 1). They are not estimated and add
#'   no parameters. These are *not* frequency weights or starting values.
#'   Zero probabilities are not allowed (the component would be
#'   unidentified).
#' @param shared optional list saying which parameters are the same in all
#'   components: `mean`, a character vector of formula terms (use
#'   `"(Intercept)"` for the intercept; a factor is shared as a whole term),
#'   and `dis1`, `dis2` (`TRUE`/`FALSE`) for the blocks of the first and
#'   second distribution parameter. Distribution parameters can only be
#'   shared between components of the same family. By default nothing is
#'   shared.
#' @param dis_param_formula_1,dis_param_formula_2 optional one-sided formulas
#'   for the first and second distribution parameter of the components, as
#'   in [countreg()]. They apply to every component whose family has the
#'   parameter (coefficients are component-specific unless shared).
#' @param offset,weights an offset (column names, log scale) and frequency
#'   weights, as in [countreg()]. With `group_var`, weights must be constant
#'   within each panel: a weight stands for that many independent copies of
#'   the panel.
#' @param group_var optional column name(s) identifying panels. All rows of a
#'   panel belong to the same (unknown) component, so the panel, not the
#'   row, is the independent unit; rows are independent given the
#'   component.
#' @param method optimizer for the direct maximization of the observed
#'   log-likelihood (see [countreg()]); default `"BFGS"`.
#' @param max.iters maximum number of iterations of the direct maximization.
#' @param start.vals optional starting values: a vector with one value per
#'   parameter (named as in `coef()`, or in that order). It is used as the
#'   first start; the others are generated.
#' @param n_starts number of starts. Starts are reproducible from `seed` and
#'   are recorded in the fitted object (`fit$model$starts`).
#' @param seed optional random seed for the generated starts (and the
#'   bootstrap). The caller's random-number stream is not changed. When
#'   `NULL`, a seed is drawn from the current stream and stored.
#' @param stderr `"normal"` (from the Hessian of the observed-data
#'   log-likelihood), `"robust"` (sandwich with unit-level scores) or
#'   `"boot"` (nonparametric bootstrap of units, refitting the full mixture
#'   with several starts and aligning component labels).
#' @param bootstraps number of bootstrap replicates for `stderr = "boot"`
#'   (default 200).
#' @param ndraws,engine numerical settings of the mixed-Poisson families
#'   (PLN, PGE, PW, PLL) within each component, as in [countreg()]. They
#'   are unrelated to the number of components.
#' @param keep_data logical; store the training data in the fitted object.
#' @param verbose logical; report the progress of each start.
#' @param control a [mixture_control()] object.
#'
#' @details
#' ## Estimation
#' Each start is a response-based or random preliminary partition of the
#' units, turned into parameters by one M-step. Starts are improved by the
#' EM algorithm (posterior class probabilities in the E-step; weighted
#' updates of the component and mixing parameters in the M-step, with
#' shared parameters updated jointly), and the best distinct solutions are
#' polished by maximizing the observed-data log-likelihood directly with
#' `maxLik`. The selected solution is the best polished one that converged.
#' The start table (`fit$model$starts`) records every start's strategy,
#' log-likelihood, EM iterations and status, and whether the best solution
#' was reached by several starts; mixtures often have several local maxima,
#' so a solution found by one start only deserves caution.
#'
#' Components of the same family are exchangeable: their numbering is
#' arbitrary. After estimation they are numbered in increasing order of
#' their average expected count over the training rows (a display
#' convention). Mixtures of different families keep the order given in
#' `family`.
#'
#' ## Parameters
#' The coefficients are named `Component<k>:Mean:<term>`,
#' `Component<k>:<parameter>` (e.g. `Component2:ln(alpha)`),
#' `Shared:...` for shared parameters, and `Mixing:Component<k>:<term>` for
#' the membership logit of component \eqn{k} against the last component.
#' Component expected counts use each family's own mean (for PLN and PLL
#' it includes the lognormal factor), so shared link-scale coefficients do
#' not imply equal component means.
#'
#' ## Inference and model choice
#' The covariance matrix is the inverse of the observed information of the
#' whole free parameter vector (not of the separate M-steps). Near a
#' boundary (an almost empty component, or two almost identical
#' components) it is unreliable; use the bootstrap and the start table.
#' The usual chi-square likelihood-ratio test does **not** apply to
#' comparing numbers of components (K versus K + 1); compare such models
#' with AIC/BIC on the same sample, start stability and held-out
#' log scores. [regCompTest()] refuses mixture models for this reason.
#'
#' ## Interpretation
#' Classes are a modeling device. A good fit with two components does not
#' show that two real subpopulations exist; the same mixture can describe a
#' skewed or heavy-tailed distribution.
#'
#' Zero-inflated, hurdle and underreporting components, component-specific
#' formulas, and confidence/prediction/tolerance intervals for mixture
#' predictions are not available yet.
#'
#' @returns An object of class `c("countreg_mix", "flexCountReg")`
#'   (`countreg_mix()`) or `c("countreg_lc", "countreg_mix",
#'   "flexCountReg")` (`countreg_lc()`). The methods `coef()`, `vcov()`,
#'   `logLik()`, `AIC()`, `BIC()`, `summary()`, `predict()`
#'   ([predict.countreg_mix()]), `fitted()`, `residuals()`, `simulate()`,
#'   `tidy()`, `glance()` and `augment()` work as for other models, and
#'   [component_probs()] returns prior or posterior class probabilities.
#'   [rootogram()], [cureplot()] and [margEffTable()] (delta-method standard
#'   errors) accept mixture models.
#'
#' @references
#' McLachlan, G. J., & Peel, D. (2000). *Finite Mixture Models*. Wiley.
#'
#' Grün, B., & Leisch, F. (2008). FlexMix version 2: Finite mixtures with
#' concomitant variables and varying and constant parameters. *Journal of
#' Statistical Software*, 28(4), 1-35.
#'
#' Wedel, M., DeSarbo, W. S., Bult, J. R., & Ramaswamy, V. (1993). A latent
#' class Poisson regression model for heterogeneous count data. *Journal of
#' Applied Econometrics*, 8(4), 397-411.
#'
#' @seealso [mixture_control()], [component_probs()],
#'   [predict.countreg_mix()], [count_families]
#'
#' @examples
#' # A two-component Poisson mixture of simulated counts
#' set.seed(1)
#' d <- data.frame(y = c(rpois(300, 1), rpois(200, 6)))
#' m <- countreg_mix(y ~ 1, data = d, family = "POISSON", components = 2,
#'                   n_starts = 5, seed = 42)
#' summary(m)
#' component_probs(m, conditioning = "prior")[1, ]
#'
#' \donttest{
#' # The same model through the latent-class interface
#' lc <- countreg_lc(y ~ 1, data = d, family = "POISSON", classes = 2,
#'                   n_starts = 5, seed = 42)
#' all.equal(logLik(m), logLik(lc))
#'
#' # Latent-class NB2 regression with a membership covariate
#' data("abalone_rings")
#' fit <- countreg_lc(rings ~ shell_weight, data = abalone_rings,
#'                    family = "NB2", classes = 2,
#'                    membership_formula = ~ sex, n_starts = 5, seed = 1)
#' summary(fit)
#' }
#' @name countreg_mix
#' @export
countreg_mix <- function(formula, data, family = "POISSON", components = 2L,
                         mixing_formula = ~ 1, fixed_mixing = NULL,
                         shared = NULL, dis_param_formula_1 = NULL,
                         dis_param_formula_2 = NULL, offset = NULL,
                         weights = NULL, group_var = NULL, method = "BFGS",
                         max.iters = 1000, start.vals = NULL, n_starts = 10L,
                         seed = NULL, stderr = c("normal", "robust", "boot"),
                         bootstraps = NULL, ndraws = 1500,
                         engine = "quadrature", keep_data = TRUE,
                         verbose = FALSE, control = mixture_control()) {
  .fcr_fm_fit(
    entry = "countreg_mix", call = match.call(), formula = formula,
    data = data, family = family, components = components,
    components_missing = missing(components), mixing_formula = mixing_formula,
    fixed_mixing = fixed_mixing, shared = shared,
    dis_param_formula_1 = dis_param_formula_1,
    dis_param_formula_2 = dis_param_formula_2, offset = offset,
    weights = weights, group_var = group_var, method = method,
    max.iters = max.iters, start.vals = start.vals, n_starts = n_starts,
    seed = seed, stderr = stderr, bootstraps = bootstraps, ndraws = ndraws,
    engine = engine, keep_data = keep_data, verbose = verbose,
    control = control)
}

#' @rdname countreg_mix
#' @export
countreg_lc <- function(formula, data, family = "POISSON", classes = 2L,
                        membership_formula = ~ 1, fixed_mixing = NULL,
                        shared = NULL, dis_param_formula_1 = NULL,
                        dis_param_formula_2 = NULL, offset = NULL,
                        weights = NULL, group_var = NULL, method = "BFGS",
                        max.iters = 1000, start.vals = NULL, n_starts = 10L,
                        seed = NULL, stderr = c("normal", "robust", "boot"),
                        bootstraps = NULL, ndraws = 1500,
                        engine = "quadrature", keep_data = TRUE,
                        verbose = FALSE, control = lc_control()) {
  .fcr_fm_fit(
    entry = "countreg_lc", call = match.call(), formula = formula,
    data = data, family = family, components = classes,
    components_missing = missing(classes),
    mixing_formula = membership_formula, fixed_mixing = fixed_mixing,
    shared = shared, dis_param_formula_1 = dis_param_formula_1,
    dis_param_formula_2 = dis_param_formula_2, offset = offset,
    weights = weights, group_var = group_var, method = method,
    max.iters = max.iters, start.vals = start.vals, n_starts = n_starts,
    seed = seed, stderr = stderr, bootstraps = bootstraps, ndraws = ndraws,
    engine = engine, keep_data = keep_data, verbose = verbose,
    control = control)
}

# Evaluate expr with a fixed seed without changing the caller's
# random-number stream.
.fcr_with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old)) {
    if (exists(".Random.seed", envir = globalenv()))
      rm(".Random.seed", envir = globalenv())
  } else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  expr
}

# ---- Specification ----------------------------------------------------------

# Validate the arguments into the settings of a fit (everything needed to
# rebuild the training problem from data, e.g. in bootstrap refits).
.fcr_fm_settings <- function(a) {
  lab_prefix <- if (a$entry == "countreg_lc") "Class" else "Component"
  fam <- .fcr_fm_families(a$family, a$components, a$components_missing,
                          lab_prefix)
  K <- fam$K
  check_formula <- function(f, nm) {
    if (is.null(f)) return(NULL)
    f <- tryCatch(stats::as.formula(f), error = function(e) NULL)
    if (is.null(f) || length(f) != 2L) {
      stop("`", nm, "` must be a one-sided formula such as ~ x.",
           call. = FALSE)
    }
    f
  }
  mix_nm <- if (a$entry == "countreg_lc") "membership_formula" else
    "mixing_formula"
  mixing_formula <- check_formula(a$mixing_formula %||% ~ 1, mix_nm)
  fixed <- a$fixed_mixing
  if (!is.null(fixed)) {
    if (length(all.vars(mixing_formula)) > 0L) {
      stop("Give either `fixed_mixing` or a `", mix_nm, "` with ",
           "covariates, not both.", call. = FALSE)
    }
    if (!is.numeric(fixed) || length(fixed) != K || any(!is.finite(fixed)) ||
        any(fixed < 0)) {
      stop("`fixed_mixing` must be ", K, " finite nonnegative ",
           "probabilities, one per component.", call. = FALSE)
    }
    if (abs(sum(fixed) - 1) > 1e-8) {
      stop("`fixed_mixing` must sum to 1 (it sums to ", format(sum(fixed)),
           "); it is not rescaled.", call. = FALSE)
    }
    if (any(fixed == 0)) {
      stop("`fixed_mixing` has a zero probability: that component would ",
           "have no data and unidentified parameters. Remove the ",
           "component instead.", call. = FALSE)
    }
    fixed <- as.numeric(fixed / sum(fixed))
  }
  mix_mode <- if (K == 1L) "none" else if (!is.null(fixed)) "fixed" else
    if (length(all.vars(mixing_formula))) "covariate" else "constant"
  params <- lapply(fam$families, get_params)
  f1 <- check_formula(a$dis_param_formula_1, "dis_param_formula_1")
  f2 <- check_formula(a$dis_param_formula_2, "dis_param_formula_2")
  if (!is.null(f1) && all(vapply(params, function(p) is.null(p[[1]]),
                                 logical(1)))) {
    stop("`dis_param_formula_1` was given, but no component family has a ",
         "distribution parameter to model.", call. = FALSE)
  }
  if (!is.null(f2) && all(vapply(params, function(p) is.null(p[[2]]),
                                 logical(1)))) {
    stop("`dis_param_formula_2` was given, but no component family has a ",
         "second distribution parameter.", call. = FALSE)
  }
  engine <- match.arg(a$engine, c("quadrature", "halton", "poilog"))
  if (engine == "poilog" && "PLN" %in% fam$families) {
    rlang::check_installed("poilog", reason = "for `engine = \"poilog\"`.")
  }
  .fcr_check_controls(a$max.iters, a$ndraws)
  if (!inherits(a$control, "fcr_mixture_control")) {
    stop("`control` must be made by mixture_control() or lc_control().",
         call. = FALSE)
  }
  n_starts <- a$n_starts
  if (!is.numeric(n_starts) || length(n_starts) != 1L ||
      !is.finite(n_starts) || n_starts < 1 || n_starts != round(n_starts)) {
    stop("`n_starts` must be a single whole number of at least 1.",
         call. = FALSE)
  }
  if (!is.null(a$seed) && (!is.numeric(a$seed) || length(a$seed) != 1L ||
                           !is.finite(a$seed))) {
    stop("`seed` must be a single number.", call. = FALSE)
  }
  if (!is.null(a$group_var) && (!is.character(a$group_var) ||
                                !length(a$group_var))) {
    stop("`group_var` must be column name(s).", call. = FALSE)
  }
  list(entry = a$entry, formula = stats::as.formula(a$formula),
       families = fam$families, labels = fam$labels, K = K,
       mixing_formula = mixing_formula, mix_mode = mix_mode,
       fixed_mixing = fixed, shared = .fcr_fm_check_shared(a$shared),
       dis1 = f1, dis2 = f2, offset = a$offset, weights = a$weights,
       group_var = a$group_var, ndraws = a$ndraws, engine = engine,
       method = .fcr_check_method(a$method), max.iters = a$max.iters,
       n_starts = as.integer(n_starts), control = a$control)
}

# Build the training problem from the settings and data.
.fcr_fm_setup <- function(st, data) {
  K <- st$K
  aux <- list(dis1 = st$dis1, dis2 = st$dis2,
              mixing = if (st$mix_mode %in% c("constant", "covariate"))
                st$mixing_formula)
  prep <- prepare_model_data(st$formula, data, aux_formulas = aux,
                             offset = st$offset, weights = st$weights,
                             group_var = st$group_var)
  n <- length(prep$y)
  ones <- matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))
  params <- lapply(st$families, get_params)
  has1 <- any(vapply(params, function(p) !is.null(p[[1]]), logical(1)))
  has2 <- any(vapply(params, function(p) !is.null(p[[2]]), logical(1)))
  des <- list(X = prep$X, offset = prep$offset,
              D1 = if (has1) prep$aux$dis1 %||% ones,
              D2 = if (has2) prep$aux$dis2 %||% ones,
              Z = if (st$mix_mode %in% c("constant", "covariate"))
                prep$aux$mixing %||% ones,
              n = n)
  w <- prep$weights
  unit <- NULL
  n_u <- n
  Zu <- des$Z
  w_unit <- w
  unit_ids <- NULL
  if (!is.null(st$group_var)) {
    unit <- as.integer(prep$group)
    unit_ids <- levels(prep$group)
    n_u <- length(unit_ids)
    first <- match(seq_len(n_u), unit)
    w_unit <- w[first]
    if (any(abs(w - w_unit[unit]) > 1e-10 * pmax(1, abs(w)))) {
      stop("With `group_var`, `weights` must be constant within each ",
           "panel (a weight is the number of copies of the whole panel).",
           call. = FALSE)
    }
    if (!is.null(des$Z)) {
      Zu <- des$Z[first, , drop = FALSE]
      dif <- abs(des$Z - Zu[unit, , drop = FALSE])
      bad <- colnames(des$Z)[colSums(dif > 1e-10 * pmax(1, abs(des$Z))) > 0]
      if (length(bad)) {
        stop("Mixing (membership) variables must be constant within each ",
             "panel, because the whole panel belongs to one class. ",
             "Not constant: ", toString(bad), ". Use a panel-level summary ",
             "(e.g. the panel mean) as a new column if that is what you ",
             "intend.", call. = FALSE)
      }
    }
  }
  if (!is.null(Zu)) .fcr_check_design(Zu, "mixing_formula")
  map <- .fcr_fm_map(st$families, des$X, prep$blueprint$terms, des$D1,
                     des$D2, !is.null(st$dis1), !is.null(st$dis2), Zu,
                     st$shared)
  haltons <- randtoolbox::halton(st$ndraws)
  list(st = st, prep = prep, des = des, y = as.numeric(prep$y), n = n,
       unit = unit, unit_ids = unit_ids, n_u = n_u, Zu = Zu, w_row = w,
       w_unit = w_unit, K = K, map = map, mix_mode = st$mix_mode,
       fixed_mixing = st$fixed_mixing, haltons = haltons,
       normed_haltons = stats::qnorm(haltons), engine = st$engine,
       analytic = .fcr_fm_analytic(st$families),
       blueprint = prep$blueprint)
}

# ---- Starting values ----------------------------------------------------------

# Pooled one-component starting values: the Poisson regression for the mean
# coefficients of every component, the family's starting value (or the
# moment estimate of the dispersion) for the distribution parameters, and
# equal mixing probabilities.
.fcr_fm_pooled <- function(spec) {
  map <- spec$map
  pm <- suppressWarnings(stats::glm.fit(spec$des$X, spec$y,
                                        weights = spec$w_row,
                                        offset = spec$des$offset,
                                        family = stats::poisson()))
  b <- pm$coefficients
  b[is.na(b)] <- 0
  theta <- stats::setNames(numeric(map$n_par), map$names)
  for (k in seq_len(spec$K)) {
    theta[map$mean_idx[[k]]] <- b
    fam <- map$families[k]
    for (s in c("d1", "d2")) {
      idx <- map[[paste0(s, "_idx")]][[k]]
      if (!length(idx)) next
      D <- if (s == "d1") spec$des$D1 else spec$des$D2
      st0 <- rep_len(.fcr_get_family(fam)$start, 1)
      if (s == "d1") {
        a0 <- .fcr_moment_disp_start(fam, spec$y, pm$fitted.values,
                                     spec$w_row)
        if (!is.null(a0)) st0 <- a0
      }
      v <- numeric(length(idx))
      int <- which(colnames(D) == "(Intercept)")
      v[if (length(int)) int else 1L] <- st0
      theta[idx] <- v
    }
  }
  list(theta = theta, mu = pm$fitted.values)
}

# Starting responsibilities (n_u x K) from a preliminary partition of the
# units: "quantile" orders the units by their observed-to-pooled-mean ratio
# and splits them into K groups; "jitter" adds noise to that score first;
# "random" assigns units at random. Each unit gets 0.8 on its group.
.fcr_fm_start_tau <- function(spec, mu0, strategy, spread) {
  K <- spec$K
  n_u <- spec$n_u
  if (K == 1L) return(matrix(1, n_u, 1))
  num <- spec$w_row * spec$y
  den <- spec$w_row * mu0
  if (!is.null(spec$unit)) {
    num <- rowsum(num, spec$unit, reorder = TRUE)[, 1]
    den <- rowsum(den, spec$unit, reorder = TRUE)[, 1]
  }
  score <- log((num + 0.5) / (den + 0.5))
  g <- if (strategy == "random") {
    sample.int(K, n_u, replace = TRUE)
  } else {
    if (strategy == "jitter") {
      s <- stats::sd(score)
      score <- score + stats::rnorm(n_u, 0, spread * (if (is.finite(s) &&
                                                         s > 0) s else 1))
    }
    r <- rank(score, ties.method = "random")
    pmin(K, ceiling(K * r / n_u))
  }
  tau <- matrix(0.2 / (K - 1), n_u, K)
  tau[cbind(seq_len(n_u), g)] <- 0.8
  tau
}

# Parameters from starting responsibilities: one mixing step and one
# component step (with more iterations) from the pooled values.
.fcr_fm_theta_from_tau <- function(spec, theta0, tau) {
  th <- .fcr_fm_mstep_mix(theta0, spec, tau, 100L)
  if (is.null(th)) th <- theta0
  .fcr_fm_mstep_comp(th, spec, tau, 100L)
}

# ---- EM and direct maximization ---------------------------------------------------

.fcr_fm_total <- function(e, spec) {
  if (is.null(e)) return(-Inf)
  pos <- spec$w_unit > 0
  sum(spec$w_unit[pos] * e$ll[pos])
}

.fcr_fm_em <- function(theta, spec, control) {
  e <- .fcr_fm_eval(theta, spec)
  ll <- .fcr_fm_total(e, spec)
  if (!is.finite(ll)) {
    return(list(theta = theta, loglik = -Inf, iterations = 0L,
                converged = FALSE, status = "invalid start"))
  }
  status <- "iteration limit"
  conv <- FALSE
  it <- 0L
  while (it < control$em_maxit) {
    it <- it + 1L
    th <- .fcr_fm_mstep_mix(theta, spec, e$tau, max(50L, control$mstep_maxit))
    if (is.null(th)) {
      status <- "empty component"
      break
    }
    th <- .fcr_fm_mstep_comp(th, spec, e$tau, control$mstep_maxit)
    e2 <- .fcr_fm_eval(th, spec)
    ll2 <- .fcr_fm_total(e2, spec)
    if (!is.finite(ll2) || ll2 < ll - 1e-8 * (abs(ll) + 1)) {
      status <- "likelihood decreased"
      break
    }
    dpar <- max(abs(th - theta))
    dll <- ll2 - ll
    theta <- th
    e <- e2
    ll <- ll2
    if (abs(dll) <= control$em_tol * (abs(ll) + control$em_tol) &&
        dpar <= sqrt(control$em_tol)) {
      conv <- TRUE
      status <- "converged"
      break
    }
  }
  list(theta = theta, loglik = ll, iterations = it, converged = conv,
       status = status)
}

# Direct maximization of the observed log-likelihood with maxLik (unit-level
# contributions, so robust standard errors and BHHH work).
.fcr_fm_maxlik <- function(theta, spec, method, max.iters, print.level = 0) {
  fn <- function(p) .fcr_fm_loglik_units(p, spec)
  gr <- if (spec$analytic) function(p) .fcr_fm_grad_units(p, spec)
  fit <- tryCatch(
    maxLik::maxLik(fn, grad = gr, start = theta, method = method,
                   control = list(iterlim = max.iters,
                                  printLevel = print.level)),
    error = function(e) NULL)
  if (is.null(fit) || !is.finite(fit$maximum)) return(NULL)
  names(fit$estimate) <- spec$map$names
  fit$converged <- .fcr_converged(fit, method)
  fit
}

# Run all starts, polish the best distinct ones and select the solution.
# Returns the selected maxLik fit and the start table.
.fcr_fm_search <- function(spec, n_starts, seed, start.vals = NULL,
                           verbose = FALSE, polish_print = 0) {
  st <- spec$st
  control <- st$control
  map <- spec$map
  pooled <- .fcr_fm_pooled(spec)
  if (spec$K == 1L) n_starts <- 1L
  seeds <- .fcr_with_seed(seed, sample.int(.Machine$integer.max, n_starts))
  strategies <- character(n_starts)
  starts <- vector("list", n_starts)
  for (s in seq_len(n_starts)) {
    if (s == 1L && !is.null(start.vals)) {
      strategies[s] <- "user"
      starts[[s]] <- start.vals
      next
    }
    strategies[s] <- if (spec$K == 1L) "pooled" else
      if (s == 1L || (s == 2L && !is.null(start.vals))) "quantile" else
        c("jitter", "random")[(s %% 2L) + 1L]
    if (spec$K == 1L) {
      starts[[s]] <- pooled$theta
      next
    }
    tau <- .fcr_with_seed(seeds[s], .fcr_fm_start_tau(
      spec, pooled$mu, strategies[s], control$start_spread))
    starts[[s]] <- .fcr_fm_theta_from_tau(spec, pooled$theta, tau)
  }
  direct <- control$algorithm == "direct" || spec$K == 1L
  res <- vector("list", n_starts)
  for (s in seq_len(n_starts)) {
    if (direct) {
      f <- .fcr_fm_maxlik(starts[[s]], spec, st$method, st$max.iters)
      res[[s]] <- list(theta = if (is.null(f)) starts[[s]] else f$estimate,
                       loglik = if (is.null(f)) -Inf else f$maximum,
                       iterations = 0L, converged = NA,
                       status = if (is.null(f)) "failed" else "direct",
                       fit = f)
    } else {
      res[[s]] <- .fcr_fm_em(starts[[s]], spec, control)
    }
    if (verbose) {
      message(sprintf("Start %d/%d (%s): log-likelihood %.4f, %s", s,
                      n_starts, strategies[s], res[[s]]$loglik,
                      res[[s]]$status))
    }
  }
  ll_start <- vapply(res, `[[`, numeric(1), "loglik")
  # Polish the best distinct solutions
  ord <- order(ll_start, decreasing = TRUE)
  ord <- ord[is.finite(ll_start[ord])]
  if (!length(ord)) {
    stop("No start gave a finite log-likelihood. Check the data and the ",
         "model specification.", call. = FALSE)
  }
  chosen <- integer(0)
  for (s in ord) {
    if (length(chosen) >= control$n_polish) break
    if (any(abs(ll_start[chosen] - ll_start[s]) <=
            1e-6 * (abs(ll_start[s]) + 1))) next
    chosen <- c(chosen, s)
  }
  fits <- vector("list", n_starts)
  final_ll <- rep(NA_real_, n_starts)
  conv <- rep(NA, n_starts)
  for (s in chosen) {
    f <- if (direct) res[[s]]$fit else
      .fcr_fm_maxlik(res[[s]]$theta, spec, st$method, st$max.iters,
                     polish_print)
    fits[[s]] <- f
    if (!is.null(f)) {
      final_ll[s] <- f$maximum
      conv[s] <- isTRUE(f$converged)
    }
  }
  ok <- which(!is.na(final_ll))
  if (!length(ok)) {
    stop("Direct maximization failed for every candidate solution.",
         call. = FALSE)
  }
  good <- ok[conv[ok] %in% TRUE]
  pick_from <- if (length(good)) good else ok
  best <- pick_from[which.max(final_ll[pick_from])]
  # With analytic scores, a few Newton-Raphson steps (numerical Hessian of
  # the analytic gradient) tighten the gradient after a quasi-Newton
  # optimizer; kept only if the log-likelihood does not decrease.
  if (spec$analytic && st$method != "NR") {
    f2 <- .fcr_fm_maxlik(fits[[best]]$estimate, spec, "NR", 20L)
    if (!is.null(f2) && isTRUE(f2$converged) &&
        f2$maximum >= fits[[best]]$maximum - 1e-8) {
      f2$refined_from <- st$method
      fits[[best]] <- f2
      final_ll[best] <- f2$maximum
      conv[best] <- TRUE
    }
  }
  best_ll <- final_ll[best]
  tol <- max(1e-3, 1e-6 * abs(best_ll))
  reached <- pmax(ll_start, final_ll, na.rm = TRUE)
  table <- data.frame(
    start = seq_len(n_starts), strategy = strategies, seed = seeds,
    em_iterations = vapply(res, function(r) as.integer(r$iterations),
                           integer(1)),
    em_status = vapply(res, `[[`, character(1), "status"),
    loglik_start = ll_start, polished = seq_len(n_starts) %in% chosen,
    loglik = final_ll, converged = conv,
    reached_best = is.finite(reached) & reached >= best_ll - tol,
    selected = seq_len(n_starts) == best, stringsAsFactors = FALSE)
  list(fit = fits[[best]], table = table, any_converged = length(good) > 0,
       better_failed = any(reached[-best] > best_ll + tol, na.rm = TRUE))
}

# Average expected count of each component over the training rows.
.fcr_fm_comp_avg_mean <- function(theta, spec) {
  pars <- .fcr_fm_unpack(theta, spec$map)
  w <- spec$w_row
  vapply(seq_len(spec$K), function(k) {
    cv <- .fcr_fm_comp_values(spec$des, pars, k)
    m <- suppressWarnings(.fcr_family_mean(spec$map$families[k], cv$mu,
                                           cv$alpha, cv$sigma))
    sum(w * m, na.rm = TRUE) / sum(w)
  }, numeric(1))
}

# Display order of exchangeable components: increasing average expected
# count. Components of different families, or with different fixed
# probabilities, keep their order.
.fcr_fm_canonical_perm <- function(theta, spec) {
  K <- spec$K
  if (K == 1L || length(unique(spec$map$families)) > 1L ||
      (!is.null(spec$fixed_mixing) &&
       length(unique(spec$fixed_mixing)) > 1L)) {
    return(seq_len(K))
  }
  order(.fcr_fm_comp_avg_mean(theta, spec))
}

# ---- Main fitting function -------------------------------------------------------

.fcr_fm_fit <- function(entry, call, formula, data, family, components,
                        components_missing, mixing_formula, fixed_mixing,
                        shared, dis_param_formula_1, dis_param_formula_2,
                        offset, weights, group_var, method, max.iters,
                        start.vals, n_starts, seed, stderr, bootstraps,
                        ndraws, engine, keep_data, verbose, control) {
  stderr <- match.arg(tolower(stderr[1]), c("normal", "robust", "boot"))
  if (stderr == "boot") {
    if (is.null(bootstraps)) bootstraps <- 200L
    if (!is.numeric(bootstraps) || length(bootstraps) != 1L ||
        !is.finite(bootstraps) || bootstraps < 2) {
      stop("`bootstraps` must be a single number of at least 2.",
           call. = FALSE)
    }
    bootstraps <- as.integer(bootstraps)
  }
  st <- .fcr_fm_settings(list(
    entry = entry, formula = formula, family = family,
    components = components, components_missing = components_missing,
    mixing_formula = mixing_formula, fixed_mixing = fixed_mixing,
    shared = shared, dis_param_formula_1 = dis_param_formula_1,
    dis_param_formula_2 = dis_param_formula_2, offset = offset,
    weights = weights, group_var = group_var, method = method,
    max.iters = max.iters, n_starts = n_starts, seed = seed,
    ndraws = ndraws, engine = engine, control = control))
  data <- as.data.frame(data)
  spec <- .fcr_fm_setup(st, data)
  prep <- spec$prep
  # Keep the formulas with `.` expanded
  st$formula <- prep$formula
  st$dis1 <- prep$aux_formulas$dis1
  st$dis2 <- prep$aux_formulas$dis2
  if (!is.null(prep$aux_formulas$mixing)) {
    st$mixing_formula <- prep$aux_formulas$mixing
  }
  spec$st <- st
  if (prep$n_dropped > 0) {
    message(prep$n_dropped, " row(s) with missing values were removed.")
  }
  status <- .fcr_fm_family_status(unique(st$families))
  if (any(status == "experimental")) {
    message("Mixtures of family ",
            toString(unique(st$families)[status == "experimental"]),
            " are experimental: the component distributions are those of ",
            "countreg(), but recovery of such mixtures has not been ",
            "validated. Check the start table and the stability of the fit.")
  }
  if (stderr == "boot" && !is.null(weights)) {
    .fcr_check_boot_freq(spec$w_unit, length(spec$w_unit))
  }
  if (!is.null(start.vals)) start.vals <- .fcr_check_start(start.vals,
                                                           spec$map$names)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)

  # ---- Search ----
  search <- .fcr_fm_search(spec, st$n_starts, seed, start.vals, verbose,
                           if (verbose) 1 else 0)
  fit <- search$fit
  if (!search$any_converged) {
    warning("The optimizer did not report convergence for any start (",
            trimws(fit$message), "). Consider more `max.iters`, more ",
            "`n_starts` or a different `method`.", call. = FALSE)
  }
  if (search$better_failed) {
    warning("A start that did not converge reached a higher ",
            "log-likelihood than the selected solution; see ",
            "fit$model$starts.", call. = FALSE)
  }

  # ---- Display order of exchangeable components ----
  theta <- fit$estimate
  perm <- .fcr_fm_canonical_perm(theta, spec)
  H <- fit$hessian
  gO <- fit$gradientObs
  g <- fit$gradient
  if (!identical(perm, seq_len(spec$K))) {
    A <- .fcr_fm_perm_matrix(spec$map, perm)
    Ainv <- solve(A)
    theta <- as.vector(A %*% theta)
    if (!is.null(H)) H <- t(Ainv) %*% H %*% Ainv
    if (!is.null(gO)) gO <- gO %*% Ainv
    if (!is.null(g)) g <- as.vector(g %*% Ainv)
  }
  names(theta) <- spec$map$names
  if (!is.null(H)) dimnames(H) <- list(spec$map$names, spec$map$names)
  .fcr_fm_finish(spec, theta, H, g, gO, fit, search, stderr, bootstraps,
                 seed, data, call, keep_data, entry)
}

.fcr_fm_finish <- function(spec, theta, H, g, gO, fit, search, stderr,
                           bootstraps, seed, data, call, keep_data, entry) {
  st <- spec$st
  map <- spec$map
  nm <- map$names
  # ---- Covariance ----
  boot <- NULL
  if (stderr == "boot") {
    boot <- .fcr_fm_bootstrap(spec, theta, H, data, bootstraps, seed)
    V <- .fcr_boot_cov(boot$draws)
    vcov_type <- "bootstrap"
  } else {
    V <- .fcr_hessian_vcov(list(hessian = H, estimate = theta))
    vcov_type <- "hessian"
    if (stderr == "robust") {
      keep <- spec$w_unit > 0
      meat <- crossprod(gO[keep, , drop = FALSE] / sqrt(spec$w_unit[keep]))
      V <- V %*% meat %*% V
      vcov_type <- "robust"
    }
  }
  dimnames(V) <- list(nm, nm)

  # ---- Posterior summaries ----
  e <- .fcr_fm_eval(theta, spec)
  tau <- e$tau
  prior <- exp(e$logpi)
  colnames(tau) <- colnames(prior) <- st$labels
  w <- spec$w_unit
  ent <- -sum(w * rowSums(ifelse(tau > 0, tau * log(tau), 0)))
  K <- spec$K
  n_bic <- sum(w)
  n_par <- length(theta)
  bic <- -2 * fit$maximum + log(n_bic) * n_par
  modal <- max.col(tau, ties.method = "first")
  classes <- data.frame(
    component = seq_len(K), label = st$labels, family = st$families,
    prior_share = colSums(w * prior) / n_bic,
    posterior_share = colSums(w * tau) / n_bic,
    effective_size = colSums(w * tau),
    modal_count = as.vector(tapply(w, factor(modal, levels = seq_len(K)),
                                   sum)),
    mean_count = .fcr_fm_comp_avg_mean(theta, spec),
    stringsAsFactors = FALSE)
  classes$modal_count[is.na(classes$modal_count)] <- 0
  classes$small <- classes$posterior_share < st$control$small_class

  # ---- Numerical checks ----
  numerical <- .fcr_fm_numerical(H, g, fit$maximum)
  identical_pairs <- .fcr_fm_identical(theta, spec)
  if (length(identical_pairs)) {
    warning("Components ", paste(identical_pairs, collapse = ", "),
            " are nearly identical (same family and almost the same ",
            "parameters): the fit is effectively one component fewer, and ",
            "their shares and means are arbitrary. Use fewer components.",
            call. = FALSE)
  }

  out <- list(
    estimate = theta, maximum = fit$maximum, hessian = H, gradient = g,
    gradientObs = gO, code = fit$code, message = fit$message,
    iterations = fit$iterations, method = st$method,
    converged = isTRUE(fit$converged), vcov = V, vcov_type = vcov_type,
    stderr = stderr, se = stats::setNames(.fcr_se(V), nm),
    bootstraps = if (stderr == "boot") bootstraps,
    boot_failures = boot$failures, n_boot_ok = boot$n_ok,
    n_boot_ambiguous = boot$n_ambiguous,
    modelType = entry, entry = entry,
    family = st$families, K = K, labels = st$labels,
    component_ids = seq_len(K),
    family_status = .fcr_fm_family_status(st$families),
    formula = st$formula, dis_param_formula_1 = st$dis1,
    dis_param_formula_2 = st$dis2, mixing_formula = st$mixing_formula,
    mix_mode = st$mix_mode, fixed_mixing = st$fixed_mixing,
    shared = st$shared, map = map, settings = st,
    blueprint = spec$blueprint, offset = spec$blueprint$offset_label,
    group_var = st$group_var, unit = spec$unit, unit_ids = spec$unit_ids,
    y = spec$y, weights = spec$w_row, unit_weights = spec$w_unit,
    nobs = spec$n, n_bic = n_bic,
    n_groups = if (!is.null(spec$unit)) spec$n_u,
    n_units = spec$n_u, rows_used = spec$prep$rows_used,
    n_dropped = spec$prep$n_dropped,
    posterior = tau, prior = prior, classes = classes, entropy = ent,
    relative_entropy = if (K > 1L) 1 - ent / (n_bic * log(K)) else NA_real_,
    ICL = bic + 2 * ent, starts = search$table, seed = seed,
    refined_from = fit$refined_from, identical_pairs = identical_pairs,
    n_starts = nrow(search$table), control = st$control,
    ndraws = st$ndraws, engine = st$engine, numerical = numerical)
  .createFlexCountReg(model = out, data = data, call = call,
                      keep_data = keep_data)
}

# Pairs of components ("1-2") of the same family whose component-specific
# parameters differ by less than `tol` (a degenerate fit: the mixture is
# effectively one component fewer).
.fcr_fm_identical <- function(theta, spec, tol = 1e-3) {
  map <- spec$map
  K <- spec$K
  out <- character(0)
  if (K < 2L) return(out)
  for (a in seq_len(K - 1L)) for (b in (a + 1L):K) {
    if (map$families[a] != map$families[b]) next
    ia <- map$comp_specific[[a]]
    ib <- map$comp_specific[[b]]
    if (length(ia) != length(ib)) next
    if (!length(ia) || max(abs(theta[ia] - theta[ib])) < tol) {
      out <- c(out, paste0(a, "-", b))
    }
  }
  out
}

# Numerical checks at the solution: gradient, definiteness and conditioning
# of the observed information (-H).
.fcr_fm_numerical <- function(H, g, ll) {
  out <- list(max_abs_gradient = if (!is.null(g)) max(abs(g)) else NA_real_,
              min_eigenvalue = NA_real_, condition = NA_real_,
              information = "unavailable")
  if (is.null(H) || any(!is.finite(H))) return(out)
  ev <- tryCatch(eigen(-(H + t(H)) / 2, symmetric = TRUE,
                       only.values = TRUE)$values, error = function(e) NULL)
  if (is.null(ev)) return(out)
  out$min_eigenvalue <- min(ev)
  tol <- sqrt(.Machine$double.eps) * max(abs(ev))
  out$information <- if (min(ev) > tol) "positive definite" else
    if (min(ev) >= -tol) "singular" else "indefinite"
  if (out$information == "positive definite") {
    out$condition <- max(ev) / min(ev)
  }
  out
}

# ---- Bootstrap ------------------------------------------------------------------------

# Nonparametric bootstrap of the independent units (rows, or whole panels
# with fresh identifiers), refitting the full mixture with the original
# estimates plus `boot_starts` generated starts, and aligning component
# labels to the original fit (among permutations that keep each family).
.fcr_fm_bootstrap <- function(spec, theta, H, data, B, seed) {
  st <- spec$st
  map <- spec$map
  d_used <- data[spec$prep$rows_used, , drop = FALSE]
  perms <- .fcr_fm_allowed_perms(map$families, spec$fixed_mixing)
  V0 <- tryCatch(solve(-H), error = function(e) NULL)
  scale <- if (!is.null(V0)) sqrt(pmax(diag(V0), 0)) else rep(1, length(theta))
  scale[!is.finite(scale) | scale <= 0] <- 1
  n_amb <- 0L
  unit_rows <- if (is.null(spec$unit)) as.list(seq_len(spec$n)) else
    split(seq_len(spec$n), spec$unit)
  gname <- .fcr_unused_name(names(d_used), ".fcr_boot_unit")
  wname <- .fcr_unused_name(names(d_used), ".fcr_boot_w")
  refit <- function(idx, w = NULL) {
    rows <- unlist(unit_rows[idx], use.names = FALSE)
    d <- d_used[rows, , drop = FALSE]
    stb <- st
    if (!is.null(spec$unit)) {
      d[[gname]] <- rep(seq_along(idx), lengths(unit_rows[idx]))
      stb$group_var <- gname
    }
    if (!is.null(w)) {
      d[[wname]] <- rep(w, lengths(unit_rows[idx]))
      stb$weights <- wname
    } else if (!is.null(st$weights)) {
      d[[wname]] <- spec$w_row[rows]
      stb$weights <- wname
    }
    sp <- .fcr_fm_setup(stb, d)
    sp$st <- stb
    s <- .fcr_fm_search(sp, st$control$boot_starts + 1L,
                        sample.int(.Machine$integer.max, 1L),
                        start.vals = theta)
    est <- s$fit$estimate
    # Align labels with the original estimates
    dist <- apply(perms, 1, function(p) {
      sum(((.fcr_fm_permute(est, map, p) - theta) / scale)^2)
    })
    o <- order(dist)
    if (length(o) > 1L && dist[o[2]] < 1.5 * dist[o[1]]) {
      n_amb <<- n_amb + 1L
    }
    est <- .fcr_fm_permute(est, map, perms[o[1], ])
    names(est) <- map$names
    .fcr_refit_result(est, s$fit$converged)
  }
  freq <- if (!is.null(st$weights)) spec$w_unit
  draws <- .fcr_with_seed(seed + 1, .fcr_bootstrap(
    n = spec$n_u, B = B, coef_names = map$names, refit = refit, freq = freq))
  if (n_amb > 0L) {
    warning(n_amb, " of ", B, " bootstrap refits matched the original ",
            "components ambiguously (two labelings almost equally close). ",
            "Class-specific standard errors may be unreliable; label-",
            "invariant quantities (overall predictions) are not affected.",
            call. = FALSE)
  }
  list(draws = draws, failures = attr(draws, "failures"),
       n_ok = sum(stats::complete.cases(draws)), n_ambiguous = n_amb)
}
