#' Random Parameters Count Regression Models
#'
#' @name countreg.rp
#' @aliases countreg_rp
#' @param formula an R formula. This formula should specify the outcome and the
#'   independent variables that have fixed parameters.
#' @param rpar_formula a symbolic description of the model related specifically
#'   to the random parameters. This should not include an outcome variable. If
#'   the intercept is random, include it in this formula. If the intercept is
#'   fixed, include it in \code{formula} but not in \code{rpar_formula}.
#' @param data a dataframe that has all of the variables in the \code{formula}
#'   and \code{rpar_formula}.
#' @param family the name of the distribution/model type to estimate. Default is
#'   "NB2". Options include "Poisson", "NB1", "NB2", "NBP", "PIG", "Sichel",
#'   etc. (see \code{\link{count_families}} for the full list). DW, DPO,
#'   HYPERPOIS and GAMMACOUNT solve their mean equation for every
#'   observation and draw, so they are much slower here than closed-form
#'   families; consider fewer \code{ndraws} while exploring.
#'   With \code{"GP2"}, the dispersion must be nonnegative: the
#'   underdispersed GP2 requires \eqn{1 + a\mu > 0}, which random
#'   parameters with unbounded support violate for some draws of every
#'   observation. Use \code{"GP1"}, whose parameter space does not depend
#'   on \eqn{\mu}, for underdispersed random-parameter models.
#' @param rpardists an optional named vector whose names are the random
#'   parameters and values the distribution. The distribution options include
#'   normal ("n"), lognormal ("ln"), triangular ("t"), uniform ("u"), and gamma
#'   ("g"). Values are matched to the random parameters by name; every random
#'   parameter must be named. A single unnamed value applies to all random
#'   parameters. If not provided, normal is used. See Details for how each
#'   distribution's \code{St.Dev} parameter is defined.
#' @param dis_param_formula_1 a symbolic description of the model for the 
#'        first parameter of the count distribution (e.g., ln(alpha) for NB2).
#' @param dis_param_formula_2 a symbolic description of the model for the 
#'        second parameter of the count distribution (if applicable).
#' @param het_mean_formula an optional one-sided formula for heterogeneity in
#'   the means of the random parameters. Each random parameter gets its own
#'   coefficient for each variable (named \code{HetMean:<parameter>:<variable>}),
#'   and the effect is additive: \eqn{\beta_{ik} = \mu_k + z_i'\delta_k}
#'   (Wood and Gayah, 2025, Eq. 2). For lognormal parameters the shift is on the
#'   log scale, and for gamma parameters it is multiplicative,
#'   \eqn{\mu_k \exp(z_i'\delta_k)}, so the mean stays positive. The
#'   intercept is dropped.
#' @param het_var_formula an optional one-sided formula for heterogeneity in
#'   the scales of the random parameters:
#'   \eqn{\sigma_{ik} = \sigma_k \exp(s_i'\phi)}, with \eqn{\phi} shared by
#'   all random parameters. The intercept is dropped.
#' @param ndraws the number of Halton draws to use for estimating the random 
#'        parameters.
#' @param scrambled if the Halton draws should be scrambled. Scrambled
#'   draws are generated from a random seed that is stored in the fit, so
#'   predictions reuse the same draws (and \code{set.seed()} before fitting
#'   makes the fit reproducible).
#' @param mix_draws with \code{engine = "halton"}, the number of Halton
#'   draws for integrating over the mixing distribution of the PLN, PGE, PLL
#'   and PW families (inside each random-parameter draw). The default is
#'   \code{min(ndraws, 200)}. Not used by the default
#'   \code{engine = "quadrature"}, which integrates adaptively.
#' @param correlated if the random parameters should be correlated. 
#'        If TRUE, only normal distributions are used.
#' @param group_var an optional variable name (string), or a vector with one
#'   value per row, defining the panel structure (repeated measures). If
#'   provided, each panel shares one draw of the random parameters, and the
#'   likelihood and standard errors account for the panel structure.
#' @param panel_id deprecated; use \code{group_var}.
#' @param stderr the type of standard errors: \code{"normal"} (Hessian-based,
#'   the default), \code{"robust"} (sandwich; clustered by panel when
#'   \code{group_var} is given), or \code{"boot"} (bootstrap that resamples
#'   whole panels and refits the model).
#' @param bootstraps the number of bootstrap samples when
#'   \code{stderr = "boot"} (default 200).
#' @param keep_data logical; if \code{FALSE}, the training data are not stored
#'   in the fitted object. \code{predict()} then needs \code{newdata}.
#' @param offset the name of a variable, or vector of variable names, to be
#'   used as an offset (summed when several are given). An offset can instead
#'   be given with \code{offset()} in \code{formula}; use only one of the two.
#' @param weights the name of a variable, or a numeric vector, of frequency
#'   weights. With \code{panel_id}, weights apply to whole panels (the panel
#'   log-likelihood is weighted), so they should be constant within a panel.
#' @param method optimization method (e.g., "BHHH", "BFGS", "NM").
#' @param max.iters maximum number of iterations.
#' @param engine the numerical method for the mixed Poisson families PLN, PGE, PW
#'   and PLL: \code{"quadrature"} (default; adaptive quadrature on the log
#'   scale, accurate for every count), \code{"halton"} (the average over
#'   Halton draws used by flexCountReg 0.1.2 and earlier, kept to reproduce
#'   earlier fits; it is badly biased for counts far in the tail) or, for PLN
#'   only, \code{"poilog"} (the poilog package).
#' @param start.vals optional vector of starting values, one per parameter.
#'   Named values are matched to the parameters by name.
#' @param verbose logical.
#' 
#' @details
#' \code{countreg_rp()} is the same function as \code{countreg.rp()}, under a
#' name that tools do not mistake for an S3 method.
#'
#' \strong{Random parameter distributions.} Each random parameter is
#' \eqn{\beta = \mu + \sigma \upsilon}, where the \code{:Mean} and
#' \code{:St.Dev} estimates are \eqn{\mu} and \eqn{\sigma}:
#' \itemize{
#'   \item normal ("n"): \eqn{\upsilon \sim N(0, 1)}, so \eqn{\sigma} is the
#'   standard deviation.
#'   \item uniform ("u") and triangular ("t"): \eqn{\upsilon} is uniform or
#'   symmetric triangular on \eqn{[-1, 1]}, so \eqn{\beta} lies in
#'   \eqn{[\mu - \sigma, \mu + \sigma]}. Here \eqn{\sigma} is a half-width;
#'   the standard deviation is \eqn{\sigma/\sqrt{3}} (uniform) or
#'   \eqn{\sigma/\sqrt{6}} (triangular). Models fitted with version 0.1.1 or earlier
#'   used \eqn{[\mu - \sigma/2, \mu + \sigma/2]} for the uniform; their
#'   \eqn{\sigma} is converted automatically in \code{predict()}.
#'   \item lognormal ("ln"): \eqn{\beta = \exp(\mu + \sigma \upsilon)},
#'   always positive. The expected value of \eqn{\exp(x\beta)} is infinite
#'   for \eqn{x > 0}, so marginal (Exact or Simulated) predictions do not exist
#'   for observations with a positive value of a lognormal random parameter's
#'   variable. Use the variable's negative (so \eqn{x \le 0}), a gamma random
#'   parameter, or \code{method = "Individual"} predictions instead.
#'   \item gamma ("g"): mean \eqn{\mu > 0} and standard deviation
#'   \eqn{\sigma}, so the coefficient is always positive.
#' }
#'
#' The mixing distribution of the PLN, PGE, PLL and PW families uses its own
#' Halton dimension, independent of the random parameters' draws.
#'
#' @references
#' Wood, J.S., Gayah, V. (2025). Out-of-sample prediction and interpretation for
#' random parameter generalized linear models. \emph{Accident Analysis and
#' Prevention}, 220, 108147.
#'
#' @returns
#' An object of class `flexCountReg` which is a list with the following components:
#' \itemize{
#'  \item model: the fitted model object.
#'  \item data: the data frame used to fit the model.
#'  \item call: the matched call.
#'  \item formula: the formula used to fit the model.
#' }
#' 
#' @examples
#' \donttest{
#' # Load data
#' data("washington_roads")
#' washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
#' 
#' # 1. Basic Random Parameters Negative Binomial (NB2)
#' rp_nb2 <- countreg.rp(Total_crashes ~ lnaadt + lnlength,
#'                       rpar_formula = ~ -1 + speed50,
#'                       data = washington_roads,
#'                       family = "NB2",
#'                       rpardists = c(speed50 = "n"),
#'                       ndraws = 100,
#'                       method = "BHHH")
#' summary(rp_nb2)
#' 
#' # 2. Random Parameters with Panel Structure (if 'site_id' exists)
#' # rp_panel <- countreg.rp(Total_crashes ~ -1 + lnaadt,
#' #                        rpar_formula = ~ speed50,
#' #                        data = washington_roads,
#' #                        group_var = "site_id",
#' #                        family = "NB2",
#' #                        ndraws = 100)
#' 
#' # 3. Generalized Random Parameters Model with Heterogeneity
#' rp_gen <- countreg.rp(Total_crashes ~ lnaadt,
#'                       rpar_formula = ~ -1 + speed50,
#'                       dis_param_formula_1 = ~ lnlength, 
#'                       het_mean_formula = ~ AADT10kplus,
#'                       data = washington_roads,
#'                       family = "NB2",
#'                       rpardists = c(speed50 = "n"),
#'                       ndraws = 100)
#' summary(rp_gen)
#' 
#' # 4. Random Parameters Poisson Model with panel specification
#' rp_poisson <- countreg.rp(Total_crashes ~ lnaadt,
#'                       rpar_formula = ~ -1 + speed50,
#'                       het_mean_formula = ~ AADT10kplus,
#'                       data = washington_roads,
#'                       family = "POISSON",
#'                       rpardists = c(speed50 = "n"),
#'                       ndraws = 100,
#'                       group_var = "ID")
#' summary(rp_poisson)
#' }
#' 
#' @importFrom maxLik maxLik
#' @importFrom randtoolbox halton
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom stats terms reformulate coef
#' @importFrom tibble as_tibble
#' 
#' @include helpers.R tri.R  createFlexCountReg.R
#' @export
countreg.rp <- function(formula, rpar_formula, data, family = "NB2",
                         rpardists = NULL,
                         dis_param_formula_1 = NULL,
                         dis_param_formula_2 = NULL,
                         het_mean_formula = NULL, 
                         het_var_formula = NULL,
                         ndraws = 500, scrambled = FALSE,
                         correlated = FALSE, group_var = NULL,
                         weights = NULL, offset = NULL,
                         method = 'BHHH', max.iters = 1000,
                         start.vals = NULL, verbose = FALSE,
                         stderr = "normal", bootstraps = NULL,
                         engine = "quadrature", keep_data = TRUE,
                         panel_id = NULL, mix_draws = NULL) {
  
  print.level <- if (verbose) 2 else 0
  if (!is.null(panel_id)) {
    .fcr_deprecated_arg("countreg.rp", "panel_id", "group_var")
    if (is.null(group_var)) group_var <- panel_id
  }
  panel_id <- group_var
  stderr <- .fcr_resolve_stderr(stderr, missing(stderr), bootstraps)
  
  # --- 1. Family, method and argument checks ---
  family <- toupper(gsub("[^[:alnum:]]", "", family))
  if (!family %in% .fcr_families) {
    stop("Unknown `family` \"", family, "\". Must be one of: ",
         toString(.fcr_families), ".", call. = FALSE)
  }
  params <- get_params(family)
  engine <- match.arg(engine, c("quadrature", "halton", "poilog"))
  if (engine == "poilog" && family == "PLN") {
    rlang::check_installed("poilog", reason = "for `engine = \"poilog\"`.")
  }
  method <- .fcr_check_method(method)
  .fcr_check_controls(max.iters, ndraws)
  .fcr_check_dis_formulas(family, params, dis_param_formula_1,
                          dis_param_formula_2)
  
  # Panel id: a column name, or a vector with one value per row
  data <- as.data.frame(data)
  group_var <- NULL
  if (!is.null(panel_id)) {
    if (is.character(panel_id) && length(panel_id) == 1L) {
      if (!panel_id %in% names(data)) {
        stop("`group_var` \"", panel_id, "\" is not a column of `data`.",
             call. = FALSE)
      }
      group_var <- panel_id
    } else if (length(panel_id) == nrow(data)) {
      data$.fcr_panel_id <- panel_id
      group_var <- ".fcr_panel_id"
    } else {
      stop("`group_var` must be a column name in `data` or a vector with ",
           "one value per row of `data`.", call. = FALSE)
    }
  }
  
  if (correlated && !is.null(rpardists) && any(rpardists != "n")) {
    warning("When correlated = TRUE, only normal distributions are used. ",
            "Resetting rpardists to 'n'.", call. = FALSE)
    rpardists <- NULL
  }
  
  # --- 2. Data: one complete-case model frame for everything ---
  prep <- prepare_model_data(formula, data,
                             aux_formulas = list(
                               rpar = rpar_formula,
                               dis1 = dis_param_formula_1,
                               dis2 = dis_param_formula_2,
                               het_mean = het_mean_formula,
                               het_var = het_var_formula),
                             offset = offset, weights = weights,
                             group_var = group_var)
  # Keep the formulas with `.` expanded, for prediction, methods and refits
  # (a refit must not expand `.` against the bootstrap panel column)
  formula <- prep$formula
  rpar_formula <- prep$aux_formulas$rpar
  dis_param_formula_1 <- prep$aux_formulas$dis1
  dis_param_formula_2 <- prep$aux_formulas$dis2
  het_mean_formula <- prep$aux_formulas$het_mean
  het_var_formula <- prep$aux_formulas$het_var
  if (prep$n_dropped > 0) {
    message(prep$n_dropped, " row(s) with missing values were removed.")
  }
  y <- prep$y
  N_obs <- length(y)
  X_Fixed <- prep$X
  X_offset <- prep$offset
  X_dis_1 <- prep$aux$dis1
  X_dis_2 <- prep$aux$dis2
  X_het_mean <- prep$aux$het_mean
  X_het_var <- prep$aux$het_var
  panel_group <- if (is.null(prep$group)) factor(seq_len(N_obs)) else
    prep$group
  
  # Panel-level frequency weights: the panel log-likelihood is weighted
  w_obs <- prep$weights
  w_panel <- tapply(w_obs, panel_group, `[`, 1)
  w_panel <- as.numeric(w_panel[levels(panel_group)])
  if (!is.null(weights) &&
      any(tapply(w_obs, panel_group, function(v) diff(range(v))) > 0)) {
    stop("`weights` vary within panels. Weights are frequency weights of ",
         "the panels (the independent units), so they must be constant ",
         "within each panel.", call. = FALSE)
  }
  if (stderr == "boot" && !is.null(weights)) {
    .fcr_check_boot_freq(w_panel, length(w_panel))
  }
  
  # --- 3. Random parameters ---
  X_rand <- prep$aux$rpar
  rand_has_int <- "(Intercept)" %in% colnames(X_rand)
  if (rand_has_int) {
    named_int <- !is.null(names(rpardists)) &&
      any(grepl("intercept", names(rpardists), ignore.case = TRUE))
    if ("(Intercept)" %in% colnames(X_Fixed)) {
      if (named_int) {
        stop("The intercept is in both `formula` and `rpar_formula`. For a ",
             "random intercept, remove it from `formula` (add `- 1`).",
             call. = FALSE)
      }
      # The intercept implied by rpar_formula is not a random parameter
      X_rand <- X_rand[, colnames(X_rand) != "(Intercept)", drop = FALSE]
    }
  }
  rpar <- colnames(X_rand)
  N_rand <- length(rpar)
  if (N_rand == 0) {
    stop("`rpar_formula` has no random parameters.", call. = FALSE)
  }
  rpardists <- .fcr_match_rpardists(rpardists, rpar, correlated)
  # A variable in both formulas: with a free-location random parameter
  # (normal, uniform, triangular) only the sum of the fixed coefficient and
  # the random parameter's mean is identified. Lognormal and gamma random
  # parameters have a fixed support, so a fixed shift is identified.
  both <- intersect(rpar, colnames(X_Fixed))
  free_loc <- both[rpardists[both] %in% c("n", "u", "t")]
  if (length(free_loc)) {
    stop("Variable(s) ", toString(free_loc), " are in both `formula` and ",
         "`rpar_formula` with a normal, uniform or triangular random ",
         "parameter. Only the sum of the fixed coefficient and the random ",
         "parameter's mean is identified; remove them from `formula` (the ",
         "random parameter's mean is the average effect).", call. = FALSE)
  }
  
  N_fixed <- ncol(X_Fixed)
  
  # Distribution parameter counts
  N_dis_1 <- if (is.null(params[[1]])) 0 else if (is.null(X_dis_1)) 1 else
    ncol(X_dis_1)
  N_dis_2 <- if (is.null(params[[2]])) 0 else if (is.null(X_dis_2)) 1 else
    ncol(X_dis_2)
  N_het_mean <- if (is.null(X_het_mean)) 0 else ncol(X_het_mean) * N_rand
  N_het_var <- if (is.null(X_het_var)) 0 else ncol(X_het_var)
  
  # Halton draws: one dimension per random parameter, plus a separate one
  # for the mixing distribution of PLN, PGE, PLL and PW families
  if (is.null(mix_draws)) mix_draws <- min(ndraws, 200)
  .fcr_check_controls(ndraws = mix_draws)
  draw_seed <- if (scrambled) sample.int(.Machine$integer.max, 1L)
  hdraws <- .fcr_halton(ndraws, N_rand + 1, scrambled,
                        draw_seed)[, seq_len(N_rand), drop = FALSE]
  dist_haltons <- .fcr_halton(mix_draws, N_rand + 1, scrambled,
                              draw_seed)[, N_rand + 1]
  normed_dist_haltons <- stats::qnorm(dist_haltons)
  
  # --- 4. Parameter names and starting values ---
  x_names <- c(colnames(X_Fixed), paste0(rpar, ":Mean"))
  if (correlated) {
    x_names <- c(x_names, .fcr_chol_names(rpar))
  } else {
    x_names <- c(x_names, paste0(rpar, ":St.Dev"))
  }
  if (N_het_mean > 0) {
    x_names <- c(x_names, paste0("HetMean:", rep(rpar, each = ncol(X_het_mean)),
                                 ":", colnames(X_het_mean)))
  }
  if (N_het_var > 0) {
    x_names <- c(x_names, paste0("HetVar:", colnames(X_het_var)))
  }
  if (N_dis_1 > 0) {
    x_names <- c(x_names, if (is.null(X_dis_1)) params[[1]] else
      paste0(params[[1]], ":", colnames(X_dis_1)))
  }
  if (N_dis_2 > 0) {
    x_names <- c(x_names, if (is.null(X_dis_2)) params[[2]] else
      paste0(params[[2]], ":", colnames(X_dis_2)))
  }
  
  if (!is.null(start.vals)) {
    start <- .fcr_check_start(start.vals, x_names)
  } else {
    # Fixed and random means from a Poisson GLM on the same rows, offset
    # and weights
    X_glm <- cbind(X_Fixed, X_rand)
    glm_mod <- tryCatch(
      stats::glm.fit(X_glm, y, weights = w_obs, offset = X_offset,
                     family = stats::poisson()),
      error = function(e) NULL)
    b <- if (is.null(glm_mod)) rep(0, ncol(X_glm)) else glm_mod$coefficients
    b[is.na(b)] <- 0
    s_fixed <- b[seq_len(N_fixed)]
    s_rand <- b[N_fixed + seq_len(N_rand)]
    # Gamma random parameters need a positive mean
    is_g <- rpardists == "g"
    s_rand[is_g] <- pmax(abs(s_rand[is_g]), 0.1)
    if (correlated) {
      # Lower-triangular Cholesky factor 0.1 * I, packed
      s_var <- .fcr_chol_pack(diag(0.1, N_rand))
    } else {
      s_var <- ifelse(is_g, pmax(0.5 * s_rand, 0.05), 0.1)
    }
    # Scalar distribution parameters start near the Poisson limit (small
    # overdispersion); GP1/GP2 start exactly at it (dispersion 0), because a
    # strongly underdispersed start gives zero probability to larger counts
    # GP2 starts slightly inside its allowed region a >= 0 (see below)
    # Families with their own starting values (e.g. at their Poisson or
    # geometric limit, ln(shape) = 0) declare them in the registry
    fam_start <- .fcr_get_family(family)$start_rp
    dis_start <- if (family == "GP1") 0 else if (family == "GP2") 0.05 else
      if (!is.null(fam_start)) fam_start[1] else -2
    dis_start_2 <- if (!is.null(fam_start)) rep_len(fam_start, 2)[2] else -2
    start <- c(s_fixed, s_rand, s_var, rep(0, N_het_mean + N_het_var),
               rep(dis_start, N_dis_1 > 0 && is.null(X_dis_1)),
               rep(0, if (!is.null(X_dis_1)) N_dis_1 else 0),
               rep(dis_start_2, N_dis_2 > 0 && is.null(X_dis_2)),
               rep(0, if (!is.null(X_dis_2)) N_dis_2 else 0))
    names(start) <- x_names
  }
  
  # --- 5. Likelihood ---
  local_probFunc <- get_probFunc(family)
  n_var_params <- if (correlated) N_rand * (N_rand + 1) / 2 else N_rand
  
  ll_countreg_rp <- function(p) {
    idx <- 0
    take <- function(k) {
      if (k == 0) return(NULL)
      out <- p[idx + seq_len(k)]
      idx <<- idx + k
      out
    }
    fixed_coefs <- take(N_fixed)
    rand_means <- take(N_rand)
    rand_vars <- take(n_var_params)
    het_mean_coefs <- take(N_het_mean)
    if (!is.null(het_mean_coefs)) {
      het_mean_coefs <- matrix(het_mean_coefs, ncol = N_rand)
    }
    het_var_coefs <- take(N_het_var)
    p_dis_1 <- take(N_dis_1)
    p_dis_2 <- take(N_dis_2)
    
    vec_1 <- if (N_dis_1 == 0) NULL else if (is.null(X_dis_1))
      rep(exp(p_dis_1), N_obs) else exp(as.vector(X_dis_1 %*% p_dis_1))
    vec_2 <- if (N_dis_2 == 0) NULL else if (is.null(X_dis_2))
      rep(exp(p_dis_2), N_obs) else exp(as.vector(X_dis_2 %*% p_dis_2))
    
    # GP2 with negative dispersion a requires 1 + a * mu > 0, which
    # unbounded random parameters violate for some draws of every
    # observation: the model is only defined for a >= 0 (alpha = a + 1)
    if (family == "GP2" && any(vec_1 < 1)) {
      return(rep(-Inf, length(w_panel)))
    }
    mu_fixed <- exp(as.vector(X_Fixed %*% fixed_coefs) + X_offset)
    
    draws_info <- generate_random_draws(hdraws = hdraws, 
                                        random_coefs_means = rand_means, 
                                        rand_var_params = rand_vars, 
                                        rpardists = rpardists, 
                                        rpar = rpar, 
                                        X_rand = X_rand,
                                        het_mean_coefs = het_mean_coefs, 
                                        X_het_mean = X_het_mean,
                                        het_var_coefs = het_var_coefs, 
                                        X_het_var = X_het_var,
                                        correlated = correlated)
    mu_mat <- exp(draws_info$xb_rand_mat) * mu_fixed
    D <- ncol(mu_mat)
    
    probs_flat <- tryCatch(
      local_probFunc(y = rep(y, times = D), 
                     predicted = as.vector(mu_mat), 
                     alpha = if (!is.null(vec_1)) rep(vec_1, times = D),
                     sigma = if (!is.null(vec_2)) rep(vec_2, times = D),
                     haltons = dist_haltons, 
                     normed_haltons = normed_dist_haltons,
                     engine = engine, log = TRUE),
      error = function(e) rep(NA_real_, N_obs * D))
    
    # Panel log-likelihood per draw, averaged over draws on the log scale
    log_probs <- matrix(probs_flat, nrow = N_obs, ncol = D)
    sum_log <- rowsum(log_probs, group = panel_group, reorder = TRUE)
    m <- apply(sum_log, 1, max)
    ll_panel <- m + log(rowMeans(exp(sum_log - m)))
    ll_panel[!is.finite(m)] <- m[!is.finite(m)]
    # Panels with zero weight contribute exactly 0 (also when -Inf)
    ll_panel[w_panel == 0] <- 0
    
    ll_panel * w_panel
  }
  
  # --- 6. Optimization ---
  fit <- maxLik::maxLik(
    ll_countreg_rp,
    start = start,
    method = method,
    control = list(iterlim = max.iters, printLevel = print.level))
  names(fit$estimate) <- x_names
  fit$converged <- .fcr_converged(fit, method)
  if (!isTRUE(fit$converged)) {
    warning("The optimizer did not report convergence (",
            trimws(fit$message), "). Consider a different `method` or more ",
            "`max.iters`.", call. = FALSE)
  }
  
  # --- 7. Result Packaging ---
  fit$x_names <- x_names
  fit$formula <- formula
  fit$rpar_formula <- rpar_formula
  fit$rpar <- rpar
  fit$het_mean_formula <- het_mean_formula
  fit$het_var_formula <- het_var_formula
  fit$dis_param_formula_1 <- dis_param_formula_1
  fit$dis_param_formula_2 <- dis_param_formula_2
  fit$family <- family
  fit$modelType <- "countreg.rp"
  fit$method <- method
  fit$ndraws <- ndraws
  fit$mix_draws <- mix_draws
  fit$draw_seed <- draw_seed
  fit$scrambled <- scrambled
  fit$correlated <- correlated
  fit$panel_id <- panel_id 
  fit$group_var <- panel_id
  fit$offset <- prep$blueprint$offset_label
  fit$engine <- engine
  fit$rpardists <- rpardists
  # Draw conventions (flexCountReg >= 0.1.2): uniform draws on
  # [mu - sigma, mu + sigma], additive heterogeneity in means with one
  # coefficient per random parameter, separate Halton dimension for the
  # mixing distribution. Fits without this field used the older conventions.
  fit$draw_version <- 2L
  fit$blueprint <- prep$blueprint
  fit$y <- y
  fit$nobs <- N_obs
  fit$n_groups <- nlevels(panel_group)
  # Independent units for BIC: the panels, counted with their frequency
  # weights
  fit$n_bic <- sum(w_panel)
  # Panel of each analysis row (the independent units), for refits
  fit$panel_index <- if (!is.null(prep$group)) as.integer(panel_group)
  # Sample evidence for likelihood-ratio comparisons
  fit$lr_sample <- .fcr_lr_sample(prep, data, group_var)
  fit$correlated <- correlated
  fit$weights <- prep$weights
  fit$rows_used <- prep$rows_used
  fit$n_dropped <- prep$n_dropped
  
  if (stderr == "boot") {
    if (is.null(bootstraps)) {
      bootstraps <- 200L
      message("`bootstraps` not given; using 200 bootstrap samples.")
    }
    # Resample whole panels (rows when there is no panel id) and refit
    # from the full-sample estimates
    boot_data <- data[prep$rows_used, , drop = FALSE]
    rows_by_panel <- split(seq_len(N_obs), panel_group)
    # Frequency weights are constant within a panel: a panel with weight w
    # stands for w independent panels, so new panel frequencies are drawn
    # (see .fcr_bootstrap) and replace the weights of the panel's rows
    boot_freq <- if (!is.null(weights)) w_panel
    draws <- .fcr_bootstrap(
      n = length(rows_by_panel), B = as.integer(bootstraps),
      coef_names = x_names, freq = boot_freq,
      refit = function(idx, w = NULL) {
        sel <- rows_by_panel[idx]
        r <- unlist(sel, use.names = FALSE)
        bd <- boot_data[r, , drop = FALSE]
        boot_id <- .fcr_unused_name(names(bd), ".fcr_boot_panel")
        bd[[boot_id]] <- rep(seq_along(sel), lengths(sel))
        w_r <- if (!is.null(w)) rep(w, lengths(sel))
        m <- countreg.rp(formula = formula, rpar_formula = rpar_formula,
                         data = bd, family = family,
                         rpardists = rpardists,
                         dis_param_formula_1 = dis_param_formula_1,
                         dis_param_formula_2 = dis_param_formula_2,
                         het_mean_formula = het_mean_formula,
                         het_var_formula = het_var_formula,
                         ndraws = ndraws, scrambled = scrambled,
                         mix_draws = mix_draws,
                         correlated = correlated,
                         group_var = boot_id,
                         weights = w_r,
                         offset = offset, method = method,
                         max.iters = max.iters, start.vals = fit$estimate,
                         engine = engine, keep_data = FALSE)
        .fcr_refit_result(m$model$estimate, m$model$converged)
      })
    V <- .fcr_boot_cov(draws)
    fit$boot_failures <- attr(draws, "failures")
    fit$n_boot_ok <- sum(stats::complete.cases(draws))
    fit$bootstraps <- as.integer(bootstraps)
    vcov_type <- "bootstrap"
  } else if (stderr == "robust") {
    # One log-likelihood contribution per panel: cluster-robust sandwich,
    # with frequency weights per panel
    V <- .fcr_robust_vcov(fit, w_panel)
    vcov_type <- "robust"
  } else {
    V <- .fcr_hessian_vcov(fit)
    vcov_type <- "hessian"
  }
  dimnames(V) <- list(x_names, x_names)
  fit$vcov <- V
  fit$vcov_type <- vcov_type
  fit$stderr <- stderr
  fit$se <- .fcr_se(V)
  names(fit$se) <- x_names
  fit$bootstrapped_se <- if (stderr == "boot") fit$se else NULL
  
  obj <- .createFlexCountReg(model = fit, 
                             data = data, keep_data = keep_data,
                             call = match.call(), 
                             formula = formula)
  return(obj)
}

# Validate random-parameter distributions and order them like `rpar`.
.fcr_match_rpardists <- function(rpardists, rpar, correlated) {
  valid <- c("n", "ln", "t", "u", "g")
  if (is.null(rpardists) || correlated) {
    return(stats::setNames(rep("n", length(rpar)), rpar))
  }
  rpardists <- unlist(rpardists)
  # Check the values before any expansion, so a scalar is validated too
  if (!is.character(rpardists) || !length(rpardists) ||
      anyNA(rpardists)) {
    stop("`rpardists` must be a character vector of distribution codes (",
         toString(valid), ") without missing values.", call. = FALSE)
  }
  bad <- setdiff(unique(rpardists), valid)
  if (length(bad)) {
    stop("Unknown random-parameter distribution(s): ", toString(bad),
         ". Use one of ", toString(valid), ".", call. = FALSE)
  }
  nm <- names(rpardists)
  if (!is.null(nm) && (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm))) {
    stop("Names of `rpardists` must be non-empty and unique.", call. = FALSE)
  }
  if (is.null(nm)) {
    if (length(rpardists) == 1L) {
      return(stats::setNames(rep(rpardists, length(rpar)), rpar))
    }
    if (length(rpardists) != length(rpar)) {
      stop("`rpardists` must be a named vector with one entry per random ",
           "parameter: ", toString(rpar), ".", call. = FALSE)
    }
    nm <- rpar
  }
  nm[grepl("^\\(?intercept\\)?$", nm, ignore.case = TRUE)] <- "(Intercept)"
  names(rpardists) <- nm
  miss <- setdiff(rpar, nm)
  extra <- setdiff(nm, rpar)
  if (length(miss) || length(extra)) {
    stop("Names of `rpardists` must match the random parameters (",
         toString(rpar), ").",
         if (length(miss)) paste0(" Missing: ", toString(miss), "."),
         if (length(extra)) paste0(" Unknown: ", toString(extra), "."),
         call. = FALSE)
  }
  rpardists[rpar]
}

#' @rdname countreg.rp
#' @usage NULL
#' @export
countreg_rp <- countreg.rp
