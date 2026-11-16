#' Predictions for flexCountReg models
#'
#' @name predict.flexCountReg
#' @param object a model object estimated using this R package.
#' @param newdata optional data frame for which to generate predictions. It
#'   must contain every variable used by the model's formulas (including
#'   \code{rpar_formula}, \code{dis_param_formula_1/2},
#'   \code{underreport_formula} and heterogeneity formulas) and any offset
#'   columns. The outcome is not required. If omitted, predictions are made for
#'   the data used to fit the model.
#' @param new_data the same as \code{newdata}, under the name used by
#'   tidymodels. Give one or the other.
#' @param ... optional arguments passed to the function. This includes `method`.
#'   Passing the data as `data =` is deprecated; use `newdata`.
#' @param type the type of prediction. \code{"response"} (default) is the
#'   expected observed count. For models with underreporting this is
#'   \eqn{\mu \times p}, the true mean times the probability of being reported.
#'   \code{"true"} is the true expected count \eqn{\mu} (excluding the
#'   underreporting probability), \code{"prob_reported"} is the probability of
#'   being reported \eqn{p}, and \code{"link"} is the linear predictor
#'   \eqn{X\beta} plus the offset. For models without underreporting,
#'   \code{"true"} equals \code{"response"} and \code{"prob_reported"} is 1.
#'   For zero-inflated and hurdle models (\code{zero_formula}),
#'   \code{"response"} includes the zero model, \code{"true"} is the mean of
#'   the count distribution alone, and \code{"zero"} is the probability from
#'   the zero model: the probability of a structural zero (zero-inflated) or
#'   of a zero count (hurdle).
#'   \code{"prob"} returns a matrix of probabilities \eqn{P(Y = k)} with one
#'   row per row of \code{newdata} and one column per count
#'   \eqn{k = 0, \dots,} \code{max_count} (passed in \code{...}; default the
#'   largest count in the training data), for rootograms and for scoring
#'   models out of sample. With \code{overflow = TRUE} (in \code{...}), a last
#'   column holds the tail mass \eqn{P(Y > } \code{max_count}\eqn{)}, so each
#'   row sums to 1. Types other than \code{"response"} are available
#'   for \code{countreg} models only.
#' 
#' @details
#' The returned vector always has one element per row of \code{newdata}, in the
#' same order. Rows with missing values in any variable the model uses get
#' \code{NA}. Factor levels and contrasts from the training data are used, so
#' \code{newdata} may contain only a subset of the levels. An error is raised
#' if an offset column is missing from \code{newdata}.
#' @note optional parameter `method`: Only valid for random parameters models
#'   (`countreg.rp`). Options include \code{Simulated} (default),
#'   \code{Individual} (or \code{EB}), or \code{Exact}.
#' @note optional parameter `eb_type`: for \code{method = "Individual"},
#'   either \code{"coefficients"} (default) or \code{"mean"}.
#'
#' @description
#' Generates predictions for the expected count (lambda) for observations.
#' 
#' For \strong{countreg.rp} (Random Parameters) models, three methods are
#' available:
#' \itemize{
#'   \item \strong{Simulated}: Uses Halton draws to simulate the random
#'   parameters and averages the outcomes. This is a simulation-based
#'   approximation.
#'   \item \strong{Individual} (alias \code{EB}): Empirical Bayes
#'   observation-specific predictions (Wood and Gayah, 2025, Eqs. 5-7). Each
#'   simulation draw is weighted by the likelihood of the unit's observed
#'   outcomes (all of a panel's rows when \code{panel_id} was used), giving
#'   conditional coefficients \eqn{E[\beta_i | y_i]}. By default these are
#'   used in the mean function, as in the paper; pass
#'   \code{eb_type = "mean"} to instead return the conditional mean
#'   \eqn{E[\mu_{it} | y_i]}. The conditional coefficients and their
#'   conditional standard deviations are returned as the attributes
#'   \code{"eb_coefficients"} and \code{"eb_sd"}. Requires the outcome (and
#'   the \code{panel_id} column, if used) in \code{newdata}, so it is meant
#'   for data with observed outcomes.
#'   \item \strong{Exact}: Uses the moment generating functions (MGFs) of
#'   the random parameter distributions (Wood and Gayah, 2025, Appendix A),
#'   including heterogeneity in means and variances. This is exact for normal,
#'   uniform, triangular and gamma parameters and has no simulation error. For
#'   lognormal parameters it uses a Laplace (saddlepoint) approximation, which
#'   is accurate for \eqn{x < 0}; for \eqn{x > 0} (and for gamma parameters
#'   with \eqn{x\sigma^2/\mu \ge 1}) the expected value is infinite and
#'   \code{Inf} is returned with a warning.
#' }
#' 
#' For \strong{countreg}, \strong{poisLindRE}, and \strong{RENB} models, the
#' function calculates the expected value \eqn{\mu = \exp(X\beta)} (with
#' appropriate adjustments for specific families like PLN, PLL or
#' underreporting).
#' 
#' @references
#' Wood, J.S., Gayah, V. (2025). Out-of-sample prediction and interpretation for
#' random parameter generalized linear models. \emph{Accident Analysis and
#' Prevention}, 220, 108147.
#' 
#' @returns A numeric vector of predicted expected counts for each observation 
#'  in the provided data. If no data is provided, the predictions for the data 
#'  used in estimating the model are provided.
#'
#' @importFrom randtoolbox halton
#' @importFrom utils head tail
#' @include corr_haltons.R halton_dists.R  helpers.R
#' 
#' @examples
#' \donttest{
#' # Load data and create a dummy variable
#' data("washington_roads")
#' washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
#' 
#' # =========================================================================
#' # 1. Fixed Parameter Model (countreg)
#' # =========================================================================
#' nb2_fixed <- countreg(Total_crashes ~ lnaadt + lnlength + speed50,
#'                       data = washington_roads, 
#'                       family = "NB2")
#' pred_fixed <- predict(nb2_fixed, newdata = washington_roads)
#' 
#' # =========================================================================
#' # 2. Random Parameters Model (countreg.rp)
#' # =========================================================================
#' rp_nb2 <- countreg.rp(Total_crashes ~ lnaadt + lnlength,
#'                       rpar_formula = ~ -1 + speed50,
#'                       data = washington_roads,
#'                       family = "NB2",
#'                       rpardists = c(speed50 = "n"),
#'                       ndraws = 100)
#' 
#' # Method A: Simulated (Default)
#' pred_sim <- predict(rp_nb2, newdata = washington_roads, method = "Simulated")
#' 
#' # Method B: Exact (Analytical MGF)
#' pred_exact <- predict(rp_nb2, newdata = washington_roads, method = "Exact")
#'
#' # Method C: Individual (Empirical Bayes, requires the observed outcome)
#' pred_eb <- predict(rp_nb2, newdata = washington_roads, method = "Individual")
#' 
#' # =========================================================================
#' # 3. Random Effects Models (poisLindRE / RENB)
#' # =========================================================================
#' pl_re <- poisLind.re(Total_crashes ~ lnaadt + lnlength,
#'                      data = washington_roads,
#'                      group_var = "ID")
#' pred_pl_re <- predict(pl_re, newdata = washington_roads)
#' }
#' @export
predict.flexCountReg <- function(object, newdata = NULL, ...,
                                 type = c("response", "true",
                                          "prob_reported", "link", "prob",
                                          "zero"),
                                 new_data = NULL){
  
  type <- match.arg(type)
  
  # Extract optional parameters from '...'
  additional_args  <- list(...)
  
  # `new_data` is the tidymodels name for `newdata`
  if (!is.null(new_data) && is.null(newdata)) newdata <- new_data
  if (!is.null(newdata)) {
    data <- as.data.frame(newdata)
  } else if (!is.null(additional_args$data)) {
    warning("The `data` argument of predict() is deprecated; use `newdata`.",
            call. = FALSE)
    data <- as.data.frame(additional_args$data)
  } else {
    # Default to training data
    data <- as.data.frame(object$data)
  }
  
  model <- object$model
  
  if (!is.null(object$modelType)) {
    modtype <- object$modelType
  } else if (!is.null(model$modelType)) {
    modtype <- model$modelType
  } else {
    modtype <- "countreg" 
  }
  
  if (type == "prob") {
    max_count <- additional_args$max_count
    if (is.null(max_count)) {
      y_tr <- object$model[["y"]]
      max_count <- if (!is.null(y_tr)) max(y_tr) else 20
    }
    if (!is.numeric(max_count) || length(max_count) != 1L ||
        !is.finite(max_count) || max_count < 0 ||
        max_count != round(max_count)) {
      stop("`max_count` must be a single non-negative whole number.",
           call. = FALSE)
    }
    P <- .fcr_prob_matrix(object, data, as.integer(max_count))
    if (isTRUE(additional_args$overflow)) {
      # Tail mass P(Y > max_count), so each row sums to 1
      P <- cbind(P, pmax(0, 1 - rowSums(P)))
      colnames(P)[ncol(P)] <- paste0(">", max_count)
    }
    return(P)
  }
  
  if (type != "response" && modtype != "countreg") {
    stop("type = \"", type, "\" is only available for countreg() models.",
         call. = FALSE)
  }
  
  # Handle method argument
  if (is.null(additional_args$method)) {
    if(modtype == "countreg.rp") {
      method <- "Simulated"
    } else {
      method <- "Standard" 
    }
  } else {
    method <- additional_args$method
  }
  
  # =========================================================================
  # 1. RANDOM PARAMETERS MODELS (countreg.rp)
  # =========================================================================
  if(modtype == "countreg.rp"){
    
    method <- match.arg(method, c("Simulated", "Exact", "Individual", "EB"))
    if (method == "EB") method <- "Individual"
    
    family <- if(!is.null(model$family)) model$family else "NB2" 
    rpardists <- model$rpardists
    correlated <- model$correlated
    scrambled <- model$scrambled
    ndraws <- max(model$ndraws, 500)
    params <- get_params(family)
    coefs <- model$estimate
    
    # Conventions of the fit: fits made with flexCountReg 0.1.1 or earlier have no
    # draw_version and used uniform draws on [mu - sigma/2, mu + sigma/2] and
    # multiplicative heterogeneity in means shared by all random parameters.
    legacy <- is.null(model$draw_version)
    het_form <- if (legacy) "legacy" else "additive"
    
    # --- Matrices, built with the training terms and factor levels ---
    bp <- .fcr_blueprint(object, aux_formulas = list(
      rpar = model$rpar_formula,
      dis1 = model$dis_param_formula_1, dis2 = model$dis_param_formula_2,
      het_mean = model$het_mean_formula, het_var = model$het_var_formula))
    fd <- forge_new_data(bp, data)
    X_Fixed <- fd$X
    X_offset <- fd$offset
    
    # Random parameters, as named in the fitted coefficients
    rpar <- if (!is.null(model$rpar)) model$rpar else
      sub(":Mean$", "", grep(":Mean$", names(coefs), value = TRUE))
    X_rand <- fd$aux$rpar[, rpar, drop = FALSE]
    N_rand <- length(rpar)
    if (is.null(rpardists)) rpardists <- stats::setNames(rep("n", N_rand), rpar)
    if (!is.null(names(rpardists)) && all(rpar %in% names(rpardists))) {
      rpardists <- rpardists[rpar]
    }
    
    # --- Coefficients, matched by name ---
    fixed_coefs <- .fcr_coefs(coefs, colnames(X_Fixed))
    random_coefs_means <- .fcr_coefs(coefs, paste0(rpar, ":Mean"))
    if (correlated) {
      rand_var_params <- coefs[grepl("^Chol:", names(coefs))]
    } else {
      rand_var_params <- .fcr_coefs(coefs, paste0(rpar, ":St.Dev"))
      if (legacy && any(rpardists == "u")) {
        rlang::inform(paste(
          "This model was fit with flexCountReg 0.1.1 or earlier, whose uniform",
          "random parameters span mu +/- St.Dev/2. Its uniform St.Dev values",
          "are halved to the current convention (mu +/- St.Dev); refit the",
          "model to update the estimates."),
          .frequency = "once", .frequency_id = "fcr_legacy_uniform")
        rand_var_params[rpardists == "u"] <- rand_var_params[rpardists == "u"] / 2
      }
    }
    
    X_het_mean <- fd$aux$het_mean
    X_het_var <- fd$aux$het_var
    het_mean_coefs <- NULL
    if (!is.null(X_het_mean)) {
      het_mean_coefs <- if (legacy) {
        .fcr_coefs(coefs, paste0("HetMean:", colnames(X_het_mean)))
      } else {
        matrix(.fcr_coefs(coefs, paste0(
          "HetMean:", rep(rpar, each = ncol(X_het_mean)), ":",
          colnames(X_het_mean))), ncol = N_rand)
      }
    }
    het_var_coefs <- if (!is.null(X_het_var)) {
      .fcr_coefs(coefs, paste0("HetVar:", colnames(X_het_var)))
    }
    
    # Distribution parameters (for the PLN/PLL mean correction)
    vec_1 <- .fcr_dist_param(coefs, params[[1]], fd$aux$dis1, nrow(data))
    vec_2 <- .fcr_dist_param(coefs, params[[2]], fd$aux$dis2, nrow(data))
    mean_adj <- .fcr_mean_adjust(family, vec_1, vec_2)
    # For GP1/GP2 with a negative dispersion E[Y] is not proportional to the
    # mean parameter, so it is averaged over the draws instead
    linear_mean <- .fcr_mean_is_linear(family, vec_1)
    if (method == "Exact" && !linear_mean) {
      message("Exact predictions need a mean proportional to exp(X beta), ",
              "which does not hold for ", family, " with a negative ",
              "dispersion; using Simulated predictions.")
      method <- "Simulated"
    }
    
    mu_fixed <- exp(as.vector(X_Fixed %*% fixed_coefs) + X_offset)
    
    # --- Exact: moment generating functions (Wood and Gayah, 2025) ---
    # E[exp(x beta)] per row; also used to check, for Simulated predictions,
    # whether the expected value exists
    if (method %in% c("Exact", "Simulated")) {
      loc <- .fcr_rp_locations(random_coefs_means, rpardists, nrow(data),
                               het_mean_coefs, X_het_mean, het_var_coefs,
                               X_het_var, het_form)
      if (correlated) {
        L <- .fcr_chol_unpack(rand_var_params, N_rand)
        W <- X_rand * loc$var_multiplier
        rp_factor <- exp(rowSums(X_rand * loc$mu_mat) +
                           0.5 * rowSums((W %*% L)^2))
      } else {
        rp_factor <- rep(1, nrow(data))
        for (k in seq_len(N_rand)) {
          rp_factor <- rp_factor *
            rp_mgf(rpardists[[k]], X_rand[, k], loc$mu_mat[, k],
                   abs(rand_var_params[k]) * loc$var_multiplier)
        }
      }
      no_mean <- is.infinite(rp_factor)
      if (any(no_mean)) {
        warning(sum(no_mean), " prediction(s) are infinite: the expected ",
                "value does not exist for a lognormal random parameter with ",
                "a positive x, or a gamma random parameter with ",
                "x * sigma^2 / mu >= 1. See ?countreg.rp.", call. = FALSE)
      }
      if (method == "Exact") {
        return(as.vector(mu_fixed * rp_factor * mean_adj))
      }
    }
    
    # --- Simulated or Individual (Empirical Bayes) ---
    hdraws <- 
      .fcr_halton(ndraws, N_rand + 1, isTRUE(scrambled),
                  model$draw_seed)[, seq_len(N_rand), drop = FALSE]
    
    draw_args <- list(
      hdraws = hdraws, 
      random_coefs_means = random_coefs_means, 
      rand_var_params = rand_var_params, 
      rpardists = rpardists, 
      rpar = rpar, 
      het_mean_coefs = het_mean_coefs, 
      X_het_mean = X_het_mean,
      het_var_coefs = het_var_coefs, 
      X_het_var = X_het_var,
      correlated = correlated,
      het_form = het_form)
    
    draws_info <- do.call(generate_random_draws,
                          c(draw_args, list(X_rand = X_rand)))
    # Conditional means (N x draws), as in the likelihood
    mu_mat <- exp(draws_info$xb_rand_mat) * mu_fixed
    
    if (method == 'Simulated') {
      out <- if (linear_mean) {
        as.vector(rowMeans(mu_mat)) * mean_adj
      } else {
        D <- ncol(mu_mat)
        m <- .fcr_family_mean(family, as.vector(mu_mat), rep(vec_1, D),
                              if (!is.null(vec_2)) rep(vec_2, D))
        rowMeans(matrix(m, nrow = nrow(mu_mat)))
      }
      # A finite average of draws is not an estimate of a mean that does
      # not exist
      out[no_mean] <- Inf
      return(out)
    }
    
    dist_haltons <- if (legacy) {
      randtoolbox::halton(model$ndraws, 1)
    } else {
      .fcr_halton(model$mix_draws %||% model$ndraws, N_rand + 1,
                  isTRUE(scrambled), model$draw_seed)[, N_rand + 1]
    }
    return(.fcr_rp_eb_predict(
      model = model, data = data, family = family, mu_mat = mu_mat,
      mu_fixed = mu_fixed, X_rand = X_rand, draw_args = draw_args,
      vec_1 = vec_1, vec_2 = vec_2, mean_adj = mean_adj,
      dist_haltons = dist_haltons,
      eb_type = if (is.null(additional_args$eb_type)) "coefficients"
                else additional_args$eb_type))
  }
  
  # =========================================================================
  # 2a. NEGATIVE MULTINOMIAL FAMILY (negmultinom)
  # =========================================================================
  else if (modtype == "NegMultinom") {
    return(.nm_predict(object, data, method = method))
  }
  
  # =========================================================================
  # 2b. RANDOM EFFECTS MODELS (poisLindRE / RENB)
  # =========================================================================
  else if (modtype %in% c("poisLindRE", "RENB", "REPoisson")) {
    
    bp <- .fcr_blueprint(object)
    fd <- forge_new_data(bp, data)
    # Coefficients by name from the estimates (so that methods that perturb
    # the estimates, such as delta-method standard errors, change the
    # predictions); fits from earlier versions may only have beta_pred
    beta_pred <- if (all(colnames(fd$X) %in% names(model$estimate))) {
      .fcr_coefs(model$estimate, colnames(fd$X))
    } else {
      model$beta_pred
    }
    if (length(beta_pred) != ncol(fd$X)) {
      stop("The model's coefficients do not match its design matrix.",
           call. = FALSE)
    }
    
    mu <- exp(as.vector(fd$X %*% beta_pred) + fd$offset)
    return(mu)
  }
  
  # =========================================================================
  # 3. FIXED PARAMETER MODELS (countreg)
  # =========================================================================
  else {
    
    parts <- .fcr_countreg_parts(object, data)
    if (type == "zero" && is.null(parts$zero_prob)) {
      stop("type = \"zero\" requires a model with `zero_formula`.",
           call. = FALSE)
    }
    switch(type,
           link = parts$link,
           prob_reported = parts$p_report,
           true = parts$mu_true,
           zero = if (parts$zero_type == "inflated") parts$zero_prob else
             1 - parts$zero_prob,
           response = parts$response)
  }
}

# Components of a countreg() prediction for the rows of `data`: the linear
# predictor (with offset), the reporting probability, the distribution
# parameters on the scale used by get_probFunc(), the mean parameter
# exp(link), and the true mean (with the PLN/PLL mean factor).
.fcr_countreg_parts <- function(object, data) {
  model <- object$model
  coefs <- model$estimate
  family <- model$family
  params <- get_params(family)
  
  bp <- .fcr_blueprint(object, aux_formulas = list(
    dis1 = model$dis_param_formula_1, dis2 = model$dis_param_formula_2,
    underreport = model$underreport_formula, zero = model$zero_formula))
  fd <- forge_new_data(bp, data)
  n <- nrow(as.data.frame(data))
  
  beta_pred <- .fcr_coefs(coefs, colnames(fd$X))
  link <- as.vector(fd$X %*% beta_pred) + fd$offset
  
  p_report <- rep(1, n)
  X_und <- fd$aux$underreport
  if (!is.null(X_und)) {
    und_pars <- .fcr_coefs(coefs, paste0("Underreporting:", colnames(X_und)))
    und_lin <- as.vector(X_und %*% und_pars)
    und_family <- model$underreport_family
    if (identical(und_family, "probit") && is.null(model$underreport_sign)) {
      # Fits from flexCountReg 0.1.1 and earlier used P(reported) = 1 - pnorm(z'g)
      p_report <- stats::pnorm(und_lin, lower.tail = FALSE)
    } else {
      p_report <- .report_prob(und_lin, und_family)
    }
  }
  
  vec_1 <- .fcr_dist_param(coefs, params[[1]], fd$aux$dis1, n)
  vec_2 <- .fcr_dist_param(coefs, params[[2]], fd$aux$dis2, n)
  mu_link <- exp(link)
  # alpha and sigma are the parameters of the reported counts' distribution
  # (equal to the true counts' parameters without underreporting)
  thin <- .fcr_thin_params(family, vec_1, vec_2,
                           if (!is.null(X_und)) p_report)
  parts <- list(link = link, mu_link = mu_link, p_report = p_report,
                alpha = thin$alpha, sigma = thin$sigma,
                alpha_true = vec_1, sigma_true = vec_2, family = family,
                mu_true = .fcr_family_mean(family, mu_link, vec_1, vec_2),
                zero_type = model$zero_type)
  parts$response <- parts$mu_true * p_report
  
  # Zero-inflated or hurdle model
  X_zero <- fd$aux$zero
  if (!is.null(X_zero)) {
    zp <- stats::plogis(as.vector(
      X_zero %*% .fcr_coefs(coefs, paste0("Zero:", colnames(X_zero)))))
    parts$zero_prob <- zp
    if (model$zero_type == "inflated") {
      parts$response <- (1 - zp) * parts$mu_true
    } else {
      lf0 <- .fcr_parts_density(object, parts, 0, log = TRUE)
      parts$f0 <- exp(lf0)
      parts$lf0 <- lf0
      # 1 - P(Y = 0), accurate when P(Y = 0) is close to 1
      parts$response <- zp * parts$mu_true / -expm1(lf0)
    }
  }
  parts
}

# Count-family density P(Y = k) for every row of a .fcr_countreg_parts()
# result, evaluated as in the likelihood (before any zero model).
.fcr_parts_density <- function(object, parts, k, log = FALSE) {
  m <- object$model
  dens <- get_probFunc(parts$family)
  haltons <- randtoolbox::halton(m$ndraws %||% 1500)
  mean_par <- parts$mu_link * parts$p_report
  out <- rep(NA_real_, length(mean_par))
  ok <- !is.na(mean_par)
  if (any(ok)) {
    out[ok] <- dens(y = rep(k, sum(ok)), predicted = mean_par[ok],
                    alpha = if (!is.null(parts$alpha)) parts$alpha[ok],
                    sigma = if (!is.null(parts$sigma)) parts$sigma[ok],
                    haltons = haltons, normed_haltons = stats::qnorm(haltons),
                    engine = m$engine %||% "halton", log = log)
  }
  out
}

# ---- Internal helpers for predict() ------------------------------------------

# Training-data blueprint (terms, factor levels, contrasts, offsets). Fits from
# flexCountReg >= 0.1.2 store it; for older fits and for estimators that do
# not yet store it, it is rebuilt from the training data kept in the object.
.fcr_blueprint <- function(object, aux_formulas = list()) {
  model <- object$model
  if (!is.null(model$blueprint)) return(model$blueprint)
  if (is.null(object$data)) {
    stop("The model object has no stored blueprint or training data, so ",
         "predictions cannot be made. Refit the model.", call. = FALSE)
  }
  offset <- model[["offset"]]   # exact match: `$` would match offset_values
  if (.fcr_has_formula_offset(model$formula)) offset <- NULL
  prepare_model_data(model$formula, object$data,
                     aux_formulas = aux_formulas, offset = offset)$blueprint
}

# Coefficients by name, with a clear error if any are missing.
.fcr_coefs <- function(coefs, nm) {
  miss <- setdiff(nm, names(coefs))
  if (length(miss)) {
    stop("Coefficient(s) not found in the model: ", toString(miss),
         call. = FALSE)
  }
  as.vector(coefs[nm])
}

# Value of a distribution parameter (exp of its linear predictor).
.fcr_dist_param <- function(coefs, pname, X_dis, n) {
  if (is.null(pname)) return(NULL)
  if (is.null(X_dis)) return(rep(exp(.fcr_coefs(coefs, pname)), n))
  exp(as.vector(X_dis %*%
                  .fcr_coefs(coefs, paste0(pname, ":", colnames(X_dis)))))
}


# Panel (cluster) identifier for the rows of `data` in a countreg.rp model.
.fcr_rp_panel <- function(model, data) {
  pid <- model$panel_id
  if (is.null(pid)) return(seq_len(nrow(data)))
  if (is.character(pid) && length(pid) == 1L) {
    if (!pid %in% names(data)) {
      stop("method = \"Individual\" requires the panel id column `", pid,
           "` in `newdata`.", call. = FALSE)
    }
    return(data[[pid]])
  }
  if (length(pid) == nrow(data)) return(pid)
  stop("The model was fit with a `panel_id` vector; method = \"Individual\" ",
       "is only available for data with the same rows. Refit with ",
       "`panel_id` given as a column name.", call. = FALSE)
}

# Empirical Bayes (observation-specific) predictions for countreg.rp models,
# following Wood and Gayah (2025), Eqs. 5-7. Each simulation draw r is weighted
# by the likelihood of unit i's observed outcomes at that draw,
#   W_ir = L_ir / sum_r L_ir,  L_ir = prod_{t in i} f(y_it | mu_itr, theta),
# and the conditional coefficients are E[beta_i | y_i] = sum_r W_ir beta_ir.
#
# eb_type = "coefficients" (default, as in the paper) plugs E[beta_i | y_i]
# into the mean function; eb_type = "mean" returns E[mu_it | y_i] =
# sum_r W_ir mu_itr. The conditional coefficients and their conditional
# standard deviations are attached as attributes "eb_coefficients" and
# "eb_sd" (one row per observation, one column per random parameter).
.fcr_rp_eb_predict <- function(model, data, family, mu_mat, mu_fixed, X_rand,
                               draw_args, vec_1, vec_2, mean_adj, eb_type,
                               dist_haltons) {
  eb_type <- match.arg(eb_type, c("coefficients", "mean"))
  
  y_name <- all.vars(model$formula)[1]
  if (!y_name %in% names(data)) {
    stop("method = \"Individual\" (Empirical Bayes) requires the outcome `",
         y_name, "` in `newdata`.", call. = FALSE)
  }
  y <- data[[y_name]]
  panel <- .fcr_rp_panel(model, data)
  
  N <- nrow(mu_mat)
  D <- ncol(mu_mat)
  rpar <- colnames(X_rand)
  out <- rep(NA_real_, N)
  eb_coef <- matrix(NA_real_, N, length(rpar), dimnames = list(NULL, rpar))
  eb_sd <- eb_coef
  
  ok <- !is.na(y) & !is.na(panel) & stats::complete.cases(mu_mat)
  if (!is.null(vec_1)) ok <- ok & !is.na(vec_1)
  if (!is.null(vec_2)) ok <- ok & !is.na(vec_2)
  if (!any(ok)) {
    return(structure(out, eb_coefficients = eb_coef, eb_sd = eb_sd))
  }
  n_ok <- sum(ok)
  
  # Likelihood of each observation at each draw, as in estimation
  prob_fun <- get_probFunc(family)
  mu_ok <- mu_mat[ok, , drop = FALSE]
  probs <- prob_fun(y = rep(y[ok], times = D),
                    predicted = as.vector(mu_ok),
                    alpha = if (!is.null(vec_1)) rep(vec_1[ok], times = D),
                    sigma = if (!is.null(vec_2)) rep(vec_2[ok], times = D),
                    haltons = dist_haltons,
                    normed_haltons = stats::qnorm(dist_haltons),
                    engine = if (is.null(model$engine)) "halton"
                             else model$engine, log = TRUE)
  log_probs <- matrix(probs, nrow = n_ok, ncol = D)
  
  # Panel log-likelihood per draw, then normalized weights (log-sum-exp)
  pg <- as.character(panel[ok])
  log_L <- rowsum(log_probs, group = pg, reorder = FALSE)
  row_max <- apply(log_L, 1, max)
  W <- exp(log_L - row_max)
  W <- W / rowSums(W)
  bad <- !is.finite(row_max) | !stats::complete.cases(W)
  if (any(bad)) {
    warning(sum(bad), " unit(s) have zero likelihood at every draw; their ",
            "Empirical Bayes weights are set to be equal across draws.",
            call. = FALSE)
    W[bad, ] <- 1 / D
  }
  W_obs <- W[match(pg, rownames(log_L)), , drop = FALSE]
  
  # Conditional coefficients: draws of each random parameter, built with the
  # same function used in estimation (unit design column for parameter k)
  for (k in seq_along(rpar)) {
    E_k <- matrix(0, nrow = n_ok, ncol = length(rpar))
    E_k[, k] <- 1
    args_k <- draw_args
    if (!is.null(args_k$X_het_mean)) {
      args_k$X_het_mean <- args_k$X_het_mean[ok, , drop = FALSE]
    }
    if (!is.null(args_k$X_het_var)) {
      args_k$X_het_var <- args_k$X_het_var[ok, , drop = FALSE]
    }
    B_k <- do.call(generate_random_draws,
                   c(args_k, list(X_rand = E_k)))$xb_rand_mat
    m1 <- rowSums(W_obs * B_k)
    m2 <- rowSums(W_obs * B_k^2)
    eb_coef[ok, k] <- m1
    eb_sd[ok, k] <- sqrt(pmax(m2 - m1^2, 0))
  }
  
  a1 <- if (!is.null(vec_1)) vec_1[ok]
  a2 <- if (!is.null(vec_2)) vec_2[ok]
  if (eb_type == "coefficients") {
    mu_eb <- mu_fixed[ok] *
      exp(rowSums(X_rand[ok, , drop = FALSE] * eb_coef[ok, , drop = FALSE]))
    out[ok] <- .fcr_family_mean(family, mu_eb, a1, a2)
  } else {
    m <- .fcr_family_mean(family, as.vector(mu_ok), rep(a1, D),
                          if (!is.null(a2)) rep(a2, D))
    out[ok] <- rowSums(W_obs * matrix(m, nrow = n_ok))
  }
  
  structure(out, eb_coefficients = eb_coef, eb_sd = eb_sd)
}
