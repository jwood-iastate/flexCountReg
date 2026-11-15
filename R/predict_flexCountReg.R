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
#'   Types other than \code{"response"} are available for \code{countreg}
#'   models only.
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
#'   \item \strong{Exact}: Uses the analytical Moment Generating Functions
#'   (MGFs) of the random parameter distributions to calculate the exact
#'   expected value. This method is faster and removes simulation error.
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
#' @import randtoolbox stats modelr rlang 
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
                                          "prob_reported", "link")){
  
  type <- match.arg(type)
  
  # Extract optional parameters from '...'
  additional_args  <- list(...)
  
  if (!is.null(additional_args$new_data) && is.null(newdata)) {
    newdata <- additional_args$new_data
  }
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
    
    # --- Matrices, built with the training terms and factor levels ---
    bp <- .fcr_blueprint(object, aux_formulas = list(
      dis1 = model$dis_param_formula_1, dis2 = model$dis_param_formula_2,
      het_mean = model$het_mean_formula, het_var = model$het_var_formula))
    fd <- forge_new_data(bp, data)
    X_Fixed <- fd$X
    X_offset <- fd$offset
    
    # Random-parameter matrix; same intercept rule as in estimation
    train <- as.data.frame(object$data)
    rpar_tt <- stats::delete.response(stats::terms(model$rpar_formula))
    rpar_xlev <- stats::.getXlevels(
      rpar_tt, stats::model.frame(rpar_tt, train, na.action = stats::na.omit))
    X_rand <- .fcr_rhs_matrix(model$rpar_formula, data, xlev = rpar_xlev)
    if("(Intercept)" %in% colnames(X_rand) && !is.null(rpardists)){
      if(!any(grepl("intercept", names(rpardists), ignore.case=TRUE))){
        X_rand <- X_rand[ , colnames(X_rand) != "(Intercept)", drop=FALSE]
      }
    }
    rpar <- colnames(X_rand)
    N_rand <- ncol(X_rand)
    
    # --- Coefficients, matched by name ---
    fixed_coefs <- .fcr_coefs(coefs, colnames(X_Fixed))
    random_coefs_means <- .fcr_coefs(coefs, paste0(rpar, ":Mean"))
    if (correlated) {
      rand_var_params <- coefs[grepl("^Chol:", names(coefs))]
    } else {
      rand_var_params <- .fcr_coefs(coefs, paste0(rpar, ":St.Dev"))
    }
    
    X_het_mean <- fd$aux$het_mean
    X_het_var <- fd$aux$het_var
    het_mean_coefs <- if (!is.null(X_het_mean)) {
      .fcr_coefs(coefs, paste0("HetMean:", colnames(X_het_mean)))
    }
    het_var_coefs <- if (!is.null(X_het_var)) {
      .fcr_coefs(coefs, paste0("HetVar:", colnames(X_het_var)))
    }
    
    # Distribution parameters (for the PLN/PLL mean correction)
    vec_1 <- .fcr_dist_param(coefs, params[[1]], fd$aux$dis1, nrow(data))
    vec_2 <- .fcr_dist_param(coefs, params[[2]], fd$aux$dis2, nrow(data))
    mean_adj <- .fcr_mean_adjust(family, vec_1, vec_2)
    
    # --- Prediction Method Implementation ---
    
    if (method == 'Exact') {
      # Base prediction
      pred_exact <- exp(as.vector(X_Fixed %*% fixed_coefs) + X_offset)
      
      # Correlated Normal
      if (correlated) {
        L <- matrix(0, N_rand, N_rand)
        L[lower.tri(L, diag = TRUE)] <- rand_var_params
        var_rand <- rowSums((X_rand %*% L)^2)
        mu_rand <- as.vector(X_rand %*% random_coefs_means)
        pred_exact <- pred_exact * exp(mu_rand + 0.5 * var_rand)
        
      } else {
        # Independent
        for (i in 1:N_rand) {
          dist_type <- rpardists[i]
          x_col <- X_rand[, i]
          mu <- random_coefs_means[i]
          sigma <- rand_var_params[i]
          
          correction <- rep(1, nrow(data))
          
          if (dist_type == "n") { # Normal
            correction <- exp(x_col * mu + 0.5 * (x_col^2) * (sigma^2))
          } else if (dist_type == "ln") { # Lognormal
            arg_w <- -exp(mu) * x_col * sigma^2
            valid_idx <- arg_w >= -1/exp(1)
            if (any(valid_idx)) {
              x_val <- x_col[valid_idx]
              w_val <- lamW::lambertW0(arg_w[valid_idx])
              num <- x_val * exp(mu - w_val) - (w_val^2)/(2 * sigma^2)
              den <- sigma * sqrt(abs(x_val * exp(mu - w_val) - (1/sigma^2)))
              correction[valid_idx] <- exp(num) / den
            }
            if(any(!valid_idx)) correction[!valid_idx] <- NA
          } else if (dist_type == "t") { # Triangular
            nz <- abs(x_col) > 1e-8
            if (any(nz)) {
              x_nz <- x_col[nz]
              term <- exp(x_nz * (mu - sigma))
              correction[nz] <- 
                (exp(2 * sigma * x_nz) - 2 * exp(sigma * x_nz) + 1) * term / 
                (sigma^2 * x_nz^2)
            }
          } else if (dist_type == "u") { # Uniform
            nz <- abs(x_col) > 1e-8
            if (any(nz)) {
              x_nz <- x_col[nz]
              correction[nz] <- 
                (exp(x_nz * (mu + sigma)) - exp(x_nz * (mu - sigma))) / 
                (2 * sigma * x_nz)
            }
          } else if (dist_type == "g") { # Gamma
            term_base <- 1 - (sigma^2 * x_col) / mu
            if(any(term_base <= 0)) 
              warning("Gamma exact prediction undefined for some observations.")
            correction <- 
              ifelse(term_base > 0, term_base^(-(mu^2)/(sigma^2)), NA)
          }
          pred_exact <- pred_exact * correction
        }
      }
      
      # PLN/PLL mean correction
      pred_exact <- pred_exact * mean_adj
      
      return(pred_exact)
      
    } else {
      # Simulated or Individual (Empirical Bayes)
      hdraws <- 
        as.matrix(randtoolbox::halton(ndraws, N_rand, mixed = scrambled))
      
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
        correlated = correlated)
      
      draws_info <- do.call(generate_random_draws,
                            c(draw_args, list(X_rand = X_rand)))
      
      rpar_mat <- exp(draws_info$xb_rand_mat)
      mu_fixed <- exp(as.vector(X_Fixed %*% fixed_coefs) + X_offset)
      # Conditional means (N x draws), as in the likelihood
      mu_mat <- sweep(rpar_mat, 1, mu_fixed, "*")
      
      if (method == 'Simulated') {
        return(as.vector(rowMeans(mu_mat)) * mean_adj)
      }
      
      return(.fcr_rp_eb_predict(
        model = model, data = data, family = family, mu_mat = mu_mat,
        mu_fixed = mu_fixed, X_rand = X_rand, draw_args = draw_args,
        vec_1 = vec_1, vec_2 = vec_2, mean_adj = mean_adj,
        eb_type = if (is.null(additional_args$eb_type)) "coefficients"
                  else additional_args$eb_type))
    }
  }
  
  # =========================================================================
  # 2a. NEGATIVE MULTINOMIAL FAMILY (negmultinom)
  # =========================================================================
  else if (modtype == "NegMultinom") {
    return(.nm_predict(model, data, method = method))
  }
  
  # =========================================================================
  # 2b. RANDOM EFFECTS MODELS (poisLindRE / RENB)
  # =========================================================================
  else if (modtype %in% c("poisLindRE", "RENB", "REPoisson")) {
    
    bp <- .fcr_blueprint(object)
    fd <- forge_new_data(bp, data)
    beta_pred <- model$beta_pred
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
    
    coefs <- model$estimate
    family <- model$family
    params <- get_params(family)
    
    bp <- .fcr_blueprint(object, aux_formulas = list(
      dis1 = model$dis_param_formula_1, dis2 = model$dis_param_formula_2,
      underreport = model$underreport_formula))
    fd <- forge_new_data(bp, data)
    
    # Coefficients matched by name
    beta_pred <- .fcr_coefs(coefs, colnames(fd$X))
    link <- as.vector(fd$X %*% beta_pred) + fd$offset
    if (type == "link") return(link)
    
    # Reporting probability
    p_report <- rep(1, nrow(data))
    X_und <- fd$aux$underreport
    if (!is.null(X_und)) {
      und_pars <- .fcr_coefs(coefs, paste0("Underreporting:", colnames(X_und)))
      und_lin <- as.vector(X_und %*% und_pars)
      und_family <- model$underreport_family
      if (identical(und_family, "probit") && is.null(model$underreport_sign)) {
        # Fits from flexCountReg <= 0.1.2 used P(reported) = 1 - pnorm(z'g)
        p_report <- stats::pnorm(und_lin, lower.tail = FALSE)
      } else {
        p_report <- .report_prob(und_lin, und_family)
      }
    }
    if (type == "prob_reported") return(p_report)
    
    # True mean, with the mean correction for PLN and PLL
    vec_1 <- .fcr_dist_param(coefs, params[[1]], fd$aux$dis1, nrow(data))
    vec_2 <- .fcr_dist_param(coefs, params[[2]], fd$aux$dis2, nrow(data))
    mu_true <- exp(link) * .fcr_mean_adjust(family, vec_1, vec_2)
    if (type == "true") return(mu_true)
    
    return(mu_true * p_report)
  }
}

# ---- Internal helpers for predict() ------------------------------------------

# Training-data blueprint (terms, factor levels, contrasts, offsets). Fits from
# flexCountReg >= 0.1.3 store it; for older fits and for estimators that do
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

# Ratio of the distribution mean to exp(X*beta). For the Poisson-lognormal
# (sigma = parameter 1) and Poisson-Lindley-lognormal (sigma = parameter 2),
# the lognormal factor has mean exp(sigma^2 / 2).
.fcr_mean_adjust <- function(family, vec_1, vec_2) {
  if (family == "PLN" && !is.null(vec_1)) return(exp(vec_1^2 / 2))
  if (family == "PLL" && !is.null(vec_2)) return(exp(vec_2^2 / 2))
  1
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
                               draw_args, vec_1, vec_2, mean_adj, eb_type) {
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
  dist_haltons <- randtoolbox::halton(model$ndraws, 1)
  mu_ok <- mu_mat[ok, , drop = FALSE]
  probs <- prob_fun(y = rep(y[ok], times = D),
                    predicted = as.vector(mu_ok),
                    alpha = if (!is.null(vec_1)) rep(vec_1[ok], times = D),
                    sigma = if (!is.null(vec_2)) rep(vec_2[ok], times = D),
                    haltons = dist_haltons,
                    normed_haltons = stats::qnorm(dist_haltons),
                    engine = if (is.null(model$engine)) "halton"
                             else model$engine)
  log_probs <- matrix(log(probs), nrow = n_ok, ncol = D)
  
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
  
  adj <- rep_len(mean_adj, N)[ok]
  if (eb_type == "coefficients") {
    out[ok] <- mu_fixed[ok] *
      exp(rowSums(X_rand[ok, , drop = FALSE] * eb_coef[ok, , drop = FALSE])) *
      adj
  } else {
    out[ok] <- rowSums(W_obs * mu_ok) * adj
  }
  
  structure(out, eb_coefficients = eb_coef, eb_sd = eb_sd)
}
