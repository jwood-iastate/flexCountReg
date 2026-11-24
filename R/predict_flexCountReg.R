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
#' @param interval \code{"none"} (default), \code{"confidence"} or
#'   \code{"prediction"} (\code{countreg} models). See the section
#'   \emph{Intervals}.
#' @param level the confidence or prediction level of the interval.
#'
#' @section Intervals:
#' For \code{countreg} models, \code{interval = "confidence"} gives a
#' confidence interval for the expected count (\code{type = "response"}, the
#' expected observed count, or \code{"true"}), and
#' \code{interval = "prediction"} a prediction interval for one future
#' observed count (\code{type = "response"}). The result is then a data
#' frame with columns \code{fit}, \code{lwr} and \code{upr}, with the
#' settings in the attribute \code{"interval"}. For content tolerance
#' intervals see \code{\link{tolerance_interval}}.
#'
#' Estimation uncertainty is represented by parameter vectors chosen by the
#' optional argument \code{uncertainty} (in \code{...}):
#' \itemize{
#'   \item \code{"normal"} (default): the asymptotic normal distribution of
#'   the estimates, \eqn{\hat\theta + R^\top z} with
#'   \eqn{R^\top R} the covariance matrix and \eqn{z} the standard normal
#'   quantiles of \code{ndraws} (default 500) Halton points. Halton
#'   sequences are deterministic, so the intervals are reproducible and have
#'   less simulation error than pseudo-random draws; for models with many
#'   parameters, \code{scrambled = TRUE} (with a \code{seed}) uses scrambled
#'   Halton points, which avoid the correlation between high dimensions of
#'   the standard sequence. When the covariance matrix is not positive
#'   definite (e.g. a dispersion parameter at its boundary) the parametric
#'   bootstrap is used instead, with a message.
#'   \item \code{"bootstrap"}: a parametric bootstrap. \code{bootstraps}
#'   (default 200) response vectors are simulated from the fitted model on
#'   the training rows and the model is refit to each with its frozen
#'   specification; refits that fail or do not converge are rejected and
#'   counted. Use \code{seed} for reproducibility. Not available for models
#'   with weights.
#' }
#' The confidence interval takes the \eqn{(1 - level)/2} and
#' \eqn{(1 + level)/2} quantiles of the expected count over these parameter
#' vectors. The prediction interval uses the predictive distribution of the
#' count, the probability mass function of the observed count averaged over
#' the parameter vectors, and returns its integer quantiles: the smallest
#' counts \eqn{L} and \eqn{U} with \eqn{P(Y \le L) \ge (1 - level)/2} and
#' \eqn{P(Y \le U) \ge (1 + level)/2}. Because counts are discrete, the
#' predictive coverage is usually above \code{level}; the limits are not
#' interpolated. Both include zero-inflated, hurdle and underreporting
#' components. \code{max_count} (in \code{...}, default 1e5) bounds the
#' counts used to resolve the upper tail.
#'
#' In a Monte Carlo study of Poisson and NB2 regressions (n = 30 and 200),
#' prediction intervals at level 0.90 covered 0.93 to 0.97 of future counts,
#' while plug-in intervals that ignore estimation uncertainty fell to 0.88
#' for the smaller sample.
#'
#' @details
#' The returned vector always has one element per row of \code{newdata}, in the
#' same order. Rows with missing values in any variable the estimate needs
#' get \code{NA}: a row with a missing distribution-parameter, zero-model or
#' reporting covariate has no count distribution, so its probabilities
#' (\code{type = "prob"}) and intervals are \code{NA}, while its expected
#' count is returned when it does not depend on that parameter (for example
#' the NB2 mean) and is \code{NA} otherwise. Factor levels and contrasts from the training data are used, so
#' \code{newdata} may contain only a subset of the levels. An error is raised
#' if an offset column is missing from \code{newdata}.
#' @note optional parameter `method`: Only valid for random parameters models
#'   (`countreg_rp`). Options include \code{Simulated} (default),
#'   \code{Individual} (or \code{EB}), or \code{Exact}.
#' @note optional parameter `eb_type`: for \code{method = "Individual"},
#'   either \code{"coefficients"} (default) or \code{"mean"}.
#'
#' @description
#' Generates predictions for the expected count (lambda) for observations.
#' 
#' For \strong{countreg_rp} (Random Parameters) models, three methods are
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
#'   lognormal parameters with \eqn{x < 0} the lognormal MGF has no closed
#'   form, and the Laplace approximation of Asmussen, Jensen and
#'   Rojas-Nandayapa (2016) is used: accurate, but approximate despite the
#'   method's name (it needs the lamW package). For \eqn{x > 0} with a
#'   lognormal parameter (and for gamma parameters with
#'   \eqn{x\sigma^2/\mu \ge 1}) the expected value is infinite, and
#'   \code{Inf} is returned with a warning, for both \code{Exact} and
#'   \code{Simulated} predictions. Whether the expected value exists is
#'   decided from the parameters alone, so \code{Simulated} predictions do
#'   not need lamW. A finite expected value too large to represent is
#'   reported as an overflow.
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
#'  used in estimating the model are provided. With \code{interval =
#'  "confidence"} or \code{"prediction"}, a data frame with columns
#'  \code{fit}, \code{lwr} and \code{upr} (see the section Intervals).
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
#' # 2. Random Parameters Model (countreg_rp)
#' # =========================================================================
#' rp_nb2 <- countreg_rp(Total_crashes ~ lnaadt + lnlength,
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
                                 new_data = NULL,
                                 interval = c("none", "confidence",
                                              "prediction"),
                                 level = 0.95){

  type <- match.arg(type)
  interval <- match.arg(interval)
  
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
  # Correlated random effects: the group-mean columns of `newdata`
  if (!is.null(model$cre)) data <- .fcr_cre_newdata(model$cre, data)
  
  if (!is.null(object$modelType)) {
    modtype <- object$modelType
  } else if (!is.null(model$modelType)) {
    modtype <- model$modelType
  } else {
    modtype <- "countreg" 
  }
  
  if (modtype == "countreg_fe") {
    if (type != "response" || interval != "none") {
      stop("countreg_fe() models give type = \"response\" predictions ",
           "without intervals.", call. = FALSE)
    }
    return(.fcr_fe_predict(object, data,
                           additional_args$newgroup %||% "none"))
  }

  if (interval != "none") {
    if (modtype != "countreg") {
      stop("Intervals are available for countreg() models.", call. = FALSE)
    }
    ia <- additional_args
    return(.fcr_predict_interval(
      object, data, type, interval, level,
      uncertainty = match.arg(ia$uncertainty %||% "normal",
                              c("normal", "bootstrap")),
      ndraws = ia$ndraws %||% 500, bootstraps = ia$bootstraps %||% 200,
      seed = ia$seed, scrambled = isTRUE(ia$scrambled),
      max_count = ia$max_count %||% 1e5))
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
  # 1. RANDOM PARAMETERS MODELS (countreg_rp)
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
    # Scrambled draws of fits before draw_version 3 ignored their seed
    draw_seed <- if (isTRUE(model$draw_version >= 3)) model$draw_seed
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
    
    # --- Does the expected value exist? (no MGF evaluation needed) ---
    # --- Exact: moment generating functions (Wood and Gayah, 2025) ---
    if (method %in% c("Exact", "Simulated")) {
      loc <- .fcr_rp_locations(random_coefs_means, rpardists, nrow(data),
                               het_mean_coefs, X_het_mean, het_var_coefs,
                               X_het_var, het_form)
      sd_k <- function(k) abs(rand_var_params[k]) * loc$var_multiplier
      no_mean <- rep(FALSE, nrow(data))
      if (!correlated) {
        for (k in seq_len(N_rand)) {
          ex <- .fcr_rp_mean_exists(rpardists[[k]], X_rand[, k],
                                    loc$mu_mat[, k], sd_k(k))
          no_mean <- no_mean | (!is.na(ex) & !ex)
        }
      }
      if (any(no_mean)) {
        warning(sum(no_mean), " prediction(s) are infinite: the expected ",
                "value does not exist for a lognormal random parameter with ",
                "a positive x, or a gamma random parameter with ",
                "x * sigma^2 / mu >= 1. See ?countreg_rp.", call. = FALSE)
      }
      if (method == "Exact") {
        if (correlated) {
          L <- .fcr_chol_unpack(rand_var_params, N_rand)
          W <- X_rand * loc$var_multiplier
          log_rp <- rowSums(X_rand * loc$mu_mat) + 0.5 * rowSums((W %*% L)^2)
        } else {
          log_rp <- rep(0, nrow(data))
          for (k in seq_len(N_rand)) {
            log_rp <- log_rp + rp_log_mgf(rpardists[[k]], X_rand[, k],
                                          loc$mu_mat[, k], sd_k(k))
          }
        }
        out <- exp(log(mu_fixed) + log_rp + log(mean_adj))
        overflow <- is.infinite(out) & !no_mean
        if (any(overflow)) {
          warning(sum(overflow), " Exact prediction(s) exist but exceed the ",
                  "largest representable number (overflow); they are Inf.",
                  call. = FALSE)
        }
        out[no_mean] <- Inf
        return(as.vector(out))
      }
    }
    
    # --- Simulated or Individual (Empirical Bayes) ---
    # Empirical Bayes for fits with per-unit draws regenerates each unit's
    # own draws (the ones used in estimation); marginal (Simulated)
    # predictions use one set of draws for every row
    hdraws <- if (method == "Individual" &&
                  identical(model$draws, "per_unit")) {
      .fcr_unit_halton(.fcr_rp_draw_units(model, data), model$ndraws,
                       N_rand, isTRUE(scrambled), draw_seed)
    } else {
      .fcr_halton(ndraws, N_rand + 1, isTRUE(scrambled),
                  draw_seed)[, seq_len(N_rand), drop = FALSE]
    }
    
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
                  isTRUE(scrambled), draw_seed)[, N_rand + 1]
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
  fd <- .fcr_countreg_fd(object, data)
  .fcr_parts_from_lp(object, .fcr_countreg_lp(object, fd))
}

# Design matrices of all model components for the rows of `data`
.fcr_countreg_fd <- function(object, data) {
  model <- object$model
  bp <- .fcr_blueprint(object, aux_formulas = list(
    dis1 = model$dis_param_formula_1, dis2 = model$dis_param_formula_2,
    underreport = model$underreport_formula, zero = model$zero_formula))
  forge_new_data(bp, data)
}

# Linear predictors of the model components at the coefficient vector
# `coefs` (default: the estimates): mean (with offset), dis1 and dis2 (the
# logarithms of the distribution parameters), und (reporting) and zero.
# Absent components are NULL.
.fcr_countreg_lp <- function(object, fd, coefs = object$model$estimate) {
  params <- get_params(object$model$family)
  n <- nrow(fd$X)
  lin <- function(X, prefix) {
    if (is.null(X)) return(NULL)
    as.vector(X %*% .fcr_coefs(coefs, paste0(prefix, colnames(X))))
  }
  dis <- function(pname, X) {
    if (is.null(pname)) return(NULL)
    if (is.null(X)) return(rep(.fcr_coefs(coefs, pname), n))
    lin(X, paste0(pname, ":"))
  }
  list(mean = as.vector(fd$X %*% .fcr_coefs(coefs, colnames(fd$X))) +
         fd$offset,
       dis1 = dis(params[[1]], fd$aux$dis1),
       dis2 = dis(params[[2]], fd$aux$dis2),
       und = lin(fd$aux$underreport, "Underreporting:"),
       zero = lin(fd$aux$zero, "Zero:"))
}

# Prediction components from linear predictors (see .fcr_countreg_lp()),
# one element per row; the rows may be data rows or points of a parameter
# set.
.fcr_parts_from_lp <- function(object, lp) {
  model <- object$model
  family <- model$family
  link <- lp$mean
  n <- length(link)

  p_report <- rep(1, n)
  if (!is.null(lp$und)) {
    und_lin <- lp$und
    und_family <- model$underreport_family
    if (identical(und_family, "probit") && is.null(model$underreport_sign)) {
      # Fits from flexCountReg 0.1.1 and earlier used P(reported) = 1 - pnorm(z'g)
      p_report <- stats::pnorm(und_lin, lower.tail = FALSE)
    } else {
      p_report <- .report_prob(und_lin, und_family)
    }
  }
  
  vec_1 <- if (!is.null(lp$dis1)) exp(lp$dis1)
  vec_2 <- if (!is.null(lp$dis2)) exp(lp$dis2)
  mu_link <- exp(link)
  # alpha and sigma are the parameters of the reported counts' distribution
  # (equal to the true counts' parameters without underreporting)
  thin <- .fcr_thin_params(family, vec_1, vec_2,
                           if (!is.null(lp$und)) p_report)
  parts <- list(link = link, mu_link = mu_link, p_report = p_report,
                alpha = thin$alpha, sigma = thin$sigma,
                alpha_true = vec_1, sigma_true = vec_2, family = family,
                mu_true = .fcr_family_mean(family, mu_link, vec_1, vec_2),
                zero_type = model$zero_type)
  parts$response <- parts$mu_true * p_report
  
  # Zero-inflated or hurdle model
  if (!is.null(lp$zero)) {
    zp <- stats::plogis(lp$zero)
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


# Panel (cluster) identifier for the rows of `data` in a countreg_rp model.
# Draw unit of each row of `data` for a fit with per-unit draws: units of
# the estimation sample get their own draws back (matched by panel id, or by
# row position without panels); new units get new draws, in order of
# appearance. Rows without a panel id get unit 1 (they are not predicted).
.fcr_rp_draw_units <- function(model, data) {
  panel <- as.character(.fcr_rp_panel(model, data))
  u <- match(panel, model$draw_units)
  new <- is.na(u) & !is.na(panel)
  if (any(new)) {
    u[new] <- length(model$draw_units) + match(panel[new], unique(panel[new]))
  }
  u[is.na(u)] <- 1L
  u
}

.fcr_rp_panel <- function(model, data) {
  pid <- model$panel_id
  if (is.null(pid)) return(seq_len(nrow(data)))
  if (is.character(pid) && length(pid) == 1L) {
    if (!pid %in% names(data)) {
      if (isTRUE(model$group_from_vector)) {
        stop("method = \"Individual\" needs each row's panel id. The model ",
             "was fit with `group_var` given as a vector, which is stored ",
             "in the training data as the column `", pid, "`. Add that ",
             "column to `newdata`, or refit with `group_var` naming a ",
             "column of `data`.", call. = FALSE)
      }
      stop("method = \"Individual\" requires the panel id column `", pid,
           "` in `newdata`.", call. = FALSE)
    }
    return(data[[pid]])
  }
  # Fits from flexCountReg 0.1.2 and earlier stored the vector itself; it
  # can only be matched to rows by position
  if (length(pid) == nrow(data)) {
    warning("The model was fit with `group_var` given as a vector, so ",
            "panels are matched to the rows of `newdata` by position. This ",
            "is only correct for the training rows in their original ",
            "order; refit to match by id.", call. = FALSE)
    return(pid)
  }
  stop("The model was fit with a `panel_id` vector; method = \"Individual\" ",
       "is only available for data with the same rows. Refit with ",
       "`panel_id` given as a column name.", call. = FALSE)
}

# Empirical Bayes (observation-specific) predictions for countreg_rp models,
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
  
  # The response as in estimation: the formula's left-hand side evaluated on
  # the data (so I(count / 2) ~ x uses count / 2, not the raw `count`)
  lhs <- model$formula[[2]]
  miss <- setdiff(all.vars(lhs), names(data))
  if (length(miss)) {
    stop("method = \"Individual\" (Empirical Bayes) requires the outcome `",
         paste(miss, collapse = "`, `"), "` in `newdata`.", call. = FALSE)
  }
  y <- eval(lhs, data, environment(model$formula) %||% parent.frame())
  if (is.logical(y) && all(is.na(y))) y <- as.numeric(y)   # all missing
  if (!is.numeric(y) || NCOL(y) != 1L || length(y) != nrow(data)) {
    stop("The model's response must evaluate to one numeric value per row ",
         "of `newdata`.", call. = FALSE)
  }
  y <- as.vector(y)
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
    args_k$hdraws <- .fcr_draw_rows(args_k$hdraws, which(ok))
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
