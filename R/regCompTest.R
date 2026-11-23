.mcfadden_r2 <- function(ll_full, ll_null) {
  if (is.finite(ll_full) && is.finite(ll_null) && ll_null != 0) {
    1 - (ll_full / ll_null)
  } else {
    NA_real_
  }
}

# Training rows of a fitted model, rebuilt with its blueprint: outcome,
# design matrix, offset and weights for the rows used in estimation.
.fcr_fit_pieces <- function(object, data = NULL) {
  m <- object$model
  if (is.null(data)) {
    d <- .fcr_training_rows(object)
  } else {
    d <- as.data.frame(data)
    if (!is.null(object$data) && nrow(d) == nrow(as.data.frame(object$data)) &&
        !is.null(m[["rows_used"]])) {
      d <- d[m[["rows_used"]], , drop = FALSE]
    }
  }
  bp <- m[["blueprint"]] %||% .fcr_blueprint(object)
  fd <- forge_new_data(bp, d)
  y <- as.numeric(stats::model.response(stats::model.frame(
    object$formula, d, na.action = stats::na.pass)))
  ok <- stats::complete.cases(fd$X) & !is.na(fd$offset) & !is.na(y)
  w <- m[["weights"]]
  if (length(w) != length(y)) w <- rep(1, length(y))
  list(y = y[ok], X = fd$X[ok, , drop = FALSE], offset = fd$offset[ok],
       weights = w[ok], data = d[ok, , drop = FALSE],
       offset_cols = bp$offset_cols)
}

# Poisson GLM log-likelihood on given pieces (intercept-only if X is NULL).
.fcr_poisson_ll <- function(pc, X = NULL) {
  if (is.null(X)) X <- matrix(1, length(pc$y), 1)
  f <- stats::glm.fit(X, pc$y, weights = pc$weights, offset = pc$offset,
                      family = stats::poisson())
  list(LL = sum(pc$weights * stats::dpois(pc$y, f$fitted.values, log = TRUE)),
       converged = isTRUE(f$converged),
       n.coef = ncol(X))
}

# countreg() fit of `family` on the pieces: intercept-only or with the
# model's mean variables. Offsets are passed as one summed column.
.fcr_countreg_base <- function(pc, formula, family, intercept_only, ...) {
  d <- pc$data
  d$.fcr_y <- pc$y
  d$.fcr_off <- pc$offset
  tt <- stats::delete.response(stats::terms(formula))
  rhs <- if (intercept_only) "1" else
    paste(attr(tt, "term.labels"), collapse = " + ")
  if (!nzchar(rhs)) rhs <- "1"
  # With the model's variables, keep its intercept choice: a model without
  # an intercept is compared with a base model without one (pc$X, the
  # Poisson base's design, has no intercept column either)
  if (!intercept_only && attr(tt, "intercept") == 0L) {
    rhs <- paste(rhs, "- 1")
  }
  # In the model formula's environment, so functions it uses (such as a
  # user's transformation) are found as in the original fit
  f <- stats::as.formula(paste(".fcr_y ~", rhs),
                         env = environment(formula) %||% parent.frame())
  args <- list(...)
  if (is.null(args$method)) args$method <- "BFGS"
  fit <- suppressWarnings(do.call(countreg, c(list(
    formula = f, data = d, family = family,
    offset = if (any(pc$offset != 0)) ".fcr_off" else NULL,
    weights = pc$weights, keep_data = FALSE), args)))
  list(LL = fit$model$maximum, n.coef = length(fit$model$estimate),
       converged = isTRUE(fit$model$converged))
}

# Log-likelihood of the baseline model for McFadden's pseudo R^2, fit on the
# same rows, offset and weights as the model. By default this is an
# intercept-only Poisson model for every model, so the pseudo R^2 is
# comparable across count distributions. `basemodel` can name another
# countreg() family, and `variables = TRUE` includes the model's mean
# covariates in the baseline.
.fit_null_ll <- function(trained_model, pc, basemodel = "Poisson",
                         variables = FALSE, ...) {
  if (identical(toupper(basemodel), "POISSON")) {
    return(tryCatch(.fcr_poisson_ll(pc, if (variables) pc$X else NULL)$LL,
                    error = function(e) NA_real_))
  }
  tryCatch(.fcr_countreg_base(pc, trained_model$formula, basemodel,
                              intercept_only = !variables, ...)$LL,
           error = function(e) NA_real_)
}

get_model_components <- function(trained_model, data, basemodel = "Poisson",
                                 variables = FALSE, ...) {
  if (!inherits(trained_model, "flexCountReg")) {
    stop("The supplied object is not a flexCountReg model.", call. = FALSE)
  }
  pc <- .fcr_fit_pieces(trained_model, data)
  LL <- as.numeric(stats::logLik(trained_model))
  n.coef <- length(stats::coef(trained_model))
  LL_null <- .fit_null_ll(trained_model, pc, basemodel = basemodel,
                         variables = variables, ...)
  list(
    fit = trained_model$model,
    object = trained_model,
    pieces = pc,
    y = pc$y,
    LL = LL,
    LL_null = LL_null,
    PseudoR2 = .mcfadden_r2(LL, LL_null),
    n.coef = n.coef,
    nobs = stats::nobs(trained_model),
    n_bic = attr(stats::logLik(trained_model), "nobs"),
    AIC = stats::AIC(trained_model),
    BIC = stats::BIC(trained_model)
  )
}

#' Compare Regression Models with Likelihood Ratio Test, AIC, and BIC
#'
#' This function compares fitted regression model objects using the
#' Likelihood Ratio (LR) test, Akaike Information Criterion (AIC), and Bayesian
#' Information Criterion (BIC).
#'
#' When only one fitted model is supplied, the function estimates a base model
#' on the same rows, offset and weights and compares the fitted model against
#' it. When two fitted models are supplied, the function compares them
#' directly; the model with fewer parameters is the restricted (base) model of
#' the LR test.
#'
#' The LR test is reported only when the restricted model is nested in the
#' larger one with a regular restriction: the same estimator, zero and
#' underreporting components and random parameters, the restricted model's
#' parameters a subset of the larger model's, and the same family or a
#' nested family (Poisson in NB1, NB2, PLN, PIG1, PIG2, PIG, GP1, GP2 and
#' COM; NB1 and NB2 in NBP; PL in PLG and PLL). When the restriction puts
#' one parameter on the boundary of its space (e.g. Poisson in NB2), the
#' p-value uses a 50:50 mixture of chi-square distributions with df - 1 and
#' df degrees of freedom. Two fitted models must use the same observations,
#' outcome, offsets and weights. Otherwise the LR statistic and p-value are
#' \code{NA} (\code{LR_type} gives the kind of restriction); AIC and BIC are
#' always reported.
#'
#' @name regCompTest
#' @param model A fitted regression model object.
#' @param data An optional data frame containing the variables in the model.
#'   If not supplied, the original data used to estimate the model will be
#'   used when available.
#' @param basemodel A character string specifying the family of base model to
#'   compare against (options include the family from
#'   \code{\link{countreg}} or "Poisson"). Default is "Poisson". It is
#'   used for the LR test of a single model and, in both workflows, as the
#'   baseline of McFadden's Pseudo R^2.
#' @param variables Logical. If \code{TRUE}, the base model will include the
#'   same variables as the provided model. If \code{FALSE}, the base model
#'   will be an intercept-only model. Default is \code{FALSE}.
#' @param model2 An optional second fitted regression model object. If
#'   supplied, the function compares \code{model} and \code{model2} directly
#'   and does not estimate a separate base model.
#' @param print Logical. If \code{TRUE}, a table of the results will be shown.
#'   If \code{FALSE}, the table of results will not be printed to the console.
#' @param ... Additional arguments to be passed to the base model fitting
#'   function. These arguments are used only when a single fitted model is
#'   supplied and the function needs to estimate the base model.
#' @returns A list containing the comparison results.
#' 
#' For the single-model workflow, the returned list includes \code{LL},
#' \code{LLbase}, \code{LR}, \code{LRdof}, \code{AIC}, \code{AICbase},
#' \code{BIC}, \code{BICbase}, \code{LR_pvalue}, \code{PseudoR2},
#' \code{statistics}, \code{gtTable}, \code{latexTable}, and \code{htmlTable}.
#' 
#' For the two-model workflow, the returned list includes \code{model1},
#' \code{model2}, \code{base_model}, \code{comparison_model},
#' \code{base_index}, \code{comparison_index}, \code{LL_base},
#' \code{LL_comparison}, \code{LR}, \code{LRdof}, \code{LR_pvalue},
#' \code{statistics}, \code{gtTable}, \code{latexTable}, and \code{htmlTable}.
#'
#' @include metrics.R
#' @importFrom knitr kable
#' @importFrom tibble tibble
#' @importFrom gt gt tab_header fmt_number
#'
#' @details The function performs the following steps for the single-model
#' workflow:
#' \enumerate{
#' \item Fits the base model, either a Poisson regression or another specified
#'   model.
#' \item Computes the log-likelihoods of the provided model and the base
#'   model.
#' \item Calculates the AIC and BIC for both models.
#' \item Conducts a Likelihood Ratio test to compare the fitted model against
#'   the estimated base model when the fitted model has more parameters than
#'   the base model.
#' \item Computes McFadden's Pseudo R^2, \eqn{1 - LL/LL_0}, where
#'   \eqn{LL_0} is the log-likelihood of the base model (by default an
#'   intercept-only Poisson model, the base count model, so that the pseudo
#'   R^2 is comparable across models and distributions).
#' }
#'
#' For the two-model workflow, the model with fewer parameters is the
#' restricted (base) model and the other is the comparison (full) model. The
#' models must be fit to the same observations (an error is raised otherwise)
#' and should be nested; a warning is given when the restricted model has the
#' higher log-likelihood (a sign that the models are not nested or that one
#' did not converge), and no LR test is reported when they have the same
#' number of parameters. The LR test assumes the restricted model is in the
#' interior of the parameter space; when it is on the boundary (for example
#' Poisson within NB2), the p-value is conservative.
#'
#' The Likelihood-Ratio test is computed as
#' \deqn{LR = 2 (LL_{comparison} - LL_{base})}
#' and is chi-squared with degrees of freedom equal to the difference in the
#' number of parameters. AIC and BIC are those of \code{AIC()} and
#' \code{BIC()}; for panel models, BIC uses the number of panels.
#'
#' @examples
#' 
#' # Example 1: Compare one fitted model against an estimated Poisson base model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT > 10000, 1, 0)
#'
#' nbp.base <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                      ShouldWidth04 + AADTover10k,
#'                      data = washington_roads, family = "NBP", method = "NM",
#'                      max.iters = 3000)
#'
#' regCompTest(nbp.base, data = washington_roads, basemodel = "Poisson",
#'             variables = TRUE, print = TRUE)
#'
#' # Example 2: Compare two fitted models directly
#' nb2.base <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                      ShouldWidth04 + AADTover10k,
#'                      data = washington_roads, family = "NB2", method = "NM",
#'                      max.iters = 3000)
#'
#' regCompTest(nbp.base, data = washington_roads, model2 = nb2.base,
#'             print = TRUE)
#'
#' @export
regCompTest <- function(
    model, data = NULL, basemodel = "Poisson",
    variables = FALSE, print = FALSE, ..., model2 = NULL) {
  if (inherits(model, "countreg_mix") || inherits(model2, "countreg_mix")) {
    stop("regCompTest() does not test finite-mixture or latent-class ",
         "models: the chi-square reference distribution of the ",
         "likelihood-ratio statistic does not hold when the number of ",
         "components changes or a component can vanish (a boundary of the ",
         "parameter space). Compare such models with AIC/BIC on the same ",
         "sample (regCompTable()), start stability and held-out log ",
         "scores.", call. = FALSE)
  }
  
  safe_round <- function(x, digits = 4) {
    if (length(x) == 0L || all(is.na(x))) {
      return(NA_real_)
    }
    round(x, digits)
  }
  
  make_single_model_stats <- function(model_res, base_res, nest) {
    test <- .fcr_lr_pvalue(model_res$LL, base_res$LL, nest)
    reason <- test$reason
    LR <- test$LR
    LRdof <- if (isTRUE(nest$nested)) nest$df else NA_integer_
    LR_pvalue <- test$p
    
    statistics <- tibble::tibble(
      Statistic = c(
        "Log-likelihood", "Number of parameters", "AIC", "BIC",
        "LR Test Statistic", "LR degrees of freedom",
        "LR p-value", "McFadden's Pseudo R^2"
      ),
      Model = c(
        safe_round(model_res$LL),
        model_res$n.coef,
        safe_round(model_res$AIC),
        safe_round(model_res$BIC),
        safe_round(LR),
        LRdof,
        safe_round(LR_pvalue),
        safe_round(model_res$PseudoR2)
      ),
      BaseModel = c(
        safe_round(base_res$LL),
        base_res$n.coef,
        safe_round(base_res$AIC),
        safe_round(base_res$BIC),
        NA,
        NA,
        NA,
        safe_round(base_res$PseudoR2)
      )
    )
    
    list(
      LR = LR,
      LRdof = LRdof,
      LR_pvalue = LR_pvalue,
      reason = reason,
      statistics = statistics
    )
  }
  
  make_pairwise_stats <- function(res1, res2, base_index, nest) {
    comparison_index <- if (base_index == 1L) 2L else 1L
    base_res <- if (base_index == 1L) res1 else res2
    comparison_res <- if (comparison_index == 1L) res1 else res2
    
    # comparison = full model, base = restricted model
    test <- .fcr_lr_pvalue(comparison_res$LL, base_res$LL, nest)
    LR <- test$LR
    LRdof <- if (isTRUE(nest$nested)) nest$df else NA_integer_
    LR_pvalue <- test$p
    
    statistics <- tibble::tibble(
      Statistic = c(
        "Log-likelihood", "Number of parameters", "AIC", "BIC",
        "LR Test Statistic", "LR degrees of freedom", "LR p-value",
        "McFadden's Pseudo R^2"
      ),
      Model = c(
        safe_round(comparison_res$LL),
        round(comparison_res$n.coef, 0),
        safe_round(comparison_res$AIC),
        safe_round(comparison_res$BIC),
        safe_round(LR),
        LRdof,
        safe_round(LR_pvalue),
        safe_round(comparison_res$PseudoR2)
      ),
      BaseModel = c(
        safe_round(base_res$LL),
        round(base_res$n.coef, 0),
        safe_round(base_res$AIC),
        safe_round(base_res$BIC),
        NA,
        NA,
        NA,
        safe_round(base_res$PseudoR2)
      )
    )
    
    list(
      base_index = base_index,
      comparison_index = comparison_index,
      base_model = base_res,
      comparison_model = comparison_res,
      LL_base = base_res$LL,
      LL_comparison = comparison_res$LL,
      LR = LR,
      LRdof = LRdof,
      LR_pvalue = LR_pvalue,
      statistics = statistics
    )
  }
  
  if (!is.null(model2)) {
    res_a <- get_model_components(model, data, basemodel, variables, ...)
    res_b <- get_model_components(model2, data, basemodel, variables, ...)
    if (res_a$nobs != res_b$nobs) {
      stop("The two models were fit to different numbers of observations (",
           res_a$nobs, " and ", res_b$nobs, "); an LR test is not valid.",
           call. = FALSE)
    }
    .fcr_check_same_sample(res_a$pieces, res_b$pieces)
    sa <- model$model$lr_sample
    sb <- model2$model$lr_sample
    if (!is.null(sa) && !is.null(sb)) {
      why <- .fcr_same_sample(sa, sb)
      if (nzchar(why)) {
        stop("The two models were not fit to the same sample (", why,
             "); an LR test is not valid.", call. = FALSE)
      }
    }
    # The restricted (base) model is the one with fewer parameters
    if (res_a$n.coef <= res_b$n.coef) {
      res1 <- res_a; res2 <- res_b; restricted_is_first <- TRUE
    } else {
      res1 <- res_b; res2 <- res_a; restricted_is_first <- FALSE
    }
  } else {
    res1 <- get_model_components(model, data, basemodel, variables, ...)
  }
  
  if (is.null(model2)) {
    # Base model on the same rows, offset and weights as the fitted model
    pc <- res1$pieces
    if (identical(toupper(basemodel), "POISSON")) {
      base_fit <- .fcr_poisson_ll(pc, if (variables) pc$X else NULL)
    } else {
      base_fit <- .fcr_countreg_base(pc, model$formula, basemodel,
                                     intercept_only = !variables, ...)
    }
    LLbase <- base_fit$LL
    n.coef.base <- base_fit$n.coef
    mean_cols <- colnames(pc$X)
    base_spec <- .fcr_base_spec(
      basemodel, if (variables) mean_cols else "(Intercept)")
    nest <- .fcr_nesting(base_spec, .fcr_model_spec(model))
    # The base model's mean (intercept only, or the model's design) must be
    # in the column space of the model's mean design
    if (nest$nested && !variables &&
        !.fcr_in_span(matrix(1, nrow(pc$X), 1), pc$X)) {
      nest <- list(nested = FALSE, reason = paste(
        "the model's mean design does not contain an intercept, so the",
        "intercept-only base model is not nested in it"))
    }
    if (nest$nested && identical(model$model$weight_type, "sampling")) {
      nest <- list(nested = FALSE, reason = paste(
        "the model uses sampling weights, so its log-likelihood is a",
        "pseudo-likelihood"))
    }
    if (nest$nested && !isTRUE(model$model$converged)) {
      nest <- list(nested = FALSE, reason = paste(
        "the model did not report convergence (or its convergence is",
        "unknown)"))
    }
    if (nest$nested && !isTRUE(base_fit$converged %||% TRUE)) {
      nest <- list(nested = FALSE,
                   reason = "the base model did not converge")
    }
    if (!nest$nested && res1$n.coef > n.coef.base) {
      message("No LR test: the ", toupper(basemodel), " base model is not ",
              "nested in the fitted model with a regular restriction (",
              nest$reason, "). Compare the models with AIC and BIC.")
    }
    AICbase <- -2 * LLbase + 2 * n.coef.base
    BICbase <- -2 * LLbase + log(res1$n_bic) * n.coef.base
    
    base_res <- list(
      LL = LLbase,
      n.coef = n.coef.base,
      AIC = AICbase,
      BIC = BICbase
    )
    
    summary_res <- make_single_model_stats(res1, base_res, nest)
    if (nest$nested && is.na(summary_res$LR)) {
      warning("No LR test: ", summary_res$reason, ".", call. = FALSE)
    }
    
    gtTable <- gt::gt(summary_res$statistics) |>
      tab_header(title = "Model Comparison Statistics") |>
      fmt_number(columns = c("Model", "BaseModel"), decimals = 4)
    
    test <- list(
      model = res1$fit,
      y = res1$y,
      LL = res1$LL,
      LLbase = LLbase,
      LR = summary_res$LR,
      LRdof = summary_res$LRdof,
      AIC = res1$AIC,
      AICbase = AICbase,
      BIC = res1$BIC,
      BICbase = BICbase,
      LR_pvalue = summary_res$LR_pvalue,
      LR_type = if (nest$nested && !is.na(summary_res$LR)) nest$type else
        NA_character_,
      LR_reason = if (nest$nested && !is.na(summary_res$LR)) NA_character_
        else nest$reason %||% summary_res$reason %||% NA_character_,
      PseudoR2 = res1$PseudoR2,
      statistics = summary_res$statistics,
      gtTable = gtTable,
      latexTable = knitr::kable(
        summary_res$statistics,
        format = "latex",
        booktabs = TRUE,
        caption = "Model Comparison Statistics"
      ),
      htmlTable = knitr::kable(
        summary_res$statistics,
        format = "html",
        table.attr = "class='table table-striped'",
        caption = "Model Comparison Statistics"
      )
    )
  } else {
    nest <- .fcr_lr_eligibility(res1$object, res2$object)
    if (res1$n.coef == res2$n.coef) {
      warning("The models have the same number of parameters, so they are ",
              "not nested; the LR test is not reported.", call. = FALSE)
    } else if (!nest$nested) {
      warning("The LR test is not reported: ", nest$reason, ". Compare the ",
              "models with AIC and BIC.", call. = FALSE)
    } else {
      chk <- .fcr_lr_pvalue(res2$LL, res1$LL, nest)
      if (is.na(chk$LR)) {
        warning("The LR test is not reported: ", chk$reason, " (check its ",
                "convergence).", call. = FALSE)
        nest <- list(nested = FALSE, reason = chk$reason)
      }
    }
    pair_res <- make_pairwise_stats(res2, res1, base_index = 2L, nest)
    pair_res$base_index <- if (restricted_is_first) 1L else 2L
    pair_res$comparison_index <- 3L - pair_res$base_index
    
    gtTable <- gt::gt(pair_res$statistics) |>
      tab_header(title = "Model Comparison Statistics") |>
      fmt_number(columns = c("Model", "BaseModel"), decimals = 4)
    
    test <- list(
      base_model = res1$fit,
      comparison_model = res2$fit,
      base_index = pair_res$base_index,
      comparison_index = pair_res$comparison_index,
      LL_base = pair_res$LL_base,
      LL_comparison = pair_res$LL_comparison,
      LR = pair_res$LR,
      LRdof = pair_res$LRdof,
      LR_pvalue = pair_res$LR_pvalue,
      LR_type = if (isTRUE(nest$nested)) nest$type else NA_character_,
      LR_reason = if (isTRUE(nest$nested)) NA_character_ else
        nest$reason %||% NA_character_,
      PseudoR2_base = pair_res$base_model$PseudoR2,
      PseudoR2_comparison = pair_res$comparison_model$PseudoR2,
      statistics = pair_res$statistics,
      gtTable = gtTable,
      latexTable = knitr::kable(
        pair_res$statistics,
        format = "latex",
        booktabs = TRUE,
        caption = "Model Comparison Statistics"
      ),
      htmlTable = knitr::kable(
        pair_res$statistics,
        format = "html",
        table.attr = "class='table table-striped'",
        caption = "Model Comparison Statistics"
      )
    )
  }
  
  if (print) {
    base::print(test$gtTable)
  }
  
  return(test)
}

# Two fits compared by an LR test must use the same observations, outcome,
# offsets and weights.
.fcr_check_same_sample <- function(pa, pb) {
  same <- length(pa$y) == length(pb$y) && isTRUE(all.equal(pa$y, pb$y)) &&
    isTRUE(all.equal(pa$weights, pb$weights)) &&
    isTRUE(all.equal(pa$offset, pb$offset))
  if (!same) {
    stop("The two models were not fit to the same observations, outcome, ",
         "offsets and weights; an LR test is not valid.", call. = FALSE)
  }
  invisible(NULL)
}

