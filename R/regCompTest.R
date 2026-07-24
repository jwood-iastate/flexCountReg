.mcfadden_r2 <- function(ll_full, ll_null) {
  if (is.finite(ll_full) && is.finite(ll_null) && ll_null != 0) {
    1 - (ll_full / ll_null)
  } else {
    NA_real_
  }
}

.fit_null_ll <- function(trained_model, data) {
  fit <- trained_model$model
  if (is.null(fit)) {
    return(NA_real_)
  }
  
  # Use the model's own formula/data to build the response
  formula <- fit$formula
  mf <- stats::model.frame(formula, data)
  y <- as.numeric(stats::model.response(mf))
  n <- length(y)
  
  # Offset if present
  off <- stats::model.offset(mf)
  if (is.null(off)) off <- rep(0, n)
  
  # Case 1: Poisson GLM-style model
  if (!is.null(fit$family)) {
    # Intercept-only Poisson with offset:
    # mu_i = exp(beta0 + offset_i)
    beta0 <- log(sum(y) / sum(exp(off)))
    mu0 <- exp(beta0 + off)
    return(sum(stats::dpois(y, mu0, log = TRUE)))
  }
  
  # Case 2: your custom countreg-style models
  # Replace this with the null-model fitting call your package uses
  # (e.g., countreg(y ~ 1, data = ..., family = ...))
  if (!is.null(fit$family)) {
    null_fit <- try(
      countreg(
        stats::as.formula(paste(deparse(formula[[2L]]), "~ 1")),
        data = data,
        family = fit$family,
        method = fit$method
      ),
      silent = TRUE
    )
    
    if (!inherits(null_fit, "try-error") && !is.null(null_fit$model$maximum)) {
      return(null_fit$model$maximum)
    }
  }
  
  NA_real_
}

get_model_components <- function(trained_model, data) {
  fit <- trained_model$model
  if (is.null(fit)) {
    stop("The supplied object does not contain a fitted model in `$model`.")
  }
  
  formula <- fit$formula
  mod_df <- stats::model.frame(formula, data)
  y <- as.numeric(stats::model.response(mod_df))
  
  LL <- fit$maximum
  n.coef <- length(fit$estimate)
  nobs <- length(y)
  
  LL_null <- .fit_null_ll(trained_model, data)
  PseudoR2 <- .mcfadden_r2(LL, LL_null)
  
  list(
    fit = fit,
    y = y,
    LL = LL,
    LL_null = LL_null,
    PseudoR2 = PseudoR2,
    n.coef = n.coef,
    nobs = nobs,
    AIC = myAIC(LL, n.coef),
    BIC = myBIC(LL, n.coef, nobs)
  )
}

#' Compare Regression Models with Likelihood Ratio Test, AIC, and BIC
#'
#' This function compares fitted regression model objects using the
#' Likelihood Ratio (LR) test, Akaike Information Criterion (AIC), and Bayesian
#' Information Criterion (BIC).
#'
#' When only one fitted model is supplied, the function estimates a base model
#' from the supplied data and compares the fitted model against that base
#' model. When two fitted models are supplied, the function compares the two
#' fitted models directly and uses the model with the larger log-likelihood as
#' the base model for the LR test.
#'
#' @name regCompTest
#' @param model A fitted regression model object.
#' @param data An optional data frame containing the variables in the model.
#'   If not supplied, the original data used to estimate the model will be
#'   used when available.
#' @param basemodel A character string specifying the family of base model to
#'   compare against (options include the family from
#'   \code{\link{countreg}} or "Poisson"). Default is "Poisson".
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
#' @import tibble knitr
#' @importFrom dplyr %>%
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
#' \item Computes McFadden's Pseudo R^2.
#' }
#'
#' For the two-model workflow, the function compares the two fitted models
#' directly. The model with the larger log-likelihood is treated as the base
#' model for the LR test, and the other model is treated as the comparison
#' model.
#'
#' The Likelihood-Ratio test is computed as \deqn{LR = -2 (LL_{base} - LL_{comparison})}.
#' The test is chi-squared with degrees of freedom
#' \deqn{dof = |N_{base \ params} - N_{comparison \ params}|}.
#' The AIC is calculated as \deqn{AIC = -2 \cdot LL + 2 \cdot nparam}, and the
#' BIC is calculated as \deqn{BIC = -2 \cdot LL + nparam \cdot \log(n)}.
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
  
  safe_round <- function(x, digits = 4) {
    if (length(x) == 0L || all(is.na(x))) {
      return(NA_real_)
    }
    round(x, digits)
  }
  
  first_non_null <- function(...) {
    vals <- list(...)
    for (v in vals) {
      if (!is.null(v)) {
        return(v)
      }
    }
    NULL
  }
  
  make_single_model_stats <- function(model_res, base_res) {
    LR <- NA_real_
    LRdof <- NA_integer_
    LR_pvalue <- NA_real_
    
    if (model_res$n.coef > base_res$n.coef) {
      LR <- -2 * (base_res$LL - model_res$LL)
      LRdof <- model_res$n.coef - base_res$n.coef
      if (is.finite(LR) && is.finite(LRdof) && LR > 0 && LRdof > 0) {
        LR_pvalue <- pchisq(LR, LRdof, lower.tail = FALSE)
      } else {
        LR_pvalue <- 1
      }
    }
    
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
      statistics = statistics
    )
  }
  
  make_pairwise_stats <- function(res1, res2, base_index) {
    comparison_index <- if (base_index == 1L) 2L else 1L
    base_res <- if (base_index == 1L) res1 else res2
    comparison_res <- if (comparison_index == 1L) res1 else res2
    
    LR <- -2 * (comparison_res$LL - base_res$LL)
    LRdof <- abs(base_res$n.coef - comparison_res$n.coef)
    LR_pvalue <- NA_real_
    if (is.finite(LR) && is.finite(LRdof) && LR > 0 && LRdof > 0) {
      LR_pvalue <- pchisq(LR, LRdof, lower.tail = FALSE)
    } else if (LR == 0) {
      LR_pvalue <- 1
    }
    
    statistics <- tibble::tibble(
      Statistic = c(
        "Log-likelihood", "Number of parameters", "AIC", "BIC",
        "LR Test Statistic", "LR degrees of freedom", "LR p-value",
        "McFadden's Pseudo R^2"
      ),
      Model = c(
        safe_round(base_res$LL),
        round(base_res$n.coef, 0),
        safe_round(base_res$AIC),
        safe_round(base_res$BIC),
        safe_round(LR),
        LRdof,
        safe_round(LR_pvalue),
        safe_round(base_res$PseudoR2)
      ),
      BaseModel = c(
        safe_round(comparison_res$LL),
        round(comparison_res$n.coef, 0),
        safe_round(comparison_res$AIC),
        safe_round(comparison_res$BIC),
        NA,
        NA,
        NA,
        safe_round(comparison_res$PseudoR2)
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
  
  if (is.null(data)) {
    data <- first_non_null(model$data, if (!is.null(model2)) model2$data else NULL)
  }
  if (is.null(data)) {
    stop("No data were supplied and the model object does not contain `$data`.")
  }
  
  if(!is.null(model2)){
    m1LL = model$model$logLik
    m2LL = model2$model$logLik
    
    if (m1LL<m2LL){
      res1 = get_model_components(model2, data)
      res2 = get_model_components(model, data)
    }else{
      res1 = get_model_components(model, data)
      res2 = get_model_components(model2, data)
    }
  }else{
    res1 <- get_model_components(model, data)
  }
  
  if (is.null(model2)) {
    if (variables) {
      if (basemodel == "Poisson") {
        base_mod <- glm(res1$fit$formula, data, family = poisson(link = "log"))
      } else {
        base_mod <- countreg(res1$fit$formula, data, family = basemodel, ...)
        base_mod <- base_mod$model
      }
    } else {
      if (basemodel == "Poisson") {
        base_mod <- glm(res1$y ~ 1, data, family = poisson(link = "log"))
      } else {
        base_mod <- countreg(res1$y ~ 1, data, family = basemodel, ...)
        base_mod <- base_mod$model
      }
    }
    
    if (basemodel == "Poisson") {
      LLbase <- sum(dpois(base_mod$y, base_mod$fitted.values, log = TRUE))
      n.coef.base <- length(base_mod$coefficients)
      AICbase <- myAIC(LLbase, n.coef.base)
      BICbase <- myBIC(LLbase, n.coef.base, res1$nobs)
    } else {
      LLbase <- base_mod$maximum
      n.coef.base <- length(base_mod$estimate)
      AICbase <- myAIC(LLbase, n.coef.base)
      BICbase <- myBIC(LLbase, n.coef.base, res1$nobs)
    }
    
    base_res <- list(
      LL = LLbase,
      n.coef = n.coef.base,
      AIC = AICbase,
      BIC = BICbase
    )
    
    summary_res <- make_single_model_stats(res1, base_res)
    
    gtTable <- gt::gt(summary_res$statistics) %>%
      tab_header(title = "Model Comparison Statistics") %>%
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
    base_index <- if (res1$LL >= res2$LL) 1L else 2L
    pair_res <- make_pairwise_stats(res1, res2, base_index)
    
    gtTable <- gt::gt(pair_res$statistics) %>%
      tab_header(title = "Model Comparison Statistics") %>%
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

# Utility functions for AIC and BIC
myAIC <- function(LL, nparam) {
  return(-2 * LL + 2 * nparam)
}

myBIC <- function(LL, nparam, n) {
  return(-2 * LL + nparam * log(n))
}

