#' Compare Regression Models with Likelihood Ratio Test, AIC, and BIC
#'
#' This function compares fitted regression model objects using the
#' Likelihood Ratio (LR) test, Akaike Information Criterion (AIC), and Bayesian
#' Information Criterion (BIC).
#'
#' @name regCompTest
#' @param model A fitted regression model object.
#' @param data An options data frame containing the variables in the model. If
#'   not supplied, the original data used to estimate the model will be used.
#' @param basemodel A character string specifying the family of base model to
#'   compare against (options include the family from \code{\link{countreg}} or
#'   "Poisson"). Default is "Poisson".
#' @param variables Logical. If \code{TRUE}, the base model will include the
#'   same variables as the provided model. If \code{FALSE}, the base model will
#'   be an intercept-only model. Default is \code{FALSE}.
#' @param print Logical. If \code{TRUE}, a table of the results will be shown.
#'   If \code{FALSE}, the table of results will not be printed to the console.
#' @param ... Additional arguments to be passed to the base model fitting
#'   function - options are any argument from the \code{\link{countreg}}
#'   function.
#' @param model2 An optional second fitted regression model object to compare.
#'   If supplied, the function returns comparison results for both trained
#'   models in addition to the base-model comparison.
#' @returns A list containing the following components.
#' \itemize{
#' \item For the original single-model workflow: \code{LL}, \code{LLbase},
#'   \code{LR}, \code{LRdof}, \code{AIC}, \code{AICbase}, \code{BIC},
#'   \code{BICbase}, \code{LR_pvalue}, \code{PseudoR2}, \code{statistics},
#'   \code{gtTable}, \code{latexTable}, and \code{htmlTable}.
#' \item If \code{model2} is supplied: \code{model1}, \code{model2}, and a
#'   combined \code{statistics} table plus the associated table objects.
#' }
#'
#' @include metrics.R
#' @import tibble knitr
#' @importFrom dplyr %>%
#' @importFrom gt gt tab_header fmt_number
#'
#' @details The function performs the following steps:
#' \enumerate{
#' \item Fits the base model, either a Poisson regression or another specified
#' model.
#' \item Computes the log-likelihoods of the provided model(s) and the base
#' model.
#' \item Calculates the AIC and BIC for the model(s).
#' \item Conducts a Likelihood Ratio test to compare the model(s) to the base
#' model (if the provided model has more parameters than the base model).
#' \item Computes McFadden's Pseudo R^2.
#' }
#'
#' The Likelihood-Ratio test is computed as \deqn{LR = -2 (LL_{base} - LL_{model})}.
#' The test is chi-squared with degrees of freedom
#' \deqn{dof=N_{model \ params}-N_{base \ params}}.
#' The AIC is calculated as \deqn{AIC = -2 \cdot LL + 2 \cdot nparam}, and the
#' BIC is calculated as \deqn{BIC = -2 \cdot LL + nparam \cdot \log(n)}.
#'
#' @examples
#'
#' # Comparing the NBP model with the NB2 model and base Poisson model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
#'
#' nbp.base <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                     ShouldWidth04 + AADTover10k,
#'                     data=washington_roads, family = 'NBP', method = 'NM',
#'                     max.iters=3000)
#' nb2.base <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                     ShouldWidth04 + AADTover10k,
#'                     data=washington_roads, family = 'NB2', method = 'NM',
#'                     max.iters=3000)                    
#'                    
#' regCompTest(nbp.base, washington_roads, basemodel="Poisson", 
#' model2=nb2.base, print=TRUE)
#'
#' @export
regCompTest <- function(
    model, data = NULL, basemodel = "Poisson",
    variables = FALSE, print = FALSE, ..., model2 = NULL){
  
  safe_round <- function(x, digits = 4) {
    if (length(x) == 0L || all(is.na(x))) {
      return(NA_real_)
    }
    round(x, digits)
  }
  
  get_model_components <- function(trained_model) {
    fit <- trained_model$model
    
    formula <- fit$formula
    mod_df <- stats::model.frame(formula, data)
    y <- as.numeric(stats::model.response(mod_df))
    
    LL <- fit$maximum
    
    if (variables) {
      if (basemodel == "Poisson") {
        base_mod <- glm(formula, data, family = poisson(link = "log"))
      } else {
        base_mod <- countreg(formula, data, family = basemodel, ...)
        base_mod <- base_mod$model
      }
    } else {
      if (basemodel == "Poisson") {
        base_mod <- glm(y ~ 1, data, family = poisson(link = "log"))
      } else {
        base_mod <- countreg(y ~ 1, data, family = basemodel, ...)
        base_mod <- base_mod$model
      }
    }
    
    if (basemodel == "Poisson") {
      LLbase <- sum(dpois(base_mod$y, base_mod$fitted.values, log = TRUE))
      n.coef.base <- length(base_mod$coefficients)
    } else {
      LLbase <- base_mod$maximum
      n.coef.base <- length(base_mod$estimate)
    }
    
    n.coef <- length(fit$estimate)
    
    LR <- NA_real_
    LRdof <- NA_integer_
    LR_pvalue <- NA_real_
    
    if (n.coef > n.coef.base) {
      LR <- -2 * (LLbase - LL)
      LRdof <- n.coef - n.coef.base
      
      if (is.finite(LR) && LR > 0) {
        LR_pvalue <- pchisq(LR, LRdof, lower.tail = FALSE)
      } else {
        LR_pvalue <- 1
      }
    }
    
    AIC <- myAIC(LL, n.coef)
    AICbase <- myAIC(LLbase, n.coef.base)
    
    BIC <- myBIC(LL, n.coef, length(y))
    BICbase <- myBIC(LLbase, n.coef.base, length(y))
    
    PseudoR2 <- 1 - LL / LLbase
    
    statistics <- tibble::tibble(
      Statistic = c(
        "AIC", "BIC", "LR Test Statistic", "LR degrees of freedom",
        "LR p-value", "McFadden's Pseudo R^2"
      ),
      Model = c(
        safe_round(AIC),
        safe_round(BIC),
        safe_round(LR),
        LRdof,
        safe_round(LR_pvalue),
        safe_round(PseudoR2)
      ),
      BaseModel = c(
        safe_round(AICbase),
        safe_round(BICbase),
        NA,
        NA,
        NA,
        NA
      )
    )
    
    gtTable <- gt::gt(statistics) %>%
      tab_header(
        title = "Model Comparison Statistics"
      ) %>%
      fmt_number(
        columns = c(rlang::sym("Model"), rlang::sym("BaseModel")),
        decimals = 4
      )
    
    latexTable <- knitr::kable(
      statistics,
      format = "latex",
      booktabs = TRUE,
      caption = "Model Comparison Statistics"
    )
    
    htmlTable <- knitr::kable(
      statistics,
      format = "html",
      table.attr = "class='table table-striped'",
      caption = "Model Comparison Statistics"
    )
    
    list(
      model = fit,
      y = y,
      LL = LL,
      LLbase = LLbase,
      LR = LR,
      LRdof = LRdof,
      AIC = AIC,
      AICbase = AICbase,
      BIC = BIC,
      BICbase = BICbase,
      LR_pvalue = LR_pvalue,
      PseudoR2 = PseudoR2,
      statistics = statistics,
      gtTable = gtTable,
      latexTable = latexTable,
      htmlTable = htmlTable
    )
  }
  
  if (is.null(data)) { # Use object data if no new data are provided
    data <- model$data
  }
  
  res1 <- get_model_components(model)
  
  if (is.null(model2)) {
    test <- res1
  } else {
    res2 <- get_model_components(model2)
    
    combined_stats <- tibble::tibble(
      Statistic = res1$statistics$Statistic,
      Model1 = res1$statistics$Model,
      BaseModel1 = res1$statistics$BaseModel,
      Model2 = res2$statistics$Model,
      BaseModel2 = res2$statistics$BaseModel
    )
    
    combined_gt <- gt::gt(combined_stats) %>%
      tab_header(
        title = "Model Comparison Statistics"
      ) %>%
      fmt_number(
        columns = c(
          rlang::sym("Model1"),
          rlang::sym("BaseModel1"),
          rlang::sym("Model2"),
          rlang::sym("BaseModel2")
        ),
        decimals = 4
      )
    
    test <- list(
      model1 = res1,
      model2 = res2,
      statistics = combined_stats,
      gtTable = combined_gt,
      latexTable = knitr::kable(
        combined_stats,
        format = "latex",
        booktabs = TRUE,
        caption = "Model Comparison Statistics"
      ),
      htmlTable = knitr::kable(
        combined_stats,
        format = "html",
        table.attr = "class='table table-striped'",
        caption = "Model Comparison Statistics"
      )
    )
  }
  
  if (print) {
    print(test$gtTable)
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