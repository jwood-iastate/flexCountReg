#' Compare Regression Models with Likelihood Ratio Test, AIC, and BIC
#'
#' This function compares a given regression model to a base model using the Likelihood Ratio (LR) test, Akaike Information Criterion (AIC), and Bayesian Information Criterion (BIC).
#'
#' @param model A fitted regression model object.
#' @param data A data frame containing the variables in the model.
#' @param basemodel A character string specifying the type of base model to compare against. Default is "Poisson". Other options include any model specified from the \code{\link{flexCountReg}} function.
#' @param variables Logical. If \code{TRUE}, the base model will include the same variables as the provided model. If \code{FALSE}, the base model will be an intercept-only model. Default is \code{FALSE}.
#' @param print Logical. If \code{TRUE}, a table of the results will be shown. If \code{FALSE}, the table of results will not be printed to the console.
#' @param ... Additional arguments to be passed to the base model fitting function - options are any argument from the \code{\link{flexCountReg}} function.
#' @return A list containing the following components:
#' \item{LL}{Log-likelihood of the provided model.}
#' \item{LLbase}{Log-likelihood of the base model.}
#' \item{LR}{Likelihood Ratio statistic.}
#' \item{LRdof}{Degrees of freedom for the Likelihood Ratio test.}
#' \item{AIC}{Akaike Information Criterion for the provided model.}
#' \item{AICbase}{Akaike Information Criterion for the base model.}
#' \item{BIC}{Bayesian Information Criterion for the provided model.}
#' \item{BICbase}{Bayesian Information Criterion for the base model.}
#' \item{LR_pvalue}{P-value for the Likelihood Ratio test.}
#' \item{PseudoR2}{McFadden's Pseudo R^2.}
#' \item{statistics}{A tibble format summary of the results.}
#' \item{gtTable}{A \link[gt]{gt} table object summarizing the results.}
#' \item{latexTable}{Latex code for a table summarizing the results.}
#' \item{htmlTable}{HTML table summarizing the results.}
#' 
#' @include metrics.R
#' @import tibble knitr
#' @importFrom dplyr %>%
#' @importFrom gt gt tab_header fmt_number
#' 
#' @details The function performs the following steps:
#' \enumerate{
#' \item Fits the base model, either a Poisson regression or another specified model.
#' \item Computes the log-likelihoods of both the provided model and the base model.
#' \item Calculates the AIC and BIC for both models.
#' \item Conducts a Likelihood Ratio test to compare the models (if the provided model has more parameters than the base model).
#' \item Computes McFadden's Pseudo R^2.
#' }
#' 
#' The Likelihood-Ratio test is computed as \deqn{LR = -2 (LL_{base \ model}-LL_{model})}. The test is chi-squared with degrees of freedom \deqn{dof=N_{model \ params}-N_{base \ mode \ params}}.
#' The AIC is calculated as \deqn{AIC = -2 \cdot LL + 2 \cdot nparam}, and the BIC is calculated as \deqn{BIC = -2 \cdot LL + nparam \cdot \log(n)}.
#' @examples
#' 
#' # Comparing the NBP model with the NB2 model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
#'
#' nbp.base <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                     ShouldWidth04 + AADTover10k,
#'                     data=washington_roads, form = 'nbp', method = 'NM',
#'                     max.iters=3000)
#' comptests <- regCompTest(nbp.base, washington_roads, basemodel="NB2", print=TRUE)
#' 
#' @export
regCompTest <- function(model, data, basemodel = "Poisson", variables = FALSE, print=FALSE, ...){
  
  test <- c()
  
  formula <- model$formula
  mod_df <- stats::model.frame(formula, data)
  y <- as.numeric(stats::model.response(mod_df))
  
  LL <- model$maximum # The log-likelihood of the original model
  
  if(variables){
    if(basemodel=="Poisson"){
      base_mod <- glm(formula, data, family = poisson(link = "log"))
    }
    else{
      base_mod <- flexCountReg(formula, data, dist=basemodel, rpar_formula=NULL, ...)
    }
  }
  else{
    if(basemodel=="Poisson"){
      base_mod <- glm(y ~ 1, data, family = poisson(link = "log"))
    }
    else{
      base_mod <- flexCountReg(y ~ 1, data, dist=basemodel, rpar_formula=NULL, ...)
    }
  }
  
  # Get Log-Likelihood for the Base Model
  if (basemodel=="Poisson"){
    LLbase <- sum(dpois(base_mod$y, base_mod$fitted.values, log=TRUE))
  }
  else{
    LLbase <- base_mod$maximum
  }
  
  # get number of coefficients in the base model
  if (basemodel=="Poisson"){
    n.coef.base <- length(base_mod$coefficients)
  }
  else{
    n.coef.base <- length(base_mod$estimate)
  }
  
  # number of coefficients in the model provided
  n.coef <- length(model$estimate)
  
  test$LL <- LL
  test$LLbase <- LLbase
  
  # Compute Likelihood Ratio Test if the provided model has more parameters than the comparison model
  if (n.coef>n.coef.base){
    test$LR <- -2 * (LLbase - LL) # LR Statistic
    test$LRdof <- n.coef - n.coef.base # LR Degrees of Freedom
    
    if (test$LR > 0) {
      test$LR_pvalue <- pchisq(test$LR, test$LRdof, lower.tail = FALSE)  # LR p-Value
    } else {
      test$LR_pvalue <- 1
    }
  }
  else{
    test$LR <- NULL 
    test$LRdof <- NULL
    test$LR_pvalue <- NULL
  }
  
  
  # Calculate AIC and BIC for the main model and the base model
  test$AIC <- myAIC(LL, n.coef)
  test$AICbase <- myAIC(LLbase, n.coef.base)
  
  test$BIC <- myBIC(LL, n.coef, length(y))
  test$BICbase <- myBIC(LLbase, n.coef.base, length(y))
  
  
  
  # Compute McFadden's Pseudo R^2, based on a Poisson intercept-only model
  test$PseudoR2 <- 1 - LL / LLbase
  
  # Generate a table of the results and values
  statistics <- tibble::tibble(
    Statistic = c("AIC", "BIC", "LR Test Statistic", "LR degrees of freedom", "LR p-value", "McFadden's Pseudo R^2"),
    Model = c(round(test$AIC, 4), 
              round(test$BIC, 4), 
              round(test$LR, 4), 
              test$LRdof, 
              ifelse(test$LR_pvalue < 0.0001, "<0.0001", round(test$LR_pvalue, 4)), 
              round(test$PseudoR2, 4)),
    BaseModel = c(round(test$AICbase, 4), 
                  round(test$BICbase, 4), 
                  NA, 
                  NA, 
                  NA, 
                  NA)
  )
  
  test$statistics <- statistics
  
  # gt Table
  gtTable <- gt::gt(statistics) %>%
    tab_header(
      title = "Model Comparison Statistics"
    ) %>%
    fmt_number(
      columns = c(rlang::sym("Model"), rlang::sym("BaseModel")),
      decimals = 4
    )
  
  test$gtTable <- gtTable
  
  if (print){
    print(gtTable)
  }
  
  # LaTeX table
  test$latexTable <- knitr::kable(statistics, format = "latex", booktabs = TRUE, caption = "Model Comparison Statistics")
  
  # HTML table
  test$htmlTable <- knitr::kable(statistics, format = "html", table.attr = "class='table table-striped'", caption = "Model Comparison Statistics")
  
  return(test)
}

# Utility functions for AIC and BIC
myAIC <- function(LL, nparam) {
  return(-2 * LL + 2 * nparam)
}

myBIC <- function(LL, nparam, n) {
  return(-2 * LL + nparam * log(n))
}
