#' Create a Table Comparing Regression Models with AIC, BIC, and McFadden's Pseudo-R-Squared
#'
#' This function creates tables comparing the flexCountReg package models supplied to the function. 
#'
#' @name regCompTable
#' @param models A named list of fitted flexCountReg model objects. This must include 2 or more models.
#' @param coefs A logical. The default value `TRUE` indicates that the coefficients from the models should be included in the table of comparisons.
#' @param AIC A logical. The default value `TRUE` indicates that AIC values for the models should be included.
#' @param BIC A logical. The default value `TRUE` indicates that BIC values for the models should be included.
#' @param RSquare A logical. The default value `TRUE` indicates that the McFadden's Pseudo-R-Squared statistic (comparing against a Poisson regression model) should be included.
#' @param tableType The type of table format to return. Options include "tibble" for returning the table as a tibble, "gt" for a \link[gt]{gt} table object, or "latex" for a latex table. The default is "tibble".
#' @param digits An integer value indicating the number of decimals to round the table values to.
#'
#' @include regCompTest.R
#' @import tibble knitr  
#' @importFrom dplyr %>% mutate sym
#' @importFrom stats as.formula model.frame model.response dpois glm
#' @importFrom gt gt tab_header fmt_number
#' 
#' @examples
#' 
#' # Comparing the NBP model with the NB2 and NB1 models
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
#'
#' nb.1 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                     ShouldWidth04 + AADTover10k,
#'                     data=washington_roads, family = 'NB1', method = 'NM')
#'                     
#' nb.2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                     ShouldWidth04 + AADTover10k,
#'                     data=washington_roads, family = 'NB2', method = 'NM')
#'                     
#'                     
#' nb.p <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                     ShouldWidth04 + AADTover10k,
#'                     data=washington_roads, family = 'NBP', method = 'NM')
#'                     
#'                     
#' comptable <- regCompTable(list("NB-1"=nb.1, "NB-2"=nb.2, "NB-P"=nb.p), tableType="latex")
#' print(comptable)
#' @export
regCompTable <- function(models, coefs=TRUE, AIC=TRUE, BIC=TRUE, RSquare=TRUE, 
                         tableType="tibble", digits=3){
  
  Nmodels <- length(models)
  if(Nmodels < 2) {
    msg <- paste0("Provide a list of 2 or more models for the `models` input ", 
                  "to call the `regCompTable` function")
    warning(msg)
  }
  
  # Extract the names of all coefficients in any of the models supplied and create an initial tibble with a column of the coefficient names and other stats
  vars <- c() # names of estimated coefficients
  modNames <- names(models)
  for (i in seq_along(models)){
    variables <- names(models[[i]]$model$estimate) # Get the variable names
    for (item in variables) { # Add any new, unique, variables/coefficient names
      if (!(item %in% vars)) {
        vars <- append(vars, item)
      }
    }
  }
  vars <- append(vars, "N Obs.") # Add the number of observed datapoints
  vars <- append(vars, "LL") # Add the Log Likelihood
  
  if(AIC){
    vars <- append(vars, "AIC")
  }
  if(BIC){
    vars <- append(vars, "BIC")
  }
  if(RSquare){
    vars <- append(vars, "Pseudo-R-Sq.")
  }
  
  summaryTibble <- tibble(`Parameter` = vars) # Initial tibble
  
  for (i in seq_along(modNames)){ 
    formula <- as.formula(models[[i]]$formula) # not including random parameters
    data <- models[[i]]$data
    mod_df <- model.frame(formula, data)
    y <- as.numeric(model.response(mod_df))
    LL <- as.numeric(models[[i]]$model$maximum) # The log-likelihood of the original model
    base_mod <- glm(y ~ 1, data, family = poisson(link = "log"))
    LLbase <- sum(dpois(base_mod$y, base_mod$fitted.values, log=TRUE))
    n.coef <- length(models[[i]]$model$estimate)
    aic <- round(myAIC(LL, n.coef), digits)
    bic <- round(myBIC(LL, n.coef, length(y)), digits)
    rsq <- round(1 - LL / LLbase, digits)
    Nobs <- length(y)
    
    coefNames <- names(models[[i]]$model$estimate)
    coefs <- models[[i]]$model$estimate
    if (!is.null(models[[i]]$model$bootstraps)){
      se <- models[[i]]$model$bootstrapped_se
    }else{
      se <- sqrt(diag(-1/(models[[i]]$model$hessian)))
    }
    t <- abs(coefs / se)
    coefsRounded <- round(coefs, digits)
    seRounded <- round(se, digits)
    tablecoefs <- c()
    
    for (j in seq_along(vars)){
      if (vars[j] %in% coefNames){
        for (k in seq_along(coefNames)){
          if (coefNames[k] == vars[j]){
            if(t[k] >= 3.29){
              coefval <- paste0(coefsRounded[k], " (", seRounded[k], ")***")
            }
            else if (t[k] >= 2.5758) {
              coefval <- paste0(coefsRounded[k], " (", seRounded[k], ")**")
            }
            else if (t[k] >= 1.96){
              coefval <- paste0(coefsRounded[k], " (", seRounded[k], ")*")
            }
            else{
              coefval <- paste0(coefsRounded[k], " (", seRounded[k], ")")
            }
            tablecoefs <- append(tablecoefs, coefval)
          }
        }
      } else if (vars[j] == "N Obs.") {
        tablecoefs <- append(tablecoefs, Nobs)
      } else if (vars[j]  == "LL") {
        tablecoefs <- append(tablecoefs, round(as.numeric(LL), digits))
      } else if (vars[j]  == "AIC") {
        tablecoefs <- append(tablecoefs, aic)
      } else if (vars[j]  == "BIC") {
        tablecoefs <- append(tablecoefs, bic)
      } else if (vars[j]  == "Pseudo-R-Sq.") {
        tablecoefs <- append(tablecoefs, rsq)
      } else{
        tablecoefs <- append(tablecoefs, "---")
      }
    }
    summaryTibble <- summaryTibble %>%
      mutate(!!sym(modNames[i]) := tablecoefs)
  }
  
  # Adding a row with the notation
  # notation_row <- tibble(
  #   Parameter = "Note",
  #   !!!setNames(rep("p-value codes: * (p<=0.05), ** (p<=0.01), *** (p<=0.001)", length(modNames)), modNames)
  # )
  
  
  if(tableType == "tibble"){
    return(summaryTibble)
  } 
  else if (tableType == "gt"){
    gtTable <- gt(summaryTibble) %>%
      tab_header(
        title = "Model Comparisons",
        subtitle = "Mean (Standard Error)"
      ) %>%
      gt::tab_source_note(
        source_note = "p-value codes: *=(p<=0.05), **=(p<=0.01), ***=(p<=0.001)"
      )
    return(gtTable)
  }
  else if (tableType == "latex"){
    latexTable <- knitr::kable(summaryTibble, format = "latex", booktabs = TRUE, caption = "Model Comparison Statistics")
    return(latexTable)
  }
}
